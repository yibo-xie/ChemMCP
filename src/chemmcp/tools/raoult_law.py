"""
拉乌尔定律计算工具
计算理想溶液中各组分的蒸气压、总蒸气压和气相组成。
"""
import logging
from typing import Dict, Any, List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RaoultLaw(BaseTool):
    """
    拉乌尔定律计算工具。

    对于理想溶液，P_i = x_i × P_i*，并使用 Dalton 分压定律计算气相组成。
    """
    __version__                 = "0.1.0"
    name                        = "RaoultLaw"
    func_name                   = "calculate_raoult_law"
    description                 = "Calculate vapor pressures of ideal solutions using Raoult's law and vapor-phase compositions using Dalton's law of partial pressures."
    implementation_description  = "Applies Raoult's law (P_i = x_i * P_i*) for each component to get partial pressures, sums them for total pressure, then uses Dalton's law (y_i = P_i / P_total) for vapor composition."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Raoult's Law", "Vapor Pressure", "Solution Chemistry", "Thermodynamics", "Ideal Solution"]
    required_envs               = []

    code_input_sig = [
        ("mole_fractions",        "list", "N/A",     "List of mole fractions of each component in the liquid phase (must sum to ~1)."),
        ("pure_vapor_pressures",  "list", "N/A",     "List of pure-component vapor pressures in same units (e.g., mmHg, kPa, atm). Must match length of mole_fractions."),
        ("component_names",       "list", '["A","B"]', "Optional list of component names."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Semicolon-separated: 'x1,x2,...;P1*,P2*,...;name1,name2,...'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing partial_pressures, total_vapor_pressure, vapor_phase_compositions, and validation info."),
    ]

    examples = [
        {
            "code_input": {
                "mole_fractions": [0.4, 0.6],
                "pure_vapor_pressures": [119.0, 37.0],
                "component_names": ["Benzene", "Toluene"],
            },
            "text_input": {
                "input_params": "0.4,0.6;119.0,37.0;Benzene,Toluene",
            },
            "output": {
                "result": {
                    "partial_pressures": {"Benzene": 47.6, "Toluene": 22.2},
                    "total_vapor_pressure": 69.8,
                    "vapor_phase_compositions": {"Benzene": 0.682, "Toluene": 0.318},
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        mole_fractions: List[float],
        pure_vapor_pressures: List[float],
        component_names: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """核心逻辑：应用拉乌尔定律和道尔顿分压定律。"""
        # ---- 输入验证 ----
        if len(mole_fractions) != len(pure_vapor_pressures):
            raise ChemMCPError("mole_fractions and pure_vapor_pressures must have the same length.")
        if len(mole_fractions) == 0:
            raise ChemMCPError("Input lists cannot be empty.")

        n = len(mole_fractions)
        total_x = sum(mole_fractions)
        if abs(total_x - 1.0) > 0.01:
            logger.warning(f"Mole fractions sum to {total_x}, normalizing.")
            mole_fractions = [x / total_x for x in mole_fractions]

        for i, x in enumerate(mole_fractions):
            if x < 0:
                raise ChemMCPError(f"Mole fraction of component {i} is negative.")

        for i, p in enumerate(pure_vapor_pressures):
            if p < 0:
                raise ChemMCPError(f"Vapor pressure of component {i} cannot be negative.")

        if component_names is None:
            component_names = [f"Component_{i+1}" for i in range(n)]
        elif len(component_names) != n:
            component_names = [f"Component_{i+1}" for i in range(n)]

        # ---- Raoult 定律：分压 P_i = x_i * P_i* ----
        partial_pressures = {}
        for i in range(n):
            p_partial = mole_fractions[i] * pure_vapor_pressures[i]
            partial_pressures[component_names[i]] = round(p_partial, 6)

        # ---- 总蒸气压 ----
        total_pressure = sum(partial_pressures.values())

        # ---- Dalton 定律：气相组成 y_i = P_i / P_total ----
        vapor_compositions = {}
        if total_pressure > 0:
            for name, p_partial in partial_pressures.items():
                vapor_compositions[name] = round(p_partial / total_pressure, 6)
        else:
            for name in component_names:
                vapor_compositions[name] = 0.0

        logger.info(f"Raoult's law: P_total = {total_pressure:.4f}, "
                     f"vapor compo = {vapor_compositions}")

        return {
            "partial_pressures": partial_pressures,
            "total_vapor_pressure": round(total_pressure, 6),
            "vapor_phase_compositions": vapor_compositions,
            "liquid_mole_fractions": {component_names[i]: round(mole_fractions[i], 6) for i in range(n)},
            "number_of_components": n,
            "component_names": component_names,
            "explanation": (
                f"Raoult's law applied to {n}-component ideal solution. "
                f"Total vapor pressure = {total_pressure:.4f}. "
                f"Vapor phase is enriched in the more volatile component(s)."
            ),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split(";")
            mole_fractions = [float(x.strip()) for x in parts[0].split(",")]
            pure_vapor_pressures = [float(x.strip()) for x in parts[1].split(",")]
            component_names = None
            if len(parts) > 2:
                component_names = [s.strip() for s in parts[2].split(",")]

            return self._run_base(mole_fractions, pure_vapor_pressures, component_names)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'x1,x2,...;P1*,P2*,...;name1,name2,...'"
            )
