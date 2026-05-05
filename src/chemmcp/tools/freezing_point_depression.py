"""
凝固点降低计算工具（依数性性质）
计算非挥发性溶质引起的溶液凝固点降低：ΔTf = i × Kf × m
"""
import logging
from typing import Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FreezingPointDepression(BaseTool):
    """
    凝固点降低计算工具。

    基于依数性公式 ΔTf = i × Kf × m 计算溶液的凝固点降低值和新凝固点。
    """
    __version__                 = "0.1.0"
    name                        = "FreezingPointDepression"
    func_name                   = "calculate_freezing_point_depression"
    description                 = "Calculate freezing point depression of a solution using the colligative property formula ΔTf = i × Kf × m."
    implementation_description  = "Uses the freezing point depression formula: ΔTf = i·Kf·m, where i is the van't Hoff factor, Kf is the cryoscopic constant, and m is molality."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Colligative Properties", "Freezing Point", "Solution Chemistry", "Physical Chemistry"]
    required_envs               = []

    code_input_sig = [
        ("molality",              "float", "N/A",     "Molality of solution in mol solute / kg solvent."),
        ("cryoscopic_constant_kf", "float", "1.86",   "Cryoscopic constant (Kf) in K·kg/mol (water = 1.86)."),
        ("vanthoff_factor_i",     "float", "1.0",     "van't Hoff factor (i=1 for nonelectrolyte, i≈2 for NaCl, etc.)."),
        ("normal_freezing_point_k","float","273.15",   "Normal freezing point of pure solvent in Kelvin."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'molality [Kf] [i] [fp_K]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing delta_tf_k, new_freezing_point_k, and explanation."),
    ]

    examples = [
        {
            "code_input": {
                "molality": 0.5,
                "cryoscopic_constant_kf": 1.86,
                "vanthoff_factor_i": 1.0,
                "normal_freezing_point_k": 273.15,
            },
            "text_input": {
                "input_params": "0.5 1.86 1.0 273.15",
            },
            "output": {
                "result": {
                    "delta_tf_k": 0.93,
                    "new_freezing_point_k": 272.22,
                }
            },
        },
        # NaCl in water example
        {
            "code_input": {
                "molality": 0.1,
                "cryoscopic_constant_kf": 1.86,
                "vanthoff_factor_i": 1.9,
                "normal_freezing_point_k": 273.15,
            },
            "text_input": {
                "input_params": "0.1 1.86 1.9 273.15",
            },
            "output": {
                "result": {
                    "delta_tf_k": 0.3534,
                    "new_freezing_point_k": 272.7966,
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
        molality: float,
        cryoscopic_constant_kf: float = 1.86,
        vanthoff_factor_i: float = 1.0,
        normal_freezing_point_k: float = 273.15,
    ) -> Dict[str, Any]:
        """核心逻辑：ΔTf = i × Kf × m"""
        if molality < 0:
            raise ChemMCPError("Molality cannot be negative.")
        if cryoscopic_constant_kf < 0:
            raise ChemMCPError("Cryoscopic constant cannot be negative.")
        if vanthoff_factor_i < 0:
            raise ChemMCPError("van't Hoff factor cannot be negative.")
        if normal_freezing_point_k <= 0:
            raise ChemMCPError("Normal freezing point must be positive.")

        delta_tf = vanthoff_factor_i * cryoscopic_constant_kf * molality
        new_fp = normal_freezing_point_k - delta_tf

        logger.info(f"Freezing point depression: ΔTf={delta_tf:.4f}K, new T_f={new_fp:.4f}K "
                     f"(m={molality}, Kf={cryoscopic_constant_kf}, i={vanthoff_factor_i})")

        return {
            "delta_tf_k": round(delta_tf, 6),
            "new_freezing_point_k": round(new_fp, 6),
            "new_freezing_point_c": round(new_fp - 273.15, 4),
            "parameters_used": {
                "molality_mol_kg": molality,
                "cryoscopic_constant_Kf": cryoscopic_constant_kf,
                "vanthoff_factor_i": vanthoff_factor_i,
                "normal_freezing_point_K": normal_freezing_point_k,
            },
            "explanation": (
                f"ΔTf = i × Kf × m = {vanthoff_factor_i} × {cryoscopic_constant_kf} × {molality} "
                f"= {delta_tf:.4f} K. New freezing point = {new_fp:.4f} K ({new_fp - 273.15:.2f} °C)."
            ),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            kwargs = {"molality": float(parts[0])}
            if len(parts) > 1: kwargs["cryoscopic_constant_kf"] = float(parts[1])
            if len(parts) > 2: kwargs["vanthoff_factor_i"] = float(parts[2])
            if len(parts) > 3: kwargs["normal_freezing_point_k"] = float(parts[3])
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'molality [Kf] [i] [fp_K]'")
