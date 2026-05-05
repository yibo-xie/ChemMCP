"""
渗透压计算工具（依数性性质）
使用 van't Hoff 公式计算溶液渗透压：π = i × c × R × T
"""
import logging
from typing import Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class OsmoticPressure(BaseTool):
    """
    渗透压计算工具。

    使用 van't Hoff 渗透压公式 π = i × c × R × T 计算渗透压。
    """
    __version__                 = "0.1.0"
    name                        = "OsmoticPressure"
    func_name                   = "calculate_osmotic_pressure"
    description                 = "Calculate osmotic pressure of a solution using the van't Hoff equation: π = i × c × R × T."
    implementation_description  = "Uses the van't Hoff formula for osmotic pressure: π = i·c·R·T, where i is van't Hoff factor, c is molarity (mol/L), R = 0.08206 L·atm/(K·mol), and T is temperature in Kelvin."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Colligative Properties", "Osmotic Pressure", "Solution Chemistry", "Physical Chemistry"]
    required_envs               = []

    code_input_sig = [
        ("molarity_c",        "float", "N/A",       "Molarity of solution in mol/L (mol solute / L solution)."),
        ("temperature_k",     "float", "298.15",     "Temperature in Kelvin."),
        ("vanthoff_factor_i", "float", "1.0",        "van't Hoff factor (i=1 for nonelectrolyte, i≈2 for NaCl, etc.)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'molarity [T(K)] [i]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing osmotic_pressure_atm, osmotic_pressure_pa, and explanation."),
    ]

    examples = [
        {
            "code_input": {
                "molarity_c": 0.1,
                "temperature_k": 298.15,
                "vanthoff_factor_i": 1.0,
            },
            "text_input": {
                "input_params": "0.1 298.15 1.0",
            },
            "output": {
                "result": {
                    "osmotic_pressure_atm": 2.447,
                    "osmotic_pressure_pa": 247900.0,
                }
            },
        },
        # Glucose example
        {
            "code_input": {
                "molarity_c": 0.05,
                "temperature_k": 310.15,
                "vanthoff_factor_i": 1.0,
            },
            "text_input": {
                "input_params": "0.05 310.15 1.0",
            },
            "output": {
                "result": {
                    "osmotic_pressure_atm": 1.273,
                    "osmotic_pressure_pa": 129000.0,
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
        """初始化气体常数。"""
        self.R_atm = 0.08206   # L·atm/(K·mol)
        self.R_J = 8.314462618  # J/(mol·K) = Pa·m³/(mol·K)

    def _run_base(
        self,
        molarity_c: float,
        temperature_k: float = 298.15,
        vanthoff_factor_i: float = 1.0,
    ) -> Dict[str, Any]:
        """核心逻辑：π = i × c × R × T"""
        if molarity_c < 0:
            raise ChemMCPError("Molarity cannot be negative.")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")
        if vanthoff_factor_i < 0:
            raise ChemMCPError("van't Hoff factor cannot be negative.")

        # π (atm) = i * c (mol/L) * R (L·atm/(K·mol)) * T (K)
        pi_atm = vanthoff_factor_i * molarity_c * self.R_atm * temperature_k

        # 转换为 Pa: 1 atm = 101325 Pa
        pi_pa = pi_atm * 101325.0

        # 也用 R_J 计算（c 需要转换为 mol/m³）
        # c (mol/L) = c (mol/m³) / 1000
        # π (Pa) = i * c_m3 * R_J * T = i * (c_L * 1000) * R_J * T
        pi_pa_check = vanthoff_factor_i * molarity_c * 1000.0 * self.R_J * temperature_k

        logger.info(f"Osmotic pressure: π={pi_atm:.4f} atm ({pi_pa:.1f} Pa) "
                     f"(c={molarity_c} M, T={temperature_k}K, i={vanthoff_factor_i})")

        return {
            "osmotic_pressure_atm": round(pi_atm, 6),
            "osmotic_pressure_pa": round(pi_pa, 2),
            "osmotic_pressure_bar": round(pi_atm * 1.01325, 4),
            "osmotic_pressure_mmHg": round(pi_atm * 760.0, 4),
            "parameters_used": {
                "molarity_M": molarity_c,
                "temperature_K": temperature_k,
                "vanthoff_factor_i": vanthoff_factor_i,
                "R_L_atm_mol_K": self.R_atm,
            },
            "explanation": (
                f"π = i × c × R × T = {vanthoff_factor_i} × {molarity_c} × {self.R_atm} × {temperature_k} "
                f"= {pi_atm:.4f} atm = {pi_pa:.1f} Pa = {pi_atm * 760.0:.2f} mmHg."
            ),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            kwargs = {"molarity_c": float(parts[0])}
            if len(parts) > 1: kwargs["temperature_k"] = float(parts[1])
            if len(parts) > 2: kwargs["vanthoff_factor_i"] = float(parts[2])
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'molarity [T(K)] [i]'")
