import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HeatCapacityLookup(BaseTool):
    """
    查询物质的定压/定容热容及温度依赖关系（Shomate多项式拟合）。
    
    Shomate方程: Cp° = a + bT + cT² + dT³ + e/T²
    单位: J/(mol·K)
    """
    __version__ = "0.1.0"
    name = "HeatCapacityLookup"
    func_name = "lookup_heat_capacity"
    description = "Look up heat capacity (Cp or Cv) for common gases using Shomate polynomial equation with temperature dependence."
    implementation_description = "Uses Shomate polynomial coefficients: Cp = a + bT + cT² + dT³ + e/T². Includes data for N2, O2, H2, CO2, H2O(g), NH3, CH4, Ar, CO. Cv estimated via Cp - R for ideal gases."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Heat Capacity", "Physical Chemistry", "Shomate Equation"]
    required_envs = []

    # Shomate coefficients (a, b, c, d, e) for T in K, Cp in J/(mol·K)
    # Valid range: 298-1200K (approximate)
    _shomate_data = {
        "N2":   {"a": 28.98,  "b": 0.01584e-2,  "c": -0.5735e-5,  "d": 2.884e-9,   "e": -0.4153e5},
        "O2":   {"a": 29.66,  "b": 0.00614e-2,   "c": -0.1460e-5,  "d": 1.202e-9,   "e": -0.8170e5},
        "H2":   {"a": 27.28,  "b": 0.03264e-2,   "c": -0.5046e-5,  "d": -0.1687e-9,  "e": 1.1288e5},
        "CO":   {"a": 28.11,  "b": 0.01675e-2,   "c": -0.5586e-5,  "d": 2.862e-9,   "e": -0.3660e5},
        "CO2":  {"a": 24.99,  "b": 0.05519e-2,   "c": -0.3383e-5,  "d": 7.954e-9,   "e": -0.7281e5},
        "H2O(g)": {"a": 30.09,"b": 0.06832e-2,   "c": 0.7794e-5,   "d": -6.583e-9,  "e": 0.3602e5},
        "NH3":  {"a": 26.49,  "b": 0.02360e-2,   "c": 0.1707e-5,    "d": 1.478e-9,   "e": -1.122e5},
        "CH4":  {"a": 18.89,  "b": 0.05202e-2,   "c": 1.1987e-5,    "d": -11.355e-9, "e": -2.824e5},
        "Ar":   {"a": 20.80,  "b": 0.0,          "c": 0.0,         "d": 0.0,        "e": 0.0},       # monatomic ideal gas
        "SO2":  {"a": 25.72,  "b": 0.05796e-2,   "c": -0.3810e-5,  "d": 8.761e-9,   "e": -0.7937e5},
        "NO":   {"a": 29.35,  "b": 0.00985e-2,   "c": -0.3460e-5,  "d": 2.682e-9,   "e": -0.4850e5},
    }

    code_input_sig = [
        ("substance", "str", "N/A", "Substance name (case-insensitive): N2, O2, H2, CO, CO2, H2O(g), NH3, CH4, Ar, SO2, NO."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("cp_type", "str", "cp", "Heat capacity type: 'cp' or 'cv'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'substance [temperature_k] [cp_type]'. Example: 'N2 298.15 cp'"),
    ]

    output_sig = [
        ("value", "float", "Heat capacity value in J/(mol·K)."),
        ("unit", "str", "Unit of the result."),
        ("method", "str", "Calculation method used."),
        ("coefficients", "str", "Shomate coefficients used (a,b,c,d,e)."),
    ]

    examples = [
        {
            "code_input": {
                "substance": "N2",
                "temperature_k": 298.15,
                "cp_type": "cp",
            },
            "text_input": {
                "input_params": "N2 298.15 cp",
            },
            "output": {
                "value": 29.125,
                "unit": "J/(mol·K)",
                "method": "Shomate polynomial at T=298.15 K",
                "coefficients": "a=28.98, b=1.584e-04, c=-5.735e-06, d=2.884e-09, e=-41530",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R_gas = 8.314  # J/(mol·K) for ideal gas correction

    def _run_base(self, substance: str, temperature_k: float = 298.15, cp_type: str = "cp") -> dict:
        """Core logic: calculate Cp/Cv using Shomate equation."""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")

        key = substance.strip()
        # Case-insensitive lookup
        if key not in self._shomate_data:
            # Try case-insensitive match
            for k in self._shomate_data:
                if k.lower() == key.lower():
                    key = k
                    break
            else:
                available = ", ".join(sorted(self._shomate_data.keys()))
                raise ChemMCPError(f"Unknown substance '{substance}'. Available: {available}")

        coeffs = self._shomate_data[key]
        a, b, c, d, e = coeffs["a"], coeffs["b"], coeffs["c"], coeffs["d"], coeffs["e"]

        T = temperature_k
        cp = a + b * T + c * T**2 + d * T**3 + e / T**2

        if cp_type.lower().strip() == "cv":
            # For ideal gas: Cv = Cp - R
            value = cp - self.R_gas
            unit = "J/(mol·K)"
            method = f"Cv via Shomate Cp - R at T={T} K"
        else:
            value = cp
            unit = "J/(mol·K)"
            method = f"Shomate polynomial at T={T} K"

        coeff_str = f"a={a}, b={b:.6e}, c={c:.6e}, d={d:.6e}, e={e:.1f}"

        logger.info(f"HeatCapacity: {substance} @ {T} K, {cp_type.upper()} = {value:.4f} {unit}")
        return {
            "value": round(value, 4),
            "unit": unit,
            "method": method,
            "coefficients": coeff_str,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            substance = parts[0]
            temp = float(parts[1]) if len(parts) > 1 else 298.15
            cp_type = parts[2] if len(parts) > 2 else "cp"
            return self._run_base(substance, temp, cp_type)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'substance [temperature_k] [cp_type]'")
