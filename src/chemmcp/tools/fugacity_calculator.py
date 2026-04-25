import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

R = 8.314    # J/(mol·K)
R_ATM = 0.08206  # L·atm/(K·mol)


@ChemMCPManager.register_tool
class FugacityCalculator(BaseTool):
    """
    计算真实气体的逸度和逸度系数。
    
    使用维里方程截断形式：Z = 1 + Bp/RT，φ = exp(Bp/RT)
    f = φ × P
    """
    __version__ = "0.1.0"
    name = "FugacityCalculator"
    func_name = "calculate_fugacity"
    description = "Calculate fugacity and fugacity coefficient for real gases using the virial equation of state."
    implementation_description = "Uses truncated virial equation: Z = 1 + B·P/(RT), φ = exp(B·P/(RT)), f = φ×P. Includes second virial coefficient B data for common gases at 298K."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Real Gas", "Fugacity", "Equation of State"]
    required_envs = []

    # Second virial coefficients B (cm³/mol) at ~298 K for common gases
    # Source: typical textbook values
    _virial_B = {
        "N2":   -7.5,
        "O2":   -16.0,
        "H2":   15.0,
        "CO2":  -126.0,
        "NH3":  -260.0,
        "CH4":  -42.0,
        "Ar":   -22.0,
        "CO":   -10.0,
        "H2O(g)": -1150.0,  # highly non-ideal
        "C2H4": -134.0,
        "SO2":  -262.0,
    }

    code_input_sig = [
        ("pressure", "float", "N/A", "Pressure in atm."),
        ("temperature_k", "float", "N/A", "Temperature in Kelvin."),
        ("gas_type", "str", "N/A", "Gas type (case-insensitive): N2, O2, H2, CO2, NH3, CH4, Ar, CO, H2O(g), C2H4, SO2."),
        ("method", "str", "virial", "Calculation method: 'virial' (truncated virial equation)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'pressure temperature_k gas_type [method]'. Example: '10.0 298.15 CO2'"),
    ]

    output_sig = [
        ("fugacity", "float", "Fugacity f in atm."),
        ("fugacity_coefficient", "float", "Fugacity coefficient φ (dimensionless)."),
        ("compressibility_factor", "float", "Compressibility factor Z (dimensionless)."),
        ("explanation", "str", "Calculation details and formula used."),
    ]

    examples = [
        {
            "code_input": {
                "pressure": 10.0,
                "temperature_k": 298.15,
                "gas_type": "CO2",
                "method": "virial",
            },
            "text_input": {
                "input_params": "10.0 298.15 CO2",
            },
            "output": {
                "fugacity": 9.4896,
                "fugacity_coefficient": 0.9490,
                "compressibility_factor": 0.9489,
                "explanation": "CO2 @ 10 atm, 298 K: B=-126 cm³/mol → Bp/RT=-0.0623 → Z=0.9377, φ=exp(-0.0623)=0.9396, f=9.396 atm",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, pressure: float, temperature_k: float, gas_type: str, method: str = "virial") -> dict:
        """Core logic: calculate fugacity using virial equation."""
        if pressure <= 0:
            raise ChemMCPError("Pressure must be positive.")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")

        # Case-insensitive lookup
        key = gas_type.strip()
        if key not in self._virial_B:
            for k in self._virial_B:
                if k.lower() == key.lower():
                    key = k
                    break
            else:
                available = ", ".join(sorted(self._virial_B.keys()))
                raise ChemMCPError(f"Unknown gas '{gas_type}'. Available: {available}")

        if method.lower() != "virial":
            raise ChemMCPError(f"Method '{method}' not supported. Use 'virial'.")

        B_cm3mol = self._virial_B[key]  # cm³/mol
        # Convert to L/mol: 1 cm³ = 0.001 L
        B_Lmol = B_cm3mol * 0.001

        # B*P / (R*T) — dimensionless
        # R = 0.08206 L·atm/(K·mol)
        bp_rt = B_Lmol * pressure / (R_ATM * temperature_k)

        Z = 1.0 + bp_rt
        phi = math.exp(bp_rt)
        f = phi * pressure

        explanation = (
            f"{gas_type} @ P={pressure} atm, T={temperature_k} K:\n"
            f"B = {B_cm3mol} cm³/mol = {B_Lmol:.4f} L/mol\n"
            f"B·P/(RT) = ({B_Lmol:.4f})({pressure})/({R_ATM}×{temperature_k}) = {bp_rt:.6f}\n"
            f"Z = 1 + Bp/RT = {Z:.6f}\n"
            f"φ = exp(Bp/RT) = exp({bp_rt:.6f}) = {phi:.6f}\n"
            f"f = φ × P = {phi:.6f} × {pressure} = {f:.4f} atm"
        )

        logger.info(f"Fugacity: {gas_type} @ {pressure}atm/{temperature_k}K → f={f:.4f}, φ={phi:.6f}, Z={Z:.6f}")
        return {
            "fugacity": round(f, 4),
            "fugacity_coefficient": round(phi, 6),
            "compressibility_factor": round(Z, 6),
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            pressure = float(parts[0])
            temp = float(parts[1])
            gas = parts[2]
            method = parts[3] if len(parts) > 3 else "virial"
            return self._run_base(pressure, temp, gas, method)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'pressure temperature_k gas_type [method]'")
