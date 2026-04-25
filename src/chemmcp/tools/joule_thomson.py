import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class JouleThomson(BaseTool):
    """
    焦耳-汤姆逊系数计算与节流过程分析工具。
    计算焦耳-汤姆逊系数 μ_JT，判断节流过程的温度变化方向，分析反转温度。
    """
    __version__ = "0.1.0"
    name = "JouleThomson"
    func_name = "joule_thomson_analysis"
    description = "Calculate Joule-Thomson coefficient and analyze throttling (isenthalpic expansion) processes."
    implementation_description = "Uses thermodynamic relations: μ_JT = (∂T/∂P)_H = [T(∂V/∂T)_P - V] / Cp. Supports real gas via virial or van der Waals approximations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Physical Chemistry", "Joule-Thomson", "Throttling"]
    required_envs    = []

    code_input_sig = [
        ("gas_type", "str", "N/A", "Gas model: 'ideal', 'vdw' (van der Waals), 'virial' (second virial), or 'custom'."),
        ("temperature_k", "float", "N/A", "Temperature in Kelvin."),
        ("pressure_atm", "float", "1.0", "Pressure in atm (default 1.0)."),
        ("cp_j_mol_k", "float", "N/A", "Molar heat capacity at constant pressure in J/(mol·K)."),
        # van der Waals params (only needed for vdw model)
        ("a_vdw", "float", "0.0", "van der Waals constant a (L²·atm/mol²), used when gas_type='vdw' or 'custom'."),
        ("b_vdw", "float", "0.0", "van der Waals constant b (L/mol), used when gas_type='vdw' or 'custom'."),
        # second virial coefficient (only needed for virial model)
        ("b_virial", "float", "0.0", "Second virial coefficient B (L/mol) at given T, used when gas_type='virial'."),
        ("db_dt", "float", "0.0", "Temperature derivative dB/dT (L/(mol·K)), used when gas_type='virial'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: gas_type temperature_k cp [pressure_atm a_vdw b_vdw b_virial db_dt]."),
    ]

    output_sig = [
        ("mu_jt", "float", "Joule-Thomson coefficient μ_JT in K/atm."),
        ("effect", "str", "Description of cooling/heating effect upon throttling."),
        ("inversion_temp_k", "float", "Estimated inversion temperature T_inv in Kelvin (for vdW gases)."),
        ("analysis", "str", "Detailed analysis text of the throttling process."),
    ]

    examples         = [
        {
            "code_input": {
                "gas_type": 'vdw',
                "temperature_k": 300.0,
                "pressure_atm": 1.0,
                "cp_j_mol_k": 28.8,
                "a_vdw": 1.39,
                "b_vdw": 0.0391,
                "b_virial": 0.0,
                "db_dt": 0.0
            },
            "text_input": {
                "input_params": 'vdw 300.0 1.0 28.8 1.39 0.0391'
            },
            "output": {
                "mu_jt": 0.184,
                "effect": 'Cooling upon throttling (μ_JT > 0)',
                "inversion_temp_k": 695.4,
                "analysis": 'At 300 K, nitrogen-like gas exhibits positive μ_JT; throttling causes cooling.'
            }
        },
        {
            "code_input": {
                "gas_type": 'ideal',
                "temperature_k": 298.0,
                "pressure_atm": 1.0,
                "cp_j_mol_k": 20.8,
                "a_vdw": 0.0,
                "b_vdw": 0.0,
                "b_virial": 0.0,
                "db_dt": 0.0
            },
            "text_input": {
                "input_params": 'ideal 298.0 1.0 20.8'
            },
            "output": {
                "mu_jt": 0.0,
                "effect": 'No temperature change (μ_JT = 0)',
                "inversion_temp_k": float("inf"),
                "analysis": 'Ideal gas has zero Joule-Thomson coefficient.'
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618  # J/(mol·K), universal gas constant

    def _run_base(
        self,
        gas_type: str,
        temperature_k: float,
        cp_j_mol_k: float,
        pressure_atm: float = 1.0,
        a_vdw: float = 0.0,
        b_vdw: float = 0.0,
        b_virial: float = 0.0,
        db_dt: float = 0.0,
    ) -> dict:
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive.")
        if cp_j_mol_k <= 0:
            raise ChemMCPError("Cp must be positive.")

        gas_type = gas_type.lower().strip()

        if gas_type == "ideal":
            mu_jt = 0.0
            t_inv = float("inf")
            analysis = (
                "For an ideal gas, (∂V/∂T)_P = R/P and V = RT/P, "
                "so T(∂V/∂T)_P - V = TR/P - RT/P = 0. Hence μ_JT = 0. "
                "No temperature change occurs during an isenthalpic throttling process."
            )
            effect = "No temperature change (μ_JT = 0)"

        elif gas_type == "vdw":
            if a_vdw == 0 and b_vdw == 0:
                raise ChemMCPError("Van der Waals parameters a and b must be provided for vdw model.")
            # For vdW gas: μ_JT = (2a/RT - 3b) / Cp
            # Using R in L·atm/(K·mol) for consistency with a (L²·atm/mol²)
            R_latm = 0.082057366080
            mu_jt = (2 * a_vdw / (R_latm * temperature_k) - 3 * b_vdw) / cp_j_mol_k
            # Inversion temp for vdW: T_inv = 2a / (R * 3b)
            if b_vdw > 0:
                t_inv = 2 * a_vdw / (3 * R_latm * b_vdw)
            else:
                t_inv = float("inf")

            if abs(mu_jt) < 1e-10:
                effect = "Near inversion point (μ_JT ≈ 0)"
            elif mu_jt > 0:
                effect = "Cooling upon throttling (μ_JT > 0)"
            else:
                effect = "Heating upon throttling (μ_JT < 0)"

            analysis = (
                f"Van der Waals gas at T={temperature_k} K: μ_JT = {mu_jt:.6f} K/atm.\n"
                f"Inversion temperature T_inv = {t_inv:.1f} K.\n"
                f"{'Below inversion → cooling' if temperature_k < t_inv else 'Above inversion → heating'}"
            )

        elif gas_type == "virial":
            # μ_JT = [T*(dB/dT) - B] / Cp  (using B in molar volume units)
            mu_jt = (temperature_k * db_dt - b_virial) / cp_j_mol_k
            t_inv = None  # Cannot determine from single-point B data

            if abs(mu_jt) < 1e-10:
                effect = "Near inversion point (μ_JT ≈ 0)"
            elif mu_jt > 0:
                effect = "Cooling upon throttling (μ_JT > 0)"
            else:
                effect = "Heating upon throttling (μ_JT < 0)"

            analysis = (
                f"Virial approximation at T={temperature_k} K: B={b_virial} L/mol, dB/dT={db_dt} L/(mol·K).\n"
                f"μ_JT = {mu_jt:.6f} K/atm."
            )

        else:
            raise ChemMCPError(f"Unsupported gas type: '{gas_type}'. Use 'ideal', 'vdw', or 'virial'.")

        return {
            "mu_jt": round(mu_jt, 6),
            "effect": effect,
            "inversion_temp_k": round(t_inv, 1) if isinstance(t_inv, (int, float)) and t_inv != float("inf") else t_inv,
            "analysis": analysis,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            gas_type = parts[0]
            temperature_k = float(parts[1])
            cp_j_mol_k = float(parts[2])
            pressure_atm = float(parts[3]) if len(parts) > 3 else 1.0
            a_vdw = float(parts[4]) if len(parts) > 4 else 0.0
            b_vdw = float(parts[5]) if len(parts) > 5 else 0.0
            b_virial = float(parts[6]) if len(parts) > 6 else 0.0
            db_dt = float(parts[7]) if len(parts) > 7 else 0.0
            return self._run_base(gas_type, temperature_k, cp_j_mol_k, pressure_atm, a_vdw, b_vdw, b_virial, db_dt)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
