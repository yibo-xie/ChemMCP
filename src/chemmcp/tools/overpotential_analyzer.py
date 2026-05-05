import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class OverpotentialAnalyzer(BaseTool):
    """
    Overpotential analysis tool.
    
    Computes total overpotential, activation/concentration breakdown,
    Tafel slope, and exchange current density estimation.
    """
    __version__ = "0.1.0"
    name = "OverpotentialAnalyzer"
    func_name = "analyze_overpotential"
    description = "Analyze electrode overpotential including Tafel kinetics, activation/concentration breakdown, and exchange current density estimation."
    implementation_description = "Uses Butler-Volmer and Tafel equations: eta_act=(RT/anF)ln(j/j0), eta_conc=(RT/nF)ln(1-j/jL)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Overpotential", "Tafel", "Kinetics"]
    required_envs = []

    code_input_sig = [
        ("applied_potential_v", "float", "N/A", "Applied electrode potential vs reference (V)."),
        ("equilibrium_potential_v", "float", "0.0", "Equilibrium (Nernst) potential (V). Default 0."),
        ("current_density_ma_cm2", "float", "N/A", "Measured current density (mA/cm2)."),
        ("transfer_coefficient", "float", "0.5", "Charge transfer coefficient alpha (0<alpha<=1). Default 0.5."),
        ("electrons_transferred", "int", "1", "Number of electrons transferred n. Default 1."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin. Default 298.15."),
        ("limiting_current_density", "float", "None", "Limiting current density j_L (mA/cm2). None ignores conc overpotential."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: applied_V eq_V j alpha n T [j_L]."),
    ]

    output_sig = [
        ("overpotential_v", "float", "Total overpotential eta = E_app - E_eq (V)."),
        ("activation_overpotential_v", "float", "Activation (kinetic) overpotential eta_act (V)."),
        ("concentration_overpotential_v", "float", "Concentration overpotential eta_conc (V), or None."),
        ("tafel_slope_mv", "float", "Tafel slope b = ln(10)*RT/anF (mV/dec)."),
        ("exchange_current_density_est", "float", "Estimated exchange current density j0 (mA/cm2)."),
        ("analysis_summary", "str", "Human-readable summary."),
    ]

    examples = [
        {
            "code_input": {
                "applied_potential_v": -0.5,
                "equilibrium_potential_v": 0.0,
                "current_density_ma_cm2": -10.0,
                "transfer_coefficient": 0.5,
                "electrons_transferred": 1,
                "temperature_k": 298.15,
                "limiting_current_density": None,
            },
            "text_input": {
                "input_params": "-0.5 0.0 -10.0 0.5 1 298.15",
            },
            "output": {
                "overpotential_v": -0.5,
                "activation_overpotential_v": -0.1183,
                "concentration_overpotential_v": None,
                "tafel_slope_mv": 59.16,
                "exchange_current_density_est": 0.001,
                "analysis_summary": "Cathodic process with total eta=-0.500V",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.F = 96485.33212
        self.R = 8.314462618

    def _run_base(self, applied_potential_v, equilibrium_potential_v=0.0,
                   current_density_ma_cm2=None, transfer_coefficient=0.5,
                   electrons_transferred=1, temperature_k=298.15,
                   limiting_current_density=None) -> dict:
        eta_total = applied_potential_v - equilibrium_potential_v
        n = electrons_transferred
        alpha = transfer_coefficient
        T = temperature_k

        if alpha <= 0 or alpha > 1:
            raise ChemMCPError("Transfer coefficient must be in (0, 1].")
        if T <= 0:
            raise ChemMCPError("Temperature must be positive.")
        if n < 1:
            raise ChemMCPError("Electrons transferred must be >= 1.")

        tafel_slope_v = (math.log(10) * self.R * T) / (alpha * n * self.F)
        tafel_slope_mv = round(tafel_slope_v * 1000, 2)

        eta_act = None
        j0_est = None
        eta_conc = None

        if current_density_ma_cm2 is not None:
            j = abs(current_density_ma_cm2)
            if j > 0:
                abs_eta = abs(eta_total)
                if abs_eta > 1e-10 and tafel_slope_mv > 1e-10:
                    j0_est = j * (10 ** (-abs_eta / tafel_slope_mv))
                    eta_act_abs = (tafel_slope_mv / 1000.0) * math.log10(max(j / max(j0_est, 1e-30), 1e-30))
                    sign = 1 if current_density_ma_cm2 >= 0 else -1
                    eta_act = sign * eta_act_abs
                    if limiting_current_density is not None:
                        jL = abs(limiting_current_density)
                        ratio = j / jL
                        if ratio < 1:
                            eta_conc_val = (self.R * T / (n * self.F)) * math.log(1 - ratio)
                            eta_conc = sign * eta_conc_val
                else:
                    eta_act = 0.0
                    j0_est = float("inf")

        result = {
            "overpotential_v": round(eta_total, 6),
            "activation_overpotential_v": round(eta_act, 6) if eta_act is not None else None,
            "concentration_overpotential_v": round(eta_conc, 6) if eta_conc is not None else None,
            "tafel_slope_mv": tafel_slope_mv,
            "exchange_current_density_est": round(j0_est, 6) if j0_est is not None else None,
        }

        direction = "anodic" if eta_total > 0 else "cathodic" if eta_total < 0 else "equilibrium"
        parts = [f"{direction} process" if direction != "equilibrium" else "at equilibrium",
                 f"total eta={eta_total:.3f}V"]
        if eta_act is not None:
            parts.append(f"eta_act={eta_act:.3f}V")
        if eta_conc is not None:
            parts.append(f"eta_conc={eta_conc:.3f}V")
        parts.append(f"b={tafel_slope_mv:.1f} mV/dec")
        if j0_est is not None:
            parts.append(f"j0~={j0_est:.4g} mA/cm2")
        result["analysis_summary"] = ", ".join(parts)
        logger.info(f"Overpotential analysis: {result['analysis_summary']}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 5:
                raise ValueError("Need at least 5 params.")
            kwargs = {
                "applied_potential_v": float(parts[0]),
                "equilibrium_potential_v": float(parts[1]),
                "current_density_ma_cm2": float(parts[2]),
                "transfer_coefficient": float(parts[3]),
                "electrons_transferred": int(parts[4]),
            }
            if len(parts) > 5:
                kwargs["temperature_k"] = float(parts[5])
            if len(parts) > 6:
                val = parts[6]
                kwargs["limiting_current_density"] = None if val.upper() == "NONE" else float(val)
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}")
