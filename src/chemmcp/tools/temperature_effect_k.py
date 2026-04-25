import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TemperatureEffectK(BaseTool):
    """
    计算温度对平衡常数的影响（van't Hoff 方程）。
    ln(K2/K1) = -ΔH°/R × (1/T2 - 1/T1)
    """
    __version__                = "0.1.0"
    name                       = "TemperatureEffectK"
    func_name                  = "temperature_effect_K"
    description                = "Calculate the effect of temperature on equilibrium constant using van't Hoff equation."
    implementation_description = "Uses integrated van't Hoff equation: ln(K2/K1) = -ΔH°/R × (1/T2 - 1/T1). Assumes ΔH° is constant over temperature range."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Thermodynamics", "van't Hoff", "Equilibrium", "Temperature"]
    required_envs              = []

    code_input_sig             = [
        ("K1",                "float", "N/A",   "Equilibrium constant at temperature T1."),
        ("T1",                "float", "N/A",   "Initial temperature in Kelvin."),
        ("T2",                "float", "N/A",   "Final temperature in Kelvin."),
        ("delta_h",           "float", "N/A",   "Standard enthalpy change of reaction in kJ/mol."),
        ("r_gas",             "float", "8.314", "Gas constant in J/(mol·K). Default: 8.314."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Space-separated: 'K1 T1 T2 delta_h [r_gas]'."),
    ]

    output_sig                 = [
        ("K2",                "float", "Equilibrium constant at temperature T2."),
        ("ln_ratio",          "float", "Natural log of K2/K1."),
        ("direction",         "str",   "'increases' if K2 > K1, 'decreases' if K2 < K1, 'unchanged' if equal."),
        ("explanation",       "str",   "Brief explanation based on Le Chatelier principle."),
    ]

    examples                   = [
        {
            "code_input": {
                "K1": 4.0,
                "T1": 300.0,
                "T2": 400.0,
                "delta_h": 40.0,
                "r_gas": 8.314,
            },
            "text_input": {
                "input_params": "4.0 300 400 40",
            },
            "output": {
                "K2": 10.19,
                "ln_ratio": 1.605,
                "direction": "increases",
                "explanation": "Endothermic reaction (ΔH>0): increasing T shifts equilibrium toward products, K increases.",
            }
        },
        {
            "code_input": {
                "K1": 100.0,
                "T1": 300.0,
                "T2": 400.0,
                "delta_h": -50.0,
                "r_gas": 8.314,
            },
            "text_input": {
                "input_params": "100 300 400 -50",
            },
            "output": {
                "K2": 7.39,
                "ln_ratio": -2.606,
                "direction": "decreases",
                "explanation": "Exothermic reaction (ΔH<0): increasing T shifts equilibrium toward reactants, K decreases.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, K1: float, T1: float, T2: float, delta_h: float,
                  r_gas: float = 8.314) -> dict:
        if T1 <= 0 or T2 <= 0:
            raise ChemMCPError("Temperatures must be positive in Kelvin.")
        if K1 <= 0:
            raise ChemMCPError("K1 must be positive.")

        # ln(K2/K1) = -ΔH°/R * (1/T2 - 1/T1); ΔH in kJ/mol → J/mol
        ln_ratio = (-delta_h * 1000.0 / r_gas) * (1.0 / T2 - 1.0 / T1)
        K2 = K1 * math.exp(ln_ratio)

        if K2 > K1 * 1.0001:
            direction = "increases"
        elif K2 < K1 * 0.9999:
            direction = "decreases"
        else:
            direction = "unchanged"

        # Generate explanation
        if delta_h > 0:
            explanation = (
                f"Endothermic reaction (ΔH={delta_h} kJ/mol > 0): "
                f"{'increasing' if T2 > T1 else 'decreasing'} T shifts equilibrium toward products, K {'increases' if T2 > T1 else 'decreases'}."
            )
        elif delta_h < 0:
            explanation = (
                f"Exothermic reaction (ΔH={delta_h} kJ/mol < 0): "
                f"{'increasing' if T2 > T1 else 'decreasing'} T shifts equilibrium toward reactants, K {'decreases' if T2 > T1 else 'increases'}."
            )
        else:
            explanation = "Thermoneutral reaction (ΔH≈0): temperature change has negligible effect on K."

        logger.info(f"van't Hoff: K1={K1} @ {T1}K → K2={K2:.4f} @ {T2}K (ΔH={delta_h})")
        return {
            "K2": round(K2, 4),
            "ln_ratio": round(ln_ratio, 4),
            "direction": direction,
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 4:
                raise ValueError("Need at least K1 T1 T2 delta_h")
            K1 = float(parts[0])
            T1 = float(parts[1])
            T2 = float(parts[2])
            dh = float(parts[3])
            rg = float(parts[4]) if len(parts) > 4 else 8.314
            return self._run_base(K1, T1, T2, dh, rg)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
