import logging
import math
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SpontaneityCheck(BaseTool):
    """
    判断反应自发性（基于 ΔG 判据）。
    ΔG < 0: 自发 (spontaneous)
    ΔG = 0: 平衡态 (equilibrium)
    ΔG > 0: 非自发 (non-spontaneous, reverse direction spontaneous)

    同时提供温度对自发性的影响分析（焓驱/熵驱反应）。
    """
    __version__ = "0.1.0"
    name = "SpontaneityCheck"
    func_name = "check_spontaneity"
    description = "Determine reaction spontaneity using Gibbs free energy criterion (ΔG). Analyzes whether a reaction is spontaneous at given temperature and provides thermodynamic driving force analysis."
    implementation_description = "Applies the Gibbs criterion: ΔG < 0 → spontaneous; ΔG = 0 → equilibrium; ΔG > 0 → non-spontaneous. Also classifies reactions as enthalpy-driven or entropy-driven and calculates the temperature threshold where spontaneity changes."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spontaneity", "Gibbs Energy", "Thermodynamics", "Reaction Feasibility"]
    required_envs = []

    code_input_sig = [
        ("delta_g", "float", "N/A", "Gibbs free energy change in kJ/mol (ΔG). Can be positive or negative."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (default: 298.15 K)."),
        ("delta_h", "float", "None", "Optional: Enthalpy change in kJ/mol (ΔH). If provided, enables deeper analysis of driving forces."),
        ("delta_s", "float", "None", "Optional: Entropy change in J/(mol·K) (ΔS). If provided with ΔH, enables temperature-dependence analysis."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: delta_g [temperature_k] [delta_h] [delta_s]. E.g., '-50.2 298 -120.5 150' or just '-17.5'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with is_spontaneous, direction, delta_g, temperature_k, criterion, driving_force_analysis, equilibrium_analysis, and optional temperature_threshold."),
    ]

    examples = [
        {
            "code_input": {
                "delta_g": -50.2,
                "temperature_k": 298.15,
                "delta_h": -120.5,
                "delta_s": 150.0,
            },
            "text_input": {
                "input_str": "-50.2 298.15 -120.5 150",
            },
            "output": {
                "result": {
                    "is_spontaneous": True,
                    "direction": "forward (products favored)",
                    "delta_g_kj_per_mol": -50.2,
                    "temperature_k": 298.15,
                    "criterion": "ΔG = -50.2 kJ/mol < 0 → SPONTANEOUS as written",
                    "driving_force": "enthalpy-driven AND entropy-driven (both ΔH < 0 and ΔS > 0 favor spontaneity)",
                    "equilibrium_analysis": "K >> 1 at this temperature; reaction strongly favors products",
                    "temperature_range": "Spontaneous at ALL temperatures (both enthalpy and entropy favor forward direction)",
                }
            },
        },
        {
            "code_input": {
                "delta_g": +25.8,
                "temperature_k": 298.15,
                "delta_h": +80.3,
                "delta_s": -120.0,
            },
            "text_input": {
                "input_str": "+25.8 298.15 +80.3 -120",
            },
            "output": {
                "result": {
                    "is_spontaneous": False,
                    "direction": "reverse (reactants favored)",
                    "delta_g_kj_per_mol": 25.8,
                    "temperature_k": 298.15,
                    "criterion": "ΔG = +25.8 kJ/mol > 0 → NON-SPONTANEOUS as written (reverse is spontaneous)",
                    "driving_force": "Both enthalpy (+ΔH) and entropy (-ΔS) oppose the forward reaction",
                    "equilibrium_analysis": "K << 1; reactants heavily favored",
                    "temperature_range": "Non-spontaneous at ALL temperatures (reverse always favored)",
                }
            },
        },
        {
            "code_input": {
                "delta_g": 0.0,
                "temperature_k": 373.15,
                "delta_h": 40.7,
                "delta_s": 109.1,
            },
            "text_input": {
                "input_str": "0.0 373.15 40.7 109.1",
            },
            "output": {
                "result": {
                    "is_spontaneous": False,  # technically at equilibrium
                    "direction": "at EQUILIBRIUM",
                    "delta_g_kj_per_mol": 0.0,
                    "temperature_k": 373.15,
                    "criterion": "ΔG ≈ 0 → SYSTEM AT EQUILIBRIUM",
                    "driving_force": "enthalpy-opposed but entropy-driven; exactly balanced at T = 373.15 K",
                    "equilibrium_analysis": "K = 1 at this temperature",
                    "temperature_threshold": {"T_eq": 373.15, "note": "Below T_eq: non-spontaneous; Above T_eq: spontaneous"},
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize constants."""
        self.R = 8.314e-3  # Gas constant in kJ/(mol·K)

    def _run_base(self, delta_g: float, temperature_k: float = 298.15,
                  delta_h: Optional[float] = None, delta_s: Optional[float] = None) -> dict:
        """Core logic: determine spontaneity from ΔG."""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")

        # Tolerance for treating as zero
        EPSILON = 0.01  # kJ/mol

        # Basic determination
        if abs(delta_g) <= EPSILON:
            verdict = "EQUILIBRIUM"
            is_spon = None  # neither strictly true nor false
            direction = "at equilibrium (no net change)"
            criterion_str = f"ΔG = {delta_g:.4f} kJ/mol ≈ 0 → SYSTEM AT EQUILIBRIUM"
        elif delta_g < 0:
            verdict = "SPONTANEOUS"
            is_spon = True
            direction = "forward (products favored)"
            criterion_str = f"ΔG = {delta_g:.2f} kJ/mol < 0 → SPONTANEOUS as written"
        else:
            verdict = "NON-SPONTANEOUS"
            is_spon = False
            direction = "reverse (reactants favored)"
            criterion_str = f"ΔG = {delta_g:+.2f} kJ/mol > 0 → NON-SPONTANEOUS as written (reverse direction is spontaneous)"

        # Equilibrium constant estimate
        K_str = self._estimate_K(delta_g, temperature_k)

        # Driving force analysis
        driving_force = self._analyze_driving_force(delta_h, delta_s)

        # Temperature range analysis (if both ΔH and ΔS provided)
        temp_analysis = self._analyze_temperature_dependence(delta_h, delta_s, temperature_k)

        result = {
            "is_spontaneous": is_spon if is_spon is not None else bool(is_spon),
            "verdict": verdict,
            "direction": direction,
            "delta_g_kj_per_mol": round(delta_g, 4),
            "temperature_k": round(temperature_k, 2),
            "criterion": criterion_str,
            "driving_force": driving_force,
            "equilibrium_analysis": f"K = {K_str}" if K_str else "K cannot be reliably estimated",
        }

        if temp_analysis:
            result["temperature_range"] = temp_analysis

        return result

    def _run_text(self, input_str: str) -> dict:
        """Parse space-separated text input."""
        try:
            parts = input_str.strip().split()
            if not parts:
                raise ValueError("Empty input")

            delta_g = float(parts[0])
            temp = float(parts[1]) if len(parts) > 1 else 298.15
            dh = float(parts[2]) if len(parts) > 2 else None
            ds = float(parts[3]) if len(parts) > 3 else None

            return self._run_base(delta_g, temp, dh, ds)
        except ValueError as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Expected format: 'delta_g [temperature_k] [delta_h] [delta_s]' (all numeric values)"
            )

    def _estimate_K(self, delta_g: float, T: float) -> str:
        """Estimate equilibrium constant from ΔG° = -RT ln K."""
        if abs(delta_g / (self.R * T)) > 500:
            return "extreme" if delta_g < 0 else "≈ 0"
        try:
            K = math.exp(-delta_g / (self.R * T))
            if K > 1e10:
                return "very large (K >> 1, products strongly favored)"
            elif K > 1000:
                return f"{K:.2e} (products favored)"
            elif K > 0.001:
                return f"{K:.4g}"
            elif K > 1e-10:
                return f"{K:.2e} (reactants favored)"
            else:
                return "very small (K ≈ 0, reactants strongly favored)"
        except OverflowError:
            return "extreme" if delta_g < 0 else "≈ 0"

    def _analyze_driving_force(self, delta_h: Optional[float], delta_s: Optional[float]) -> str:
        """Analyze what drives the reaction thermodynamically."""
        if delta_h is None and delta_s is None:
            return "No ΔH/ΔS provided — cannot analyze driving forces."

        parts = []
        h_favor_forward = delta_h is not None and delta_h < 0
        s_favor_forward = delta_s is not None and delta_s > 0
        h_oppose_forward = delta_h is not None and delta_h > 0
        s_oppose_forward = delta_s is not None and delta_s < 0

        if h_favor_forward and s_favor_forward:
            return "Enthalpy-driven AND entropy-driven (both ΔH < 0 and ΔS > 0 favor spontaneity). Reaction is spontaneous at ALL temperatures."
        elif h_oppose_forward and s_oppose_forward:
            return "Both enthalpy (+ΔH) and entropy (-ΔS) oppose the forward reaction. Reverse direction is spontaneous at ALL temperatures."
        elif h_favor_forward and s_oppose_forward:
            return "Enthalpy-driven (ΔH < 0 favors spontaneity), but entropy opposes (-ΔS). Reaction is spontaneous only at LOW temperatures where |TΔS| < |ΔH|."
        elif h_oppose_forward and s_favor_forward:
            return "Entropy-driven (ΔS > 0 favors spontaneity), but enthalpy opposes (+ΔH). Reaction is spontaneous only at HIGH temperatures where TΔS > |ΔH|."
        elif delta_h is not None:
            return f"Only ΔH provided ({delta_h:+.2f} kJ/mol): {'exothermic (heat-releasing)' if delta_h < 0 else 'endothermic (heat-absorbing)'}. Need ΔS for complete analysis."
        elif delta_s is not None:
            return f"Only ΔS provided ({delta_s:+.2f} J/(mol·K)): {'entropy-increasing' if delta_s > 0 else 'entropy-decreasing'}. Need ΔH for complete analysis."
        return "Unable to classify."

    def _analyze_temperature_dependence(self, delta_h: Optional[float],
                                         delta_s: Optional[float],
                                         current_T: float) -> Optional[dict]:
        """Analyze how spontaneity changes with temperature."""
        if delta_h is None or delta_s is None or abs(delta_s) < 1e-6:
            return None

        # T_equilibrium where ΔG = 0: T_eq = ΔH / ΔS (with unit conversion: ΔS in J→kJ)
        T_eq = delta_h / (delta_s / 1000.0)  # K

        analysis = {}
        h_sign = "negative (exothermic)" if delta_h < 0 else "positive (endothermic)"
        s_sign = "positive (disorder increases)" if delta_s > 0 else "negative (order increases)"

        if T_eq > 0:
            analysis["T_equilibrium_K"] = round(T_eq, 2)
            analysis["T_equilibrium_C"] = round(T_eq - 273.15, 2)

            if delta_h < 0 and delta_s < 0:
                # Enthalpy-favored, entropy-opposed
                analysis["rule"] = f"Spontaneous when T < {round(T_eq, 2)} K ({round(T_eq - 273.15, 2)} °C); Non-spontaneous when T > T_eq"
                if current_T < T_eq:
                    analysis["current_status"] = f"T = {round(current_T, 2)} K < T_eq → SPONTANEOUS region ✓"
                else:
                    analysis["current_status"] = f"T = {round(current_T, 2)} K > T_eq → NON-SPONTANEOUS region ✗"

            elif delta_h > 0 and delta_s > 0:
                # Entropy-favored, enthalpy-opposed
                analysis["rule"] = f"Spontaneous when T > {round(T_eq, 2)} K ({round(T_eq - 273.15, 2)} °C); Non-spontaneous when T < T_eq"
                if current_T > T_eq:
                    analysis["current_status"] = f"T = {round(current_T, 2)} K > T_eq → SPONTANEOUS region ✓"
                else:
                    analysis["current_status"] = f"T = {round(current_T, 2)} K < T_eq → NON-SPONTANEOUS region ✗"

            else:
                analysis["note"] = f"T_eq = {round(T_eq, 2)} K (theoretical crossover point)"

            analysis["thermodynamic_profile"] = (
                f"ΔH = {delta_h:+.2f} kJ/mol ({h_sign}), "
                f"ΔS = {delta_s:+.2f} J/(mol·K) ({s_sign})"
            )
        else:
            # T_eq ≤ 0 or negative — no meaningful crossover in physical range
            if delta_h < 0 and delta_s > 0:
                analysis["rule"] = "Spontaneous at ALL temperatures (always favorable)"
            elif delta_h > 0 and delta_s < 0:
                analysis["rule"] = "Non-spontaneous at ALL temperatures (reverse always favorable)"
            else:
                analysis["rule"] = "Edge case — consult detailed phase diagram"

        return analysis
