import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CapacityFactorCalculator(BaseTool):
    """
    容量因子 k' 计算工具。
    计算色谱容量因子（保留因子）k'，评估化合物在固定相上的保留行为，优化分离条件。
    """
    __version__ = "0.1.0"
    name = "CapacityFactorCalculator"
    func_name = "calculate_capacity_factor"
    description = "Calculate chromatographic capacity factor (retention factor) k' and related retention parameters for method development and optimization."
    implementation_description = (
        "Calculates capacity factor k' = (tR - t₀)/t₀ from retention time and dead time. "
        "Also computes: adjusted retention t'R, retention ratio, phase ratio β, "
        "distribution constant K = k'·β, and resolution prediction. "
        "Provides retention optimization guidance based on k' range analysis."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "Capacity Factor", "Retention Factor", "HPLC", "Method Development", "k prime"]
    required_envs = []

    code_input_sig = [
        ("retention_time_min", "float", "N/A", "Analyte retention time in minutes."),
        ("dead_time_min", "float", "N/A", "Column dead time (void time) t₀ in minutes."),
        ("peak_width_half_height_min", "None", "None", "Peak width at half height in minutes (for N and Rs estimation)."),
        ("column_length_mm", "None", "None", "Column length in mm (for L/u calculation)."),
        ("flow_rate_mL_min", "None", "None", "Flow rate in mL/min (for linear velocity)."),
        ("column_internal_diameter_mm", "None", "None", "Column ID in mm (for phase ratio)."),
        ("target_k_range", "tuple", "None", "Desired k' range as (min, max) for optimization advice."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'tR=5.2 t0=0.8 Wh=0.15' or 'tR_list=3.1,5.2,8.7 t0=0.8'."),
    ]

    output_sig = [
        ("retention_analysis", "dict", "Complete retention analysis including k', N, resolution potential, and optimization recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "retention_time_min": 5.2,
                "dead_time_min": 0.82,
                "peak_width_half_height_min": 0.15,
            },
            "text_input": {"input_params": "tR=5.2 t0=0.82 Wh=0.15"},
            "output": {"retention_analysis": {"k_prime": 5.34, "status": "optimal"}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_k_prime(self, tR: float, t0: float) -> float:
        """k' = (tR - t0) / t0."""
        if t0 <= 0:
            raise ChemMCPError("Dead time must be positive.")
        return (tR - t0) / t0

    def _assess_k_range(self, k: float) -> dict:
        """Assess k' value against optimal range."""
        if k < 0:
            return {"range": "invalid", "status": "🚨 ERROR — elutes before void peak!", "color": "red"}
        elif k < 0.5:
            return {"range": "too_low", "status": "⚠ Too low (k'<0.5) — poor separation from solvent front",
                    "color": "orange", "action": "Weaken mobile phase or lower temperature"}
        elif k < 1.0:
            return {"range": "low", "status": "Low (0.5-1) — minimal retention",
                    "color": "yellow", "action": "Consider slightly weaker mobile phase"}
        elif k < 2.0:
            return {"range": "acceptable_lower", "status": "Acceptable lower (1-2)",
                    "color": "lightgreen", "action": "OK for fast screening"}
        elif k <= 10.0:
            return {"range": "optimal", "status": "✅ Optimal (2-10) — ideal separation range",
                    "color": "green", "action": "Well within optimal range"}
        elif k <= 20.0:
            return {"range": "high", "status": "High (10-20) — long retention but acceptable",
                    "color": "yellow", "action": "Consider stronger mobile phase to reduce run time"}
        else:
            return {"range": "too_high", "status": "🔶 Very high (k'>20) — excessive retention",
                    "color": "orange", "action": "Strengthen mobile phase significantly or use gradient"}

    def _calc_resolution_prediction(self, k1: float, k2: float, alpha: float,
                                     N: Optional[int]) -> Optional[dict]:
        """Predict Rs between two peaks."""
        if N is None or alpha is None or alpha <= 1:
            return None
        # Rs = (√N/4) · (α-1/α) · (k2_avg/(1+k2_avg))
        k_avg = (k1 + k2) / 2
        rs = (math.sqrt(N) / 4) * ((alpha - 1) / alpha) * (k_avg / (1 + k_avg))
        return {
            "predicted_Rs": round(rs, 3),
            "baseline_resolution": round(rs, 1),
            "baseline_separation": "Baseline" if rs >= 1.5 else "Partial" if rs >= 1.0 else "Not resolved",
        }

    def _run_base(self, retention_time_min: float,
                  dead_time_min: float,
                  peak_width_half_height_min: Optional[float] = None,
                  column_length_mm: Optional[float] = None,
                  flow_rate_mL_min: Optional[float] = None,
                  column_internal_diameter_mm: Optional[float] = None,
                  target_k_range: Optional[tuple] = None) -> dict:
        """Core logic."""

        # Primary calculation
        k_prime = self._calc_k_prime(retention_time_min, dead_time_min)
        tR_adjusted = retention_time_min - dead_time_min  # t'R

        # Derived quantities
        N = None
        if peak_width_half_height_min:
            N = int(round(5.54 * (retention_time_min / peak_width_half_height_min) ** 2))

        # Phase ratio (approximate for cylindrical column)
        beta = None
        if column_internal_diameter_mm:
            # For packed columns, β ≈ Vm/Vs ≈ (void fraction)/(1-void fraction) ≈ 1.5-3 typical
            # More precisely: β depends on stationary phase film thickness
            beta_approx = 1.65  # typical for fully porous C18, 5μm
            K_dist = k_prime * beta_approx
            beta = {"phase_ratio_beta": round(beta_approx, 2),
                     "distribution_constant_K": round(K_dist, 1)}

        # Assessment
        assessment = self._assess_k_range(k_prime)

        # Retention breakdown
        retention_breakdown = {
            "retention_time_tR_min": retention_time_min,
            "dead_time_t0_min": dead_time_min,
            "adjusted_retention_tR_prime_min": round(tR_adjusted, 4),
            "capacity_factor_k_prime": round(k_prime, 4),
            "retention_on_column_fraction": round(tR_adjusted / retention_time_min * 100, 1),
            "void_fraction": round(dead_time_min / retention_time_min * 100, 1),
        }

        result = {
            "retention_analysis": {
                "retention_parameters": retention_breakdown,
                "efficiency": {"theoretical_plates_N": N} if N is not None else None,
                "phase_distribution": beta,
                "k_prime_assessment": assessment,
                "optimization_recommendations": self._get_recommendations(
                    k_prime, assessment, peak_width_half_height_min, target_k_range),
            }
        }
        return result

    def _get_recommendations(self, k: float, assessment: dict,
                              Wh: Optional[float], target: Optional[tuple]) -> List[str]:
        recs = [assessment["status"]]
        action = assessment.get("action")
        if action:
            recs.append(f"💡 Suggestion: {action}")

        if target:
            tmin, tmax = target
            if k < tmin:
                recs.append(f"Target k' range is [{tmin}, {tmax}] — current value ({k:.2f}) below target.")
                recs.append("  → Weaken mobile phase (increase water/aqueous % for RP-HPLC)")
            elif k > tmax:
                recs.append(f"Target k' range is [{tmin}, {tmax}] — current value ({k:.2f}) above target.")
                recs.append("  → Strengthen mobile phase (increase organic % for RP-HPLC)")

        if Wh and k < 1:
            recs.append("⚠ Low k' compounds are most affected by extra-column volume — ensure minimal system void.")

        if k > 30:
            recs.append("Consider gradient elution for this strongly retained compound.")

        return recs[:6]

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                km = {"tR": "retention_time_min", "t0": "dead_time_min",
                      "Wh": "peak_width_half_height_min", "L": "column_length_mm",
                      "F": "flow_rate_mL_min", "ID": "column_internal_diameter_mm"}.get(k, k)
                try: kwargs[km] = float(v)
                except ValueError: kwargs[km] = v
        return self._run_base(**kwargs)
