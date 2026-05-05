import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SelectivityFactorCalculator(BaseTool):
    """
    选择性因子 α 计算工具。
    计算色谱分离的选择性因子（分离因子）α，评估相邻峰的分离能力，指导方法优化。
    """
    __version__ = "0.1.0"
    name = "SelectivityFactorCalculator"
    func_name = "calculate_selectivity_factor"
    description = "Calculate chromatographic selectivity factor (separation factor) α between analyte pairs and predict resolution for method optimization."
    implementation_description = (
        "Calculates selectivity factor α = k₂'/k₁' = (tR2-t₀)/(tR1-t₀) between adjacent or specified peak pairs. "
        "Also computes: resolution Rs using the fundamental resolution equation, "
        "required N for baseline resolution at given α and k', "
        "and retention window analysis. Provides optimization guidance based on thermodynamic "
        "(mobile phase composition, temperature) vs kinetic (N) approaches."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "Selectivity Factor", "Alpha", "Resolution", "HPLC", "Method Development"]
    required_envs = []

    code_input_sig = [
        ("peak_retention_times", "list", "N/A", "List of retention times in minutes [tR1, tR2, ...]."),
        ("dead_time_min", "float", "N/A", "Column dead time t₀ in minutes."),
        ("target_resolution_Rs", "float", "1.5", "Desired resolution (default 1.5 = baseline)."),
        ("column_efficiency_N", "None", "None", "Current column efficiency (theoretical plates)."),
        ("temperature_K", "None", "None", "Column temperature in Kelvin (for van't Hoff analysis)."),
        ("temperature2_K", "None", "None", "Second temperature for van't Hoff selectivity prediction."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'tRs=3.5,5.2,8.1 t0=0.8 N=10000' or 'tR1=3.5 tR2=5.2 t0=0.8'."),
    ]

    output_sig = [
        ("selectivity_analysis", "dict", "Complete selectivity analysis including α values, resolution predictions, required N, and optimization strategies."),
    ]

    examples = [
        {
            "code_input": {
                "peak_retention_times": [3.52, 5.18, 8.33],
                "dead_time_min": 0.82,
                "column_efficiency_N": 12000,
            },
            "text_input": {"input_params": "tRs=3.52,5.18,8.33 t0=0.82 N=12000"},
            "output": {"selectivity_analysis": {"alpha_values": [...], "resolution_predictions": [...]}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_alpha(self, tR1: float, tR2: float, t0: float) -> float:
        """α = k₂'/k₁' = (tR2-t0)/(tR1-t0)."""
        if tR1 <= t0:
            raise ChemMCPError(f"Peak 1 tR ({tR1}) must be greater than dead time t0 ({t0}).")
        return (tR2 - t0) / (tR1 - t0)

    def _calc_resolution(self, alpha: float, k_avg: float, N: int) -> float:
        """Fundamental resolution equation: Rs = (√N/4)·(α-1/α)·(k_avg/(1+k_avg))."""
        if alpha <= 1 or k_avg < 0:
            return 0.0
        return (math.sqrt(N) / 4.0) * ((alpha - 1) / alpha) * (k_avg / (1 + k_avg))

    def _calc_required_N(self, alpha: float, k_avg: float, target_Rs: float) -> Optional[int]:
        """Calculate N needed for target resolution."""
        if alpha <= 1 or k_avg < 0 or target_Rs <= 0:
            return None
        # Rearrange: N = [4·Rs·α/(α-1) · (1+k_avg)/k_avg]²
        term = (4 * target_Rs * alpha / (alpha - 1)) * ((1 + k_avg) / k_avg)
        return int(math.ceil(term ** 2))

    def _assess_alpha(self, alpha: float) -> dict:
        """Interpret α value."""
        if alpha < 1.01:
            return {"level": "co-elution", "emoji": "🔴",
                    "description": "Co-elution risk — α ≈ 1 means no selectivity difference",
                    "difficulty": "Very difficult to resolve"}
        elif alpha < 1.05:
            return {"level": "very_poor", "emoji": "🔶",
                    "description": "Very poor selectivity — requires very high N",
                    "difficulty": "Difficult — needs N > 25000"}
        elif alpha < 1.10:
            return {"level": "poor", "emoji": "⚠️",
                    "description": "Poor selectivity — challenging separation",
                    "difficulty": "Moderately difficult"}
        elif alpha < 1.20:
            return {"level": "moderate", "emoji": "🟡",
                    "description": "Moderate selectivity — achievable with good column",
                    "difficulty": "Achievable with standard HPLC"}
        elif alpha < 1.50:
            return {"level": "good", "emoji": "🟢",
                    "description": "Good selectivity — straightforward separation",
                    "difficulty": "Easy — routine HPLC"}
        else:
            return {"level": "excellent", "emoji": "✅",
                    "description": "Excellent selectivity — well-resolved peaks",
                    "difficulty": "Trivially separable"}

    def _calc_thermodynamic_contributions(self, alpha: float, T1_K: Optional[float],
                                           T2_K: Optional[float]) -> Optional[dict]:
        """
        Estimate thermodynamic contributions to selectivity.
        Uses van't Hoff approximation: ln(k) = -ΔH°/RT + ΔS°/R
        Selectivity change with temperature depends on enthalpy difference.
        """
        if T1_K is None:
            return None

        result = {"temperature_T1_K": T1_K}

        if T2_K is not None and T2_K != T1_K:
            # Approximate: d(ln α)/d(1/T) = -Δ(ΔH°)/R
            # Without specific data, estimate typical effect
            # For most RP-HPLC: α changes ~1-3% per 10°C near room temperature
            delta_T = abs(T2_K - T1_K)
            typical_change_per_10K = 0.02  # 2% change in α per 10K
            estimated_change = typical_change_per_10K * (delta_T / 10)

            if T2_K > T1_K:
                # Higher T usually reduces α slightly in RP-HPLC
                alpha_est = alpha * (1 - estimated_change)
            else:
                alpha_est = alpha * (1 + estimated_change)

            result["temperature_T2_K"] = T2_K
            result["estimated_alpha_at_T2"] = round(alpha_est, 4)
            result["change_percent"] = round((alpha_est - alpha) / alpha * 100, 2)
            result["note"] = (
                "Higher temperature typically reduces selectivity in reversed-phase HPLC. "
                "Lower temperature may improve α but increases backpressure."
            )

        return result

    def _run_base(self, peak_retention_times: list,
                  dead_time_min: float,
                  target_resolution_Rs: float = 1.5,
                  column_efficiency_N: Optional[int] = None,
                  temperature_K: Optional[float] = None,
                  temperature2_K: Optional[float] = None) -> dict:
        """Core logic."""

        if len(peak_retention_times) < 2:
            raise ChemMCPError("Need at least 2 retention times.")

        sorted_tRs = sorted(peak_retention_times)
        n_peaks = len(sorted_tRs)

        # Calculate all adjacent-pair α values
        pair_analyses = []
        for i in range(n_peaks - 1):
            tR1, tR2 = sorted_tRs[i], sorted_tRs[i + 1]
            alpha = self._calc_alpha(tR1, tR2, dead_time_min)
            k1 = (tR1 - dead_time_min) / dead_time_min
            k2 = (tR2 - dead_time_min) / dead_time_min
            k_avg = (k1 + k2) / 2

            # Resolution (actual or predicted)
            Rs_actual = None
            if column_efficiency_N:
                Rs_actual = self._calc_resolution(alpha, k_avg, column_efficiency_N)

            # Required N for target Rs
            N_required = self._calc_required_N(alpha, k_avg, target_resolution_Rs)

            assessment = self._assess_alpha(alpha)

            pair_analyses.append({
                "peak_pair": f"{i+1}/{i+2}",
                "tR1_min": tR1,
                "tR2_min": tR2,
                "k1_prime": round(k1, 4),
                "k2_prime": round(k2, 4),
                "alpha": round(alpha, 4),
                "k_average": round(k_avg, 4),
                "assessment": assessment,
                "resolution": {
                    "predicted_Rs": round(Rs_actual, 3) if Rs_actual is not None else None,
                    "baseline_separated": (Rs_actual or 0) >= target_resolution_Rs if Rs_actual is not None else None,
                    "required_N_for_baseline": N_required,
                    "current_N": column_efficiency_N,
                    "N_sufficient": (column_efficiency_N or 0) >= (N_required or float('inf')) if column_efficiency_N and N_required else None,
                } if column_efficiency_N or N_required else None,
            })

        # Overall worst-case (critical pair)
        critical_pair = min(pair_analyses, key=lambda p: p["alpha"])
        min_alpha = critical_pair["alpha"]

        # Thermodynamic analysis
        thermo = self._calc_thermodynamic_contributions(min_alpha, temperature_K, temperature2_K)

        # Optimization strategy
        strategy = self._get_optimization_strategy(min_alpha, critical_pair, column_efficiency_N)

        result = {
            "selectivity_analysis": {
                "summary": {
                    "n_peaks": n_peaks,
                    "dead_time_t0_min": dead_time_min,
                    "minimum_alpha": round(min_alpha, 4),
                    "critical_pair": f"Peak {pair_analyses.index(critical_pair)+1}/{pair_analyses.index(critical_pair)+2}",
                    "overall_selectivity": self._assess_alpha(min_alpha)["description"],
                },
                "pairwise_analysis": pair_analyses,
                "thermodynamic_analysis": thermo,
                "optimization_strategy": strategy,
            }
        }
        return result

    def _get_optimization_strategy(self, min_alpha: float, critical_pair: dict,
                                    N: Optional[int]) -> List[str]:
        recs = []

        if min_alpha >= 1.5:
            recs.append("✅ Excellent selectivity — no special optimization needed.")
            return recs

        if min_alpha < 1.05:
            recs.append("🔴 Critical pair has very low α — consider:")
            recs.append("  1. Change stationary phase (different selectivity)")
            recs.append("  2. Use different mobile phase pH/ion-pair reagent")
            recs.append("  3. Try HILIC/SFC mode instead of RPLC")
        elif min_alpha < 1.10:
            recs.append("⚠ Low α — try these approaches in order:")
            recs.append("  1. Optimize mobile phase organic modifier type (MeOH vs ACN)")
            recs.append("  2. Fine-tune pH (±0.5-1 unit around pKa)")
            recs.append("  3. Adjust column temperature (lower often helps)")
            recs.append("  4. Consider phenyl-hexyl, polar-embedded, or chiral phases")
        elif min_alpha < 1.20:
            recs.append("🟡 Moderate α — standard optimizations should suffice:")
            recs.append("  1. Optimize %B (organic strength) gradient")
            recs.append("  2. Temperature adjustment (±10°C)")
            recs.append("  3. Use longer column or smaller particles for higher N")

        if N and min_alpha > 1.05:
            req_N = self._calc_required_N(
                min_alpha, critical_pair.get("k_average", 2), 1.5)
            if req_N and N < req_N:
                recs.append(f"📊 Current N ({N}) below required ({req_N}) for baseline resolution.")
                recs.append(f"  → Need {req_N//N + 1}x more plates: use longer/smaller-particle column.")

        return recs[:7]

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                if k == "tRs":
                    kwargs["peak_retention_times"] = [float(x) for x in v.split(",")]
                elif k == "tR1":
                    kwargs.setdefault("peak_retention_times", []).append(float(v))
                elif k == "tR2":
                    kwargs.setdefault("peak_retention_times", []).append(float(v))
                elif k == "t0":
                    kwargs["dead_time_min"] = float(v)
                elif k == "N":
                    kwargs["column_efficiency_N"] = int(float(v))
                elif k == "Rs_target":
                    kwargs["target_resolution_Rs"] = float(v)
                elif k == "T":
                    kwargs["temperature_K"] = float(v) + 273.15
        # Ensure we have a list
        if isinstance(kwargs.get("peak_retention_times"), list) and len(kwargs["peak_retention_times"]) == 1:
            pass  # need at least 2
        return self._run_base(**kwargs)
