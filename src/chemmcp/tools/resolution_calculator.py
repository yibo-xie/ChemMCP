import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ResolutionCalculator(BaseTool):
    """
    色谱分离度计算与优化建议工具。
    计算色谱峰之间的分离度(Rs)，预测分离效果，并提供优化建议以达到目标分离度。
    """
    __version__ = "0.1.0"
    name = "ResolutionCalculator"
    func_name = "calculate_resolution"
    description = "Calculate chromatographic resolution (Rs) between peak pairs and provide optimization suggestions to achieve target separation."
    implementation_description = (
        "Implements fundamental resolution equation: Rs = (√N/4)·(α-1/α)·(k₂/(1+k_avg)). "
        "Also supports direct calculation from retention times and peak widths, "
        "and provides optimization pathways for N (efficiency), α (selectivity), and k (retention)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "Resolution", "Separation", "Method Optimization", "HPLC"]
    required_envs = []

    code_input_sig = [
        ("calculation_mode", "str", "direct", "Mode: 'direct' (from tR/W), 'fundamental' (from N, α, k), or 'predict' (estimate from structures)."),
        # Direct mode inputs
        ("tR1_min", "None", "None", "Retention time of peak 1 (min) [direct mode]."),
        ("tR2_min", "None", "None", "Retention time of peak 2 (min) [direct mode]."),
        ("W1_min", "None", "None", "Baseline width of peak 1 (min) [direct mode]."),
        ("W2_min", "None", "None", "Baseline width of peak 2 (min) [direct mode]."),
        ("Wh1_min", "None", "None", "Half-height width of peak 1 (min) [direct mode]."),
        ("Wh2_min", "None", "None", "Half-height width of peak 2 (min) [direct mode]."),
        # Fundamental mode inputs
        ("N", "None", "None", "Theoretical plate count [fundamental mode]."),
        ("alpha", "None", "None", "Separation factor α = k₂/k₁ [fundamental mode]."),
        ("k1", "None", "None", "Retention factor of first peak [fundamental mode]."),
        ("k2", "None", "None", "Retention factor of second peak [fundamental mode]."),
        # General
        ("target_rs", "float", "1.5", "Target resolution (baseline: Rs=1.5)."),
        ("column_length_mm", "None", "None", "Column length for optimization suggestions."),
        ("particle_size_um", "None", "None", "Particle size for efficiency estimation."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'tR1=3.2 tR2=3.5 Wb1=0.20 Wb2=0.22' or 'N=10000 alpha=1.10 k2=3.0'."),
    ]

    output_sig = [
        ("resolution", "dict", "Complete resolution analysis including Rs value, baseline assessment, and step-by-step optimization plan if needed."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_mode": "direct",
                "tR1_min": 3.2,
                "tR2_min": 3.5,
                "W1_min": 0.20,
                "W2_min": 0.22,
            },
            "text_input": {"input_params": "tR1=3.2 tR2=3.5 Wb1=0.20 Wb2=0.22"},
            "output": {"resolution": {"Rs": 1.43, "baseline_separated": False}},
        },
        {
            "code_input": {
                "calculation_mode": "fundamental",
                "N": 10000,
                "alpha": 1.10,
                "k2": 3.0,
                "k1": 2.73,
            },
            "text_input": {"input_params": "N=10000 alpha=1.10 k2=3.0"},
            "output": {"resolution": {"Rs": 1.52, "baseline_separated": True}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _rs_direct(self, tR1: float, tR2: float, W1: float, W2: float) -> dict:
        """Rs = 2(tR2 - tR1) / (W1 + W2)."""
        rs = 2 * abs(tR2 - tR1) / (W1 + W2)
        return {"value": round(rs, 3), "method": "Direct (baseline widths)"}

    def _rs_wh(self, tR1: float, tR2: float, Wh1: float, Wh2: float) -> dict:
        """Rs from half-height widths."""
        # Approximate Wb ≈ 1.7 * Wh (for Gaussian peaks)
        Wb1_est = 1.699 * Wh1
        Wb2_est = 1.699 * Wh2
        rs = 2 * abs(tR2 - tR1) / (Wb1_est + Wb2_est)
        return {"value": round(rs, 3), "method": "From half-height widths (Gaussian approximation)"}

    def _rs_fundamental(self, N: float, alpha: float, k1: Optional[float] = None,
                         k2: Optional[float] = None) -> dict:
        """Rs = (√N/4) · ((α-1)/α) · (k₂/(1+k_avg))."""
        if k1 is not None and k2 is not None:
            k_avg = (k1 + k2) / 2
            k_term = k2 / (1 + k_avg)
        elif k2 is not None:
            k_term = k2 / (1 + k2)
        else:
            k_term = 1.0

        term_selectivity = (alpha - 1) / alpha
        term_efficiency = math.sqrt(N) / 4
        rs = term_efficiency * term_selectivity * k_term

        return {
            "value": round(rs, 3),
            "method": "Fundamental resolution equation",
            "breakdown": {
                "efficiency_term": round(term_efficiency, 4),
                "selectivity_term": round(term_selectivity, 4),
                "retention_term": round(k_term, 4),
                "product": f"{term_efficiency:.4f} × {term_selectivity:.4f} × {k_term:.4f}",
            }
        }

    def _optimization_plan(self, current_rs: float, target_rs: float,
                            mode: str, **kwargs) -> List[Dict]:
        """Generate step-by-step optimization plan."""
        if current_rs >= target_rs:
            return [{"step": 1, "action": "✓ Target achieved", "detail": f"Rs={current_rs:.2f} ≥ {target_rs}", "category": "none"}]

        ratio_needed = (target_rs / max(current_rs, 0.01)) ** 2  # squared since Rs ∝ √N
        plan = []
        step = 1

        # Strategy 1: Increase N (longer column, smaller particles)
        L = kwargs.get("column_length_mm")
        dp = kwargs.get("particle_size_um")
        if mode == "fundamental":
            N_current = kwargs.get("N", 5000)
            N_needed = N_current * ratio_needed

            if L and dp:
                # Option A: Longer column
                L_new = min(L * ratio_needed, 500)
                N_from_L = N_current * (L_new / L)
                rs_from_L = current_rs * math.sqrt(L_new / L)

                # Option B: Smaller particles
                dp_options = [5.0, 3.5, 2.6, 2.1, 1.7, 1.3]
                best_dp = dp
                best_rs_dp = current_rs
                for dp_opt in dp_options:
                    if dp_opt < dp:
                        n_ratio = dp / dp_opt  # N inversely proportional to dp
                        rs_new = current_rs * math.sqrt(n_ratio)
                        if rs_new > best_rs_dp:
                            best_rs_dp = rs_new
                            best_dp = dp_opt

                plan.append({
                    "step": step, "category": "efficiency",
                    "action": f"Increase column length to {int(L_new)} mm",
                    "detail": f"Estimated Rs ≈ {rs_from_L:.2f}. Current: {L} mm.",
                    "feasibility": "easy" if L_new <= 250 else "moderate",
                })
                step += 1

                if best_dp < dp:
                    plan.append({
                        "step": step, "category": "efficiency",
                        "action": f"Switch to smaller particles ({dp} → {best_dp} μm)",
                        "detail": f"Estimated Rs ≈ {best_rs_dp:.2f}. May require UHPLC system.",
                        "feasibility": "moderate" if best_dp <= 2.6 else "advanced",
                    })
                    step += 1

        # Strategy 2: Improve selectivity (most powerful lever)
        alpha_current = kwargs.get("alpha", 1.05)
        if mode == "fundamental" and alpha_current < 1.5:
            # How much α improvement needed?
            # Rs ∝ (α-1)/α, so we need new_α such that (new_α-1)/new_α = target/current * (α-1)/α
            current_sel = (alpha_current - 1) / alpha_current
            needed_sel = current_sel * (target_rs / max(current_rs, 0.01))
            # Solve: (α'-1)/α' = needed_sel → α' = 1/(1-needed_sel)
            if needed_sel < 1:
                alpha_needed = 1 / (1 - needed_sel)
                plan.append({
                    "step": step, "category": "selectivity",
                    "action": f"Improve selectivity (α: {alpha_current:.3f} → {alpha_needed:.3f})",
                    "detail": (
                        "Try: different stationary phase (C8→Phenyl→PFP), "
                        "change organic modifier (ACN↔MeOH↔THF), "
                        "adjust pH ±1-2 units, change temperature ±15°C"
                    ),
                    "feasibility": "high_impact",
                    "note": "Selectivity is the most powerful lever — small α changes give large Rs gains",
                })
                step += 1

        # Strategy 3: Optimize retention (k range)
        plan.append({
            "step": step, "category": "retention",
            "action": "Optimize retention factors into k=2-10 range",
            "detail": "Weaken initial %B if k<2; strengthen if k>10. Gradient: adjust start/end %B.",
            "feasibility": "easy",
        })

        return plan

    def _run_base(self, calculation_mode: str = "direct",
                  tR1_min: Optional[float] = None, tR2_min: Optional[float] = None,
                  W1_min: Optional[float] = None, W2_min: Optional[float] = None,
                  Wh1_min: Optional[float] = None, Wh2_min: Optional[float] = None,
                  N: Optional[float] = None, alpha: Optional[float] = None,
                  k1: Optional[float] = None, k2: Optional[float] = None,
                  target_rs: float = 1.5,
                  column_length_mm: Optional[float] = None,
                  particle_size_um: Optional[float] = None) -> dict:

        # Calculate Rs based on mode
        if calculation_mode == "fundamental" and N is not None and alpha is not None:
            rs_result = self._rs_fundamental(N, alpha, k1, k2)
        elif calculation_mode in ("direct", "predict") and tR1_min and tR2_min:
            if W1_min and W2_min:
                rs_result = self._rs_direct(tR1_min, tR2_min, W1_min, W2_min)
            elif Wh1_min and Wh2_min:
                rs_result = self._rs_wh(tR1_min, tR2_min, Wh1_min, Wh2_min)
            else:
                raise ChemMCPError("For direct mode, need either baseline widths (W1, W2) or half-height widths (Wh1, Wh2)")
        else:
            raise ChemMCPError(
                f"Insufficient parameters for '{calculation_mode}' mode. "
                f"For 'direct': need tR1, tR2, and widths. "
                f"For 'fundamental': need N, alpha, and at least one k value."
            )

        rs_val = rs_result["value"]

        # Assessment
        baseline_sep = rs_val >= 1.5
        assessment = {
            "rs_value": rs_val,
            "baseline_separated": baseline_sep,
            "peak_purity_estimate": min(99.9, 50 + 50 * math.tanh((rs_val - 0.8) * 2)),
            "grade": (
                "Excellent (Rs ≥ 2.0)" if rs_val >= 2.0 else
                "Baseline (1.5 ≤ Rs < 2.0)" if rs_val >= 1.5 else
                "Partial (1.0 ≤ Rs < 1.5)" if rs_val >= 1.0 else
                "Poor co-elution (Rs < 1.0)"
            ),
        }

        # Optimization plan
        opt_plan = self._optimization_plan(
            rs_val, target_rs, calculation_mode,
            column_length_mm=column_length_mm,
            particle_size_um=particle_size_um,
            N=N, alpha=alpha,
        )

        result = {
            "resolution": {
                "result": rs_result,
                "assessment": assessment,
                "target_resolution": target_rs,
                "gap_to_target": round(max(0, target_rs - rs_val), 3),
                "optimization_plan": opt_plan,
                "quick_reference": {
                    "Rs < 0.5": "Severe co-elution — completely overlapping peaks",
                    "Rs = 0.5-1.0": "Partial separation — shoulders visible",
                    "Rs = 1.0-1.5": "Partial baseline — ~94% purity at Rs=1.0, ~98.7% at Rs=1.5",
                    "Rs = 1.5-2.0": "Baseline separation — standard acceptance criterion",
                    "Rs > 2.0": "Excellent separation — comfortable margin",
                },
            }
        }
        return result

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        for part in input_params.strip().split():
            if "=" in part:
                key, val = part.split("=", 1)
                try:
                    kwargs[key] = float(val)
                except ValueError:
                    kwargs[key] = val
        return self._run_base(**kwargs)
