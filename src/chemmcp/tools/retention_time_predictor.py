import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ─── Solvent Strength Parameters (Snyder) ───
_SNYDER_P = {
    "acetonitrile": 0.652,
    "methanol": 0.729,
    "ethanol": 0.688,
    "tetrahydrofuran": 0.450,
    "isopropanol": 0.606,
    "dioxane": 0.561,
}

# Approximate linear solvent strength relationship parameters
# log(k) = log(k_w) - S * φ, where φ is volume fraction of organic modifier
# Typical S values for small molecules on C18
_S_VALUES = {
    "very_small_nonpolar": 2.5,   # MW < 150, non-polar
    "small_moderate": 3.5,         # MW 150-300, moderate polarity
    "medium_polar": 4.5,           # MW 200-500, polar
    "large_polar": 5.5,            # MW > 400, polar/ionic
    "peptide": 8.0,                # Peptides
}


@ChemMCPManager.register_tool
class RetentionTimePredictor(BaseTool):
    """
    保留时间预测与方法转移工具。
    基于线性溶剂强度(LSS)模型预测HPLC保留时间，支持不同色谱柱和仪器间的方法转移。
    """
    __version__ = "0.1.0"
    name = "RetentionTimePredictor"
    func_name = "predict_retention_time"
    description = "Predict HPLC retention times using Linear Solvent Strength (LSS) model and perform method transfer between different columns/instruments."
    implementation_description = (
        "Implements the Linear Solvent Strength (LSS) model: log(k) = log(kw) - S·φ, "
        "where k is retention factor, φ is organic fraction, and S is solvent strength parameter. "
        "Supports method transfer via column geometry scaling and gradient compression factors."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["HPLC", "Retention Time", "Method Transfer", "LSS Model", "Chromatography"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "predict", "Mode: 'predict' (new prediction), 'transfer' (method transfer), or 'calibrate' (estimate k_w and S)."),
        ("analyte_type", "str", "small_moderate", "Analyte category for S value estimation."),
        ("column_length_mm", "float", "100", "Column length in mm."),
        ("column_id_mm", "float", "2.1", "Column inner diameter in mm."),
        ("flow_rate_ml_min", "float", "0.3", "Flow rate in mL/min."),
        ("organic_start_pct", "float", "5.0", "Initial organic percentage in gradient."),
        ("organic_end_pct", "float", "95.0", "Final organic percentage."),
        ("gradient_time_min", "float", "10.0", "Gradient time in minutes."),
        ("dead_time_min", "None", "None", "Column dead time t₀ (optional; calculated if not provided)."),
        ("known_retention_min", "None", "None", "Known retention time at known conditions (for calibration/transfer)."),
        ("target_column_length_mm", "None", "None", "Target column length for method transfer (mm)."),
        ("target_column_id_mm", "None", "None", "Target column ID for method transfer (mm)."),
        ("target_flow_rate_ml_min", "None", "None", "Target flow rate for method transfer (mL/min)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Description of retention problem. Example: 'predict retention for drug-like compound on 100x2.1mm column 0.3 mL/min 5-95%ACN in 10min'"),
    ]

    output_sig = [
        ("prediction", "dict", "Retention time prediction with k values, method transfer parameters, and scaled conditions if applicable."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "predict",
                "analyte_type": "small_moderate",
                "column_length_mm": 100,
                "column_id_mm": 2.1,
                "flow_rate_ml_min": 0.3,
                "organic_start_pct": 5,
                "organic_end_pct": 95,
                "gradient_time_min": 10,
            },
            "text_input": {"input_params": "predict retention drug-like compound 100x2.1mm 0.3mL/min 5-95% ACN 10min"},
            "output": {"prediction": {"retention_time_min": 7.23, "k_factor": 4.2}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_dead_time(self, length_mm: float, id_mm: float, flow_rate: float) -> float:
        """Estimate column dead time from geometry."""
        # Volumetric dead time = column volume / flow rate
        vol_ml = math.pi * (id_mm / 2) ** 2 * length_mm / 1000  # mm³ → mL
        total_porosity = 0.68  # typical for fully porous particles
        return vol_ml * total_porosity / flow_rate

    def _predict_isocratic_k(self, phi: float, kw: float, S: float) -> float:
        """Calculate retention factor at given organic fraction."""
        log_k = math.log10(kw) - S * phi
        return max(10 ** log_k, 0.01)

    def _predict_gradient_retention(self, t0: float, kw: float, S: float,
                                     phi_start: float, phi_end: float,
                                     tg: float) -> dict:
        """Predict retention time under gradient conditions (simplified LSS)."""
        beta = (phi_end - phi_start) * t0 / tg  # Gradient steepness

        log_k_end = math.log10(max(kw, 1.01)) - S * phi_end
        k_end = 10 ** log_k_end

        if k_end < 1e-3:
            # Elutes before gradient effectively starts
            return {"tR_min": t0 * 1.05, "k": 0.05, "note": "Elutes near void"}

        # Simplified gradient retention approximation
        # For linear gradient: average k ≈ (1/b) * ln(2.3 * kw * S * b + 1)
        b = S * beta
        if abs(b) < 0.001:
            # Nearly isocratic
            phi_avg = (phi_start + phi_end) / 2
            k_avg = self._predict_isocratic_k(phi_avg, kw, S)
            tR = t0 * (1 + k_avg)
        else:
            try:
                k_avg = (1 / b) * math.log(2.3 * kw * S * b * (t0 / (t0 * (1 + kw ** 0.5))) + 1)
                k_avg = max(min(k_avg, 500), 0.02)
            except (ValueError, OverflowError):
                k_avg = kw ** 0.5
            tR = t0 * (1 + k_avg)

        # More practical estimation based on empirical rules
        # Compounds typically elute between 20%-80% of gradient time for well-designed methods
        phi_elution = (math.log10(max(kw, 1.1)) / S) if S > 0 else 0.5
        phi_elution = max(0, min(1, phi_elution))

        if phi_start <= phi_elution <= phi_end:
            frac = (phi_elution - phi_start) / (phi_end - phi_start)
            tR_est = t0 + frac * tg
        elif phi_elution < phi_start:
            tR_est = t0 * (1 + self._predict_isocratic_k(phi_start, kw, S))
        else:
            tR_est = t0 + tg * 0.9  # Elutes late in gradient

        tR_est = max(t0 * 1.02, min(tR_est, t0 + tg + 2))
        k_final = (tR_est - t0) / t0

        return {
            "tR_min": round(tR_est, 2),
            "k": round(k_final, 2),
            "estimated_elution_phi": round(phi_elution, 3),
            "elution_position_pct_gradient": round(frac * 100 if phi_start <= phi_elution <= phi_end else 90, 1),
        }

    def _method_transfer(self, tR_known: float, t0_source: float, t0_target: float,
                         L_source: float, L_target: float,
                         id_source: float, id_target: float,
                         F_source: float, F_target: float,
                         tg_source: float) -> dict:
        """Calculate method transfer scaling factors."""
        # Geometric scaling factor
        geom_scale = (L_target / L_source) * (id_source / id_target) ** 2 * (F_source / F_target)

        # Gradient compression/extension factor
        t0_ratio = t0_target / t0_source
        grad_scale = t0_ratio  # Scale gradient time proportionally to t0 change

        tR_predicted = t0_target + (tR_known - t0_source) * geom_scale
        tg_new = tg_source * grad_scale

        return {
            "source_conditions": {
                "t0_min": round(t0_source, 3),
                "L_mm": L_source,
                "ID_mm": id_source,
                "F_ml_min": F_source,
            },
            "target_conditions": {
                "t0_min": round(t0_target, 3),
                "L_mm": L_target,
                "ID_mm": id_target,
                "F_ml_min": F_target,
            },
            "scaling_factors": {
                "geometric_scale_factor": round(geom_scale, 3),
                "gradient_scale_factor": round(grad_scale, 3),
            },
            "predicted_retention_target_min": round(tR_predicted, 2),
            "recommended_gradient_time_min": round(tg_new, 2),
            "transfer_quality_note": (
                "Good transfer" if 0.8 < geom_scale < 1.5 else
                "Significant scale change — verify selectivity"
            ),
        }

    def _run_base(self, mode: str = "predict", analyte_type: str = "small_moderate",
                  column_length_mm: float = 100, column_id_mm: float = 2.1,
                  flow_rate_ml_min: float = 0.3,
                  organic_start_pct: float = 5.0, organic_end_pct: float = 95.0,
                  gradient_time_min: float = 10.0,
                  dead_time_min: Optional[float] = None,
                  known_retention_min: Optional[float] = None,
                  target_column_length_mm: Optional[float] = None,
                  target_column_id_mm: Optional[float] = None,
                  target_flow_rate_ml_min: Optional[float] = None) -> dict:
        """Core logic."""

        # Dead time
        t0 = dead_time_min if dead_time_min else self._calc_dead_time(column_length_mm, column_id_mm, flow_rate_ml_min)

        # Get S and estimate kw
        S = _S_VALUES.get(analyte_type, _S_VALUES["small_moderate"])
        kw = 10 ** (S * 0.3)  # Assume k≈5 at 30% B as reference point

        phi_s = organic_start_pct / 100.0
        phi_e = organic_end_pct / 100.0

        result_data = {}

        if mode == "transfer" and known_retention_min and target_column_length_mm:
            t0_target = self._calc_dead_time(
                target_column_length_mm,
                target_column_id_mm or column_id_mm,
                target_flow_rate_ml_min or flow_rate_ml_min
            )
            transfer_result = self._method_transfer(
                known_retention_min, t0, t0_target,
                column_length_mm, target_column_length_mm,
                column_id_mm, target_column_id_mm or column_id_mm,
                flow_rate_ml_min, target_flow_rate_ml_min or flow_rate_ml_min,
                gradient_time_min
            )
            result_data["transfer"] = transfer_result
            result_data["mode"] = "method_transfer"

        elif mode == "calibrate":
            if known_retention_min:
                k_observed = (known_retention_min - t0) / t0
                phi_avg = (phi_s + phi_e) / 2
                # Back-calculate kw: log(kw) = log(k) + S*φ
                log_kw = math.log10(max(k_observed, 0.01)) + S * phi_avg
                kw_calc = 10 ** log_kw
                result_data["calibration"] = {
                    "observed_tR_min": known_retention_min,
                    "observed_k": round(k_observed, 3),
                    "estimated_kw": round(kw_calc, 2),
                    "S_value_used": S,
                    "note": f"Based on observed k={k_observed:.2f} at avg φ={phi_avg:.2f}",
                }
            result_data["mode"] = "calibration"
            result_data["dead_time"] = round(t0, 3)

        else:  # predict
            pred = self._predict_gradient_retention(t0, kw, S, phi_s, phi_e, gradient_time_min)

            # Also predict at nearby conditions for robustness check
            pred_early = self._predict_gradient_retention(t0, kw, S, phi_s, phi_e, gradient_time_min * 0.8)
            pred_late = self._predict_gradient_retention(t0, kw, S, phi_s, phi_e, gradient_time_min * 1.2)

            result_data["prediction"] = {
                "primary": pred,
                "robustness_check": {
                    "if_gradient_80pct": pred_early,
                    "if_gradient_120pct": pred_late,
                },
                "model_parameters": {
                    "kw_estimate": round(kw, 2),
                    "S_value": S,
                    "analyte_type": analyte_type,
                    "dead_time_t0_min": round(t0, 3),
                },
                "interpretation": self._interpret_retention(pred["tR_min"], t0, gradient_time_min),
            }
            result_data["mode"] = "prediction"
            result_data["column_info"] = {
                "length_mm": column_length_mm,
                "inner_diameter_mm": column_id_mm,
                "flow_rate_ml_min": flow_rate_ml_min,
                "linear_velocity_mm_min": round(column_length_mm / t0, 1) if t0 > 0 else 0,
            }

        return {"prediction": result_data}

    def _interpret_retention(self, tR: float, t0: float, tg: float) -> str:
        k = (tR - t0) / t0
        pos = (tR - t0) / tg * 100 if tg > 0 else 0
        if k < 0.5:
            return f"Very early eluting (k={k:.1f}). Consider starting at lower %B or using weaker initial mobile phase."
        elif k < 2:
            return f"Early-mid eluting (k={k:.1f}). Acceptable range but could increase retention for better resolution."
        elif k < 15:
            return f"Good retention range (k={k:.1f}). Well-positioned within the gradient ({pos:.0f}%)."
        elif k < 30:
            return f"Late eluting (k={k:.1f}). Consider steeper gradient or higher %B start."
        else:
            return f"Very strongly retained (k={k:.1f}). May need stronger eluent or longer gradient."

    def _run_text(self, input_params: str) -> dict:
        text = input_params.lower()
        mode = "predict"
        if "transfer" in text:
            mode = "transfer"
        elif "calibrate" in text:
            mode = "calibrate"

        # Extract numbers
        import re
        numbers = [float(x) for x in re.findall(r'[\d.]+', text)]

        atype = "small_moderate"
        if any(w in text for w in ["peptide", "protein"]):
            atype = "peptide"
        elif any(w in text for w in ["large", "high_mw"]):
            atype = "large_polar"
        elif "non.polar" in text or "hydrophob" in text:
            atype = "very_small_nonpolar"

        kwargs = {"mode": mode, "analyte_type": atype}
        if len(numbers) >= 2:
            kwargs["column_length_mm"] = numbers[0]
            kwargs["column_id_mm"] = numbers[1]
        if len(numbers) >= 3:
            kwargs["flow_rate_ml_min"] = numbers[2]
        if len(numbers) >= 5:
            kwargs["organic_start_pct"] = numbers[3]
            kwargs["organic_end_pct"] = numbers[4]
        if len(numbers) >= 6:
            kwargs["gradient_time_min"] = numbers[5]

        return self._run_base(**kwargs)
