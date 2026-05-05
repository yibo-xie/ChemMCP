"""
Recovery Calculator — 加标回收率计算与评估
支持单值与平行样，自定义判定区间
"""
import logging
import math
from typing import List, Optional, Dict, Any, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RecoveryCalculator(BaseTool):
    """
    加标回收率计算器：计算单次或多次平行测定的回收率，
    统计平均值、RSD，并根据自定义判定区间进行放行判定。
    """
    __version__ = "0.1.0"
    name = "RecoveryCalculator"
    func_name = "calculate_recovery"
    description = "Calculate spike recovery rate with statistical evaluation (mean, RSD, pass/fail assessment). Supports single values and replicate measurements."
    implementation_description = "Computes individual recoveries as (measured - baseline) / spike_added × 100%, then calculates mean recovery, RSD%, and compares against user-defined acceptance criteria for method validation/QC release."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Recovery", "Spike Recovery", "Method Validation", "QC", "Quality Control", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("measured_values", "float | list[float]", "required", "Measured value(s) after spiking. Single float or list of replicates."),
        ("spike_added", "float", "required", "Amount of spike added (same units as measured)."),
        ("baseline", "float", "0.0", "Background/baseline value in unspiked sample."),
        ("acceptance_min", "float", "80.0", "Lower acceptance bound (%). Recovery below this = FAIL."),
        ("acceptance_max", "float", "120.0", "Upper acceptance bound (%). Recovery above this = FAIL."),
        ("unit", "str", "", "Unit of measurement (for reporting; e.g., 'mg/L', 'μg/g', 'ppb')."),
        ("sample_id", "str", "", "Sample identifier for report labeling."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., '95.2 98.7 102.3 100.0' or single: '95.5 100 0 80 120'."),
    ]

    output_sig = [
        ("result", "dict", "JSON containing: recoveries[], mean_recovery, rsd_percent, acceptance_range, is_pass, assessment, and detailed statistics."),
    ]

    examples = [
        {
            "code_input": {
                "measured_values": [95.2, 98.7, 102.3],
                "spike_added": 100.0,
                "baseline": 0.0,
                "acceptance_min": 80.0,
                "acceptance_max": 120.0,
                "unit": "μg/L",
                "sample_id": "QC-2026-001",
            },
            "text_input": {
                "input_params": "95.2 98.7 102.3 100.0",
            },
            "output": {
                "result": {
                    "mode": "recovery_assessment",
                    "note": "Example output structure.",
                    "is_pass": True,
                }
            }
        },
        {
            "code_input": {
                "measured_values": 45.0,
                "spike_added": 100.0,
                "baseline": 5.0,
                "acceptance_min": 90.0,
                "acceptance_max": 110.0,
                "unit": "mg/kg",
                "sample_id": "SOIL-0042",
            },
            "text_input": {
                "input_params": "45.0 100.0 5.0 90 110",
            },
            "output": {
                "result": {
                    "mode": "recovery_assessment",
                    "note": "Single value with baseline subtraction.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _to_list(values: Union[float, int, List]) -> List[float]:
        if isinstance(values, (list, tuple)):
            return [float(v) for v in values]
        return [float(values)]

    @staticmethod
    def _calc_rsd(values: List[float]) -> float:
        """Calculate relative standard deviation (%) from a list of values."""
        n = len(values)
        if n < 2:
            return 0.0
        mean_val = sum(values) / n
        if mean_val == 0:
            return 0.0
        variance = sum((v - mean_val) ** 2 for v in values) / (n - 1)
        std_dev = math.sqrt(variance)
        return (std_dev / abs(mean_val)) * 100

    @staticmethod
    def _calc_confidence_interval(values: List[float], confidence: float = 0.95) -> tuple:
        """Calculate confidence interval for the mean using t-distribution approximation."""
        import math
        n = len(values)
        if n < 2:
            return (values[0], values[0])
        mean_val = sum(values) / n
        std_dev = math.sqrt(sum((v - mean_val) ** 2 for v in values) / (n - 1))
        # Approximate t-value (for 95% CI, reasonable approximations)
        t_approx = {2: 12.71, 3: 4.30, 4: 3.18, 5: 2.78, 6: 2.57,
                     7: 2.45, 8: 2.36, 9: 2.31, 10: 2.26}
        t = t_approx.get(n, 1.96)  # fallback to z-value for large n
        margin = t * std_dev / math.sqrt(n)
        return (round(mean_val - margin, 3), round(mean_val + margin, 3))

    def _assess_precision(self, rsd: float) -> dict:
        """Assess RSD against typical analytical chemistry benchmarks."""
        if rsd <= 2:
            level = "Excellent"
            verdict = "✓ Precision is excellent — suitable for regulatory analysis."
        elif rsd <= 5:
            level = "Good"
            verdict = "✓ Good precision — acceptable for most applications."
        elif rsd <= 10:
            level = "Acceptable"
            verdict = "~ Acceptable precision — may need investigation at trace levels."
        elif rsd <= 15:
            level = "Marginal"
            verdict = "⚠ Marginal precision — review method; consider additional replicates."
        else:
            level = "Poor"
            verdict = "✗ Poor precision — significant variability; investigate source."

        return {"rsd_percent": round(rsd, 3), "level": level, "verdict": verdict}

    def _run_base(self, measured_values: Union[float, List] = None,
                  spike_added: float = 0.0, baseline: float = 0.0,
                  acceptance_min: float = 80.0, acceptance_max: float = 120.0,
                  unit: str = "", sample_id: str = "") -> dict:

        # Validate inputs
        if measured_values is None:
            raise ChemMCPError("'measured_values' is required.")
        if spike_added <= 0:
            raise ChemMCPError(f"'spike_added' must be positive. Got: {spike_added}")

        values = self._to_list(measured_values)
        n = len(values)

        # Calculate individual recoveries
        recoveries = []
        for i, val in enumerate(values):
            net_measured = val - baseline
            rec = (net_measured / spike_added) * 100
            recoveries.append({
                "replicate": i + 1,
                "measured_value": val,
                "net_measured": round(net_measured, 4),
                "recovery_percent": round(rec, 3),
            })

        rec_values = [r["recovery_percent"] for r in recoveries]
        mean_rec = sum(rec_values) / n
        rsd = self._calc_rsd(rec_values)

        ci_low, ci_high = self._calc_confidence_interval(rec_values)

        # Pass/fail assessment
        all_in_range = all(acceptance_min <= r <= acceptance_max for r in rec_values)
        mean_in_range = acceptance_min <= mean_rec <= acceptance_max
        is_pass = all_in_range and mean_in_range

        # Detailed assessment
        low_recs = [r for r in rec_values if r < acceptance_min]
        high_recs = [r for r in rec_values if r > acceptance_max]

        if is_pass:
            assessment = f"PASS ✓ All {n} replicates within [{acceptance_min}%, {acceptance_max}%]. Mean={mean_rec:.1f}%, RSD={rsd:.2f}%."
        elif not mean_in_range:
            assessment = f"FAIL ✗ Mean recovery ({mean_rec:.1f}%) outside acceptance range [{acceptance_min}%, {acceptance_max}%]."
        else:
            outliers = len(low_recs) + len(high_recs)
            assessment = f"FAIL ✗ {outliers}/{n} replicate(s) outside range. Low: {len(low_recs)}, High: {len(high_recs)}."

        # Regulatory guidance reference
        guidance = self._get_guidance(mean_rec, rsd, n)

        return {"result": {
            "mode": "recovery_assessment",
            "sample_id": sample_id or "unspecified",
            "unit": unit or "not specified",
            "parameters": {
                "n_replicates": n,
                "spike_added": spike_added,
                "baseline": baseline,
                "acceptance_range_pct": [acceptance_min, acceptance_max],
            },
            "individual_recoveries": recoveries,
            "statistics": {
                "mean_recovery_percent": round(mean_rec, 3),
                "median_recovery_percent": round(sorted(rec_values)[n // 2], 3),
                "std_deviation_percent": round(math.sqrt(
                    sum((r - mean_rec) ** 2 for r in rec_values) / max(n - 1, 1)), 3),
                "rsd_percent": round(rsd, 3),
                "range_percent": [round(min(rec_values), 3), round(max(rec_values), 3)],
                "confidence_interval_95_pct": [ci_low, ci_high],
                "standard_error": round(
                    math.sqrt(sum((r - mean_rec) ** 2 for r in rec_values) / max(n - 1, 1))
                    / math.sqrt(n), 4),
            },
            "precision_assessment": self._assess_precision(rsd),
            "accuracy_assessment": {
                "bias_from_100_pct": round(mean_rec - 100, 3),
                "bias_interpretation": self._interpret_bias(mean_rec),
            },
            "acceptance": {
                "is_pass": is_pass,
                "all_replicates_in_range": all_in_range,
                "mean_in_range": mean_in_range,
                "acceptance_criteria": f"{acceptance_min}% ≤ Recovery ≤ {acceptance_max}%",
                "assessment": assessment,
            },
            "regulatory_guidance": guidance,
            "recommendations": self._make_recommendations(is_pass, mean_rec, rsd, n, low_recs, high_recs),
        }}

    @staticmethod
    def _interpret_bias(mean_rec: float) -> str:
        bias = mean_rec - 100
        abs_bias = abs(bias)
        if abs_bias <= 2:
            return "Negligible bias — excellent accuracy."
        elif abs_bias <= 5:
            return f"Slight {'positive' if bias > 0 else ''}bias ({bias:+.1f}%). Generally acceptable."
        elif abs_bias <= 10:
            return f"Moderate {'positive' if bias > 0 else ''}bias ({bias:+.1f}%). Review calibration/matrix effects."
        elif abs_bias <= 20:
            return f"Significant {'positive' if bias > 0 else ''}bias ({bias:+.1f}%). Investigation needed."
        else:
            return f"Severe bias ({bias:+.1f}%). Major systematic error — do NOT use data without correction."

    @staticmethod
    def _get_guidance(mean_rec: float, rsd: float, n: int) -> dict:
        """Reference to common regulatory guidelines."""
        guidelines = {
            "ICH Q2(R1)": {
                "typical_acceptance": "80-120% (may vary by concentration level)",
                "precision": "RSD ≤ 5% for repeatability (≤15% at LLOQ)",
                "note": "Pharmaceutical validation guideline.",
            },
            "EPA/ISO (environmental)": {
                "typical_acceptance": "70-130% (matrix-dependent)",
                "precision": "RSD ≤ 20% for environmental samples",
                "note": "Environmental method validation.",
            },
            "AOAC/SFDA (food)": {
                "typical_acceptance": "80-120%",
                "precision": "RSD varies by analyte/matrix",
                "note": "Food safety analysis standards.",
            },
            "General best practice": {
                "typical_acceptance": f"User-defined (current: custom range)",
                "precision": f"Current RSD={rsd:.2f}%",
                "note": f"n={n} replicates; more replicates → tighter CI.",
            },
        }
        return guidelines

    @staticmethod
    def _make_recommendations(is_pass: bool, mean_rec: float, rsd: float,
                              n: int, lows: list, highs: list) -> List[str]:
        recs = []
        if is_pass:
            recs.append("✓ Recovery test PASSED. Data can be used for reporting.")
            if rsd > 5:
                recs.append("→ Consider investigating source of variability for improved precision.")
        else:
            recs.append("✗ Recovery test FAILED. Do not use data until root cause identified.")

        if mean_rec < 85:
            recs.append("→ Low recovery: check extraction efficiency, analyte degradation, or matrix binding.")
        elif mean_rec > 115:
            recs.append("→ High recovery: check for contamination, interference, or over-correction.")

        if len(lows) > 0:
            recs.append(f"→ {len(lows)} replicate(s) below minimum — re-test or exclude with justification.")
        if len(highs) > 0:
            recs.append(f"→ {len(highs)} replicate(s) above maximum — check for pipetting error or contamination.")

        if n < 3:
            recs.append("→ Consider running ≥3 replicates for meaningful statistical evaluation.")

        recs.append("→ Document all results in lab notebook/LIMS with raw data attachment.")
        return recs

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) >= 4:
                # Format: value1 [value2 ...] spike_added [baseline [min max]]
                vals = [float(p) for p in parts[:-3]] if len(parts) > 4 else [float(parts[0])]
                spike = float(parts[-3])
                base = float(parts[-2]) if len(parts) > 2 else 0.0
                acc_min = float(parts[-1]) if len(parts) > 3 else 80.0
                acc_max = 120.0  # default
                if len(parts) > 5:
                    acc_max = float(parts[-1])
                    acc_min = float(parts[-2])
                return self._run_base(measured_values=vals, spike_added=spike,
                                       baseline=base, acceptance_min=acc_min,
                                       acceptance_max=acc_max)
            elif len(parts) == 2:
                return self._run_base(measured_values=float(parts[0]),
                                       spike_added=float(parts[1]))
            else:
                raise ValueError("Need at least: measured_value spike_added")
        except (IndexError, ValueError, ZeroDivisionError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
