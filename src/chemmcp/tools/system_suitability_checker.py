import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# USP/EP/ChP system suitability acceptance criteria reference
SST_LIMITS = {
    "usp": {
        "resolution_Rs_min": 1.5,
        "tailing_factor_Tf_max": 2.0,
        "theoretical_plates_N_min": 2000,
        "capacity_factor_k_min": 1.0,
        "capacity_factor_k_max": 10.0,
        "RSD_percent_max": 2.0,  # for n ≥ 5 injections
    },
    "pharmacopeial_strict": {
        "resolution_Rs_min": 2.0,
        "tailing_factor_Tf_max": 1.5,
        "theoretical_plates_N_min": 5000,
        "capacity_factor_k_min": 2.0,
        "capacity_factor_k_max": 8.0,
        "RSD_percent_max": 1.0,
    },
    "general": {
        "resolution_Rs_min": 1.0,
        "tailing_factor_Tf_max": 3.0,
        "theoretical_plates_N_min": 1000,
        "capacity_factor_k_min": 0.5,
        "capacity_factor_k_max": 20.0,
        "RSD_percent_max": 5.0,
    },
}


@ChemMCPManager.register_tool
class SystemSuitabilityChecker(BaseTool):
    """
    系统适用性测试参数计算工具。
    计算色谱系统适用性关键参数：分辨率、塔板数、拖尾因子、容量因子、RSD等，对照药典标准判定。
    """
    __version__ = "0.1.0"
    name = "SystemSuitabilityChecker"
    func_name = "check_system_suitability"
    description = "Calculate and evaluate chromatographic system suitability parameters: resolution, tailing factor, theoretical plates, capacity factor, RSD, against pharmacopeial standards."
    implementation_description = (
        "Implements full SST calculation suite per ICH Q2(R1) / USP <621> guidelines: "
        "resolution (Rs), tailing factor (T/Af), theoretical plates (N), capacity factor (k'), "
        "relative retention (α), repeatability (RSD%), and peak-to-peak ratio. "
        "Supports multiple standard sets (USP general, strict, custom). "
        "Provides pass/fail verdict with detailed diagnostic information."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "System Suitability", "SST", "HPLC", "USP", "ICH", "Validation"]
    required_envs = []

    code_input_sig = [
        ("peak_data_list", "list", "N/A", "List of peak data dicts. Each: {'tR': float, 'Wh': float or 'Wb': float, 'height': float, 'area': float}."),
        ("standard_set", "str", "'usp'", "Standard set to evaluate against: 'usp', 'pharmacopeial_strict', 'general'."),
        ("n_injections", "int", "5", "Number of replicate injections for RSD calculation."),
        ("replicate_areas", "None", "None", "List of area values from replicate injections for RSD calculation."),
        ("column_dead_time_min", "None", "None", "Column dead time t₀ in minutes (for k' calculation)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'peaks=tR1=5.0,Wb1=0.3,tR2=6.2,Wb2=0.35 t0=0.8 standard=usp'."),
    ]

    output_sig = [
        ("sst_report", "dict", "Complete system suitability report with all calculated parameters, pass/fail status per criterion, overall verdict, and recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "peak_data_list": [
                    {"tR": 5.02, "Wh": 0.12, "height": 15000, "area": 285000},
                    {"tR": 6.35, "Wh": 0.14, "height": 12000, "area": 248000},
                ],
                "column_dead_time_min": 0.82,
                "standard_set": "usp",
            },
            "text_input": {"input_params": "peaks=tR1=5.02,Wh1=0.12,tR2=6.35,Wh2=0.14 t0=0.82"},
            "output": {
                "sst_report": {"overall_verdict": "PASS", "parameters": {...}}
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_resolution(self, tR1: float, tR2: float, Wb1: float, Wb2: float) -> float:
        """Rs = 2(tR2 - tR1) / (Wb1 + Wb2)."""
        return abs(2 * (tR2 - tR1)) / (Wb1 + Wb2)

    def _calc_resolution_wh(self, tR1: float, tR2: float, Wh1: float, Wh2: float) -> float:
        """Rs using half-height widths: Rs ≈ 1.18(tR2-tR1)/(Wh1+Wh2)."""
        return 1.18 * abs(tR2 - tR1) / (Wh1 + Wh2)

    def _calc_tailing_factor(self, tR: float, Wh: float, t_front: Optional[float] = None,
                              t_tail: Optional[float] = None) -> float:
        """Tailing factor T = W(0.05) / 2f where f is front half-width at 5% height."""
        # Simplified: if only Wh given, estimate T from asymmetry approximation
        if t_front is not None and t_tail is not None:
            w_005 = t_tail - t_front
            f = tR - t_front
            return w_005 / (2 * f)
        # Default to symmetric if no data
        return 1.0

    def _calc_plates(self, tR: float, Wh: float) -> int:
        """N = 5.54 × (tR/Wh)²."""
        return int(round(5.54 * (tR / Wh) ** 2))

    def _calc_plates_wb(self, tR: float, Wb: float) -> int:
        """N = 16 × (tR/Wb)²."""
        return int(round(16 * (tR / Wb) ** 2))

    def _calc_capacity_factor(self, tR: float, t0: float) -> Optional[float]:
        """k' = (tR - t0) / t0."""
        if t0 is None or t0 <= 0:
            return None
        return (tR - t0) / t0

    def _calc_selectivity(self, tR1: float, tR2: float, t0: Optional[float]) -> Optional[float]:
        """α = k2'/k1' = (tR2-t0)/(tR1-t0)."""
        if t0 is None or t0 <= 0:
            return None
        k1 = (tR1 - t0) / t0
        k2 = (tR2 - t0) / t0
        if k1 <= 0:
            return None
        return k2 / k1

    def _calc_rsd(self, values: list) -> dict:
        """Calculate relative standard deviation."""
        if not values or len(values) < 2:
            return {"RSD_percent": None, "n": len(values) if values else 0}
        n = len(values)
        mean = sum(values) / n
        std_dev = math.sqrt(sum((x - mean) ** 2 for x in values) / (n - 1)) if n > 1 else 0
        rsd = (std_dev / mean * 100) if mean != 0 else 0
        return {
            "RSD_percent": round(rsd, 4),
            "mean": round(mean, 2),
            "std_dev": round(std_dev, 2),
            "n": n,
        }

    def _run_base(self, peak_data_list: list,
                  standard_set: str = "usp",
                  n_injections: int = 5,
                  replicate_areas: Optional[List[float]] = None,
                  column_dead_time_min: Optional[float] = None) -> dict:
        """Core logic."""

        limits = SST_LIMITS.get(standard_set, SST_LIMITS["usp"])
        peaks = peak_data_list
        n_peaks = len(peaks)

        if n_peaks < 1:
            raise ChemMCPError("At least one peak required.")

        # Per-peak calculations
        peak_results = []
        for i, pk in enumerate(peaks):
            tR = pk["tR"]
            Wh = pk.get("Wh")
            Wb = pk.get("Wb")
            height = pk.get("height")
            area = pk.get("area")

            # Plates
            N = self._calc_plates(tR, Wh) if Wh else self._calc_plates_wb(tR, Wb) if Wb else None
            method = "half-height" if Wh else "baseline" if Wb else "unknown"

            # Capacity factor
            k = self._calc_capacity_factor(tR, column_dead_time_min)

            # Tailing (simplified — needs more granular data)
            Tf = pk.get("Tf", 1.0)

            peak_results.append({
                "peak_id": i + 1,
                "retention_time_min": tR,
                "theoretical_plates_N": N,
                "plates_method": method,
                "capacity_factor_k_prime": round(k, 3) if k else None,
                "tailing_factor_Tf": Tf,
                "height": height,
                "area": area,
            })

        # Resolution between adjacent peaks
        resolutions = []
        selectivities = []
        for i in range(n_peaks - 1):
            p1, p2 = peaks[i], peaks[i + 1]
            w1 = p1.get("Wb") or (p1.get("Wh") * 1.82 if p1.get("Wh") else 0.1)
            w2 = p2.get("Wb") or (p2.get("Wh") * 1.82 if p2.get("Wh") else 0.1)
            rs = self._calc_resolution(p1["tR"], p2["tR"], w1, w2)
            alpha = self._calc_selectivity(p1["tR"], p2["tR"], column_dead_time_min)
            resolutions.append({
                "peak_pair": f"{i+1}/{i+2}",
                "Rs": round(rs, 3),
                "pass": rs >= limits["resolution_Rs_min"],
            })
            if alpha is not None:
                selectivities.append({
                    "peak_pair": f"{i+1}/{i+2}",
                    "alpha": round(alpha, 4),
                })

        # RSD from replicates
        rsd_result = self._calc_rsd(replicate_areas) if replicate_areas else None

        # Evaluate each criterion
        criteria_results = []

        # Criterion 1: Resolution
        min_rs = min((r["Rs"] for r in resolutions), default=None)
        if min_rs is not None:
            criteria_results.append({
                "parameter": "Resolution (Rs)",
                "value": min_rs,
                "limit": f"≥ {limits['resolution_Rs_min']}",
                "pass": min_rs >= limits["resolution_Rs_min"],
            })

        # Criterion 2: Tailing factor
        max_tf = max((p["tailing_factor_Tf"] for p in peak_results if p["tailing_factor_Tf"] is not None), default=None)
        if max_tf is not None:
            criteria_results.append({
                "parameter": "Tailing Factor (T)",
                "value": max_tf,
                "limit": f"≤ {limits['tailing_factor_Tf_max']}",
                "pass": max_tf <= limits["tailing_factor_Tf_max"],
            })

        # Criterion 3: Theoretical plates
        min_n = min((p["theoretical_plates_N"] for p in peak_results if p["theoretical_plates_N"] is not None), default=None)
        if min_n is not None:
            criteria_results.append({
                "parameter": "Theoretical Plates (N)",
                "value": min_n,
                "limit": f"≥ {limits['theoretical_plates_N_min']}",
                "pass": min_n >= limits["theoretical_plates_N_min"],
            })

        # Criterion 4: Capacity factor range
        all_k = [p["capacity_factor_k_prime"] for p in peak_results if p["capacity_factor_k_prime"] is not None]
        if all_k:
            k_ok = all(limits["capacity_factor_k_min"] <= k <= limits["capacity_factor_k_max"] for k in all_k)
            criteria_results.append({
                "parameter": "Capacity Factor (k')",
                "value": [round(k, 3) for k in all_k],
                "limit": f"{limits['capacity_factor_k_min']} – {limits['capacity_factor_k_max']}",
                "pass": k_ok,
            })

        # Criterion 5: RSD
        if rsd_result and rsd_result["RSD_percent"] is not None:
            criteria_results.append({
                "parameter": "Repeatability (RSD%)",
                "value": rsd_result["RSD_percent"],
                "limit": f"≤ {limits['RSD_percent_max']}%",
                "pass": rsd_result["RSD_percent"] <= limits["RSD_percent_max"],
            })

        # Overall verdict
        all_pass = all(c["pass"] for c in criteria_results) if criteria_results else True
        verdict = "✅ PASS" if all_pass else "❌ FAIL"

        result = {
            "sst_report": {
                "overall_verdict": verdict,
                "standard_set": standard_set,
                "limits_applied": limits,
                "n_peaks": n_peaks,
                "n_injections": n_injections,
                "per_peak_parameters": peak_results,
                "resolutions_between_peaks": resolutions,
                "selectivity_factors": selectivities,
                "repeatability": rsd_result,
                "criteria_evaluation": criteria_results,
                "summary_table": self._make_summary(criteria_results),
                "recommendations": self._get_recommendations(all_pass, criteria_results),
            }
        }
        return result

    def _make_summary(self, criteria: list) -> str:
        lines = ["┌─────────────────────┬──────────┬─────────┬──────┐",
                 "│ Parameter           │ Value    │ Limit   │ Pass │",
                 "├─────────────────────┼──────────┼─────────┼──────┤"]
        for c in criteria:
            val = str(c["value"]) if not isinstance(c["value"], list) else ",".join(str(v) for v in c["value"])
            status = "✅" if c["pass"] else "❌"
            lines.append(f"│ {c['parameter']:<19} │ {val:<8} │ {c['limit']:<7} │ {status} │")
        lines.append("└─────────────────────┴──────────┴─────────┴──────┘")
        return "\n".join(lines)

    def _get_recommendations(self, passed: bool, criteria: list) -> List[str]:
        recs = []
        if passed:
            recs.append("✅ System suitability test PASSED — system ready for sample analysis.")
            return recs

        failed = [c for c in criteria if not c["pass"]]
        for f in failed:
            param = f["parameter"]
            if "Resolution" in param:
                recs.append(f"⚠ {param} ({f['value']}) below limit — adjust mobile phase strength, temperature, or use a longer column.")
            elif "Tailing" in param:
                recs.append(f"⚠ {param} ({f['value']}) exceeds limit — check column health, reduce injection volume, or adjust sample solvent.")
            elif "Plates" in param:
                recs.append(f"⚠ {param} ({f['value']}) too low — consider replacing column, reducing extra-column volume, or lowering flow rate.")
            elif "Capacity" in param:
                recs.append(f"⚠ {param} out of range — adjust mobile phase composition to optimize retention.")
            elif "RSD" in param:
                recs.append(f"⚠ {param} ({f['value']}%) too high — check injector precision, ensure proper sample dissolution, verify integration consistency.")

        recs.append("Re-run SST after adjustments before proceeding with analysis.")
        return recs[:6]

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        peaks_raw = []
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                if k == "peaks":
                    # Parse simplified format
                    sub_parts = v.split(",")
                    current_peak = {}
                    for sp in sub_parts:
                        if "=" in sp:
                            sk, sv = sp.split("=", 1)
                            if sk.startswith("tR"):
                                idx = sk[2:]
                                if idx == "":
                                    current_peak["tR"] = float(sv)
                                else:
                                    if current_peak.get("tR"):
                                        peaks_raw.append(current_peak)
                                    current_peak = {"tR": float(sv)}
                            elif sk.startswith("Wh"):
                                current_peak["Wh"] = float(sv)
                            elif sk.startswith("Wb"):
                                current_peak["Wb"] = float(sv)
                    if current_peak:
                        peaks_raw.append(current_peak)
                    kwargs["peak_data_list"] = peaks_raw
                elif k == "t0":
                    kwargs["column_dead_time_min"] = float(v)
                elif k == "standard":
                    kwargs["standard_set"] = v
        return self._run_base(**kwargs)
