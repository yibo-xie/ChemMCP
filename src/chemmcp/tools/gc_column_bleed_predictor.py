import logging
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Reference data: maximum isothermal operating temperatures for common GC stationary phases
# Source: Agilent/Restek/Supelec column specifications (publicly available)
STATIONARY_PHASE_TEMP_LIMITS = {
    # Non-polar / Low-polarity phases
    "DB-1": 260, "DB-1ms": 260, "HP-1": 260, "SE-30": 280,
    "DB-5": 325, "DB-5ms": 325, "HP-5": 325, "SE-54": 300,
    "DB-35": 340, "DB-35ms": 340, "HP-35": 340,
    "DB-624": 260, "DB-624ms": 260,
    "DB-1301": 290, "DB-1301ms": 290,

    # Mid-polarity phases
    "DB-17": 270, "DB-17ms": 270, "HP-50": 270, "OV-17": 275,
    "DB-1701": 270, "DB-1701ms": 270,
    "DB-225": 240, "DB-225ms": 240,
    "DB-WAX": 250, "DB-WAXms": 250, "HP-INNOWax": 250,
    "DB-FFAP": 250, "HP-FFAP": 245,

    # PEG / Wax phases
    "PEG-20M": 220, "Carbowax 20M": 220,
    "Stabilwax": 260,

    # Specialty columns
    "DB-23": 260,  # FAME analysis
    "DB-210": 260,
    "PLOT Q": 250,  # Poraplot Q
    "PLOT S": 200,  # Poraplot S

    # Generic categories (fallback)
    "polydimethylsiloxane": 260,
    "5%-phenyl-methylpolysiloxane": 325,
    "14%-cyanopropylphenyl-methylpolysiloxane": 270,
    "50%-phenyl-dimethylpolysiloxane": 370,
    "polyethylene_glycol": 250,
    "wax": 250,
    "ffap": 250,
}

# Bleed severity thresholds (relative to max temp)
BLEED_THRESHOLDS = {
    "safe": 0.85,      # < 85% of max temp — negligible bleed
    "moderate": 0.95,   # 85-95% — noticeable baseline rise
    "high": 1.0,        # 95-100% — significant bleed
    "critical": 1.05,   # > 100% — severe bleed risk, column damage likely
}


@ChemMCPManager.register_tool
class GcColumnBleedPredictor(BaseTool):
    """
    色谱柱流失温度限制提醒工具。
    根据GC柱固定相类型、柱温程序参数，预测柱流失风险，给出温度安全建议。
    """
    __version__ = "0.1.0"
    name = "GcColumnBleedPredictor"
    func_name = "predict_column_bleed"
    description = "Predict GC column stationary phase bleed risk based on column type, temperature program, and operating conditions."
    implementation_description = (
        "Uses reference temperature limits for common GC stationary phases to assess "
        "bleed risk at given oven temperatures. Provides safety margins, "
        "baseline drift warnings, and recommendations for temperature programming "
        "to minimize column degradation and maintain detection sensitivity."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "GC", "Column Bleed", "Temperature", "Stationary Phase"]
    required_envs = []

    code_input_sig = [
        ("stationary_phase", "str", "N/A", "Name of the stationary phase (e.g., 'DB-5', 'DB-WAX', 'HP-5')."),
        ("oven_temp_c", "float", "N/A", "Current or planned oven temperature in °C."),
        ("isothermal_hold_time_min", "float", "0.0", "Isothermal hold time in minutes."),
        ("ramp_rate_c_per_min", "None", "None", "Temperature ramp rate in °C/min for gradient mode."),
        ("final_temp_c", "None", "None", "Final temperature in °C for gradient mode."),
        ("column_length_m", "float", "30.0", "Column length in meters (affects total thermal exposure)."),
        ("film_thickness_um", "float", "0.25", "Film thickness in μm (thicker films bleed more)."),
        ("detector_type", "str", "'MS'", "Detector type: 'MS', 'FID', 'ECD', 'TCD', 'NPD', 'FPD'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters string like 'phase=DB-5 T=300 hold=10 ramp=10 Tfinal=320 detector=MS'."),
    ]

    output_sig = [
        ("bleed_assessment", "dict", "Complete bleed risk assessment including severity level, safety margin, temperature recommendations, and detector-specific notes."),
    ]

    examples = [
        {
            "code_input": {
                "stationary_phase": "DB-5",
                "oven_temp_c": 310,
                "isothermal_hold_time_min": 15,
            },
            "text_input": {"input_params": "phase=DB-5 T=310 hold=15"},
            "output": {
                "bleed_assessment": {"severity": "moderate", "max_temp": 325, "safety_margin_pct": 4.7}
            },
        },
        {
            "code_input": {
                "stationary_phase": "DB-WAX",
                "oven_temp_c": 255,
                "detector_type": "MS",
            },
            "text_input": {"input_params": "phase=DB-WAX T=255 detector=MS"},
            "output": {
                "bleed_assessment": {"severity": "critical", "max_temp": 250}
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _lookup_max_temp(self, phase: str) -> int:
        """Look up max isothermal temperature for a given stationary phase."""
        phase_upper = phase.strip()
        if phase_upper in STATIONARY_PHASE_TEMP_LIMITS:
            return STATIONARY_PHASE_TEMP_LIMITS[phase_upper]
        # Fuzzy match
        for key, val in STATIONARY_PHASE_TEMP_LIMITS.items():
            if phase_upper.lower() == key.lower():
                return val
            if phase_upper.lower() in key.lower() or key.lower() in phase_upper.lower():
                return val
        raise ChemMCPError(
            f"Unknown stationary phase '{phase}'. "
            f"Known phases include: {', '.join(sorted(set(STATIONARY_PHASE_TEMP_LIMITS.keys()))[:20])}..."
        )

    def _assess_severity(self, ratio: float) -> tuple:
        """Return (severity_level, color_emoji, description)."""
        if ratio <= BLEED_THRESHOLDS["safe"]:
            return ("safe", "✅", "Negligible bleed — well within safe operating range")
        elif ratio <= BLEED_THRESHOLDS["moderate"]:
            return ("moderate", "⚠️", "Moderate bleed — expect slight baseline drift")
        elif ratio <= BLEED_THRESHOLDS["high"]:
            return ("high", "🔶", "High bleed — significant baseline rise, reduced sensitivity")
        else:
            return ("critical", "🚨", "CRITICAL — exceeds max temp! Column damage likely")

    def _detector_sensitivity(self, detector: str) -> dict:
        """Detector-specific bleed impact assessment."""
        info = {
            "MS": {
                "impact": "HIGH",
                "note": "Mass spectrometer is extremely sensitive to column bleed; "
                        "siloxane ions (m/z 73, 147, 207, 281, 355) will dominate low-mass range. "
                        "Consider using MS-grade columns and keeping temp ≥20°C below limit.",
            },
            "FID": {
                "impact": "LOW-MODERATE",
                "note": "FID tolerates moderate bleed well; baseline may rise slightly. "
                        "Generally not a concern unless doing trace analysis.",
            },
            "ECD": {
                "impact": "VERY HIGH",
                "note": "ECD is highly sensitive to halogenated bleed products. "
                        "Keep temperature well below limit for trace-level work.",
            },
            "TCD": {
                "impact": "MODERATE",
                "note": "TCD responds to all compounds; bleed increases baseline noise.",
            },
            "NPD": {
                "impact": "HIGH",
                "note": "N-P detector sensitive to nitrogen-containing bleed fragments. "
                        "Use lower temperatures for N/P trace analysis.",
            },
            "FPD": {
                "impact": "MODERATE-HIGH",
                "note": "Sulfur/phosphorus mode can be affected by S/P containing bleed species.",
            },
        }
        return info.get(detector.upper(), {"impact": "UNKNOWN", "note": f"No specific guidance for {detector}."})

    def _calc_thermal_stress(self, temp: float, hold_time: float, length: float,
                              film: float, ramp_rate: Optional[float],
                              final_temp: Optional[float]) -> dict:
        """Estimate cumulative thermal stress."""
        base_stress = (temp / 200.0) ** 2 * (hold_time / 10.0) * (film / 0.25)
        stress_factor = round(base_stress * (length / 30.0), 2)

        ramp_stress = None
        if ramp_rate and final_temp and final_temp > temp:
            avg_ramp_temp = (temp + final_temp) / 2.0
            ramp_duration = (final_temp - temp) / max(ramp_rate, 0.1)
            ramp_stress = round((avg_ramp_temp / 200.0) ** 2 * (ramp_duration / 10.0) * (film / 0.25), 2)

        return {
            "isothermal_stress_index": stress_factor,
            "ramp_stress_index": ramp_stress,
            "interpretation": (
                "Low (<1)" if stress_factor < 1 else
                "Moderate (1-5)" if stress_factor < 5 else
                "High (5-20)" if stress_factor < 20 else
                "Very high (>20)"
            ),
        }

    def _run_base(self, stationary_phase: str, oven_temp_c: float,
                  isothermal_hold_time_min: float = 0.0,
                  ramp_rate_c_per_min: Optional[float] = None,
                  final_temp_c: Optional[float] = None,
                  column_length_m: float = 30.0,
                  film_thickness_um: float = 0.25,
                  detector_type: str = "MS") -> dict:
        """Core logic: assess GC column bleed risk."""

        max_temp = self._lookup_max_temp(stationary_phase)
        ratio = oven_temp_c / max_temp
        severity, emoji, desc = self._assess_severity(ratio)

        # Final temp check during ramp
        final_severity = None
        if final_temp_c:
            final_ratio = final_temp_c / max_temp
            final_severity = self._assess_severity(final_ratio)[0]

        # Thermal stress
        thermal = self._calc_thermal_stress(
            oven_temp_c, isothermal_hold_time_min, column_length_m,
            film_thickness_um, ramp_rate_c_per_min, final_temp_c,
        )

        # Detector-specific advice
        det_info = self._detector_sensitivity(detector_type)

        # Recommendations
        recommendations = self._generate_recommendations(
            severity, ratio, max_temp, detector_type,
            film_thickness_um, ramp_rate_c_per_min, final_temp_c,
        )

        result = {
            "bleed_assessment": {
                "column_info": {
                    "stationary_phase": stationary_phase,
                    "max_isothermal_temperature_c": max_temp,
                    "current_operating_temp_c": oven_temp_c,
                    "temperature_utilization_pct": round(ratio * 100, 1),
                },
                "risk_evaluation": {
                    "severity_level": severity,
                    "status_emoji": emoji,
                    "description": desc,
                    "safety_margin_c": round(max_temp - oven_temp_c, 1),
                    "safety_margin_pct": round((1 - ratio) * 100, 1),
                },
                "gradient_analysis": {
                    "has_ramp": ramp_rate_c_per_min is not None,
                    "ramp_rate_c_per_min": ramp_rate_c_per_min,
                    "final_temp_c": final_temp_c,
                    "final_temp_severity": final_severity,
                } if ramp_rate_c_per_min else None,
                "thermal_stress": thermal,
                "detector_impact": {
                    "detector_type": detector_type.upper(),
                    **det_info,
                },
                "recommendations": recommendations,
                "siloxane_ions_to_monitor": [73, 147, 207, 281, 355] if detector_type.upper() == "MS" else None,
            }
        }
        return result

    def _generate_recommendations(self, severity: str, ratio: float, max_temp: float,
                                   detector: str, film: float,
                                   ramp_rate: Optional[float], final_temp: Optional[float]) -> List[str]:
        recs = []

        if severity == "critical":
            recs.append(f"🚨 IMMEDIATE: Reduce oven temperature below {max_temp}°C to prevent permanent column damage!")
            recs.append("Consider replacing the column if operated above limit for extended periods.")
        elif severity == "high":
            recs.append(f"Reduce temperature by at least {round((ratio - 0.90) * max_temp)}°C for safer operation.")
            recs.append("Limit isothermal hold time to minimize cumulative thermal exposure.")
        elif severity == "moderate":
            recs.append(f"Acceptable for short-term runs; consider reducing by 10-20°C for long sequences.")

        if film > 0.5:
            recs.append("Thick film (>0.5μm) increases bleed — consider thin-film alternative for high-temp work.")
        if detector.upper() == "MS" and ratio > 0.80:
            recs.append("For MS detection: use MS-certified column, set ion source to exclude m/z < 40 if possible.")
        if detector.upper() == "ECD" and ratio > 0.75:
            recs.append("For ECD: keep temperature at least 30°C below limit due to high sensitivity to bleed.")
        if final_temp and final_temp > max_temp * 0.95:
            recs.append(f"Gradient final temp ({final_temp}°C) approaches/exceeds limit — reduce final temp or shorten high-temp hold.")

        if severity in ("safe", "moderate"):
            recs.append("✓ Current conditions are within acceptable range for routine analysis.")

        return recs[:8]

    def _run_text(self, input_params: str) -> dict:
        """Parse key=value format input."""
        kwargs = {}
        for part in input_params.strip().split():
            if "=" in part:
                key, val = part.split("=", 1)
                key_map = {
                    "phase": "stationary_phase", "T": "oven_temp_c",
                    "hold": "isothermal_hold_time_min", "ramp": "ramp_rate_c_per_min",
                    "Tfinal": "final_temp_c", "L": "column_length_m",
                    "df": "film_thickness_um", "detector": "detector_type",
                }
                mapped_key = key_map.get(key, key)
                try:
                    kwargs[mapped_key] = float(val)
                except ValueError:
                    kwargs[mapped_key] = val
        return self._run_base(**kwargs)
