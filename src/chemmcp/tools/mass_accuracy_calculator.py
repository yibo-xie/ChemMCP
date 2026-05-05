"""
Mass Accuracy Calculator - Calculates mass error in ppm, mmu, Da units
between theoretical and observed masses for MS quality control.
"""

import logging
import math
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Mass accuracy grade thresholds (ppm)
_ACCURACY_GRADES = [
    ("exceptional", 1.0,   "Orbitrap/FT-ICR at high resolution — excellent calibration"),
    ("excellent",   3.0,   "High-resolution Q-TOF with good calibration"),
    ("very_good",   5.0,   "Standard Q-TOF performance"),
    ("good",        10.0,  "Acceptable for most LC-MS applications"),
    ("moderate",    20.0,  "Unit resolution instrument or low-mass region"),
    ("poor",        50.0,  "Needs recalibration or investigation"),
    ("unacceptable", float("inf"), "Instrument requires immediate maintenance"),
]


@ChemMCPManager.register_tool
class MassAccuracyCalculator(BaseTool):
    """
    质量精度计算器 — 计算理论质量与观测质量之间的误差（ppm、mmu、Da）。
    
    用于质谱质量控制，评估仪器精度状态和校准需求。
    """
    __version__      = "0.1.0"
    name             = "MassAccuracyCalculator"
    func_name        = "calculate_mass_accuracy"
    description      = "Calculate mass accuracy in ppm, mmu, and Da between theoretical and observed mass values."
    implementation_description = "Computes mass error as: error_ppm = (observed - theoretical) / theoretical × 10⁶; error_mmu = (observed - theoretical) × 1000; error_da = observed - theoretical. Provides accuracy grading and calibration suggestions based on typical instrument performance benchmarks."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Mass Spectrometry", "Accuracy", "Calibration", "Quality Control", "ppm"]
    required_envs    = []

    code_input_sig   = [
        ("theoretical_mass", "float", "N/A", "Theoretical exact mass (in Daltons)."),
        ("observed_mass", "float", "N/A", "Observed m/z value from the instrument."),
        ("charge_state", "int", "1", "Charge state of the ion (used to convert m/z to neutral mass if needed)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'theoretical_mass observed_mass [charge_state]'. Example: '286.1438 286.1445 1'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with error_ppm, error_mmu, error_da, accuracy_grade, grade_description, calibration_suggestion, and charge-corrected values."),
    ]

    examples         = [
        {
            "code_input": {
                "theoretical_mass": 286.1438,
                "observed_mass": 286.1445,
                "charge_state": 1,
            },
            "text_input": {
                "input_params": "286.1438 286.1445 1"
            },
            "output": {
                "result": {
                    "theoretical_mass": 286.1438,
                    "observed_mass": 286.1445,
                    "error_da": 0.0007,
                    "error_mmu": 0.70,
                    "error_ppm": 2.45,
                    "accuracy_grade": "excellent",
                    "grade_description": "High-resolution Q-TOF with good calibration",
                    "calibration_suggestion": "Within acceptable range; continue routine monitoring",
                }
            },
        },
        {
            "code_input": {
                "theoretical_mass": 553.2770,
                "observed_mass": 553.2900,
                "charge_state": 2,
            },
            "text_input": {
                "input_params": "553.2770 553.2900 2"
            },
            "output": {
                "result": {
                    "error_ppm": 23.49,
                    "accuracy_grade": "moderate",
                    "calibration_suggestion": "Consider recalibration; check lock mass",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, theoretical_mass: float, observed_mass: float, charge_state: int = 1) -> dict:
        """Core logic."""
        if theoretical_mass <= 0:
            raise ChemMCPError("Theoretical mass must be positive.")
        if observed_mass <= 0:
            raise ChemMCPError("Observed mass must be positive.")
        if charge_state == 0:
            raise ChemMCPError("Charge state cannot be zero.")

        # If comparing m/z values, convert to neutral mass equivalent
        # For ppm calculation on m/z level (common practice):
        theo_mz = theoretical_mass / abs(charge_state)
        obs_mz = observed_mass / abs(charge_state)

        # Calculate errors at m/z level (standard practice)
        error_da = obs_mz - theo_mz
        error_mmu = error_da * 1000  # millimass units
        error_ppm = (error_da / theo_mz) * 1e6 if theo_mz != 0 else 0.0

        # Also calculate at neutral mass level
        error_da_neutral = observed_mass - theoretical_mass
        error_mmu_neutral = error_da_neutral * 1000
        error_ppm_neutral = (error_da_neutral / theoretical_mass) * 1e6 if theoretical_mass != 0 else 0.0

        # Determine accuracy grade
        grade, threshold, description = self._get_grade(abs(error_ppm))

        # Calibration suggestion
        cal_suggestion = self._calibration_suggestion(abs(error_ppm), grade)

        return {
            "result": {
                "input_theoretical_mass": round(theoretical_mass, 6),
                "input_observed_mass": round(observed_mass, 6),
                "charge_state": charge_state,
                # m/z-level errors (standard reporting)
                "error_ppm_mz_level": round(error_ppm, 4),
                "error_mmu_mz_level": round(error_mmu, 4),
                "error_da_mz_level": round(error_da, 8),
                # Neutral mass-level errors
                "error_ppm_neutral_level": round(error_ppm_neutral, 4),
                "error_mmu_neutral_level": round(error_mmu_neutral, 4),
                "error_da_neutral_level": round(error_da_neutral, 8),
                # Grade and recommendation
                "accuracy_grade": grade,
                "grade_threshold_ppm": threshold,
                "grade_description": description,
                "calibration_suggestion": cal_suggestion,
                # Sign information
                "sign": "positive (heavier than expected)" if error_da > 0 else "negative (lighter than expected)",
                # Reference info
                "notes": (
                    "Mass Accuracy Notes:\n"
                    "• < 3 ppm: High-resolution instrument (Orbitrap/FT-ICR/Q-TOF)\n"
                    "• 3–10 ppm: Standard Q-TOF acceptable range\n"
                    "• 10–20 ppm: Unit-res quadrupole or low-m/z region\n"
                    f"• Current measurement: {abs(error_ppm):.2f} ppm ({grade})\n"
                    "• Always use lock mass or internal calibration for best accuracy\n"
                    "• Error sign indicates systematic bias direction"
                ),
            }
        }

    def _get_grade(self, abs_ppm: float) -> Tuple[str, float, str]:
        """Get accuracy grade from ppm value."""
        for grade, threshold, desc in _ACCURACY_GRADES:
            if abs_ppm <= threshold:
                return grade, threshold, desc
        return _ACCURACY_GRADES[-1]

    def _calibration_suggestion(self, abs_ppm: float, grade: str) -> str:
        """Generate calibration/maintenance suggestion."""
        if grade in ("exceptional", "excellent"):
            return "Excellent accuracy. Continue routine QC monitoring."
        elif grade == "very_good":
            return "Good accuracy. No action needed; include in trend monitoring."
        elif grade == "good":
            return "Acceptable for most applications. Monitor for drift trends."
        elif grade == "moderate":
            return "Consider checking calibration status. Verify lock mass performance."
        elif grade == "poor":
            return "Recommend recalibration. Check mass axis calibration solution and ion source condition."
        else:
            return "URGENT: Instrument needs immediate recalibration and possibly service. Do not run critical samples."

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            theo = float(parts[0])
            obs = float(parts[1])
            chg = int(parts[2]) if len(parts) > 2 else 1
            return self._run_base(theo, obs, chg)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'theoretical_mass observed_mass [charge_state]'")
