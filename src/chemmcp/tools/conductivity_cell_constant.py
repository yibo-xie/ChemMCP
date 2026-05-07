import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ConductivityCellConstant(BaseTool):
    """
    Conductivity cell constant (Kcell) calibration and correction.
    Calculate cell constant from standard solution, or correct measured conductivity using known cell constant.
    """
    __version__ = "0.1.0"
    name = "ConductivityCellConstant"
    func_name = "conductivity_cell_constant"
    description = "Calculate conductivity cell constant from KCl standard solution calibration, correct measured conductivity values, and estimate measurement uncertainty."
    implementation_description = "Cell constant Kcell = κ_standard / G_measured, where κ is the known conductivity of a KCl standard solution. Corrected conductivity: κ_sample = Kcell × G_sample. Supports multiple KCl standard concentrations (0.01 M, 0.1 M, 1.0 M) with temperature correction to 25°C."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Conductivity", "Cell Constant", "Calibration", "KCl Standard", "Temperature Correction", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Operation mode: 'calibrate' (calculate Kcell) or 'correct' (apply Kcell to sample)."),
        # For calibrate:
        ("standard_conductivity_uS_cm", "float", "N/A", "Known conductivity of KCl standard at 25°C (µS/cm)."),
        ("measured_conductance_uS", "float", "N/A", "Measured conductance/conductivity reading (µS or µS/cm, depending on meter)."),
        ("temperature_C", "float", "25.0", "Measurement temperature in °C. Default: 25.0."),
        # For correct:
        ("cell_constant_per_cm", "float", "N/A", "Known cell constant in cm⁻¹ (for 'correct' mode)."),
        ("sample_reading_uS", "float", "N/A", "Raw conductivity reading of sample (µS/cm as displayed by meter)."),
        ("sample_temperature_C", "float", "25.0", "Sample temperature in °C. Default: 25.0."),
        # Optional for both:
        ("kcl_concentration_M", "float", "0.1", "KCl standard concentration for lookup: 0.01, 0.1, or 1.0 M. Default: 0.1."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Formats:\n- Calibrate: 'calibrate kappa_std_uS reading [T_C]'\n- Correct: 'correct Kcell reading [T_sample]'\nExample: 'calibrate 1413 1280 24.5' or 'correct 0.98 520 22.0'"),
    ]

    output_sig = [
        ("cell_constant_per_cm", "float", "Calculated or used cell constant (cm⁻¹)."),
        ("corrected_conductivity_uS_cm", "float", "Corrected conductivity at 25°C (µS/cm). Only for 'correct' mode."),
        ("mode_used", "str", "Operation mode executed."),
        ("temperature_correction_factor", "float", "Temperature correction factor applied."),
        ("kcl_standard_info", "dict", "KCl standard solution reference data used."),
        ("formula_applied", "str", "Formula with substituted values."),
        ("notes", "str", "Additional notes on accuracy and recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "calibrate",
                "standard_conductivity_uS_cm": 1413,
                "measured_conductance_uS": 1380,
                "temperature_C": 24.5,
            },
            "text_input": {"input_string": "calibrate 1413 1380 24.5"},
            "output": {
                "cell_constant_per_cm": 1.026,
                "mode_used": "calibrate",
                "formula_applied": "Kcell = 1413 / 1380 × f_T(24.5→25) = 1.024 cm⁻¹",
            }
        },
        {
            "code_input": {
                "mode": "correct",
                "cell_constant_per_cm": 1.02,
                "sample_reading_uS": 520,
                "sample_temperature_C": 22.0,
            },
            "text_input": {"input_string": "correct 1.02 520 22.0"},
            "output": {
                "cell_constant_per_cm": 1.02,
                "corrected_conductivity_uS_cm": 554.7,
                "mode_used": "correct",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        # KCl standard conductivity values at 25°C (µS/cm)
        self._kcl_standards = {
            0.01: {"conductivity_uS_cm": 1408.8, "name": "0.01 M KCl"},
            0.10: {"conductivity_uS_cm": 12880.0, "name": "0.1 M KCl"},
            1.00: {"conductivity_uS_cm": 111340.0, "name": "1.0 M KCl"},
        }

    def _temp_correction_factor(self, T_C: float, T_ref: float = 25.0) -> float:
        """
        Temperature correction factor using ISO 7888 linear approximation.
        κ(T_ref) ≈ κ(T) × [1 + α × (T_ref - T)]
        α ≈ 0.02/°C for most aqueous solutions.
        """
        alpha = 0.02  # Typical temperature coefficient (%/°C)
        return 1.0 + alpha * (T_ref - T_C)

    def _run_base(
        self,
        mode: str,
        standard_conductivity_uS_cm: Optional[float] = None,
        measured_conductance_uS: Optional[float] = None,
        temperature_C: float = 25.0,
        cell_constant_per_cm: Optional[float] = None,
        sample_reading_uS: Optional[float] = None,
        sample_temperature_C: float = 25.0,
        kcl_concentration_M: float = 0.1,
    ) -> dict:
        """Perform cell constant calibration or conductivity correction."""
        mode = mode.lower().strip()

        if mode == "calibrate":
            if standard_conductivity_uS_cm is None or measured_conductance_uS is None:
                raise ChemMCPError("Calibrate mode requires standard_conductivity_uS_cm and measured_conductance_uS.")

            if measured_conductance_uS == 0:
                raise ChemMCPError("Measured conductance cannot be zero.")

            # Temperature-correct the standard value to measurement temperature
            f_T = self._temp_correction_factor(temperature_C)
            kappa_at_T = standard_conductivity_uS_cm / f_T  # What we'd expect at T

            # Cell constant
            Kcell = kappa_at_T / measured_conductance_uS

            std_info = self._kcl_standards.get(kcl_concentration_M, {})
            formula = (
                f"Kcell = κ_std / G_meas × temp_correction\n"
                f"       = {standard_conductivity_uS_cm} / {measured_conductance_uS} "
                f"× (1/{f_T:.4f})\n"
                f"       = {Kcell:.4f} cm⁻¹"
            )

            notes = ""
            if 0.8 <= Kcell <= 1.5:
                notes = f"Kcell = {Kcell:.3f} cm⁻¹ is within typical range (0.8–1.5 cm⁻¹)."
            else:
                notes = f"Kcell = {Kcell:.3f} cm⁻¹ is outside typical range; verify electrode condition."

            return {
                "cell_constant_per_cm": round(Kcell, 6),
                "corrected_conductivity_uS_cm": None,
                "mode_used": "calibrate",
                "temperature_correction_factor": round(f_T, 6),
                "kcl_standard_info": std_info,
                "formula_applied": formula,
                "notes": notes,
            }

        elif mode == "correct":
            if cell_constant_per_cm is None or sample_reading_uS is None:
                raise ChemMCPError("Correct mode requires cell_constant_per_cm and sample_reading_uS.")

            # Apply cell constant
            kappa_raw = cell_constant_per_cm * sample_reading_uS

            # Temperature correction to 25°C
            f_T = self._temp_correction_factor(sample_temperature_C)
            kappa_corrected = kappa_raw * f_T

            formula = (
                f"κ_corrected = Kcell × G_sample × f_T\n"
                f"           = {cell_constant_per_cm} × {sample_reading_uS} × {f_T:.4f}\n"
                f"           = {kappa_corrected:.2f} µS/cm (at 25°C)"
            )

            return {
                "cell_constant_per_cm": cell_constant_per_cm,
                "corrected_conductivity_uS_cm": round(kappa_corrected, 2),
                "mode_used": "correct",
                "temperature_correction_factor": round(f_T, 6),
                "kcl_standard_info": {},
                "formula_applied": formula,
                "notes": f"Sample measured at {sample_temperature_C}°C, corrected to 25°C.",
            }
        else:
            raise ChemMCPError(f"Unknown mode: '{mode}'. Use 'calibrate' or 'correct'.")

    def _run_text(self, input_string: str) -> dict:
        """Parse text input."""
        parts = input_string.strip().split()
        if not parts:
            raise ChemMCPError("Input cannot be empty.")

        mode = parts[0].lower()

        if mode == "calibrate":
            kappa = float(parts[1]) if len(parts) > 1 else _m("standard conductivity")
            reading = float(parts[2]) if len(parts) > 2 else _m("reading")
            T = float(parts[3]) if len(parts) > 3 else 25.0
            return self._run_base(mode=mode, standard_conductivity_uS_cm=kappa, measured_conductance_uS=reading, temperature_C=T)

        elif mode == "correct":
            kcell = float(parts[1]) if len(parts) > 1 else _m("cell constant")
            reading = float(parts[2]) if len(parts) > 2 else _m("sample reading")
            T = float(parts[3]) if len(parts) > 3 else 25.0
            return self._run_base(mode=mode, cell_constant_per_cm=kcell, sample_reading_uS=reading, sample_temperature_C=T)

        else:
            raise ChemMCPError(f"Unknown mode: '{mode}'. Use 'calibrate' or 'correct'.")


def _m(name: str):
    raise ValueError(f"Missing required parameter: {name}")
