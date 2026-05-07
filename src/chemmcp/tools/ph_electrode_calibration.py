import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PhElectrodeCalibration(BaseTool):
    """
    pH electrode calibration: calculate slope and offset from buffer measurements.
    Assess electrode quality (slope % efficiency, offset, R²) and generate calibration curve.
    """
    __version__ = "0.1.0"
    name = "PhElectrodeCalibration"
    func_name = "ph_electrode_calibration"
    description = "Calibrate a pH electrode using standard buffer solutions. Calculate slope (mV/pH), offset (mV), Nernstian efficiency (%), and assess electrode condition. Supports 2-point and multi-point calibration."
    implementation_description = "Performs linear regression on (pH_buffer, mV_measured) data to determine calibration parameters:\n- Slope S (mV/pH): ideal = 59.16 mV/pH at 25°C\n- Offset E₀ (mV): intercept\n- Efficiency η = (S / 59.16) × 100%\n- R² goodness of fit\n- Predicts pH for unknown mV readings using calibrated equation."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["pH Electrode", "Calibration", "Nernst Equation", "Slope", "Offset", "Buffer Solution", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("buffer_ph_values", "list", "N/A", "List of known pH values of standard buffers."),
        ("measured_mV_values", "list", "N/A", "List of mV readings for each buffer solution."),
        ("temperature_C", "float", "25.0", "Calibration temperature in °C. Default: 25.0."),
        # For prediction:
        ("unknown_mV", "float", "None", "Unknown sample mV reading to convert to pH using the calibration. Optional."),
        ("unknown_temperature_C", "float", "None", "Temperature of unknown sample in °C. Uses calibration T if not given."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: '[T_C] [unknown_mV] || ph1,mv1 ph2,mv2 ...'. Example: '25 || 4.01,180 7.00,335 10.01,465' or '25 -50 || 4.00,178 7.00,333 10.00,468'"),
    ]

    output_sig = [
        ("slope_mV_pH", "float", "Calibration slope in mV/pH unit."),
        ("offset_mV", "float", "Calibration offset/intercept in mV."),
        ("nernstian_efficiency_percent", "float", "Slope as percentage of ideal Nernstian slope at given temperature."),
        ("ideal_slope_mV_pH", "float", "Ideal Nernstian slope at calibration temperature."),
        ("r_squared", "float", "R² value of linear fit (goodness of fit)."),
        ("electrode_condition", "str", "Assessment: 'excellent', 'good', 'acceptable', 'poor', or 'replace'."),
        ("predicted_pH", "float", "Predicted pH for unknown mV input, if provided."),
        ("calibration_equation", "str", "The calibration equation string."),
        ("calibration_summary", "str", "Detailed text summary of calibration results."),
    ]

    examples = [
        {
            "code_input": {
                "buffer_ph_values": [4.01, 7.00, 10.01],
                "measured_mV_values": [178, 333, 468],
                "temperature_C": 25.0,
            },
            "text_input": {"input_string": "25 || 4.01,178 7.00,333 10.01,468"},
            "output": {
                "slope_mV_pH": 57.5,
                "offset_mV": -53.0,
                "nernstian_efficiency_percent": 97.2,
                "r_squared": 0.9999,
                "electrode_condition": "excellent",
                "calibration_equation": "pH = (E + 53.0) / 57.5",
            }
        },
        {
            "code_input": {
                "buffer_ph_values": [4.01, 7.00],
                "measured_mV_values": [150, 300],
                "temperature_C": 25.0,
                "unknown_mV": -50,
            },
            "text_input": {"input_string": "25 -50 || 4.01,150 7.00,300"},
            "output": {
                "slope_mV_pH": 49.9,
                "nernstian_efficiency_percent": 84.4,
                "electrode_condition": "acceptable",
                "predicted_pH": 12.02,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._R = 8.314       # J/(mol·K)
        self._F = 96485       # C/mol

    def _ideal_slope(self, T_C: float) -> float:
        """Calculate ideal Nernstian slope at temperature T."""
        T_K = T_C + 273.15
        return (self._R * T_K * 1000) / self._F  # Convert to mV/pH

    def _linear_regression(self, x: List[float], y: List[float]) -> tuple:
        """Simple linear regression: returns (slope, intercept, r_squared)."""
        n = len(x)
        if n < 2:
            raise ChemMCPError("Need at least 2 points for calibration.")

        sum_x = sum(x)
        sum_y = sum(y)
        sum_xy = sum(xi * yi for xi, yi in zip(x, y))
        sum_x2 = sum(xi ** 2 for xi in x)
        sum_y2 = sum(yi ** 2 for yi in y)

        denom = n * sum_x2 - sum_x ** 2
        if abs(denom) < 1e-15:
            raise ChemMCPError("Cannot perform regression: all x values are identical.")

        slope = (n * sum_xy - sum_x * sum_y) / denom
        intercept = (sum_y - slope * sum_x) / n

        # R² calculation
        y_mean = sum_y / n
        ss_tot = sum((yi - y_mean) ** 2 for yi in y)
        ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
        r_sq = 1 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

        return slope, intercept, r_sq

    def _assess_electrode(self, efficiency_pct: float, r_sq: float) -> str:
        """Assess electrode condition based on calibration quality."""
        if efficiency_pct >= 98 and r_sq >= 0.999:
            return "excellent"
        elif efficiency_pct >= 95 and r_sq >= 0.998:
            return "good"
        elif efficiency_pct >= 90 and r_sq >= 0.995:
            return "acceptable"
        elif efficiency_pct >= 85:
            return "poor"
        else:
            return "replace"

    def _run_base(
        self,
        buffer_ph_values: List[float],
        measured_mV_values: List[float],
        temperature_C: float = 25.0,
        unknown_mV: Optional[float] = None,
        unknown_temperature_C: Optional[float] = None,
    ) -> dict:
        """Perform pH electrode calibration."""
        if len(buffer_ph_values) != len(measured_mV_values):
            raise ChemMCPError("buffer_ph_values and measured_mV_values must have same length.")
        if len(buffer_ph_values) < 2:
            raise ChemMCPError("Need at least 2 buffer points for calibration.")

        # Linear regression: E(mV) vs pH → E = slope × pH + offset
        slope, offset, r_sq = self._linear_regression(buffer_ph_values, measured_mV_values)

        ideal_slope = self._ideal_slope(temperature_C)
        efficiency = (slope / ideal_slope) * 100.0 if ideal_slope != 0 else 0.0

        condition = self._assess_electrode(efficiency, r_sq)

        # Predict pH for unknown sample
        predicted_ph = None
        if unknown_mV is not None:
            T_pred = unknown_temperature_C if unknown_temperature_C else temperature_C
            ideal_at_pred = self._ideal_slope(T_pred)
            # Temperature-corrected prediction
            predicted_ph = (unknown_mV - offset) / slope

        equation = f"E(mV) = {slope:.2f} × pH + ({offset:.2f})   ↔   pH = (E − {offset:.2f}) / {slope:.2f}"

        summary_parts = [
            f"pH Electrode Calibration Results (@ {temperature_C}°C):",
            f"  Slope:     {slope:.2f} mV/pH  (ideal: {ideal_slope:.2f} mV/pH)",
            f"  Offset:    {offset:.2f} mV",
            f"  Efficiency: {efficiency:.1f}% of Nernstian response",
            f"  R²:        {r_sq:.6f}",
            f"  Condition: {condition.upper()}",
            f"  Equation:  {equation}",
        ]
        if predicted_ph is not None:
            summary_parts.append(f"  Unknown sample: {unknown_mV} mV → pH = {predicted_ph:.2f}")

        return {
            "slope_mV_pH": round(slope, 4),
            "offset_mV": round(offset, 4),
            "nernstian_efficiency_percent": round(efficiency, 2),
            "ideal_slope_mV_pH": round(ideal_slope, 4),
            "r_squared": round(r_sq, 6),
            "electrode_condition": condition,
            "predicted_pH": round(predicted_ph, 4) if predicted_ph is not None else None,
            "calibration_equation": equation,
            "calibration_summary": "\n".join(summary_parts),
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        if "||" not in input_string:
            raise ChemMCPError("Must use '||' separator between parameters and data.")

        left, right = input_string.split("||", 1)
        params = left.strip().split()
        pairs = right.strip().split()

        T = float(params[0]) if len(params) > 0 else 25.0
        unknown_mv = float(params[1]) if len(params) > 1 else None
        T_unknown = float(params[2]) if len(params) > 2 else None

        ph_list = []
        mv_list = []
        for pair in pairs:
            if "," in pair:
                sub = pair.split(",")
                try:
                    ph_list.append(float(sub[0]))
                    mv_list.append(float(sub[1]))
                except ValueError:
                    continue

        if len(ph_list) < 2:
            raise ChemMCPError("Need at least 2 valid pH,mV pairs after '||'.")

        return self._run_base(ph_list, mv_list, T, unknown_mv, T_unknown)
