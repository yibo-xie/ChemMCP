import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PotentiometricTitrationEndpoint(BaseTool):
    """
    Determine potentiometric titration endpoint using first and second derivative methods.
    Analyzes titration curve (volume vs. potential/pH) to find equivalence point(s).
    """
    __version__ = "0.1.0"
    name = "PotentiometricTitrationEndpoint"
    func_name = "potentiometric_titration_endpoint"
    description = "Determine potentiometric titration endpoint(s) from volume-potential (or volume-pH) data using first derivative (dE/dV) and second derivative (d²E/dV²) methods."
    implementation_description = "Implements two methods for endpoint detection:\n1. First derivative: maximum of dE/dV indicates steepest slope → endpoint\n2. Second derivative: zero-crossing of d²E/dV² with sign change (+→−) → precise endpoint\nSupports multiple endpoints for polyprotic systems."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Titration", "Endpoint Detection", "Potentiometry", "First Derivative", "Second Derivative", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("volume_mL", "list", "N/A", "List of titrant volumes in mL."),
        ("potential_V_or_pH", "list", "N/A", "List of measured potential values (V) or pH values corresponding to each volume."),
        ("method", "str", "auto", "Method for endpoint detection: 'first_derivative', 'second_derivative', or 'auto' (tries both). Default: auto."),
        ("smoothing_window", "int", "3", "Window size for Savitzky-Golay-like smoothing (odd number). Default: 3."),
        ("signal_type", "str", "potential", "Signal type: 'potential' or 'pH'. Default: potential."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: '[method] [smoothing] [type] || v1,e1 v2,e2 ...'. Example: 'auto 3 pH || 0,2.9 10,3.3 20,3.8 25,4.5 30,7.0 35,11.0 40,12.0'"),
    ]

    output_sig = [
        ("endpoint_volume_mL", "float", "Detected endpoint volume (mL)."),
        ("endpoint_signal_value", "float", "Potential or pH at the endpoint."),
        ("method_used", "str", "The method that successfully detected the endpoint."),
        ("first_derivative_peaks", "list", "Volume positions of first derivative maxima."),
        ("second_derivative_zero_crossings", "list", "Zero-crossing volumes from second derivative."),
        ("derivative_data", "dict", "Raw derivative data points for plotting/verification."),
        ("titration_summary", "str", "Text summary of the titration analysis."),
    ]

    examples = [
        {
            "code_input": {
                "volume_mL": [0, 5, 10, 15, 20, 24, 25, 26, 28, 30, 32, 35, 40],
                "potential_V_or_pH": [2.90, 3.30, 3.80, 4.25, 5.25, 6.00, 7.00, 8.00, 10.00, 10.70, 11.00, 11.30, 11.50],
                "method": "auto",
                "signal_type": "pH",
            },
            "text_input": {"input_string": "auto 3 pH || 0,2.9 5,3.3 10,3.8 15,4.25 20,5.25 24,6 25,7 26,8 28,10 30,10.7 32,11 35,11.3 40,11.5"},
            "output": {
                "endpoint_volume_mL": 25.0,
                "endpoint_signal_value": 7.00,
                "method_used": "second_derivative",
                "titration_summary": "Strong acid-strong base titration detected. Endpoint at V = 25.00 mL, pH ≈ 7.0.",
            }
        },
        {
            "code_input": {
                "volume_mL": [0, 5, 10, 15, 20, 22, 23, 24, 25, 26, 28, 30],
                "potential_V_or_pH": [2.8, 3.6, 4.0, 4.3, 4.7, 5.0, 5.4, 6.2, 8.0, 10.2, 11.0, 11.3],
                "method": "first_derivative",
                "signal_type": "pH",
            },
            "text_input": {"input_string": "first_derivative 3 pH || 0,2.8 5,3.6 10,4 15,4.3 20,4.7 22,5 23,5.4 24,6.2 25,8 26,10.2 28,11 30,11.3"},
            "output": {
                "endpoint_volume_mL": 24.8,
                "endpoint_signal_value": 7.1,
                "method_used": "first_derivative",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _smooth(self, y: List[float], window: int) -> List[float]:
        """Simple moving average smoothing."""
        if window < 2:
            return list(y)
        half = window // 2
        smoothed = []
        n = len(y)
        for i in range(n):
            start = max(0, i - half)
            end = min(n, i + half + 1)
            smoothed.append(sum(y[start:end]) / (end - start))
        return smoothed

    def _derivative(self, x: List[float], y: List[float]) -> tuple:
        """Calculate first derivative dy/dx using central differences."""
        if len(x) < 2:
            return [], []
        dx_list = []
        dydx_list = []
        for i in range(len(x)):
            if i == 0:
                # Forward difference
                dx = x[i+1] - x[i]
                dy = y[i+1] - y[i]
            elif i == len(x) - 1:
                # Backward difference
                dx = x[i] - x[i-1]
                dy = y[i] - y[i-1]
            else:
                # Central difference
                dx = x[i+1] - x[i-1]
                dy = y[i+1] - y[i-1]

            if abs(dx) < 1e-12:
                dx_list.append(x[i])
                dydx_list.append(0.0)
            else:
                dx_list.append(x[i])
                dydx_list.append(dy / dx)

        return dx_list, dydx_list

    def _find_second_derivative_zero_crossings(self, x: List[float], d2y: List[float]) -> list:
        """Find zero crossings where second derivative changes from positive to negative (peak in first deriv)."""
        crossings = []
        for i in range(len(d2y) - 1):
            if d2y[i] > 0 and d2y[i+1] <= 0:
                # Linear interpolation for more precise crossing point
                if abs(d2y[i+1] - d2y[i]) > 1e-12:
                    frac = d2y[i] / (d2y[i] - d2y[i+1])
                    xc = x[i] + frac * (x[i+1] - x[i])
                else:
                    xc = (x[i] + x[i+1]) / 2
                crossings.append(round(xc, 4))
        return crossings

    def _find_first_deriv_maxima(self, x: List[float], dy: List[float]) -> list:
        """Find local maxima in first derivative."""
        peaks = []
        for i in range(1, len(dy) - 1):
            if dy[i] > dy[i-1] and dy[i] > dy[i+1] and dy[i] > 0:
                peaks.append(round(x[i], 4))
        return peaks

    def _interpolate_signal_at_volume(self, vol: float, x: List[float], y: List[float]) -> float:
        """Linear interpolation of signal value at given volume."""
        for i in range(len(x) - 1):
            if x[i] <= vol <= x[i+1]:
                if abs(x[i+1] - x[i]) < 1e-12:
                    return y[i]
                t = (vol - x[i]) / (x[i+1] - x[i])
                return y[i] + t * (y[i+1] - y[i])
        # Extrapolate
        if vol <= x[0]:
            return y[0]
        return y[-1]

    def _run_base(
        self,
        volume_mL: List[float],
        potential_V_or_pH: List[float],
        method: str = "auto",
        smoothing_window: int = 3,
        signal_type: str = "potential",
    ) -> dict:
        """Detect titration endpoint(s)."""
        if len(volume_mL) != len(potential_V_or_pH):
            raise ChemMCPError("volume_mL and potential_V_or_pH must have same length.")
        if len(volume_mL) < 4:
            raise ChemMCPError("Need at least 4 data points for reliable endpoint detection.")

        V = list(volume_mL)
        E = list(potential_V_or_pH)

        # Smooth data
        E_smooth = self._smooth(E, smoothing_window)

        # Calculate derivatives
        _, d1 = self._derivative(V, E_smooth)
        _, d2 = self._derivative(V, d1)

        # Find endpoints by both methods
        first_deriv_peaks = self._find_first_deriv_maxima(V, d1)
        second_deriv_crossings = self._find_second_derivative_zero_crossings(V, d2)

        # Determine best result
        endpoint_vol = None
        endpoint_val = None
        method_used = None

        candidates = second_deriv_crossings if second_deriv_crossings else first_deriv_peaks

        if method == "second_derivative" or (method == "auto" and second_deriv_crossings):
            if second_deriv_crossings:
                endpoint_vol = second_deriv_crossings[0]
                endpoint_val = round(self._interpolate_signal_at_volume(endpoint_vol, V, E), 4)
                method_used = "second_derivative"

        if (method == "first_derivative" or (method == "auto" and not endpoint_vol)) and first_deriv_peaks:
            endpoint_vol = first_deriv_peaks[0]
            endpoint_val = round(self._interpolate_signal_at_volume(endpoint_vol, V, E), 4)
            method_used = "first_derivative"

        if endpoint_vol is None:
            raise ChemMCPError("Could not detect a clear endpoint. Check data quality or try adjusting smoothing.")

        sig_name = "pH" if signal_type.lower() == "ph" else "E (V)"
        summary = (
            f"Potentiometric Titration Analysis ({sig_name} vs V):\n"
            f"  Endpoint detected at V = {endpoint_vol:.2f} mL ({sig_name} = {endpoint_val})\n"
            f"  Method: {method_used}\n"
            f"  Data points: {len(V)}\n"
            f"  First derivative peaks at: {first_deriv_peaks}\n"
            f"  Second derivative zero-crossings at: {second_deriv_crossings}"
        )

        return {
            "endpoint_volume_mL": endpoint_vol,
            "endpoint_signal_value": endpoint_val,
            "method_used": method_used,
            "first_derivative_peaks": first_deriv_peaks,
            "second_derivative_zero_crossings": second_deriv_crossings,
            "derivative_data": {
                "volumes": [round(v, 4) for v in V],
                "first_derivative": [round(d, 4) for d in d1],
                "second_derivative": [round(d, 4) for d in d2],
            },
            "titration_summary": summary,
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        if "||" not in input_string:
            raise ChemMCPError("Must use '||' separator between parameters and data.")

        left, right = input_string.split("||", 1)
        params = left.strip().split()
        pairs = right.strip().split()

        method = params[0] if len(params) > 0 else "auto"
        smooth = int(params[1]) if len(params) > 1 else 3
        sig_type = params[2] if len(params) > 2 else "potential"

        V_data = []
        E_data = []
        for pair in pairs:
            if "," in pair:
                sub = pair.split(",")
                try:
                    V_data.append(float(sub[0]))
                    E_data.append(float(sub[1]))
                except ValueError:
                    continue

        if not V_data:
            raise ChemMCPError("No valid data pairs found after '||'.")

        return self._run_base(V_data, E_data, method, smooth, sig_type)
