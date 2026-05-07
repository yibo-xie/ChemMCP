import logging
from typing import Dict, List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _central_diff(x_data: List[float], y_data: List[float]) -> List[float]:
    """Central difference for interior points, forward/backward at boundaries."""
    n = len(x_data)
    if n < 2:
        raise ChemMCPError("Need at least 2 data points.")
    derivs = [0.0] * n

    # Forward difference for first point
    derivs[0] = (y_data[1] - y_data[0]) / (x_data[1] - x_data[0])

    # Central difference for interior
    for i in range(1, n - 1):
        dx = x_data[i + 1] - x_data[i - 1]
        if abs(dx) < 1e-15:
            derivs[i] = 0.0
        else:
            derivs[i] = (y_data[i + 1] - y_data[i - 1]) / dx

    # Backward difference for last point
    derivs[n - 1] = (y_data[n - 1] - y_data[n - 2]) / (x_data[n - 1] - x_data[n - 2])

    return [round(d, 8) for d in derivs]


def _forward_diff(x_data: List[float], y_data: List[float]) -> List[float]:
    """Forward difference: f'(xi) ≈ (f(xi+1) - f(xi)) / (xi+1 - xi)."""
    n = len(x_data)
    derivs = []
    for i in range(n - 1):
        dx = x_data[i + 1] - x_data[i]
        if abs(dx) < 1e-15:
            derivs.append(0.0)
        else:
            derivs.append(round((y_data[i + 1] - y_data[i]) / dx, 8))
    # Last point: use same as previous
    derivs.append(derivs[-1] if derivs else 0.0)
    return derivs


def _backward_diff(x_data: List[float], y_data: List[float]) -> List[float]:
    """Backward difference: f'(xi) ≈ (f(xi) - f(xi-1)) / (xi - xi-1)."""
    n = len(x_data)
    derivs = [0.0]  # First point: use next
    for i in range(1, n):
        dx = x_data[i] - x_data[i - 1]
        if abs(dx) < 1e-15:
            derivs.append(0.0)
        else:
            derivs.append(round((y_data[i] - y_data[i - 1]) / dx, 8))
    derivs[0] = derivs[1] if n > 1 else 0.0
    return derivs


def _savitzky_golay(x_data: List[float], y_data: List[float],
                    window_size: int = 5, poly_order: int = 2) -> tuple:
    """
    Savitzky-Golay filter for smoothing and differentiation.
    Returns (first_derivatives, second_derivatives).
    """
    try:
        from scipy.signal import savgol_filter
        import numpy as np
        first = savgol_filter(y_data, window_size, poly_order, deriv=1, delta=x_data[1]-x_data[0] if len(x_data)>1 else 1.0)
        second = savgol_filter(y_data, window_size, poly_order, deriv=2, delta=x_data[1]-x_data[0] if len(x_data)>1 else 1.0)
        return ([round(float(d), 8) for d in first],
                [round(float(d), 8) for d in second])
    except ImportError:
        pass

    # Pure Python fallback: polynomial fitting in sliding window
    return _savitzky_golay_pure(x_data, y_data, window_size, poly_order)


def _savitzky_golay_pure(x_data: List[float], y_data: List[float],
                         window_size: int, poly_order: int) -> tuple:
    """Pure Python Savitzky-Golay implementation."""
    n = len(y_data)
    half_w = window_size // 2
    first_deriv = [0.0] * n
    second_deriv = [0.0] * n

    try:
        import numpy as np
        has_numpy = True
    except ImportError:
        has_numpy = False

    for i in range(n):
        left = max(0, i - half_w)
        right = min(n - 1, i + half_w)
        x_win = [x_data[j] - x_data[i] for j in range(left, right + 1)]
        y_win = [y_data[j] for j in range(left, right + 1)]

        if has_numpy:
            coeffs = np.polyfit(x_win, y_win, min(poly_order, len(x_win) - 1))
            # First derivative of polynomial at x=0 is the linear coeff
            pder1 = np.polyder(coeffs)
            pder2 = np.polyder(pder1)
            first_deriv[i] = round(float(np.polyval(pder1, 0)), 8)
            second_deriv[i] = round(float(np.polyval(pder2, 0)), 8)
        else:
            # Simple central diff within the window
            wlen = right - left
            if wlen >= 2 and i > left and i < right:
                dx = x_data[i + 1] - x_data[i - 1]
                if abs(dx) > 1e-15:
                    first_deriv[i] = round((y_data[i + 1] - y_data[i - 1]) / dx, 8)

    return first_deriv, second_deriv


@ChemMCPManager.register_tool
class NumericalDifferentiation(BaseTool):
    """
    数值微分工具 —— 实验数据斜率提取。
    对离散实验数据进行数值微分，支持多种差分方法。
    """
    __version__ = "0.1.0"
    name = "NumericalDifferentiation"
    func_name = "numerical_differentiation"
    description = (
        "Compute numerical derivatives of discrete experimental data. "
        "Supports multiple methods: central/forward/backward differences "
        "and Savitzky-Golay filtering. Essential for kinetic data analysis."
    )
    implementation_description = (
        "Implements finite difference schemes (central, forward, backward) "
        "and Savitzky-Golay smooth differentiation with optional scipy backend."
    )
    oss_dependencies = [
        ("scipy", "https://scipy.org/", "BSD-3-Clause"),
        ("numpy", "https://numpy.org/", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Differentiation", "Experimental Data", "Kinetics", "Derivative"]
    required_envs = []

    code_input_sig = [
        ("x_data", "list", "N/A", "List of x values (independent variable)."),
        ("y_data", "list", "N/A", "List of y values (dependent variable)."),
        ("method", "str", "central", "Method: 'central', 'forward', 'backward', or 'savitzky_golay'."),
        ("window_size", "int", "5", "Window size for savitzky_golay method (default: 5, must be odd)."),
        ("poly_order", "int", "2", "Polynomial order for savitzky_golay (default: 2)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Format: 'x1,x2,x3; y1,y2,y3; [method]; [window_size]; [poly_order]' "
         "Example: '0,1,2,3,4; 0,1,4,9,16; central'"),
    ]

    output_sig = [
        ("derivatives", "list", "First derivative dy/dx at each data point."),
        ("x_at_derivative", "list", "x values corresponding to derivatives."),
        ("second_derivatives", "list or None", "Second derivative (only for savitzky_golay method)."),
    ]

    examples = [
        {
            "code_input": {
                "x_data": [0.0, 1.0, 2.0, 3.0, 4.0],
                "y_data": [0.0, 1.0, 4.0, 9.0, 16.0],
                "method": "central",
            },
            "text_input": {
                "input_str": "0,1,2,3,4; 0,1,4,9,16; central",
            },
            "output": {
                "derivatives": [1.0, 2.0, 4.0, 6.0, 7.0],
                "x_at_derivative": [0.0, 1.0, 2.0, 3.0, 4.0],
                "second_derivatives": None,
            },
        },
        {
            "code_input": {
                "x_data": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
                "y_data": [0.0, 0.25, 1.0, 2.25, 4.0, 6.25, 9.0],
                "method": "savitzky_golay",
                "window_size": 5,
                "poly_order": 2,
            },
            "text_input": {
                "input_str": "0,0.5,1,1.5,2,2.5,3; 0,0.25,1,2.25,4,6.25,9; savitzky_golay; 5; 2",
            },
            "output": {
                "derivatives": [0.0, 0.5, 1.0, 2.0, 3.0, 4.5, 6.0],
                "x_at_derivative": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
                "second_derivatives": [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        x_data: List[float],
        y_data: List[float],
        method: str = "central",
        window_size: int = 5,
        poly_order: int = 2,
    ) -> Dict:
        if len(x_data) != len(y_data):
            raise ChemMCPError(f"x_data ({len(x_data)}) != y_data ({len(y_data)})")
        if len(x_data) < 2:
            raise ChemMCPError("Need at least 2 data points.")

        method = method.lower().strip()

        if method == "central":
            derivs = _central_diff(x_data, y_data)
            second = None
        elif method == "forward":
            derivs = _forward_diff(x_data, y_data)
            second = None
        elif method == "backward":
            derivs = _backward_diff(x_data, y_data)
            second = None
        elif method == "savitzky_golay":
            if window_size % 2 == 0:
                window_size += 1  # must be odd
            if window_size < 3:
                window_size = 3
            if poly_order >= window_size:
                poly_order = window_size - 1
            derivs, second = _savitzky_golay(x_data, y_data, window_size, poly_order)
        else:
            raise ChemMCPError(f"Unknown method: {method}. Use: central/forward/backward/savitzky_golay")

        result = {
            "derivatives": derivs,
            "x_at_derivative": list(x_data),
            "second_derivatives": second,
        }
        logger.info(f"Numerical differentiation ({method}): {len(derivs)} points computed")
        return result

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = [p.strip() for p in input_str.split(";")]
            if len(parts) < 2:
                raise ValueError("Need at least 2 parts: x_data; y_data")

            x_data = [float(v) for v in parts[0].split(",")]
            y_data = [float(v) for v in parts[1].split(",")]
            method = parts[2].strip().lower() if len(parts) > 2 else "central"
            window_size = int(parts[3]) if len(parts) > 3 else 5
            poly_order = int(parts[4]) if len(parts) > 4 else 2

            return self._run_base(x_data, y_data, method, window_size, poly_order)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
