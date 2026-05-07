import logging
from typing import Dict, List, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _cubic_spline_coefficients(x: List[float], y: List[float]):
    """
    Compute natural cubic spline coefficients.
    Returns: (a_coeffs, b_coeffs, c_coeffs, d_coeffs, x_nodes)
    where spline on [x[i], x[i+1]] = a_i + b_i*(x-xi) + c_i*(x-xi)^2 + d_i*(x-xi)^3
    """
    n = len(x) - 1  # number of intervals

    try:
        import numpy as np
        return _cubic_spline_numpy(x, y)
    except ImportError:
        pass

    # Pure Python tridiagonal solver for natural spline (S''(x0)=S''(xn)=0)
    h = [x[i + 1] - x[i] for i in range(n)]

    # Build tridiagonal system for c coefficients (second derivatives at nodes)
    # h_{i-1}*c_{i-1} + 2*(h_{i-1}+h_i)*c_i + h_i*c_{i+1} = 3*((y_{i+1}-y_i)/h_i - (y_i-y_{i-1})/h_{i-1})
    alpha = [0.0] * (n + 1)
    l_list = [1.0] * (n + 1)
    mu_list = [0.0] * (n + 1)
    z_list = [0.0] * (n + 1)

    for i in range(1, n):
        alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1])

    for i in range(1, n):
        l_list[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu_list[i - 1]
        if abs(l_list[i]) < 1e-15:
            l_list[i] = 1e-15
        mu_list[i] = h[i] / l_list[i]
        z_list[i] = (alpha[i] - h[i - 1] * z_list[i - 1]) / l_list[i]

    l_list[n] = 1.0
    z_list[n] = 0.0
    c = [0.0] * (n + 1)
    b = [0.0] * n
    d = [0.0] * n
    a = list(y)

    for j in range(n - 1, -1, -1):
        c[j] = z_list[j] - mu_list[j] * c[j + 1]
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    return a, b, c[:n], d, x


def _cubic_spline_numpy(x, y):
    """Use scipy/numpy for cubic spline."""
    try:
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(x, y, bc_type='natural')
        return ("numpy_scipy", cs, x, y)
    except ImportError:
        pass
    try:
        import numpy as np
        # Manual numpy implementation
        n = len(x) - 1
        h = np.diff(x)
        # Build tridiagonal matrix
        A = np.zeros((n + 1, n + 1))
        rhs = np.zeros(n + 1)
        for i in range(1, n):
            A[i, i - 1] = h[i - 1]
            A[i, i] = 2 * (h[i - 1] + h[i])
            A[i, i + 1] = h[i]
            rhs[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

        # Natural boundary conditions
        A[0, 0] = 1.0
        A[n, n] = 1.0

        c = np.linalg.solve(A, rhs)
        a = np.array(y, dtype=float)
        b = np.zeros(n)
        d_vec = np.zeros(n)

        for j in range(n):
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (2 * c[j] + c[j + 1]) / 3
            d_vec[j] = (c[j + 1] - c[j]) / (3 * h[j])

        return (list(a), list(b), list(c[:n]), list(d_vec), x)
    except ImportError:
        return _cubic_spline_coefficients(x, y)


def _eval_spline(spline_data, x_query: float) -> tuple:
    """Evaluate cubic spline at a point, return (value, slope)."""
    if isinstance(spline_data[0], str) and spline_data[0] == "numpy_scipy":
        _, cs, _, _ = spline_data
        val = float(cs(x_query))
        # Derivative
        try:
            slope = float(cs(x_query, 1))
        except Exception:
            slope = None
        return val, slope

    a, b, c, d, x_nodes = spline_data
    n = len(x_nodes) - 1

    # Find interval
    if x_query <= x_nodes[0]:
        idx = 0
    elif x_query >= x_nodes[-1]:
        idx = n - 1
    else:
        idx = 0
        for i in range(n):
            if x_nodes[i] <= x_query <= x_nodes[i + 1]:
                idx = i
                break

    dx = x_query - x_nodes[idx]
    val = a[idx] + b[idx] * dx + c[idx] * dx ** 2 + d[idx] * dx ** 3
    slope = b[idx] + 2 * c[idx] * dx + 3 * d[idx] * dx ** 2
    return round(val, 8), round(slope, 8)


@ChemMCPManager.register_tool
class InterpolationSpline(BaseTool):
    """
    样条插值工具 —— 热力学数据表内插。
    使用三次样条插值对表格数据进行内插和外推。
    """
    __version__ = "0.1.0"
    name = "InterpolationSpline"
    func_name = "interpolation_spline"
    description = (
        "Perform cubic spline interpolation on tabulated data. "
        "Essential for interpolating thermodynamic property tables."
    )
    implementation_description = (
        "Implements natural cubic spline interpolation with optional scipy backend. "
        "Returns interpolated values and slopes at query points."
    )
    oss_dependencies = [
        ("scipy", "https://scipy.org/", "BSD-3-Clause"),
        ("numpy", "https://numpy.org/", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Interpolation", "Thermodynamics", "Data Processing", "Spline"]
    required_envs = []

    code_input_sig = [
        ("x_known", "list", "N/A", "Known x values (must be sorted ascending)."),
        ("y_known", "list", "N/A", "Known y values corresponding to each x."),
        ("x_query", "float_or_list", "N/A", "Single float or list of floats to interpolate at."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Format: 'x1,x2,x3; y1,y2,y3; xq1,xq2,...' "
         "Example: '0,1,2,4; 0,1,8,16; 0.5,1.5,3'"),
    ]

    output_sig = [
        ("y_interpolated", "float or list", "Interpolated y value(s) at query point(s)."),
        ("slope_at_point", "float or list or None", "First derivative dy/dx at query point(s)."),
    ]

    examples = [
        {
            "code_input": {
                "x_known": [0.0, 1.0, 2.0, 4.0],
                "y_known": [0.0, 1.0, 8.0, 16.0],
                "x_query": 1.5,
            },
            "text_input": {
                "input_str": "0,1,2,4; 0,1,8,16; 1.5",
            },
            "output": {
                "y_interpolated": 3.375,
                "slope_at_point": 6.5,
            },
        },
        {
            "code_input": {
                "x_known": [0.0, 1.0, 2.0, 4.0],
                "y_known": [0.0, 1.0, 8.0, 16.0],
                "x_query": [0.5, 1.5, 3.0],
            },
            "text_input": {
                "input_str": "0,1,2,4; 0,1,8,16; 0.5,1.5,3",
            },
            "output": {
                "y_interpolated": [0.25, 3.375, 12.25],
                "slope_at_point": [0.5, 6.5, 4.5],
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        x_known: List[float],
        y_known: List[float],
        x_query: Union[float, List[float]],
    ) -> Dict:
        if len(x_known) != len(y_known):
            raise ChemMCPError(f"x_known ({len(x_known)}) and y_known ({len(y_known)}) must have same length.")
        if len(x_known) < 3:
            raise ChemMCPError("Need at least 3 points for cubic spline interpolation.")

        # Sort by x
        paired = sorted(zip(x_known, y_known), key=lambda p: p[0])
        xs = [p[0] for p in paired]
        ys = [p[1] for p in paired]

        spline_data = _cubic_spline_coefficients(xs, ys)

        is_single = isinstance(x_query, (int, float))
        queries = [float(x_query)] if is_single else [float(q) for q in x_query]

        y_interp = []
        slopes = []
        for q in queries:
            val, slope = _eval_spline(spline_data, q)
            y_interp.append(val)
            slopes.append(slope)

        result = {
            "y_interpolated": y_interp[0] if is_single else y_interp,
            "slope_at_point": slopes[0] if is_single else slopes,
        }
        logger.info(f"Spline interpolation at {x_query}: {result}")
        return result

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = [p.strip() for p in input_str.split(";")]
            if len(parts) < 3:
                raise ValueError("Need 3 semicolon-separated parts: x_known; y_known; x_query")

            x_known = [float(v) for v in parts[0].split(",")]
            y_known = [float(v) for v in parts[1].split(",")]
            q_parts = parts[2].strip().split(",")
            if len(q_parts) == 1 and "." in parts[2] or parts[2].strip().lstrip("-").replace(".", "").isdigit():
                x_query = float(parts[2].strip())
            else:
                x_query = [float(v) for v in q_parts]

            return self._run_base(x_known, y_known, x_query)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
