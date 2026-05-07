import logging
import json
import math
import numpy as np
from scipy import integrate

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NumericalIntegrator(BaseTool):
    """
    数值积分工具，用于轨道重叠积分、热力学函数计算等。
    支持多种数值积分方法：梯形法、辛普森法、高斯求积、自适应积分。
    """
    __version__ = "0.1.0"
    name = "NumericalIntegrator"
    func_name = "numerical_integrate"
    description = "Numerical integration for orbital overlap integrals, thermodynamic functions, partition functions, and general chemistry calculations."
    implementation_description = "Supports trapezoidal rule, Simpson's rule, Gaussian quadrature (via scipy), and adaptive quadrature. Accepts function expressions as strings or discrete data points."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
        ("scipy", "https://scipy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Integration", "Thermodynamics", "Orbital Overlap", "Quadrature"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Integration mode: 'function' (from expression) or 'data' (from discrete points)."),
        ("method", "str", "simpson", "Method: 'trapezoidal', 'simpson', 'gaussian', or 'adaptive'."),
        # For mode=function:
        ("func_expr", "str", "x**2", "Function expression in variable x (e.g., 'x**2', 'sin(x)', 'exp(-x**2)')."),
        ("lower_bound", "float", "0.0", "Lower integration limit."),
        ("upper_bound", "float", "1.0", "Upper integration limit."),
        ("n_points", "int", "1000", "Number of points for trapezoidal/Simpson methods."),
        # For mode=data:
        ("x_data", "list", "null", "X coordinates for data mode (list of floats)."),
        ("y_data", "list", "null", "Y values for data mode (list of floats)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("integral_value", "float", "The computed integral value."),
        ("method_used", "str", "The method actually used."),
        ("n_evaluations", "int", "Number of function evaluations performed."),
        ("estimated_error", "float or null", "Estimated error for adaptive/adaptive methods."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "function",
                "method": "simpson",
                "func_expr": "x**2",
                "lower_bound": 0.0,
                "upper_bound": 1.0,
                "n_points": 100,
            },
            "text_input": {
                "input_str": '{"mode":"function","method":"simpson","func_expr":"x**2","lower_bound":0,"upper_bound":1,"n_points":100}',
            },
            "output": {
                "integral_value": 0.333333,
                "method_used": "simpson",
                "n_evaluations": 101,
                "estimated_error": None,
            },
        },
        {
            "code_input": {
                "mode": "function",
                "method": "gaussian",
                "func_expr": "sin(x)",
                "lower_bound": 0.0,
                "upper_bound": math.pi,
                "n_points": 50,
            },
            "text_input": {
                "input_str": '{"mode":"function","method":"gaussian","func_expr":"sin(x)","lower_bound":0,"upper_bound":3.14159}',
            },
            "output": {
                "integral_value": 2.0,
                "method_used": "gaussian",
                "n_evaluations": 50,
                "estimated_error": None,
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _eval_func(self, expr: str, x: float) -> float:
        """Safely evaluate a mathematical expression."""
        # Allow only safe mathematical functions
        safe_dict = {
            "sin": math.sin, "cos": math.cos, "tan": math.tan,
            "exp": math.exp, "log": math.log, "log10": math.log10,
            "sqrt": math.sqrt, "abs": abs, "pi": math.pi, "e": math.e,
            "asin": math.asin, "acos": math.acos, "atan": math.atan,
            "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
            "**": pow,
        }
        try:
            return eval(expr, {"__builtins__": {}}, {**safe_dict, "x": x})
        except Exception as e:
            raise ChemMCPError(f"Failed to evaluate expression '{expr}' at x={x}: {e}")

    def _run_base(
        self,
        mode: str = "function",
        method: str = "simpson",
        func_expr: str = "x**2",
        lower_bound: float = 0.0,
        upper_bound: float = 1.0,
        n_points: int = 1000,
        x_data: list = None,
        y_data: list = None,
    ) -> dict:
        """Core logic: numerical integration."""

        if mode == "data":
            if not x_data or not y_data:
                raise ChemMCPError("Data mode requires both x_data and y_data.")
            if len(x_data) != len(y_data):
                raise ChemMCPError("x_data and y_data must have the same length.")
            x_arr = np.array(x_data, dtype=float)
            y_arr = np.array(y_data, dtype=float)
            n_eval = len(x_arr)

            if method == "trapezoidal":
                result = float(np.trapz(y_arr, x_arr))
            elif method == "simpson":
                result = float(integrate.simpson(y_arr, x=x_arr))
            else:
                result = float(np.trapz(y_arr, x_arr))
                method = "trapezoidal"

            return {
                "integral_value": round(result, 6),
                "method_used": method,
                "n_evaluations": n_eval,
                "estimated_error": None,
            }

        elif mode == "function":
            a, b = lower_bound, upper_bound
            if a >= b:
                raise ChemMCPError("lower_bound must be less than upper_bound.")

            f_vec = np.vectorize(lambda x: self._eval_func(func_expr, x))

            if method == "trapezoidal":
                x_vals = np.linspace(a, b, max(n_points, 2))
                y_vals = f_vec(x_vals)
                result = float(np.trapz(y_vals, x_vals))
                n_eval = len(x_vals)

            elif method == "simpson":
                x_vals = np.linspace(a, b, max(n_points + 1, 3))  # Simpson needs odd number of points
                if len(x_vals) % 2 == 0:
                    x_vals = np.linspace(a, b, len(x_vals) + 1)
                y_vals = f_vec(x_vals)
                result = float(integrate.simpson(y_vals, x=x_vals))
                n_eval = len(x_vals)

            elif method == "gaussian":
                result, abs_err = integrate.fixed_quad(f_vec, a, b, n=max(n_points, 5))
                result = float(result)
                n_eval = n_points
                return {
                    "integral_value": round(result, 6),
                    "method_used": "gaussian",
                    "n_evaluations": n_eval,
                    "estimated_error": round(float(abs_err) if abs_err is not None else 0.0, 10),
                }

            elif method == "adaptive":
                result, abs_err = integrate.quad(f_vec, a, b)
                result = float(result)
                return {
                    "integral_value": round(result, 6),
                    "method_used": "adaptive",
                    "n_evaluations": -1,  # quad doesn't report this simply
                    "estimated_error": round(abs_err, 10),
                }

            else:
                raise ChemMCPError(f"Unknown method: {method}. Use 'trapezoidal', 'simpson', 'gaussian', or 'adaptive'.")

            logger.info(f"Integrated '{func_expr}' from {a} to {b} using {method}: {result}")
            return {
                "integral_value": round(result, 6),
                "method_used": method,
                "n_evaluations": n_eval,
                "estimated_error": None,
            }
        else:
            raise ChemMCPError(f"Unknown mode: {mode}. Use 'function' or 'data'.")

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            return self._run_base(**params)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
