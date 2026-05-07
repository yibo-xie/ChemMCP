import logging
import math
from typing import Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

from .partial_derivative import _safe_eval

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NewtonRaphsonSolver(BaseTool):
    """
    牛顿-拉夫森法求根工具 —— 非线性方程求根。
    使用牛顿迭代法求解非线性方程 f(x) = 0 的根。
    """
    __version__ = "0.1.0"
    name = "NewtonRaphsonSolver"
    func_name = "newton_raphson_solver"
    description = (
        "Find roots of nonlinear equations using Newton-Raphson iteration method. "
        "Widely used for solving equilibrium equations and pH calculations."
    )
    implementation_description = (
        "Implements Newton-Raphson iteration: x_{n+1} = x_n - f(x_n)/f'(x_n). "
        "Uses central difference for derivative evaluation."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Methods", "Root Finding", "Nonlinear Equations", "Newton-Raphson"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Function expression f(x) to find root of (set equal to zero)."),
        ("variable", "str", "N/A", "Variable name in the expression, e.g., 'x'."),
        ("initial_guess", "float", "N/A", "Initial guess for the root."),
        ("tolerance", "float", "1e-8", "Convergence tolerance (default: 1e-8)."),
        ("max_iterations", "int", "100", "Maximum number of iterations (default: 100)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Space-separated: 'func_expr variable initial_guess [tolerance] [max_iter]'. "
         "Example: 'x**3-x-2 x 1.5'"),
    ]

    output_sig = [
        ("root", "float", "The approximated root value."),
        ("iterations", "int", "Number of iterations performed."),
        ("converged", "bool", "Whether the method converged within tolerance."),
        ("f_at_root", "float", "Function value at the found root (should be ≈0)."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "x**3 - x - 2",
                "variable": "x",
                "initial_guess": 1.5,
                "tolerance": 1e-8,
                "max_iterations": 100,
            },
            "text_input": {
                "input_str": "x**3-x-2 x 1.5",
            },
            "output": {
                "root": 1.5213797,
                "iterations": 4,
                "converged": True,
                "f_at_root": 0.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        func_expr: str,
        variable: str,
        initial_guess: float,
        tolerance: float = 1e-8,
        max_iterations: int = 100,
    ) -> Dict:
        x = float(initial_guess)
        h = 1e-7

        for i in range(1, max_iterations + 1):
            try:
                f_x = _safe_eval(func_expr, {variable: x})
            except Exception as e:
                raise ChemMCPError(f"Failed to evaluate function at x={x}: {e}")

            if abs(f_x) < tolerance:
                return {
                    "root": round(x, 10),
                    "iterations": i - 1,
                    "converged": True,
                    "f_at_root": round(f_x, 12),
                }

            # Derivative via central difference
            try:
                fp_x = (_safe_eval(func_expr, {variable: x + h}) -
                        _safe_eval(func_expr, {variable: x - h})) / (2 * h)
            except Exception:
                raise ChemMCPError("Failed to evaluate derivative.")

            if abs(fp_x) < 1e-15:
                raise ChemMCPError(f"Derivative near zero at x={x}. Cannot continue Newton-Raphson.")

            x_new = x - f_x / fp_x

            if abs(x_new - x) < tolerance:
                f_final = _safe_eval(func_expr, {variable: x_new})
                return {
                    "root": round(x_new, 10),
                    "iterations": i,
                    "converged": True,
                    "f_at_root": round(f_final, 12),
                }
            x = x_new

        # Did not converge
        f_final = _safe_eval(func_expr, {variable: x})
        return {
            "root": round(x, 10),
            "iterations": max_iterations,
            "converged": False,
            "f_at_root": round(f_final, 12),
        }

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = input_str.strip().split()
            if len(parts) < 3:
                raise ValueError("Need at least 3 parts: func_expr variable initial_guess")

            func_expr = parts[0]
            variable = parts[1]
            initial_guess = float(parts[2])
            tolerance = float(parts[3]) if len(parts) > 3 else 1e-8
            max_iterations = int(parts[4]) if len(parts) > 4 else 100

            return self._run_base(func_expr, variable, initial_guess, tolerance, max_iterations)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
