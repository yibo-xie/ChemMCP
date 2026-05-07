import logging
from typing import Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

from .partial_derivative import _safe_eval

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BisectionMethod(BaseTool):
    """
    二分法求根工具 —— pH计算、平衡常数求解。
    使用二分法在给定区间内寻找连续函数的根。
    """
    __version__ = "0.1.0"
    name = "BisectionMethod"
    func_name = "bisection_method"
    description = (
        "Find roots of continuous functions using bisection method. "
        "Robust and guaranteed to converge for continuous functions on an interval "
        "where f(a) and f(b) have opposite signs. Useful for pH calculations "
        "and equilibrium constant problems."
    )
    implementation_description = (
        "Implements bisection method: repeatedly halves interval [a,b], selecting "
        "the subinterval where sign change occurs. Convergence is linear but guaranteed."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Methods", "Root Finding", "Equilibrium", "pH Calculation"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Function expression f(x) to find root of."),
        ("variable", "str", "N/A", "Variable name in the expression, e.g., 'x'."),
        ("a", "float", "N/A", "Left bound of the search interval (f(a) should be nonzero)."),
        ("b", "float", "N/A", "Right bound of the search interval (f(b) should be nonzero)."),
        ("tolerance", "float", "1e-8", "Convergence tolerance for interval width (default: 1e-8)."),
        ("max_iterations", "int", "100", "Maximum number of iterations (default: 100)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Space-separated: 'func_expr variable a b [tolerance] [max_iter]'. "
         "Example: 'x**3-x-2 x 1 2'"),
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
                "a": 1.0,
                "b": 2.0,
                "tolerance": 1e-8,
                "max_iterations": 100,
            },
            "text_input": {
                "input_str": "x**3-x-2 x 1 2",
            },
            "output": {
                "root": 1.5213797,
                "iterations": 27,
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
        a: float,
        b: float,
        tolerance: float = 1e-8,
        max_iterations: int = 100,
    ) -> Dict:
        fa = _safe_eval(func_expr, {variable: a})
        fb = _safe_eval(func_expr, {variable: b})

        if fa * fb > 0:
            raise ChemMCPError(
                f"f(a)={fa} and f(b)={fb} have same sign. Bisection requires opposite signs."
            )

        for i in range(1, max_iterations + 1):
            c = (a + b) / 2.0
            fc = _safe_eval(func_expr, {variable: c})

            if abs(fc) < tolerance or abs(b - a) / 2 < tolerance:
                return {
                    "root": round(c, 10),
                    "iterations": i,
                    "converged": True,
                    "f_at_root": round(fc, 12),
                }

            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc

        # Max iterations reached
        root = (a + b) / 2.0
        f_final = _safe_eval(func_expr, {variable: root})
        return {
            "root": round(root, 10),
            "iterations": max_iterations,
            "converged": False,
            "f_at_root": round(f_final, 12),
        }

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = input_str.strip().split()
            if len(parts) < 4:
                raise ValueError("Need at least 4 parts: func_expr variable a b")

            func_expr = parts[0]
            variable = parts[1]
            a_val = float(parts[2])
            b_val = float(parts[3])
            tolerance = float(parts[4]) if len(parts) > 4 else 1e-8
            max_iterations = int(parts[5]) if len(parts) > 5 else 100

            return self._run_base(func_expr, variable, a_val, b_val, tolerance, max_iterations)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
