import logging
import math
from typing import Dict, List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _safe_eval(expr: str, variables_dict: Dict[str, float]) -> float:
    """Safely evaluate a math expression with given variable bindings."""
    allowed_names = {
        "sin": math.sin, "cos": math.cos, "tan": math.tan,
        "exp": math.exp, "log": math.log, "log10": math.log10,
        "log2": math.log2, "sqrt": math.sqrt, "abs": abs,
        "pi": math.pi, "e": math.e,
        "asin": math.asin, "acos": math.acos, "atan": math.atan,
        "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
        "ceil": math.ceil, "floor": math.floor,
        "pow": pow, "max": max, "min": min,
    }
    allowed_names.update(variables_dict)
    try:
        return eval(expr, {"__builtins__": {}}, allowed_names)
    except Exception as e:
        raise ChemMCPError(f"Failed to evaluate expression '{expr}': {e}")


@ChemMCPManager.register_tool
class PartialDerivative(BaseTool):
    """
    偏导数计算工具 —— 热力学偏导关系。
    使用中心差分格式数值计算多元函数的偏导数。
    """
    __version__ = "0.1.0"
    name = "PartialDerivative"
    func_name = "partial_derivative"
    description = (
        "Compute partial derivatives of multivariable functions numerically "
        "using central difference scheme. Useful for thermodynamic partial "
        "derivative relations."
    )
    implementation_description = (
        "Uses central difference formula: ∂f/∂x ≈ (f(x+h)-f(x-h))/(2h). "
        "Supports arbitrary math expressions via safe eval."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Methods", "Thermodynamics", "Calculus", "Partial Derivative"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Mathematical expression string, e.g., 'x**2*y + sin(y)'."),
        ("variables", "list", "N/A", "List of variable names, e.g., ['x', 'y']."),
        ("eval_point", "dict", "N/A", "Dictionary mapping variable names to float values, e.g., {'x': 2.0, 'y': 3.14}."),
        ("var_to_diff", "str", "N/A", "Name of the variable to differentiate with respect to."),
        ("step_size", "float", "1e-5", "Step size for numerical differentiation (default: 1e-5)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Space-separated: 'func_expr var1,var2,... val1,val2,... var_to_diff [step_size]'. "
         "Example: 'x**2*y+sin(y) x,y 2,3.14 x 1e-5'"),
    ]

    output_sig = [
        ("result", "float", "The computed partial derivative value."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "x**2 * y + sin(y)",
                "variables": ["x", "y"],
                "eval_point": {"x": 2.0, "y": 3.14},
                "var_to_diff": "x",
                "step_size": 1e-5,
            },
            "text_input": {
                "input_str": "x**2*y+sin(y) x,y 2,3.14 x 1e-5",
            },
            "output": {
                "result": 12.56,
            },
        },
        {
            "code_input": {
                "func_expr": "x**2 * y + sin(y)",
                "variables": ["x", "y"],
                "eval_point": {"x": 2.0, "y": 3.14},
                "var_to_diff": "y",
                "step_size": 1e-5,
            },
            "text_input": {
                "input_str": "x**2*y+sin(y) x,y 2,3.14 y 1e-5",
            },
            "output": {
                "result": 4.9999967,  # ≈ 4 - cos(3.14) + 4 = ~5.0
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
        variables: List[str],
        eval_point: Dict[str, float],
        var_to_diff: str,
        step_size: float = 1e-5,
    ) -> float:
        if var_to_diff not in variables:
            raise ChemMCPError(f"Variable '{var_to_diff}' not found in variables list: {variables}")

        h = step_size
        point_plus = dict(eval_point)
        point_minus = dict(eval_point)
        point_plus[var_to_diff] += h
        point_minus[var_to_diff] -= h

        f_plus = _safe_eval(func_expr, point_plus)
        f_minus = _safe_eval(func_expr, point_minus)

        result = (f_plus - f_minus) / (2.0 * h)
        logger.info(f"Partial derivative of '{func_expr}' w.r.t. {var_to_diff} at {eval_point} = {result}")
        return round(result, 8)

    def _run_text(self, input_str: str) -> float:
        try:
            parts = input_str.strip().split()
            if len(parts) < 4:
                raise ValueError("Need at least 4 parts: func_expr variables values var_to_diff")

            func_expr = parts[0]
            variables = parts[1].split(",")
            values = [float(v) for v in parts[2].split(",")]
            var_to_diff = parts[3]
            step_size = float(parts[4]) if len(parts) > 4 else 1e-5

            if len(variables) != len(values):
                raise ValueError(f"Mismatch: {len(variables)} variables but {len(values)} values")

            eval_point = dict(zip(variables, values))
            return self._run_base(func_expr, variables, eval_point, var_to_diff, step_size)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
