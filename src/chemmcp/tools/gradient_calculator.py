import logging
import math
from typing import Dict, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

from .partial_derivative import _safe_eval

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GradientCalculator(BaseTool):
    """
    梯度计算工具 —— 势能面分析、几何优化。
    计算多元函数的梯度向量（所有一阶偏导数）。
    """
    __version__ = "0.1.0"
    name = "GradientCalculator"
    func_name = "gradient_calculator"
    description = (
        "Compute gradient vector of a multivariable function numerically. "
        "Essential for potential energy surface analysis and geometry optimization."
    )
    implementation_description = (
        "Computes ∇f = (∂f/∂x₁, ∂f/∂x₂, ..., ∂f/∂xₙ) using central differences "
        "for each variable component."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Methods", "Optimization", "Potential Energy Surface", "Gradient"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Mathematical expression string, e.g., 'x**2 + y**2'."),
        ("variables", "list", "N/A", "List of variable names, e.g., ['x', 'y']."),
        ("eval_point", "dict", "N/A", "Dictionary mapping variable names to float values."),
        ("step_size", "float", "1e-5", "Step size for numerical differentiation (default: 1e-5)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Space-separated: 'func_expr var1,var2,... val1,val2,... [step_size]'. "
         "Example: 'x**2+y**2 x,y 3,4'"),
    ]

    output_sig = [
        ("gradient", "dict", "Dictionary mapping each variable name to its partial derivative value."),
        ("magnitude", "float", "Euclidean norm (magnitude) of the gradient vector."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "x**2 + y**2",
                "variables": ["x", "y"],
                "eval_point": {"x": 3.0, "y": 4.0},
                "step_size": 1e-5,
            },
            "text_input": {
                "input_str": "x**2+y**2 x,y 3,4",
            },
            "output": {
                "gradient": {"x": 6.0, "y": 8.0},
                "magnitude": 10.0,
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
        step_size: float = 1e-5,
    ) -> Dict:
        h = step_size
        grad = {}
        for var in variables:
            point_plus = dict(eval_point)
            point_minus = dict(eval_point)
            point_plus[var] += h
            point_minus[var] -= h
            f_plus = _safe_eval(func_expr, point_plus)
            f_minus = _safe_eval(func_expr, point_minus)
            grad[var] = round((f_plus - f_minus) / (2.0 * h), 8)

        magnitude = round(math.sqrt(sum(v ** 2 for v in grad.values())), 8)
        logger.info(f"Gradient at {eval_point}: {grad}, magnitude={magnitude}")
        return {
            "gradient": grad,
            "magnitude": magnitude,
        }

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = input_str.strip().split()
            if len(parts) < 3:
                raise ValueError("Need at least 3 parts: func_expr variables values")

            func_expr = parts[0]
            variables = parts[1].split(",")
            values = [float(v) for v in parts[2].split(",")]
            step_size = float(parts[3]) if len(parts) > 3 else 1e-5

            if len(variables) != len(values):
                raise ValueError(f"Mismatch: {len(variables)} variables but {len(values)} values")

            eval_point = dict(zip(variables, values))
            return self._run_base(func_expr, variables, eval_point, step_size)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
