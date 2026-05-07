import logging
import math
from typing import Dict, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

from .partial_derivative import _safe_eval

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LaplacianOperator(BaseTool):
    """
    拉普拉斯算符工具 —— 薛定谔方程、电荷密度分析。
    计算标量场的拉普拉斯算子 ∇²f = Σ ∂²f/∂xᵢ²。
    """
    __version__ = "0.1.0"
    name = "LaplacianOperator"
    func_name = "laplacian_operator"
    description = (
        "Compute the Laplacian (∇²f) of a scalar field function. "
        "Used in Schrödinger equation and charge density analysis."
    )
    implementation_description = (
        "Computes ∇²f = Σᵢ ∂²f/∂xᵢ² using central finite differences for each variable. "
        "Sum of second pure partial derivatives."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Methods", "Quantum Mechanics", "PDE", "Laplacian"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Mathematical expression string."),
        ("variables", "list", "N/A", "List of variable names."),
        ("eval_point", "dict", "N/A", "Dictionary mapping variable names to float values."),
        ("step_size", "float", "1e-5", "Step size for numerical differentiation (default: 1e-5)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Space-separated: 'func_expr var1,var2,... val1,val2,... [step_size]'. "
         "Example: 'x**2+y**2+z**2 x,y,z 1,2,3'"),
    ]

    output_sig = [
        ("result", "float", "The computed Laplacian value (∇²f)."),
        ("second_partials", "dict", "Dictionary of each variable's second partial derivative contribution."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "x**2 + y**2 + z**2",
                "variables": ["x", "y", "z"],
                "eval_point": {"x": 1.0, "y": 2.0, "z": 3.0},
                "step_size": 1e-5,
            },
            "text_input": {
                "input_str": "x**2+y**2+z**2 x,y,z 1,2,3",
            },
            "output": {
                "result": 6.0,
                "second_partials": {"x": 2.0, "y": 2.0, "z": 2.0},
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
        second_partials = {}
        total = 0.0

        for var in variables:
            p_plus = dict(eval_point); p_plus[var] += h
            p_minus = dict(eval_point); p_minus[var] -= h

            f_pp = _safe_eval(func_expr, p_plus)
            f_0 = _safe_eval(func_expr, eval_point)
            f_mm = _safe_eval(func_expr, p_minus)
            d2 = (f_pp - 2 * f_0 + f_mm) / (h * h)
            second_partials[var] = round(d2, 8)
            total += d2

        logger.info(f"Laplacian of '{func_expr}' at {eval_point} = {total}")
        return {
            "result": round(total, 8),
            "second_partials": second_partials,
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
