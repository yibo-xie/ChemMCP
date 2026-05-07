import json
import logging
import math
import re
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ErrorPropagation(BaseTool):
    """
    误差传播计算工具。
    基于一阶泰勒展开（偏导数法）进行不确定度传播分析，支持基本运算和常见函数。
    """
    __version__ = "0.1.0"
    name = "ErrorPropagation"
    func_name = "propagate_error"
    description = "Error/uncertainty propagation calculation for experimental data using first-order Taylor expansion (partial derivative method)."
    implementation_description = "Implements general error propagation via partial derivatives. Parses mathematical expressions and computes propagated uncertainty for independent (and optionally correlated) variables. Supports +, -, *, /, ^, exp, log, ln, sin, cos, tan, sqrt, abs."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Uncertainty", "Error Analysis", "Experimental Data", "Statistics", "Chemometrics"]
    required_envs = []

    code_input_sig = [
        ("expression", "str", "N/A", "Mathematical expression relating variables (e.g., 'w * h', 'm / V', 'x**2 * y')."),
        ("variables", "dict", "N/A", "Dict of variable name -> {'value': float, 'uncertainty': float}. E.g., {'w': {'value':10, 'uncertainty':0.1}, ...}"),
        ("correlation_matrix", "dict", "{}", "Optional correlation matrix for correlated variables. Keys like '(var1,var2)': corr_coeff."),
        ("absolute_tolerance", "float", "1e-8", "Numerical step size for central difference partial derivatives."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with expression, variables, and optional correlation_matrix."),
    ]

    output_sig = [
        ("value", "float", "Computed value of the expression."),
        ("uncertainty", "float", "Propagated absolute uncertainty (standard uncertainty)."),
        ("relative_uncertainty", "float", "Relative uncertainty (Δy/|y|) as fraction."),
        ("percent_uncertainty", "float", "Relative uncertainty as percentage."),
        ("partial_derivatives", "dict", "Partial derivative ∂f/∂xi for each variable at the nominal point."),
        ("contribution_dict", "dict", "Variance contribution (∂f/∂xi · σi)² from each variable."),
        ("sensitivity_coefficients", "dict", "Sensitivity coefficient (∂f/∂xi · σi / σf) normalized."),
        ("expression", "str", "The input expression."),
        ("n_variables", "int", "Number of variables in the expression."),
    ]

    examples = [
        {
            "code_input": {
                "expression": "w * h",
                "variables": {
                    "w": {"value": 10.0, "uncertainty": 0.1},
                    "h": {"value": 5.0, "uncertainty": 0.05},
                },
            },
            "text_input": {
                "params_str": '{"expression":"w*h","variables":{"w":{"value":10,"uncertainty":0.1},"h":{"value":5,"uncertainty":0.05}}}'
            },
            "output": {
                "value": 50.0,
                "uncertainty": 1.118033988749895,
                "relative_uncertainty": 0.02236,
            },
        },
        {
            "code_input": {
                "expression": "m / V",
                "variables": {
                    "m": {"value": 100.0, "uncertainty": 0.5},
                    "V": {"value": 50.0, "uncertainty": 0.2},
                },
            },
            "text_input": {
                "params_str": '{"expression":"m/V","variables":{"m":{"value":100,"uncertainty":0.5},"V":{"value":50,"uncertainty":0.2}}}'
            },
            "output": {
                "value": 2.0,
                "uncertainty": 0.014142135623730951,
            },
        },
        {
            "code_input": {
                "expression": "x ** 2",
                "variables": {
                    "x": {"value": 5.0, "uncertainty": 0.1},
                },
            },
            "text_input": {"params_str": '{"expression":"x**2","variables":{"x":{"value":5,"uncertainty":0.1}}}'},
            "output": {
                "value": 25.0,
                "uncertainty": 1.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- safe expression evaluator ----

    _SAFE_NAMES = {
        "sin": math.sin, "cos": math.cos, "tan": math.tan,
        "exp": math.exp, "log": math.log10, "ln": math.log,
        "sqrt": math.sqrt, "abs": abs,
        "pi": math.pi, "e": math.e,
        "asin": math.asin, "acos": math.acos, "atan": math.atan,
        "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    }

    @staticmethod
    def _eval_expr(expr: str, var_values: Dict[str, float]) -> float:
        """Safely evaluate a math expression with given variable values."""
        ns = dict(ErrorPropagation._SAFE_NAMES)
        ns.update(var_values)
        # Only allow alphanumeric, underscore, operators, dots, parentheses
        expr_clean = expr.strip()
        if re.match(r'^[\w\s\+\-\*\/\^\(\)\.\,\[\]]+$', expr_clean.replace("**", "")):
            try:
                # Replace ^ with ** for Python eval
                expr_py = expr_clean.replace("^", "**")
                return eval(expr_py, {"__builtins__": {}}, ns)
            except Exception as e:
                raise ChemMCPError(f"Failed to evaluate expression '{expr}': {e}")
        else:
            raise ChemMCPError(f"Expression contains disallowed characters: '{expr}'")

    @staticmethod
    def _partial_derivative(expr: str, var: str, var_values: Dict[str, float], h: float = 1e-8) -> float:
        """Central difference partial derivative ∂f/∂var."""
        values_up = dict(var_values)
        values_dn = dict(var_values)
        x0 = var_values.get(var, 0.0)
        step = max(h * abs(x0) if x0 != 0 else h, 1e-15)
        values_up[var] = x0 + step
        values_dn[var] = x0 - step
        f_plus = ErrorPropagation._eval_expr(expr, values_up)
        f_minus = ErrorPropagation._eval_expr(expr, values_dn)
        return (f_plus - f_minus) / (2.0 * step)

    def _run_base(
        self,
        expression: str,
        variables: Dict[str, Dict[str, float]],
        correlation_matrix: Optional[Dict[str, float]] = None,
        absolute_tolerance: float = 1e-8,
    ) -> dict:
        """Core logic: propagate uncertainties through an expression."""
        if not expression or not expression.strip():
            raise ChemMCPError("Expression cannot be empty.")
        if not variables:
            raise ChemMCPError("Variables dict cannot be empty.")

        corrs = correlation_matrix or {}

        # Build nominal value dict
        nominals = {}
        for vname, vinfo in variables.items():
            if not isinstance(vinfo, dict):
                raise ChemMCPError(f"Variable '{vname}' must be a dict with 'value' and 'uncertainty' keys.")
            nominals[vname] = vinfo["value"]

        # Compute nominal function value
        f0 = self._eval_expr(expression, nominals)

        # Compute partial derivatives and variance contributions
        partials = {}
        contributions = {}
        var_names = list(variables.keys())

        for vname in var_names:
            sigma = variables[vname].get("uncertainty", 0.0)
            df_dx = self._partial_derivative(expression, vname, nominals, h=absolute_tolerance)
            partials[vname] = round(df_dx, 12)
            contributions[vname] = (df_dx * sigma) ** 2

        # Total variance (including covariance terms if correlations provided)
        total_var = sum(contributions.values())
        for key, rho in corrs.items():
            match = re.match(r'\((\w+),(\w+)\)', key)
            if match:
                vi, vj = match.group(1), match.group(2)
                if vi in partials and vj in partials:
                    si = variables[vi].get("uncertainty", 0.0)
                    sj = variables[vj].get("uncertainty", 0.0)
                    total_var += 2.0 * rho * partials[vi] * si * partials[vj] * sj

        sigma_f = math.sqrt(max(total_var, 0.0))
        rel_sigma = abs(sigma_f / f0) if abs(f0) > 1e-30 else (float('inf') if sigma_f > 0 else 0.0)

        # Sensitivity coefficients (normalized)
        sensitivity = {}
        if sigma_f > 1e-30:
            for vname in var_names:
                sensitivity[vname] = round(contributions[vname] ** 0.5 / sigma_f, 6)

        logger.info(f"ErrorPropagation: f={f0} ± {sigma_f} ({rel_sigma*100:.4f}%)")
        return {
            "value": f0,
            "uncertainty": round(sigma_f, 12),
            "relative_uncertainty": round(rel_sigma, 12),
            "percent_uncertainty": round(rel_sigma * 100, 8),
            "partial_derivatives": {k: round(v, 12) for k, v in partials.items()},
            "contribution_dict": {k: round(v, 14) for k, v in contributions.items()},
            "sensitivity_coefficients": sensitivity,
            "expression": expression.strip(),
            "n_variables": len(var_names),
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input. Expected JSON string with keys: expression, variables, ...")
        return self._run_base(**kwargs)
