import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class UncertaintyPropagator(BaseTool):
    """
    测量不确定度传递计算工具。
    基于 GUM（测量不确定度表示指南），通过偏导数法传递不确定度。
    """
    __version__ = "0.1.0"
    name = "UncertaintyPropagator"
    func_name = "propagate_uncertainty"
    description = "Propagate measurement uncertainty through mathematical expressions using GUM (Guide to Uncertainty in Measurement) methodology."
    implementation_description = "Uses numerical differentiation (central difference) to compute sensitivity coefficients (∂f/∂xi), then combines uncertainties: uc² = Σ(∂f/∂xi)²·u²(xi). Supports expanded uncertainty U=k·uc."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Uncertainty", "GUM", "Metrology", "QA/QC", "Statistics"]
    required_envs = []

    code_input_sig = [
        ("variables", "dict", "N/A", "Dict of variable names → {'value': float, 'uncertainty': float}. E.g., {'mass': {'value': 1.234, 'u': 0.002}, 'volume': {'value': 50.0, 'u': 0.05}}."),
        ("expression", "str", "N/A", "Mathematical expression as string using variable names. Supported ops: +, -, *, /, **, sqrt(), sin(), cos(), log(), log10(), exp(). Example: 'mass / volume * 1000'."),
        ("coverage_factor", "float", "2.0", "Coverage factor k for expanded uncertainty (k=2 for ~95% confidence)."),
        ("correlations", "dict", "", "Correlation matrix between variables: {('var1','var2'): r}. Default: all independent (r=0)."),
        ("perturbation", "float", "1e-8", "Perturbation size for numerical differentiation."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("result_value", "float", "Computed result value."),
        ("standard_uncertainty", "float", "Combined standard uncertainty uc."),
        ("expanded_uncertainty", "float", "Expanded uncertainty U = k × uc."),
        ("sensitivity_coefficients", "dict", "Sensitivity coefficients ∂f/∂xi for each variable."),
        ("uncertainty_contributions", "dict", "Each variable's contribution (% of uc²) and absolute contribution."),
        ("uncertainty_budget", "list", "Sorted list of contributions by magnitude."),
        ("relative_uncertainty", "float", "Relative standard uncertainty uc/|y| (%)."),
        ("coverage_factor_used", "float", "k used."),
        ("formula_detail", "str", "Detailed formula breakdown."),
    ]

    examples = [
        {
            "code_input": {
                "variables": {"mass": {"value": 0.1523, "u": 0.0002}, "volume": {"value": 25.00, "u": 0.03}},
                "expression": "mass / volume * 1000",
                "coverage_factor": 2.0,
            },
            "text_input": {"params_str": "see code input"},
            "output": {
                "result_value": 6.092,
                "expanded_uncertainty": 0.015,
            },
        },
    ]

    # Safe math environment for eval
    _SAFE_NAMES = {
        "sqrt": math.sqrt, "sin": math.sin, "cos": math.cos,
        "tan": math.tan, "exp": math.exp, "log": math.log,
        "log10": math.log10, "abs": abs, "pi": math.pi,
        "e": math.e,
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _eval_expr(self, expr: str, var_values: dict) -> float:
        """Safely evaluate a mathematical expression."""
        env = dict(self._SAFE_NAMES)
        env.update(var_values)
        try:
            result = eval(expr, {"__builtins__": {}}, env)
        except Exception as e:
            raise ChemMCPError(f"Failed to evaluate expression '{expr}': {e}")
        return float(result)

    def _run_base(
        self,
        variables: Dict[str, Dict[str, float]],
        expression: str,
        coverage_factor: float = 2.0,
        correlations: Optional[Dict[tuple, float]] = None,
        perturbation: float = 1e-8,
    ) -> dict:
        """Core logic: propagate uncertainty through expression."""
        if not variables:
            raise ChemMCPError("No variables provided.")
        if not expression or not expression.strip():
            raise ChemMCPError("Expression is empty.")

        corr = correlations or {}

        # Compute nominal value
        var_values = {k: v["value"] for k, v in variables.items()}
        y0 = self._eval_expr(expression.strip(), var_values)

        # Numerical partial derivatives (central difference)
        sens_coeffs: Dict[str, float] = {}
        for name, vinfo in variables.items():
            x0 = vinfo["value"]
            u = vinfo.get("u", 0.0)
            h = max(abs(x0) * perturbation, perturbation) if x0 != 0 else perturbation

            var_plus = dict(var_values)
            var_plus[name] = x0 + h
            var_minus = dict(var_values)
            var_minus[name] = x0 - h

            f_plus = self._eval_expr(expression.strip(), var_plus)
            f_minus = self._eval_expr(expression.strip(), var_minus)
            sens_coeffs[name] = (f_plus - f_minus) / (2 * h)

        # Combined variance: uc² = Σ(∂f/∂xi)²ui² + 2ΣΣ(∂f/∂xi)(∂f/∂xj)r_ij ui uj
        var_names = list(variables.keys())
        uc_sq = 0.0
        contrib: Dict[str, Dict[str, float]] = {}

        # Initialize all contrib entries first
        for name in var_names:
            ci = sens_coeffs[name]
            ui = variables[name].get("u", 0.0)
            contrib[name] = {"absolute": ci * ui, "squared": (ci * ui) ** 2, "percent": 0.0}

        # Calculate combined variance with correlations
        for i, name_i in enumerate(var_names):
            ci = sens_coeffs[name_i]
            ui = variables[name_i].get("u", 0.0)
            for j, name_j in enumerate(var_names):
                cj = sens_coeffs[name_j]
                uj = variables[name_j].get("u", 0.0)
                key = (name_i, name_j)
                r = corr.get(key, 1.0 if i == j else 0.0)
                if i != j:
                    # Cross terms: add half to each
                    cross_term = r * ci * ui * cj * uj
                    contrib[name_i]["squared"] += cross_term
                # when i==j, diagonal already initialized above

        for name in var_names:
            uc_sq += contrib[name]["squared"]

        uc = math.sqrt(max(uc_sq, 0.0))
        U = coverage_factor * uc
        rel_u = (uc / abs(y0) * 100) if abs(y0) > 1e-30 else float('inf')

        # Percentage contributions
        total_uc_sq = uc_sq if uc_sq > 0 else 1e-30
        for name in contrib:
            # recalc individual squared contribution (diagonal only for % display)
            ci = sens_coeffs[name]
            ui = variables[name].get("u", 0.0)
            diag_contrib = (ci * ui) ** 2
            contrib[name]["percent"] = (diag_contrib / total_uc_sq) * 100 if total_uc_sq > 0 else 0.0

        # Budget sorted by |contribution|
        budget = sorted(
            [{"variable": n, "sensitivity_coefficient": round(sens_coeffs[n], 6),
              "uncertainty_u": variables[n].get("u", 0),
              "contribution_percent": round(contrib[n]["percent"], 4),
              "absolute_contribution": round(abs(sens_coeffs[n] * variables[n].get("u", 0)), 8)}
             for n in var_names],
            key=lambda x: abs(x["contribution_percent"]),
            reverse=True,
        )

        # Formula detail
        terms = [f"({c:.4g}×{variables[n]['u']:.4g})²" for n, c in sens_coeffs.items()]
        formula_detail = (
            f"f = {expression}\n"
            f"y = {y0:.6g}\n"
            f"uc² = {' + '.join(terms)}\n"
            f"uc = √{uc_sq:.6g} = {uc:.6g}\n"
            f"U ({coverage_factor:.1f}σ) = {U:.6g}"
        )

        logger.info(f"Uncertainty propagation: y={y0:.6g}, uc={uc:.6g}, U={U:.6g} (k={coverage_factor})")
        return {
            "result_value": round(y0, 10),
            "standard_uncertainty": round(uc, 10),
            "expanded_uncertainty": round(U, 10),
            "sensitivity_coefficients": {n: round(c, 10) for n, c in sens_coeffs.items()},
            "uncertainty_contributions": {n: {
                "absolute": round(abs(c * variables[n].get("u", 0)), 10),
                "percent": round(contrib[n]["percent"], 4),
            } for n, c in sens_coeffs.items()},
            "uncertainty_budget": budget,
            "relative_uncertainty": round(rel_u, 6),
            "coverage_factor_used": coverage_factor,
            "formula_detail": formula_detail,
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
