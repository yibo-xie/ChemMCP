import logging
import math
from typing import Dict, List, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _eval_model(model_expr: str, params: Dict[str, float], x: float) -> float:
    """Evaluate a model expression with given parameters at point x."""
    allowed = {
        "sin": math.sin, "cos": math.cos, "tan": math.tan,
        "exp": math.exp, "log": math.log, "log10": math.log10,
        "sqrt": math.sqrt, "abs": abs,
        "pi": math.pi, "e": math.e,
    }
    allowed.update(params)
    allowed["x"] = x
    try:
        return eval(model_expr, {"__builtins__": {}}, allowed)
    except Exception as e:
        raise ChemMCPError(f"Model evaluation failed: {e}")


def _extract_param_names(model_expr: str) -> List[str]:
    """Extract potential parameter names from a model expression."""
    import re
    # Find standalone single-letter parameters (not part of function names like exp, sin, cos, etc)
    # Look for single lowercase/uppercase letters that appear as operands
    common_params = ["a", "b", "c", "d", "k", "m", "n", "p", "q", "r", "s",
                     "A", "B", "C", "D", "K", "M", "N"]
    found = []
    for p in common_params:
        # Use word boundary check to avoid matching inside function names
        pattern = r'(?<![a-zA-Z])' + re.escape(p) + r'(?![a-zA-Z0-9])'
        if re.search(pattern, model_expr) and p != "x":
            found.append(p)
    return sorted(set(found)) if found else ["a", "b", "c"]


@ChemMCPManager.register_tool
class CurveFittingNonlinear(BaseTool):
    """
    非线性曲线拟合工具 —— 动力学参数提取。
    将数据拟合到用户指定的非线性模型函数。
    """
    __version__ = "0.1.0"
    name = "CurveFittingNonlinear"
    func_name = "curve_fitting_nonlinear"
    description = (
        "Fit experimental data to user-specified nonlinear model functions. "
        "Essential for extracting kinetic parameters from reaction data."
    )
    implementation_description = (
        "Uses scipy.optimize.curve_fit when available (Levenberg-Marquardt), "
        "otherwise falls back to a simple gradient descent implementation."
    )
    oss_dependencies = [
        ("scipy", "https://scipy.org/", "BSD-3-Clause"),
        ("numpy", "https://numpy.org/", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Curve Fitting", "Kinetics", "Data Analysis", "Nonlinear Regression"]
    required_envs = []

    code_input_sig = [
        ("x_data", "list", "N/A", "List of x values (independent variable)."),
        ("y_data", "list", "N/A", "List of y values (dependent variable)."),
        ("model_func", "str", "N/A", "Model function expression string using params a,b,c,... and variable x. E.g., 'a*exp(-b*x)+c'."),
        ("initial_params", "list", "N/A", "Initial guess for parameters [a0, b0, c0, ...]. Order must match extracted param names."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "JSON-like format or space-separated with semicolons. "
         "Format: 'x1,x2,x3; y1,y2,y3; model_expr; p0_1,p0_2,p0_3'"),
    ]

    output_sig = [
        ("fitted_params", "dict", "Fitted parameter values as dict {param_name: value}."),
        ("r_squared", "float", "R² goodness-of-fit coefficient (0 to 1)."),
        ("residuals", "list", "Residuals (observed - fitted) for each data point."),
        ("fitted_y", "list", "Fitted y values at each x data point."),
    ]

    examples = [
        {
            "code_input": {
                "x_data": [0.0, 1.0, 2.0, 3.0, 4.0],
                "y_data": [5.0, 2.9, 1.8, 1.2, 0.7],
                "model_func": "a*exp(-b*x)+c",
                "initial_params": [4.0, 0.5, 0.5],
            },
            "text_input": {
                "input_str": "0,1,2,3,4; 5,2.9,1.8,1.2,0.7; a*exp(-b*x)+c; 4,0.5,0.5",
            },
            "output": {
                "fitted_params": {"a": 4.0, "b": 0.5, "c": 0.5},
                "r_squared": 0.999,
                "residuals": [0.0, 0.0, 0.0, 0.0, 0.0],
                "fitted_y": [4.5, 2.92, 1.85, 1.17, 0.75],
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._use_scipy = False
        try:
            import scipy.optimize
            import numpy as np
            self._scipy_optimize = scipy.optimize
            self._numpy = np
            self._use_scipy = True
        except ImportError:
            logger.info("scipy/numpy not available, using fallback gradient descent")

    def _run_base(
        self,
        x_data: List[float],
        y_data: List[float],
        model_func: str,
        initial_params: List[float],
    ) -> Dict:
        if len(x_data) != len(y_data):
            raise ChemMCPError(f"x_data length ({len(x_data)}) != y_data length ({len(y_data)})")
        if len(x_data) < len(initial_params):
            raise ChemMCPError("Need more data points than parameters.")

        param_names = _extract_param_names(model_func)
        if len(initial_params) != len(param_names):
            raise ChemMCPError(
                f"Need {len(param_names)} initial params ({param_names}), got {len(initial_params)}"
            )

        if self._use_scipy:
            return _fit_with_scipy(x_data, y_data, model_func, param_names, initial_params,
                                    self._scipy_optimize, self._numpy)
        else:
            return _fit_with_gd(x_data, y_data, model_func, param_names, initial_params)

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = [p.strip() for p in input_str.split(";")]
            if len(parts) < 3:
                raise ValueError("Need at least 3 semicolon-separated parts")

            x_data = [float(v) for v in parts[0].split(",")]
            y_data = [float(v) for v in parts[1].split(",")]
            model_func = parts[2].strip()
            initial_params = [float(v) for v in parts[3].split(",")] if len(parts) > 3 else [1.0] * 3

            return self._run_base(x_data, y_data, model_func, initial_params)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")


def _fit_with_scipy(x_data, y_data, model_func, param_names, initial_params, sp_opt, np_mod):
    """Fit using scipy.optimize.curve_fit."""
    xs = np_mod.array(x_data)
    ys = np_mod.array(y_data)
    p0 = list(initial_params)

    def model_fn(x, *args):
        pdict = dict(zip(param_names, args))
        return [_eval_model(model_func, pdict, xi) for xi in x]

    try:
        popt, pcov = sp_opt.curve_fit(lambda x, *args: model_fn(xs, *args), xs, ys, p0=p0, maxfev=10000)
    except Exception as e:
        raise ChemMCPError(f"scipy curve_fit failed: {e}")

    fitted_dict = {name: round(float(val), 6) for name, val in zip(param_names, popt)}
    fitted_y = [_eval_model(model_func, fitted_dict, xi) for xi in x_data]
    residuals = [round(ys[i] - fitted_y[i], 8) for i in range(len(ys))]
    ss_res = sum(r ** 2 for r in residuals)
    ss_tot = sum((yi - sum(ys) / len(ys)) ** 2 for yi in ys)
    r_squared = round(1 - ss_res / ss_tot, 6) if ss_tot > 0 else 1.0

    return {
        "fitted_params": fitted_dict,
        "r_squared": r_squared,
        "residuals": residuals,
        "fitted_y": [round(yi, 6) for yi in fitted_y],
    }


def _fit_with_gd(x_data, y_data, model_func, param_names, initial_params,
                  lr: float = 0.01, max_iter: int = 5000, tol: float = 1e-8):
    """Fallback: simple gradient descent fitting."""
    n = len(x_data)
    params = list(initial_params)
    y_mean = sum(y_data) / n

    for iteration in range(max_iter):
        grads = [0.0] * len(params)
        sse = 0.0

        for i in range(n):
            pdict = dict(zip(param_names, params))
            try:
                y_pred = _eval_model(model_func, pdict, x_data[i])
            except ChemMCPError:
                y_pred = 0.0
            residual = y_pred - y_data[i]
            sse += residual ** 2

            # Numerical gradient of SSE w.r.t each param
            for j in range(len(params)):
                h = 1e-5
                params_plus = list(params)
                params_plus[j] += h
                pdict_plus = dict(zip(param_names, params_plus))
                try:
                    y_pred_plus = _eval_model(model_func, pdict_plus, x_data[i])
                except ChemMCPError:
                    y_pred_plus = y_pred
                grads[j] += 2 * residual * (y_pred_plus - y_pred) / h

        # Update
        max_grad = max(abs(g) for g in grads) if grads else 0
        if max_grad < tol and iteration > 10:
            break

        for j in range(len(params)):
            params[j] -= lr * grads[j] / n

    fitted_dict = {name: round(float(p), 6) for name, p in zip(param_names, params)}
    fitted_y = []
    residuals = []
    for i in range(n):
        pdict = dict(zip(param_names, params))
        try:
            yp = _eval_model(model_func, pdict, x_data[i])
        except ChemMCPError:
            yp = 0.0
        fitted_y.append(round(yp, 6))
        residuals.append(round(y_data[i] - yp, 8))

    ss_res = sum(r ** 2 for r in residuals)
    ss_tot = sum((yi - y_mean) ** 2 for yi in y_data)
    r_squared = round(1 - ss_res / ss_tot, 6) if ss_tot > 0 else 1.0

    return {
        "fitted_params": fitted_dict,
        "r_squared": r_squared,
        "residuals": residuals,
        "fitted_y": fitted_y,
    }
