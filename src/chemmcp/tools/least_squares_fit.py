import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _beta_reg_inc(a: float, b: float, x: float, max_iter: int = 200) -> float:
    """
    Regularized incomplete beta function I_x(a, b) using continued fraction expansion.
    Used for computing t-distribution CDF/p-values.
    """
    if x < 0 or x > 1:
        return 0.0 if x < 0.5 else 1.0
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0

    # Use symmetry for efficiency
    if x > (a + 1) / (a + b + 2):
        return 1.0 - _beta_reg_inc(b, a, 1.0 - x, max_iter)

    # Continued fraction for I_x(a,b)
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)

    # Front factor
    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lbeta) / a

    # Modified Lentz's method for continued fraction
    f = 1.0
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1)
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d

    for m in range(1, max_iter + 1):
        # Even step
        m2 = 2 * m
        aa = m * (b - m) * x / ((a + m2 - 1) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        f *= delta
        if abs(delta - 1.0) < 1e-12:
            break

        # Odd step
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        f *= delta
        if abs(delta - 1.0) < 1e-12:
            break

    return front * f


def _t_cdf(t_val: float, df: int) -> float:
    """CDF of Student's t-distribution via regularized incomplete beta."""
    x = df / (df + t_val ** 2)
    ib = _beta_reg_inc(df / 2.0, 0.5, x)
    if t_val >= 0:
        return 1.0 - 0.5 * ib
    else:
        return 0.5 * ib


def _t_ppf(p: float, df: int) -> float:
    """
    Approximate inverse CDF (percent point function) of t-distribution.
    Uses Newton's method with normal approximation as starting guess.
    """
    if p <= 0 or p >= 1:
        raise ChemMCPError("p must be strictly between 0 and 1.")
    # Starting guess from normal approximation
    import random as _r
    try:
        # Approximate inverse normal
        if abs(p - 0.5) < 1e-10:
            t = 0.0
        elif p < 0.5:
            t = -_approx_inverse_normal(1.0 - p)
        else:
            t = _approx_inverse_normal(p)
    except Exception:
        t = 0.0

    # Refine with Newton-Raphson on CDF
    for _ in range(50):
        cdf_val = _t_cdf(t, df)
        # PDF of t-distribution
        pdf_val = math.exp(
            math.lgamma((df + 1) / 2.0) - math.lgamma(df / 2.0)
            - 0.5 * math.log(df * math.pi)
            - (df + 1) / 2.0 * math.log(1.0 + t ** 2 / df)
        )
        if pdf_val < 1e-30:
            break
        step = (cdf_val - p) / pdf_val
        t -= step
        if abs(step) < 1e-12:
            break
    return t


def _approx_inverse_normal(p: float) -> float:
    """Approximate inverse standard normal CDF (Abramowitz & Stegun)."""
    if p < 1e-10:
        return -6.0
    if p > 1 - 1e-10:
        return 6.0
    if abs(p - 0.5) < 1e-10:
        return 0.0

    q = p - 0.5
    if abs(q) <= 0.42:
        r = q * q
        t = q * (((a4 * r + a3) * r + a2) * r + a1) / ((((b4 * r + b3) * r + b2) * r + b1) * r + 1.0)
    else:
        r = p if q > 0 else 1.0 - p
        r = math.sqrt(-math.log(r))
        t = c0 + c1 * r + c2 * r * r + c3 * r ** 3 + c4 * r ** 4 + c5 * r ** 5
        if q < 0:
            t = -t
    return t


# Rational approximation constants for inverse normal
a1 = -3.969683028665376e+01
a2 = 2.209460984245205e+02
a3 = -2.759285104469687e+02
a4 = 1.383577518672690e+02
a5 = -3.066479806614716e+01
a6 = 2.506628277459239e+00
b1 = -5.447609879822406e+01
b2 = 1.615858368580409e+02
b3 = -1.556989798598866e+02
b4 = 6.680131188771972e+01
b5 = -1.328068155288572e+01
c0 = -7.784894002430293e-03
c1 = -3.223964580411365e-01
c2 = -2.400758277161838e+00
c3 = -2.549732539343734e+00
c4 = 4.374664141464968e+00
c5 = 2.938163982698783e+00


@ChemMCPManager.register_tool
class LeastSquaresFit(BaseTool):
    """
    最小二乘拟合工具。
    线性回归、标准曲线拟合，输出完整统计量：方程、R²、Pearson r、p值、置信区间等。
    """
    __version__ = "0.1.0"
    name = "LeastSquaresFit"
    func_name = "least_squares_fit"
    description = "Least squares linear regression and standard curve fitting with full statistics: equation, R², Pearson r, p-value, confidence intervals, residuals."
    implementation_description = "Implements ordinary least squares (OLS) linear regression with comprehensive statistical diagnostics. Includes slope/intercept with confidence intervals at user-specified alpha level, t-test for slope significance, ANOVA table components (SSR/SSE/SST), RMSE, MAE, standard error of estimate s_yx, predicted values and residuals. Supports forced-through-origin regression."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Regression", "Linear Fit", "Standard Curve", "Chemometrics", "Statistics", "p-value", "Confidence Interval"]
    required_envs = []

    code_input_sig = [
        ("x", "List[float]", "N/A", "Independent variable values (x data points)."),
        ("y", "List[float]", "N/A", "Dependent variable values (y data points)."),
        ("alpha", "float", "0.95", "Confidence level for CI (e.g., 0.95 for 95% CI)."),
        ("force_through_origin", "bool", "False", "Force regression line through origin (intercept=0)."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with x, y, alpha, force_through_origin."),
    ]

    output_sig = [
        ("equation", "str", "Fitted equation string (y = ax + b)."),
        ("slope", "float", "Regression slope (a in y=ax+b)."),
        ("intercept", "float", "Regression y-intercept (b in y=ax+b)."),
        ("r_squared", "float", "Coefficient of determination R²."),
        ("pearson_r", "float", "Pearson correlation coefficient r."),
        ("pearson_r_pvalue", "float", "Two-tailed p-value for Pearson correlation significance."),
        ("slope_pvalue", "float", "Two-tailed p-value testing H₀: slope=0."),
        ("sse", "float", "Sum of squared errors (residuals)."),
        ("rmse", "float", "Root mean square error."),
        ("mae", "float", "Mean absolute error."),
        ("s_yx", "float", "Standard error of the estimate (residual standard deviation)."),
        ("slope_ci", "list", "Confidence interval for slope [lower, upper] at given alpha."),
        ("intercept_ci", "list", "Confidence interval for intercept [lower, upper] at given alpha."),
        ("y_pred", "list", "Predicted/fitted y values."),
        ("residuals", "list", "Residuals (y_observed - y_predicted)."),
        ("n", "int", "Number of data points."),
        ("df_residual", "int", "Degrees of freedom for residuals."),
        ("ssr", "float", "Sum of squares due to regression."),
        ("sst", "float", "Total sum of squares."),
        ("adjusted_r_squared", "float", "Adjusted R²."),
    ]

    examples = [
        {
            "code_input": {
                "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                "y": [0.1, 1.1, 2.0, 2.9, 4.1, 5.0],
                "alpha": 0.95,
            },
            "text_input": {"params_str": '{"x":[0,1,2,3,4,5],"y":[0.1,1.1,2,2.9,4.1,5],"alpha":0.95}'},
            "output": {
                "equation": "y = 0.9714x + 0.0967",
                "r_squared": 0.9986,
                "pearson_r": 0.9993,
                "slope_pvalue": 0.0001,
            },
        },
        {
            "code_input": {
                "x": [0.0, 1.0, 2.0, 5.0, 10.0],
                "y": [0.02, 0.105, 0.205, 0.502, 1.01],
                "alpha": 0.95,
            },
            "text_input": {"params_str": '{"x":[0,1,2,5,10],"y":[0.02,0.105,0.205,0.502,1.01],"alpha":0.95}'},
            "output": {
                "r_squared": 0.9998,
                "rmse": 0.0103,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _mean(data: List[float]) -> float:
        return sum(data) / len(data) if data else 0.0

    @staticmethod
    def _stats_1d(data: List[float]) -> tuple:
        """Return (mean, variance_sample, std_sample)."""
        n = len(data)
        m = sum(data) / n
        var = sum((x - m) ** 2 for x in data) / (n - 1) if n > 1 else 0.0
        return m, var, math.sqrt(var)

    def _run_base(
        self,
        x: List[float],
        y: List[float],
        alpha: float = 0.95,
        force_through_origin: bool = False,
    ) -> dict:
        """Core logic: OLS linear regression with full diagnostics."""
        if not x or not y:
            raise ChemMCPError("x and y cannot be empty.")
        n = len(x)
        if n != len(y):
            raise ChemMCPError(f"Length mismatch: {n} x values vs {len(y)} y values.")
        if n < 2:
            raise ChemMCPError(f"Need at least 2 data points, got {n}.")
        if not (0 < alpha < 1):
            raise ChemMCPError(f"alpha must be between 0 and 1 (exclusive), got {alpha}.")

        # ---- OLS computation ----
        if force_through_origin:
            sxx = sum(xi ** 2 for xi in x)
            sxy = sum(xi * yi for xi, yi in zip(x, y))
            if sxx < 1e-30:
                raise ChemMCPError("All x values are zero; cannot fit through origin.")
            slope = sxy / sxx
            intercept = 0.0
            p = 1  # one parameter (slope only)
        else:
            mx = self._mean(x)
            my = self._mean(y)
            sxx = sum((xi - mx) ** 2 for xi in x)
            sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
            if sxx < 1e-30:
                raise ChemMCPError("Variance in x is zero; cannot fit.")
            slope = sxy / sxx
            intercept = my - slope * mx
            p = 2  # two parameters

        # Predicted values & residuals
        y_pred = [slope * xi + intercept for xi in x]
        residuals = [yi - ypi for yi, ypi in zip(y, y_pred)]

        # Sum of squares
        my_mean = self._mean(y)
        sst = sum((yi - my_mean) ** 2 for yi in y)
        sse = sum(r ** 2 for r in residuals)
        ssr = sst - sse

        # R²
        r_sq = 1.0 - sse / sst if sst > 1e-30 else 1.0

        # Adjusted R²
        adj_r2 = 1.0 - (1.0 - r_sq) * (n - 1) / (n - p) if n > p else 0.0

        # Standard error of estimate
        df_res = n - p
        s_yx = math.sqrt(sse / df_res) if df_res > 0 else 0.0

        # RMSE, MAE
        rmse = math.sqrt(sse / n) if n > 0 else 0.0
        mae = sum(abs(r) for r in residuals) / n

        # Pearson r
        sx = math.sqrt(sxx / (n - 1)) if n > 1 else 0.0
        sy_vec = [(yi - my_mean) for yi in y]
        syy = sum(v ** 2 for v in sy_vec)
        sy = math.sqrt(syy / (n - 1)) if n > 1 else 0.0
        pearson_r = sxy / ((n - 1) * sx * sy) if sx * sy > 1e-30 else (1.0 if sxy >= 0 else -1.0)

        # t-statistic for slope (if not forced through origin)
        slope_se = s_yx / math.sqrt(sxx) if (not force_through_origin and sxx > 1e-30) else (s_yx / math.sqrt(sxx) if sxx > 1e-30 else float('inf'))
        t_slope = slope / slope_se if slope_se > 1e-30 else (float('inf') if abs(slope) > 1e-30 else 0.0)

        # Two-tailed p-value for slope
        if df_res > 0 and abs(t_slope) < 1e8:
            slope_pval = 2.0 * (1.0 - _t_cdf(abs(t_slope), df_res))
        elif abs(t_slope) >= 1e8:
            slope_pval = 0.0
        else:
            slope_pval = 1.0

        # Two-tailed p-value for Pearson r (same t-statistic essentially)
        if df_res > 0 and abs(pearson_r) < 1.0 - 1e-10:
            t_r = pearson_r * math.sqrt((n - 2) / (1.0 - min(pearson_r ** 2, 0.999999)))
            r_pval = 2.0 * (1.0 - _t_cdf(abs(t_r), df_res))
        else:
            r_pval = 0.0 if abs(pearson_r) > 0.9999 else 1.0

        # Confidence intervals for slope and intercept
        tail_prob = (1.0 - alpha) / 2.0
        try:
            t_crit = _t_ppf(1.0 - tail_prob, df_res)
        except Exception:
            t_crit = 1.96  # fallback to normal approx

        if not force_through_origin:
            slope_ci_lower = slope - t_crit * slope_se
            slope_ci_upper = slope + t_crit * slope_se
            # SE of intercept
            x_bar = self._mean(x)
            intercept_se = s_yx * math.sqrt(1.0 / n + x_bar ** 2 / sxx)
            int_ci_lower = intercept - t_crit * intercept_se
            int_ci_upper = intercept + t_crit * intercept_se
            slope_ci = [round(slope_ci_lower, 10), round(slope_ci_upper, 10)]
            intercept_ci = [round(int_ci_lower, 10), round(int_ci_upper, 10)]
        else:
            slope_se_orig = s_yx / math.sqrt(sxx) if sxx > 1e-30 else 0.0
            slope_ci = [round(slope - t_crit * slope_se_orig, 10), round(slope + t_crit * slope_se_orig, 10)]
            intercept_ci = [0.0, 0.0]

        logger.info(
            f"LeastSquaresFit: y={slope:.6g}x+{intercept:.6g}, "
            f"R²={r_sq:.6f}, r={pearson_r:.6f}, p(slope)={slope_pval:.2e}"
        )

        return {
            "equation": f"y = {slope:.6g}x + {intercept:.6g}",
            "slope": round(slope, 12),
            "intercept": round(intercept, 12),
            "r_squared": round(r_sq, 6),
            "pearson_r": round(pearson_r, 6),
            "pearson_r_pvalue": round(min(r_pval, 1.0), 8),
            "slope_pvalue": round(min(slope_pval, 1.0), 8),
            "sse": round(sse, 12),
            "rmse": round(rmse, 8),
            "mae": round(mae, 8),
            "s_yx": round(s_yx, 8),
            "slope_ci": slope_ci,
            "intercept_ci": intercept_ci,
            "y_pred": [round(v, 10) for v in y_pred],
            "residuals": [round(r, 10) for r in residuals],
            "n": n,
            "df_residual": df_res,
            "ssr": round(ssr, 12),
            "sst": round(sst, 12),
            "adjusted_r_squared": round(adj_r2, 6),
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input. Expected JSON string with keys: x, y, alpha, ...")
        return self._run_base(**kwargs)
