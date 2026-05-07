import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RegressionDiagnostics(BaseTool):
    """
    回归诊断工具。
    执行残差分析、杠杆值计算、Cook距离、影响点识别、Durbin-Watson检验等回归诊断。
    """
    __version__ = "0.1.0"
    name = "RegressionDiagnostics"
    func_name = "diagnose_regression"
    description = "Regression diagnostics: residual analysis, leverage values, Cook's distance, influential points, standardized residuals, DFFITS, Durbin-Watson test."
    implementation_description = "Performs comprehensive OLS regression diagnostics including: residuals (raw/standardized/studentized), leverage (hat matrix diagonal), Cook's D, DFFITS, Durbin-Watson statistic for autocorrelation, normality tests of residuals."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Regression", "Diagnostics", "Residuals", "Statistics", "QA/QC", "Leverage"]
    required_envs = []

    code_input_sig = [
        ("x_data", "list", "N/A", "Independent variable data (list of numbers)."),
        ("y_data", "list", "N/A", "Dependent variable data (list of numbers)."),
        ("residuals", "list", "", "Pre-computed residuals (optional; if provided with fitted_values, skips regression)."),
        ("fitted_values", "list", "", "Pre-computed fitted/predicted values (optional)."),
        ("influence_threshold_cook", "float", "4/n", "Cook's distance threshold for influential points (default 4/n or 1)."),
        ("leverage_threshold", "float", "2p/n", "Leverage threshold (default 2(p+1)/n for simple linear)."),
        ("confidence_level", "float", "0.95", "Confidence level for intervals."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("regression_summary", "dict", "Basic regression info (slope, intercept, R², etc.)."),
        ("residuals", "dict", "Raw, standardized, and studentized residuals."),
        ("leverage", "dict", "Leverage values per point and threshold analysis."),
        ("cooks_distance", "dict", "Cook's D values and influential point identification."),
        ("dffits", "list", "DFFITS values for each observation."),
        ("durbin_watson", "dict", "Durbin-Watson statistic and autocorrelation assessment."),
        ("normality_test", "dict", "Normality assessment of residuals (skewness, kurtosis, Shapiro-like proxy)."),
        ("influential_points", "list", "List of identified influential/high-leverage points."),
        ("diagnostics_summary", "str", "Overall diagnostic summary and recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "x_data": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                "y_data": [1.2, 2.1, 2.9, 4.2, 5.1, 5.8, 6.9, 8.1, 9.0, 10.1],
            },
            "text_input": {"params_str": "see code input"},
            "output": {"regression_summary": {"r_squared": 0.998}},
        },
        {
            "code_input": {
                "x_data": [1, 2, 3, 4, 5, 6, 7, 8],
                "y_data": [2.1, 3.9, 6.2, 8.0, 10.1, 11.8, 14.5, 15.9],
            },
            "text_input": {"params_str": "outlier data"},
            "output": {},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _mean(data: List[float]) -> float:
        return sum(data) / len(data)

    @staticmethod
    def _std(data: List[float], ddof: int = 1) -> float:
        n = len(data)
        m = sum(data) / n
        return math.sqrt(sum((x - m)**2 for x in data) / max(n - ddof, 1))

    def _simple_linear_regression(self, x: List[float], y: List[float]) -> dict:
        """Simple linear regression y = a + bx."""
        n = len(x)
        mx = self._mean(x)
        my = self._mean(y)
        ss_xx = sum((xi - mx)**2 for xi in x)
        ss_xy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
        ss_yy = sum((yi - my)**2 for yi in y)

        slope = ss_xy / ss_xx if ss_xx > 0 else 0
        intercept = my - slope * mx

        fitted = [intercept + slope * xi for xi in x]
        residuals = [yi - yh for yi, yh in zip(y, fitted)]

        sst = ss_yy
        sse = sum(r**2 for r in residuals)
        ssr = sst - sse
        r_sq = 1 - sse / sst if sst > 0 else 1.0

        p = 2  # parameters (slope + intercept)
        mse = sse / (n - p) if n > p else 0
        std_err = math.sqrt(mse) if mse >= 0 else 0

        # Standard error of slope and intercept
        se_slope = math.sqrt(mse / ss_xx) if ss_xx > 0 else 0
        se_intercept = math.sqrt(mse * (1/n + mx**2/ss_xx)) if ss_xx > 0 else 0

        # t-statistics for coefficients
        t_slope = slope / se_slope if se_slope > 0 else float('inf')
        t_int = intercept / se_intercept if se_intercept > 0 else float('inf')

        return {
            "slope": slope, "intercept": intercept,
            "fitted_values": fitted, "residuals": residuals,
            "r_squared": r_sq, "r_squared_adj": 1 - (1-r_sq)*(n-1)/(n-p) if n > p else 0,
            "sse": sse, "ssr": ssr, "sst": sst,
            "mse": mse, "rmse": std_err,
            "se_slope": se_slope, "se_intercept": se_intercept,
            "t_slope": t_slope, "t_intercept": t_int,
            "n": n, "p": p,
            "mx": mx, "my": my,
            "ss_xx": ss_xx,
        }

    def _run_base(
        self,
        x_data: List[float],
        y_data: List[float],
        residuals: Optional[List[float]] = None,
        fitted_values: Optional[List[float]] = None,
        influence_threshold_cook: Optional[float] = None,
        leverage_threshold: Optional[float] = None,
        confidence_level: float = 0.95,
    ) -> dict:
        """Core logic: perform full regression diagnostics."""
        if not x_data or not y_data:
            raise ChemMCPError("x_data and y_data are required.")
        if len(x_data) != len(y_data):
            raise ChemMCPError(f"Length mismatch: {len(x_data)} x vs {len(y_data)} y.")

        n = len(x_data)
        if n < 3:
            raise ChemMCPError(f"Need at least 3 points for meaningful diagnostics. Got {n}.")

        # --- Step 1: Regression (or use provided values) ---
        if residuals is not None and fitted_values is not None:
            reg = {}
            reg["residuals"] = list(residuals)
            reg["fitted_values"] = list(fitted_values)
            reg["n"] = n
            # Recompute basic stats from provided data
            reg["my"] = self._mean(y_data)
            reg["mse"] = sum(r**2 for r in residuals) / max(n - 2, 1)
            reg["rmse"] = math.sqrt(max(reg["mse"], 0))
            reg["mx"] = self._mean(x_data)
            reg["ss_xx"] = sum((xi - reg["mx"])**2 for xi in x_data)
        else:
            reg = self._simple_linear_regression(x_data, y_data)

        e = reg["residuals"]
        y_hat = reg["fitted_values"]
        mse = reg.get("mse", sum(r**2 for r in e) / max(n-2, 1))
        rmse = reg.get("rmse", math.sqrt(max(mse, 0)))
        mx = reg.get("mx", self._mean(x_data))
        ss_xx = reg.get("ss_xx", sum((xi - mx)**2 for xi in x_data))
        p_params = 2  # simple linear regression: intercept + slope

        # --- Step 2: Leverage (hat matrix diagonal) ---
        # For simple linear regression: h_ii = 1/n + (xi - x̄)²/Σ(xj-x̄)²
        h_vals = []
        for i, xi in enumerate(x_data):
            hi = 1.0/n + ((xi - mx)**2 / ss_xx) if ss_xx > 0 else 1.0/n
            h_vals.append(hi)

        avg_h = p_params / n  # average leverage
        lev_thresh = leverage_threshold if leverage_threshold is not None else (2.0 * p_params / n)
        high_leverage = [{"index": i, "x": x_data[i], "y": y_data[i], "leverage": round(h_vals[i], 6),
                          "exceeds_threshold": h_vals[i] > lev_thresh}
                         for i in range(n) if h_vals[i] > lev_thresh]

        # --- Step 3: Residuals analysis ---
        # Standardized residuals: ri = ei / (RMSE × √(1-hii))
        std_resids = []
        for i in range(n):
            denom = rmse * math.sqrt(max(1 - h_vals[i], 1e-12))
            std_resids.append(e[i] / denom if denom > 0 else 0)

        # Studentized residuals (approximate using deleted residuals formula)
        # t_i = ri × √((n-p-1)/(n-p-ri²)) — approximate
        stud_resids = []
        dof_del = max(n - p_params - 1, 1)
        for i in range(n):
            ri = std_resids[i]
            factor = math.sqrt(dof_del / max(dof_del + 0, 1e-12))
            # More accurate approximation
            si_squared = mse * (1 - h_vals[i]) / max(1 - h_vals[i], 1e-12)
            si = math.sqrt(max(si_squared, 0))
            stud_resids.append(e[i] / si if si > 0 else 0)

        # --- Step 4: Cook's Distance ---
        cook_thresh = influence_threshold_cook if influence_threshold_cook is not None else (4.0 / n)
        cooks_d = []
        for i in range(n):
            di = (std_resids[i]**2 / p_params) * (h_vals[i] / max(1 - h_vals[i], 1e-12))
            cooks_d.append(di)

        influential = [
            {"index": i, "x": x_data[i], "y": y_data[i],
             "cooks_d": round(cooks_d[i], 6), "leverage": round(h_vals[i], 6),
             "standardized_residual": round(std_resids[i], 4)}
            for i in range(n) if cooks_d[i] > cook_thresh
        ]

        # --- Step 5: DFFITS ---
        dffits_vals = []
        for i in range(n):
            dfi = std_resids[i] * math.sqrt(h_vals[i] / max(1 - h_vals[i], 1e-12))
            dffits_vals.append(round(dfi, 6))

        # Threshold: 2√((p+1)/n)
        dffits_thresh = 2 * math.sqrt(p_params / n)
        high_dffits = [{"index": i, "dffits": dffits_vals[i]} for i in range(n) if abs(dffits_vals[i]) > dffits_thresh]

        # --- Step 6: Durbin-Watson Statistic ---
        dw_num = sum((e[i] - e[i-1])**2 for i in range(1, n))
        dw_den = sum(ei**2 for ei in e)
        dw_stat = dw_num / dw_den if dw_den > 0 else 2.0
        # DW ≈ 2: no autocorrelation; DW → 0: positive AC; DW → 4: negative AC
        if dw_stat < 1.5:
            dw_assessment = f"Positive autocorrelation suspected (DW={dw_stat:.3f} < 1.5)"
        elif dw_stat > 2.5:
            dw_assessment = f"Negative autocorrelation suspected (DW={dw_stat:.3f} > 2.5)"
        else:
            dw_assessment = f"No significant autocorrelation (DW={dw_stat:.3f} ≈ 2.0)"

        # --- Step 7: Normality of residuals ---
        res_mean = self._mean(e)
        res_std = self._std(e, ddof=1)
        skewness = sum(((ri - res_mean) / res_std)**3 for ri in e) / n if res_std > 0 else 0
        kurtosis = sum(((ri - res_mean) / res_std)**4 for ri in e) / n - 3 if res_std > 0 else 0
        # Simple normality check: |skewness| < 2 and |kurtosis| < 3 roughly normal
        is_normal_approx = abs(skewness) < 2 and abs(kurtosis) < 3
        norm_assessment = (
            f"Skewness={skewness:.4f}, Kurtosis={kurtosis:.4f}. "
            f"{'Approximately normal' if is_normal_approx else 'Potential non-normality detected'} "
            f"(check Q-Q plot for confirmation)."
        )

        # --- Step 8: Influential points summary ---
        all_influential_indices = set()
        for pt in high_leverage:
            all_influential_indices.add(pt["index"])
        for pt in influential:
            all_influential_indices.add(pt["index"])
        for pt in high_dffits:
            all_influential_indices.add(pt["index"])

        influencels = sorted(all_influential_indices)

        # --- Summary ---
        issues = []
        if high_leverage:
            issues.append(f"⚠️ {len(high_leverage)} high-leverage point(s) (threshold={lev_thresh:.4f}).")
        if influential:
            issues.append(f"⚠️ {len(influential)} influential point(s) by Cook's D (threshold={cook_thresh:.4f}).")
        if high_dffits:
            issues.append(f"⚠️ {len(high_dffits)} high-DFFITS point(s) (threshold={dffits_thresh:.4f}).")
        if dw_stat < 1.5 or dw_stat > 2.5:
            issues.append(f"⚠️ Autocorrelation in residuals ({dw_assessment}).")
        if not is_normal_approx:
            issues.append(f"⚠️ Non-normal residuals ({norm_assessment}).")

        if issues:
            diag_summary = "Regression diagnostics:\n" + "\n".join(issues) + (
                f"\n\nTotal influential points: {len(influencels)} at indices: {influencels}."
                if influencels else ""
            )
        else:
            diag_summary = (
                f"✅ No major issues detected.\n"
                f"   R²={reg.get('r_squared', 0):.4f}, RMSE={rmse:.6g}, DW={dw_stat:.3f}, "
                f"no high-leverage/influential points, residuals approximately normal."
            )

        logger.info(f"Regression diagnostics: R²={reg.get('r_squared',0):.4f}, DW={dw_stat:.3f}, influential={len(influential)}")
        return {
            "regression_summary": {
                "slope": round(reg.get("slope", 0), 8),
                "intercept": round(reg.get("intercept", 0), 8),
                "r_squared": round(reg.get("r_squared", 0), 6),
                "r_squared_adj": round(reg.get("r_squared_adj", 0), 6),
                "rmse": round(rmse, 8),
                "n": n,
                "equation": f"y = {reg.get('intercept', 0):.6g} + {reg.get('slope', 0):.6g}·x",
            },
            "residuals": {
                "raw": [round(ei, 8) for ei in e],
                "standardized": [round(sr, 6) for sr in std_resids],
                "studentized": [round(st, 6) for st in stud_resids],
                "mean": round(res_mean, 8),
                "std": round(res_std, 8),
            },
            "leverage": {
                "values": [round(h, 6) for h in h_vals],
                "average": round(avg_h, 6),
                "threshold_used": round(lev_thresh, 6),
                "high_leverage_points": high_leverage,
            },
            "cooks_distance": {
                "values": [round(cd, 6) for cd in cooks_d],
                "threshold_used": round(cook_thresh, 6),
                "influential_points": influential,
            },
            "dffits": dffits_vals,
            "durbin_watson": {
                "statistic": round(dw_stat, 6),
                "assessment": dw_assessment,
                "range_interpretation": "0=positive AC, 2=no AC, 4=negative AC",
            },
            "normality_test": {
                "skewness": round(skewness, 6),
                "kurtosis": round(kurtosis, 6),
                "approximately_normal": is_normal_approx,
                "assessment": norm_assessment,
            },
            "influential_points": influencels,
            "diagnostics_summary": diag_summary,
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
