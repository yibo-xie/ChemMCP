import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ---- Approximate t-distribution critical values ----
# Keys: (df_rounded, alpha_two_tailed) -> t_critical
# For one-tailed: use alpha/2
def _t_critical_approx(df: float, alpha: float, two_tail: bool = True) -> float:
    """Approximate t-critical value using inverse CDF approximation."""
    if df <= 0:
        return float('inf')
    # Use approximation via normal distribution for large df
    if df > 100:
        import math as _m
        a = 1 - (alpha / 2) if two_tail else 1 - alpha
        z = _sqrt(2) * _erfinv(2 * a - 1)
        return z + (z**3 + z) / (4 * df) + (5*z**5 + 16*z**3 + 3*z) / (96 * df**2)

    # Approximate using Wilson-Hilferty for small df
    # Use precomputed table for common values and interpolate
    _t_table = {
        # df: {alpha_2tail}
        1: {0.10:6.314, 0.05:12.706, 0.02:31.821, 0.01:63.657},
        2: {0.10:2.920, 0.05:4.303, 0.02:6.965, 0.01:9.925},
        3: {0.10:2.353, 0.05:3.182, 0.02:4.541, 0.01:5.841},
        4: {0.10:2.132, 0.05:2.776, 0.02:3.747, 0.01:4.604},
        5: {0.10:2.015, 0.05:2.571, 0.02:3.365, 0.01:4.032},
        6: {0.10:1.943, 0.05:2.447, 0.02:3.143, 0.01:3.707},
        7: {0.10:1.895, 0.05:2.365, 0.02:2.998, 0.01:3.499},
        8: {0.10:1.860, 0.05:2.306, 0.02:2.896, 0.01:3.355},
        9: {0.10:1.833, 0.05:2.262, 0.02:2.821, 0.01:3.250},
        10:{0.10:1.812, 0.05:2.228, 0.02:2.764, 0.01:3.169},
        11:{0.10:1.796, 0.05:2.201, 0.02:2.718, 0.01:3.106},
        12:{0.10:1.782, 0.05:2.179, 0.02:2.681, 0.01:3.055},
        13:{0.10:1.771, 0.05:2.160, 0.02:2.650, 0.01:3.012},
        14:{0.10:1.761, 0.05:2.145, 0.02:2.624, 0.01:2.977},
        15:{0.10:1.753, 0.05:2.131, 0.02:2.602, 0.01:2.947},
        16:{0.10:1.746, 0.05:2.120, 0.02:2.583, 0.01:2.921},
        17:{0.10:1.740, 0.05:2.110, 0.02:2.567, 0.01:2.898},
        18:{0.10:1.734, 0.05:2.101, 0.02:2.552, 0.01:2.878},
        19:{0.10:1.728, 0.05:2.093, 0.02:2.539, 0.01:2.861},
        20:{0.10:1.725, 0.05:2.086, 0.02:2.528, 0.01:2.845},
        25:{0.10:1.708, 0.05:2.060, 0.02:2.485, 0.01:2.787},
        30:{0.10:1.697, 0.05:2.042, 0.02:2.457, 0.01:2.750},
        40:{0.10:1.684, 0.05:2.021, 0.02:2.423, 0.01:2.704},
        50:{0.10:1.676, 0.05:2.009, 0.02:2.403, 0.01:2.678},
        60:{0.10:1.671, 0.05:2.000, 0.02:2.390, 0.01:2.660},
        80:{0.10:1.664, 0.05:1.990, 0.02:2.374, 0.01:2.639},
        100:{0.10:1.660, 0.05:1.984, 0.02:2.364, 0.01:2.626},
    }

    a_key = alpha if two_tail else alpha * 2  # store as two-tailed in table
    int_df = int(round(df))

    if int_df in _t_table:
        alphas = sorted(_t_table[int_df].keys())
        if a_key in _t_table[int_df]:
            t_val = _t_table[int_df][a_key]
        elif a_key < alphas[0]:
            t_val = _t_table[int_df][alphas[0]]
        elif a_key > alphas[-1]:
            t_val = _t_table[int_df][alphas[-1]]
        else:
            # Interpolate
            lo = max(a for a in alphas if a <= a_key)
            hi = min(a for a in alphas if a >= a_key)
            frac = (a_key - lo) / (hi - lo) if hi != lo else 0
            t_val = _t_table[int_df][lo] + frac * (_t_table[int_df][hi] - _t_table[int_df][lo])
    else:
        # Extrapolate from nearest
        keys = sorted(k for k in _t_table if k <= int_df)
        if keys:
            base = keys[-1]
            if a_key in _t_table[base]:
                base_t = _t_table[base][a_key]
            else:
                base_t = list(_t_table[base].values())[len(_t_table[base]) // 2]
            # Rough adjustment for higher df
            ratio = math.sqrt(base / df) if df > base else 1.0
            t_val = base_t * max(ratio, 0.95)
        else:
            t_val = 2.0  # fallback

    return t_val


def _p_value_from_t(t_stat: float, df: float) -> float:
    """Approximate two-tailed p-value from t-statistic.
    Uses Abramowitz & Stegun normal CDF + Fisher correction for t-distribution."""
    import math as _m
    if df <= 0:
        return 0.0
    abs_t = abs(t_stat)

    # Fisher's approximation: convert t to z
    if df >= 3:
        z = abs_t * (1 - 1/(4*df)) / _m.sqrt(1 + abs_t**2 / (2*df))
    else:
        z = abs_t  # fallback for very small df

    # Standard normal CDF (Abramowitz & Stegun 26.2.17)
    p = 2.0 * (1.0 - _norm_cdf(z))
    return max(0.0, min(1.0, p))


def _norm_cdf(z):
    """Standard normal CDF using math.erf (exact)."""
    import math as _m
    return 0.5 * (1.0 + _m.erf(z / _m.sqrt(2.0)))


def _beta_reg(x: float, a: float, b: float, max_iter: int = 200) -> float:
    """Approximate regularized incomplete beta function I_x(a,b) using continued fraction."""
    import math as _m
    if x < 0 or x > 1:
        raise ValueError("x must be in [0,1]")
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0

    # Use continued fraction representation for I_x(a,b)
    # Based on Lentz's method
    bt = _m.exp(
        _m.lgamma(a + b) - _m.lgamma(a) - _m.lgamma(b)
        + a * _m.log(x) + b * _m.log(max(1 - x, 1e-30))
    )

    if x < (a + 1) / (a + b + 2):
        return bt * _beta_cf(x, a, b, max_iter) / a
    else:
        return 1.0 - bt * _beta_cf(1 - x, b, a, max_iter) / b


def _beta_cf(x: float, a: float, b: float, max_iter: int = 200) -> float:
    """Continued fraction for incomplete beta function."""
    import math as _m
    EPS = 1e-12
    f = 1.0
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1)
    if abs(d) < EPS:
        d = EPS
    d = 1.0 / d
    f = d

    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa_m = m * (b - m) * x / ((a + m2 - 1) * (a + m2))
        d = 1.0 + aa_m * d
        if abs(d) < EPS:
            d = EPS
        c = 1.0 + aa_m / c
        if abs(c) < EPS:
            c = EPS
        d = 1.0 / d
        f *= d * c

        aa_m2 = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1))
        d = 1.0 + aa_m2 * d
        if abs(d) < EPS:
            d = EPS
        c = 1.0 + aa_m2 / c
        if abs(c) < EPS:
            c = EPS
        d = 1.0 / d
        delta = d * c
        f *= delta
        if abs(delta - 1.0) < EPS:
            break

    return f


_sqrt = math.sqrt


@ChemMCPManager.register_tool
class TTestCalculator(BaseTool):
    """
    t 检验显著性分析工具。
    支持单样本、配对、独立样本（等方差/异方差）t检验，自动选择或手动指定类型。
    """
    __version__ = "0.1.0"
    name = "TTestCalculator"
    func_name = "perform_t_test"
    description = "Perform t-test significance analysis: one-sample, paired, independent (equal/unequal variance). Auto-detect test type."
    implementation_description = "Implements Student's t-test with Cohen's d effect size, confidence intervals, and full descriptive statistics. Supports auto mode that selects appropriate test based on input."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["t-test", "Statistics", "Hypothesis Testing", "Significance", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("sample_a", "list", "N/A", "First sample data (list of numbers)."),
        ("sample_b", "list", "", "Second sample data (for independent/paired tests)."),
        ("test_type", "str", "auto", "Test type: 'auto', 'one-sample', 'paired', 'independent'."),
        ("population_mean", "float", "0.0", "Population mean for one-sample test."),
        ("alpha", "float", "0.05", "Significance level."),
        ("equal_var", "bool", "True", "Assume equal variances for independent test (Welch's t-test if False)."),
        ("alternative", "str", "two-sided", "Alternative hypothesis: 'two-sided', 'greater', 'less'."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("test_type", "str", "Final test type performed."),
        ("t_statistic", "float", "t-statistic value."),
        ("p_value", "float", "p-value (two-tailed)."),
        ("degrees_of_freedom", "float", "Degrees of freedom."),
        ("mean_group_a", "float", "Mean of group A."),
        ("mean_group_b", "float", "Mean of group B (or population mean)."),
        ("mean_difference", "float", "Difference between means."),
        ("confidence_interval", "list", "(1−α) confidence interval for the difference."),
        ("effect_size", "float", "Cohen's d effect size."),
        ("reject_null", "bool", "Whether to reject H₀ at α level."),
        ("summary", "str", "Natural language summary."),
        ("descriptive_stats", "dict", "Descriptive statistics per group (n, mean, std)."),
    ]

    examples = [
        {
            "code_input": {
                "sample_a": [23.5, 24.1, 22.8, 23.9, 24.2, 23.7],
                "sample_b": [21.3, 22.1, 21.8, 22.5, 21.9, 22.3],
                "test_type": "independent",
                "alpha": 0.05,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"reject_null": True, "p_value": 0.002},
        },
        {
            "code_input": {
                "sample_a": [100.2, 99.8, 100.5, 99.7, 100.1, 99.9],
                "population_mean": 100.0,
                "test_type": "one-sample",
                "alpha": 0.05,
            },
            "output": {"reject_null": False, "p_value": 0.75},
            "text_input": {"params_str": "one-sample test"},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _desc(data: List[float]) -> dict:
        n = len(data)
        m = sum(data) / n
        s = math.sqrt(sum((x - m) ** 2 for x in data) / max(n - 1, 1)) if n > 1 else 0.0
        se = s / math.sqrt(n) if n > 0 else 0.0
        return {"n": n, "mean": round(m, 6), "std": round(s, 6), "se": round(se, 6),
                "min": round(min(data), 6), "max": round(max(data), 6)}

    def _run_base(
        self,
        sample_a: List[float],
        sample_b: Optional[List[float]] = None,
        test_type: str = "auto",
        population_mean: float = 0.0,
        alpha: float = 0.05,
        equal_var: bool = True,
        alternative: str = "two-sided",
    ) -> dict:
        """Core logic: perform t-test."""
        alt = alternative.lower().strip()
        tt = test_type.lower().strip()

        # Auto-detect test type
        if tt == "auto":
            if sample_b is None or not sample_b:
                tt = "one-sample"
            else:
                tt = "independent"

        desc_a = self._desc(sample_a)

        if tt == "one-sample":
            n = desc_a["n"]
            mean_a = desc_a["mean"]
            std_a = desc_a["std"]
            se_a = desc_a["se"]

            t_stat = (mean_a - population_mean) / se_a if se_a > 0 else 0.0
            df = n - 1
            p_val = _p_value_from_t(t_stat, df)
            mean_diff = mean_a - population_mean
            desc_b = {"n": 0, "mean": population_mean, "std": 0}

            # CI for mean
            t_crit = _t_critical_approx(df, alpha, two_tail=True)
            ci_low = mean_a - t_crit * se_a
            ci_high = mean_a + t_crit * se_a
            ci = [round(ci_low, 6), round(ci_high, 6)]

            # Cohen's d (one-sample)
            cohens_d = mean_diff / std_a if std_a > 0 else 0.0

        elif tt == "paired":
            if not sample_b or len(sample_b) != len(sample_a):
                raise ChemMCPError("Paired test requires two samples of equal length.")
            diffs = [a - b for a, b in zip(sample_a, sample_b)]
            desc_d = self._desc(diffs)
            n = desc_d["n"]
            t_stat = desc_d["mean"] / desc_d["se"] if desc_d["se"] > 0 else 0.0
            df = n - 1
            p_val = _p_value_from_t(t_stat, df)
            mean_diff = desc_d["mean"]
            desc_b = self._desc(sample_b)

            t_crit = _t_critical_approx(df, alpha, two_tail=True)
            ci_low = desc_d["mean"] - t_crit * desc_d["se"]
            ci_high = desc_d["mean"] + t_crit * desc_d["se"]
            ci = [round(ci_low, 6), round(ci_high, 6)]
            cohens_d = desc_d["mean"] / desc_d["std"] if desc_d["std"] > 0 else 0.0

        elif tt == "independent":
            if not sample_b:
                raise ChemMCPError("Independent t-test requires sample_b.")
            desc_b = self._desc(sample_b)
            na, nb = desc_a["n"], desc_b["n"]
            ma, mb = desc_a["mean"], desc_b["mean"]
            sa, sb = desc_a["std"], desc_b["std"]
            mean_diff = ma - mb

            if equal_var:
                # Pooled variance (Student's t)
                sp_sq = ((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2)
                se_diff = math.sqrt(sp_sq * (1/na + 1/nb))
                df = na + nb - 2
            else:
                # Welch's t-test
                se_diff = math.sqrt(sa**2/na + sb**2/nb)
                num = (sa**2/na + sb**2/nb)**2
                denom = ((sa**2/na)**2/(na-1) + (sb**2/nb)**2/(nb-1)) if (na > 1 and nb > 1) else 1
                df = num / denom if denom > 0 else max(na + nb - 2, 1)

            t_stat = mean_diff / se_diff if se_diff > 0 else 0.0
            p_val = _p_value_from_t(t_stat, df)

            t_crit = _t_critical_approx(df, alpha, two_tail=True)
            ci_low = mean_diff - t_crit * se_diff
            ci_high = mean_diff + t_crit * se_diff
            ci = [round(ci_low, 6), round(ci_high, 6)]

            # Cohen's d (pooled SD)
            pooled_sd = math.sqrt(sp_sq) if equal_var else math.sqrt((sa**2 + sb**2) / 2)
            cohens_d = mean_diff / pooled_sd if pooled_sd > 0 else 0.0

        else:
            raise ChemMCPError(f"Unknown test_type '{test_type}'. Use: auto, one-sample, paired, independent.")

        reject = p_val < alpha

        # Summary
        if reject:
            summary = (
                f"{tt.replace('_', '-').title()} t-test: t({df:.1f})={t_stat:.4f}, "
                f"p={p_val:.4f} < α={alpha} → Reject H₀. "
                f"Significant difference detected (mean diff = {mean_diff:.4g})."
            )
        else:
            summary = (
                f"{tt.replace('_', '-').title()} t-test: t({df:.1f})={t_stat:.4f}, "
                f"p={p_val:.4f} ≥ α={alpha} → Fail to reject H₀. "
                f"No significant difference (mean diff = {mean_diff:.4g})."
            )

        logger.info(f"t-test ({tt}): t={t_stat:.4f}, p={p_val:.4f}, df={df:.1f}, reject={reject}")
        return {
            "test_type": tt,
            "t_statistic": round(t_stat, 6),
            "p_value": round(p_val, 6),
            "degrees_of_freedom": round(df, 4),
            "mean_group_a": desc_a["mean"],
            "mean_group_b": desc_b.get("mean", population_mean),
            "mean_difference": round(mean_diff, 6),
            "confidence_interval": ci,
            "effect_size": round(cohens_d, 6),
            "reject_null": reject,
            "summary": summary,
            "descriptive_stats": {"group_a": desc_a, "group_b": desc_b},
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
