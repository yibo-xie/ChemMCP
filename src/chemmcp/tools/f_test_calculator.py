import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _f_critical_approx(df1: int, df2: int, alpha: float) -> float:
    """
    Approximate F-critical value using inverse CDF approximation.
    For common values, uses lookup table; otherwise approximates.
    """
    # Simplified lookup for common alpha values
    if alpha == 0.05:
        _f_table_005 = {
            (1, 5): 6.61, (1, 10): 4.96, (1, 20): 4.35, (1, 30): 4.17, (1, 60): 4.00, (1, 120): 3.92,
            (2, 5): 5.79, (2, 10): 4.10, (2, 20): 3.49, (2, 30): 3.32, (2, 60): 3.15, (2, 120): 3.07,
            (3, 5): 5.41, (3, 10): 3.71, (3, 20): 3.10, (3, 30): 2.92, (3, 60): 2.76, (3, 120): 2.68,
            (4, 5): 5.19, (4, 10): 3.48, (4, 20): 2.87, (4, 30): 2.69, (4, 60): 2.53, (4, 120): 2.45,
            (5, 5): 5.05, (5, 10): 3.33, (5, 20): 2.71, (5, 30): 2.53, (5, 60): 2.37, (5, 120): 2.29,
            (6, 5): 4.95, (6, 10): 3.22, (6, 20): 2.60, (6, 30): 2.42, (6, 60): 2.25, (6, 120): 2.18,
            (7, 5): 4.88, (7, 10): 3.14, (7, 20): 2.51, (7, 30): 2.33, (7, 60): 2.17, (7, 120): 2.09,
            (8, 5): 4.82, (8, 10): 3.07, (8, 20): 2.45, (8, 30): 2.27, (8, 60): 2.10, (8, 120): 2.02,
            (9, 5): 4.77, (9, 10): 3.02, (9, 20): 2.39, (9, 30): 2.21, (9, 60): 2.04, (9, 120): 1.97,
            (10, 5): 4.74, (10, 10): 2.98, (10, 20): 2.35, (10, 30): 2.16, (10, 60): 2.00, (10, 120): 1.92,
        }
    elif alpha == 0.01:
        _f_table_001 = {
            (1, 5): 16.26, (1, 10): 10.04, (1, 20): 8.10, (1, 30): 7.56, (1, 60): 7.08, (1, 120): 6.85,
            (2, 5): 13.27, (2, 10): 7.56, (2, 20): 5.85, (2, 30): 5.39, (2, 60): 4.98, (2, 120): 4.79,
            (3, 5): 12.06, (3, 10): 6.55, (3, 20): 4.94, (3, 30): 4.51, (3, 60): 4.13, (3, 120): 3.95,
            (4, 5): 11.39, (4, 10): 5.99, (4, 20): 4.43, (4, 30): 4.02, (4, 60): 3.65, (4, 120): 3.48,
            (5, 5): 10.97, (5, 10): 5.64, (5, 20): 4.10, (5, 30): 3.70, (5, 60): 3.34, (5, 120): 3.17,
            (6, 5): 10.67, (6, 10): 5.39, (6, 20): 3.87, (6, 30): 3.47, (6, 60): 3.12, (6, 120): 2.95,
            (7, 5): 10.46, (7, 10): 5.20, (7, 20): 3.70, (7, 30): 3.30, (7, 60): 2.95, (7, 120): 2.78,
            (8, 5): 10.29, (8, 10): 5.06, (8, 20): 3.56, (8, 30): 3.17, (8, 60): 2.82, (8, 120): 2.65,
            (9, 5): 10.16, (9, 10): 4.94, (9, 20): 3.45, (9, 30): 3.07, (9, 60): 2.72, (9, 120): 2.55,
            (10, 5): 10.05, (10, 10): 4.85, (10, 20): 3.37, (10, 30): 2.98, (10, 60): 2.64, (10, 120): 2.47,
        }
    elif alpha == 0.10:
        _f_table_010 = {
            (1, 5): 4.06, (1, 10): 3.29, (1, 20): 2.97, (1, 30): 2.84, (1, 60): 2.72, (1, 120): 2.66,
            (2, 5): 3.78, (2, 10): 2.92, (2, 20): 2.59, (2, 30): 2.46, (2, 60): 2.34, (2, 120): 2.28,
            (3, 5): 3.62, (3, 10): 2.73, (3, 20): 2.38, (3, 30): 2.25, (3, 60): 2.13, (3, 120): 2.07,
            (4, 5): 3.52, (4, 10): 2.61, (4, 20): 2.25, (4, 30): 2.12, (4, 60): 2.00, (4, 120): 1.93,
            (5, 5): 3.45, (5, 10): 2.52, (5, 20): 2.16, (5, 30): 2.03, (5, 60): 1.91, (5, 120): 1.84,
        }

    table = locals().get(f"_f_table_{str(alpha).replace('.', '')}", {})
    
    # Try exact match first
    key = (df1, df2)
    if key in table:
        return float(table[key])

    # Try swapped (F distribution is not symmetric but we can try nearest)
    # Find closest df pair and interpolate
    if table:
        keys = list(table.keys())
        best_key = min(keys, key=lambda k: abs(k[0]-df1)*0.5 + abs(k[1]-df2))
        base_val = table[best_key]
        # Rough adjustment
        ratio = math.sqrt(best_key[0] / df1) * math.sqrt(best_key[1] / df2) if df1 > 0 and df2 > 0 else 1.0
        return max(base_val * min(max(ratio, 0.8), 1.3), 1.01)

    # Fallback approximation using chi-squared / normal
    import math as _m
    z = 1.645 if alpha == 0.05 else (2.326 if alpha == 0.01 else 1.282)
    return (_m.sqrt(2.0 / (max(df1, 1) - 1)) * z + 1) ** 3


def _p_value_from_f(F_stat: float, df1: int, df2: int) -> float:
    """Approximate p-value from F-statistic."""
    import math as _m
    if F_stat <= 0 or df1 <= 0 or df2 <= 0:
        return 1.0
    x = df2 / (df2 + df1 * F_stat)
    p = _beta_reg(x, df2 / 2, df1 / 2)
    return 1 - p


# Reuse beta_reg from t_test module logic (inline it here)
def _beta_reg(x: float, a: float, b: float, max_iter: int = 200) -> float:
    import math as _m
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = _m.exp(
        _m.lgamma(a+b) - _m.lgamma(a) - _m.lgamma(b)
        + a*_m.log(x) + b*_m.log(max(1-x, 1e-30))
    )
    if x < (a+1)/(a+b+2):
        return bt * _beta_cf(x, a, b, max_iter) / a
    else:
        return 1 - bt * _beta_cf(1-x, b, a, max_iter) / b

def _beta_cf(x: float, a: float, b: float, max_iter: int = 200) -> float:
    import math as _m
    EPS = 1e-12
    f = c = d = 1.0
    d = 1.0 - (a+b)*x/(a+1)
    if abs(d) < EPS:
        d = EPS
    d = 1.0 / d
    f = d
    for m in range(1, max_iter+1):
        m2 = 2 * m
        aa_m = m*(b-m)*x / ((a+m2-1)*(a+m2))
        d = 1 + aa_m * d
        if abs(d) < EPS:
            d = EPS
        c = 1 + aa_m / c
        if abs(c) < EPS:
            c = EPS
        d = 1 / d
        f *= d * c
        aa_m2 = -(a+m)*(a+b+m)*x / ((a+m2)*(a+m2+1))
        d = 1 + aa_m2 * d
        if abs(d) < EPS:
            d = EPS
        c = 1 + aa_m2 / c
        if abs(c) < EPS:
            c = EPS
        d = 1 / d
        delta = d * c
        f *= delta
        if abs(delta-1) < EPS:
            break
    return f


@ChemMCPManager.register_tool
class FTestCalculator(BaseTool):
    """
    F 检验方差齐性检验工具。
    比较两组样本的方差是否相等，支持原始样本和摘要统计两种输入模式。
    """
    __version__ = "0.1.0"
    name = "FTestCalculator"
    func_name = "perform_f_test"
    description = "F-test for equality of variances between two groups. Supports raw samples or summary statistics input."
    implementation_description = "Implements two-sample F-test for variance homogeneity. F = s₁²/s₂² with larger variance in numerator. Supports both raw data input and summary statistic (variance/size) mode."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["F-test", "Variance", "Statistics", "Hypothesis Testing", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("sample1", "list", "", "First sample data (list of floats)."),
        ("sample2", "list", "", "Second sample data (list of floats)."),
        ("variance1", "float", "", "Variance of group 1 (alternative to sample1)."),
        ("variance2", "float", "", "Variance of group 2 (alternative to sample2)."),
        ("size1", "int", "", "Sample size of group 1 (required when using variance input)."),
        ("size2", "int", "", "Sample size of group 2 (required when using variance input)."),
        ("alpha", "float", "0.05", "Significance level."),
        ("alternative", "str", "two-sided", "Alternative hypothesis: 'two-sided', 'greater', 'less'."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("f_statistic", "float", "F-statistic value."),
        ("p_value", "float", "p-value."),
        ("df", "list", "Degrees of freedom [df1, df2]."),
        ("critical_values", "dict", "Critical F values at given alpha."),
        ("reject_null", "bool", "Whether to reject H₀ (equal variances)."),
        ("variance_summary", "dict", "Summary of each group's variance."),
        ("summary", "str", "Natural language interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "sample1": [12.5, 13.1, 12.8, 12.9, 13.0],
                "sample2": [14.2, 15.1, 14.8, 14.5, 14.9, 15.0],
                "alpha": 0.05,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"reject_null": False},
        },
        {
            "code_input": {
                "variance1": 0.25, "variance2": 1.44,
                "size1": 10, "size2": 10,
                "alpha": 0.05,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"reject_null": True, "f_statistic": 5.76},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _var(data: List[float], ddof: int = 1) -> float:
        n = len(data)
        m = sum(data) / n
        return sum((x - m) ** 2 for x in data) / (n - ddof)

    @staticmethod
    def _mean(data: List[float]) -> float:
        return sum(data) / len(data)

    def _run_base(
        self,
        sample1: Optional[List[float]] = None,
        sample2: Optional[List[float]] = None,
        variance1: Optional[float] = None,
        variance2: Optional[float] = None,
        size1: Optional[int] = None,
        size2: Optional[int] = None,
        alpha: float = 0.05,
        alternative: str = "two-sided",
    ) -> dict:
        """Core logic: perform F-test for variance equality."""
        alt = alternative.lower().strip()

        # Determine input mode
        use_samples = sample1 is not None and sample2 is not None

        if use_samples:
            v1 = self._var(sample1)
            v2 = self._var(sample2)
            n1 = len(sample1)
            n2 = len(sample2)
            mean1 = self._mean(sample1)
            mean2 = self._mean(sample2)
        elif variance1 is not None and variance2 is not None:
            if size1 is None or size2 is None:
                raise ChemMCPError("When using variance input, both size1 and size2 are required.")
            v1 = float(variance1)
            v2 = float(variance2)
            n1 = int(size1)
            n2 = int(size2)
            mean1 = mean2 = None
        else:
            raise ChemMCPError("Provide either (sample1, sample2) or (variance1, variance2, size1, size2).")

        if n1 < 2 or n2 < 2:
            raise ChemMCPError("Each group must have at least 2 observations (df ≥ 1).")
        if v1 <= 0 or v2 <= 0:
            raise ChemMCPError("Variances must be positive.")

        # F = larger variance / smaller variance
        if v1 >= v2:
            F_stat = v1 / v2
            df1 = n1 - 1
            df2 = n2 - 1
            larger_group = 1
        else:
            F_stat = v2 / v1
            df1 = n2 - 1
            df2 = n1 - 1
            larger_group = 2

        p_val = _p_value_from_f(F_stat, df1, df2)
        crit = _f_critical_approx(df1, df2, alpha)

        # Decision
        if alt == "two-sided":
            reject = p_val < alpha
        elif alt == "greater":
            reject = F_stat > crit
        else:  # less
            reject = F_stat < (1 / crit) if crit > 0 else False

        var_summary = {
            "group_1": {"variance": round(v1, 8), "std": round(math.sqrt(v1), 8), "n": n1, "mean": round(mean1, 6) if mean1 else None},
            "group_2": {"variance": round(v2, 8), "std": round(math.sqrt(v2), 8), "n": n2, "mean": round(mean2, 6) if mean2 else None},
        }

        if reject:
            summary = (
                f"F-test: F({df1},{df2})={F_stat:.4f}, p={p_val:.4f}, "
                f"F_crit({alpha})={crit:.4f} → Reject H₀. "
                f"Variances are significantly different (unequal). "
                f"Group {larger_group} has significantly larger variance."
            )
        else:
            summary = (
                f"F-test: F({df1},{df2})={F_stat:.4f}, p={p_val:.4f}, "
                f"F_crit({alpha})={crit:.4f} → Fail to reject H₀. "
                f"No significant evidence of unequal variances. Variance homogeneity assumption holds."
            )

        logger.info(f"F-test: F({df1},{df2})={F_stat:.4f}, p={p_val:.4f}, reject={reject}")
        return {
            "f_statistic": round(F_stat, 6),
            "p_value": round(p_val, 6),
            "df": [df1, df2],
            "critical_values": {"alpha": alpha, "F_critical": round(crit, 6)},
            "reject_null": reject,
            "variance_summary": var_summary,
            "summary": summary,
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
