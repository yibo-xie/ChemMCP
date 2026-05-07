import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _f_critical_anova(df1: int, df2: int, alpha: float = 0.05) -> float:
    """Approximate F-critical value for ANOVA."""
    # Reuse approximation
    return _f_critical_approx(df1, df2, alpha)


def _f_critical_approx(df1: int, df2: int, alpha: float = 0.05) -> float:
    """Approximate F-critical value."""
    if alpha == 0.05:
        table = {
            (1,5):6.61,(1,10):4.96,(1,20):4.35,(1,30):4.17,(1,60):4.00,(1,120):3.92,
            (2,5):5.79,(2,10):4.10,(2,20):3.49,(2,30):3.32,(2,60):3.15,(2,120):3.07,
            (3,5):5.41,(3,10):3.71,(3,20):3.10,(3,30):2.92,(3,60):2.76,(3,120):2.68,
            (4,5):5.19,(4,10):3.48,(4,20):2.87,(4,30):2.69,(4,60):2.53,(4,120):2.45,
            (5,5):5.05,(5,10):3.33,(5,20):2.71,(5,30):2.53,(5,60):2.37,(5,120):2.29,
            (6,5):4.95,(6,10):3.22,(6,20):2.60,(6,30):2.42,(6,60):2.25,(6,120):2.18,
            (7,5):4.88,(7,10):3.14,(7,20):2.51,(7,30):2.33,(7,60):2.17,(7,120):2.09,
            (8,5):4.82,(8,10):3.07,(8,20):2.45,(8,30):2.27,(8,60):2.10,(8,120):2.02,
            (9,5):4.77,(9,10):3.02,(9,20):2.39,(9,30):2.21,(9,60):2.04,(9,120):1.97,
            (10,5):4.74,(10,10):2.98,(10,20):2.35,(10,30):2.16,(10,60):2.00,(10,120):1.92,
            (15,5):4.62,(15,10):2.86,(15,20):2.21,(15,30):2.03,(15,60):1.86,(15,120):1.79,
            (20,5):4.56,(20,10):2.77,(20,20):2.12,(20,30):1.93,(20,60):1.75,(20,120):1.68,
            (30,5):4.50,(30,10):2.70,(30,20):2.04,(30,30):1.85,(30,60):1.67,(30,120):1.59,
            (50,5):4.44,(50,10):2.64,(50,20):1.98,(50,30):1.78,(50,60):1.60,(50,120):1.52,
            (100,5):4.39,(100,10):2.58,(100,20):1.91,(100,30):1.72,(100,60):1.54,(100,120):1.46,
        }
    elif alpha == 0.01:
        table = {
            (1,5):16.26,(1,10):10.04,(1,20):8.10,(1,30):7.56,(1,60):7.08,(1,120):6.85,
            (2,5):13.27,(2,10):7.56,(2,20):5.85,(2,30):5.39,(2,60):4.98,(2,120):4.79,
            (3,5):12.06,(3,10):6.55,(3,20):4.94,(3,30):4.51,(3,60):4.13,(3,120):3.95,
            (4,5):11.39,(4,10):5.99,(4,20):4.43,(4,30):4.02,(4,60):3.65,(4,120):3.48,
            (5,5):10.97,(5,10):5.64,(5,20):4.10,(5,30):3.70,(5,60):3.34,(5,120):3.17,
            (10,5):10.05,(10,10):4.85,(10,20):3.37,(10,30):2.98,(10,60):2.64,(10,120):2.47,
            (20,5):9.44,(20,10):4.43,(20,20):2.94,(20,30):2.55,(20,60):2.20,(20,120):2.02,
            (30,5):9.18,(30,10):4.21,(30,20):2.73,(30,30):2.34,(30,60):1.99,(30,120):1.81,
        }
    else:  # 0.10
        table = {
            (1,5):4.06,(1,10):3.29,(1,20):2.97,(1,30):2.84,(1,60):2.72,(1,120):2.66,
            (2,5):3.78,(2,10):2.92,(2,20):2.59,(2,30):2.46,(2,60):2.34,(2,120):2.28,
            (3,5):3.62,(3,10):2.73,(3,20):2.38,(3,30):2.25,(3,60):2.13,(3,120):2.07,
            (5,5):3.45,(5,10):2.52,(5,20):2.16,(5,30):2.03,(5,60):1.91,(5,120):1.84,
            (10,5):4.74,(10,10):2.98,(10,20):2.35,(10,30):2.16,(10,60):2.00,(10,120):1.92,
        }

    key = (df1, df2)
    if key in table:
        return float(table[key])
    
    keys = list(table.keys())
    best = min(keys, key=lambda k: abs(k[0]-df1)*0.5 + abs(k[1]-df2))
    base_val = table[best]
    ratio = math.sqrt(best[0]/max(df1,1)) * math.sqrt(best[1]/max(df2,1))
    return max(base_val * min(max(ratio, 0.8), 1.3), 1.01)


def _p_value_from_f(F_stat: float, df1: int, df2: int) -> float:
    """Approximate p-value from F-statistic."""
    if F_stat <= 0 or df1 <= 0 or df2 <= 0:
        return 1.0
    x = df2 / (df2 + df1 * F_stat)
    p = _beta_reg(x, df2/2, df1/2)
    return 1 - p


def _beta_reg(x: float, a: float, b: float, max_iter: int = 200) -> float:
    import math as _m
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = _m.exp(_m.lgamma(a+b)-_m.lgamma(a)-_m.lgamma(b)+a*_m.log(x)+b*_m.log(max(1-x,1e-30)))
    if x < (a+1)/(a+b+2):
        return bt * _beta_cf(x, a, b, max_iter) / a
    return 1 - bt * _beta_cf(1-x, b, a, max_iter) / b

def _beta_cf(x: float, a: float, b: float, max_iter: int = 200) -> float:
    import math as _m
    EPS = 1e-12; f = c = d = 1.0
    d = 1-(a+b)*x/(a+1)
    if abs(d)<EPS: d=EPS
    d = 1 / d
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
class AnovaAnalyzer(BaseTool):
    """
    方差分析工具。
    支持单因素方差分析（one-way ANOVA）和双因素方差分析（two-way ANOVA，含交互项）。
    """
    __version__ = "0.1.0"
    name = "AnovaAnalyzer"
    func_name = "perform_anova"
    description = "Analysis of variance (ANOVA): one-way and two-way (with interaction term)."
    implementation_description = "Implements one-way ANOVA for comparing means across multiple groups, and two-way ANOVA for analyzing effects of two factors with interaction. Outputs full ANOVA table with SS, MS, F, p-values."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["ANOVA", "Variance Analysis", "Statistics", "Experimental Design", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "one_way", "ANOVA mode: 'one_way' or 'two_way'."),
        ("groups", "list", "", "List of groups (each group is a list of numbers). Required for one_way."),
        ("factor_a", "list", "", "Factor A level labels for each observation (for two_way)."),
        ("factor_b", "list", "", "Factor B level labels for each observation (for two_way)."),
        ("values", "list", "", "Observed values corresponding to factor_a/factor_b (for two_way)."),
        ("alpha", "float", "0.05", "Significance level."),
        ("group_labels", "list", "", "Optional labels for groups (one_way mode)."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("mode", "str", "ANOVA mode performed."),
        ("anova_table", "dict", "Full ANOVA table with Source, DF, SS, MS, F, p-value, significant at α."),
        ("f_critical", "dict", "Critical F values at given α."),
        ("conclusion", "str", "Natural language conclusion."),
        ("summary_stats", "dict", "Group-level summary statistics (n, mean, std)."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "one_way",
                "groups": [[23.5, 24.1, 22.8], [21.3, 22.1, 21.8, 22.5], [25.0, 24.8, 25.2, 24.9]],
                "alpha": 0.05,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"anova_table": {}},
        },
        {
            "code_input": {
                "mode": "two_way",
                "factor_a": ["A1","A1","A1","A1","A2","A2","A2","A2"],
                "factor_b": ["B1","B1","B2","B2","B1","B1","B2","B2"],
                "values": [10,12,15,18,14,16,20,23],
                "alpha": 0.05,
            },
            "text_input": {"params_str": "two_way input"},
            "output": {"anova_table": {}},
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
    def _var(data: List[float], ddof: int = 0) -> float:
        n = len(data)
        m = sum(data) / n
        return sum((x - m) ** 2 for x in data) / max(n - ddof, 1)

    def _one_way_anova(self, groups: List[List[float]], alpha: float, labels: Optional[List[str]] = None) -> dict:
        """Perform one-way ANOVA."""
        k = len(groups)
        if k < 2:
            raise ChemMCPError("One-way ANOVA requires at least 2 groups.")

        all_data = []
        group_means = []
        group_sizes = []
        group_stds = []

        for g in groups:
            all_data.extend(g)
            group_means.append(self._mean(g))
            group_sizes.append(len(g))
            group_stds.append(math.sqrt(self._var(g, ddof=1)))

        N = len(all_data)
        grand_mean = self._mean(all_data)

        # Sum of Squares
        ss_between = sum(ni * (mi - grand_mean)**2 for ni, mi in zip(group_sizes, group_means))
        ss_total = sum((x - grand_mean)**2 for x in all_data)
        ss_within = ss_total - ss_between

        # Degrees of freedom
        df_between = k - 1
        df_within = N - k
        df_total = N - 1

        # Mean squares
        ms_between = ss_between / df_between if df_between > 0 else 0
        ms_within = ss_within / df_within if df_within > 0 else 0

        # F-statistic
        f_stat = ms_between / ms_within if ms_within > 0 else float('inf')

        # p-value
        p_val = _p_value_from_f(f_stat, df_between, df_within)

        # Critical F
        f_crit = _f_critical_anova(df_between, df_within, alpha)

        sig = p_val < alpha

        anova_table = {
            "between_groups": {
                "source": "Between Groups",
                "df": df_between,
                "ss": round(ss_between, 6),
                "ms": round(ms_between, 6),
                "f": round(f_stat, 4),
                "p_value": round(p_val, 6),
                "significant_at_0_05": sig,
            },
            "within_groups": {
                "source": "Within Groups (Error)",
                "df": df_within,
                "ss": round(ss_within, 6),
                "ms": round(ms_within, 6),
                "f": None,
                "p_value": None,
                "significant_at_0_05": None,
            },
            "total": {
                "source": "Total",
                "df": df_total,
                "ss": round(ss_total, 6),
                "ms": None,
                "f": None,
                "p_value": None,
                "significant_at_0_05": None,
            },
        }

        f_critical_dict = {"F_critical(α={})".format(alpha): round(f_crit, 4)}

        label_list = labels or [f"Group_{i+1}" for i in range(k)]
        summary_stats = {label_list[i]: {"n": group_sizes[i], "mean": round(group_means[i], 6), "std": round(group_stds[i], 6)}
                         for i in range(k)}

        conclusion_parts = []
        if sig:
            conclusion_parts.append(
                f"Between-group difference is statistically significant "
                f"(F({df_between},{df_within})={f_stat:.3f}, p={p_val:.4f} < {alpha})."
            )
            # Find which groups differ most
            max_idx = group_means.index(max(group_means))
            min_idx = group_means.index(min(group_means))
            conclusion_parts.append(
                f"Group '{label_list[max_idx]}' (mean={group_means[max_idx]:.3g}) differs from "
                f"Group '{label_list[min_idx]}' (mean={group_means[min_idx]:.3g})."
            )
        else:
            conclusion_parts.append(
                f"No significant between-group difference detected "
                f"(F({df_between},{df_within})={f_stat:.3f}, p={p_val:.4f} ≥ {alpha})."
            )

        return {
            "anova_table": anova_table,
            "f_critical": f_critical_dict,
            "conclusion": "\n".join(conclusion_parts),
            "summary_stats": summary_stats,
        }

    def _two_way_anova(self, factor_a: List[str], factor_b: List[str], values: List[float], alpha: float) -> dict:
        """Perform two-way ANOVA with interaction."""
        n = len(values)
        if n != len(factor_a) or n != len(factor_b):
            raise ChemMCPError("Lengths of factor_a, factor_b, and values must all be equal.")

        levels_a = sorted(set(factor_a))
        levels_b = sorted(set(factor_b))
        a = len(levels_a)
        b = len(levels_b)

        # Build cell matrix
        cell_sums = {}
        cell_counts = {}
        row_sums = {}
        row_counts = {}
        col_sums = {}
        col_counts = {}

        for i in range(n):
            fa, fb, v = factor_a[i], factor_b[i], values[i]
            cell_key = (fa, fb)
            cell_sums[cell_key] = cell_sums.get(cell_key, 0) + v
            cell_counts[cell_key] = cell_counts.get(cell_key, 0) + 1
            row_sums[fa] = row_sums.get(fa, 0) + v
            row_counts[fa] = row_counts.get(fa, 0) + 1
            col_sums[fb] = col_sums.get(fb, 0) + v
            col_counts[fb] = col_counts.get(fb, 0) + 1

        grand_sum = sum(values)
        grand_mean = grand_sum / n

        # Calculate sums of squares
        ss_total = sum((v - grand_mean)**2 for v in values)

        # SSA (Factor A main effect)
        ssa = sum(row_sums[la]**2 / row_counts[la] for la in levels_a) - grand_sum**2 / n

        # SSB (Factor B main effect)
        ssb = sum(col_sums[lb]**2 / col_counts[lb] for lb in levels_b) - grand_sum**2 / n

        # SSAB (Interaction)
        ssab = 0
        for la in levels_a:
            for lb in levels_b:
                ck = (la, lb)
                if ck in cell_counts and cell_counts[ck] > 0:
                    ssab += cell_sums[ck]**2 / cell_counts[ck]
        ssab -= grand_sum**2 / n
        ssab -= ssa
        ssab -= ssb

        # SSE (Error)
        sse = 0
        for i in range(n):
            fa, fb, v = factor_a[i], factor_b[i], values[i]
            ck = (fa, fb)
            cell_mean = cell_sums[ck] / cell_counts[ck]
            sse += (v - cell_mean)**2

        # Degrees of freedom
        df_a = a - 1
        df_b = b - 1
        df_ab = (a - 1) * (b - 1)
        df_e = n - a * b
        df_t = n - 1

        # Mean squares
        ms_a = ssa / df_a if df_a > 0 else 0
        ms_b = ssb / df_b if df_b > 0 else 0
        ms_ab = ssab / df_ab if df_ab > 0 else 0
        ms_e = sse / df_e if df_e > 0 else 0

        # F statistics
        f_a = ms_a / ms_e if ms_e > 0 else float('inf')
        f_b = ms_b / ms_e if ms_e > 0 else float('inf')
        f_ab = ms_ab / ms_e if ms_e > 0 else float('inf')

        # p-values
        p_a = _p_value_from_f(f_a, df_a, df_e)
        p_b = _p_value_from_f(f_b, df_b, df_e)
        p_ab = _p_value_from_f(f_ab, df_ab, df_e)

        # Critical values
        fc_a = _f_critical_anova(df_a, df_e, alpha)
        fc_b = _f_critical_anova(df_b, df_e, alpha)
        fc_ab = _f_critical_anova(df_ab, df_e, alpha)

        def make_row(source, df, ss, ms, f, pv, fc):
            return {
                "source": source, "df": df, "ss": round(ss, 6),
                "ms": round(ms, 6) if ms is not None else None,
                "f": round(f, 4) if f is not None else None,
                "p_value": round(pv, 6) if pv is not None else None,
                "significant_at_0_05": pv < alpha if pv is not None else None,
            }

        anova_table = {
            "factor_A": make_row("Factor A", df_a, ssa, ms_a, f_a, p_a, fc_a),
            "factor_B": make_row("Factor B", df_b, ssb, ms_b, f_b, p_b, fc_b),
            "interaction_AB": make_row("Interaction A×B", df_ab, ssab, ms_ab, f_ab, p_ab, fc_ab),
            "error": make_row("Error", df_e, sse, ms_e, None, None, None),
            "total": make_row("Total", df_t, ss_total, None, None, None, None),
        }

        f_critical = {
            f"F_A(α={alpha})": round(fc_a, 4),
            f"F_B(α={alpha})": round(fc_b, 4),
            f"F_AB(α={alpha})": round(fc_ab, 4),
        }

        parts = []
        if p_a < alpha:
            parts.append(f"Factor A has a significant effect (F={f_a:.3f}, p={p_a:.4f}).")
        else:
            parts.append(f"Factor A has NO significant effect (F={f_a:.3f}, p={p_a:.4f}).")
        if p_b < alpha:
            parts.append(f"Factor B has a significant effect (F={f_b:.3f}, p={p_b:.4f}).")
        else:
            parts.append(f"Factor B has NO significant effect (F={f_b:.3f}, p={p_b:.4f}).")
        if p_ab < alpha:
            parts.append(f"Interaction A×B is significant (F={f_ab:.3f}, p={p_ab:.4f}).")
        else:
            parts.append(f"Interaction A×B is NOT significant (F={f_ab:.3f}, p={p_ab:.4f}).")

        return {
            "anova_table": anova_table,
            "f_critical": f_critical,
            "conclusion": "\n".join(parts),
            "summary_stats": {
                "levels_A": levels_a, "levels_B": levels_b,
                "n_per_cell": {ck: cell_counts.get(ck, 0) for ck in [(la, lb) for la in levels_a for lb in levels_b]},
            },
        }

    def _run_base(
        self,
        mode: str = "one_way",
        groups: Optional[List[List[float]]] = None,
        factor_a: Optional[List[str]] = None,
        factor_b: Optional[List[str]] = None,
        values: Optional[List[float]] = None,
        alpha: float = 0.05,
        group_labels: Optional[List[str]] = None,
    ) -> dict:
        """Core logic: perform ANOVA."""
        m = mode.lower().strip()

        if m == "one_way":
            if not groups or len(groups) < 2:
                raise ChemMCPError("One-way ANOVA requires 'groups' with at least 2 groups.")
            result = self._one_way_anova(groups, alpha, labels=group_labels)
        elif m == "two_way":
            if not factor_a or not factor_b or not values:
                raise ChemMCPError("Two-way ANOVA requires 'factor_a', 'factor_b', and 'values'.")
            result = self._two_way_anova(factor_a, factor_b, values, alpha)
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use 'one_way' or 'two_way'.")

        result["mode"] = m
        logger.info(f"ANOVA ({m}): completed")
        return result

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
