import logging
import math
from typing import List, Optional

import numpy as np
from scipy import stats

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MethodComparisonEvaluator(BaseTool):
    """
    方法比对评估工具：配对t检验 + Bland-Altman分析。
    
    用于比较两种分析方法（如新方法vs参比方法）的测量结果一致性。
    """
    __version__ = "0.1.0"
    name             = "MethodComparisonEvaluator"
    func_name        = "compare_methods"
    description      = "Evaluate agreement between two analytical methods using paired t-test and Bland-Altman analysis."
    implementation_description = "Implements paired Student's t-test (scipy.stats.ttest_rel) and Bland-Altman analysis (mean difference, LOA, bias) for method comparison in analytical chemistry validation."
    oss_dependencies = [
        ("NumPy", "https://numpy.org", "BSD"),
        ("SciPy", "https://scipy.org", "BSD"),
    ]
    services_and_software = []
    categories       = ["General"]
    tags             = ["Statistics", "Method Validation", "Analytical Chemistry", "QA/QC"]
    required_envs    = []

    code_input_sig   = [
        ("method_a", "list", "N/A", "List of measurements from Method A (reference method)."),
        ("method_b", "list", "N/A", "List of measurements from Method B (new method)."),
        ("confidence_level", "float", "0.95", "Confidence level for intervals (e.g., 0.95 for 95%)."),
        ("alpha", "float", "0.05", "Significance level for t-test (e.g., 0.05)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'method_a_values method_b_values [confidence_level] [alpha]'. Lists use commas."),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary containing paired_t_test, bland_altman, and conclusion."),
    ]

    examples         = [
        {
            "code_input": {
                "method_a": [10.2, 10.5, 9.8, 10.1, 10.3, 9.9, 10.4, 10.0],
                "method_b": [10.1, 10.6, 9.7, 10.2, 10.4, 9.8, 10.5, 10.1],
                "confidence_level": 0.95,
                "alpha": 0.05,
            },
            "text_input": {
                "input_params": "10.2,10.5,9.8,10.1,10.3,9.9,10.4,10.0 10.1,10.6,9.7,10.2,10.4,9.8,10.5,10.1 0.95 0.05",
            },
            "output": {
                "result": {
                    "paired_t_test": {"t_statistic": "...", "p_value": "...", "df": "...", "significant": "..."},
                    "bland_altman": {"bias": "...", "mean_difference": "...", "sd_diff": "...", "loa_upper": "...", "loa_lower": "..."},
                    "conclusion": "..."
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        method_a: List[float],
        method_b: List[float],
        confidence_level: float = 0.95,
        alpha: float = 0.05,
    ) -> dict:
        """
        核心逻辑：执行配对t检验和Bland-Altman分析
        """
        # 输入验证
        if len(method_a) != len(method_b):
            raise ChemMCPError("Method A and Method B must have the same number of measurements.")
        if len(method_a) < 3:
            raise ChemMCPError("At least 3 pairs of measurements are required.")

        a = np.array(method_a, dtype=float)
        b = np.array(method_b, dtype=float)
        n = len(a)

        # ---- 1. 配对 t 检验 ----
        t_stat, p_value = stats.ttest_rel(a, b)
        df = n - 1
        significant = p_value < alpha

        # 置信区间（均值差）
        diff = a - b
        mean_diff = float(np.mean(diff))
        se = float(np.std(diff, ddof=1) / math.sqrt(n))
        t_crit = float(stats.t.ppf((1 + confidence_level) / 2, df))
        ci_low = round(mean_diff - t_crit * se, 6)
        ci_high = round(mean_diff + t_crit * se, 6)

        # ---- 2. Bland-Altman 分析 ----
        mean_ab = (a + b) / 2.0  # 两方法均值
        sd_diff = float(np.std(diff, ddof=1))

        # Limits of Agreement (LOA): bias ± 1.96 * SD
        loa_upper = round(mean_diff + 1.96 * sd_diff, 6)
        loa_lower = round(mean_diff - 1.96 * sd_diff, 6)

        # LOA 的置信区间
        se_loa = math.sqrt(1.0 / n + 1.96**2 / (n - 1)) * sd_diff if n > 1 else 0
        loa_ci_upper = round(loa_upper + t_crit * se_loa, 6)
        loa_ci_lower = round(loa_lower - t_crit * se_loa, 6)

        # ---- 3. 回归分析（Passing-Bablok 可选：简单线性回归）----
        slope, intercept, r_value, p_reg, std_err = stats.linregress(a, b)

        # ---- 4. 结论 ----
        if not significant and abs(mean_diff) < 0.1 * abs(np.mean(a)):
            conclusion = "No significant systematic difference between methods; good agreement."
        elif not significant:
            conclusion = "No statistically significant difference, but bias may be clinically/analytically relevant."
        else:
            conclusion = f"Significant difference detected (p={p_value:.4g}<{alpha}). Methods may not be interchangeable."

        result = {
            "sample_size": n,
            "paired_t_test": {
                "t_statistic": round(float(t_stat), 6),
                "p_value": round(float(p_value), 6),
                "df": int(df),
                "alpha": alpha,
                "significant": bool(significant),
                "mean_diff_ci": [ci_low, ci_high],
                "confidence_level": confidence_level,
            },
            "bland_altman": {
                "bias": round(mean_diff, 6),
                "sd_difference": round(sd_diff, 6),
                "loa_upper": loa_upper,
                "loa_lower": loa_lower,
                "loa_ci_upper": loa_ci_upper,
                "loa_ci_lower": loa_ci_lower,
            },
            "regression": {
                "slope": round(float(slope), 6),
                "intercept": round(float(intercept), 6),
                "r_squared": round(float(r_value ** 2), 6),
                "pearson_r": round(float(r_value), 6),
            },
            "conclusion": conclusion,
        }

        logger.info(f"Method comparison completed: n={n}, p={p_value:.4g}, bias={mean_diff:.4f}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need at least two lists of values.")

            a_vals = [float(x) for x in parts[0].split(",")]
            b_vals = [float(x) for x in parts[1].split(",")]
            conf = float(parts[2]) if len(parts) > 2 else 0.95
            alp = float(parts[3]) if len(parts) > 3 else 0.05

            return self._run_base(a_vals, b_vals, conf, alp)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'a1,a2,a3 b1,b2,b3 [conf] [alpha]'")
