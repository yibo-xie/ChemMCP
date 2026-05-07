import logging
import math
from typing import List, Dict, Optional

import numpy as np
from scipy import stats as sp_stats

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MeasurementPrecisionCalculator(BaseTool):
    """
    精密度计算工具：重复性（repeatability）和再现性（reproducibility）。
    
    支持单因素方差分析(ANOVA)分解变异来源，计算不同浓度水平的精密度。
    """
    __version__ = "0.1.0"
    name             = "MeasurementPrecisionCalculator"
    func_name        = "calculate_precision"
    description      = "Calculate analytical precision: repeatability (intra-day) and reproducibility (inter-day) with ANOVA-based variance decomposition."
    implementation_description = "Computes SD, RSD% for each concentration level, performs one-way ANOVA to separate within-run and between-run variance components, following ICH Q2(R1) guidelines."
    oss_dependencies = [
        ("NumPy", "https://numpy.org", "BSD"),
        ("SciPy", "https://scipy.org", "BSD"),
    ]
    services_and_software = []
    categories       = ["General"]
    tags             = ["Precision", "Repeatability", "Reproducibility", "ICH Q2", "Method Validation"]
    required_envs    = []

    code_input_sig   = [
        ("measurements_by_run", "list_of_lists", "N/A", "List of runs, each run is a list of measurement values. e.g., [[run1_vals], [run2_vals], ...]"),
        ("concentration_levels_or_None", "list_or_None", "None", "Optional nominal concentrations for each run group."),
        ("acceptance_rsd", "float", "15.0", "Maximum acceptable RSD% for passing criteria (ICH typical: 15%, LLOQ: 20%)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated runs, comma-separated values: 'v1,v2,v3; v4,v5,v6; v7,v8,v9 [rsd_limit]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with precision statistics by level, overall precision, ANOVA table, and acceptance judgment."),
    ]

    examples         = [
        {
            "code_input": {
                "measurements_by_run": [
                    [10.1, 10.3, 9.9, 10.2, 10.0, 10.1],   # Day 1 / Run 1
                    [10.4, 10.1, 10.3, 10.0, 10.2, 9.98],   # Day 2 / Run 2
                    [9.95, 10.25, 10.05, 9.85, 10.15, 10.1], # Day 3 / Run 3
                ],
                "concentration_levels_or_None": None,
                "acceptance_rsd": 15.0,
            },
            "text_input": {
                "input_params": "10.1,10.3,9.9,10.2,10.0,10.1; 10.4,10.1,10.3,10.0,10.2,9.98; 9.95,10.25,10.05,9.85,10.15,10.1 15",
            },
            "output": {
                "result": {
                    "overall_precision": {"sd": "...", "rsd_pct": "..."},
                    "repeatability": {...},
                    "reproducibility": {...},
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
        measurements_by_run: List[List[float]],
        concentration_levels_or_None: Optional[List[float]] = None,
        acceptance_rsd: float = 15.0,
    ) -> dict:
        """核心逻辑：精密度计算 + ANOVA"""
        if len(measurements_by_run) < 2:
            raise ChemMCPError("At least 2 runs are required for precision evaluation.")

        runs = [np.array(r, dtype=float) for r in measurements_by_run]
        n_runs = len(runs)
        all_data = np.concatenate(runs)
        grand_mean = float(np.mean(all_data))
        n_total = len(all_data)

        # ---- 各运行统计 ----
        run_stats = []
        all_sds = []
        all_means = []
        for i, r in enumerate(runs):
            m = float(np.mean(r))
            s = float(np.std(r, ddof=1))
            rsd = (s / m * 100.0) if m != 0 else float('inf')
            accepted = rsd <= acceptance_rsd
            run_stats.append({
                "run": int(i + 1),
                "n": len(r),
                "mean": round(m, 6),
                "sd": round(s, 6),
                "rsd_pct": round(rsd, 4),
                "accepted": accepted,
            })
            all_sds.append(s)
            all_means.append(m)

        # ---- 复合精密度（所有数据合并） ----
        overall_sd = float(np.std(all_data, ddof=1))
        overall_mean = float(np.mean(all_data))
        overall_rsd = (overall_sd / overall_mean * 100.0) if overall_mean != 0 else float('inf')

        # ---- 单因素 ANOVA ----
        ss_between = sum(len(r) * (float(np.mean(r)) - grand_mean) ** 2 for r in runs)
        ss_within = sum(float(np.sum((r - np.mean(r)) ** 2)) for r in runs)
        df_between = n_runs - 1
        df_within = n_total - n_runs
        ms_between = ss_between / df_between if df_between > 0 else 0
        ms_within = ss_within / df_within if df_within > 0 else 0
        f_stat = ms_between / ms_within if ms_within > 0 else float('inf')
        p_value_anova = 1.0 - sp_stats.f.cdf(f_stat, df_between, df_within) if f_stat != float('inf') else 0.0

        anova_table = {
            "source_between": {"ss": round(ss_between, 8), "df": int(df_between), "ms": round(ms_between, 8)},
            "source_within": {"ss": round(ss_within, 8), "df": int(df_within), "ms": round(ms_within, 8)},
            "f_statistic": round(f_stat, 6),
            "p_value": round(p_value_anova, 6),
        }

        # ---- 变异分量 ----
        # σ²_within (repeatability variance)
        var_within = ms_within
        # σ²_between (between-run variance component)
        n_per_group = np.mean([len(r) for r in runs])
        var_between = max((ms_within - var_within) / n_per_group, 0) if n_per_group > 0 else 0

        sd_repeatability = math.sqrt(var_within) if var_within >= 0 else 0
        sd_reproducibility = math.sqrt(var_within + var_between) if (var_within + var_between) >= 0 else 0
        rsd_reproducibility = (sd_reproducibility / overall_mean * 100.0) if overall_mean != 0 else float('inf')
        rsd_repeatability = (sd_repeatability / overall_mean * 100.0) if overall_mean != 0 else float('inf')

        # ---- 判断 ----
        overall_accepted = overall_rsd <= acceptance_rsd

        result = {
            "summary": {
                "n_runs": n_runs,
                "n_total_measurements": n_total,
                "grand_mean": round(grand_mean, 6),
                "overall_sd": round(overall_sd, 6),
                "overall_rsd_pct": round(overall_rsd, 4),
                "accepted": overall_accepted,
                "acceptance_criteria_rsd_pct": acceptance_rsd,
            },
            "repeatability": {
                "description": "Within-run (intra-day) precision",
                "sd": round(sd_repeatability, 6),
                "rsd_pct": round(rsd_repeatability, 4),
            },
            "reproducibility": {
                "description": "Between-run (inter-day) precision",
                "sd": round(sd_reproducibility, 6),
                "rsd_pct": round(rsd_reproducibility, 4),
            },
            "variance_components": {
                "within_run_variance": round(var_within, 8),
                "between_run_variance_component": round(var_between, 8),
            },
            "anova_table": anova_table,
            "by_run_statistics": run_stats,
            "conclusion": (
                f"Overall RSD = {overall_rsd:.2f}%. "
                f"Repeatability RSD = {rsd_repeatability:.2f}%, Reproducibility RSD = {rsd_reproducibility:.2f}%. "
                f"{'PASS' if overall_accepted else 'FAIL'} (criteria: RSD ≤ {acceptance_rsd}%)."
            ),
        }

        logger.info(f"Precision calculation: overall RSD={overall_rsd:.2f}%, F={f_stat:.3f}, p={p_value_anova:.4g}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            # Last token (if numeric) is rsd_limit; rest is data
            tokens = input_params.strip().split()
            rsd_lim = 15.0
            data_str = input_params.strip()
            if tokens and tokens[-1].replace(".", "").replace("-", "").isdigit():
                rsd_lim = float(tokens[-1])
                data_str = " ".join(tokens[:-1])
            runs_raw = data_str.split(";")
            measurements = []
            for r in runs_raw:
                r = r.strip()
                if not r:
                    continue
                measurements.append([float(x) for x in r.split(",")])
            return self._run_base(measurements, None, rsd_lim)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'v1,v2,v3; v4,v5,v6; ... [rsd_limit]'")
