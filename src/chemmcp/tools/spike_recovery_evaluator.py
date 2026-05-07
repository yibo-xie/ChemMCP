import logging
import math
from typing import List, Optional

import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SpikeRecoveryEvaluator(BaseTool):
    """
    加标回收率统计评估工具。
    
    用于分析方法验证中的准确度评估，计算加标回收率和相对标准偏差。
    """
    __version__ = "0.1.0"
    name             = "SpikeRecoveryEvaluator"
    func_name        = "evaluate_spike_recovery"
    description      = "Calculate spike recovery rates and assess accuracy for analytical method validation."
    implementation_description = "Computes individual recovery percentages (C_spiked - C_unspiked) / C_added × 100%, mean recovery, RSD%, and judges against acceptance criteria (typically 80-120%)."
    oss_dependencies = [
        ("NumPy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories       = ["General"]
    tags             = ["Spike Recovery", "Method Validation", "Accuracy", "Analytical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("unspiked_values", "list", "N/A", "Measured concentrations of unspiked samples (same units as spiked_conc)."),
        ("spiked_values", "list", "N/A", "Measured concentrations of spiked samples."),
        ("spiked_concentration", "float", "N/A", "The concentration of spike added to each sample."),
        ("acceptance_low", "float", "80.0", "Lower acceptance limit for recovery (%)."),
        ("acceptance_high", "float", "120.0", "Upper acceptance limit for recovery (%)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'u1,u2,u3 s1,s2,s3 spiked_conc [acc_low] [acc_high]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with recovery statistics, acceptance judgment, and detailed results."),
    ]

    examples         = [
        {
            "code_input": {
                "unspiked_values": [0.52, 0.48, 0.51, 0.49, 0.50, 0.53],
                "spiked_values": [1.48, 1.45, 1.50, 1.47, 1.49, 1.52],
                "spiked_concentration": 1.0,
                "acceptance_low": 80.0,
                "acceptance_high": 120.0,
            },
            "text_input": {
                "input_params": "0.52,0.48,0.51,0.49,0.50,0.53 1.48,1.45,1.50,1.47,1.49,1.52 1.0 80 120",
            },
            "output": {
                "result": {
                    "mean_recovery_pct": "...",
                    "rsd_pct": "...",
                    "accepted": True,
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
        unspiked_values: List[float],
        spiked_values: List[float],
        spiked_concentration: float,
        acceptance_low: float = 80.0,
        acceptance_high: float = 120.0,
    ) -> dict:
        """核心逻辑：加标回收率评估"""
        if len(unspiked_values) != len(spiked_values):
            raise ChemMCPError("Unspiked and spiked value lists must have the same length.")
        if len(unspiked_values) < 2:
            raise ChemMCPError("At least 2 pairs of measurements are required.")
        if spiked_concentration <= 0:
            raise ChemMCPError("Spiked concentration must be positive.")

        u = np.array(unspiked_values, dtype=float)
        s = np.array(spiked_values, dtype=float)
        c_added = float(spiked_concentration)

        # 计算每个样品的回收率
        recoveries = ((s - u) / c_added) * 100.0

        # 统计量
        mean_rec = float(np.mean(recoveries))
        sd_rec = float(np.std(recoveries, ddof=1))
        rsd_rec = (sd_rec / mean_rec * 100.0) if mean_rec != 0 else float('inf')

        n = len(recoveries)
        se = sd_rec / math.sqrt(n)
        t_crit = 2.571  # approximate t(0.975, df=5), will compute properly
        from scipy import stats as sp_stats
        try:
            t_crit = float(sp_stats.t.ppf(0.975, n - 1))
        except Exception:
            pass

        ci_low = round(mean_rec - t_crit * se, 4)
        ci_high = round(mean_rec + t_crit * se, 4)

        # 判断是否合格
        accepted = acceptance_low <= mean_rec <= acceptance_high

        # 详细结果
        details = []
        for i in range(n):
            r = round(float(recoveries[i]), 4)
            in_range = acceptance_low <= r <= acceptance_high
            details.append({
                "sample": int(i + 1),
                "unspiked": round(float(u[i]), 6),
                "spiked": round(float(s[i]), 6),
                "recovery_pct": r,
                "accepted": in_range,
            })

        result = {
            "n_samples": n,
            "spiked_concentration": c_added,
            "recovery_statistics": {
                "mean_recovery_pct": round(mean_rec, 4),
                "sd_recovery_pct": round(sd_rec, 4),
                "rsd_pct": round(rsd_rec, 4),
                "min_recovery_pct": round(float(np.min(recoveries)), 4),
                "max_recovery_pct": round(float(np.max(recoveries)), 4),
                "median_recovery_pct": round(float(np.median(recoveries)), 4),
                "ci_95_pct": [ci_low, ci_high],
            },
            "acceptance_criteria": {
                "lower_limit": acceptance_low,
                "upper_limit": acceptance_high,
                "accepted": accepted,
            },
            "sample_details": details,
            "conclusion": (
                f"Mean recovery = {mean_rec:.2f}%, RSD = {rsd_rec:.2f}%. "
                f"{'PASS' if accepted else 'FAIL'}: recovery {'within' if accepted else 'outside'} [{acceptance_low}%, {acceptance_high}%]."
            ),
        }

        logger.info(f"Spike recovery evaluation: mean={mean_rec:.2f}%, RSD={rsd_rec:.2f}%, accepted={accepted}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            u_vals = [float(x) for x in parts[0].split(",")]
            s_vals = [float(x) for x in parts[1].split(",")]
            sp_conc = float(parts[2])
            acc_lo = float(parts[3]) if len(parts) > 3 else 80.0
            acc_hi = float(parts[4]) if len(parts) > 4 else 120.0
            return self._run_base(u_vals, s_vals, sp_conc, acc_lo, acc_hi)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'u1,u2,u3 s1,s2,s3 spiked_conc [low] [high]'")
