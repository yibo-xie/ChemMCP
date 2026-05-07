import logging
import math
from typing import List, Optional

import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class QcChartGenerator(BaseTool):
    """
    质控图（Shewhart Control Chart）生成与Westgard规则判断。
    
    生成 X-bar 图、R 图、S 图的控制限数据，并应用 Westgard 多规则判断系统。
    """
    __version__ = "0.1.0"
    name             = "QcChartGenerator"
    func_name        = "generate_qc_chart"
    description      = "Generate Shewhart QC control charts (X-bar, R, S) with Westgard multi-rule evaluation."
    implementation_description = "Calculates control limits (UCL/CL/LCL) for X-bar, R, and S charts, then applies Westgard rules (1-2s, 1-3s, 2-2s, R-4s, 4-1s, 10-x) to detect out-of-control conditions."
    oss_dependencies = [
        ("NumPy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories       = ["General"]
    tags             = ["QC", "Control Chart", "Quality Assurance", "Westgard Rules"]
    required_envs    = []

    code_input_sig   = [
        ("values", "list", "N/A", "List of QC measurement values."),
        ("chart_type", "str", "xbar", "Chart type: 'xbar', 'r', 's', or 'all'."),
        ("subgroup_size", "int", "1", "Subgroup size for R/S charts (default=1 for individual measurements)."),
        ("target_mean", "float_or_None", "None", "Target mean value; if None, calculated from data."),
        ("target_sd", "float_or_None", "None", "Target SD value; if None, calculated from data."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'v1,v2,v3,... [chart_type] [subgroup_size] [target_mean] [target_sd]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with chart data, control limits, Westgard rule violations, and status."),
    ]

    examples         = [
        {
            "code_input": {
                "values": [100.2, 100.5, 99.8, 100.1, 100.3, 99.9, 100.4, 100.0,
                           100.1, 99.7, 100.2, 100.6, 99.8, 100.0, 100.3, 99.95,
                           100.15, 102.5, 100.05, 99.85],
                "chart_type": "all",
                "subgroup_size": 1,
                "target_mean": None,
                "target_sd": None,
            },
            "text_input": {
                "input_params": "100.2,100.5,99.8,100.1,100.3,99.9,100.4,100.0,100.1,99.7,100.2,100.6,99.8,100.0,100.3,99.95,100.15,102.5,100.05,99.85 all 1",
            },
            "output": {
                "result": {
                    "xbar_chart": { "ucl": "...", "cl": "...", "lcl": "...", "points": [...] },
                    "westgard_violations": [...],
                    "status": "..."
                }
            },
        },
    ]

    # 常数表：A2, D3, D4, B3, B4 for different subgroup sizes
    _CONSTANTS = {
        2: {"A2": 1.880, "D3": 0,      "D4": 3.267, "B3": 0,      "B4": 3.267},
        3: {"A2": 1.023, "D3": 0,      "D4": 2.574, "B3": 0,      "B4": 2.568},
        4: {"A2": 0.729, "D3": 0,      "D4": 2.282, "B3": 0,      "B4": 2.266},
        5: {"A2": 0.577, "D3": 0,      "D4": 2.114, "B3": 0,      "B4": 2.089},
        6: {"A2": 0.483, "D3": 0,      "D4": 2.004, "B3": 0.030,   "B4": 1.970},
        7: {"A2": 0.419, "D3": 0.076,  "D4": 1.924, "B3": 0.118,   "B4": 1.882},
        8: {"A2": 0.373, "D3": 0.136,  "D4": 1.864, "B3": 0.185,   "B4": 1.815},
        9: {"A2": 0.337, "D3": 0.184,  "D4": 1.816, "B3": 0.239,   "B4": 1.761},
        10:{"A2": 0.308, "D3": 0.223,  "D4": 1.777, "B3": 0.284,   "B4": 1.716},
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        values: List[float],
        chart_type: str = "xbar",
        subgroup_size: int = 1,
        target_mean: Optional[float] = None,
        target_sd: Optional[float] = None,
    ) -> dict:
        """
        核心逻辑：生成质控图数据 + Westgard规则判断
        """
        if len(values) < 2:
            raise ChemMCPError("At least 2 data points are required.")

        data = np.array(values, dtype=float)
        n_total = len(data)
        ct = chart_type.lower()

        # 使用目标值或计算统计量
        mean_val = target_mean if target_mean is not None else float(np.mean(data))
        sd_val = target_sd if target_sd is not None else float(np.std(data, ddof=1))

        result = {}
        all_points_with_status = []

        if ct in ("xbar", "all"):
            ucl_x = round(mean_val + 3 * sd_val, 6)
            lcl_x = round(mean_val - 3 * sd_val, 6)
            ucl_2s = round(mean_val + 2 * sd_val, 6)
            lcl_2s = round(mean_val - 2 * sd_val, 6)

            points_x = []
            for i, v in enumerate(data):
                z_score = (v - mean_val) / sd_val if sd_val > 0 else 0
                status = "in_control"
                if abs(z_score) > 3:
                    status = "out_of_control_3s"
                elif abs(z_score) > 2:
                    status = "warning_2s"
                points_x.append({
                    "index": int(i),
                    "value": round(float(v), 6),
                    "z_score": round(z_score, 4),
                    "status": status,
                })
                all_points_with_status.append(status)

            result["xbar_chart"] = {
                "ucl": ucl_x,
                "cl": round(mean_val, 6),
                "lcl": lcl_x,
                "ucl_2s": ucl_2s,
                "lcl_2s": lcl_2s,
                "mean": round(mean_val, 6),
                "sd": round(sd_val, 6),
                "n_points": n_total,
                "points": points_x,
            }

        if ct in ("r", "s", "all") and subgroup_size > 1:
            # 分组数据计算 R/S 图
            n_groups = n_total // subgroup_size
            subgroups = data[:n_groups * subgroup_size].reshape(n_groups, subgroup_size)
            group_means = np.mean(subgroups, axis=1)
            group_ranges = np.ptp(subgroups, axis=1)
            group_sds = np.std(subgroups, ddof=1, axis=1)

            r_bar = float(np.mean(group_ranges))
            s_bar = float(np.mean(group_sds))
            c = self._CONSTANTS.get(subgroup_size, self._CONSTANTS[5])

            if ct in ("r", "all"):
                result["r_chart"] = {
                    "ucl": round(r_bar * c["D4"], 6),
                    "cl": round(r_bar, 6),
                    "lcl": round(r_bar * c["D3"], 6),
                    "points": [round(float(v), 6) for v in group_ranges],
                }

            if ct in ("s", "all"):
                result["s_chart"] = {
                    "ucl": round(s_bar * c["B4"], 6),
                    "cl": round(s_bar, 6),
                    "lcl": round(s_bar * max(c["B3"], 0), 6),
                    "points": [round(float(v), 6) for v in group_sds],
                }
        elif ct in ("r", "s"):
            result["note"] = "R/S charts require subgroup_size > 1."

        # ---- Westgard 规则判断 ----
        violations = self._apply_westgard_rules(data, mean_val, sd_val)
        result["westgard_violations"] = violations
        result["in_control"] = len(violations) == 0
        result["status"] = "In Control" if len(violations) == 0 else f"Out of Control ({len(violations)} violation(s))"

        logger.info(f"QC chart generated: n={n_total}, mean={mean_val:.4f}, sd={sd_val:.4f}, violations={len(violations)}")
        return result

    def _apply_westgard_rules(self, data: np.ndarray, mean: float, sd: float) -> list:
        """应用 Westgard 多规则"""
        violations = []
        n = len(data)
        if sd == 0 or n < 2:
            return violations

        z_scores = [(float(d) - mean) / sd for d in data]

        # Rule 1_2s: 一个点超出 ±2s（警告）
        warning_indices = [i for i, z in enumerate(z_scores) if abs(z) > 2]

        # Rule 1_3s: 一个点超出 ±3s → 失控
        for i, z in enumerate(z_scores):
            if abs(z) > 3:
                violations.append({
                    "rule": "1_3s",
                    "point_index": int(i),
                    "value": round(float(data[i]), 6),
                    "z_score": round(z, 4),
                    "description": f"Point {i+1} exceeds ±3s (z={z:.2f})",
                })

        # Rule 2_2s: 连续2个点同侧超过 +2s 或 -2s
        for i in range(n - 1):
            if z_scores[i] > 2 and z_scores[i+1] > 2:
                violations.append({"rule": "2_2s", "range": f"{i+1}-{i+2}", "side": "+2s",
                                   "description": f"Points {i+1}-{i+2} both exceed +2s"})
            if z_scores[i] < -2 and z_scores[i+1] < -2:
                violations.append({"rule": "2_2s", "range": f"{i+1}-{i+2}", "side": "-2s",
                                   "description": f"Points {i+1}-{i+2} both exceed -2s"})

        # Rule R_4s: 连续2个点相差 > 4s（一个高一个低）
        for i in range(n - 1):
            if abs(z_scores[i] - z_scores[i+1]) >= 4:
                violations.append({"rule": "R_4s", "range": f"{i+1}-{i+2}",
                                   "description": f"Points {i+1}-{i+2} span > 4s ({z_scores[i]:.1f} to {z_scores[i+1]:.1f})"})

        # Rule 4_1s: 连续4个点同侧超过 +1s 或 -1s
        for i in range(n - 3):
            if all(z > 1 for z in z_scores[i:i+4]):
                violations.append({"rule": "4_1s", "range": f"{i+1}-{i+4}", "side": "+1s",
                                   "description": f"Points {i+1}-{i+4} all exceed +1s"})
            if all(z < -1 for z in z_scores[i:i+4]):
                violations.append({"rule": "4_1s", "range": f"{i+1}-{i+4}", "side": "-1s",
                                   "description": f"Points {i+1}-{i+4} all below -1s"})

        # Rule 10_x: 连续10个点在均值同一侧
        for i in range(n - 9):
            if all(z > 0 for z in z_scores[i:i+10]):
                violations.append({"rule": "10_x", "range": f"{i+1}-{i+10}", "side": "above_mean",
                                   "description": f"Points {i+1}-{i+10} all above mean"})
            if all(z < 0 for z in z_scores[i:i+10]):
                violations.append({"rule": "10_x", "range": f"{i+1}-{i+10}", "side": "below_mean",
                                   "description": f"Points {i+1}-{i+10} all below mean"})

        return violations

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            vals = [float(x) for x in parts[0].split(",")]
            ct = parts[1] if len(parts) > 1 else "xbar"
            sg = int(parts[2]) if len(parts) > 2 else 1
            tm = float(parts[3]) if len(parts) > 3 and parts[3].lower() != "none" else None
            tsd = float(parts[4]) if len(parts) > 4 and parts[4].lower() != "none" else None
            return self._run_base(vals, ct, sg, tm, tsd)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'v1,v2,... [chart_type] [sg_size] [target_mean] [target_sd]'")
