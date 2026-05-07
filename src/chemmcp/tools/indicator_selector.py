import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IndicatorSelector(BaseTool):
    """
    酸碱指示剂选择工具：基于 pKa 匹配和当量点 pH 选择最合适的指示剂。
    
    内置常用酸碱指示剂数据库，支持自动匹配和评分。
    """
    __version__ = "0.1.0"
    name             = "IndicatorSelector"
    func_name        = "select_indicator"
    description      = "Select appropriate acid-base titration indicator based on expected pH at equivalence point and pKa matching."
    implementation_description = "Uses a built-in database of common acid-base indicators with pH transition ranges. Scores and ranks indicators by how well their transition range covers the expected equivalence point pH."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Indicator", "Acid-Base Titration", "Analytical Chemistry", "pKa"]
    required_envs    = []

    # 内置指示剂数据库
    INDICATOR_DB = [
        {"name": "Methyl Violet (1st)",     "range_low": 0.0,  "range_high": 1.6,  "acid_color": "Yellow",   "base_color": "Blue",      "pKa_approx": None, "pK_type": "multi"},
        {"name": "Thymol Blue (1st)",        "range_low": 1.2,  "range_high": 2.8,  "acid_color": "Red",     "base_color": "Yellow",    "pKa_approx": 1.5,  "pK_type": "1st"},
        {"name": "Methyl Orange",            "range_low": 3.1,  "range_high": 4.4,  "acid_color": "Red",     "base_color": "Yellow",    "pKa_approx": 3.4,  "pK_type": "sulfonphthalein"},
        {"name": "Bromophenol Blue",         "range_low": 3.0,  "range_high": 4.6,  "acid_color": "Yellow",  "base_color": "Purple",   "pKa_approx": 3.85, "pK_type": "sulfonphthalein"},
        {"name": "Bromocresol Green",        "range_low": 3.8,  "range_high": 5.4,  "acid_color": "Yellow",  "base_color": "Blue",      "pKa_approx": 4.68, "pK_type": "sulfonphthalein"},
        {"name": "Methyl Red",               "range_low": 4.4,  "range_high": 6.2,  "acid_color": "Red",     "base_color": "Yellow",    "pKa_approx": 5.0,  "pK_type": "azo"},
        {"name": "Bromocresol Purple",       "range_low": 5.2,  "range_high": 6.8,  "acid_color": "Yellow",  "base_color": "Purple",   "pKa_approx": 6.3,  "pK_type": "sulfonphthalein"},
        {"name": "Chlorophenol Red",         "range_low": 5.0,  "range_high": 6.6,  "acid_color": "Yellow",  "base_color": "Red",       "pKa_approx": 5.8,  "pK_type": "sulfonphthalein"},
        {"name": "Bromothymol Blue",         "range_low": 6.0,  "range_high": 7.6,  "acid_color": "Yellow",  "base_color": "Blue",      "pKa_approx": 7.1,  "pK_type": "sulfonphthalein"},
        {"name": "Phenol Red",               "range_low": 6.4,  "range_high": 8.0,  "acid_color": "Yellow",  "base_color": "Red",       "pKa_approx": 7.54, "pK_type": "sulfonphthalein"},
        {"name": "Neutral Red",              "range_low": 6.8,  "range_high": 8.0,  "acid_color": "Red",     "base_color": "Yellow",    "pKa_approx": 7.4,  "pK_type": "phenazine"},
        {"name": "Cresol Purple (1st)",      "range_low": 1.2,  "range_high": 2.8,  "acid_color": "Red",     "base_color": "Yellow",    "pKa_approx": 1.9,  "pK_type": "sulfonphthalein"},
        {"name": "Cresol Purple (2nd)",      "range_low": 7.2,  "range_high": 8.8,  "acid_color": "Yellow",  "base_color": "Purple",   "pKa_approx": 8.0,  "pK_type": "sulfonphthalein"},
        {"name": "Meta-Cresol Purple",       "range_low": 7.6,  "range_high": 9.2,  "acid_color": "Yellow",  "base_color": "Purple",   "pKa_approx": 8.32, "pK_type": "sulfonphthalein"},
        {"name": "Phenolphthalein",          "range_low": 8.2,  "range_high": 10.0, "acid_color": "Colorless","base_color": "Pink/Red", "pKa_approx": 9.3,  "pK_type": "phthalein"},
        {"name": "Thymolphthalein",          "range_low": 9.3,  "range_high": 10.5, "acid_color": "Colorless","base_color": "Blue",      "pKa_approx": 9.7,  "pK_type": "phthalein"},
        {"name": "Alizarin Yellow R",        "range_low": 10.1, "range_high": 12.0, "acid_color": "Yellow",  "base_color": "Red",       "pKa_approx": 11.0, "pK_type": "azo"},
        {"name": "Malachite Green (2nd)",    "range_low": 11.4, "range_high": 13.6, "acid_color": "Yellow",  "base_color": "Colorless", "pKa_approx": 12.5, "pK_type": "triarylmethane"},
        {"name": "Tropeolin OOO",            "range_low": 11.0, "range_high": 13.0, "acid_color": "Brown",   "base_color": "Orange",    "pKa_approx": 12.0, "pK_type": "azo"},
    ]

    code_input_sig   = [
        ("equivalence_ph", "float", "N/A", "Expected pH at the equivalence point."),
        ("titration_type", "str", "general", "Titration type: 'strong_acid_strong_base', 'weak_acid_strong_base', 'strong_acid_weak_base', or 'general'."),
        ("tolerance", "float", "1.0", "Maximum acceptable deviation of indicator midpoint from equivalence pH (pH units)."),
        ("top_n", "int", "5", "Number of top recommendations to return."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'equivalence_ph [titration_type] [tolerance] [top_n]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with best indicator, ranked alternatives, matching details, and selection rationale."),
    ]

    examples         = [
        {
            "code_input": {
                "equivalence_ph": 8.72,
                "titration_type": "weak_acid_strong_base",
                "tolerance": 1.0,
                "top_n": 5,
            },
            "text_input": {
                "input_params": "8.72 weak_acid_strong_base 1.0 5",
            },
            "output": {
                "result": {
                    "best_indicator": "...",
                    "ranked_indicators": [...],
                    "selection_rationale": "..."
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
        equivalence_ph: float,
        titration_type: str = "general",
        tolerance: float = 1.0,
        top_n: int = 5,
    ) -> dict:
        """核心逻辑：指示剂选择与评分"""
        if equivalence_ph < 0 or equivalence_ph > 14:
            raise ChemMCPError("Equivalence pH must be between 0 and 14.")

        scored = []
        for ind in self.INDICATOR_DB:
            mid = (ind["range_low"] + ind["range_high"]) / 2.0
            span = ind["range_high"] - ind["range_low"]

            # 判断是否在变色范围内
            in_range = ind["range_low"] <= equivalence_ph <= ind["range_high"]
            distance = abs(mid - equivalence_ph)

            # 评分（0-100）
            if in_range:
                score = 100.0 - (distance / (span / 2.0) * 30) if span > 0 else 100.0
            elif distance <= tolerance:
                score = 70.0 - (distance / tolerance * 30)
            else:
                score = max(0, 40.0 - (distance - tolerance) * 10)

            # 窄范围加分（更敏锐的终点检测）
            if span < 1.5:
                score += 5

            scored.append({
                **ind,
                "midpoint_ph": round(mid, 3),
                "distance_from_eq_ph": round(distance, 3),
                "in_transition_range": in_range,
                "score": round(min(score, 100), 2),
            })

        # 按分数排序
        scored.sort(key=lambda x: x["score"], reverse=True)
        ranked = scored[:top_n]
        best = ranked[0]

        # 选择理由
        if best["in_transition_range"]:
            reason = f"'{best['name']}' has transition range {best['range_low']}-{best['range_high']} which contains the expected pH {equivalence_ph:.2f}. Color change: {best['acid_color']} → {best['base_color']}."
        else:
            reason = f"'{best['name']}' is closest available (midpoint={best['midpoint_ph']}, ΔpH={best['distance_from_eq_ph']}). Consider using potentiometric detection for better accuracy."

        result = {
            "equivalence_ph": round(equivalence_ph, 4),
            "titration_type": titration_type,
            "tolerance_pH_units": tolerance,
            "best_indicator": {
                "name": best["name"],
                "transition_range": f"{best['range_low']}-{best['range_high']}",
                "midpoint_ph": best["midpoint_ph"],
                "color_change": f"{best['acid_color']} → {best['base_color']}",
                "score": best["score"],
                "in_range": best["in_transition_range"],
            },
            "ranked_alternatives": [
                {
                    "rank": i + 1,
                    "name": r["name"],
                    "range": f"{r['range_low']}-{r['range_high']}",
                    "color_change": f"{r['acid_color']} → {r['base_color']}",
                    "score": r["score"],
                    "in_range": r["in_transition_range"],
                } for i, r in enumerate(ranked)
            ],
            "all_in_range": [s["name"] for s in scored if s["in_transition_range"]],
            "selection_rationale": reason,
        }

        logger.info(f"Indicator selected: '{best['name']}' for eq_pH={equivalence_ph:.2f}, score={best['score']}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            eq_ph = float(parts[0])
            ttype = parts[1] if len(parts) > 1 else "general"
            tol = float(parts[2]) if len(parts) > 2 else 1.0
            tn = int(parts[3]) if len(parts) > 3 else 5
            return self._run_base(eq_ph, ttype, tol, tn)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'equivalence_ph [titration_type] [tolerance] [top_n]'")
