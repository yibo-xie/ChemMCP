import logging
from typing import List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 内置 pKb 数据库：常见碱的 pKb 值
# 数据来源: CRC Handbook of Chemistry & Physics, 标准化学数据
_PKB_DATABASE = {
    # 强碱
    "sodium hydroxide":       {"formula": "NaOH",    "pKb": None,     "type": "strong",   "category": "strong"},
    "naoh":                   {"formula": "NaOH",    "pKb": None,     "type": "strong",   "category": "strong"},
    "potassium hydroxide":    {"formula": "KOH",     "pKb": None,     "type": "strong",   "category": "strong"},
    "koh":                    {"formula": "KOH",     "pKb": None,     "type": "strong",   "category": "strong"},
    "calcium hydroxide":      {"formula": "Ca(OH)2", "pKb": None,     "type": "strong",   "category": "strong"},
    "ca(oh)2":                {"formula": "Ca(OH)2", "pKb": None,     "type": "strong",   "category": "strong"},
    "barium hydroxide":       {"formula": "Ba(OH)2", "pKb": None,     "type": "strong",   "category": "strong"},

    # 弱碱（分子碱）
    "ammonia":                {"formula": "NH3",     "pKb": 4.75,     "type": "monoprotic", "category": "weak"},
    "nh3":                    {"formula": "NH3",     "pKb": 4.75,     "type": "monoprotic", "category": "weak"},
    "methylamine":            {"formula": "CH3NH2",  "pKb": 3.36,     "type": "monoprotic", "category": "weak"},
    "ethylamine":             {"formula": "C2H5NH2", "pKb": 3.25,     "type": "monoprotic", "category": "weak"},
    "dimethylamine":          {"formula": "(CH3)2NH","pKb": 3.27,     "type": "monoprotic", "category": "weak"},
    "aniline":                {"formula": "C6H5NH2","pKb": 9.38,      "type": "monoprotic", "category": "weak"},
    "pyridine":               {"formula": "C5H5N",   "pKb": 8.75,     "type": "monoprotic", "category": "weak"},
    "hydrazine (1st)":        {"formula": "N2H4",    "pKb": [6.07, 15.0], "type": "diprotic", "category": "weak"},
    "ethylenediamine":        {"formula": "C2H8N2",  "pKb": [3.29, 6.44], "type": "diprotic", "category": "weak"},

    # 共轭碱（弱酸的共轭碱）
    "acetate":                {"formula": "CH3COO-","pKb": 9.24,      "type": "monoprotic", "category": "conjugate_base"},
    "sodium acetate":         {"formula": "CH3COONa","pKb": 9.24,     "type": "monoprotic", "category": "conjugate_base"},
    "fluoride":               {"formula": "F-",      "pKb": 10.83,    "type": "monoprotic", "category": "conjugate_base"},
    "carbonate":              {"formula": "CO3^2-","pKb": [3.67, 7.65], "type": "diprotic", "category": "conjugate_base"},
    "sodium carbonate":       {"formula": "Na2CO3", "pKb": [3.67, 7.65], "type": "diprotic", "category": "conjugate_base"},
    "bicarbonate":            {"formula": "HCO3-",  "pKb": 7.65,      "type": "monoprotic", "category": "conjugate_base"},
    "sodium bicarbonate":     {"formula": "NaHCO3", "pKb": 7.65,      "type": "monoprotic", "category": "conjugate_base"},
    "phosphate":              {"formula": "PO4^3-","pKb": [1.74, 6.80, 11.65], "type": "triprotic", "category": "conjugate_base"},
    "hydrogen phosphate":     {"formula": "HPO4^2-","pKb": [1.74, 6.80],"type": "diprotic", "category": "conjugate_base"},
}


@ChemMCPManager.register_tool
class GetPkb(BaseTool):
    """
    查询常见碱的 pKb 值。
    包含内置数据库，支持强碱、分子弱碱、共轭碱。
    """
    __version__ = "0.1.0"
    name = "GetPkb"
    func_name = "get_pkb"
    description = "Look up pKb values for common bases from a built-in database."
    implementation_description = "Uses a built-in database of ~20 common bases with pKb values from CRC Handbook / standard references."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Acid-Base", "pKb", "Database", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("base_name", "str", "N/A", "Name of the base to look up (case-insensitive)."),
    ]

    text_input_sig = [
        ("base_name", "str", "N/A", "Name of the base to look up."),
    ]

    output_sig = [
        ("base_name", "str", "The matched base name."),
        ("formula", "str", "Formula of the base."),
        ("pKb", "float or list", "pKb value(s). Single float for monoprotic, list for polyprotic bases."),
        ("type", "str", "Base type: strong, monoprotic, diprotic, triprotic."),
        ("category", "str", "Category: strong, weak, conjugate_base."),
        ("n_basic", "int", "Number of basic sites (None for strong bases)."),
    ]

    examples = [
        {
            "code_input": {
                "base_name": "ammonia"
            },
            "text_input": {
                "base_name": "ammonia"
            },
            "output": {
                "base_name": "ammonia",
                "formula": "NH3",
                "pKb": 4.75,
                "type": "monoprotic",
                "category": "weak",
                "n_basic": 1,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """加载 pKb 数据库"""
        self.database = _PKB_DATABASE

    def _run_base(self, base_name: str) -> dict:
        """查询 pKb"""
        if not base_name or not base_name.strip():
            raise ChemMCPError("Base name cannot be empty.")

        key = base_name.strip().lower()

        # 精确匹配
        if key in self.database:
            entry = self.database[key]
            return self._format_result(key, entry)

        # 模糊匹配
        matches = [k for k in self.database if key in k or k in key]
        if len(matches) == 1:
            entry = self.database[matches[0]]
            return self._format_result(matches[0], entry)
        elif len(matches) > 1:
            raise ChemMCPError(
                f"Multiple matches found for '{base_name}': {matches}. Please be more specific."
            )

        available = sorted(self.database.keys())
        raise ChemMCPError(
            f"Base '{base_name}' not found in database. "
            f"Available bases: {available[:20]}{'...' if len(available) > 20 else ''}"
        )

    def _format_result(self, name: str, entry: dict) -> dict:
        pKb_val = entry["pKb"]
        b_type = entry["type"]

        n_basic = None
        if b_type == "monoprotic":
            n_basic = 1
        elif b_type == "diprotic":
            n_basic = 2
        elif b_type == "triprotic":
            n_basic = 3

        return {
            "base_name": name,
            "formula": entry["formula"],
            "pKb": pKb_val,
            "type": b_type,
            "category": entry.get("category", "unknown"),
            "n_basic": n_basic,
        }

    def _run_text(self, base_name: str) -> dict:
        return self._run_base(base_name)

    def list_available_bases(self) -> List[str]:
        return sorted(self.database.keys())
