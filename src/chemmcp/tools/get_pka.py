import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 内置 pKa 数据库：常见酸的 pKa 值
# 数据来源: CRC Handbook of Chemistry & Physics, Wikipedia (标准化学数据)
_PKA_DATABASE = {
    # 一元酸 (Monoprotic acids)
    "hydrochloric acid":      {"formula": "HCl",     "pKa": -7.0,   "type": "monoprotic", "category": "strong"},
    "hcl":                    {"formula": "HCl",     "pKa": -7.0,   "type": "monoprotic", "category": "strong"},
    "nitric acid":            {"formula": "HNO3",    "pKa": -1.4,   "type": "monoprotic", "category": "strong"},
    "hno3":                   {"formula": "HNO3",    "pKa": -1.4,   "type": "monoprotic", "category": "strong"},
    "sulfuric acid (1st)":    {"formula": "H2SO4",   "pKa": -3.0,   "type": "diprotic_1st", "category": "strong"},
    "perchloric acid":        {"formula": "HClO4",   "pKa": -10.0,  "type": "monoprotic", "category": "strong"},
    "hclo4":                  {"formula": "HClO4",   "pKa": -10.0,  "type": "monoprotic", "category": "strong"},
    "acetic acid":            {"formula": "CH3COOH","pKa": 4.76,   "type": "monoprotic", "category": "weak"},
    "ch3cooh":                {"formula": "CH3COOH","pKa": 4.76,   "type": "monoprotic", "category": "weak"},
    "formic acid":            {"formula": "HCOOH",   "pKa": 3.75,   "type": "monoprotic", "category": "weak"},
    "hcooh":                  {"formula": "HCOOH",   "pKa": 3.75,   "type": "monoprotic", "category": "weak"},
    "benzoic acid":           {"formula": "C6H5COOH","pKa": 4.20,  "type": "monoprotic", "category": "weak"},
    "hydrofluoric acid":      {"formula": "HF",      "pKa": 3.17,   "type": "monoprotic", "category": "weak"},
    "hf":                     {"formula": "HF",      "pKa": 3.17,   "type": "monoprotic", "category": "weak"},
    "hypochlorous acid":      {"formula": "HClO",    "pKa": 7.53,   "type": "monoprotic", "category": "weak"},
    "hclo":                   {"formula": "HClO",    "pKa": 7.53,   "type": "monoprotic", "category": "weak"},
    "cyanic acid":            {"formula": "HCN",     "pKa": 9.31,   "type": "monoprotic", "category": "weak"},
    "hcn":                    {"formula": "HCN",     "pKa": 9.31,   "type": "monoprotic", "category": "weak"},
    "boric acid":             {"formula": "H3BO3",   "pKa": 9.24,   "type": "monoprotic", "category": "weak"},  # 实际为路易斯酸，表观 pKa
    "phenol":                 {"formula": "C6H5OH", "pKa": 9.99,   "type": "monoprotic", "category": "weak"},
    "bicarbonate":            {"formula": "HCO3-",   "pKa": 10.33,  "type": "monoprotic", "category": "weak"},
    "hydrogen carbonate":     {"formula": "HCO3-",   "pKa": 10.33,  "type": "monoprotic", "category": "weak"},
    "ammonium":               {"formula": "NH4+",    "pKa": 9.25,   "type": "monoprotic", "category": "weak"},

    # 二元酸 (Diprotic acids)
    "carbonic acid":          {"formula": "H2CO3",   "pKa": [6.35, 10.33], "type": "diprotic", "category": "weak"},
    "h2co3":                  {"formula": "H2CO3",   "pKa": [6.35, 10.33], "type": "diprotic", "category": "weak"},
    "oxalic acid":            {"formula": "H2C2O4",  "pKa": [1.25, 4.27],  "type": "diprotic", "category": "weak"},
    "h2c2o4":                 {"formula": "H2C2O4",  "pKa": [1.25, 4.27],  "type": "diprotic", "category": "weak"},
    "sulfuric acid":          {"formula": "H2SO4",   "pKa": [-3.0, 1.99],  "type": "diprotic", "category": "mixed"},
    "h2so4":                  {"formula": "H2SO4",   "pKa": [-3.0, 1.99],  "type": "diprotic", "category": "mixed"},
    "hydrosulfuric acid":     {"formula": "H2S",     "pKa": [7.04, 19.0],  "type": "diprotic", "category": "weak"},
    "h2s":                    {"formula": "H2S",     "pKa": [7.04, 19.0],  "type": "diprotic", "category": "weak"},
    "chromic acid":           {"formula": "H2CrO4",  "pKa": [0.74, 6.49],  "type": "diprotic", "category": "weak"},
    "hydrazoic acid":         {"formula": "HN3",     "pKa": [4.65, 14.0],  "type": "diprotic", "category": "weak"},

    # 三元酸 (Triprotic acids)
    "phosphoric acid":        {"formula": "H3PO4",   "pKa": [2.15, 7.20, 12.35], "type": "triprotic", "category": "weak"},
    "h3po4":                  {"formula": "H3PO4",   "pKa": [2.15, 7.20, 12.35], "type": "triprotic", "category": "weak"},
    "citric acid":            {"formula": "C6H8O7",  "pKa": [3.13, 4.76, 6.40],  "type": "triprotic", "category": "weak"},
    "arsenic acid":           {"formula": "H3AsO4",  "pKa": [2.26, 6.76, 11.29], "type": "triprotic", "category": "weak"},
    "boric acid (full)":      {"formula": "H3BO3",   "pKa": [9.24, 12.74, 13.80], "type": "triprotic", "category": "weak"},
}


@ChemMCPManager.register_tool
class GetPka(BaseTool):
    """
    查询常见酸的 pKa 值。
    包含内置数据库，支持一元酸、二元酸、三元酸。
    """
    __version__ = "0.1.0"
    name = "GetPka"
    func_name = "get_pka"
    description = "Look up pKa values for common acids from a built-in database."
    implementation_description = "Uses a built-in database of ~30 common acids with pKa values from CRC Handbook / standard references."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Acid-Base", "pKa", "Database", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("acid_name", "str", "N/A", "Name of the acid to look up (case-insensitive, supports common names and formulas)."),
    ]

    text_input_sig = [
        ("acid_name", "str", "N/A", "Name of the acid to look up."),
    ]

    output_sig = [
        ("acid_name", "str", "The matched acid name."),
        ("formula", "str", "Molecular formula of the acid."),
        ("pKa", "float or list", "pKa value(s). Single float for monoprotic, list for polyprotic."),
        ("type", "str", "Acid type: monoprotic, diprotic, or triprotic."),
        ("category", "str", "Strength category: strong, weak, or mixed."),
        ("n_protonic", "int", "Number of dissociable protons."),
    ]

    examples = [
        {
            "code_input": {
                "acid_name": "acetic acid"
            },
            "text_input": {
                "acid_name": "acetic acid"
            },
            "output": {
                "acid_name": "acetic acid",
                "formula": "CH3COOH",
                "pKa": 4.76,
                "type": "monoprotic",
                "category": "weak",
                "n_protonic": 1,
            }
        },
        {
            "code_input": {
                "acid_name": "phosphoric acid"
            },
            "text_input": {
                "acid_name": "phosphoric acid"
            },
            "output": {
                "acid_name": "phosphoric acid",
                "formula": "H3PO4",
                "pKa": [2.15, 7.20, 12.35],
                "type": "triprotic",
                "category": "weak",
                "n_protonic": 3,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """加载 pKa 数据库"""
        self.database = _PKA_DATABASE

    def _run_base(self, acid_name: str) -> dict:
        """查询 pKa"""
        if not acid_name or not acid_name.strip():
            raise ChemMCPError("Acid name cannot be empty.")

        key = acid_name.strip().lower()

        # 精确匹配
        if key in self.database:
            entry = self.database[key]
            return self._format_result(key, entry)

        # 模糊匹配（包含搜索）
        matches = [k for k in self.database if key in k or k in key]
        if len(matches) == 1:
            entry = self.database[matches[0]]
            return self._format_result(matches[0], entry)
        elif len(matches) > 1:
            raise ChemMCPError(
                f"Multiple matches found for '{acid_name}': {matches}. Please be more specific."
            )

        # 未找到
        available = sorted(set(
            k for k, v in self.database.items()
            if isinstance(v.get("pKa"), (int, float)) or True
        ))
        raise ChemMCPError(
            f"Acid '{acid_name}' not found in database. "
            f"Available acids: {available[:20]}{'...' if len(available) > 20 else ''}"
        )

    def _format_result(self, name: str, entry: dict) -> dict:
        """格式化输出结果"""
        pKa_val = entry["pKa"]
        acid_type = entry["type"]

        n_protonic = 1
        if "tri" in acid_type:
            n_protonic = 3
        elif "di" in acid_type:
            n_protonic = 2

        return {
            "acid_name": name,
            "formula": entry["formula"],
            "pKa": pKa_val,
            "type": acid_type,
            "category": entry.get("category", "unknown"),
            "n_protonic": n_protonic,
        }

    def _run_text(self, acid_name: str) -> dict:
        """文本接口直接调用核心逻辑"""
        return self._run_base(acid_name)

    def list_available_acids(self) -> List[str]:
        """列出所有可用的酸"""
        return sorted(self.database.keys())
