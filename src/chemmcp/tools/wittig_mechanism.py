import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Wittig 反应机理知识库
WITTIG_MECHANISM_DATA = {
    "general": {
        "name": "Wittig Reaction",
        "overview": "Wittig反应是醛或酮与磷叶立德（ylide）反应生成烯烃和氧化三苯基膦的反应。该反应是合成烯烃的重要方法，具有高度的区域和立体选择性。",
        "steps": [
            {
                "step": 1,
                "name": "成盐（Nucleophilic Addition / Betaine Formation）",
                "description": "磷叶立德的碳负离子进攻羰基碳，形成偶极中间体（betaine）。这是速率决定步骤。",
                "electron_flow": "叶立德碳负离子 → 羰基碳（亲核加成）",
                "key_features": ["碳负离子亲核进攻", "C=O π键断裂", "形成四配位氧负离子"]
            },
            {
                "step": 2,
                "name": "Betaine 中间体",
                "description": "形成的 betaine 是一个两性离子中间体，同时带有正电荷（磷上）和负电荷（氧上）。",
                "electron_flow": "无净电子流动，中间体重排",
                "key_features": ["两性离子结构", "P-C 和 C-O 键共存", "可旋转"]
            },
            {
                "step": 3,
                "name": "Oxaphosphetane 形成（环化）",
                "description": "Betaine 通过分子内环化形成四元环 oxaphosphetane 中间体。这一步是可逆的。",
                "electron_flow": "氧负离子进攻磷原子（分子内 SN2）",
                "key_features": ["四元环结构", "P-O-C-C 环", "立体化学在此步骤确定"]
            },
            {
                "step": 4,
                "name": "消除（Elimination）",
                "description": "Oxaphosphetane 经协同消除分解为烯烃和氧化三苯基膦。这一步是不可逆的，驱动反应完成。",
                "electron_flow": "P-O 键断裂，C=C π键形成（协同过程）",
                "key_features": ["强 P=O 键形成提供驱动力（ΔG < 0）", "生成稳定烯烃", "Ph3P=O 为副产物"]
            }
        ],
        "stereochemistry": {
            "non_stabilized_ylide": {
                "type": "非稳定叶立德（alkyl substituents）",
                "selectivity": "Z-烯烃为主（动力学控制）",
                "explanation": "由于 oxaphosphetane 形成过程中的立体电子效应，非稳定叶立德倾向于给出 Z-烯烃。betaine 在消除前不能充分旋转。"
            },
            "stabilized_ylide": {
                "type": "稳定叶立德（EWG 如 COR, COOR, CN 取代）",
                "selectivity": "E-烯烃为主（热力学控制）",
                "explanation": "稳定叶立德反应较慢，betaine 有时间达到平衡，优先生成热力学更稳定的 E-烯烃。"
            },
            "semi_stabilized_ylide": {
                "type": "半稳定叶立德（aryl, vinyl 取代）",
                "selectivity": "E/Z 混合物，通常 E 优势",
                "explanation": "中等活性，立体选择性介于两者之间。"
            }
        }
    },
    "variants": [
        {
            "name": "Standard Wittig Reaction",
            "substrate": "醛/酮 + Ph3P=CHR",
            "product": "R'CH=CR2 + Ph3P=O",
            "conditions": "无水条件，THF 或 DCM 溶剂，0°C 至回流"
        },
        {
            "name": "Horner-Wadsworth-Emmons (HWE) Reaction",
            "substrate": "醛/酮 + phosphonate ester",
            "product": "E-烯烃（高选择性）",
            "advantages": "副产物水溶性磷酸盐易分离，E 选择性高"
        },
        {
            "name": "Wittig-Horner Reaction",
            "substrate": "醛/酮 + phosphonate anion",
            "product": "α,β-不饱和酯/酮",
            "note": "HWE 的变体"
        },
        {
            "name": "Schlosser Modification",
            "substrate": "标准 Wittig 条件 + LiBr/Li盐",
            "product": "将 Z-选择翻转为 E-选择",
            "mechanism": "通过锂盐促进 betaine 可逆化和 E-oxaphosphetane 优先形成"
        }
    ],
    "common_examples": [
        {
            "reactants": "benzaldehyde + Ph3P=CH2 (methylenetriphenylphosphorane)",
            "product": "styrene (PhCH=CH2)",
            "note": "最简单的 Wittig 反应示例"
        },
        {
            "reactants": "acetone + Ph3P=CHCH3",
            "product": "2-methylpropene",
            "note": "酮类底物的 Wittig 反应"
        },
        {
            "reactants": "cyclohexanone + Ph3P=CHCO2Et (stabilized ylide)",
            "product": "ethyl cyclohexenecarboxylate (E-major)",
            "note": "稳定叶立德，E 选择性"
        }
    ],
    "limitations": [
        "对空间位阻大的酮效果较差",
        "需要严格无水无氧条件（叶立德对水和空气敏感）",
        "三苯基膦衍生物分子量大，原子经济性差",
        "某些情况下 E/Z 选择性不够理想"
    ]
}


@ChemMCPManager.register_tool
class WittigMechanism(BaseTool):
    """
    Wittig 反应机理分析工具。
    提供完整的 Wittig 反应机理步骤、立体化学、变体及示例。
    """
    __version__ = "0.1.0"
    name = "WittigMechanism"
    func_name = "wittig_mechanism_analysis"
    description = "Analyze and explain the Wittig reaction mechanism including salt formation, betaine intermediate, oxaphosphetane formation, and elimination steps."
    implementation_description = "Knowledge-based tool using a built-in database of Wittig reaction mechanisms, stereochemical rules, variants, and examples."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Mechanism", "Wittig", "Ylide", "Alkene Synthesis", "Stereochemistry"]
    required_envs = []

    code_input_sig = [
        ("query", "str", "N/A", "Query about Wittig reaction: 'general', 'steps', 'stereochemistry', 'variants', 'examples', or a specific question."),
        ("detail_level", "str", "standard", "Detail level: 'brief', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated query and optional detail level, e.g., 'stereochemistry detailed'."),
    ]

    output_sig = [
        ("result", "str", "Detailed analysis result of the Wittig reaction mechanism."),
    ]

    examples = [
        {
            "code_input": {"query": "steps", "detail_level": "detailed"},
            "text_input": {"query_text": "steps detailed"},
            "output": {"result": "1. 成盐..."}
        },
        {
            "code_input": {"query": "stereochemistry", "detail_level": "standard"},
            "text_input": {"query_text": "stereochemistry"},
            "output": {"result": "..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, query: str, detail_level: str = "standard") -> str:
        """Core logic for Wittig mechanism analysis."""
        q = query.strip().lower()
        dl = detail_level.lower()

        if q in ("general", "overview", "intro", ""):
            return self._format_general(dl)
        elif q in ("steps", "mechanism", "step-by-step"):
            return self._format_steps(dl)
        elif q in ("stereochemistry", "stereoselectivity", "e/z"):
            return self._format_stereochemistry(dl)
        elif q in ("variants", "modifications", "hwe"):
            return self._format_variants(dl)
        elif q in ("examples", "example"):
            return self._format_examples(dl)
        elif q in ("limitations", "disadvantages"):
            return "\n".join(f"- {item}" for item in WITTIG_MECHANISM_DATA["limitations"])
        else:
            # Try to answer as a general query
            return self._search_query(query, dl)

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split(None, 1)
        query = parts[0]
        detail_level = parts[1] if len(parts) > 1 else "standard"
        return self._run_base(query, detail_level)

    def _format_general(self, detail_level: str) -> str:
        data = WITTIG_MECHANISM_DATA["general"]
        lines = [
            f"## {data['name']} 反应机理概述",
            "",
            data["overview"],
            "",
            f"**总反应式**: R₁R₂C=O + Ph₃P=CR₃R₄ → R₁R₂C=CR₃R₄ + Ph₃P=O",
            "",
            f"### 反应阶段 ({len(data['steps'])} 步)",
        ]
        for step in data["steps"]:
            lines.append(f"**{step['step']}. {step['name']}**")
            if detail_level != "brief":
                lines.append(f"   {step['description']}")
            if detail_level == "detailed":
                lines.append(f"   电子流向: {step['electron_flow']}")
                lines.append(f"   关键特征: {', '.join(step['key_features'])}")
            lines.append("")
        lines.append("### 立体化学要点")
        sc = data["stereochemistry"]
        for key, val in sc.items():
            lines.append(f"- **{val['type']}**: {val['selectivity']}")
            if detail_level != "brief":
                lines.append(f"  {val['explanation']}")
        return "\n".join(lines)

    def _format_steps(self, detail_level: str) -> str:
        data = WITTIG_MECHANISM_DATA["general"]["steps"]
        lines = ["## Wittig 反应分步机理", ""]
        for step in data:
            lines.append(f"### 第 {step['step']} 步: {step['name']}")
            lines.append(f"{step['description']}")
            if detail_level in ("standard", "detailed"):
                lines.append(f"\n**电子流向**: {step['electron_flow']}")
            if detail_level == "detailed":
                lines.append(f"\n**关键特征**:")
                for feat in step["key_features"]:
                    lines.append(f"  • {feat}")
            lines.append("")
        return "\n".join(lines)

    def _format_stereochemistry(self, detail_level: str) -> str:
        data = WITTIG_MECHANISM_DATA["general"]["stereochemistry"]
        lines = ["## Wittig 反应立体化学", ""]
        for key, val in data.items():
            lines.append(f"### {val['type']}")
            lines.append(f"- **选择性**: {val['selectivity']}")
            lines.append(f"- **解释**: {val['explanation']}")
            lines.append("")
        lines.append("### 总结规律")
        lines.append("| 叶立德类型 | 取代基 | 主要产物 | 控制类型 |")
        lines.append("|-----------|--------|---------|---------|")
        lines.append("| 非稳定型 | alkyl | **Z-烯烃** | 动力学 |")
        lines.append("| 半稳定型 | aryl, vinyl | **E-烯烃** (优势) | 混合 |")
        lines.append("| 稳定型 | EWG (COR, COOR, CN) | **E-烯烃** | 热力学 |")
        return "\n".join(lines)

    def _format_variants(self, detail_level: str) -> str:
        lines = ["## Wittig 反应主要变体", ""]
        for v in WITTIG_MECHANISM_DATA["variants"]:
            lines.append(f"### {v['name']}")
            lines.append(f"- **底物**: {v['substrate']}")
            lines.append(f"- **产物**: {v['product']}")
            if "conditions" in v:
                lines.append(f"- **条件**: {v['conditions']}")
            if "advantages" in v:
                lines.append(f"- **优点**: {v['advantages']}")
            if "note" in v:
                lines.append(f"- **备注**: {v['note']}")
            lines.append("")
        return "\n".join(lines)

    def _format_examples(self, detail_level: str) -> str:
        lines = ["## Wittig 反应典型实例", ""]
        for i, ex in enumerate(WITTIG_MECHANISM_DATA["common_examples"], 1):
            lines.append(f"### 例 {i}")
            lines.append(f"- **反应物**: {ex['reactants']}")
            lines.append(f"- **产物**: {ex['product']}")
            lines.append(f"- **说明**: {ex['note']}")
            lines.append("")
        return "\n".join(lines)

    def _search_query(self, query: str, detail_level: str) -> str:
        """Search through knowledge base for relevant info."""
        results = []
        # Search in various sections
        if any(kw in query for kw in ["betaine", "成盐", "salt", "step 1", "第一步"]):
            results.append(self._format_steps(detail_level))
        if any(kw in query for kw in ["stereo", "立体", "e/z", "e z", "selectivity"]):
            results.append(self._format_stereochemistry(detail_level))
        if any(kw in query for kw in ["variant", "变体", "hwe", "horner"]):
            results.append(self._format_variants(detail_level))
        if any(kw in query for kw in ["oxaphosphetane", "消除", "elimination"]):
            results.append(self._format_steps(detail_level))

        if results:
            return "\n\n---\n\n".join(results)

        # Default: return overview
        return self._format_general(detail_level)
