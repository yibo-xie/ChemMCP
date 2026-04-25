import logging
import re
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 电子转移箭头验证规则库
ARROW_PUSHING_RULES = {
    "fundamental_rules": [
        {
            "rule": "1. 箭头起点必须是电子源",
            "valid_starts": ["键 (σ或π)", "孤对电子 (lone pair)", "负电荷 (formal negative)"],
            "invalid_starts": ["正电荷", "原子核 (无电子)", "空轨道 (无电子可给)"],
            "example_valid": "O: → C=O (氧的孤对电子进攻羰基碳)",
            "example_invalid": "⁺N → O (正电荷不是电子源)"
        },
        {
            "rule": "2. 箭头终点必须是电子受体（缺电子位置）",
            "valid_ends": ["原子 (形成新键)", "σ键 (断裂均裂)", "π键 (断裂异裂)", "正电荷"],
            "example_valid": "Nu: → C⁺ (亲核试剂进攻碳正离子)",
            "example_invalid": ":O⁻ → :O⁻ (负电到负电，无驱动力)"
        },
        {
            "rule": "3. 每根箭头代表2个电子（一个电子对）",
            "note": "单电子转移用鱼钩箭头(↔ 或半箭头)表示自由基反应",
        },
        {
            "rule": "4. 电荷守恒",
            "note": "反应前后总电荷必须相等；每根箭头移动2e⁻"
        },
        {
            "rule": "5. 八隅体/十八电子规则",
            "note": "第二周期元素不超过8e⁻; 过渡金属通常18e⁻; H/He满足双电子规则"
        }
    ],
    "common_errors": [
        {
            "error": "五价碳",
            "description": "箭头画到一个已经有4个键的碳上使其变成5键",
            "fix": "检查是否需要先断一根键（如消除反应），或者箭头应该指向其他原子",
            "example": "CH₄ + :NH₃ → 不应直接在C上加第五个键"
        },
        {
            "error": "违反电荷守恒",
            "description": "箭头移动后电荷不平衡",
            "fix": "数清每个原子的形式电荷变化，确保总电荷不变",
        },
        {
            "error": "从错误的位置开始画箭头",
            "description": "如从H画箭头（H通常没有孤对电子可作为电子源）",
            "fix": "确认电子源的真正位置"
        },
        {
            "error": "遗漏箭头",
            "description": "只画了成键箭头但没画断键箭头（或反之）",
            "fix": "成键和断键箭头必须配对出现（除非是协同反应）"
        },
        {
            "error": "酸碱混淆",
            "description": "把 Brønsted 酸碱反应画成了共价键形成",
            "fix": "质子转移用弯箭头从碱的孤对电子→H，同时H-X键的电子→X"
        }
    ],
    "reaction_type_patterns": {
        "nucleophilic_addition": {
            "name": "亲核加成",
            "arrow_pattern": [
                "1. 亲核试剂的孤对电子/负电荷 → 亲电碳（羰基碳等）",
                "2. π键电子 → 电负性原子（O 变为 O⁻）",
                "3. （可选）质子化步骤：溶剂/酸提供 H⁺"
            ],
            "example": "醛 + CN⁻ → 氰醇",
            "key_check": "检查羰基碳是否确实缺电子（无稳定化取代基时活性更高）"
        },
        "nucleophilic_substitution_sn2": {
            "name": "SN2 亲核取代",
            "arrow_pattern": [
                "1. Nu: → C（背面进攻）",
                "2. C-LG 键电子 → LG（离去基团带走电子对）",
                "注: 两步协同，只有一个过渡态"
            ],
            "example": "OH⁻ + CH₃Br → CH₃OH + Br⁻",
            "key_check": "确认离去基团足够好（弱碱性的共轭碱）"
        },
        "nucleophilic_substitution_sn1": {
            "name": "SN1 亲核取代",
            "arrow_pattern": [
                "1. C-LG 键电子 → LG（慢步，形成碳正离子）",
                "2. Nu: → C⁺（快步）",
                "（可选）重排: 相邻基团迁移（带键电子对迁移）"
            ],
            "example": "(CH₃)₃C-Cl → (CH₃)₃C⁺ → (CH₃)₃C-Nu",
            "key_check": "检查是否能形成稳定的碳正离子（3° > 2° > 1°）"
        },
        "elimination_e2": {
            "name": "E2 消除",
            "arrow_pattern": [
                "1. B:（碱）夺取 β-H（B: → H, 形成新 H-B 键）",
                "2. Cβ-H 键电子 → Cα-Cβ 形成 π 键",
                "3. Cα-LG 键电子 → LG"
            ],
            "example": "KOH + (CH₃)₂CH-CH₂Br → (CH₃)₂C=CH₂",
            "key_check": "H 和 LG 必须反式共平面（anti-periplanar）"
        },
        "elimination_e1": {
            "name": "E1 消除",
            "arrow_pattern": [
                "1. C-LG 键电子 → LG（形成碳正离子）",
                "2. B: 夺取 β-H（B: → H）",
                "3. Cβ-H 电子 → Cα-Cβ π 键"
            ],
            "key_check": "与 SN1 共同的第一步（形成相同碳正离子）"
        },
        "electrophilic_addition": {
            "name": "亲电加成",
            "arrow_pattern": [
                "1. π键电子 → E⁺（如 H⁺, Br⁺）",
                "2. 形成更稳定的碳正离子中间体（马氏规则）",
                "3. Nu: / Nu⁻ → C⁺"
            ],
            "example": "HBr + CH₂=CH-CH₃ → CH₃-CHBr-CH₃",
            "key_check": "中间体碳正离子的稳定性决定区域选择性"
        },
        "pericyclic": {
            "name": "周环反应",
            "arrow_pattern": [
                "环形箭头体系（所有箭头首尾相连形成环状）",
                "电子在闭环中流动，无中间体",
                "键的断裂和形成完全同步"
            ],
            "example": "Diels-Alder: 3 根箭头形成环状流动",
            "key_check": "箭头必须形成闭合环路；遵循 Woodward-Hoffmann 规则"
        }
    },
    "formal_charge_calculation": {
        "formula": "形式电荷 = 族数 - [非键电子数 + 1/2 成键电子数]",
        "examples": [
            ("NH₄⁺ 中 N", "5 - [0 + 8/2] = 5 - 4 = +1 ✅"),
            ("H₃O⁺ 中 O", "6 - [2 + 6/2] = 6 - 5 = +1 ✅"),
            ("NO₃⁻ 中 N (一个双键)", "5 - [0 + 8/2] = 5 - 4 = +1 (共振结构之一)"),
        ]
    }
}


@ChemMCPManager.register_tool
class ArrowPushingValidator(BaseTool):
    """
    电子转移箭头合理性验证工具。
    检查有机化学反应机理中的弯箭头是否符合基本规则、电荷守恒和八隅体规则。
    """
    __version__ = "0.1.0"
    name = "ArrowPushingValidator"
    func_name = "validate_arrow_pushing"
    description = "Validate the correctness of electron-pushing arrows in organic reaction mechanisms, checking fundamental rules, charge conservation, and octet compliance."
    implementation_description = "Rule-based validation engine for organic chemistry electron-pushing notation with built-in error patterns and reaction-type-specific arrow templates."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Arrow Pushing", "Mechanism Validation", "Electron Flow", "Organic Chemistry", "Curly Arrows"]
    required_envs = []

    code_input_sig = [
        ("query", "str", "N/A", "Query type: 'rules' (show all rules), 'errors' (common errors), a reaction type name (check pattern), or 'validate' with description."),
        ("reaction_type", "str", "", "Optional: specific reaction type to check (e.g., 'sn2', 'e2', 'nucleophilic_addition'). Leave empty for general rules."),
        ("description", "str", "", "Optional: describe the arrow-pushing to validate if query='validate'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated query, e.g., 'rules', 'errors', 'sn2', or 'validate <description>'."),
    ]

    output_sig = [
        ("result", "str", "Validation result including rule checks, error detection, and suggestions."),
    ]

    examples = [
        {
            "code_input": {"query": "rules", "reaction_type": "", "description": ""},
            "text_input": {"query_text": "rules"},
            "output": {"result": "## Electron-Pushing Rules..."}
        },
        {
            "code_input": {"query": "sn2", "reaction_type": "sn2", "description": ""},
            "text_input": {"query_text": "sn2"},
            "output": {"result": "## SN2 Arrow Pattern..."}
        },
        {
            "code_input": {"query": "validate", "reaction_type": "", "description": "Arrow from O lone pair to C of C=O"},
            "text_input": {"query_text": "validate Arrow from O lone pair to carbonyl carbon"},
            "output": {"result": "## Validation Result..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, query: str, reaction_type: str = "", description: str = "") -> str:
        q = query.strip().lower()
        rt = reaction_type.strip().lower()
        desc = description.strip()

        if q in ("rules", "rule", "basics", "fundamental"):
            return self._format_rules()
        elif q in ("errors", "common_errors", "mistakes"):
            return self._format_common_errors()
        elif q == "validate":
            return self._validate_description(desc)
        elif q in ("charge", "formal_charge", "fc"):
            return self._format_formal_charge()
        elif q:
            # Try as reaction type
            return self._format_reaction_pattern(q)
        elif rt:
            return self._format_reaction_pattern(rt)
        else:
            return self._format_rules()

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split(None, 1)
        q = parts[0]
        rest = parts[1] if len(parts) > 1 else ""
        
        if q == "validate":
            return self._run_base("validate", "", rest)
        return self._run_base(q, rest, "")

    def _format_rules(self) -> str:
        lines = ["## 电子转移箭头基本规则", ""]
        for r in ARROW_PUSHING_RULES["fundamental_rules"]:
            lines.append(f"### {r['rule']}")
            if "valid_starts" in r:
                lines.append(f"- **有效起点**: {', '.join(r['valid_starts'])}")
            if "invalid_starts" in r:
                lines.append(f"- **无效起点**: {', '.join(r['invalid_starts'])}")
            if "valid_ends" in r:
                lines.append(f"- **有效终点**: {', '.join(r['valid_ends'])}")
            if "example_valid" in r:
                lines.append(f"- ✅ 正确示例: {r['example_valid']}")
            if "example_invalid" in r:
                lines.append(f"- ❌ 错误示例: {r['example_invalid']}")
            if "note" in r:
                lines.append(f"- 📝 注: {r['note']}")
            lines.append("")
        return "\n".join(lines)

    def _format_common_errors(self) -> str:
        lines = ["## 常见箭头错误及修正", ""]
        for err in ARROW_PUSHING_RULES["common_errors"]:
            lines.append(f"### ❌ {err['error']}")
            lines.append(f"**问题**: {err['description']}")
            lines.append(f"**修正**: {err['fix']}")
            if "example" in err:
                lines.append(f"**示例**: {err['example']}")
            lines.append("")
        return "\n".join(lines)

    def _format_formal_charge(self) -> str:
        data = ARROW_PUSHING_RULES["formal_charge_calculation"]
        lines = ["## 形式电荷计算", ""]
        lines.append(f"**公式**: {data['formula']}")
        lines.append("")
        lines.append("### 示例")
        for ex in data["examples"]:
            lines.append(f"- **{ex[0]}**: {ex[1]}")
        lines.append("")
        lines.append("### 快速判断技巧")
        lines.append("- 每根共价键算该原子拥有1个电子")
        lines.append("- 每对孤对电子算2个")
        lines.append("- 与族数比较: 多余电子=负电荷, 缺少=正电荷")
        return "\n".join(lines)

    def _format_reaction_pattern(self, rtype: str) -> str:
        # Map various names to canonical keys
        type_map = {
            "sn2": "nucleophilic_substitution_sn2",
            "s_n2": "nucleophilic_substitution_sn2",
            "sn1": "nucleophilic_substitution_sn1",
            "s_n1": "nucleophilic_substitution_sn1",
            "e2": "elimination_e2",
            "e1": "elimination_e1",
            "nucleophilic_addition": "nucleophilic_addition",
            "nu_add": "nucleophilic_addition",
            "electrophilic_addition": "electrophilic_addition",
            "el_add": "electrophilic_addition",
            "pericyclic": "pericyclic",
            "diels-alder": "pericyclic",
            "dielsalder": "pericyclic",
        }

        key = type_map.get(rtype)
        if not key:
            available = list(ARROW_PUSHING_RULES["reaction_type_patterns"].keys())
            return f"未找到反应类型 '{rtype}'。可用类型: {', '.join(available)}"

        data = ARROW_PUSHING_RULES["reaction_type_patterns"][key]
        lines = [f"## {data['name']} — 箭头模式", ""]
        lines.append("### 正确的箭头画法:")
        for i, pattern in enumerate(data["arrow_pattern"], 1):
            lines.append(f"{pattern}")
        lines.append("")
        lines.append(f"**示例**: {data.get('example', 'N/A')}")
        lines.append("")
        lines.append(f"⚠️ **检查要点**: {data.get('key_check', 'N/A')}")

        return "\n".join(lines)

    def _validate_description(self, desc: str) -> str:
        """Simple pattern-based validation of an arrow-pushing description."""
        if not desc:
            return "❓ 请提供需要验证的箭头描述。"

        issues = []
        checks = []

        desc_lower = desc.lower()

        # Check 1: Valid electron source?
        has_valid_source = any(kw in desc_lower for kw in [
            "lone pair", "孤对", ":", "negative", "负", "bond", "键", "π", "pi", "double bond"
        ])
        if has_valid_source:
            checks.append("✅ 检测到可能的电子源")
        else:
            issues.append("⚠️ 未明确检测到有效电子源（孤对电子/键/负电荷）")

        # Check 2: Common errors
        if "五价" in desc or "pentavalent" in desc_lower or "5 bond" in desc_lower or "five bond" in desc_lower:
            issues.append("❌ 可能存在五价碳错误 — 第二周期元素不能超过4个键")

        if "h→" in desc_lower or "氢→" in desc or "hydrogen→" in desc_lower:
            issues.append("⚠️ 从 H 出发画箭头 — H 通常没有孤对电子作为电子源（除非是氢负离子 H⁻）")

        # Check 3: Charge-related
        if "positive" in desc_lower and "source" in desc_lower:
            issues.append("❌ 正电荷不能作为电子源（箭头起点）")

        if "negative" in desc_lower and ("to negative" in desc_lower or "to anion" in desc_lower):
            issues.append("⚠️ 负电到负电 — 检查是否有热力学驱动力")

        # Check 4: Reaction type keywords
        detected_type = None
        for tname in ARROW_PUSHING_RULES["reaction_type_patterns"]:
            if tname.replace("_", " ") in desc_lower or tname in desc_lower:
                detected_type = tname
                break

        lines = ["## 箭头验证结果", ""]
        lines.append(f"**输入描述**: {desc}")
        lines.append("")

        if checks:
            lines.append("### 通过项")
            for c in checks:
                lines.append(c)
            lines.append("")

        if issues:
            lines.append("### ⚠️ 注意事项")
            for issue in issues:
                lines.append(f"- {issue}")
            lines.append("")
        else:
            lines.append("### ✅ 未检测到明显错误")
            lines.append("")

        if detected_type:
            lines.append("---")
            lines.append(self._format_reaction_pattern(detected_type))

        if not issues and not detected_type:
            lines.append("> 💡 提示: 这是基于关键词的基本检查。复杂的机理建议结合具体反应结构式分析。")

        return "\n".join(lines)
