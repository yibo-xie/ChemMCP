import logging
from typing import Optional, List, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 周环反应知识库
PERICYCLIC_DATA = {
    "electrocyclic": {
        "name": "电环化反应 (Electrocyclic Reactions",
        "description": "共轭多烯末端碳原子之间形成 σ 键（或其逆过程），π 电子重新排列为 σ 键。",
        "rules": [
            "**Woodward-Hoffmann 规则** (基于轨道对称性):",
            "• 热反应: 4n π电子 → 对旋(disrotatory); 4n+2 π电子 → 顺旋(conrotatory)",
            "• 光化学反应: 4n π电子 → 顺旋(conrotatory); 4n+2 π电子 → 对旋(disrotatory)",
            "",
            "**Dewar-Zimmerman 规则** (基于芳香性过渡态):",
            "• 芳香性过渡态(4n+2 e⁻) → 热允许",
            "• 反芳香性过渡态(4n e⁻) → 光化学允许"
        ],
        "examples": [
            {"reactants": "1,3,5-己三烯 (6π e⁻)", "thermal": "顺旋 → 环己二烯", "photochemical": "对旋 → 环己二烯"},
            {"reactants": "丁二烯 (4π e⁻)", "thermal": "对旋 → 环丁烯", "photochemical": "顺旋 → 环丁烯"},
            {"reactants": "环丁烯开环", "thermal": "对旋 → 丁二烯 (4e⁻)", "photochemical": "顺旋 → 丁二烯"},
        ],
        "stereochemistry": "立体专一性 — 旋转方式决定产物构型"
    },
    "cycloaddition": {
        "name": "环加成反应 (Cycloadditions",
        "description": "两个或多个 π 体系通过形成两个新的 σ 键而连接成环的反应。",
        "classification": [
            "[i+j] 命名法: i 和 j 分别为两个组分参与反应的 π 电子数",
            "[4+2]: Diels-Alder 反应 — 最重要",
            "[2+2]: 酮/烯光化学环加成",
            "[3+2]: 1,3-偶极环加成",
            "[4+1], [5+2], [6+3] 等: 较少见但已知"
        ],
        "selection_rules": [
            "**热反应**: π电子总数 = 4n+2 → 允许 (suprafacial-suprafacial)",
            "**光化学反应**: π电子总数 = 4n → 允许 (suprafacial-suprafacial)",
            "",
            "| 反应类型 | π电子总数 | 热反应 | 光反应 |",
            "|---------|-----------|--------|--------|",
            "| Diels-Alder [4+2] | 6 (4n+2) | ✅ 允许 | ❌ 禁阻 |",
            "| [2+2] 烯烃二聚 | 4 (4n) | ❌ 禁阻 | ✅ 允许 |",
            "| [4+4] 二聚 | 8 (4n) | ❌ 禁阻 | ✅ 允许 |",
            "| [6+4] | 10 (4n+2) | ✅ 允许 | ❌ 禁阻 |"
        ],
        "examples": [
            {"name": "Diels-Alder", "description": "diene + dienophile -> cyclohexene", "stereo": "endo 选择性, suprafacial"},
            {"name": "[2+2] Photo", "description": "ketene + alkene -> cyclobutanone", "stereo": "通常非立体专一"},
            {"name": "1,3-Dipolar", "description": "azoalkyne + alkene -> pyrazole", "stereo": "取决于具体体系"},
        ]
    },
    "sigmatropic": {
        "name": "σ迁移反应 (Sigmatropic Rearrangements)",
        "description": "σ键沿着共轭 π 体系迁移到新位置，同时 π 键重新排列。",
        "notation": "[i,j]-迁移: σ键从位置 i 迁移到位置 j (编号方向一致)",
        "common_types": [
            "[1,3]-H shift: 少见（HOMO对称性不匹配，除非有特殊取代）",
            "[1,5]-H shift: 常见，热允许（如 1,3-戊二烯）",
            "[1,7]-H shift: 可在链状体系中发生",
            "[3,3]-Cope rearrangement: 1,5-二烯 → 另一个1,5-二烯",
            "[3,3]-Claisen rearrangement: 烯丙基乙烯基醚 → γ,δ-不饱和羰基",
            "[2,3]-Wittig rearrangement: 醚 α-碳负离子重排",
            "[5,5]: 在某些天然产物合成中观察到"
        ],
        "selection_rules": [
            "**[i,j]-σ迁移选择规则**:",
            "• 热反应: (i+j) = 4n+2 → suprafacial 允许; (i+j) = 4n → antarafacial 允许",
            "• H 迁移: [1,n]-H shift 中 n 为奇数时热允许(suprafacial), n 为偶数时需 antarafacial",
            "• C 迁移: 类似规则，但 C 可以翻转构型(antarafacial)",
            "",
            "| 迁移类型 | i+j | 热条件 | 常见程度 |",
            "|---------|-----|--------|---------|",
            "| [1,3]-H | 4 (4n) | ❌ 需antarafacial | 稀少 |",
            "| [1,5]-H | 6 (4n+2) | ✅ suprafacial | **常见** |",
            "| [1,7]-H | 8 (4n) | ❌ 需antarafacial | 链状可发生 |",
            "| [3,3]-Cope | 6 (4n+2) | ✅ suprafacial | **非常常见** |",
            "| [3,3]-Claisen | 6 (4n+2) | ✅ suprafacial | **非常常见** |",
            "| [2,3] | 5 (4n+1) | 特殊处理 | 常见 |"
        ],
        "key_examples": [
            ("Cope 重排", "1,5-己二烯 → 另一个 1,5-己二烯", "[3,3]-迁移"),
            ("Claisen 重排", "烯丙基乙烯基醚 → 4-戊烯醛", "[3,3]-迁移, 合成上有用"),
            ("Cope 消除", "胺氧化物 → 烯烃 + N,N-二甲羟胺", "[3,3] 同面, syn消除"),
            ("Fischer 吲哚合成", "苯腙 → 吲哚", "[3,3]-σ迁移关键步骤"),
        ]
    },
    "chirality_transfer": {
        "name": "手性转移 (Chirality Transfer in Pericyclic Reactions",
        "concepts": [
            "电环化反应中，对旋(conrotation)保持轴手性关系；对旋(disrotation)反转",
            "[1,n]-σ迁移中的手性: suprafacial H迁移保持构型; C迁移可翻转",
            "Cope/Claisen [3,3]迁移: 过渡态类似椅式/船式构象，影响立体化学",
            "周环反应的立体专一性使其成为不对称合成的有力工具"
        ]
    }
}


@ChemMCPManager.register_tool
class PericyclicAnalyzer(BaseTool):
    """
    周环反应分析工具。
    支持电环化、环加成、σ迁移等周环反应的选择规则、机理和立体化学分析。
    """
    __version__ = "0.1.0"
    name = "PericyclicAnalyzer"
    func_name = "analyze_pericyclic_reaction"
    description = "Analyze pericyclic reactions including electrocyclic reactions, cycloadditions (Diels-Alder, etc.), and sigmatropic rearrangements with Woodward-Hoffmann rules."
    implementation_description = "Knowledge-based tool implementing Woodward-Hoffmann and Dewar-Zimmerman selection rules for pericyclic reactions."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Pericyclic", "Electrocyclic", "Cycloaddition", "Sigmatropic", "Woodward-Hoffmann", "Orbital Symmetry"]
    required_envs = []

    code_input_sig = [
        ("reaction_type", "str", "N/A", "Type of pericyclic reaction: 'electrocyclic', 'cycloaddition', 'sigmatropic', 'all', or a specific query."),
        ("pi_electrons", "int", "N/A", "Number of π electrons involved (for electrocyclic/cycloaddition analysis). Use -1 if not applicable."),
        ("condition", "str", "thermal", "Reaction condition: 'thermal' or 'photochemical'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated parameters: 'reaction_type pi_electrons condition', e.g., 'electrocyclic 4 thermal' or 'sigmatropic detailed'."),
    ]

    output_sig = [
        ("result", "str", "Analysis result including allowed/forbidden prediction, mechanism details, and stereochemistry."),
    ]

    examples = [
        {
            "code_input": {"reaction_type": "electrocyclic", "pi_electrons": 4, "condition": "thermal"},
            "text_input": {"query_text": "electrocyclic 4 thermal"},
            "output": {"result": "## 电环化反应分析 (4π e⁻, 热反应)..."}
        },
        {
            "code_input": {"reaction_type": "cycloaddition", "pi_electrons": 6, "condition": "thermal"},
            "text_input": {"query_text": "cycloaddition 6 thermal"},
            "output": {"result": "## [4+2] Diels-Alder 反应..."}
        },
        {
            "code_input": {"reaction_type": "sigmatropic", "pi_electrons": -1, "condition": "thermal"},
            "text_input": {"query_text": "sigmatropic"},
            "output": {"result": "## σ迁移反应..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, pi_electrons: int = -1, condition: str = "thermal") -> str:
        rt = reaction_type.strip().lower()
        cond = condition.strip().lower()

        if rt in ("electrocyclic", "electro"):
            return self._analyze_electrocyclic(pi_electrons, cond)
        elif rt in ("cycloaddition", "cycloadd", "diels-alder", "dielsalder", "da"):
            return self._analyze_cycloaddition(pi_electrons, cond)
        elif rt in ("sigmatropic", "sigma", "sigmap", "cope", "claisen"):
            return self._analyze_sigmatropic(cond)
        elif rt in ("all", "overview", "list", ""):
            return self._format_overview()
        else:
            # Try to find best match
            if any(kw in rt for kw in ["electro"]):
                return self._analyze_electrocyclic(pi_electrons, cond)
            elif any(kw in rt for kw in ["cyclo", "diels"]):
                return self._analyze_cycloaddition(pi_electrons, cond)
            elif any(kw in rt for kw in ["sigma", "migrat", "cope", "claisen"]):
                return self._analyze_sigmatropic(cond)
            return self._format_overview()

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split()
        rtype = parts[0] if len(parts) > 0 else "all"
        pi_e = int(parts[1]) if len(parts) > 1 and parts[1].lstrip("-").isdigit() else -1
        cond = parts[2] if len(parts) > 2 else "thermal"
        return self._run_base(rtype, pi_e, cond)

    def _is_4n_plus_2(self, n: int) -> bool:
        """Check if n is of form 4k+2."""
        return (n - 2) % 4 == 0 and n >= 2

    def _is_4n(self, n: int) -> bool:
        """Check if n is of form 4k."""
        return n % 4 == 0 and n > 0

    def _analyze_electrocyclic(self, pi_e: int, cond: str) -> str:
        lines = [f"## 电环化反应分析 ({pi_e}π 电子, {cond}反应)", ""]
        
        if pi_e < 0:
            data = PERICYCLIC_DATA["electrocyclic"]
            lines.append(data["description"])
            lines.append("")
            for rule in data["rules"]:
                lines.append(rule)
                lines.append("")
            lines.append("### 示例")
            for ex in data["examples"]:
                lines.append(f"- {ex['reactants']}:")
                lines.append(f"  - 热: {ex['thermal']}")
                lines.append(f"  - 光: {ex['photochemical']}")
            return "\n".join(lines)

        is_4np2 = self._is_4n_plus_2(pi_e)
        is_4n = self._is_4n(pi_e)

        lines.append(f"### Woodward-Hoffmann 规则判断")
        lines.append(f"")
        lines.append(f"**π 电子数**: {pi_e}")
        lines.append(f"**分类**: {'4n+2 型' if is_4np2 else '4n 型' if is_4n else '其他'}")
        lines.append(f"")

        if cond in ("thermal", "heat", "Δ"):
            if is_4np2:
                mode = "顺旋 (conrotatory)"
                allowed = "✅ **对称性允许**"
            elif is_4n:
                mode = "对旋 (disrotatory)"
                allowed = "✅ **对称性允许**"
            else:
                mode = "需进一步分析"
                allowed = "⚠️ 非标准 π 电子数"
        elif cond in ("photochemical", "photo", "hν", "hv"):
            if is_4n:
                mode = "顺旋 (conrotatory)"
                allowed = "✅ **对称性允许**"
            elif is_4np2:
                mode = "对旋 (disrotatory)"
                allowed = "✅ **对称性允许**"
            else:
                mode = "需进一步分析"
                allowed = "⚠️ 非标准 π 电子数"
        else:
            mode = "未知"
            allowed = f"❓ 未知条件 '{cond}'"

        lines.append(f"| 条件 | 旋转方式 | 判断 |")
        lines.append(f"|------|---------|------|")
        lines.append(f"| {cond} | {mode} | {allowed} |")

        lines.append("")
        lines.append("### Dewar-Zimmerman 分析")
        if cond in ("thermal", "heat", "Δ"):
            ts_aromatic = is_4np2
            lines.append(f"- 过渡态 {'具有芳香性 (4n+2 e⁻)' if ts_aromatic else '具有反芳香性 (4n e⁻)'}")
            lines.append(f"- 芳香性过渡态 → 热反应{'允许 ✅' if ts_aromatic else '禁阻 ❌'}")
        else:
            ts_aromatic = is_4n
            lines.append(f"- 过渡态 {'具有反芳香性 (4n e⁻)' if not ts_aromatic else '具有芳香性 (4n+2 e⁻)'}")
            lines.append(f"- 光化学反应通过激发态改变轨道对称性")

        lines.append("")
        lines.append("### 立体化学结果")
        if "顺旋" in mode:
            lines.append("- 顺旋: 两端原子同向旋转（都顺时针或都逆时针）")
            lines.append("- 对于链状体系: 决定产物的相对构型（E/Z 或 R/S）")
        elif "对旋" in mode:
            lines.append("- 对旋: 两端原子反向旋转（一个顺时针一个逆时针）")
            lines.append("- 对于链状体系: 决定产物的相对构型")

        return "\n".join(lines)

    def _analyze_cycloaddition(self, pi_e: int, cond: str) -> str:
        lines = [f"## 环加成反应分析 (总 π 电子数={pi_e}, {cond})", ""]
        
        if pi_e < 0:
            data = PERICYCLIC_DATA["cycloaddition"]
            lines.append(data["description"])
            lines.append("")
            lines.append("**分类体系:**")
            for c in data["classification"]:
                lines.append(f"- {c}")
            lines.append("")
            lines.append("### 选择规则")
            for rule in data["selection_rules"]:
                lines.append(rule)
                lines.append("")
            return "\n".join(lines)

        is_4np2 = self._is_4n_plus_2(pi_e)
        is_4n = self._is_4n(pi_e)

        lines.append("### Woodward-Hoffmann 规则判断 (suprafacial-suprafacial)")
        lines.append("")
        lines.append(f"**总 π 电子数**: {pi_e}")

        if cond in ("thermal", "heat", "Δ"):
            if is_4np2:
                allowed = "✅ **对称性允许** (suprafacial-suprafacial)"
            elif is_4n:
                allowed = "❌ **对称性禁阻** (suprafacial-suprafacial); 需要 antarafacial 组分"
            else:
                allowed = "⚠️ 非标准电子数"
        elif cond in ("photochemical", "photo", "hν"):
            if is_4n:
                allowed = "✅ **对称性允许** (suprafacial-suprafacial)"
            elif is_4np2:
                allowed = "❌ **对称性禁阻** (suprafacial-suprafacial)"
            else:
                allowed = "⚠️ 非标准电子数"
        else:
            allowed = f"❓ 未知条件"

        lines.append(f"**判断**: {allowed}")
        lines.append("")

        # Classify by component sizes
        classifications = {
            4: "[2+2] 烯酮/烯烃环加成 或 [4+0](电环化)",
            6: "[4+2] Diels-Alder 反应 (最常见)",
            8: "[4+4] 烯烃二聚 或 [6+2]",
            10: "[6+4] 或 [8+2]",
            2: "[2+0] (实际为电环化)",
            3: "[2+1] (卡宾加成, 非周环)",
            5: "[4+1] 或 [3+2] (1,3-偶极)",
        }
        classification = classifications.get(pi_e, "未知组合")
        lines.append(f"**可能类型**: {classification}")

        if pi_e == 6 and cond in ("thermal", "heat", "Δ"):
            lines.append("")
            lines.append("### 🎯 Diels-Alder 反应要点")
            lines.append("- 双烯体(HOMO) + 亲双烯体(LUMO) 轨道相互作用")
            lines.append("- 正常电子需求: 富电子双烯体 + 缺电子亲双烯体")
            lines.append("- 反向电子需求: 缺电子双烯体 + 富电子亲双烯体")
            lines.append("- 内型(endo)选择性: 次级轨道相互作用驱动")
            lines.append("- 立体专一性: 相对构型完全保留")

        return "\n".join(lines)

    def _analyze_sigmatropic(self, cond: str) -> str:
        data = PERICYCLIC_DATA["sigmatropic"]
        lines = ["## σ迁移反应分析", ""]
        lines.append(data.get("description", "σ键沿共轭π体系迁移到新位置"))
        lines.append("")
        lines.append("### 命名法: [i,j]-迁移")
        lines.append("- i: 断裂的σ键起始位置")
        lines.append("- j: 新形成的σ键目标位置")
        lines.append("")
        lines.append("### 常见类型与选择规则")
        for rule in data.get("selection_rules", []):
            lines.append(rule)
            lines.append("")
        lines.append("### 重要实例")
        for ex in data.get("key_examples", []):
            lines.append(f"- **{ex[0]}**: {ex[1]} ({ex[2]})")
        return "\n".join(lines)

    def _format_overview(self) -> str:
        lines = ["## 周环反应总览", ""]
        lines.append("| 类型 | 描述 | 代表反应 | 关键规则 |")
        lines.append("|------|------|---------|---------|")
        lines.append("| 电环化 | π→σ 相互转化 | 丁二烯↔环丁烯 | W-H 旋转模式 |")
        lines.append("| 环加成 | 两分子π→两根σ | Diels-Alder [4+2] | W-H 电子计数 |")
        lines.append("| σ迁移 | σ键沿π体系移动 | Cope, Claisen [3,3] | W-H 面/迁移距离 |")
        lines.append("")
        lines.append("### 核心规则: Woodward-Hoffmann")
        lines.append("- 基于分子轨道对称性守恒原理")
        lines.append("- 热反应: 轨道相位匹配决定允许/禁阻")
        lines.append("- 光化学反应: 激发态改变轨道相位，规则翻转")
        lines.append("- Dewar-Zimmerman: 过渡态芳香性判据（等价表述）")
        return "\n".join(lines)
