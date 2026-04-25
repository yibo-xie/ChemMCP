import logging
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 立体选择性知识库
STEREOSELECTIVITY_DATA = {
    "syn_anti_addition": {
        "name": "Syn vs Anti 加成立体化学",
        "description": "不同反应中两个基团从同侧(syn)还是异侧(anti)加到底物上。",
        "reactions": [
            {"reaction": "催化氢化 (H₂/Pd, Pt, Ni)", "mode": "**syn** 加成", "reason": "两个 H 从金属表面同侧转移"},
            {"reaction": "硼氢化-氧化 (BH₃ then H₂O₂/NaOH)", "mode": "**syn** 加成", "reason": "B 和 H 通过环状过渡态同时加成"},
            {"reaction": "OsO₄ 二羟基化", "mode": "**syn** 加成", "reason": "[3+2] 环加成形成酯环，顺式开环"},
            {"reaction": "KMnO₄ 冷稀溶液二羟基化", "mode": "**syn** 加成", "reason": "类似 OsO₄ 的 [3+2] 机理"},
            {"reaction": "X₂ (Br₂/Cl₂) 对烯烃卤化", "mode": "**anti** 加成", "reason": "卤鎓离子中间体迫使 X⁻ 从背面进攻"},
            {"reaction": "卤醇形成 (X₂/H₂O)", "mode": "**anti** 加成", "reason": "同样经过卤鎓离子中间体"},
            {"reaction": "环氧化物酸性开环", "mode": "**anti** 加成", "reason": "SN2-like 背面进攻"},
            {"reaction": "环氧化物碱性开环", "mode": "**anti** 加成", "reason": "SN2 机理"},
            {"reaction": "溴代醇内醚化 (Br₂)", "mode": "**anti** 加成", "reason": "溴鎓离子机理"},
        ]
    },
    "e2_stereochemistry": {
        "name": "E2 消除的立体化学",
        "requirement": "H-Cβ-Cα-LG 必须呈 **反式共平面 (anti-periplanar)** 构象",
        "dihedral_angle": "二面角 ≈ 180°",
        "consequences": [
            "在环己烷体系中: H 和 LG 必须 **双竖键 (diaxial)** 关系",
            "这决定了消除的立体化学结果",
            "某些构型下无法满足 anti-periplanar → 不能发生 E2（或极慢）"
        ],
        "examples": [
            {"substrate": "menthyl chloride (新薄荷基氯)", "result": "E2 极慢 — Cl 为平伏键(equatorial), 无反式共平面 H"},
            {"substrate": "neomenthyl chloride (新异薄荷基氯)", "result": "E2 快速 — Cl 为直立键(axial), 有反式共平面 H"},
            {"substrate": "trans-1,2-二溴环己烷 E2", "result": "生成环己炔 (2× E2, 两对 diaxial H/Br)"},
            {"substrate": "cis-1,2-二溴环己烷 E2", "result": "生成 1,3-环己二烯 (非炔烃, 因无 diaxial 排列)"},
        ]
    },
    "sn2_stereochemistry": {
        "name": "SN2 的立体化学",
        "configuration": "**Walden 翻转** (完全构型翻转)",
        "mechanism": "Nu 从离去基团背面进攻 → 过渡态为三角双锥 → LG 离去 → 构型反转",
        "example": "(R)-2-溴丁烷 + OH⁻ → (S)-2-丁醇 (100% 翻转)",
        "note": "这是 SN2 反应的特征性标志 — 如果观察到外消旋化则不是纯 SN2"
    },
    "alkene hydrogenation_stereo": {
        "name": "烯烃氢化的立体选择性",
        "syn_addition": "H₂ 在催化剂表面 **同侧** 加成",
        "factors": [
            "顺式(z/ cis)烯烃 → meso 产物 (内消旋) 或特定非对映体",
            "反式(E/trans)烯烃 → racemic 或 dl 对",
            "催化剂表面吸附可影响选择性（均相催化如 Wilkinson's 催化剂有更高选择性）"
        ],
        "lindlar": "Lindlar Pd/BaSO₄ + 喹啉: 炔烃 → **顺式(Z)烯烃** (syn 加成停在烯烃)",
        "na_nh3": "Na/液NH₃: 炔烃 → **反式(E)烯烃** (反式电子转移)"
    },
    "diels_alder_stereo": {
        "name": "Diels-Alder 反应的立体化学",
        "features": [
            "**endo 选择性**: 内型产物优先 (Alder endo rule)",
            "原因: 次级轨道相互作用(SOI)降低 endo 过渡态能量",
            "",
            "**立体专一性**: 双烯体和亲双烯体的相对构型完全保留到产物中",
            "- 顺式-反式关系: 双烯体的顺式在产物中仍为顺式",
            "- E/Z 烯烃: 亲双烯体的 E/Z 在产物中保留",
            "",
            "**对映选择性**: 手性亲双烯体或手性催化剂可诱导对映选择性"
        ],
        "example": "cyclopentadiene + maleic anhydride → endo-norbornene adduct (动力学控制)"
    },
    "chiral_synthesis": {
        "name": "不对称合成中的立体选择性",
        "methods": [
            ("手性池 (Chiral Pool)", "使用天然手性原料 (氨基酸、糖类、萜类)"),
            ("手性辅助剂 (Chiral Auxiliary)", "如 Evans oxazolidinone, 可回收"),
            ("不对称催化 (Asymmetric Catalysis)", "Noyori, Sharpless, Knowles 等诺贝尔奖工作"),
            ("酶催化 (Biocatalysis)", "酶的高对映选择性"),
            ("结晶诱导动态拆分 (CIDR)", "外消旋化+选择性结晶"),
        ],
        "key_reactions": [
            ("Sharpless 不对称环氧化", " allylic alcohol → epoxide, >90% ee"),
            ("Sharpless 不对称双羟基化", " alkene → diol, AD-mix α/β"),
            ("Jacobsen-Katsuki 环氧化", " unfunctionalized alkene → epoxide"),
            ("Noyori 不对称氢化", " C=C/C=O → chiral product, Ru-BINAP 催化"),
            ("Corey-Bakshi-Shibata (CBS) 还原", " ketone → chiral alcohol, >99% ee"),
        ]
    }
}


@ChemMCPManager.register_tool
class StereoselectivityPredictor(BaseTool):
    """
    立体选择性预测工具。
    预测有机化学反应的立体化学结果，包括 syn/anti 加成、E2/SN2 立体化学、Diels-Alder 等。
    """
    __version__ = "0.1.0"
    name = "StereoselectivityPredictor"
    func_name = "predict_stereoselectivity"
    description = "Predict stereoselectivity in organic reactions including syn/anti addition, E2/SN2 stereochemistry, Diels-Alder endo/exo selectivity, and asymmetric synthesis."
    implementation_description = "Knowledge-based tool covering stereochemical outcomes of major organic reaction classes with detailed mechanistic explanations."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Stereoselectivity", "Stereochemistry", "Syn/Anti", "E2", "SN2", "Diels-Alder", "Asymmetric Synthesis"]
    required_envs = []

    code_input_sig = [
        ("reaction_type", "str", "N/A", "Type of reaction: 'syn_anti', 'e2', 'sn2', 'hydrogenation', 'diels_alder', 'asymmetric', or a specific query."),
        ("substrate", "str", "", "Substrate description (optional, for context-specific prediction)."),
        ("detail_level", "str", "standard", "Detail level: 'brief', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated query, e.g., 'syn_anti bromination' or 'e2 cyclohexane detailed'."),
    ]

    output_sig = [
        ("result", "str", "Stereoselectivity prediction with detailed reasoning."),
    ]

    examples = [
        {
            "code_input": {"reaction_type": "syn_anti", "substrate": "bromination of alkene", "detail_level": "standard"},
            "text_input": {"query_text": "syn_anti bromination"},
            "output": {"result": "## Syn/Anti 立体化学..."}
        },
        {
            "code_input": {"reaction_type": "e2", "substrate": "cyclohexyl halide", "detail_level": "detailed"},
            "text_input": {"query_text": "e2 cyclohexyl detailed"},
            "output": {"result": "## E2 立体化学..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, substrate: str = "", detail_level: str = "standard") -> str:
        rt = reaction_type.strip().lower()
        sub = substrate.strip().lower()
        dl = detail_level.lower()

        if rt in ("syn_anti", "syn/anti", "syn anti", "addition stereo", "syn", "anti"):
            return self._predict_syn_anti(sub, dl)
        elif rt in ("e2", "elimination stereo", "e2 stereo"):
            return self._predict_e2(dl)
        elif rt in ("sn2", "substitution stereo", "walden"):
            return self._predict_sn2(dl)
        elif rt in ("hydrogenation", "h2", "catalytic hydrogenation", "lindlar", "na/nh3"):
            return self._predict_hydrogenation(dl)
        elif rt in ("diels-alder", "dielsalder", "da", "pericyclic stereo", "endo", "exo"):
            return self._predict_diels_alder(dl)
        elif rt in ("asymmetric", "chiral", "enantioselective", "ee", "不对称"):
            return self._predict_asymmetric(dl)
        elif rt in ("all", "overview", ""):
            return self._format_overview()
        else:
            # Keyword matching
            combined = f"{rt} {sub}".lower()
            if any(kw in combined for kw in ["bromin", "halogen", "br2", "cl2", "dihydroxy"]):
                return self._predict_syn_anti(sub, dl)
            if any(kw in combined for kw in ["eliminat", "e2"]):
                return self._predict_e2(dl)
            if any(kw in combined for kw in ["substitut", "sn2", "nucleophilic"]):
                return self._predict_sn2(dl)
            if any(kw in combined for kw in ["hydrogen", "h2", "lindlar", "reduction"]):
                return self._predict_hydrogenation(dl)
            return self._format_overview()

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split(None, 1)
        rtype = parts[0] if parts else ""
        rest = parts[1] if len(parts) > 1 else ""
        
        # Check for detail level
        dl = "standard"
        for token in rest.split():
            if token in ("brief", "standard", "detailed"):
                dl = token
                rest = rest.replace(token, "").strip()
                break
        
        return self._run_base(rtype, rest, dl)

    def _predict_syn_anti(self, substrate: str, dl: str) -> str:
        data = STEREOSELECTIVITY_DATA["syn_anti_addition"]
        lines = ["## Syn vs Anti 加成立体化学", ""]
        lines.append(data["description"])
        lines.append("")

        lines.append("### 各反应的加成模式")
        lines.append("")
        lines.append("| 反应 | 模式 | 原因 |")
        lines.append("|------|------|------|")
        for r in data["reactions"]:
            lines.append(f"| {r['reaction']} | {r['mode']} | {r['reason']} |")

        if dl == "detailed":
            lines.append("")
            lines.append("### 记忆技巧")
            lines.append("- **Syn (同侧)**: 催化氢化、硼氢化、OsO₄/KMnO₄ 二羟基化 — 都涉及**协同**或**表面**过程")
            lines.append("- **Anti (异侧)**: 卤素加成(X₂)、卤醇 — 都经过**三元环卤鎓离子**中间体")

        return "\n".join(lines)

    def _predict_e2(self, dl: str) -> str:
        data = STEREOSELECTIVITY_DATA["e2_stereochemistry"]
        lines = ["## E2 消除的立体化学", ""]
        lines.append(f"### 核心要求")
        lines.append(f"**{data['requirement']}**")
        lines.append(f"- 二面角: {data['dihedral_angle']}")
        lines.append("")
        
        lines.append("### 关键后果")
        for c in data["consequences"]:
            lines.append(f"- {c}")
        
        if dl != "brief":
            lines.append("")
            lines.append("### 典型示例")
            for ex in data["examples"]:
                lines.append(f"- **{ex['substrate']}**: {ex['result']}")

        if dl == "detailed":
            lines.append("")
            lines.append("### 环己烷体系详解")
            lines.append("在椅式环己烷中:")
            lines.append("- 直立键(axial): 上下方向，与相邻 axial 键呈 anti-periplanar (180°)")
            lines.append("- 平伏键(equatorial): 斜向外展，与任何键都不呈 180°")
            lines.append("")
            lines.append("→ 因此: 只有当 LG 在 axial 位时才能发生 E2 消除")
            lines.append("→ 如果 LG 在 equatorial 位: 需先 ring-flip 到 axial 构象")

        return "\n".join(lines)

    def _predict_sn2(self, dl: str) -> str:
        data = STEREOSELECTIVITY_DATA["sn2_stereochemistry"]
        lines = ["## SN2 亲核取代的立体化学", ""]
        lines.append(f"### 核心特征: **{data['configuration']}**")
        lines.append("")
        lines.append(data["mechanism"])
        lines.append("")
        lines.append(f"**示例**: {data['example']}")
        lines.append("")
        lines.append(f"> {data['note']}")

        if dl == "detailed":
            lines.append("")
            lines.append("### 立体化学证据")
            lines.append("- 如果底物是光学纯的 (R), SN2 产物应为光学纯的 (S)")
            lines.append("- 观察到的旋光度应大小相等、方向相反")
            lines.append("- 外消旋化(racemization) → 可能是 SN1 或其他机理")
            lines.append("- 部分外消旋 → 可能有 SN1/SN2 竞争或邻基参与")

        return "\n".join(lines)

    def _predict_hydrogenation(self, dl: str) -> str:
        data = STEREOSELECTIVITY_DATA["alkene hydrogenation_stereo"]
        lines = ["## 氢化反应的立体选择性", ""]
        lines.append(f"### 催化氢化: **{data['syn_addition']}**")
        lines.append("")
        for f in data["factors"]:
            lines.append(f"- {f}")
        lines.append("")
        lines.append(f"### Lindlar 催化: {data['lindlar']}")
        lines.append(f"### Na/液NH₃: {data['na_nh3']}")
        
        if dl == "detailed":
            lines.append("")
            lines.append("### 炔烃还原路径对比")
            lines.append("| 方法 | 产物构型 | 机理 |")
            lines.append("|------|---------|------|")
            lines.append("| H₂/Lindlar Pd-BaSO₄ | **Z (顺式)** | syn 加成, 停在烯烃 |")
            lines.append("| Na / liquid NH₃ | **E (反式)** | 溶解电子, 反式电子转移 |")
            lines.append("| H₂ / excess Pd-C | 烷烃 | 完全氢化 |")
            lines.append("| Birch (Na/NH₃, aromatic) | 1,4-环己二烯 | 非共轭还原 |")

        return "\n".join(lines)

    def _predict_diels_alder(self, dl: str) -> str:
        data = STEREOSELECTIVITY_DATA["diels_alder_stereo"]
        lines = ["## Diels-Alder 反应的立体化学", ""]
        for f in data["features"]:
            lines.append(f"- {f}" if not f == "" else "")
        lines.append("")
        lines.append(f"**示例**: {data['example']}")

        if dl == "detailed":
            lines.append("")
            lines.append("### Endo vs Exo 能量图")
            lines.append("- Endo 过渡态能量更低 (SOI stabilization)")
            lines.append("- 但 Exo 产物通常更稳定 (较少位阻)")
            lines.append("- 动力学控制 → endo; 热力学控制 → exo")
            lines.append("")
            lines.append("### 立体专一性实例")
            lines.append("- (E,E)-双烯体 + 亲双烯体 → 产物中取代基互为 trans")
            lines.append("- (Z,Z)-双烯体 + 亲双烯体 → 产物中取代基互为 cis")

        return "\n".join(lines)

    def _predict_asymmetric(self, dl: str) -> str:
        data = STEREOSELECTIVITY_DATA["chiral_synthesis"]
        lines = ["## 不对称合成的立体选择性", ""]
        lines.append("### 主要方法")
        for method, desc in data["methods"]:
            lines.append(f"- **{method}**: {desc}")
        lines.append("")
        lines.append("### 经典不对称反应")
        for name, result in data["key_reactions"]:
            lines.append(f"- **{name}**: {result}")

        return "\n".join(lines)

    def _format_overview(self) -> str:
        lines = ["## 立体选择性总览", ""]
        lines.append("| 反应类型 | 立体化学结果 | 决定因素 |")
        lines.append("|---------|------------|---------|")
        lines.append("| X₂ 卤化烯烃 | Anti 加成 | 卤鎓离子中间体 |")
        lines.append("| 催化氢化 | Syn 加成 | 金属表面协同转移 |")
        lines.append("| OsO₄/KMnO₄ 二羟基化 | Syn 加成 | [3+2] 环加成 |")
        lines.append("| BH₃ 硼氢化 | Syn 加成 | 四元环过渡态 |")
        lines.append("| SN2 取代 | Walden 翻转 | 背面进攻 |")
        lines.append("| E2 消除 | Anti-periplanar 要求 | 二面角 180° |")
        lines.append("| Diels-Alder | Endo 优先 + 立体专一 | SOI + 轨道对称 |")
        lines.append("| Lindlar 还原 | Z-烯烃 | Syn 加成停在烯烃 |")
        lines.append("| Na/液NH₃ 还原 | E-烯烃 | 反式电子转移 |")
        return "\n".join(lines)
