import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 还原反应机理知识库
REDUCTION_MECHANISM_DATA = {
    "nabh4_reduction": {
        "name": "NaBH₄ (Sodium Borohydride) Reduction",
        "reagents": "NaBH₄ in MeOH or EtOH",
        "substrate": "醛 → 伯醇 / 酮 → 仲醇 / 酰氯 →醛 / 环氧化物开环",
        "mechanism_steps": [
            "1. BH₄⁻ 作为氢负离子供体（亲核试剂），羰基碳作为亲电中心",
            "2. BH₄⁻ 的 H⁻ 进攻羰基碳，π电子转移到氧上形成烷氧基硼中间体",
            "3. 溶剂(ROH)质子化烷氧基，释放产物醇和 BH₃（进一步反应）",
            "4. 总计每分子 NaBH₄ 可还原 4 分子羰基化合物"
        ],
        "reactivity_order": "醛 > 酮 > 环氧化物 > 酰氯 >> 酯/酰胺 (不反应)",
        "features": ["温和条件 (0°C ~ rt)", "选择性还原 C=O", "在质子溶剂中稳定", "安全易操作"],
        "limitations": ["不能还原酯、酰胺、羧酸", "对水敏感但比 LiAlH₄ 好得多", "通常不还原孤立 C=C"],
        "chemoselectivity": [
            "可还原醛酮而不影响酯、酰胺、硝基、卤素",
            "不还原孤立碳碳双键（C=C）",
            "可还原 α,β-不饱和羰基化合物中的羰基（1,2-还原为主）或共轭加成（1,4-还原，需 CeCl₃ 等 Lewis 酸催化）",
            "Luche 还原: NaBH₄/CeCl₃ 选择性 1,2-还原 α,β-不饱和羰基"
        ]
    },
    "lialh4_reduction": {
        "name": "LiAlH₄ (Lithium Aluminium Hydride) Reduction",
        "reagents": "LiAlH₄ in dry Et₂O or THF",
        "substrate": "醛→醇 / 酮→醇 / 酯→伯醇 / 酰胺→胺 / 酸→醇 / 酰氯→醇 / 环氧化物→醇 / 腈→伯胺",
        "mechanism_steps": [
            "1. AlH₄⁻ 中 H⁻ 进攻缺电子羰基碳（强亲核性）",
            "2. 形成烷氧基铝氢化物中间体",
            "3. 进一步的 H⁻ 转移可能发生（取决于底物类型）",
            "4. 水解步骤（后处理）：加入 H₂O 或稀酸分解 Al-O 键，释放醇",
            "5. 水解剧烈放热！必须小心操作"
        ],
        "reactivity_order": "酰氯 > 醛 > 酮 > 酯 ≈ 羧酸 > 酰胺 > 腈",
        "features": ["最强常用的还原剂", "几乎还原所有含羰基官能团", "还原腈为伯胺"],
        "limitations": ["极度水/空气敏感", "无水无氧操作", "水解剧烈放热", "无化学选择性（会还原几乎所有羰基）", "不能还原孤立 C=C"],
        "safety_notes": "与水反应剧烈放出 H₂ 气体；必须在惰性气氛中操作；使用时需严格干燥"
    },
    "li_in_nh3_reduction": {
        "name": "Li / liquid NH₃ (Dissolving Metal Reduction)",
        "reagents": "Li metal in liquid NH₃ (-33°C), 通常加 t-BuOH 作质子源",
        "substrate": "炔烃 → 反式烯烃 / 芳环 → 1,4-环己二烯(Birch还原) / α,β-不饱和酮 → 饱和酮",
        "mechanism_steps": [
            "1. Li 在液氨中溶解生成溶剂化电子 e⁻(solv)",
            "2. 底物接受单电子形成自由基阴离子",
            "3. 自由基阴离子从质子源(t-BuOH/NH₃)夺取质子",
            "4. 再接受一个电子形成阴离子",
            "5. 第二次质子化得到最终产物",
            "6. 对于炔烃: 反式加成（anti addition）→ 反式烯烃"
        ],
        "features": ["Birch 还原: 芳香族 → 非共轭二烯", "炔烃 → 反式烯烃（高立体选择性）", "共轭还原（1,4-还原）"],
        "limitations": ["需要液氨（-33°C）", "操作复杂", "某些官能团不兼容"]
    },
    "catalytic_hydrogenation": {
        "name": "Catalytic Hydrogenation (Pd, Pt, Ni)",
        "reagents": "H₂ + Pd/C, PtO₂ (Adams catalyst), or Raney Ni",
        "substrate": "烯烃/炔烃 → 烷烃 / 硝基 → 胺 / 腈 → 胺 / 醛/酮 → 醇 / 苄基脱保护",
        "mechanism_steps": [
            "1. H₂ 在金属表面解离吸附为两个 M-H 键（Heterolytic 或 homolytic）",
            "2. 不饱和键配位到金属表面（π-络合物形成）",
            "3. 逐步转移 H 原子到底物上（顺式加成 syn addition）",
            "4. 饱和产物从金属表面脱附",
            "5. **关键**: 顺式加成（syn addition）— 两个 H 从同侧加上去"
        ],
        "selectivity_rules": [
            "Pd/BaSO₄ (Lindlar催化剂): 炔烃 → 顺式烯烃（停住）",
            "Na/液NH₃: 炔烃 → 反式烯烃",
            "Pd/C: 芐基醚/酯容易氢解（脱保护常用）",
            "PtO₂: 还原活性最高"
        ],
        "features": ["原子经济性好（H₂）", "操作简单", "可放大生产"],
        "limitations": ["需要压力设备", "催化剂昂贵（Pd, Pt）", "可能过度还原", "中毒失活（S, P 化合物）"]
    },
    "other_reductions": {
        "name": "Other Important Reduction Methods",
        "methods": [
            {
                "name": "DIBAL-H (Diisobutylaluminium hydride)",
                "scope": "酯 → 醛 (低温控制); 腈 → 亚胺; 环氧化物 → 醇",
                "note": "体积大导致位阻控制，可在 -78°C 将酯停在醛阶段"
            },
            {
                "name": "B₂H₆ / BH₃·THF (Hydroboration)",
                "scope": "烯烃 → 硼烷 → 氧化 → 醇 (反马氏)",
                "note": "实际是硼氢化-氧化两步，但常归在还原类"
            },
            {
                "name": "Wolff-Kishner / Huang-Minlon",
                "scope": "醛/酮 → 亚甲基 (C=O → CH₂)",
                "note": "碱性条件 (KOH/肼/乙二醇)，与 Clemmensen 互补"
            },
            {
                "name": "Clemmensen Reduction",
                "scope": "醛/酮 → 亚甲基 (C=O → CH₂)",
                "note": "酸性条件 (Zn(Hg)/HCl)，酸敏感底物用 Wolff-Kishner"
            },
            {
                "name": "McFadyen-Stevens",
                "scope": "酰肼 → 醛",
                "note": "温和条件下将羧酸衍生物转化为醛"
            }
        ]
    }
}


@ChemMCPManager.register_tool
class ReductionMechanism(BaseTool):
    """
    各类还原反应机理分析工具。
    支持 NaBH₄、LiAlH₄、Li/NH₃、催化氢化等常见还原反应的详细机理分析。
    """
    __version__ = "0.1.0"
    name = "ReductionMechanism"
    func_name = "reduction_mechanism_analysis"
    description = "Analyze and explain reduction reaction mechanisms including NaBH4, LiAlH4, dissolving metal reduction, and catalytic hydrogenation."
    implementation_description = "Knowledge-based tool with built-in mechanism database for common organic reduction reactions."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Reduction", "Mechanism", "NaBH4", "LiAlH4", "Hydrogenation", "Organic Chemistry"]
    required_envs = []

    code_input_sig = [
        ("reaction_name", "str", "N/A", "Name of the reduction reaction: 'nabh4', 'lialh4', 'li_nh3', 'hydrogenation', 'dibal-h', 'wolff-kishner', 'clemmensen', or 'all' for overview."),
        ("detail_level", "str", "standard", "Detail level: 'brief', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated reaction name and detail level, e.g., 'nabh4 detailed' or 'all'."),
    ]

    output_sig = [
        ("result", "str", "Detailed analysis of the reduction reaction mechanism."),
    ]

    examples = [
        {
            "code_input": {"reaction_name": "nabh4", "detail_level": "detailed"},
            "text_input": {"query_text": "nabh4 detailed"},
            "output": {"result": "## NaBH4 Reduction..."}
        },
        {
            "code_input": {"reaction_name": "all", "detail_level": "brief"},
            "text_input": {"query_text": "all brief"},
            "output": {"result": "..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_name: str, detail_level: str = "standard") -> str:
        r = reaction_name.strip().lower().replace("-", "_").replace(" ", "_")
        dl = detail_level.lower()

        mapping = {
            "nabh4": "nabh4_reduction",
            "nabhd4": "nabh4_reduction",
            "sodium_borohydride": "nabh4_reduction",
            "lialh4": "lialh4_reduction",
            "lithium_aluminium_hydride": "lialh4_reduction",
            "lithium_aluminum_hydride": "lialh4_reduction",
            "laih4": "lialh4_reduction",
            "li_nh3": "li_in_nh3_reduction",
            "linh3": "li_in_nh3_reduction",
            "dissolving_metal": "li_in_nh3_reduction",
            "birch": "li_in_nh3_reduction",
            "hydrogenation": "catalytic_hydrogenation",
            "catalytic": "catalytic_hydrogenation",
            "pdc": "catalytic_hydrogenation",
            "pd/c": "catalytic_hydrogenation",
            "dibal": "other_reductions",
            "dibal-h": "other_reductions",
            "wolff-kishner": "other_reductions",
            "wolffkishner": "other_reductions",
            "huang-minlon": "other_reductions",
            "clemmensen": "other_reductions",
            "all": None,
            "overview": None,
            "list": None,
        }

        key = mapping.get(r)
        if key is None:
            if r in mapping.values():
                key = r
            else:
                return self._format_overview(dl)

        if key:
            return self._format_single(key, dl)
        return self._format_overview(dl)

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split(None, 1)
        name = parts[0]
        detail = parts[1] if len(parts) > 1 else "standard"
        return self._run_base(name, detail)

    def _format_overview(self, detail_level: str) -> str:
        lines = ["## 有机还原反应机理总览", ""]
        lines.append("| 还原剂 | 强度 | 适用范围 | 条件 | 选择性 |")
        lines.append("|--------|------|---------|------|--------|")
        items = [
            ("NaBH₄", "中等", "醛→醇, 酮→醇", "MeOH, 0°C~rt", "高（不影响酯/酸）"),
            ("LiAlH₄", "极强", "所有羰基化合物", "THF, 无水", "低（全还原）"),
            ("Li/NH₃", "强", "炔烃→反式烯烃; Birch还原", "液NH₃, -33°C", "专一"),
            ("H₂/Pd-C", "中等偏强", "C=C, C≡C, NO₂, CN等", "常压/加压", "取决于催化剂"),
            ("DIBAL-H", "强", "酯→醛(-78°C)", "甲苯, -78°C", "温度控制"),
            ("Wolff-Kishner", "特殊", "C=O → CH₂", "KOH/肼/Δ", "碱性条件专用"),
            ("Clemmensen", "特殊", "C=O → CH₂", "Zn(Hg)/HCl", "酸性条件专用"),
        ]
        for item in items:
            lines.append(f"| {item[0]} | {item[1]} | {item[2]} | {item[3]} | {item[4]} |")

        if detail_level != "brief":
            lines.append("")
            lines.append("### 选择指南")
            lines.append("- **只需还原醛酮**: NaBH₄（首选，安全简单）")
            lines.append("- **需要还原酯/酸/酰胺**: LiAlH₄")
            lines.append("- **炔烃→反式烯烃**: Li/液NH₃")
            lines.append("- **炔烃→顺式烯烃**: Lindlar Pd/BaSO₄")
            lines.append("- **羰基变亚甲基**: 酸性用 Clemmensen，碱性用 Wolff-Kishner")
            lines.append("- **酯停在醛**: DIBAL-H, -78°C")
        return "\n".join(lines)

    def _format_single(self, key: str, detail_level: str) -> str:
        data = REDUCTION_MECHANISM_DATA[key]
        lines = [f"## {data['name']} 反应机理", ""]
        lines.append(f"**还原剂**: {data.get('reagents', 'N/A')}")
        lines.append(f"**适用范围**: {data.get('substrate', 'N/A')}")
        lines.append("")

        if detail_level != "brief":
            lines.append("### 分步机理")
            for step in data["mechanism_steps"]:
                lines.append(f"- {step}")
            lines.append("")

        if "reactivity_order" in data:
            lines.append(f"**活性顺序**: {data['reactivity_order']}")
            lines.append("")

        lines.append("### 特点")
        for feat in data["features"]:
            lines.append(f"- ✅ {feat}")

        lines.append("")
        lines.append("### 局限性")
        for lim in data["limitations"]:
            lines.append(f"- ⚠️ {lim}")

        extra_fields = ["chemoselectivity", "selectivity_rules", "safety_notes"]
        for field in extra_fields:
            if field in data and detail_level == "detailed":
                lines.append(f"\n### {field.replace('_', ' ').title()}")
                if isinstance(data[field], list):
                    for item in data[field]:
                        lines.append(f"- {item}")
                else:
                    lines.append(str(data[field]))

        if "methods" in data:
            lines.append("\n### 其他重要还原方法")
            for m in data["methods"]:
                lines.append(f"\n**{m['name']}**")
                lines.append(f"- 适用: {m['scope']}")
                lines.append(f"- 备注: {m['note']}")

        return "\n".join(lines)
