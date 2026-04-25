import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 氧化反应机理知识库
OXIDATION_MECHANISM_DATA = {
    "swern_oxidation": {
        "name": "Swern Oxidation",
        "reagents": "(COCl)₂ (草酰氯), DMSO, Et₃N",
        "substrate": "伯醇/仲醇 → 醛/酮",
        "mechanism_steps": [
            "1. 活化: DMSO 与草酰氯在 -78°C 下反应生成活性氯代锍盐中间体 [ClCH₂SMe]⁺[OCOCOCl]⁻",
            "2. 烷氧基锍盐形成: 醇进攻活化的 DMSO，脱去 HCl 形成烷氧基锍盐中间体",
            "3. 脱质子: Et₃N 夺取 α-氢，通过五元环过渡态消除生成羰基化合物和二甲硫醚(CH₃SCH₃)",
            "4. 副产物: CO, CO₂, CH₃SCH₃ (气体逸出驱动反应)"
        ],
        "features": ["温和条件 (-78°C)", "不过度氧化（伯醇→醛）", "无重金属", "副产物为气体易分离"],
        "limitations": ["低温操作", "产生恶臭二甲硫醚", "草酰剧毒"]
    },
    "pcc_oxidation": {
        "name": "PCC (Pyridinium Chlorochromate) Oxidation",
        "reagents": "PCC in CH₂Cl₂",
        "substrate": "伯醇 → 醛 / 仲醇 → 酮",
        "mechanism_steps": [
            "1. 铬酸酯形成: 醇羟基氧原子配位到 Cr(VI) 中心，取代一个氯离子",
            "2. 形成 chromate ester 中间体（铬酸酯）",
            "3. β-消除: 碱（常为吡啶或另一个醇分子）夺取 α-H，Cr(VI) 同时被还原为 Cr(IV)",
            "4. Cr(IV) 歧化为 Cr(III) 和 Cr(V)，最终 Cr(III) 为稳定产物"
        ],
        "features": ["中性条件", "伯醇停在醛阶段", "在无水 CH₂Cl₂ 中进行"],
        "limitations": ["Cr(VI) 有毒致癌", "环境不友好", "酸性条件下可能过度氧化"]
    },
    "jones_oxidation": {
        "name": "Jones Oxidation",
        "reagents": "CrO₃ / H₂SO₄ / 丙酮",
        "substrate": "伯醇 → 羧酸 / 仲醇 → 酮",
        "mechanism_steps": [
            "1. 在酸性介质中形成铬酸 H₂CrO₄ 或 HCrO₄⁻",
            "2. 醇与铬酸形成 chromate ester",
            "3. 类似 PCC 的消除机理，但水相条件允许醛水合物进一步氧化为羧酸"
        ],
        "features": ["反应迅速", "便宜", "氧化能力强"],
        "limitations": ["伯醇过度氧化至羧酸", "酸敏感官能团不兼容", "Cr(VI) 剧毒"]
    },
    "dess_martin_oxidation": {
        "name": "Dess-Martin Periodinane (DMP) Oxidation",
        "reagents": "Dess-Martin 试剂 (C₁₃H₁₃IO₈) in CH₂Cl₂",
        "substrate": "伯醇 → 醛 / 仲醇 → 酮",
        "mechanism_steps": [
            "1. 醇对碘(V) 中心进行亲核取代，置换一个乙酸根",
            "2. 形成 alkoxyperiodinane 中间体",
            "3. 碱性消除: 通过四元环过渡态，乙酸根夺取 α-H，同时 I(V) 还原为 I(III)",
            "4. 生成羰基化合物和碘(III)副产物"
        ],
        "features": ["极其温和（室温）", "高化学选择性", "官能团耐受性好", "无毒"],
        "limitations": ["试剂较贵", "对潮湿敏感", "可能发生过度氧化（α-氧化）"]
    },
    "tpap_oxidation": {
        "name": "TPAP/NMO (Ley-Griffith) Oxidation",
        "reagents": "TPAP (tetrapropylammonium perruthenate), NMO (氧化剂)",
        "substrate": "伯醇 → 醛 / 仲醇 → 酮",
        "mechanism_steps": [
            "1. Ru(VII) 与醇形成 Ru-烷氧化物",
            "2. β-消除: 通过氢化物转移机理，Ru(VII) 还原为 Ru(V)",
            "3. NMO 将 Ru(V) 重新氧化为 Ru(VII)（催化循环）"
        ],
        "features": ["催化量 Ru (1-5 mol%)", "室温快速反应", "高选择性"],
        "limitations": ["TPAP 价格高且可能爆炸性", "需要化学计量的共氧化剂"]
    },
    "collins_oxidation": {
        "name": "Collins Oxidation",
        "reagents": "CrO₃·pyridine complex (Collins 试剂) in CH₂Cl₂",
        "substrate": "伯醇 → 醛 / 仲醇 → 酮",
        "mechanism_steps": [
            "1. 吡啶-CrO₃ 配合物中的 Cr(VI) 与醇形成 chromate ester",
            "2. 无水条件下的 β-消除（类似 PCC 但吡啶作为碱和配体）",
            "3. 生成醛/酮和 Cr(IV) 物种"
        ],
        "features": ["无水条件，伯醇不停在醛", "比 Jones 温和"],
        "limitations": ["试剂需新鲜制备", "吸湿性强", "Cr 毒性"]
    },
    "activated_dmsO_oxidations": {
        "name": "Other Activated DMSO Oxidations",
        "variants": [
            {"name": "Parikh-Doering", "reagent": "SO₃·Py + DMSO, Et₃N", "note": "Swern 的变体"},
            {"name": "Albright-Onodera", "reagent": "TFAA + DMSO", "note": "用 TFAA 替代草酰氯"},
            {"name": "Moffatt", "reagent": "DCC + DMSO", "note": "最早发现的 DMSO 氧化法"}
        ]
    }
}


@ChemMCPManager.register_tool
class OxidationMechanism(BaseTool):
    """
    各类氧化反应机理分析工具。
    支持 Swern、PCC、Jones、Dess-Martin、TPAP、Collins 等常见氧化反应的详细机理分析。
    """
    __version__ = "0.1.0"
    name = "OxidationMechanism"
    func_name = "oxidation_mechanism_analysis"
    description = "Analyze and explain oxidation reaction mechanisms including Swern, PCC, Jones, Dess-Martin, TPAP, and Collins oxidations."
    implementation_description = "Knowledge-based tool with built-in mechanism database for common organic oxidation reactions."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Oxidation", "Mechanism", "Swern", "PCC", "Dess-Martin", "Organic Chemistry"]
    required_envs = []

    code_input_sig = [
        ("reaction_name", "str", "N/A", "Name of the oxidation reaction: 'swern', 'pcc', 'jones', 'dess_martin', 'tpap', 'collins', or 'all' for overview."),
        ("detail_level", "str", "standard", "Detail level: 'brief', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated reaction name and detail level, e.g., 'swern detailed' or 'all'."),
    ]

    output_sig = [
        ("result", "str", "Detailed analysis of the oxidation reaction mechanism."),
    ]

    examples = [
        {
            "code_input": {"reaction_name": "swern", "detail_level": "detailed"},
            "text_input": {"query_text": "swern detailed"},
            "output": {"result": "## Swern Oxidation..."}
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
        """Core logic for oxidation mechanism analysis."""
        r = reaction_name.strip().lower().replace("-", "_").replace(" ", "_")
        dl = detail_level.lower()

        mapping = {
            "swern": "swern_oxidation",
            "pcc": "pcc_oxidation",
            "jones": "jones_oxidation",
            "dess-martin": "dess_martin_oxidation",
            "dessmartin": "dess_martin_oxidation",
            "dmp": "dess_martin_oxidation",
            "tpap": "tpap_oxidation",
            "ley-griffith": "tpap_oxidation",
            "leygriffith": "tpap_oxidation",
            "collins": "collins_oxidation",
            "parikh-doering": "activated_dmso_oxidations",
            "parikhdoering": "activated_dmso_oxidations",
            "moffatt": "activated_dmso_oxidations",
            "dmso": "activated_dmso_oxidations",
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
        lines = ["## 有机氧化反应机理总览", ""]
        lines.append("| 反应名称 | 氧化剂 | 底物范围 | 温度 | 特点 |")
        lines.append("|---------|--------|---------|------|------|")
        summary_items = [
            ("Swern", "(COCl)₂/DMSO/Et₃N", "1°→醛, 2°→酮", "-78°C", "温和, 无金属"),
            ("PCC", "CrO₃·py·HCl", "1°→醛, 2°→酮", "rt", "中性条件"),
            ("Jones", "CrO₃/H₂SO₄", "1°→酸, 2°→酮", "0°C→rt", "强氧化"),
            ("Dess-Martin", "Periodinane", "1°→醛, 2°→酮", "rt", "温和, 无毒"),
            ("TPAP/NMO", "Ru(VII)/NMO", "1°→醛, 2°→酮", "rt", "催化量"),
            ("Collins", "CrO₃·py₂", "1°→醛, 2°→酮", "rt", "无水"),
        ]
        for item in summary_items:
            lines.append(f"| {item[0]} | {item[1]} | {item[2]} | {item[3]} | {item[4]} |")

        if detail_level != "brief":
            lines.append("")
            lines.append("### 选择指南")
            lines.append("- **需要温和无金属**: Swern 或 Dess-Martin")
            lines.append("- **大量底物/成本敏感**: Jones (注意安全)")
            lines.append("- **酸敏感底物**: Dess-Martin > Swern > TPAP")
            lines.append("- **避免 Cr 毒性**: Dess-Martin, TPAP, 或 Swern")
        return "\n".join(lines)

    def _format_single(self, key: str, detail_level: str) -> str:
        data = OXIDATION_MECHANISM_DATA[key]
        name = data["name"]

        lines = [f"## {name} 反应机理", ""]
        lines.append(f"**氧化剂**: {data['reagents']}")
        lines.append(f"**适用范围**: {data['substrate']}")
        lines.append("")

        if detail_level != "brief":
            lines.append("### 分步机理")
            for step in data["mechanism_steps"]:
                lines.append(f"- {step}")
            lines.append("")

        lines.append("### 特点")
        for feat in data["features"]:
            lines.append(f"- ✅ {feat}")

        lines.append("")
        lines.append("### 局限性")
        for lim in data["limitations"]:
            lines.append(f"- ⚠️ {lim}")

        if "variants" in data and detail_level == "detailed":
            lines.append("")
            lines.append("### 相关变体")
            for v in data["variants"]:
                lines.append(f"- **{v['name']}**: {v['reagent']} — {v['note']}")

        return "\n".join(lines)
