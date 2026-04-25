import logging
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 离去基团知识库
LEAVING_GROUP_DATA = {
    "ranking": [
        # (rank, leaving_group, pKa_of_conjugate_acid, type, notes)
        (1, "N₂ (diazonium)", -10, "excellent", "极好的离去基团，自发离去生成稳定 N₂ 气体"),
        (2, "TsO⁻ / tosylate (p-toluenesulfonate)", -2, "excellent", "常用最佳离去基团之一，共振稳定化"),
        (3, "MsO⁻ / mesylate (methanesulfonate)", -1, "excellent", "类似 tosylate 但体积更小"),
        (4, "TfO⁻ / triflate (trifluoromethanesulfonate)", -14, "exceptional", "最强非超氧化物离去基团，几乎任何 SN2/SN1"),
        (5, "I⁻ (iodide)", -7, "very good", "最大原子半径，电荷分散好；也是好的亲核试剂"),
        (6, "Br⁻ (bromide)", -9, "good", "常用离去基团"),
        (7, "Cl⁻ (chloride)", -7, "moderate", "中等离去能力"),
        (8, "H₂O (water from -OH₂⁺)", -1.7, "moderate", "-OH 本身很差，质子化后变好"),
        (9, "N₃⁻ (azide)", 4.6, "moderate", "本身是好的亲核试剂，但也可作为离去基团"),
        (10, "CH₃COO⁻ (acetate)", 4.76, "poor", "弱酸共轭碱，离去能力差"),
        (11, "F⁻ (fluoride)", 3.2, "very poor", "HF 键强，F⁻ 碱性强，极难离去（但在某些特殊情况下可）"),
        (12, "-OH (hydroxide)", 15.7, "terrible", "极差的离去基团，强碱，通常需质子化或转化"),
        (13, "-OR (alkoxide)", 16-18, "terrible", "比 -OH 更差，强碱"),
        (14, "-NH₂ (amide)", 38, "extremely terrible", "几乎不能作为离去基团"),
        (15, "-CR₃⁻ (carbanion)", 40, "worst", "极强的碱，绝不可能作为离去基团"),
        (16, "NO₃⁻ (nitrate)", -1.4, "good", "硝酸根，较好的离去基团"),
        (17, "POCl₃-derived (from -OH)", 0, "good", "醇经 POCl₃/Py 转化为好的离去基团"),
        (18, "RSO₂- (sulfonate esters general)", 0, "excellent", "磺酸酯类都是优秀的离去基团"),
    ],
    "principles": [
        "**基本规则**: 离去基团能力 ∝ 共轭酸的酸性 (pKa 越小 → 离去越好)",
        "**好的离去基团**: 弱碱 (共轭酸 pKa < 5)，能稳定负电荷",
        "**差的离去基团**: 强碱 (共轭酸 pKa > 10)，难以承受额外电子对",
        "**关键因素**:",
        "- 电荷稳定性: 共振、诱导效应、极化性",
        "- 溶剂化: 大的、可极化的离子在质子溶剂中更稳定",
        "- 键强度: H-LG 键越弱 → LG 越容易离去",
    ],
    "activation_methods": {
        "-OH (alcohol)": [
            ("质子化", "用 HX (HCl, HBr, HI) 或 H₂SO₄ 将 -OH 质子化为 -OH₂⁺"),
            ("转化为磺酸酯", "TsCl/Py 或 MsCl/Py → OTs 或 OMs (极佳 LG)"),
            ("转化为无机酯", "SOCl₂, PCl₃, PCl₅, POCl₃ → Cl 为 LG"),
            ("Appel 反应", "CBr₄/PPh₃ → Br 为 LG (温和条件)"),
            ("Mitsunobu 反应", "DEAD/PPh₃ → 构型翻转的取代"),
        ],
        "-NH₂ (amine)": [
            ("转化为重氮盐", "NaNO₂/HCl (0-5°C) → N₂⁺ 作为 LG (Sandmeyer)"),
            ("转化为铵盐", "彻底甲基化后消除 (Hofmann/Cope)"),
        ],
        "-COOH (carboxylic acid)": [
            ("质子化/活化", "转化为酰氯 (SOCl₂), 酸酐, 酯等"),
            ("脱羧反应", "特定条件下 -COOH 以 CO₂ 形式离去"),
        ]
    },
    "context_dependence": {
        "SN2_vs_SN1": {
            "SN2": "LG 能力影响速率但不改变机理; 好 LG 加速 SN2",
            "SN1": "LG 离去是决速步; LG 能力对 SN1 至关重要 — 差的 LG 基本不发生 SN1"
        },
        "E1_E2": "与取代类似: E1 需要 LG 先离去形成碳正离子; E2 中 LG 与 β-H 同步离去",
        "solvent_effects": "极性质子溶剂(水, 醇)通过氢键稳定 LG⁻ → 增强 LG 能力; 极性非质子溶剂(DMSO, DMF) 不稳定化 LG⁻ → 相对减弱 LG 效应但增强 Nu 活性"
    }
}


@ChemMCPManager.register_tool
class LeavingGroupRanker(BaseTool):
    """
    离去基团离去能力比较工具。
    排序和比较常见离去基团的离去能力，并提供活化方法。
    """
    __version__ = "0.1.0"
    name = "LeavingGroupRanker"
    func_name = "rank_leaving_groups"
    description = "Compare and rank the leaving group ability of common groups in organic substitution and elimination reactions."
    implementation_description = "Knowledge-based tool with comprehensive leaving group database including ranking by ability, pKa correlation, activation methods for poor leaving groups, and context-dependent behavior."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Leaving Group", "Substitution", "Elimination", "Organic Chemistry", "pKa"]
    required_envs = []

    code_input_sig = [
        ("query", "str", "ranking", "Query type: 'ranking' (full ranked list), 'compare' (specific groups), 'activate' (how to activate poor LG), or 'principles'."),
        ("groups", "str", "", "Comma-separated list of leaving groups to compare (for 'compare' mode). e.g., 'I, Br, Cl, OH, OTs'."),
        ("detail_level", "str", "standard", "Detail level: 'brief', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated query, e.g., 'ranking', 'compare I Br Cl OH', 'activate OH', or 'principles detailed'."),
    ]

    output_sig = [
        ("result", "str", "Ranked comparison of leaving groups with explanations."),
    ]

    examples = [
        {
            "code_input": {"query": "ranking", "groups": "", "detail_level": "standard"},
            "text_input": {"query_text": "ranking"},
            "output": {"result": "## Leaving Group Ranking..."}
        },
        {
            "code_input": {"query": "compare", "groups": "I, Br, Cl, F, OH, OTs", "detail_level": "standard"},
            "text_input": {"query_text": "compare I Br Cl F OH OTs"},
            "output": {"result": "## Comparison..."}
        },
        {
            "code_input": {"query": "activate", "groups": "OH", "detail_level": "standard"},
            "text_input": {"query_text": "activate OH"},
            "output": {"result": "## Activating -OH..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, query: str = "ranking", groups: str = "", detail_level: str = "standard") -> str:
        q = query.strip().lower()
        gs = groups.strip()
        dl = detail_level.lower()

        if q in ("rank", "ranking", "list"):
            return self._format_ranking(dl)
        elif q in ("compare", "comparison", "vs"):
            return self._compare_groups(gs, dl)
        elif q in ("activate", "activation", "how to"):
            return self._activation_methods(gs, dl)
        elif q in ("principle", "principles", "theory", "basis"):
            return self._format_principles(dl)
        elif q in ("context", "solvent", "sn1", "sn2"):
            return self._format_context(dl)
        else:
            # Try to interpret as a compare query
            if q and not gs:
                return self._compare_groups(q, dl)
            return self._format_ranking(dl)

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split(None, 1)
        q = parts[0] if parts else "ranking"
        rest = parts[1] if len(parts) > 1 else ""
        
        # Check for detail level
        dl = "standard"
        for token in rest.split():
            if token in ("brief", "standard", "detailed"):
                dl = token
                rest = rest.replace(token, "").strip()
                break
        
        return self._run_base(q, rest, dl)

    def _format_ranking(self, dl: str) -> str:
        lines = ["## 离去基团离去能力排名", ""]
        
        if dl != "brief":
            lines.append("**核心原则**: 离去基团能力 ≈ 其共轭酸的酸性 (pKa 越小越好)")
            lines.append("")
        
        lines.append("| 排名 | 离去基团 | 共轭酸 pKa | 等级 | 备注 |")
        lines.append("|------|---------|-----------|------|------|")

        for rank, lg, pka, grade, notes in LEAVING_GROUP_DATA["ranking"]:
            # Clean up pKa display
            if isinstance(pka, int) and pka < -10:
                pka_str = "~ -14" if pka == -14 else str(pka)
            elif pka == "-":
                pka_str = "—"
            else:
                pka_str = str(pka)
            
            emoji = {"exceptional": "🌟", "excellent": "✅✅", "very good": "✅", 
                    "good": "🟢", "moderate": "🟡", "poor": "🟠", 
                    "very poor": "🔴", "terrible": "❌", "extremely terrible": "💀", "worst": "☠️"}
            lines.append(f"| {rank} | {lg} | {pka_str} | {emoji.get(grade, '')} {grade} | {notes} |")

        if dl == "detailed":
            lines.append("")
            lines.append("### 快速记忆")
            lines.append("- **最佳 LG**: TfO⁻ > TsO⁻ > MsO⁻ > I⁻ > Br⁻ > Cl⁻ >> F⁻ > ⁻OH > ⁻OR > ⁻NH₂")
            lines.append("- **记住**: 磺酸酯 > 卤素 > 含氧阴离子 (除质子化的)")
            lines.append("- **最差**: ⁻NH₂ 和碳负离子几乎不可能作为离去基团")

        return "\n".join(lines)

    def _compare_groups(self, groups_str: str, dl: str) -> str:
        if not groups_str:
            return self._format_ranking(dl)

        # Parse groups
        requested = [g.strip().lower().replace(",", "") for g in groups_str.replace(",", " ").split()]
        requested = [g for g in requested if g]

        # Build lookup from our data
        lookup = {}
        for rank, lg, pka, grade, notes in LEAVING_GROUP_DATA["ranking"]:
            lg_lower = lg.lower().split("/")[0].split(" ")[0].replace("⁻","").replace("-","").replace("(","").replace(")","")
            lookup[lg_lower] = (rank, lg, pka, grade, notes)
        # Also add aliases
        alias_map = {
            "tosylate": "tso", "mesylate": "mso", "triflate": "tfo", "azide": "n3",
            "acetate": "oac", "water": "h2o", "hydroxide": "oh", "alkoxide": "or",
            "amide": "nh2", "carbanion": "cr3", "nitrate": "no3",
            "iodide": "i", "bromide": "br", "chloride": "cl", "fluoride": "f",
            "ts": "tso", "ms": "mso", "tf": "tfo",
        }

        lines = [f"## 离去基团比较: {groups_str}", ""]
        results = []
        not_found = []

        for g in requested:
            key = g.replace("⁻","").replace("-","").replace(" ","").lower()
            match = lookup.get(key) or lookup.get(alias_map.get(key, ""))
            if match:
                results.append(match)
            else:
                not_found.append(g)

        # Sort by rank
        results.sort(key=lambda x: x[0])

        if results:
            lines.append("| 离去基团 | 排名 | 等级 | pKa | 备注 |")
            lines.append("|---------|------|------|-----|------|")
            for rank, lg, pka, grade, notes in results:
                lines.append(f"| {lg} | #{rank} | {grade} | {pka} | {notes} |")
            lines.append("")
            lines.append(f"**结论**: {' > '.join(lg for _, lg, _, _, _ in results)} (从好到差排序)")

        if not_found:
            lines.append(f"\n⚠️ 未找到: {', '.join(not_found)}")
            lines.append("> 可用的离去基团包括: I⁻, Br⁻, Cl⁻, F⁻, ⁻OH, OTs⁻, OMs⁻, OTf⁻, N₃⁻, CH₃COO⁻, H₂O, NO₃⁻ 等")

        return "\n".join(lines)

    def _activation_methods(self, target: str, dl: str) -> str:
        data = LEAVING_GROUP_DATA["activation_methods"]
        target_clean = target.strip().lower().replace("-","").replace(" ","") if target else ""

        lines = ["## 差离去基团的活化方法", ""]
        
        if target_clean:
            # Find matching activation methods
            found = False
            for key, methods in data.items():
                key_clean = key.lower().replace("-","").replace(" ","")
                if target_clean in key_clean or key_clean in target_clean:
                    found = True
                    lines.append(f"### 活化 {key}")
                    for i, (method, desc) in enumerate(methods, 1):
                        lines.append(f"{i}. **{method}**: {desc}")
                    break
            
            if not found:
                lines.append(f"> 未找到 '{target}' 的具体活化方法。以下是常见情况:")
                for key, methods in data.items():
                    lines.append(f"\n**{key}**:")
                    for method, desc in methods:
                        lines.append(f"- {method}: {desc}")
        else:
            for key, methods in data.items():
                lines.append(f"### {key}")
                for method, desc in methods:
                    lines.append(f"- **{method}**: {desc}")
                lines.append("")

        return "\n".join(lines)

    def _format_principles(self, dl: str) -> str:
        lines = ["## 离去基团能力的理论基础", ""]
        for p in LEAVING_GROUP_DATA["principles"]:
            lines.append(p)
            lines.append("")
        return "\n".join(lines)

    def _format_context(self, dl: str) -> str:
        data = LEAVING_GROUP_DATA.get("context_dependent", {})
        lines = ["## 离去基团在不同反应环境中的行为", ""]
        for ctx, desc in data.items():
            lines.append(f"### {ctx}")
            lines.append(desc)
            lines.append("")
        return "\n".join(lines)
