import logging
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 亲核性知识库
NUCLEOPHILICITY_DATA = {
    "ranking": [
        # General ranking in protic solvent (polar protic, e.g., H2O, ROH)
        (1, "RS⁻ / thiolate", "very strong", "大原子 + 可极化 + 弱碱"),
        (2, "CN⁻ (cyanide)", "very strong", "碳亲核试剂，sp 杂化，轨道暴露"),
        (3, "I⁻ (iodide)", "very strong", "最大可极化性，极好的 Nu 但弱碱"),
        (4, "SH⁻ (hydrosulfide)", "very strong", "类似 RS⁻"),
        (5, "Br⁻ (bromide)", "strong", "良好可极化性"),
        (6, "HO⁻ (hydroxide)", "strong", "小体积强碱，在质子溶剂中活性中等偏上"),
        (7, "CH₃O⁻ (methoxide)", "strong", "强碱但比 HO⁻ 体积大"),
        (8, "Cl⁻ (chloride)", "moderate", "中等可极化性"),
        (9, "CH₃COO⁻ (acetate)", "weak-moderate", "共振稳定化降低亲核性"),
        (10, "F⁻ (fluoride)", "weak in protic / strong in aprotic", "高电负性+强溶剂化→质子溶剂中差; 非质子溶剂中极佳"),
        (11, "H₂O (water)", "weak", "中性分子，亲核性弱"),
        (12, "ROH (alcohol)", "weak", "中性分子"),
        (13, "CH₃COOH (acetic acid)", "very weak", "几乎不表现亲核性"),
        (14, "NH₃ (ammonia)", "moderate", "中性分子但有孤对电子"),
        (15, "RNH₂ (amine)", "moderate-good", "有机胺，常用亲核试剂"),
        (16, "N₃⁻ (azide)", "good", "线性、可极化，也是好的离去基团"),
        (17, "Ph₃P (triphenylphosphine)", "special", "用于 Mitsunobu, Wittig 等特殊反应"),
        (18, "enolate ion", "variable", "取决于结构; 通常为 strong"),
        (19, "organometallic (RLi, RMgX)", "very strong (hard)", "硬亲核试剂，与羰基反应极佳"),
    ],
    "principles": [
        "**基本矛盾**: 亲核性 ≠ 碱性",
        "- 碱性: 对 H⁺ 的亲和力（热力学性质）",
        "- 亲核性: 对缺电子碳(或其他原子)的进攻能力（动力学性质）",
        "",
        "**影响亲核性的关键因素**:",
        "",
        "1. **电荷**: 负离子 > 中性分子 (RO⁻ > ROH > ROH₂⁺)",
        "",
        "2. **可极化性(Polarizability)**:",
        "   - 周期表中同族: 从上到下 ↑ 可极化性 → ↑ 亲核性",
        "   - I⁻ > Br⁻ > Cl⁻ > F⁻ (在**质子溶剂**中)",
        "   - RSe⁻ > RS⁻ > RO⁻ (同族趋势)",
        "",
        "3. **溶剂效应 (关键!)**:",
        "   - **极性质子溶剂** (H₂O, ROH): 通过氢键强烈溶剂化小阴离子",
        "     → 大的/可极化的阴离子被较少溶剂化 → 更活泼",
        "     → 排序: I⁻ > Br⁻ > Cl⁻ > F⁻ (与碱性相反!)",
        "   - **极性非质子溶剂** (DMSO, DMF, acetone, acetonitrile, HMPA):",
        "     → 只溶剂化阳离子，不溶剂化阴离子（或很弱）",
        "     → 小而强的碱（如 F⁻）变得非常活泼（\"裸露\"阴离子）",
        "     → 排序: F⁻ > Cl⁻ > Br⁻ > I⁻ (与碱性一致)",
        "",
        "4. **立体位阻**:",
        "   - 位阻越小 → 亲核性越强 (SN2 需要)",
        "   - CH₃O⁻ > t-BuO⁻ (虽然 t-BuO⁻ 碱性更强)",
        "   - 这就是为什么 t-BuK 是消除试剂而非取代试剂",
        "",
        "5. **轨道特性 (HSAB)**:",
        "   - 硬亲核试剂 (small, low polarizability): HO⁻, RO⁻, NH₃, RLi, RMgX",
        "     → 偏好进攻硬亲电中心 (羰基碳, H⁺)",
        "   - 软亲核试剂 (large, high polarizability): I⁻, RS⁻, Ph₃P, CN⁻",
        "     → 偏好进攻软亲电中心 (软碳如饱和碳, Pd 等)",
        "",
        "6. **共振效应**:",
        "   - 电荷离域 → 亲核性下降",
        "   - CH₃COO⁻ < HO⁻ (乙酸根共振稳定化)",
        "   - enolate 的 O-端 vs C-端: C-端更亲核（硬/软差异）",
    ],
    "solvent_comparison": {
        "protic_solvent": {
            "examples": "H₂O, ROH (MeOH, EtOH), carboxylic acids",
            "halide_ranking": "I⁻ > Br⁻ > Cl⁻ >> F⁻",
            "reason": "小阴离子被强氢键溶剂化，难以脱离溶剂笼进攻底物",
            "note": "F⁻ 在水中几乎无亲核性（形成稳定的 HF₂⁻ 或强氢键络合物）"
        },
        "aprotic_solvent": {
            "examples": "DMSO, DMF, acetone, MeCN, HMPA, THF (弱)",
            "halide_ranking": "F⁻ > Cl⁻ > Br⁻ > I⁻",
            "reason": "阳离子被溶剂化但阴离子基本 \"裸露\"，碱性/亲核性一致",
            "note": "这就是为什么 F⁻ 在 DMSO 中是优秀的 SN2 试剂（如 TBAF 用于脱硅基）"
        }
    },
    "special_cases": {
        "ambident_nucleophiles": [
            ("CN⁻", "C-端进攻 (主要) vs N-端进攻", "C-端较软，对饱和碳更有利; N-端对羰基有利"),
            ("NO₂⁻", "N-端 vs O-端", "取决于溶剂和底物"),
            ("Enolate", "C-端 (烷基化) vs O-端 (酰化)", "Hard electrophile → O; Soft electrophile → C (Irreversible-Trombe Model)"),
        ],
        "intramolecular": ["分子内亲核性远大于分子间 (邻近效应/Chelate effect / 有效浓度)"],
        "alpha_effect": ["具有邻位孤对电子的亲核试剂活性异常高 (e.g., NH₂OH > NH₃; NH₂NH₂ > NH₃) — 原因: 邻基孤对电子通过 n→σ* 相互作用稳定过渡态"],
    },
    "nucleophile_vs_base": {
        "title": "何时是亲核试剂? 何时是碱?",
        "rules": [
            "高温 → 有利于消除(E2/E1) — 消除有更高的活化熵 (ΔS‡ > 0)",
            "位阻大的底物 → 有利于消除 (E2)",
            "位阻大的碱 (t-BuOK, LDA) → 有利于消除 (无法 SN2 进攻)",
            "强碱 + 伯底物 → SN2 与 E2 竞争",
            "弱碱 + 叔底物 → SN1/E1",
            "亲核性强但碱性弱的试剂 (I⁻, RS⁻, CN⁻) → 有利于取代",
            "可极化性高的试剂 → 有利于 SN2 (soft-soft 匹配)"
        ]
    }
}


@ChemMCPManager.register_tool
class NucleophilicityRanker(BaseTool):
    """
    亲核试剂亲核性比较工具。
    排序和比较常见亲核试剂的亲核性，分析溶剂效应和影响因素。
    """
    __version__ = "0.1.0"
    name = "NucleophilicityRanker"
    func_name = "rank_nucleophiles"
    description = "Compare and rank the nucleophilicity of common nucleophiles, analyze solvent effects, and distinguish between nucleophilicity and basicity."
    implementation_description = "Knowledge-based tool with comprehensive nucleophilicity database including rankings in different solvents, HSAB theory, ambident nucleophiles, and nucleophile-vs-base selectivity rules."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Nucleophilicity", "Basicity", "Solvent Effects", "SN2", "Organic Chemistry", "HSAB"]
    required_envs = []

    code_input_sig = [
        ("query", "str", "ranking", "Query type: 'ranking' (full ranked list), 'compare' (specific groups), 'solvent' (solvent effects), 'vs_base' (Nu vs base), or 'principles'."),
        ("nucleophiles", "str", "", "Comma-separated list of nucleophiles to compare (for 'compare' mode). e.g., 'I, Br, Cl, F, OH, OMe, SH, CN'."),
        ("detail_level", "str", "standard", "Detail level: 'brief', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated query, e.g., 'ranking', 'compare I Br Cl F OH', 'solvent detailed', or 'vs_base'."),
    ]

    output_sig = [
        ("result", "str", "Ranked comparison of nucleophiles with explanations."),
    ]

    examples = [
        {
            "code_input": {"query": "ranking", "nucleophiles": "", "detail_level": "standard"},
            "text_input": {"query_text": "ranking"},
            "output": {"result": "## Nucleophilicity Ranking..."}
        },
        {
            "code_input": {"query": "compare", "nucleophiles": "I, Br, Cl, F, OH, OCH3, SH, CN", "detail_level": "standard"},
            "text_input": {"query_text": "compare I Br Cl F OH OCH3 SH CN"},
            "output": {"result": "## Comparison..."}
        },
        {
            "code_input": {"query": "solvent", "nucleophiles": "", "detail_level": "detailed"},
            "text_input": {"query_text": "solvent detailed"},
            "output": {"result": "## Solvent Effects..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, query: str = "ranking", nucleophiles: str = "", detail_level: str = "standard") -> str:
        q = query.strip().lower()
        nucs = nucleophiles.strip()
        dl = detail_level.lower()

        if q in ("rank", "ranking", "list"):
            return self._format_ranking(dl)
        elif q in ("compare", "comparison", "vs"):
            return self._compare_nucleophiles(nucs, dl)
        elif q in ("solvent", "solvent effect", "solvent_effects"):
            return self._format_solvent_effects(dl)
        elif q in ("vs_base", "vsbase", "nu_vs_base", "nucleophile_vs_base"):
            return self._format_nu_vs_base(dl)
        elif q in ("principle", "principles", "theory", "basis"):
            return self._format_principles(dl)
        elif q in ("ambident", "special", "alpha_effect"):
            return self._format_special_cases(dl)
        else:
            if q and not nucs:
                return self._compare_nucleophiles(q, dl)
            return self._format_ranking(dl)

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split(None, 1)
        q = parts[0] if parts else "ranking"
        rest = parts[1] if len(parts) > 1 else ""
        
        dl = "standard"
        for token in rest.split():
            if token in ("brief", "standard", "detailed"):
                dl = token
                rest = rest.replace(token, "").strip()
                break
        
        return self._run_base(q, rest, dl)

    def _format_ranking(self, dl: str) -> str:
        lines = ["## 亲核试剂亲核性排名 (极性质子溶剂中)", ""]
        
        if dl != "brief":
            lines.append("**默认条件**: 极性质子溶剂 (如 H₂O, EtOH)")
            lines.append("> ⚠️ 注意: 在极性非质子溶剂中排序可能完全不同!")
            lines.append("")

        lines.append("| 排名 | 亲核试剂 | 强度 | 特点 |")
        lines.append("|------|---------|------|------|")

        for rank, nu, strength, notes in NUCLEOPHILICITY_DATA["ranking"]:
            emoji_map = {
                "very strong": "🔥🔥", "strong": "🔥", "moderate-good": "✅",
                "moderate": "🟢", "weak-moderate": "🟡", "weak": "🟠",
                "very weak": "🔴", "variable": "🔄", "special": "⭐"
            }
            lines.append(f"| {rank} | {nu} | {emoji_map.get(strength, '')} {strength} | {notes} |")

        if dl == "detailed":
            lines.append("")
            lines.append("### 快速记忆口诀")
            lines.append("**质子溶剂中**: RS⁻ ≈ CN⁻ > I⁻ > Br⁻ > Cl⁻ > HO⁻ > F⁻")
            lines.append("**非质子溶剂中**: F⁻ > Cl⁻ > Br⁻ > I⁻ (顺序翻转!)")

        return "\n".join(lines)

    def _compare_nucleophiles(self, nucs_str: str, dl: str) -> str:
        if not nucs_str:
            return self._format_ranking(dl)

        requested = [n.strip() for n in nucs_str.replace(",", " ").split()]
        requested = [n for n in requested if n]

        # Build lookup
        lookup = {}
        for rank, nu, strength, notes in NUCLEOPHILICITY_DATA["ranking"]:
            nu_lower = nu.lower()
            key = nu_lower.split("/")[0].split(" ")[0].replace("⁻","").replace("-","").replace("(","").replace(")","").replace("_","")
            lookup[key] = (rank, nu, strength, notes)
            # Also add short aliases
            for alias in [nu_lower.replace("⁻",""), nu_lower.replace("-","")]:
                simple = alias.split("(")[0].replace(" ","").replace("/","").lower()
                if simple and simple not in lookup:
                    lookup[simple] = (rank, nu, strength, notes)

        alias_map = {
            "i": "iodide", "br": "bromide", "cl": "chloride", "f": "fluoride",
            "oh": "hydroxide", "ome": "methoxide", "oet": "ethoxide", "otbu": "tert-butoxide",
            "sh": "hydrosulfide", "rs": "thiolate", "cn": "cyanide", "n3": "azide",
            "oac": "acetate", "h2o": "water", "roh": "alcohol", "nh3": "ammonia",
            "rnh2": "amine", "pPh3": "triphenylphosphine", "rli": "organolithium",
            "rmgx": "grignard", "enolate": "enolate ion",
        }

        lines = [f"## 亲核试剂比较: {nucs_str}", ""]
        results = []
        not_found = []

        for n in requested:
            key = n.replace("⁻","").replace("-","").replace(" ","").lower()
            match = lookup.get(key) or lookup.get(alias_map.get(key, ""))
            if match:
                results.append(match)
            else:
                not_found.append(n)

        results.sort(key=lambda x: x[0])

        if results:
            lines.append("| 亲核试剂 | 排名 | 强度 | 特点 |")
            lines.append("|---------|------|------|------|")
            for rank, nu, strength, notes in results:
                lines.append(f"| {nu} | #{rank} | {strength} | {notes} |")
            lines.append("")
            lines.append(f"**结论 (质子溶剂中)**: {' > '.join(nu for _, nu, _, _ in results)}")

        if not_found:
            lines.append(f"\n⚠️ 未找到: {', '.join(not_found)}")
            lines.append("> 常用亲核试剂包括: I⁻, Br⁻, Cl⁻, F⁻, ⁻OH, CH₃O⁻, t-BuO⁻, RS⁻, CN⁻, N₃⁻, CH₃COO⁻, H₂O, NH₃, RNH₂, Ph₃P, enolate, RLi, RMgX")

        return "\n".join(lines)

    def _format_solvent_effects(self, dl: str) -> str:
        data = NUCLEOPHILICITY_DATA["solvent_comparison"]
        lines = ["## 溶剂效应对亲核性的影响", ""]
        lines.append("### ⚠️ 这是理解亲核性的最关键概念之一!")
        lines.append("")

        for stype, sdata in data.items():
            title = "极性质子溶剂" if "protic" in stype else "极性非质子溶剂"
            lines.append(f"### {title}")
            lines.append(f"- **示例**: {sdata['examples']}")
            lines.append(f"- **卤素离子排序**: **{sdata['halide_ranking']}**")
            lines.append(f"- **原因**: {sdata['reason']}")
            lines.append(f"- 📝 {sdata['note']}")
            lines.append("")

        if dl == "detailed":
            lines.append("### 为什么会有这种差异?")
            lines.append("")
            lines.append("**质子溶剂中的氢键网络:**")
            lines.append("- 小阴离子(F⁻, Cl⁻)电荷密度高 → 形成强氢键 → 被 \"包裹\" → 难以进攻")
            lines.append("- 大阴离子(I⁻, Br⁻)电荷密度低 → 弱氢键 → 较自由 → 容易进攻")
            lines.append("")
            lines.append("**非质子溶剂的情况:**")
            lines.append("- 阳离子(Na⁺, K⁺)被深度溶剂化（氧原子配位）")
            lines.append("- 阴离子基本不受溶剂束缚 → \"裸露\"状态")
            lines.append("- 此时碱性越强 = 亲核性越强 (F⁻ > Cl⁻ > Br⁻ > I⁻)")

        return "\n".join(lines)

    def _format_nu_vs_base(self, dl: str) -> str:
        data = NUCLEOPHILICITY_DATA["nucleophile_vs_base"]
        lines = [f"## {data['title']}", ""]
        for rule in data["rules"]:
            lines.append(f"- {rule}")
        lines.append("")
        lines.append("### 经验总结")
        lines.append("| 试剂类型 | 倾向 | 典型代表 |")
        lines.append("|---------|------|---------|")
        lines.append("| 强碱 + 小体积 | SN2 > E2 | HO⁻, CH₃O⁻, CN⁻ |")
        lines.append("| 强碱 + 大体积 | E2 >> SN2 | t-BuOK, LDA |")
        lines.append("| 弱碱 + 高可极化性 | SN2 >> E1/E2 | I⁻, RS⁻, Ph₃P, CN⁻ |")
        lines.append("| 中性分子 | 取代为主 (需催化) | H₂O, ROH, RNH₂ |")
        lines.append("| 可极化阴离子 | SN2 优异 | I⁻, Br⁻, RS⁻, SeR⁻ |")

        return "\n".join(lines)

    def _format_principles(self, dl: str) -> str:
        lines = ["## 亲核性的理论基础", ""]
        for p in NUCLEOPHILICITY_DATA["principles"]:
            lines.append(p)
        return "\n".join(lines)

    def _format_special_cases(self, dl: str) -> str:
        data = NUCLEOPHILICITY_DATA["special_cases"]
        lines = ["## 特殊情况: 两可亲核剂与其他效应", ""]
        
        lines.append("### 两可亲核试剂 (Ambident Nucleophiles)")
        for nu, desc, note in data.get("ambident_nucleophiles", []):
            lines.append(f"- **{nu}**: {desc}")
            lines.append(f"  {note}")
        
        lines.append("")
        for k, v in data.items():
            if k != "ambident_nucleophiles":
                lines.append(f"### {k.replace('_', ' ').title()}")
                if isinstance(v, list):
                    for item in v:
                        lines.append(f"- {item}")

        return "\n".join(lines)
