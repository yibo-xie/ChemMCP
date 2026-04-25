import logging
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 区域选择性知识库
REGIOSELECTIVITY_DATA = {
    "markovnikov_rule": {
        "name": "Markovnikov 规则 (马氏规则)",
        "statement": "在 HX 对不对称烯烃的亲电加成中，H 加在含氢较多的碳上，X 加在含氢较少（取代较多）的碳上。",
        "modern_version": "H 加成到能形成更稳定碳正离子的那个碳原子上",
        "carbocation_stability": "3° > 2° > 1° > methyl ≈ vinyl ≈ aryl",
        "examples": [
            {"reactant": "propene (CH₃-CH=CH₂) + HBr", "major": "2-bromopropane (CH₃-CHBr-CH₃)", "minor": "1-bromopropane", "ratio": ">99:1"},
            {"reactant": "2-methylpropene + HI", "major": "tert-butyl iodide (3° carbocation)", "minor": "isobutyl iodide"},
            {"reactant": "styrene (Ph-CH=CH₂) + HCl", "major": "1-chloro-1-phenylethane (benzylic cation stabilized)"},
        ],
        "exceptions": [
            {"case": "Kharasch effect / 过氧化物效应", "description": "HBr + ROOR → Anti-Markovnikov (自由基机理)", "note": "仅对 HBr 有效; HCl/HI 不受影响"},
            {"case": "F-C 烷基化重排", "description": "可能发生 1,2-H 或 1,2-烷基迁移形成更稳定的碳正离子"},
            {"case": "硼氢化反应", "description": "BH₃ 加成 → Anti-Markovnikov (syn addition, 立体因素控制)"},
        ]
    },
    "addition_to_conjugated_systems": {
        "name": "共轭体系加成的区域选择性",
        "rules": [
            "**1,2-加成 vs 1,4-加成** (α,β-不饱和羰基化合物):",
            "- 低温/动力学控制: **1,2-加成** (直接进攻羰基碳)",
            "- 高温/热力学控制: **1,4-加成** (共轭加成, Michael-type)",
            "",
            "**影响因素**: ",
            "- 硬亲核试剂(如 RLi, RMgX, LiAlH₄): 倾向 1,2-加成",
            "- 软亲核试剂(如 RS⁻, CN⁻, 胺类): 倾向 1,4-加成 (Michael addition)",
            "- CeCl₃ (Luche 还原): NaBH₄/CeCl₃ 强制 1,2-选择性"
        ],
        "examples": [
            {"substrate": "butenone (CH₂=CH-CO-CH₃) + HBr", "product": "1,2- 和 1,4- 混合物", "control": "温度决定比例"},
            {"substrate": "butenone + NH₃", "product": "1,4-加成产物 (Michael addition)", "note": "软亲核试剂"},
            {"substrate": "butenone + MeMgBr (-78°C)", "product": "1,2-加成产物为主", "note": "硬亲核试剂+低温"},
        ]
    },
    "aromatic_substitution_directing_effects": {
        "name": "芳香族取代的定位效应",
        "ortho_para_directors": [
            ("强活化", ["-NH₂", "-NHR", "-NR₂"], "孤对电子共振给电子"),
            ("中等活化", ["-OH", "-OR"], "共振给电子"),
            ("弱活化", ["-NHCOR", "-OCOR", "-R (alkyl)", "-Ph (aryl)"], "超共轭或弱共振"),
            ("弱钝化但 o/p 定位", ["-F", "--Cl", "--Br", "--I"], "诱导吸电子 > 共振给电子"),
        ],
        "meta_directors": [
            ("强钝化", ["-NO₂", "-NR₃⁺", "-CF₃", "-CCl₃"], "强吸电子诱导和共振"),
            ("中等钝化", ["-CN", "-SO₃H", "-CHO", "-COR", "-COOH", "-COOR"], "极性共振吸电子"),
        ],
        "rules": [
            "1. 活化基团加速反应; 钝化基团减慢反应",
            "2. 卤素特殊: 钝化环但 o/p 定位 (共振给电子 > 诱导吸电子)",
            "3. 多取代苯: 最强活化基团主导定位",
            "4. 位阻: o-位有大的取代基时 p-位产物增加"
        ]
    },
    "elimination_regioselectivity": {
        "name": "消除反应的区域选择性 (Zaitsev vs Hofmann)",
        "zaitsev_rule": "消除反应优先生成取代较多的烯烃（热力学更稳定）",
        "hofmann_rule": "大体积碱导致消除生成取代较少的烯烃（动力学控制，位阻最小的 β-H 更易被夺取）",
        "factors": {
            "zaitsev_favored": ["小体积碱 (NaOEt, KOH, EtONa)", "高温", "底物无特别位阻"],
            "hofmann_favored": ["大体积碱 (t-BuOK, LDA)", "季铵盐 (Hofmann elimination)", "锍盐", "β-碳有支链"],
        },
        "special_cases": [
            {"case": "Hofmann 消除 (季铵盐)", "description": "R₃N⁺-CH₂-CH₂-R' + OH⁻ → R'CH=CH₂ + R₃N + H₂O", "always_hofmann": True},
            {"case": "Cope 消除", "description": "胺氧化物热解 → 烯烃 + N,N-二甲羟胺", "stereochemistry": "syn 消除, 通常 Hofmann"},
            {"case": "Chugaev 消除", "description": "黄原酸酯热解 → 烯烃", "stereochemistry": "syn 消除, 类似 Cope"},
        ]
    },
    "other_rules": {
        "name": "其他区域选择性规则",
        "rules": [
            {
                "name": "卤醇形成 (Halohydrin formation)",
                "rule": "X₂/H₂O 对烯烃加成 → OH 加在更取代的碳上（类似马氏规则，因为卤鎓离子中间体中更取代的碳带更多 δ⁺）",
            },
            {
                "name": "炔烃水合 (Hydration of alkynes)",
                "rule": "Hg²⁺催化水合 → Markovnikov 规则适用 → 末端炔烃生成甲基酮（烯醇→酮互变异构）",
            },
            {
                "name": "环氧开环 (Epoxide ring opening)",
                "rule": "酸性条件: 在更取代的碳上开环（类似 SN1，部分碳正离子特征）; 碱性条件: 在位阻较小的碳上开环（SN2 特征）",
            },
            {
                "name": "Diels-Alder 内型/外型选择",
                "rule": "内型(endo)产物通常优势（次级轨道相互作用）— Alder endo rule",
            },
        ]
    }
}


@ChemMCPManager.register_tool
class RegioselectivityPredictor(BaseTool):
    """
    区域选择性预测工具。
    基于 Markovnikov 规则、定位效应等有机化学原理预测反应的区域选择性。
    """
    __version__ = "0.1.0"
    name = "RegioselectivityPredictor"
    func_name = "predict_regioselectivity"
    description = "Predict regioselectivity in organic reactions including Markovnikov rule, aromatic directing effects, Zaitsev/Hofmann elimination, and conjugate addition."
    implementation_description = "Knowledge-based tool implementing organic chemistry regioselectivity rules with detailed explanations and exception handling."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Regioselectivity", "Markovnikov", "Directing Effects", "Zaitsev", "Organic Chemistry"]
    required_envs = []

    code_input_sig = [
        ("reaction_type", "str", "N/A", "Type of reaction: 'markovnikov', 'anti-markovnikov', 'aromatic', 'elimination', 'conjugated', or a specific query."),
        ("substrate", "str", "", "Substrate molecule description (e.g., 'propene', 'toluene', '2-methyl-2-butanol')."),
        ("reagent", "str", "", "Reagent or condition (e.g., 'HBr', 'HNO3/H2SO4', 'NaOEt/heat')."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated query, e.g., 'markovnikov propene HBr' or 'aromatic toluene nitration'."),
    ]

    output_sig = [
        ("result", "str", "Regioselectivity prediction with detailed reasoning."),
    ]

    examples = [
        {
            "code_input": {"reaction_type": "markovnikov", "substrate": "propene", "reagent": "HBr"},
            "text_input": {"query_text": "markovnikov propene HBr"},
            "output": {"result": "## Markovnikov 预测..."}
        },
        {
            "code_input": {"reaction_type": "elimination", "substrate": "2-bromo-2-methylbutane", "reagent": "KOEt/EtOH/heat"},
            "text_input": {"query_text": "elimination 2-bromo-2-methylbutane KOEt heat"},
            "output": {"result": "## 消除区域选择性..."}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, substrate: str = "", reagent: str = "") -> str:
        rt = reaction_type.strip().lower()
        sub = substrate.strip().lower()
        rea = reagent.strip().lower()

        if rt in ("markovnikov", "markovnikoff", "马氏", "markov"):
            return self._predict_markovnikov(sub, rea)
        elif rt in ("anti-markovnikov", "anti_markovnikov", "anti-m", "反马氏", "antimarkovnikov"):
            return self._predict_anti_markovnikov(sub, rea)
        elif rt in ("aromatic", "aromatic_se", "定位", "directing", "sear", "benzene"):
            return self._predict_aromatic(sub, rea)
        elif rt in ("elimination", "zaitsev", "hofmann", "消除", "e2", "e1"):
            return self._predict_elimination(sub, rea)
        elif rt in ("conjugated", "michael", "1,2-", "1,4-", "conjugate", "共轭"):
            return self._predict_conjugated(sub, rea)
        elif rt in ("all", "overview", ""):
            return self._format_overview()
        else:
            # Try keyword matching
            combined = f"{rt} {sub} {rea}".lower()
            if any(kw in combined for kw in ["hbr", "hi", "hcl", "hx ", "addition"]):
                if "peroxide" in combined or "roor" in combined:
                    return self._predict_anti_markovnikov(sub, rea)
                return self._predict_markovnikov(sub, rea)
            if any(kw in combined for kw in ["benzene", "nitration", "sulfonation", "halogenation", "fcr"]):
                return self._predict_aromatic(sub, rea)
            if any(kw in combined for kw in ["naoet", "koh/heat", "tbuok", "eliminat"]):
                return self._predict_elimination(sub, rea)
            return self._format_overview()

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split()
        rtype = parts[0] if parts else ""
        sub = parts[1] if len(parts) > 1 else ""
        rea = " ".join(parts[2:]) if len(parts) > 2 else ""
        return self._run_base(rtype, sub, rea)

    def _predict_markovnikov(self, substrate: str, reagent: str) -> str:
        data = REGIOSELECTIVITY_DATA["markovnikov_rule"]
        lines = [f"## Markovnikov 区域选择性预测", ""]
        
        lines.append(f"**规则**: {data['statement']}")
        lines.append(f"**现代表述**: {data['modern_version']}")
        lines.append("")
        
        if substrate and reagent:
            lines.append(f"### 反应分析")
            lines.append(f"- **底物**: {substrate}")
            lines.append(f"- **试剂**: {reagent}")
            
            # Simple heuristic for common substrates
            sub_lower = substrate.lower()
            rea_lower = reagent.lower()
            
            # Check for peroxide/Kharasch
            if "peroxide" in rea_lower or "roor" in rea_lower:
                lines.append("")
                lines.append("⚠️ **检测到过氧化物 — 此反应遵循 **Anti-Markovnikov** 规则 (自由基机理)")
                lines.append("- 仅对 HBr 有效")
                return "\n".join(lines)
            
            lines.append("")
            lines.append("### 🎯 预测结果: **遵循 Markovnikov 规则**")
            lines.append(f"- H 加在含 H 较多的碳上")
            lines.append(f"- X({reagent.replace('h','').strip() or '亲电试剂'}) 加在取代较多的碳上")
            lines.append(f"- **原因**: 形成更稳定的碳正离子中间体 ({data['carbocation_stability']})")

        lines.append("")
        lines.append("### 典型示例")
        for ex in data["examples"]:
            lines.append(f"- {ex['reactant']}")
            lines.append(f"  - 主要产物: **{ex['major']}**")
            if "minor" in ex:
                lines.append(f"  - 次要产物: {ex['minor']} ({ex.get('ratio', '')})")
        lines.append("")
        lines.append("### ⚠️ 重要例外")
        for exc in data["exceptions"]:
            lines.append(f"- **{exc['case']}**: {exc['description']}")
            if "note" in exc:
                lines.append(f"  📝 {exc['note']}")
        
        return "\n".join(lines)

    def _predict_anti_markovnikov(self, substrate: str, reagent: str) -> str:
        lines = ["## Anti-Markovnikov 区域选择性预测", ""]
        lines.append("### Kharasch 效应 (过氧化物效应)")
        lines.append("")
        lines.append("**适用范围**: 仅 **HBr** + 过氧化物 (ROOR) 或 UV 光照")
        lines.append("")
        lines.append("#### 自由基机理:")
        steps = [
            "1. 引发: ROOR → 2RO· (均裂)",
            "2. RO· + HBr → ROH + Br·",
            "3. **Br· 进攻烯烃**: Br 加在含 H 较多的碳上（形成更稳定的自由基中间体）",
            "   - 自由基稳定性: 3° > 2° > 1° (与碳正离子类似)",
            "4. R· + HBr → RH + Br· (链传递)"
        ]
        for s in steps:
            lines.append(s)
        
        lines.append("")
        lines.append("### 🎯 预测结果: **Anti-Markovnikov**")
        if substrate:
            lines.append(f"- 底物: {substrate} + HBr/ROOR")
            lines.append("- Br 加在含 H 较多（取代较少）的碳上")
            lines.append("- H 加在取代较多的碳上")
        lines.append("")
        lines.append("⚠️ **重要限制**:")
        lines.append("- 仅对 **HBr** 有效 (HCl 的 H-Cl 键太强, H-I 键太弱)")
        lines.append("- HCl/HI 即使有过氧化物也仍按 Markovnikov 进行")
        lines.append("")
        lines.append("### 其他 Anti-Markovnikov 反应:")
        lines.append("- **硼氢化-氧化**: BH₃ → Anti-Markovnikov 醇 (立体因素控制)")
        lines.append("- **羰基还原(DIBAL-H)**: 酯 → 醛 (非典型区域选择性)")

        return "\n".join(lines)

    def _predict_aromatic(self, substrate: str, reagent: str) -> str:
        data = REGIOSELECTIVITY_DATA["aromatic_substitution_directing_effects"]
        lines = [f"## 芳香族取代定位效应预测", ""]
        
        if substrate:
            lines.append(f"**底物**: {substrate}")
            default_reagent = reagent or '(未指定)'
            lines.append(f"**试剂**: {default_reagent}")
            lines.append("")
            
            # Try to identify substituent on benzene
            sub_lower = substrate.lower()
            detected_group = None
            
            group_map = {
                "toluene/methyl/ch3": ("alkyl", "o/p", "弱活化"),
                "phenol/oh/hydroxy": ("oh/or", "o/p", "强活化"),
                "anisole/och3/methoxy": ("or", "o/p", "强活化"),
                "aniline/nh2/amino": ("nh2", "o/p", "最强活化"),
                "nitro/no2": ("no2", "m", "最强钝化"),
                "benzaldehyde/cho/aldehyde": ("cho", "m", "中等钝化"),
                "benzoic acid/cooh/carboxyl": ("cooh", "m", "中等钝化"),
                "acetophenone/cor/ketone": ("cor", "m", "中等钝化"),
                "benzonitrile/cn": ("cn", "m", "中等钝化"),
                "chlorobenzene/cl/chloro": ("halogen", "o/p", "弱钝化"),
                "bromobenzene/br/bromo": ("halogen", "o/p", "弱钝化"),
            }
            
            best_match = None
            for key, (group, pos, act) in group_map.items():
                if any(kw in sub_lower for kw in key.split("/")):
                    best_match = (group, pos, act)
                    break
            
            if best_match:
                group, position, activation = best_match
                lines.append(f"### 🎯 检测到取代基: **{group.upper()}** 类型")
                lines.append(f"- **定位效应**: **{position}-定位**")
                lines.append(f"- **活化/钝化**: {activation}")
                
                if position == "o/p":
                    lines.append(f"\n**预测**: 取代主要发生在 **邻位(o-)和对位(p-)**")
                    lines.append("- 邻位和对位电子云密度更高（共振给电子或超共轭）")
                    if group == "halogen":
                        lines.append("  ⚠️ 卤素特殊: 虽然钝化环但仍是 o/p 定位（共振给电子 > 诱导吸电子）")
                else:
                    lines.append(f"\n**预测**: 取代主要发生在 **间位(m-)**")
                    lines.append("- 间位相对电子云密度更高（吸电子基团使 o/p 位更缺电子）")
            else:
                lines.append("> 未识别出具体取代基。请参考以下完整定位规则:")
        else:
            lines.append("> 请指定底物以获得具体预测。\n")
        
        lines.append("")
        lines.append("### 完整定位规则速查")
        lines.append("")
        lines.append("| 取代基类型 | 代表基团 | 定位 | 活化/钝化 |")
        lines.append("|-----------|---------|------|----------|")
        for strength, groups, reason in data["ortho_para_directors"]:
            grp_str = ", ".join(g.replace("-","") for g in groups)
            lines.append(f"| {strength} | {grp_str} | **o/p** | 活化 |")
        for strength, groups, _ in data["meta_directors"]:
            grp_str = ", ".join(g.replace("-","") for g in groups)
            lines.append(f"| {strength} | {grp_str} | **m** | 钝化 |")

        return "\n".join(lines)

    def _predict_elimination(self, substrate: str, reagent: str) -> str:
        data = REGIOSELECTIVITY_DATA["elimination_regioselectivity"]
        lines = ["## 消除反应区域选择性 (Zaitsev vs Hofmann)", ""]
        
        lines.append(f"**Zaitsev 规则**: {data['zaitsev_rule']}")
        lines.append(f"**Hofmann 规则**: {data['hofmann_rule']}")
        lines.append("")
        
        if reagent:
            rea_lower = reagent.lower()
            # Determine which rule applies
            is_bulk_base = any(kw in rea_lower for kw in ["tbuok", "t-buok", "lda", "bulk", "大体积"])
            is_quaternary = "quaternary" in rea_lower or "nr3" in rea_lower or "铵" in reagent
            is_small_base = any(kw in rea_lower for kw in ["naoet", "koh", "naoh", "etona"])

            if is_bulk_base or is_quaternary:
                rule = "Hofmann"
                lines.append("### 🎯 预测: **Hofmann 产物为主** (取代较少的烯烃)")
                if is_bulk_base:
                    lines.append(f"- 原因: 大体积碱 ({reagent}) 位阻大，倾向于夺取空间位阻最小的 β-H")
                if is_quaternary:
                    lines.append("- 原因: 季铵盐的 Hofmann 消除总是给出 Hofmann 产物")
            elif is_small_base:
                rule = "Zaitsev"
                lines.append("### 🎯 预测: **Zaitsev 产物为主** (取代较多的烯烃)")
                lines.append(f"- 原因: 小体积碱 ({reagent}) 无明显位阻偏好，热力学稳定产物占优")
            else:
                rule = "unknown"
                lines.append("### ⚠️ 无法确定具体规则")
                lines.append(f"- 试剂 '{reagent}' 的信息不足以判断")
                lines.append("- 一般情况下默认 Zaitsev 产物为主")
        else:
            lines.append("### 因素判断")
            lines.append("")
            lines.append("**Zaitsef 产物 favored when:**")
            for f in data["factors"]["zaitsev_favored"]:
                lines.append(f"- {f}")
            lines.append("")
            lines.append("**Hofmann 产物 favored when:**")
            for f in data["factors"]["hofmann_favored"]:
                lines.append(f"- {f}")

        lines.append("")
        lines.append("### 特殊情况")
        for case in data.get("special_cases", []):
            lines.append(f"- **{case['case']}**: {case['description']}")
            if "stereochemistry" in case:
                lines.append(f"  立体化学: {case['stereochemistry']}")

        return "\n".join(lines)

    def _predict_conjugated(self, substrate: str, reagent: str) -> str:
        data = REGIOSELECTIVITY_DATA["addition_to_conjugated_systems"]
        lines = ["## 共轭体系加成区域选择性 (1,2- vs 1,4-)", ""]
        
        for rule in data["rules"]:
            lines.append(rule)
            lines.append("")
        
        if data.get("examples"):
            lines.append("### 示例")
            for ex in data["examples"]:
                lines.append(f"- {ex['substrate']} + {ex.get('reagent','Nu')}")
                lines.append(f"  → {ex['product']}")
                if "note" in ex:
                    lines.append(f"  📝 {ex['note']}")
        
        return "\n".join(lines)

    def _format_overview(self) -> str:
        lines = ["## 区域选择性总览", ""]
        lines.append("| 反应类型 | 控制规则 | 关键因素 |")
        lines.append("|---------|---------|---------|")
        lines.append("| 烯烃亲电加成(HX) | Markovnikov | 碳正离子稳定性 |")
        lines.append("| HBr/过氧化物 | Anti-Markovnikov | 自由基稳定性 |")
        lines.append("| 硼氢化-氧化 | Anti-Markovnikov | 立体因素 |")
        lines.append("| 芳香族取代 | 定位效应(o/p vs m) | 取代基电子效应 |")
        lines.append("| E2 消除 | Zaitsev/Hofmann | 碱体积 + 温度 |")
        lines.append("| α,β-不饱和羰基加成 | 1,2- vs 1,4- | 亲核试剂硬度 + 温度 |")
        return "\n".join(lines)
