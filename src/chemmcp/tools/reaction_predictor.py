import logging
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 反应预测知识库 — 基于有机化学规则
REACTION_PREDICTION_DATA = {
    "nucleophilic_substitution": {
        "name": "亲核取代反应",
        "patterns": [
            {"substrate": "R-X (alkyl halide) + Nu⁻", "product": "R-Nu + X⁻", "conditions": "SN1(3°/allylic/benzylic + polar protic), SN2(1°/methyl + polar aprotic)"},
            {"substrate": "R-OTs (tosylate) + Nu⁻", "product": "R-Nu", "note": "OTs 是极好的离去基团，类似 SN2"},
            {"substrate": "R-OH + HX", "product": "R-X + H₂O", "note": "需质子化 -OH 变成更好的离去基团 (-OH₂⁺)"},
        ],
        "sn1_vs_sn2": {
            "sn2_favored": ["甲基卤代烷", "伯卤代烷", "无位阻的仲卤代烷", "强亲核试剂", "极性非质子溶剂(DMF, DMSO, acetone)", "低温"],
            "sn1_favored": ["叔卤代烷", "烯丙型/苄型卤代烷", "弱亲核试剂/中性亲核试剂", "极性质子溶剂(H₂O, ROH)", "高温", "底物能形成稳定碳正离子"],
        }
    },
    "elimination": {
        "name": "消除反应",
        "patterns": [
            {"substrate": "R-CH₂-CH₂-X + strong base", "product": "R-CH=CH₂ + HX", "type": "E2 (anti-periplanar)"},
            {"substrate": "R₃C-X + weak base/heat", "product": "R₂C=CHR + HX", "type": "E1 (via carbocation)"},
            {"substrate": "R-CH(OH)-CH₃ + Δ", "product": "R-CH=CH₂ + H₂O", "type": "E1 (acid-catalyzed dehydration)"},
        ],
        "zaitsev_vs_hofmann": {
            "zaitsef": ["小体积强碱 (NaOEt, KOH)", "高温", "产物: 取代多的烯烃（更稳定）"],
            "hofmann": ["大体积强碱 (t-BuOK, LDA)", "产物: 取代少的烯烃（位阻控制）"],
        }
    },
    "addition_to_carbonyl": {
        "name": "羰基加成",
        "patterns": [
            {"reagent": "Grignard (RMgX) + aldehyde/ketone", "product": "醇 (1° from aldehyde, 2° from ketone)", "note": "不能与酯中的C=O停止（继续加成）"},
            {"reagent": "Organolithium (RLi) + aldehyde/ketone", "product": "醇", "note": "比 Grignard 更活泼"},
            {"reagent": "NaBH₄ / LiAlH₄ + aldehyde", "product": "伯醇"},
            {"reagent": "NaBH₄ / LiAlH₄ + ketone", "product": "仲醇"},
            {"reagent": "HCN + aldehyde/ketone", "product": "氰醇 (cyanohydrin)"},
            {"reagent": "hydrazine + ketone", "product": "腙 (hydrazone)", "note": "Wolff-Kishner 第一步"},
            {"reagent": "primary amine + aldehyde/ketone", "product": "亚胺 (imine/Schiff base)"},
            {"reagent": "secondary amine + aldehyde/ketone", "product": "烯胺 (enamine)"},
            {"reagent": "H₂O/H₃O⁺ + aldehyde/ketone (hydrate formation)", "product": "geminal diol (水合物的平衡偏向羰基)"},
            {"reagent": "ROH/H⁺ + aldehyde/ketone (acetal formation)", "product": "缩醛 (aldehyde) / 缩酮 (ketone)"},
        ]
    },
    "addition_to_alkene": {
        "name": "烯烃加成",
        "patterns": [
            {"reagent": "H-X (HBr, HI, HCl)", "product": "Markovnikov 卤代烷", "mechanism": "electrophilic addition via carbocation"},
            {"reagent": "HBr + peroxides (ROOR)", "product": "Anti-Markovnikov 溴代烷", "mechanism": "radical chain (Kharasch effect)"},
            {"reagent": "X₂ (Br₂, Cl₂) in inert solvent", "product": "邻二卤化物 (vicinal dihalide)", "mechanism": "halonium ion intermediate → anti addition"},
            {"reagent": "HO-X (X₂ in H₂O)", "product": "卤代醇 (halohydrin)", "regioselectivity": "OH 加在更取代的碳上（更稳定碳正离子特征）"},
            {"reagent": "H₂O/H₂SO₄ (acid-catalyzed hydration)", "product": "醇 (Markovnikov)", "note": "可逆反应"},
            {"reagent": "Hg(OAc)₂ / H₂O then NaBH₄ (oxymercuration-demercuration)", "product": "醇 (Markovnikov, no rearrangement)", "advantage": "避免碳正离子重排"},
            {"reagent": "BH₃·THF then H₂O₂/NaOH (hydroboration-oxidation)", "product": "醇 (Anti-Markovnikov, syn addition)", "note": "B 加在取代少的一端"},
            {"reagent": "H₂ / metal catalyst (Pd, Pt, Ni)", "product": "烷烃 (syn addition)", "note": "催化氢化"},
            {"reagent": "OsO₄ / NMO or KMnO₄ (cold, dilute)", "product": "顺式邻二醇 (cis-vicinal diol)", "mechanism": "[3+2] 环加成"},
            {"reagent": "KMnO₄ (hot, conc.) or O₃ then Zn/H₂O", "product": "羰基化合物 (cleavage)", "note": "臭氧氧化后还原水解得醛/酮; 强氧化得羧酸/酮"},
            {"reagent": "peroxyacid (mCPBA)", "product": "环氧化物 (epoxide)", "stereochemistry": "立体专一性保留"},
        ]
    },
    "addition_to_alkyne": {
        "name": "炔烃加成",
        "patterns": [
            {"reagent": "2 eq. H-X", "product": "偕二卤代烷 (geminal dihalide)", "note": "两步 Markovnikov 加成"},
            {"reagent": "1 eq. H-X", "product": "乙烯基卤化物 (vinyl halide)"},
            {"reagent": "H₂ / Lindlar Pd-BaSO₄", "product": "顺式烯烃 (Z-alkene)", "note": "停在烯烃阶段"},
            {"reagent": "H₂ / Na / liquid NH₃", "product": "反式烯烃 (E-alkene)", "note": "溶解金属还原"},
            {"reagent": "H₂O / HgSO₄ / H₂SO₄", "product": "烯醇互变异构为醛(Markovnikov) 或 酮", "note": "末端炔烃→甲基酮"},
            {"reagent": "2 eq. Br₂", "product": "四卤代烷"},
            {"reagent": "1 eq. Br₂", "product": "(E)-二溴烯烃 (trans-dibromoalkene)"},
        ]
    },
    "oxidation_reactions": {
        "name": "氧化反应",
        "patterns": [
            {"substrate": "1° alcohol", "reagent": "PCC / Swern / Dess-Martin", "product": "醛 (停在醛)"},
            {"substrate": "1° alcohol", "reagent": "Jones (CrO₃/H₂SO₄) / KMnO₄", "product": "羧酸"},
            {"substrate": "2° alcohol", "reagent": "任何常用氧化剂", "product": "酮"},
            {"substrate": "3° alcohol", "reagent": "-", "product": "不反应（无α-H）"},
            {"substrate": "aldehyde", "reagent": "Tollens' [Ag(NH₃)₂]⁺", "product": "羧酸 + Ag镜(定性检测)"},
            {"substrate": "aldehyde/ketone", "reagent": "Fehling's / Benedict's", "product": "醛阳性(红色Cu₂O); 酮阴性(除α-羟基酮)"},
            {"substrate": "alkyl side-chain of aromatics", "reagent": "KMnO₄ / Na₂Cr₂O₇ / HNO₃", "product": "苯甲酸衍生物 (只要有 α-H 的侧链都能被氧化为-COOH)"},
        ]
    },
    "aromatic_electrophilic_substitution": {
        "name": "芳香族亲电取代 (SEAr)",
        "patterns": [
            {"reagent": "X₂ / FeX₃ or AlX₃", "product": "卤代芳烃 (halogenation)", "director": "o/p directing but deactivating"},
            {"reagent": "HNO₃ / H₂SO₄ (mixed acid)", "product": "硝基化合物 (nitration)", "director": "m-directing, strongly deactivating"},
            {"reagent": "R-Cl / AlCl₃ (Friedel-Crafts alkylation)", "product": "烷基苯", "issues": "易重排、多烷基化、无效定位基失活环"},
            {"reagent": "RCOCl / AlCl₃ (Friedel-Crafts acylation)", "product": "芳基酮 (aryl ketone)", "advantages": "不重排、单酰基化、产物是间位定位基"},
            {"reagent": "ROR / H⁺ (酸催化) or SO₃ / H₂SO₄", "product": "烷氧基苯 (alkylation) 或 磺酸 (sulfonation)", "note": "磺化为可逆反应，可用于占位策略"},
        ],
        "directing_effects": {
            "ortho/para_activating": ["-NH₂, -NHR, -NR₂ (最强)", "-OH, -OR", "-NHCOR", "- alkyl (-R)", "-Ph (苯基)", "-X (卤素, 弱钝化但 o/p 定位)"],
            "meta_deactivating": ["-NO₂ (最强)", "-CN", "-SO₃H", "-CHO, -COR", "-COOH, -COOR", "-NR₃⁺"],
        }
    }
}


@ChemMCPManager.register_tool
class ReactionPredictor(BaseTool):
    """
    反应产物预测工具。
    基于有机化学反应规则和机理知识，给定反应物和试剂预测主要产物。
    """
    __version__ = "0.1.0"
    name = "ReactionPredictor"
    func_name = "predict_reaction_product"
    description = "Predict the major organic reaction product(s) given reactants and reagents, based on established organic chemistry rules and mechanisms."
    implementation_description = "Knowledge-based reaction prediction engine covering nucleophilic substitution, elimination, carbonyl addition, alkene/alkyne addition, oxidation, and aromatic electrophilic substitution."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Reaction Prediction", "Organic Chemistry", "Product Prediction", "Mechanism-Based"]
    required_envs = []

    code_input_sig = [
        ("reactants", "str", "N/A", "Reactant(s): e.g., '2-methyl-2-butanol', '1-propene', 'benzaldehyde', 'cyclohexanone'."),
        ("reagents", "str", "N/A", "Reagent(s)/conditions: e.g., 'HBr', 'NaOEt/heat', 'PCC', 'H2/Pd-C', 'excess NH3'."),
        ("detail_level", "str", "standard", "Detail level: 'brief' (product only), 'standard' (product+reasoning), or 'detailed' (full mechanism)."),
    ]

    text_input_sig = [
        ("query_text", "str", "N/A", "Space-separated reactants and reagents, e.g., 'cyclohexene Br2 in H2O standard'."),
    ]

    output_sig = [
        ("result", "str", "Predicted product with reasoning and mechanism notes."),
    ]

    examples = [
        {
            "code_input": {"reactants": "1-methylcyclohexanol", "reagents": "conc. H2SO4 / heat", "detail_level": "standard"},
            "text_input": {"query_text": "1-methylcyclohexanol conc.H2SO4 heat standard"},
            "output": {"result": "## 预测结果: methylenecyclohexane (major)..."}
        },
        {
            "code_input": {"reactants": "styrene", "reagents": "HBr", "detail_level": "brief"},
            "text_input": {"query_text": "styrene HBr brief"},
            "output": {"result": "1-bromo-1-phenylethane (Markovnikov product)"}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reactants: str, reagents: str, detail_level: str = "standard") -> str:
        r = reactants.strip().lower()
        rea = reagents.strip().lower()
        dl = detail_level.lower()

        # Classify the reaction type based on keywords
        rtype = self._classify_reaction(r, rea)

        if rtype:
            return self._predict_by_type(r, rea, rtype, dl)
        
        return self._general_prediction(r, rea, dl)

    def _run_text(self, query_text: str) -> str:
        parts = query_text.strip().split()
        # Try to find detail level at end
        dl = "standard"
        if parts and parts[-1] in ("brief", "standard", "detailed"):
            dl = parts.pop()

        if len(parts) >= 2:
            # Find where reagents start (look for common reagent keywords)
            reagent_keywords = ["hbr", "hi", "hcl", "h2so4", "naoh", "koh", "naoet", "tbuok", 
                               "lda", "pcc", "nacn", "nh3", "br2", "cl2", "nabh4", "lialh4",
                               "h2", "pd", "pt", "ni", "kmno4", "o3", "mcpba", "bh3", "hg",
                               "conc.", "heat", "Δ", "roor", "peroxides", "hv", "light",
                               "alc.", "alcohol", "aq.", "water", "cold", "dilute", "hot",
                               "conc", "dil", "toluene", "thf", "dmso", "ether", "acid", "base",
                               "tscl", "pyridine", "socl2", "px3", "px5", "sox2", "dess-martin",
                               "swern", "jones", "tollens", "febr3", "alcl3", "hno3", "fcr",
                               "friedel-crafts", "acylation", "alkylation", "grignard", "rmgx",
                               "rli", "organolithium", "h2o2", "naoh", "nmo", "os04", "lindlar",
                               "na/nh3", "nano2", "hcl", "cu cn", "cun", "sn1", "sn2", "e1", "e2"]

            # Simple heuristic: last part before dl is likely reagent
            # For now just split: first word = reactant, rest = reagents
            reactants = parts[0]
            reagents_str = " ".join(parts[1:]) if len(parts) > 1 else ""
            return self._run_base(reactants, reagents_str, dl)

        return self._run_base(query_text, "", dl)

    def _classify_reaction(self, reactants: str, reagents: str) -> Optional[str]:
        """Classify reaction type from reactant and reagent descriptions."""
        r_combined = f"{reactants} {reagents}".lower()

        # Substitution patterns
        sub_kw = ["halide", "bromide", "chloride", "iodide", "tosylate", "mesylate", "alkyl ",
                  "-x", "r-x", "substitut", " sn1", " sn2", " naoh", " nacn", " nan3", " sh",
                  " ch3coona", " acetate", " i-", " br-", " cl-"]
        elim_kw = ["eliminat", "e1 ", "e2 ", "naoet", "tbuok", " koh/heat", "naoh/heat", 
                   "alc.koh", "alcoholic", "dehydrat", " h2so4/heat", "conc.h2so4", "pocl3"]
        add_carbonyl_kw = ["aldehyde", "ketone", "carbonyl", "c=o", "formaldehyde", "acetone",
                          "cyclohexanone", "benzaldehyde", "acetaldehyde", "butanone",
                          "grignard", "rmgx", "rli", "nabh4", "lialh4", "hcn", "hydrazine",
                          "nh2oh", "amine", "imine", "enamine", "acetal", "hemiacetal", "hydrate"]
        add_alkene_kw = ["alkene", "alkene", "ene ", "= ", "ethylene", "propene", "butene",
                        "styrene", "cyclohexene", "cyclopentene", "octene", "hexene",
                        "hbr", "hi", "hcl", "br2", "cl2", "h2o/h+", "hydration", "hydroboration",
                        "oxymercuration", "halogenation", "halohydrin", "hydrogenation", "h2/",
                        "kmno4", "o3", "ozonolysis", "mcpba", "epoxidation", "os04", "dihydroxylation"]
        add_alkyne_kw = ["alkyne", "alkyne", "≡ ", "ethyne", "acetylene", "propyne", "phenylacetylene",
                        "lindlar", "na/nh3", "liquid nh3", "hgso4", "hydration of alkyne"]
        oxid_kw = ["oxid", " alcohol", "ol ", "1°", "2°", "primary alc", "secondary alc",
                  "pcc", "swern", "dess-martin", "jones", "kmno4", "cr", "tollens", "fehling",
                  "benedict", "chromic", "pyridinium"]
        arom_kw = ["benzene", "toluene", "phenol", "anisole", "nitrobenzene", "aniline",
                 "acetophenone", "benzoic", "aromatic", "sear", "friedel", "febr3", "alcl3",
                 "nitration", "sulfonation", "halogenation (aromatic)", "fcr", "acylat", "alkylat"]

        # Score each category
        scores = {}
        for kw_list, cat in [(sub_kw, "substitution"), (elim_kw, "elimination"),
                              (add_carbonyl_kw, "carbonyl_addition"), (add_alkene_kw, "alkene_addition"),
                              (add_alkyne_kw, "alkyne_addition"), (oxid_kw, "oxidation"),
                              (arom_kw, "aromatic_se")]:
            score = sum(1 for kw in kw_list if kw in r_combined)
            if score > 0:
                scores[cat] = score

        if scores:
            return max(scores, key=scores.get)
        return None

    def _predict_by_type(self, reactants: str, reagents: str, rtype: str, dl: str) -> str:
        data = REACTION_PREDICTION_DATA.get(rtype.replace(" ", "_"))
        if not data:
            return self._general_prediction(reactants, reagents, dl)

        lines = [f"## 反应预测: {data.get('name', rtype)}", ""]
        lines.append(f"**反应物**: {reactants}")
        lines.append(f"**试剂**: {reagents}")
        lines.append("")

        # Find matching pattern
        best_match = None
        best_score = 0
        for pat in data.get("patterns", []):
            score = 0
            search_in = f"{pat.get('substrate', '')} {pat.get('reagent', '')} ".lower()
            for word in f"{reactants} {reagents}".lower().split():
                if len(word) > 2 and word in search_in:
                    score += 1
            if score > best_score:
                best_score = score
                best_match = pat

        if best_match:
            lines.append(f"### 🎯 预测主要产物")
            lines.append(f"**{best_match.get('product', '未知')}**")
            
            if dl != "brief":
                if "conditions" in best_match:
                    lines.append(f"\n**条件**: {best_match['conditions']}")
                if "mechanism" in best_match:
                    lines.append(f"\n**机理**: {best_match['mechanism']}")
                if "note" in best_match:
                    lines.append(f"\n📝 **注**: {best_match['note']}")
                if "regioselectivity" in best_match:
                    lines.append(f"\n**区域选择性**: {best_match['regioselectivity']}")
                if "stereochemistry" in best_match:
                    lines.append(f"**立体化学**: {best_match['stereochemistry']}")
                if "advantage" in best_match:
                    lines.append(f"✅ **优势**: {best_match['advantage']}")
                if "issues" in best_match:
                    lines.append(f"⚠️ **注意**: {best_match['issues']}")

        # Show additional info for detailed mode
        if dl == "detailed":
            extra_keys = [k for k in data.keys() if k not in ("name", "patterns")]
            for ek in extra_keys:
                edata = data[ek]
                if isinstance(edata, dict):
                    lines.append(f"\n### {ek.replace('_', ' ').title()}")
                    for k, v in edata.items():
                        if isinstance(v, list):
                            lines.append(f"- **{k}**: {', '.join(v)}")
                        else:
                            lines.append(f"- **{k}**: {v}")

        if not best_match:
            lines.append("\n未找到精确匹配。以下是该类反应的一般规律:")
            for pat in data.get("patterns", [])[:3]:
                lines.append(f"- {pat.get('substrate', '?')} + {pat.get('reagent', '?')}")
                lines.append(f"  → {pat.get('product', '?')}")

        return "\n".join(lines)

    def _general_prediction(self, reactants: str, reagents: str, dl: str) -> str:
        lines = [f"## 反应预测", ""]
        lines.append(f"**反应物**: {reactants}")
        lines.append(f"**试剂**: {reagents}")
        lines.append("")
        lines.append("> ⚠️ 无法精确分类此反应。以下是基于关键词的一般分析:\n")
        
        # Keyword-based heuristics
        hints = []
        combined = f"{reactants} {reagents}".lower()

        if any(w in combined for w in ["hbr", "hi", "hcl"]):
            if any(w in combined for w in ["alkene", "ene", "propene", "styrene", "cyclohex"]):
                hints.append("可能是 **烯烃亲电加成** (Markovnikov 规则)")
                if "peroxide" in combined or "roor" in combined:
                    hints.append("  → 注意: 过氧化物存在时为 **反马氏** (自由基加成)")
            elif any(w in combined for w in ["alcohol", "ol "]):
                hints.append("可能是 **醇的卤化** (需先质子化)")

        if any(w in combined for w in ["naoet", "koh", "naoh"]) and ("heat" in combined or "Δ" in combined):
            hints.append("可能是 **消除反应 (E2)** — 检查 β-H 是否存在")

        if "strong base" in combined or "tbuok" in combined or "lda" in combined:
            hints.append("强碱存在 → 可能发生 **消除 (E2)** 而非取代")

        if "pcc" in combined or "swern" in combined or "dess-martin" in combined:
            hints.append("温和氧化剂 → **醛/酮** (伯醇→醛, 仲醇→酮, 不进一步氧化)")

        if "grignard" in combined or "rmgx" in combined or "rli" in combined:
            hints.append("有机金属试剂 → **对羰基的亲核加成** → 生成醇")

        if "h2" in combined and ("pd" in combined or "pt" in combined or "ni" in combined):
            hints.append("**催化氢化** → 饱和化合物 (syn addition)")

        if hints:
            for h in hints:
                lines.append(f"- {h}")
        else:
            lines.append("- 请提供更具体的反应物或试剂信息以便精确预测")
            lines.append("")
            lines.append("### 支持的反应类型:")
            for key, val in REACTION_PREDICTION_DATA.items():
                lines.append(f"- **{val.get('name', key)}**")

        return "\n".join(lines)
