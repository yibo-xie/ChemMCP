import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AsymmetricSynthesisGuide(BaseTool):
    """
    不对称合成方法选择指南 - 根据目标产物类型推荐合适的不对称合成策略。
    涵盖手性池、不对称催化、酶法拆分、结晶诱导等方法。
    """
    __version__ = "0.1.0"
    name             = "AsymmetricSynthesisGuide"
    func_name        = "asymmetric_synthesis_guide"
    description      = "Guide for selecting asymmetric synthesis methods based on target molecule type, substrate, and constraints."
    implementation_description = "Knowledge-based system covering chiral pool, asymmetric catalysis (hydrogenation, epoxidation, dihydroxylation), enzymatic resolution, crystallization-induced dynamic resolution, and auxiliary-based approaches."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Asymmetric Synthesis", "Stereocontrol", "Chirality", "Catalysis", "Enantioselectivity"]
    required_envs    = []

    code_input_sig   = [
        ("target_type", "str", "N/A", "Target product type: 'chiral alcohol', 'chiral amine', 'chiral acid', 'amino acid', 'epoxide', 'chiral ketone', 'chiral alkene', 'allene', 'axial chirality'."),
        ("substrate_hint", "str", "", "Optional hint about starting material (e.g., 'ketone', 'alkene', 'imine')."),
        ("constraints", "str", "", "Optional constraints: 'industrial scale', 'high ee (>99%)', 'metal-free', 'mild conditions'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'target_type [substrate_hint] [constraints]'. Example: 'chiral alcohol ketone high ee'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed asymmetric synthesis guide with ranked method recommendations."),
    ]

    examples         = [
        {
            "code_input": {"target_type": "chiral alcohol", "substrate_hint": "ketone", "constraints": ""},
            "text_input": {"input_params": "chiral alcohol ketone"},
            "output": {"result": "Recommended: CBS Reduction..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_knowledge_base()

    def _build_knowledge_base(self):
        """Build comprehensive asymmetric synthesis knowledge base."""
        self.methods = {
            # === Chiral Alcohol Methods ===
            "chiral alcohol": [
                {
                    "name": "CBS (Corey-Bakshi-Shibata) Reduction",
                    "catalyst": "CBS catalyst (oxazaborolidine) + BH3·THF or catecholborane",
                    "substrate": "Prochiral ketone → chiral secondary alcohol",
                    "typical_ee": "90-99%+",
                    "conditions": "0°C to rt, inert atmosphere, THF solvent",
                    "pros": ["Very high enantioselectivity", "Broad substrate scope", "Well-understood mechanism"],
                    "cons": ["CBS catalyst is moisture-sensitive/expensive", "BH3 is pyrophoric", "Stoichiometric or sub-stoichiometric boron reagent needed"],
                    "scale": "Lab to multi-kg",
                    "industrial_examples": "Merck HIV protease inhibitor synthesis; Taxol side chain",
                    "metal_free": True,
                    "cost": "Medium-High",
                    "ranking_score": 95,
                },
                {
                    "name": "Asymmetric Hydrogenation (Noyori-type)",
                    "catalyst": "Ru-BINAP or Ru-DuPHOS complexes",
                    "substrate": "β-Keto ester / β-ketoamide / functionalized ketone → chiral alcohol",
                    "typical_ee": "95-99%+",
                    "conditions": "H2 (1-100 atm), 25-50°C, alcoholic solvent",
                    "pros": ["Atom economical (H2)", "Catalytic (low loading possible)", "Industrial proven", "High throughput"],
                    "cons": ["Requires specialized equipment (pressure)", "Expensive Ru catalysts", "Substrate scope limited to certain ketone types"],
                    "scale": "Lab to >100 ton industrial",
                    "industrial_examples": "Monsanto L-DOPA process; Takasago (-)-menthol; Metolachlor (Syngenta)",
                    "metal_free": False,
                    "cost": "Medium (catalyst recycling possible)",
                    "ranking_score": 98,
                },
                {
                    "name": "Enzymatic Ketone Reduction (KRED)",
                    "catalyst": "Ketoreductase enzyme + NAD(P)H cofactor recycling system",
                    "substrate": "Prochiral ketone → chiral alcohol",
                    "typical_ee": "95-99%+",
                    "conditions": "Aqueous buffer or biphasic, 25-37°C, pH 6-8",
                    "pros": ["Excellent ee", "Green chemistry (water as solvent)", "Mild conditions", "Catalyst = protein (renewable)"],
                    "cons": ["Enzyme engineering may be needed for novel substrates", "Cofactor recycling adds complexity", "Limited substrate scope per enzyme"],
                    "scale": "Lab to industrial (Codexis, etc.)",
                    "industrial_examples": "Statins (atorvastatin intermediate); Montelukast intermediate",
                    "metal_free": True,
                    "cost": "Low-Medium (at scale)",
                    "ranking_score": 92,
                },
                {
                    "name": "Baker's Yeast Reduction",
                    "catalyst": "Saccharomyces cerevisiae (whole-cell biocatalyst)",
                    "substrate": "β-Keto ester / simple ketone → chiral alcohol",
                    "typical_ee": "70-98%",
                    "conditions": "Aqueous glucose medium, 25-30°C, aerobic/anaerobic",
                    "pros": ["Very cheap", "No special equipment", "Environmentally benign"],
                    "cons": ["Variable selectivity", "Over-reduction possible", "Product isolation from biomass can be tricky"],
                    "scale": "Lab scale typically",
                    "industrial_examples": "Historical: (S)-β-hydroxy esters",
                    "metal_free": True,
                    "cost": "Very Low",
                    "ranking_score": 70,
                },
                {
                    "name": "Chiral Auxiliary-Based Reduction",
                    "catalyst": "Chiral auxiliary (Evans oxazolidinone, etc.) + standard reducing agent",
                    "substrate": "Auxiliary-bound ketone/iminium → chiral product after auxiliary removal",
                    "typical_ee": "95-99%+",
                    "conditions": "Standard reduction conditions (NaBH4, LiAlH4, etc.), then cleavage",
                    "pros": ["Predictable stereochemistry", "Broad scope", "Well-established protocols"],
                    "cons": ["Additional steps (attach/remove auxiliary)", "Stoichiometric auxiliary (waste)", "Lower atom economy"],
                    "scale": "Lab to kg",
                    "industrial_examples": "Evans aldol methodology widely used in pharma",
                    "metal_free": True,
                    "cost": "Medium",
                    "ranking_score": 80,
                },
            ],
            # === Chiral Amine Methods ===
            "chiral amine": [
                {
                    "name": "Asymmetric Reductive Amination",
                    "catalyst": "Ir/phosphine or Rh/diphosphine complex, or organocatalyst",
                    "substrate": "Ketone/aldehyde + ammonia/amine source → chiral amine",
                    "typical_ee": "90-99%",
                    "conditions": "H2 pressure or silane reducing agent, 25-60°C",
                    "pros": ["Direct access to chiral amines", "One-pot from carbonyl", "Good functional group tolerance"],
                    "cons": ["Metal catalyst cost", "Competing imine formation issues", "Need for H2 source"],
                    "scale": "Lab to pilot",
                    "industrial_examples": "Sitagliptin (Merck) via asymmetric hydrogenation of enamine",
                    "metal_free": False,
                    "cost": "Medium-High",
                    "ranking_score": 88,
                },
                {
                    "name": "Enzymatic Transamination / Amine Dehydrogenase",
                    "catalyst": "Transaminase (TA) or amine dehydrogenase (AmDH)",
                    "substrate": "Ketone + amine donor (isopropylamine, alanine) → chiral amine",
                    "typical_ee": "95-99%+",
                    "conditions": "Aqueous buffer, PLP cofactor, 30-50°C",
                    "pros": ["Excellent ee", "Green", "Catalytic", "Broadening scope via protein engineering"],
                    "cons": ["Equilibrium limitation (need to drive)", "Cofactor required", "Some substrates not accepted"],
                    "scale": "Lab to industrial",
                    "industrial_examples": "Sitagliptin (Codexis transaminase process); API intermediates",
                    "metal_free": True,
                    "cost": "Low-Medium",
                    "ranking_score": 90,
                },
                {
                    "name": "Chiral Resolution (Classical)",
                    "catalyst": "Chiral resolving agent (tartaric acid, camphorsulfonic acid, etc.)",
                    "substrate": "Racemic amine → diastereomeric salts → separation",
                    "typical_ee": ">99% after recrystallization",
                    "conditions": "Salt formation in solvent, crystallization, basification",
                    "pros": ["Simple equipment", "Very high purity achievable", "Max yield = 50% per cycle (with racemization can approach 100%)"],
                    "cons": ["Max 50% yield without racemization", "Resolving agent cost/waste", "Many crystallizations needed"],
                    "scale": "Industrial (very common)",
                    "industrial_examples": "(S)-ibuprofen; (R)-(−)-3-chloro-1-phenylpropanamine; Many APIs",
                    "metal_free": True,
                    "cost": "Low-Medium",
                    "ranking_score": 75,
                },
            ],
            # === Amino Acid Methods ===
            "amino acid": [
                {
                    "name": "Chiral Pool (Natural Amino Acids)",
                    "catalyst": "N/A - use natural L-amino acids directly",
                    "substrate": "Natural L-amino acids as starting materials",
                    "typical_ee": "100% (naturally occurring)",
                    "conditions": "Standard peptide/protection chemistry",
                    "pros": ["Perfect chirality", "Cheap at scale", "Wide variety available (20 proteinogenic + many non-natural)", "Well-established protection strategies"],
                    "cons": ["Limited to L-configuration naturally", "May need derivatization", "Not all structures accessible"],
                    "scale": "Any scale",
                    "industrial_examples": "Penicillin; Peptide drugs; Semi-synthetic antibiotics",
                    "metal_free": True,
                    "cost": "Low",
                    "ranking_score": 85,
                },
                {
                    "name": "Strecker Synthesis (Asymmetric)",
                    "catalyst": "Chiral Lewis acid or phase-transfer catalyst",
                    "substrate": "Aldehyde + amine + cyanide → α-amino nitrile → amino acid",
                    "typical_ee": "80-97%",
                    "conditions": "Various: PTC (NaOH/toluene), Lewis acid (−40°C), or organocatalyst",
                    "pros": ["Direct amino acid synthesis", "Flexible aldehyde component", "Can make both D and L forms"],
                    "cons": ["Cyanide handling", "Hydrolysis step needed", "ee can vary with substrate"],
                    "scale": "Lab to pilot",
                    "industrial_examples": "Non-proteinogenic amino acid synthesis",
                    "metal_free": "Varies by variant",
                    "cost": "Medium",
                    "ranking_score": 78,
                },
                {
                    "name": "Enzymatic Resolution (Aminoacylase / Hydrolase)",
                    "catalyst": "Aminoacylase, lipase, or esterase",
                    "substrate": "N-Acetyl-DL-amino acid → L-amino acid + D-N-acetyl (recyclable)",
                    "typical_ee": ">99%",
                    "conditions": "Aqueous, pH 7-8, 37°C, immobilized enzyme preferred",
                    "pros": ["Industrial proven (Tanabe Seiyaku since 1954)", "Continuous process possible", "High yields with racemization loop"],
                    "cons": ["Only 50% per pass without racemization", "Enzyme immobilization needed for stability"],
                    "scale": "Multi-ton industrial",
                    "industrial_examples": "Tanabe L-Met, L-Val, L-Ala process (50+ years running)",
                    "metal_free": True,
                    "cost": "Low (at industrial scale)",
                    "ranking_score": 87,
                },
            ],
            # === Epoxide Methods ===
            "epoxide": [
                {
                    "name": "Sharpless Epoxidation",
                    "catalyst": "Ti(OiPr)4 + chiral tartrate ester (DET or DIET)",
                    "substrate": "Primary allylic alcohol → 2,3-epoxy alcohol",
                    "typical_ee": "90-99%+",
                    "conditions": "CH2Cl2, −20°C to 0°C, TBHP as oxidant, molecular sieves (4Å)",
                    "pros": ["Reliable and predictable (enantioface determined by tartrate)", "High ee", "Broad allylic alcohol scope"],
                    "cons": ["Only works for allylic alcohols", "Ti reagent sensitive to water", "TBHP handling"],
                    "scale": "Lab to multi-kg",
                    "industrial_examples": "(-)-Disparlure; Glycidol; (+)-disparlure pheromone; Total synthesis intermediates",
                    "metal_free": False,
                    "cost": "Medium",
                    "ranking_score": 93,
                },
                {
                    "name": "Sharpless Asymmetric Dihydroxylation (AD)",
                    "catalyst": "OsO4 (catalytic) + chiral ligand (DHQD/DHQ PHAL, etc.) + co-oxidant",
                    "substrate": "Alkene → vicinal diol (can be converted to epoxide)",
                    "typical_ee": "90-99%+",
                    "conditions": "t-BuOH/H2O, 0°C, K3[Fe(CN)6] or NMO as co-oxidant",
                    "pros": ["Works on most alkenes (not just allylic alcohols)", "Predictable face selectivity (table exists)", "Catalytic Os"],
                    "cons": ["Osmium toxicity/cost", "Slow reaction for some substrates", "Over-oxidation possible"],
                    "scale": "Lab to kg",
                    "industrial_examples": "Taxol side chain; Chloramphenicol; Indinavir intermediate",
                    "metal_free": False,
                    "cost": "Medium-High",
                    "ranking_score": 91,
                },
                {
                    "name": "Jacobsen/Katsuki Epoxidation (Mn-Salen)",
                    "catalyst": "Chiral Mn(III)-salen complex",
                    "substrate": "Unfunctionalized cis-alkenes (especially aryl-substituted) → epoxide",
                    "typical_ee": "85-98%",
                    "conditions": "CH2Cl2, 0°C to rt, oxidant (NaOCl, m-CPBA, PhIO)",
                    "pros": ["Works on unfunctionalized alkenes", "No directing group needed", "Good for cis-disubstituted and terminal alkenes"],
                    "cons": ["Best for cis-alkenes (trans lower ee)", "Mn-salen synthesis can be tedious", "Oxidant choice matters"],
                    "scale": "Lab to kg",
                    "industrial_examples": "Indirubin derivatives; Pharmacologically active epoxides",
                    "metal_free": False,
                    "cost": "Medium",
                    "ranking_score": 86,
                },
            ],
            # === General / Broad Methods ===
            "general": [
                {
                    "name": "Organocatalysis (Proline & Derivatives)",
                    "catalyst": "L-Proline or Jorgensen-Hayashi catalyst, MacMillan catalyst, etc.",
                    "substrate": "Aldehydes/ketones (aldol, Mannich, α-amination, α-oxyamination)",
                    "typical_ee": "90-99%+",
                    "conditions": "Often solvent-free or simple solvents (DMSO, CHCl3), rt to 0°C",
                    "pros": ["Metal-free", "Cheap catalysts", "Green", "Nobel Prize recognized (2021)"],
                    "cons": ["High catalyst loading often needed (10-30 mol%)", "Limited to activated substrates (aldehydes mostly)", "Reaction times can be long"],
                    "scale": "Lab to pilot",
                    "industrial_examples": "List et al. organocatalytic processes under development",
                    "metal_free": True,
                    "cost": "Very Low-Low",
                    "ranking_score": 84,
                },
                {
                    "name": "Phase-Transfer Catalysis (PTC) Asymmetric Alkylation",
                    "catalyst": "Chiral quaternary ammonium salt (Maruoka, Shibasaki, etc.)",
                    "substrate": "Glycinate Schiff base / phenol / other acidic C-H → alkylated product",
                    "typical_ee": "90-99%",
                    "conditions": "Two-phase (toluene/aq. NaOH or KOH), 0°C to rt",
                    "pros": ["Simple setup", "Operates under strong basic conditions", "Scalable", "Metal-free"],
                    "cons": ["Limited substrate classes well-developed", "Catalyst synthesis can be complex", "Emulsion formation issues"],
                    "scale": "Lab to multi-ton industrial",
                    "industrial_examples": "(R)- or (S)-phenylalanine (via glycinate Schiff base); Merck HIV drug intermediate",
                    "metal_free": True,
                    "cost": "Medium",
                    "ranking_score": 82,
                },
                {
                    "name": "Dynamic Kinetic Resolution (DKR)",
                    "catalyst": "Enzyme (resolution) + metal catalyst (racemization) or organocatalyst",
                    "substrate": "Racemate → single enantiomer (theoretical 100% yield)",
                    "typical_ee": "95-99%, yield up to 90%+",
                    "conditions": "Combines resolution and racemization conditions",
                    "pros": ["Can exceed 50% yield limit of classical resolution", "One-pot transformation", "High efficiency"],
                    "cons": ["Matching rates of resolution and racemization is tricky", "More complex optimization", "Catalyst compatibility issues"],
                    "scale": "Lab to pilot/industrial",
                    "industrial_examples": "Alcohol DKR (lipase + Ru); Amine DKR (enzyme + Pd); β-Keto ester DKR",
                    "metal_free": "Usually requires metal for racemization",
                    "cost": "Medium-High",
                    "ranking_score": 89,
                },
            ],
        }

        # Cross-reference map for related target types
        self.target_aliases = {
            "chiral ketone": "chiral alcohol",  # Often made via reduction of prochiral ketone
            "secondary alcohol": "chiral alcohol",
            "primary alcohol": "chiral alcohol",
            "α-amino acid": "amino acid",
            "non-natural amino acid": "amino acid",
            "oxirane": "epoxide",
            "chirality": "general",
            "enantiopure": "general",
            "stereoselective": "general",
        }

    def _run_base(self, target_type: str, substrate_hint: str = "", constraints: str = "") -> str:
        """Generate asymmetric synthesis guide."""
        target_key = target_type.lower().strip()
        if target_key in self.target_aliases:
            target_key = self.target_aliases[target_key]

        methods = self.methods.get(target_key)
        if not methods:
            # Try general methods as fallback
            methods = self.methods.get("general", [])
            if not methods:
                raise ChemMCPError(f"No methods found for target type '{target_type}'. Try: {list(self.methods.keys())}")

        # Apply constraint-based ranking adjustments
        constraint_lower = constraints.lower() if constraints else ""
        substrate_lower = substrate_hint.lower() if substrate_hint else ""

        for m in methods:
            score = float(m.get("ranking_score", 50))

            if constraint_lower:
                if "industrial" in constraint_lower and m["scale"] in ("Industrial", "multi-ton industrial", "Lab to >100 ton industrial"):
                    score += 15
                if "high ee" in constraint_lower or ">99%" in constraint_lower:
                    ee = m.get("typical_ee", "")
                    if "99%" in ee:
                        score += 10
                    elif "95" in ee or "97" in ee:
                        score += 5
                if "metal free" in constraint_lower or "no metal" in constraint_lower:
                    if m.get("metal_free") is True:
                        score += 15
                    else:
                        score -= 20
                if "cost-effective" in constraint_lower or "cheap" in constraint_lower:
                    cost = m.get("cost", "")
                    if "Low" in cost or "Very Low" in cost:
                        score += 10
                    elif "High" in cost:
                        score -= 5

            if substrate_lower:
                sub = m.get("substrate", "").lower()
                for kw in substrate_lower.split():
                    if kw in sub and len(kw) > 2:
                        score += 5

            m["_adjusted_score"] = score

        methods.sort(key=lambda x: x.get("_adjusted_score", 0), reverse=True)

        parts = [f"## Asymmetric Synthesis Guide: {target_type.title()}\n"]
        parts.append(f"Found **{len(methods)}** recommended methods (ranked by suitability):\n")

        for i, m in enumerate(methods, 1):
            score = m.get("_adjusted_score", 0)
            star = " ⭐⭐⭐" if score >= 100 else " ⭐⭐" if score >= 90 else " ⭐"
            parts.append(f"### Method {i}: {m['name']}{star} *(Score: {score:.0f})*")
            parts.append(f"- **Catalyst/System:** {m['catalyst']}")
            parts.append(f"- **Substrate:** {m['substrate']}")
            parts.append(f"- **Typical ee:** {m['typical_ee']}")
            parts.append(f"- **Conditions:** {m['conditions']}")
            parts.append(f"- **Pros:** {'; '.join(m['pros'])}")
            parts.append(f"- **Cons:** {'; '.join(m['cons'])}")
            parts.append(f"- **Scale:** {m['scale']}")
            parts.append(f"- **Industrial Examples:** {m['industrial_examples']}")
            parts.append(f"- **Metal-free:** {'Yes ✅' if m['metal_free'] else 'No ❌'}")
            parts.append(f"- **Cost:** {m['cost']}")
            parts.append("")

        # Add decision flowchart summary
        parts.append("\n---\n### 🧭 Quick Decision Guide\n")
        parts.append("| Your Priority | Recommended Approach |")
        parts.append("|---|---|")
        parts.append("| Highest ee (>99%) | CBS reduction / Asymmetric hydrogenation / Enzymatic |")
        parts.append("| Industrial scale | Asymmetric hydrogenation / Enzymatic / Classical resolution |")
        parts.append("| Metal-free | Organocatalysis / Biocatalysis / Chiral pool |")
        parts.append("| Lowest cost | Baker's yeast / Chiral pool / Simple resolution |")
        parts.append("| Novel substrate | Chiral auxiliary / Develop new enzymatic route |")
        parts.append("| Green chemistry | Enzymatic / Organocatalysis / Water-based |")

        return "\n".join(parts)

    def _run_text(self, input_params: str) -> str:
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Input must include target_type. Format: 'target_type [substrate_hint] [constraints]'")
        target = parts[0]
        substrate = parts[1] if len(parts) > 1 else ""
        constraints = " ".join(parts[2:]) if len(parts) > 2 else ""
        return self._run_base(target, substrate, constraints)
