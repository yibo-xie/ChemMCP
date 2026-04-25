import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive catalyst database organized by type and reaction
_CATALYST_DB = {
    "acid": {
        "description": "Brønsted and Lewis acid catalysts for various organic transformations.",
        "catalysts": {
            "H2SO4 (conc.)": {
                "type": "Brønsted strong acid", "strength": "very strong", "pKa": -3,
                "uses": ["esterification", "dehydration", "hydrolysis", "nitration (mixed acid)", "Friedel-Crafts (with Lewis acid)", "polymerization initiation"],
                "conditions": "0-180°C depending on application. Often used as solvent or co-solvent.",
                "pros": ["Very cheap", "Strongest common acid", "Versatile", "Dehydrating agent"],
                "cons": ["Highly corrosive", "Causes charring/sulfonation side reactions", "Hard to remove", "Oxidizing at high concentration"],
                "safety": "☠️ Causes severe burns. Add acid to water (never reverse). Use face shield for conc. work.",
            },
            "HCl (conc.)": {
                "type": "Brønsted strong acid", "strength": "strong", "pKa": -7,
                "uses": ["ester hydrolysis", "deprotection (Boc, t-butyl)", "chlorination", "acidic workup", "diazotization"],
                "conditions": "Reflux in dioxane or aqueous/organic mixtures. Gas evolution possible.",
                "pros": ["Volatile (can be removed)", "Doesn't oxidize like HNO3/H2SO4", "Cheap"],
                "cons": ["Corrosive fumes", "Aqueous — not always compatible with water-sensitive substrates", "Can cause chlorination side reactions"],
                "safety": "Corrosive fumes. Use in fume hood. Can damage glass over time.",
            },
            "p-Toluenesulfonic acid (p-TsOH)": {
                "type": "Brønsted strong acid (organic)", "strength": "strong", "pKa": -2.8,
                "uses": ["esterification", "acetal formation/protection", "Pechmann condensation", "Mannich reaction", "mild acid catalysis"],
                "conditions": "RT to reflux in toluene/benzene with Dean-Stark trap. Often used with molecular sieves.",
                "pros": ["Solid (easy to handle)", "Milder than H2SO4", "Non-oxidizing", "Soluble in organic solvents", "Easy to remove (wash with NaHCO3)"],
                "cons": ["More expensive than mineral acids", "Still corrosive", "Can be hard to separate from product"],
                "safety": "Corrosive solid. Avoid dust inhalation. Less hazardous than mineral acids.",
            },
            "AlCl3": {
                "type": "Lewis acid", "strength": "very strong (hard)", "pKa": "N/A",
                "uses": ["Friedel-Crafts alkylation", "Friedel-Crafts acylation", "isomerization", "Diels-Alder catalysis"],
                "conditions": "Anhydrous DCM or nitrobenzene as solvent. 0°C → RT. Stoichiometric (1+ eq typically).",
                "pros": ["Very powerful FC catalyst", "Well-established protocols", "Relatively cheap"],
                "cons": ["Moisture sensitive (quenches violently with H2O)", "Often needed stoichiometrically", "Difficult workup (aqueous quench messy)", "Over-alkylation common"],
                "safety": "React violently with water (exothermic HCl evolution). Handle under N2. Corrosive.",
            },
            "BF3·OEt2": {
                "type": "Lewis acid", "strength": "strong", "pKa": "N/A",
                "uses": ["epoxide opening", "carbocation rearrangements", "Mukaiyama aldol", "polymerization initiator"],
                "conditions": "Anhydrous DCM, 0°C → RT. Catalytic (5-20 mol%).",
                "pros": ["Catalytic amounts sufficient", "Liquid (easy to measure)", "Mild compared to AlCl3", "Selective for certain transformations"],
                "cons": ["Expensive", "Moisture sensitive", "Gaseous BF3 byproduct concerns", "Can cause rearrangement side reactions"],
                "safety": "Releases HF if exposed to moisture. Toxic. Use in fume hood.",
            },
            "FeCl3": {
                "type": "Lewis acid", "strength": "moderate-strong", "pKa": "N/A",
                "uses": ["Friedel-Crafts (milder than AlCl3)", "aromatic halogenation catalyst", "phenol coupling", "oxidation"],
                "conditions": "DCM or nitromethane. RT to 80°C. Can be catalytic or stoichiometric.",
                "pros": ["Cheaper than many Lewis acids", "Less moisture-sensitive than AlCl3", "Can be used catalytically", "Green profile (iron-based)"],
                "cons": ["Weaker than AlCl3 for difficult FC reactions", "Colored (yellow-green), can complicate purification", "Can act as one-electron oxidizer"],
                "safety": "Corrosive, environmental toxin. Moderate hazard level.",
            },
            "Amberlyst-15": {
                "type": "Solid acid (ion exchange resin)", "strength": "moderate-strong (like p-TsOH)", "pKa": "~-2",
                "uses": ["esterification", "protection/deprotection", "alkylation", "continuous flow chemistry"],
                "conditions": "Pack in column or stir in flask. RT to reflux. Can be filtered off and reused.",
                "pros": ["Reusable", "Easy separation (filtration)", "Non-toxic", "Works in green chemistry", "No aqueous workup needed"],
                "cons": ["Limited surface area/slow diffusion", "Not suitable for all substrates", "Mechanical degradation over time", "Swelling in solvents"],
                "safety": "Very safe — solid resin. Minimal hazards.",
            },
        }
    },
    "base": {
        "description": "Base catalysts for elimination, condensation, deprotonation, and substitution reactions.",
        "catalysts": {
            "NaOH/KOH": {
                "type": "Inorganic hydroxide (strong base)", "strength": "very strong", "pKb": "-1.7 / -0.6",
                "uses": ["saponification", "E2 elimination", "aldol condensation", "Claisen-Schmidt", "haloform reaction"],
                "conditions": "Aqueous or alcoholic solution. RT to reflux. Concentration: 1-10 M typical.",
                "pros": ["Very cheap", "Very strong base", "Versatile", "High solubility in water/alcohol"],
                "cons": ["Not compatible with water-sensitive substrates", "Can cause hydrolysis of esters/amides", "Promotes E2 over SN2 with secondary/tertiary substrates"],
                "safety": "Caustic — causes severe burns. Dissolves glass slowly. Use PPE.",
            },
            "NaOEt/KOtBu": {
                "type": "Alkoxide base", "strength": "very strong", "pKb": "-~2",
                "uses": ["E2 elimination", "Claisen condensation", "Dieckmann condensation", "intramolecular alkylation"],
                "conditions": "Anhydrous ethanol (NaOEt) or t-butanol (KOtBu). RT to reflux under N2.",
                "pros": ["Strong, non-nucleophilic bases", "KOtBu is bulky → promotes E2 selectively", "Solubility in organic solvents"],
                "cons": ["Moisture sensitive", "NaOEt can act as nucleophile (SN2 side reaction)", "Require anhydrous conditions"],
                "safety": "Flammable solutions. Moisture-sensitive. Standard base precautions.",
            },
            "t-BuOK": {
                "type": "Bulky alkoxide (sterically hindered)", "strength": "very strong", "pKb": "~-2",
                "uses": ["E2 elimination (especially Hofmann-oriented)", "deprotonation of weak C-H acids", "isomerization"],
                "conditions": "t-BuOH or THF as solvent. RT to 80°C. Often used neat or concentrated.",
                "pros": ["★ Best E2 promoter — steric bulk prevents SN2", "Deprotonates very weak acids", "Commercially available"],
                "cons": ["Expensive", "Very hygroscopic", "t-BuOH has high BP (hard to remove)", "Too bulky for some applications"],
                "safety": "Pyrophoric when dry! Store under N2. Handle with care.",
            },
            "LDA (Lithium Diisopropylamide)": {
                "type": "Sterically hindered amide base", "strength": "very strong (pKa conj ~36)", "pKb": "N/A",
                "uses": ["kinetic enolate formation", "regioselective deprotonation", "elimination to form alkynes"],
                "conditions": "THF as solvent, -78°C (dry ice/acetone). Generated from i-Pr2NH + n-BuLi.",
                "pros": ["★ Gold standard for kinetic enolates", "Excellent regioselectivity", "Non-nucleophilic (bulky)", "Predictable behavior"],
                "cons": ["Must be prepared fresh (or purchased as solution)", "Requires cryogenic temperatures (-78°C)", "Expensive", "Moisture/air sensitive"],
                "safety": "⚠️ Pyrophoric! Strict inert atmosphere. Cryogenic handling. Professional use only.",
            },
            "DBU": {
                "type": "Amidine base (non-nucleophilic)", "strength": "strong", "pKa": "~12",
                "uses": ["elimination (E2)", "isomerization", "dehydrohalogenation", "carbene formation", "catalytic base"],
                "conditions": "Various organic solvents. RT to reflux. Often used catalytically (10-30 mol%).",
                "pros": ["Non-nucleophilic", "Liquid (easy to handle)", "Mild enough for sensitive substrates", "Can be used catalytically"],
                "cons": ["Expensive", "Not strong enough for some deprotonations", "Basic enough to epimerize stereocenters", "Has odor"],
                "safety": "Corrosive. Causes skin/eye irritation. Use gloves and fume hood.",
            },
            "DABCO": {
                "type": "Tertiary amine base", "strength": "moderate", "pKa": "8.8",
                "uses": ["nucleophilic catalyst (Baylis-Hillman)", "base for eliminations", "Schotten-Baumann conditions"],
                "conditions": "DCM or other aprotic solvents. RT.",
                "pros": ["Inexpensive", "Low toxicity", "Good nucleophilicity catalyst", "Solid (stable storage)"],
                "cons": ["Weak base (not for demanding deprotonations)", "Can act as nucleophile (may not be desired)", "Limited scope"],
                "safety": "Generally safe. Irritant. Low hazard profile.",
            },
            "NaH": {
                "type": "Hydride base", "strength": "very strong (H2 pKa ~35)", "pKb": "N/A",
                "uses": ["enolate formation", "alkoxide generation", "deprotonation of C-H, O-H, N-H", " Williamson ether synthesis"],
                "conditions": "DMF, THF, or DMSO. 0°C → RT under N2. H2 gas evolution!",
                "pros": ["Very strong base", "Byproduct is only H2 gas (clean)", "Available as mineral oil dispersion (easier handling)", "Cheap"],
                "cons": ["☠️ H2 gas evolution — fire/explosion risk", "Pyrophoric when finely divided", "Incompatible with everything protic", "Cannot be used catalytically"],
                "safety": "⚠️ ⚠️ **EXTREME CAUTION** — reacts VIOLENTLY with water releasing H2. Fire hazard. Use inert atmosphere, no sparks/flames nearby.",
            },
            "K2CO3/Cs2CO3": {
                "type": "Carbonate base (mild)", "strength": "weak-moderate", "pKb": "N/A",
                "uses": ["SN2 alkylation (mild)", "alkylation of phenols/amines", "Suzuki coupling (base component)", "O/N-alkylations"],
                "conditions": "Acetone, DMF, acetonitrile, or DMSO. RT to reflux. Cs2CO3 more soluble.",
                "pros": ["Mild — doesn't promote elimination", "Inexpensive", "Cs2CO3 soluble in organic media", "Compatible with many functional groups"],
                "cons": ["Too weak for enolate formation or E2", "Cs2CO3 expensive (but often worth it)", "Slow reaction rates sometimes"],
                "safety": "Generally safe. Dust irritant. Standard PPE sufficient.",
            },
        }
    },
    "metal": {
        "description": "Transition metal catalysts for cross-coupling, hydrogenation, oxidation, and C-H activation.",
        "catalysts": {
            "Pd(PPh3)4 (Tetrakis(triphenylphosphine)palladium(0))": {
                "type": "Pd(0) catalyst", "class": "cross-coupling",
                "uses": ["Suzuki coupling", "Heck reaction", "Sonogashira", "Stille coupling", "Buchwald-Hartwig amination", "Negishi coupling"],
                "conditions": "Toluene/EtOH/H2O or dioxane. 80-110°C under N2. 1-5 mol% loading.",
                "pros": ["Most versatile Pd(0) source", "Air-stable as solid (though Pd(0) degrades)", "Well-precedented", "Commercially available in high purity"],
                "cons": ["Phosphine ligands can be cumbersome", "Sensitive to oxygen in solution", "Expensive (Pd cost)", "Can form Pd black precipitate"],
                "safety": "Palladium compound — toxic. Phosphine ligands have odor. Use in fume hood.",
            },
            "Pd(dppf)Cl2": {
                "type": "Pd(II) pre-catalyst", "class": "cross-coupling",
                "uses": ["Suzuki-Miyaura coupling", "Buchwald-Hartwig amination", "C-N bond formation", "C-O bond formation"],
                "conditions": "Dioxane or 1,4-dioxane/water. 80-100°C. 1-3 mol%. Base required (K2CO3, Cs2CO3).",
                "pros": ["Excellent for Suzuki couplings", "dppf bidentate ligand stabilizes Pd", "Good air stability", "Reliable performance"],
                "cons": ["Expensive", "dppf ligand synthesis adds cost", "May need optimization for challenging substrates"],
                "safety": "Standard Pd catalyst precautions.",
            },
            "Pd/C (Palladium on Carbon)": {
                "type": "Heterogeneous Pd catalyst", "class": "reduction/hydrogenation",
                "uses": ["Hydrogenation (alkenes, alkynes, nitro, benzyl)", "transfer hydrogenation", "Carbonyl reduction (with H2)", "debenzylation", "dehalogenation"],
                "conditions": "EtOH, EtOAc, MeOH, or hexane. H2 balloon (1 atm) or Parr shaker (30-60 psi). RT to 50°C. 5-10 wt%.",
                "pros": ["Simple workup (filter off)", "Very versatile reductions", "Can be reused (sometimes)", "Well-understood mechanism"],
                "cons": ["PYROPHORIC WHEN DRY (!!!)", "Poisoned by S compounds", "Over-reduction possible", "Not selective for partial hydrogenation without modification (Lindlar's for cis)"],
                "safety": "🔥🔥 **ALWAYS KEEP WET!** Filter through Celite, never let dry on filter. Spontaneous combustion risk.",
            },
            "Grubbs Catalyst (2nd Gen)": {
                "type": "Ru carbene catalyst", "class": "olefin metathesis",
                "uses": ["Ring-closing metathesis (RCM)", "Cross metathesis (CM)", "Ring-opening metathesis (ROMP)", "Ring-opening cross metathesis (ROCM)"],
                "conditions": "DCM (reflux) or tolene (reflux). Degassed solvent. N2 atmosphere. 1-5 mol%.",
                "pros": ["Revolutionary for ring formation", "Tolerates many functional groups", "Well-established protocols", "2nd gen: more active, tolerant"],
                "cons": ["Very expensive (Ru complex)", "Air-sensitive (some versions)", "Ethylene gas evolution (must vent)", "Ruthenium removal can be tricky"],
                "safety": "Ruthenium compound — moderate toxicity. Use standard precautions.",
            },
            "CuI / CuCN": {
                "type": "Copper(I) catalyst", "class": "coupling/cycloaddition",
                "uses": ["Sonogashira coupling (with Pd)", "Glaser coupling", "Click chemistry (CuAAC)", "Ullmann coupling", "Goldberg reaction"],
                "conditions": "Varies widely. Sonogashira: Pd/Cu in amine solvent. Click: CuSO4 + sodium ascorbate in t-BuOH/H2O. RT.",
                "pros": ["Inexpensive", "Essential for Sonogashira and click chemistry", "Click reaction: bioorthogonal, reliable", "Copper-catalyzed reactions are atom-economical"],
                "cons": ["Cu salts colored (complicate purification)", "Homocoupling (Glaser) side reaction", "Copper removal can be problematic", "Some Cu salts are toxic"],
                "safety": "Copper salts are toxic to aquatic life. Avoid environmental release. Standard PPE.",
            },
            "Ni(dppp)Cl2": {
                "type": "Ni(II) catalyst", "class": "cross-coupling (Pd alternative)",
                "uses": ["Kumada coupling", "Negishi coupling", "Suzuki-type (Ni)", "reductive coupling", "C-O activation"],
                "conditions": "Ether or THF. RT to 60°C. Often with Grignard or organozinc reagents.",
                "pros": ["Much cheaper than Pd!", "Good for C(sp²)-C(sp³) couplings", "Can activate traditionally inert bonds (C-O)", "Earth-abundant metal"],
                "cons": ["More sensitive to air/moisture than Pd systems", "Narrower substrate scope in some cases", "Functional group tolerance lower than Pd", "Less developed (fewer precedents)"],
                "safety": "Nickel compound — allergen and suspected carcinogen. Handle with care.",
            },
            "Ti(OiPr)4 (Titanium tetraisopropoxide)": {
                "type": "Lewis acid / reductant", "class": "carbonyl addition / McMurry",
                "uses": ["McMurry coupling", "Sharpless epoxidation (with tartrate)", "Reductive coupling", "Pinacol coupling"],
                "conditions": "Varies: McMurry needs reducing agent (Zn). Sharpless: Ti(OiPr)4 + TBHP + tartrate, DCM, -20°C.",
                "pros": ["Sharpless epoxidation: ★ enantioselective!", "McMurry: forms alkenes from carbonyls", "Relatively inexpensive Ti source"],
                "cons": ["Moisture sensitive", "Sharpless requires careful optimization", "McMurry: low-yielding for many substrates", "Ti residues hard to remove"],
                "safety": "Irritant. Reacts with water releasing iPrOH. Standard PPE.",
            },
        }
    },
    "enzyme": {
        "description": "Biocatalysts for asymmetric synthesis, resolutions, and green chemistry applications.",
        "catalysts": {
            "Lipase (CAL-B, Candida antarctica)": {
                "type": "Hydrolase enzyme", "class": "resolution / esterification",
                "uses": ["Kinetic resolution of racemic alcohols/esters", "Regioselective acylation", "Transesterification", "Desymmetrization"],
                "conditions": "Aqueous buffer or organic solvent (MTBE, toluene). 25-40°C. pH 7-8 optimal. Vinyl acetate as acyl donor common.",
                "pros": ["★ Excellent enantioselectivity (often >99% ee)", "Mild conditions (aqueous, RT)", "Reusable (immobilized form)", "Broad substrate acceptance", "Green chemistry"],
                "cons": ["Slow (hours to days)", "Substrate size limitations", "Can be inhibited by products", "Immobilized enzyme expensive initially"],
                "safety": "Very safe — biological catalyst. No special hazards beyond normal lab practice.",
            },
            "KRED (Ketoreductase)": {
                "type": "Oxidoreductase enzyme", "class": "asymmetric reduction",
                "uses": ["Asymmetric reduction of ketones → chiral alcohols", "Prochiral ketone desymmetrization", "Deracemization of alcohols"],
                "conditions": "Aqueous/organic biphasic or cosolvent. 25-30°C. pH 6-8. NAD(P)H cofactor recycling system (GDH/glucose).",
                "pros": ["★ Outstanding enantio- and diastereoselectivity", "Direct access to chiral alcohols (no resolution waste)", "Cofactor recycling makes it practical", "Growing commercial availability (Codex, etc.)"],
                "cons": ["Requires cofactor (NADPH) system", "Narrow substrate scope per enzyme variant", "Optimization needed per substrate", "Enzyme cost (but improving)"],
                "safety": "Safe biological catalyst. Aqueous conditions.",
            },
            "P450 BM3 (Cytochrome P450)": {
                "type": "Monooxygenase enzyme", "class": "C-H oxidation",
                "uses": ["Regioselective C-H hydroxylation", "Epoxidation of alkenes", "Heteroatom oxidation"],
                "conditions": "Aqueous buffer, 25-30°C. Requires NADPH and O2 supply. Whole-cell or purified enzyme.",
                "pros": ["Site-selective oxidation (unmatched by chemical methods)", "Mild conditions (room temp, water)", "Catalytic in nature's oxidant (O2)", "Potential for late-stage functionalization"],
                "cons": ["Very slow (TOF often <1 min⁻¹)", "NADPH requirement (costly cofactor)", "Limited substrate scope (natural enzymes)", "Engineering required for non-native substrates"],
                "safety": "Biological system. Safe handling. Cell cultures require biosafety level 1 practices.",
            },
            "Nitrile hydratase / nitrilase": {
                "type": "Hydratase / lyase enzyme", "class": "nitrile conversion",
                "uses": ["Nitrile → amide (hydratase)", "Nitrile → carboxylic acid (nitrilase)", "Industrial production of acrylamide, nicotinamide"],
                "conditions": "Aqueous buffer, pH 7-8, 25-40°C. Whole cells often used.",
                "pros": ["Industrial-scale proven (acrylamide >10k tons/year)", "100% atom economy possible", "No protecting groups needed", "Mild conditions"],
                "cons": ["Product inhibition possible", "pH control critical", "Substrate/product toxicity to enzyme possible"],
                "safety": "Industrial biocatalyst. Safe under standard biochemical handling.",
            },
            "Penicillin G acylase (PGA)": {
                "type": "Amidase enzyme", "class": "antibiotic intermediate synthesis",
                "uses": ["Penicillin G → 6-APA (semisynthetic β-lactam precursor)", "Kinetic resolution of amines/esters", "amide synthesis (reverse hydrolysis)"],
                "conditions": "Aqueous buffer, pH 7.5-8.0, 35-37°C. Immobilized form available.",
                "pros": ["★ Industrial workhorse for β-lactam antibiotics", "Highly regioselective", "Immobilized — reusable", "Mild, aqueous conditions"],
                "cons": ["Specific to phenylacetyl (and similar) groups", "Product inhibition", "Narrow substrate scope outside penicillins"],
                "safety": "Safe biological catalyst. Standard biochemistry practices.",
            },
        }
    },
}

# Reaction type → catalyst type mapping
_REACTION_CATALYST_MAP = {
    "substitution_SN2": ["base:NaOEt", "base:K2CO3", "base:Cs2CO3"],
    "substitution_SN1": ["acid:H2SO4", "acid:p-TsOH", "acid:BF3·OEt2"],
    "elimination_E2": ["base:t-BuOK", "base:NaOEt/KOtBu", "base:DBU", "base:LDA"],
    "elimination_E1": ["acid:H2SO4", "acid:p-TsOH"],
    "cross_coupling": ["metal:Pd(PPh3)4", "metal:Pd(dppf)Cl2", "metal:Ni(dppp)Cl2", "metal:CuI"],
    "hydrogenation": ["metal:Pd/C"],
    "oxidation": ["metal:Ti(OiPr)4", "enzyme:P450 BM3"],
    "reduction_carbonyl": ["enzyme:KRED"],
    "resolution": ["enzyme:Lipase (CAL-B)"],
    "esterification": ["acid:H2SO4", "acid:p-TsOH", "acid:Amberlyst-15"],
    "hydrolysis": ["base:NaOH/KOH", "enzyme:Penicillin G acylase (PGA)"],
    "Friedel-Crafts": ["acid:AlCl3", "acid:FeCl3", "acid:BF3·OEt2"],
    "enolate_chemistry": ["base:LDA", "base:NaH", "base:KHMDS"],
    "alkene_metathesis": ["metal:Grubbs Catalyst (2nd Gen)"],
    "click_chemistry": ["metal:CuI"],
    "asymmetric_synthesis": ["enzyme:KRED", "enzyme:Lipase (CAL-B)", "metal:Ti(OiPr)4"],
    "nitrile_conversion": ["enzyme:Nitrile hydratase / nitrilase"],
    "C-H_activation": ["enzyme:P450 BM3", "metal:Pd catalysts"],
}


@ChemMCPManager.register_tool
class CatalystRecommender(BaseTool):
    """
    推荐催化剂的工具（酸、碱、金属、酶）。
    内置全面的催化剂数据库，涵盖酸、碱、过渡金属和生物酶四大类，支持按反应类型和约束条件推荐。
    """
    __version__      = "0.1.0"
    name             = "CatalystRecommender"
    func_name        = "recommend_catalyst"
    description      = "Recommend appropriate catalysts (acid, base, metal, enzyme) for a given reaction type, with detailed pros/cons, conditions, and safety information."
    implementation_description = "Uses embedded database of 26+ catalysts across 4 categories (acid, base, metal, enzyme) with detailed properties, uses, conditions, pros/cons, and safety data."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Catalysis", "Synthesis", "Organic Chemistry", "Green Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("reaction_type", "str", "N/A", "Type of reaction (e.g., 'cross_coupling', 'hydrogenation', 'esterification', 'E2', 'resolution')."),
        ("catalyst_type", "str", "any", "Preferred catalyst category: 'acid', 'base', 'metal', 'enzyme', or 'any' for all types."),
        ("constraints", "str", "None", "Constraints: 'cheap', 'green', 'mild_conditions', 'high_selectivity'. Use 'None' for default."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'reaction_type [catalyst_type] [constraint]'. Example: 'cross_coupling metal cheap'"),
    ]

    output_sig       = [
        ("result", "str", "Recommended catalyst(s) with full details including conditions, pros/cons, and safety."),
    ]

    examples         = [
        {
            "code_input": {"reaction_type": "E2", "catalyst_type": "base", "constraints": "None"},
            "text_input": {"input_text": "E2 base"},
            "output": {
                "result": """## Catalyst Recommendation for: E2 Elimination

### ⭐ Recommended Bases

#### 1. t-BuPotassium tert-butoxide)
| Property | Detail |
|----------|--------|
| Type | Bulky alkoxide (sterically hindered) |
| Strength | Very strong |
| Conditions | t-BuOH or THF, RT to 80°C |

**Why it's great for E2:** Steric bulk prevents SN2 substitution entirely → clean elimination.

✅ **Pros:**
- ★ Best E2 promoter — steric bulk prevents SN2
- Deprotes very weak acids
- Commercially available

❌ **Cons:**
- Expensive
- Very hygroscopic
- t-BuOH has high BP (hard to remove)

⚠️ **Safety:** Pyrophoric when dry! Store under N2.

---

#### 2. DBU (1,8-Diazabicyclo[5.4.0]undec-7-ene)
| Property | Detail |
|----------|--------|
| Type | Amidine base (non-nucleophilic) |
| Strength | Strong (pKa ~12) |
| Conditions | Various organic solvents, RT to reflux |

**Why it works:** Non-nucleophilic → won't compete via SN2. Mild enough for sensitive substrates.

✅ **Pros:**
- Non-nucleophilic
- Liquid (easy to handle)
- Can be used catalytically

❌ **Cons:**
- Expensive
- Not strong enough for some deprotonations

---

### Also Consider
- **LDA**: For E2 where you also want kinetic control (-78°C)
- **NaOEt/KOtBu**: Cheaper alternative but less selective"""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, catalyst_type: str = "any", constraints: str = "None") -> str:
        """Core logic: recommend catalysts."""
        rt = reaction_type.strip().lower()
        ct = catalyst_type.strip().lower() if catalyst_type else "any"
        cons = constraints.strip() if constraints and constraints.upper() != "NONE" else None

        lines = []
        lines.append(f"## Catalyst Recommendation for: {rt}\n")

        # Map reaction type to recommended catalyst keys
        rec_keys = self._get_recommendations(rt, ct)

        if not rec_keys:
            lines.append(f"⚠️ No recommendations found for '{rt}' with catalyst type '{ct}'.")
            lines.append(f"\n**Available reaction types:**\n")
            for k in sorted(_REACTION_CATALYST_MAP.keys()):
                lines.append(f"- `{k}`")
            return "\n".join(lines)

        # Display recommendations
        lines.append("### ⭐ Recommended Catalysts\n")
        for i, (cat_type, cat_name) in enumerate(rec_keys[:5]):
            details = self._get_catalyst_details(cat_type, cat_name)
            if details:
                lines.append(f"#### {i+1}. {cat_name}")
                lines.append(f"| Property | Detail |")
                lines.append(f"|----------|--------|")
                for k, v in details.get("summary", {}).items():
                    lines.append(f"| {k} | {v} |")
                lines.append("")
                if details.get("why"):
                    lines.append(f"**Why it's great:** {details['why']}\n")
                lines.append("✅ **Pros:**")
                for p in details["pros"]:
                    lines.append(f"- {p}")
                lines.append("\n❌ **Cons:**")
                for c in details["cons"]:
                    lines.append(f"- {c}")
                lines.append(f"\n⚠️ **Safety:** {details['safety']}\n")
                lines.append("---\n")

        # Constraint filter note
        if cons:
            lines.append(f"### 📋 Constraint Note: `{cons}`\n")
            filtered = self._apply_constraint(rec_keys, cons)
            if filtered:
                lines.append("**Best matches for this constraint:**")
                for ct_name, cn in filtered[:3]:
                    lines.append(f"- ✅ {cn}")
            else:
                lines.append("Consider relaxing constraint or exploring alternatives.")

        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        rt = parts[0] if parts else "cross_coupling"
        ct = parts[1] if len(parts) > 1 else "any"
        cons = parts[2] if len(parts) > 2 else "None"
        return self._run_base(rt, ct, cons)

    def _get_recommendations(self, rt, ct):
        # Try exact match first
        rec_list = _REACTION_CATALYST_MAP.get(rt)
        if not rec_list:
            # Fuzzy match
            for key in _REACTION_CATALYST_MAP:
                if rt in key or key in rt:
                    rec_list = _REACTION_CATALYST_MAP[key]
                    break
        if not rec_list:
            return []

        result = []
        for item in rec_list:
            item_ct, item_name = item.split(":", 1)
            if ct == "any" or item_ct == ct:
                result.append((item_ct, item_name))
        return result

    def _get_catalyst_details(self, cat_type, cat_name):
        cat_db = _CATALYST_DB.get(cat_type, {})
        catalysts = cat_db.get("catalysts", {})
        cat = catalysts.get(cat_name)
        if not cat:
            # Try fuzzy match
            for k in catalysts:
                if cat_name.lower() in k.lower() or k.lower() in cat_name.lower():
                    cat = catalysts[k]
                    cat_name = k
                    break
        if not cat:
            return None

        return {
            "summary": {
                "Type": cat["type"],
                "Strength": cat.get("strength", "N/A"),
                "Conditions": cat["conditions"].split(".")[0] + ".",
            },
            "why": cat.get("uses", [])[0] if cat.get("uses") else "",
            "pros": cat.get("pros", []),
            "cons": cat.get("cons", []),
            "safety": cat.get("safety", "See SDS."),
        }

    def _apply_constraint(self, rec_keys, constraint):
        cl = constraint.lower()
        results = []
        for ct, cn in rec_keys:
            details = self._get_catalyst_details(ct, cn)
            if not details:
                continue
            if cl == "cheap":
                if any(x in details["safety"].lower() for x in ["cheap", "inexpensive", "low cost"]) or \
                   any("cheap" in p.lower() or "inexpensive" in p.lower() for p in details["pros"]):
                    results.append((ct, cn))
            elif cl == "green":
                if "enzyme" in ct or "iron" in cn.lower() or "biological" in details["safety"].lower():
                    results.append((ct, cn))
            elif cl == "mild_conditions":
                if "RT" in details["summary"]["Conditions"] and "☠️" not in details["safety"][:3]:
                    results.append((ct, cn))
            elif cl == "high_selectivity":
                if any("selective" in p.lower() or "enantio" in p.lower() for p in details["pros"]):
                    results.append((ct, cn))
        return results
