import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive solvent database
_SOLVENT_DB = {
    # Common solvents for organic reactions: (dielectric_const, bp_C, polarity, protic/aprotic, donor_num, acceptor_num, safety, typical_uses)
    "water": {
        "formula": "H2O", "dielectric": 78.4, "bp": 100.0, "polarity": "very high",
        "type": "protic", "donor_num": 33, "acceptor_num": 47.2,
        "safety": "Non-toxic, non-flammable", "cost": "Very low",
        "miscible_with": ["alcohols", "acetone", "acetonitrile", "THF"],
        "not_miscible_with": ["hexane", "toluene", "ether"],
        "uses": ["SN1", "E1", "acid-base", "ionic reactions", "hydrolysis", "aqueous workup"],
        "notes": "Excellent for ionic/polar reactions; promotes carbocation formation (SN1/E1). Poor for organometallics.",
    },
    "methanol": {
        "formula": "CH3OH", "dielectric": 32.6, "bp": 64.7, "polarity": "high",
        "type": "protic", "donor_num": 19, "acceptor_num": 41.3,
        "safety": "Toxic (can cause blindness), flammable", "cost": "Low",
        "miscible_with": ["water", "ethanol", "acetone", "chloroform"],
        "not_miscible_with": ["hexane", "alkanes"],
        "uses": ["SN1", "E1", "solvolysis", "esterification", "Grignard quenching"],
        "notes": "Protic — promotes ionization. Good for SN1/E1. Reacts with strong bases and organometallics.",
    },
    "ethanol": {
        "formula": "C2H5OH", "dielectric": 24.3, "bp": 78.4, "polarity": "moderate-high",
        "type": "protic", "donor_num": 20, "acceptor_num": 39,
        "safety": "Flammable", "cost": "Low",
        "miscible_with": ["water", "organic solvents"],
        "not_miscible_with": ["alkanes"],
        "uses": ["SN1", "E1", "fermentation", "reductions (NaBH4)"],
        "notes": "Similar to methanol but less polar. Good general-purpose protic solvent.",
    },
    "t-butanol": {
        "formula": "(CH3)3COH", "dielectric": 12.5, "bp": 82.0, "polarity": "low-moderate",
        "type": "protic", "donor_num": 25, "acceptor_num": 34,
        "safety": "Flammable", "cost": "Low-Moderate",
        "miscible_with": ["water", "ethanol"],
        "not_miscible_with": ["alkanes"],
        "uses": ["E2 elimination", "SN1 (with acid)"],
        "notes": "Bulky alcohol — favors E2 over SN2 due to steric effects on the solvent itself.",
    },
    "acetic_acid": {
        "formula": "CH3COOH", "dielectric": 6.2, "bp": 118.0, "polarity": "high",
        "type": "protic", "donor_num": 12, "acceptor_num": 44.8,
        "safety": "Corrosive, pungent odor", "cost": "Low",
        "miscible_with": ["water", "ethanol"],
        "not_miscible_with": ["alkanes", "ether"],
        "uses": ["SN1", "E1", "acylation", "electrophilic aromatic substitution", "nitration"],
        "notes": "Acidic medium — excellent for reactions requiring acidic conditions.",
    },
    "DMSO": {
        "formula": "(CH3)2SO", "dielectric": 46.7, "bp": 189.0, "polarity": "very high",
        "type": "aprotic_polar", "donor_num": 29.8, "acceptor_num": 17.7,
        "safety": "Penetrates skin easily (carries dissolved chemicals)", "cost": "Moderate",
        "miscible_with": ["water", "alcohols", "chloroform", "acetone"],
        "not_miscible_with": ["hydrocarbons"],
        "uses": ["SN2", "oxidations (Swern, Pfitzner-Moffatt)", "organometallic reactions", "biological"],
        "notes": "★ BEST SN2 SOLVENT — polar aprotic doesn't H-bond to anions → naked anions are highly reactive. High boiling point allows heated reactions.",
    },
    "DMF": {
        "formula": "HCON(CH3)2", "dielectric": 36.7, "bp": 153.0, "polarity": "high",
        "type": "aprotic_polar", "donor_num": 26.6, "acceptor_num": 16.0,
        "safety": "Reproductive hazard, teratogenic", "cost": "Moderate",
        "miscible_with": ["water", "most organic solvents"],
        "not_miscible_with": ["alkanes", "ether"],
        "uses": ["SN2", "nucleophilic substitutions", "coupling reactions", "Vilsmeier-Haack"],
        "notes": "Excellent polar aprotic for SN2. Can be hard to remove (high BP). Toxicity concerns.",
    },
    "acetone": {
        "formula": "(CH3)2CO", "dielectric": 20.7, "bp": 56.0, "polarity": "moderate-high",
        "type": "aprotic_polar", "donor_num": 17, "acceptor_num": 12.5,
        "safety": "Highly flammable", "cost": "Low",
        "miscible_with": ["water", "alcohols", "chloroform", "ether"],
        "not_miscible_with": ["alkanes"],
        "uses": ["SN2", "oxidation product removal", "halogenation", "nucleophilic addition"],
        "notes": "Good polar aprotic for SN2. Low BP makes it easy to remove. Miscible with water.",
    },
    "acetonitrile": {
        "formula": "CH3CN", "dielectric": 37.5, "bp": 81.6, "polarity": "high",
        "type": "aprotic_polar", "donor_num": 14.1, "acceptor_num": 22.7,
        "safety": "Toxic, releases cyanide if burned", "cost": "Low-Moderate",
        "miscible_with": ["water", "most organic solvents"],
        "not_miscible_with": ["alkanes", "fats"],
        "uses": ["SN2", "SN1", "electrophilic reactions", "HPLC mobile phase", "ion pairing"],
        "notes": "Polar aprotic — good for SN2. Inert to many reagents. Easy to remove (moderate BP).",
    },
    "THF": {
        "formula": "C4H8O", "dielectric": 7.5, "bp": 66.0, "polarity": "moderate",
        "type": "aprotic_polar", "donor_num": 20, "acceptor_num": 8.0,
        "safety": "Forms peroxides on storage (explosive!)", "cost": "Moderate",
        "miscible_with": ["water", "ether", "alcohols", "benzene"],
        "not_miscible_with": ["alkanes"],
        "uses": ["Grignard", "organolithium", "hydride reductions", "polymerization", "SN2"],
        "notes": "Essential for organometallic chemistry. Ethers solvate cations well. MUST check for peroxides before distilling!",
    },
    "diethyl_ether": {
        "formula": "(C2H5)2O", "dielectric": 4.3, "bp": 34.6, "polarity": "low",
        "type": "aprotic_polar", "donor_num": 19.2, "acceptor_num": 3.9,
        "safety": "Highly flammable, forms peroxides", "cost": "Low",
        "miscible_with": ["most organic solvents"],
        "not_miscible_with": ["water"],
        "uses": ["Grignard", "Wittig", "extractions", "liquid-liquid extraction"],
        "notes": "Classic ether for Grignard/Wittig. Very low BP. Forms explosive peroxides — use BHT-stabilized or fresh.",
    },
    "dichloromethane": {
        "formula": "CH2Cl2", "dielectric": 8.9, "bp": 40.0, "polarity": "moderate",
        "type": "aprotic_polar", "donor_num": 0, "acceptor_num": 8.9,
        "safety": "Suspected carcinogen, volatile", "cost": "Low",
        "miscible_with": ["most organic solvents"],
        "not_miscible_with": ["water"],
        "uses": ["extraction", "chromatography", "many organic reactions as inert medium"],
        "notes": "Inert to most conditions. Excellent extraction solvent (immiscible with water). Low BP for easy removal.",
    },
    "chloroform": {
        "formula": "CHCl3", "dielectric": 4.8, "bp": 61.2, "polarity": "low-moderate",
        "type": "aprotic_polar", "donor_num": 0, "acceptor_num": 8.4,
        "safety": "Suspected carcinogen, liver toxin", "cost": "Low",
        "miscible_with": ["most organic solvents", "ethanol"],
        "not_miscible_with": ["water"],
        "uses": ["extraction", "chromatography", "phase-transfer catalysis"],
        "notes": "Slightly more polar than DCM. Often used interchangeably with DCM.",
    },
    "toluene": {
        "formula": "C6H5CH3", "dielectric": 2.38, "bp": 110.6, "polarity": "very low",
        "type": "aprotic_nonpolar", "donor_num": 0, "acceptor_num": 1.4,
        "safety": "Flammable, neurotoxin", "cost": "Low",
        "miscible_with": ["organic solvents"],
        "not_miscible_with": ["water"],
        "uses": ["Diels-Alder", "[3+2] cycloadditions", "high-temp reactions", "reflux"],
        "notes": "Nonpolar — good for pericyclic reactions and thermal reactions. Higher BP than benzene (safer alternative).",
    },
    "hexane": {
        "formula": "C6H14", "dielectric": 1.88, "bp": 69.0, "polarity": "very low",
        "type": "aprotic_nonpolar", "donor_num": 0, "acceptor_num": 0,
        "safety": "Highly flammable, neurotoxic", "cost": "Low",
        "miscible_with": ["nonpolar organics"],
        "not_miscible_with": ["water", "alcohols", "DMSO"],
        "uses": ["chromatography (nonpolar eluent)", "extraction of nonpolar compounds", "recrystallization"],
        "notes": "Very nonpolar. Used in column chromatography and for extracting nonpolar compounds.",
    },
    "pyridine": {
        "formula": "C5H5N", "dielectric": 12.4, "bp": 115.3, "polarity": "moderate",
        "type": "aprotic_polar", "donor_num": 33.1, "acceptor_num": 6.7,
        "safety": "Foul odor, toxic", "cost": "Moderate",
        "miscible_with": ["water", "most organic solvents"],
        "not_miscible_with": ["alkanes"],
        "uses": ["base catalyst", "acylation", "silyl deprotection", "as base + solvent dual role"],
        "notes": "Basic aromatic amine — acts as both solvent AND base. Essential for acylation reactions.",
    },
    "HMPA": {
        "formula": "[(CH3)2N]3PO (HMPA)", "dielectric": 30.0, "bp": 235.0, "polarity": "high",
        "type": "aprotic_polar", "donor_num": 38.8, "acceptor_num": 10.6,
        "safety": "Probable carcinogen, suspected reproductive toxin", "cost": "High",
        "miscible_with": ["water", "most organic solvents"],
        "not_miscible_with": ["alkanes"],
        "uses": ["SN2 enhancement", "anion activation", "enhancing nucleophilicity"],
        "notes": "Most powerful polar aprotic for 'naked' anions. TOXIC — avoid when possible. DMPU is safer alternative.",
    },
}

# Reaction type → recommended solvents mapping
_REACTION_SOLVENT_MAP = {
    "SN2": {
        "primary": ["DMSO", "acetone", "acetonitrile", "DMF"],
        "secondary": ["THF", "DMF", "DMSO"],
        "avoid": ["protic solvents (slow down nucleophile)"],
        "reason": "Polar aprotic solvents don't H-bond to nucleophiles, making them highly reactive ('naked anions').",
    },
    "SN1": {
        "primary": ["water", "methanol", "ethanol", "acetic acid"],
        "secondary": ["acetone/water mixtures", "TFA"],
        "avoid": ["strong bases in aprotic solvents (promote SN2 instead)"],
        "reason": "Protic solvents stabilize carbocation intermediates and promote ionization of the leaving group.",
    },
    "E2": {
        "primary": ["t-butanol", "ethanol", "DMSO", "THF"],
        "secondary": ["t-BuOK/t-BuOH system"],
        "avoid": ["polar aprotic (may favor SN2)"],
        "reason": "Bulky bases in alcohols favor elimination. Heat promotes E2.",
    },
    "E1": {
        "primary": ["ethanol", "water", "acetic acid"],
        "secondary": ["formic acid", "TFA"],
        "avoid": ["strong bases"],
        "reason": "Requires weak base/nucleophile + heat in protic solvent to allow carbocation formation.",
    },
    "Grignard": {
        "primary": ["diethyl_ether", "THF"],
        "secondary": [],
        "avoid": ["protic solvents (quench Grignard!)", "anything with acidic protons"],
        "reason": "Ether oxygen coordinates to Mg, stabilizing the Grignard reagent. Must be absolutely dry.",
    },
    "organolithium": {
        "primary": ["THF", "diethyl_ether", "hexane/ether mixtures"],
        "secondary": [],
        "avoid": ["any protic or electrophilic solvent"],
        "reason": "Even more reactive than Grignard. THF solvates Li+ well. Must be anhydrous, inert atmosphere.",
    },
    "reduction_NaBH4": {
        "primary": ["methanol", "ethanol", "THF", "water"],
        "secondary": [],
        "avoid": ["strong acids"],
        "reason": "NaBH4 is stable in alcohols and protic solvents. Methanol/EtOH are common.",
    },
    "reduction_LiAlH4": {
        "primary": ["diethyl_ether", "THF"],
        "secondary": [],
        "avoid": ["ANY protic solvent (explosive reaction!)", "chlorinated solvents"],
        "reason": "LiAlH4 reacts violently with any proton source. Must be strictly anhydrous ether/THF.",
    },
    "oxidation_Swern": {
        "primary": ["DMSO", "dichloromethane"],
        "secondary": [],
        "avoid": ["protic solvents"],
        "reason": "DMSO is both reactant and solvent. DCM as co-solvent. Must be anhydrous, cold (-78°C to RT).",
    },
    "oxidation_PCC": {
        "primary": ["dichloromethane"],
        "secondary": [],
        "avoid": ["protic solvents (decompose PCC)"],
        "reason": "PDC/PCC oxidants are used in DCM. Non-aqueous conditions prevent over-oxidation to carboxylic acids.",
    },
    "oxidation_Jones": {
        "primary": ["acetone"],
        "secondary": [],
        "avoid": ["solvents that react with Cr(VI)"],
        "reason": "Jones reagent (CrO3/H2SO4) in acetone. Acetone is resistant to oxidation.",
    },
    "Wittig": {
        "primary": ["THF", "diethyl_ether"],
        "secondary": ["DMF", "DMSO"],
        "avoid": ["protic solvents (destroy ylide)"],
        "reason": "Phosphorus ylide requires anhydrous conditions. THF is standard.",
    },
    "Diels-Alder": {
        "primary": ["toluene", "benzene", "xylene", "dichloromethane"],
        "secondary": ["without solvent (neat)"],
        "avoid": ["polar protic (interfere with orbital alignment)"],
        "reason": "Nonpolar or moderately polar solvents. Often run neat at high temperature.",
    },
    "Friedel-Crafts_alkylation": {
        "primary": ["dichloromethane", "nitrobenzene", "CS2"],
        "secondary": ["nitromethane"],
        "avoid": ["protic solvents (deactivate Lewis acid)"],
        "reason": "Lewis acid catalyst (AlCl3, FeCl3) requires non-coordinating, anhydrous solvent.",
    },
    "Friedel-Crafts_acylation": {
        "primary": ["dichloromethane", "nitrobenzene", "CS2"],
        "secondary": [],
        "avoid": ["coordinating solvents (compete with substrate for Lewis acid)"],
        "reason": "Same as alylation — Lewis acid catalysis needs non-basic, anhydrous media.",
    },
    "nitration": {
        "primary": ["sulfuric_acid (as solvent)", "acetic_acid"],
        "secondary": ["dichloromethane (mixed acid)"],
        "avoid": ["bases (neutralize nitrating agent)"],
        "reason": "Conc. HNO3/H2SO4 mixture. Sulfuric acid acts as solvent and dehydrating agent.",
    },
    "sulfonation": {
        "primary": ["sulfuric_acid", "oleum"],
        "secondary": ["SO3 in dioxane"],
        "avoid": ["bases, water"],
        "reason": "Fuming sulfuric acid is both reagent and solvent.",
    },
    "amide_coupling": {
        "primary": ["DMF", "dichloromethane", "THF"],
        "secondary": ["acetonitrile", "NMP"],
        "avoid": ["protic solvents (interfere with coupling reagents)"],
        "reason": "EDC/HOBt, DCC, HATU couplings need aprotic conditions. DMF dissolves everything.",
    },
    "hydrogenation": {
        "primary": ["ethanol", "methanol", "ethyl_acetate", "hexane"],
        "secondary": ["THF", "acetic_acid"],
        "avoid": ["compounds with catalyst poisons (S, P, halides)"],
        "reason": "Pd/C, PtO2 catalyst. Solvent must dissolve substrate and be H2-compatible.",
    },
    "halogenation_addition": {
        "primary": ["dichloromethane", "chloroform", "CCl4"],
        "secondary": [],
        "avoid": ["protic solvents (may cause side reactions)"],
        "reason": "Br2, Cl2 additions to alkenes. Inert halogenated solvents are ideal.",
    },
    "hydroboration": {
        "primary": ["THF", "diglyme"],
        "secondary": ["diglyme"],
        "avoid": ["protic solvents (destroy borane)"],
        "reason": "BH3·THF complex is standard. THF stabilizes borane adduct.",
    },
    "suzuki_coupling": {
        "primary": ["toluene/ethanol/water mixture", "dioxane/water", "DMF/water"],
        "secondary": ["THF/water"],
        "avoid": ["strictly anhydrous conditions (needs water for transmetalation)"],
        "reason": "Pd-catalyzed cross-coupling. Requires aqueous base (K2CO3, Cs2CO3) for transmetalation.",
    },
    "Heck_reaction": {
        "primary": ["DMF", "acetonitrile", "NMP", "dioxane"],
        "secondary": [],
        "avoid": ["protic solvents (in some cases)"],
        "reason": "Pd-catalyzed coupling. Polar aprotic at elevated temperature (80-120°C).",
    },
    "alkyne_chemistry": {
        "primary": ["THF", "diethyl_ether", "liquid ammonia", "hexane"],
        "secondary": [],
        "avoid": ["protic unless intentional (e.g., metal-ammonia reduction)"],
        "reason": "Depends on specific transformation. Terminal alkynes need basic conditions.",
    },
    "enolate_chemistry": {
        "primary": ["THF", "DMF", "DMSO", "diethyl_ether"],
        "secondary": ["toluene (for LDA at low temp)"],
        "avoid": ["protic solvents (quench enolate)"],
        "reason": "LDA, NaH, KHMDS generate enolates. Must be anhydrous. THF/DMF common.",
    },
    "peptide_synthesis": {
        "primary": ["DMF", "dichloromethane", "NMP"],
        "secondary": ["acetonitrile"],
        "avoid": ["protic solvents (cause side reactions)"],
        "reason": "SPPS uses DMF/NMP for swelling resin and coupling. DCM for washes.",
    },
}


@ChemMCPManager.register_tool
class SolventSelector(BaseTool):
    """
    根据反应类型推荐合适溶剂的工具。
    内置全面的溶剂数据库和反应类型映射，基于极性、质子性、配位能力、安全性等维度进行推荐。
    """
    __version__      = "0.1.0"
    name             = "SolventSelector"
    func_name        = "select_solvent"
    description      = "Recommend appropriate solvents based on reaction type, with detailed properties, safety notes, and reasoning."
    implementation_description = "Uses comprehensive embedded database of 18+ common organic solvents with properties (dielectric constant, BP, polarity, type, safety) and maps 28+ reaction types to optimal solvent choices."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Solvent", "Reaction Conditions", "Organic Chemistry", "Lab Safety"]
    required_envs    = []

    code_input_sig   = [
        ("reaction_type", "str", "N/A", "Type of reaction (e.g., 'SN2', 'SN1', 'Grignard', 'reduction', 'oxidation', 'Diels-Alder', 'suzuki_coupling')."),
        ("constraints", "str", "None", "Optional constraints: 'low_toxicity', 'high_bp', 'easy_removal', 'water_miscible', 'anhydrous_only'. Use 'None' for no constraints."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'reaction_type [constraint]'. Example: 'SN2 low_toxicity'"),
    ]

    output_sig       = [
        ("result", "str", "Recommended solvent(s) with full details table, reasoning, and safety warnings."),
    ]

    examples         = [
        {
            "code_input": {"reaction_type": "SN2", "constraints": "None"},
            "text_input": {"input_text": "SN2"},
            "output": {
                "result": """## Solvent Recommendation for: SN2

### ⭐ Primary Recommendations

| # | Solvent | Formula | Dielectric | BP (°C) | Type | Safety |
|---|---------|---------|------------|---------|------|--------|
| 1 | **DMSO** | (CH₃)₂SO | 46.7 | 189.0 | Polar Aprotic | ⚠️ Penetrates skin |
| 2 | **Acetone** | (CH₃)₂CO | 20.7 | 56.0 | Polar Aprotic | 🔥 Flammable |
| 3 | **Acetonitrile** | CH₃CN | 37.5 | 81.6 | Polar Aprotic | ☠️ Toxic |

### Why These Solvents?
Polar aprotic solvents don't hydrogen-bond to nucleophiles, creating 'naked' anions that are highly reactive toward SN2 displacement.

### 🚫 Avoid
- Protic solvents (water, alcohols) — H-bond to nucleophile, slowing it dramatically
- Nonpolar solvents — don't dissolve ionic intermediates

### Key Properties Comparison
- **DMSO**: Best SN2 performance, but high BP (hard to remove), skin-penetrating hazard
- **Acetone**: Easy to remove (BP 56°C), cheap, miscible with water
- **Acetonitrile**: Good balance, moderate BP, inert to many reagents"""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, constraints: str = "None") -> str:
        """Core logic: recommend solvents for given reaction type."""
        rt = reaction_type.strip().lower()
        cons = constraints.strip() if constraints and constraints.upper() != "NONE" else None

        lines = []
        lines.append(f"## Solvent Recommendation for: {reaction_type}\n")

        # Find matching reaction type (fuzzy match)
        match_key = self._match_reaction(rt)
        rec = _REACTION_SOLVENT_MAP.get(match_key)

        if rec is None:
            lines.append(f"⚠️ Reaction type '{reaction_type}' not found in database.")
            lines.append(f"\n**Available reaction types:**\n")
            for k in sorted(_REACTION_SOLVENT_MAP.keys()):
                lines.append(f"- `{k}`")
            return "\n".join(lines)

        # Primary recommendations
        lines.append("### ⭐ Primary Recommendations\n")
        primary = rec["primary"]
        table_rows = []
        for i, solv_name in enumerate(primary):
            s = _SOLVENT_DB.get(solv_name)
            if s:
                safety_icon = "✅" if "non-toxic" in s["safety"].lower() or "Low" in s.get("cost","") else \
                             ("⚠️" if any(x in s["safety"].lower() for x in ["toxic", "penetrates", "carcinogen"]) else "🔥")
                table_rows.append(
                    f"| {i+1} | **{solv_name.replace('_', ' ').title()}** | {s['formula']} | "
                    f"{s['dielectric']} | {s['bp']} | {s['type'].replace('_', ' ').title()} | {safety_icon} {s['safety']} |"
                )

        if table_rows:
            lines.append("| # | Solvent | Formula | Dielectric | BP (°C) | Type | Safety |")
            lines.append("|---|---------|---------|------------|---------|------|--------|")
            lines.extend(table_rows)
            lines.append("")

        # Reasoning
        lines.append(f"### Why These Solvents?\n{rec['reason']}\n")

        # Secondary recommendations
        if rec.get("secondary"):
            lines.append("### Alternative Options\n")
            for s in rec["secondary"]:
                lines.append(f"- **{s.replace('_', ' ').title()}**")
            lines.append("")

        # Avoid section
        if rec.get("avoid"):
            lines.append("### 🚫 Avoid\n")
            for a in rec["avoid"]:
                lines.append(f"- {a}")
            lines.append("")

        # Detailed property comparison for top picks
        lines.append("### Key Properties Comparison\n")
        for solv_name in primary[:3]:
            s = _SOLVENT_DB.get(solv_name)
            if s:
                notes_key = "notes"
                lines.append(f"**{solv_name.replace('_', ' ').title()}**: {s[notes_key]}")
                lines.append("")

        # Apply constraints filter
        if cons:
            lines.append("---\n### Constraint Filter: `{cons}`\n")
            filtered = self._apply_constraints(primary + rec.get("secondary", []), cons)
            if filtered:
                lines.append("**Filtered recommendations:**")
                for f in filtered:
                    lines.append(f"- ✅ {f.replace('_', ' ').title()}")
            else:
                lines.append("⚠️ No solvents fully meet this constraint. Consider relaxing requirements.")

        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        rt = parts[0] if parts else "SN2"
        cons = parts[1] if len(parts) > 1 else "None"
        return self._run_base(rt, cons)

    def _match_reaction(self, rt):
        exact = _REACTION_SOLVENT_MAP.get(rt)
        if exact:
            return rt
        # Fuzzy match
        for key in _REACTION_SOLVENT_MAP:
            if rt in key or key in rt:
                return key
        # Partial keyword match
        keywords = {
            "sn2": "SN2", "sn1": "SN1", "e2": "E2", "e1": "E1",
            "grignard": "Grignard", "organolithium": "organolithium",
            "reduction": "reduction_NaBH4", "oxidation": "oxidation_PCC",
            "wittig": "Wittig", "diels": "Diels-Alder", "diels-alder": "Diels-Alder",
            "friedel": "Friedel-Crafts_alkylation", "nitration": "nitration",
            "suzuki": "suzuki_coupling", "heck": "Heck_reaction",
            "hydrogenation": "hydrogenation", "hydroboration": "hydroboration",
            "enolate": "enolate_chemistry", "peptide": "peptide_synthesis",
            "amide": "amide_coupling", "halogenation": "halogenation_addition",
            "alkyne": "alkyne_chemistry", "sulfonation": "sulfonation",
        }
        for kw, val in keywords.items():
            if kw in rt:
                return val
        return rt

    def _apply_constraints(self, solvent_list, constraint):
        result = []
        cl = constraint.lower()
        for sn in solvent_list:
            s = _SOLVENT_DB.get(sn)
            if not s:
                continue
            ok = True
            if cl == "low_toxicity":
                if any(x in s["safety"].lower() for x in ["carcinogen", "teratogenic", "reproductive"]):
                    ok = False
            elif cl == "high_bp":
                if s["bp"] < 100:
                    ok = False
            elif cl == "easy_removal":
                if s["bp"] > 120:
                    ok = False
            elif cl == "water_miscible":
                if "water" not in s.get("miscible_with", []):
                    ok = False
            elif cl == "anhydrous_only":
                if s["type"] == "protic":
                    ok = False
            if ok:
                result.append(sn)
        return result
