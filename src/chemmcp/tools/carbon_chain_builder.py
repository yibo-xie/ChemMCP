"""
Carbon Chain Builder - suggests carbon chain elongation and shortening strategies
with specific reaction sequences, reagents, and conditions.
"""

import logging
import re
from typing import Dict, List, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# === CHAIN ELONGATION METHODS (carbon count increase) ===
ELONGATION_METHODS: List[Dict[str, Any]] = [
    {
        "name": "Wittig Reaction (+1 C from aldehyde)",
        "carbon_change": +1,
        "starting_fg": ["aldehyde"],
        "product_fg": ["alkene"],
        "reagents": "Ph₃P=CHR (ylide), THF",
        "conditions": "0°C → RT, under N₂, 2-12 h",
        "mechanism": "[2+2] cycloaddition → oxaphosphetane → Ph₃PO + alkene",
        "yield": "50-90%",
        "stereochemistry": "Z (non-stabilized ylide) / E (stabilized ylide)",
        "example": "R-CHO + Ph₃P=CH₂ → R-CH=CH₂",
        "notes": "One of the most reliable alkene-forming reactions; byproduct Ph₃PO easily removed",
    },
    {
        "name": "Grignard Addition (+n C from carbonyl)",
        "carbon_change": None,  # depends on R group size
        "starting_fg": ["aldehyde", "ketone", "ester", "CO2", "epoxide"],
        "product_fg": ["alcohol (1°/2°/3°)", "carboxylic acid", "alcohol"],
        "reagents": "R'-MgX (Grignard) or R'-Li (organolithium), ether/THF",
        "conditions": "-10°C → RT, anhydrous, under N₂/Ar",
        "mechanism": "Nucleophilic addition to carbonyl carbon (or epoxide ring-opening)",
        "yield": "60-90%",
        "stereochemistry": "New stereocenter formed (racemic unless chiral auxiliary used)",
        "example": "R-CHO + CH₃MgBr → R-CH(OH)-CH₃ (+1 C); R-CHO + PhMgBr → R-CH(OH)-Ph (+6 C)",
        "notes": "Most versatile chain-building method; R' can be 1°-3° alkyl, vinyl, arynyl; sensitive to protic solvents and active H",
    },
    {
        "name": "Cyanohydrin Formation / Cyanation (+1 C)",
        "carbon_change": +1,
        "starting_fg": ["aldehyde", "ketone", "alkyl_halide"],
        "product_fg": ["cyanohydrin", "cyanohydrin", "nitrile"],
        "reagents": "KCN/NaCN (for carbonyls) or NaCN/DMSO (for halides)",
        "conditions": "RT (carbonyl) or reflux (SN2 cyanation)",
        "mechanism": "Nucleophilic addition of CN⁻ to C=O; SN2 displacement for halides",
        "yield": "65-90% (carbonyl); 70-90% (halide, 1° only)",
        "stereochemistry": "New chiral center at cyanohydrin carbon (racemic)",
        "example": "R-CHO + KCN → R-CH(OH)-CN → (hydrolysis) → R-CH(OH)-COOH",
        "notes": "CN can be further converted to COOH (hydrolysis), CH₂NH₂ (reduction), CHO (DIBAL-H), etc.",
    },
    {
        "name": "Aldol Condensation (+2 to +n C)",
        "carbon_change": +2,
        "starting_fg": ["aldehyde", "ketone"],
        "product_fg": ["α,β-unsaturated carbonyl", "β-hydroxy carbonyl"],
        "reagents": "Base (NaOH, LDA) or acid catalyst",
        "conditions": "0°C → RT (kinetic) or heat (thermodynamic dehydration)",
        "mechanism": "Enolate formation → nucleophilic addition to carbonyl → (dehydration)",
        "yield": "50-85% (crossed-aldol can be lower)",
        "stereochemistry": "Possible syn/anti selectivity with chiral auxiliaries or Evans oxazolidinones",
        "example": "Ph-CHO + CH₃CHO (LDA) → Ph-CH(OH)-CH(CH₃)-OH → (-H₂O) → Ph-CH=CH-CHO",
        "notes": "Crossed-aldol requires careful control (one component non-enolizable ideal); directed aldol with boron enolates gives high stereocontrol",
    },
    {
        "name": "Arndt-Eistert Homologation (+1 C, acid → acid)",
        "carbon_change": +1,
        "starting_fg": ["carboxylic_acid"],
        "product_fg": ["carboxylic_acid (homologated)"],
        "reagents": "1) SOCl₂ 2) CH₂N₂ 3) Ag₂O/H₂O or AgBF₄/H₂O or hv/heat",
        "conditions": "Multi-step: acid chloride → diazoketone → Wolff rearrangement → homologated acid",
        "mechanism": "Acid chloride + CH₂N₂ → α-diazoketone → Wolff rearrangement (1,2-shift) → ketene hydrolysis",
        "yield": "50-75% overall",
        "stereochemistry": "Retention at migrating carbon (concerted rearrangement)",
        "example": "R-COOH → R-COCl → R-CO-CHN₂ → R-CH₂-COOH",
        "notes": "Classic carboxylic acid homologation; CH₂N₂ is hazardous (explosive) — use with care; safer alternatives exist (e.g., Rosemund-von Braun variant)",
    },
    {
        "name": "Homologation via Sulfur Ylide (+1 C)",
        "carbon_change": +1,
        "starting_fg": ["aldehyde", "ketone"],
        "product_fg": ["epoxide (from aldehyde)", "epoxide/oxirane (from ketone)"],
        "reagents": "Dimethylsulfonium methylide (Me₃S⁺I⁻, base) or dimethylsulfoxonium methylide",
        "conditions": "DMSO or THF, NaH or KOtBu, 0°C → RT",
        "mechanism": "Sulfur ylide nucleophilic addition to carbonyl → epoxide formation",
        "yield": "55-80%",
        "stereochemistry": "Trans-selective from acyclic aldehydes (with appropriate ylide)",
        "example": "R-CHO + Me₃S⁺I⁻/base → R-CH---CH₂ (epoxide) → (ring-open) → R-CH(OH)-CH₂OH",
        "notes": "Corey-Chaykovsky reaction; epoxides can be opened regioselectively to give 1,2-difunctional products",
    },
    {
        "name": "Malonic Ester Synthesis (+2 C)",
        "carbon_change": +2,
        "starting_fg": ["alkyl_halide"],
        "product_fg": ["carboxylic_acid (substituted acetic acid)"],
        "reagents": "CH₂(COOEt)₂, EtONa, then RX, then hydrolysis + decarboxylation",
        "conditions": "1) EtONa/EtOH 2) RX 3) H₃O⁺, heat (decarboxylation)",
        "mechanism": "Enolate alkylation of malonic ester → hydrolysis → β-keto acid → decarboxylation",
        "yield": "50-80% overall (3 steps)",
        "stereochemistry": "None (prochiral center may become chiral if R ≠ H)",
        "example": "RX + CH₂(COOEt)₂ → R-CH(COOEt)₂ → R-CH₂-COOH",
        "notes": "Classic method for α-substituted carboxylic acids; double alkylation possible; acetoacetic ester variant gives ketones",
    },
    {
        "name": "Suzuki-Miyaura Coupling (aryl/vinyl +n C)",
        "carbon_change": None,  # variable
        "starting_fg": ["vinyl_halide", "aryl_halide", "triflate"],
        "product_fg": ["biaryl", "conjugated diene/stilbene"],
        "reagents": "R'-B(OH)₂ (boronic acid), Pd(PPh₃)₄ or Pd(dppf)Cl₂, base (K₂CO₃, Cs₂CO₃)",
        "conditions": "Reflux, dioxane/H₂O or toluene/EtOH/H₂O, under N₂, 60-100°C",
        "mechanism": "Oxidative addition → transmetallation → reductive elimination (Pd catalytic cycle)",
        "yield": "70-95%",
        "stereochemistry": "Retains configuration of vinyl groups",
        "example": "Ar-Br + Ph-B(OH)₂ → Ar-Ph (biaryl); vinyl-Br + Ph-B(OH)₂ → styrene derivative",
        "notes": "Extremely reliable and tolerant of many functional groups; widely used in pharma and materials; boronic acids are stable and commercially available",
    },
    {
        "name": "Heck Reaction (+n C, alkene coupling)",
        "carbon_change": None,
        "starting_fg": ["vinyl_halide", "aryl_halide"],
        "product_fg": ["substituted_alkene"],
        "reagents": "Alkene (terminal or internal), Pd(OAc)₂, base (Et₃N, K₂CO₃), possibly phosphine ligand",
        "conditions": "80-130°C, DMF or acetonitrile or dioxane, under N₂",
        "mechanism": "Pd(0) oxidative addition → alkene coordination/migratory insertion → β-H elimination",
        "yield": "50-85%",
        "stereochemistry": "Typically trans (E) selective for terminal alkenes",
        "example": "Ar-I + CH₂=CH-R → Ar-CH=CH-R (E-isomer major)",
        "notes": "Forms new C-C bond between aryl/vinyl halide and alkene; no need for organometallic reagent; intramolecular version powerful for cyclization",
    },
    {
        "name": "Diels-Alder Cycloaddition (+4 C from diene)",
        "carbon_change": +4,
        "starting_fg": ["diene", "dienophile (alkene/alkyne)"],
        "product_fg": ["cyclohexene", "cyclohexadiene"],
        "reagents": "Heat (or Lewis acid catalyst like AlCl₃, EtAlCl₂)",
        "conditions": "Heat (reflux, toluene/xylene) or RT-Lewis acid catalysis",
        "mechanism": "[4+2] concerted cycloaddition (pericyclic, suprafacial on both components)",
        "yield": "40-95% (highly substrate-dependent)",
        "stereochemistry": "Endo rule (kinetic product) vs exo (thermodynamic); endo usually favored under kinetic control",
        "example": "Butadiene + ethylene → cyclohexene; Butadiene + maleic anhydride → endo-bicyclo[2.2.1]hept-5-ene-2,3-dicarboxylic anhydride",
        "notes": "Powerful method for building 6-membered rings with up to 4 new stereocenters in one step; inverse electron demand DA also possible",
    },
]

# === CHAIN SHORTENING METHODS ===
SHORTENING_METHODS: List[Dict[str, Any]] = [
    {
        "name": "Hofmann Elimination (-1 C per cycle)",
        "carbon_change": -1,
        "starting_fg": ["amine (quaternary ammonium)"],
        "product_fg": ["alkene"],
        "reagents": "Excess CH₃I (to form quat salt), then Ag₂O/H₂O, heat",
        "conditions": "1) Quaternization: RT 2) Elimination: 100-150°C",
        "mechanism": "E2 elimination (Hofmann rule: less substituted alkene favored due to sterics of bulky NMe₃⁺)",
        "yield": "30-70% per degradation cycle",
        "stereochemistry": "Hofmann selectivity (less substituted alkene)",
        "example": "R-CH₂-CH₂-NH₂ → (excess CH₃I) → [R-CH₂-CH₂-NMe₃]⁺I⁻ → (Ag₂O, Δ) → R-CH=CH₂ + NMe₃",
        "notes": "Classical degradation method; each cycle removes one carbon; useful for structure determination but low-yielding for synthesis",
    },
    {
        "name": "Oxidative Cleavage of Alkenes (variable C loss)",
        "carbon_change": -2,  # typical for terminal alkene
        "starting_fg": ["alkene"],
        "product_fg": ["ketone/aldehyde", "carboxylic_acid/ketone", "carboxylic_acid"],
        "reagents": "Ozonolysis (O₃ then Zn/HOAc or Me₂S) or KMnO₄ (hot) or OsO₄/NaIO₄",
        "conditions": "O₃: -78°C in CH₂Cl₂ then workup; KMnO₄: reflux, basic or acidic",
        "mechanism": "1,3-dipolar cycloaddition of O₃ → molozide → ozonide → reductive/oxidative workup",
        "yield": "65-95% (ozonolysis); 30-80% (KMnO₄)",
        "stereochemistry": "Destroys stereochemistry at the double bond",
        "example": "R-CH=CH₂ → (O₃/Zn) → R-CHO + HCHO; R-CH=CH-R' → (O₃) → R-CHO + R'-CHO; (KMnO₄ hot) → R-COOH + R'-COOH",
        "notes": "Ozonolysis is cleanest; choice of reductant (Zn vs Me₂S vs PPh₃) determines whether aldehydes or alcohols are obtained; hot KMnO₄ over-oxidizes to acids",
    },
    {
        "name": "Decarboxylation (-1 C)",
        "carbon_change": -1,
        "starting_fg": ["carboxylic_acid (β-keto or β-diacid or malonic type)"],
        "product_fg": ["ketone", "acid", "alkane"],
        "reagents": "Heat (with or without acid/base catalyst)",
        "conditions": "150-200°C (thermal) or milder with Cu salts or lead tetraacetate (HTAD/ Barton decarboxylation)",
        "mechanism": "Six-membered cyclic transition state (β-keto acids) or radical pathway (Barton)",
        "yield": "50-90% (β-keto); 40-80% (Barton)",
        "stereochemistry": "Racemization at α-position (via enol intermediate for thermal)",
        "example": "R-CO-CH₂-COOH → (Δ) → R-CO-CH₃ + CO₂; R-CH₂-COOH (malonic) → R-CH₃ + CO₂",
        "notes": "β-Keto acids decarboxylate readily upon heating; malonic ester synthesis relies on this; Hunsdiecker reaction (Ag salt + Br₂) gives R-Br (loss of 1 C)",
    },
    {
        "name": "Hunsdiecker Reaction (-1 C, acid → halide)",
        "carbon_change": -1,
        "starting_fg": ["carboxylic_acid (as silver salt)"],
        "product_fg": ["alkyl_halide (bromide)"],
        "reagents": "Ag-salt of RCOOH + Br₂ (or I₂, or (PhIO)₂/I₂ — Kochi modification)",
        "conditions": "CCl₄ or other inert solvent, reflux or hv (light)",
        "mechanism": "Radical chain: RCOO-Ag + Br₂ → RCOOBr → R• + CO₂ + Br• → RBr",
        "yield": "30-70%",
        "stereochemistry": "Racemization (radical process)",
        "example": "R-COOH → (AgNO₃/NH₃) → RCOOAg → (+ Br₂) → R-Br + CO₂ + AgBr",
        "notes": "Classical radical decarboxylative halogenation; Kochi modification uses Pb(OAc)₄/LiCl or (PhIO)₂/I₂ for better yields; Cristol–Firth modification uses Hg salts",
    },
    {
        "name": "Kochi Decarboxylation (-1 C, acid → alkane)",
        "carbon_change": -1,
        "starting_fg": ["carboxylic_acid"],
        "product_fg": ["alkane"],
        "reagents": "(Pb(OAc)₄ + LiCl) or (PhI(OCOCF₃)₂ + Cu catalyst) + Hünig's base",
        "conditions": "Benzene or arene solvent, reflux, under N₂",
        "mechanism": "Lead tetraacetate or hypervalent iodine-promoted radical decarboxylation",
        "yield": "40-75%",
        "stereochemistry": "Racemization (radical)",
        "example": "R-COOH → (Pb(OAc)₄, LiCl, hν) → R-Cl + CO₂ (or R-H if no halide source)",
        "notes": "Improved over classical Hunsdiecker; modern variants use photoredox catalysis for milder conditions",
    },
    {
        "name": "Periodate Cleavage of Vicinal Diols (-2 C total, diol → 2 carbonyls)",
        "carbon_change": -2,
        "starting_fg": ["vicinal_diol (1,2-diol)"],
        "product_fg": ["aldehyde/ketone (two fragments)"],
        "reagents": "NaIO₄ (sodium periodate), H₂O/THF or aq. acetone",
        "conditions": "RT, 30 min - 6 h, pH ~7",
        "mechanism": "Cyclic periodate ester intermediate → C-C bond cleavage",
        "yield": "80-98%",
        "stereochemistry": "Non-stereospecific (cleaves bond regardless of stereochemistry)",
        "example": "R-CH(OH)-CH(OH)-R' → (NaIO₄) → R-CHO + R'-CHO (or ketone if R/R' = alkyl)",
        "notes": "Very clean and high-yielding; complementary to Pb(OAc)₄ oxidation; often used after Upjohn dihydroxylation of alkenes for two-step cleavage",
    },
    {
        "name": "Ozone Workup to Carboxylic Acids (oxidative)",
        "carbon_change": -2,  # terminal alkene loses 1 C as CO₂
        "starting_fg": ["alkene"],
        "product_fg": ["carboxylic_acid", "ketone"],
        "reagents": "O₃ then H₂O₂ (oxidative workup)",
        "conditions": "-78°C (ozonolysis) then add H₂O₂, warm to RT",
        "mechanism": "Ozonolysis → ozonide → oxidative cleavage to acid(s)/ketone(s)",
        "yield": "60-85%",
        "stereochemistry": "Bond cleaved, stereochemistry lost",
        "example": "R-CH=CH₂ → (O₃, then H₂O₂) → R-COOH + CO₂ + H₂O",
        "notes": "Terminal alkene gives carboxylic acid (loss of 1 C as CO₂); internal gives 2 acids or keto-acid",
    },
]


@ChemMCPManager.register_tool
class CarbonChainBuilder(BaseTool):
    __version__      = "0.1.0"
    name             = "CarbonChainBuilder"
    func_name        = "build_carbon_chain"
    description      = "Suggest carbon chain elongation and shortening strategies with specific reaction sequences, reagents, conditions, and expected outcomes."
    implementation_description = "Database-driven approach covering Wittig, Grignard, aldol, cyanation, Suzuki, Heck, Diels-Alder for elongation; Hofmann, oxidative cleavage, decarboxylation, Hunsdiecker for shortening."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Carbon Chain", "Chain Elongation", "Chain Shortening", "Synthesis"]
    required_envs    = []

    code_input_sig   = [
        ("current_smiles", "str", "N/A", "SMILES string of the current molecule."),
        ("operation", "str", "elongate", "Operation: 'elongate' or 'shorten'."),
        ("target_carbon_count", "int", "0", "Target number of carbons (0 = suggest all applicable methods)."),
        ("method_preference", "str", "auto", "Method preference: 'auto', 'named_reaction', 'organometallic', 'pericyclic', 'degradation', 'oxidative'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'SMILES operation [target_C] [preference]'. Example: 'CCO elongate 5 auto'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: current_smiles, operation, current_carbon_count, suggested_methods (list of {name, carbon_change, starting_fg, reagents, conditions, yield, example}), summary."),
    ]

    examples         = [
        {
            "code_input": {
                "current_smiles": "CCO",
                "operation": "elongate",
                "target_carbon_count": 0,
                "method_preference": "auto",
            },
            "text_input": {"input_params": "CCO elongate 0 auto"},
            "output": {
                "result": {
                    "current_smiles": "CCO",
                    "operation": "elongate",
                    "current_carbon_count": 2,
                    "suggested_methods": [
                        {
                            "name": "Grignard Addition",
                            "carbon_change": "+1 to +n",
                            "reagents": "R'MgX (or R'Li), ether/THF",
                            "yield": "60-90%",
                            "example": "Convert alcohol to aldehyde first, then add Grignard reagent",
                        }
                    ],
                    "summary": "Ethanol (C2): oxidize to acetaldehyde, then apply Grignard/Wittig/cyanation for chain growth",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, current_smiles: str, operation: str = "elongate", target_carbon_count: int = 0, method_preference: str = "auto") -> dict:
        """Core logic: suggest chain modification methods."""
        if not current_smiles:
            raise ChemMCPError("Current SMILES is required.")

        op = operation.lower().strip()
        if op not in ("elongate", "shorten"):
            raise ChemMCPError("Operation must be 'elongate' or 'shorten'")

        # Count carbons in current molecule
        c_count = self._count_carbons(current_smiles)

        # Determine delta needed
        if target_carbon_count > 0:
            delta = target_carbon_count - c_count
            if (delta > 0 and op == "shorten") or (delta < 0 and op == "elongate"):
                return {
                    "current_smiles": current_smiles,
                    "operation": op,
                    "current_carbon_count": c_count,
                    "target_carbon_count": target_carbon_count,
                    "note": f"⚠ Target ({target_carbon_count} C) conflicts with operation ({op}). Showing {op} methods anyway.",
                    "suggested_methods": [],
                }

        # Select method database
        all_methods = ELONGATION_METHODS if op == "elongate" else SHORTENING_METHODS

        # Filter by preference
        if method_preference != "auto":
            pref = method_preference.lower()
            filtered = []
            for m in all_methods:
                name_lower = m["name"].lower()
                if pref == "named_reaction" and any(x in name_lower for x in ["wittig", "grignard", "aldol", "heck", "suzuki", "diels", "hofmann", "hunsdiecker"]):
                    filtered.append(m)
                elif pref == "organometallic" and any(x in name_lower for x in ["grignard", "suzuki", "heck", "ylide"]):
                    filtered.append(m)
                elif pref == "pericyclic" and "diels" in name_lower:
                    filtered.append(m)
                elif pref == "degradation" and any(x in name_lower for x in ["hofmann", "decarboxylat", "hunsdiecker", "kochi"]):
                    filtered.append(m)
                elif pref == "oxidative" and any(x in name_lower for x in ["oxidative", "periodate", "ozone", "kmno4"]):
                    filtered.append(m)
                else:
                    filtered.append(m)
            methods = filtered if filtered else all_methods
        else:
            methods = all_methods

        # Analyze current FG to prioritize relevant methods
        current_fgs = self._detect_fg(current_smiles)

        # Score and rank methods
        scored_methods = []
        for m in methods:
            relevance = 0
            for fg in current_fgs:
                if fg in m.get("starting_fg", []):
                    relevance += 2
                if fg in m.get("product_fg", []):
                    relevance += 1

            scored_methods.append({**m, "_relevance_score": relevance})

        scored_methods.sort(key=lambda x: x["_relevance_score"], reverse=True)

        # Format output
        suggested = []
        for m in scored_methods[:8]:  # top 8
            suggested.append({
                "name": m["name"],
                "carbon_change": m.get("carbon_change", "variable"),
                "starting_functional_group_required": m.get("starting_fg", []),
                "reagents": m.get("reagents", ""),
                "conditions": m.get("conditions", ""),
                "expected_yield": m.get("yield", "?"),
                "example_reaction": m.get("example", ""),
                "key_notes": m.get("notes", ""),
                "relevance_to_current_molecule": "High" if m["_relevance_score"] >= 2 else ("Medium" if m["_relevance_score"] > 0 else "General"),
            })

        # Build summary
        if op == "elongate":
            summary = f"Molecule has {c_count} carbons. "
            if target_carbon_count > c_count:
                summary += f"Need +{target_carbon_count - c_count} carbons. "
            summary += f"Suggested {len(suggested)} chain elongation method(s)."
        else:
            summary = f"Molecule has {c_count} carbons. "
            if target_carbon_count > 0 and target_carbon_count < c_count:
                summary += f"Need -{c_count - target_carbon_count} carbons. "
            summary += f"Suggested {len(suggested)} chain shortening method(s)."

        return {
            "current_smiles": current_smiles,
            "operation": op,
            "current_carbon_count": c_count,
            "target_carbon_count": target_carbon_count if target_carbon_count > 0 else None,
            "detected_functional_groups": current_fgs,
            "method_preference_used": method_preference,
            "suggested_methods": suggested,
            "summary": summary,
        }

    @staticmethod
    def _count_carbons(smiles: str) -> int:
        """Count carbon atoms in SMILES."""
        s = smiles
        count = 0
        # Uppercase C not followed by lowercase (which would be Cl, etc.)
        i = 0
        while i < len(s):
            if s[i] == 'C':
                if i + 1 < len(s) and s[i+1].islower():
                    # This is Cl, Ca, etc. — skip
                    if s[i+1] == 'l':  # Cl
                        i += 2
                        continue
                count += 1
            i += 1
        # Also count aromatic carbons (lowercase c)
        count += s.count('c')
        return count

    @staticmethod
    def _detect_fg(smiles: str) -> List[str]:
        """Detect functional groups in SMILES."""
        s = smiles
        fgs = []
        if re.search(r'C\(=O\)[OH,O]', s) or re.search(r'[OH]C\(=O\)', s):
            fgs.append("carboxylic_acid")
        if re.search(r'C\(=O\)OC|OC\(=O\)', s):
            fgs.append("ester")
        if re.search(r'C\(=O\)N|NC\(=O\)', s):
            fgs.append("amide")
        if re.search(r'(?:^|[^\(])C(=O)(?![a-z])', s) and not re.search(r'C\(=O\)[ON]', s):
            fgs.append("ketone")
        if re.search(r'C=O$', s) or re.search(r'\(=O\)[^a-z)]', s):
            fgs.append("aldehyde")
        if re.search(r'CO(?![a-z])', s) and not re.search(r'C\(=O\)', s):
            fgs.append("alcohol")
        if re.search(r'C=C', s):
            fgs.append("alkene")
        if re.search(r'C#C', s):
            fgs.append("alkyne")
        if re.search(r'NC|N\(', s) and not re.search(r'C\(=O\)N', s):
            fgs.append("amine")
        if re.search(r'[BClFI](?=[a-z]|[0-9()#=\[\]]|$)', s):
            fgs.append("alkyl_halide")
        if re.search(r'c[1-9]', s) or re.search(r'c1.*c1', s):
            fgs.append("aromatic_ring")
        if not fgs:
            fgs.append("alkane")
        return fgs

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            smiles_str = parts[0]
            op = parts[1] if len(parts) > 1 else "elongate"
            tgt_c = int(parts[2]) if len(parts) > 2 else 0
            pref = parts[3] if len(parts) > 3 else "auto"
            return self._run_base(smiles_str, op, tgt_c, pref)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'SMILES operation [target_C] [preference]'")
