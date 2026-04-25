import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RingFormationStrategy(BaseTool):
    """
    成环反应策略工具 - 推荐合适的成环方法。
    支持3-8元环的常见成环策略，包括分子内缩合、关环复分解、内酰胺化等。
    """
    __version__ = "0.1.0"
    name             = "RingFormationStrategy"
    func_name        = "ring_formation_strategy"
    description      = "Recommend ring-forming reaction strategies based on target ring size, starting materials, and constraints."
    implementation_description = "Uses a knowledge base of named ring-forming reactions (aldol, Dieckmann, lactamization, RCM, Nazarov, etc.) to recommend strategies for 3-8 membered rings."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Ring Formation", "Cyclization", "Synthetic Strategy", "Heterocycles"]
    required_envs    = []

    code_input_sig   = [
        ("target_ring_size", "int", "N/A", "Target ring size (3-8 membered rings supported)."),
        ("starting_material_hint", "str", "", "Optional hint about starting material type (e.g., 'diacid', 'dihalide', 'enone', 'amino alcohol')."),
        ("constraints", "str", "", "Optional constraints (e.g., 'stereocontrol', 'mild conditions', 'no metal')."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'ring_size [starting_material] [constraints]'. Example: '6 diacid stereocontrol'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed strategy recommendation including reaction name, mechanism type, typical conditions, pros/cons, and examples."),
    ]

    examples         = [
        {
            "code_input": {
                "target_ring_size": 6,
                "starting_material_hint": "diacid",
                "constraints": ""
            },
            "text_input": {
                "input_params": "6 diacid"
            },
            "output": {
                "result": "Recommended: Dieckmann Condensation..."
            },
        },
        {
            "code_input": {
                "target_ring_size": 5,
                "starting_material_hint": "",
                "constraints": "stereocontrol"
            },
            "text_input": {
                "input_params": "5 '' stereocontrol"
            },
            "output": {
                "result": "Recommended: Intramolecular Aldol..."
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_knowledge_base()

    def _build_knowledge_base(self):
        """Build comprehensive ring formation knowledge base."""
        self.strategies = {
            # === 3-membered rings ===
            3: [
                {
                    "name": "Intramolecular Nucleophilic Substitution (S_N2)",
                    "mechanism": "Halohydrin cyclization / Epoxide formation from halohydrins under base",
                    "substrate_requirement": "β-halo alcohol or similar 1,2-difunctionalized compound",
                    "conditions": "Base (NaOH, NaH), mild temperature (0-25°C)",
                    "pros": ["High yield", "Mild conditions", "Stereospecific"],
                    "cons": ["Limited to 3-membered rings", "Requires specific functional groups"],
                    "examples": "Epoxide from halohydrin; Aziridine from β-halo amine",
                    "stereocontrol": "High (inversion at carbon)",
                    "metal_free": True,
                    "key_reference": "Corey-Chaykovsky epoxidation (for carbonyl → epoxide)",
                    "category": "nucleophilic_substitution",
                },
                {
                    "name": "Simmons-Smith Cyclopropanation",
                    "mechanism": "Carbenoid addition to alkene (CH2I2 + Zn(Cu) → :CH2)",
                    "substrate_requirement": "Alkene with appropriate tether to form cyclopropane",
                    "conditions": "Zn(Cu) couple, CH2I2, Et2O, 0°C to rt",
                    "pros": ["Stereospecific (cis-addition)", "Tolerates many functional groups"],
                    "cons": ["Requires Zn reagent", "Diastereoselectivity depends on substrate"],
                    "examples": "Cyclopropane from terminal alkene; Fused cyclopropanes",
                    "stereocontrol": "High (stereospecific cis addition)",
                    "metal_free": False,
                    "key_reference": "Simmons-Smith, JACS 1958",
                    "category": "carbene_addition",
                },
            ],
            # === 4-membered rings ===
            4: [
                {
                    "name": "[2+2] Photocycloaddition",
                    "mechanism": "Photochemical [π2s + π2s] cycloaddition between two alkenes or enone + alkene",
                    "substrate_requirement": "Two π-systems in proximity (often intramolecular)",
                    "conditions": "hv (UV light, ~300 nm), often photosensitizer, low temperature",
                    "pros": ["Atom economical", "Unique access to strained rings", "Can form complex polycycles"],
                    "cons": ["Requires UV equipment", "Regio- and stereoselectivity can be tricky"],
                    "examples": "Oxetane from carbonyl + alkene (Paternò–Büchi); Cyclobutane from two alkenes",
                    "stereocontrol": "Moderate to high (suprafacial addition)",
                    "metal_free": True,
                    "key_reference": "de Mayo reaction; Paternò–Büchi",
                    "category": "pericyclic",
                },
                {
                    "name": "Intramolecular Aldol (4-exo-trig)",
                    "mechanism": "Enolate attack on ketone/aldehyde via 4-exo-trig cyclization",
                    "substrate_requirement": "1,4-Diketone or keto-aldehyde with enolizable position",
                    "conditions": "Base (LDA, KOH), -78°C to 0°C",
                    "pros": ["Readily available substrates", "Forms oxetanol/cyclobutanone"],
                    "cons": ["4-exo-trig is disfavored (Baldwin rules)", "Competing pathways common"],
                    "examples": "Cyclobutanone from 1,4-diketide",
                    "stereocontrol": "Moderate",
                    "metal_free": True,
                    "key_reference": "Baldwin rules for ring closure",
                    "category": "aldol_condensation",
                },
                {
                    "name": "Ketene Dimerization / [2+2]",
                    "mechanism": "Ketene dimerization or ketene + imine → β-lactam (Staudinger synthesis)",
                    "substrate_requirement": "Acid chloride + Et3N (generates ketene) + imine",
                    "conditions": "Acid chloride, Et3N, imine, CH2Cl2, 0°C → rt",
                    "pros": ["Direct β-lactam synthesis", "Well-established for antibiotics"],
                    "cons": ["Ketenes are reactive/hazardous", "Limited scope"],
                    "examples": "β-Lactam core of penicillins/cephalosporins",
                    "stereocontrol": "Moderate (cis-selective typically)",
                    "metal_free": True,
                    "key_reference": "Staudinger β-lactam synthesis",
                    "category": "cycloaddition",
                },
            ],
            # === 5-membered rings ===
            5: [
                {
                    "name": "Intramolecular Aldol Condensation (5-endo-trig)",
                    "mechanism": "Enolate attacks carbonyl within same molecule to form 5-membered ring",
                    "substrate_requirement": "1,5-Dicarbonyl compound or equivalent",
                    "conditions": "Base (NaOH, KOH, LDA), often heat (reflux)",
                    "pros": ["Very reliable", "Substrates readily available", "High yields common"],
                    "cons": ["May give mixture of stereoisomers", "Dehydration competes"],
                    "examples": "Cyclopentenone synthesis; Hydrindanone systems; Jasmonoid natural products",
                    "stereocontrol": "Moderate (can be controlled with chiral auxiliaries/catalysis)",
                    "metal_free": True,
                    "key_reference": "Robinson annulation (related); Nazarov cyclization variant",
                    "category": "aldol_condensation",
                },
                {
                    "name": "Intramolecular Michael Addition (5-exo-trig)",
                    "mechanism": "Nucleophile adds to α,β-unsaturated carbonyl in 5-exo-trig fashion",
                    "substrate_requirement": "Nucleophile tethered to Michael acceptor (1,4-relationship)",
                    "conditions": "Base (DBU, K2CO3, amine), rt to reflux",
                    "pros": ["Favored by Baldwin rules (exo-trig)", "Mild conditions", "Broad substrate scope"],
                    "cons": ["1,4 vs 1,2-addition competition possible", "Polymerization side reactions"],
                    "examples": "Cyclopentane fused systems; Bicyclic natural product cores",
                    "stereocontrol": "Moderate to high with chiral catalysts",
                    "metal_free": True,
                    "key_reference": "Michael-initiated ring closures",
                    "category": "michael_addition",
                },
                {
                    "name": "Intramolecular Williamson Ether Synthesis",
                    "mechanism": "S_N2 displacement of halide by alkoxide to form cyclic ether",
                    "substrate_requirement": " Halo alcohol with halide and OH separated by 3 carbons (for 5-ring)",
                    "conditions": "Base (NaOH, NaH, K2CO3), polar aprotic solvent (DMF, THF), heat",
                    "pros": ["Simple, reliable", "High yields", "Tetrahydrofuran (THF) ring is stable"],
                    "cons": ["Requires good leaving group", "S_N2 limits (no tertiary centers)", "Competition with elimination"],
                    "examples": "THF derivatives; Sugar-derived tetrahydrofurans; Oxygen heterocycles",
                    "stereocontrol": "High (inversion at electrophilic carbon)",
                    "metal_free": True,
                    "key_reference": "Williamson ether synthesis (classical)",
                    "category": "nucleophilic_substitution",
                },
                {
                    "name": "Paal-Knorr Pyrrole/Furan/Thiophene Synthesis",
                    "mechanism": "1,4-Diketone condensation with primary amine (pyrrole) or acid (furan) or P4S10/Lawesson's (thiophene)",
                    "substrate_requirement": "1,4-Diketone + amine (for pyrrole) / acid catalyst (for furan)",
                    "conditions": "Amine (for pyrrole): AcOH, reflux; Acid (for furan): H2SO4 or p-TsOH, heat",
                    "pros": ["One-step heterocycle formation", "Versatile (N/O/S heterocycles)", "High yields"],
                    "cons": ["Requires 1,4-diketone precursor", "Limited substitution pattern control"],
                    "examples": "Pyrrole alkaloids; Furan natural products; Thiophene pharmaceuticals",
                    "stereocontrol": "N/A (aromatic products)",
                    "metal_free": True,
                    "key_reference": "Paal-Knorr synthesis (1884, classical)",
                    "category": "condensation",
                },
                {
                    "name": "1,3-Dipolar Cycloaddition (Intramolecular)",
                    "mechanism": "[3+2] cycloaddition between dipole (nitrone, azomethine ylide, diazo) and dipolarophile",
                    "substrate_requirement": "Dipole precursor + tethered alkene/alkyne as dipolarophile",
                    "conditions": "Thermal or Lewis acid catalyzed, solvent varies",
                    "pros": ["High atom economy", "Excellent stereocontrol", "Multiple bonds formed in one step"],
                    "cons": ["Requires specific dipole/dipolarophile pair", "Regioselectivity can be an issue"],
                    "examples": "Pyrrolidines (natural product cores); Isoxazolines; Proline derivatives",
                    "stereocontrol": "High (concerted, suprafacial)",
                    "metal_free": True,
                    "key_reference": "Huisgen cycloaddition; Intramolecular nitrone cycloaddition",
                    "category": "cycloaddition",
                },
            ],
            # === 6-membered rings ===
            6: [
                {
                    "name": "Dieckmann Condensation",
                    "mechanism": "Intramolecular Claisen condensation of diester → β-keto ester (cyclic)",
                    "substrate_requirement": "Diester with ester groups at 1,6-positions",
                    "conditions": "Strong base (NaOEt, NaH), ethanol or THF, 0°C → rt, then acid workup",
                    "pros": ["Classic, well-understood", "High yields", "Product is versatile β-keto ester"],
                    "cons": ["Only works for diesters", "Requires strong base", "Acidic α-H needed"],
                    "examples": "Cyclohexenone after decarboxylation; Wieland-Miescher ketone precursor",
                    "stereocontrol": "N/A (planar enolate intermediate)",
                    "metal_free": True,
                    "key_reference": "Dieckmann, Ber. 1894; Classical method",
                    "category": "claisen_condensation",
                },
                {
                    "name": "Robinson Annulation",
                    "mechanism": "Michael addition followed by intramolecular aldol condensation → cyclohexenone",
                    "substrate_requirement": "Enolizable ketone + methyl vinyl ketone (MVK) or equivalent",
                    "conditions": "Base (KOH, NaOH, pyrrolidine), ethanol/water, heat (reflux)",
                    "pros": ["Powerful C-C bond forming sequence", "Forms bicyclic systems easily", "Widely used in total synthesis"],
                    "cons": ["Overalkylation possible", "Stereochemistry can be hard to control", "Basic conditions may not tolerate all FGs"],
                    "examples": "Steroid skeleton construction; Wieland-Miescher ketone; Carvone synthesis",
                    "stereocontrol": "Moderate (can be controlled with chiral auxiliaries)",
                    "metal_free": True,
                    "key_reference": "Raphael, Robinson; JCS 1935",
                    "category": "annulation",
                },
                {
                    "name": "Diels-Alder Cycloaddition (Intramolecular)",
                    "mechanism": "[4+2] cycloaddition between conjugated diene and dienophile (tethered)",
                    "substrate_requirement": "Diene and dienophile in same molecule (through tether)",
                    "conditions": "Thermal (heat, often >100°C) or Lewis acid catalyzed (EtAlCl2, etc.)",
                    "pros": ["Extremely powerful", "Up to 4 new stereocenters in one step", "Atom economical", "Predictable endo/exo selectivity"],
                    "cons": ["May require high temperature", "Cis-trans diene geometry matters", "Electron-poor dienophile usually needed"],
                    "examples": "Polycyclic natural products; Decalin systems; Endiandric acids (cascading iDA)",
                    "stereocontrol": "Very high (concerted, endo rule)",
                    "metal_free": True,
                    "key_reference": "Diels & Alder, 1928; Nobel Prize 1950",
                    "category": "cycloaddition",
                },
                {
                    "name": "Intramolecular Heck Reaction",
                    "mechanism": "Pd-catalyzed aryl/vinyl halide addition to tethered alkene → cyclization",
                    "substrate_requirement": "Aryl/vinyl halide (or triflate) + tethered alkene",
                    "conditions": "Pd(OAc)2 or Pd(PPh3)4, base (Et3N, K2CO3), polar solvent (DMF, MeCN), 80-120°C",
                    "pros": ["Forms C-C bond to sp2 carbon", "Tolerates many functional groups", "Good for medium/large rings too"],
                    "cons": ["Requires palladium catalyst", "β-Hydride elimination competing", "Expensive"],
                    "examples": "Benzofused rings; Indole synthesis; Lactam-containing macrocycles",
                    "stereocontrol": "Moderate (syn addition, then elimination)",
                    "metal_free": False,
                    "key_reference": "Heck, 1970s; Nobel Prize 2010",
                    "category": "cross_coupling",
                },
                {
                    "name": "Intramolecular Mukaiyama Aldol / Prins Cyclization",
                    "mechanism": "Lewis acid-promoted cyclization of homoallylic alcohol with aldehyde → tetrahydropyran",
                    "substrate_requirement": "Homoallylic alcohol + aldehyde (or equivalent)",
                    "conditions": "Lewis acid (SnCl4, BF3·OEt2, TMSOTf), CH2Cl2, -78°C to rt",
                    "pros": ["Excellent for oxygen heterocycles", "High stereocontrol (Prins cyclization)", "Mild conditions"],
                    "cons": ["Moisture sensitive", "Lewis acid required", "Side reactions possible"],
                    "examples": "Tetrahydropyran sugars; Polyether ionophore fragments; Prins cyclization products",
                    "stereocontrol": "High (chelation-controlled or Felkin-Anh model)",
                    "metal_free": False,
                    "key_reference": "Mukaiyama aldol; Prins cyclization",
                    "category": "lewis_acid_catalysis",
                },
                {
                    "name": "Biginelli Dihydropyrimidine Synthesis",
                    "mechanism": "Three-component condensation of aldehyde, β-keto ester, and urea/thiourea",
                    "substrate_requirement": "Aldehyde + ethyl acetoacetate + urea (or thiourea)",
                    "conditions": "Acid catalyst (HCl, p-TsOH) or Lewis acid, ethanol, reflux",
                    "pros": ["Multicomponent (high efficiency)", "Diverse pharmacophore (DHPMs have bioactivity)", "Simple setup"],
                    "cons": ["Limited to dihydropyrimidinones", "Yield can vary widely", "Acid-sensitive groups problematic"],
                    "examples": "Calcium channel modulators (monastrol); Antiviral agents; Antihypertensives",
                    "stereocontrol": "Low (achiral unless using chiral catalyst)",
                    "metal_free": True,
                    "key_reference": "Biginelli, 1893; Modern variants with Lewis acids",
                    "category": "multicomponent",
                },
            ],
            # === 7-membered rings ===
            7: [
                {
                    "name": "Ring-Closing Metathesis (RCM, 7-membered)",
                    "mechanism": "Ru-carbene catalyzed [2+2+2] cycloreversion/reformation → cyclic alkene",
                    "substrate_requirement": "Diene with terminal alkenes at 1,7-positions",
                    "conditions": "Grubbs 2nd gen catalyst (G-II) or Hoveyda-Grubbs, CH2Cl2 or toluene, reflux (40-80°C)",
                    "pros": ["Most general method for 7-rings", "Highly functional group tolerant", "Mild conditions", "E-selectivity controllable"],
                    "cons": ["Ruthenium catalyst expensive", "E/Z selectivity can be mixed", "Removal of Ru residues needed for pharma"],
                    "examples": "Azepanes; Benzazepines; Cycloheptene derivatives; Macrocyclic precursors",
                    "stereocontrol": "Moderate (E/Z mixture, can be controlled with catalyst choice)",
                    "metal_free": False,
                    "key_reference": "Grubbs & Schrock, 1990s; Nobel Prize 2005",
                    "category": "metathesis",
                },
                {
                    "name": "Intramolecular Aldol (7-endo-trig, disfavored)",
                    "mechanism": "Enolate attack on distant carbonyl (thermodynamically controlled)",
                    "substrate_requirement": "1,7-Dicarbonyl compound",
                    "conditions": "Strong base, high dilution, thermodynamic control",
                    "pros": ["Simple conceptually", "No transition metals"],
                    "cons": ["7-endo-trig disfavored (Baldwin)", "Low yield", "Many side products"],
                    "examples": "Cycloheptanone derivatives (low yield routes)",
                    "stereocontrol": "Low",
                    "metal_free": True,
                    "key_reference": "Baldwin rules favor exo over endo",
                    "category": "aldol_condensation",
                },
                {
                    "name": "Azepine/Benzazepine Formation via Beckmann/Beckmann Fragmentation",
                    "mechanism": "Beckmann rearrangement of cyclohexanone oxime → caprolactam (7-membered lactam)",
                    "substrate_requirement": "Cyclohexanone oxime",
                    "conditions": "Acid (H2SO4, PPA, PCl5, SOCl2, TsCl), heat",
                    "pros": ["Industrial process (nylon-6 precursor)", "Reliable", "Scalable"],
                    "cons": ["Limited to caprolactam-type products", "Harsh acidic conditions"],
                    "examples": "ε-Caprolactam (nylon-6 monomer); Benzazepine pharmaceutical intermediates",
                    "stereocontrol": "N/A (lactam product)",
                    "metal_free": True,
                    "key_reference": "Beckmann, 1886; Industrial nylon process",
                    "category": "rearrangement",
                },
            ],
            # === 8-membered rings ===
            8: [
                {
                    "name": "Ring-Closing Metathesis (RCM, 8-membered)",
                    "mechanism": "Ru-catalyzed RCM for medium-sized rings (8-membered)",
                    "substrate_requirement": "Diene with terminal alkenes at 1,8-positions",
                    "conditions": "Grubbs 2nd gen or Z-selective catalyst, high dilution (0.001-0.01 M), CH2Cl2, reflux",
                    "pros": ["Best method for 8-rings", "Functional group tolerance excellent", "High dilution minimizes oligomerization"],
                    "cons": ["Medium rings are challenging (transannular strain)", "High dilution needed (wasteful)", "Expensive catalyst"],
                    "examples": "Cyclooctene; Oxocanes; Medium-ring natural product fragments (e.g., taxol fragments)",
                    "stereocontrol": "Moderate (E/Z mixture)",
                    "metal_free": False,
                    "key_reference": "Grubbs RCM for medium rings; Z-catalysts for E-control",
                    "category": "metathesis",
                },
                {
                    "name": "Transannular Reactions (from larger rings)",
                    "mechanism": "Form larger ring first (RCM/macrolactonization), then transannular cyclization",
                    "substrate_requirement": "Larger ring (12-16) with appropriately positioned reactive groups",
                    "conditions": "Depends on transannular step (could be aldol, Michael, SN2, etc.)",
                    "pros": ["Avoids direct 8-membered ring formation difficulties", "Can build complex polycyclic systems"],
                    "cons": ["Multi-step approach", "Lower overall yield", "Complex planning needed"],
                    "examples": "Taxol core; Complex alkaloids; Polycyclic ether natural products",
                    "stereocontrol": "Varies by transannular step",
                    "metal_free": "Varies",
                    "key_reference": "Nicolaou's transannular strategies",
                    "category": "transannular",
                },
                {
                    "name": "Lactonization / Macrolactonization (8-membered lactone)",
                    "mechanism": "Intramolecular esterification of ω-hydroxy acid → medium-ring lactone",
                    "substrate_requirement": "ω-Hydroxy acid with OH and COOH at 1,8-positions",
                    "conditions": "High dilution, coupling agent (DCC, EDCI, Yamaguchi, Keck macrolactonization), mild",
                    "pros": ["Well-developed for macrolides", "Many coupling methods available", "Stereochemistry retained"],
                    "cons": ["High dilution required", "Competing dimer/oligomer formation", "Coupling agents expensive"],
                    "examples": "8-Membered ring lactone natural products; Medium-ring segments of macrolides",
                    "stereocontrol": "High (retention of configuration)",
                    "metal_free": True,
                    "key_reference": "Yamaguchi macrolactonization; Corey-Nicolaou; Keck",
                    "category": "lactonization",
                },
            ],
        }

        # Additional special strategies that cross-cut ring sizes
        self.special_strategies = {
            "nitrogen_heterocycle": {
                "name": "Nitrogen Heterocycle Formation Strategies",
                "methods": [
                    {"method": "Lactamization (amide bond formation)", "rings": [4, 5, 6, 7], "notes": "Amino acid → lactam; peptide coupling reagents"},
                    {"method": "Intramolecular Mitsunobu", "rings": [5, 6], "notes": "Alcohol + amine/nucleophile → N-heterocycle"},
                    {"method": "Bischler-Napieralski / Pictet-Spengler", "rings": [6], "notes": "Tetrahydroisoquinoline formation"},
                    {"method": "Intramolecular Amination (Buchwald-Hartwig)", "rings": [5, 6, 7, 8], "notes": "Pd-catalyzed C-N bond formation for ring closure"},
                    {"method": "Schmidt Reaction / Beckmann", "rings": [5, 6, 7], "notes": "Ketone + HN3 → lactam (Schmidt); Oxime → lactam (Beckmann)"},
                ]
            },
            "oxygen_heterocycle": {
                "name": "Oxygen Heterocycle Formation Strategies",
                "methods": [
                    {"method": "Williamson Ether Synthesis", "rings": [3, 4, 5, 6], "notes": "Classic S_N2 cyclization"},
                    {"method": "Prins Cyclization", "rings": [6], "notes": "Homoallylic alcohol + aldehyde → tetrahydropyran"},
                    {"method": "Epoxide Ring Opening (intramolecular)", "rings": [5, 6], "notes": "Epoxide opened by internal nucleophile"},
                    {"method": "Tandem Epoxidation-Cyclization", "rings": [5, 6], "notes": "Form epoxide then intramolecular open"},
                ]
            }
        }

    def _run_base(self, target_ring_size: int, starting_material_hint: str = "", constraints: str = "") -> str:
        """Generate ring formation strategy recommendations."""
        if target_ring_size < 3:
            raise ChemMCPError("Ring size must be >= 3.")
        if target_ring_size > 8:
            raise ChemMCPError("Ring sizes > 8 are considered macrocycles. Use macrocyclization strategies (RCM, macrolactonization, etc.).")

        strategies = self.strategies.get(target_ring_size, [])
        if not strategies:
            raise ChemMCPError(f"No strategies available for {target_ring_size}-membered rings.")

        # Filter by hints/constraints if provided
        hint_lower = (starting_material_hint.lower() + " " + constraints.lower()).strip() if starting_material_hint or constraints else ""

        result_parts = [f"## Ring Formation Strategies for {target_ring_size}-Membered Rings\n"]
        result_parts.append(f"Found **{len(strategies)}** recommended strategies:\n")

        for i, s in enumerate(strategies, 1):
            # Apply hint-based filtering/ranking
            relevance_score = 1.0
            if hint_lower:
                # Check substrate match
                sub_req = s.get("substrate_requirement", "").lower()
                for kw in hint_lower.split():
                    if kw in sub_req:
                        relevance_score += 0.5

                # Check constraint matches
                if "stereo" in hint_lower and s.get("stereocontrol") == "High":
                    relevance_score += 1.0
                if "no metal" in hint_lower or "metal free" in hint_lower:
                    if s.get("metal_free"):
                        relevance_score += 1.0
                    else:
                        relevance_score -= 0.5
                if "mild" in hint_lower:
                    if any(x in s.get("conditions", "").lower() for x in ["rt", "0°", "25°", "mild"]):
                        relevance_score += 0.5

            s["_relevance"] = relevance_score

        # Sort by relevance
        strategies.sort(key=lambda x: x.get("_relevance", 0), reverse=True)

        for i, s in enumerate(strategies, 1):
            rel_tag = " ⭐" if s.get("_relevance", 0) > 2 else ""
            result_parts.append(f"### Strategy {i}: {s['name']}{rel_tag}")
            result_parts.append(f"- **Mechanism:** {s['mechanism']}")
            result_parts.append(f"- **Substrate:** {s['substrate_requirement']}")
            result_parts.append(f"- **Conditions:** {s['conditions']}")
            result_parts.append(f"- **Pros:** {', '.join(s['pros'])}")
            result_parts.append(f"- **Cons:** {', '.join(s['cons'])}")
            result_parts.append(f"- **Examples:** {s['examples']}")
            result_parts.append(f"- **Stereocontrol:** {s['stereocontrol']}")
            result_parts.append(f"- **Metal-free:** {'Yes' if s['metal_free'] else 'No'}")
            result_parts.append(f"- **Category:** {s['category']}")
            result_parts.append("")

        # Add special strategies if relevant
        if hint_lower:
            for key, spec in self.special_strategies.items():
                if any(kw in hint_lower for kw in key.split("_")):
                    result_parts.append(f"\n### Additional: {spec['name']}\n")
                    for m in spec["methods"]:
                        rings_info = ", ".join(str(r) for r in m["rings"])
                        result_parts.append(f"- **{m['method']}** (rings: {rings_info}): {m['notes']}")
                    result_parts.append("")

        # Add Baldwin rules summary
        result_parts.append("\n---\n### 📐 Baldwin's Rules Reference\n")
        baldwin_data = {
            3: ["tet (favored)", "trig (disfavored)"],
            4: ["tet (favored)", "trig (disfavored)", "dig (favored)"],
            5: ["exo-tet (favored)", "exo-trig (favored)", "exo-dig (favored)", "endo-trig (disfavored)"],
            6: ["exo-tet (favored)", "exo-trig (favored)", "exo-dig (favored)", "endo-trig (allowed but slow)"],
            7: ["exo-tet (favored)", "exo-trig (favored)", "endo-trig (disfavored)"],
            8: ["exo-tet (favored)", "exo-trig (favored)", "transannular strain significant"],
        }
        if target_ring_size in baldwin_data:
            result_parts.append(f"For {target_ring_size}-membered rings: {', '.join(baldwin_data[target_ring_size])}")

        return "\n".join(result_parts)

    def _run_text(self, input_params: str) -> str:
        """Parse text input."""
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Input must include at least ring_size. Format: 'ring_size [starting_material] [constraints]'")
        try:
            ring_size = int(parts[0])
        except ValueError:
            raise ChemMCPError(f"First parameter must be integer ring size, got '{parts[0]}'")
        start_hint = parts[1] if len(parts) > 1 else ""
        constraints = " ".join(parts[2:]) if len(parts) > 2 else ""
        return self._run_base(ring_size, start_hint, constraints)
