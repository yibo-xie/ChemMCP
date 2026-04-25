"""
Functional Group Interconversion - suggests pathways for interconverting
functional groups with reagents, conditions, mechanisms, and expected yields.
"""

import logging
from typing import Dict, List, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Comprehensive functional group interconversion database
# Each entry: (source_fg, target_fg) -> [pathway options]
FG_INTERCONVERSION_DB: Dict[str, Dict[str, List[Dict[str, Any]]]] = {
    # === ALCOHOL conversions ===
    "alcohol → aldehyde": [
        {
            "method": "Oxidation (mild)",
            "reagents": "PCC / PDC in CH₂Cl₂",
            "conditions": "RT, 1-2 h, anhydrous",
            "mechanism": "Chromate ester formation + β-elimination",
            "yield": "70-85%",
            "notes": "Stops at aldehyde; over-oxidation to acid is minimal",
            "selectivity": "Primary alcohol only; secondary → ketone",
        },
        {
            "method": "Swern oxidation",
            "reagents": "(COCl)₂, DMSO, Et₃N",
            "conditions": "-78°C → RT, CH₂Cl₂",
            "mechanism": "Activation of DMSO with chloroalkane + base-mediated elimination",
            "yield": "75-90%",
            "notes": "Mild, acid-free; good for sensitive substrates",
            "selectivity": "Primary → aldehyde; secondary → ketone",
        },
        {
            "method": "Dess-Martin periodinane",
            "reagents": "DMP in CH₂Cl₂",
            "conditions": "RT, 30 min - 2 h",
            "mechanism": "Hypervalent iodine-mediated oxidation via alkoxyperiodinane",
            "yield": "85-95%",
            "notes": "Very mild, high-yielding, operationally simple",
            "selectivity": "Excellent for primary/secondary alcohols",
        },
    ],
    "alcohol → ketone": [
        {
            "method": "Oxidation (strong)",
            "reagents": "Jones reagent (CrO₃/H₂SO₄/acetone)",
            "conditions": "0°C → RT, aqueous acetone",
            "mechanism": "Chromate ester formation + elimination",
            "yield": "70-95%",
            "notes": "Classic method; acidic conditions may affect acid-sensitive groups",
            "selectivity": "Secondary alcohol → ketone; primary → carboxylic acid",
        },
    ],
    "alcohol → carboxylic_acid": [
        {
            "method": "Strong oxidation",
            "reagents": "KMnO₄ (aq.) or Na₂Cr₂O₇/H₂SO₄ or HNO₃",
            "conditions": "Heat or reflux depending on oxidant",
            "mechanism": "Multi-step oxidation through aldehyde intermediate",
            "yield": "60-90%",
            "notes": "Jones oxidation also gives acid from 1° alcohols",
            "selectivity": "Primary alcohols only",
        },
    ],
    "alcohol → alkene": [
        {
            "method": "Acid-catalyzed dehydration",
            "reagents": "H₂SO₄ or H₃PO₄ or TsOH / heat",
            "conditions": "80-180°C (depends on substrate)",
            "mechanism": "E1 or E2 depending on conditions",
            "yield": "60-80%",
            "notes": "Follows Zaitsev rule (more substituted alkene favored); possible rearrangement",
            "selectivity": "May give mixture of isomers",
        },
        {
            "method": "POCl₃ / pyridine (for sensitive substrates)",
            "reagents": "POCl₃, pyridine, 0°C → RT",
            "conditions": "0°C → RT, 1-4 h",
            "mechanism": "E2 elimination via phosphoryl ester intermediate",
            "yield": "70-90%",
            "notes": "Mild, non-acidic; good for acid-sensitive substrates",
            "selectivity": "Often gives less substituted (Hofmann) product with POCl₃/py",
        },
    ],
    "alcohol → alkyl_halide": [
        {
            "method": "Nucleophilic substitution (HX)",
            "reagents": "HX (HCl, HBr, HI) or PBr₃ / SOCl₂",
            "conditions": "RT or gentle heating",
            "mechanism": "SN1 (3°) or SN2 (1°) with inversion",
            "yield": "65-95%",
            "notes": "PBr₃ for bromides; SOCl₂ for chlorides (with pyridine)",
            "selectivity": "May involve rearrangement for tertiary alcohols (SN1)",
        },
    ],

    # === ALDEHYDE conversions ===
    "aldehyde → alcohol": [
        {
            "method": "Reduction",
            "reagents": "NaBH₄ in MeOH/EtOH or LiAlH₄ in ether",
            "conditions": "0°C → RT, 30 min - 2 h",
            "mechanism": "Hydride transfer to carbonyl carbon",
            "yield": "80-99%",
            "notes": "NaBH₄ is milder and selective; LiAlH₄ is more reactive",
            "selectivity": "Chemoselective: reduces aldehydes/ketones but not esters/amides (NaBH₄)",
        },
    ],
    "aldehyde → carboxylic_acid": [
        {
            "method": "Oxidation",
            "reagents": "Ag₂O (Tollens') / KMnO₄ / CrO₃ / NaClO₂ (Pinnick)",
            "conditions": "RT to reflux depending on method",
            "mechanism": "Hydride transfer or nucleophilic addition of oxidant",
            "yield": "70-95%",
            "notes": "Pinnick oxidation (NaClO₂) is highly chemoselective",
            "selectivity": "Pinnick: tolerates alkenes, alkynes; Ag₂O is specific for α-hydroxy aldehydes",
        },
    ],
    "aldehyde → alkene": [
        {
            "method": "Wittig reaction",
            "reagents": "Ph₃P=CHR (phosphonium ylide), THF",
            "conditions": "0°C → RT, under N₂",
            "mechanism": "[2+2] cycloaddition → oxaphosphetane → elimination",
            "yield": "50-90%",
            "notes": "Non-stabilized ylides give Z-alkene; stabilized ylides give E-alkene",
            "selectivity": "E/Z depends on ylide type and conditions",
        },
    ],

    # === KETONE conversions ===
    "ketone → alcohol": [
        {
            "method": "Reduction",
            "reagents": "NaBH₄ / LiAlH₄ / selective reducing agents",
            "conditions": "RT, protic solvent (NaBH₄) or aprotic (LiAlH₄)",
            "mechanism": "Hydride delivery to carbonyl C",
            "yield": "80-99%",
            "notes": "Can be stereoselective (chiral reduction with CBS catalyst)",
            "selectivity": "Gives secondary alcohol",
        },
    ],
    "ketone → alkane": [
        {
            "method": "Clemmensen reduction",
            "reagents": "Zn(Hg), HCl (conc.), reflux",
            "conditions": "Reflux, strongly acidic",
            "mechanism": "Carbenoid mechanism on zinc surface",
            "yield": "50-85%",
            "notes": "Acid-sensitive substrates not suitable",
            "selectivity": "Reduces C=O to CH₂",
        },
        {
            "method": "Wolff-Kishner reduction",
            "reagents": "NH₂NH₂, KOH, ethylene glycol, 180-200°C",
            "conditions": "High temperature, basic",
            "mechanism": "Formation of hydrazone → loss of N₂",
            "yield": "40-80%",
            "notes": "Base-sensitive substrates not suitable; complementary to Clemmensen",
            "selectivity": "Reduces C=O to CH₂",
        },
    ],
    "ketone → alkene": [
        {
            "method": "Wittig reaction",
            "reagents": "Ph₃P=CHR (stabilized ylide preferred)",
            "conditions": "THF or DCM, RT-reflux",
            "mechanism": "Standard Wittig olefination",
            "yield": "60-90%",
            "notes": "Tetrasubstituted alkenes from ketones",
            "selectivity": "Usually E-selective with stabilized ylides",
        },
    ],

    # === CARBOXYYLIC ACID conversions ===
    "carboxylic_acid → alcohol": [
        {
            "method": "Reduction (via activated derivative)",
            "reagents": "LiAlH₄ (direct) or 1) SOCl₂ 2) LiAlH₄ or B₂H₆",
            "conditions": "Reflux in ether/THF (LiAlH₄)",
            "mechanism": "Nucleophilic acyl substitution then hydride addition",
            "yield": "70-95%",
            "notes": "LiAlH₄ directly reduces acids; BH₃·THF is more selective",
            "selectivity": "Reduces COOH to CH₂OH",
        },
    ],
    "carboxylic_acid → ester": [
        {
            "method": "Fischer esterification",
            "reagents": "R'OH, H⁺ (cat. H₂SO₄ or TsOH)",
            "conditions": "Reflux with water removal (Dean-Stark)",
            "mechanism": "Protonation, nucleophilic attack, dehydration",
            "yield": "60-85%",
            "notes": "Equilibrium reaction; use excess alcohol or remove water",
            "selectivity": "Works for most acids; phenols need special conditions",
        },
        {
            "method": "Steglich esterification (coupling)",
            "reagents": "R'OH, DCC, DMAP, CH₂Cl₂",
            "conditions": "0°C → RT, 12-24 h",
            "mechanism": "DCC activates acid to O-acylisourea; DMAP catalyzes acyl transfer",
            "yield": "80-95%",
            "notes": "Mild, high-yielding; good for sensitive/macrocycle substrates",
            "selectivity": "Excellent functional group tolerance",
        },
    ],
    "carboxylic_acid → amide": [
        {
            "method": "Coupling activation",
            "reagents": "R'NH₂, EDCI/HOBt or DCC or (COCl)₂ then amine",
            "conditions": "RT, DMF or CH₂Cl₂",
            "mechanism": "Acyl activation followed by aminolysis",
            "yield": "70-95%",
            "notes": "EDCI/HOBt minimizes racemization for amino acids",
            "selectivity": "Widely used in peptide synthesis",
        },
    ],

    # === ALKENE conversions ===
    "alkene → alkane": [
        {
            "method": "Catalytic hydrogenation",
            "reagents": "H₂, Pd/C or PtO₂ or Raney Ni",
            "conditions": "RT, 1-4 atm H₂", "solvent": "EtOH, EtOAc, or hexanes",
            "mechanism": "Syn-addition of H₂ on metal surface",
            "yield": "90-99%",
            "notes": "Stereospecific syn addition; cis-alkene → meso/achiral product",
            "selectivity": "Selective for C=C over C=O (with Pd/C); Lindlar's for alkynes → cis-alkenes",
        },
    ],
    "alkene → alcohol": [
        {
            "method": "Oxymercuration-demercuration",
            "reagents": "Hg(OAc)₂, H₂O/THF, then NaBH₄",
            "conditions": "RT, 1-2 h total",
            "mechanism": "Markovnikov addition of HgOAc, then anti-reduction",
            "yield": "80-95%",
            "notes": "Markovnikov regioselectivity; no rearrangement",
            "selectivity": "Gives Markovnikov alcohol (anti-Markovnikov use hydroboration-oxidation)",
        },
        {
            "method": "Hydroboration-oxidation",
            "reagents": "1) BH₃·THF 2) H₂O₂, NaOH",
            "conditions": "0°C → RT, then basic workup",
            "mechanism": "Syn addition of B-H (anti-Markovnikov), then oxidation",
            "yield": "75-95%",
            "notes": "Anti-Markovnikov regioselectivity; syn stereochemistry",
            "selectivity": "Gives less substituted (primary) alcohol from terminal alkenes",
        },
    ],
    "alkene → diol": [
        {
            "method": "OsO₄-catalyzed dihydroxylation",
            "reagents": "OsO₄ (cat.), NMO or K₃Fe(CN)₆, acetone/H₂O",
            "conditions": "RT, 4-12 h",
            "mechanism": "[3+2] cycloaddition of OsO₄ across double bond",
            "yield": "70-95%",
            "notes": "Syn stereoselectivity; Sharpless AD gives enantioselective version",
            "selectivity": "Vicinal syn-diol formation",
        },
        {
            "method": "Cold KMnO₄ (dilute)",
            "reagents": "KMnO₄ (cold, dilute), NaOH (aq.)",
            "conditions": "0°C, pH > 10",
            "mechanism": "Cyclic manganate ester hydrolysis",
            "yield": "50-80%",
            "notes": "Older method; OsO₄ is preferred for reliability",
            "selectivity": "syn-Diol; over-oxidation if too harsh",
        },
    ],
    "alkene → epoxide": [
        {
            "method": "mCPBA epoxidation",
            "reagents": "mCPBA (meta-chloroperoxybenzoic acid), CH₂Cl₂",
            "conditions": "0°C → RT, 2-8 h, dark",
            "mechanism": "Concerted one-step oxygen transfer (stereospecific)",
            "yield": "70-90%",
            "notes": "Stereospecific retention of alkene geometry",
            "selectivity": "Electron-rich alkenes react faster",
        },
    ],

    # === ALKYNE conversions ===
    "alkyne → cis_alkene": [
        {
            "method": "Lindlar hydrogenation",
            "reagents": "H₂, Lindlar's catalyst (Pd/CaCO₃/Pb/quinoline)",
            "conditions": "RT, 1 atm H₂, hexane or EtOAc",
            "mechanism": "Syn addition of H₂ on poisoned Pd surface",
            "yield": "80-95%",
            "notes": "Poisoned catalyst stops at cis-alkene; quinoline prevents over-reduction",
            "selectivity": "cis-Alkene specifically; trans requires Na/NH₃",
        },
    ],
    "alkyne → trans_alkene": [
        {
            "method": "Dissolving metal reduction",
            "reagents": "Na or Li, NH₃(l), -78°C",
            "conditions": "-78°C, liquid ammonia, THF co-solvent",
            "mechanism": "Single-electron transfer → radical anion → trans-alkene",
            "yield": "70-90%",
            "notes": "Gives trans-alkene stereoselectively",
            "selectivity": "trans-Alkene specifically",
        },
    ],
    "alkyne → ketone": [
        {
            "method": "Hydration (Markovnikov)",
            "reagents": "HgSO₄, H₂SO₄ (aq.), heat",
            "conditions": "Reflux, aqueous acidic",
            "mechanism": "Hg²⁺-catalyzed keto-enol tautomerization",
            "yield": "65-85%",
            "notes": "Markovnikov addition gives methyl ketone from terminal alkyne",
            "selectivity": "Terminal alkyne → methyl ketone",
        },
        {
            "method": "Hydroboration-oxidation (anti-Markovnikov)",
            "reagents": "1) Disiamylborane or R₂BH 2) H₂O₂, NaOH",
            "conditions": "0°C → RT, then basic workup",
            "mechanism": "Anti-Markovnikov addition of B-H, then oxidation",
            "yield": "60-80%",
            "notes": "Terminal alkyne → aldehyde (not ketone!)",
            "selectivity": "Anti-Markovnikov: terminal alkyne → aldehyde",
        },
    ],

    # === HALIDE conversions ===
    "alkyl_halide → alcohol": [
        {
            "method": "SN1/SN2 hydrolysis",
            "reagents": "Aqueous NaOH or AgNO₃(aq.) or moist Ag₂O",
            "conditions": "Heat or RT depending on substrate",
            "mechanism": "Nucleophilic substitution (SN2 for 1°, SN1 for 3°)",
            "yield": "60-90%",
            "notes": "Ag⁺ promotes ionization for unreactive halides",
            "selectivity": "Possible elimination as side reaction (especially with strong base)",
        },
    ],
    "alkyl_halide → alkene": [
        {
            "method": "Elimination (E2)",
            "reagents": "KOtBu or DBU or NaOEt, heat",
            "conditions": "Reflux, aprotic solvent (THF, DMSO)",
            "mechanism": "Concerted base-promoted β-H abstraction + X leaving",
            "yield": "60-90%",
            "notes": "Zaitsev product favored with small bases; Hofmann with bulky bases",
            "selectivity": "E2 requires antiperiplanar H-X arrangement",
        },
    ],
    "alkyl_halide → nitrile": [
        {
            "method": "SN2 cyanation",
            "reagents": "NaCN or KCN, DMSO or EtOH/H₂O",
            "conditions": "Reflux (for 1° halides primarily)",
            "mechanism": "CN⁻ as nucleophile in SN2 displacement",
            "yield": "70-90%",
            "notes": "Only works well for primary halides (steric hindrance for 2°/3°)",
            "selectivity": "SN2 with inversion of configuration",
        },
    ],

    # === AMINE conversions ===
    "amine → amide": [
        {
            "method": "Acylation",
            "reagents": "Acid chloride / anhydride, base (Et₃N, pyridine)",
            "conditions": "0°C → RT, CH₂Cl₂ or THF",
            "mechanism": "Nucleophilic acyl substitution (addition-elimination)",
            "yield": "70-95%",
            "notes": "Use Schotten-Baumann conditions for water-soluble amines",
            "selectivity": "Overacylation possible with primary amines (use 1 equiv)",
        },
    ],
}


# Normalized FG name mapping
FG_NAME_MAP: Dict[str, str] = {
    "alcohol": "alcohol", "oh": "alcohol", "r-oh": "alcohol",
    "aldehyde": "aldehyde", "rcho": "aldehyde", "cho": "aldehyde",
    "ketone": "ketone", "rcor": "ketone", "cor": "ketone", "c=o": "ketone",
    "carboxylic_acid": "carboxylic_acid", "cooh": "carboxylic_acid", "acid": "carboxylic_acid",
    "ester": "ester", "coor": "ester",
    "amide": "amide", "conh2": "amide", "conhr": "amide",
    "alkene": "alkene", "c=c": "alkene", "olefin": "alkene",
    "alkane": "alkane", "c-c": "alkane",
    "alkyne": "alkyne", "c#c": "alkyne",
    "alkyl_halide": "alkyl_halide", "halide": "alkyl_halide", "rx": "alkyl_halide",
    "nitrile": "nitrile", "cn": "nitrile",
    "amine": "amine", "nh2": "amine", "r-nh2": "amine",
    "diol": "diol", "epoxide": "epoxide",
    "aniline": "aniline", "phenol": "phenol",
}


@ChemMCPManager.register_tool
class FunctionalGroupInterconversion(BaseTool):
    __version__      = "0.1.0"
    name             = "FunctionalGroupInterconversion"
    func_name        = "convert_functional_group"
    description      = "Suggest step-by-step pathways for interconverting functional groups with reagents, conditions, mechanisms, and yields."
    implementation_description = "Queries a comprehensive database of named organic reactions covering major functional group transformations. Returns ranked pathways with practical details."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Functional Group", "Organic Synthesis", "Reaction Pathway", "Interconversion"]
    required_envs    = []

    code_input_sig   = [
        ("source_fg", "str", "N/A", "Source functional group name (e.g., 'alcohol', 'aldehyde', 'carboxylic_acid', 'alkene', 'alkyl_halide', 'amine')."),
        ("target_fg", "str", "N/A", "Target functional group name."),
        ("context_smiles", "str", "", "Optional SMILES string providing molecular context (affects selectivity notes)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'source_fg target_fg [context_SMILES]'. Example: 'alcohol aldehyde'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: source_fg, target_fg, pathways (list of {method, reagents, conditions, mechanism, yield, notes, selectivity}), conversion_possible, alternative_routes."),
    ]

    examples         = [
        {
            "code_input": {"source_fg": "alcohol", "target_fg": "aldehyde", "context_smiles": ""},
            "text_input": {"input_params": "alcohol aldehyde"},
            "output": {
                "result": {
                    "source_fg": "alcohol",
                    "target_fg": "aldehyde",
                    "conversion_possible": True,
                    "pathway_count": 3,
                    "best_method": "Dess-Martin periodinane (85-95% yield, very mild)",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, source_fg: str, target_fg: str, context_smiles: str = "") -> dict:
        """Core logic: look up FG interconversion pathways."""
        # Normalize names
        src_norm = FG_NAME_MAP.get(source_fg.lower().strip(), source_fg.lower().strip())
        tgt_norm = FG_NAME_MAP.get(target_fg.lower().strip(), target_fg.lower().strip())

        key = f"{src_norm} → {tgt_norm}"

        pathways = FG_INTERCONVERSION_DB.get(key)

        if not pathways:
            # Try reverse direction
            rev_key = f"{tgt_norm} → {src_norm}"
            if FG_INTERCONVERSION_DB.get(rev_key):
                return {
                    "source_fg": source_fg,
                    "target_fg": target_fg,
                    "conversion_possible": False,
                    "note": f"Direct {source_fg} → {target_fg} not found, but reverse ({target_fg} → {source_fg}) exists.",
                    "alternative_suggestion": f"Consider multi-step route or check if {target_fg} → {source_fg} → ... → {target_fg} is viable",
                    "pathways": [],
                }

            return {
                "source_fg": source_fg,
                "target_fg": target_fg,
                "conversion_possible": False,
                "note": f"No direct conversion pathway found for {source_fg} → {target_fg}",
                "available_source_groups": sorted(set(k.split(" → ")[0] for k in FG_INTERCONVERSION_DB.keys())),
                "available_target_groups": sorted(set(k.split(" → ")[1] for k in FG_INTERCONVERSION_DB.keys())),
                "pathways": [],
            }

        # Context-aware adjustments
        context_notes = []
        if context_smiles:
            import re as _re
            s = context_smiles
            if _re.search(r'C=C', s) and tgt_norm == "alkene":
                context_notes.append("⚠ Substrate already contains alkene: check chemoselectivity")
            if _re.search(r'[ONSFClBrI]', s) and src_norm == "alkyl_halide":
                context_notes.append("Note: Other heteroatoms present — consider protecting groups")

        best = max(pathways, key=lambda p: self._parse_yield_range(p.get("yield", "0%"))[1])

        return {
            "source_fg": source_fg,
            "target_fg": target_fg,
            "conversion_possible": True,
            "pathway_count": len(pathways),
            "pathways": pathways,
            "best_method": best["method"],
            "best_yield": best["yield"],
            "context_notes": context_notes,
            "alternative_routes": self._find_alternative_routes(src_norm, tgt_norm),
        }

    @staticmethod
    def _parse_yield_range(yield_str: str) -> Tuple[int, int]:
        """Parse yield range like '70-85%' into (low, high)."""
        try:
            nums = [int(x) for x in yield_str.replace("%", "").split("-")]
            return (nums[0], nums[-1])
        except Exception:
            return (0, 100)

    def _find_alternative_routes(self, src: str, tgt: str) -> List[str]:
        """Find indirect routes via intermediate functional groups."""
        intermediates = set()
        for k in FG_INTERCONVERSION_DB.keys():
            parts = k.split(" → ")
            if parts[0] == src:
                intermediates.add(parts[1])

        routes = []
        for inter in intermediates:
            if f"{inter} → {tgt}" in FG_INTERCONVERSION_DB:
                routes.append(f"{src} → {inter} → {tgt}")

        return routes[:5]

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            src = parts[0]
            tgt = parts[1] if len(parts) > 1 else ""
            ctx = parts[2] if len(parts) > 2 else ""
            return self._run_base(src, tgt, ctx)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'source_fg target_fg [context_SMILES]'")
