"""
Retrosynthesis Analyzer - performs retrosynthetic analysis with bond disconnection strategies
to break down a target molecule into simpler precursors.
"""

import logging
import re
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Strategic bond disconnection rules
# Each rule: (priority, name, SMILES_pattern_hint, rationale)
DISCONNECTION_RULES: List[Dict[str, Any]] = [
    {"priority": 1, "name": "Adjacent to heteroatom (C-X)", "pattern_hint": "C-O,C-N,C-S", "rationale": "Heteroatom bonds are often formed in the last step"},
    {"priority": 2, "name": "α to carbonyl (C-C=O)", "pattern_hint": "CC(=O),C(=O)C", "rationale": "Carbonyl α-position is highly reactive for C-C bond formation"},
    {"priority": 3, "name": "Branched position (tertiary carbon)", "pattern_hint": "C(C)(C)C", "rationale": "Branch points often arise from alkylation or addition reactions"},
    {"priority": 4, "name": "Ring fusion / bridging bond", "pattern_hint": "C12CCCCC1.*C2", "rationale": "Ring-forming reactions are convergent and high-yielding"},
    {"priority": 5, "name": "Conjugated system (diene/ene-yne)", "pattern_hint": "C=CC=C,C#CC#C", "rationale": "Aldol, Wittig, or coupling reactions build conjugation"},
    {"priority": 6, "name": "Aromatic side-chain attachment", "pattern_hint": "c1ccccc1C", "rationale": "Friedel-Crafts, Suzuki, etc. attach side chains late"},
    {"priority": 7, "name": "Ester / amide bond", "pattern_hint": "C(=O)OC,C(=O)N", "rationale": "Condensation reactions are typically final steps"},
]

# Common synthons and their synthetic equivalents
SYNTHON_EQUIVALENTS: Dict[str, List[Dict[str, str]]] = {
    "-CH₂⁻ (nucleophilic)": [
        {"reagent": "CH₃I (methylation)", "reaction_type": "Alkylation"},
        {"reagent": "formaldehyde + base (homologation)", "reaction_type": "Aldol-type addition"},
        {"reagent": "epoxide ring opening", "reaction_type": "Nucleophilic substitution"},
    ],
    "-CHO (electrophilic carbonyl)": [
        {"reagent": "DMF (Vilsmeier formylation)", "reaction_type": "Formylation"},
        {"reagent": "CO/HCl, AlCl₃ (Gattermann-Koch)", "reaction_type": "Formylation"},
        {"reagent": "PCC/PDC oxidation of alcohol", "reaction_type": "Oxidation"},
    ],
    "-COCH₃ (acyl electrophile)": [
        {"reagent": "acetyl chloride / AlCl₃ (FC acylation)", "reaction_type": "Friedel-Crafts acylation"},
        {"reagent": "acetic anhydride", "reaction_type": "Acetylation"},
    ],
    "-C≡N (cyanide nucleophile)": [
        {"reagent": "KCN / NaCN", "reaction_type": "SN2 cyanation"},
        {"reagent": "TMSCN (trimethylsilyl cyanide)", "reaction_type": "Cyanation"},
    ],
    "=CH₂ (carbene equivalent)": [
        {"reagent": "CH₂I₂/Zn(Cu) (Simmons-Smith)", "reaction_type": "Cyclopropanation"},
        {"reagent": "CH₂N₂ (diazomethane)", "reaction_type": "Methylation / cycloaddition"},
    ],
    "+CH₂-CH₂⁺ (1,2-dication synthon)": [
        {"reagent": "BrCH₂CH₂Br (1,2-dibromoethane)", "reaction_type": "Double alkylation"},
        {"reagent": "ethylene oxide opening", "reaction_type": "Epoxide ring-opening"},
    ],
}

# Common building block molecules with SMILES
BUILDING_BLOCKS: List[Dict[str, Any]] = [
    {"name": "Ethylene glycol", "smiles": "OCCO", "role": "1,2-diol synthon"},
    {"name": "Benzaldehyde", "smiles": "O=Cc1ccccc1", "role": "aromatic aldehyde donor"},
    {"name": "Acetone", "smiles": "CC(=O)C", "role": "methyl ketone / 3-carbon unit"},
    {"name": "Diethyl malonate", "smiles": "CCOC(=O)CC(=O)OCC", "role": "active methylene source"},
    {"name": "Cyclohexanone", "smiles": "O=C1CCCCC1", "role": "6-membered cyclic ketone"},
    {"name": "Phenylmagnesium bromide", "smiles": "c1ccccc[Mg]Br", "role": "aryl nucleophile (Grignard)"},
    {"name": "Ethyl acetoacetate", "smiles": "CC(=O)CC(=O)OCC", "role": "β-keto ester (active methylene)"},
    {"name": "Nitrobenzene", "smiles": "O=[N+]([O-])c1ccccc1", "role": "nitroarene → aniline precursor"},
]


@ChemMCPManager.register_tool
class RetrosynthesisAnalyzer(BaseTool):
    __version__      = "0.1.0"
    name             = "RetrosynthesisAnalyzer"
    func_name        = "analyze_retrosynthesis"
    description      = "Perform retrosynthetic analysis on a target molecule using strategic bond disconnection rules."
    implementation_description = "Applies Corey's disconnection principles: prioritize bonds adjacent to functional groups, then branched positions, rings, and conjugated systems. Builds a retrosynthetic tree with suggested synthons and starting materials."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Retrosynthesis", "Disconnection", "Synthesis Planning", "Organic Synthesis"]
    required_envs    = []

    code_input_sig   = [
        ("target_smiles", "str", "N/A", "SMILES string of the target molecule."),
        ("max_depth", "int", "3", "Maximum depth of retrosynthetic analysis (number of disconnection rounds)."),
        ("focus_strategy", "str", "auto", "Focus strategy: 'auto', 'heteroatom', 'carbonyl', 'ring', 'conjugation'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'target_SMILES [max_depth] [strategy]'. Example: 'CC(=O)c1ccccc1 3 auto'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: target_smiles, retrosynthetic_tree (nested dict per level), suggested_disconnections (list of {level, bond, rationale, synthons, equivalents}), total_steps_estimate."),
    ]

    examples         = [
        {
            "code_input": {
                "target_smiles": "CC(=O)c1ccccc1",
                "max_depth": 2,
                "focus_strategy": "auto",
            },
            "text_input": {"input_params": "CC(=O)c1ccccc1 2 auto"},
            "output": {
                "result": {
                    "target_smiles": "CC(=O)c1ccccc1",
                    "total_steps_estimate": 2,
                    "suggested_disconnections": [
                        {
                            "level": 1,
                            "disconnection": "C(arom)-C(=O) bond",
                            "rationale": "Aromatic side-chain attachment: Friedel-Crafts acylation of benzene",
                            "synthons": ["Acyl cation [CH₃CO⁺]", "Aromatic π-system"],
                            "equivalents": ["Acetyl chloride / AlCl₃", "Benzene"],
                        },
                        {
                            "level": 2,
                            "disconnection": "Acetyl chloride from carboxylic acid",
                            "rationale": "Acid chloride preparation from acetic acid",
                            "synthones": ["Acetyl group", "Chloride"],
                            "equivalents": ["Acetic acid + SOCl₂ / PCl₃", ""],
                        }
                    ],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, target_smiles: str, max_depth: int = 3, focus_strategy: str = "auto") -> dict:
        """Core logic: retrosynthetic analysis."""
        if not target_smiles:
            raise ChemMCPError("Target SMILES is required.")

        max_depth = min(max(1, max_depth), 5)  # limit to reasonable depth

        # Analyze target structure features
        features = self._analyze_features(target_smiles)

        # Determine strategy
        if focus_strategy == "auto":
            focus_strategy = self._select_strategy(features)

        # Build retrosynthetic tree
        tree = self._build_tree(target_smiles, max_depth, focus_strategy, features)

        return {
            "target_smiles": target_smiles,
            "target_features": features,
            "strategy_used": focus_strategy,
            "retrosynthetic_tree": tree["tree"],
            "suggested_disconnections": tree["disconnections"],
            "total_steps_estimate": len(tree["disconnections"]),
        }

    def _analyze_features(self, smiles: str) -> Dict[str, Any]:
        """Identify structural features in the target molecule."""
        s = smiles
        return {
            "has_aromatic_ring": bool(re.search(r'c[1-9]', s)) or bool(re.search(r'c1.*c1', s)),
            "has_carbonyl": bool(re.search(r'C(=O)', s)),
            "has_alcohol": bool(re.search(r'CO(?![a-z])', s)) or ('CO' in s),
            "has_amine": bool(re.search(r'NC|N\(', s)),
            "has_ester": bool(re.search(r'C(=O)OC|OC(=O)', s)),
            "has_amide": bool(re.search(r'C(=O)N|NC(=O)', s)),
            "has_double_bond": bool(re.search(r'C=C', s)),
            "has_triple_bond": bool(re.search(r'C#C', s)),
            "has_ring": bool(re.search(r'[0-9]', s)),
            "ring_count": len(re.findall(r'[0-9]', s)) if re.search(r'[0-9]', s) else 0,
            "carbon_count": len(re.findall(r'(?<![a-z])C(?![a-z])', s)) + len(re.findall(r'c', s)),
            "heteroatoms": list(set(re.findall(r'[ONSPFClBrI](?![a-z])', s))),
        }

    def _select_strategy(self, features: dict) -> str:
        """Auto-select best disconnection strategy based on features."""
        if features.get("has_ester") or features.get("has_amide"):
            return "heteroatom"
        if features.get("has_carbonyl") and not features.get("has_aromatic_ring"):
            return "carbonyl"
        if features.get("has_aromatic_ring") and features.get("carbon_count", 0) > 8:
            return "ring"
        if features.get("has_double_bond") or features.get("has_triple_bond"):
            return "conjugation"
        return "heteroatom"  # default

    def _build_tree(self, smiles: str, depth: int, strategy: str, features: dict) -> dict:
        """Build retrosynthetic tree recursively."""
        disconnections: List[Dict[str, Any]] = []
        current_target = smiles

        for level in range(1, depth + 1):
            disc = self._suggest_disconnection(current_target, level, strategy, features)
            if not disc:
                break
            disconnections.append(disc)
            current_target = disc.get("precursor_smiles", current_target)

        tree = {
            "target": smiles,
            "levels": [],
        }
        prev = smiles
        for i, d in enumerate(disconnections):
            tree["levels"].append({
                "level": i + 1,
                "target": prev,
                "disconnection": d["bond"],
                "rationale": d["rationale"],
                "synthons": d["synthons"],
                "synthetic_equivalents": d["equivalents"],
                "precursor": d.get("precursor", "starting material"),
            })
            prev = d.get("precursor", prev)

        return {"tree": tree, "disconnections": disconnections}

    def _suggest_disconnection(self, smiles: str, level: int, strategy: str, features: dict) -> Optional[Dict[str, Any]]:
        """Suggest one disconnection at given level."""
        s = smiles

        # Strategy-based suggestions
        suggestions: List[Dict[str, Any]] = []

        if strategy in ("heteroatom", "auto"):
            # Check ester/amide
            if re.search(r'C(=O)O[C]', s):
                suggestions.append({
                    "bond": "Ester C-O bond disconnection",
                    "rationale": "Disconnect ester → carboxylic acid + alcohol (or acyl halide + alcohol)",
                    "synthons": ["R-COOH (or R-COCl)", "R'-OH"],
                    "equivalents": ["Carboxylic acid (→ acyl halide)", "Alcohol"],
                    "precursor": "Carboxylic acid + Alcohol",
                })

            if re.search(r'C(=O)N', s):
                suggestions.append({
                    "bond": "Amide C-N bond disconnection",
                    "rationale": "Disconnect amide → amine + carboxylic acid derivative",
                    "synthons": ["R-COOH (or R-COCl)", "R'-NH₂"],
                    "equivalents": ["Acyl halide / activated ester", "Amine"],
                    "precursor": "Amine + Acylating agent",
                })

            # Check ether/alcohol C-O
            if re.search(r'COC(?![a-z])', s) or re.search(r'COc', s):
                suggestions.append({
                    "bond": "Ether/alcohol C-O bond disconnection",
                    "rationale": "Williamson ether synthesis (reverse): alkoxide + alkyl halide",
                    "synthons": ["R-O⁻ (alkoxide)", "R'-X (alkyl halide)"],
                    "equivalents": ["Alcohol (deprotonated)", "Primary alkyl halide"],
                    "precursor": "Alcohol + Alkyl halide",
                })

        if strategy in ("carbonyl", "auto"):
            # Aryl-C=O (ketone attached to aromatic)
            if re.search(r'c[1-9].*C\(=O\)|C\(=O\)c[1-9]', s):
                suggestions.append({
                    "bond": "Ar-C(=O) bond disconnection (Friedel-Crafts)",
                    "rationale": "Aryl ketone via Friedel-Crafts acylation of arene",
                    "synthons": ["Acyl cation [R-CO⁺]", "Arene (π-nucleophile)"],
                    "equivalents": ["Acyl chloride / AlCl₃", "Arene (benzene derivative)"],
                    "precursor": "Arene + Acyl chloride (Lewis acid)",
                })

            # Aliphatic ketone α-position
            if re.search(r'CC\(=O\)C|C\(=O\)CC', s) and not re.search(r'c[1-9]', s):
                suggestions.append({
                    "bond": "Ketone α-C-C bond disconnection",
                    "rationale": "Ketone formation via alkylation of enolate or Grignard addition to nitrile/acid",
                    "synthons": ["Enolate (nucleophile)", "Alkyl halide (electrophile)"],
                    "equivalents": ["Ketone enolate (LDA base)", "Alkyl halide"],
                    "precursor": "Smaller ketone + Alkyl halide (enolate alkylation)",
                })

        if strategy in ("ring", "auto") and features.get("has_ring"):
            suggestions.append({
                "bond": "Ring-forming bond disconnection",
                "rationale": "Open ring at strategic position; consider intramolecular aldol, Dieckmann, or cyclization",
                "synthons": ["Dinucleophile fragment", "Dielectrophile fragment"],
                "equivalents": ["Dicarbonyl compound", "Dihalide / diacid derivative"],
                "precursor": "Linear precursor suitable for ring-closing reaction",
            })

        if strategy in ("conjugation", "auto"):
            if re.search(r'C=CC=C', s):
                suggestions.append({
                    "bond": "Conjugated diene central bond disconnection",
                    "rationale": "Build conjugation via aldol condensation or Wittig reaction",
                    "synthons": ["Enone/Enal synthon", "Wittig ylide / Enolate"],
                    "equivalents": ["α,β-Unsaturated carbonyl compound", "Phosphonium ylide or enolate"],
                    "precursor": "Two smaller carbonyl compounds (aldol/Wittig)",
                })

            if re.search(r'C#C', s):
                suggestions.append({
                    "bond": "Alkyne bond disconnection",
                    "rationale": "Form alkyne via alkylation of acetylide or elimination",
                    "synthons": ["Acetylide anion", "Alkyl halide"],
                    "equivalents": ["Terminal alkyne (NaNH₂)", "Alkyl halide"],
                    "precursor": "Terminal alkyne + Alkyl halide (SN2)",
                })

        # Fallback: generic suggestion
        if not suggestions:
            suggestions.append({
                "bond": "General C-C bond disconnection",
                "rationale": "Apply standard disconnection: identify most reactive/best-stabilized synthons",
                "synthons": ["Nucleophilic synthon", "Electrophilic synthon"],
                "equivalents": ["Organometallic reagent (Grignard, organolithium)", "Carbonyl / alkyl halide"],
                "precursor": "Two simpler fragments joined by C-C bond forming reaction",
            })

        # Return top suggestion for this level (rotate through options)
        idx = (level - 1) % len(suggestions)
        return suggestions[idx]

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            smiles_str = parts[0]
            depth = int(parts[1]) if len(parts) > 1 else 3
            strategy = parts[2] if len(parts) > 2 else "auto"
            return self._run_base(smiles_str, depth, strategy)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'SMILES [depth] [strategy]'")
