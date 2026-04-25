"""
Disconnection Suggester - suggests optimal bond disconnection positions
in a target molecule with rationale and difficulty ratings.
"""

import logging
import re
from typing import Dict, List, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Disconnection scoring criteria
# Each criterion: (name, pattern, score_weight)
DISCONNECTION_CRITERIA: List[Dict[str, Any]] = [
    {
        "name": "Adjacent to heteroatom (C-O, C-N, C-S)",
        "patterns": [r'C(?=[ONSFClBrI])', r'(?<=[CONS])C', r'COC', r'CNC'],
        "score": 10,
        "rationale": "Heteroatom bonds are often formed in the final step (high reliability)",
        "difficulty": "Easy",
        "typical_reaction": "Nucleophilic substitution / Condensation",
    },
    {
        "name": "α to carbonyl (C-C=O)",
        "patterns": [r'CC\(=O\)', r'C\(=O\)CC', r'C\(=O\)C1', r'1C\(=O\)'],
        "score": 9,
        "rationale": "Carbonyl α-position is highly reactive for enolate chemistry",
        "difficulty": "Moderate",
        "typical_reaction": "Enolate alkylation / Aldol / Claisen condensation",
    },
    {
        "name": "β to carbonyl (conjugated position)",
        "patterns": [r'C=CC\(=O\)', r'C\(=O\)C=C', r'C=CC(=O)'],
        "score": 8,
        "rationale": "Conjugation allows Michael addition or conjugate addition strategies",
        "difficulty": "Moderate",
        "typical_reaction": "Michael addition / Robinson annulation / Wittig-Horner",
    },
    {
        "name": "Benzylic position (Ar-CH₂)",
        "patterns": [r'[cC][1-9].*C[^0-9a-z]|[cC]c[1-9]C'],
        "score": 8,
        "rationale": "Benzylic carbocations/stabilized anions are accessible",
        "difficulty": "Moderate-Easy",
        "typical_reaction": "Friedel-Crafts / Benzylic halogenation → substitution",
    },
    {
        "name": "Allylic position (C=C-C)",
        "patterns": [r'C=CC', r'C=C[C]'],
        "score": 7,
        "rationale": "Allylic systems allow SN2' or rearrangement pathways",
        "difficulty": "Moderate",
        "typical_reaction": "Allylic halogenation / SN2' / Rearrangement",
    },
    {
        "name": "Tertiary carbon (branched)",
        "patterns": [r'C\(C\)\(C\)C', r'C\(C\)C\(C\)'],
        "score": 6,
        "rationale": "Branch points suggest convergent assembly from smaller fragments",
        "difficulty": "Harder",
        "typical_reaction": "Tertiary alkylation / Multi-step fragment coupling",
    },
    {
        "name": "Ring bond (cyclization point)",
        "patterns": [r'C[0-9]', r'[0-9]C'],
        "score": 5,
        "rationale": "Ring-opening reveals linear precursor; consider ring-closing strategy",
        "difficulty": "Variable",
        "typical_reaction": "Intramolecular aldol / Dieckmann / Lactamization / RCM",
    },
    {
        "name": "Ester/Amide bond (condensation)",
        "patterns": [r'C\(=O\)OC', r'C\(=O\)N', r'OC\(=O\)', r'NC\(=O\)'],
        "score": 10,
        "rationale": "Condensation bonds are reliable and high-yielding to disconnect",
        "difficulty": "Easy",
        "typical_reaction": "Esterification / Amidation / DCC coupling",
    },
    {
        "name": "Aromatic C-X bond",
        "patterns": [r'[cC][1-9][ONSFClBrI]', r'[ONSFClBrI][cC][1-9]'],
        "score": 4,
        "rationale": "Direct aromatic substitution may be challenging; check directing groups",
        "difficulty": "Harder",
        "typical_reaction": "SEAr (requires activating group) / Metal-catalyzed cross-coupling",
    },
]

# Difficulty rating scale
DIFFICULTY_DETAILS: Dict[str, Dict[str, str]] = {
    "Easy": {"yield_range": "70-95%", "steps_needed": "1-2", "reliability": "High", "note": "Standard textbook reactions"},
    "Moderate": {"yield_range": "50-80%", "steps_needed": "2-4", "reliability": "Medium-High", "note": "Well-established but needs optimization"},
    "Moderate-Easy": {"yield_range": "60-85%", "steps_needed": "2-3", "reliability": "High", "note": "Reliable with proper conditions"},
    "Harder": {"yield_range": "30-70%", "steps_needed": "3-6", "reliability": "Medium", "note": "May require protecting groups or special conditions"},
    "Variable": {"yield_range": "20-80%", "steps_needed": "3-8", "reliability": "Variable", "note": "Highly substrate-dependent; test multiple conditions"},
}


@ChemMCPManager.register_tool
class DisconnectionSuggester(BaseTool):
    __version__      = "0.1.0"
    name             = "DisconnectionSuggester"
    func_name        = "suggest_disconnections"
    description      = "Suggest and rank optimal bond disconnection positions in a target molecule with rationale, expected synthons, and difficulty ratings."
    implementation_description = "Scans the target molecule against strategic disconnection criteria, scores each potential disconnection site, ranks by synthetic feasibility, and provides reaction type suggestions."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Disconnection", "Retrosynthesis", "Synthetic Planning", "Bond Analysis"]
    required_envs    = []

    code_input_sig   = [
        ("target_smiles", "str", "N/A", "SMILES string of the target molecule."),
        ("strategy", "str", "auto", "Strategy: 'auto', 'functional_group', 'heteroatom', 'ring', 'carbonyl'."),
        ("max_suggestions", "int", "5", "Maximum number of disconnection suggestions to return."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'target_SMILES [strategy] [max_suggestions]'. Example: 'CC(=O)OCC auto 5'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: target_smiles, ranked_disconnections (list of {rank, position, bond_type, score, rationale, difficulty, synthons, typical_reaction}), strategy_used."),
    ]

    examples         = [
        {
            "code_input": {
                "target_smiles": "CC(=O)Oc1ccccc1",
                "strategy": "auto",
                "max_suggestions": 3,
            },
            "text_input": {"input_params": "CC(=O)Oc1ccccc1 auto 3"},
            "output": {
                "result": {
                    "target_smiles": "CC(=O)Oc1ccccc1",
                    "strategy_used": "auto",
                    "ranked_disconnections": [
                        {
                            "rank": 1,
                            "position": "Ester bond C(=O)-O",
                            "bond_type": "Ester/Amide bond (condensation)",
                            "score": 10,
                            "rationale": "Condensation bonds are reliable and high-yielding to disconnect",
                            "difficulty": "Easy",
                            "synthons": ["Acetate synthon [CH₃COO⁻]", "Phenol synthon [Ph-O⁻]"],
                            "typical_reaction": "Esterification (Fischer) / Steglich esterification",
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

    def _run_base(self, target_smiles: str, strategy: str = "auto", max_suggestions: int = 5) -> dict:
        """Core logic: suggest and rank disconnections."""
        if not target_smiles:
            raise ChemMCPError("Target SMILES is required.")

        s = target_smiles
        max_suggestions = min(max(1, max_suggestions), 15)

        # Score all applicable criteria
        scored: List[Dict[str, Any]] = []
        for crit in DISCONNECTION_CRITERIA:
            # Filter by strategy
            strat_name = crit["name"].lower()
            if strategy != "auto":
                if strategy == "functional_group" and "carbonyl" not in strat_name and "heteroatom" not in strat_name and "ester" not in strat_name:
                    # Still include FG-related ones
                    pass
                if strategy == "heteroatom" and not any(x in strat_name for x in ["heteroatom", "ester", "amide"]):
                    if "heteroatom" not in crit["name"].lower() and "ester" not in crit["name"].lower() and "amide" not in crit["name"].lower():
                        continue
                if strategy == "ring" and "ring" not in strat_name.lower():
                    continue
                if strategy == "carbonyl" and "carbonyl" not in strat_name.lower() and "α" not in crit["name"] and "β" not in crit["name"]:
                    if "carbonyl" not in crit["name"].lower():
                        continue

            matched_patterns = []
            for pat in crit["patterns"]:
                try:
                    if re.search(pat, s):
                        matched_patterns.append(pat)
                except Exception:
                    pass

            if matched_patterns:
                diff_info = DIFFICULTY_DETAILS.get(crit["difficulty"], {})
                scored.append({
                    "name": crit["name"],
                    "score": crit["score"],
                    "rationale": crit["rationale"],
                    "difficulty": crit["difficulty"],
                    "difficulty_details": diff_info,
                    "typical_reaction": crit["typical_reaction"],
                    "matched_patterns": len(matched_patterns),
                })

        # Sort by score descending
        scored.sort(key=lambda x: x["score"], reverse=True)

        # Take top N suggestions
        top = scored[:max_suggestions]

        # Generate synthons for each suggestion
        ranked = []
        for i, item in enumerate(top):
            synthons = self._generate_synthons(s, item["name"])
            ranked.append({
                "rank": i + 1,
                "position": item["name"],
                "bond_type": item["name"],
                "score": item["score"],
                "rationale": item["rationale"],
                "difficulty": item["difficulty"],
                "yield_range": item["difficulty_details"].get("yield_range", "Unknown"),
                "steps_estimate": item["difficulty_details"].get("steps_needed", "?"),
                "synthons": synthons,
                "typical_reaction": item["typical_reaction"],
            })

        return {
            "target_smiles": target_smiles,
            "strategy_used": strategy,
            "total_sites_found": len(scored),
            "ranked_disconnections": ranked,
        }

    @staticmethod
    def _generate_synthons(smiles: str, bond_type: str) -> List[str]:
        """Generate suggested synthons for a given disconnection type."""
        bt = bond_type.lower()

        if "heteroatom" in bt or "ester" in bt or "amide" in bt:
            return ["Electrophilic fragment (R-X or R-CO-X)", "Nucleophilic heteroatom fragment (HO⁻, RO⁻, H₂N⁻)"]
        if "carbonyl" in bt and "α" in bt:
            return ["Enolate nucleophile (stabilized carbanion)", "Electrophile (alkyl halide, carbonyl)"]
        if "β" in bt or "conjugat" in bt:
            return ["Michael donor (enolate / cuprate)", "Michael acceptor (α,β-unsaturated carbonyl)"]
        if "benzylic" in bt:
            return ["Benzylic cation equivalent / benzylic anion", "Arene + electrophile/nucleophile"]
        if "allylic" in bt:
            return ["Allylic electrophile (halide / ester)", "Nucleophile (soft/hard depending on conditions)"]
        if "tertiary" in bt or "branch" in bt:
            return ["Tertiary carbocation equivalent", "Tertiary carbanion / organometallic"]
        if "ring" in bt:
            return ["Dinucleophile fragment", "Dielectrophile fragment (for ring closure)"]
        if "aromatic" in bt:
            return ["Arene π-system (activated/deactivated)", "Side-chain electrophile"]

        return ["Electrophilic synthon (+)", "Nucleophilic synthon (−)"]

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            smiles_str = parts[0]
            strategy = parts[1] if len(parts) > 1 else "auto"
            max_n = int(parts[2]) if len(parts) > 2 else 5
            return self._run_base(smiles_str, strategy, max_n)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'SMILES [strategy] [max_suggestions]'")
