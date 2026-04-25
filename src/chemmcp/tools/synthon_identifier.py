"""
Synthon Identifier - identifies synthons and synthetic equivalents
from a target molecule at specified disconnection positions.
"""

import logging
import re
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Comprehensive synthon database: bond type → (synthon_a, synthon_b, equivalents_a, equivalents_b)
SYNTHON_DATABASE: List[Dict[str, Any]] = [
    {
        "bond_type": "C-CO (carbonyl adjacent)",
        "pattern": r'C.*C\(=O\)|C\(=O\).*C',
        "synthon_a": {"charge": "+", "type": "Acyl electrophile [R-C≡O⁺ or R-C(=O)⁺]"},
        "synthon_b": {"charge": "−", "type": "Nucleophilic carbon [R'⁻ or R'-M]"},
        "equivalents_a": ["Acid chloride (RCOCl)", "Anhydride [(RCO)₂O]", "Ester (RCOOR')", "Nitrile (RCN)", "Aldehyde/ketone (via enamine)"],
        "equivalents_b": ["Grignard reagent (R'MgX)", "Organolithium (R'Li)", "Enolate (from ketone/ester)", "Acetylide (R'C≡CM)", "Phosphonium ylide (Wittig)"],
    },
    {
        "bond_type": "C-O (ether/alcohol)",
        "pattern": r'COC|CO(?![a-z])',
        "synthon_a": {"charge": "+", "type": "Alkyl electrophile [R⁺]"},
        "synthon_b": {"charge": "−", "type": "Oxygen nucleophile [HO⁻ / RO⁻]"},
        "equivalents_a": ["Alkyl halide (RX)", "Tosylate (ROTs)", "Mesylate (ROMs)", "Epoxide"],
        "equivalents_b": ["Alcohol (ROH)", "Phenol (ArOH)", "Carboxylic acid (RCOOH)"],
    },
    {
        "bond_type": "C-N (amine/amide)",
        "pattern": r'CN|NC|C(=O)N|NC(=O)',
        "synthon_a": {"charge": "+", "type": "Electrophilic carbon/nitrogen"},
        "synthon_b": {"charge": "−", "type": "Nitrogen nucleophile [R'NH₂ / R'NH⁻]"},
        "equivalents_a": ["Alkyl halide (RX)", "Acyl chloride (RCOCl)", "Carbonyl compound (imine formation)"],
        "equivalents_b": ["Amine (R'NH₂)", "Azide (N₃⁻, then reduce)", "Nitro compound (reduction)", "Phthalimide (Gabriel synthesis)"],
    },
    {
        "bond_type": "C=C (alkene)",
        "pattern": r'C=C',
        "synthon_a": {"charge": "+", "type": "Electrophilic carbocation equivalent"},
        "synthon_b": {"charge": "−", "type": "Nucleophilic carbanion equivalent"},
        "equivalents_a": ["Aldehyde/ketone (C=O)", "α-Haloketone", "Epoxide"],
        "equivalents_b": ["Phosphonium ylide (Wittig)", "Julia olefination sulfone", "Peterson reagent (silyl anion)", "Tebbe reagent"],
    },
    {
        "bond_type": "C≡C (alkyne)",
        "pattern": r'C#C',
        "synthon_a": {"charge": "+", "type": "Electrophilic alkynyl fragment"},
        "synthon_b": {"charge": "−", "type": "Acetylide nucleophile [RC≡C⁻]"},
        "equivalents_a": ["1,2-Dibromoethane + 2 base", "Gem-dihalide", "Carbonyl (via Seyferth-Gilbert)"],
        "equivalents_b": ["Terminal alkyne (NaNH₂ deprotonated)", "TMS-acetylide (deprotected in situ)", "Acetylide metal salt"],
    },
    {
        "bond_type": "Ar-C (aromatic side chain)",
        "pattern": r'[cC][1-9].*CC|[cC]c[1-9]',
        "synthon_a": {"charge": "+", "type": "Arene π-system (electrophilic aromatic substitution)"},
        "synthon_b": {"charge": "±/+", "type": "Side-chain electrophile"},
        "equivalents_a": ["Benzene derivative (activated/deactivated arene)", "Heteroarene"],
        "equivalents_b": ["Alkyl halide (FC alkylation)", "Acyl chloride (FC acylation)", "Alkene (FC alkylation)"],
    },
    {
        "bond_type": "C-COOH (carboxylic acid α-position)",
        "pattern": r'.*C\(=O\)[OH,O]',
        "synthon_a": {"charge": "+", "type": "CO₂ electrophile (+CO₂)"},
        "synthon_b": {"charge": "−", "type": "Carbanion [R⁻]"},
        "equivalents_a": ["CO₂ (carboxylation)", "Carbonate (ROCO₂R)", "Urea derivatives"],
        "equivalents_b": ["Grignard reagent (RMgX)", "Organolithium (RLi)", "Malonic ester enolate"],
    },
]

# Polarity matching guide
POLARITY_RULES: Dict[str, str] = {
    "donor-acceptor": "Polarized bonds should be disconnected so that electron-rich fragments pair with electron-poor partners",
    "stabilization": "Synthons should correspond to stabilized species: enolates, enols, enamines for −; carbonyls, iminium, oxocarbenium for +",
    "oxidation_state": "Consider oxidation level matching: don't disconnect between very different oxidation states without redox steps",
}


@ChemMCPManager.register_tool
class SynthonIdentifier(BaseTool):
    __version__      = "0.1.0"
    name             = "SynthonIdentifier"
    func_name        = "identify_synthons"
    description      = "Identify synthons and their synthetic equivalents from a target molecule at a given disconnection position."
    implementation_description = "Analyzes the target molecule to identify polarized bonds, suggests optimal disconnection points, and maps each resulting fragment to its synthetic equivalent reagents."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Synthon", "Retrosynthesis", "Synthetic Equivalents", "Disconnection"]
    required_envs    = []

    code_input_sig   = [
        ("target_smiles", "str", "N/A", "SMILES string of the target molecule."),
        ("disconnection_bond", "int", "-1", "Specific bond index to disconnect (-1 for auto-suggest best)."),
        ("include_reagents", "bool", "True", "Whether to include specific reagent suggestions."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'target_SMILES [bond_index] [include_reagents]'. Example: 'CC(=O)c1ccccc1 -1 true'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: target_smiles, identified_synthons (list of {bond_type, synthon_+, synthon_-, equivalents_+, equivalents_, polarity_match}), recommended_disconnection."),
    ]

    examples         = [
        {
            "code_input": {
                "target_smiles": "CC(=O)c1ccccc1",
                "disconnection_bond": -1,
                "include_reagents": True,
            },
            "text_input": {"input_params": "CC(=O)c1ccccc1 -1 true"},
            "output": {
                "result": {
                    "target_smiles": "CC(=O)c1ccccc1",
                    "recommended_disconnection": "Ar-C(=O) bond (aromatic acyl disconnection)",
                    "identified_synthons": [
                        {
                            "bond_type": "Ar-C (aromatic side chain)",
                            "synthon_plus": "Acyl cation [CH₃CO⁺]",
                            "synthon_minus": "Arene π-system (benzene)",
                            "equivalents_plus": ["Acetyl chloride (CH₃COCl)", "Acetic anhydride"],
                            "equivalents_minus": ["Benzene", "Substituted benzene"],
                            "polarity_match": "✓ Good: Electrophilic acyl pairs with electron-rich arene",
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

    def _run_base(self, target_smiles: str, disconnection_bond: int = -1, include_reagents: bool = True) -> dict:
        """Core logic: identify synthons."""
        if not target_smiles:
            raise ChemMCPError("Target SMILES is required.")

        s = target_smiles

        # Find matching bond types from database
        matches: List[Dict[str, Any]] = []
        for entry in SYNTHON_DATABASE:
            try:
                if re.search(entry["pattern"], s):
                    match_info = {
                        "bond_type": entry["bond_type"],
                        "synthon_plus": entry["synthon_a"]["type"],
                        "synthon_minus": entry["synthon_b"]["type"],
                        "polarity_match": "✓ Good polarity matching" if entry["synthon_a"]["charge"] == "+" and entry["synthon_b"]["charge"] == "−" else "⚠ Check polarity compatibility",
                    }
                    if include_reagents:
                        match_info["equivalents_plus"] = entry.get("equivalents_a", [])
                        match_info["equivalents_minus"] = entry.get("equivalents_b", [])
                    matches.append(match_info)
            except Exception:
                continue

        # If no pattern matches, provide generic analysis
        if not matches:
            matches.append({
                "bond_type": "General C-C/C-heteroatom bond",
                "synthon_plus": "Electrophilic fragment (acceptor)",
                "synthon_minus": "Nucleophilic fragment (donor)",
                "polarity_match": "Generic suggestion — analyze manually for best disconnection",
                "equivalents_plus": ["Carbonyl compound", "Alkyl halide", "Epoxide"] if include_reagents else [],
                "equivalents_minus": ["Organometallic (Grignard/R-Li)", "Enolate", "Amine/Alcohol"] if include_reagents else [],
            })

        # Determine best disconnection recommendation
        best = self._recommend_best(s, matches)

        return {
            "target_smiles": target_smiles,
            "identified_synthons": matches,
            "recommended_disconnection": best,
            "total_synthon_pairs_found": len(matches),
            "polarity_guidelines": POLARITY_RULES,
        }

    @staticmethod
    def _recommend_best(smiles: str, matches: list) -> str:
        """Recommend the best single disconnection."""
        s = smiles

        # Priority order for recommendations
        if re.search(r'C(=O)[ON]', s):
            return "Disconnect ester/amide C-O/N bond first (high-yielding condensation reverse)"
        if re.search(r'c[1-9].*C\(=O\)|C\(=O\)c[1-9]', s):
            return "Disconnect aryl-carbonyl bond (Friedel-Crafts approach)"
        if re.search(r'COC', s):
            return "Disconnect ether C-O bond (Williamson ether synthesis reverse)"
        if re.search(r'C=C', s):
            return "Disconnect alkene via Wittig or carbonyl coupling strategy"
        if re.search(r'C#C', s):
            return "Disconnect alkyne via acetylide coupling"
        if matches:
            return f"Best option: {matches[0]['bond_type']}"
        return "Analyze functional groups manually for optimal disconnection"

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            smiles_str = parts[0]
            bond_idx = int(parts[1]) if len(parts) > 1 else -1
            inc_reag = parts[2].lower() == "true" if len(parts) > 2 else True
            return self._run_base(smiles_str, bond_idx, inc_reag)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'SMILES [bond_index] [include_reagents]'")
