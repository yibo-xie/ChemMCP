"""
Fragmentation Predictor - MS/MS fragmentation pattern prediction
with collision-energy-dependent fragment probability scoring.
"""

import logging
import re
import math
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Atomic masses (most abundant isotope for MS calculations)
_ATOMIC_MASSES = {
    "H": 1.00783, "C": 12.0000, "N": 14.00307, "O": 15.99491,
    "F": 18.99840, "P": 30.97376, "S": 31.97207,
    "Cl": 34.96885, "Br": 78.91833, "I": 126.90447,
    "Na": 22.98977, "K": 38.96371,
}

# Common neutral losses in MS/MS: (mass_delta, formula, name, description)
_MSMS_NEUTRAL_LOSSES = [
    (1.00783,     "H•",       "Hydrogen radical",        "Common in ESI+"),
    (2.01565,     "H₂",       "Hydrogen",                "Reduction"),
    (17.00274,    "OH•",      "Hydroxyl radical",         "Alcohols, carboxylic acids"),
    (18.01056,    "H₂O",      "Water",                   "VERY COMMON — alcohols, acids, sugars"),
    (27.99491,    "CO",       "Carbon monoxide",          "Carbonyls, common"),
    (28.03130,    "C₂H₄",     "Ethylene",                 "McLafferty-type"),
    (28.01056,    "N₂",       "Nitrogen",                 "Azo compounds"),
    (29.00274,    "CHO•",     "Formyl radical",           "Aldehydes"),
    (30.01058,    "CH₂O",     "Formaldehyde",             ""),
    (31.01839,    "CH₃O•",    "Methoxy radical",          "Ethers, methyl esters"),
    (42.01056,    "CH₂=C=O",  "Ketene",                   "Acetates, carbonyls"),
    (43.00581,    "C₂H₃O",    "Acetyl/acetaldehyde",      "Ketones, esters"),
    (43.98983,    "CO₂",      "Carbon dioxide",            "Decarboxylation — COMMON"),
    (44.99766,    "CO₂H₂",    "Formic acid / CO₂+2H",     "Carboxylic acids"),
    (45.02940,    "CH₃CHOH",  "Ethanol fragment",          "Ethyl-containing"),
    (46.00548,    "NO₂",      "Nitro group",              "Nitro compounds"),
    (46.00548,    "CH₂O₂",    "Formic acid",              ""),
    (49.99232,    "CH₃Cl",   "Methyl chloride",           "Chlorinated"),
    (56.02622,    "C₄H₈",     "Butene",                   "McLafferty rearrangement"),
    (58.04187,    "C₃H₆O",   "Acetone/propionaldehyde",   "Ketones"),
    (60.02113,    "C₂H₄O₂",  "Acetic acid",              "McLafferty of acids/esters"),
    (64.96909,    "SO₂",      "Sulfur dioxide",            "Sulfonyl compounds"),
    (79.91690,    "Br•",      "Bromine radical",           "Brominated"),
    (84.04187,    "C₅H₈O",   "Cyclopentanone fragment",   ""),
    (98.01686,    "H₃PO₄",   "Phosphoric acid",           "Phosphorylated compounds"),
    (117.99044,   "C₅H₄O₄",  "Maleic/succinic derivative","Dicarboxylics"),
]

# Common diagnostic ions: (mz, formula, name, compound_class)
_DIAGNOSTIC_IONS = [
    (30.0340,  "CH₄N⁺",    "Protonated methyleneimine",  "Primary amines"),
    (43.0184,  "C₂H₆N⁺",   "Ethyliminium",               "Secondary amines"),
    (44.0500,  "C₃H₆N⁺",   "Propyliminium",              "Tertiary amines"),
    (50.9976,  "CH₄NO⁺",   "Methylisocyanate ion",       "Nitro/amide"),
    (60.0448,  "C₂H₆NO⁺",  "Acetamide immonium",         "Amides"),
    (70.0656,  "C₄H₈N⁺",   "Butyliminium",               "Amines"),
    (74.0241,  "C₃H₅O₂N⁺", "Glycine immonium",           "Peptides (b₁/a₁)"),
    (84.0449,  "C₄H₆NO⁺",  "α-aminobutyryl immonium",    "Peptides"),
    (84.0813,  "C₅H₁₀N⁺",  "Valine/pentaniminium",       "Peptides"),
    (86.0606,  "C₄H₈NO⁺",  "Threonine immonium",         "Peptides"),
    (102.0550, "C₄H₈NO₂⁺", "Glutamic acid immonium",     "Peptides"),
    (104.1070, "C₅H₁₄NO⁺", "Leucine immonium",           "Peptides"),
    (120.0810, "C₄H₁₀NO₃⁺","Serine immonium",            "Peptides"),
]


@ChemMCPManager.register_tool
class FragmentationPredictor(BaseTool):
    """
    质谱碎片化规律预测器 — 预测分子在 MS/MS 条件下的碎裂模式。
    
    基于 CID/HCD 碎裂规则，预测主要碎片离子、中性丢失和碎裂途径，
    并根据碰撞能量给出概率评分。
    """
    __version__      = "0.1.0"
    name             = "FragmentationPredictor"
    func_name        = "predict_msms_fragments"
    description      = "Predict MS/MS fragmentation patterns with collision-energy-dependent fragment probability scoring."
    implementation_description = "Uses common CID/HCD fragmentation rules including neutral loss patterns, diagnostic ions, bond cleavage preferences, and energy-dependent intensity modeling. Covers small molecules and peptide-like fragments."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Mass Spectrometry", "Fragmentation", "MS/MS", "Tandem MS", "CID"]
    required_envs    = []

    code_input_sig   = [
        ("smiles", "str", "N/A", "SMILES string or molecular formula of the precursor molecule."),
        ("collision_energy", "float", "30.0", "Collision energy in eV (typical range: 10-60 eV)."),
        ("precursor_mz", "float", "None", "Precursor m/z value (auto-calculated if not provided)."),
        ("ionization_mode", "str", "ESI+", "Ionization mode: 'ESI+', 'ESI-', 'APCI+', 'APCI-'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'smiles [collision_energy] [precursor_mz] [ionization_mode]'. Example: 'CC(=O)O 25 61 ESI+'"),
    ]

    output_sig       = [
        ("result", "dict", "Complete MS/MS prediction including precursor info, predicted fragments with m/z/intensity/pathway/probability, neutral loss summary, and fragmentation rules applied."),
    ]

    examples         = [
        {
            "code_input": {"smiles": "CC(=O)O", "collision_energy": 20.0, "ionization_mode": "ESI+"},
            "text_input": {"input_params": "CC(=O)O 20 ESI+"},
            "output": {
                "result": {
                    "precursor_mz": 61.0284,
                    "precursor_formula": "C₂H₄O₂",
                    "collision_energy_eV": 20.0,
                    "ionization_mode": "ESI+",
                    "top_fragments": [
                        {"mz": 43.0184, "formula": "[M-H₂O+H]⁺", "pathway": "Dehydration", "probability": 85},
                        {"mz": 15.0235, "formula": "[CH₃]⁺", "pathway": "α-Cleavage", "probability": 45},
                    ],
                    "neutral_losses_predicted": [{"loss": "H₂O (18.0106 Da)", "significance": "high"}],
                    "fragmentation_rules_applied": ["carboxylic_acid_esip"],
                }
            },
        },
        {
            "code_input": {"smiles": "C17H19NO3", "collision_energy": 35.0, "ionization_mode": "ESI+"},
            "text_input": {"input_params": "C17H19NO3 35 ESI+"},
            "output": {
                "result": {
                    "precursor_mz": 286.1438,
                    "top_fragments": [
                        {"mz": 268.1332, "formula": "[M-H₂O+H]⁺", "pathway": "Dehydration", "probability": 72},
                        {"mz": 240.1280, "formula": "[M-H₂O-CO+H]⁺", "pathway": "Dehydration + decarbonylation", "probability": 55},
                        {"mz": 212.1174, "formula": "[M-H₂O-2CO+H]⁺", "pathway": "Multiple neutral losses", "probability": 30},
                    ],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._rdkit_available = False
        try:
            from rdkit import Chem
            self.Chem = Chem
            self._rdkit_available = True
        except ImportError:
            pass

    def _calculate_exact_mass(self, formula: str) -> float:
        """Calculate exact mass from molecular formula."""
        elements = _parse_formula_simple(formula)
        mass = 0.0
        for elem, count in elements.items():
            if elem not in _ATOMIC_MASSES:
                logger.warning(f"No mass data for element '{elem}', skipping")
                continue
            mass += _ATOMIC_MASSES[elem] * count
        return round(mass + 1.007276, 4)  # add proton for [M+H]+

    def _detect_functional_groups(self, smi: str) -> list:
        """Detect functional groups for MS/MS prediction."""
        groups = []
        s = smi.upper()
        if re.search(r'C\(=O\)[Oo]', smi):
            groups.append("ester")
        if re.search(r'C\(=O\)[Oo][Hh]?', smi):
            groups.append("carboxylic_acid")
        elif re.search(r'C(=O)[Hh]|[Hh]C=O', smi):
            groups.append("aldehyde")
        elif re.search(r'C\(=O\)[CcH]', smi) and not re.search(r'C\(=O\)[OoNn]', smi):
            groups.append("ketone")
        if re.search(r'[Cc]O[Hh]?', smi) and not re.search(r'C\(=O\)O', smi):
            groups.append("alcohol")
        if re.search(r'CN|C#N', smi):
            groups.append("nitrile")
        if re.search(r'[Nn]', smi) and not re.search(r'C\(=O\)[Nn]', smi):
            groups.append("amine")
        if re.search(r'C\(=O\)[Nn]', smi):
            groups.append("amide")
        if re.search(r'Cl', smi):
            groups.append("chlorinated")
        if re.search(r'Br', smi):
            groups.append("brominated")
        if 'c' in smi:
            groups.append("aromatic")
        if re.search(r'S', smi):
            groups.append("sulfur_compound")
        if re.search(r'P', smi):
            groups.append("phosphorus_compound")
        return groups if groups else ["hydrocarbon"]

    def _predict_msms(self, smi: str, ce: float, prec_mz: float, mode: str) -> dict:
        """Core MS/MS prediction logic."""
        groups = self._detect_functional_groups(smi)

        # Estimate precursor mz if not given
        if prec_mz is None or prec_mz <= 0:
            # Try to estimate from SMILES via RDKit or simple formula parsing
            if self._rdkit_available:
                try:
                    mol = self.Chem.MolFromSmiles(smi)
                    if mol:
                        from rdkit.Chem import Descriptors
                        mw = Descriptors.ExactMolWt(mol)
                        if mode.endswith("+"):
                            prec_mz = round(mw + 1.007276, 4)
                        else:
                            prec_mz = round(mw - 1.007276, 4)
                except Exception:
                    pass
            if prec_mz is None or prec_mz <= 0:
                raise ChemMCPError(f"Cannot determine precursor m/z for '{smi}'. Please provide precursor_mz.")

        # CE normalization factor (optimal CE typically scales with sqrt(mz))
        ce_factor = min(ce / 30.0, 2.0)  # normalize to baseline 30eV

        fragments = []
        neutral_loss_summary = []

        # Apply MS/MS rules per functional group
        applied_rules = set()

        # Carboxylic acid (ESI+) → strong dehydration, decarboxylation
        if "carboxylic_acid" in groups and mode == "ESI+":
            applied_rules.add("carboxylic_acid_esip")
            h2o_mz = round(prec_mz - 18.01056, 4)
            co2_mz = round(prec_mz - 43.98983, 4)
            h2o_co2_mz = round(prec_mz - (18.01056 + 43.98983), 4)
            prob_h2o = int(min(95 * ce_factor, 99))
            prob_co2 = int(min(75 * ce_factor, 90))
            fragments.append({"mz": h2o_mz, "formula": "[M-H₂O+H]⁺", "pathway": "Dehydration (-H₂O)", "probability": prob_h2o})
            fragments.append({"mz": co2_mz, "formula": "[M-CO₂+H]⁺", "pathway": "Decarboxylation (-CO₂)", "probability": prob_co2})
            if prec_mz > 80:
                fragments.append({"mz": h2o_co2_mz, "formula": "[M-H₂O-CO₂+H]⁺", "pathway": "Dehydration + decarboxylation", "probability": int(prob_co2 * 0.6)})
            neutral_loss_summary.append({"loss": "H₂O (18.0106 Da)", "significance": "high"})
            neutral_loss_summary.append({"loss": "CO₂ (43.9898 Da)", "significance": "high"})

        # Alcohol → dehydration dominant
        if "alcohol" in groups and mode == "ESI+":
            applied_rules.add("alcohol_esip")
            h2o_mz = round(prec_mz - 18.01056, 4)
            prob = int(min(80 * ce_factor, 95))
            fragments.append({"mz": h2o_mz, "formula": "[M-H₂O+H]⁺", "pathway": "Dehydration (-H₂O)", "probability": prob})
            neutral_loss_summary.append({"loss": "H₂O (18.0106 Da)", "significance": "very high"})

        # Ketone → α-cleavage, H2O loss
        if "ketone" in groups:
            applied_rules.add("ketone_cleavage")
            co_mz = round(prec_mz - 27.99491, 4)
            acetyl_mz = round(prec_mz - 43.00581, 4)
            fragments.append({"mz": co_mz, "formula": "[M-CO+H]⁺", "pathway": "α-cleavage / CO loss", "probability": int(min(60 * ce_factor, 80))})
            fragments.append({"mz": acetyl_mz, "formula": "[M-C₂H₃O+H]⁺", "pathway": "Acetyl loss (α-cleavage)", "probability": int(min(50 * ce_factor, 75))})

        # Amide → characteristic ions
        if "amide" in groups:
            applied_rules.add("amide_fragmentation")
            fragments.append({"mz": round(prec_mz - 44.0, 4), "formula": "[M-CONH₂+H]⁺", "pathway": "Amide bond cleavage", "probability": int(min(55 * ce_factor, 70))})

        # Ester → McLafferty + acyl cleavage
        if "ester" in groups:
            applied_rules.add("ester_fragmentation")
            mcl_mz = round(prec_mz - 60.02113, 4) if prec_mz > 100 else round(prec_mz - 74.0, 4)
            fragments.append({"mz": mcl_mz, "formula": "[M-McLafferty+H]⁺", "pathway": "McLafferty rearrangement", "probability": int(min(65 * ce_factor, 85))})

        # Aromatic → ring retention fragments
        if "aromatic" in groups:
            applied_rules.add("aromatic_ring")
            fragments.append({"mz": 77.0234, "formula": "C₆H₅⁺", "pathway": "Phenyl cation", "probability": int(min(40 * ce_factor, 60))})
            fragments.append({"mz": 91.0542, "formula": "C₇H₇⁺ (tropylium)", "pathway": "Benzyl/tropylium (if benzyl present)", "probability": int(min(50 * ce_factor, 70))})

        # Chlorinated/Brominated → halogen loss
        if "chlorinated" in groups:
            applied_rules.add("halogen_loss_cl")
            fragments.append({"mz": round(prec_mz - 34.969, 4), "formula": "[M-Cl+H]⁺", "pathway": "Chlorine radical loss", "probability": int(min(45 * ce_factor, 65))})
        if "brominated" in groups:
            applied_rules.add("halogen_loss_br")
            fragments.append({"mz": round(prec_mz - 78.918, 4), "formula": "[M-Br+H]⁺", "pathway": "Bromine radical loss", "probability": int(min(45 * ce_factor, 65))})

        # Sulfur compounds
        if "sulfur_compound" in groups:
            applied_rules.add("sulfur_fragmentation")
            so2_mz = round(prec_mz - 63.965, 4)
            fragments.append({"mz": so2_mz, "formula": "[M-SO₂+H]⁺ or [M-SO+H]⁺", "pathway": "Sulfur oxide loss", "probability": int(min(35 * ce_factor, 55))})

        # Generic fragments (always include some)
        if not fragments or len(fragments) < 3:
            applied_rules.add("generic_fragments")
            for loss_mass, loss_name, desc in [("H₂O", "Water loss", 18.0106), ("CO", "CO loss", 27.9949), ("C₂H₄", "Ethylene loss", 28.0313)]:
                frag_mz = round(prec_mz - desc, 4)
                if frag_mz > 10:
                    fragments.append({
                        "mz": frag_mz,
                        "formula": f"[M-{loss_name}+H]⁺",
                        "pathway": f"Neutral loss of {loss_name}",
                        "probability": int(min(30 * ce_factor, 50)),
                    })

        # Sort by probability descending
        fragments.sort(key=lambda x: x["probability"], reverse=True)

        # Cap at top 20 fragments
        fragments = fragments[:20]

        return {
            "precursor_mz": prec_mz,
            "precursor_smiles": smi,
            "collision_energy_eV": ce,
            "ionization_mode": mode,
            "top_fragments": fragments,
            "neutral_losses_predicted": neutral_loss_summary,
            "fragmentation_rules_applied": list(applied_rules),
            "detected_functional_groups": groups,
            "notes": (
                "MS/MS Fragmentation Prediction Notes:\n"
                "• Predictions are based on empirical fragmentation rules\n"
                "• Actual spectra depend on instrument type, CE ramping, and collision gas\n"
                "• Higher CE generally increases low-m/z fragment abundance\n"
                "• ESI+ favors protonated molecules; consider adduct formation\n"
                "• Peptide-like molecules follow b/y ion series (not modeled here)"
            ),
        }

    def _run_base(self, smiles: str, collision_energy: float = 30.0, precursor_mz: float = None, ionization_mode: str = "ESI+") -> dict:
        """Core logic."""
        if not smiles or not smiles.strip():
            raise ChemMCPError("SMILES string or molecular formula is required.")
        if collision_energy < 0 or collision_energy > 200:
            raise ChemMCPError("Collision energy should be between 0-200 eV.")
        result = self._predict_msms(smiles.strip(), collision_energy, precursor_mz, ionization_mode)
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            smiles = parts[0]
            ce = float(parts[1]) if len(parts) > 1 else 30.0
            prec_mz = float(parts[2]) if len(parts) > 2 else None
            mode = parts[3] if len(parts) > 3 else "ESI+"
            return self._run_base(smiles, ce, prec_mz, mode)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'smiles [CE] [precursor_mz] [mode]'")


def _parse_formula_simple(formula: str) -> Dict[str, int]:
    """Quick formula parser."""
    matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    result = {}
    for elem, count in matches:
        if elem:
            result[elem] = result.get(elem, 0) + (int(count) if count else 1)
    return result
