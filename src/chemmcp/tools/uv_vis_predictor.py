"""
UV-Vis Absorption Predictor - predicts UV-Vis absorption based on chromophores
and auxochromes using Woodward-Fieser rules and empirical data.
"""

import logging
import re
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Woodward-Fieser base values for dienes (nm)
DIENE_BASE_VALUES: Dict[str, int] = {
    "homoannular": 253,   # cyclic diene in same ring (cisoid)
    "heteroannular": 214, # cyclic diene in different rings (transoid)
    "acyclic": 217,       # acyclic diene
}

# Woodward-Fieser substituent increments for dienes (nm)
DIENE_SUBSTITUENTS: Dict[str, int] = {
    "alkyl_ring_residue": 5,
    "exo_cyclic_c_double": 5,
    "double_bond_extending": 30,
    "halogen": 5,
    "auxochrome_OR": 6,
    "auxochrome_NR2": 60,
    "auxochrome_SR": 30,
}

# Woodward-Fieser base values for enones (α,β-unsaturated ketones) (nm)
ENONE_BASE_VALUES: Dict[str, int] = {
    "acyclic_enone": 215,
    "6-membered_ring_enone": 215,
    "5-membered_ring_enone": 202,
    "αβ_unsaturated_aldehyde": 207,
}

# Enone substituent increments (nm)
ENONE_SUBSTITUENTS: Dict[str, int] = {
    "alpha_alkyl": 10,
    "beta_alkyl_trans": 12,
    "beta_alkyl_cis": 0,
    "gamma_alkyl_or_longer": 5,
    "alpha_OH_OR_NHR": 35,
    "alpha_OAc": 10,
    "alpha_Cl_Br": 0,
    "beta_OH_OR": 31,
    "beta_OR": 85,  # actually beta-OR is special
    "beta_NR2": 95,
    "beta_Sr": 85,
    "gamma_OH_OR": 50,
    "delta_OH_OR": 50,
    "exocyclic_double_bond": 5,
    "homodiene_ext": 39,
    "extending_double_bond": 30,
}

# Common chromophore reference data (chromophore -> lambda_max nm, epsilon estimate)
CHROMOPHORE_DATA: List[Dict[str, Any]] = [
    {"name": "Isolated C=C", "lambda_max": 175, "epsilon": 10000, "smiles_pattern": r"C=C"},
    {"name": "Conjugated Diene", "lambda_max": 220, "epsilon": 20000, "smiles_pattern": r"C=CC=C"},
    {"name": "Triene", "lambda_max": 260, "epsilon": 35000, "smiles_pattern": r"C=CC=CC=C"},
    {"name": "Carbonyl (ketone)", "lambda_max": 280, "epsilon": 15, "smiles_pattern": r"C(=O)"},
    {"name": "Carbonyl (aldehyde)", "lambda_max": 290, "epsilon": 16, "smiles_pattern": r"C(=O)"},
    {"name": "Carboxylic acid", "lambda_max": 204, "epsilon": 41, "smiles_pattern": r"C(=O)[OH]"},
    {"name": "Ester", "lambda_max": 205, "epsilon": 50, "smiles_pattern": r"C(=O)OC"},
    {"name": "Amide", "lambda_max": 208, "epsilon": 32, "smiles_pattern": r"C(=O)N"},
    {"name": "Nitro", "lambda_max": 270, "epsilon": 14, "smiles_pattern": r"N(=O)=O"},
    {"name": "Nitroso", "lambda_max": 300, "epsilon": 100, "smiles_pattern": r"N=O"},
    {"name": "Azo (N=N)", "lambda_max": 340, "epsilon": 5, "smiles_pattern": r"N=N"},
    {"name": "Benzene ring", "lambda_max": 255, "epsilon": 215, "smiles_pattern": r"c1ccccc1"},
    {"name": "Phenol", "lambda_max": 270, "epsilon": 1450, "smiles_pattern": r"Oc1ccccc1"},
    {"name": "Aniline", "lambda_max": 230, "epsilon": 8600, "smiles_pattern": r"Nc1ccccc1"},
    {"name": "Stilbene", "lambda_max": 295, "epsilon": 25000, "smiles_pattern": r"C=CC=Cc1ccccc1"},
    {"name": "Enone (α,β)", "lambda_max": 220, "epsilon": 10000, "smiles_pattern": r"C=CC(=O)"},
    {"name": "α,β-Unsaturated aldehyde", "lambda_max": 210, "epsilon": 11500, "smiles_pattern": r"C=CC=O"},
    {"name": "Quinone", "lambda_max": 245, "epsilon": 20000, "smiles_pattern": r"O=C1C=CC(=O)C=CC1"},
]

# Auxochromes and their bathochromic shifts (nm)
AUXOCHROME_SHIFTS: Dict[str, int] = {
    "-OH, -OR": 7,
    "-NH2, -NHR, -NR2": 20,
    "-Cl, -Br": 5,
    "-SH, -SR": 10,
    "-SO3H": 15,
}

# Solvent correction factors for UV-Vis (nm)
SOLVENT_CORRECTIONS: Dict[str, int] = {
    "water": -8,
    "methanol": 0,
    "ethanol": 0,
    "hexane": 11,
    "cyclohexane": 11,
    "ether": 7,
    "chloroform": 1,
    "acetonitrile": 0,
    "dioxane": 5,
    "dichloromethane": 1,
    "dmf": 0,
    "dmso": -7,
    "thf": 0,
}


@ChemMCPManager.register_tool
class UvVisPredictor(BaseTool):
    __version__      = "0.1.0"
    name             = "UvVisPredictor"
    func_name        = "predict_uv_vis"
    description      = "Predict UV-Vis absorption maxima (λmax) based on chromophores, auxochromes, and solvent effects using Woodward-Fieser rules."
    implementation_description = "Uses Woodward-Fieser rules for conjugated dienes and α,β-unsaturated carbonyls, plus empirical chromophore database with auxochrome shifts and solvent corrections."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["UV-Vis", "Spectroscopy", "Chromophore", "Woodward-Fieser"]
    required_envs    = []

    code_input_sig   = [
        ("smiles", "str", "N/A", "SMILES string of the molecule to analyze."),
        ("solvent", "str", "methanol", "Solvent name (e.g., 'methanol', 'hexane', 'water')."),
        ("analysis_mode", "str", "auto", "Analysis mode: 'auto' (detect automatically), 'diene', 'enone', or 'general'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'SMILES [solvent] [mode]'. Example: 'C=CC=CC(C)=C methanol auto'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: detected_chromophores, predicted_absorptions (list of {lambda_max_nm, epsilon_estimate, chromophore, notes}), solvent_correction, overall_assessment."),
    ]

    examples         = [
        {
            "code_input": {"smiles": "C=CC=CC(C)=C", "solvent": "methanol", "analysis_mode": "auto"},
            "text_input": {"input_params": "C=CC=CC(C)=C methanol auto"},
            "output": {
                "result": {
                    "detected_chromophores": ["Conjugated Diene"],
                    "predicted_absorptions": [{"lambda_max_nm": 235, "epsilon_estimate": 25000, "chromophore": "Conjugated Diene (substituted)", "notes": "Alkyl-substituted acyclic diene"}],
                    "solvent_correction": 0,
                    "overall_assessment": "Expected λmax ~235 nm in methanol",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, smiles: str, solvent: str = "methanol", analysis_mode: str = "auto") -> dict:
        """Core logic: predict UV-Vis absorption."""
        if not smiles:
            raise ChemMCPError("SMILES string cannot be empty.")

        solv_corr = SOLVENT_CORRECTIONS.get(solvent.lower(), 0)

        absorptions: List[Dict[str, Any]] = []
        detected: List[str] = []

        # Detect patterns in SMILES
        s_lower = smiles

        # Count double bonds in conjugation
        conj_diene_match = re.findall(r'C=C', s_lower)
        conj_count = len(conj_diene_match)

        # Check for enone pattern
        has_enone = bool(re.search(r'C=CC(=O)|C=CC=O|C\(=O\)C=C', s_lower))
        has_carbonyl = bool(re.search(r'C(=O)', s_lower))
        has_aromatic = bool(re.search(r'c[1-9]', s_lower)) or bool(re.search(r'c1.*c1', s_lower))
        has_nitro = bool(re.search(r'N\(?=\)(=O)', s_lower)) or ('N(=O)' in s_lower or '[N+](=O)' in s_lower)
        has_phenol = bool(re.search(r'Oc[1-9]', s_lower)) or ('Oc1' in s_lower)
        has_aniline = bool(re.search(r'Nc[1-9]', s_lower)) or ('Nc1' in s_lower)
        has_ester = bool(re.search(r'C(=O)OC|OC(=O)', s_lower))
        has_amide = bool(re.search(r'C(=O)N|NC(=O)', s_lower))
        has_acid = bool(re.search(r'C(=O)[OH]|C(=O)O', s_lower))
        has_azo = bool(re.search(r'N=N', s_lower))

        if analysis_mode == "auto":
            if has_enone:
                analysis_mode = "enone"
            elif conj_count >= 2:
                analysis_mode = "diene"
            else:
                analysis_mode = "general"

        if analysis_mode == "diene":
            base_val = DIENE_BASE_VALUES["acyclic"]
            increment = 0

            # Check for homo/heteroannular
            if re.search(r'[C]1.*C=C.*C=C.*1', s_lower):
                base_val = DIENE_BASE_VALUES["homoannular"]
            elif re.search(r'[C]1.*C=C.*1.*[C]1.*C=C.*1', s_lower):
                base_val = DIENE_BASE_VALUES["heteroannular"]

            # Alkyl substituents: count non-H atoms attached roughly
            alkyl_count = len(re.findall(r'[C](?=[^0-9=()#\[\]])(?![a-z])', s_lower)) - conj_count * 2
            alkyl_count = max(0, alkyl_count)
            increment += alkyl_count * DIENE_SUBSTITUENTS["alkyl_ring_residue"]

            # Exocyclic double bond check
            if re.search(r'=C([C])', s_lower) or re.search(r'C1.*=C1', s_lower):
                increment += DIENE_SUBSTITUENTS["exo_cyclic_c_double"]

            # Extending double bonds beyond diene
            extra_db = max(0, conj_count - 2)
            increment += extra_db * DIENE_SUBSTITUENTS["double_bond_extending"]

            # Auxochrome detection
            if re.search(r'OC|O\(', s_lower):
                increment += DIENE_SUBSTITUENTS["auxochrome_OR"]
            if re.search(r'NC|N\(', s_lower):
                increment += DIENE_SUBSTITUENTS["auxochrome_NR2"]
            if re.search(r'SC|S\(', s_lower):
                increment += DIENE_SUBSTITUENTS["auxochrome_SR"]

            lam_max = base_val + increment + solv_corr
            detected.append(f"Conjugated Diene system ({conj_count} double bonds)")
            absorptions.append({
                "lambda_max_nm": lam_max,
                "epsilon_estimate": 10000 + conj_count * 8000,
                "chromophore": f"Diene (Woodward-Fieser)",
                "notes": f"Base={base_val}nm + increments={increment}nm + solvent({solvent})={solv_corr}nm",
            })

        elif analysis_mode == "enone":
            base_val = ENONE_BASE_VALUES["acyclic_enone"]
            if re.search(r'[C]1.*C=C.*C(=O).*1', s_lower):
                base_val = ENONE_BASE_VALUES["6-membered_ring_enone"]
            elif re.search(r'[C]1.*C=C.*C(=O).*1', s_lower) and re.search(r'.{0,4}C1.{0,3}', s_lower):
                pass  # keep default

            increment = 0
            # Alpha substituents
            if re.search(r'\(C\)C(=O)|C\(C\)C(=O)', s_lower):
                increment += ENONE_SUBSTITUENTS["alpha_alkyl"]
            # Beta substituents
            if re.search(r'C=C\(C\)C(=O)', s_lower):
                increment += ENONE_SUBSTITUENTS["beta_alkyl_trans"]
            # Gamma or longer
            if re.search(r'C=CC\(C\)C(=O)', s_lower):
                increment += ENONE_SUBSTITUENTS["gamma_alkyl_or_longer"]
            # Extending conjugation
            extra_db = max(0, conj_count - 2)
            increment += extra_db * ENONE_SUBSTITUENTS["extending_double_bond"]

            lam_max = base_val + increment + solv_corr
            detected.append("α,β-Unsaturated Carbonyl (Enone)")
            absorptions.append({
                "lambda_max_nm": lam_max,
                "epsilon_estimate": 12000,
                "chromophore": "Enone (Woodward-Fieser)",
                "notes": f"Base={base_val}nm + increments={increment}nm + solvent({solvent})={solv_corr}nm",
            })

        # General mode: detect all chromophores
        general_chromo = []
        if has_aromatic:
            general_chromo.append({"name": "Benzene/Aromatic ring", "base_lambda": 255})
        if has_carbonyl and not has_enone:
            if has_acid:
                general_chromo.append({"name": "Carboxylic acid", "base_lambda": 204})
            elif has_ester:
                general_chromo.append({"name": "Ester", "base_lambda": 205})
            elif has_amide:
                general_chromo.append({"name": "Amide", "base_lambda": 208})
            else:
                general_chromo.append({"name": "Carbonyl (ketone/aldehyde)", "base_lambda": 280})
        if has_nitro:
            general_chromo.append({"name": "Nitro group", "base_lambda": 270})
        if has_azo:
            general_chromo.append({"name": "Azo group", "base_lambda": 340})
        if has_phenol:
            general_chromo.append({"name": "Phenol derivative", "base_lambda": 270})
        if has_aniline:
            general_chromo.append({"name": "Aniline derivative", "base_lambda": 230})

        if analysis_mode == "general" or (analysis_mode != "diene" and analysis_mode != "enone"):
            for gc in general_chromo:
                lam = gc["base_lambda"] + solv_corr
                detected.append(gc["name"])
                absorptions.append({
                    "lambda_max_nm": lam,
                    "epsilon_estimate": 500 if "Carbonyl" in gc["name"] else (1000 if gc["name"] == "Nitro group" else 200),
                    "chromophore": gc["name"],
                    "notes": f"Empirical value, corrected by solvent ({solv_corr:+d} nm)",
                })

        if not absorptions:
            # Fallback: basic chromophore scan from database
            for ch in CHROMOPHORE_DATA:
                try:
                    pat = ch["smiles_pattern"].replace('c', '[cC]')
                    if re.search(pat, s_lower):
                        lam = ch["lambda_max"] + solv_corr
                        detected.append(ch["name"])
                        absorptions.append({
                            "lambda_max_nm": lam,
                            "epsilon_estimate": ch["epsilon"],
                            "chromophore": ch["name"],
                            "notes": f"Database match, solvent correction: {solv_corr:+d} nm",
                        })
                except Exception:
                    pass

        if not absorptions:
            raise ChemMCPError(f"No recognizable chromophores found in SMILES: '{smiles}'")

        assessment = f"Predicted {len(absorptions)} absorption band(s)"
        if absorptions:
            max_lam = max(a["lambda_max_nm"] for a in absorptions)
            min_lam = min(a["lambda_max_nm"] for a in absorptions)
            assessment += f" ranging from {min_lam}-{max_lam} nm"

        return {
            "detected_chromophores": detected,
            "predicted_absorptions": absorptions,
            "solvent_correction": solv_corr,
            "overall_assessment": assessment,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            smiles_str = parts[0]
            solvent = parts[1] if len(parts) > 1 else "methanol"
            mode = parts[2] if len(parts) > 2 else "auto"
            return self._run_base(smiles_str, solvent, mode)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'SMILES [solvent] [mode]'")
