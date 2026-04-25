import logging
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# NMR chemical shift reference data (¹H and ¹³C)
# Shifts in ppm relative to TMS = 0 ppm
# Based on standard NMR reference tables (Silverstein, Pretsch, etc.)

# ¹H NMR chemical shift ranges (ppm)
H1_SHIFT_DATA = {
    # Aliphatic protons
    "H-C(sp³) (alkane)": {"range": (0.7, 1.3), "typical": 0.9, "notes": "Terminal CH₃"},
    "H-C(sp³) (CH₂ in chain)": {"range": (1.2, 1.4), "typical": 1.3, "notes": "Methylene"},
    "H-C(sp³) (allylic)": {"range": (1.6, 2.1), "typical": 1.8, "notes": "CH₂ adjacent to C=C"},
    "H-C(sp³) (α to carbonyl)": {"range": (2.0, 2.6), "typical": 2.3, "notes": "α to C=O"},
    "H-C(sp³) (α to aromatic)": {"range": (2.3, 2.8), "typical": 2.5, "notes": "Benzyl position"},
    "H-C(sp³) (α to ester/acid)": {"range": (2.0, 2.4), "typical": 2.2, "notes": ""},
    "H-C(sp³) (α to nitrile)": {"range": (2.1, 2.5), "typical": 2.3, "notes": "-CH₂CN"},
    "H-C(sp³) (α to halogen Cl/Br)": {"range": (3.3, 4.0), "typical": 3.5, "notes": "Electronegative substituent effect"},
    "H-C(sp³) (α to oxygen, ether)": {"range": (3.2, 3.9), "typical": 3.5, "notes": "-O-CH₂-"},
    "H-C(sp³) (α to oxygen, alcohol)": {"range": (3.3, 4.0), "typical": 3.6, "notes": "-OH / -O-CH"},
    "H-C(sp³) (α to nitrogen)": {"range": (2.5, 3.5), "typical": 3.0, "notes": "Amine α-protons"},
    # Olefinic protons
    "H-C(sp²) (alkene, monosubst.)": {"range": (4.9, 5.4), "typical": 5.1, "notes": "=CH₂ terminal"},
    "H-C(sp²) (alkene, disubst. trans)": {"range": (5.2, 5.7), "typical": 5.4, "notes": "trans alkene"},
    "H-C(sp²) (alkene, disubst. cis)": {"range": (5.5, 6.3), "typical": 5.9, "notes": "cis alkene"},
    "H-C(sp²) (vinyl, trisubst.)": {"range": (5.3, 5.9), "typical": 5.6, "notes": ""},
    "H-C(sp²) (aromatic)": {"range": (6.5, 8.5), "typical": 7.27, "notes": "Benzene: 7.27 ppm"},
    "H-C(sp²) (aromatic, electron-withdrawing)": {"range": (7.5, 9.0), "typical": 8.2, "notes": "Nitrobenzene etc."},
    "H-C(sp²) (heteroaromatic, furan)": {"range": (6.5, 7.8), "typical": 7.3, "notes": ""},
    "H-C(sp²) (heteroaromatic, pyridine)": {"range": (7.2, 8.6), "typical": 8.5, "notes": "α-H ~8.5 ppm"},
    "H-C(sp²) (heteroaromatic, pyrrole)": {"range": (6.3, 7.0), "typical": 6.5, "notes": "N-H proton ~7-8 ppm"},
    # Aldehyde
    "H-C(=O) (aldehyde)": {"range": (9.0, 10.2), "typical": 9.8, "notes": "Very deshielded"},
    # Carboxylic acid
    "COOH (carboxylic acid)": {"range": (10.0, 13.0), "typical": 11.5, "notes": "Broad, variable"},
    # Hydroxyl
    "O-H (alcohol)": {"range": (0.5, 5.5), "typical": 2.0, "notes": "Variable; concentration-dependent"},
    "O-H (phenol)": {"range": (4.0, 12.0), "typical": 7.0, "notes": "Broad, H-bond dependent"},
    # Amine/amide
    "N-H (amine)": {"range": (0.5, 5.5), "typical": 2.0, "notes": "Broad; exchangeable"},
    "N-H (amide)": {"range": (5.5, 9.0), "typical": 7.5, "notes": "Broad; often coupled"},
    # Alkyne
    "H-C(sp) (alkyne)": {"range": (1.8, 3.1), "typical": 2.5, "notes": "Shielded by anisotropy"},
}

# ¹³C NMR chemical shift ranges (ppm)
C13_SHIFT_DATA = {
    "C (alkane, CH₃)": {"range": (0, 35), "typical": 10, "notes": "Primary carbon"},
    "C (alkane, CH₂)": {"range": (15, 55), "typical": 30, "notes": "Secondary carbon"},
    "C (alkane, CH)": {"range": (25, 60), "typical": 35, "notes": "Tertiary carbon"},
    "C (alkane, quaternary)": {"range": (25, 50), "typical": 35, "notes": "Quaternary sp³ C"},
    "C (allylic/C=C adjacent)": {"range": (20, 45), "typical": 30, "notes": ""},
    "C (α to carbonyl)": {"range": (30, 65), "typical": 45, "notes": "α-C of ketone/aldehyde"},
    "C (α to oxygen, ether/alcohol)": {"range": (50, 90), "typical": 65, "notes": "-O-CH₂-, -O-CH<"},
    "C (α to nitrogen)": {"range": (35, 70), "typical": 50, "notes": ""},
    "C (alkene, =CH₂)": {"range": (100, 125), "typical": 115, "notes": "Terminal alkene C"},
    "C (alkene, =CH-)": {"range": (120, 145), "typical": 135, "notes": "Internal alkene C"},
    "C (alkene, =C<)": {"range": (130, 155), "typical": 145, "notes": "Trisubstituted alkene C"},
    "C (aromatic)": {"range": (115, 160), "typical": 128, "notes": "Benzene: 128.5 ppm"},
    "C (aromatic, ipso to EWG)": {"range": (130, 165), "typical": 148, "notes": "Electron-withdrawing group effect"},
    "C (aromatic, ipso to EDG)": {"range": (120, 150), "typical": 138, "notes": "Electron-donating group effect"},
    "C≡C (alkyne)": {"range": (65, 95), "typical": 80, "notes": "Shielded by anisotropy"},
    "C≡N (nitrile)": {"range": (112, 126), "typical": 119, "notes": ""},
    "C=O (ketone)": {"range": (195, 220), "typical": 210, "notes": ""},
    "C=O (aldehyde)": {"range": (190, 205), "typical": 200, "notes": ""},
    "C=O (carboxylic acid)": {"range": (170, 185), "typical": 178, "notes": ""},
    "C=O (ester)": {"range": (160, 180), "typical": 172, "notes": ""},
    "C=O (amide)": {"range": (160, 180), "typical": 172, "notes": ""},
    "C=O (anhydride)": {"range": (165, 180), "typical": 170, "notes": "Two signals"},
    "C=O (acid chloride)": {"range": (175, 195), "typical": 185, "notes": "Desheilded"},
}


@ChemMCPManager.register_tool
class NmrChemicalShift(BaseTool):
    """
    NMR化学位移预测工具。
    预测分子中不同类型氢原子和碳原子的NMR化学位移范围，基于取代基效应和经验数据库。
    """
    __version__ = "0.1.0"
    name = "NmrChemicalShift"
    func_name = "predict_nmr_chemical_shift"
    description = "Predict NMR (¹H and ¹³C) chemical shifts for various atomic environments based on empirical databases."
    implementation_description = "Uses comprehensive empirical chemical shift correlation tables for ¹H and ¹³C NMR spectroscopy, covering aliphatic, olefinic, aromatic, heteroaromatic, and functionalized environments with substituent effects."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "NMR", "Chemical Shift", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("nucleus", "str", "1H", "Nucleus type: '1H' for proton or '13C' for carbon-13."),
        ("environments", "list", "N/A", "List of atomic environment descriptions (e.g., ['methyl_alkane', 'aromatic', 'aldehyde'])."),
        ("smiles", "str", "None", "Optional SMILES for additional context."),
        ("solvent", "str", "CDCl3", "NMR solvent (affects shifts slightly)."),
        ("include_coupling_info", "bool", "False", "Whether to include typical coupling constant ranges."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: nucleus env1 env2 ... e.g., '1H methyl_aromatic aldehyde'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing predicted chemical shifts for each environment, including range, typical value, multiplicity, and notes."),
    ]

    examples = [
        {
            "code_input": {
                "nucleus": "1H",
                "environments": ["methyl_alkane", "aromatic", "aldehyde"],
                "smiles": None,
                "solvent": "CDCl3",
                "include_coupling_info": False,
            },
            "text_input": {
                "input_params": "1H methyl_alkane aromatic aldehyde",
            },
            "output": {
                "result": {
                    "nucleus": "¹H",
                    "solvent": "CDCl₃",
                    "predictions": [
                        {"environment": "H-C(sp³) (alkane)", "shift_range_ppm": "0.7-1.3", "typical_ppm": 0.9, "multiplicity": "t (if CH₃)", "integration": "3H", "notes": "Terminal CH₃"},
                        {"environment": "H-C(sp²) (aromatic)", "shift_range_ppm": "6.5-8.5", "typical_ppm": 7.27, "multiplicity": "m", "integration": "~5H", "notes": "Benzene ring"},
                        {"environment": "H-C(=O) (aldehyde)", "shift_range_ppm": "9.0-10.2", "typical_ppm": 9.8, "multiplicity": "s", "integration": "1H", "notes": "Very deshielded"},
                    ],
                    "reference": "TMS = 0.0 ppm",
                }
            }
        },
        {
            "code_input": {
                "nucleus": "13C",
                "environments": ["ketone_carbonyl", "aromatic", "alcohol_alpha"],
                "smiles": None,
                "solvent": "CDCl3",
                "include_coupling_info": False,
            },
            "text_input": {
                "input_params": "13C ketone_carbonyl aromatic alcohol_alpha",
            },
            "output": {
                "result": {
                    "nucleus": "¹³C",
                    "solvent": "CDCl₃",
                    "predictions": [
                        {"environment": "C=O (ketone)", "shift_range_ppm": "195-220", "typical_ppm": 210, "notes": "Carbonyl carbon"},
                        {"environment": "C (aromatic)", "shift_range_ppm": "115-160", "typical_ppm": 128, "notes": "Aromatic ring carbons"},
                        {"environment": "C (α to oxygen, ether/alcohol)", "shift_range_ppm": "50-90", "typical_ppm": 65, "notes": "Carbon bonded to oxygen"},
                    ],
                    "reference": "TMS = 0.0 ppm",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h1_db = dict(H1_SHIFT_DATA)
        self.c13_db = dict(C13_SHIFT_DATA)

    def _run_base(self, nucleus: str, environments: List[str], smiles: str = None,
                  solvent: str = "CDCl3", include_coupling_info: bool = False) -> dict:
        """Core logic."""
        nuc = nucleus.upper().replace(" ", "")
        if nuc in ("1H", "H1", "PROTON"):
            db = self.h1_db
            nuc_label = "¹H"
        elif nuc in ("13C", "C13", "CARBON"):
            db = self.c13_db
            nuc_label = "¹³C"
        else:
            raise ChemMCPError(f"Unsupported nucleus '{nucleus}'. Use '1H' or '13C'.")

        if not environments:
            raise ChemMCPError("At least one environment must be provided.")

        predictions = []
        for env in environments:
            matched = self._match_env(env.lower(), db)
            if matched is None:
                predictions.append({
                    "environment": env,
                    "status": "not_found",
                    "note": f"No match found for '{env}' in {nuc_label} database.",
                })
                continue

            key, data = matched
            lo, hi = data["range"]
            pred = {
                "environment": key,
                "shift_range_ppm": f"{lo}-{hi}",
                "typical_ppm": data["typical"],
                "midpoint_ppm": round((lo + hi) / 2, 1),
                "width_ppm": round(hi - lo, 1),
                "notes": data.get("notes", ""),
            }

            if include_coupling_info and nuc == "1H":
                pred.update(self._coupling_info(key))

            predictions.append(pred)

        result = {
            "nucleus": nuc_label,
            "solvent": solvent,
            "predictions": predictions,
            "reference": "TMS = 0.0 ppm",
        }

        if include_coupling_info:
            result["coupling_note"] = "Coupling constants are approximate; actual values depend on molecular geometry and solvent."

        return {"result": result}

    def _match_env(self, env: str, db: dict):
        """Match user input to database key."""
        direct_map_h1 = {
            "methyl_alkane": "H-C(sp³) (alkane)",
            "methylene": "H-C(sp³) (CH₂ in chain)",
            "allylic": "H-C(sp³) (allylic)",
            "alpha_carbonyl": "H-C(sp³) (α to carbonyl)",
            "benzyl": "H-C(sp³) (α to aromatic)",
            "ether": "H-C(sp³) (α to oxygen, ether)",
            "alcohol": "H-C(sp³) (α to oxygen, alcohol)",
            "amine_alpha": "H-C(sp³) (α to nitrogen)",
            "halide_alpha": "H-C(sp³) (α to halogen Cl/Br)",
            "alkene_terminal": "H-C(sp²) (alkene, monosubst.)",
            "alkene_trans": "H-C(sp²) (alkene, disubst. trans)",
            "alkene_cis": "H-C(sp²) (alkene, disubst. cis)",
            "aromatic": "H-C(sp²) (aromatic)",
            "aldehyde": "H-C(=O) (aldehyde)",
            "carboxylic_acid": "COOH (carboxylic acid)",
            "oh_alcohol": "O-H (alcohol)",
            "oh_phenol": "O-H (phenol)",
            "nh_amine": "N-H (amine)",
            "nh_amide": "N-H (amide)",
            "alkyne": "H-C(sp) (alkyne)",
        }
        direct_map_c13 = {
            "methyl_alkane": "C (alkane, CH₃)",
            "methylene": "C (alkane, CH₂)",
            "methine": "C (alkane, CH)",
            "quaternary_sp3": "C (alkane, quaternary)",
            "ether": "C (α to oxygen, ether/alcohol)",
            "alcohol": "C (α to oxygen, ether/alcohol)",
            "alkene": "C (alkene, =CH-)",
            "aromatic": "C (aromatic)",
            "alkyne": "C≡C (alkyne)",
            "nitrile": "C≡N (nitrile)",
            "ketone": "C=O (ketone)",
            "aldehyde": "C=O (aldehyde)",
            "carboxylic_acid": "C=O (carboxylic acid)",
            "ester": "C=O (ester)",
            "amide": "C=O (amide)",
        }

        dmap = direct_map_h1 if db is self.h1_db else direct_map_c13
        if env in dmap:
            key = dmap[env]
            if key in db:
                return key, db[key]

        # Fuzzy match
        for key in db:
            clean_key = key.lower().replace(" ", "_").replace("-", "").replace("(", "").replace(")", "").replace(",", "")
            clean_env = env.replace(" ", "_").replace("-", "").replace("(", "").replace(")", "")
            if clean_env in clean_key or clean_key in clean_env or any(w in clean_key for w in clean_env.split("_")):
                return key, db[key]

        return None

    @staticmethod
    def _coupling_info(env_key: str) -> dict:
        """Add typical coupling information for ¹H."""
        coupling_data = {
            "H-C(sp³) (alkane)": {"multiplicity": "t (triplet) if CH₃", "J_typical_Hz": "7", "integration": "3H"},
            "H-C(sp³) (CH₂ in chain)": {"multiplicity": "q (quartet) if next to CH₃", "J_typical_Hz": "7", "integration": "2H"},
            "H-C(sp²) (aromatic)": {"multiplicity": "m (multiplet)", "J_typical_Hz": "7-9 (ortho)", "integration": "variable"},
            "H-C(=O) (aldehyde)": {"multiplicity": "s (singlet) or small J", "J_typical_Hz": "1-3 (allylic)", "integration": "1H"},
            "H-C(sp²) (alkene, disubst. trans)": {"multiplicity": "d (doublet) or q", "J_typical_Hz": "12-18 (trans)", "integration": "1H"},
            "H-C(sp²) (alkene, disubst. cis)": {"multiplicity": "d (doublet) or q", "J_typical_Hz": "6-12 (cis)", "integration": "1H"},
            "O-H (alcohol)": {"multiplicity": "s (broad singlet)", "J_typical_Hz": "exchange-broadened", "integration": "1H (exchangeable)"},
            "N-H (amine)": {"multiplicity": "s (broad singlet)", "J_typical_Hz": "exchange-broadened", "integration": "1-2H (exchangeable)"},
        }
        return coupling_data.get(env_key, {"multiplicity": "variable", "J_typical_Hz": "depends on structure"})

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            nucleus = parts[0]
            envs = parts[1:]
            return self._run_base(nucleus, envs)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
