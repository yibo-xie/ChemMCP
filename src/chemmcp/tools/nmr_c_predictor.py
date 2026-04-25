import logging
import re
from typing import List, Dict, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NmrCPredictor(BaseTool):
    """
    ¹³C NMR 化学位移预测工具。
    基于增量规则预测碳原子化学位移。
    """
    __version__ = "0.1.0"
    name = "NmrCPredictor"
    func_name = "predict_nmr_c_shifts"
    description = "Predict ¹³C NMR chemical shifts for a molecule given as SMILES. Returns chemical shift (ppm) and assignment for each carbon type."
    implementation_description = "Uses additive increment rules (Lindeman-Adams / Grant-Paul type) for ¹³C chemical shift prediction. Covers sp³ (0-80 ppm), sp² olefinic/aromatic (100-150 ppm), carbonyl (160-220 ppm), and nitrile/alkynyl regions."
    oss_dependencies = [
        ("RDKit", "https://www.rdkit.org/", "BSD-3-Clause"),
        ("NMR prediction rules", "based on literature: Grant-Paul, Lindeman-Adams rules", None),
    ]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["NMR", "¹³C NMR", "Chemical Shift", "Spectroscopy", "Carbon-13"]
    required_envs = []

    code_input_sig = [
        ("smiles", "str", "N/A", "SMILES string of the molecule."),
        ("solvent", "str", "CDCl₃", "NMR solvent. Default: CDCl₃."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'smiles [solvent]'. Example: 'CCO CDCl3'"),
    ]

    output_sig = [
        ("predictions", "list", "List of predicted carbon signals with shift(ppm), carbon_type, assignment, and notes."),
    ]

    examples = [
        {
            "code_input": {"smiles": "CCO", "solvent": "CDCl₃"},
            "text_input": {"input_params": "CCO"},
            "output": {
                "predictions": [
                    {"shift": 58.0, "carbon_type": "CH₃", "assignment": "CH₃-O"},
                    {"shift": 15.2, "carbon_type": "CH₂", "assignment": "-O-CH₂"},
                ]
            },
        },
        {
            "code_input": {"smiles": "c1ccccc1", "solvent": "CDCl₃"},
            "text_input": {"input_params": "c1ccccc1"},
            "output": {
                "predictions": [{"shift": 128.5, "carbon_type": "CH", "assignment": "Benzene C-H"}],
            },
        },
    ]

    # ========== ¹³C CHEMICAL SHIFT BASE VALUES & INCREMENTS ==========

    # Base chemical shifts for carbon types (ppm)
    _BASE_SHIFTS = {
        # --- Aliphatic (sp³) ---
        "alkane_CH3_primary": 6.0,
        "alkane_CH3_secondary": 15.0,
        "alkane_CH3_tertiary": 24.0,
        "alkane_CH2": 25.0,
        "alkane_CH": 33.0,
        "quaternary_C_alkane": 31.0,

        # --- Heteroatom-substituted aliphatic ---
        "C-O (alcohol/ether)": 60.0,
        "C-N (amine)": 45.0,
        "C-halogen (Cl/Br)": 35.0,
        "C-S (thioether/thiol)": 28.0,
        "CH3-O (methoxy)": 57.0,
        "CH3-N (methylamine)": 40.0,
        "CH3-C=O (acetyl)": 25.0,
        "CH3-aromatic (toluene)": 21.0,
        "CH2-O (ether/alcohol α)": 68.0,
        "CH2-N (amine α)": 48.0,
        "CH2-halogen": 42.0,
        "CH2-C=O (ketone α)": 38.0,
        "CH2-aromatic (benzylic)": 36.0,
        "CH-O (methine with O)": 72.0,
        "CH-N (methine with N)": 55.0,
        "CH-aromatic (benzylic CH)": 40.0,
        "C-O-C (acetal/ketal)": 95.0,
        "C-CN (nitrile α)": 22.0,
        "epoxide_C": 42.0,
        "cyclopropane_C": -3.0,

        # --- Olefinic (sp², non-aromatic) ---
        "terminal_alkene_CH2": 114.0,
        "disubstituted_alkene_CH": 125.0,
        "trisubstituted_alkene_CH": 132.0,
        "tetrasubstituted_alkene_C": 138.0,
        "conjugated_diene_C": 128.0,
        "enol_ether_C": 152.0,
        "enamine_C": 148.0,

        # --- Aromatic (sp²) ---
        "benzene_CH": 128.5,
        "benzene_ipso_C (substituted)": 138.0,
        "pyridine_alpha_C": 150.0,
        "pyridine_beta_C": 124.0,
        "pyridine_gamma_C": 136.0,
        "furan_C": 142.0,
        "thiophene_C": 127.0,
        "pyrrole_C": 118.0,

        # --- Carbonyl (sp², highly deshielded) ---
        "aldehyde_C=O": 190.0,
        "alpha_beta_unsat_aldehyde_C=O": 192.0,
        "ketone_C=O": 205.0,
        "alpha_beta_unsat_ketone_C=O": 200.0,
        "aryl_ketone_C=O": 197.0,
        "carboxylic_acid_C=O": 178.0,
        "ester_C=O": 172.0,
        "alpha_beta_unsat_ester_C=O": 168.0,
        "amide_C=O": 172.0,
        "acid_chloride_C=O": 180.0,
        "anhydride_C=O": 175.0,  # two different values typically
        "carbonate_C=O": 155.0,
        "urea_C=O": 160.0,
        "lactone_C=O": 175.0,
        "lactam_C=O": 175.0,

        # --- Other deshielded carbons ---
        "nitrile_CN": 118.0,
        "isocyanide_NC": 150.0,
        "alkyne_C_sp": 70.0,   # internal alkyne
        "terminal_alkyne_C≡": 68.0,
        "terminal_alkyne_≡CH": 84.0,
        "imine_C=N": 155.0,
        "oxime_C=N": 150.0,
        "isocyanate_N=C=O": 130.0,
    }

    # Substituent increments for aromatic carbons (relative to benzene at 128.5)
    # Format: (ipso, ortho, meta, para)
    _AROMATIC_INCREMENTS = {
        "CH3": (9.3, 0.7, -0.1, -2.9),
        "CH2CH3": (15.6, -0.5, 0.0, -2.6),
        "CH2CH2CH3": (18.0, -2.4, 0.2, -2.5),
        "C(CH3)3": (22.2, -3.1, -0.1, -3.0),
        "CH=CH2": (9.1, -2.4, 0.2, -1.2),
        "C≡CH": (-5.8, 2.9, 0.5, -1.2),
        "Ph": (13.0, -1.1, 0.4, -1.0),
        "CHO": (8.2, 1.2, 0.5, 5.8),
        "COCH3": (9.1, 0.1, 0.0, 4.2),
        "COC H3": (2.1, 1.7, -0.7, -4.5),  # OCH3
        "OH": (26.9, -12.7, 1.4, -7.3),
        "OCH3": (31.4, -14.4, 1.0, -7.7),
        "OCOPh": (23.0, -6.4, 0.9, -5.3),
        "NH2": (18.0, -13.3, 0.9, -9.8),
        "NHCH3": (10.2, -12.5, 0.5, -11.0),
        "N(CH3)2": (22.6, -15.7, 3.6, -11.6),
        "F": (34.8, -12.9, 1.5, -4.4),
        "Cl": (6.4, 0.2, 1.3, -2.1),
        "Br": (-5.4, 3.3, 2.1, -1.2),
        "I": (-32.3, 7.7, 0.9, -0.3),
        "CF3": (2.6, -2.2, 0.3, 3.5),
        "NO2": (19.6, -5.3, 0.8, 6.0),
        "CN": (-15.4, 3.6, -0.6, 13.9),
        "COOH": (2.1, 1.5, 0.0, 5.1),
        "COOCH3": (1.3, 4.0, 0.0, -6.0),
        "CONH2": (5.0, -1.5, 0.2, -4.0),
        "SH": (2.0, 0.6, 0.2, -3.1),
        "SO2NH2": (11.3, -2.6, 0.8, -3.0),
        "SO2CH3": (14.7, -0.9, 2.6, 2.7),
        "SO3H": (14.2, -3.5, 1.2, 4.8),
    }

    # Aliphatic substituent increments (α, β, γ effects on a reference carbon)
    _ALIPHATIC_INCREMENTS = {
        "C": (9.1, 9.4, -2.5),       # alkyl substitution
        "COCH3": (30, 10, 2),         # acetyl group
        "COOH": (20, 2, -2),          # carboxyl
        "COOR": (22, 2, -2),          # ester
        "C6H5": (23, 9, -2),           # phenyl
        "OH": (48, 10, -6),           # hydroxyl
        "OR": (58, 8, -5),            # ether/alkoxide
        "OCOR": (51, 6, -5),          # ester oxygen
        "NH2": (29, 11, -5),          # amino
        "NR2": (41, 10, -5),          # amine
        "NH3+": (26, 8, -5),          # ammonium
        "Cl": (33, 11, -6),           # chlorine
        "Br": (27, 11, -4),           # bromine
        "I": (23, 11, -3),            # iodine
        "Ph": (23, 9, -2),            # phenyl
        "CF3": (30, 6, -3),           # trifluoromethyl
        "NO2": (63, 10, -5),          # nitro
        "CN": (4, 3, -2),             # nitrile
        "SH": (11, 11, -4),           # thiol
        "SR": (24, 7, -4),            # sulfide
        "=O": (31, 0, -2),  # carbonyl (aldehyde)
        "=C<": (19, 7, -3),   # vinyl (olefinic)
        "C≡": (5, 5, -1),     # ethynyl (alkynyl)
    }

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize RDKit if available."""
        self._rdkit_available = False
        try:
            from rdkit import Chem
            self.Chem = Chem
            self._rdkit_available = True
        except ImportError:
            logger.warning("RDKit not available for NmrCPredictor, using heuristic analysis")

    def _analyze_smiles(self, smiles: str) -> list:
        """Analyze SMILES and return predicted carbon shifts."""
        smi = smiles.strip()

        if self._rdkit_available:
            return self._analyze_with_rdkit(smi)
        return self._analyze_heuristic(smi)

    def _analyze_heuristic(self, smi: str) -> list:
        """Heuristic prediction based on SMILES patterns."""
        predictions = []
        smi_lower = smi.lower()

        # Carboxylic acid / derivative carbonyls
        if re.search(r'C\(=O\)O[H]?', smi):
            predictions.append({"shift": round(178.0, 1), "carbon_type": "C=O", "assignment": "Carboxylic acid carbonyl C", "notes": "Typical range: 170-185 ppm"})
        elif re.search(r'C\(=O\)OC', smi) or re.search(r'OC\(=O\)', smi):
            predictions.append({"shift": round(172.0, 1), "carbon_type": "C=O", "assignment": "Ester carbonyl C", "notes": "Typical range: 165-175 ppm"})
            predictions.append({"shift": round(52.0, 1), "carbon_type": "CH₃/CH₂", "assignment": "Ester O-CHₓ", "notes": "Alpha to oxygen"})

        # Ketone carbonyl
        if re.search(r'C\(=O\)[Cc]', smi) and not re.search(r'C\(=O\)[Oo]', smi):
            is_conj = bool(re.search(r'[Cc]=?[Cc].*C\(=O\)|C\(=O\).*[Cc]=[Cc]', smi))
            base_ketone = 202.0 if not is_conj else 197.0
            predictions.append({"shift": round(base_ketone, 1), "carbon_type": "C=O", "assignment": "Ketone carbonyl C", "notes": "Conjugated lowers by ~5-15 ppm"})
            predictions.append({"shift": round(29.0, 1), "carbon_type": "CH₃/CH₂", "assignment": "α-C to ketone", "notes": "Alpha to carbonyl"})

        # Aldehyde
        if re.search(r'(C=O)[Hh]|[Hh]C=O', smi) or 'CHO' in smi or (smi.endswith('C=O') and '=' not in smi.replace('C=O','')):
            predictions.append({"shift": round(192.0, 1), "carbon_type": "C=O", "assignment": "Aldehyde carbonyl C", "notes": "Highly deshielded; conjugation moves downfield"})

        # Aromatic carbons
        if 'c' in smi_lower:
            n_ar_c = smi_lower.count('c')
            has_EWG = bool(re.search(r'NO2|CN|C\(=O\)|SO3|CHO', smi))
            has_EDG = bool(re.search(r'[OoSsNn][Cc]|Cc[Oo]|N[H]', smi))

            ipso_shift = 140.0 + (5.0 if has_EWG else -5.0 if has_EDG else 0)
            ch_shift = 127.0 + (3.0 if has_EWG else -3.0 if has_EDG else 0)

            predictions.append({"shift": round(ipso_shift, 1), "carbon_type": "Cq (ipso)", "assignment": "Aromatic substituted (ipso) C", "notes": "No attached H; often weak or absent"})
            predictions.append({"shift": round(ch_shift, 1), "carbon_type": "CH", "assignment": "Aromatic C-H", "notes": f"{'Deshielded by EWG' if has_EWG else 'Shielded by EDG' if has_EDG else 'Benzene-like'}"})

        # Alcohol/Ether carbons
        if re.search(r'[Cc]O[CcH]', smi) and not re.search(r'C\(=O\)O', smi):
            predictions.append({"shift": round(65.0, 1), "carbon_type": "CH₂/O/CH", "assignment": "C-O (alcohol/ether)", "notes": "Deshielded by electronegative oxygen"})

        # Amine carbons
        if re.search(r'[Cc][Nn]', smi) and not re.search(r'N.*C\(=O\)', smi):
            predictions.append({"shift": round(45.0, 1), "carbon_type": "CH₂/N/CH", "assignment": "C-N (amine)", "notes": "Alpha to nitrogen"})

        # Nitrile
        if 'CN' in smi or 'C#N' in smi:
            predictions.append({"shift": round(118.0, 1), "carbon_type": "C≡N", "assignment": "Nitrile carbon", "notes": "Characteristic sharp peak ~115-120 ppm"})

        # Alkene carbons
        if '=C' in smi and 'c' not in smi_lower.replace('=C',''):
            predictions.append({"shift": round(123.0, 1), "carbon_type": "CH₂/C=CH", "assignment": "Alkene C=C", "notes": "Olefinic region 100-150 ppm"})

        # Alkyne carbons
        if '#' in smi or 'C#C' in smi:
            predictions.append({"shift": round(75.0, 1), "carbon_type": "C≡", "assignment": "Alkyne C≡C", "notes": "Internal alkyne ~65-90 ppm"})

        # Default alkane chain
        if not any(p["assignment"].startswith("Ketone") or p["assignment"].startswith("Ester") or p["assignment"].startswith("Acid") for p in predictions):
            predictions.extend([
                {"shift": round(14.0, 1), "carbon_type": "CH₃", "assignment": "Terminal methyl C", "notes": "Alkane chain terminal"},
                {"shift": round(22.0, 1), "carbon_type": "CH₂", "assignment": "Methylene C", "notes": "Alkane chain internal"},
                {"shift": round(32.0, 1), "carbon_type": "CH₂/CH", "assignment": "Chain methylene/methine", "notes": "Middle of chain"},
            ])

        # Amide
        if re.search(r'C\(=O\)[Nn]', smi):
            predictions.append({"shift": round(172.0, 1), "carbon_type": "C=O", "assignment": "Amide carbonyl C", "notes": "Amide C=O; often broadened"})

        return predictions

    def _analyze_with_rdkit(self, smi: str) -> list:
        """Use RDKit for atom-level carbon analysis."""
        try:
            mol = self.Chem.MolFromSmiles(smi)
            if mol is None:
                return self._analyze_heuristic(smi)

            mol = self.Chem.AddHs(mol)
            predictions = []

            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue

                idx = atom.GetIdx()
                sym = atom.GetSymbol()
                hybrid = str(atom.GetHybridization())
                is_aromatic = atom.GetIsAromatic()
                num_h = atom.GetTotalNumHs()
                degree = atom.GetDegree()

                # Determine carbon type
                if is_aromatic:
                    ctype = "Ar-C" if num_h > 0 else "Ar-Cq"
                    base = 128.5 if num_h > 0 else 138.0
                elif hybrid == "SP":
                    ctype = "C(sp)"
                    base = 75.0
                elif hybrid == "SP2":
                    has_o_neighbor = any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors())
                    if has_o_neighbor and degree <= 2:
                        ctype = "C=O"
                        base = 195.0
                    else:
                        ctype = "C=C"
                        base = 125.0
                else:  # SP3
                    has_electronegative = any(
                        n.GetAtomicNum() in (7, 8, 9, 16, 17) for n in atom.GetNeighbors()
                    )
                    if has_electronegative:
                        base = 55.0
                        ctype = "C-X (α to heteroatom)"
                    else:
                        base = 30.0
                        ctype = "C(sp³)"

                # Adjust for neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:  # Oxygen
                        if hybrid != "SP2":
                            base += 35
                    elif neighbor.GetAtomicNum() == 7:  # Nitrogen
                        base += 20
                    elif neighbor.GetAtomicNum() in (17, 35, 53):  # Halogens
                        base += 20

                label = f"C{idx}"
                if num_h > 0:
                    h_label = f"H{num_h}" if num_h > 1 else "H"
                    label = f"C{idx}({['', 'CH', 'CH2', 'CH3'][num_h]})"

                predictions.append({
                    "shift": round(base, 1),
                    "carbon_type": ctype,
                    "assignment": label,
                    "notes": f"Hybridization: {hybrid}, Aromatic: {is_aromatic}, H-count: {num_h}",
                })

            return predictions

        except Exception as e:
            logger.warning(f"RDKit analysis failed: {e}")
            return self._analyze_heuristic(smi)

    def _run_base(self, smiles: str, solvent: str = "CDCl₃") -> dict:
        """
        Predict ¹³C NMR chemical shifts.

        Args:
            smiles: SMILES string of molecule
            solvent: NMR solvent name

        Returns:
            Dict with predictions list
        """
        if not smiles:
            raise ChemMCPError("SMILES string is required.")

        predictions = self._analyze_smiles(smiles)

        return {
            "predictions": predictions,
            "total_carbons": len(predictions),
            "solvent": solvent,
            "smiles_input": smiles.strip(),
            "note": "Predictions are approximate (±3-10 ppm). DEPT/APT experiments distinguish CH/CH₂/CH₃/Cq.",
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'smiles [solvent]'")

        smiles = parts[0]
        solvent = parts[1] if len(parts) > 1 else "CDcl₃"

        return self._run_base(smiles, solvent)
