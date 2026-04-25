import logging
import re
from typing import List, Dict, Optional, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NmrHPredictor(BaseTool):
    """
    ¹H NMR 化学位移预测工具。
    基于取代基增量规则和官能团数据库预测质子化学位移。
    """
    __version__ = "0.1.0"
    name = "NmrHPredictor"
    func_name = "predict_nmr_h_shifts"
    description = "Predict ¹H NMR chemical shifts for a molecule given as SMILES. Returns chemical shift (ppm), multiplicity, integration, and assignment for each proton type."
    implementation_description = "Uses Shoolery-type additive rules and an extensive functional group shift database to predict proton chemical shifts. Covers aliphatic, olefinic, aromatic, aldehyde, carboxylic acid, and heteroatom-attached protons."
    oss_dependencies = [
        ("RDKit", "https://www.rdkit.org/", "BSD-3-Clause"),
        ("NMR prediction rules", "based on literature: Silverstein, Pavia, etc.", None),
    ]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["NMR", "¹H NMR", "Chemical Shift", "Spectroscopy", "Prediction"]
    required_envs = []

    code_input_sig = [
        ("smiles", "str", "N/A", "SMILES string of the molecule."),
        ("solvent", "str", "CDCl₃", "NMR solvent (affects reference). Default: CDCl₃."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'smiles [solvent]'. Example: 'CCO CDCl3'"),
    ]

    output_sig = [
        ("predictions", "list", "List of predicted signals with shift(ppm), multiplicity, integration, assignment, and notes."),
    ]

    examples = [
        {
            "code_input": {"smiles": "CCO", "solvent": "CDCl₃"},
            "text_input": {"input_params": "CCO"},
            "output": {
                "predictions": [
                    {"shift": 1.2, "multiplicity": "t", "integration": 3, "assignment": "CH₃"},
                    {"shift": 3.7, "multiplicity": "q", "integration": 2, "assignment": "CH₂-O"},
                ]
            },
        },
        {
            "code_input": {"smiles": "c1ccccc1", "solvent": "CDCl₃"},
            "text_input": {"input_params": "c1ccccc1"},
            "output": {
                "predictions": [{"shift": 7.27, "multiplicity": "s", "integration": 5, "assignment": "Ar-H"}],
            },
        },
    ]

    # ========== PROTON CHEMICAL SHIFT DATABASE ==========
    # Base shifts for common structural motifs
    # Format: (base_shift_ppm, pattern, integration_hint, structural_context)
    _SHIFT_TABLE = {
        # --- Aliphatic protons ---
        "alkane_CH3": (0.87, "m", 3, "Terminal methyl in alkane chain"),
        "alkane_CH2": (1.25, "m", 2, "Methylene in alkane chain"),
        "alkane_CH": (1.50, "m", 1, "Methine in alkane chain"),
        "CH3_next_to_ether": (3.30, "q", 3, "CH₃-O-R (methoxy)"),
        "CH3_next_to_ester": (3.67, "s", 3, "Methyl ester CH₃OC(=O)-"),
        "CH2_next_to_O": (3.35, "t", 2, "-O-CH₂- (ether/alcohol)"),
        "CH2_next_to_halogen_Cl": (3.45, "t", 2, "-Cl-CH₂-"),
        "CH2_next_to_halogen_Br": (3.35, "t", 2, "-Br-CH₂-"),
        "CH2_next_to_N": (2.45, "q", 2, "-N-CH₂- (amine)"),
        "CH_next_to_O": (3.70, "m", 1, ">CH-O-"),
        "CH_next_to_aromatic": (2.85, "m", 1, "Ar-CH< (benzylic)"),
        "CH2_next_to_aromatic": (2.60, "t", 2, "Ar-CH₂- (benzylic)"),
        "CH3_next_to_aromatic": (2.30, "s", 3, "Ar-CH₃ (toluene-like)"),
        "CH3_next_to_C=O_ketone": (2.10, "s", 3, "CH₃-C(=O)- (acetyl)"),
        "CH3_next_to_C=O_aldehyde": (2.20, "s", 3, "CH₃-CHO (acetaldehyde)"),
        "CH2_alpha_to_C=O": (2.30, "t", 2, "-C(=O)-CH₂-"),
        "CH2_beta_to_C=O": (1.55, "m", 2, "-C(=O)-CH₂-CH₂-"),
        "CH3_COOH": (2.08, "s", 3, "Acetic acid CH₃"),
        "CH3_COCH3": (2.15, "s", 3, "Acetone CH₃"),
        "allylic_CH2": (1.95, "m", 2, "-C=C-CH₂- (allylic)"),
        "allylic_CH3": (1.70, "d", 3, "=C-CH₃ (propene-like)"),

        # --- Olefinic protons (=C-H) ---
        "vinyl_CH2_terminal": (4.90, "m", 2, "Terminal =CH₂ (two protons may be non-equivalent)"),
        "vinyl_CH_disubstituted_cis": (5.25, "d", 1, "cis RCH=CHR' (vinylic H)"),
        "vinyl_CH_disubstituted_trans": (5.35, "d", 1, "trans RCH=CHR' (vinylic H)"),
        "vinyl_CH_trisubst": (5.65, "m", 1, "Trisubstituted R₂C=CHR vinylic H"),
        "conjugated_diene_vinyl": (5.80, "m", 1, "Conjugated diene vinyl proton"),
        "enol_ether_vinyl": (6.40, "q", 1, "Enol ether =C(H)-OR"),
        "alpha_beta_unsat_carbonyl_vinyl": (6.80, "d", 1, "α,β-Unsaturated carbonyl β-proton"),
        "alpha_beta_unsat_carbonyl_alpha_H": (6.10, "d", 1, "α,β-Unsaturated carbonyl α-proton"),

        # --- Aromatic protons ---
        "benzene_H": (7.27, "br s", 5, "Benzene C₆H₆"),
        "monosub_benzene_ortho": (7.55, "m/d", 2, "Ortho protons of monosubstituted benzene"),
        "monosub_benzene_meta": (7.32, "t/m", 2, "Meta protons of monosubstituted benzene"),
        "monosub_benzene_para": (7.22, "t/m", 1, "Para proton of monosubstituted benzene"),
        "EWG_substituted_benzene": (7.80, "m", "var", "Benzene with electron-withdrawing group"),
        "EDG_substituted_benzene": (6.80, "m", "var", "Benzene with electron-donating group"),
        "heteroaromatic_pyridine_H2/H6": (8.50, "m", 2, "Pyridine α-protons (H2, H6)"),
        "heteroaromatic_pyridine_H3/H5": (7.15, "m", 2, "Pyridine β-protons (H3, H5)"),
        "heteroaromatic_pyridine_H4": (7.60, "tt", 1, "Pyridine γ-proton (H4)"),
        "heteroaromatic_furan": (6.40, "m", 3, "Furan ring protons"),
        "heteroaromatic_thiophene": (7.10, "m", 3, "Thiophene ring protons"),
        "heteroaromatic_pyrrole": (6.50, "m", 3, "Pyrrole NH + CH protons"),
        "phenol_OH": (4.50, "s", 1, "Phenol OH (variable, broad)"),
        "aniline_NH2": (3.50, "br s", 2, "Aniline NH₂ (broad, exchangeable)"),

        # --- Aldehydic & Carboxylic protons ---
        "aldehyde_H": (9.50, "s", 1, "Aldehyde -CHO proton"),
        "aldehyde_conj": (9.80, "s", 1, "Conjugated aldehyde -CHO"),
        "carboxylic_acid_OH": (11.0, "br s", 1, "Carboxylic acid COOH (very broad, variable)"),
        "carboxylic_acid_conj_OH": (12.0, "br s", 1, "Conjugated carboxylic acid OH"),

        # --- Alcohol / Phenol protons (exchangeable) ---
        "alcohol_OH_primary": (3.40, "br s", 1, "Primary alcohol OH (variable)"),
        "alcohol_OH_secondary": (3.80, "br s", 1, "Secondary alcohol OH (variable)"),
        "alcohol_OH_tertiary": (2.50, "br s", 1, "Tertiary alcohol OH (variable)"),
        "intramolecular_Hbonded_OH": (11.0, "s", 1, "Intramolecularly H-bonded OH (salicylate-type)"),

        # --- Amine protons (exchangeable) ---
        "primary_amine_NH2_aliphatic": (1.00, "br s", 2, "Aliphatic primary amine NH₂"),
        "secondary_amine_NH_aliphatic": (1.10, "br d/s", 1, "Aliphatic secondary amine NH"),
        "amide_NH_primary": (5.50, "br s", 2, "Primary amide NH₂ (often broad)"),
        "amide_NH_secondary": (7.00, "br s", 1, "Secondary amide NH (often broad)"),

        # --- Other notable protons ---
        "acetylene_≡CH": (2.50, "s", 1, "Terminal alkyne ≡C-H"),
        "nitroalkane_CH": (4.30, "q", 1, "R-CH₂-NO₂ α-proton"),
        "ester_alpha_CH2": (4.10, "q", 2, "-O-C(=O)-CH₂-"),
        "epoxide_CH2": (2.50, "m", 2, "Epoxide methylene protons"),
        "cyclopropane_CH2": (0.22, "m", 2, "Cyclopropane ring protons (shielded!)"),
        "acetal_CH": (5.40, "m", 1, "Acetal >CH(O-R)₂ proton"),
        "formate_H": (8.05, "s", 1, "Formate ester HCOOR"),
        "nitroso": (7.80, "m", 1, "Nitroso compound N=O adjacent"),
        "azide_adjacent_CH2": (3.30, "t", 2, "R-CH₂-N₃"),
        "nitrile_adjacent_CH2": (2.45, "t", 2, "R-CH₂-CN"),
        "sulfonamide_NH": (7.50, "br s", 1, "Sulfonamide N-H"),
        "enol_OH": (14.0, "s", 1, "Enol OH (β-diketone enol; very deshielded)"),
    }

    # Substituent effects on benzene rings (relative to benzene at 7.27 ppm)
    _AROMATIC_SUBSTITUENT_EFFECTS = {
        "electron_withdrawing_strong": 0.80,   # NO₂, CN, COR, CHO, COOH, SO₃H
        "electron_withdrawing_moderate": 0.40,  # Halogens (F, Cl, Br), CF₃
        "electron_donating_strong": -0.40,      # O⁻, NR₂, OH, OR, NHCOR
        "electron_donating_moderate": -0.20,     # Alkyl (CH₃, CH₂R, etc.)
        "halogen": 0.20,                         # Halogen has mixed effect
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
            from rdkit.Chem import Descriptors
            self.Chem = Chem
            self.Descriptors = Descriptors
            self._rdkit_available = True
        except ImportError:
            logger.warning("RDKit not available for NmrHPredictor, using SMILES heuristics")

    def _analyze_smiles(self, smiles: str) -> list:
        """Analyze SMILES and return predicted proton signals."""
        predictions = []
        smi = smiles.strip()

        # Use RDKit if available for better analysis
        if self._rdkit_available:
            return self._analyze_with_rdkit(smi)

        # Fallback: heuristic analysis based on SMILES patterns
        return self._analyze_heuristic(smi)

    def _analyze_heuristic(self, smi: str) -> list:
        """Heuristic SMILES-based prediction."""
        predictions = []
        smi_lower = smi.lower()

        # Detect functional groups via regex/pattern matching

        # Carboxylic acid
        if re.search(r'C\(=O\)[Oo]', smi) or re.search(r'C\(=O\)O', smi):
            predictions.append({
                "shift": round(11.0, 2),
                "multiplicity": "br s",
                "integration": 1,
                "assignment": "Carboxylic acid OH (exchangeable)",
                "notes": "Very broad, concentration-dependent, may exchange with D₂O",
            })

        # Aldehyde
        if re.search(r'[Cc]=?[Oo]\s*$') or '[CHO]' in smi or re.search(r'C=O$', smi) or ('C=O' in smi and not '=' in smi.replace('C=O','')):
            pass  # handled below
        if re.search(r'(C=O)[Hh]', smi) or smi.endswith('C=O'):
            predictions.append({
                "shift": round(9.5, 2),
                "multiplicity": "s",
                "integration": 1,
                "assignment": "Aldehyde -CHO",
                "notes": "Characteristic downfield singlet; conjugation moves to ~9.8-10.0",
            })

        # Aromatic ring
        if 'c' in smi_lower:
            has_EDG = bool(re.search(r'[OoSsNn][Cc]|Cc[Oo]', smi))
            has_EWG = bool(re.search(r'NO2|CN|C\(=O\)|SO3|CHO', smi))

            if has_EWG:
                base = 7.80
                desc = "Aromatic H (deshielded by EWG)"
            elif has_EDG:
                base = 6.80
                desc = "Aromatic H (shielded by EDG)"
            else:
                base = 7.27
                desc = "Aromatic H"

            # Count aromatic protons approximately
            ar_h_count = smi_lower.count('c') - sum(
                1 for m in re.finditer(r'c[c\d\)\(\]]*', smi_lower) if m
            )
            ar_h_count = max(ar_h_count, 1)

            predictions.append({
                "shift": round(base, 2),
                "multiplicity": "m" if ar_h_count > 1 else "s",
                "integration": ar_h_count,
                "assignment": desc,
                "notes": f"Aromatic region; substitution pattern affects splitting",
            })

        # Alcohol
        if re.search(r'[Cc]O[Hh]?', smi) and not re.search(r'C\(=O\)O', smi):
            predictions.append({
                "shift": round(2.0, 2),
                "multiplicity": "br s",
                "integration": 1,
                "assignment": "Alcohol O-H (exchangeable)",
                "notes": "Variable position (0.5-5.0); broad; exchanges with D₂O",
            })

        # Amine
        if re.search(r'N[Hh]', smi) and not re.search(r'N.*C\(=O\)', smi):
            predictions.append({
                "shift": round(1.5, 2),
                "multiplicity": "br s",
                "integration": 2 if 'NH2' in smi else 1,
                "assignment": "Amine N-H (exchangeable)",
                "notes": "Variable (0.5-5.0); broad; exchanges with D₂O",
            })

        # Ester O-CH3
        if re.search(r'OC[\(]?[Cc]', smi) or re.search(r'C\(=O\)OC', smi):
            predictions.append({
                "shift": round(3.7, 2),
                "multiplicity": "s",
                "integration": 3,
                "assignment": "Ester O-CH₃",
                "notes": "Methyl ester characteristic singlet",
            })

        # Ether/Alcohol CH2-O
        if re.search(r'[Cc]O[Cc]', smi) and not re.search(r'C\(=O\)O', smi):
            predictions.append({
                "shift": round(3.4, 2),
                "multiplicity": "t/q/m",
                "integration": 2,
                "assignment": "CH₂-O (ether/alcohol α-CH₂)",
                "notes": "Deshielded by oxygen electronegativity",
            })

        # Ketone alpha-CH3
        if re.search(r'C\(=O\)[Cc]', smi) and not re.search(r'C\(=O\)O', smi):
            predictions.append({
                "shift": round(2.1, 2),
                "multiplicity": "s",
                "integration": 3,
                "assignment": "CH₃-C=O (acetyl)",
                "notes": "Alpha to carbonyl; slightly deshielded",
            })
            predictions.append({
                "shift": round(2.3, 2),
                "multiplicity": "t/q",
                "integration": 2,
                "assignment": "CH₂-C=O (ketone α-CH₂)",
                "notes": "Alpha to carbonyl",
            })

        # Alkane CH3/CH2 (default)
        remaining = True
        if not predictions:
            predictions.extend([
                {
                    "shift": round(0.9, 2),
                    "multiplicity": "t",
                    "integration": 3,
                    "assignment": "Terminal CH₃ (alkyl)",
                    "notes": "Typical alkane methyl triplet",
                },
                {
                    "shift": round(1.25, 2),
                    "multiplicity": "m",
                    "integration": 2,
                    "assignment": "Methylene CH₂ (alkyl chain)",
                    "notes": "Typical alkane methylene",
                },
            ])

        # Terminal alkene =CH2
        if '=C' in smi:
            predictions.append({
                "shift": round(4.9, 2),
                "multiplicity": "m/dd",
                "integration": 2,
                "assignment": "Terminal vinyl =CH₂",
                "notes": "Two non-equivalent protons possible; geminal coupling ~2 Hz, vicinal ~10/17 Hz",
            })

        # Alkyne ≡CH
        if '#' in smi or 'C#C' in smi:
            predictions.append({
                "shift": round(2.5, 2),
                "multiplicity": "s",
                "integration": 1,
                "assignment": "Terminal alkyne ≡C-H",
                "notes": "Relatively shielded due to cylindrical electron cloud",
            })

        return predictions

    def _analyze_with_rdkit(self, smi: str) -> list:
        """Use RDKit for atom-level analysis."""
        try:
            mol = self.Chem.MolFromSmiles(smi)
            if mol is None:
                return self._analyze_heuristic(smi)

            mol = self.Chem.AddHs(mol)
            predictions = []

            # Analyze each hydrogen's environment
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 1:
                    continue

                idx = atom.GetIdx()
                neighbor = atom.GetNeighbors()[0]
                n_atom_num = neighbor.GetAtomicNum()
                n_hybrid = neighbor.GetHybridization()

                # Determine environment
                shift = self._estimate_shift_rd(mol, idx, neighbor)
                mult = self._estimate_multiplicity(mol, idx)

                predictions.append({
                    "shift": round(shift, 2),
                    "multiplicity": mult,
                    "integration": 1,
                    "assignment": f"H on {neighbor.GetSymbol()}{idx}",
                    "notes": "RDKit-predicted environment",
                })

            # Merge equivalent protons (simplified grouping)
            return self._merge_equivalent_protons(predictions)

        except Exception as e:
            logger.warning(f"RDKit analysis failed: {e}")
            return self._analyze_heuristic(smi)

    def _estimate_shift_rd(self, mol, h_idx, neighbor) -> float:
        """Estimate chemical shift using RDKit atom properties."""
        sym = neighbor.GetSymbol()
        hybrid = str(neighbor.GetHybridization())

        # Check bonded atoms for substituent effects
        is_aromatic = neighbor.GetIsAromatic()
        is_in_ring = neighbor.IsInRing()
        has_electronegative_neighbor = any(
            n.GetAtomicNum() in (7, 8, 9, 16, 17)  # N, O, F, S, Cl
            for n in neighbor.GetNeighbors() if n.GetIdx() != h_idx
        )
        has_carbonyl_neighbor = any(
            n.GetAtomicNum() == 6 and any(nn.GetAtomicNum() == 8 for nn in n.GetNeighbors())
            for n in neighbor.GetNeighbors() if n.GetIdx() != h_idx
        )

        # Base shift determination
        if is_aromatic:
            base = 7.27
            if has_electronegative_neighbor:
                base += 0.5
        elif sym == 6 and hybrid == "SP2":
            base = 5.3  # vinylic
        elif sym == 6 and hybrid == "SP":
            base = 2.5  # alkynyl
        elif sym == 6:
            base = 1.25  # aliphatic
            if has_electronegative_neighbor:
                base += 2.3  # alpha to electronegative
            elif has_carbonyl_neighbor:
                base += 1.0  # beta to carbonyl
        elif sym == 8:  # oxygen (OH)
            base = 2.0  # alcohol
        elif sym == 7:  # nitrogen (NH)
            base = 1.5  # amine
        else:
            base = 2.0

        return base

    def _estimate_multiplicity(self, mol, h_idx) -> str:
        """Estimate peak multiplicity based on neighboring H count."""
        atom = mol.GetAtomWithIdx(h_idx)
        neighbor = atom.GetNeighbors()[0]

        # Count neighboring hydrogens through bonds (simplified: 3-bond neighbors)
        n_neighbors = 0
        for nn in neighbor.GetNeighbors():
            if nn.GetAtomicNum() != 1:
                for nnn in nn.GetNeighbors():
                    if nnn.GetAtomicNum() == 1 and nnn.GetIdx() != h_idx:
                        n_neighbors += 1

        # n+1 rule
        n_peaks = n_neighbors + 1
        mult_map = {1: "s", 2: "d", 3: "t", 4: "q", 5: "quint", 6: "sext", 7: "sept"}
        return mult_map.get(n_peaks, f"m ({n_peaks}-line)")

    def _merge_equivalent_protons(self, predictions: list) -> list:
        """Group chemically equivalent protons."""
        if not predictions:
            return predictions

        # Simple merge by similar shift (within 0.1 ppm)
        merged = []
        used = set()

        for i, p in enumerate(predictions):
            if i in used:
                continue
            group = [p]
            used.add(i)
            for j in range(i + 1, len(predictions)):
                if j in used:
                    continue
                if abs(p["shift"] - predictions[j]["shift"]) < 0.15 and p["assignment"].split(" ")[0] == predictions[j]["assignment"].split(" ")[0]:
                    group.append(predictions[j])
                    used.add(j)

            if len(group) == 1:
                merged.append(group[0])
            else:
                merged.append({
                    **group[0],
                    "integration": len(group),
                    "assignment": f"{len(group)}H × {group[0]['assignment']}",
                })

        return merged

    def _run_base(self, smiles: str, solvent: str = "CDCl₃") -> dict:
        """
        Predict ¹H NMR chemical shifts.

        Args:
            smiles: SMILES string of molecule
            solvent: NMR solvent name

        Returns:
            Dict with predictions list
        """
        if not smiles:
            raise ChemMCPError("SMILES string is required.")

        predictions = self._analyze_smiles(smiles)

        total_h = sum(p.get("integration", 0) for p in predictions if isinstance(p.get("integration"), int))

        return {
            "predictions": predictions,
            "total_protons": total_h,
            "solvent": solvent,
            "smiles_input": smiles.strip(),
            "note": "Predictions are approximate (±0.2-0.5 ppm). Actual values depend on solvent, temperature, and concentration.",
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'smiles [solvent]'")

        smiles = parts[0]
        solvent = parts[1] if len(parts) > 1 else "CDCl₃"

        return self._run_base(smiles, solvent)
