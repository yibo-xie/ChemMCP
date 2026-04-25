import logging
import re
from typing import List, Dict, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MassSpecFragmenter(BaseTool):
    """
    质谱碎片化模式预测工具。
    预测分子的主要碎片离子（m/z 值）及碎裂途径。
    """
    __version__ = "0.1.0"
    name = "MassSpecFragmenter"
    func_name = "predict_mass_spec_fragments"
    description = "Predict mass spectrometry fragmentation patterns for a molecule given as SMILES or molecular formula. Returns molecular ion, base peak candidates, and fragment ions with m/z values and fragmentation pathways."
    implementation_description = "Uses common EI-MS fragmentation rules (McLafferty rearrangement, alpha-cleavage, inductive cleavage, loss of neutral molecules) to predict major fragment ions. Covers common neutral losses and characteristic fragmentation patterns for alcohols, ketones, esters, aromatics, amines, halogenated compounds, etc."
    oss_dependencies = [
        ("RDKit", "https://www.rdkit.org/", "BSD-3-Clause"),
        ("Fragmentation rules", "based on McLafferty, Tureček, and standard organic MS textbooks", None),
    ]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Mass Spectrometry", "Fragmentation", "m/z", "EI-MS", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("smiles", "str", "N/A", "SMILES string of the molecule (or molecular formula)."),
        ("ionization_mode", "str", "EI+", "Ionization mode: 'EI+' (electron impact), 'CI+' (chemical ionization), 'ESI+' (electrospray), 'ESI-' (negative ESI)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'smiles [ionization_mode]'. Example: 'CCO EI+'"),
    ]

    output_sig = [
        ("fragmentation_data", "dict", "Complete fragmentation analysis including molecular ion, fragments with m/z, formulas, and pathway descriptions."),
    ]

    examples = [
        {
            "code_input": {"smiles": "CCO", "ionization_mode": "EI+"},
            "text_input": {"input_params": "CCO"},
            "output": {
                "fragmentation_data": {
                    "molecular_ion_mz": 46,
                    "base_peak_mz": 31,
                    "fragments": [{"mz": 31, "formula": "CH₂O⁺", "pathway": "α-cleavage next to oxygen"}],
                }
            },
        },
        {
            "code_input": {"smiles": "c1ccccc1", "ionization_mode": "EI+"},
            "text_input": {"input_params": "c1ccccc1"},
            "output": {
                "fragmentation_data": {
                    "molecular_ion_mz": 78,
                    "base_peak_mz": 77,
                    "fragments": [{"mz": 77, "formula": "C₆H₅⁺", "pathway": "Loss of H•"}],
                }
            },
        },
    ]

    # ========== ATOMIC WEIGHTS (most abundant isotope) ==========
    _ATOMIC_MASSES = {
        "H": 1.00783, "C": 12.0000, "N": 14.00307, "O": 15.99491,
        "F": 18.99840, "Na": 22.98977, "Si": 27.97693, "P": 30.97376,
        "S": 31.97207, "Cl": 34.96885, "Br": 78.91833, "I": 126.90447,
        "B": 11.00931, "Al": 26.98154, "Ca": 39.96259, "Fe": 55.93494,
    }

    # ========== COMMON NEUTRAL LOSSES ==========
    _NEUTRAL_LOSSES = [
        (1.0078, "H•", "Hydrogen radical loss"),
        (2.0157, "H₂", "Hydrogen molecule loss"),
        (15.9949, "O", "Oxygen atom loss (rare)"),
        (17.0027, "OH•", "Hydroxyl radical loss (alcohols)"),
        (18.0106, "H₂O", "Water loss (alcohols, carboxylic acids) — VERY COMMON"),
        (26.0015, "C₂H₂", "Acetylene loss (aromatics)"),
        (27.9949, "CO", "Carbon monoxide loss (carbonyl compounds) — COMMON"),
        (28.0313, "C₂H₄", "Ethylene loss ( McLafferty-type)"),
        (28.0313, "CH₂=CH₂", "Ethylene elimination"),
        (29.0027, "CHO•", "Formyl radical loss (aldehydes)"),
        (30.0106, "CH₂O", "Formaldehyde loss"),
        (31.0184, "CH₃O•", "Methoxy radical loss (ethers/methyl esters)"),
        (42.0106, "CH₂=C=O", "Ketene loss (common in carbonyls)"),
        (43.0058, "C₂H₃O", "Acetyl radical / acetaldehyde loss"),
        (44.9977, "CO₂", "Carbon dioxide loss (decarboxylation) — COMMON"),
        (45.0058, "CH₃CHO", "Acetaldehyde loss"),
        (46.0058, "NO₂", "Nitro group loss (nitro compounds)"),
        (46.0058, "CH₂CH₂OH", "Ethanol loss"),
        (47.0132, "CH₃OCH₂", "Methyl ether fragment"),
        (48.0000, "CH₃SH", "Methanethiol loss"),
        (49.9923, "CH₃Cl", "Methyl chloride loss"),
        (56.0262, "C₄H₈", "Butene loss (McLafferty)"),
        (57.0146, "C₄H₉•", "Butyl radical loss"),
        (58.0419, "C₃H₆O", "Acetone/propionaldehyde loss"),
        (59.0133, "C₂H₃O₂", "Carboxyl radical loss"),
        (60.0212, "C₂H₄O₂", "Acetic acid loss"),
        (64.9621, "SO₂", "Sulfur dioxide loss (sulfonyl compounds)"),
        (73.0468, "C₃H₇O", "Propoxy/propionyl loss"),
        (79.9169, "Br•", "Bromine radical loss"),
        (91.0548, "C₇H₇•", "Tropylium/benzyl loss"),
        (105.0704, "C₈H₉•", "Styryl/tropyl loss"),
        (127.0012, "I•", "Iodine radical loss"),
    ]

    # ========== CHARACTERISTIC FRAGMENTATION PATTERNS ==========
    _FRAGMENTATION_RULES = {
        # Alcohol patterns
        "alcohol_primary": {
            "key_fragments": "M-18 (H₂O loss), M-29 (CHO loss), M-46 (2×H₂O)",
            "base_peak_candidate": "m/z 31 (CH₂=OH⁺)",
            "notes": "α-cleavage gives CH₂=OH⁺ at m/z 31; dehydration prominent",
        },
        "alcohol_secondary": {
            "key_fragments": "M-18, α-cleavage both sides of OH-bearing C",
            "base_peak_candidate": "m/z 45 (CH₃CH=OH⁺) or larger RCH=OH⁺",
            "notes": "α-cleavage gives RCH=OH⁺; dehydration common",
        },
        "alcohol_tertiary": {
            "key_fragments": "M-15 (CH₃ loss), M-18 (H₂O), M-33 (CH₃+H₂O)",
            "base_peak_candidate": "R₃C⁺ (tertiary carbocation)",
            "notes": "Carbocation stability drives fragmentation; often lose CH₃ then H₂O",
        },

        # Carbonyl patterns
        "aldehyde": {
            "key_fragments": "M-1 (H loss → acylium), M-29 (CHO loss), McLafferty if γ-H present",
            "base_peak_candidate": "m/z 29 (CHO⁺) or M-1 (R-C≡O⁺ acylium)",
            "notes": "β-cleavage gives acylium ion R-C≡O⁺ (strong peak); McLafferty if chain ≥3C",
        },
        "ketone": {
            "key_fragments": "α-cleavage on both sides of C=O, McLafferty rearrangement",
            "base_peak_candidate": "Acyl ion R-C≡O⁺ or enol oxonium ion",
            "notes": "α-cleavage gives acylium ions; McLafferty if γ-H on either side; loss of CO (28) from acyl ions",
        },
        "carboxylic_acid": {
            "key_fragments": "M-17 (OH loss), M-18 (H₂O), M-45 (COOH loss), M-44 (CO₂)",
            "base_peak_candidate": "m/z 45 (COOH⁺) for small acids; m/z 60 (McLafferty)",
            "notes": "McLafferty rearrangement gives m/z 60 for carboxylic acids with γ-H; decarboxylation (M-44) common",
        },
        "ester": {
            "key_fragments": "α-cleavage next to C=O and O, McLafferty, loss of OR•",
            "base_peak_candidate": "R-C≡O⁺ (acyl) or ⁺O=CH-R (oxonium) or m/z 74 (McLafferty for methyl esters)",
            "notes": "McLafferty gives m/z 74 for ethyl/methyl esters with γ-H; two α-cleavage sites",
        },

        # Amine patterns
        "amine_primary_aliphatic": {
            "key_fragments": "α-cleavage giving iminium ion CH₂=NH₂⁺ (m/z 30)",
            "base_peak_candidate": "m/z 30 (CH₂=NH₂⁺) — very characteristic!",
            "notes": "α-cleavage gives CH₂=NH₂⁺ at m/z 30 (diagnostic for primary amines); also M-1",
        },
        "amine_secondary_aliphatic": {
            "key_fragments": "α-cleavage giving RCH=NH-CH₃⁺ type ions",
            "base_peak_candidate": "m/z 44, 58, 72... (iminium series, +14 per CH₂)",
            "notes": "Iminium ion series: m/z 30, 44, 58, 72, 86... (CnH2n+2N⁺)",
        },
        "aniline_aromatic": {
            "key_fragments": "M-HCN (M-27), M-HCN-H (M-28), M-1",
            "base_peak_candidate": "m/z 66 (C₅H₆⁺, cyclopentadiene-like after HCN loss)",
            "notes": "Aromatic amine loses HCN readily; strong M peak usually present",
        },

        # Ether patterns
        "ether_aliphatic": {
            "key_fragments": "α-cleavage on both sides of O, loss of alkyl radicals",
            "base_peak_candidate": "m/z 45 (CH₂=OH-CH₃⁺), 59, 73... (oxonium series)",
            "notes": "α-cleavage gives oxonium ions: CnH2n+1O⁺ series (m/z 31, 45, 59, 73...)",
        },

        # Halogenated compound patterns
        "chlorinated": {
            "key_fragments": "M, M+2 (³⁷Cl, ~3:1 ratio); Cl• loss (M-35)",
            "base_peak_candidate": "Aryl-Cl: M-Cl (benzoyl-type); Alkyl-Cl: M-Cl/HCl",
            "notes": "Characteristic 3:1 isotope pattern for single Cl; α-cleavage next to Cl",
        },
        "brominated": {
            "key_fragments": "M, M+2 (⁸¹Br, ~1:1 ratio); Br• loss (M-79)",
            "base_peak_candidate": "Often M-Br or aryl cation",
            "notes": "Characteristic 1:1 isotope pattern for single Br; easy to identify",
        },

        # Aromatic patterns
        "benzene_homolog": {
            "key_fragments": "M (often strong), M-1 (H loss), sequential C₂H₂ (26) losses",
            "base_peak_candidate": "m/z 77 (C₆H₅⁺ phenyl cation), m/z 51 (C₄H₃⁺)",
            "notes": "Strong M+ peak; tropylium m/z 91 for benzyl; phenyl m/z 77; stepwise C₂H₂ loss",
        },
        "alkyl_benzene": {
            "key_fragments": "Benzylic cleavage → tropylium ion m/z 91 (BASE PEAK!)",
            "base_peak_candidate": "m/z 91 (C₇H₇⁺ tropylium) — diagnostic for benzyl groups!",
            "notes": "Benzylic cleavage gives resonance-stabilized tropylium ion at m/z 91; very intense",
        },

        # Nitro compound patterns
        "nitro_aliphatic": {
            "key_fragments": "M-NO (M-30), M-NO₂ (M-46), M-NO₂-H (M-47)",
            "base_peak_candidate": "Often NO₂⁺ (m/z 46) or hydrocarbon fragment",
            "notes": "Nitro group losses dominate; weak M+ peak",
        },
        "nitro_aromatic": {
            "key_fragments": "M-NO (M-30), M-NO₂ (M-46), NO₂⁺ (m/z 46)",
            "base_peak_candidate": "m/z 46 (NO₂⁺) or M-NO (phenoxide-type)",
            "notes": "Strong NO₂⁺ peak at m/z 46; characteristic of nitroaromatics",
        },

        # Nitrile patterns
        "nitrile_aliphatic": {
            "key_fragments": "M-1 (H loss), α-cleavage giving CN-containing fragments",
            "base_peak_candidate": "m/z 41 (CH₂=C=N⁺ or C₃H₅⁺) for propionitrile-type",
            "notes": "McLafferty possible if γ-H; cyano-containing fragments common",
        },
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
            logger.warning("RDKit not available for MassSpecFragmenter")

    def _parse_formula(self, smiles_or_formula: str) -> dict:
        """Parse SMILES or formula into element counts."""
        s = smiles_or_formula.strip()

        # Try as SMILES first
        if self._rdkit_available:
            try:
                mol = self.Chem.MolFromSmiles(s)
                if mol:
                    from rdkit.Chem import Descriptors
                    return {
                        "type": "SMILES",
                        "smiles": self.Chem.MolToSmiles(mol),
                        "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
                        "mol_wt": round(Descriptors.MolWt(mol), 2),
                    }
            except Exception:
                pass

        # Try as molecular formula
        return self._parse_molecular_formula(s)

    def _parse_molecular_formula(self, formula: str) -> dict:
        """Parse molecular formula like C6H12O6 into element counts and mass."""
        import re
        elements = {}
        matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        total_mass = 0.0

        for elem, count_str in matches:
            if not elem:
                continue
            count = int(count_str) if count_str else 1
            elements[elem] = elements.get(elem, 0) + count
            if elem in self._ATOMIC_MASSES:
                total_mass += self._ATOMIC_MASSES[elem] * count

        return {
            "type": "Formula",
            "formula": formula,
            "elements": elements,
            "exact_mass": round(total_mass, 4),
        }

    def _detect_functional_groups(self, smi: str) -> list:
        """Detect functional groups from SMILES/heuristic analysis."""
        groups = []

        if re.search(r'C\(=O\)O[H]?', smi):
            groups.append("carboxylic_acid")
        elif re.search(r'C\(=O\)[Oc]', smi):
            groups.append("ester")
        elif re.search(r'(C=O)[Hh]|[Hh]C=O', smi) or smi.endswith('C=O'):
            groups.append("aldehyde")
        elif re.search(r'C\(=O\)[Cc]', smi) and not re.search(r'C\(=O\)[OoNn]', smi):
            groups.append("ketone")
        elif re.search(r'[Cc]O[H]?', smi) and not re.search(r'C\(=O\)O', smi):
            groups.append("alcohol")
        if re.search(r'CN|C#N', smi):
            groups.append("nitrile")
        if re.search(r'[Nn][Hh]2|[Nn][Hh](?![Cc]=O)', smi) and not re.search(r'C\(=O\)[Nn]', smi):
            groups.append("amine")
        if re.search(r'C\(=O\)[Nn]', smi):
            groups.append("amide")
        if re.search(r'Cl', smi):
            groups.append("chlorinated")
        if re.search(r'Br', smi):
            groups.append("brominated")
        if 'c' in smi.lower():
            groups.append("aromatic")
            if re.search(r'c1ccccc1[Cc]', smi):
                groups.append("alkyl_benzene")
        if re.search(r'NO2|N.*O=O', smi):
            groups.append("nitro")
        if re.search(r'[Cc]O[Cc]', smi) and not re.search(r'C\(=O\)O', smi):
            groups.append("ether")
        if re.search(r'S', smi):
            groups.append("sulfur_compound")

        return groups if groups else ["hydrocarbon"]

    def _predict_fragments(self, smi: str, mode: str = "EI+") -> dict:
        """Predict fragmentation pattern based on functional groups."""
        groups = self._detect_functional_groups(smi)
        parsed = self._parse_formula(smi)
        M = parsed.get("exact_mass", 0)

        if M <= 0:
            raise ChemMCPError(f"Could not determine molecular weight for '{smi}'")

        fragments = []
        base_peak_candidates = []

        # Apply fragmentation rules for each detected functional group
        applied_rules = set()
        for grp in groups:
            if grp in self._FRAGMENTATION_RULES and grp not in applied_rules:
                rule = self._FRAGMENTATION_RULES[grp]
                applied_rules.add(grp)

                # Generate specific fragments based on rules
                if grp == "alcohol_primary":
                    fragments.append({"mz": round(M - 18.0106, 1), "formula": "[M-H₂O]⁺", "pathway": "Dehydration (H₂O loss)", "intensity": "moderate"})
                    fragments.append({"mz": 31.0, "formula": "CH₂=OH⁺", "pathway": "α-Cleavage (primary alcohol)", "intensity": "strong (base peak candidate)"})
                    base_peak_candidates.append((31.0, "CH₂=OH⁺"))

                elif grp == "alcohol_secondary":
                    fragments.append({"mz": round(M - 18.0106, 1), "formula": "[M-H₂O]⁺", "pathway": "Dehydration", "intensity": "moderate"})
                    fragments.append({"mz": 45.0, "formula": "CH₃CH=OH⁺", "pathway": "α-Cleavage (secondary alcohol)", "intensity": "strong"})
                    base_peak_candidates.append((45.0, "CH₃CH=OH⁺"))

                elif grp == "alcohol_tertiary":
                    fragments.append({"mz": round(M - 15.0, 1), "formula": "[M-CH₃]⁺", "pathway": "Methyl loss from tertiary C", "intensity": "moderate"})
                    fragments.append({"mz": round(M - 18.0106, 1), "formula": "[M-H₂O]⁺", "pathway": "Dehydration", "intensity": "moderate"})

                elif grp == "aldehyde":
                    fragments.append({"mz": round(M - 1.0078, 1), "formula": "[M-H]⁺ (acylium)", "pathway": "β-Cleavage → acylium ion", "intensity": "strong"})
                    fragments.append({"mz": 29.0, "formula": "CHO⁺", "pathway": "Formyl cation (small aldehydes)", "intensity": "moderate"})
                    if M > 50:  # McLafferty possible
                        fragments.append({"mz": round(M - 44.0, 1), "formula": "[M-CH₃CHO]⁺", "pathway": "McLafferty rearrangement", "intensity": "strong"})
                    base_peak_candidates.append((round(M - 1.0078, 1), "acylium [M-H]⁺"))
                    base_peak_candidates.append((29.0, "CHO⁺"))

                elif grp == "ketone":
                    fragments.append({"mz": round(M - 27.9949, 1), "formula": "[M-CO]⁺", "pathway": "CO loss from molecular ion", "intensity": "moderate"})
                    fragments.append({"mz": round(M - 43.0, 1), "formula": "[M-CH₃CO]⁺", "pathway": "Acetyl loss (α-cleavage)", "intensity": "strong"})
                    if M > 60:
                        fragments.append({"mz": round(M - 58.0, 1), "formula": "[M-McLafferty]⁺", "pathway": "McLafferty rearrangement", "intensity": "strong"})

                elif grp == "carboxylic_acid":
                    fragments.append({"mz": round(M - 44.9977, 1), "formula": "[M-CO₂]⁺", "pathway": "Decarboxylation", "intensity": "strong"})
                    fragments.append({"mz": round(M - 18.0106, 1), "formula": "[M-H₂O]⁺", "pathway": "Dehydration", "intensity": "moderate"})
                    fragments.append({"mz": 60.0, "formula": "C₂H₄O₂⁺", "pathway": "McLafferty rearrangement (characteristic!)", "intensity": "strong (base peak)"})
                    base_peak_candidates.append((60.0, "C₂H₄O₂⁺ (McLafferty)"))

                elif grp == "ester":
                    fragments.append({"mz": round(M - 31.0184, 1), "formula": "[M-OCH₃]⁺", "pathway": "Alkoxy loss (α-cleavage)", "intensity": "strong"})
                    fragments.append({"mz": round(M - 43.0, 1), "formula": "[M-CH₃CO]⁺", "pathway": "Acyl loss (α-cleavage)", "intensity": "strong"})
                    fragments.append({"mz": 74.0, "formula": "C₃H₆O₂⁺? or C₂H₆O₂N?", "pathway": "McLafferty rearrangement (methyl/ethyl ester)", "intensity": "strong (base peak)"})
                    base_peak_candidates.append((74.0, "McLafferty ester"))

                elif grp == "amine_primary_aliphatic":
                    fragments.append({"mz": 30.0, "formula": "CH₂=NH₂⁺", "pathway": "α-Cleavage (DIAGNOSTIC for primary amine!)", "intensity": "very strong (base peak)"})
                    base_peak_candidates.append((30.0, "CH₂=NH₂⁺ (diagnostic)"))

                elif grp == "amine_secondary_aliphatic":
                    fragments.append({"mz": 44.0, "formula": "CH₂=NH-CH₃⁺", "pathway": "α-Cleavage (iminium ion)", "intensity": "strong"})
                    fragments.append({"mz": 58.0, "formula": "(CH₃)₂C=NH₂⁺? or C₃H₈N⁺", "pathway": "Higher homolog iminium", "intensity": "moderate"})

                elif grp == "amide":
                    fragments.append({"mz": round(M - 44.9977, 1), "formula": "[M-CO₂]⁺", "pathway": "Loss of CONH (isocyanate) or related", "intensity": "moderate"})
                    fragments.append({"mz": 44.0, "formula": "CONH₂⁺ or CH₂=NH-OH?", "pathway": "Amide-specific fragment", "intensity": "moderate"})

                elif grp == "ether_aliphatic":
                    fragments.append({"mz": 45.0, "formula": "CH₂=OH-CH₃⁺", "pathway": "α-Cleavage (oxonium ion)", "intensity": "strong"})
                    fragments.append({"mz": 59.0, "formula": "C₂H₅O=CH₂⁺ or CH₃O=CH-CH₃⁺", "pathway": "Higher oxonium ion", "intensity": "moderate"})

                elif grp == "chlorinated":
                    fragments.append({"mz": round(M - 35.0, 1), "formula": "[M-Cl]⁺", "pathway": "Chlorine radical loss", "intensity": "moderate"})
                    fragments.append({"mz": round(M - 36.0, 1), "formula": "[M-HCl]⁺", "pathway": "HCl elimination", "intensity": "moderate"})

                elif grp == "brominated":
                    fragments.append({"mz": round(M - 79.0, 1), "formula": "[M-Br]⁺", "pathway": "Bromine radical loss", "intensity": "moderate"})
                    fragments.append({"mz": round(M - 80.0, 1), "formula": "[M-HBr]⁺", "pathway": "HBr elimination", "intensity": "moderate"})

                elif grp == "aromatic":
                    fragments.append({"mz": round(M - 1.0078, 1), "formula": "[M-H]⁺", "pathway": "H loss from aromatic ring", "intensity": "moderate"})
                    fragments.append({"mz": 77.0, "formula": "C₆H₅⁺", "pathway": "Phenyl cation", "intensity": "strong"})

                elif grp == "alkyl_benzene":
                    fragments.append({"mz": 91.0, "formula": "C₇H₇⁺ (tropylium)", "pathway": "Benzylic cleavage → tropylium ion (DIAGNOSTIC!)", "intensity": "very strong (base peak)"})
                    fragments.append({"mz": 77.0, "formula": "C₆H₅⁺", "pathway": "Phenyl cation (tropylium - CH₂)", "intensity": "moderate"})
                    fragments.append({"mz": 65.0, "formula": "C₅H₅⁺", "pathway": "Cyclopentadienyl cation", "intensity": "weak"})
                    base_peak_candidates.append((91.0, "C₇H₇⁺ tropylium (diagnostic!)"))

                elif grp == "nitro":
                    fragments.append({"mz": 46.0, "formula": "NO₂⁺", "pathway": "Nitro cation (DIAGNOSTIC for nitro compounds!)", "intensity": "strong"})
                    fragments.append({"mz": round(M - 30.0, 1), "formula": "[M-NO]⁺", "pathway": "NO loss", "intensity": "moderate"})
                    fragments.append({"mz": round(M - 46.0, 1), "formula": "[M-NO₂]⁺", "pathway": "NO₂ loss", "intensity": "moderate"})
                    base_peak_candidates.append((46.0, "NO₂⁺ (diagnostic)"))

                elif grp == "nitrile":
                    fragments.append({"mz": 41.0, "formula": "CH₂=C=N⁺ or C₃H₅⁺", "pathway": "Cyano-containing fragment (common)", "intensity": "strong"})
                    fragments.append({"mz": round(M - 27.0, 1), "formula": "[M-HCN]⁺", "pathway": "HCN loss", "intensity": "moderate"})

        # If no specific rules matched, apply generic hydrocarbon fragmentation
        if not fragments:
            fragments = self._generic_hydrocarbon_fragmentation(M, smi)

        # Determine most likely base peak
        if base_peak_candidates:
            base_mz, base_formula = base_peak_candidates[0]
        else:
            base_mz = M
            base_formula = "M⁺•"

        return {
            "molecular_ion_mz": round(M, 1),
            "molecular_ion_intensity": "variable (often moderate for EI)",
            "base_peak_mz": base_mz,
            "base_peak_formula": base_formula,
            "fragments": fragments[:15],  # Top 15 fragments
            "detected_functional_groups": groups,
            "fragmentation_rules_applied": list(applied_rules),
        }

    def _generic_hydrocarbon_fragmentation(self, M: float, smi: str) -> list:
        """Generic fragmentation for unrecognized/hydrocarbon structures."""
        fragments = [
            {"mz": round(M, 1), "formula": "M⁺•", "pathway": "Molecular ion", "intensity": "weak-moderate"},
            {"mz": round(M - 15.0, 1), "formula": "[M-CH₃]⁺", "pathway": "Methyl loss", "intensity": "weak"},
            {"mz": round(M - 29.0, 1), "formula": "[M-CHO/C₂H₅]⁺", "pathway": "Ethyl/formyl loss", "intensity": "weak"},
            {"mz": 43.0, "formula": "C₃H₇⁺", "pathway": "Propyl cation (common hydrocarbon fragment)", "intensity": "moderate"},
            {"mz": 57.0, "formula": "C₄H₉⁺", "pathway": "Butyl cation", "intensity": "moderate"},
            {"mz": 41.0, "formula": "C₃H₅⁺ (allyl)", "pathway": "Allyl cation (resonance stabilized)", "intensity": "moderate-strong"},
            {"mz": 29.0, "formula": "C₂H₅⁺/CHO⁺", "pathway": "Ethyl/formyl cation", "intensity": "moderate"},
            {"mz": 27.0, "formula": "C₂H₃⁺", "pathway": "Vinyl/acetylene cation", "intensity": "weak"},
        ]
        return fragments

    def _run_base(self, smiles: str, ionization_mode: str = "EI+") -> dict:
        """
        Predict MS fragmentation.

        Args:
            smiles: SMILES string or molecular formula
            ionization_mode: Ionization mode

        Returns:
            Dict with fragmentation data
        """
        if not smiles:
            raise ChemMCPError("SMILES string or molecular formula is required.")

        result = self._predict_fragments(smiles.strip(), ionization_mode)
        result["ionization_mode"] = ionization_mode
        result["input"] = smiles.strip()

        # Add general notes
        result["general_notes"] = (
            "Electron Impact (EI) Mass Spectrometry:\n"
            "• 70 eV electron beam generates M⁺• (radical cation)\n"
            "• Fragmentation occurs via unimolecular dissociation\n"
            "• Base peak = most intense peak (100% relative intensity)\n"
            "• Neutral fragments are NOT detected (only charged species)\n"
            "• Isotope patterns: Cl shows 3:1 (³⁵Cl:³⁷Cl), Br shows 1:1 (⁷⁹Br:⁸¹Br)\n"
            "• Nitrogen rule: odd nominal mass → odd number of N atoms (and vice versa)"
        )

        return {"fragmentation_data": result}

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'smiles [mode]'")

        smiles = parts[0]
        mode = parts[1] if len(parts) > 1 else "EI+"

        return self._run_base(smiles, mode)
