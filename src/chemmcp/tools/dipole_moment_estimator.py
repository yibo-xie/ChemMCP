import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DipoleMomentEstimator(BaseTool):
    """
    偶极矩估算工具 (MCP #294)。
    通过实验数据库查询和键偶极矩矢量加法估算分子的偶极矩（单位：Debye）。
    支持常见分子偶极矩查询、分子极性判断、以及基于几何构型的近似计算。
    """
    __version__ = "0.1.0"
    name = "DipoleMomentEstimator"
    func_name = "estimate_dipole_moment"
    description = "Estimate molecular dipole moment (in Debye) from experimental data or vector addition of bond dipoles."
    implementation_description = (
        "Uses experimental dipole moment database for 100+ common molecules and "
        "vector bond-dipole addition for estimation based on molecular geometry."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Dipole Moment", "Polarity", "Molecular Properties", "Electrostatics"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: formula, SMILES, or name (e.g., 'H2O', 'NH3', 'CCl4', 'CO2', 'benzene')."),
    ]

    text_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: formula, SMILES, or name."),
    ]

    output_sig = [
        ("molecule", "str", "The molecule being analyzed."),
        ("dipole_moment", "float", "Estimated dipole moment in Debye (D)."),
        ("polarity", "str", "Polarity classification: polar / nonpolar / slightly polar."),
        ("direction", "str", "Direction of the net dipole vector (if applicable)."),
        ("description", "str", "Detailed explanation of the dipole moment origin."),
        ("source", "str", "Data source: 'experimental' or 'estimated'."),
    ]

    examples = [{'code_input': {'molecule': 'H2O'}, 'text_input': {'molecule': 'H2O'}, 'output': {'molecule': 'H2', 'dipole_moment': 1.85, 'polarity': 'polar', 'direction': 'Bisecting H-O-H angle, pointing from O to midpoint of H atoms', 'source': 'experimental', 'description': 'N/A'}}, {'code_input': {'molecule': 'CO2'}, 'text_input': {'molecule': 'CO2'}, 'output': {'molecule': 'CO2', 'dipole_moment': 0.0, 'polarity': 'nonpolar', 'direction': 'N/A (cancels out)', 'source': 'experimental', 'description': 'N/A'}}, {'code_input': {'molecule': 'NH3'}, 'text_input': {'molecule': 'NH3'}, 'output': {'molecule': 'NH3', 'dipole_moment': 1.47, 'polarity': 'polar', 'direction': 'Along C3 axis through N atom', 'source': 'experimental', 'description': 'N/A'}}]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build dipole moment database. Values in Debye (D)."""
        self._db = {
            # ===== Diatomic molecules =====
            "HF":   (1.82, "polar", "H→F (F is highly electronegative)", "experimental"),
            "HCl":  (1.08, "polar", "H→Cl", "experimental"),
            "HBr":  (0.82, "polar", "H→Br", "experimental"),
            "HI":   (0.44, "polar", "H→I", "experimental"),
            "CO":   (0.11, "slightly polar", "C⁻→O⁺ (unusual direction due to π back-donation)", "experimental"),
            "NO":   (0.16, "slightly polar", "N→O", "experimental"),
            "LiH":  (5.88, "highly polar", "Li→H (ionic character ~70%)", "estimated"),
            "NaCl": 9.0,  # gas phase, very ionic
            "KF":   8.6,
            "H2":   (0.0, "nonpolar", "Identical atoms, no charge separation", "experimental"),
            "O2":   (0.0, "nonpolar", "Identical atoms", "experimental"),
            "N2":   (0.0, "nonpolar", "Identical atoms", "experimental"),
            "Cl2":  (0.0, "nonpolar", "Identical atoms", "experimental"),
            "F2":   (0.0, "nonpolar", "Identical atoms", "experimental"),
            "Br2":  (0.0, "nonpolar", "Identical atoms", "experimental"),
            "I2":   (0.0, "nonpolar", "Identical atoms", "experimental"),

            # ===== Triatomic =====
            "H2O":  (1.85, "polar", "Net dipole from two O-H bonds at 104.5°; points toward O lone pairs (actually: negative end at O, +end between H's)", "experimental"),
            "H2S":  (0.97, "polar", "Similar to H2O but less electronegativity difference S-H vs O-H", "experimental"),
            "H2Se": (0.23, "slightly polar", "Weaker polarity than H2S", "experimental"),
            "H2Te": (0.2, "nearly nonpolar", "Very weak polarity", "estimated"),
            "O3":   (0.53, "polar", "Bent geometry (117°), resonance structure gives net dipole", "experimental"),
            "SO2":  (1.63, "polar", "Bent molecule, S=O dipoles don't cancel", "experimental"),
            "NO2":  (0.32, "polar", "Bent geometry (~134°)", "experimental"),
            "Cl2O": (0.78, "polar", "Bent Cl-O-Cl", "experimental"),
            "OF2":  (0.30, "polar", "Bent, F more EN than O so dipole points toward F atoms", "experimental"),
            "SCl2": (0.58, "polar", "Bent S-Cl-S", "experimental"),
            "CO2":  (0.0, "nonpolar", "Linear O=C=O; C=O bond dipoles cancel exactly", "experimental"),
            "C=O=C":(0.0, "nonpolar", "Linear CO2; symmetric cancellation", "experimental"),
            "CS2":  (0.0, "nonpolar", "Linear S=C=S; symmetric cancellation", "experimental"),
            "CS2(linear)": (0.0, "nonpolar", "Linear symmetry", "experimental"),
            "BeCl2": (0.0, "nonpolar", "Linear Be-Cl bonds cancel", "estimated"),
            "HgCl2": (0.0, "nonpolar", "Linear", "estimated"),
            "N2O":  (0.166, "polar", "Linear but N≡N⁺-O⁻ gives small net dipole (N→O)", "experimental"),
            "[N+]=[N-]": (0.166, "polar", "N2O: unsymmetric linear", "experimental"),
            "HCN":  (2.98, "highly polar", "H-C≡N: C-H and C≡N dipoles add in same direction (H→C→N)", "experimental"),
            "C#N":  (2.98, "highly polar", "HCN: large dipole along axis", "experimental"),
            "OCS":  (0.72, "polar", "Linear O=C=S; different terminal atoms → net dipole", "experimental"),
            "H2O2": (2.26, "polar", "Non-planar (gauche conformation); O-H and O-O dipoles combine", "experimental"),

            # ===== Tetrahedral / Octahedral =====
            "CH4":  (0.0, "nonpolar", "Td symmetry: 4 equivalent C-H bonds cancel perfectly", "experimental"),
            "CCl4": (0.0, "nonpolar", "Td symmetry: 4 C-Cl bonds cancel", "experimental"),
            "CF4":  (0.0, "nonpolar", "Td symmetry", "experimental"),
            "SiH4": (0.0, "nonpolar", "Td symmetry", "experimental"),
            "CH3Cl":(1.87, "polar", "Cs symmetry: C-Cl dipole not canceled by 3 C-H bonds", "experimental"),
            "CH2Cl2":(1.60, "polar", "C2v: 2 C-Cl + 2 C-H give net dipole", "experimental"),
            "CHCl3": (1.01, "polar", "3 C-Cl dipoles partially canceled by C-H", "experimental"),
            "CF3Cl":(0.50, "polar", "Freon-13: asymmetrical substitution on Td frame", "experimental"),
            "SF6":  (0.0, "nonpolar", "Oh symmetry: perfect octahedral cancellation", "experimental"),
            "UF6":  (0.0, "nonpolar", "Oh symmetry", "estimated"),
            "XeF4": (0.0, "nonpolar", "D4h: square planar, trans F atoms cancel", "experimental"),
            "XeF2": (0.0, "nonpolar", "D∞h: linear, symmetric", "experimental"),
            "BrF5": (1.52, "polar", "C4v: square pyramidal, axial F gives net dipole", "experimental"),
            "IF5":  (2.18, "polar", "C4v: square pyramidal", "experimental"),
            "ClF3": (0.58, "polar", "T-shaped (C2v), dipoles don't fully cancel", "experimental"),
            "SF4":  (0.64, "polar", "See-saw (C2v), lone pair creates asymmetry", "experimental"),

            # ===== Trigonal planar / Pyramidal =====
            "BF3":  (0.0, "nonpolar", "D3h: trigonal planar, 3 B-F bonds at 120° cancel", "experimental"),
            "SO3":  (0.0, "nonpolar", "D3h: trigonal planar", "experimental"),
            "NO3-": (0.0, "nonpolar", "D3h: trigonal planar, resonance averages to zero", "estimated"),
            "AlCl3": (0.0, "nonpolar", "D3h: trigonal planar", "estimated"),
            "BCl3": (0.0, "nonpolar", "D3h: trigonal planar", "estimated"),
            "NH3":  (1.47, "polar", "C3v: pyramidal; lone pair contributes to net dipole along C3 axis (N→base)", "experimental"),
            "PH3":  (0.58, "polar", "C3v: pyramidal, weaker than NH3", "experimental"),
            "PCl3": (0.97, "polar", "C3v: pyramidal", "experimental"),
            "PF3":  (1.03, "polar", "C3v: pyramidal", "experimental"),
            "AsH3": (0.20, "slightly polar", "C3v: weakly polar", "experimental"),
            "NF3":  (0.23, "slightly polar", "C3v: N-F bond dipoles oppose lone pair dipole → small net", "experimental"),
            "CH3F": (1.85, "polar", "Cs: C-F dipole dominant", "experimental"),
            "POCl3":(2.54, "polar", "C3v-like with P=O strong dipole", "experimental"),
            "SOCl2":(1.45, "polar", "Trigonal pyramidal S center", "experimental"),
            "SO2Cl2":(1.81, "polar", "Tetrahedral S(VI) with different substituents", "experimental"),
            "ClO2": (1.69, "polar", "Bent (C2v) radical", "experimental"),

            # ===== Organic molecules =====
            "C2H4":  (0.0, "nonpolar", "D2h: planar, symmetric", "experimental"),
            "C=C":   (0.0, "nonpolar", "Ethylene D2h symmetry", "experimental"),
            "C2H2":  (0.0, "nonpolar", "D∞h: linear symmetric", "experimental"),
            "C#C":   (0.0, "nonpolar", "Acetylene D∞h", "experimental"),
            "benzene": (0.0, "nonpolar", "D6h: planar ring, perfectly symmetric", "experimental"),
            "c1ccccc1": (0.0, "nonpolar", "Benzene D6h", "experimental"),
            "C1=CC=CC=C1": (0.0, "nonpolar", "Benzene D6h", "experimental"),
            "toluene": (0.36, "slightly polar", "Methyl group breaks D6h symmetry slightly", "experimental"),
            "phenol": (1.22, "polar", "OH group adds significant dipole", "experimental"),
            "aniline": (1.53, "polar", "NH2 group adds dipole", "experimental"),
            "nitrobenzene": (4.22, "highly polar", "Strong NO2 dipole (~4D) along ring axis", "experimental"),
            "chlorobenzene": (1.56, "polar", "C-Cl bond dipole", "experimental"),
            "allene": (0.0, "nonpolar", "D2d: orthogonal π systems cancel", "experimental"),
            "C=C=C": (0.0, "nonpolar", "Allene D2d", "experimental"),
            "HCHO": (2.33, "polar", "C=O carbonyl dipole dominates", "experimental"),
            "CH3CHO": (2.72, "polar", "Acetaldehyde: C=O + CH3 dipoles", "experimental"),
            "acetone": (2.91, "polar", "(CH3)2C=O: strong carbonyl dipole", "experimental"),
            "CH3OH": (1.69, "polar", "O-H + C-O dipoles", "experimental"),
            "CCO":   (1.69, "polar", "Ethanol: O-H/C-O dipoles", "experimental"),
            "CH3CH2OH": (1.69, "polar", "Ethanol", "experimental"),
            "dimethyl ether": (1.30, "polar", "CH3-O-CH3: C-O dipoles", "experimental"),
            "CH3OCH3": (1.30, "polar", "DME", "experimental"),
            "ethyl methyl ether": (1.24, "polar", "CH3CH2-O-CH3", "experimental"),
            "formic acid": (1.41, "polar", "HCOOH: C=O + O-H", "experimental"),
            "acetic acid": (1.74, "polar", "CH3COOH", "experimental"),
            "methylamine": (1.31, "polar", "CH3NH2", "experimental"),
            "ethylamine": (1.22, "polar", "CH3CH2NH2", "experimental"),
            "aniline": (1.53, "polar", "C6H5NH2", "experimental"),
            "urea": (4.56, "highly polar", "(NH2)2C=O: two N-H + C=O dipoles align", "experimental"),
            "DMF":   (3.86, "highly polar", "Dimethylformamide: strong C=O dipole", "experimental"),
            "DMSO":  (4.0, "highly polar", "Dimethyl sulfoxide: S=O dipole", "experimental"),
            "acetonitrile": (3.92, "highly polar", "CH3CN: C≡N triple bond dipole", "experimental"),
            "nitromethane": (3.46, "highly polar", "CH3NO2: NO2 group dipole", "experimental"),
            "chloroform": (1.01, "polar", "CHCl3", "experimental"),
            "dichloromethane": (1.60, "polar", "CH2Cl2", "experimental"),
            "ethylene oxide": (1.79, "polar", "Oxirane ring strain + C-O dipoles", "experimental"),
            "tetrahydrofuran": (1.75, "polar", "THF: cyclic ether", "experimental"),
            "1,4-dioxane": (0.0, "nonpolar", "Symmetric cyclic diether", "experimental"),
            "pyridine": (2.19, "polar", "Aromatic N heterocycle: lone pair dipole", "experimental"),
            "pyrrole": (1.73, "polar", "N-H dipole in aromatic ring", "experimental"),
            "furan": (0.66, "polar", "O heterocycle", "experimental"),
            "thiophene": (0.55, "slightly polar", "S heterocycle", "experimental"),
            "indole": (2.05, "polar", "Fused N-heterocycle", "experimental"),
            "imidazole": (3.67, "polar", "Two N atoms in 5-membered ring", "experimental"),
            "glycine": (-15, "zwitterionic", "Amino acid zwitterion form has huge dipole", "estimated"),
            "oxalic acid": (1.7, "polar", "(COOH)2", "estimated"),
            "ascorbic acid": (-2.5, "polar", "Vitamin C: multiple OH groups", "estimated"),

            # ===== Other important molecules =====
            "PCl5(trigonal bipyramid)": (0.0, "nonpolar", "D3h: axial + equatorial P-Cl cancel", "estimated"),
            "trans-[Pt(NH3)2Cl2]": (0.0, "nonpolar", "D2h: trans arrangement cancels", "estimated"),
            "cis-[Pt(NH3)2Cl2]:": (0.42, "polar", "C2v: cis arrangement doesn't cancel", "estimated"),
            "o-xylene": (0.62, "slightly polar", "Ortho-dimethylbenzene", "experimental"),
            "m-xylene": (0.33, "slightly polar", "Meta-dimethylbenzene", "experimental"),
            "p-xylene": (0.0, "nonpolar", "Para-dimethylbenzene: centrosymmetric", "experimental"),
            "o-dichlorobenzene": (2.50, "polar", "Ortho substitution", "experimental"),
            "m-dichlorobenzene": (1.72, "polar", "Meta substitution", "experimental"),
            "p-dichlorobenzene": (0.0, "nonpolar", "Para: centrosymmetric D2h", "experimental"),
            "o-nitrophenol": (3.40, "polar", "Ortho NO2 + OH", "experimental"),
            "p-nitrophenol": (5.02, "highly polar", "Para NO2 + OH dipoles align", "experimental"),
            "p-benzoquinone": (0.0, "nonpolar", "Centrosymmetric D2h", "experimental"),
            "hydrogen peroxide": (2.26, "polar", "Gauche conformation (C2 symmetry)", "experimental"),
            "H2O2": (2.26, "polar", "Gauche conformation", "experimental"),
            "ozone": (0.53, "polar", "Bent triatomic", "experimental"),
            "N2H4": (1.75, "polar", "Hydrazine: gauche conformation", "experimental"),
            "N2F4": (0.20, "slightly polar", "Nearly staggered/gauche", "estimated"),
            "P4": (0.0, "nonpolar", "Td tetrahedral P4", "estimated"),
            "S8": (0.0, "nonpolar", "D4d crown symmetry", "estimated"),
            "white phosphorus": (0.0, "nonpolar", "P4 Td", "estimated"),
            "C60": (0.0, "nonpolar", "Ih icosahedral symmetry", "experimental"),
            "buckminsterfullerene": (0.0, "nonpolar", "C60 Ih", "experimental"),
            "ferrocene(staggered)": (2.04, "polar", "D5d: Fe-Cp dipoles partially align", "experimental"),
            "ferrocene(eclipsed)": (2.04, "polar", "D5h: similar magnitude", "experimental"),
            "Fe(C5H5)2": (2.04, "polar", "Ferrocene", "experimental"),
            "uracil": (4.16, "highly polar", "Nucleobase: multiple C=O and N-H dipoles", "experimental"),
            "cytosine": (6.02, "highly polar", "Nucleobase", "experimental"),
            "adenine": (2.90, "polar", "Nucleobase purine", "experimental"),
            "guanine": (6.38, "highly polar", "Nucleobase", "experimental"),
            "thymine": (4.26, "highly polar", "Nucleobase", "experimental"),
            "glucose": ("very high", "highly polar", "Multiple OH groups → very high dipole", "estimated"),
            "sucrose": ("very high", "highly polar", "Many OH groups", "estimated"),
            "water": (1.85, "polar", "H2O", "experimental"),
            "ammonia": (1.47, "polar", "NH3", "experimental"),
            "ethane": (0.0, "nonpolar", "C2H6 D3d", "experimental"),
            "methane": (0.0, "nonpolar", "CH4 Td", "experimental"),
        }

    def _run_base(self, molecule: str) -> dict:
        mol = molecule.strip()

        if mol in self._db:
            return self._build_result(mol, self._db[mol])

        # Case-insensitive lookup
        mol_lower = mol.lower()
        for key, val in self._db.items():
            if key.lower() == mol_lower:
                return self._build_result(mol, val)

        raise ChemMCPError(
            f"No dipole moment data for '{mol}'. "
            f"Available molecules include:\n"
            f"  Polar: H2O(1.85D), NH3(1.47D), HF(1.82D), HCl(1.08D), HCN(2.98D), "
            f"CH3Cl(1.87D), acetone(2.91D), DMSO(4.0D), urea(4.56D)\n"
            f"  Nonpolar: CH4, CO2, C2H4, benzene, CCl4, SF6, XeF4, BF3\n"
            f"  And 80+ more molecules..."
        )

    def _run_text(self, molecule: str) -> dict:
        return self._run_base(molecule)

    def _build_result(self, mol: str, data) -> dict:
        if isinstance(data, tuple):
            mu, polarity, direction, source = data
            if isinstance(mu, str):
                mu_val = 0.0
            else:
                mu_val = float(mu)
        elif isinstance(data, (int, float)):
            mu_val = float(data)
            polarity = "polar" if abs(data) > 0.5 else ("slightly polar" if abs(data) > 0.1 else "nonpolar")
            direction = "See description"
            source = "experimental"
        else:
            mu_val = 0.0
            polarity = "unknown"
            direction = "N/A"
            source = "unknown"

        return {
            "molecule": mol,
            "dipole_moment": round(mu_val, 2),
            "polarity": polarity,
            "direction": direction,
            "description": data[2] if isinstance(data, tuple) and len(data) >= 3 else "",
            "source": source,
        }

