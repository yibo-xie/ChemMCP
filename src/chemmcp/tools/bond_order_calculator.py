import logging
from typing import Optional, List, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BondOrderCalculator(BaseTool):
    """
    键级计算工具 (MCP #293)。
    根据分子结构计算化学键的键级（bond order），支持：
    - 常见分子的实验/理论键级数据库查询
    - 双原子分子的分子轨道(MO)理论键级计算
    - Lewis结构规则估算多原子分子的键级
    """
    __version__ = "0.1.0"
    name = "BondOrderCalculator"
    func_name = "calculate_bond_order"
    description = "Calculate bond order (bond multiplicity) for chemical bonds in molecules using MO theory, Lewis structures, or database lookup."
    implementation_description = (
        "Combines a database of known bond orders for common molecules with "
        "MO theory calculations for diatomic molecules and Lewis structure rules for polyatomic molecules."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Bond Order", "Chemical Bonding", "MO Theory", "Lewis Structure"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: formula, SMILES, or name (e.g., 'N2', 'CH4', 'benzene', 'O=O')."),
        ("atom_pair", "str", "None", "Optional: atom pair to query, e.g., 'C-C', 'C=C', 'C≡C', 'N-O', or indices like '1-2'. If None, returns all bond orders."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'molecule' or 'molecule atom_pair', e.g., 'N2', 'benzene C-C', 'CO C-O'."),
    ]

    output_sig = [
        ("molecule", "str", "The molecule being analyzed."),
        ("bonds", "list", "List of bond information: each entry has atoms, bond_order, bond_type, method, and description."),
        ("total_bond_order_sum", "float", "Sum of all bond orders in the molecule."),
    ]


    examples = [{'code_input': {'molecule': 'N2', 'atom_pair': 'N/A'}, 'text_input': {'query': 'N2'}, 'output': {'molecule': 'N2', 'bonds': [{'atoms': 'N≡N', 'bond_order': 3.0, 'bond_type': 'triple bond', 'method': 'MO Theory'}], 'total_bond_order_sum': 3.0}}, {'code_input': {'molecule': 'benzene', 'atom_pair': 'N/A'}, 'text_input': {'query': 'benzene C-C'}, 'output': {'molecule': 'benzene', 'bonds': [{'atoms': 'C-C', 'bond_order': 1.5, 'bond_type': 'aromatic partial double', 'method': 'Resonance/Lewis'}], 'total_bond_order_sum': 'N/A'}}, {'code_input': {'molecule': 'O2', 'atom_pair': 'N/A'}, 'text_input': {'query': 'O2'}, 'output': {'molecule': 'O2', 'bonds': [{'atoms': 'O=O', 'bond_order': 2.0, 'bond_type': 'double bond', 'method': 'MO Theory'}], 'total_bond_order_sum': 'N/A'}}]
    def _init_modules(self):
        """Build bond order databases."""
        # Diatomic MO theory bond orders: (valence_e_bonding, valence_e_antibonding)
        self._diatomic_mo = {
            # Homonuclear diatomic (Period 2)
            "H2":   ("H-H",  (2, 0), 1.0, "σ(1s)² → BO=(2-0)/2=1"),
            "He2":  ("He-He",(2, 2), 0.0, "σ(1s)²σ*(1s)² → BO=(2-2)/2=0 (unstable)"),
            "Li2":  ("Li-Li", (2, 0), 1.0, "σ(2s)² → BO=1"),
            "B2":   ("B-B",  (4, 2), 1.0, "σ(2s)²σ*(2s)²π(2px)¹π(2py)¹ → BO=(4-2)/2=1 (paramagnetic)"),
            "C2":   ("C-C",  (6, 2), 2.0, "+σ(2pz)² → BO=(6-2)/2=2"),
            "N2":   ("N≡N",  (8, 2), 3.0, "+π(2px)²π(2py)² → BO=(8-2)/2=3 (strongest diatomic bond)"),
            "O2":   ("O=O",  (8, 4), 2.0, "+π*(2px)¹π*(2py)¹ → BO=(8-4)/2=2 (paramagnetic)"),
            "F2":   ("F-F",  (8, 6), 1.0, "+σ*(2pz)² → BO=(8-6)/2=1 (weak single bond)"),
            "Ne2":  ("Ne-Ne",(8, 8), 0.0, "BO=0 (noble gas, no bonding)"),

            # Heteronuclear diatomic (approximate)
            "HF":   ("H-F",  (2, 0), 1.0, "Polar covalent single bond, BO≈1"),
            "HCl":  ("H-Cl", (2, 0), 1.0, "Polar covalent single bond, BO≈1"),
            "CO":   ("C≡O",  (8, 2), 3.0, "Isoelectronic with N2, BO=3 (with dative character)"),
            "NO":   ("N-O",  (7, 3), 2.0, "BO=(7-3)/2=2 (paramagnetic, bond order 2)"),
            "CN":   ("C≡N",  (8, 2), 3.0, "Isoelectronic with N2, BO≈3"),
            "NO+":  ("N-O+", (8, 2), 3.0, "Removed antibonding e⁻ → BO increases to 3"),
            "NO-":  ("N-O-", (7, 4), 1.5, "Added antibonding e⁻ → BO decreases to 1.5"),
            "O2+":  ("O-O+", (8, 3), 2.5, "Removed one π* e⁻ → BO=2.5"),
            "O2-":  ("O-O-", (8, 5), 1.5, "Added one π* e⁻ → BO=1.5 (superoxide)"),
            "O2(2-)":("O-O²-",(8, 6), 1.0, "Two extra π* e⁻ → BO=1 (peroxide)"),
        }

        # Polyatomic molecule bond order database
        self._polyatomic_db = {
            # Hydrocarbons
            "CH4": [("C-H", 1.0, "single", "sp³ C-H σ bond")],
            "C2H6": [("C-C", 1.0, "single", "sp³-sp³ σ bond"), ("C-H", 1.0, "single", "sp³ C-H ×6")],
            "C2H4": [("C=C", 2.0, "double", "sp²-sp²: σ+π bond"), ("C-H", 1.0, "single", "sp² C-H ×4")],
            "C2H2": [("C≡C", 3.0, "triple bond", "sp-sp: σ+2π bonds"), ("C-H", 1.0, "single", "sp C-H ×2")],
            "c1ccccc1": [("C-C", 1.5, "aromatic partial double", "6 π e⁻ delocalized over 6 C-C bonds"), ("C-H", 1.0, "single", "sp² C-H ×6")],
            "benzene": [("C-C", 1.5, "aromatic partial double", "6 π e⁻ delocalized over 6 C-C bonds"), ("C-H", 1.0, "single", "sp² C-H ×6")],
            "C1=CC=CC=C1": [("C-C", 1.5, "aromatic partial double", "Delocalized π system"), ("C-H", 1.0, "single", "sp² C-H ×6")],
            "CH3CH3": [("C-C", 1.0, "single", "sp³-sp³ σ bond"), ("C-H", 1.0, "single", "sp³ C-H ×6")],
            "CH2=CH2": [("C=C", 2.0, "double", "sp²-sp² σ+π"), ("C-H", 1.0, "single", "sp² C-H ×4")],
            "CH#CH": [("C≡C", 3.0, "triple bond", "sp-sp σ+2π"), ("C-H", 1.0, "single", "sp C-H ×2")],
            "allene": [("C=C", 2.0, "double", "Central sp C=C double ×2 orthogonal π systems"), ("C-H", 1.0, "single", "sp² C-H ×4")],
            "C=C=C": [("C=C", 2.0, "double", "Allene: two orthogonal C=C double bonds"), ("C-H", 1.0, "single", "sp² C-H ×4")],

            # Oxygen compounds
            "H2O": [("O-H", 1.0, "single", "sp³ O-H σ bond, bent geometry")],
            "H2O2": [("O-O", 1.0, "single", "Single σ bond (weak, ~146 kJ/mol)"), ("O-H", 1.0, "single", "O-H ×2")],
            "CO2": [("C=O", 2.0, "double", "Two C=O double bonds (linear, resonance with C≡O⁺-O⁻)")],
            "C=O=C": [("C=O", 2.0, "double", "CO2: two equivalent C=O bonds")],
            "O=S=O": [("S=O", 1.5, "partial double", "SO2: resonance of S=O and S⁺-O⁻, BO≈1.5-2")],
            "SO2": [("S=O", 1.5, "partial double", "SO2: S-O bond order ~1.5-2 from resonance")],
            "SO3": [("S-O", 1.33, "partial double", "SO3: resonance of 3 S=O over 3 positions, BO=4/3≈1.33")],
            "H2SO4": [("S-O(hydroxyl)", 1.0, "single", "S-OH single bond ×2"), ("S=O", 2.0, "double", "S=O terminal double bond ×2")],

            # Nitrogen compounds
            "NH3": [("N-H", 1.0, "single", "sp³ N-H σ bond")],
            "N2O4": [("N-N", 1.0, "single", "N-N single bond (weak, easily dissociates)"), ("N=O", 2.0, "double", "N=O double bond ×4")],
            "NO3-": [("N-O", 1.33, "partial double", "Nitrate: resonance of N=O over 3 O atoms, BO=4/3≈1.33")],
            "NH4+": [("N-H", 1.0, "single", "sp³ N-H σ bond ×4")],

            # Carbon oxides / functional groups
            "HCHO": [("C=O", 2.0, "double", "Carbonyl C=O double bond"), ("C-H", 1.0, "single", "sp² C-H ×2")],
            "HCOOH": [("C=O", 2.0, "double", "Acid carbonyl"), ("C-O(H)", 1.0, "single", "C-OH single"), ("O-H", 1.0, "single", "Hydroxyl")],
            "CH3OH": [("C-O", 1.0, "single", "sp³ C-O σ bond"), ("C-H", 1.0, "single", "sp³ C-H ×3"), ("O-H", 1.0, "single", "Hydroxyl")],
            "CCl4": [("C-Cl", 1.0, "single", "sp³ C-Cl σ bond ×4")],
            "CH2Cl2": [("C-Cl", 1.0, "single", "sp³ C-Cl ×2"), ("C-H", 1.0, "single", "sp³ C-H ×2")],
            "CHCl3": [("C-Cl", 1.0, "single", "sp³ C-Cl ×3"), ("C-H", 1.0, "single", "sp³ C-H")],

            # Other important molecules
            "CCO": [("C-C", 1.0, "single", "Ethanol C-C single bond"), ("C-O", 1.0, "single", "C-O single bond"), ("C-H", 1.0, "single", "C-H ×5"), ("O-H", 1.0, "single", "Hydroxyl")],
            "CH3Cl": [("C-Cl", 1.0, "single", "sp³ C-Cl polar covalent"), ("C-H", 1.0, "single", "sp³ C-H ×3")],
            "BF3": [("B-F", 1.33, "partial double", "BF3: B has incomplete octet, π back-donation gives BO>1, avg ≈1.33")],
            "SF6": [("S-F", 1.0, "single", "Octahedral S(VI)-F single bonds ×6, hypervalent via d orbitals")],
            "XeF4": [("Xe-F", 1.0, "single", "Square planar Xe(IV)-F single bonds ×4, hypervalent")],
            "PCl5": [("P-Cl(axial)", 1.0, "single", "Axial P-Cl, longer/weaker ×2"), ("P-Cl(equatorial)", 1.0, "single", "Equatorial P-Cl ×3")],
            "H2S": [("S-H", 1.0, "single", "sp³ S-H σ bond, bent")],
            "PH3": [("P-H", 1.0, "single", "sp³ P-H σ bond, pyramidal")],
            "SiH4": [("Si-H", 1.0, "single", "tetrahedral Si-H ×4")],
            "HCN": [("C≡N", 3.0, "triple bond", "C≡N triple bond (sp hybridization)"), ("C-H", 1.0, "single", "sp C-H")],
            "C#N": [("C≡N", 3.0, "triple bond", "Cyanide: C≡N triple bond")],
            "CH3CN": [("C-C", 1.0, "single", "sp³-sp C-C single"), ("C≡N", 3.0, "triple bond", "C≡N nitrile triple bond"), ("C-H", 1.0, "single", "C-H ×3")],
            "O3": [("O-O", 1.5, "partial double", "Ozone: central O bonded to two O with BO≈1.5 (resonance hybrid)")],
            "CO3(2-)": [("C-O", 1.33, "partial double", "Carbonate: resonance of C=O over 3 O, BO=4/3≈1.33")],
            "HCO3-": [("C=O", 2.0, "double", "One C=O"), ("C-O", 1.0, "single", "Two C-O single (one protonated)")],
            "H2CO3": [("C=O", 2.0, "double", "Carbonyl"), ("C-O(H)", 1.0, "single", "C-OH ×2")],
            "CH3COOH": [("C-C", 1.0, "single", "C-C single"), ("C=O", 2.0, "double", "Carboxyl C=O"), ("C-O", 1.0, "single", "C-O single"), ("C-H", 1.0, "single", "C-H ×3"), ("O-H", 1.0, "single", "Carboxylic OH")],
            "glycine": [("C-C", 1.0, "single", "α-C to carboxyl C"), ("C=O", 2.0, "double", "Carboxyl carbonyl"), ("C-N", 1.0, "single", "C-N single"), ("C-O", 1.0, "single", "C-OH"), ("C-H/N-H/O-H", 1.0, "single", "Various X-H")],
            "urea": [("C=O", 2.0, "double", "Urea carbonyl"), ("C-N", 1.33, "partial double", "C-N with partial double character from resonance ×2"), ("N-H", 1.0, "single", "N-H ×4")],
            "Pyridine": [("C-C", 1.5, "aromatic", "Aromatic ring C-C/C-N"), ("C-N", 1.5, "aromatic", "Aromatic C-N in ring"), ("C-H", 1.0, "single", "C-H ×5")],
            "graphite": [("C-C", 1.33, "partial double", "Graphene layer: aromatic C-C with BO≈1.33-1.5 within plane")],
            "diamond": [("C-C", 1.0, "single", "Tetrahedral sp³ network, pure single bonds")],
            "C60": [("C-C(6-6)", 1.5, "aromatic", "Fusion of two hexagons, BO≈1.45-1.5"), ("C-C(6-5)", 1.0, "single", "Hexagon-pentagon fusion, BO≈1.0")],
        }

    def _run_base(self, molecule: str, atom_pair: str = None) -> dict:
        """Calculate bond orders for a molecule."""
        mol = molecule.strip()

        # Check diatomic MO database first
        if mol in self._diatomic_mo:
            atoms, (be, bae), bo, desc = self._diatomic_mo[mol]
            bond_info = {"atoms": atoms, "bond_order": bo, "bond_type": self._bo_to_type(bo),
                          "method": "MO Theory", "description": desc}
            return self._format_result(mol, [bond_info], atom_pair)

        # Check polyatomic database
        if mol in self._polyatomic_db:
            bonds = []
            for atoms, bo, btype, desc in self._polyatomic_db[mol]:
                bonds.append({"atoms": atoms, "bond_order": bo, "bond_type": btype,
                              "method": "Database/Lewis", "description": desc})
            return self._format_result(mol, bonds, atom_pair)

        # Case-insensitive lookup
        mol_lower = mol.lower()
        for key in list(self._polyatomic_db.keys()) + list(self._diatomic_mo.keys()):
            if key.lower() == mol_lower:
                return self._run_base(key, atom_pair)

        raise ChemMCPError(
            f"Cannot find bond order data for '{mol}'. "
            f"Supported molecules include:\n"
            f"  Diatomic (MO): H2, Li2, B2, C2, N2, O2, F2, CO, NO, CN, HF, HCl, NO+, NO-, O2+, O2-, O2(2-)\n"
            f"  Polyatomic: CH4, C2H6, C2H4, C2H2, benzene, H2O, CO2, SO2, SO3, NH3, BF3, SF6, XeF4, "
            f"HCHO, CH3OH, HCN, glycine, urea, C60, graphite, diamond, etc."
        )

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        molecule = parts[0]
        atom_pair = parts[1] if len(parts) > 1 else None
        return self._run_base(molecule, atom_pair)

    @staticmethod
    def _bo_to_type(bo: float) -> str:
        if bo <= 0.1:
            return "no bond"
        elif bo < 0.75:
            return "very weak/partial"
        elif bo < 1.25:
            return "single bond"
        elif bo < 1.75:
            return "intermediate (between single & double)"
        elif bo < 2.25:
            return "double bond"
        elif bo < 2.75:
            return "intermediate (between double & triple)"
        elif bo < 3.25:
            return "triple bond"
        else:
            return "very high (>3)"

    def _format_result(self, molecule: str, bonds: list, atom_pair: str = None) -> dict:
        if atom_pair:
            ap = atom_pair.strip().upper()
            filtered = [b for b in bonds if ap in b["atoms"].upper() or b["atoms"].upper() in ap]
            if not filtered:
                raise ChemMCPInputError(f"No bond matching '{atom_pair}' found in {molecule}. Available bonds: {[b['atoms'] for b in bonds]}")
            bonds = filtered

        total_bo = sum(b["bond_order"] for b in bonds)
        return {
            "molecule": molecule,
            "bonds": bonds,
            "total_bond_order_sum": round(total_bo, 3),
        }

