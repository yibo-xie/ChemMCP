import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

# Valence electrons for elements (group number for main group)
VALENCE_ELECTRONS = {
    "H": 1, "He": 2,
    "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
    "K": 1, "Ca": 2, "Ga": 3, "Ge": 4, "As": 5, "Se": 6, "Br": 7, "Kr": 8,
    "I": 7, "Xe": 8,
}

# Predefined Lewis structures for common molecules
LEWIS_STRUCTURES: dict = {
    # Format: central_atom, bonds, lone_pairs_per_atom, total_valence_e, structure_text
    "H2O": {
        "formula": "H2O", "central": "O", "total_ve": 8,
        "bonds": [("O", "H", 1), ("O", "H", 1)],
        "lone_pairs": {"O": 2},
        "structure": """
    ·· 
  H — O — H
    ·· 
""",
        "description": "Bent molecular geometry (~104.5°), tetrahedral electron geometry.",
        "octet_satisfied": True,
    },
    "CO2": {
        "formula": "CO2", "central": "C", "total_ve": 16,
        "bonds": [("C", "O", 2), ("C", "O", 2)],
        "lone_pairs": {"O": 2, "O_2": 2},
        "structure": """
  ··     ··
  O = C = O
  ··     ··
""",
        "description": "Linear molecular geometry (180°), sp hybridized carbon.",
        "octet_satisfied": True,
    },
    "NH3": {
        "formula": "NH3", "central": "N", "total_ve": 8,
        "bonds": [("N", "H", 1), ("N", "H", 1), ("N", "H", 1)],
        "lone_pairs": {"N": 1},
        "structure": """
      ··
    H — N — H
      |
      H
      ··
""",
        "description": "Trigonal pyramidal molecular geometry (~107°), tetrahedral electron geometry.",
        "octet_satisfied": True,
    },
    "CH4": {
        "formula": "CH4", "central": "C", "total_ve": 8,
        "bonds": [("C", "H", 1), ("C", "H", 1), ("C", "H", 1), ("C", "H", 1)],
        "lone_pairs": {},
        "structure": """
      H
      |
  H — C — H
      |
      H
""",
        "description": "Tetrahedral molecular geometry (109.5°), sp³ hybridized carbon.",
        "octet_satisfied": True,
    },
    "SO2": {
        "formula": "SO2", "central": "S", "total_ve": 18,
        "bonds": [("S", "O", 2), ("S", "O", 1)],
        "lone_pairs": {"S": 1, "O": 2, "O_2": 3},
        "structure": """(Resonance structures exist)
  ··           ··
  O = S — O:
  ··   ··   ··
↔
  ··       ··
:O — S = O
··  ··   ··
""",
        "description": "Bent molecular geometry (~119°), trigonal planar electron geometry. Resonance hybrid of two structures.",
        "octet_satisfied": True,
    },
    "PCl5": {
        "formula": "PCl5", "central": "P", "total_ve": 40,
        "bonds": [("P", "Cl", 1)] * 5,
        "lone_pairs": {},
        "structure": """
        Cl
        |
    Cl — P — Cl
       / \\
      Cl   Cl
""",
        "description": "Trigonal bipyramidal molecular geometry (90° and 120°), sp³d hybridized phosphorus. Expanded octet (10 electrons around P).",
        "octet_satisfied": False,  # expanded octet
    },
    "SF6": {
        "formula": "SF6", "central": "S", "total_ve": 48,
        "bonds": [("S", "F", 1)] * 6,
        "lone_pairs": {},
        "structure": """
      F   F
      |   |
  F — S — F
      |   |
      F   F
""",
        "description": "Octahedral molecular geometry (90°), sp³d² hybridized sulfur. Expanded octet (12 electrons around S).",
        "octet_satisfied": False,
    },
    "BF3": {
        "formula": "BF3", "central": "B", "total_ve": 24,
        "bonds": [("B", "F", 1), ("B", "F", 1), ("B", "F", 1)],
        "lone_pairs": {"F": 3, "F_2": 3, "F_3": 3},
        "structure": """
      ··   ··   ··
    F — B — F
      ··   ··
      F
      ··
""",
        "description": "Trigonal planar molecular geometry (120°), sp² hybridized boron. Boron is electron-deficient (only 6 valence electrons).",
        "octet_satisfied": False,  # electron deficient
    },
    "N2": {
        "formula": "N2", "central": None, "total_ve": 10,
        "bonds": [("N", "N", 3)],
        "lone_pairs": {"N": 1, "N_2": 1},
        "structure": """
  ··     ··
  N ≡ N
  ··     ··
""",
        "description": "Diatomic molecule with triple bond. Linear (180°). Very strong bond (945 kJ/mol).",
        "octet_satisfied": True,
    },
    "O2": {
        "formula": "O2", "central": None, "total_ve": 12,
        "bonds": [("O", "O", 2)],
        "lone_pairs": {"O": 2, "O_2": 2},
        "structure": """
  ··   ··
  O = O
  ··   ··
""",
        "description": "Diatomic with double bond. Paramagnetic (two unpaired electrons in π* antibonding MOs). Bond order = 2.",
        "octet_satisfied": True,
    },
    "HCl": {
        "formula": "HCl", "central": None, "total_ve": 8,
        "bonds": [("H", "Cl", 1)],
        "lone_pairs": {"Cl": 3},
        "structure": """
  ·· ·· ··
  H — Cl:
  ·· ·· ··
""",
        "description": "Diatomic polar molecule. Linear. Dipole moment due to electronegativity difference.",
        "octet_satisfied": True,
    },
    "H2S": {
        "formula": "H2S", "central": "S", "total_ve": 8,
        "bonds": [("S", "H", 1), ("S", "H", 1)],
        "lone_pairs": {"S": 2},
        "structure": """
    ··
  H — S — H
    ··
""",
        "description": "Bent molecular geometry (~92°), more acute than H2O because S has larger atomic size.",
        "octet_satisfied": True,
    },
    "PH3": {
        "formula": "PH3", "central": "P", "total_ve": 8,
        "bonds": [("P", "H", 1), ("P", "H", 1), ("P", "H", 1)],
        "lone_pairs": {"P": 1},
        "structure": """
      ··
    H — P — H
      |
      H
      ··
""",
        "description": "Trigonal pyramidal molecular geometry (~93.5°). Less basic than NH3.",
        "octet_satisfied": True,
    },
    "CCl4": {
        "formula": "CCl4", "central": "C", "total_ve": 32,
        "bonds": [("C", "Cl", 1)] * 4,
        "lone_pairs": {"Cl": 3, "Cl_2": 3, "Cl_3": 3, "Cl_4": 3},
        "structure": """
      Cl
      |
  Cl — C — Cl
      |
      Cl
""",
        "description": "Tetrahedral molecular geometry (109.5°). Nonpolar despite polar bonds (symmetric).",
        "octet_satisfied": True,
    },
    "NO2": {
        "formula": "NO2", "central": "N", "total_ve": 17,
        "bonds": [("N", "O", 2), ("N", "O", 1)],
        "lone_pairs": {"N": 0, "O": 2, "O_2": 3},  # odd-electron species
        "structure": """(Resonance; radical with unpaired electron)
  ··         ··
  O = N — O·
  ··   ··   ··
""",
        "description": "Bent molecular geometry (~134°). Odd-electron species (radical, brown gas). Resonance stabilized.",
        "octet_satisfied": False,  # radical
    },
    "SO3": {
        "formula": "SO3", "central": "S", "total_ve": 24,
        "bonds": [("S", "O", 2), ("S", "O", 2), ("S", "O", 2)],
        "lone_pairs": {"O": 2, "O_2": 2, "O_3": 2},
        "structure": """(Resonance)
     ··   ··   ··
   O = S = O
     ··   ··
     O
     ··
""",
        "description": "Trigonal planar molecular geometry (120°). Resonance hybrid of three equivalent structures.",
        "octet_satisfied": False,  # expanded octet on S
    },
    "CH2O": {
        "formula": "CH2O", "central": "C", "total_ve": 12,
        "bonds": [("C", "O", 2), ("C", "H", 1), ("C", "H", 1)],
        "lone_pairs": {"O": 2},
        "structure": """
    ··
    H
    |
    C = O
    |
    H
    ··
""",
        "description": "Trigonal planar molecular geometry (~120°). Polar molecule (carbonyl group).",
        "octet_satisfied": True,
    },
    "HCN": {
        "formula": "HCN", "central": "C", "total_ve": 10,
        "bonds": [("H", "C", 1), ("C", "N", 3)],
        "lone_pairs": {"N": 1},
        "structure": """
  ··
  H — C ≡ N:
  ··
""",
        "description": "Linear molecular geometry (180°). Weak acid (hydrocyanic acid). Toxic.",
        "octet_satisfied": True,
    },
}


@ChemMCPManager.register_tool
class DrawLewisStructure(BaseTool):
    __version__ = "0.1.0"
    name = "DrawLewisStructure"
    func_name = 'draw_lewis_structure'
    description = "Draw the Lewis structure of a molecule showing valence electrons, bonding pairs, and lone pairs."
    implementation_description = "Uses a rule-based algorithm with predefined Lewis structures for common molecules. Shows electron dot notation, bond orders, lone pairs, octet status, and molecular geometry description."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Lewis Structure", "Valence Electrons", "Molecular Geometry", "Chemical Bonding"]
    required_envs = []

    code_input_sig = [
        ('formula', 'str', 'N/A', 'Molecular formula (e.g., H2O, CO2, NH3, CH4, SO2, PCl5, SF6)'),
    ]
    text_input_sig = [
        ('formula', 'str', 'N/A', 'Molecular formula'),
    ]
    output_sig = [
        ('formula', 'str', 'Molecular formula'),
        ('lewis_diagram', 'str', 'ASCII/text Lewis structure diagram'),
        ('total_valence_electrons', 'int', 'Total count of valence electrons'),
        ('bond_info', 'list', 'List of bonds between atoms'),
        ('lone_pair_info', 'dict', 'Lone pair distribution per atom'),
        ('geometry_description', 'str', 'Molecular geometry description'),
        ('octet_status', 'str', 'Whether octet rule is satisfied'),
    ]
    
        
    examples = [
        {'code_input': {'formula': 'H2O'}, 'text_input': {'formula': 'H2O'}, 'output': {'formula': 'H2O', 'lewis_diagram': '...', 'total_valence_electrons': 8, 'bond_info': '...', 'lone_pair_info': '...', 'geometry_description': 'bent (104.5)', 'octet_status': '...'}},
        {'code_input': {'formula': 'CO2'}, 'text_input': {'formula': 'CO2'}, 'output': {'formula': 'CO2', 'lewis_diagram': '...', 'total_valence_electrons': 16, 'bond_info': '...', 'lone_pair_info': '...', 'geometry_description': 'linear (180)', 'octet_status': '...'}},
        {'code_input': {'formula': 'SF6'}, 'text_input': {'formula': 'SF6'}, 'output': {'formula': 'SF6', 'lewis_diagram': '...', 'total_valence_electrons': 48, 'bond_info': '...', 'lone_pair_info': '...', 'geometry_description': 'octahedral (90)', 'octet_status': '...'}},
    ]
    def _run_base(self, formula: str) -> dict:
        # Normalize formula
        f = formula.strip().upper()
        
        if f not in LEWIS_STRUCTURES:
            available = sorted(list(LEWIS_STRUCTURES.keys()))
            raise ChemMCPInputError(
                f"Lewis structure not available for '{formula}'. "
                f"Available molecules: {available}. "
                f"Note: This tool supports common molecules with predefined structures."
            )
        
        ls = LEWIS_STRUCTURES[f]
        return {
            "formula": ls["formula"],
            "lewis_diagram": ls["structure"].strip(),
            "total_valence_electrons": ls["total_ve"],
            "bond_info": [{"atom1": b[0], "atom2": b[1], "order": b[2]} for b in ls["bonds"]],
            "lone_pair_info": ls["lone_pairs"],
            "geometry_description": ls["description"],
            "octet_status": "satisfied" if ls.get("octet_satisfied") else "not satisfied / exception",
        }


if __name__ == "__main__":
    run_mcp_server()
