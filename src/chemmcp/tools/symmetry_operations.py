import logging
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SymmetryOperations(BaseTool):
    """
    对称操作演示与特征标表查询工具 (MCP #292)。
    输入点群符号，返回该点群的完整对称操作列表、特征标表（含所有不可约表示的χ值）、
    以及各不可约表示的对称性说明。
    覆盖常见点群：C1, Cs, Ci, C2, C2v, C3v, C4v, D3, D3d, D3h, D4h, D6h, Td, Oh, D∞h, C∞v 等
    """
    __version__ = "0.1.0"
    name = "SymmetryOperations"
    func_name = "query_symmetry_operations"
    description = "Query symmetry operations and character table for a given Schoenflies point group."
    implementation_description = (
        "Contains built-in character tables for 20+ common point groups. "
        "Returns full symmetry operation list, character table with all irreducible representations, "
        "and interpretation of each irrep's symmetry behavior."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Symmetry", "Character Table", "Group Theory", "Point Group"]
    required_envs = []

    code_input_sig = [
        ("point_group", "str", "N/A", "Schoenflies point group symbol (e.g., 'C2v', 'D3h', 'Td', 'Oh', 'D∞h')."),
    ]

    text_input_sig = [
        ("point_group", "str", "N/A", "Schoenflies point group symbol."),
    ]

    output_sig = [
        ("point_group", "str", "The queried point group."),
        ("order", "int", "Order of the group (number of symmetry operations)."),
        ("symmetry_operations", "list", "List of all symmetry operations with descriptions."),
        ("character_table", "dict", "Full character table: {irrep: {class: χ_value}}."),
        ("irreps_info", "dict", "Description of each irreducible representation's symmetry properties."),
        ("description", "str", "Human-readable description of the point group."),
    ]

    examples = [{'code_input': {'point_group': 'C2v'}, 'text_input': {'point_group': 'C2v'}, 'output': {'point_group': 'C2v', 'order': 4, 'symmetry_operations': ['E', 'C2(z)', 'σv(xz)', "σv'(yz)"], 'character_table': {'A1': {'E': 1, 'C2(z)': 1, 'σv(xz)': 1, "σv'(yz)": 1}}, 'description': 'C2v is common for bent molecules like H2O.', 'irreps_info': 'N/A'}}, {'code_input': {'point_group': 'Td'}, 'text_input': {'point_group': 'Td'}, 'output': {'point_group': 'Td', 'order': 24, 'symmetry_operations': ['E', '4C3', '3C2', '6S4', '6σd'], 'description': 'Td is the tetrahedral group for molecules like CH4.', 'character_table': 'N/A', 'irreps_info': 'N/A'}}]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build character table database."""
        self._char_tables = {
            # ===== C1 =====
            "C1": {
                "order": 1,
                "classes": ["E"],
                "operations": [("E", "Identity")],
                "table": {
                    "A": {"E": 1},
                },
                "irreps_info": {
                    "A": "Totally symmetric. All functions transform as A in C1.",
                },
                "description": "C1 is the trivial group with only identity. Asymmetric molecules like CHClFBr belong here.",
            },

            # ===== Cs =====
            "Cs": {
                "order": 2,
                "classes": ["E", "σ"],
                "operations": [("E", "Identity"), ("σ", "Reflection in mirror plane")],
                "table": {
                    "A'": {"E": 1, "σ": 1},
                    "A''": {"E": 1, "σ": -1},
                },
                "irreps_info": {
                    "A'": "Symmetric under reflection (in-plane vibrations: x, y, Rx², Ry², Rz², xy).",
                    "A''": "Antisymmetric under reflection (out-of-plane: z, xz, yz; Rz).",
                },
                "description": "Cs has only identity and one mirror plane. Molecules like CH3Cl belong to Cs.",
            },

            # ===== Ci =====
            "Ci": {
                "order": 2,
                "classes": ["E", "i"],
                "operations": [("E", "Identity"), ("i", "Inversion through center")],
                "table": {
                    "Ag": {"E": 1, "i": 1},
                    "Au": {"E": 1, "i": -1},
                },
                "irreps_info": {
                    "Ag": "Symmetric under inversion (gerade: s, d orbitals; x², y², z², xy).",
                    "Au": "Antisymmetric under inversion (ungerade: p, f orbitals; x, y, z).",
                },
                "description": "Ci has only identity and inversion center. meso-Tartaric acid belongs to Ci.",
            },

            # ===== C2 =====
            "C2": {
                "order": 2,
                "classes": ["E", "C2"],
                "operations": [("E", "Identity"), ("C2", "180° rotation about principal axis")],
                "table": {
                    "A": {"E": 1, "C2": 1},
                    "B": {"E": 1, "C2": -1},
                },
                "irreps_info": {
                    "A": "Symmetric under C2 rotation (z, Rz, x², y², z², xy).",
                    "B": "Antisymmetric under C2 rotation (x, y, Rx, Ry, xz, yz).",
                },
                "description": "C2 has only a 2-fold axis. Twisted H2O2 (cis) belongs to C2.",
            },

            # ===== C2v =====
            "C2v": {
                "order": 4,
                "classes": ["E", "C2(z)", "σv(xz)", "σv'(yz)"],
                "operations": [
                    ("E", "Identity"),
                    ("C2(z)", "180° rotation about z-axis"),
                    ("σv(xz)", "Reflection in xz plane"),
                    ("σv'(yz)", "Reflection in yz plane"),
                ],
                "table": {
                    "A1": {"E": 1, "C2(z)": 1, "σv(xz)": 1, "σv'(yz)": 1},
                    "A2": {"E": 1, "C2(z)": 1, "σv(xz)": -1, "σv'(yz)": -1},
                    "B1": {"E": 1, "C2(z)": -1, "σv(xz)": 1, "σv'(yz)": -1},
                    "B2": {"E": 1, "C2(z)": -1, "σv(xz)": -1, "σv'(yz)": 1},
                },
                "irreps_info": {
                    "A1": "Totally symmetric (z, x², y², z²). IR + Raman active.",
                    "A2": "Symmetric to C2, antisymmetric to σv (Rz, xy). Raman only.",
                    "B1": "Antisymmetric to C2, symmetric to σv(xz) (x, Rx, xz). IR + Raman active.",
                    "B2": "Antisymmetric to C2, symmetric to σv'(yz) (y, Ry, yz). IR + Raman active.",
                },
                "description": "C2v: Very common! Bent molecules (H2O, H2S, SO2, NO2), CH2Cl2, cis-complexes. All 4 irreps are 1D.",
            },

            # ===== C3v =====
            "C3v": {
                "order": 6,
                "classes": ["E", "2C3", "3σv"],
                "operations": [
                    ("E", "Identity"),
                    ("2C3", "120° and 240° rotations about z-axis"),
                    ("3σv", "Three vertical mirror planes containing z-axis"),
                ],
                "table": {
                    "A1": {"E": 1, "2C3": 1, "3σv": 1},
                    "A2": {"E": 1, "2C3": 1, "3σv": -1},
                    "E": {"E": 2, "2C3": -1, "3σv": 0},
                },
                "irreps_info": {
                    "A1": "Totally symmetric (z, x²+y², z²). IR + Raman active.",
                    "A2": "Symmetric to rotation, antisymmetric to σv (Rz). Inactive.",
                    "E": "2D representation (x,y), (Rx,Ry), (x²-y²,xy), (xz,yz). IR + Raman active (doubly degenerate).",
                },
                "description": "C3v: Trigonal pyramidal molecules (NH3, PH3, PCl3, CH3Cl). Has 2D E irrep.",
            },

            # ===== C4v =====
            "C4v": {
                "order": 8,
                "classes": ["E", "C4", "C4³", "2C2", "4σv"],
                "operations": [
                    ("E", "Identity"),
                    ("C4", "90° rotation about z-axis"),
                    ("C4³", "270° rotation about z-axis"),
                    ("2C2", "Two 180° rotations about axes ⊥ z"),
                    ("4σv", "Four vertical mirror planes"),
                ],
                "table": {
                    "A1": {"E": 1, "C4": 1, "C4³": 1, "2C2": 1, "4σv": 1},
                    "A2": {"E": 1, "C4": 1, "C4³": 1, "2C2": 1, "4σv": -1},
                    "B1": {"E": 1, "C4": -1, "C4³": -1, "2C2": 1, "4σv": 1},
                    "B2": {"E": 1, "C4": -1, "C4³": -1, "2C2": -1, "4σv": 1},
                    "E": {"E": 2, "C4": 0, "C4³": 0, "2C2": -2, "4σv": 0},
                },
                "irreps_info": {
                    "A1": "Totally symmetric (z, x²+y², z²). IR + Raman.",
                    "A2": "(Rz). Inactive.",
                    "B1": "(x²-y²). Raman only.",
                    "B2": "(xy). Raman only.",
                    "E": "2D: (x,y), (Rx,Ry), (xz,yz). IR + Raman (doubly degenerate).",
                },
                "description": "C4v: Square pyramidal (BrF5, XeOF4, IF5) or 4-coordinate planar with axial ligands.",
            },

            # ===== C2h =====
            "C2h": {
                "order": 4,
                "classes": ["E", "C2", "i", "σh"],
                "operations": [
                    ("E", "Identity"),
                    ("C2", "180° rotation about principal axis"),
                    ("i", "Inversion through center"),
                    ("σh", "Horizontal mirror plane"),
                ],
                "table": {
                    "Ag": {"E": 1, "C2": 1, "i": 1, "σh": 1},
                    "Bg": {"E": 1, "C2": -1, "i": 1, "σh": -1},
                    "Au": {"E": 1, "C2": 1, "i": -1, "σh": -1},
                    "Bu": {"E": 1, "C2": -1, "i": -1, "σh": 1},
                },
                "irreps_info": {
                    "Ag": "g+symmetric (Rz, x², y², z², xy). Raman only.",
                    "Bg": "g+antisymmetric (Rx, Ry, xz, yz). Raman only.",
                    "Au": "u+symmetric (z). IR only.",
                    "Bu": "u+antisymmetric (x, y). IR only.",
                },
                "description": "C2h: Planar trans-conformations (trans-HOOH, planar N2F2). g/u labels from inversion.",
            },

            # ===== C3h =====
            "C3h": {
                "order": 6,
                "classes": ["E", "C3", "C3²", "σh", "2S3"],
                "operations": [
                    ("E", "Identity"), ("C3", "120° rotation"), ("C3²", "240° rotation"),
                    ("σh", "Horizontal mirror"), ("2S3", "Improper rotations S3⁺/S3⁻"),
                ],
                "table": {
                    "A'": {"E": 1, "C3": 1, "C3²": 1, "σh": 1, "2S3": 1},
                    "E'": {"E": 2, "C3": -1, "C3²": -1, "σh": 2, "2S3": -1},
                    "A''": {"E": 1, "C3": 1, "C3²": 1, "σh": -1, "2S3": -1},
                    "E''": {"E": 2, "C3": -1, "C3²": -1, "σh": -2, "2S3": 1},
                },
                "irreps_info": {
                    "A'": "Totally symmetric (x²+y², z²). Raman.",
                    "E'": "2D: (x,y), (x²-y², xy). IR + Raman.",
                    "A''": "(z). IR only.",
                    "E''": "2D: (Rx,Ry), (xz, yz). Raman.",
                },
                "description": "C3h: B(OH)3 planar, certain metal complexes with horizontal mirror.",
            },

            # ===== D3 =====
            "D3": {
                "order": 6,
                "classes": ["E", "2C3", "3C2"],
                "operations": [
                    ("E", "Identity"),
                    ("2C3", "120° and 240° rotations about principal axis"),
                    ("3C2", "Three 2-fold axes perpendicular to principal axis"),
                ],
                "table": {
                    "A1": {"E": 1, "2C3": 1, "3C2": 1},
                    "A2": {"E": 1, "2C3": 1, "3C2": -1},
                    "E": {"E": 2, "2C3": -1, "3C2": 0},
                },
                "irreps_info": {
                    "A1": "Totally symmetric (x²+y², z²). Raman only.",
                    "A2": "(Rz). Inactive.",
                    "E": "2D: (x,y), (Rx,Ry), (x²-y²,xy), (xz,yz). IR + Raman.",
                },
                "description": "D3: Chiral tris-chelate complexes (Co(en)₃³⁺ twisted). No mirror planes → optically active.",
            },

            # ===== D3d =====
            "D3d": {
                "order": 12,
                "classes": ["E", "2C3", "3C2", "i", "2S6", "3σd"],
                "operations": [
                    ("E", "Identity"), ("2C3", "C3±120°"), ("3C2", "3 C2⊥C3"),
                    ("i", "Inversion"), ("2S6", "Improper S6±60°"), ("3σd", "3 dihedral mirrors"),
                ],
                "table": {
                    "A1g": {"E": 1, "2C3": 1, "3C2": 1, "i": 1, "2S6": 1, "3σd": 1},
                    "A2g": {"E": 1, "2C3": 1, "3C2": -1, "i": 1, "2S6": 1, "3σd": -1},
                    "Eg": {"E": 2, "2C3": -1, "3C2": 0, "i": 2, "2S6": -1, "3σd": 0},
                    "A1u": {"E": 1, "2C3": 1, "3C2": 1, "i": -1, "2S6": -1, "3σd": -1},
                    "A2u": {"E": 1, "2C3": 1, "3C2": -1, "i": -1, "2S6": -1, "3σd": 1},
                    "Eu": {"E": 2, "2C3": -1, "3C2": 0, "i": -2, "2S6": 1, "3σd": 0},
                },
                "irreps_info": {
                    "A1g": "Totally symmetric g (x²+y², z²). Raman only.",
                    "A2g": "(Rz). Raman only.",
                    "Eg": "2D g: (xz,yz), (x²-y²,xy). Raman only (doubly degenerate).",
                    "A1u": "u-symmetric. Inactive.",
                    "A2u": "(z). IR only.",
                    "Eu": "2D u: (x,y), (Rx,Ry). IR only (doubly degenerate).",
                },
                "description": "D3d: Staggered ethane (CH3-CH3), metal acetylacetonates. Anti-conformation of ethane.",
            },

            # ===== D3h =====
            "D3h": {
                "order": 12,
                "classes": ["E", "C3", "3C2'", "σh", "2S3", "3σv"],
                "operations": [
                    ("E", "Identity"), ("C3", "120° rotation about z"),
                    ("3C2'", "Three C2 axes ⊥ to z in molecular plane"),
                    ("σh", "Horizontal mirror plane (molecular plane)"),
                    ("2S3", "Improper rotations S3±60°"),
                    ('3σv', 'Three vertical mirrors containing C3'),
                ],
                "table": {
                    "A1'": {"E": 1, "C3": 1, "3C2'": 1, "σh": 1, "2S3": 1, "3σv": 1},
                    "A2'": {"E": 1, "C3": 1, "3C2'": -1, "σh": 1, "2S3": 1, "3σv": -1},
                    "E'": {"E": 2, "C3": -1, "3C2'": 0, "σh": 2, "2S3": -1, "3σv": 1},
                    "A1''": {"E": 1, "C3": 1, "3C2'": 1, "σh": -1, "2S3": -1, "3σv": -1},
                    "A2''": {"E": 1, "C3": 1, "3C2'": -1, "σh": -1, "2S3": -1, "3σv": 1},
                    "E''": {"E": 2, "C3": -1, "3C2'": 0, "σh": -2, "2S3": 1, "3σv": 0},
                },
                "irreps_info": {
                    "A1'": "Totally symmetric (x²+y², z²). Raman only.",
                    "A2'": "(Rz). Inactive.",
                    "E'": "2D: (x,y), (x²-y²,xy). IR + Raman (doubly degenerate).",
                    "A1''": "Inactive.",
                    "A2''": "(z). IR only.",
                    "E''": "2D: (Rx,Ry), (xz,yz). Raman only (doubly degenerate).",
                },
                "description": "D3h: Trigonal planar (BF3, SO3, NO3⁻, CO3²⁻), trigonal bipyramidal (PCl5). Very important!",
            },

            # ===== D2h =====
            "D2h": {
                "order": 8,
                "classes": ["E", "C2(z)", "C2(y)", "C2(x)", "i", "σ(xy)", "σ(xz)", "σ(yz)"],
                "operations": [
                    ("E", "Identity"), ("C2(z)", "C2 about z"), ("C2(y)", "C2 about y"),
                    ("C2(x)", "C2 about x"), ("i", "Inversion"), ("σ(xy)", "Horizontal mirror"),
                    ("σ(xz)", "xz mirror"), ("σ(yz)", "yz mirror"),
                ],
                "table": {
                    "Ag": {"E": 1, "C2(z)": 1, "C2(y)": 1, "C2(x)": 1, "i": 1, "σ(xy)": 1, "σ(xz)": 1, "σ(yz)": 1},
                    "B1g": {"E": 1, "C2(z)": 1, "C2(y)": -1, "C2(x)": -1, "i": 1, "σ(xy)": 1, "σ(xz)": 1, "σ(yz)": -1},
                    "B2g": {"E": 1, "C2(z)": -1, "C2(y)": 1, "C2(x)": -1, "i": 1, "σ(xy)": 1, "σ(xz)": -1, "σ(yz)": 1},
                    "B3g": {"E": 1, "C2(z)": -1, "C2(y)": -1, "C2(x)": 1, "i": 1, "σ(xy)": 1, "σ(xz)": -1, "σ(yz)": -1},
                    "Au": {"E": 1, "C2(z)": 1, "C2(y)": 1, "C2(x)": 1, "i": -1, "σ(xy)": -1, "σ(xz)": -1, "σ(yz)": -1},
                    "B1u": {"E": 1, "C2(z)": 1, "C2(y)": -1, "C2(x)": -1, "i": -1, "σ(xy)": -1, "σ(xz)": -1, "σ(yz)": 1},
                    "B2u": {"E": 1, "C2(z)": -1, "C2(y)": 1, "C2(x)": -1, "i": -1, "σ(xy)": -1, "σ(xz)": 1, "σ(yz)": -1},
                    "B3u": {"E": 1, "C2(z)": -1, "C2(y)": -1, "C2(x)": 1, "i": -1, "σ(xy)": -1, "σ(xz)": 1, "σ(yz)": 1},
                },
                "irreps_info": {
                    "Ag": "Totally symmetric g (x², y², z², xy). Raman only.",
                    "B1g": "(Rz, xz). Raman only.",
                    "B2g": "(Rx, yz). Raman only.",
                    "B3g": "(Ry, xy). Raman only.",
                    "Au": "Inactive.",
                    "B1u": "IR inactive (no translation component).",
                    "B2u": "IR inactive.",
                    "B3u": "(z). IR only.",
                },
                "description": "D2h: Ethylene (C2H4), N2O4(planar), B2H6. All 8 irreps are 1D. Highest abelian group.",
            },

            # ===== D2d =====
            "D2d": {
                "order": 8,
                "classes": ["E", "S4", "C2(z)", "2C2'", "2σd"],
                "operations": [
                    ("E", "Identity"), ("S4", "90° improper rotation (S4/S4³)"),
                    ("C2(z)", "180° rotation about z"), ("2C2'", "Two C2 axes in xy plane"),
                    ("2σd", "Two dihedral mirror planes"),
                ],
                "table": {
                    "A1": {"E": 1, "S4": 1, "C2(z)": 1, "2C2'": 1, "2σd": 1},
                    "A2": {"E": 1, "S4": 1, "C2(z)": 1, "2C2'": -1, "2σd": -1},
                    "B1": {"E": 1, "S4": -1, "C2(z)": 1, "2C2'": 1, "2σd": -1},
                    "B2": {"E": 1, "S4": -1, "C2(z)": 1, "2C2'": -1, "2σd": 1},
                    "E": {"E": 2, "S4": 0, "C2(z)": -2, "2C2'": 0, "2σd": 0},
                },
                "irreps_info": {
                    "A1": "Totally symmetric (x²+y², z²). Raman only.",
                    "A2": "(Rz). Inactive.",
                    "B1": "(x²-y²). Raman only.",
                    "B2": "(z). IR only.",
                    "E": "2D: (x,y), (Rx,Ry), (xz,yz). IR + Raman (doubly degenerate).",
                },
                "description": "D2d: Allene (H2C=C=CH2), puckered cyclobutane, staggered conformation of some molecules.",
            },

            # ===== D4h =====
            "D4h": {
                "order": 16,
                "classes": ["E", "C4", "C4³", "C2", "2C2'", "2C2''", "i", "S4", "σh", "2σv", "2σd"],
                "operations": [
                    ("E", "Identity"), ("C4", "90° rotation"), ("C4³", "270° rotation"),
                    ("C2", "180° rotation about z"), ("2C2'", "C2 through atoms"),
                    ("2C2''", "C2 between atoms"), ("i", "Inversion"),
                    ("S4", "Improper rotation"), ("σh", "Horizontal mirror"),
                    ("2σv", "Vertical mirrors through atoms"), ("2σd", "Diagonal mirrors"),
                ],
                "table": {
                    "A1g": {"E": 1, "C4": 1, "C4³": 1, "C2": 1, "2C2'": 1, "2C2''": 1, "i": 1, "S4": 1, "σh": 1, "2σv": 1, "2σd": 1},
                    "A2g": {"E": 1, "C4": 1, "C4³": 1, "C2": 1, "2C2'": -1, "2C2''": -1, "i": 1, "S4": 1, "σh": 1, "2σv": -1, "2σd": -1},
                    "B1g": {"E": 1, "C4": -1, "C4³": -1, "C2": 1, "2C2'": 1, "2C2''": -1, "i": 1, "S4": -1, "σh": 1, "2σv": 1, "2σd": -1},
                    "B2g": {"E": 1, "C4": -1, "C4³": -1, "C2": 1, "2C2'": -1, "2C2''": 1, "i": 1, "S4": -1, "σh": 1, "2σv": -1, "2σd": 1},
                    "Eg": {"E": 2, "C4": 0, "C4³": 0, "C2": -2, "2C2'": 0, "2C2''": 0, "i": 2, "S4": 0, "σh": 2, "2σv": 0, "2σd": 0},
                    "A1u": {"E": 1, "C4": 1, "C4³": 1, "C2": 1, "2C2'": 1, "2C2''": 1, "i": -1, "S4": -1, "σh": -1, "2σv": -1, "2σd": -1},
                    "A2u": {"E": 1, "C4": 1, "C4³": 1, "C2": 1, "2C2'": -1, "2C2''": -1, "i": -1, "S4": -1, "σh": -1, "2σv": 1, "2σd": 1},
                    "B1u": {"E": 1, "C4": -1, "C4³": -1, "C2": 1, "2C2'": 1, "2C2''": -1, "i": -1, "S4": 1, "σh": -1, "2σv": -1, "2σd": 1},
                    "B2u": {"E": 1, "C4": -1, "C4³": -1, "C2": 1, "2C2'": -1, "2C2''": 1, "i": -1, "S4": 1, "σh": -1, "2σv": 1, "2σd": -1},
                    "Eu": {"E": 2, "C4": 0, "C4³": 0, "C2": -2, "2C2'": 0, "2C2''": 0, "i": -2, "S4": 0, "σh": -2, "2σv": 0, "2σd": 0},
                },
                "irreps_info": {
                    "A1g": "Totally symmetric g (x²+y², z²). Raman only.",
                    "A2g": "(Rz). Raman only.",
                    "B1g": "(x²-y²). Raman only.",
                    "B2g": "(xy). Raman only.",
                    "Eg": "2D g: (xz, yz), (Rx, Ry). Raman only (doubly degenerate).",
                    "A1u-A2u-B1u-B2u": "Mostly inactive except:",
                    "A2u": "(z). IR only.",
                    "Eu": "2D u: (x, y). IR only (doubly degenerate).",
                },
                "description": "D4h: Square planar complexes (XeF4, PtCl4²⁻, Ni(CN)4²⁻). One of the most important groups in coordination chemistry!",
            },

            # ===== D5d =====
            "D5d": {
                "order": 20,
                "classes": ["E", "2C5", "2C5²", "5C2", "i", "2S10", "2S10³", "5σd"],
                "operations": [
                    ("E", "Identity"), ("2C5", "C5±72°"), ("2C5²", "C5±144°"),
                    ("5C2", "Five C2⊥C5"), ("i", "Inversion"),
                    ("2S10", "S10±36°"), ("2S10³", "S10±108°"), ("5σd", "5 dihedral mirrors"),
                ],
                "table": {
                    "A1g": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": 1, "i": 1, "2S10": 1, "2S10³": 1, "5σd": 1},
                    "A2g": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": -1, "i": 1, "2S10": 1, "2S10³": 1, "5σd": -1},
                    "E1g": {"E": 2, "2C5": 0.618, "2C5²": -1.618, "5C2": 0, "i": 2, "2S10": 0.618, "2S10³": -1.618, "5σd": 0},
                    "E2g": {"E": 2, "2C5": -1.618, "2C5²": 0.618, "5C2": 0, "i": 2, "2S10": -1.618, "2S10³": 0.618, "5σd": 0},
                    "A1u": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": 1, "i": -1, "2S10": -1, "2S10³": -1, "5σd": -1},
                    "A2u": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": -1, "i": -1, "2S10": -1, "2S10³": -1, "5σd": 1},
                    "E1u": {"E": 2, "2C5": 0.618, "2C5²": -1.618, "5C2": 0, "i": -2, "2S10": -0.618, "2S10³": 1.618, "5σd": 0},
                    "E2u": {"E": 2, "2C5": -1.618, "2C5²": 0.618, "5C2": 0, "i": -2, "2S10": 1.618, "2S10³": -0.618, "5σd": 0},
                },
                "irreps_info": {
                    "A1g": "Totally symmetric g (x²+y², z²). Raman only.",
                    "A2g": "(Rz). Raman only.",
                    "E1g-E2g": "2D g representations. Raman only (doubly degenerate).",
                    "A1u": "Inactive.",
                    "A2u": "(z). IR only.",
                    "E1u-E2u": "2D u representations. IR only (doubly degenerate).",
                },
                "description": "D5d: Staggered ferrocene (Fe(C5H5)2 staggered). Important in organometallic chemistry.",
            },

            # ===== D5h =====
            "D5h": {
                "order": 20,
                "classes": ["E", "2C5", "2C5²", "5C2", "σh", "2S5", "2S5³", "5σv"],
                "operations": [
                    ("E", "Identity"), ("2C5", "C5±72°"), ("2C5²", "C5±144°"),
                    ("5C2", "Five C2⊥C5"), ("σh", "Horizontal mirror"),
                    ("2S5", "S5±36°"), ("2S5³", "S5±108°"), ("5σv", "5 vertical mirrors"),
                ],
                "table": {
                    "A1'": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": 1, "σh": 1, "2S5": 1, "2S5³": 1, "5σv": 1},
                    "A2'": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": -1, "σh": 1, "2S5": 1, "2S5³": 1, "5σv": -1},
                    "E1'": {"E": 2, "2C5": 0.618, "2C5²": -1.618, "5C2": 0, "σh": 2, "2S5": 0.618, "2S5³": -1.618, "5σv": 1},
                    "E2'": {"E": 2, "2C5": -1.618, "2C5²": 0.618, "5C2": 0, "σh": 2, "2S5": -1.618, "2S5³": 0.618, "5σv": 1},
                    "A1''": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": 1, "σh": -1, "2S5": -1, "2S5³": -1, "5σv": -1},
                    "A2''": {"E": 1, "2C5": 1, "2C5²": 1, "5C2": -1, "σh": -1, "2S5": -1, "2S5³": -1, "5σv": 1},
                    "E1''": {"E": 2, "2C5": 0.618, "2C5²": -1.618, "5C2": 0, "σh": -2, "2S5": -0.618, "2S5³": 1.618, "5σv": 0},
                    "E2''": {"E": 2, "2C5": -1.618, "2C5²": 0.618, "5C2": 0, "σh": -2, "2S5": 1.618, "2S5³": -0.618, "5σv": 0},
                },
                "irreps_info": {
                    "A1'": "Totally symmetric (x²+y², z²). Raman only.",
                    "A2'": "(Rz). Inactive.",
                    "E1'": "2D: (x,y). IR + Raman (doubly degenerate).",
                    "E2'": "2D: (x²-y²,xy). Raman only.",
                    "A1''": "Inactive.",
                    "A2''": "(z). IR only.",
                    "E1''-E2''": "2D''. Raman only.",
                },
                "description": "D5h: Eclipsed ferrocene, pentagonal planar complexes.",
            },

            # ===== D6h =====
            "D6h": {
                "order": 24,
                "classes": ["E", "C6", "C6⁵", "C3=C6²", "C3²=C6⁴", "C2=C6³", "2C2'", "2C2''", "i", "2S3", "2S6", "2S6⁵", "σh", "2σv", "2σd", "2σ'"],
                "operations": [
                    ("E", "Identity"), ("C6", "60° rotation"), ("C6⁵", "300° rotation"),
                    ("C3", "120° rotation"), ("C3²", "240° rotation"), ("C2", "180° rotation"),
                    ("2C2'", "C2 through opposite atoms"), ("2C2''", "C2 between atoms"),
                    ("i", "Inversion"), ("2S3", "S3±60°"), ("2S6", "S6±30°"),
                    ("2S6⁵", "S6⁵"), ("σh", "Horizontal mirror (molecular plane)"),
                    ("2σv", "Mirrors through atoms"), ("2σd", "Diagonal mirrors"),
                    ("2σ'", "Additional mirrors"),
                ],
                "table": {
                    "A1g": {"E": 1, "C6": 1, "C6⁵": 1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": 1, "2C2''": 1, "i": 1, "2S3": 1, "2S6": 1, "2S6⁵": 1, "σh": 1, "2σv": 1, "2σd": 1, "2σ'": 1},
                    "A2g": {"E": 1, "C6": 1, "C6⁵": 1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": -1, "2C2''": -1, "i": 1, "2S3": 1, "2S6": 1, "2S6⁵": 1, "σh": 1, "2σv": -1, "2σd": -1, "2σ'": -1},
                    "B1g": {"E": 1, "C6": -1, "C6⁵": -1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": 1, "2C2''": -1, "i": 1, "2S3": 1, "2S6": -1, "2S6⁵": -1, "σh": 1, "2σv": 1, "2σd": -1, "2σ'": -1},
                    "B2g": {"E": 1, "C6": -1, "C6⁵": -1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": -1, "2C2''": 1, "i": 1, "2S3": 1, "2S6": -1, "2S6⁵": -1, "σh": 1, "2σv": -1, "2σd": 1, "2σ'": -1},
                    "E1g": {"E": 2, "C6": 1, "C6⁵": 1, "C3": -1, "C3²": -1, "C2": -2, "2C2'": 0, "2C2''": 0, "i": 2, "2S3": -1, "2S6": 1, "2S6⁵": 1, "σh": 2, "2σv": 0, "2σd": 0, "2σ'": 0},
                    "E2g": {"E": 2, "C6": -1, "C6⁵": -1, "C3": -1, "C3²": -1, "C2": 2, "2C2'": 0, "2C2''": 0, "i": 2, "2S3": -1, "2S6": -1, "2S6⁵": -1, "σh": 2, "2σv": 0, "2σd": 0, "2σ'": 0},
                    "A1u": {"E": 1, "C6": 1, "C6⁵": 1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": 1, "2C2''": 1, "i": -1, "2S3": -1, "2S6": -1, "2S6⁵": -1, "σh": -1, "2σv": -1, "2σd": -1, "2σ'": -1},
                    "A2u": {"E": 1, "C6": 1, "C6⁵": 1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": -1, "2C2''": -1, "i": -1, "2S3": -1, "2S6": -1, "2S6⁵": -1, "σh": -1, "2σv": 1, "2σd": 1, "2σ'": 1},
                    "B1u": {"E": 1, "C6": -1, "C6⁵": -1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": 1, "2C2''": -1, "i": -1, "2S3": -1, "2S6": 1, "2S6⁵": 1, "σh": -1, "2σv": -1, "2σd": 1, "2σ'": 1},
                    "B2u": {"E": 1, "C6": -1, "C6⁵": -1, "C3": 1, "C3²": 1, "C2": 1, "2C2'": -1, "2C2''": 1, "i": -1, "2S3": -1, "2S6": 1, "2S6⁵": 1, "σh": -1, "2σv": 1, "2σd": -1, "2σ'": 1},
                    "E1u": {"E": 2, "C6": 1, "C6⁵": 1, "C3": -1, "C3²": -1, "C2": -2, "2C2'": 0, "2C2''": 0, "i": -2, "2S3": 1, "2S6": -1, "2S6⁵": -1, "σh": -2, "2σv": 0, "2σd": 0, "2σ'": 0},
                    "E2u": {"E": 2, "C6": -1, "C6⁵": -1, "C3": -1, "C3²": -1, "C2": 2, "2C2'": 0, "2C2''": 0, "i": -2, "2S3": 1, "2S6": 1, "2S6⁵": 1, "σh": -2, "2σv": 0, "2σd": 0, "2σ'": 0},
                },
                "irreps_info": {
                    "A1g": "Totally symmetric g (x²+y², z²). Raman only.",
                    "A2g": "(Rz). Raman only.",
                    "B1g-B2g": "Raman only.",
                    "E1g": "2D g: (xz, yz). Raman only (doubly degenerate).",
                    "E2g": "2D g: (x²-y², xy). Raman only (doubly degenerate).",
                    "A1u-A2u-B1u-B2u": "Mostly inactive.",
                    "E1u": "2D u: (x, y). IR only (doubly degenerate).",
                    "E2u": "2D u. Inactive.",
                },
                "description": "D6h: Benzene (C6H6)! The most important aromatic molecule. 12 classes, 12 irreps.",
            },

            # ===== Td =====
            "Td": {
                "order": 24,
                "classes": ["E", "8C3", "3C2", "6S4", "6σd"],
                "operations": [
                    ("E", "Identity"),
                    ("8C3", "Eight 120° rotations about 4 body diagonals (C3 and C3² each)"),
                    ("3C2", "Three 180° rotations about axes through midpoints of opposite edges"),
                    ("6S4", "Six improper 90° rotations (S4 and S4³ about each of 3 C4-like axes)"),
                    ("6σd", "Six dihedral mirror planes"),
                ],
                "table": {
                    "A1": {"E": 1, "8C3": 1, "3C2": 1, "6S4": 1, "6σd": 1},
                    "A2": {"E": 1, "8C3": 1, "3C2": 1, "6S4": -1, "6σd": -1},
                    "E": {"E": 2, "8C3": -1, "3C2": 2, "6S4": 0, "6σd": 0},
                    "T1": {"E": 3, "8C3": 0, "3C2": -1, "6S4": 1, "6σd": -1},
                    "T2": {"E": 3, "8C3": 0, "3C2": -1, "6S4": -1, "6σd": 1},
                },
                "irreps_info": {
                    "A1": "Totally symmetric (x²+y²+z²). Raman only.",
                    "A2": "Inactive.",
                    "E": "2D: (2z²-x²-y², x²-y²). Raman only (doubly degenerate).",
                    "T1": "3D: (Rx, Ry, Rz). Inactive (magnetic dipole allowed).",
                    "T2": "3D: (x, y, z), (xy, xz, yz). IR + Raman active (triply degenerate).",
                },
                "description": "Td: Tetrahedral! CH4, CCl4, SiH4, P4, CF4. No center of inversion → no mutual exclusion rule. 5 irreps up to 3D.",
            },

            # ===== Oh =====
            "Oh": {
                "order": 48,
                "classes": ["E", "8C3", "6C2", "6C4", "3C2(=C4²)", "i", "6S4", "8S6", "3σh", "6σd"],
                "operations": [
                    ("E", "Identity"),
                    ("8C3", "Eight C3 rotations about body diagonals"),
                    ("6C2", "Six C2 rotations about face-center axes"),
                    ("6C4", "Six C4 rotations about Cartesian axes"),
                    ("3C2(=C4²)", "Three C2 (= C4²) about Cartesian axes"),
                    ("i", "Inversion through center"),
                    ("6S4", "Six S4 improper rotations"),
                    ("8S6", "Eight S6 improper rotations"),
                    ("3σh", "Three horizontal mirrors (xy, xz, yz)"),
                    ("6σd", "Six diagonal/dihedral mirror planes"),
                ],
                "table": {
                    "A1g": {"E": 1, "8C3": 1, "6C2": 1, "6C4": 1, "3C2": 1, "i": 1, "6S4": 1, "8S6": 1, "3σh": 1, "6σd": 1},
                    "A2g": {"E": 1, "8C3": 1, "6C2": -1, "6C4": -1, "3C2": 1, "i": 1, "6S4": -1, "8S6": 1, "3σh": 1, "6σd": -1},
                    "Eg": {"E": 2, "8C3": -1, "6C2": 0, "6C4": 0, "3C2": 2, "i": 2, "6S4": 0, "8S6": -1, "3σh": 2, "6σd": 0},
                    "T1g": {"E": 3, "8C3": 0, "6C2": -1, "6C4": 1, "3C2": -1, "i": 3, "6S4": 1, "8S6": 0, "3σh": -1, "6σd": -1},
                    "T2g": {"E": 3, "8C3": 0, "6C2": 1, "6C4": -1, "3C2": -1, "i": 3, "6S4": -1, "8S6": 0, "3σh": -1, "6σd": 1},
                    "A1u": {"E": 1, "8C3": 1, "6C2": 1, "6C4": 1, "3C2": 1, "i": -1, "6S4": -1, "8S6": -1, "3σh": -1, "6σd": -1},
                    "A2u": {"E": 1, "8C3": 1, "6C2": -1, "6C4": -1, "3C2": 1, "i": -1, "6S4": 1, "8S6": -1, "3σh": -1, "6σd": 1},
                    "Eu": {"E": 2, "8C3": -1, "6C2": 0, "6C4": 0, "3C2": 2, "i": -2, "6S4": 0, "8S6": 1, "3σh": -2, "6σd": 0},
                    "T1u": {"E": 3, "8C3": 0, "6C2": -1, "6C4": 1, "3C2": -1, "i": -3, "6S4": -1, "8S6": 0, "3σh": 1, "6σd": 1},
                    "T2u": {"E": 3, "8C3": 0, "6C2": 1, "6C4": -1, "3C2": -1, "i": -3, "6S4": 1, "8S6": 0, "3σh": 1, "6σd": -1},
                },
                "irreps_info": {
                    "A1g": "Totally symmetric g (x²+y²+z²). Raman only.",
                    "A2g": "Inactive.",
                    "Eg": "2D g: (2z²-x²-y², x²-y²). Raman only (doubly degenerate).",
                    "T1g": "3D g: (Rx, Ry, Rz). Magnetic dipole / Raman (triply degenerate).",
                    "T2g": "3D g: (xy, xz, yz). Raman only (triply degenerate).",
                    "A1u-A2u": "Inactive.",
                    "Eu": "2D u. Inactive.",
                    "T1u": "3D u: (x, y, z). IR active! Translations transform as T1u (triply degenerate).",
                    "T2u": "3D u. Inactive.",
                },
                "description": "Oh: Octahedral! SF6, UF6, Mo(CO)6, [Fe(CN)6]³⁻. The most important coordination geometry. 10 irreps up to 3D.",
            },

            # ===== Ih =====
            "Ih": {
                "order": 120,
                "classes": ["E", "12C5", "12C5²", "20C3", "15C2", "i", "12S10", "12S10³", "20S6", "15σ"],
                "operations": [
                    ("E", "Identity"), ("12C5", "Twelve C5 rotations"), ("12C5²", "Twelve C5² rotations"),
                    ("20C3", "Twenty C3 rotations"), ("15C2", "Fifteen C2 rotations"),
                    ("i", "Inversion"), ("12S10", "Twelve S10"), ("12S10³", "Twelve S10³"),
                    ("20S6", "Twenty S6"), ("15σ", "Fifteen mirror planes"),
                ],
                "table": {
                    "Ag": {"E": 1, "12C5": 1, "12C5²": 1, "20C3": 1, "15C2": 1, "i": 1, "12S10": 1, "12S10³": 1, "20S6": 1, "15σ": 1},
                    "T1g": {"E": 3, "12C5": 0.618*3, "12C5²": -1.618*3, "20C3": 0, "15C2": -1, "i": 3, "12S10": 0.618*3, "12S10³": -1.618*3, "20S6": 0, "15σ": -1},
                    "T2g": {"E": 3, "12C5": -1.618*3, "12C5²": 0.618*3, "20C3": 0, "15C2": -1, "i": 3, "12S10": -1.618*3, "12S10³": 0.618*3, "20S6": 0, "15σ": -1},
                    "Gg": {"E": 4, "12C5": 1, "12C5²": 1, "20C3": 1, "15C2": 0, "i": 4, "12S10": 1, "12S10³": 1, "20S6": 1, "15σ": 0},
                    "Hg": {"E": 5, "12C5": 0, "12C5²": 0, "20C3": -1, "15C2": 1, "i": 5, "12S10": 0, "12S10³": 0, "20S6": -1, "15σ": 1},
                    "Au": {"E": 1, "12C5": 1, "12C5²": 1, "20C3": 1, "15C2": 1, "i": -1, "12S10": -1, "12S10³": -1, "20S6": -1, "15σ": -1},
                    "T1u": {"E": 3, "12C5": 0.618*3, "12C5²": -1.618*3, "20C3": 0, "15C2": -1, "i": -3, "12S10": -0.618*3, "12S10³": 1.618*3, "20S6": 0, "15σ": 1},
                    "T2u": {"E": 3, "12C5": -1.618*3, "12C5²": 0.618*3, "20C3": 0, "15C2": -1, "i": -3, "12S10": 1.618*3, "12S10³": -0.618*3, "20S6": 0, "15σ": 1},
                    "Gu": {"E": 4, "12C5": 1, "12C5²": 1, "20C3": 1, "15C2": 0, "i": -4, "12S10": -1, "12S10³": -1, "20S6": -1, "15σ": 0},
                    "Hu": {"E": 5, "12C5": 0, "12C5²": 0, "20C3": -1, "15C2": 1, "i": -5, "12S10": 0, "12S10³": 0, "20S6": 1, "15σ": -1},
                },
                "irreps_info": {
                    "Ag": "Totally symmetric g. Raman only.",
                    "T1g-T2g-Gg-Hg": "g-type representations. Raman only.",
                    "Au": "Inactive.",
                    "T1u": "3D u: (x, y, z). IR active (translations).",
                    "T2u-Gu-Hu": "u-type. Mostly inactive.",
                },
                "description": "Ih: Icosahedral! C60 (buckminsterfullerene), B12H12²⁻. The highest symmetry point group. Order=120, 10 irreps up to 5D (H)!",
            },

            # ===== D∞h =====
            "D∞h": {
                "order": -1,  # infinite
                "classes": ["E", "C∞+", "C∞-", "∞C2⊥", "i", "S∞+", "∞σv", "σh"],
                "operations": [
                    ("E", "Identity"),
                    ("C∞+", "Rotation by +α about molecular axis"),
                    ("C∞-", "Rotation by -α about molecular axis"),
                    ("∞C2⊥", "Infinite C2 axes perpendicular to molecular axis"),
                    ("i", "Inversion through center"),
                    ("S∞+", "Improper rotation by +α"),
                    ("∞σv", "Infinite vertical mirror planes"),
                    ("σh", "Horizontal mirror plane perpendicular to axis"),
                ],
                "table": {
                    "Σg+": {"E": 1, "C∞+": 1, "C∞-": 1, "∞C2⊥": 1, "i": 1, "S∞+": 1, "∞σv": 1, "σh": 1},
                    "Σg-": {"E": 1, "C∞+": 1, "C∞-": 1, "∞C2⊥": 1, "i": 1, "S∞+": 1, "∞σv": -1, "σh": -1},
                    "Πg": {"E": 2, "C∞+": "2*cos(φ)", "C∞-": "2*cos(-φ)", "∞C2⊥": 0, "i": 2, "S∞+": "-2*cos(φ)", "∞σv": 0, "σh": 0},
                    "Δg": {"E": 2, "C∞+": "2*cos(2φ)", "C∞-": "2*cos(-2φ)", "∞C2⊥": 0, "i": 2, "S∞+": "2*cos(2φ)", "∞σv": 0, "σh": 0},
                    "Σu+": {"E": 1, "C∞+": 1, "C∞-": 1, "∞C2⊥": 1, "i": -1, "S∞+": -1, "∞σv": 1, "σh": -1},
                    "Σu-": {"E": 1, "C∞+": 1, "C∞-": 1, "∞C2⊥": 1, "i": -1, "S∞+": -1, "∞σv": -1, "σh": 1},
                    "Πu": {"E": 2, "C∞+": "2*cos(φ)", "C∞-": "2*cos(-φ)", "∞C2⊥": 0, "i": -2, "S∞+": "2*cos(φ)", "∞σv": 0, "σh": 0},
                    "Δu": {"E": 2, "C∞+": "2*cos(2φ)", "C∞-": "2*cos(-2φ)", "∞C2⊥": 0, "i": -2, "S∞+": "-2*cos(2φ)", "∞σv": 0, "σh": 0},
                },
                "irreps_info": {
                    "Σg+": "Totally symmetric gerade (z²). Raman only.",
                    "Σg-": "Gerade, antisymmetric to σv (Rz). Inactive.",
                    "Πg": "2D g: (Rx, Ry). Raman (magnetic dipole).",
                    "Δg+": "Higher angular momentum g states. Raman.",
                    "Σu+": "Ungerade (z along axis). IR active.",
                    "Σu-": "Ungerade. Inactive.",
                    "Πu": "2D u: (x, y). IR active (bond-bending modes).",
                    "Δu": "Higher u states. Mostly inactive.",
                },
                "description": "D∞h: Linear homonuclear diatomic (H2, O2, N2, Cl2) or centrosymmetric linear (CO2, C2H2, [N+]=[N-]). Infinite order group. g/u labeling from inversion center.",
            },

            # ===== C∞v =====
            "C∞v": {
                "order": -1,
                "classes": ["E", "C∞+", "C∞-", "∞σv"],
                "operations": [
                    ("E", "Identity"),
                    ("C∞+", "Rotation by +α about molecular axis"),
                    ("C∞-", "Rotation by -α about molecular axis"),
                    ("∞σv", "Infinite vertical mirror planes containing axis"),
                ],
                "table": {
                    "Σ+": {"E": 1, "C∞+": 1, "C∞-": 1, "∞σv": 1},
                    "Σ-": {"E": 1, "C∞+": 1, "C∞-": 1, "∞σv": -1},
                    "Π": {"E": 2, "C∞+": "2*cos(φ)", "C∞-": "2*cos(-φ)", "∞σv": 0},
                    "Δ": {"E": 2, "C∞+": "2*cos(2φ)", "C∞-": "2*cos(-2φ)", "∞σv": 0},
                    "Φ": {"E": 2, "C∞+": "2*cos(3φ)", "C∞-": "2*cos(-3φ)", "∞σv": 0},
                },
                "irreps_info": {
                    "Σ+": "Totally symmetric (z, x²+y², z²). IR + Raman active.",
                    "Σ-": "Antisymmetric to σv (Rz). Inactive.",
                    "Π": "2D: (x, y), (Rx, Ry). IR + Raman active (doubly degenerate).",
                    "Δ": "2D higher. Raman active.",
                    "Φ": "2D even higher. Raman active.",
                },
                "description": "C∞v: Linear heteronuclear diatomic (HCl, CO, HF, NO, HCN) or any polar linear molecule. No inversion center → no g/u labels. Infinite order.",
            },
        }

    def _run_base(self, point_group: str) -> dict:
        """Query symmetry operations and character table for a point group."""
        pg = point_group.strip()

        # Normalize common variants
        pg_normalized = self._normalize(pg)

        if pg_normalized not in self._char_tables:
            available = sorted(self._char_tables.keys())
            raise ChemMCPError(
                f"Unknown point group '{pg}'. Available point groups: {available}\n"
                f"Note: Use subscripts like C2v, C3v, D3h, D4h, D6h, Td, Oh, D∞h, C∞v, etc."
            )

        data = self._char_tables[pg_normalized]
        return {
            "point_group": pg_normalized,
            "order": data["order"],
            "symmetry_operations": data["operations"],
            "character_table": data["table"],
            "irreps_info": data["irreps_info"],
            "description": data["description"],
            "num_classes": len(data["classes"]),
            "num_irreps": len(data["table"]),
        }

    def _run_text(self, point_group: str) -> dict:
        return self._run_base(point_group)

    @staticmethod
    def _normalize(pg: str) -> str:
        """Normalize point group name variations."""
        p = pg.strip()
        # Common aliases
        mapping = {
            "c2v": "C2v", "c3v": "C3v", "c4v": "C4v", "c5v": "C5v",
            "c2h": "C2h", "c3h": "C3h",
            "d3": "D3", "d3d": "D3d", "d3h": "D3h", "d2h": "D2h",
            "d2d": "D2d", "d4h": "D4h", "d5d": "D5d", "d5h": "D5h",
            "d6h": "D6h", "d4d": "D4d", "d6d": "D6d",
            "td": "Td", "oh": "Oh", "ih": "Ih",
            "c1": "C1", "cs": "Cs", "ci": "Ci", "c2": "C2", "c3": "C3",
            "dinfh": "D∞h", "d infinity h": "D∞h", "doo_h": "D∞h",
            "cinfv": "C∞v", "c infinity v": "C∞v", "coo_v": "C∞v",
            "dinf": "D∞h", "cinf": "C∞v",
        }
        lower_key = p.lower().replace(" ", "")
        if lower_key in mapping:
            return mapping[lower_key]
        if p in mapping:
            return mapping[p]
        return p

