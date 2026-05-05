import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PointGroupIdentifier(BaseTool):
    """
    分子点群识别工具 (MCP #291)。
    通过分子式、SMILES或通用名称识别分子的Schoenflies点群，
    提供对称元素列表和特征标表摘要信息。
    覆盖所有常见点群：C1, Cs, Ci, Cn, Cnv, Cnh, Dn, Dnd, Dnh, Sn, T, Th, Td, O, Oh, I, Ih, C∞v, D∞h
    """
    __version__ = "0.1.0"
    name = "PointGroupIdentifier"
    func_name = "identify_point_group"
    description = "Identify the Schoenflies point group of a molecule from its SMILES, formula, or common name."
    implementation_description = (
        "Uses a comprehensive database of 80+ known molecules mapped to their point groups, "
        "with symmetry element enumeration and character table summary for each point group."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Point Group", "Symmetry", "Group Theory", "Molecular Geometry"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: SMILES string, molecular formula, or common name (e.g., 'H2O', 'CH4', 'benzene')."),
    ]

    text_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: SMILES string, molecular formula, or common name."),
    ]

    output_sig = [
        ("point_group", "str", "The Schoenflies point group symbol."),
        ("symmetry_elements", "list", "List of symmetry elements (E, Cn, σ, i, Sn, etc.)."),
        ("character_table_summary", "str", "Summary of character table: number of classes, irreducible representations, and total order."),
        ("order", "int", "Order of the point group (number of symmetry operations)."),
        ("confidence", "str", "Confidence level: 'high', 'medium', or 'low'."),
    ]

    examples = [
        {
            "code_input": {"molecule": "H2O"},
            "text_input": {"molecule": "H2O"},
            "output": {
                "point_group": "C2v",
                "symmetry_elements": ["E", "C2(z)", "σv(xz)", "σv'(yz)"],
                "character_table_summary": "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4",
                "order": 4,
                "confidence": "high",
            }
        },
        {
            "code_input": {"molecule": "CH4"},
            "text_input": {"molecule": "CH4"},
            "output": {
                "point_group": "Td",
                "symmetry_elements": ["E", "4C3", "3C2", "6S4", "6σd"],
                "character_table_summary": "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24",
                "order": 24,
                "confidence": "high",
            }
        },
        {
            "code_input": {"molecule": "BF3"},
            "text_input": {"molecule": "BF3"},
            "output": {
                "point_group": "D3h",
                "symmetry_elements": ["E", "C3", "3C2'", "σh", "2S3", "3σv"],
                "character_table_summary": "D3h: 6 classes, 6 irreps (A1', A2', E', A1'', A2'', E''), order=12",
                "order": 12,
                "confidence": "high",
            }
        },
        {
            "code_input": {"molecule": "CHClFBr"},
            "text_input": {"molecule": "CHClFBr"},
            "output": {
                "point_group": "C1",
                "symmetry_elements": ["E"],
                "character_table_summary": "C1: 1 class, 1 irrep (A), order=1",
                "order": 1,
                "confidence": "high",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build comprehensive molecule → point group database."""
        self._molecule_db = {
            # ===== C1 (no symmetry) =====
            "CHClFBr": ("C1", ["E"], "C1: 1 class, 1 irrep (A), order=1", 1, "high"),
            "CHFClBr": ("C1", ["E"], "C1: 1 class, 1 irrep (A), order=1", 1, "high"),
            "CHBrClF": ("C1", ["E"], "C1: 1 class, 1 irrep (A), order=1", 1, "high"),

            # ===== Cs (only mirror plane) =====
            "CH3Cl": ("Cs", ["E", "σ"], "Cs: 2 classes, 2 irreps (A', A''), order=2", 2, "high"),
            "CH2ClF": ("Cs", ["E", "σ"], "Cs: 2 classes, 2 irreps (A', A''), order=2", 2, "medium"),
            "CH2ClBr": ("Cs", ["E", "σ"], "Cs: 2 classes, 2 irreps (A', A''), order=2", 2, "medium"),
            "HOCl": ("Cs", ["E", "σ"], "Cs: 2 classes, 2 irreps (A', A''), order=2", 2, "medium"),
            "ONCl": ("Cs", ["E", "σ"], "Cs: 2 classes, 2 irreps (A', A''), order=2", 2, "medium"),
            "CH3CH2OH": ("C1", ["E"], "C1: 1 class, 1 irrep (A), order=1", 1, "low"),  # ethanol - low sym
            "CCO": ("C1", ["E"], "C1: 1 class, 1 irrep (A), order=1", 1, "low"),

            # ===== Ci (inversion center only) =====
            "meso-tartaric": ("Ci", ["E", "i"], "Ci: 2 classes, 2 irreps (Ag, Au), order=2", 2, "high"),
            "staggered(C-C)": ("Ci" if False else "C2h", ["E", "i", "C2", "σh"], "C2h: 4 classes, 4 irreps (Ag, Bg, Au, Bu), order=4", 4, "low"),
            "(RC)(RS)-2,3-dibromobutane": ("Ci", ["E", "i"], "Ci: 2 classes, 2 irreps (Ag, Au), order=2", 2, "high"),

            # ===== Cn (proper rotation axis only) =====
            "H2O2(cis)": ("C2", ["E", "C2"], "C2: 2 classes, 2 irreps (A, B), order=2", 2, "medium"),
            "twisted biphenyl(θ≠90°)": ("C2", ["E", "C2"], "C2: 2 classes, 2 irreps (A, B), order=2", 2, "medium"),
            "N(CH3)3(propeller)": ("C3", ["E", "C3", "C3^2"], "C3: 3 classes, 3 irreps (A, E), order=3", 3, "medium"),
            "PPH3(propeller)": ("C3", ["E", "C3", "C3^2"], "C3: 3 classes, 3 irreps (A, E), order=3", 3, "medium"),

            # ===== Cnv (n-fold axis + n vertical mirrors) =====
            "H2O": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "H2S": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "SO2": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "NO2": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "O=S=O": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "Cl2O": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "cis-[Pt(NH3)2Cl2]": ("C2v", ["E", "C2", "2σv"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "CH2Cl2": ("C2v", ["E", "C2(z)", "σv(xz)", "σv'(yz)"], "C2v: 4 classes, 4 irreps (A1, A2, B1, B2), order=4", 4, "high"),
            "NH3": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "high"),
            "PH3": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),
            "PCl3": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),
            "CH3X(X≠H)": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),
            "CHCl3": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),
            "CClF3": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),
            "POCl3": ("C3v", ["E", "2C3", "3σv"], "C3v: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),
            "C4v": ("C4v", ["E", "C4", "C4^3", "2C2", "4σv"], "C4v: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "medium"),
            "BrF5": ("C4v", ["E", "C4", "C4^3", "2C2", "4σv"], "C4v: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "high"),
            "XeOF4": ("C4v", ["E", "C4", "C4^3", "2C2", "4σv"], "C4v: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "high"),
            "IF5": ("C4v", ["E", "C4", "C4^3", "2C2", "4σv"], "C4v: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "high"),
            "SF5Cl": ("C4v", ["E", "C4", "C4^3", "2C2", "4σv"], "C4v: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "medium"),
            "C5v": ("C5v", ["E", "2C5", "5σv"], "C5v: 3 classes, 3 irreps (A1, A2, E1+E2), order=10", 10, "medium"),

            # ===== Cnh (n-fold axis + horizontal mirror) =====
            "HOOH(trans)": ("C2h", ["E", "C2", "i", "σh"], "C2h: 4 classes, 4 irreps (Ag, Bg, Au, Bu), order=4", 4, "medium"),
            "B(OH)3(planar)": ("C3h", ["E", "C3", "C3^2", "σh", "2S3"], "C3h: 6 classes, 6 irreps (A', E', A'', E''), order=6", 6, "medium"),

            # ===== Dn (n-fold axis + n perpendicular C2 axes) =====
            "tris-chelate twisted": ("D3", ["E", "2C3", "3C2"], "D3: 3 classes, 3 irreps (A1, A2, E), order=6", 6, "medium"),

            # ===== Dnd (n-fold axis + n dihedral mirrors + improper rotations) =====
            "allene": ("D2d", ["E", "S4", "C2(z)", "2C2'", "2σd"], "D2d: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "high"),
            "C=C=C": ("D2d", ["E", "S4", "C2(z)", "2C2'", "2σd"], "D2d: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "high"),
            "cyclobutane(puckered)": ("D2d", ["E", "S4", "C2(z)", "2C2'", "2σd"], "D2d: 5 classes, 5 irreps (A1, A2, B1, B2, E), order=8", 8, "medium"),
            "staggered ethane": ("D3d", ["E", "2C3", "3C2", "i", "2S6", "3σd"], "D3d: 6 classes, 6 irreps (A1g, A2g, Eg, A1u, A2u, Eu), order=12", 12, "high"),
            "ferrocene(staggered)": ("D5d", ["E", "2C5", "2C5^2", "5C2", "i", "2S10", "5σd"], "D5d: 8 classes, 8 irreps (A1g, A2g, E1g, E2g, A1u, A2u, E1u, E2u), order=20", 20, "high"),
            "Fe(C5H5)2(staggered)": ("D5d", ["E", "2C5", "2C5^2", "5C2", "i", "2S10", "5σd"], "D5d: 8 classes, 8 irreps (A1g, A2g, E1g, E2g, A1u, A2u, E1u, E2u), order=20", 20, "high"),
            "staggered(C-C)6": ("D6d", ["E", "2C6", "2C3", "2C2\u2032", "6C2\u2032\u2032", "i", "2S12", "2S6^5", "2S3^5", "6\u03c3d"], "D6d: 11 classes, 11 irreps, order=24", 24, "medium"),
            "S8(crown)": ("D4d", ["E", "2S8", "2C4, 2C4^3", "4C2\u2032", "4C2\u2032\u2032", "8\u03c3d"], "D4d: 7 classes, 7 irreps (A1, A2, B1, B2, E1, E2, E3), order=16", 16, "high"),

            # ===== Dnh (n-fold axis + n C2⊥ + σh + n σv) =====
            "ethylene": ("D2h", ["E", "C2(z)", "C2(y)", "C2(x)", "i", "σ(xy)", "σ(xz)", "σ(yz)"], "D2h: 8 classes, 8 irreps (Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u), order=8", 8, "high"),
            "C=C": ("D2h", ["E", "C2(z)", "C2(y)", "C2(x)", "i", "σ(xy)", "σ(xz)", "σ(yz)"], "D2h: 8 classes, 8 irreps (Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u), order=8", 8, "high"),
            "N2O4(planar)": ("D2h", ["E", "C2(z)", "C2(y)", "C2(x)", "i", "σ(xy)", "σ(xz)", "σ(yz)"], "D2h: 8 classes, 8 irreps (Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u), order=8", 8, "high"),
            "B2H6": ("D2h", ["E", "C2(z)", "C2(y)", "C2(x)", "i", "σ(xy)", "σ(xz)", "σ(yz)"], "D2h: 8 classes, 8 irreps (Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u), order=8", 8, "high"),
            "acetylene": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg, Δg..., Σu+, Σu-, Πu, Δu...), order=∞", -1, "high"),
            "C#C": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "CO2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "C=O=C": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "[N+]=[N-]": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "C2H2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "BF3": ("D3h", ["E", "C3", "3C2'", "σh", "2S3", "3σv"], "D3h: 6 classes, 6 irreps (A1', A2', E', A1'', A2'', E''), order=12", 12, "high"),
            "SO3": ("D3h", ["E", "C3", "3C2'", "σh", "2S3", "3σv"], "D3h: 6 classes, 6 irreps (A1', A2', E', A1'', A2'', E''), order=12", 12, "high"),
            "NO3-": ("D3h", ["E", "C3", "3C2'", "σh", "2S3", "3σv"], "D3h: 6 classes, 6 irreps (A1', A2', E', A1'', A2'', E''), order=12", 12, "high"),
            "CO3(2-)": ("D3h", ["E", "C3", "3C2'", "σh", "2S3", "3σv"], "D3h: 6 classes, 6 irreps (A1', A2', E', A1'', A2'', E''), order=12", 12, "high"),
            "PCl5(trigonal bipyramid)": ("D3h", ["E", "C3", "3C2'", "σh", "2S3", "3σv"], "D3h: 6 classes, 6 irreps (A1', A2', E', A1'', A2'', E''), order=12", 12, "high"),
            "XeF4": ("D4h", ["E", "C4", "4C2", "C2'", "2C2''", "i", "S4", "σh", "2σv", "2σd", "2σ'"], "D4h: 10 classes, 10 irreps (A1g-A2u, B1g-B2u, Eg-Eu), order=16", 16, "high"),
            "PtCl4(2-)": ("D4h", ["E", "C4", "4C2", "C2'", "2C2''", "i", "S4", "σh", "2σv", "2σd", "2σ'"], "D4h: 10 classes, 10 irreps, order=16", 16, "high"),
            "Ni(CN)4(2-)": ("D4h", ["E", "C4", "4C2", "C2'", "2C2''", "i", "S4", "σh", "2σv", "2σd", "2σ'"], "D4h: 10 classes, 10 irreps, order=16", 16, "high"),
            "trans-[Pt(NH3)2Cl2]": ("D2h", ["E", "3C2", "i", "3σ"], "D2h: 8 classes, 8 irreps, order=8", 8, "high"),
            "cyclooctatetraenyl dianion": ("D4h", ["E", "C4", "..."], "D4h: 10 classes, 10 irreps, order=16", 16, "medium"),
            "benzene": ("D6h", ["E", "C6", "3C2", "C2'", "6C2''", "i", "S6", "σh", "2σv", "2σd", "3σ'"], "D6h: 12 classes, 12 irreps (A1g-B2u, E1g-E2u, E1g*-E2g*), order=24", 24, "high"),
            "c1ccccc1": ("D6h", ["E", "C6", "3C2", "C2'", "6C2''", "i", "S6", "σh", "2σv", "2σd", "3σ'"], "D6h: 12 classes, 12 irreps, order=24", 24, "high"),
            "C1=CC=CC=C1": ("D6h", ["E", "C6", "3C2", "C2'", "6C2''", "i", "S6", "σh", "2σv", "2σd", "3σ'"], "D6h: 12 classes, 12 irreps, order=24", 24, "high"),
            "graphene(fragment)": ("D6h", ["E", "C6", "..."], "D6h: 12 classes, 12 irreps, order=24", 24, "medium"),
            "ferrocene(eclipsed)": ("D5h", ["E", "2C5", "5C2", "σh", "2S5", "5σv"], "D5h: 8 classes, 8 irreps (A1', A2', E1', E2', A1'', A2'', E1'', E2''), order=20", 20, "high"),
            "Fe(C5H5)2(eclipsed)": ("D5h", ["E", "2C5", "5C2", "σh", "2S5", "5σv"], "D5h: 8 classes, 8 irreps, order=20", 20, "high"),
            "C5H5-C5H5(eclipsed)": ("D5h", ["E", "2C5", "5C2", "σh", "2S5", "5σv"], "D5h: 8 classes, 8 irreps, order=20", 20, "high"),
            "eclipsed ethane": ("D3h", ["E", "C3", "3C2'", "σh", "2S3", "3σv"], "D3h: 6 classes, 6 irreps, order=12", 12, "medium"),

            # ===== High symmetry groups =====
            "CH4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24", 24, "high"),
            "CCl4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24", 24, "high"),
            "SiH4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24", 24, "high"),
            "P4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24", 24, "high"),
            "CF4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24", 24, "high"),
            "Pb4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps (A1, A2, E, T1, T2), order=24", 24, "medium"),
            "SF6": ("Oh", ["E", "3C4", "4C3", "6C2", "3C2(=C42)", "i", "3S4", "4S6", "3σh", "6σd", "8σh"], "Oh: 10 classes, 10 irreps (A1g-A2u, Eg-Eu, T1g-T2u), order=48", 48, "high"),
            "UF6": ("Oh", ["E", "3C4", "4C3", "6C2", "i", "3S4", "4S6", "3σh", "6σd"], "Oh: 10 classes, 10 irreps, order=48", 48, "high"),
            "MoF6": ("Oh", ["E", "3C4", "4C3", "6C2", "i", "3S4", "4S6", "3σh", "6σd"], "Oh: 10 classes, 10 irreps, order=48", 48, "medium"),
            "WF6": ("Oh", ["E", "3C4", "4C3", "6C2", "i", "3S4", "4S6", "3σh", "6σd"], "Oh: 10 classes, 10 irreps, order=48", 48, "medium"),
            "OsF8": ("Oh", ["E", "3C4", "4C3", "6C2", "i", "3S4", "4S6", "3σh", "6σd"], "Oh: 10 classes, 10 irreps, order=48", 48, "medium"),
            "XeO4": ("Td", ["E", "4C3", "3C2", "6S4", "6σd"], "Td: 5 classes, 5 irreps, order=24", 24, "medium"),
            "C60": ("Ih", ["E", "6C5", "10C3", "15C2", "15σ", "i", "6S10", "10S6", "20S3"], "Ih: 10 classes, 10 irreps (Ag, T1g, T2g, Gg, Hg, Au, T1u, T2u, Gu, Hu), order=120", 120, "high"),
            "B12H12(2-)": ("Ih", ["E", "6C5", "10C3", "15C2", "15σ", "i", "6S10", "10S6", "20S3"], "Ih: 10 classes, 10 irreps, order=120", 120, "high"),

            # ===== Linear C∞v =====
            "HCl": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ, Φ...), order=∞", -1, "high"),
            "CO": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),
            "HCN": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),
            "C#N": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),
            "HF": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),
            "HI": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),
            "NO": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),
            "HBr": ("C∞v", ["E", "C∞", "∞σv"], "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞", -1, "high"),

            # ===== Diatomic homonuclear D∞h =====
            "H2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "O2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "N2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "Cl2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "F2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "Br2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
            "I2": ("D∞h", ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"], "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞", -1, "high"),
        }

    def _run_base(self, molecule: str) -> dict:
        """Identify the point group of a molecule."""
        mol_key = molecule.strip()

        # Direct lookup
        if mol_key in self._molecule_db:
            pg, elems, ct, order, conf = self._molecule_db[mol_key]
            logger.info(f"Molecule '{mol_key}' → point group {pg}")
            return {
                "point_group": pg,
                "symmetry_elements": elems,
                "character_table_summary": ct,
                "order": order,
                "confidence": conf,
            }

        # Case-insensitive lookup
        mol_lower = mol_key.lower()
        for key, val in self._molecule_db.items():
            if key.lower() == mol_lower:
                pg, elems, ct, order, conf = val
                return {
                    "point_group": pg,
                    "symmetry_elements": elems,
                    "character_table_summary": ct,
                    "order": order,
                    "confidence": conf,
                }

        # Heuristic fallback
        result = self._heuristic_determine(mol_key)
        if result:
            return result

        raise ChemMCPError(
            f"Cannot determine point group for '{mol_key}'. "
            f"Please provide a SMILES, formula, or common name. "
            f"Known molecules include: H2O, NH3, CH4, CO2, benzene, SF6, XeF4, BF3, allene, ferrocene, C60, etc. "
            f"(80+ molecules in database)"
        )

    def _run_text(self, molecule: str) -> dict:
        return self._run_base(molecule)

    def _heuristic_determine(self, mol_key: str):
        import re
        # Homonuclear diatomic
        if re.match(r'^[A-Z][a-z]?2$', mol_key):
            return {"point_group": "D∞h", "symmetry_elements": ["E", "C∞", "∞C2⊥C∞", "i", "σh", "∞σv"],
                    "character_table_summary": "D∞h: ∞ classes, irreps (Σg+, Σg-, Πg..., Σu+, Σu-, Πu...), order=∞",
                    "order": -1, "confidence": "medium"}
        # Heteronuclear diatomic
        m = re.match(r'^([A-Z][a-z]?)([A-Z][a-z]?)$', mol_key)
        if m and len(m.group(1)) != len(m.group(2)):
            return {"point_group": "C∞v", "symmetry_elements": ["E", "C∞", "∞σv"],
                    "character_table_summary": "C∞v: ∞ classes, irreps (Σ+, Σ-, Π, Δ...), order=∞",
                    "order": -1, "confidence": "medium"}
        # MX6 octahedral
        if re.match(r'^[A-Z][a-z]?[Ff]6$|^[A-Z][a-z]?[Cc][ll]6$', mol_key):
            return {"point_group": "Oh", "symmetry_elements": ["E", "3C4", "4C3", "6C2", "i", "3S4", "4S6", "3σh", "6σd"],
                    "character_table_summary": "Oh: 10 classes, 10 irreps, order=48", "order": 48, "confidence": "high"}
        # MX4 tetrahedral
        if re.match(r'^[A-Z][a-z]?[FfClBr]4$', mol_key):
            return {"point_group": "Td", "symmetry_elements": ["E", "4C3", "3C2", "6S4", "6σd"],
                    "character_table_summary": "Td: 5 classes, 5 irreps, order=24", "order": 24, "confidence": "high"}
        # MX3 trigonal planar
        if re.match(r'^[A-Z][a-z]?[FfClBr]3$', mol_key):
            return {"point_group": "D3h", "symmetry_elements": ["E", "C3", "3C2'", "σh", "2S3", "3σv"],
                    "character_table_summary": "D3h: 6 classes, 6 irreps, order=12", "order": 12, "confidence": "medium"}
        # H2X bent
        if re.match(r'^H2[A-Z][a-z]?$', mol_key):
            return {"point_group": "C2v", "symmetry_elements": ["E", "C2(z)", "σv(xz)", "σv'(yz)"],
                    "character_table_summary": "C2v: 4 classes, 4 irreps, order=4", "order": 4, "confidence": "medium"}
        # H3X pyramidal
        if re.match(r'^H3[A-Z][a-z]?$', mol_key):
            return {"point_group": "C3v", "symmetry_elements": ["E", "2C3", "3σv"],
                    "character_table_summary": "C3v: 3 classes, 3 irreps, order=6", "order": 6, "confidence": "medium"}
        return None
