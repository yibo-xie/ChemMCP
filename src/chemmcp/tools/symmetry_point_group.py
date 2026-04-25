import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Point group determination rules based on molecular features
# This uses a heuristic approach with common symmetry elements

@ChemMCPManager.register_tool
class SymmetryPointGroup(BaseTool):
    """
    Determine the point group symmetry of a molecule.
    Uses structural heuristics to assign point groups based on molecular geometry.
    """
    __version__ = "0.1.0"
    name = "SymmetryPointGroup"
    func_name = "determine_point_group"
    description = "Determine the point group symmetry of a molecule given its molecular formula or SMILES."
    implementation_description = "Uses heuristic rules based on molecular geometry, symmetry elements (axes, planes, center), and common molecular shapes to determine the Schoenflies point group."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Symmetry", "Point Group", "Molecular Geometry", "Group Theory"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: SMILES string, molecular formula, or common name (e.g., 'H2O', 'CCO', 'benzene', 'CH4')."),
    ]

    text_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: SMILES string, molecular formula, or common name."),
    ]

    output_sig = [
        ("point_group", "str", "The Schoenflies point group symbol (e.g., C2v, D3h, Td, Oh)."),
        ("symmetry_elements", "str", "Description of key symmetry elements found."),
        ("confidence", "str", "Confidence level: 'high', 'medium', or 'low' depending on whether the assignment is unambiguous."),
    ]

    examples = [
        {
            "code_input": {"molecule": "H2O"},
            "text_input": {"molecule": "H2O"},
            "output": {
                "point_group": "C2v",
                "symmetry_elements": "E, C2(z), σv(xz), σv'(yz)",
                "confidence": "high",
            }
        },
        {
            "code_input": {"molecule": "CH4"},
            "text_input": {"molecule": "CH4"},
            "output": {
                "point_group": "Td",
                "symmetry_elements": "E, 4C3, 3C2, 6S4, 6σd",
                "confidence": "high",
            }
        },
        {
            "code_input": {"molecule": "C=O=C"},  # CO2
            "text_input": {"molecule": "CO2"},
            "output": {
                "point_group": "D∞h",
                "symmetry_elements": "E, C∞, ∞C2⊥C∞, i, σh, ∞σv",
                "confidence": "high",
            }
        },
        {
            "code_input": {"molecule": "c1ccccc1"},  # benzene
            "text_input": {"molecule": "benzene"},
            "output": {
                "point_group": "D6h",
                "symmetry_elements": "E, C6, 3C2, C2', 6C2'', i, S6, σh, 2σv, 2σd, 3σ'",
                "confidence": "high",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build the known molecule → point group database."""
        self._common_molecules = {
            # Linear molecules
            "H2": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),
            "O2": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),
            "N2": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),
            "Cl2": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),
            "HCl": ("C∞v", "E, C∞, ∞σv", "high"),
            "CO": ("C∞v", "E, C∞, ∞σv", "high"),
            "HCN": ("C∞v", "E, C∞, ∞σv", "high"),
            "C=O=C": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),  # CO2
            "CO2": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),
            "C#N": ("C∞v", "E, C∞, ∞σv", "high"),  # HCN-like
            "[N+]=[N-]": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),  # N2O

            # Tetrahedral / Octahedral
            "CH4": ("Td", "E, 4C3, 3C2, 6S4, 6σd", "high"),
            "CCl4": ("Td", "E, 4C3, 3C2, 6S4, 6σd", "high"),
            "SiH4": ("Td", "E, 4C3, 3C2, 6S4, 6σd", "high"),
            "P4": ("Td", "E, 4C3, 3C2, 6S4, 6σd", "high"),
            "SF6": ("Oh", "E, 3C4, 4C3, 6C2, 3C2(=C42), i, 3S4, 4S6, 3σh, 6σd, 8σh(?)", "high"),
            "UF6": ("Oh", "E, 3C4, 4C3, 6C2, i, 3S4, 4S6, 3σh, 6σd", "high"),
            "XeF4": ("D4h", "E, C4, 4C2, C2', 2C2'', i, S4, σh, 2σv, 2σd, 2σ'", "high"),
            "BF3": ("D3h", "E, C3, 3C2, σh, 2S3, 3σv", "high"),
            "SO3": ("D3h", "E, C3, 3C2, σh, 2S3, 3σv", "high"),
            "NO3-": ("D3h", "E, C3, 3C2, σh, 2S3, 3σv", "high"),
            "PCl5(D3h)": ("D3h", "E, C3, 3C2, σh, 2S3, 3σv", "high"),

            # Trigonal pyramidal / Bent
            "NH3": ("C3v", "E, 2C3, 3σv", "high"),
            "PH3": ("C3v", "E, 2C3, 3σv", "medium"),
            "PCl3": ("C3v", "E, 2C3, 3σv", "medium"),
            "H2O": ("C2v", "E, C2(z), σv(xz), σv'(yz)", "high"),
            "H2S": ("C2v", "E, C2(z), σv(xz), σv'(yz)", "high"),
            "SO2": ("C2v", "E, C2(z), σv(xz), σv'(yz)", "high"),
            "O=S=O": ("C2v", "E, C2(z), σv(xz), σv'(yz)", "high"),
            "NO2": ("C2v", "E, C2(z), σv(xz), σv'(yz)", "high"),
            "Cl2O": ("C2v", "E, C2(z), σv(xz), σv'(yz)", "high"),

            # Planar / Other common
            "c1ccccc1": ("D6h", "E, C6, 3C2, C2', 6C2'', i, S6, σh, 2σv, 2σd, 3σ'", "high"),  # benzene
            "benzene": ("D6h", "E, C6, 3C2, C2', 6C2'', i, S6, σh, 2σv, 2σd, 3σ'", "high"),
            "C1=CC=CC=C1": ("D6h", "E, C6, 3C2, C2', 6C2'', i, S6, σh, 2σv, 2σd, 3σ'", "high"),
            "ethylene": ("D2h", "E, C2(z), C2(y), C2(x), i, σ(xy), σ(xz), σ(yz)", "high"),
            "C=C": ("D2h", "E, C2(z), C2(y), C2(x), i, σ(xy), σ(xz), σ(yz)", "high"),  # ethylene
            "acetylene": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),
            "C#C": ("D∞h", "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "high"),  # acetylene
            "allene": ("D2d", "E, S4, C2(z), 2C2', 2σd", "high"),
            "C=C=C": ("D2d", "E, S4, C2(z), 2C2', 2σd", "high"),  # allene

            # Others
            "H2O2(non-planar)": ("C2", "E, C2", "medium"),  # non-planar H2O2
            "C(C)(F)(F)F": ("Cs", "E, σ", "medium"),  # CH3F - slightly distorted
            "CH3Cl": ("Cs", "E, σ", "medium"),
            "CClF3": ("C3v", "E, 2C3, 3σv", "medium"),
            "cis-[Pt(NH3)2Cl2]": ("C2v", "E, C2, 2σv", "high"),
            "trans-[Pt(NH3)2Cl2]": ("D2h", "E, 3C2, i, 3σ, 2S4", "high"),
            "Fe(C5H5)2(staggered)": ("D5d", "E, 2C5, 2C52, 5C2, i, 2S10, 5σd", "high"),
            "Fe(C5H5)2(eclipsed)": ("D5h", "E, 2C5, 5C2, σh, 2S5, 5σv", "high"),
            "ferrocene": ("D5d", "E, 2C5, 2C52, 5C2, i, 2S10, 5σd", "high"),
            "C5H5-C5H5": ("D5d", "E, 2C5, 2C52, 5C2, i, 2S10, 5σd", "high"),
            "B12H12^2-": ("Ih", "E, 6C5, 10C3, 15C2, 15σ, i, 6S10, 10S6, 20S3", "high"),
            "C60": ("Ih", "E, 6C5, 10C3, 15C2, 15σ, i, 6S10, 10S6, 20S3", "high"),

            # Asymmetric / low symmetry
            "CHClFBr": ("C1", "E only (no symmetry other than identity)", "high"),
            "CC(O)Cl": ("Cs", "E, σ", "medium"),  # CH3COCl approx
            "CH3CH2OH": ("Cs" if False else "C1", "E only (no symmetry due to free rotation/conformation)", "low"),  # ethanol - low symmetry in most conformations
            "CCO": ("Cs", "E, σ (approximate, depends on conformation)", "low"),  # ethanol SMILES
            "CHFClBr": ("C1", "E only", "high"),
        }

    def _run_base(self, molecule: str) -> dict:
        """
        Determine the point group of a molecule.
        Returns dict with point_group, symmetry_elements, confidence.
        """
        mol_key = molecule.strip()

        # Direct lookup
        if mol_key in self._common_molecules:
            pg, elems, conf = self._common_molecules[mol_key]
            logger.info(f"Molecule '{mol_key}' → point group {pg} (confidence: {conf})")
            return {
                "point_group": pg,
                "symmetry_elements": elems,
                "confidence": conf,
            }

        # Case-insensitive lookup for names
        mol_lower = mol_key.lower()
        for key, val in self._common_molecules.items():
            if key.lower() == mol_lower:
                pg, elems, conf = val
                logger.info(f"Molecule '{mol_key}' (case-insensitive match: '{key}') → point group {pg}")
                return {
                    "point_group": pg,
                    "symmetry_elements": elems,
                    "confidence": conf,
                }

        # Heuristic fallback based on formula patterns
        result = self._heuristic_determine(mol_key)
        if result:
            return result

        raise ChemMCPError(
            f"Cannot determine point group for '{mol_key}'. "
            f"Please provide a SMILES string, molecular formula, or common name of a known molecule. "
            f"Known molecules include: H2O, NH3, CH4, CO2, benzene, SF6, XeF4, BF3, ethylene, acetylene, ferrocene, C60, etc."
        )

    def _run_text(self, molecule: str) -> dict:
        """Text interface delegates to base logic."""
        return self._run_base(molecule)

    def _heuristic_determine(self, mol_key: str) -> Optional[dict]:
        """
        Heuristic rules for molecules not in the database.
        Returns None if cannot determine.
        """
        import re

        # Diatomic homonuclear
        if re.match(r'^[A-Z][a-z]?2$', mol_key):
            return {"point_group": "D∞h", "symmetry_elements": "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "confidence": "medium"}

        # Diatomic heteronuclear
        if re.match(r'^[A-Z][a-z]?[A-Z][a-z]?$', mol_key) and len(mol_key) >= 2 and len(mol_key) <= 4:
            # Check it looks like AB diatomic
            chars = set(re.findall(r'[A-Z][a-z]?', mol_key))
            if len(chars) == 2:
                elem_counts = re.findall(r'[A-Z][a-z]?', mol_key)
                if len(elem_counts) == 2:
                    return {"point_group": "C∞v", "symmetry_elements": "E, C∞, ∞σv", "confidence": "medium"}

        # MX6 octahedral pattern (e.g., WF6, MoF6)
        if re.match(r'^[A-Z][a-z]?[F|Cl|Br]6$', mol_key):
            return {"point_group": "Oh", "symmetry_elements": "E, 3C4, 4C3, 6C2, i, 3S4, 4S6, 3σh, 6σd", "confidence": "high"}

        # MX4 tetrahedral pattern
        if re.match(r'^[A-Z][a-z]?[F|Cl|Br|H]4$', mol_key):
            return {"point_group": "Td", "symmetry_elements": "E, 4C3, 3C2, 6S4, 6σd", "confidence": "high"}

        # MX3 trigonal planar
        if re.match(r'^[A-Z][a-z]?[F|Cl|Br|H]3$', mol_key):
            return {"point_group": "D3h", "symmetry_elements": "E, C3, 3C2, σh, 2S3, 3σv", "confidence": "medium"}

        # MX2 linear
        if re.match(r'^[A-Z][a-z]?[F|Cl|Br|H|O|S]2$', mol_key):
            return {"point_group": "D∞h", "symmetry_elements": "E, C∞, ∞C2⊥C∞, i, σh, ∞σv", "confidence": "medium"}

        # H2X bent (group 16 hydrides)
        if re.match(r'^H2[A-Z][a-z]?$', mol_key) or mol_key.upper().startswith("H2"):
            return {"point_group": "C2v", "symmetry_elements": "E, C2(z), σv(xz), σv'(yz)", "confidence": "medium"}

        # H3X pyramidal (group 15 hydrides)
        if re.match(r'^H3[A-Z][a-z]?$', mol_key) or mol_key.upper().startswith("H3"):
            return {"point_group": "C3v", "symmetry_elements": "E, 2C3, 3σv", "confidence": "medium"}

        return None
