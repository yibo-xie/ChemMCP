import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive VSEPR geometry data: (bp, lp) → full analysis
VSEPR_DATA = {
    (2, 0): {
        "steric_number": 2, "hybridization": "sp", "electron_geometry": "linear", "molecular_geometry": "linear",
        "ideal_angle": "180°",
        "3d_description": "Two bonding pairs repel to maximum separation (180°). Linear arrangement minimizes electron pair repulsion.",
        "examples": "CO2, BeCl2, HCN, C2H2 (acetylene), NO2⁻, XeF2 (with 3 LP)",
        "lone_pair_arrangement": "N/A",
        "deviation_note": "No deviation from ideal angle; linear molecules have exactly 180° bond angles.",
        "point_group": "D∞h (homonuclear) or C∞v (heteronuclear)",
    },
    (3, 0): {
        "steric_number": 3, "hybridization": "sp²", "electron_geometry": "trigonal planar", "molecular_geometry": "trigonal planar",
        "ideal_angle": "120°",
        "3d_description": "Three bonding pairs in a plane, 120° apart. All atoms lie in the same plane with the central atom at the center of an equilateral triangle.",
        "examples": "BF3, SO3, NO3⁻, CO3²⁻, AlCl3, GaCl3, HCHO (carbonyl C), C2H4 (each C)",
        "lone_pair_arrangement": "N/A",
        "deviation_note": "Planar trigonal geometry is rigid; bond angles are exactly 120° for symmetric substituents.",
        "point_group": "D3h",
    },
    (2, 1): {
        "steric_number": 3, "electron_geometry": "trigonal planar", "molecular_geometry": "bent/angular",
        "ideal_angle": "<120° (~119°)",
        "3d_description": "Three domains total: 2 BP + 1 LP. Lone pair occupies one vertex of trigonal plane, pushing two bonding pairs together to an angle <120°.",
        "examples": "SO2 (119°), NO2 (134°), SnCl2 (~95°), O3 (117°)",
        "lone_pair_arrangement": "Lone pair in sp² orbital, occupying one trigonal vertex",
        "deviation_note": "LP-BP repulsion > BP-BP repulsion, compressing bond angle below 120°. The exact angle depends on electronegativity differences.",
        "point_group": "C2v",
    },
    (4, 0): {
        "steric_number": 4, "hybridization": "sp³", "electron_geometry": "tetrahedral", "molecular_geometry": "tetrahedral",
        "ideal_angle": "109.5°",
        "3d_description": "Four bonding pairs directed toward vertices of a regular tetrahedron. Central atom at center, four substituents at corners. All bond angles equal at 109.5°.",
        "examples": "CH4, CCl4, SiH4, CF4, NH4⁺, PO₄³⁻, SO₄²⁻, ClO₄⁻, MnO₄⁻",
        "lone_pair_arrangement": "N/A",
        "deviation_note": "Perfect tetrahedral symmetry when all substituents identical. Slight deviations occur with different substituents (e.g., CH3Cl: angles ~109° but not all equal).",
        "point_group": "Td",
    },
    (3, 1): {
        "steric_number": 4, "hybridization": "sp³", "electron_geometry": "tetrahedral", "molecular_geometry": "trigonal pyramidal",
        "ideal_angle": "<109.5° (~107°)",
        "3d_description": "Four domains: 3 BP + 1 LP arranged tetrahedrally. Three atoms form a pyramid base, central atom slightly above base plane. Lone pair occupies the fourth tetrahedral position.",
        "examples": "NH3 (107°), PH3 (93.5°), PCl3 (~100°), AsH3 (91.8°), NF3 (102.5°), SO₃²⁻, ClO₃⁻, H₃O⁺",
        "lone_pair_arrangement": "Lone pair in one sp³ orbital, pointing away from the three bonded atoms (the 'apex' of the inverted tetrahedron)",
        "deviation_note": "LP-BP repulsion > BP-BP, compressing angles from 109.5° toward ~107° (NH3). For heavier group 15 elements (PH3, AsH3), angles approach 90° indicating nearly pure p-orbital bonding.",
        "point_group": "C3v",
    },
    (2, 2): {
        "steric_number": 4, "hybridization": "sp³", "electron_geometry": "tetrahedral", "molecular_geometry": "bent/angular",
        "ideal_angle": "<109.5° (~104.5°)",
        "3d_description": "Four domains: 2 BP + 2 LP in tetrahedral arrangement. Two lone pairs and two bonding pairs. The molecule is bent/V-shaped as both LPs occupy tetrahedral positions.",
        "examples": "H2O (104.5°), H2S (92.3°), H2Se (91°), H2Te (~90°), OF2 (103.8°), SCl2 (~103°), Cl2O (110.9°)",
        "lone_pair_arrangement": "Two lone pairs in two sp³ orbitals, occupying ~2/3 of tetrahedral volume. Strong LP-LP repulsion further compresses bond angle beyond pyramidal case.",
        "deviation_note": "LP-LP repulsion >> LP-BP >> BP-BP. In H2O, this gives 104.5° (vs 107° for NH3). For H2S/H2Se/H2Te, angles approach 90° as s-character in bonding orbitals decreases.",
        "point_group": "C2v",
    },
    (5, 0): {
        "steric_number": 5, "hybridization": "sp³d", "electron_geometry": "trigonal bipyramidal", "molecular_geometry": "trigonal bipyramidal",
        "ideal_angle": "90°, 120°, 180°",
        "3d_description": "Five bonding pairs: 3 equatorial (120° apart in a plane) + 2 axial (180°, perpendicular to equatorial plane). Axial bonds are typically longer than equatorial bonds due to more 90° repulsions (3 vs 2).",
        "examples": "PCl5, PF5, AsF5, SbF5, PCl5(gas phase)",
        "lone_pair_arrangement": "N/A",
        "deviation_note": "Equatorial-equatorial angles are 120°, axial-equatorial are 90°, axial-axial is 180°. Axial bonds longer than equatorial (e.g., PCl5: axial 240pm, equatorial 202pm).",
        "point_group": "D3h",
    },
    (4, 1): {
        "steric_number": 5, "electron_geometry": "trigonal bipyramidal", "molecular_geometry": "see-saw / distorted tetrahedron",
        "ideal_angle": "~90°, ~120°, <180°",
        "3d_description": "Five domains: 4 BP + 1 LP. Lone pair occupies an EQUATORIAL position (minimizes 90° repulsions: only 2 axial BPs at 90° vs 3 if axial). Result is a see-saw shape.",
        "examples": "SF4 (axial F-S-F = 173°, equatorial F-S-F = 101.6°), TeCl4, IF4⁺",
        "lone_pair_arrangement": "Lone pair in equatorial position (preferred over axial because it has only 2 neighbors at 90° instead of 3)",
        "deviation_note": "Strong LP-BP repulsion bends axial bonds away from 180° (SF4: 173°). Equatorial bonds also bent from ideal 120° by LP presence.",
        "point_group": "C2v",
    },
    (3, 2): {
        "steric_number": 5, "electron_geometry": "trigonal bipyramidal", "molecular_geometry": "T-shaped",
        "ideal_angle": "~90°, 180°",
        "3d_description": "Five domains: 3 BP + 2 LP. Both LPs occupy equatorial positions (trans to each other, minimizing LP-LP repulsion at 120°). Three atoms form a T-shape: 2 axial + 1 equatorial.",
        "examples": "ClF3 (F-Cl-F(axial)=175°, equatorial=87.5°), BrF3 (~86°, 172°), IF3",
        "lone_pair_arrangement": "Two lone pairs in equatorial positions, 120° apart (maximizing LP-LP separation)",
        "deviation_note": "Two strong LP-BP repulsions bend axial bonds from 180°. ClF3 axial angle compressed to ~175°.",
        "point_group": "C2v",
    },
    (2, 3): {
        "steric_number": 5, "electron_geometry": "trigonal bipyramidal", "molecular_geometry": "linear",
        "ideal_angle": "180°",
        "3d_description": "Five domains: 2 BP + 3 LP. All 3 LPs occupy equatorial positions (120° apart, maximizing LP-LP separation). Two axial bonding pairs give linear molecular geometry.",
        "examples": "XeF2, I3⁻, ICl₂⁻, KrF₂",
        "lone_pair_arrangement": "Three lone pairs in all three equatorial positions (120° apart); bonding pairs occupy axial positions",
        "deviation_note": "Despite having only 2 bonded atoms, steric number is 5! Linear shape arises from TBP e-domain geometry with 3 equatorial LPs.",
        "point_group": "D∞h",
    },
    (6, 0): {
        "steric_number": 6, "hybridization": "sp³d²", "electron_geometry": "octahedral", "molecular_geometry": "octahedral",
        "ideal_angle": "90°, 180°",
        "3d_description": "Six bonding pairs directed toward vertices of a regular octahedron. Central atom at body center, 6 substituents at face centers of a cube. All adjacent angles 90°, opposite atoms 180°. All positions equivalent.",
        "examples": "SF6, UF6, MoF6, WF6, [Fe(CN)6]³⁻, [Co(NH3)6]³⁺, [Fe(H2O)6]³⁺",
        "lone_pair_arrangement": "N/A",
        "deviation note": "Perfect octahedral symmetry. All bond lengths equal, all 90° angles identical. Highest coordination number without lone pairs among common geometries.",
        "point_group": "Oh",
    },
    (5, 1): {
        "steric_number": 6, "electron_geometry": "octahedral", "molecular_geometry": "square pyramidal",
        "ideal_angle": "~90° (basal), <90° (apex-basal)",
        "3d_description": "Six domains: 5 BP + 1 LP. Lone pair occupies one apex position. Five atoms form a square pyramid: 4 basal (square) + 1 apical. LP pushes basal atoms slightly downward.",
        "examples": "BrF5 (basal F-Br-F = 90°, apical F-Br-F(basal) = 84.8°), IF5 (82°, 84°), XeOF4",
        "lone_pair_arrangement": "Lone pair at one octahedral apex; trans basal atoms bent slightly toward each other",
        "deviation_note": "LP-BP repulsion compresses basal-apex angles below 90° (BrF5: 84.8°). Basal square remains near 90°.",
        "point_group": "C4v",
    },
    (4, 2): {
        "steric_number": 6, "electron_geometry": "octahedral", "molecular_geometry": "square planar",
        "ideal_angle": "90°, 180°",
        "3d_description": "Six domains: 4 BP + 2 LP. Two LPs occupy TRANS (opposite) AXIAL positions — this minimizes LP-LP repulsion (180° apart, no 90° LP-LP interaction). Four bonding pairs form a perfect square in the equatorial plane.",
        "examples": "XeF4, PtCl4²⁻, Ni(CN)4²⁻, PdCl4²⁻, AuCl4⁻, ICl4⁻",
        "lone_pair_arrangement": "Two lone pairs trans to each other at axial positions (180° apart, maximally separated)",
        "deviation_note": "Trans arrangement of LPs is crucial — cis would place them at 90° (much higher repulsion). Square planar complexes often involve d⁸ metal ions (Ni²⁺, Pd²⁺, Pt²⁺, Au³⁺).",
        "point_group": "D4h",
    },
}

# Known molecule → (bonding_pairs, lone_pairs) lookup
MOLECULE_VSEPR_DB = {
    # Steric number 2
    "CO2": (2, 0), "BeCl2": (2, 0), "HCN": (2, 0), "C2H2": (2, 0), "acetylene": (2, 0),
    "NO2+": (2, 0),
    # Steric number 3 - trigonal planar
    "BF3": (3, 0), "SO3": (3, 0), "NO3-": (3, 0), "CO3(2-)": (3, 0), "AlCl3": (3, 0),
    "GaCl3": (3, 0), "HCHO": (3, 0), "formaldehyde": (3, 0), "C2H4": (3, 0),
    "ethylene": (3, 0), "benzene": (3, 0), "c1ccccc1": (3, 0), "graphite": (3, 0),
    # Steric number 3 - bent
    "SO2": (2, 1), "NO2": (2, 1), "SnCl2": (2, 1), "O3": (2, 1), "PbCl2": (2, 1),
    # Steric number 4 - tetrahedral
    "CH4": (4, 0), "CCl4": (4, 0), "CF4": (4, 0), "SiH4": (4, 0), "GeH4": (4, 0),
    "SnH4": (4, 0), "PbH4": (4, 0), "NH4+": (4, 0), "PO4(3-)": (4, 0),
    "SO4(2-)": (4, 0), "ClO4-": (4, 0), "MnO4-": (4, 0), "CrO4(2-)": (4, 0),
    "SiO4(4-)": (4, 0), "diamond": (4, 0), "CH3Cl": (4, 0), "CH2Cl2": (4, 0),
    "CHCl3": (4, 0), "CCl3F": (4, 0), "POCl3": (4, 0), "OsO4": (4, 0),
    # Steric number 4 - pyramidal
    "NH3": (3, 1), "PH3": (3, 1), "AsH3": (3, 1), "SbH3": (3, 1), "NF3": (3, 1),
    "PCl3": (3, 1), "PBr3": (3, 1), "PI3": (3, 1), "PF3": (3, 1), "AsCl3": (3, 1),
    "SO3(2-)": (3, 1), "ClO3-": (3, 1), "H3O+": (3, 1), "P(CH3)3": (3, 1),
    # Steric number 4 - bent
    "H2O": (2, 2), "H2S": (2, 2), "H2Se": (2, 2), "H2Te": (2, 2), "OF2": (2, 2),
    "SCl2": (2, 2), "SBr2": (2, 2), "Cl2O": (2, 2), "H2O2": (2, 2),
    # Steric number 5 - trigonal bipyramidal
    "PCl5": (5, 0), "PF5": (5, 0), "AsF5": (5, 0), "SbF5": (5, 0),
    # Steric number 5 - see-saw
    "SF4": (4, 1), "TeCl4": (4, 1), "IF4+": (4, 1),
    # Steric number 5 - T-shaped
    "ClF3": (3, 2), "BrF3": (3, 2), "IF3": (3, 2), "ICl3": (3, 2),
    # Steric number 5 - linear (TBP)
    "XeF2": (2, 3), "I3-": (2, 3), "ICl2-": (2, 3), "KrF2": (2, 3),
    # Steric number 6 - octahedral
    "SF6": (6, 0), "UF6": (6, 0), "MoF6": (6, 0), "WF6": (6, 0), "ReF6": (6, 0),
    "[Fe(CN)6](3-)": (6, 0), "[Fe(H2O)6](3+)": (6, 0), "[Co(NH3)6](3+)": (6, 0),
    # Steric number 6 - square pyramidal
    "BrF5": (5, 1), "IF5": (5, 1), "XeOF4": (5, 1),
    # Steric number 6 - square planar
    "XeF4": (4, 2), "PtCl4(2-)": (4, 2), "Ni(CN)4(2-)": (4, 2), "PdCl4(2-)": (4, 2),
    "AuCl4(-)": (4, 2), "ICl4(-)": (4, 2),
}


@ChemMCPManager.register_tool
class VseprGeometry(BaseTool):
    """
    VSEPR理论预测分子几何构型工具 (MCP #296)。
    基于价层电子对互斥理论(VSEPR)，全面预测分子的：
    - 电子域几何构型和分子几何构型
    - 理想键角与实际键角（含偏差说明）
    - 杂化轨道类型
    - 孤对电子排布及其影响
    - 3D结构描述与点群
    比现有 PredictVseprGeometry 更全面详细。
    """
    __version__ = "0.1.0"
    name = "VseprGeometry"
    func_name = 'predict_vsepr_geometry'
    description = "Comprehensive VSEPR-based molecular geometry prediction with detailed analysis of electron domain geometry, bond angles, hybridization, lone pair effects, and 3D structure."
    implementation_description = (
        "Uses VSEPR theory with full (bp, lp) combinations for steric numbers 2-6. "
        "Provides idealized and actual bond angles, lone pair arrangement details, "
        "deviation explanations, point group assignment, and 3D structural descriptions."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["VSEPR", "Molecular Geometry", "Hybridization", "Chemical Bonding", "Steric Number"]
    required_envs = []

    code_input_sig = [
        ('molecule', 'str', 'N/A', 'Molecule identifier or use bonding_pairs+lone_pairs mode'),
        ('bonding_pairs', 'int', 'None', 'Number of bonding pairs around central atom (overrides molecule lookup)'),
        ('lone_pairs', 'int', 'None', 'Number of lone pairs on central atom (overrides molecule lookup)'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Query: molecule name/formula or "bp lp" e.g., "H2O", "SF6", "4 2".'),
    ]
    output_sig = [
        ('molecule', 'str', 'Molecule identifier'),
        ('steric_number', 'int', 'Total electron domains (BP + LP)'),
        ('bonding_pairs', 'int', 'Number of bonding pairs'),
        ('lone_pairs', 'int', 'Number of lone pairs'),
        ('electron_domain_geometry', 'str', 'Arrangement of ALL electron domains'),
        ('molecular_geometry', 'str', 'Arrangement of ATOMS only'),
        ('ideal_bond_angle', 'str', 'Ideal bond angle(s)'),
        ('hybridization', 'str', 'Central atom hybridization'),
        ('coordination_number', 'int', 'Number of bonded atoms'),
        ('lone_pair_arrangement', 'str', 'Description of lone pair placement'),
        ('deviation_explanation', 'str', 'Why actual angles differ from ideal'),
        ('3d_structure_description', 'str', 'Detailed 3D spatial description'),
        ('point_group', 'str', 'Molecular point group (symmetry)'),
        ('examples', 'str', 'Example molecules with same geometry'),
    ]

    examples = [{'code_input': {'molecule': 'H2O', 'bonding_pairs': 'N/A', 'lone_pairs': 'N/A'}, 'text_input': {'query': 'H2O'}, 'output': {'molecule': 'H2O', 'steric_number': 4, 'bonding_pairs': 2, 'lone_pairs': 2, 'molecular_geometry': 'bent/angular', 'ideal_bond_angle': '104.5°', 'hybridization': 'sp³', 'point_group': 'C2v', '3d_structure_description': 'N/A', 'coordination_number': 'N/A', 'deviation_explanation': 'N/A', 'electron_domain_geometry': 'N/A', 'examples': 'N/A', 'lone_pair_arrangement': 'N/A'}}, {'code_input': {'molecule': 'XeF4', 'bonding_pairs': 'N/A', 'lone_pairs': 'N/A'}, 'text_input': {'query': 'XeF4'}, 'output': {'molecule': 'XeF4', 'steric_number': 6, 'bonding_pairs': 4, 'lone_pairs': 2, 'molecular_geometry': 'square planar', 'ideal_bond_angle': '90°, 180°', 'hybridization': 'sp³d²', 'point_group': 'D4h', '3d_structure_description': 'N/A', 'coordination_number': 'N/A', 'deviation_explanation': 'N/A', 'electron_domain_geometry': 'N/A', 'examples': 'N/A', 'lone_pair_arrangement': 'N/A'}}, {'code_input': {'molecule': 'custom', 'bonding_pairs': 5, 'lone_pairs': 1}, 'text_input': {'query': '5 1'}, 'output': {'molecular_geometry': 'square pyramidal', 'hybridization': 'sp³d', '3d_structure_description': 'N/A', 'bonding_pairs': 'N/A', 'coordination_number': 'N/A', 'deviation_explanation': 'N/A', 'electron_domain_geometry': 'N/A', 'examples': 'N/A', 'ideal_bond_angle': 'N/A', 'lone_pair_arrangement': 'N/A', 'lone_pairs': 'N/A', 'molecule': 'N/A', 'point_group': 'N/A', 'steric_number': 'N/A'}}]

    def __init__(self, init: bool = True, interface: str = 'code'):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, molecule: str, bonding_pairs: int = None, lone_pairs: int = None) -> dict:
        mol = molecule.strip()

        if bonding_pairs is not None and lone_pairs is not None:
            bp, lp = int(bonding_pairs), int(lone_pairs)
            sn = bp + lp
            mol_label = f"custom ({bp} BP, {lp} LP)"
        else:
            # Look up from database
            mol_upper = mol.upper().replace(" ", "")
            result = MOLECULE_VSEPR_DB.get(mol_upper)
            if result is None:
                # Try case-insensitive
                for k, v in MOLECULE_VSEPR_DB.items():
                    if k.upper().replace(" ", "") == mol_upper:
                        result = v
                        break
            if result is None:
                avail = sorted(set(MOLECULE_VSEPR_DB.keys()))
                raise ChemMCPInputError(
                    f"Cannot determine VSEPR for '{mol}'. Provide bonding_pairs & lone_pairs explicitly, "
                    f"or use known molecule. Available: {avail[:40]}..."
                )
            bp, lp = result
            sn = bp + lp
            mol_label = mol

        key = (bp, lp)
        if key not in VSEPR_DATA:
            raise ChemMCPInputError(
                f"No VSEPR data for bonding_pairs={bp}, lone_pairs={lp}. "
                f"Supported: steric numbers 2-6 with all valid (bp,lp) combinations."
            )

        g = VSEPR_DATA[key]
        return {
            "molecule": mol_label,
            "steric_number": bp + lp,
            "hybridization": g.get("hybridization", g.get("hyb", "sp³")),
            "bonding_pairs": bp,
            "lone_pairs": lp,
            "electron_domain_geometry": g["electron_geometry"],
            "molecular_geometry": g["molecular_geometry"],
            "lone_pair_arrangement": g["lone_pair_arrangement"],
            "deviation_explanation": g.get("deviation_note", g.get("deviation note", "")),
            "3d_structure_description": g["3d_description"],
            "point_group": g.get("point_group", "varies"),
            "example_molecules": g.get("examples", ""),
            "ideal_bond_angle": g.get("angles", "?"),
            "coordination_number": bp + lp,
        }

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            return self._run_base(molecule="custom", bonding_pairs=int(parts[0]), lone_pairs=int(parts[1]))
        return self._run_base(molecule=query)
