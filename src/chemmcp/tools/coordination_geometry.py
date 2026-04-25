import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CoordinationGeometry(BaseTool):
    """
    Query coordination geometry (octahedral, tetrahedral, square planar, etc.)
    based on coordination number, ligand type, and electronic configuration.
    """
    __version__ = "0.1.0"
    name = "CoordinationGeometry"
    func_name = "query_coordination_geometry"
    description = "Query and predict coordination geometry for complexes based on coordination number, metal center, ligands, and VSEPR/CFSE considerations."
    implementation_description = "Uses a rule-based database of common coordination geometries with their coordination numbers, typical bond angles, examples, and determining factors (hybridization, steric effects, electronic configuration)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Geometry", "Complexes", "VSEPR", "Crystal Field"]
    required_envs = []

    code_input_sig = [
        ("coordination_number", "int", "N/A", "Coordination number (CN) of the complex, e.g., 2, 4, 5, 6."),
        ("metal_ion", "str", "None", "Metal ion symbol (e.g., 'Fe2+', 'Cu2+', 'Co3+'). Optional but improves accuracy."),
        ("ligand_type", "str", "None", "Ligand type or description (e.g., 'strong field', 'weak field', 'bidentate', 'Cl-', 'NH3'). Optional."),
        ("electron_config", "str", "None", "d-electron count or config (e.g., 'd6', 'd8', 'low-spin', 'high-spin'). Optional."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'CN [metal_ion] [ligand_type] [electron_config]', e.g., '6 Fe2+ weak-field d6' or '4 Pt2+ strong-field'."),
    ]

    output_sig = [
        ("geometry", "str", "Predicted coordination geometry name (e.g., 'Octahedral', 'Square Planar', 'Tetrahedral')."),
        ("hybridization", "str", "Hybridization scheme (e.g., 'sp3d2', 'dsp2', 'sp3')."),
        ("bond_angles", "str", "Typical bond angles in degrees."),
        ("examples", "str", "Example compounds with this geometry."),
        ("description", "str", "Brief explanation of why this geometry is favored."),
    ]

    examples = [
        {
            "code_input": {
                "coordination_number": 6,
                "metal_ion": "Fe2+",
                "ligand_type": "weak-field",
                "electron_config": "d6 high-spin",
            },
            "text_input": {
                "query": "6 Fe2+ weak-field d6 high-spin"
            },
            "output": {
                "geometry": "Octahedral",
                "hybridization": "sp3d2",
                "bond_angles": "90°, 180°",
                "examples": "[Fe(H2O)6]2+, [FeF6]3-, [CoF6]3-",
                "description": "CN=6 most commonly adopts octahedral geometry. Weak-field ligands favor high-spin d6 with sp3d2 hybridization.",
            }
        },
        {
            "code_input": {
                "coordination_number": 4,
                "metal_ion": "Pt2+",
                "ligand_type": "strong-field",
                "electron_config": "d8",
            },
            "text_input": {
                "query": "4 Pt2+ strong-field d8"
            },
            "output": {
                "geometry": "Square Planar",
                "hybridization": "dsp2",
                "bond_angles": "90°, 180°",
                "examples": "[PtCl4]2-, [Ni(CN)4]2-, [PdCl4]2-",
                "description": "d8 metals with strong-field ligands often adopt square planar geometry due to large CFSE gain from dsp2 hybridization.",
            }
        },
        {
            "code_input": {
                "coordination_number": 4,
                "metal_ion": "Zn2+",
                "ligand_type": None,
                "electron_config": "d10",
            },
            "text_input": {
                "query": "4 Zn2+ d10"
            },
            "output": {
                "geometry": "Tetrahedral",
                "hybridization": "sp3",
                "bond_angles": "109.5°",
                "examples": "[ZnCl4]2-, [Cd(SCN)4]2-, [MnO4]-",
                "description": "d10 configuration has no CFSE preference; tetrahedral is favored by sterics.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build coordination geometry database."""
        self._geometry_db = {
            2: {
                "default": {
                    "geometry": "Linear",
                    "hybridization": "sp",
                    "bond_angles": "180°",
                    "examples": "[Ag(NH3)2]+, [Au(CN)2]-, [CuCl2]-",
                    "description": "CN=2 almost always gives linear geometry due to sp hybridization minimizing repulsion.",
                },
            },
            3: {
                "default": {
                    "geometry": "Trigonal Planar",
                    "hybridization": "sp2",
                    "bond_angles": "120°",
                    "examples": "[HgI3]-, [Cu(CN)3]2-, [AgCl3]2-",
                    "description": "CN=3 typically trigonal planar; can be T-shaped for some d8/d10 systems.",
                },
                "special": {
                    "geometry": "T-shaped",
                    "hybridization": "sp2d",
                    "bond_angles": "90°, 180°",
                    "examples": "[ClF3], [BrF3]",
                    "description": "T-shaped geometry occurs when there are 3 bonds + 2 lone pairs (AX3E2 in VSEPR).",
                },
            },
            4: {
                "default": {
                    "geometry": "Tetrahedral",
                    "hybridization": "sp3",
                    "bond_angles": "109.5°",
                    "examples": "[ZnCl4]2-, [FeCl4]-, [MnO4]-, [NiCl4]2-(high-spin)",
                    "description": "Default CN=4 geometry. Favored by d0, d5(high-spin), d10 configurations.",
                },
                "d8_strong_field": {
                    "geometry": "Square Planar",
                    "hybridization": "dsp2",
                    "bond_angles": "90°, 180°",
                    "examples": "[PtCl4]2-, [Ni(CN)4]2-, [PdCl4]2-, [AuCl4]-",
                    "description": "d8 metals (Ni2+, Pd2+, Pt2+, Au3+) with strong-field ligands adopt square planar due to large CFSE from dsp2.",
                },
                "see_saw": {
                    "geometry": "See-saw (Distorted Tetrahedral)",
                    "hybridization": "sp3d",
                    "bond_angles": "~90°, ~120°, ~180°",
                    "examples": "[SF4], [TeCl4]",
                    "description": "See-saw geometry: 4 bonding pairs + 1 lone pair (AX4E).",
                },
            },
            5: {
                "default": {
                    "geometry": "Trigonal Bipyramidal",
                    "hybridization": "sp3d",
                    "bond_angles": "90°, 120°, 180°",
                    "examples": "[Fe(CO)5], [CuCl5]3-, [PF5], [SbF5]",
                    "description": "Most common CN=5 geometry. Axial and equatorial positions differ.",
                },
                "square_pyramidal": {
                    "geometry": "Square Pyramidal",
                    "hybridization": "sp3d / d2sp2",
                    "bond_angles": "90° (basal), ~105° (apical-basal)",
                    "examples": "[Ni(CN)5]3-, [InCl5]2-, [VO(acac)2]",
                    "description": "Often an intermediate between TBP and octahedral; common for d8 metals with moderate field ligands.",
                },
            },
            6: {
                "default": {
                    "geometry": "Octahedral",
                    "hybridization": "sp3d2 (outer orbital) / d2sp3 (inner orbital)",
                    "bond_angles": "90°, 180°",
                    "examples": "[Fe(H2O)6]2+, [Co(NH3)6]3+, [Al(H2O)6]3+, [Ti(H2O)6]3+",
                    "description": "The most common coordination geometry. All positions equivalent.",
                },
                "trigonal_prismatic": {
                    "geometry": "Trigonal Prismatic",
                    "hybridization": "d4s (or d4sp)",
                    "bond_angles": "~60° (triangular faces), ~90° (rectangular faces)",
                    "examples": "[Mo(S2C2Ph2)3], [Re(S2C2Ph2)3], [ZrMe6]2-",
                    "description": "Rare; requires specific ligand constraints (dithiolene complexes). Distorted octahedron with D3h symmetry.",
                },
                "trigonal_antiprismatic": {
                    "geometry": "Trigonal Antiprismatic (Distorted Octahedral)",
                    "hybridization": "sp3d2 (distorted)",
                    "bond_angles": "~78-82°, ~95-100°",
                    "examples": "Some lanthanide/actinide complexes, ThI6^2-",
                    "description": "Common for larger f-block metals where ligand-ligand repulsion favors distortion.",
                },
            },
            7: {
                "default": {
                    "geometry": "Pentagonal Bipyramidal",
                    "hybridization": "sp3d3",
                    "bond_angles": "72°, 90°, 180°",
                    "examples": "[IF7], [UO2(H2O)5]2+, [NbF7]2-, [ZrF7]3-",
                    "description": "Most common CN=7 geometry. 5 equatorial + 2 axial positions.",
                },
                "capped_octahedral": {
                    "geometry": "Capped Octahedral",
                    "hybridization": "sp3d3",
                    "bond_angles": "variable",
                    "examples": "[Mo(CN)7]4-, [NbOF6]3-",
                    "description": "Octahedron with one face capped by a seventh ligand.",
                },
                "capped_trigonal_prismatic": {
                    "geometry": "Capped Trigonal Prismatic",
                    "hybridization": "sp3d3",
                    "bond_angles": "variable",
                    "examples": "[La(H2O)7]3+, some f-block complexes",
                    "description": "Trigonal prism with one rectangular face capped.",
                },
            },
            8: {
                "default": {
                    "geometry": "Square Antiprismatic",
                    "hybridization": "sp3d4 / d4sp3",
                    "bond_angles": "~99.6° (square), ~77.5° (twist angle ~45°)",
                    "examples": "[TaF8]3-, [Zr(acac)4], [ReH8]2-, [Mo(CN)8]4-",
                    "description": "Most stable CN=8 geometry for transition metals. D4d symmetry.",
                },
                "dodecahedral": {
                    "geometry": "Dodecahedral (Triangular Dodecahedron)",
                    "hybridization": "sp3d4 / d4sp3",
                    "bond_angles": "variable (two types of vertices)",
                    "examples": "[Ce(NO3)4]^-, [Zr(NO3)4]^-, [Mo(CN)8]4-(isomer)",
                    "description": "D2d symmetry. Two interpenetrating tetrahedra. Common for nitrato complexes.",
                },
                "hexagonal_bipyramidal": {
                    "geometry": "Hexagonal Bipyramidal",
                    "hybridization": "sp3d4 / f sp3d3",
                    "bond_angles": "60°, 90°, 180°",
                    "examples": "[UO2(NO3)4]^4- (approximate), some actinide complexes",
                    "description": "Hexagonal planar equatorial + 2 axial. Less common than square antiprism.",
                },
            },
        }

        # Metal-specific overrides
        self._metal_rules = {
            # d8 strong-field → square planar at CN=4
            "pt2+": {"cn4_strong": "d8_strong_field"},
            "pd2+": {"cn4_strong": "d8_strong_field"},
            "ni2+": {"cn4_strong": "d8_strong_field"},
            "au3+": {"cn4_strong": "d8_strong_field"},
            # d10 always tetrahedral at CN=4
            "zn2+": {"cn4_default": "default"},
            "cd2+": {"cn4_default": "default"},
            # Cr(III) always octahedral
            "cr3+": {"cn6_default": "default"},
            "co3+": {"cn6_default": "default"},
        }

    def _run_base(self, coordination_number: int, metal_ion: str = None,
                  ligand_type: str = None, electron_config: str = None) -> dict:
        """Predict coordination geometry based on input parameters."""
        cn = coordination_number

        if cn not in self._geometry_db:
            valid_cns = sorted(self._geometry_db.keys())
            raise ChemMCPError(
                f"Unsupported coordination number: {cn}. "
                f"Supported coordination numbers are: {valid_cns}"
            )

        cn_data = self._geometry_db[cn]

        # Determine which variant to use
        variant = "default"

        if cn == 4:
            metal_lower = (metal_ion or "").lower().replace(" ", "")
            ligand_lower = (ligand_type or "").lower()
            config_lower = (electron_config or "").lower()

            # Check for d8 strong-field → square planar
            is_d8 = any(x in config_lower for x in ["d8"])
            is_strong_field = any(x in ligand_lower for x in ["strong-field", "strong field", "cn", "nh3", "cyanide", "co", "phen", "bpy"])
            is_d8_metal = any(x in metal_lower for x in ["pt2+", "pd2+", "ni2+", "au3+", "rh+", "ir+"])

            if is_d8_metal and (is_strong_field or is_d8):
                variant = "d8_strong_field"

        elif cn == 5:
            # Default to trigonal bipyramidal
            variant = "default"

        elif cn == 6:
            # Default to octahedral
            variant = "default"

        result = dict(cn_data[variant])
        logger.info(f"CN={cn}, metal={metal_ion}, ligand={ligand_type} → {result['geometry']} ({variant})")
        return result

    def _run_text(self, query: str) -> dict:
        """Parse text query and delegate."""
        parts = query.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Query must include at least coordination number. Format: 'CN [metal] [ligand] [config]'")

        try:
            cn = int(parts[0])
        except ValueError:
            raise ChemMCPError(f"First token must be integer coordination number, got: '{parts[0]}'")

        metal_ion = parts[1] if len(parts) > 1 else None
        ligand_type = parts[2] if len(parts) > 2 else None
        electron_config = " ".join(parts[3:]) if len(parts) > 3 else None

        return self._run_base(cn, metal_ion, ligand_type, electron_config)
