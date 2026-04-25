import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

# Crystal structure database for elements and common compounds
CRYSTAL_STRUCTURES: dict = {
    # Elements
    "H":  {"structure": "hexagonal close-packed (hcp, solid H2)", "cn": None, "packing": None, "a_pm": 374, "c_pm": 603},
    "He": {"structure": "hcp (at low T, high P)", "cn": None, "packing": None, "note": "No crystal structure at STP (gas)"},
    "Li": {"structure": "body-centered cubic (bcc)", "cn": 8, "packing": 0.68, "a_pm": 351},
    "Be": {"structure": "hexagonal close-packed (hcp)", "cn": 12, "packing": 0.74, "a_pm": 229, "c_pm": 358},
    "B":  {"structure": "rhombohedral (complex)", "cn": None, "packing": None, "note": "Complex rhombohedral structure with B12 icosahedra"},
    "C":  {"structure": "diamond cubic (diamond) / hexagonal (graphite) / fullerene / nanotube", "cn": 4, "packing": 0.34, "note": "Multiple allotropes: diamond (fcc), graphite (layered hcp-like)"},
    "N":  {"structure": "cubic (solid N2, Pa=3)", "cn": None, "packing": None, "note": "Gas at STP; solid N2 has cubic Pa3 structure"},
    "O":  {"structure": "monoclinic (solid O2, γ-phase)", "cn": None, "packing": None, "note": "Gas at STP; multiple solid phases"},
    "F":  {"structure": "cubic (solid F2)", "cn": None, "packing": None, "note": "Gas at STP"},
    "Ne": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "note": "At very low temperature"},
    "Na": {"structure": "body-centered cubic (bcc)", "cn": 8, "packing": 0.68, "a_pm": 429},
    "Mg": {"structure": "hexagonal close-packed (hcp)", "cn": 12, "packing": 0.74, "a_pm": 321, "c_pm": 521},
    "Al": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 405},
    "Si": {"structure": "diamond cubic", "cn": 4, "packing": 0.34, "a_pm": 543},
    "Ge": {"structure": "diamond cubic", "cn": 4, "packing": 0.34, "a_pm": 566},
    "Sn": {"structure": "β-Sn: tetragonal (white tin) / α-Sn: diamond (gray tin, <13°C)", "cn": [6, 4], "note": "Allotropic transformation at 13.2°C"},
    "Fe": {"structure": "α-Fe: bcc (ferrite) / γ-Fe: fcc (austenite, 912-1394°C) / δ-Fe: bcc", "cn": [8, 12], "note": "Allotropes: bcc at RT, fcc at high T"},
    "Cu": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 361},
    "Ag": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 409},
    "Au": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 408},
    "Zn": {"structure": "hexagonal close-packed (hcp)", "cn": 12, "packing": 0.74, "a_pm": 266, "c_pm": 495},
    "Ca": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 558},
    "Ni": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 352},
    "Pt": {"structure": "face-centered cubic (fcc)", "cn": 12, "packing": 0.74, "a_pm": 392},
    # Ionic compounds
    "NaCl": {"structure": "rock salt (NaCl-type, fcc)", "cn": [6, 6], "madelung": 1.74756, "a_pm": 564, "space_group": "Fm-3m",
             "description": "Each Na⁺ octahedrally coordinated by 6 Cl⁻ and vice versa. Two interpenetrating fcc lattices."},
    "CsCl": {"structure": "cesium chloride (CsCl-type, simple cubic)", "cn": [8, 8], "madelung": 1.76267, "a_pm": 412, "space_group": "Pm-3m",
             "description": "Each Cs⁺ at cube center coordinated by 8 Cl⁻ at corners."},
    "ZnS": {"structure": "zinc blende (sphalerite, fcc) / wurtzite (hcp)", "cn": [4, 4], "madelung": 1.63806, "a_pm": 541,
            "description": "Tetrahedral coordination of Zn²⁺ by S²⁻. Zinc blende: fcc-based; wurtzite: hcp-based."},
    "CaF2": {"structure": "fluorite (CaF₂-type, fcc)", "cn": [8, 4], "madelung": 2.51939, "a_pm": 546,
              "description": "Each Ca²⁺ coordinated by 8 F⁻ in a cubic arrangement; each F⁻ tetrahedrally coordinated by 4 Ca²⁺."},
    "TiO2": {"structure": "rutile (tetragonal)", "cn": [6, 3], "a_pm": 459, "c_pm": 296,
              "description": "Each Ti⁴⁺ octahedrally coordinated by 6 O²⁻; each O²⁻ bonded to 3 Ti⁴⁺ in a planar arrangement."},
    "MgO": {"structure": "rock salt (NaCl-type, fcc)", "cn": [6, 6], "madelung": 1.74756, "a_pm": 421,
            "description": "Ionic oxide with rock salt structure. High melting point (2852°C)."},
    "Al2O3": {"structure": "corundum (α-Al₂O₃, trigonal)", "cn": [6, 4], "madelung": 4.1719,
               "description": "Each Al³⁺ octahedrally coordinated by 6 O²⁻; each O²⁻ surrounded by 4 Al³⁺. Very hard (Mohs 9)."},
}


@ChemMCPManager.register_tool
class GetCrystalStructure(BaseTool):
    __version__ = "0.1.0"
    name = "GetCrystalStructure"
    func_name = 'get_crystal_structure'
    description = "Query crystal structure type for elements or common ionic compounds."
    implementation_description = "Uses a built-in crystallographic database with structure type, coordination number, packing efficiency, lattice parameters, space group information for elements and common compounds."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Crystal Structure", "Solid State Chemistry", "Crystallography", "Materials Science"]
    required_envs = []

    code_input_sig = [
        ('material', 'str', 'N/A', 'Element symbol or compound formula (e.g., Fe, NaCl, CsCl, ZnS, CaF2)'),
    ]
    text_input_sig = [
        ('material', 'str', 'N/A', 'Element symbol or compound formula'),
    ]
    output_sig = [
        ('material', 'str', 'Material identifier'),
        ('crystal_system', 'str', 'Crystal system/structure type'),
        ('coordination_number', 'int/list', 'Coordination number(s)'),
        ('packing_efficiency', 'float', 'Atomic packing factor (if applicable)'),
        ('lattice_parameters', 'dict', 'Lattice constants'),
        ('description', 'str', 'Detailed structural description'),
    ]
    
        
    examples = [
        {'code_input': {'material': 'NaCl'}, 'text_input': {'material': 'NaCl'}, 'output': {'material': 'NaCl', 'crystal_system': 'rock salt (fcc)', 'coordination_number': [6, 6], 'packing_efficiency': '74%', 'lattice_parameters': {'a': 5.64}, 'description': '...'}},
        {'code_input': {'material': 'Fe'}, 'text_input': {'material': 'Fe'}, 'output': {'material': 'Fe', 'crystal_system': 'bcc', 'coordination_number': [8], 'packing_efficiency': '68%', 'lattice_parameters': {'a': 2.87}, 'description': '...'}},
        {'code_input': {'material': 'CsCl'}, 'text_input': {'material': 'CsCl'}, 'output': {'material': 'CsCl', 'crystal_system': 'CsCl (simple cubic)', 'coordination_number': [8, 8], 'packing_efficiency': '68%', 'lattice_parameters': {'a': 4.12}, 'description': '...'}},
    ]
    def _run_base(self, material: str) -> dict:
        mat = material.strip()
        
        # Try direct lookup
        if mat in CRYSTAL_STRUCTURES:
            d = CRYSTAL_STRUCTURES[mat]
        else:
            # Try capitalized version
            cap = mat[0].upper() + mat[1:] if len(mat) > 1 else mat.upper()
            if cap in CRYSTAL_STRUCTURES:
                d = CRYSTAL_STRUCTURES[cap]
            else:
                available_elements = sorted([k for k in CRYSTAL_STRUCTURES.keys() if len(k) <= 3])
                available_compounds = sorted([k for k in CRYSTAL_STRUCTURES.keys() if len(k) > 3])
                raise ChemMCPInputError(
                    f"Crystal structure data not found for '{mat}'. "
                    f"Available elements: {available_elements}. "
                    f"Available compounds: {available_compounds}."
                )

        result = {
            "material": mat,
            "crystal_system": d.get("structure"),
            "coordination_number": d.get("cn"),
            "packing_efficiency": d.get("packing"),
            "lattice_parameters": {k: v for k, v in d.items() if k.startswith('a_') or k.startswith('c_')},
            "space_group": d.get("space_group"),
            "madelung_constant": d.get("madelung"),
            "description": d.get("description", ""),
        }
        if "note" in d:
            result["note"] = d["note"]
        return result


if __name__ == "__main__":
    run_mcp_server()
