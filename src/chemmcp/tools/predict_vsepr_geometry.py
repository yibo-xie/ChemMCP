import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

# VSEPR geometry mapping: steric_number -> (electron_geometry, molecular_geometry, bond_angles, hybridization)
VSEPR_GEOMETRIES = {
    # (bonding_pairs, lone_pairs) -> info
    (2, 0): {"eg": "linear", "mg": "linear", "angles": "180°", "hybrid": "sp", "example": "CO2, BeCl2"},
    (3, 0): {"eg": "trigonal planar", "mg": "trigonal planar", "angles": "120°", "hybrid": "sp²", "example": "BF3, SO3"},
    (2, 1): {"eg": "trigonal planar", "mg": "bent/angular", "angles": "<120° (~119°)", "hybrid": "sp²", "example": "SO2, SnCl2"},
    (4, 0): {"eg": "tetrahedral", "mg": "tetrahedral", "angles": "109.5°", "hybrid": "sp³", "example": "CH4, CCl4"},
    (3, 1): {"eg": "tetrahedral", "mg": "trigonal pyramidal", "angles": "<109.5° (~107°)", "hybrid": "sp³", "example": "NH3, PCl3"},
    (2, 2): {"eg": "tetrahedral", "mg": "bent/angular", "angles": "<109.5° (~104.5°)", "hybrid": "sp³", "example": "H2O, H2S"},
    (5, 0): {"eg": "trigonal bipyramidal", "mg": "trigonal bipyramidal", "angles": "90°, 120°, 180°", "hybrid": "sp³d", "example": "PCl5"},
    (4, 1): {"eg": "trigonal bipyramidal", "mg": "see-saw", "angles": "~90°, ~120°", "hybrid": "sp³d", "example": "SF4"},
    (3, 2): {"eg": "trigonal bipyramidal", "mg": "T-shaped", "angles": "~90°, 180°", "hybrid": "sp³d", "example": "ClF3"},
    (2, 3): {"eg": "trigonal bipyramidal", "mg": "linear", "angles": "180°", "hybrid": "sp³d", "example": "XeF2, I3⁻"},
    (6, 0): {"eg": "octahedral", "mg": "octahedral", "angles": "90°, 180°", "hybrid": "sp³d²", "example": "SF6"},
    (5, 1): {"eg": "octahedral", "mg": "square pyramidal", "angles": "~90°", "hybrid": "sp³d²", "example": "BrF5"},
    (4, 2): {"eg": "octahedral", "mg": "square planar", "angles": "90°, 180°", "hybrid": "sp³d²", "example": "XeF4"},
}

# SMILES-based molecule geometry database
MOLECULE_VSEPR = {
    # SMILES -> (bp, lp) for central atom analysis
    "O": ("H2O", 2, 2),  # water
    "C(=O)=O": ("CO2", 2, 0),  # CO2
    "N": ("NH3", 3, 1),  # ammonia
    "C": ("CH4", 4, 0),  # methane
    "S(=O)(=O)": ("SO2", 2, 1),  # SO2 (simplified)
    "S(F)(F)(F)(F)(F)F": ("SF6", 6, 0),
    "P(Cl)(Cl)(Cl)(Cl)Cl": ("PCl5", 5, 0),
    "B(F)(F)F": ("BF3", 3, 0),
    "N#N": ("N2", None, None),  # diatomic
    "O=O": ("O2", None, None),
    "Cl": ("HCl", None, None),
    "C(C)(C)C": ("CCl4" if False else "neopentane-like", 4, 0),
    "N=C=O": ("NCO/HCN", 2, 0),  # linear
}


@ChemMCPManager.register_tool
class PredictVseprGeometry(BaseTool):
    __version__ = "0.1.0"
    name = "PredictVseprGeometry"
    func_name = 'predict_vsepr_geometry'
    description = "Predict molecular geometry using VSEPR theory based on bonding pairs and lone pairs around a central atom."
    implementation_description = "Uses VSEPR (Valence Shell Electron Pair Repulsion) theory to predict electron geometry, molecular geometry, bond angles, and hybridization from steric number (bonding pairs + lone pairs)."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["VSEPR", "Molecular Geometry", "Hybridization", "Chemical Bonding"]
    required_envs = []

    code_input_sig = [
        ('molecule', 'str', 'N/A', 'Molecular formula or SMILES string (e.g., H2O, NH3, SF6, C=O)'),
        ('bonding_pairs', 'int', 'N/A', 'Optional: number of bonding pairs around central atom'),
        ('lone_pairs', 'int', 'N/A', 'Optional: number of lone pairs on central atom'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Molecule name or formula, e.g., \"water\" or \"H2O\" or \"SF6\"'),
    ]
    output_sig = [
        ('molecule', 'str', 'Molecule identifier'),
        ('steric_number', 'int', 'Total number of electron domains (BP + LP)'),
        ('electron_geometry', 'str', 'Electron pair geometry'),
        ('molecular_geometry', 'str', 'Molecular geometry (atom positions only)'),
        ('bond_angles', 'str', 'Approximate bond angles'),
        ('hybridization', 'str', 'Predicted hybridization of central atom'),
        ('description', 'str', 'Detailed geometry description'),
    ]
    
        
    
    examples = [
        {'code_input': {'molecule': 'H2O', 'bonding_pairs': 2, 'lone_pairs': 2}, 'text_input': {'query': 'H2O'}, 'output': {'molecule': 'H2O', 'steric_number': 4, 'electron_geometry': 'tetrahedral', 'molecular_geometry': 'bent/angular', 'bond_angles': '104.5', 'hybridization': 'sp3', 'description': '...'}},
        {'code_input': {'molecule': 'SF6', 'bonding_pairs': 6, 'lone_pairs': 0}, 'text_input': {'query': 'SF6'}, 'output': {'molecule': 'SF6', 'steric_number': 6, 'electron_geometry': 'octahedral', 'molecular_geometry': 'octahedral', 'bond_angles': '90', 'hybridization': 'sp3d2', 'description': '...'}},
    ]
    def _run_base(self, molecule: str, bonding_pairs: int = None, lone_pairs: int = None) -> dict:
        mol = molecule.strip()
        
        # If bp/lp explicitly provided, use them directly
        if bonding_pairs is not None and lone_pairs is not None:
            bp, lp = int(bonding_pairs), int(lone_pairs)
            sn = bp + lp
        else:
            # Try to look up from known molecules
            mol_upper = mol.upper().replace(" ", "")
            lookup = {
                "H2O": (2, 2), "WATER": (2, 2), "H2S": (2, 2), "H2SE": (2, 2),
                "NH3": (3, 1), "AMMONIA": (3, 1), "PH3": (3, 1), "PHOSPHINE": (3, 1),
                "CH4": (4, 0), "METHANE": (4, 0), "CF4": (4, 0), "CCl4": (4, 0),
                "BF3": (3, 0), "SO3": (3, 0), "CO2": (2, 0), "CO2(LINEAR)": (2, 0),
                "BECL2": (2, 0), "Hgcl2": (2, 0), "CS2": (2, 0),
                "SO2": (2, 1), "NO2": (2, 1), "SNCL2": (2, 1),
                "PCl5": (5, 0), "PF5": (5, 0), "ASF5": (5, 0),
                "SF6": (6, 0), "SEF6": (6, 0), "TEF6": (6, 0),
                "XEF4": (4, 2), "BRF5": (5, 1), "XEOF4": (4, 2),
                "XEF2": (2, 3), "IF3": (3, 2), "IF5": (5, 1),
                "ICL3": (3, 2), "ICL5": (5, 1),
                "SF4": (4, 1), "CLF3": (3, 2),
                "HCN": (2, 0), "BEH2": (2, 0),
                "XEO3": (3, 1), "CLF2": (2, 2),
            }
            result = lookup.get(mol_upper)
            if result is None:
                avail = sorted([k for k in lookup.keys() if len(k) <= 10])
                raise ChemMCPInputError(
                    f"Cannot determine VSEPR for '{mol}'. Please provide bonding_pairs and lone_pairs explicitly, "
                    f"or use one of the known molecules: {avail[:30]}. "
                    f"Example: molecule='H2O' or molecule='custom', bonding_pairs=4, lone_pairs=0"
                )
            bp, lp = result
            sn = bp + lp

        key = (bp, lp)
        if key not in VSEPR_GEOMETRIES:
            raise ChemMCPInputError(
                f"No VSEPR data for bonding_pairs={bp}, lone_pairs={lp}. "
                f"Supported combinations cover steric numbers 2-6."
            )

        g = VSEPR_GEOMETRIES[key]
        return {
            "molecule": mol,
            "bonding_pairs": bp,
            "lone_pairs": lp,
            "steric_number": sn,
            "electron_geometry": g["eg"],
            "molecular_geometry": g["mg"],
            "bond_angles": g["angles"],
            "hybridization": g["hybrid"],
            "example_compounds": g["example"],
            "description": (
                f"With {sn} electron domains ({bp} bonding pairs + {lp} lone pairs), "
                f"the {mol} adopts {g['eg']} electron geometry and {g['mg']} molecular geometry. "
                f"Bond angles are approximately {g['angles']}. The central atom is {g['hybrid']} hybridized."
            ),
        }


if __name__ == "__main__":
    run_mcp_server()
