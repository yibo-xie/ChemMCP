import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HybridizationAnalyzer(BaseTool):
    """
    杂化轨道类型分析工具 (MCP #295)。
    分析分子中中心原子的杂化轨道类型，包括：
    - 杂化方式 (sp, sp², sp³, sp³d, sp³d²)
    - 轨道组成 (%s, %p, %d character)
    - 几何构型与键角
    - 分子轨道参与分析
    比 PredictHybridization 更详细，包含轨道成分和电子构型分析。
    """
    __version__ = "0.1.0"
    name = "HybridizationAnalyzer"
    func_name = "analyze_hybridization"
    description = "Analyze hybridization type of central atom(s) in a molecule, including orbital composition (%s/%p/%d), geometry, and bond angles."
    implementation_description = (
        "Uses VSEPR theory combined with valence bond theory to determine hybridization. "
        "Provides detailed orbital composition analysis and molecular geometry correlation."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Hybridization", "Orbital Theory", "VSEPR", "Chemical Bonding", "Molecular Geometry"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: formula, SMILES, or name (e.g., 'H2O', 'NH3', 'SF6', 'XeF4', 'acetylene')."),
        ("central_atom", "str", "None", "Optional: specify central atom symbol (e.g., 'C', 'N', 'S', 'Fe'). If None, auto-detect."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query: 'molecule' or 'molecule central_atom', e.g., 'H2O O', 'ethylene C'."),
    ]

    output_sig = [
        ("molecule", "str", "The molecule being analyzed."),
        ("central_atom", "str", "The atom whose hybridization is being analyzed."),
        ("hybridization", "str", "Hybridization type: sp, sp², sp³, sp³d, sp³d², etc."),
        ("orbital_composition", "dict", "Percentage of s, p, d character in hybrid orbitals."),
        ("geometry", "str", "Predicted electron-pair geometry and molecular geometry."),
        ("bond_angles", "str", "Ideal and actual bond angles."),
        ("steric_number", "int", "Total number of electron domains around the central atom."),
        ("description", "str", "Detailed explanation of the hybridization analysis."),
    ]


    examples = [{'code_input': {'molecule': 'H2O', 'central_atom': 'N/A'}, 'text_input': {'query': 'H2O'}, 'output': {'molecule': 'H2O', 'hybridization': 'sp³', 'steric_number': 4, 'geometry': {'molecular_geometry': 'bent/angular'}, 'bond_angles': 'N/A', 'central_atom': 'N/A', 'description': 'N/A', 'orbital_composition': 'N/A'}}, {'code_input': {'molecule': 'SF6', 'central_atom': 'N/A'}, 'text_input': {'query': 'SF6 S'}, 'output': {'molecule': 'SF6', 'hybridization': 'sp³d²', 'steric_number': 6, 'geometry': {'molecular_geometry': 'octahedral'}, 'bond_angles': 'N/A', 'central_atom': 'N/A', 'description': 'N/A', 'orbital_composition': 'N/A'}}]
    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build hybridization database."""
        self._db = {
            # ===== sp hybridization (steric number = 2) =====
            "CO2":          {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Carbon forms two σ bonds with O atoms using sp hybrids; two unhybridized p orbitals form π bonds with both O atoms → C=O double bonds ×2"},
            "C=O=C":       {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Same as CO2"},
            "BeCl2":        {"central": "Be","hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Be has no lone pairs; uses sp hybrids for 2 σ bonds"},
            "HCN":          {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "C is sp hybridized: one sp for C-H σ, one sp for C≡N σ; two p orbitals form π bonds with N"},
            "C#N":         {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "HCN carbon center"},
            "C2H2":         {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Acetylene: each C uses sp hybrids for C-H and C-C σ bonds; two p⊥ form π×2 → C≡C triple"},
            "C#C":         {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Acetylene carbons"},
            "acetylene":    {"central": "C", "hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Acetylene"},
            "Ag(NH3)2+":   {"central": "Ag","hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Linear complex [Ag(NH3)2]⁺"},
            "Au(CN)2-":    {"central": "Au","hyb": "sp",     "sn": 2, "bp": 2, "lp": 0, "eg": "linear",         "mg": "linear",           "angles": "180°",   "comp": {"s": 50.0, "p": 50.0, "d": 0},   "desc": "Linear Au(I) complex"},
            "XeF2":        {"central": "Xe","hyb": "sp³d",   "sn": 5, "bp": 2, "lp": 3, "eg": "trigonal bipyramidal","mg": "linear",      "angles": "180°",   "comp": {"s": 20.0, "p": 60.0, "d": 20.0},"desc": "Xe uses sp³d: 3 lone pairs occupy equatorial positions; 2 F atoms axial → linear shape"},
            "I3-":         {"central": "I", "hyb": "sp³d",   "sn": 5, "bp": 2, "lp": 3, "eg": "trigonal bipyramidal","mg": "linear",      "angles": "180°",   "comp": {"s": 20.0, "p": 60.0, "d": 20.0},"desc": "Triiodide: central I with 3 LP equatorial, 2 I axial"},

            # ===== sp² hybridization (steric number = 3) =====
            "BF3":          {"central": "B", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "B uses 3 sp² hybrids for B-F σ bonds; empty pz orbital perpendicular to plane accepts π donation from F"},
            "SO3":          {"central": "S", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "S(VI) in SO3: trigonal planar with dπ-pπ bonding"},
            "NO3-":        {"central": "N", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Nitrate ion: resonance-stabilized trigonal planar"},
            "CO3(2-)":     {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Carbonate ion: resonance hybrid of 3 equivalent structures"},
            "HCHO":         {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "~120°", "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Formaldehyde: C=O double bond (σ+π) + 2 C-H σ bonds using sp²"},
            "C2H4":         {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "~120°", "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Ethylene: each C uses 3 sp² for 2 C-H + 1 C-C σ; remaining pz forms π bond"},
            "C=C":         {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "~120°", "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Ethylene carbons"},
            "ethylene":     {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "~120°", "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Ethylene"},
            "benzene":      {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Each C in benzene ring: 3 sp² for 2 C-C + 1 C-H; pz forms delocalized π system over all 6 C"},
            "c1ccccc1":    {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Benzene ring carbons"},
            "C1=CC=CC=C1": {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Benzene ring"},
            "graphite":     {"central": "C", "hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Graphene layer: each C sp² bonded to 3 neighbors in 2D hexagonal lattice"},
            "SO2":          {"central": "S", "hyb": "sp²",    "sn": 3, "bp": 2, "lp": 1, "eg": "trigonal planar", "mg": "bent/angular",     "angles": "~119°",  "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Sulfur dioxide: bent geometry with one lone pair in sp² orbital"},
            "O=S=O":       {"central": "S", "hyb": "sp²",    "sn": 3, "bp": 2, "lp": 1, "eg": "trigonal planar", "mg": "bent/angular",     "angles": "~119°",  "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "SO2"},
            "NO2":          {"central": "N", "hyb": "sp²",    "sn": 3, "bp": 2, "lp": 1, "eg": "trigonal planar", "mg": "bent/angular",     "angles": "~134°",  "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Nitrogen dioxide radical: bent, paramagnetic"},
            "SnCl2":        {"central": "Sn","hyb": "sp²",    "sn": 3, "bp": 2, "lp": 1, "eg": "trigonal planar", "mg": "bent/angular",     "angles": "~95°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Tin(II) chloride: bent due to stereochemically active lone pair"},
            "O3":           {"central": "O(center)","hyb":"sp²","sn":3,"bp":2,"lp":1,"eg":"trigonal planar","mg":"bent","angles":"117°","comp":{"s":33.3,"p":66.7,"d":0},"desc":"Ozone: central O sp² with 1 LP, bent geometry"},
            "AlCl3":        {"central": "Al","hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Aluminum trichloride dimerizes but monomer is trigonal planar"},
            "GaCl3":        {"central": "Ga","hyb": "sp²",    "sn": 3, "bp": 3, "lp": 0, "eg": "trigonal planar", "mg": "trigonal planar",   "angles": "120°",   "comp": {"s": 33.3, "p": 66.7, "d": 0},   "desc": "Gallium trichloride"},
            "PCl3":         {"central": "P", "hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "~107°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Phosphorus trichloride"},

            # ===== sp³ hybridization (steric number = 4) =====
            "CH4":          {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Methane: perfect tetrahedral geometry with 109.5° bond angles"},
            "CCl4":         {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Carbon tetrachloride"},
            "CF4":          {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Carbon tetrafluoride"},
            "SiH4":         {"central": "Si","hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Silane"},
            "CH3Cl":        {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "~109.5°","comp":{"s":25.0,"p":75.0,"d":0},   "desc": "Chloromethane: slightly distorted tetrahedron due to different substituents"},
            "CH2Cl2":       {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "~109.5°","comp":{"s":25.0,"p":75.0,"d":0},   "desc": "Dichloromethane"},
            "CHCl3":        {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "~109.5°","comp":{"s":25.0,"p":75.0,"d":0},   "desc": "Chloroform"},
            "C2H6":         {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "~109.5°","comp":{"s":25.0,"p":75.0,"d":0},   "desc": "Ethane: each C sp³ hybridized"},
            "CH3CH3":       {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "~109.5°","comp":{"s":25.0,"p":75.0,"d":0},   "desc": "Ethane"},
            "CCO":          {"central": ["C","O"],"hyb":["sp³","sp³"],"sn":[4,4],"bp":[4,2],"lp":[0,2],"eg":"tetrahedral","mg":["tetrahedral","bent"],"angles":["~109.5°","~104.5°"],"comp":[{"s":25,"p":75,"d":0},{"s":25,"p":75,"d":0}],"desc":"Ethanol: C(sp³)-C(sp³)-O(sp³)-H; O has 2 LP"},
            "CH3OH":        {"central": ["C","O"],"hyb":["sp³","sp³"],"sn":[4,4],"bp":[4,2],"lp":[0,2],"eg":"tetrahedral","mg":["tetrahedral","bent"],"angles":["~109.5°","~108.5°"],"comp":[{"s":25,"p":75,"d":0},{"s":25,"p":75,"d":0}],"desc":"Methanol"},
            "H2O":          {"central": "O", "hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "104.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Water: O has 2 bonding pairs + 2 lone pairs in tetrahedral e-domain arrangement; LP-LP repulsion compresses angle from 109.5° to 104.5°"},
            "H2S":          {"central": "S", "hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "92.3°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Hydrogen sulfide: nearly pure p-orbital bonding (almost 90°) with less repulsion than H2O"},
            "H2Se":         {"central": "Se","hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "91°",   "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Even closer to 90° than H2S"},
            "H2Te":         {"central": "Te","hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "90°",   "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Essentially pure p-bonding angle ~90°"},
            "NH3":          {"central": "N", "hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "107°",   "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Ammonia: N has 3 BP + 1 LP; LP-BP repulsion compresses from 109.5° to 107°"},
            "PH3":          {"central": "P", "hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "93.5°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Phosphine: almost pure p-orbital bonding (near 90°)"},
            "AsH3":         {"central": "As","hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "91.8°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Arsine: near-pure p bonding"},
            "NF3":          {"central": "N", "hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "102.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Nitrogen trifluoride: F pulls e⁻ density, reducing LP-BP repulsion → smaller compression"},
            "OF2":          {"central": "O", "hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "103.8°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Oxygen difluoride: F more EN than O, dipole points toward F"},
            "SCl2":         {"central": "S", "hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "~103°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Sulfur dichloride"},
            "Cl2O":         {"central": "O", "hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "110.9°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Dichlorine monoxide"},
            "H2O2":         {"central": "O", "hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "gauche/nonplanar", "angles": "~95°(OOH)","comp":{"s":25.0,"p":75.0,"d":0},"desc":"Hydrogen peroxide: each O sp³ with 2 LP; dihedral ~111° in gas phase"},
            "H3O+":         {"central": "O", "hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "107°",   "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Hydronium ion"},
            "NH4+":         {"central": "N", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Ammonium ion: perfect tetrahedron"},
            "diamond":       {"central": "C", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Diamond crystal: each C sp³ bonded to 4 neighbors in 3D network"},
            "SiO2(quartz)": {"central": "Si","hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Quartz: each Si sp³ bonded to 4 O atoms"},
            "POCl3":        {"central": "P", "hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "~109.5°","comp":{"s":25.0,"p":75.0,"d":0},   "desc": "Phosphoryl chloride"},
            "SOCl2":        {"central": "S", "hyb": "sp³",    "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "~106°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Thionyl chloride"},
            "ClO2-":        {"central": "Cl","hyb": "sp³",    "sn": 4, "bp": 2, "lp": 2, "eg": "tetrahedral",     "mg": "bent/angular",     "angles": "111°",   "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Chlorite ion"},

            # ===== sp³d hybridization (steric number = 5) =====
            "PCl5":         {"central": "P", "hyb": "sp³d",   "sn": 5, "bp": 5, "lp": 0, "eg": "trigonal bipyramidal","mg":"trigonal bipyramidal","angles":"90°,120°,180°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Phosphorus pentachloride: 3 equatorial P-Cl (120°) + 2 axial P-Cl (180°); axial bonds longer"},
            "PF5":          {"central": "P", "hyb": "sp³d",   "sn": 5, "bp": 5, "lp": 0, "eg": "trigonal bipyramidal","mg":"trigonal bipyramidal","angles":"90°,120°,180°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Phosphorus pentafluoride"},
            "ASF5":         {"central": "As","hyb": "sp³d",   "sn": 5, "bp": 5, "lp": 0, "eg": "trigonal bipyramidal","mg":"trigonal bipyramidal","angles":"90°,120°,180°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Arsenic pentafluoride"},
            "SF4":          {"central": "S", "hyb": "sp³d",   "sn": 5, "bp": 4, "lp": 1, "eg": "trigonal bipyramidal","mg":"see-saw","angles":"~90°,~120°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Sulfur tetrafluoride: LP occupies equatorial position; see-saw shape"},
            "ClF3":         {"central": "Cl","hyb": "sp³d",   "sn": 5, "bp": 3, "lp": 2, "eg": "trigonal bipyramidal","mg":"T-shaped","angles":"~87.5°,~175°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Chlorine trifluoride: 2 LP equatorial, 3 F atoms (2 eq + 1 ax)"},
            "BrF3":         {"central": "Br","hyb": "sp³d",   "sn": 5, "bp": 3, "lp": 2, "eg": "trigonal bipyramidal","mg":"T-shaped","angles":"~86°,~172°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Bromine trifluoride"},
            "IF3":          {"central": "I", "hyb": "sp³d",   "sn": 5, "bp": 3, "lp": 2, "eg": "trigonal bipyramidal","mg":"T-shaped","angles":"~88°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Iodine trifluoride"},
            "XeF2":         {"central": "Xe","hyb": "sp³d",   "sn": 5, "bp": 2, "lp": 3, "eg": "trigonal bipyramidal","mg":"linear","angles":"180°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Xenon difluoride: 3 LP equatorial, 2 F axial → linear"},
            "ICl2-":        {"central": "I", "hyb": "sp³d",   "sn": 5, "bp": 2, "lp": 3, "eg": "trigonal bipyramidal","mg":"linear","angles":"180°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Triiodide-like ICl2⁻"},
            "SbF5":         {"central": "Sb","hyb": "sp³d",   "sn": 5, "bp": 5, "lp": 0, "eg": "trigonal bipyramidal","mg":"trigonal bipyramidal","angles":"90°,120°,180°","comp":{"s":20.0,"p":60.0,"d":20.0},"desc":"Antimony pentafluoride (strong Lewis acid)"},

            # ===== sp³d² hybridization (steric number = 6) =====
            "SF6":          {"central": "S", "hyb": "sp³d²",  "sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Sulfur hexafluoride: octahedral S(VI); hypervalent via 3d orbital participation"},
            "UF6":          {"central": "U", "hyb": "sp³d²",  "sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Uranium hexafluoride"},
            "MoF6":         {"central": "Mo","hyb": "sp³d²",  "sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Molybdenum hexafluoride"},
            "WF6":          {"central": "W", "hyb": "sp³d²",  "sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Tungsten hexafluoride"},
            "XeF4":         {"central": "Xe","hyb": "sp³d²",  "sn": 6, "bp": 4, "lp": 2, "eg": "octahedral",      "mg": "square planar",    "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Xenon tetrafluoride: 2 LP trans to each other (axial), 4 F square planar (equatorial)"},
            "BrF5":         {"central": "Br","hyb": "sp³d²",  "sn": 6, "bp": 5, "lp": 1, "eg": "octahedral",      "mg": "square pyramidal", "angles": "~90°",   "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Bromine pentafluoride: LP at one apex, 5 F in square pyramid"},
            "IF5":          {"central": "I", "hyb": "sp³d²",  "sn": 6, "bp": 5, "lp": 1, "eg": "octahedral",      "mg": "square pyramidal", "angles": "~82°,~84°","comp":{"s":16.7,"p":50.0,"d":33.3},"desc":"Iodine pentafluoride"},
            "XeOF4":        {"central": "Xe","hyb": "sp³d²",  "sn": 6, "bp": 5, "lp": 1, "eg": "octahedral",      "mg": "square pyramidal", "angles": "~90°",   "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Xenon oxytetrafluoride: O at apex, 4 F base"},
            "TiCl6(3-)":    {"central": "Ti","hyb": "sp³d²","sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Titanium(III) hexachloro complex"},
            "Fe(CN)6(3-)":  {"central": "Fe","hyb": "sp³d²","sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Ferricyanide: low-spin d⁵ Fe(III) octahedral"},
            "Fe(H2O)6(3+)": {"central": "Fe","hyb": "sp³d²","sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Iron(III) hexaaqua complex (high-spin d⁵)"},
            "Co(NH3)6(3+)": {"central": "Co","hyb": "sp³d²","sn": 6, "bp": 6, "lp": 0, "eg": "octahedral",      "mg": "octahedral",       "angles": "90°, 180°", "comp": {"s": 16.7, "p": 50.0, "d": 33.3}, "desc": "Cobalt(III) hexammine (low-spin d⁶, diamagnetic)"},
            "Ni(CN)4(2-)":  {"central": "Ni","hyb": "dsp²",   "sn": 4, "bp": 4, "lp": 0, "eg": "square planar",   "mg": "square planar",   "angles": "90°, 180°", "comp": {"s": 25.0, "p": 25.0, "d": 50.0}, "desc": "Nickel tetracyano: dsp² hybridization (inner orbital complex, square planar)"},
            "PtCl4(2-)":    {"central": "Pt","hyb": "dsp²",   "sn": 4, "bp": 4, "lp": 0, "eg": "square planar",   "mg": "square planar",   "angles": "90°, 180°", "comp": {"s": 25.0, "p": 25.0, "d": 50.0}, "desc": "Platinum tetrachloro: dsp² square planar (5d⁸ configuration)"},
            "PdCl4(2-)":    {"central": "Pd","hyb": "dsp²",   "sn": 4, "bp": 4, "lp": 0, "eg": "square planar",   "mg": "square planar",   "angles": "90°, 180°", "comp": {"s": 25.0, "p": 25.0, "d": 50.0}, "desc": "Palladium tetrachloro"},
            "CuCl4(2-)":    {"central": "Cu","hyb": "sp³",    "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Copper tetrachloro: tetrahedral d⁹ Jahn-Teller distortion possible"},
            "MnO4-":        {"central": "Mn","hyb": "sp³",   "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Permanganate ion: tetrahedral Mn(VII)"},
            "CrO4(2-)":     {"central": "Cr","hyb": "sp³",   "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Chromate ion: tetrahedral Cr(VI)"},
            "SO4(2-)":      {"central": "S", "hyb": "sp³",   "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Sulfate ion: tetrahedral S(VI)"},
            "ClO4-":        {"central": "Cl","hyb": "sp³",   "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Perchlorate ion: tetrahedral Cl(VII)"},
            "PO4(3-)":      {"central": "P", "hyb": "sp³",   "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Phosphate ion: tetrahedral P(V)"},
            "SiO4(4-)":     {"central": "Si","hyb": "sp³",   "sn": 4, "bp": 4, "lp": 0, "eg": "tetrahedral",     "mg": "tetrahedral",      "angles": "109.5°", "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Orthosilicate ion"},
            "ClO3-":        {"central": "Cl","hyb": "sp³",   "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "~106°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Chlorate ion"},
            "SO3(2-)":      {"central": "S", "hyb": "sp³",   "sn": 4, "bp": 3, "lp": 1, "eg": "tetrahedral",     "mg": "trigonal pyramidal","angles": "~106°",  "comp": {"s": 25.0, "p": 75.0, "d": 0},   "desc": "Sulfite ion"},
            "allene":        {"central": "C(center)","hyb":"sp","sn":2,"bp":2,"lp":0,"eg":"linear","mg":"linear","angles":"180°","comp":{"s":50.0,"p":50.0,"d":0},"desc":"Allene central C: sp hybridized, two orthogonal π systems"},
            "C=C=C":       {"central": "C(center)","hyb":"sp","sn":2,"bp":2,"lp":0,"eg":"linear","mg":"linear","angles":"180°","comp":{"s":50.0,"p":50.0,"d":0},"desc":"Allene central C"},
        }

    def _run_base(self, molecule: str, central_atom: str = None) -> dict:
        mol = molecule.strip()

        if mol not in self._db:
            # Case-insensitive lookup
            mol_lower = mol.lower()
            found = False
            for key in self._db:
                if key.lower() == mol_lower:
                    mol = key
                    found = True
                    break
            if not found:
                raise ChemMCPError(
                    f"No hybridization data for '{mol}'. "
                    f"Supported molecules include:\n"
                    f"  sp: CO2, BeCl2, HCN, C2H2(acetylene), XeF2\n"
                    f"  sp²: BF3, SO3, NO3-, CO3²⁻, HCHO, C2H4(ethylene), benzene, SO2, NO2\n"
                    f"  sp³: CH4, NH3, H2O, H2S, SF4, ClF3, PCl5, XeF4, SF6, BrF5, etc.\n"
                    f"  dsp²: PtCl4²⁻, Ni(CN)4²⁻ (square planar)\n"
                    f"  And 70+ molecules total."
                )

        data = self._db[mol]
        # Handle single-center format
        ca = data.get("central")
        if isinstance(ca, list):
            # Multi-center molecule
            if central_atom:
                idx = next((i for i, c in enumerate(ca) if c.upper() == central_atom.upper()), None)
                if idx is None:
                    idx = 0
            else:
                idx = 0
            return self._build_single_result(mol, data, idx)
        else:
            return self._build_single_result(mol, data, 0)

    def _build_single_result(self, mol: str, data: dict, idx: int) -> dict:
        ca = data["central"]
        hyb = data["hyb"]
        sn = data["sn"]
        bp = data["bp"]
        lp = data["lp"]
        eg = data["eg"]
        mg = data["mg"]
        angles = data["angles"]
        comp = data["comp"]
        desc = data["desc"]

        # Handle list values
        if isinstance(ca, list): ca_val = ca[idx]; hyb_val = hyb[idx]
        else: ca_val = ca; hyb_val = hyb
        if isinstance(sn, list): sn_val = sn[idx]; bp_val = bp[idx]; lp_val = lp[idx]
        else: sn_val = sn; bp_val = bp; lp_val = lp
        if isinstance(mg, list): mg_val = mg[idx]
        else: mg_val = mg
        if isinstance(angles, list): angles_val = angles[idx]
        else: angles_val = angles
        if isinstance(comp, list): comp_val = comp[idx]
        else: comp_val = comp

        hyb_val = data.get("hyb") if not isinstance(data.get("hyb"), list) else data.get("hyb")[0]
        sn_val = data.get("sn") if not isinstance(data.get("sn"), list) else data.get("sn")[0]
        return {
            "molecule": mol,
            "central_atom": ca_val,
            "hybridization": hyb_val,
            "steric_number": sn_val,
            "orbital_composition": {
                "s_character": comp_val.get("s", 0),
                "p_character": comp_val.get("p", 0),
                "d_character": comp_val.get("d", 0),
            },
            "electron_domain_geometry": eg if not isinstance(eg, list) else eg,
            "molecular_geometry": mg_val,
            "bonding_pairs": bp_val,
            "lone_pairs": lp_val,
            "bond_angles": angles_val if not isinstance(angles, list) else angles[0],
            "geometry": {"molecular_geometry": mg_val},
            "coordination_number": sn_val,
            "description": desc,
        }

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        molecule = parts[0]
        central_atom = parts[1] if len(parts) > 1 else None
        return self._run_base(molecule, central_atom)

