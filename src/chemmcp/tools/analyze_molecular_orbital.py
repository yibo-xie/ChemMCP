import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

# Molecular orbital data for diatomic molecules (period 2 homonuclear + key heteronuclear)
# Based on MO theory: σ1s < σ*1s < σ2s < σ*2s < π2p = π2p < σ2p < π*2p = π*2p < σ*2p
MO_DATA: dict = {
    "H2": {
        "formula": "H₂", "electrons": 2,
        "configuration": "(σ1s)²",
        "bond_order": 1,
        "magnetic": "diamagnetic",
        "stability": "stable",
        "description": "Simplest molecule. Both electrons in bonding σ1s orbital.",
        "orbitals": [
            {"name": "σ1s", "type": "bonding", "electrons": 2, "energy": "low"},
        ],
    },
    "He2": {
        "formula": "He₂ (hypothetical)", "electrons": 4,
        "configuration": "(σ1s)²(σ*1s)²",
        "bond_order": 0,
        "magnetic": "diamagnetic",
        "stability": "does not exist (BO=0)",
        "description": "Bonding and antibonding cancel out; He2 is not a stable molecule.",
        "orbitals": [
            {"name": "σ1s", "type": "bonding", "electrons": 2},
            {"name": "σ*1s", "type": "antibonding", "electrons": 2},
        ],
    },
    "Li2": {
        "formula": "Li₂", "electrons": 6,
        "configuration": "(σ1s)²(σ*1s)²(σ2s)²",
        "bond_order": 1,
        "magnetic": "diamagnetic",
        "stability": "stable (gas phase)",
        "description": "Valence electrons fill σ2s bonding only. Weak bond.",
        "orbitals": [
            {"name": "σ2s (valence)", "type": "bonding", "electrons": 2, "note": "core (1s) orbitals filled and cancel"},
        ],
    },
    "B2": {
        "formula": "B₂", "electrons": 10,
        "configuration": "(σ2s)²(σ*2s)²(π2p_x)¹(π2p_y)¹",
        "bond_order": 1,
        "magnetic": "paramagnetic (2 unpaired e⁻)",
        "stability": "stable",
        "description": "π2p orbitals are lower in energy than σ2p for B2 (Z<8). Two unpaired electrons make it paramagnetic.",
        "orbitals": [
            {"name": "σ2s", "type": "bonding", "electrons": 2},
            {"name": "σ*2s", "type": "antibonding", "electrons": 2},
            {"name": "π2p_x = π2p_y", "type": "bonding", "electrons": 2, "unpaired": 2},
        ],
    },
    "C2": {
        "formula": "C₂", "electrons": 12,
        "configuration": "(σ2s)²(σ*2s)²(π2p)⁴",
        "bond_order": 2,
        "magnetic": "diamagnetic",
        "stability": "stable",
        "description": "All four valence p-electrons in π bonding orbitals. Double bond character.",
        "orbitals": [
            {"name": "π2p (x+y)", "type": "bonding", "electrons": 4},
        ],
    },
    "N2": {
        "formula": "N₂", "electrons": 14,
        "configuration": "(σ2s)²(σ*2s)²(π2p)⁴(σ2p)²",
        "bond_order": 3,
        "magnetic": "diamagnetic",
        "stability": "very stable (triple bond)",
        "description": "Triple bond with all bonding orbitals filled. One of the strongest bonds known (945 kJ/mol). Very inert.",
        "orbitals": [
            {"name": "σ2s", "e": 2}, {"name": "σ*2s", "e": 2},
            {"name": "π2p (x,y)", "type": "bonding", "e": 4},
            {"name": "σ2p_z", "type": "bonding", "e": 2},
        ],
    },
    "O2": {
        "formula": "O₂", "electrons": 16,
        "configuration": "(σ2s)²(σ*2s)²(σ2p)²(π2p)⁴(π*2p)²",
        "bond_order": 2,
        "magnetic": "paramagnetic (2 unpaired e⁻ in degenerate π* orbitals)",
        "stability": "stable",
        "description": "Two unpaired electrons in π* antibonding orbitals explain O2's paramagnetism (attracted to magnetic field). Bond order 2 (double bond).",
        "orbitals": [
            {"name": "σ2p_z", "type": "bonding", "e": 2},
            {"name": "π2p (x,y)", "type": "bonding", "e": 4},
            {"name": "π*2p (x,y)", "type": "antibonding", "e": 2, "unpaired": 2},
        ],
    },
    "F2": {
        "formula": "F₂", "electrons": 18,
        "configuration": "(σ2s)²(σ*2s)²(σ2p)²(π2p)⁴(π*2p)⁴",
        "bond_order": 1,
        "magnetic": "diamagnetic",
        "stability": "stable but weak bond (F-F is weak due to lone pair repulsion)",
        "description": "Single bond. Weak F-F bond (158 kJ/mol) due to repulsion between small F atoms' lone pairs.",
        "orbitals": [
            {"name": "σ2p_z", "e": 2}, {"name": "π2p", "e": 4},
            {"name": "π*2p", "type": "antibonding", "e": 4},
        ],
    },
    "Ne2": {
        "formula": "Ne₂ (hypothetical)", "electrons": 20,
        "configuration": "(σ2s)²(σ*2s)²(σ2p)²(π2p)⁴(π*2p)⁴(σ*2p)²",
        "bond_order": 0,
        "magnetic": "diamagnetic",
        "stability": "does not exist (BO=0)",
        "description": "All bonding and antibonding orbitals cancel. Noble gases do not form stable Ne2.",
        "orbitals": [],
    },
    # Heteronuclear diatomics
    "CO": {
        "formula": "CO", "electrons": 14,
        "configuration": "similar to N₂: (σ2s)²(σ*2s)²(σ)²(π)⁴(σ_nonbonding)² or (3σ)²(4σ)²(1π)⁴(5σ)²",
        "bond_order": 3,
        "magnetic": "diamagnetic",
        "stability": "very stable (triple bond), toxic gas",
        "description": "Isoelectronic with N₂ (same electron count). Triple bond between C and O. Small dipole moment (Cδ-—Oδ+) because O donates more electron density than expected.",
        "isoelectronic_with": "N₂",
    },
    "NO": {
        "formula": "NO", "electrons": 15,
        "configuration": "similar to O₂+ : (σ2s)²(σ*2s)²(σ)²(π)⁴(π*)¹",
        "bond_order": 2.5,
        "magnetic": "paramagnetic (1 unpaired e⁻)",
        "stability": "stable free radical",
        "description": "Odd-electron species (15 electrons). Bond order 2.5. Paramagnetic free radical important in biology (signaling molecule).",
        "isoelectronic_with": "O₂⁺",
    },
    "HF": {
        "formula": "HF", "electrons": 10,
        "configuration": "(σ1s)²(σ2s)²(nonbonding on F: 3 pairs)",
        "bond_order": 1,
        "magnetic": "diamagnetic",
        "stability": "stable, polar molecule",
        "description": "Polar covalent bond. Large electronegativity difference (F=3.98, H=2.20). Strong hydrogen bonding in liquid/solid state.",
        "dipole_direction": "Hδ+—Fδ−",
    },
    "CN": {
        "formula": "CN (cyanide radical)", "electrons": 13,
        "configuration": "similar to N₂+: (σ2s)²(σ*2s)²(π)⁴(σ)¹",
        "bond_order": 2.5,
        "magnetic": "paramagnetic (1 unpaired e⁻)",
        "stability": "reactive radical",
        "description": "13-electron radical. CN⁻ ion has 14 electrons (like N₂/CO) with BO=3.",
    },
}


@ChemMCPManager.register_tool
class AnalyzeMolecularOrbital(BaseTool):
    __version__ = "0.1.0"
    name = "AnalyzeMolecularOrbital"
    func_name = 'analyze_molecular_orbital'
    description = "Analyze molecular orbital diagram for simple diatomic molecules: bond order, magnetic property, MO configuration."
    implementation_description = "Implements molecular orbital theory for period 2 homonuclear diatomic molecules (H2 through Ne2) and common heteronuclear diatomics (CO, NO, HF, CN). Returns electron configuration, bond order, magnetic behavior, and stability analysis."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Molecular Orbital Theory", "Quantum Chemistry", "Bond Order", "Magnetism"]
    required_envs = []

    code_input_sig = [
        ('molecule', 'str', 'N/A', 'Diatomic molecule formula (e.g., N2, O2, CO, NO, HF)'),
    ]
    text_input_sig = [
        ('molecule', 'str', 'N/A', 'Diatomic molecule formula'),
    ]
    output_sig = [
        ('molecule', 'str', 'Molecule identifier'),
        ('electron_configuration', 'str', 'MO electron configuration'),
        ('bond_order', 'float', 'Bond order (number of chemical bonds)'),
        ('magnetic_property', 'str', 'Diamagnetic or paramagnetic'),
        ('stability', 'str', 'Stability assessment'),
        ('description', 'str', 'Detailed explanation of MO analysis'),
    ]
    
        
    examples = [
        {'code_input': {'molecule': 'N2'}, 'text_input': {'molecule': 'N2'}, 'output': {'molecule': 'N2', 'electron_configuration': '...', 'bond_order': 3, 'magnetic_property': 'diamagnetic', 'stability': 'very stable', 'description': '...'}},
        {'code_input': {'molecule': 'O2'}, 'text_input': {'molecule': 'O2'}, 'output': {'molecule': 'O2', 'electron_configuration': '...', 'bond_order': 2, 'magnetic_property': 'paramagnetic', 'stability': 'stable', 'description': '...'}},
        {'code_input': {'molecule': 'CO'}, 'text_input': {'molecule': 'CO'}, 'output': {'molecule': 'CO', 'electron_configuration': '...', 'bond_order': 3, 'magnetic_property': 'diamagnetic', 'stability': 'stable', 'description': '...'}},
    ]
    def _run_base(self, molecule: str) -> dict:
        mol = molecule.strip().upper()
        
        # Normalize input
        normalize_map = {
            "N2": "N2", "N≡N": "N2", "N=N": "N2",
            "O2": "O2", "O=O": "O2",
            "H2": "H2", "H-H": "H2",
            "F2": "F2", "F-F": "F2",
            "B2": "B2", "C2": "C2", "LI2": "Li2", "HE2": "He2", "NE2": "Ne2",
            "CO": "CO", "C≡O": "CO", "C=O": "CO",
            "NO": "NO", "HF": "HF", "H-F": "HF", "CN": "CN",
        }
        
        key = normalize_map.get(mol, mol)
        
        if key not in MO_DATA:
            available = sorted(MO_DATA.keys())
            raise ChemMCPInputError(
                f"MO data not available for '{molecule}'. "
                f"Available molecules: {available}. "
                f"This tool covers period 2 homonuclear diatomics (H2, Li2, B2, C2, N2, O2, F2, Ne2) "
                f"and heteronuclear diatomics (CO, NO, HF, CN)."
            )

        d = MO_DATA[key]
        return {
            "molecule": d["formula"],
            "total_valence_electrons": d["electrons"],
            "electron_configuration": d["configuration"],
            "bond_order": d["bond_order"],
            "magnetic_property": d["magnetic"],
            "stability": d["stability"],
            "description": d["description"],
            "isoelectronic_with": d.get("isoelectronic_with", None),
            "dipole_direction": d.get("dipole_direction", None),
        }


if __name__ == "__main__":
    run_mcp_server()
