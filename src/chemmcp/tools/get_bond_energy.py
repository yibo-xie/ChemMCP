import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.bonding_data import get_bond_energy, BOND_ENERGIES

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetBondEnergy(BaseTool):
    __version__ = "0.1.0"
    name = "GetBondEnergy"
    func_name = 'get_bond_energy'
    description = "Query bond dissociation energies (BDE) in kJ/mol for common chemical bonds."
    implementation_description = "Uses a database of average bond dissociation energies at 298 K from standard references (CRC Handbook, NIST). Returns energy required to homolytically cleave a bond in the gas phase."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Bond Energy", "Bond Dissociation Energy", "Thermochemistry", "Chemical Bonding"]
    required_envs = []

    code_input_sig = [
        ('bond', 'str', 'N/A', 'Bond specification (e.g., C-C, C=C, C≡C, C-H, O-H, N≡N)'),
    ]
    text_input_sig = [
        ('bond', 'str', 'N/A', 'Bond specification'),
    ]
    output_sig = [
        ('bond', 'str', 'Bond specification'),
        ('energy_kj_mol', 'float', 'Bond dissociation energy in kJ/mol'),
        ('energy_kcal_mol', 'float', 'Bond dissociation energy in kcal/mol'),
        ('note', 'str', 'Interpretation note'),
    ]
    
    examples = [
        {'code_input': {'bond': 'C-C'}, 'text_input': {'bond': 'C-C'}, 'output': {'energy_kj_mol': 347, 'energy_kcal_mol': 83, 'note': 'Average C-C single bond', 'bond': 'C-C'}},
        {'code_input': {'bond': 'C=O'}, 'text_input': {'bond': 'C=O'}, 'output': {'energy_kj_mol': 799, 'energy_kcal_mol': 191, 'note': 'Carbonyl (formaldehyde/ketone)', 'bond': 'C=O'}},
        {'code_input': {'bond': 'N≡N'}, 'text_input': {'bond': 'N≡N'}, 'output': {'energy_kj_mol': 945, 'energy_kcal_mol': 226, 'note': 'N≡N triple bond (very strong)', 'bond': 'N≡N'}},
    ]
    def _run_base(self, bond: str) -> dict:
        b = bond.strip()
        
        # Normalize bond notation
        # Handle various input formats: "C-C", "CC", "C=C", "C=C", "C#C", "C≡C", etc.
        lookup_key = b
        
        if b not in BOND_ENERGIES:
            # Try to find by normalizing
            replacements = {"−": "-", "–": "-", "_": "-", " ": "", "#": "≡"}
            normalized = b
            for old, new in replacements.items():
                normalized = normalized.replace(old, new)
            
            if normalized in BOND_ENERGIES:
                lookup_key = normalized
            else:
                available = sorted(set(BOND_ENERGIES.keys()))
                raise ChemMCPInputError(
                    f"Bond energy data not found for '{b}'. "
                    f"Available bonds: {available}. "
                    f"Format: element-element for single bonds, element=element for double, element≡element for triple."
                )

        energy = BOND_ENERGIES[lookup_key]
        
        # Interpretation notes
        strength_notes = {
            "H-H": "Very strong single bond; H2 is stable.",
            "F-F": "Weak bond due to lone pair repulsion; F2 is highly reactive.",
            "N≡N": "Triple bond makes N2 extremely inert (strongest diatomic bond).",
            "C=O": "Strong double bond; carbonyl group is very stable.",
            "C≡N": "Very strong triple bond; nitriles are thermally stable.",
            "O=O": "Relatively weak bond; O2 is reactive and paramagnetic.",
            "O-O": "Very weak bond; peroxides are unstable and reactive.",
            "Cl-Cl": "Moderate bond strength; Cl2 is moderately reactive.",
            "Si-Si": "Weak bond compared to C-C; Si-Si bonds are easily broken.",
        }
        
        note = strength_notes.get(lookup_key, "")
        if not note:
            if energy > 800:
                note = "Very strong bond."
            elif energy > 500:
                note = "Strong bond."
            elif energy > 300:
                note = "Moderate bond strength."
            elif energy > 150:
                note = "Relatively weak bond."
            else:
                note = "Weak bond."

        return {
            "bond": lookup_key,
            "energy_kj_mol": round(energy, 1),
            "energy_kcal_mol": round(energy / 4.184, 2),
            "unit": "kJ/mol at 298 K (gas phase, average value)",
            "note": note,
            "available_bonds_count": len(BOND_ENERGIES),
        }


if __name__ == "__main__":
    run_mcp_server()
