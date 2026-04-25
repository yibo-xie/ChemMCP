import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)

# Electron affinity in kJ/mol (NIST CRC data)
ELECTRON_AFFINITY: dict = {
    "H": 72.8, "He": None, "Li": 59.6, "Be": None, "B": 26.7,
    "C": 121.9, "N": -6.8, "O": 141.0, "F": 328.0, "Ne": None,
    "Na": 52.9, "Mg": None, "Al": 42.5, "Si": 133.6, "P": 72.0,
    "S": 200.4, "Cl": 349.0, "Ar": None, "K": 48.4, "Ca": 2.37,
    "Sc": 18.1, "Ti": 7.6, "V": 50.6, "Cr": 64.3, "Mn": None,
    "Fe": 14.5, "Co": 63.9, "Ni": 112.0, "Cu": 118.4, "Zn": None,
    "Ga": 28.9, "Ge": 119.0, "As": 78.2, "Se": 195.0, "Br": 324.5,
    "Kr": None, "Rb": 46.9, "Sr": 5.03, "Y": 29.6, "Zr": 41.1,
    "Nb": 86.2, "Mo": 71.9, "Tc": 53.0, "Ru": 101.3, "Rh": 109.7,
    "Pd": 53.7, "Ag": 125.6, "Cd": None, "In": 28.9, "Sn": 107.3,
    "Sb": 101.2, "Te": 190.2, "I": 295.2, "Xe": None, "Cs": 45.5,
    "Ba": 13.95, "La": 48.5, "Ce": 50.0, "Pr": 50.0, "Nd": 50.0,
    "Pm": 50.0, "Sm": 50.0, "Eu": 50.0, "Gd": 50.0, "Tb": 50.0,
    "Dy": 50.0, "Ho": 50.0, "Er": 50.0, "Tm": 50.0, "Yb": 50.0,
    "Lu": 33.0, "Hf": None, "Ta": 31.0, "W": 78.6, "Re": 14.2,
    "Os": 106.3, "Ir": 151.0, "Pt": 205.3, "Au": 222.8, "Hg": None,
    "Tl": 19.2, "Pb": 35.2, "Bi": 91.2, "Po": 183.0, "At": 270.2,
}


@ChemMCPManager.register_tool
class GetElectronAffinity(BaseTool):
    __version__ = "0.1.0"
    name = "GetElectronAffinity"
    func_name = 'get_electron_affinity'
    description = "Query electron affinity of an element in kJ/mol."
    implementation_description = "Returns the first electron affinity value from NIST/CRC Handbook data. Negative values indicate energy release (exothermic). Returns None for elements that do not readily accept an electron."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Electron Affinity", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., Cl, O)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('element', 'str', 'Element symbol'),
        ('electron_affinity_kj_mol', 'float', 'Electron affinity in kJ/mol'),
        ('note', 'str', 'Interpretation note'),
    ]
    
    examples = [
        {'code_input': {'element': 'Cl'}, 'text_input': {'element': 'Cl'}, 'output': {'element': 'Cl', 'electron_affinity_kj_mol': 349.0, 'note': '...'}},
        {'code_input': {'element': 'N'}, 'text_input': {'element': 'N'}, 'output': {'element': 'N', 'electron_affinity_kj_mol': -6.8, 'note': '...'}},
    ]
    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        sym = data["symbol"]
        ea = ELECTRON_AFFINITY.get(sym)
        if ea is None:
            return {
                "element": sym,
                "electron_affinity_kj_mol": None,
                "note": f"Electron affinity data not available for {sym} (noble gases and some metals do not have stable negative ions).",
            }
        note = (
            f"Exothermic: {sym} releases {abs(ea):.1f} kJ/mol when gaining an electron."
            if ea > 0 else
            f"Endothermic: {sym} requires {-ea:.1f} kJ/mol to add an electron (unfavorable)."
            if ea < 0 else
            "Value is approximately zero."
        )
        return {
            "element": sym,
            "electron_affinity_kj_mol": round(ea, 1) if ea is not None else None,
            "unit": "kJ/mol",
            "note": note,
        }


if __name__ == "__main__":
    run_mcp_server()
