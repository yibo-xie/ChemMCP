import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.bonding_data import get_bond_length, estimate_bond_length, BOND_LENGTHS

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetBondLength(BaseTool):
    __version__ = "0.1.0"
    name = "GetBondLength"
    func_name = 'get_bond_length'
    description = "Query standard bond lengths in picometers (pm) or Angstroms."
    implementation_description = "Uses a database of average experimental bond lengths from CRC Handbook and crystallographic data. Returns bond length for a given element pair and bond type (single, double, triple, aromatic). Can also estimate from covalent radii sum if exact data is unavailable."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Bond Length", "Chemical Bonding", "Structural Chemistry"]
    required_envs = []

    code_input_sig = [
        ('element1', 'str', 'N/A', 'First element symbol (e.g., C)'),
        ('element2', 'str', 'N/A', 'Second element symbol (e.g., O)'),
        ('bond_type', 'str', 'single', 'Bond type: single, double, triple, or aromatic'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Bond specification, e.g., \"C-O single\" or \"C≡N\"'),
    ]
    output_sig = [
        ('bond', 'str', 'Bond specification'),
        ('length_pm', 'float', 'Bond length in picometers'),
        ('length_angstrom', 'float', 'Bond length in Angstroms'),
        ('source', 'str', 'Data source (experimental/covalent radii estimate)'),
    ]
    
    examples = [
        {'code_input': {'element1': 'C', 'element2': 'C', 'bond_type': 'single'}, 'text_input': {'query': 'C-C single'}, 'output': {'length_pm': 154, 'length_angstrom': 1.54, 'bond': 'C-C single', 'source': 'CRC Handbook'}},
        {'code_input': {'element1': 'C', 'element2': 'O', 'bond_type': 'double'}, 'text_input': {'query': 'C=O double'}, 'output': {'length_pm': 123, 'length_angstrom': 1.23, 'bond': 'C=O double', 'source': 'CRC Handbook'}},
        {'code_input': {'element1': 'C', 'element2': 'N', 'bond_type': 'triple'}, 'text_input': {'query': 'C≡N triple'}, 'output': {'length_pm': 116, 'length_angstrom': 1.16, 'bond': 'C≡N triple', 'source': 'CRC Handbook'}},
    ]
    def _run_base(self, element1: str, element2: str, bond_type: str = "single") -> dict:
        # Normalize bond type
        bt_map = {"single": "single", "double": "double", "triple": "triple", 
                  "aromatic": "aromatic", "1": "single", "2": "double", "3": "triple",
                  "=": "double", "#": "triple", "-": "single"}
        bt = bt_map.get(bond_type.lower(), bond_type.lower())
        
        e1 = element1.strip().capitalize()
        e2 = element2.strip().capitalize()
        bond_spec = f"{e1}-{e2} ({bt})"
        
        # Try exact lookup
        length = get_bond_length(e1, e2, bt)
        
        if length is not None:
            source = "experimental average (CRC Handbook / standard reference)"
        else:
            try:
                length = estimate_bond_length(e1, e2)
                source = "estimated from covalent radii sum"
            except ValueError:
                raise ChemMCPInputError(
                    f"No bond length data available for {e1}-{e2} ({bt}). "
                    f"Available elements in database: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, "
                    f"Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, "
                    f"Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, "
                    f"Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi."
                )

        return {
            "bond": bond_spec,
            "element1": e1,
            "element2": e2,
            "bond_type": bt,
            "length_pm": round(length, 1),
            "length_angstrom": round(length / 100, 3),
            "source": source,
        }


if __name__ == "__main__":
    run_mcp_server()
