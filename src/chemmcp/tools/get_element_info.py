import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetElementInfo(BaseTool):
    __version__ = "0.1.0"
    name = "GetElementInfo"
    func_name = 'get_element_info'
    description = "Get complete element information by element symbol or atomic number."
    implementation_description = "Uses a built-in periodic table database with IUPAC data for all 118 elements. Returns atomic number, symbol, name, atomic weight, Pauling electronegativity, electron configuration, category, group, period, and block."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Elements", "Properties", "Atomic Data"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., Fe) or atomic number (e.g., 26)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol or atomic number'),
    ]
    output_sig = [
        ('info', 'dict', 'Complete element information dictionary'),
    ]
    examples = [
        {'code_input': {'element': 'O'}, 'text_input': {'element': 'O'}, 'output': {'info': {'atomic_number': 8, 'symbol': 'O'}}},
        {'code_input': {'element': 26}, 'text_input': {'element': '26'}, 'output': {'info': {'atomic_number': 26, 'symbol': 'Fe'}}},
    ]

    def _run_base(self, element):
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}. Please provide a valid element symbol (e.g., O, Fe, U) or atomic number (1-118).")
        return {
            "atomic_number": data["atomic_number"],
            "symbol": data["symbol"],
            "name": data["name"],
            "atomic_weight": data["atomic_weight"],
            "electronegativity_pauling": data["en_pauling"],
            "electron_configuration": data["electron_config"],
            "category": data["category"],
            "group": data["group"],
            "period": data["period"],
            "block": data["block"],
        }


if __name__ == "__main__":
    run_mcp_server()
