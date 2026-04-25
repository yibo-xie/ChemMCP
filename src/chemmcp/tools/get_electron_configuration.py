import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetElectronConfiguration(BaseTool):
    __version__ = "0.1.0"
    name = "GetElectronConfiguration"
    func_name = 'get_electron_configuration'
    description = "Get the electron configuration of an element, including full and noble-gas shorthand forms."
    implementation_description = "Uses the built-in periodic table database to return both the full electron configuration (e.g., '1s² 2s² 2p⁶ 3s² 3p⁴') and the noble-gas shorthand (e.g., '[Ne] 3s² 3p⁴')."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Electron Configuration", "Quantum Chemistry"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., S, Fe)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('full_config', 'str', 'Full electron configuration'),
        ('noble_config', 'str', 'Noble gas shorthand notation'),
        ('symbol', 'str', 'Element symbol'),
    ]
    examples = [
        {'code_input': {'element': 'S'}, 'text_input': {'element': 'S'}, 'output': {'full_config': '[Ne] 3s² 3p⁴', 'noble_config': '[Ne] 3s² 3p⁴', 'symbol': 'S'}},
        {'code_input': {'element': 'Fe'}, 'text_input': {'element': 'Fe'}, 'output': {'full_config': '[Ar] 3d⁶ 4s²', 'noble_config': '[Ar] 3d⁶ 4s²', 'symbol': 'Fe'}},
    ]

    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        return {
            "full_config": data["electron_config"],
            "noble_config": data["noble_config"],
            "symbol": data["symbol"],
        }


if __name__ == "__main__":
    run_mcp_server()
