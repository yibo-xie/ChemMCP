import logging
import math
from typing import Optional
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GibbsFreeEnergy(BaseTool):
    """
    计算给定温度下的吉布斯自由能变。
    公式: ΔG = ΔH - TΔS
    """
    __version__ = "0.1.0"
    name = "GibbsFreeEnergy"
    func_name = "calculate_gibbs_free_energy"
    description = "Calculate Gibbs free energy change at a given temperature using ΔG = ΔH - TΔS."
    implementation_description = "Uses the fundamental thermodynamic relation ΔG = ΔH - TΔS. Accepts enthalpy change (kJ/mol), entropy change (J/(mol·K)), and temperature (K)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Gibbs Energy", "Physical Chemistry", "Spontaneity"]
    required_envs = []

    code_input_sig = [
        ("delta_h", "float", "N/A", "Enthalpy change ΔH in kJ/mol."),
        ("delta_s", "float", "N/A", "Entropy change ΔS in J/(mol·K)."),
        ("temperature_k", "float", "N/A", "Absolute temperature in Kelvin."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'delta_h delta_s temperature_k'. ΔH in kJ/mol, ΔS in J/(mol·K), T in K."),
    ]

    output_sig = [
        ("delta_g", "float", "Gibbs free energy change ΔG in kJ/mol."),
        ("spontaneity", "str", "'spontaneous' if ΔG<0, 'non-spontaneous' if ΔG>0, 'at_equilibrium' if ΔG≈0."),
    ]

    examples = [
        {
            "code_input": {
                "delta_h": -100.0,
                "delta_s": -200.0,
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_params": "-100.0 -200.0 298.15",
            },
            "output": {
                "delta_g": -40.37,
                "spontaneity": "spontaneous",
            },
        },
        {
            "code_input": {
                "delta_h": 50.0,
                "delta_s": 100.0,
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_params": "50.0 100.0 298.15",
            },
            "output": {
                "delta_g": 20.185,
                "spontaneity": "non-spontaneous",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, delta_h: float, delta_s: float, temperature_k: float) -> dict:
        """Core logic: ΔG = ΔH - TΔS, with ΔH in kJ/mol and ΔS in J/(mol·K)."""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")

        # Convert ΔS from J/(mol·K) to kJ/(mol·K) for consistent units
        delta_s_kj = delta_s / 1000.0
        delta_g = delta_h - temperature_k * delta_s_kj

        # Determine spontaneity
        eps = 1e-6
        if delta_g < -eps:
            spontaneity = "spontaneous"
        elif delta_g > eps:
            spontaneity = "non-spontaneous"
        else:
            spontaneity = "at_equilibrium"

        logger.info(f"ΔG calculation: ΔH={delta_h} kJ/mol, ΔS={delta_s} J/(mol·K), T={temperature_k} K → ΔG={round(delta_g, 4)} kJ/mol ({spontaneity})")

        return {
            "delta_g": round(delta_g, 4),
            "spontaneity": spontaneity,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need at least 3 parameters: delta_h delta_s temperature_k")
            delta_h = float(parts[0])
            delta_s = float(parts[1])
            temperature_k = float(parts[2])
            return self._run_base(delta_h, delta_s, temperature_k)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'delta_h delta_s temperature_k'")
