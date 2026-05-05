import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CellEmfCalculator(BaseTool):
    """
    原电池电动势计算工具。
    E°cell = E°cathode (reduction) - E°anode (reduction)
    ΔG° = -nFE°cell
    支持从半反应标准电极电势计算电池电动势和热力学量。
    """
    __version__      = "0.1.0"
    name             = "CellEmfCalculator"
    func_name        = "cell_emf_calculator"
    description      = "Calculate cell EMF (electromotive force) from standard reduction potentials of cathode and anode half-cells."
    implementation_description = "Computes E°cell = E°cathode - E°anode (both as reduction potentials), then derives ΔG° = -nFE°cell, equilibrium constant K, and spontaneity of the cell reaction."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Cell EMF", "Electrochemistry", "Galvanic Cell", "Standard Potential", "Thermodynamics"]
    required_envs    = []

    code_input_sig   = [
        ("E0_cathode_v", "float", "N/A", "Standard reduction potential of the cathode in Volts."),
        ("E0_anode_v", "float", "N/A", "Standard reduction potential of the anode in Volts."),
        ("n_electrons", "int", "2", "Number of electrons transferred in the balanced reaction. Default: 2."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin. Default: 298.15 K (25°C)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated string: 'E0_cathode E0_anode [n_electrons] [temperature_K]', e.g., '0.34 -0.76' or '1.507 0.771 5'."),
    ]

    output_sig       = [
        ("E0_cathode_V", "float", "Cathode standard reduction potential used (V)."),
        ("E0_anode_V", "float", "Anode standard reduction potential used (V)."),
        ("E0_cell_V", "float", "Standard cell EMF E°cell (V)."),
        ("n_electrons", "int", "Number of electrons transferred."),
        ("delta_G0_kJ_mol", "float", "Standard Gibbs free energy change ΔG° in kJ/mol."),
        ("delta_G0_J_mol", "float", "Standard Gibbs free energy change ΔG° in J/mol."),
        ("equilibrium_constant_K", "float", "Equilibrium constant K (dimensionless)."),
        ("spontaneous", "bool", "Whether the cell reaction is spontaneous under standard conditions."),
        ("reaction_type", "str", "'galvanic' if E°cell > 0, 'electrolytic' if E°cell < 0."),
        ("summary", "str", "Human-readable summary of the cell calculation."),
    ]

    examples         = [
        {
            "code_input": {
                "E0_cathode_v": 0.34,
                "E0_anode_v": -0.76,
                "n_electrons": 2,
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_params": "0.34 -0.76 2"
            },
            "output": {
                "E0_cathode_V": 0.34,
                "E0_anode_V": -0.76,
                "E0_cell_V": 1.10,
                "n_electrons": 2,
                "delta_G0_kJ_mol": -212.289,
                "delta_G0_J_mol": -212289.237,
                "equilibrium_constant_K": 3.882e+37,
                "spontaneous": True,
                "reaction_type": "galvanic",
                "summary": "Daniell cell (Zn|Cu): E°cell = 1.10 V, spontaneous galvanic cell.",
            }
        },
        {
            "code_input": {
                "E0_cathode_v": -0.14,
                "E0_anode_v": 0.771,
                "n_electrons": 3,
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_params": "-0.14 0.771 3"
            },
            "output": {
                "E0_cathode_V": -0.14,
                "E0_anode_V": 0.771,
                "E0_cell_V": -0.911,
                "n_electrons": 3,
                "delta_G0_kJ_mol": 263.699,
                "delta_G0_J_mol": 263698.868,
                "equilibrium_constant_K": 4.856e-47,
                "spontaneous": False,
                "reaction_type": "electrolytic",
                "summary": "Non-spontaneous: requires external voltage > 0.911 V to drive electrolysis.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618     # J/(mol·K)
        self.F = 96485.33212     # C/mol, Faraday constant

    def _run_base(self, E0_cathode_v: float, E0_anode_v: float,
                  n_electrons: int = 2, temperature_k: float = 298.15) -> dict:
        """Calculate cell EMF and thermodynamic quantities."""
        if n_electrons <= 0:
            raise ChemMCPError("Number of electrons must be positive.")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")

        R = self.R
        F = self.F
        T = temperature_k
        n = n_electrons

        # Cell EMF: E°cell = E°cathode - E°anode (both reduction potentials)
        E0_cell = E0_cathode_v - E0_anode_v

        # Gibbs free energy: ΔG° = -nFE°cell
        delta_G0_J = -n * F * E0_cell
        delta_G0_kJ = delta_G0_J / 1000.0

        # Equilibrium constant: ΔG° = -RT ln K => K = exp(-ΔG°/RT)
        try:
            import math
            if abs(T) < 1e-10 or abs(R) < 1e-10:
                K = float('inf') if delta_G0_J < 0 else 0.0
            else:
                exponent = -delta_G0_J / (R * T)
                if exponent > 700:
                    K = float('inf')
                elif exponent < -700:
                    K = 0.0
                else:
                    K = math.exp(exponent)
        except (OverflowError, ValueError):
            K = float('inf') if delta_G0_J < 0 else 0.0

        spontaneous = E0_cell > 0
        rtype = "galvanic" if spontaneous else "electrolytic"

        summary = (
            f"Cell: E°cathode = {E0_cathode_v:.3f} V, E°anode = {E0_anode_v:.3f} V\n"
            f"E°cell = {E0_cathode_v:.3f} - ({E0_anode_v:.3f}) = {E0_cell:.3f} V\n"
            f"ΔG° = -nFE°cell = -({n})({F:.0f})({E0_cell:.3f}) = {delta_G0_kJ:.3f} kJ/mol\n"
            f"K = exp(-ΔG°/RT) = {self._fmt_K(K)}\n"
            f"Reaction is {'spontaneous (galvanic)' if spontaneous else 'non-spontaneous (electrolytic)'} under standard conditions."
        )

        return {
            "E0_cathode_V": round(E0_cathode_v, 6),
            "E0_anode_V": round(E0_anode_v, 6),
            "E0_cell_V": round(E0_cell, 6),
            "n_electrons": n,
            "delta_G0_kJ_mol": round(delta_G0_kJ, 3),
            "delta_G0_J_mol": round(delta_G0_J, 3),
            "equilibrium_constant_K": K,
            "spontaneous": spontaneous,
            "reaction_type": rtype,
            "summary": summary,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse space-separated text input."""
        parts = input_params.strip().split()
        if len(parts) < 2:
            raise ChemMCPError(
                "Text input requires at least E0_cathode and E0_anode. "
                "Format: 'E0_cathode E0_anode [n] [T]'"
            )

        try:
            E0_cat = float(parts[0])
            E0_ano = float(parts[1])
            n = int(parts[2]) if len(parts) > 2 else 2
            T = float(parts[3]) if len(parts) > 3 else 298.15
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse numeric values from '{input_params}': {e}")

        return self._run_base(E0_cat, E0_ano, n, T)

    @staticmethod
    def _fmt_K(K) -> str:
        """Format equilibrium constant for display."""
        if K == 0:
            return "~0"
        elif K == float('inf'):
            return "very large (∞)"
        elif K >= 1e10 or K <= 1e-10:
            return f"{K:.4e}"
        else:
            return f"{K:.6f}"
