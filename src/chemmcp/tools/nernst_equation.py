import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NernstEquation(BaseTool):
    """
    Calculate non-standard electrode potential using the Nernst equation.
    E = E° - (RT/nF) ln(Q)  or  E = E° - (0.05916/n) log10(Q) at 25°C.
    """
    __version__ = "0.1.0"
    name = "NernstEquation"
    func_name = "nernst_equation"
    description = "Calculate non-standard cell/electrode potential using the Nernst equation. Supports both code (numeric parameters) and text (string) interfaces."
    implementation_description = "Implements the Nernst equation: E = E° − (RT/nF)·ln(Q). At 25°C: E = E° − (0.05916 V/n)·log₁₀(Q). Computes reaction quotient Q from activities/concentrations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Nernst Equation", "Electrochemistry", "Non-standard Potential", "Physical Chemistry", "Thermodynamics"]
    required_envs = []

    code_input_sig = [
        ("E0", "float", "N/A", "Standard electrode potential E° in Volts."),
        ("n", "int", "N/A", "Number of electrons transferred in the half-reaction or overall reaction."),
        ("T", "float", "298.15", "Temperature in Kelvin (K). Default: 298.15 K (25°C)."),
        ("Q", "float", "1.0", "Reaction quotient Q (dimensionless). Default: 1.0 (standard conditions)."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Space-separated values: 'E0 n [T] [Q]', e.g., '0.771 2' or '0.34 2 310 0.01'."),
    ]

    output_sig = [
        ("E_V", "float", "Calculated non-standard electrode potential in Volts."),
        ("E0_V", "float", "Standard potential used (V)."),
        ("T_K", "float", "Temperature used (K)."),
        ("Q_value", "float", "Reaction quotient used."),
        ("n_electrons", "int", "Number of electrons transferred."),
        ("RT_over_nF_V", "float", "The Nernst factor (RT/nF) in Volts at given T."),
        ("formula_applied", "str", "The Nernst equation with substituted values."),
    ]

    examples = [
        {
            "code_input": {"E0": 0.337, "n": 2, "T": 298.15, "Q": 0.01},
            "text_input": {"input_string": "0.337 2 298.15 0.01"},
            "output": {
                "E_V": 0.396,
                "E0_V": 0.337,
                "T_K": 298.15,
                "Q_value": 0.01,
                "n_electrons": 2,
                "RT_over_nF_V": 0.01284,
                "formula_applied": "E = 0.337 − (0.05916/2) × log₁₀(0.01) = 0.337 − (−0.05916) = 0.396 V",
            }
        },
        {
            "code_input": {"E0": 1.507, "n": 5, "T": 298.15, "Q": 1000},
            "text_input": {"input_string": "1.507 5 1000"},
            "output": {
                "E_V": 1.472,
                "E0_V": 1.507,
                "T_K": 298.15,
                "Q_value": 1000,
                "n_electrons": 5,
                "RT_over_nF_V": 0.00486,
                "formula_applied": "E = 1.507 − (0.05916/5) × log₁₀(1000) = 1.507 − 0.036 = 1.471 V",
            }
        },
        {
            "code_input": {"E0": 0.00, "n": 2, "T": 298.15, "Q": 1e-7},
            "text_input": {"input_string": "0 2 1e-7"},
            "output": {
                "E_V": -0.207,
                "E0_V": 0.000,
                "T_K": 298.15,
                "Q_value": 1e-7,
                "n_electrons": 2,
                "RT_over_nF_V": 0.02958,
                "formula_applied": "E = 0.000 − (0.05916/2) × log₁₀(10⁻⁷) = 0 − (−0.207) = +0.207 V",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._R = 8.314       # J/(mol·K), universal gas constant
        self._F = 96485       # C/mol, Faraday constant

    def _run_base(self, E0: float, n: int, T: float = 298.15, Q: float = 1.0) -> dict:
        """Calculate non-standard potential using Nernst equation."""
        if n == 0:
            raise ChemMCPError("Number of electrons n cannot be zero.")
        if T <= 0:
            raise ChemMCPError("Temperature T must be positive (in Kelvin).")
        if Q <= 0:
            raise ChemMCPError("Reaction quotient Q must be positive.")

        # Nernst equation: E = E° - (RT / nF) * ln(Q)
        RT_over_nF = (self._R * T) / (n * self._F)

        try:
            ln_Q = math.log(Q)
        except ValueError:
            raise ChemMCPError(f"Cannot compute ln(Q) for Q = {Q}")

        E = E0 - RT_over_nF * ln_Q

        # Also compute base-10 form for readability at 25°C
        if abs(T - 298.15) < 0.01:
            factor_05916 = 0.05916 / n  # V at 25°C
            log_Q = math.log10(Q) if Q > 0 else 0
            formula = f"E = {E0:.3f} − ({factor_05916:.5f}) × log₁₀({self._fmt_q(Q)}) = {E0:.3f} − ({log_Q:+.3f}) × {factor_05916:.5f} = {E:.3f} V"
        else:
            formula = f"E = {E0:.3f} − ({RT_over_nF:.5f}) × ln({self._fmt_q(Q)}) = {E0:.3f} − ({ln_Q:+.4f}) × {RT_over_nF:.5f} = {E:.3f} V"

        return {
            "E_V": round(E, 6),
            "E0_V": E0,
            "T_K": T,
            "Q_value": Q,
            "n_electrons": n,
            "RT_over_nF_V": round(RT_over_nF, 6),
            "formula_applied": formula,
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse space-separated input string."""
        parts = input_string.strip().split()
        if len(parts) < 2:
            raise ChemMCPError(
                f"Text input requires at least E0 and n. Format: 'E0 n [T] [Q]'. "
                f"Example: '0.771 2' or '0.34 2 310 0.01'"
            )

        try:
            E0 = float(parts[0])
            n = int(parts[1])
            T = float(parts[2]) if len(parts) > 2 else 298.15
            Q = float(parts[3]) if len(parts) > 3 else 1.0
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse numeric values from '{input_string}': {e}")

        return self._run_base(E0, n, T, Q)

    @staticmethod
    def _fmt_q(q: float) -> str:
        """Format Q value for display."""
        if q == int(q) and abs(q) < 1e6:
            return str(int(q))
        if abs(q) >= 1000 or abs(q) < 0.01:
            return f"{q:.2e}"
        return str(q)
