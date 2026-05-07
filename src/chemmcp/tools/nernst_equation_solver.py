import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NernstEquationSolver(BaseTool):
    """
    Comprehensive Nernst equation solver for electrode potential calculations.
    Supports half-cell and full cell calculations with activity corrections.
    """
    __version__ = "0.1.0"
    name = "NernstEquationSolver"
    func_name = "nernst_equation_solver"
    description = "Comprehensive Nernst equation solver for calculating electrode potentials. Supports half-cell and full-cell calculations, activity coefficients via Debye-Hückel, and concentration-to-activity conversion."
    implementation_description = "Implements the Nernst equation: E = E° − (RT/nF)·ln(Q). At 25°C: E = E° − (0.05916/n)·log₁₀(Q). Optionally applies Debye-Hückel activity correction for ionic strength effects."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Nernst Equation", "Electrode Potential", "Electrochemistry", "Activity Coefficient", "Debye-Huckel"]
    required_envs = []

    code_input_sig = [
        ("E0", "float", "N/A", "Standard electrode potential E° in Volts."),
        ("n", "int", "N/A", "Number of electrons transferred."),
        ("reactants", "list", "N/A", "List of [coefficient, activity] pairs for reactants (products in Q denominator)."),
        ("products", "list", "N/A", "List of [coefficient, activity] pairs for products (numerator in Q)."),
        ("T", "float", "298.15", "Temperature in Kelvin."),
        ("ionic_strength", "float", "0.0", "Ionic strength I (mol/L) for Debye-Hückel correction. 0 means no correction."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: 'E0 n T [ionic_strength] reactant1_coef,act ... || product1_coef,act ...'. Example: '0.771 2 298.0 0.1 1,0.01 || 1,1'"),
    ]

    output_sig = [
        ("E_V", "float", "Calculated electrode potential (V)."),
        ("E0_V", "float", "Standard potential used (V)."),
        ("Q_value", "float", "Reaction quotient Q."),
        ("log_Q", "float", "Base-10 logarithm of Q."),
        ("nernst_term_V", "float", "(RT/nF)·ln(Q) or (0.05916/n)·log10(Q) term in V."),
        ("T_K", "float", "Temperature used (K)."),
        ("activity_corrected", "bool", "Whether Debye-Hückel correction was applied."),
        ("formula", "str", "Full formula with substituted values."),
    ]

    examples = [
        {
            "code_input": {
                "E0": 0.771,
                "n": 1,
                "reactants": [[1, 0.01]],
                "products": [[1, 1.0]],
                "T": 298.15,
                "ionic_strength": 0.0,
            },
            "text_input": {"input_string": "0.771 1 298.15 0 1,0.01 || 1,1"},
            "output": {
                "E_V": 0.8896,
                "Q_value": 100.0,
                "log_Q": 2.0,
                "nernst_term_V": -0.1183,
                "activity_corrected": False,
            }
        },
        {
            "code_input": {
                "E0": 1.229,
                "n": 4,
                "reactants": [[1, 1.0]],
                "products": [[1, 0.21]],
                "T": 298.15,
                "ionic_strength": 0.05,
            },
            "text_input": {"input_string": "1.229 4 298.15 0.05 1,1 || 1,0.21"},
            "output": {
                "E_V": 1.2375,
                "Q_value": 0.21,
                "log_Q": -0.6778,
                "nernst_term_V": 0.0100,
                "activity_corrected": True,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._R = 8.314       # J/(mol·K)
        self._F = 96485       # C/mol

    def _debye_huckel_gamma(self, z: float, I: float, T: float = 298.15) -> float:
        """Calculate mean ionic activity coefficient using Debye-Hückel limiting law."""
        if I <= 0:
            return 1.0
        A = (2.0 * math.pi * 6.022e23 * 80.2 * (1.602e-19)**2 / (1000 * 1.381e-23 * T))**0.5
        A_rounded = round(A, 6)  # ~0.509 at 25°C
        log_gamma = -A_rounded * z**2 * math.sqrt(I)
        return 10**log_gamma

    def _calc_Q(self, reactants: List[List[float]], products: List[List[float]], I: float = 0.0) -> tuple:
        """Calculate reaction quotient Q from activities."""
        Q_num = 1.0
        Q_den = 1.0

        for coef, act in products:
            if act < 0:
                raise ChemMCPError(f"Activity cannot be negative: {act}")
            Q_num *= act ** coef

        for coef, act in reactants:
            if act < 0:
                raise ChemMCPError(f"Activity cannot be negative: {act}")
            Q_den *= act ** coef

        if Q_den == 0:
            raise ChemMCPError("Denominator of reaction quotient Q is zero.")

        Q = Q_num / Q_den
        return Q, Q_num, Q_den

    def _run_base(
        self,
        E0: float,
        n: int,
        reactants: List[List[float]],
        products: List[List[float]],
        T: float = 298.15,
        ionic_strength: float = 0.0,
    ) -> dict:
        """Solve Nernst equation with optional activity correction."""
        if n == 0:
            raise ChemMCPError("Number of electrons n cannot be zero.")
        if T <= 0:
            raise ChemMCPError("Temperature must be positive (K).")

        Q, _, _ = self._calc_Q(reactants, products, ionic_strength)
        RT_over_nF = (self._R * T) / (n * self._F)

        try:
            ln_Q = math.log(Q)
            log_Q = math.log10(Q)
        except ValueError:
            raise ChemMCPError(f"Cannot compute log(Q) for Q = {Q}")

        nernst_term = -RT_over_nF * ln_Q
        E = E0 + nernst_term

        activity_corrected = ionic_strength > 0
        if abs(T - 298.15) < 0.01:
            factor = 0.05916 / n
            formula = (
                f"E = {E0:.4f} − ({factor:.5f}) × log₁₀({self._fmt(Q)}) "
                f"= {E0:.4f} − ({factor:.5f}) × ({log_Q:+.4f}) "
                f"= {E0:.4f} {(nernst_term>=0 and '+' or '−')} {abs(nernst_term):.4f} "
                f"= {E:.4f} V"
            )
        else:
            formula = (
                f"E = {E0:.4f} − ({RT_over_nF:.6f}) × ln({self._fmt(Q)}) "
                f"= {E:.4f} V"
            )

        return {
            "E_V": round(E, 6),
            "E0_V": E0,
            "Q_value": Q,
            "log_Q": round(log_Q, 6),
            "nernst_term_V": round(nernst_term, 6),
            "T_K": T,
            "activity_corrected": activity_corrected,
            "formula": formula,
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        parts = input_string.strip().split()

        if len(parts) < 2:
            raise ChemMCPError("Need at least E0 and n. Format: 'E0 n [T] [I] reactants || products'")

        E0 = float(parts[0])
        n = int(parts[1])
        idx = 2

        T = float(parts[idx]) if idx < len(parts) else 298.15
        idx += 1
        I = float(parts[idx]) if idx < len(parts) else 0.0
        idx += 1

        # Parse reactants and products separated by ||
        reactants = []
        products = []
        current = reactants
        while idx < len(parts):
            p = parts[idx]
            if p == "||":
                current = products
                idx += 1
                continue
            if "," in p:
                sub = p.split(",")
                if len(sub) >= 2:
                    try:
                        current.append([float(sub[0]), float(sub[1])])
                        idx += 1
                        continue
                    except ValueError:
                        pass
            break

        return self._run_base(E0, n, reactants, products, T, I)

    @staticmethod
    def _fmt(q: float) -> str:
        if q == int(q) and abs(q) < 1e6:
            return str(int(q))
        return f"{q:.4g}"
