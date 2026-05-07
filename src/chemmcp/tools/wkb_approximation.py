import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class WKBApproximation(BaseTool):
    """
    WKB (Wentzel-Kramers-Brillouin) 近似工具。
    计算半经典隧穿概率、量子隧穿系数、束缚态能级（Bohr-Sommerfeld量子化条件）。
    """
    __version__ = "0.1.0"
    name = "WKBApproximation"
    func_name = "wkb_calculate"
    description = "Calculate WKB semi-classical approximation: tunneling probability, penetration depth, and bound state energies via Bohr-Sommerfeld quantization."
    implementation_description = "Implements WKB tunneling integral ∫√(2m(V-E))dx/ℏ, transmission coefficient T≈exp(-2γ), and Bohr-Sommerfeld quantization condition ∮p dx = (n+½)h for potential wells."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "WKB", "Tunneling", "Semi-classical", "Barrier Penetration"]
    required_envs = []

    code_input_sig = [
        ("calculation_type", "str", "N/A", "Type: 'tunneling' (barrier penetration), 'transmission' (T coefficient), 'bound_state' (energy levels), 'connection' (connection formulas)."),
        ("potential_type", "str", "'square'", "Potential shape: 'square', 'triangular', 'parabolic', 'coulomb', 'delta', 'general'."),
        ("E_eV", "float", "1.0", "Particle energy in eV."),
        ("V0_eV", "float", "5.0", "Barrier height in eV (for barrier problems) or well depth (for bound states)."),
        ("width_nm", "float", "0.5", "Barrier/well width in nanometers."),
        ("mass_kg", "float", "9.109e-31", "Particle mass in kg (default: electron mass)."),
        ("n_state", "int", "0", "Quantum number n for bound state calculation."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: calculation_type potential_type E_eV V0_eV width_nm [mass_kg n_state]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing WKB result: tunneling probability, gamma integral, connection formulas, energy levels, etc."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_type": "tunneling",
                "potential_type": "square",
                "E_eV": 1.0,
                "V0_eV": 5.0,
                "width_nm": 0.5,
                "mass_kg": 9.109e-31,
                "n_state": 0,
            },
            "text_input": {
                "input_str": "tunneling square 1.0 5.0 0.5",
            },
            "output": {
                "result": {
                    "calculation_type": "tunneling",
                    "potential": "square_barrier",
                    "gamma_wkb": "...",
                    "tunneling_probability": "...",
                    "formula": "T ≈ exp(-2γ), γ = ∫√(2m(V-E))/ℏ dx",
                }
            }
        },
        {
            "code_input": {
                "calculation_type": "bound_state",
                "potential_type": "parabolic",
                "E_eV": 0,
                "V0_eV": 10.0,
                "width_nm": 1.0,
                "n_state": 2,
            },
            "text_input": {
                "input_str": "bound_state parabolic 0 10 1.0 9.109e-31 2",
            },
            "output": {
                "result": {
                    "calculation_type": "bound_state",
                    "energy_level_n": 2,
                    "energy_eV": "...",
                    "bohr_sommerfeld_condition": "∮p dx = (n+½)h",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34   # J·s
        self.eV = 1.602176634e-19      # J per eV
        self.m_e = 9.1093837015e-31    # electron mass, kg

    def _run_base(self, calculation_type: str, potential_type: str = "square",
                  E_eV: float = 1.0, V0_eV: float = 5.0, width_nm: float = 0.5,
                  mass_kg: float = 9.109e-31, n_state: int = 0) -> dict:
        """Core logic: WKB calculation."""
        calc = calculation_type.lower().strip()
        pot = potential_type.lower().strip()
        E_J = E_eV * self.eV
        V0_J = V0_eV * self.eV
        width_m = width_nm * 1e-9

        if calc == "tunneling" or calc == "transmission":
            return self._tunneling(pot, E_J, V0_J, width_m, mass_kg)
        elif calc == "bound_state":
            return self._bound_state(pot, V0_J, width_m, mass_kg, n_state)
        elif calc == "connection":
            return self._connection_formulas(pot, E_J, V0_J, width_m, mass_kg)
        else:
            raise ChemMCPError(
                f"Unknown calculation type '{calculation_type}'. "
                f"Use: 'tunneling', 'transmission', 'bound_state', or 'connection'."
            )

    # ── Tunneling / Transmission ────────────────────────────────────
    def _tunneling(self, pot: str, E_J: float, V0_J: float, width: float, mass: float) -> dict:
        """Compute WKB tunneling probability through a barrier."""
        if E_J >= V0_J:
            # Above barrier — classical transmission with some reflection
            T_classical = 1.0
            return self._make_result("tunneling", pot, T_classical, None,
                                     f"E ({E_J/self.eV:.2f} eV) ≥ V₀ ({V0_J/self.eV:.2f} eV): particle is above barrier. Classical transmission T ≈ 1.",
                                     extra={"above_barrier": True})

        if pot == "square":
            gamma = self._gamma_square(E_J, V0_J, width, mass)
        elif pot == "triangular":
            gamma = self._gamma_triangular(E_J, V0_J, width, mass)
        elif pot == "parabolic":
            gamma = self._gamma_parabolic(E_J, V0_J, width, mass)
        elif pot == "coulomb":
            gamma = self._gamma_coulomb(E_J, V0_J, width, mass)
        elif pot == "delta":
            gamma = self._gamma_delta(V0_J, width, mass)
        else:
            raise ChemMCPError(f"Unknown potential type '{pot}'. Use: square, triangular, parabolic, coulomb, delta.")

        T = math.exp(-2.0 * gamma)

        return self._make_result("tunneling", pot, T, gamma,
                                 f"WKB transmission through {pot} barrier: T ≈ exp(-2γ)",
                                 extra={
                                     "gamma": round(gamma, 6),
                                     "E_eV": round(E_J / self.eV, 6),
                                     "V0_eV": round(V0_J / self.eV, 6),
                                     "width_nm": width * 1e9,
                                 })

    def _gamma_square(self, E: float, V0: float, a: float, m: float) -> float:
        """γ = ∫₀ᵃ √(2m(V(x)-E))/ℏ dx  for square barrier V=V0 constant."""
        kappa = math.sqrt(2.0 * m * (V0 - E)) / self.hbar
        return kappa * a

    def _gamma_triangular(self, E: float, V0: float, a: float, m: float) -> float:
        """Triangular barrier: V(x) = V0(1 - x/a) for 0<x<a."""
        # γ = (2/3) · √(2m)/ℏ · (V0-E)^{3/2} / |dV/dx|
        dVdx = V0 / a  # slope magnitude
        gamma = (2.0 / 3.0) * math.sqrt(2.0 * m) * ((V0 - E) ** 1.5) / (self.hbar * dVdx)
        return gamma

    def _gamma_parabolic(self, E: float, V0: float, a: float, m: float) -> float:
        """Parabolic barrier near top: V(x) = V0 - ½mω²x²."""
        # Approximate using curvature at barrier top
        omega = math.sqrt(2.0 * V0 / (m * a * a))
        delta_E = V0 - E
        gamma = math.pi * delta_E / (self.hbar * omega)
        return gamma

    def _gamma_coulomb(self, E: float, V0: float, a: float, m: float) -> float:
        """Coulomb-like barrier (alpha decay style)."""
        # Simplified Gamow factor approximation
        Z_eff = 50  # effective atomic number (approximate)
        e_charge = 1.602176634e-19
        # Use simplified formula
        eta = Z_eff * e_charge * e_charge / (4 * math.pi * 8.8541878128e-12 * self.hbar) * math.sqrt(m / (2 * (V0 - E)))
        gamma = eta * math.acos(math.sqrt(E / V0)) - math.sqrt((V0 / E) - 1)
        return abs(gamma)

    def _gamma_delta(self, V0: float, a: float, m: float) -> float:
        """Delta-function barrier: V(x) = V0·a·δ(x)."""
        # Transmission for delta barrier: T = 1/(1 + (mV0a/ℏ)²)
        g = m * V0 * a / (self.hbar ** 2)
        return 0.5 * math.log(1 + g * g)  # such that exp(-2γ) matches

    # ── Bound States (Bohr-Sommerfeld) ──────────────────────────────
    def _bound_state(self, pot: str, V0_J: float, width: float, mass: float, n: int) -> dict:
        """Compute bound state energies via Bohr-Sommerfeld quantization."""
        if n < 0:
            raise ChemMCPError("Quantum number n must be >= 0.")

        if pot == "parabolic":
            # Harmonic oscillator — exact result from B-S: E_n = (n+½)ℏω
            omega = math.sqrt(2.0 * V0_J / (mass * width * width))
            E_n = (n + 0.5) * self.hbar * omega
            E_n_eV = E_n / self.eV
            return self._make_result("bound_state", pot, E_n_eV, None,
                                     f"Parabolic well (harmonic oscillator) E_{n} = (n+½)ℏω = {E_n_eV:.6e} eV",
                                     extra={
                                         "quantum_number": n,
                                         "omega_rad_s": round(omega, 6),
                                         "bohr_sommerfeld": "∮p dx = (n+½)h → E_n = (n+½)ℏω",
                                     })
        elif pot == "square":
            # Infinite square well (approximate — B-S gives exact answer here)
            E_n = ((n + 1) ** 2 * math.pi ** 2 * self.hbar ** 2) / (2.0 * mass * width ** 2)
            E_n_eV = E_n / self.eV
            return self._make_result("bound_state", pot, E_n_eV, None,
                                     f"Square well E_{n+1} = (n+1)²π²ℏ²/(2mL²) = {E_n_eV:.6e} eV",
                                     extra={
                                         "quantum_number": n,
                                         "bohr_sommerfeld": "∮p dx = (n+1)h → exact for infinite well",
                                     })
        elif pot == "coulomb":
            # Hydrogen-like atom: E_n = -Z²Ry/n²
            # Use V0 as reference and width as approximate Bohr radius scale
            Ry = mass * (1.602176634e-19) ** 4 * (4 * math.pi * 8.8541878128e-12) ** 2 / (2 * self.hbar ** 2)  # rough
            Ry_eV = 13.6  # known value
            n_eff = n + 1
            E_n_eV = -Ry_eV / (n_eff * n_eff)
            return self._make_result("bound_state", pot, E_n_eV, None,
                                     f"Coulomb (hydrogen-like) E_{n_eff} = -Ry/n² = {E_n_eV:.6f} eV",
                                     extra={
                                         "quantum_number": n,
                                         "n_effective": n_eff,
                                         "rydberg_eV": Ry_eV,
                                     })
        else:
            raise ChemMCPError(f"Bound state for potential type '{pot}' not implemented. Use: parabolic, square, coulomb.")

    # ── Connection Formulas ─────────────────────────────────────────
    def _connection_formulas(self, pot: str, E_J: float, V0_J: float, width: float, mass: float) -> dict:
        """Return WKB connection formulas linking oscillatory and exponential regions."""
        formulas = {
            "turning_point_right": (
                "ψ(x) ≈ (2/κ)^(1/2) · exp(-∫ₓˣ² κ dx)    for x > x₂ (classically forbidden)"
            ),
            "turning_point_left": (
                "ψ(x) ≈ (2/κ)^(1/2) · exp(+∫ˣ¹ˣ κ dx)    for x < x₁ (classically forbidden)"
            ),
            "classical_region": (
                "ψ(x) ≈ (2/p)^(1/2) · sin(∫ˣˣ₁ p dx/ℏ + π/4)    for x₁ < x < x₂ (allowed)"
            ),
            "matching_condition": (
                "Connection rule: ψ oscillatory ↔ ψ exponential at classical turning points where E = V(x)"
            ),
            "kappa_definition": (
                "κ(x) = √(2m[V(x)-E])/ℏ    (decay constant in forbidden region)"
            ),
            "momentum_definition": (
                "p(x) = √(2m[E-V(x)])         (momentum in allowed region)"
            ),
        }
        return self._make_result("connection", pot, None, None,
                                 "WKB connection formulas at classical turning points.",
                                 extra=formulas)

    def _make_result(self, calc: str, pot: str, value, gamma, explanation: str,
                     extra: Optional[dict] = None) -> dict:
        result = {
            "calculation_type": calc,
            "potential_type": pot,
            "explanation": explanation,
        }
        if value is not None:
            result["value"] = value
        if gamma is not None:
            result["gamma"] = gamma
        if extra:
            result.update(extra)
        return {"result": result}

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split()
            calc_type = parts[0]
            pot_type = parts[1] if len(parts) > 1 else "square"
            E_ev = float(parts[2]) if len(parts) > 2 else 1.0
            V0_ev = float(parts[3]) if len(parts) > 3 else 5.0
            w_nm = float(parts[4]) if len(parts) > 4 else 0.5
            mass = float(parts[5]) if len(parts) > 5 else 9.109e-31
            n_st = int(parts[6]) if len(parts) > 6 else 0
            return self._run_base(calc_type, pot_type, E_ev, V0_ev, w_nm, mass, n_st)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
