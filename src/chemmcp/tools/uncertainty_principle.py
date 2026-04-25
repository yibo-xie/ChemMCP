import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class UncertaintyPrinciple(BaseTool):
    """
    海森堡不确定性原理数值演示工具。
    Δx · Δp ≥ ℏ/2,  ΔE · Δt ≥ ℏ/2
    """
    __version__ = "0.1.0"
    name = "UncertaintyPrinciple"
    func_name = "uncertainty_principle_demo"
    description = "Numerical demonstration of Heisenberg's Uncertainty Principle for position-momentum and energy-time."
    implementation_description = "Uses fundamental constants (h-bar) to compute minimum uncertainty products and demonstrate quantum limits for various physical systems."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Uncertainty Principle", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("system_type", "str", "N/A", "Type of system: 'position_momentum', 'energy_time', or 'electron'."),
        ("delta_x", "float", "N/A", "Position uncertainty in meters (for position_momentum type)."),
        ("mass_kg", "float", "N/A", "Particle mass in kg (for position_momentum type)."),
        ("delta_e", "float", "N/A", "Energy uncertainty in Joules (for energy_time type)."),
        ("delta_t", "float", "N/A", "Time uncertainty in seconds (for energy_time type)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated parameters depending on system_type."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing uncertainty calculation results including minimum uncertainty, actual product, and whether Heisenberg limit is satisfied."),
    ]

    examples = [
        {
            "code_input": {
                "system_type": "position_momentum",
                "delta_x": 1e-10,
                "mass_kg": 9.109e-31,
                "delta_e": 0,
                "delta_t": 0,
            },
            "text_input": {
                "input_params": "position_momentum 1e-10 9.109e-31",
            },
            "output": {
                "result": {
                    "system": "position_momentum",
                    "delta_x_m": 1e-10,
                    "min_delta_p_kg_m_s": 5.273e-25,
                    "actual_delta_p_kg_m_s": 5.273e-25,
                    "product_J_s": 5.273e-34,
                    "hbar_over_2_J_s": 5.273e-34,
                    "at_limit": True,
                    "interpretation": "The uncertainty product equals the quantum limit.",
                }
            }
        },
        {
            "code_input": {
                "system_type": "energy_time",
                "delta_x": 0,
                "mass_kg": 0,
                "delta_e": 1e-20,
                "delta_t": 1e-15,
            },
            "text_input": {
                "input_params": "energy_time 0 0 1e-20 1e-15",
            },
            "output": {
                "result": {
                    "system": "energy_time",
                    "delta_e_J": 1e-20,
                    "delta_t_s": 1e-15,
                    "product_J_s": 1e-35,
                    "hbar_over_2_J_s": 5.273e-35,
                    "above_limit": True,
                    "interpretation": "Product exceeds minimum; satisfies uncertainty principle.",
                }
            }
        },
        {
            "code_input": {
                "system_type": "electron",
                "delta_x": 1e-9,
                "mass_kg": 9.109e-31,
                "delta_e": 0,
                "delta_t": 0,
            },
            "text_input": {
                "input_params": "electron 1e-9 9.109e-31",
            },
            "output": {
                "result": {
                    "system": "electron_position_momentum",
                    "particle": "electron",
                    "mass_kg": 9.109e-31,
                    "delta_x_m": 1e-9,
                    "min_delta_v_m_s": 578700.7,
                    "interpretation": "Electron confined to 1 nm has minimum velocity uncertainty ~579 km/s.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34  # J·s, reduced Planck constant
        self.h = 6.62607015e-34       # J·s, Planck constant
        self.electron_mass = 9.1093837015e-31  # kg

    def _run_base(self, system_type: str, delta_x: float = 0.0, mass_kg: float = 0.0,
                  delta_e: float = 0.0, delta_t: float = 0.0) -> dict:
        """Core logic for uncertainty principle calculations."""
        st = system_type.lower().replace("-", "_")

        if st in ("position_momentum", "position_momentum"):
            result = self._calc_position_momentum(delta_x, mass_kg)
        elif st in ("energy_time", "energy_time"):
            result = self._calc_energy_time(delta_e, delta_t)
        elif st == "electron":
            result = self._electron_demo(delta_x, mass_kg)
        else:
            raise ChemMCPError(
                f"Unknown system_type '{system_type}'. "
                "Use 'position_momentum', 'energy_time', or 'electron'."
            )
        return {"result": result}

    def _calc_position_momentum(self, delta_x: float, mass_kg: float) -> dict:
        """Δx · Δp ≥ ℏ/2"""
        if delta_x <= 0:
            raise ChemMCPError("delta_x must be positive.")
        if mass_kg <= 0:
            raise ChemMCPError("mass_kg must be positive.")

        min_dp = self.hbar / (2 * delta_x)
        product = delta_x * min_dp
        limit = self.hbar / 2

        # Also compute minimum velocity uncertainty
        min_dv = min_dp / mass_kg

        return {
            "system": "position_momentum",
            "delta_x_m": delta_x,
            "mass_kg": mass_kg,
            "min_delta_p_kg_m_s": f"{min_dp:.4e}",
            "min_delta_v_m_s": f"{min_dv:.4e}",
            "product_J_s": f"{product:.4e}",
            "hbar_over_2_J_s": f"{limit:.4e}",
            "at_limit": abs(product - limit) < limit * 1e-6,
            "formula": "Δx · Δp ≥ ℏ/2",
            "interpretation": (
                f"When confined to Δx = {delta_x:.2e} m, "
                f"minimum momentum uncertainty is {min_dp:.3e} kg·m/s "
                f"(velocity uncertainty ≥ {min_dv:.3e} m/s for m={mass_kg:.3e} kg)."
            ),
        }

    def _calc_energy_time(self, delta_e: float, delta_t: float) -> dict:
        """ΔE · Δt ≥ ℏ/2"""
        if delta_e < 0:
            raise ChemMCPError("delta_e must be non-negative.")
        if delta_t <= 0:
            raise ChemMCPError("delta_t must be positive.")

        product = delta_e * delta_t
        limit = self.hbar / 2

        return {
            "system": "energy_time",
            "delta_e_J": delta_e,
            "delta_t_s": delta_t,
            "product_J_s": f"{product:.4e}",
            "hbar_over_2_J_s": f"{limit:.4e}",
            "satisfies_principle": product >= limit - limit * 1e-10,
            "fraction_of_minimum": round(product / limit, 6) if limit > 0 else float("inf"),
            "formula": "ΔE · Δt ≥ ℏ/2",
            "interpretation": (
                f"ΔE·Δt = {product:.3e} J·s, minimum required = {limit:.3e} J·s. "
                f"Satisfies principle: {product >= limit - limit * 1e-10}"
            ),
        }

    def _electron_demo(self, delta_x: float, mass_kg: float = None) -> dict:
        """Specialized electron uncertainty demo."""
        if mass_kg is None or mass_kg <= 0:
            mass_kg = self.electron_mass
        if delta_x <= 0:
            raise ChemMCPError("delta_x must be positive.")

        min_dp = self.hbar / (2 * delta_x)
        min_dv = min_dp / mass_kg
        min_ke = 0.5 * mass_kg * min_dv ** 2  # Joules
        min_ke_ev = min_ke / 1.602176634e-19   # eV

        return {
            "system": "electron_position_momentum",
            "particle": "electron",
            "mass_kg": mass_kg,
            "delta_x_m": delta_x,
            "min_delta_p_kg_m_s": f"{min_dp:.4e}",
            "min_delta_v_m_s": f"{min_dv:.4e}",
            "min_kinetic_energy_J": f"{min_ke:.4e}",
            "min_kinetic_energy_eV": round(min_ke_ev, 6),
            "formula": "Δx · Δp ≥ ℏ/2 → Δv ≥ ℏ/(2m·Δx)",
            "interpretation": (
                f"Electron confined to Δx = {delta_x:.2e} m has:\n"
                f"  • Minimum momentum uncertainty: {min_dp:.3e} kg·m/s\n"
                f"  • Minimum velocity uncertainty: {min_dv:.3e} m/s\n"
                f"  • Minimum kinetic energy: {min_ke_ev:.3f} eV"
            ),
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            st = parts[0]
            if st in ("position_momentum", "position_momentum"):
                return self._run_base(st, float(parts[1]), float(parts[2]))
            elif st in ("energy_time", "energy_time"):
                return self._run_base(st, 0, 0, float(parts[3]), float(parts[4]))
            elif st == "electron":
                dx = float(parts[1])
                m = float(parts[2]) if len(parts) > 2 else None
                return self._run_base("electron", dx, m or 9.109e-31)
            else:
                raise ValueError(f"Unknown system type: {st}")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
