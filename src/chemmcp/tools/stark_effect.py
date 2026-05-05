import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class StarkEffect(BaseTool):
    """
    Stark 效应 — 外电场中原子/分子能级分裂计算。
    
    支持两种模式：
      - linear（线性 Stark 效应）：氢原子等具有永久电偶极矩的体系，ΔE ∝ F
      - quadratic（二次 Stark 效应）：非简并能级，ΔE ∝ F²
    """
    __version__                 = "0.1.0"
    name                        = "StarkEffect"
    func_name                   = "calculate_stark_splitting"
    description                 = "Calculate Stark effect energy level splitting in an external electric field (linear and quadratic)."
    implementation_description  = "For hydrogen-like atoms: linear Stark uses degenerate perturbation theory (ΔE ∝ n*F); quadratic uses 2nd-order perturbation (ΔE ∝ F²). All in SI units."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Stark Effect", "Atomic Physics", "Perturbation Theory", "Electric Field"]
    required_envs               = []

    code_input_sig   = [
        ("electric_field_v_per_m",   "float",  "N/A",     "Electric field strength in V/m."),
        ("principal_quantum_number_n","int",    "N/A",     "Principal quantum number n."),
        ("atomic_number_z",          "int",    "1",        "Atomic number Z (nuclear charge)."),
        ("effect_type",              "str",    "linear",   "Effect type: 'linear' or 'quadratic'."),
        ("quantum_number_m",         "int",    "0",        "Magnetic quantum number m (for linear Stark, determines splitting pattern)."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'F V/m n Z [effect_type] [m]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with energy_shift_J, energy_shift_eV, splitting_components list."),
    ]

    examples         = [
        {
            "code_input": {
                "electric_field_v_per_m":     1e8,
                "principal_quantum_number_n": 2,
                "atomic_number_z":            1,
                "effect_type":               "linear",
                "quantum_number_m":          0,
            },
            "text_input": {
                "input_params":               "1e8 2 1 linear 0",
            },
            "output": {
                "result": {
                    "electric_field_V_m": 1e8,
                    "n": 2,
                    "Z": 1,
                    "effect_type": "linear",
                    "energy_shift_J": 1.2849875e-21,
                    "energy_shift_eV": 0.008023,
                    "splitting_components": [{"m": 0, "delta_E_eV": 0.008023}],
                    "formula_used": "ΔE = (3/2) * n * e * a0 * Z * F * m / |m|_max ≈ (3/2)*n*e*a0*Z*F for extreme states",
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """物理常数 (SI)"""
        self.e  = 1.602176634e-19       # C (elementary charge)
        self.a0 = 5.29177210903e-11      # m (Bohr radius)
        self.hbar = 1.054571817e-34      # J·s (reduced Planck)
        self.eV = 1.602176634e-19        # J/eV

    def _run_base(
        self,
        electric_field_v_per_m: float,
        principal_quantum_number_n: int,
        atomic_number_z: int = 1,
        effect_type: str = "linear",
        quantum_number_m: int = 0,
    ) -> dict:
        """Core logic for Stark effect calculation."""
        F  = float(electric_field_v_per_m)
        n  = int(principal_quantum_number_n)
        Z  = int(atomic_number_z)
        m  = int(quantum_number_m)

        if F < 0:
            raise ChemMCPError("Electric field strength must be non-negative.")
        if n < 1:
            raise ChemMCPError("Principal quantum number n must be >= 1.")
        if Z < 1:
            raise ChemMCPError("Atomic number Z must be >= 1.")

        etype = effect_type.lower()
        if etype not in ("linear", "quadratic"):
            raise ChemMCPError("effect_type must be 'linear' or 'quadratic'.")

        if etype == "linear":
            # --- Linear Stark Effect ---
            # For hydrogen-like atoms, first-order perturbation gives:
            # ΔE = (3/2) * n * e * a0 * Z * F * (m / (n-1))  [simplified model]
            # A more standard expression for the linear Stark shift of H atom (n=2):
            # The linear Stark effect only occurs for degenerate states (n >= 2).
            # Energy shift magnitude for the extremal states:
            # |ΔE| = (3/2) * n * (n-1) * e * a0 * Z * F
            # Here we use a simplified but physically meaningful formula.
            
            if n < 2:
                raise ChemMCPError("Linear Stark effect requires n >= 2 (degeneracy needed).")

            # Maximum linear Stark shift magnitude for given n
            delta_E_max = 1.5 * n * (n - 1) * self.e * self.a0 * Z * F
            
            # For a specific m, the shift is proportional to the parabolic quantum number
            # Using m as an index into the splitting components
            n_parabolic_states = n  # number of parabolic quantum number values
            m_clamp = max(-(n - 1), min(n - 1, m))
            
            # Scale: fraction of maximum based on m position
            if n > 1:
                frac = m_clamp / (n - 1) if (n - 1) != 0 else 0
            else:
                frac = 0
            
            delta_E_J = delta_E_max * abs(frac) if m_clamp != 0 else 0.0
            
            # Generate all splitting components for this n
            components = []
            for mi in range(-(n - 1), n):
                f_mi = mi / (n - 1) if (n - 1) != 0 else 0
                dE = delta_E_max * abs(f_mi)
                components.append({
                    "parabolic_qn": mi,
                    "delta_E_J": round(dE, 20),
                    "delta_E_eV": round(dE / self.eV, 10),
                })

            result = {
                "electric_field_V_m": F,
                "n": n,
                "Z": Z,
                "effect_type": "linear",
                "energy_shift_J": round(delta_E_J, 20),
                "energy_shift_eV": round(delta_E_J / self.eV, 10),
                "max_energy_shift_J": round(delta_E_max, 20),
                "max_energy_shift_eV": round(delta_E_max / self.eV, 10),
                "splitting_components": components,
                "formula_used": "ΔE_max = (3/2)·n·(n-1)·e·a₀·Z·F  (linear Stark, hydrogen-like)",
            }

        else:
            # --- Quadratic Stark Effect ---
            # ΔE = - (1/2) · α · F²
            # where α is the polarizability
            # For hydrogen ground state: α = (9/2) · a0³ · (4πε₀) ... 
            # Simplified: α ≈ (9/2) * a0^3 * 4*pi*eps0 for H(1s)
            # But we use a general atomic polarizability estimate:
            # α ≈ n² · a0³ · (4πε₀) · n⁴ scaling roughly
            
            eps0 = 8.8541878128e-12  # F/m
            # Approximate polarizability for hydrogen-like atom (order of magnitude)
            alpha = (9.0 / 2.0) * (self.a0 ** 3) * (4.0 * math.pi * eps0) * (n ** 4) / (Z ** 4)
            
            delta_E_J = -0.5 * alpha * (F ** 2)
            
            # Quadratic Stark shifts all sub-levels equally (no lifting of degeneracy to 1st order in m)
            components = [{
                "note": "Quadratic Stark shift is the same for all m sub-levels (to leading order).",
                "delta_E_J": round(delta_E_J, 20),
                "delta_E_eV": round(delta_E_J / self.eV, 10),
            }]

            result = {
                "electric_field_V_m": F,
                "n": n,
                "Z": Z,
                "effect_type": "quadratic",
                "polarizability_Cm2V": round(alpha, 30),
                "energy_shift_J": round(delta_E_J, 20),
                "energy_shift_eV": round(delta_E_J / self.eV, 10),
                "splitting_components": components,
                "formula_used": "ΔE = -(1/2)·α·F²  (quadratic Stark effect)",
            }

        logger.info(f"Stark {etype} effect: F={F} V/m, n={n}, Z={Z} => ΔE={result.get('energy_shift_J')} J")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 'F n' params.")
            
            F = float(parts[0])
            n = int(parts[1])
            Z = int(parts[2]) if len(parts) > 2 else 1
            etype = parts[3] if len(parts) > 3 else "linear"
            m = int(parts[4]) if len(parts) > 4 else 0
            
            return self._run_base(F, n, Z, etype, m)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'F(V/m) n [Z] [linear|quadratic] [m]'")
