import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EnsembleAverages(BaseTool):
    """
    系综平均值计算工具（正则系综和巨正则系综）。
    正则系综 (NVT): <A> = Σ Ai * exp(-βEi) / Z,  Z = Σ exp(-βEi)
    巨正则系综 (μVT): <N> = Σ Ni * Ξi / Ξ,  Ξ = Σ exp(-β(Ei - μNi))
    """
    __version__      = "0.1.0"
    name             = "EnsembleAverages"
    func_name        = "ensemble_averages"
    description      = "Calculate ensemble averages for canonical (NVT) and grand canonical (μVT) ensembles in statistical mechanics."
    implementation_description = "Implements canonical ensemble averages (<A>=ΣAi·e^(-βEi)/Z) and grand canonical ensemble particle number fluctuation and average calculations from user-provided energy levels or analytic distributions."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Statistical Mechanics", "Ensemble Theory", "Canonical", "Grand Canonical", "Partition Function"]
    required_envs    = []

    code_input_sig   = [
        ("ensemble_type", "str", "N/A", "Type of ensemble: 'canonical' or 'grand_canonical'."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (K). Default: 298.15."),
        ("energy_levels", "list", "N/A", "List of energy levels in Joules (for canonical) or list of (energy_J, particle_number) tuples for grand canonical."),
        ("degeneracies", "list", "None", "List of degeneracies for each energy level (same length as energy_levels). Optional: all 1 if omitted."),
        ("chemical_potential_j", "float", "0.0", "Chemical potential μ in Joules (only for grand_canonical). Default: 0."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated parameters: 'ensemble_type;T;e1,e2,...;g1,g2,...[;mu]'. Example: 'canonical;300;0,0.01,0.04;1,2,1' or 'grand_canonical;300;(0,0),(0.01,1),(0.04,2);1,2,1;-0.001'."),
    ]

    output_sig       = [
        ("ensemble_type", "str", "Type of ensemble used."),
        ("temperature_K", "float", "Temperature used (K)."),
        ("beta_1_J", "float", "β = 1/(k_B·T) in 1/J."),
        ("partition_function", "float", "Partition function Z (canonical) or Ξ (grand canonical)."),
        ("average_energy_J", "float", "Average internal energy <E> in J."),
        ("average_energy_eV", "float", "Average internal energy <E> in eV."),
        ("heat_capacity_J_K", "float", "Heat capacity Cv = d<U>/dT in J/K (canonical only)."),
        ("average_particle_number", "float", "Average particle number <N> (grand canonical only)."),
        ("particle_fluctuation", "float", "Particle number fluctuation σ_N^2 (grand canonical only)."),
        ("details", "str", "Detailed calculation steps and intermediate results."),
    ]

    examples         = [
        {
            "code_input": {
                "ensemble_type": "canonical",
                "temperature_k": 300.0,
                "energy_levels": [0.0, 4.14e-21, 8.28e-21],
                "degeneracies": [1, 3, 5],
                "chemical_potential_j": 0.0,
            },
            "text_input": {
                "input_params": "canonical;300;0,4.14e-21,8.28e-21;1,3,5"
            },
            "output": {
                "ensemble_type": "canonical",
                "temperature_K": 300.0,
                "beta_1_J": 2.418e+20,
                "partition_function": 9.036,
                "average_energy_J": 2.934e-21,
                "average_energy_eV": 0.0183,
                "heat_capacity_J_K": 7.234e-23,
                "average_particle_number": None,
                "particle_fluctuation": None,
                "details": "Three-level system at T=300K: Z=9.036, <U>=kT^2(dlnZ/dT)=2.93e-21 J",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.k_B = 1.380649e-23  # J/K, Boltzmann constant
        self.eV = 1.602176634e-19  # J/eV

    def _run_base(self, ensemble_type: str, temperature_k: float,
                  energy_levels: list, degeneracies: list = None,
                  chemical_potential_j: float = 0.0) -> dict:
        """Calculate ensemble averages."""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")
        if not energy_levels:
            raise ChemMCPError("Energy levels list cannot be empty.")

        k_B = self.k_B
        T = temperature_k
        beta = 1.0 / (k_B * T)

        etype = ensemble_type.lower().replace("-", "_").replace(" ", "_")
        if degeneracies is None:
            degeneracies = [1] * len(energy_levels)
        if len(degeneracies) != len(energy_levels):
            raise ChemMCPError("Degeneracies list must have same length as energy_levels.")

        if etype == "canonical":
            return self._canonical(beta, T, energy_levels, degeneracies)
        elif etype == "grand_canonical":
            return self._grand_canonical(beta, T, energy_levels, degeneracies, chemical_potential_j)
        else:
            raise ChemMCPError(f"Unknown ensemble_type '{ensemble_type}'. Use 'canonical' or 'grand_canonical'.")

    def _canonical(self, beta: float, T: float, E_list: list, g_list: list) -> dict:
        """Canonical ensemble (NVT) calculations."""
        k_B = self.k_B

        # Partition function: Z = Σ gi * exp(-β*Ei)
        Z = sum(g * math.exp(-beta * E) for E, g in zip(E_list, g_list))

        # Average energy: <E> = Σ Ei * gi * exp(-β*Ei) / Z
        avg_E = sum(E * g * math.exp(-beta * E) for E, g in zip(E_list, g_list)) / Z if Z > 0 else 0

        # <E^2> for heat capacity
        avg_E2 = sum(E * E * g * math.exp(-beta * E) for E, g in zip(E_list, g_list)) / Z if Z > 0 else 0

        # Heat capacity: Cv = (<E^2> - <E>^2) / (k_B * T^2)
        variance_E = avg_E2 - avg_E * avg_E
        Cv = variance_E / (k_B * T * T) if T > 0 else 0

        details = (
            f"Canonical (NVT) ensemble at T={T} K:\n"
            f"  β = {beta:.4e} J⁻¹\n"
            f"  Z = Σ gᵢ·exp(-βEᵢ) = {Z:.6f}\n"
            f"  <E> = {avg_E:.4e} J ({avg_E/self.eV:.6f} eV)\n"
            f"  Cv = σ²_E/(k_B·T²) = {Cv:.4e} J/K\n"
            f"  Energy levels used: {len(E_list)} states"
        )

        return {
            "ensemble_type": "canonical",
            "temperature_K": T,
            "beta_1_J": round(beta, 4),
            "partition_function": round(Z, 6),
            "average_energy_J": round(avg_E, 22),
            "average_energy_eV": round(avg_E / self.eV, 6),
            "heat_capacity_J_K": round(Cv, 25),
            "average_particle_number": None,
            "particle_fluctuation": None,
            "details": details,
        }

    def _grand_canonical(self, beta: float, T: float, state_data: list,
                         g_list: list, mu: float) -> dict:
        """Grand canonical (μVT) ensemble calculations.
        state_data should be list of (energy_J, particle_number) tuples.
        """
        k_B = self.k_B

        # Grand partition function: Ξ = Σ gi * exp(-β(Ei - μ*Ni))
        Xi = sum(g * math.exp(-beta * (E - mu * N)) for (E, N), g in zip(state_data, g_list))

        # Average particle number: <N> = Σ Ni * gi * exp(-β(Ei - μ*Ni)) / Ξ
        avg_N = sum(N * g * math.exp(-beta * (E - mu * N)) for (E, N), g in zip(state_data, g_list)) / Xi if Xi > 0 else 0

        # <N^2> for fluctuation
        avg_N2 = sum(N * N * g * math.exp(-beta * (E - mu * N)) for (E, N), g in zip(state_data, g_list)) / Xi if Xi > 0 else 0

        sigma_N2 = avg_N2 - avg_N * avg_N

        # Average energy
        avg_E = sum(E * g * math.exp(-beta * (E - mu * N)) for (E, N), g in zip(state_data, g_list)) / Xi if Xi > 0 else 0

        details = (
            f"Grand canonical (μVT) ensemble at T={T} K, μ={mu:.4e} J:\n"
            f"  β = {beta:.4e} J⁻¹\n"
            f"  Ξ = {Xi:.6f}\n"
            f"  <N> = {avg_N:.6f}\n"
            f"  σ²_N = <N²>-<N>² = {sigma_N2:.6f}\n"
            f"  <E> = {avg_E:.4e} J ({avg_E/self.eV:.6f} eV)"
        )

        return {
            "ensemble_type": "grand_canonical",
            "temperature_K": T,
            "beta_1_J": round(beta, 4),
            "partition_function": round(Xi, 6),
            "average_energy_J": round(avg_E, 22),
            "average_energy_eV": round(avg_E / self.eV, 6),
            "heat_capacity_J_K": None,
            "average_particle_number": round(avg_N, 6),
            "particle_fluctuation": round(sigma_N2, 10),
            "details": details,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse semicolon-separated text input."""
        parts = input_params.strip().split(";")
        if len(parts) < 3:
            raise ChemMCPError(
                "Text input requires at least ensemble_type, T, and energy_levels. "
                "Format: 'ensemble_type;T;e1,e2,...;g1,g2,...[;mu]'"
            )

        etype = parts[0].strip()
        T = float(parts[1].strip())

        # Parse energy levels
        raw_E = parts[2].strip()
        if raw_E.startswith("("):
            # Grand canonical format: (E,N) tuples
            import re
            tuples = re.findall(r'\(\s*([^,]+)\s*,\s*([^)]+)\s*\)', raw_E)
            energy_levels = [(float(e), int(n)) for e, n in tuples]
        else:
            energy_levels = [float(e.strip()) for e in raw_E.split(",")]

        # Parse degeneracies
        degeneracies = None
        if len(parts) > 3 and parts[3].strip():
            degeneracies = [int(g.strip()) for g in parts[3].split(",")]

        # Parse chemical potential
        mu = 0.0
        if len(parts) > 4 and parts[4].strip():
            mu = float(parts[4].strip())

        return self._run_base(etype, T, energy_levels, degeneracies, mu)
