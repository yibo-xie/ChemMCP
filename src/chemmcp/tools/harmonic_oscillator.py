import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HarmonicOscillator(BaseTool):
    """
    量子谐振子模型的能级和波函数计算。
    
    能级: E_v = (v + 1/2) ℏω,  where ω = sqrt(k/μ)
    波函数: ψ_v(ξ) = N_v · H_v(ξ) · exp(-ξ²/2),  ξ = sqrt(μω/ℏ)·x
    
    其中 H_v 是 Hermite 多项式，N_v 是归一化常数。
    """
    __version__ = "0.1.0"
    name = "HarmonicOscillator"
    func_name = "harmonic_oscillator"
    description = "Calculate energy eigenvalues, wavefunction properties, and classical parameters for quantum harmonic oscillator."
    implementation_description = "Solves the quantum harmonic oscillator model. Computes energy levels (E_v = (v+½)ℏω), turning points, expectation values, Hermite polynomial order, and zero-point energy. Supports both force constant and angular frequency inputs."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Harmonic Oscillator", "Vibration", "Energy Levels"]
    required_envs = []

    code_input_sig = [
        ("mass_kg", "float", "N/A", "Reduced mass of the system in kg."),
        ("quantum_number_v", "int", "N/A", "Vibrational quantum number v (v >= 0)."),
        ("force_constant_N_m", "float", "N/A", "Force constant k in N/m. Alternative: provide angular_frequency_rad_s instead."),
        ("angular_frequency_rad_s", "float", "None", "Angular frequency ω in rad/s (optional; if given, overrides force_constant)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'mass_kg k_or_omega v [omega_flag]'. Example: '1.63e-27 480 0' for reduced mass, k=480N/m, v=0"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with energy level, angular frequency, turning points, zero-point energy, wavefunction info, and spectroscopic constants."),
    ]

    examples = [
        {
            "code_input": {
                "mass_kg": 1.627e-27,
                "force_constant_N_m": 480.0,
                "quantum_number_v": 0,
                "angular_frequency_rad_s": None,
            },
            "text_input": {
                "input_params": "1.627e-27 480 0",
            },
            "output": {
                "result": {
                    "energy_J": 3.265e-20,
                    "energy_eV": 0.2038,
                    "angular_frequency_rad_s": 5.434e14,
                    "zero_point_energy_eV": 0.2038,
                    "turning_points_pm": 1.089e-11,
                    "hermite_polynomial_order": 0,
                    "vibrational_frequency_cm-1": 2874,
                }
            },
        },
        {
            "code_input": {
                "mass_kg": 1.627e-27,
                "force_constant_N_m": 480.0,
                "quantum_number_v": 1,
                "angular_frequency_rad_s": None,
            },
            "text_input": {
                "input_params": "1.627e-27 480 1",
            },
            "output": {
                "result": {
                    "energy_eV": 0.611,
                    "transition_energy_0to1_eV": 0.407,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34  # J·s (reduced Planck constant)
        self.h = 6.62607015e-34  # J·s
        self.c = 2.99792458e8  # m/s
        self.eV_per_J = 6.241509e18

    def _hermite(self, n, x):
        """Compute Hermite polynomial H_n(x) using recursion relation."""
        if n == 0:
            return 1.0
        elif n == 1:
            return 2.0 * x
        else:
            H_prev2 = 1.0  # H_0
            H_prev1 = 2.0 * x  # H_1
            for i in range(2, n + 1):
                H_curr = 2.0 * x * H_prev1 - 2.0 * (i - 1) * H_prev2
                H_prev2 = H_prev1
                H_prev1 = H_curr
            return H_prev1

    def _run_base(self, mass_kg: float, quantum_number_v: int,
                  force_constant_N_m: float = None,
                  angular_frequency_rad_s: float = None) -> dict:
        
        if quantum_number_v < 0:
            raise ChemMCPError("Quantum number v must be >= 0.")
        
        # Determine omega
        if angular_frequency_rad_s is not None and angular_frequency_rad_s > 0:
            omega = angular_frequency_rad_s
            k_calc = mass_kg * omega ** 2
        elif force_constant_N_m is not None and force_constant_N_m > 0:
            k_calc = force_constant_N_m
            omega = math.sqrt(force_constant_N_m / mass_kg)
        else:
            raise ChemMCPError("Must provide either force_constant_N_m or angular_frequency_rad_s.")

        hbar = self.hbar
        
        # Energy: E_v = (v + 1/2) * hbar * omega
        E = (quantum_number_v + 0.5) * hbar * omega
        E_eV = E * self.eV_per_J
        
        # Zero-point energy
        E_zp = 0.5 * hbar * omega
        E_zp_eV = E_zp * self.eV_per_J
        
        # Turning points: classically, E = (1/2)kx² => x_max = sqrt(2E/k)
        x_turn = math.sqrt(2 * E / k_calc) if k_calc > 0 else 0
        
        # Characteristic length scale: alpha = sqrt(m*omega/hbar) => xi = alpha*x
        alpha = math.sqrt(mass_kg * omega / hbar)
        
        # Vibrational frequency in wavenumbers (cm⁻¹): ν̃ = ω/(2πc)
        wavenumber = omega / (2 * math.pi * self.c) * 100  # convert m⁻¹ to cm⁻¹
        
        # Transition energy ΔE = hbar*omega (between adjacent levels)
        delta_E_eV = hbar * omega * self.eV_per_J
        
        # Expectation values for QHO: <x> = 0 always, <x²> = (v+1/2)*hbar/(m*omega)
        exp_x2 = (quantum_number_v + 0.5) * hbar / (mass_kg * omega)
        
        # Evaluate wavefunction at center (xi=0)
        # psi_v(0) = N_v * H_v(0) * exp(0), N_v = (alpha/sqrt(pi))^(1/2) * 1/sqrt(2^v * v!)
        N_v = (alpha / math.sqrt(math.pi)) ** 0.5 / math.sqrt(2 ** int(quantum_number_v) * math.factorial(int(quantum_number_v)))
        H_at_0 = self._hermite(quantum_number_v, 0.0)
        psi_center = N_v * H_at_0
        prob_center = psi_center ** 2

        result = {
            "quantum_number_v": quantum_number_v,
            "mass_kg": mass_kg,
            "force_constant_N_m": round(k_calc, 6),
            "angular_frequency_omega_rad_s": round(omega, 6),
            "frequency_Hz": round(omega / (2 * math.pi), 4),
            "vibrational_frequency_cm-1": round(wavenumber, 4),
            "energy_J": round(E, 25),
            "energy_eV": round(E_eV, 10),
            "zero_point_energy_J": round(E_zp, 25),
            "zero_point_energy_eV": round(E_zp_eV, 10),
            "turning_points_pm_m": round(x_turn, 25),
            "characteristic_length_alpha_1/m": round(alpha, 6),
            "expectation_x2_m2": round(exp_x2, 35),
            "wavefunction_hermite_order": quantum_number_v,
            "psi_at_center": round(psi_center, 20),
            "probability_density_at_center": round(prob_center, 30),
            "transition_energy_adjacent_eV": round(delta_E_eV, 10),
            "parity": "even" if quantum_number_v % 2 == 0 else "odd",
        }

        logger.info(f"HarmonicOscillator: v={quantum_number_v}, E={E_eV:.6f}eV, ν̃={wavenumber:.1f}cm⁻¹")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mass = float(parts[0])
            k_or_omega = float(parts[1])
            v = int(parts[2])
            omega_flag = parts[3].lower() if len(parts) > 3 else None
            
            if omega_flag == "omega" or omega_flag == "1":
                return self._run_base(mass, v, angular_frequency_rad_s=k_or_omega)
            else:
                return self._run_base(mass, v, force_constant_N_m=k_or_omega)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'mass k_or_omega v [omega_flag]'")
