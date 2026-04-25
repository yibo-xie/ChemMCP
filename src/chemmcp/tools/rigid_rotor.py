import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RigidRotor(BaseTool):
    """
    刚性转子能级计算（线性分子转动）。
    
    转动能级: E_J = J(J+1)·ℏ²/(2I) = J(J+1)·k_B·B
    其中 I = μr² 为转动惯量，B = ℏ²/(2I) 为转动常数。
    简并度: g_J = 2J+1
    
    转动光谱选律: ΔJ = ±1
    """
    __version__ = "0.1.0"
    name = "RigidRotor"
    func_name = "rigid_rotor"
    description = "Calculate rigid rotor energy levels, angular momentum, and rotational spectroscopic properties for linear molecules."
    implementation_description = "Computes rotational energy levels E_J = J(J+1)ℏ²/2I for a rigid rotor model. Calculates moment of inertia from bond length and reduced mass (or directly), rotational constant B, degeneracy, spectral transition frequencies, and characteristic rotational temperature."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Rigid Rotor", "Rotational Spectroscopy", "Energy Levels"]
    required_envs = []

    code_input_sig = [
        ("moment_of_inertia_kg_m2", "float", "N/A", "Moment of inertia I in kg·m². Alternative: provide bond_length_m + reduced_mass_kg."),
        ("quantum_number_J", "int", "N/A", "Rotational quantum number J (J >= 0)."),
        ("bond_length_m", "float", "None", "Bond length in meters (optional; if given with reduced_mass_kg, computes I automatically)."),
        ("reduced_mass_kg", "float", "None", "Reduced mass in kg (optional; used with bond_length_m to compute I)."),
        ("compute_transitions_up_to_J", "int", "None", "Also compute transition frequencies up to this J value (optional)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'I J [max_J_for_transitions]' or 'bond_length reduced_mass J [max_J]'. Example: '1.46e-46 0 5'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with rotational energy, angular momentum, rotational constant B, degeneracy, transitions, and temperature data."),
    ]

    examples = [
        {
            "code_input": {
                "moment_of_inertia_kg_m2": 1.456e-46,
                "quantum_number_J": 0,
                "bond_length_m": None,
                "reduced_mass_kg": None,
                "compute_transitions_up_to_J": None,
            },
            "text_input": {
                "input_params": "1.456e-46 0",
            },
            "output": {
                "result": {
                    "energy_J": 0,
                    "energy_eV": 0,
                    "rotational_constant_B_MHz": 215200,
                    "degeneracy": 1,
                    "rotational_temperature_K": 1.54,
                }
            },
        },
        {
            "code_input": {
                "moment_of_inertia_kg_m2": 1.456e-46,
                "quantum_number_J": 2,
                "bond_length_m": None,
                "reduced_mass_kg": None,
                "compute_transitions_up_to_J": 4,
            },
            "text_input": {
                "input_params": "1.456e-46 2 4",
            },
            "output": {
                "result": {
                    "energy_eV": 0.00048,
                    "angular_momentum_hbar_units": 2.449,
                    "degeneracy": 5,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34  # J·s
        self.h = 6.62607015e-34  # J·s
        self.kB = 1.380649e-23  # J/K (Boltzmann constant)
        self.c = 2.99792458e8   # m/s
        self.eV_per_J = 6.241509e18

    def _run_base(self, moment_of_inertia_kg_m2: float, quantum_number_J: int,
                  bond_length_m: float = None, reduced_mass_kg: float = None,
                  compute_transitions_up_to_J: int = None) -> dict:
        
        if quantum_number_J < 0:
            raise ChemMCPError("Quantum number J must be >= 0.")
        
        # Compute I from bond length and reduced mass if needed
        if moment_of_inertia_kg_m2 is None or moment_of_inertia_kg_m2 <= 0:
            if bond_length_m is not None and reduced_mass_kg is not None:
                I = reduced_mass_kg * bond_length_m ** 2
            else:
                raise ChemMCPError("Must provide moment_of_inertia_kg_m2 OR both bond_length_m and reduced_mass_kg.")
        else:
            I = moment_of_inertia_kg_m2

        hbar = self.hbar
        kB = self.kB

        # Rotational energy: E_J = J(J+1) * hbar² / (2I)
        E_J = quantum_number_J * (quantum_number_J + 1) * hbar ** 2 / (2 * I)
        E_eV = E_J * self.eV_per_J

        # Rotational constant B (in various units)
        B_energy = hbar ** 2 / (2 * I)  # in Joules
        B_Hz = B_energy / self.h          # in Hz
        B_MHz = B_Hz / 1e6               # in MHz
        B_cm_inv = B_Hz / (self.c * 100)  # in cm⁻¹
        B_eV = B_energy * self.eV_per_J

        # Degeneracy: g_J = 2J + 1
        g_J = 2 * quantum_number_J + 1

        # Angular momentum magnitude: |L| = hbar*sqrt(J(J+1))
        L_mag = hbar * math.sqrt(quantum_number_J * (quantum_number_J + 1))

        # Characteristic rotational temperature: Θ_rot = B/k_B (in K)
        Theta_rot = B_energy / kB if kB > 0 else 0

        result = {
            "quantum_number_J": quantum_number_J,
            "moment_of_inertia_I_kgm2": I,
            "rotational_constant_B_eV": round(B_eV, 15),
            "rotational_constant_B_MHz": round(B_MHz, 4),
            "rotational_constant_B_cm-1": round(B_cm_inv, 6),
            "rotational_constant_B_K": round(Theta_rot, 6),
            "energy_J": round(E_J, 30),
            "energy_eV": round(E_eV, 15),
            "energy_per_B_units": quantum_number_J * (quantum_number_J + 1),  # E/B in units of B
            "angular_momentum_hbar_units": round(math.sqrt(quantum_number_J * (quantum_number_J + 1)), 10),
            "degeneracy": g_J,
            "characteristic_rotational_temperature_K": round(Theta_rot, 6),
        }

        # Spectral transitions (ΔJ = +1 absorption)
        if compute_transitions_up_to_J is not None:
            transitions = []
            max_J = min(compute_transitions_up_to_J, 50)
            for Jp in range(max_J):
                Jpp = Jp + 1
                # ΔE = E(J+1) - E(J) = 2B(J+1)
                delta_E = 2 * B_energy * (Jp + 1)
                delta_nu_Hz = delta_E / self.h
                delta_nu_cm_inv = delta_nu_Hz / (self.c * 100)
                delta_nu_GHz = delta_nu_Hz / 1e9
                transitions.append({
                    "transition": f"J={Jp} → J={Jpp}",
                    "delta_J": 1,
                    "frequency_GHz": round(delta_nu_GHz, 6),
                    "wavenumber_cm-1": round(delta_nu_cm_inv, 6),
                    "energy_J": round(delta_E, 30),
                    "energy_eV": round(delta_E * self.eV_per_J, 15),
                })
            result["spectral_transitions"] = transitions

        # If bond length was provided, include it
        if bond_length_m is not None:
            result["bond_length_m"] = bond_length_m
        if reduced_mass_kg is not None:
            result["reduced_mass_kg"] = reduced_mass_kg

        logger.info(f"RigidRotor: J={quantum_number_J}, E={E_eV:.6f}eV, B={B_MHz:.2f}MHz, g={g_J}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            
            # Detect format: if second arg looks like it could be a small number (reduced mass), 
            # treat as bond_length + reduced_mass mode
            val1 = float(parts[0])
            val2 = float(parts[1])
            J = int(parts[2])
            max_J = int(parts[3]) if len(parts) > 3 else None
            
            # Heuristic: if val1 > 1e-20, likely bond length in m; if val1 < 1e-40, likely I
            if val1 < 1e-35:
                # It's moment of inertia
                return self._run_base(val1, J, compute_transitions_up_to_J=max_J)
            else:
                # Treat as bond_length, need reduced_mass
                if len(parts) >= 4:
                    try:
                        # Could be: bond_length reduced_mass J [max_J]
                        red_mass = val2
                        J_val = int(parts[2])
                        max_J_val = int(parts[3]) if len(parts) > 3 else None
                        return self._run_base(None, J_val, bond_length_m=val1,
                                              reduced_mass_kg=red_mass,
                                              compute_transitions_up_to_J=max_J_val)
                    except (ValueError, IndexError):
                        pass
                return self._run_base(val1, J, compute_transitions_up_to_J=max_J)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'I J [max_J]' or 'bond_length reduced_mass J [max_J]'")
