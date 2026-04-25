import logging
import math
from typing import List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
H = 6.62607015e-34        # J·s
HBAR = 1.054571817e-34     # J·s
C = 2.99792458e8           # m/s
KB = 1.380649e-23          # J/K (Boltzmann constant)
AMU = 1.66053906660e-27    # kg (atomic mass unit)


@ChemMCPManager.register_tool
class RotationalSpectrum(BaseTool):
    """
    转动光谱线位置和强度计算工具。
    计算刚性转子模型的转动跃迁能级、谱线位置（频率/波数/波长）和相对强度。
    支持线型分子和对称陀螺。
    """
    __version__ = "0.1.0"
    name = "RotationalSpectrum"
    func_name = "calculate_rotational_spectrum"
    description = "Calculate rotational spectroscopy transition energies, line positions (frequency/wavenumber/wavelength), and relative intensities for rigid rotor models."
    implementation_description = "Uses quantum mechanical rigid rotor model: E_J = B·J(J+1) for linear molecules, with selection rule ΔJ = ±1. Computes rotational constants from molecular structure and Boltzmann-distributed line intensities."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "Rotational", "Microwave", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("molecule_type", "str", "N/A", "Molecule type: 'linear' or 'symmetric_top'."),
        ("reduced_mass_amu", "float", "N/A", "Reduced mass of the molecule in atomic mass units (amu)."),
        ("bond_length_angstrom", "float", "N/A", "Bond length in Angstroms (for linear molecules)."),
        ("moment_of_inertia_kgm2", "float", "0", "Direct moment of inertia in kg·m² (overrides reduced_mass + bond_length if > 0)."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin for intensity calculation."),
        ("max_j", "int", "10", "Maximum rotational quantum number J to compute."),
        ("isotope_masses_amu", "list", "[]", "Atomic masses [m1, m2] in amu for computing reduced mass directly."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: molecule_type reduced_mass_amu bond_length_A [T_K max_j]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing rotational constant, energy levels, transition lines with positions/intensities, and spectrum summary."),
    ]

    examples = [
        {
            "code_input": {
                "molecule_type": "linear",
                "reduced_mass_amu": 0,
                "bond_length_angstrom": 1.128,
                "moment_of_inertia_kgm2": 0,
                "temperature_k": 298.15,
                "max_j": 10,
                "isotope_masses_amu": [12.0, 1.008],
            },
            "text_input": {
                "input_params": "linear 0 1.128 298.15 10 12.0 1.008",
            },
            "output": {
                "result": {
                    "molecule": "¹²C-¹H diatomic (linear)",
                    "rotational_constant_B_cm-1": round(H / (8 * math.pi**2 * C * (12.0*1.008/(12.0+1.008)) * AMU * (1.128e-10)**2) / 100, 3),
                    "moment_of_inertia_kg_m2": round((12.0*1.008/(12.0+1.008)) * AMU * (1.128e-10)**2, 35),
                    "transitions": [],
                    "temperature_k": 298.15,
                    "notes": "CO-like parameters used as example.",
                }
            }
        },
        {
            "code_input": {
                "molecule_type": "linear",
                "reduced_mass_amu": 1.0,
                "bond_length_angstrom": 1.0,
                "moment_of_inertia_kgm2": 0,
                "temperature_k": 300.0,
                "max_j": 6,
                "isotope_masses_amu": [],
            },
            "text_input": {
                "input_params": "linear 1.0 1.0 300 6",
            },
            "output": {
                "result": {
                    "rotational_constant_B_cm-1": round(H / (8 * math.pi**2 * C * 1.0 * AMU * (1.0e-10)**2) / 100, 3),
                    "max_j_computed": 6,
                    "num_transitions": 6,
                    "spectrum_range_cm-1": "...",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = H
        self.hbar = HBAR
        self.c = C
        self.kb = KB
        self.amu = AMU

    def _run_base(self, molecule_type: str, reduced_mass_amu: float = 0.0,
                  bond_length_angstrom: float = 0.0, moment_of_inertia_kgm2: float = 0.0,
                  temperature_k: float = 298.15, max_j: int = 10,
                  isotope_masses_amu: List[float] = None) -> dict:
        """Core logic."""
        mt = molecule_type.lower().replace(" ", "_")
        if mt not in ("linear", "symmetric_top"):
            raise ChemMCPError(f"Unsupported molecule_type '{molecule_type}'. Use 'linear' or 'symmetric_top'.")

        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive.")
        if max_j < 1:
            raise ChemMCPError("max_j must be >= 1.")

        # Compute moment of inertia
        if moment_of_inertia_kgm2 > 0:
            I = moment_of_inertia_kgm2
        elif isotope_masses_amu and len(isotope_masses_amu) >= 2:
            mu = (isotope_masses_amu[0] * isotope_masses_amu[1]) / sum(isotope_masses_amu)
            r = bond_length_angstrom * 1e-10
            I = mu * self.amu * r ** 2
        elif reduced_mass_amu > 0 and bond_length_angstrom > 0:
            mu = reduced_mass_amu
            r = bond_length_angstrom * 1e-10
            I = mu * self.amu * r ** 2
        else:
            raise ChemMCPError(
                "Must provide either moment_of_inertia_kgm2, "
                "or both reduced_mass_amu + bond_length_angstrom, "
                "or isotope_masses_amu + bond_length_angstrom."
            )

        # Rotational constant B (in cm⁻¹): B = h/(8π²cI)
        B_cm = self.h / (8 * math.pi ** 2 * self.c * I)   # in m⁻¹
        B_cm /= 100                                          # convert to cm⁻¹
        B_Hz = self.h * B_cm * 100 * self.c                 # B in Hz
        B_J = self.h * B_Hz                                  # B in Joules

        # Energy levels E_J = B · J(J+1) (in cm⁻¹)
        levels = []
        for j in range(max_j + 1):
            E_j = B_cm * j * (j + 1)
            levels.append({"J": j, "energy_cm-1": round(E_j, 4), "energy_J": round(E_j * 100 * self.h * self.c, 20)})

        # Transitions: ΔJ = +1 (absorption), J → J+1
        # ν̃(J→J+1) = 2B(J+1)
        transitions = []
        total_pop = 0
        populations = []

        for j in range(max_j):
            # Boltzmann population of level J
            g_j = 2 * j + 1  # degeneracy
            E_j_J = B_J * j * (j + 1)
            pop = g_j * math.exp(-E_j_J / (self.kb * temperature_k))
            populations.append(pop)
            total_pop += pop

        # Normalize populations
        if total_pop > 0:
            populations = [p / total_pop for p in populations]

        for j in range(min(max_j, len(populations))):
            nu_tilde = 2 * B_cm * (j + 1)  # wavenumber of transition
            freq_hz = nu_tilde * 100 * self.c
            wl_m = self.c / freq_hz if freq_hz > 0 else float("inf")
            intensity = populations[j] * (j + 1)  # proportional to population × degeneracy factor

            transitions.append({
                "transition": f"J={j} → J={j+1}",
                "J_lower": j,
                "J_upper": j + 1,
                "wavenumber_cm-1": round(nu_tilde, 4),
                "frequency_GHz": round(freq_hz / 1e9, 4),
                "wavelength_mm": round(wl_m * 1000, 6),
                "relative_intensity": round(intensity, 6),
                "population_J": round(populations[j], 6),
                "degeneracy_g": 2 * j + 1,
            })

        # Find most intense transition
        max_trans = max(transitions, key=lambda t: t["relative_intensity"]) if transitions else None

        result = {
            "molecule_type": molecule_type,
            "moment_of_inertia_kg_m2": f"{I:.4e}",
            "rotational_constant_B_cm-1": round(B_cm, 6),
            "rotational_constant_B_GHz": round(B_Hz / 1e9, 4),
            "rotational_constant_B_J": f"{B_J:.4e}",
            "temperature_K": temperature_k,
            "energy_levels": levels,
            "transitions": transitions,
            "num_transitions": len(transitions),
            "spectrum_range_cm-1": f"{transitions[-1]['wavenumber_cm-1']}-{transitions[0]['wavenumber_cm-1']}" if transitions else "N/A",
            "most_intense_transition": max_trans["transition"] if max_trans else None,
            "selection_rule": "ΔJ = ±1 (absorption: J → J+1)",
            "spacing": f"Equally spaced by 2B = {round(2*B_cm, 4)} cm⁻¹ ({round(2*B_Hz/1e9, 4)} GHz)",
        }

        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            mol_type = parts[0]
            mu = float(parts[1]) if len(parts) > 1 else 0.0
            bl = float(parts[2]) if len(parts) > 2 else 0.0
            T = float(parts[3]) if len(parts) > 3 else 298.15
            mj = int(parts[4]) if len(parts) > 4 else 10
            masses = [float(p) for p in parts[5:]] if len(parts) > 5 else []
            return self._run_base(mol_type, mu, bl, 0, T, mj, masses)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
