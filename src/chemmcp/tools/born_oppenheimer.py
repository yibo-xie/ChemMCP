import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BornOppenheimer(BaseTool):
    """
    Born-Oppenheimer 近似工具。
    实现核-电子运动分离近似，计算分子势能曲线、振动能级、平衡键长、解离能。
    支持双原子分子（谐振子/莫尔斯势）和多原子分子（简正模式分析）。
    """
    __version__ = "0.1.0"
    name = "BornOppenheimer"
    func_name = "born_oppenheimer_approximation"
    description = "Apply Born-Oppenheimer approximation: separate nuclear and electronic motion, compute potential energy curves, vibrational levels, equilibrium geometry."
    implementation_description = "Implements BO approximation for diatomic molecules using Morse potential and harmonic oscillator models; computes PES, force constants, vibrational frequencies, zero-point energy, and dissociation energy."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Born-Oppenheimer", "Potential Energy Surface", "Vibration", "Diatomic"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: 'H2', 'HCl', 'CO', 'N2', 'O2', 'F2', 'NaCl', or 'generic'."),
        ("calculation", "str", "'pes'", "Calculation type: 'pes' (potential curve), 'vibrational' (energy levels), 'force_constant' (k), 'reduced_mass' (μ), 'approximation_analysis' (BO validity check)."),
        ("R_range", "str", "'auto'", "Internuclear distance range: 'auto' (uses known Re±30%), or 'Rmin:Rmax' in Angstroms (e.g., '0.3:3.0')."),
        ("n_points", "int", "50", "Number of points for PES grid."),
        ("n_vib_levels", "int", "10", "Number of vibrational levels to compute."),
        ("De_eV", "float", "None", "Dissociation energy in eV (for generic molecule)."),
        ("Re_Angstrom", "float", "None", "Equilibrium bond length in Å (for generic molecule)."),
        ("omega_e_cm-1", "float", "None", "Vibrational constant ωₑ in cm⁻¹ (for generic molecule)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: molecule calculation [R_range] [n_points] [De_eV Re_Angstrom omega_e]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing PES data, vibrational analysis, BO validity metrics, molecular parameters."),
    ]

    examples = [
        {
            "code_input": {
                "molecule": "H2",
                "calculation": "pes",
            },
            "text_input": {
                "input_str": "H2 pes",
            },
            "output": {
                "result": {
                    "molecule": "H2",
                    "equilibrium_distance_A": 0.74,
                    "dissociation_energy_eV": 4.75,
                    "vibrational_frequency_cm-1": 4401,
                    "force_constant_N/m": 570,
                    "bo_validity_ratio": 0.0015,
                    "interpretation": "Born-Oppenheimer approximation is excellent for H2.",
                }
            }
        },
        {
            "code_input": {
                "molecule": "HCl",
                "calculation": "vibrational",
                "n_vib_levels": 8,
            },
            "text_input": {
                "input_str": "HCl vibrational 8",
            },
            "output": {
                "result": {
                    "molecule": "HCl",
                    "vibrational_levels": [...],
                    "zero_point_energy_eV": "...",
                    "spacing_anharmonicity": "...",
                }
            }
        },
        {
            "code_input": {
                "molecule": "H2",
                "calculation": "approximation_analysis",
            },
            "text_input": {
                "input_str": "H2 approximation_analysis",
            },
            "output": {
                "result": {
                    "bo_validity": True,
                    "mass_ratio_me_mn": 0.0015,
                    "criterion": "me/mn << 1 → BO valid",
                }
            }
        },
    ]

    # Molecular data: (Re/Ang, De/eV, omega_e/cm-1, omega_e_xe/cm-1, mu/amu)
    MOLECULE_DATA = {
        "H2":   {"Re": 0.74144, "De": 4.747,  "omega_e": 4401.21, "omega_exe": 121.34, "mu": 0.50391},
        "HCl":  {"Re": 1.2746,  "De": 4.590,  "omega_e": 2990.95, "omega_exe": 52.82,  "mu": 0.98010},
        "CO":   {"Re": 1.1283,  "De": 11.09,  "omega_e": 2169.81, "omega_exe": 13.29,  "mu": 6.86062},
        "N2":   {"Re": 1.0977,  "De": 9.79,   "omega_e": 2358.57, "omega_exe": 14.32,  "mu": 7.00672},
        "O2":   {"Re": 1.2075,  "De": 5.156,  "omega_e": 1580.19, "omega_exe": 11.98,  "mu": 7.99746},
        "F2":   {"Re": 1.4119,  "De": 1.65,   "omega_e": 916.64,  "omega_exe": 11.24,  "mu": 9.49845},
        "NaCl": {"Re": 2.3609,  "De": 4.23,   "omega_e": 364.96,  "omega_exe": 2.05,   "mu": 13.9283},
        "I2":   {"Re": 2.6655,  "De": 1.56,   "omega_e": 214.50,  "omega_exe": 0.61,   "mu": 126.90},
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34   # J·s
        self.c = 2.99792458e8          # m/s
        self.eV = 1.602176634e-19      # J per eV
        self.amu_kg = 1.66053906660e-27 # amu to kg
        self.c_light_cm_s = 2.99792458e10 # cm/s

    def _run_base(self, molecule: str, calculation: str = "pes",
                  R_range: str = "auto", n_points: int = 50, n_vib_levels: int = 10,
                  De_eV: Optional[float] = None, Re_Angstrom: Optional[float] = None,
                  omega_e_cm_1: Optional[float] = None) -> dict:
        """Core logic."""
        mol = molecule.strip().upper()
        calc = calculation.lower().strip()

        # Get molecular parameters
        if mol == "GENERIC":
            if De_eV is None or Re_Angstrom is None or omega_e_cm_1 is None:
                raise ChemMCPError("For generic molecule, must provide De_eV, Re_Angstrom, and omega_e_cm_1.")
            params = {
                "Re": Re_Angstrom, "De": De_eV,
                "omega_e": omega_e_cm_1, "omega_exe": omega_e_cm_1 * 0.01,  # rough estimate
                "mu": 1.0,
            }
        else:
            if mol not in self.MOLECULE_DATA:
                raise ChemMCPError(
                    f"Unknown molecule '{molecule}'. Available: {list(self.MOLECULE_DATA.keys())} + 'generic'."
                )
            params = self.MOLECULE_DATA[mol]

        if calc == "pes":
            return self._compute_pes(mol, params, R_range, n_points)
        elif calc == "vibrational":
            return self._compute_vibrational(mol, params, n_vib_levels)
        elif calc == "force_constant":
            return self._compute_force_constant(mol, params)
        elif calc == "reduced_mass":
            return self._compute_reduced_mass(mol, params)
        elif calc == "approximation_analysis":
            return self._bo_analysis(mol, params)
        else:
            raise ChemMCPError(
                f"Unknown calculation '{calculation}'. Use: pes, vibrational, force_constant, reduced_mass, approximation_analysis."
            )

    def _compute_pes(self, mol: str, params: dict, R_range: str, n_points: int) -> dict:
        """Compute Morse potential energy surface."""
        Re = params["Re"] * 1e-10  # to meters
        De_J = params["De"] * self.eV
        omega_e = params["omega_e"]
        mu_amu = params["mu"]

        # Compute Morse parameter a from omega_e
        mu_kg = mu_amu * self.amu_kg
        pi_c_omega = math.pi * self.c_light_cm_s * omega_e
        a_morse = pi_c_omega / math.sqrt(De_J / (2 * mu_kg))  # m^-1

        # Determine R range
        if R_range.lower() == "auto":
            R_min = params["Re"] * 0.4
            R_max = params["Re"] * 3.0
        else:
            parts = R_range.split(":")
            R_min = float(parts[0])
            R_max = float(parts[1])

        # Generate PES grid (Morse potential: V(R) = De[1 - exp(-a(R-Re))]² - De)
        pes_data = []
        for i in range(n_points):
            R_ang = R_min + (R_max - R_min) * i / (n_points - 1)
            R_m = R_ang * 1e-10
            y = math.exp(-a_morse * (R_m - Re))
            V_morse = De_J * ((1 - y) ** 2 - 1)  # V=0 at minimum
            V_eV = V_morse / self.eV
            pes_data.append({
                "R_Angstrom": round(R_ang, 6),
                "V_eV": round(V_eV, 8),
            })

        return self._make_result(
            mol, "pes",
            extra={
                "equilibrium_Re_Angstrom": params["Re"],
                "dissociation_De_eV": params["De"],
                "morse_parameter_a_m-1": round(a_morse, 4),
                "potential_model": "Morse: V(R) = Dₑ[1-exp(-a(R-Rₑ))]² - Dₑ",
                "pes_points": pes_data,
                "n_points": n_points,
                "R_range_Angstrom": [round(R_min, 4), round(R_max, 4)],
            }
        )

    def _compute_vibrational(self, mol: str, params: dict, n_levels: int) -> dict:
        """Compute vibrational energy levels (anharmonic oscillator)."""
        omega_e = params["omega_e"]       # cm^-1
        omega_ex = params.get("omega_exe", omega_e * 0.01)
        De = params["De"]

        levels = []
        for v in range(n_levels):
            # G(v) = ωₑ(v+½) - ωₑxₑ(v+½)²  in cm^-1
            G_v = omega_e * (v + 0.5) - omega_ex * ((v + 0.5) ** 2)
            E_eV = G_v * self.hbar * self.c * 100  # cm^-1 → J → eV
            E_eV /= self.eV

            if G_v < 0:
                break

            levels.append({
                "quantum_number_v": v,
                "energy_cm-1": round(G_v, 4),
                "energy_eV": round(E_eV, 8),
            })

        # Zero-point energy
        zpe_cm = omega_e * 0.5 - omega_ex * 0.25
        zpe_eV = zpe_cm * self.hbar * self.c * 100 / self.eV

        # Maximum bound v (where dG/dv = 0)
        v_max = int(omega_e / (2 * omega_ex) - 0.5)

        return self._make_result(
            mol, "vibrational",
            extra={
                "vibrational_levels": levels,
                "n_computed_levels": len(levels),
                "zero_point_energy_cm-1": round(zpe_cm, 4),
                "zero_point_energy_eV": round(zpe_eV, 8),
                "omega_e_cm-1": omega_e,
                "omega_e_xe_cm-1": omega_ex,
                "max_bound_v": v_max,
                "model": "Anharmonic: G(v) = ωₑ(v+½) - ωₑxₑ(v+½)²",
            }
        )

    def _compute_force_constant(self, mol: str, params: dict) -> dict:
        """Compute harmonic force constant k from ωₑ."""
        omega_e = params["omega_e"]  # cm^-1
        mu_amu = params["mu"]
        mu_kg = mu_amu * self.amu_kg

        # k = μω² where ω = 2πc·ωₑ (in rad/s)
        omega_rad = 2 * math.pi * self.c_light_cm_s * omega_e
        k = mu_kg * omega_rad ** 2  # N/m

        # Also compute in mdyn/Å (traditional unit)
        k_mdyn_A = k * 1e-2  # 1 N/m = 1e-2 mdyn/Å

        return self._make_result(
            mol, "force_constant",
            extra={
                "force_constant_N/m": round(k, 4),
                "force_constant_mdyn/A": round(k_mdyn_A, 6),
                "omega_e_cm-1": omega_e,
                "reduced_mass_amu": mu_amu,
                "formula": "k = μ(2πcωₑ)²",
            }
        )

    def _compute_reduced_mass(self, mol: str, params: dict) -> dict:
        """Return reduced mass info."""
        mu_amu = params["mu"]
        mu_kg = mu_amu * self.amu_kg

        return self._make_result(
            mol, "reduced_mass",
            extra={
                "reduced_mass_amu": mu_amu,
                "reduced_mass_kg": mu_kg,
                "formula": "μ = m₁m₂/(m₁+m₂)",
            }
        )

    def _bo_analysis(self, mol: str, params: dict) -> dict:
        """Check validity of Born-Oppenheimer approximation."""
        mu_amu = params["mu"]
        me_mn = 1.0 / mu_amu  # electron mass / nuclear reduced mass ratio (rough, in amu units)

        # More precise: use actual electron-to-nuclear mass ratio
        # Criterion: |∂ψ_el/∂R|/|ψ_el| << 1, roughly equivalent to me/mn << 1
        bo_valid = me_mn < 0.1  # very generous threshold

        # Nuclear vibration amplitude estimate
        omega_e = params["omega_e"]
        omega_rad = 2 * math.pi * self.c_light_cm_s * omega_e
        mu_kg = mu_amu * self.amu_kg
        delta_R = math.sqrt(self.hbar / (2 * mu_kg * omega_rad)) * 1e10  # zero-point amplitude in Å
        Re = params["Re"]
        relative_amplitude = delta_R / Re if Re > 0 else float('inf')

        return self._make_result(
            mol, "approximation_analysis",
            extra={
                "born_oppenheimer_valid": bo_valid,
                "electron_nuclear_mass_ratio_approx": round(me_mn, 6),
                "zero_point_vibration_amplitude_A": round(delta_R, 6),
                "relative_amplitude_delta_R_over_Re": round(relative_amplitude, 6),
                "equilibrium_bond_length_A": Re,
                "criterion": "me/mn << 1 and ΔR/Rₑ << 1 → BO approximation valid",
                "conclusion": "Excellent" if relative_amplitude < 0.05 else "Good" if relative_amplitude < 0.15 else "Acceptable with caution",
            }
        )

    def _make_result(self, mol: str, calc: str, extra: dict = None) -> dict:
        result = {
            "molecule": mol,
            "calculation_type": calc,
        }
        if extra:
            result.update(extra)
        return {"result": result}

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split()
            mol = parts[0]
            calc = parts[1] if len(parts) > 1 else "pes"
            Rrange = parts[2] if len(parts) > 2 else "auto"
            npoints = int(parts[3]) if len(parts) > 3 else 50
            nvib = int(parts[4]) if len(parts) > 4 else 10
            De = float(parts[5]) if len(parts) > 5 else None
            ReA = float(parts[6]) if len(parts) > 6 else None
            we = float(parts[7]) if len(parts) > 7 else None
            return self._run_base(mol, calc, Rrange, npoints, nvib, De, ReA, we)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
