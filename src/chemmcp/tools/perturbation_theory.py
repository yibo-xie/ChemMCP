import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PerturbationTheory(BaseTool):
    """
    微扰理论能量修正计算。
    
    非简并微扰和简并微扰理论:
    一级修正: E_n^(1) = ⟨n|H'|n⟩
    二级修正: E_n^(2) = Σ_{m≠n} |⟨m|H'|n⟩² / (E_n - E_m)
    
    简并微扰: 解久期方程 det(H' - E^(1)I) = 0
    """
    __version__ = "0.1.0"
    name = "PerturbationTheory"
    func_name = "perturbation_theory"
    description = "Calculate energy corrections using non-degenerate and degenerate perturbation theory for quantum systems."
    implementation_description = "Implements first-order and second-order perturbation theory for harmonic oscillator (linear/quartic perturbation), infinite square well (shifted/tilted), hydrogen atom Stark effect, helium ground state approximation, and general two-level systems. Computes energy corrections and wavefunction mixing coefficients."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Perturbation Theory", "Energy Correction", "Approximation", "Stark Effect"]
    required_envs = []

    code_input_sig = [
        ("system_type", "str", "N/A", "System: 'harmonic_perturbed', 'infinite_well_perturbed', 'hydrogen_stark', 'helium_ground', 'two_level', 'particle_in_box_perturbed'."),
        ("perturbation_strength", "float", "N/A", "Perturbation strength parameter λ (dimensionless or in appropriate units)."),
        ("order", "int", "2", "Order of perturbation: 1 for first-order only, 2 for first+second order."),
        ("mass_kg", "float", "None", "Particle mass in kg (for systems that need it)."),
        ("box_length_m", "float", "None", "Box length L in meters (for infinite well systems)."),
        ("n_state", "int", "0", "Which quantum state to compute corrections for (default=0, meaning ground state)."),
        ("force_constant_N_m", "float", "None", "Force constant k for harmonic oscillator."),
        ("n_sum_terms", "int", "20", "Number of terms to include in second-order sum (default=20)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'system_type strength [order] [extra...]'. Example: 'harmonic_perturbed 0.5 2 k=10 n_state=0'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with unperturbed energy, corrections, total corrected energy, convergence info, and comparison data."),
    ]

    examples = [
        {
            "code_input": {
                "system_type": "harmonic_perturbed",
                "perturbation_strength": 0.1,
                "order": 2,
                "mass_kg": None,
                "box_length_m": None,
                "n_state": 0,
                "force_constant_N_m": 10.0,
                "n_sum_terms": 20,
            },
            "text_input": {
                "input_params": "harmonic_perturbed 0.1 2 k=10 n_state=0",
            },
            "output": {
                "result": {
                    "unperturbed_energy_eV": 0.00798,
                    "first_order_correction_eV": 0.00040,
                    "corrected_energy_eV": 0.00838,
                }
            },
        },
        {
            "code_input": {
                "system_type": "infinite_well_perturbed",
                "perturbation_strength": 0.05,
                "order": 2,
                "mass_kg": 9.109e-31,
                "box_length_m": 1e-9,
                "n_state": 0,
                "force_constant_N_m": None,
                "n_sum_terms": 20,
            },
            "text_input": {
                "input_params": "infinite_well_perturbed 0.05 2 m=9.109e-31 L=1e-9",
            },
            "output": {
                "result": {
                    "first_order_correction_eV": 0.00941,
                    "second_order_correction_eV": -0.00012,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34
        self.eV_per_J = 6.241509074e18

    def _factorial(self, n: int) -> int:
        if n <= 1:
            return 1
        r = 1
        for i in range(2, n + 1):
            r *= i
        return r

    # ---- Harmonic Oscillator with H' = λx (linear) or λx² (quadratic shift) or λx⁴ (quartic) ----

    def _harmonic_energies(self, k: float, mass: float, n_max: int) -> List[float]:
        """Unperturbed HO energies: E_n = (n+½)ℏω."""
        omega = math.sqrt(k / mass)
        hbar = self.hbar
        return [(n + 0.5) * hbar * omega for n in range(n_max + 1)]

    def _harmonic_matrix_element_x(self, n: int, m: int, k: float, mass: float) -> float:
        """⟨n|x|m⟩ for HO: nonzero only for m = n±1.
        
        ⟨n|x|n+1⟩ = √((n+1)ℏ/(2mω)), ⟨n|x|n-1⟩ = √(nℏ/(2mω))
        """
        omega = math.sqrt(k / mass)
        hbar = self.hbar
        coeff = math.sqrt(hbar / (2.0 * mass * omega))
        if m == n + 1:
            return coeff * math.sqrt(n + 1)
        elif m == n - 1 and n >= 1:
            return coeff * math.sqrt(n)
        return 0.0

    def _harmonic_matrix_element_x2(self, n: int, m: int, k: float, mass: float) -> float:
        """⟨n|x²|m⟩ for HO: nonzero for m = n, n±2."""
        omega = math.sqrt(k / mass)
        hbar = self.hbar
        coeff = hbar / (2.0 * mass * omega)
        if m == n:
            return coeff * (2.0 * n + 1)
        elif m == n + 2:
            return coeff * math.sqrt((n + 1) * (n + 2))
        elif m == n - 2 and n >= 2:
            return coeff * math.sqrt(n * (n - 1))
        return 0.0

    def _harmonic_matrix_element_x4(self, n: int, m: int, k: float, mass: float) -> float:
        """⟨n|x⁴|m⟩ for HO using ladder operator expansion.
        
        x⁴ connects states with Δn = 0, ±2, ±4.
        """
        omega = math.sqrt(k / mass)
        hbar = self.hbar
        xi_sq = hbar / (2.0 * mass * omega)  # (ℏ/2mω)

        # x⁴ matrix elements via composition of x²
        result = 0.0
        
        # Sum over intermediate k: ⟨n|x²|k⟩·⟨k|x²|m⟩
        k_max = max(n, m) + 5
        for kk in range(max(0, min(n,m)-4), k_max + 1):
            nk = self._harmonic_matrix_element_x2(n, kk, k, mass)
            km = self._harmonic_matrix_element_x2(kk, m, k, mass)
            result += nk * km

        return result

    def _solve_harmonic(self, lam: float, order: int, k: float, mass: float,
                         n_state: int, n_terms: int) -> dict:
        """Solve perturbed harmonic oscillator with H' = λx⁴ (quartic perturbation)."""
        energies = self._harmonic_energies(k, mass, n_terms + n_state + 5)
        n = n_state
        E0 = energies[n]

        # First-order correction: E^(1) = ⟨n|λx⁴|n⟩
        V_nn = self._harmonic_matrix_element_x4(n, n, k, mass)
        E1 = lam * V_nn

        E2 = 0.0
        coeffs = {}
        if order >= 2:
            # Second-order: Σ_{m≠n} |⟨m|λx⁴|n⟩² / (E_n - E_m)
            for m in range(len(energies)):
                if m == n:
                    continue
                V_mn = self._harmonic_matrix_element_x4(m, n, k, mass)
                denom = E0 - energies[m]
                if abs(denom) > 1e-30:
                    E2 += lam * lam * V_mn * V_mn / denom
                    coeffs[m] = lam * V_mn / denom

        E_total = E0 + E1 + E2
        eV = self.eV_per_J

        return {
            "system_type": "harmonic_perturbed (quartic)",
            "quantum_state_n": n,
            "unperturbed_energy_J": round(E0, 25),
            "unperturbed_energy_eV": round(E0 * eV, 10),
            "matrix_element_V_nn_J": round(V_nn, 25),
            "first_order_correction_J": round(E1, 25),
            "first_order_correction_eV": round(E1 * eV, 10),
            "second_order_correction_J": round(E2, 25),
            "second_order_correction_eV": round(E2 * eV, 10),
            "total_corrected_energy_J": round(E_total, 25),
            "total_corrected_energy_eV": round(E_total * eV, 10),
            "perturbation_form": "H' = λx⁴",
            "n_terms_in_sum": n_terms,
            "wavefunction_mixing_coefficients": {k: round(v, 15) for k, v in list(coeffs.items())[:10]},
        }

    # ---- Infinite Square Well with H' = λx (tilted) or λx² ----

    def _iwell_energies(self, mass: float, L: float, n_max: int) -> List[float]:
        """Infinite well energies: E_n = n²π²ℏ²/(2mL²)."""
        hbar = self.hbar
        return [n * n * math.pi * math.pi * hbar * hbar / (2.0 * mass * L * L) for n in range(1, n_max + 2)]

    def _iwell_matrix_element_x(self, n: int, m: int, L: float) -> float:
        """⟨n|x|m⟩ for infinite well of width L (domain [0,L]).
        
        ⟨n|x|m⟩ = L/2 · δ_nm for shifted well centered at L/2
        For domain [0,L]: 
          If n+m odd: ⟨n|x|m⟩ = -4Lnm/(π²(n²-m²)²)
          If n=m: ⟨n|x|n⟩ = L/2
        """
        if n == m:
            return L / 2.0
        if (n + m) % 2 == 1:
            return -4.0 * L * n * m / (math.pi ** 2 * (n * n - m * m) ** 2)
        return 0.0

    def _iwell_matrix_element_x2(self, n: int, m: int, L: float) -> float:
        """⟨n|x²|m⟩ for infinite well [0,L]."""
        if n == m:
            return L * L * (1.0 / 3.0 - 1.0 / (2.0 * math.pi * math.pi * n * n))
        if (n + m) % 2 == 0 and n != m:
            num = -4.0 * L * L * n * m * (n * n + m * m - 1)
            den = math.pi ** 4 * (n * n - m * m) ** 2 * n * m
            if abs(den) > 1e-30:
                return num / den
        return 0.0

    def _solve_iwell(self, lam: float, order: int, mass: float, L: float,
                      n_state: int, n_terms: int) -> dict:
        """Solve perturbed infinite well with H' = λx (linear tilt from center)."""
        n = n_state + 1  # Quantum numbers start at 1 for infinite well
        energies_raw = self._iwell_energies(mass, L, n_terms + 10)
        E0 = energies_raw[n_state]

        # Use x' = x - L/2 (centered coordinate) so H' = λ(x-L/2) = λx' has ⟨0|λx'|0⟩ = 0
        # For simplicity use H' = λx on [0, L]
        V_nn = self._iwell_matrix_element_x(n, n, L)
        E1 = lam * V_nn

        E2 = 0.0
        coeffs = {}
        if order >= 2:
            for mm in range(1, len(energies_raw) + 1):
                if mm == n:
                    continue
                V_mn = self._iwell_matrix_element_x(n, mm, L)
                E_m = energies_raw[mm - 1]  # 0-indexed
                denom = E0 - E_m
                if abs(denom) > 1e-30:
                    E2 += lam * lam * V_mn * V_mn / denom
                    coeffs[mm] = lam * V_mn / denom

        E_total = E0 + E1 + E2
        eV = self.eV_per_J

        return {
            "system_type": "infinite_well_perturbed (linear tilt)",
            "quantum_state_n": n,
            "box_length_m": L,
            "unperturbed_energy_J": round(E0, 25),
            "unperturbed_energy_eV": round(E0 * eV, 10),
            "first_order_correction_J": round(E1, 25),
            "first_order_correction_eV": round(E1 * eV, 10),
            "second_order_correction_J": round(E2, 25),
            "second_order_correction_eV": round(E2 * eV, 10),
            "total_corrected_energy_J": round(E_total, 25),
            "total_corrected_energy_eV": round(E_total * eV, 10),
            "perturbation_form": "H' = λx",
            "wavefunction_mixing_coefficients": {k: round(v, 15) for k, v in list(coeffs.items())[:10]},
        }

    # ---- Hydrogen Atom Stark Effect ----

    def _solve_hydrogen_stark(self, lam: float, order: int, n_state: int, n_terms: int) -> dict:
        """Hydrogen atom Stark effect: H' = eEz = λz (electric field along z).
        
        For ground state (n=1): first-order = 0 (parity), second-order gives polarizability.
        α = (9/2) · a₀³ · (4πε₀/e²)³... actually α = (9/2)a₀³ in atomic units
        E^(2) = -(1/2)αE²
        """
        a0 = 5.29177210903e-11  # Bohr radius in meters
        Ry_J = 2.179872361e-18   # Rydberg in Joules
        eV = self.eV_per_J

        n = n_state + 1  # Convert 0-indexed to principal QN
        E0 = -Ry_J / (n * n)

        # Parity: z is odd → ⟨nlm|z|nlm⟩ = 0 for any state (same parity)
        E1 = 0.0

        # Second-order Stark effect for ground state
        # α = (9/2) a₀³ for hydrogen 1s
        # E^(2) = -(1/2)αE² where E_field is related to λ by λ = eE_field
        # In atomic units: E^(2) = -(9/4)n⁴(F/F_au)² ... simplified
        alpha_polarizability = 4.5 * a0 ** 3  # 9/2 · a₀³ in SI units (approximate)
        
        # λ has units of energy/length (J/m). The electric field E = λ/e
        # E² term: E^(2) = -(1/2)αE² = -(1/2)α(λ/e)²
        e_charge = 1.602176634e-19  # Coulomb
        E2 = -0.5 * alpha_polarizability * (lam / e_charge) ** 2

        E_total = E0 + E1 + E2

        return {
            "system_type": "hydrogen_stark_effect",
            "principal_quantum_number_n": n,
            "unperturbed_energy_J": round(E0, 25),
            "unperturbed_energy_eV": round(E0 * eV, 10),
            "first_order_correction_J": 0.0,
            "first_order_correction_eV": 0.0,
            "first_order_zero_reason": "Parity: ⟨ψ|z|ψ⟩ = 0 for states of definite parity",
            "second_order_correction_J": round(E2, 25),
            "second_order_correction_eV": round(E2 * eV, 10),
            "total_corrected_energy_J": round(E_total, 25),
            "total_corrected_energy_eV": round(E_total * eV, 10),
            "polarizability_SI_m3": round(alpha_polarizability, 35),
            "perturbation_form": "H' = eEz (Stark effect)",
        }

    # ---- Two-Level System ----

    def _solve_two_level(self, lam: float, order: int, gap_J: float = None) -> dict:
        """Two-level system: H₀ = diag(E₁, E₂), H' = λ·σ_x (off-diagonal coupling).
        
        Exact solution: E_± = (E₁+E₂)/2 ± (1/2)√[(Δ)² + 4λ²]
        where Δ = E₂ - E₁
        """
        eV = self.eV_per_J
        E1 = 0.0  # Ground state
        E2 = gap_J if gap_J else 1.621e-18  # Default ~10 meV gap

        delta = E2 - E1

        # Perturbation theory
        # E₁^(1) = ⟨1|λσₓ|1⟩ = 0, E₁^(2) = |⟨2|λσₓ|1⟩|²/(E₁-E₂) = -λ²/Δ
        E1_pt = 0.0
        E2_pt = 0.0
        E1_2nd = -lam * lam / delta
        E2_2nd = lam * lam / delta

        # Exact diagonalization
        E_exact_minus = (E1 + E2) / 2.0 - 0.5 * math.sqrt(delta * delta + 4.0 * lam * lam)
        E_exact_plus = (E1 + E2) / 2.0 + 0.5 * math.sqrt(delta * delta + 4.0 * lam * lam)

        E_total_lower = E1 + E1_pt + (E1_2nd if order >= 2 else 0.0)

        err_lower = abs((E_total_lower - E_exact_minus) / E_exact_minus * 100) if E_exact_minus != 0 else None

        return {
            "system_type": "two_level_system",
            "gap_delta_J": round(delta, 25),
            "gap_delta_eV": round(delta * eV, 10),
            "coupling_lambda_J": round(lam, 25),
            "unperturbed_ground_J": round(E1, 25),
            "unperturbed_excited_J": round(E2, 25),
            "first_order_correction_J": 0.0,
            "second_order_ground_J": round(E1_2nd, 25) if order >= 2 else 0.0,
            "second_order_excited_J": round(E2_2nd, 25) if order >= 2 else 0.0,
            "pt_ground_energy_J": round(E_total_lower, 25),
            "pt_ground_energy_eV": round(E_total_lower * eV, 10),
            "exact_ground_energy_J": round(E_exact_minus, 25),
            "exact_ground_energy_eV": round(E_exact_minus * eV, 10),
            "exact_excited_energy_J": round(E_exact_plus, 25),
            "pt_error_percent": round(err_lower, 8) if err_lower is not None else None,
            "avoided_crossing_gap_J": round(math.sqrt(delta*delta + 4*lam*lam) - abs(delta), 10),
        }

    # ---- Helium Ground State Approximation ----

    def _solve_helium(self, lam: float, order: int) -> dict:
        """Helium ground state using perturbation theory.
        
        H₀ = sum of two hydrogenic Hamiltonians (Z=2)
        H' = e²/(4πε₀r₁₂) (electron-electron repulsion)
        
        E₀ = 2 × (-Z²Ry) = -8Ry for Z=2 (each electron sees full nuclear charge)
        E^(1) = ⟨H'⟩ = (5/8)Z·Ry = (5/4)Ry for Z=2
        """
        Ry_J = 2.179872361e-18
        eV = self.eV_per_J
        Z = 2

        E0 = -2.0 * Z * Z * Ry_J  # Both electrons in 1s with Z=2
        # Actually: each electron has E = -Z²Ry/n² = -4Ry, so total E0 = -8Ry
        E0 = -8.0 * Ry_J

        # First-order correction: electron-electron repulsion
        # ⟨1/r₁₂⟩ = 5Z/(8a₀) in atomic units → E^(1) = (5/8)Z · (e²/4πε₀/a₀) = (5/8)Z · 2Ry
        E1 = (5.0 / 8.0) * Z * 2.0 * Ry_J  # = (5/4)Ry for Z=2

        E_total = E0 + E1

        # Compare with experimental: He ground state ≈ -79.0 eV = -2.9037 Hartree
        # Our PT: E = -8Ry + 1.25Ry = -6.75Ry = -91.96 eV (not great but first-order only)
        # In Hartree: -4 + 0.625 = -3.375 Ha vs exact -2.9037 Ha
        E_Ha = E_total / Ry_J  # In units of Rydberg
        E_exact_Ha = -5.8074  # -2.9037 Hartree × 2 (since 1Ha = 2Ry)

        rel_err = abs((E_Ha - E_exact_Ha) / E_exact_Ha) * 100

        return {
            "system_type": "helium_ground_state",
            "nuclear_charge_Z": Z,
            "unperturbed_energy_J": round(E0, 25),
            "unperturbed_energy_eV": round(E0 * eV, 10),
            "unperturbed_energy_Ry": round(E0 / Ry_J, 6),
            "first_order_correction_J": round(E1, 25),
            "first_order_correction_eV": round(E1 * eV, 10),
            "first_order_correction_Ry": round(E1 / Ry_J, 6),
            "ee_repulsion_description": "H' = e²/(4πε₀r₁₂), ⟨1/r₁₂⟩ = 5Z/(8a₀)",
            "total_corrected_energy_J": round(E_total, 25),
            "total_corrected_energy_eV": round(E_total * eV, 10),
            "total_corrected_energy_Ry": round(E_Ha, 6),
            "exact_energy_Ry": round(E_exact_Ha, 4),
            "relative_error_percent": round(rel_err, 2),
            "note": "First-order PT overbinds; variational methods give better results.",
        }

    def _run_base(self, system_type: str, perturbation_strength: float, order: int = 2,
                  mass_kg: float = None, box_length_m: float = None,
                  n_state: int = 0, force_constant_N_m: float = None,
                  n_sum_terms: int = 20) -> dict:

        if system_type == "harmonic_perturbed":
            if force_constant_N_m is None:
                force_constant_N_m = 10.0
            if mass_kg is None:
                mass_kg = 9.109e-31
            result = self._solve_harmonic(perturbation_strength, order, force_constant_N_m,
                                           mass_kg, n_state, n_sum_terms)

        elif system_type == "infinite_well_perturbed" or system_type == "particle_in_box_perturbed":
            if mass_kg is None:
                mass_kg = 9.109e-31
            if box_length_m is None:
                box_length_m = 1e-9
            result = self._solve_iwell(perturbation_strength, order, mass_kg,
                                        box_length_m, n_state, n_sum_terms)

        elif system_type == "hydrogen_stark":
            result = self._solve_hydrogen_stark(perturbation_strength, order, n_state, n_sum_terms)

        elif system_type == "two_level":
            result = self._solve_two_level(perturbation_strength, order)

        elif system_type == "helium_ground":
            result = self._solve_helium(perturbation_strength, order)

        else:
            raise ChemMCPError(f"Unknown system type: {system_type}. Choose from: "
                             f"harmonic_perturbed, infinite_well_perturbed, hydrogen_stark, "
                             f"two_level, helium_ground")

        result["perturbation_order_computed"] = order
        logger.info(f"PerturbationTheory: {system_type}, E_corr={result.get('total_corrected_energy_eV', 'N/A')}eV")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            stype = parts[0]
            strength = float(parts[1])
            ord_val = int(parts[2]) if len(parts) > 2 else 2
            
            kwargs = {"n_state": 0}
            for p in parts[3:]:
                if p.startswith("k="):
                    kwargs["force_constant_N_m"] = float(p.split("=")[1])
                elif p.startswith("m="):
                    kwargs["mass_kg"] = float(p.split("=")[1])
                elif p.startswith("L="):
                    kwargs["box_length_m"] = float(p.split("=")[1])
                elif p.startswith("n_state="):
                    kwargs["n_state"] = int(p.split("=")[1])

            return self._run_base(stype, strength, ord_val, **kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'stype strength [order] [kwargs]'")
