import logging
import math
from typing import Optional, List, Tuple
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HydrogenAtomOrbitals(BaseTool):
    """
    氢原子轨道可视化和能级计算。
    
    能级: E_n = -13.6 eV / n²
    径向波函数: R_nl(r) 使用关联拉盖尔多项式
    完整波函数: ψ_nlm(r,θ,φ) = R_nl(r) · Y_lm(θ,φ)
    """
    __version__ = "0.1.0"
    name = "HydrogenAtomOrbitals"
    func_name = "hydrogen_atom_orbitals"
    description = "Calculate hydrogen atom orbital properties: energy levels, radial wavefunction, probability density, orbital shapes, and node structure."
    implementation_description = "Computes hydrogen-like atomic orbital properties using analytical solutions of the Schrödinger equation. Includes energy eigenvalues (E_n = -13.6Z²/n² eV), radial wavefunctions via associated Laguerre polynomials, probability densities, radial distribution functions, and orbital shape classification."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Hydrogen Atom", "Orbitals", "Energy Levels", "Wavefunction"]
    required_envs = []

    code_input_sig = [
        ("principal_quantum_number_n", "int", "N/A", "Principal quantum number n (n >= 1)."),
        ("orbital_quantum_number_l", "int", "N/A", "Orbital angular momentum quantum number l (0 <= l < n)."),
        ("magnetic_quantum_number_m", "int", "0", "Magnetic quantum number m_l (-l <= m <= l, default=0)."),
        ("position_r_bohr", "float", "1.0", "Radial position in Bohr radii at which to evaluate the wavefunction (default=1.0)."),
        ("nuclear_charge_Z", "int", "1", "Nuclear charge Z for hydrogen-like ions (default=1)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'n l [m] [r_bohr] [Z]'. Example: '2 1 0 1.0 1' for 2p orbital at r=1a₀."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with energy, wavefunction values, probability density, node counts, radius info, degeneracy, and orbital classification."),
    ]

    examples = [
        {
            "code_input": {
                "principal_quantum_number_n": 1,
                "orbital_quantum_number_l": 0,
                "magnetic_quantum_number_m": 0,
                "position_r_bohr": 1.0,
                "nuclear_charge_Z": 1,
            },
            "text_input": {
                "input_params": "1 0 0 1.0 1",
            },
            "output": {
                "result": {
                    "energy_eV": -13.6057,
                    "orbital_type": "1s",
                    "n_radial_nodes": 0,
                    "total_nodes": 0,
                }
            },
        },
        {
            "code_input": {
                "principal_quantum_number_n": 2,
                "orbital_quantum_number_l": 1,
                "magnetic_quantum_number_m": 0,
                "position_r_bohr": 2.0,
                "nuclear_charge_Z": 1,
            },
            "text_input": {
                "input_params": "2 1 0 2.0 1",
            },
            "output": {
                "result": {
                    "energy_eV": -3.4014,
                    "orbital_type": "2p",
                    "n_radial_nodes": 0,
                    "angular_nodes": 1,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.a0_Bohr = 5.29177210903e-11  # Bohr radius in meters
        self.Ry_eV = 13.605693122994  # Rydberg constant in eV
        self.eV_per_J = 6.241509e18

    def _factorial(self, n: int) -> int:
        if n <= 1:
            return 1
        result = 1
        for i in range(2, n + 1):
            result *= i
        return result

    def _double_factorial(self, n: int) -> int:
        """Compute double factorial n!!."""
        if n <= 0:
            return 1
        result = 1
        while n > 0:
            result *= n
            n -= 2
        return result

    def _associated_laguerre(self, n: int, k: int, x: float) -> float:
        """Compute associated Laguerre polynomial L_n^k(x) using recursion.
        
        L_0^k(x) = 1
        L_1^k(x) = k + 1 - x
        (n+1)L_{n+1}^k(x) = (2n+k+1-x)L_n^k(x) - (n+k)L_{n-1}^k(x)
        """
        if n == 0:
            return 1.0
        elif n == 1:
            return float(k + 1 - x)

        L_prev2 = 1.0  # L_0^k
        L_prev1 = float(k + 1 - x)  # L_1^k

        for i in range(2, n + 1):
            L_curr = ((2 * i + k - 1 - x) * L_prev1 - (i + k - 1) * L_prev2) / i
            L_prev2 = L_prev1
            L_prev1 = L_curr

        return L_prev1

    def _spherical_harmonic_Y00(self, theta: float, phi: float) -> complex:
        """Y_0^0 = 1/sqrt(4π)"""
        return 1.0 / math.sqrt(4.0 * math.pi)

    def _spherical_harmonic_Y10(self, theta: float, phi: float) -> complex:
        """Y_1^0 = sqrt(3/(4π))·cos(θ)"""
        return math.sqrt(3.0 / (4.0 * math.pi)) * math.cos(theta)

    def _spherical_harmonic_Y1pm1(self, theta: float, phi: float, m: int) -> complex:
        """Y_1^{±1} = ∓sqrt(3/(8π))·sin(θ)·e^{±iφ}"""
        sign = -1.0 if m == 1 else 1.0
        amp = math.sqrt(3.0 / (8.0 * math.pi)) * math.sin(theta)
        if m == 1:
            return complex(-amp * math.cos(phi), -amp * math.sin(phi))
        else:
            return complex(amp * math.cos(phi), -amp * math.sin(phi))

    def _spherical_harmonic_Y20(self, theta: float, phi: float) -> complex:
        """Y_2^0 = sqrt(5/(16π))·(3cos²θ-1)"""
        return math.sqrt(5.0 / (16.0 * math.pi)) * (3.0 * math.cos(theta)**2 - 1.0)

    def _spherical_harmonic_Y2pm1(self, theta: float, phi: float, m: int) -> complex:
        """Y_2^{±1} = ∓sqrt(15/(8π))·sinθcosθ·e^{±iφ}"""
        sign = -1.0 if m == 1 else 1.0
        amp = math.sqrt(15.0 / (8.0 * math.pi)) * math.sin(theta) * math.cos(theta)
        if m == 1:
            return complex(-amp * math.cos(phi), -amp * math.sin(phi))
        else:
            return complex(amp * math.cos(phi), -amp * math.sin(phi))

    def _spherical_harmonic_Y2pm2(self, theta: float, phi: float, m: int) -> complex:
        """Y_2^{±2} = sqrt(15/(32π))·sin²θ·e^{±2iφ}"""
        amp = math.sqrt(15.0 / (32.0 * math.pi)) * math.sin(theta)**2
        angle = 2.0 * phi if m == 2 else -2.0 * phi
        return complex(amp * math.cos(angle), amp * math.sin(angle))

    def _get_spherical_harmonic(self, l: int, m: int, theta: float, phi: float) -> complex:
        """Get spherical harmonic Y_l^m(θ, φ)."""
        if l == 0 and m == 0:
            return self._spherical_harmonic_Y00(theta, phi)
        elif l == 1:
            if m == 0:
                return self._spherical_harmonic_Y10(theta, phi)
            elif abs(m) == 1:
                return self._spherical_harmonic_Y1pm1(theta, phi, m)
        elif l == 2:
            if m == 0:
                return self._spherical_harmonic_Y20(theta, phi)
            elif abs(m) == 1:
                return self._spherical_harmonic_Y2pm1(theta, phi, m)
            elif abs(m) == 2:
                return self._spherical_harmonic_Y2pm2(theta, phi, m)
        
        # For higher l, return simplified approximation
        raise ChemMCPError(f"Spherical harmonic Y_{l}^{m} not implemented for l > 2.")

    def _radial_wavefunction(self, n: int, l: int, rho: float, Z: int) -> float:
        """
        Compute radial wavefunction R_nl(r) where ρ = 2Zr/(na₀).
        
        R_nl(r) = sqrt((2Z/na₀)³ · (n-l-1)!/(2n[(n+l)!]³)) · ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2)
        """
        rho_val = rho
        
        # Prefactor computation
        nf = self._factorial(n - l - 1)
        nlf = self._factorial(n + l)
        prefactor = math.sqrt(
            (2.0 * Z / (n ** 3)) ** 3 * nf / (2.0 * n * (nlf ** 3))
        )
        
        # ρ^l term
        rho_power = rho_val ** l
        
        # Associated Laguerre polynomial
        L = self._associated_laguerre(n - l - 1, 2 * l + 1, rho_val)
        
        # Exponential decay
        exp_term = math.exp(-rho_val / 2.0)
        
        return prefactor * rho_power * L * exp_term

    def _most_probable_radius(self, n: int, l: int, Z: int) -> float:
        """Most probable radius r_mp for hydrogenic orbital (in Bohr radii).
        
        For ns orbitals: r_mp = n²/Z a₀
        General: solve d/dr[r²R²] = 0 → approximately n²/Z for s orbitals
        """
        # Approximate: most probable radius ≈ n²/Z for l=0, increases with l
        # More precisely from dP/dr = 0 where P(r) = r²|R_nl|²
        if l == 0:
            return n * n / Z
        elif l == n - 1:
            return n * (n + 0.5) / Z  # Circular orbits
        else:
            # Interpolation between extremes
            frac = l / (n - 1) if n > 1 else 0
            return (n * n / Z) * (1.0 + frac * 0.5)

    def _mean_radius(self, n: int, l: int, Z: int) -> float:
        """Expectation value <r> = (a₀/2Z)[3n² - l(l+1)] in units of a₀."""
        return (3.0 * n * n - l * (l + 1)) / (2.0 * Z)

    def _mean_radius_squared(self, n: int, l: int, Z: int) -> float:
        """Expectation value <r²> = (n²a₀²/2Z²)[5n²+1-3l(l+1)] in units of a₀²."""
        return (n * n / (2.0 * Z * Z)) * (5.0 * n * n + 1 - 3.0 * l * (l + 1))

    def _run_base(self, principal_quantum_number_n: int, orbital_quantum_number_l: int,
                  magnetic_quantum_number_m: int = 0, position_r_bohr: float = 1.0,
                  nuclear_charge_Z: int = 1) -> dict:

        n = principal_quantum_number_n
        l = orbital_quantum_number_l
        m = magnetic_quantum_number_m
        Z = nuclear_charge_Z

        # Validation
        if n < 1:
            raise ChemMCPError("Principal quantum number n must be >= 1.")
        if l < 0 or l >= n:
            raise ChemMCPError(f"Orbital quantum number l must satisfy 0 <= l < n (got l={l}, n={n}).")
        if abs(m) > l:
            raise ChemMCPError(f"Magnetic quantum number |m| must be <= l (got m={m}, l={l}).")
        if position_r_bohr < 0:
            raise ChemMCPError("Radial position must be >= 0.")
        if Z < 1:
            raise ChemMCPError("Nuclear charge Z must be >= 1.")

        # Orbital type label
        orbital_labels = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}
        orbital_label = orbital_labels.get(l, f"l={l}")
        orbital_type = f"{n}{orbital_label}"

        # Energy: E_n = -Ry · Z²/n² (in eV)
        energy_eV = -self.Ry_eV * Z * Z / (n * n)
        energy_J = energy_eV / self.eV_per_J

        # Radial coordinate conversion: ρ = 2Zr/(na₀), where r is in units of a₀
        rho = 2.0 * Z * position_r_bohr / n

        # Radial wavefunction value
        R_val = self._radial_wavefunction(n, l, rho, Z)

        # Angular part (evaluate at θ=π/2, φ=0 as default "equatorial" direction)
        theta = math.pi / 2.0
        phi = 0.0
        try:
            Y_val = self._get_spherical_harmonic(l, m, theta, phi)
        except ChemMCPError:
            Y_val = complex(1.0, 0.0)

        # Full wavefunction: ψ = R(r) · Y(θ,φ)
        psi_val = R_val * Y_val
        prob_density = abs(psi_val) ** 2

        # Radial distribution function: D(r) = r² · |R(r)|²
        radial_dist = position_r_bohr * position_r_bohr * R_val * R_val

        # Node counts
        n_radial_nodes = n - l - 1
        n_angular_nodes = l
        total_nodes = n - 1

        # Degeneracy (ignoring fine structure): g_n = n²
        degeneracy = n * n

        # Radii
        most_probable_r = self._most_probable_radius(n, l, Z)
        mean_r = self._mean_radius(n, l, Z)
        mean_r2 = self._mean_radius_squared(n, l, Z)
        delta_r = math.sqrt(mean_r2 - mean_r * mean_r)

        # Bohr radius for this state: a₀_eff = a₀ · n / Z
        bohr_effective = n / Z

        result = {
            "quantum_numbers": {"n": n, "l": l, "m": m},
            "orbital_type": orbital_type,
            "nuclear_charge_Z": Z,
            "energy_eV": round(energy_eV, 10),
            "energy_J": round(energy_J, 30),
            "radial_position_bohr": position_r_bohr,
            "rho_parameter": round(rho, 10),
            "radial_wavefunction_R": round(R_val, 15),
            "wavefunction_psi_real": round(psi_val.real, 15),
            "wavefunction_psi_imag": round(psi_val.imag, 15),
            "probability_density": round(prob_density, 15),
            "radial_distribution_D_r": round(radial_dist, 15),
            "n_radial_nodes": n_radial_nodes,
            "angular_nodes": n_angular_nodes,
            "total_nodes": total_nodes,
            "degeneracy": degeneracy,
            "most_probable_radius_bohr": round(most_probable_r, 6),
            "mean_radius_bohr": round(mean_r, 6),
            "mean_radius_sq_bohr2": round(mean_r2, 6),
            "uncertainty_delta_r_bohr": round(delta_r, 6),
            "effective_bohr_radius_n_by_Z": round(bohr_effective, 6),
            "ionization_energy_eV": round(abs(energy_eV), 10),
        }

        logger.info(f"HydrogenAtomOrbitals: {orbital_type}, E={energy_eV:.6f}eV, r_mp={most_probable_r:.3f}a₀")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            n = int(parts[0])
            l = int(parts[1])
            m = int(parts[2]) if len(parts) > 2 else 0
            r = float(parts[3]) if len(parts) > 3 else 1.0
            Z = int(parts[4]) if len(parts) > 4 else 1
            return self._run_base(n, l, m, r, Z)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'n l [m] [r_bohr] [Z]'")
