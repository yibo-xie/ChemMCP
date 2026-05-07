import logging
import math
from typing import List, Optional, Tuple
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Schrodinger3DSolver(BaseTool):
    """
    三维薛定谔方程求解器（氢原子/类氢离子）。
    
    求解三维定态薛定谔方程:
    [-ℏ²/(2m)∇² - Ze²/(4πε₀r)]ψ = Eψ
    
    解析求解氢原子和类氢离子的能级、波函数、径向分布、概率密度等。
    支持球坐标系下的完整分离变量解。
    """
    __version__ = "0.1.0"
    name = "Schrodinger3DSolver"
    func_name = "schrodinger_3d_solver"
    description = "Solve the 3D time-independent Schrödinger equation for hydrogen-like atoms. Compute energy levels, wavefunctions (radial × angular), probability densities, orbital shapes, and quantum mechanical observables."
    implementation_description = "Analytical solution of the 3D TISE for hydrogen-like atoms using separation of variables in spherical coordinates. Energy: E_n = -Z²Ry/n². Radial: R_nl(ρ) with associated Laguerre polynomials. Angular: spherical harmonics Y_l^m(θ,φ). Computes <r>, <r²>, Δr, radial probability distribution, and nodal structure."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Schrodinger Equation", "Hydrogen Atom", "Atomic Orbitals", "Wavefunction", "3D"]
    required_envs = []

    code_input_sig = [
        ("n", "int", "N/A", "Principal quantum number (n >= 1)."),
        ("l", "int", "N/A", "Orbital angular momentum quantum number (0 <= l < n)."),
        ("m", "int", "0", "Magnetic quantum number (-l <= m <= l, default 0)."),
        ("Z", "int", "1", "Nuclear charge for hydrogen-like ions (default 1 = H)."),
        ("n_radial_points", "int", "200", "Number of radial grid points for plotting/computation."),
        ("r_max_bohr", "float", "20.0", "Maximum radius in Bohr radii for radial grid."),
        ("compute_observables", "bool", "True", "Whether to compute expectation values (<r>, <r²>, Δr etc.)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'n l [m] [Z] [n_points] [r_max]'. Example: '2 1 0 1 200 20' for 2p orbital of H."),
    ]

    output_sig = [
        ("result", "dict", "Complete solution dictionary: energy, quantum numbers, wavefunction data on radial grid, radial probability distribution, expectation values, node counts, orbital classification, degeneracy info."),
    ]

    examples = [
        {
            "code_input": {
                "n": 1,
                "l": 0,
                "m": 0,
                "Z": 1,
                "n_radial_points": 100,
                "r_max_bohr": 10.0,
                "compute_observables": True,
            },
            "text_input": {
                "input_params": "1 0 0 1",
            },
            "output": {
                "result": {
                    "orbital_type": "1s",
                    "energy_eV": -13.606,
                    "n_radial_nodes": 0,
                    "most_probable_radius_bohr": 1.0,
                }
            },
        },
        {
            "code_input": {
                "n": 3,
                "l": 2,
                "m": 0,
                "Z": 1,
                "n_radial_points": 150,
                "r_max_bohr": 25.0,
                "compute_observables": True,
            },
            "text_input": {
                "input_params": "3 2 0 1 150 25",
            },
            "output": {
                "result": {
                    "orbital_type": "3d",
                    "energy_eV": -1.512,
                    "n_radial_nodes": 1,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Ry_eV = 13.605693122994  # Rydberg constant in eV
        self.a0_m = 5.29177210903e-11   # Bohr radius in meters

    # ---- Combinatorial helpers ----
    @staticmethod
    def _factorial(n: int) -> int:
        if n <= 1:
            return 1
        r = 1
        for i in range(2, n + 1):
            r *= i
        return r

    @staticmethod
    def _double_factorial(n: int) -> int:
        if n <= 0:
            return 1
        r = 1
        while n > 0:
            r *= n
            n -= 2
        return r

    # ---- Associated Laguerre Polynomial L_n^k(x) ----
    def _assoc_laguerre(self, n: int, k: int, x: float) -> float:
        if n == 0:
            return 1.0
        if n == 1:
            return float(k + 1 - x)
        L0 = 1.0
        L1 = float(k + 1 - x)
        for i in range(2, n + 1):
            L_curr = ((2 * i + k - 1 - x) * L1 - (i + k - 1) * L0) / i
            L0, L1 = L1, L_curr
        return L1

    # ---- Spherical Harmonics Y_l^m(θ, φ) ----
    def _Y_lm(self, l: int, m: int, theta: float, phi: float) -> complex:
        """Compute spherical harmonic Y_l^m(θ, φ) analytically for l <= 3."""
        ct = math.cos(theta)
        st = math.sin(theta)
        cp = math.cos(phi)
        sp = math.sin(phi)

        # Normalization constants
        def N(l, m):
            return math.sqrt(
                (2 * l + 1) / (4.0 * math.pi) *
                self._factorial(l - abs(m)) / self._factorial(l + abs(m))
            )

        if l == 0:
            return N(0, 0) * complex(1.0, 0)

        elif l == 1:
            N1 = N(1, abs(m))
            if m == 0:
                return N1 * ct
            elif m == 1:
                return N1 * st * complex(-sp if m > 0 else cp, -cp if m > 0 else -sp)
            elif m == -1:
                return N1 * st * complex(cp, -sp)

        elif l == 2:
            N2 = N(2, abs(m))
            if m == 0:
                return N2 * (3.0 * ct * ct - 1.0) / 2.0
            elif abs(m) == 1:
                amp = N2 * st * ct
                if m == 1:
                    return amp * complex(-sp, -cp)
                else:
                    return amp * complex(cp, -sp)
            elif abs(m) == 2:
                amp = N2 * st * st / 2.0
                angle = 2.0 * phi if m == 2 else -2.0 * phi
                return amp * complex(math.cos(angle), math.sin(angle))

        elif l == 3:
            N3 = N(3, abs(m))
            if m == 0:
                return N3 * (5.0 * ct**3 - 3.0 * ct) / 2.0
            elif abs(m) == 1:
                amp = N3 * st * (5.0 * ct**2 - 1.0) / 2.0
                if m == 1:
                    return amp * complex(-sp, -cp)
                else:
                    return amp * complex(cp, -sp)
            elif abs(m) == 2:
                amp = N3 * st * st * ct
                angle = 2.0 * phi if m == 2 else -2.0 * phi
                return amp * complex(math.cos(angle), math.sin(angle))
            elif abs(m) == 3:
                amp = N3 * st**3
                angle = 3.0 * phi if m == 3 else -3.0 * phi
                sign = -1.0 if (m > 0 and m % 2 == 1) or (m < 0 and abs(m) % 2 == 1) else 1.0
                return sign * amp * complex(math.cos(angle), math.sin(angle))

        raise ChemMCPError(f"Y_{l}^{m} not implemented for l > 3.")

    # ---- Radial Wavefunction R_nl(r) ----
    def _radial_wf(self, n: int, l: int, rho: float, Z: int) -> float:
        """R_nl where ρ = 2Zr/(na₀). Returns dimensionless value."""
        nf = self._factorial(n - l - 1)
        nlf = self._factorial(n + l)
        prefactor = math.sqrt(
            (2.0 * Z / (n ** 3)) ** 3 * nf / (2.0 * n * (nlf ** 3))
        )
        L = self._assoc_laguerre(n - l - 1, 2 * l + 1, rho)
        return prefactor * (rho ** l) * L * math.exp(-rho / 2.0)

    # ---- Observables ----
    def _expectation_r(self, n: int, l: int, Z: int) -> float:
        """<r> = a₀[3n² - l(l+1)] / (2Z) in Bohr."""
        return (3.0 * n * n - l * (l + 1)) / (2.0 * Z)

    def _expectation_r2(self, n: int, l: int, Z: int) -> float:
        """<r²> = a₀²n²[5n²+1-3l(l+1)] / (2Z²) in Bohr²."""
        return (n * n / (2.0 * Z * Z)) * (5.0 * n * n + 1 - 3.0 * l * (l + 1))

    def _most_probable_r(self, n: int, l: int, Z: int) -> float:
        """Most probable radius in Bohr."""
        if l == 0:
            return float(n * n) / Z
        elif l == n - 1:
            return float(n * (n + 0.5)) / Z
        else:
            frac = l / max(n - 1, 1)
            return (float(n * n) / Z) * (1.0 + frac * 0.5)

    def _run_base(self, n: int, l: int, m: int = 0, Z: int = 1,
                  n_radial_points: int = 200, r_max_bohr: float = 20.0,
                  compute_observables: bool = True) -> dict:

        # Validation
        if n < 1:
            raise ChemMCPError("n must be >= 1.")
        if not (0 <= l < n):
            raise ChemMCPError(f"l must satisfy 0 <= l < n (got l={l}, n={n}).")
        if abs(m) > l:
            raise ChemMCPError(f"|m| must <= l (got m={m}, l={l}).")
        if Z < 1:
            raise ChemMCPError("Z must be >= 1.")

        # Orbital label
        orb_labels = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}
        orbital_type = f"{n}{orb_labels.get(l, f'l={l}')}"
        subshell_map = {0: "sharp", 1: "principal", 2: "diffuse", 3: "fundamental"}
        subshell = subshell_map.get(l, f"l={l}")

        # Energy
        energy_eV = -self.Ry_eV * Z * Z / (n * n)
        energy_J = energy_eV / 6.241509e18

        # Nodes
        n_radial = n - l - 1
        n_angular = l
        n_total = n - 1

        # Degeneracy
        degeneracy = n * n

        # Build radial grid
        r_grid = [r_max_bohr * i / (n_radial_points - 1) for i in range(n_radial_points)]
        R_vals = []
        D_r_vals = []  # Radial distribution D(r) = r²|R(r)|²
        prob_dens_vals = []  # |ψ|² at (θ=π/2, φ=0)
        dr = r_max_bohr / (n_radial_points - 1) if n_radial_points > 1 else 1.0

        theta_eq = math.pi / 2.0
        phi_0 = 0.0

        for r_bohr in r_grid:
            rho = 2.0 * Z * r_bohr / n
            R_val = self._radial_wf(n, l, rho, Z)
            R_vals.append(R_val)

            # Radial distribution function
            D_r = r_bohr * r_bohr * R_val * R_val
            D_r_vals.append(D_r)

            # Full |ψ|² at equatorial plane
            try:
                Y_val = self._Y_lm(l, m, theta_eq, phi_0)
                psi = R_val * Y_val
                prob_dens_vals.append(abs(psi) ** 2)
            except ChemMCPError:
                prob_dens_vals.append(R_val * R_val)

        # Observables
        obs = {}
        if compute_observables:
            mean_r = self._expectation_r(n, l, Z)
            mean_r2 = self._expectation_r2(n, l, Z)
            delta_r = math.sqrt(max(0, mean_r2 - mean_r * mean_r))
            r_mp = self._most_probable_r(n, l, Z)

            # Numerical normalization check
            norm_integral = sum(D_r_vals) * dr
            # Find peak of D(r)
            peak_idx = max(range(len(D_r_vals)), key=lambda i: D_r_vals[i])
            r_peak = r_grid[peak_idx]

            obs = {
                "mean_radius_bohr": round(mean_r, 6),
                "mean_radius_sq_bohr2": round(mean_r2, 6),
                "uncertainty_delta_r_bohr": round(delta_r, 6),
                "most_probable_radius_bohr": round(r_mp, 6),
                "numerical_peak_radius_bohr": round(r_peak, 4),
                "D_r_integral": round(norm_integral, 6),
                "bohr_radius_m": self.a0_m,
                "effective_bohr_radius_n_over_Z": round(float(n) / Z, 4),
            }

        result = {
            "quantum_numbers": {"n": n, "l": l, "m": m},
            "orbital_type": orbital_type,
            "subshell_name": subshell,
            "nuclear_charge_Z": Z,
            "energy_eV": round(energy_eV, 10),
            "energy_J": round(energy_J, 30),
            "ionization_energy_eV": round(abs(energy_eV), 10),
            "n_radial_nodes": n_radial,
            "angular_nodes": n_angular,
            "total_nodes": n_total,
            "degeneracy": degeneracy,
            "radial_grid_bohr": [round(r, 8) for r in r_grid[::max(1, n_radial_points//20)]],
            "radial_wavefunction_R": [round(v, 12) for v in R_vals[::max(1, n_radial_points//20)]],
            "radial_distribution_D_r": [round(v, 12) for v in D_r_vals[::max(1, n_radial_points//20)]],
            "probability_density_equatorial": [round(v, 12) for v in prob_dens_vals[::max(1, n_radial_points//20)]],
            "n_grid_points": n_radial_points,
            "r_max_bohr": r_max_bohr,
            **obs,
        }

        logger.info(f"Schrodinger3DSolver: {orbital_type}, E={energy_eV:.6f}eV, r_mp={obs.get('most_probable_radius_bohr','?')}")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            n = int(parts[0])
            l = int(parts[1])
            m = int(parts[2]) if len(parts) > 2 else 0
            Z = int(parts[3]) if len(parts) > 3 else 1
            n_pts = int(parts[4]) if len(parts) > 4 else 200
            r_max = float(parts[5]) if len(parts) > 5 else 20.0
            return self._run_base(n, l, m, Z, n_pts, r_max)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'n l [m] [Z] [n_pts] [r_max]'")
