import logging
import math
from typing import List, Optional, Tuple
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HydrogenWavefunction(BaseTool):
    """
    氢原子波函数计算工具。
    
    计算氢原子/类氢离子的完整波函数 ψ_nlm(r,θ,φ) 及相关性质：
    - 径向部分 R_nl(r)（关联拉盖尔多项式）
    - 角向部分 Y_l^m(θ,φ)（球谐函数）
    - 完整波函数及其概率密度 |ψ|²
    - 各阶矩、期望值、节点结构
    
    氢原子波函数是量子力学中少数有解析解的系统之一，
    是理解原子结构、化学键和光谱学的基础。
    """
    __version__ = "0.1.0"
    name = "HydrogenWavefunction"
    func_name = "hydrogen_wavefunction"
    description = "Compute hydrogen atom wavefunctions ψ_nlm(r,θ,φ): radial part R_nl, angular part Y_l^m, full wavefunction, probability density |ψ|², expectation values, nodal structure, and orbital visualization data for hydrogen-like atoms."
    implementation_description = "Analytical computation of hydrogenic wavefunctions via separation of variables: ψ_nl m = R_nl(r)·Y_l^m(θ,φ). Radial part uses associated Laguerre polynomials; angular part uses spherical harmonics (associated Legendre × azimuthal). Computes energy, normalization, <r>, <r²>, Δr, radial/angular nodes, probability in shells, and grid data for visualization."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Hydrogen Atom", "Wavefunction", "Atomic Orbitals", "Radial Function", "Spherical Harmonics", "Probability Density"]
    required_envs = []

    code_input_sig = [
        ("n", "int", "N/A", "Principal quantum number (n >= 1)."),
        ("l", "int", "N/A", "Orbital angular momentum quantum number (0 <= l < n)."),
        ("m", "int", "0", "Magnetic quantum number (-l <= m <= l, default 0)."),
        ("Z", "int", "1", "Nuclear charge for hydrogen-like ions (default 1 = H)."),
        ("r_bohr", "float", "1.0", "Radial coordinate in Bohr radii (for single-point evaluation)."),
        ("theta_rad", "float", "1.5708", "Polar angle θ in radians (default π/2 = equatorial plane)."),
        ("phi_rad", "float", "0.0", "Azimuthal angle φ in radians (default 0)."),
        ("n_grid_points", "int", "100", "Number of radial grid points for distribution computation."),
        ("r_max_bohr", "float", "20.0", "Maximum radius in Bohr for radial grid."),
        ("output_mode", "str", "full", "Output mode: 'full', 'single_point', 'radial_only', 'angular_only'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'n l [m] [Z] [r] [theta] [phi] [n_pts] [r_max] [mode]'. Example: '2 1 0 1 1.0 1.5708 0' for 2p at r=a₀, equatorial plane."),
    ]

    output_sig = [
        ("result", "dict", "Complete wavefunction data: energy, R_nl value, Y_lm value, full ψ, |ψ|², quantum numbers, node info, expectation values, grid data (if requested)."),
    ]

    examples = [
        {
            "code_input": {
                "n": 1,
                "l": 0,
                "m": 0,
                "Z": 1,
                "r_bohr": 1.0,
                "theta_rad": math.pi / 2,
                "phi_rad": 0.0,
                "n_grid_points": 100,
                "r_max_bohr": 10.0,
                "output_mode": "full",
            },
            "text_input": {
                "input_params": "1 0 0 1 1.0 1.5708 0",
            },
            "output": {
                "result": {
                    "orbital_type": "1s",
                    "energy_eV": -13.606,
                    "R_nl_at_r": round(0.3679, 4),  # e^{-1} ≈ 0.368 at r=a₀ for 1s
                    "psi_modulus_squared": round(0.0429, 4),  # |ψ|²/(a₀³) at r=a₀, θ=π/2
                    "most_probable_radius_bohr": 1.0,
                }
            },
        },
        {
            "code_input": {
                "n": 2,
                "l": 1,
                "m": 0,
                "Z": 1,
                "r_bohr": 2.0,
                "theta_rad": 0.0,
                "phi_rad": 0.0,
                "output_mode": "single_point",
            },
            "text_input": {
                "input_params": "2 1 0 1 2.0 0",
            },
            "output": {
                "result": {
                    "orbital_type": "2p_z",
                    "energy_eV": -3.4015,
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
        """Compute associated Laguerre polynomial L_n^k(x)."""
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

    # ---- Radial Wavefunction R_nl(r) ----
    def _R_nl(self, n: int, l: int, rho: float, Z: int) -> float:
        """
        Radial wavefunction R_nl where ρ = 2Zr/(na₀).
        Returns dimensionless value (units: a₀^{-3/2}).
        
        Formula: R_nl = sqrt((2Z/n³)² · (n-l-1)! / (2n·[(n+l)!]³)) · ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2)
        """
        nf = self._factorial(n - l - 1)
        nlf = self._factorial(n + l)
        prefactor = math.sqrt(
            (2.0 * Z / (n ** 3)) ** 3 * nf / (2.0 * n * (nlf ** 3))
        )
        L = self._assoc_laguerre(n - l - 1, 2 * l + 1, rho)
        return prefactor * (rho ** l) * L * math.exp(-rho / 2.0)

    # ---- Spherical Harmonics Y_l^m(θ, φ) ----
    def _Y_lm(self, l: int, m: int, theta: float, phi: float) -> complex:
        """Compute spherical harmonic Y_l^m(θ, φ) analytically for l <= 5."""
        ct = math.cos(theta)
        st = math.sin(theta)
        cp = math.cos(phi)
        sp = math.sin(phi)

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

        elif l == 4:
            N4 = N(4, abs(m))
            if m == 0:
                return N4 * (35.0*ct**4 - 30.0*ct**2 + 3.0) / 8.0
            elif abs(m) == 1:
                return N4 * st * (7.0*ct**3 - 3.0*ct) * (complex(-sp, -cp) if m == 1 else complex(cp, -sp))
            elif abs(m) == 2:
                amp = N4 * st * st * (7.0*ct**2 - 1.0) / 2.0
                angle = 2.0 * phi if m == 2 else -2.0 * phi
                return amp * complex(math.cos(angle), math.sin(angle))
            elif abs(m) == 3:
                amp = N4 * st**3 * ct
                angle = 3.0 * phi if m == 3 else -3.0 * phi
                sgn = -1 if (abs(m) % 2 == 1 and m > 0) or (abs(m) % 2 == 1 and m < 0) else 1
                return sgn * amp * complex(math.cos(angle), math.sin(angle))
            elif abs(m) == 4:
                amp = N4 * st**4 / 24.0
                angle = 4.0 * phi if m == 4 else -4.0 * phi
                return amp * complex(math.cos(angle), math.sin(angle))

        elif l == 5:
            N5 = N(5, abs(m))
            if m == 0:
                return N5 * (63.0*ct**5 - 70.0*ct**3 + 15.0*ct) / 8.0
            # For higher l with m≠0, use generic approach below

        # Fallback: use associated Legendre for any remaining cases
        P = self._associated_legendre(l, abs(m), ct)
        norm = N(l, abs(m))
        azim_real = math.cos(m * phi)
        azim_imag = math.sin(m * phi)
        if m >= 0:
            return norm * P * complex(azim_real, azim_imag)
        else:
            phase = 1.0 if (abs(m) % 2 == 0) else -1.0
            return phase * norm * P * complex(azim_real, -azim_imag)

    def _associated_legendre(self, l: int, m: int, x: float) -> float:
        """Associated Legendre polynomial P_l^|m|(x) via recursion."""
        abs_m = abs(m)
        sign = 1.0 if (abs_m % 2 == 0) else -1.0
        p_mm = sign * self._double_factorial(2 * abs_m - 1)
        one_minus_x2 = max(0.0, 1.0 - x * x)
        p_mm *= one_minus_x2 ** (abs_m / 2.0)
        if l == abs_m:
            return p_mm
        p_mp1 = x * (2 * abs_m + 1) * p_mm
        if l == abs_m + 1:
            return p_mp1
        p_lm2 = p_mm
        p_lm1 = p_mp1
        for ll in range(abs_m + 2, l + 1):
            p_ll = (x * (2 * ll - 1) * p_lm1 - (ll + abs_m - 1) * p_lm2) / (ll - abs_m)
            p_lm2, p_lm1 = p_lm1, p_ll
        return p_lm1

    # ---- Analytical Expectation Values ----
    def _expectation_r(self, n: int, l: int, Z: int) -> float:
        """⟨r⟩ = a₀[3n² - l(l+1)] / (2Z) in Bohr."""
        return (3.0 * n * n - l * (l + 1)) / (2.0 * Z)

    def _expectation_r2(self, n: int, l: int, Z: int) -> float:
        """⟨r²⟩ = a₀²n²[5n²+1-3l(l+1)] / (2Z²) in Bohr²."""
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

    # ---- Full Wavefunction ψ = R × Y ----
    def _full_psi(self, n: int, l: int, m: int, r_bohr: float,
                   theta: float, phi: float, Z: int) -> complex:
        """Full wavefunction ψ_nlm(r,θ,φ) = R_nl(r) · Y_l^m(θ,φ). Units: a₀^{-3/2}."""
        rho = 2.0 * Z * r_bohr / n
        R_val = self._R_nl(n, l, rho, Z)
        Y_val = self._Y_lm(l, m, theta, phi)
        return R_val * Y_val

    # ---- Radial Probability Density D(r) = r²|R|² ----
    def _D_r(self, n: int, l: int, r_bohr: float, Z: int) -> float:
        """Radial distribution function D(r) = r²|R_nl(r)|²."""
        rho = 2.0 * Z * r_bohr / n
        R = self._R_nl(n, l, rho, Z)
        return r_bohr * r_bohr * R * R

    def _run_base(self, n: int, l: int, m: int = 0, Z: int = 1,
                  r_bohr: float = 1.0, theta_rad: float = math.pi / 2.0,
                  phi_rad: float = 0.0, n_grid_points: int = 100,
                  r_max_bohr: float = 20.0,
                  output_mode: str = "full") -> dict:

        # --- Validation ---
        if n < 1:
            raise ChemMCPError("n must be >= 1.")
        if not (0 <= l < n):
            raise ChemMCPError(f"l must satisfy 0 <= l < n (got l={l}, n={n}).")
        if abs(m) > l:
            raise ChemMCPError(f"|m| must <= l (got m={m}, l={l}).")
        if Z < 1:
            raise ChemMCPError("Z must be >= 1.")
        if r_bohr < 0:
            raise ChemMCPError("r_bohr must be non-negative.")

        mode = output_mode.lower().strip()

        # --- Orbital classification ---
        orb_labels = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}
        orbital_type = f"{n}{orb_labels.get(l, f'l={l}')}"
        subshell_map = {0: "sharp", 1: "principal", 2: "diffuse", 3: "fundamental"}
        subshell = subshell_map.get(l, f"l={l}")

        # --- Energy ---
        energy_eV = -self.Ry_eV * Z * Z / (n * n)
        energy_J = energy_eV / 6.241509e18

        # --- Nodes ---
        n_radial_nodes = n - l - 1
        n_angular_nodes = l
        n_total_nodes = n - 1
        degeneracy = n * n

        # --- Single-point evaluation ---
        rho_at_r = 2.0 * Z * r_bohr / n
        R_val = self._R_nl(n, l, rho_at_r, Z)
        Y_val = self._Y_lm(l, m, theta_rad, phi_rad)
        psi_val = R_val * Y_val
        psi_mod_sq = abs(psi_val) ** 2

        # --- Observables (analytical) ---
        mean_r = self._expectation_r(n, l, Z)
        mean_r2 = self._expectation_r2(n, l, Z)
        delta_r = math.sqrt(max(0, mean_r2 - mean_r * mean_r))
        r_mp = self._most_probable_r(n, l, Z)

        result = {
            "quantum_numbers": {"n": n, "l": l, "m": m},
            "orbital_type": orbital_type,
            "subshell_name": subshell,
            "nuclear_charge_Z": Z,
            "energy_eV": round(energy_eV, 10),
            "energy_J": round(energy_J, 30),
            "ionization_energy_eV": round(abs(energy_eV), 10),
            "degeneracy": degeneracy,
            # Single-point values
            "evaluation_point": {
                "r_bohr": r_bohr,
                "theta_rad": round(theta_rad, 10),
                "phi_rad": round(phi_rad, 10),
                "rho": round(rho_at_r, 10),
            },
            "R_nl_value": round(R_val, 14),
            "Y_lm_complex_real": round(Y_val.real, 14),
            "Y_lm_complex_imag": round(Y_val.imag, 14),
            "psi_complex_real": round(psi_val.real, 14),
            "psi_complex_imag": round(psi_val.imag, 14),
            "psi_modulus": round(abs(psi_val), 14),
            "psi_modulus_squared": round(psi_mod_sq, 14),
            # Nodes
            "n_radial_nodes": n_radial_nodes,
            "n_angular_nodes": n_angular_nodes,
            "total_nodes": n_total_nodes,
            # Observables
            "mean_radius_bohr": round(mean_r, 8),
            "mean_radius_sq_bohr2": round(mean_r2, 8),
            "uncertainty_delta_r_bohr": round(delta_r, 8),
            "most_probable_radius_bohr": round(r_mp, 6),
            "bohr_radius_m": self.a0_m,
        }

        # --- Grid computation (for 'full' and 'radial_only' modes) ---
        if mode in ("full", "radial_only"):
            dr = r_max_bohr / (n_grid_points - 1) if n_grid_points > 1 else r_max_bohr
            r_grid = [r_max_bohr * i / (n_grid_points - 1) for i in range(n_grid_points)]

            R_grid = []
            D_grid = []       # D(r) = r²R²
            psi_eq_grid = []  # |ψ|² at (θ=π/2, φ=0)
            theta_eq = math.pi / 2.0

            for r in r_grid:
                rho = 2.0 * Z * r / n
                Rv = self._R_nl(n, l, rho, Z)
                R_grid.append(Rv)
                D_grid.append(r * r * Rv * Rv)
                try:
                    Yeq = self._Y_lm(l, m, theta_eq, 0.0)
                    psieq = Rv * Yeq
                    psi_eq_grid.append(abs(psieq) ** 2)
                except (ChemMCPError, ValueError):
                    psi_eq_grid.append(Rv * Rv)

            # Normalization check
            norm_D = sum(D_grid) * dr

            # Find peaks of D(r)
            peaks = []
            for i in range(1, len(D_grid) - 1):
                if D_grid[i] > D_grid[i-1] and D_grid[i] > D_grid[i+1] and D_grid[i] > 1e-12:
                    peaks.append((round(r_grid[i], 4), round(D_grid[i], 12)))

            # Find radial nodes (where R ≈ 0)
            radial_node_radii = []
            for i in range(1, len(R_grid) - 1):
                if R_grid[i] * R_grid[i-1] < 0 and abs(R_grid[i]) < abs(R_grid[i-1]):
                    # Linear interpolation for node position
                    t = abs(R_grid[i-1]) / (abs(R_grid[i-1]) + abs(R_grid[i]))
                    r_node = r_grid[i-1] + t * dr
                    radial_node_radii.append(round(r_node, 6))

            step = max(1, n_grid_points // 25)
            result["radial_grid_data"] = {
                "n_points": n_grid_points,
                "r_max_bohr": r_max_bohr,
                "dr_bohr": round(dr, 8),
                "normalization_D_integral": round(norm_D, 8),
                "peaks": peaks[:10],  # Top 10 peaks
                "radial_node_radii": radial_node_radii,
                "r_sample": [round(r, 6) for r in r_grid[::step]],
                "R_nl_sample": [round(v, 12) for v in R_grid[::step]],
                "D_r_sample": [round(d, 12) for d in D_grid[::step]],
                "probability_density_equatorial_sample": [round(v, 12) for v in psi_eq_grid[::step]],
            }

        # --- Angular grid (for 'full' mode) ---
        if mode == "full":
            n_theta = min(45, max(10, n_grid_points // 2))
            theta_grid = [math.pi * i / (n_theta - 1) for i in range(n_theta)]
            ang_data = []
            for th in theta_grid[::max(1, n_theta // 15)]:
                try:
                    Ya = self._Y_lm(l, m, th, phi_rad)
                    ang_data.append({
                        "theta_rad": round(th, 6),
                        "theta_deg": round(math.degrees(th), 4),
                        "Y_real": round(Ya.real, 12),
                        "Y_imag": round(Ya.imag, 12),
                        "|Y|^2": round(abs(Ya)**2, 12),
                    })
                except ChemMCPError:
                    ang_data.append({"theta_rad": round(th, 6), "error": True})
            result["angular_grid_data"] = {
                "n_theta": n_theta,
                "phi_evaluated_rad": round(phi_rad, 6),
                "data": ang_data,
            }

        # --- Mode filtering ---
        if mode == "single_point":
            return {"result": {
                "orbital_type": orbital_type,
                "quantum_numbers": {"n": n, "l": l, "m": m},
                "energy_eV": round(energy_eV, 10),
                "R_nl": round(R_val, 14),
                "Y_lm": {"real": round(Y_val.real, 14), "imag": round(Y_val.imag, 14)},
                "psi": {"real": round(psi_val.real, 14), "imag": round(psi_val.imag, 14)},
                "|psi|^2": round(psi_mod_sq, 14),
                "|psi|": round(abs(psi_val), 14),
            }}
        elif mode == "radial_only":
            radial_keys = ["quantum_numbers", "orbital_type", "nuclear_charge_Z",
                           "energy_eV", "n_radial_nodes", "mean_radius_bohr",
                           "most_probable_radius_bohr", "radial_grid_data"]
            return {k: result[k] for k in radial_keys if k in result}
        elif mode == "angular_only":
            ang_keys = ["quantum_numbers", "orbital_type", "Y_lm_complex_real",
                        "Y_lm_complex_imag", "n_angular_nodes", "total_nodes",
                        "angular_grid_data"]
            return {k: result[k] for k in ang_keys if k in result}
        else:
            logger.info(f"HydrogenWavefunction: {orbital_type}, E={energy_eV:.4f}eV, |ψ|²={psi_mod_sq:.6g}")
            return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            n = int(parts[0])
            l = int(parts[1])
            m = int(parts[2]) if len(parts) > 2 else 0
            Z = int(parts[3]) if len(parts) > 3 else 1
            r = float(parts[4]) if len(parts) > 4 else 1.0
            th = float(parts[5]) if len(parts) > 5 else math.pi / 2
            ph = float(parts[6]) if len(parts) > 6 else 0.0
            n_pts = int(parts[7]) if len(parts) > 7 else 100
            r_max = float(parts[8]) if len(parts) > 8 else 20.0
            mode = parts[9] if len(parts) > 9 else "full"
            return self._run_base(n, l, m, Z, r, th, ph, n_pts, r_max, mode)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'n l [m] [Z] [r] [theta] [phi] [n_pts] [r_max] [mode]'"
            )
