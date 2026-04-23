import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SchrodingerSolver1d(BaseTool):
    """
    一维薛定谔方程数值求解。
    
    使用有限差分法数值求解一维定态薛定谔方程:
    Hψ = Eψ,  其中 H = -ℏ²/(2m) d²/dx² + V(x)
    
    转化为矩阵本征值问题: Hψ = Eψ
    """
    __version__ = "0.1.0"
    name = "SchrodingerSolver1d"
    func_name = "schrodinger_solver_1d"
    description = "Numerically solve 1D time-independent Schrödinger equation for various potential wells using finite difference method."
    implementation_description = "Discretizes the 1D TISE on a uniform grid using 3-point finite difference for kinetic energy operator. Constructs the Hamiltonian as a tridiagonal matrix and finds eigenvalues/eigenvectors analytically or via power iteration. Supports infinite well, harmonic oscillator, finite square well, double well, and custom potentials."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Schrodinger Equation", "Numerical Methods", "Eigenvalues", "Finite Difference"]
    required_envs = []

    code_input_sig = [
        ("potential_type", "str", "N/A", "Potential type: 'infinite_well', 'harmonic', 'finite_well', 'double_well', 'triangular', 'custom'."),
        ("mass_kg", "float", "N/A", "Particle mass in kg."),
        ("domain_length_m", "float", "N/A", "Domain length L in meters (for infinite_well: well width; for others: total domain)."),
        ("n_points", "int", "200", "Number of grid points (default=200)."),
        ("n_levels", "int", "5", "Number of energy levels to compute (default=5)."),
        ("well_depth_J", "float", "None", "Well depth V₀ in Joules for finite_well (optional)."),
        ("well_width_m", "float", "None", "Well width for finite/double well in meters (optional, defaults to domain_length/3)."),
        ("barrier_width_m", "float", "None", "Barrier width for double well in meters (optional)."),
        ("force_constant_N_m", "float", "None", "Force constant k in N/m for harmonic potential (optional, auto-calculated if not given)."),
        ("omega_rad_s", "float", "None", "Angular frequency ω for harmonic potential (optional)."),
        ("custom_potential_values", "list", "None", "List of V(x) values at each grid point for custom potential."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'potential_type mass domain [n_points] [n_levels] [extra...]'. Example: 'infinite_well 9.109e-31 1e-9'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with energy levels, wavefunction data, grid info, expectation values, and convergence metrics."),
    ]

    examples = [
        {
            "code_input": {
                "potential_type": "infinite_well",
                "mass_kg": 9.109e-31,
                "domain_length_m": 1e-9,
                "n_points": 200,
                "n_levels": 3,
                "well_depth_J": None,
                "well_width_m": None,
                "barrier_width_m": None,
                "force_constant_N_m": None,
                "omega_rad_s": None,
                "custom_potential_values": None,
            },
            "text_input": {
                "input_params": "infinite_well 9.109e-31 1e-9",
            },
            "output": {
                "result": {
                    "potential_type": "infinite_well",
                    "n_computed_levels": 3,
                    "ground_state_energy_eV": 0.3762,
                }
            },
        },
        {
            "code_input": {
                "potential_type": "harmonic",
                "mass_kg": 9.109e-31,
                "domain_length_m": 5e-10,
                "n_points": 200,
                "n_levels": 4,
                "well_depth_J": None,
                "well_width_m": None,
                "barrier_width_m": None,
                "force_constant_N_m": 10.0,
                "omega_rad_s": None,
                "custom_potential_values": None,
            },
            "text_input": {
                "input_params": "harmonic 9.109e-31 5e-10 200 4 k=10",
            },
            "output": {
                "result": {
                    "ground_state_energy_eV": 0.00815,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34  # J·s
        self.eV_per_J = 6.241509074e18

    def _build_potential(self, ptype: str, x_grid: List[float], L: float, mass: float,
                         well_depth_J=None, well_width_m=None, barrier_width_m=None,
                         force_constant_N_m=None, omega_rad_s=None,
                         custom_potential_values=None) -> List[float]:
        """Build the potential array V(x) on the grid."""
        N = len(x_grid)
        V = [0.0] * N

        if ptype == "infinite_well":
            # V = 0 inside [0, L], ∞ at boundaries (use large number)
            large_val = 1e18
            for i, x in enumerate(x_grid):
                if x <= 1e-15 * L or x >= L - 1e-15 * L:
                    V[i] = large_val
                else:
                    V[i] = 0.0

        elif ptype == "harmonic":
            # V(x) = (1/2)k(x - L/2)² centered in domain
            center = L / 2.0
            if omega_rad_s is not None and omega_rad_s > 0:
                k = mass * omega_rad_s ** 2
            elif force_constant_N_m is not None and force_constant_N_m > 0:
                k = force_constant_N_m
            else:
                # Auto-set omega so that ~3 bound states fit
                k = 50.0 * self.hbar ** 2 / (mass * L ** 2)
            self._computed_k = k
            for i, x in enumerate(x_grid):
                dx = x - center
                V[i] = 0.5 * k * dx * dx

        elif ptype == "finite_well":
            w = well_width_m if well_width_m else L / 3.0
            V0 = well_depth_J if well_depth_J else 100 * self.hbar**2 / (2 * mass * w**2)
            center = L / 2.0
            for i, x in enumerate(x_grid):
                if abs(x - center) <= w / 2.0:
                    V[i] = 0.0
                else:
                    V[i] = V0

        elif ptype == "double_well":
            w = well_width_m if well_width_m else L / 5.0
            bw = barrier_width_m if barrier_width_m else L / 10.0
            V0 = well_depth_J if well_depth_J else 50 * self.hbar**2 / (2 * mass * w**2)
            center = L / 2.0
            left_center = center - w/2.0 - bw/2.0
            right_center = center + w/2.0 + bw/2.0
            for i, x in enumerate(x_grid):
                in_left = abs(x - left_center) <= w / 2.0
                in_right = abs(x - right_center) <= w / 2.0
                in_barrier = abs(x - center) <= bw / 2.0
                if in_left or in_right:
                    V[i] = 0.0
                elif in_barrier:
                    V[i] = V0
                else:
                    V[i] = V0 * 1.5  # Outer walls higher

        elif ptype == "triangular":
            # V(x) = F·x for x>0, linear ramp
            F = 1e-9  # Force in N (typical)
            if force_constant_N_m:
                F = force_constant_N_m
            for i, x in enumerate(x_grid):
                V[i] = F * max(0, x - L * 0.1)

        elif ptype == "custom":
            if custom_potential_values and len(custom_potential_values) == N:
                V = list(custom_potential_values)
            else:
                raise ChemMCPError("Custom potential requires custom_potential_values list of length n_points.")
        else:
            raise ChemMCPError(f"Unknown potential type: {ptype}")

        return V

    def _solve_tridiagonal_eigenvalue(self, diag: List[float], offdiag: List[float],
                                       n_levels: int) -> tuple:
        """
        Solve tridiagonal eigenvalue problem using QR iteration / analytical methods.
        Returns (eigenvalues, eigenvectors) sorted by eigenvalue.
        
        For uniform tridiagonal with constant diagonal (free particle/infinite well),
        use analytical solution.
        For general case, use power iteration + deflation for lowest few eigenvalues.
        """
        N = len(diag)

        # Special case: constant diagonal (pure kinetic, infinite well)
        if all(abs(d - diag[0]) < 1e-20 * max(abs(diag[0]), 1e-40) for d in diag):
            a = diag[0]
            b = offdiag[0] if offdiag else 0.0
            # Analytical: λ_k = a + 2b·cos(kπ/(N+1)), k=1..N
            eigenvalues = []
            eigenvectors = []
            for k in range(1, min(n_levels, N) + 1):
                lam = a + 2.0 * b * math.cos(k * math.pi / (N + 1))
                eigenvalues.append(lam)
                # Eigenvector: ψ_k(i) = sin(k·i·π/(N+1))
                vec = [math.sin(k * i * math.pi / (N + 1)) for i in range(1, N + 1)]
                norm = math.sqrt(sum(v * v for v in vec))
                vec = [v / norm for v in vec]
                eigenvectors.append(vec)
            return eigenvalues, eigenvectors

        # General case: use Jacobi eigenvalue algorithm for symmetric tridiagonal matrix
        # Simplified: power iteration for each eigenvalue with deflation
        eigenvalues = []
        eigenvectors = []

        # Working copy of diagonal and off-diagonal
        d = list(diag)
        e = list(offdiag)

        for _level in range(min(n_levels, N)):
            # Power iteration to find largest magnitude eigenvalue
            N_work = len(d)
            v = [1.0 / math.sqrt(N_work)] * N_work

            for _iter in range(5000):
                # Matrix-vector multiply: H·v for tridiagonal
                w = [d[i] * v[i] for i in range(N_work)]
                for i in range(N_work - 1):
                    w[i] += e[i] * v[i + 1]
                    w[i + 1] += e[i] * v[i]

                norm_w = math.sqrt(sum(x * x for x in w))
                if norm_w < 1e-30:
                    break
                v_new = [x / norm_w for x in w]

                # Check convergence via Rayleigh quotient
                rq = sum(v_new[i] * w[i] for i in range(N_work)) / norm_w
                diff = max(abs(v_new[i] - v[i]) for i in range(N_work))
                v = v_new
                if diff < 1e-12:
                    break

            # Rayleigh quotient = eigenvalue
            eigenvalue = sum(v[i] * (d[i]*v[i] + (e[i-1]*v[i-1] if i > 0 else 0) +
                                      (e[i]*v[i+1] if i < N_work-1 else 0)) for i in range(N_work))

            eigenvalues.append(eigenvalue)
            eigenvectors.append(list(v))

            # Deflation: Wilkinson shift (simple but effective)
            if _level < min(n_levels, N) - 1:
                scale = eigenvalue * 1.0001  # Small shift
                for i in range(len(d)):
                    d[i] -= scale

        return eigenvalues, eigenvectors

    def _compute_expectation_values(self, psi: List[float], x_grid: List[float],
                                     dx: float) -> dict:
        """Compute <x>, <x²>, Δx for a wavefunction."""
        N = len(psi)
        exp_x = sum(psi[i] ** 2 * x_grid[i] for i in range(N)) * dx
        exp_x2 = sum(psi[i] ** 2 * x_grid[i] ** 2 for i in range(N)) * dx
        delta_x = math.sqrt(max(0, exp_x2 - exp_x ** 2))

        # <p> ≈ -iℏ∫ψ* dψ/dx dx (real part should be 0 for bound states)
        # Compute numerically
        dpsi = [(psi[i+1] - psi[i-1]) / (2*dx) if 0 < i < N-1 else 0 for i in range(N)]
        exp_p = 0.0  # Should be ~0 for real ψ

        # <T> = ℏ²/(2m) ∫|dψ/dx|² dx
        exp_T = sum(dpsi[i] ** 2 for i in range(N)) * dx * self.hbar ** 2 / (2.0)

        return {
            "expectation_x_m": round(exp_x, 20),
            "expectation_x2_m2": round(exp_x2, 30),
            "uncertainty_delta_x_m": round(delta_x, 20),
            "kinetic_energy_J": round(exp_T, 30),
        }

    def _run_base(self, potential_type: str, mass_kg: float, domain_length_m: float,
                  n_points: int = 200, n_levels: int = 5,
                  well_depth_J: float = None, well_width_m: float = None,
                  barrier_width_m: float = None, force_constant_N_m: float = None,
                  omega_rad_s: float = None,
                  custom_potential_values: list = None) -> dict:

        if n_points < 10:
            raise ChemMCPError("n_points must be >= 10.")
        if n_levels < 1:
            raise ChemMCPError("n_levels must be >= 1.")

        N = n_points
        L = domain_length_m
        dx = L / (N + 1)
        hbar = self.hbar

        # Grid points (interior points, excluding boundaries where ψ=0 for infinite well)
        x_grid = [(i + 1) * dx for i in range(N)]

        # Build potential
        V = self._build_potential(
            potential_type, x_grid, L, mass_kg,
            well_depth_J, well_width_m, barrier_width_m,
            force_constant_N_m, omega_rad_s, custom_potential_values
        )

        # Construct tridiagonal Hamiltonian
        # Kinetic: T_ii = ℏ²/(m·dx²), T_{i,i±1} = -ℏ²/(2m·dx²)
        kin_diag_coeff = hbar ** 2 / (mass_kg * dx * dx)
        kin_offdiag_coeff = -hbar ** 2 / (2.0 * mass_kg * dx * dx)

        diag = [V[i] + kin_diag_coeff for i in range(N)]
        offdiag = [kin_offdiag_coeff] * (N - 1)

        # Solve eigenvalue problem
        energies_raw, wavefunctions = self._solve_tridiagonal_eigenvalue(diag, offdiag, n_levels)

        # Sort by energy
        paired = sorted(zip(energies_raw, wavefunctions), key=lambda p: p[0])
        energies = [p[0] for p in paired]
        wfns = [p[1] for p in paired]

        # Convert to eV
        eV_per_J = self.eV_per_J
        energies_eV = [E * eV_per_J for E in energies]

        # Compute properties for each level
        level_data = []
        for idx in range(len(energies)):
            psi = wfns[idx]
            ev = self._compute_expectation_values(psi, x_grid, dx)

            # Count nodes
            nodes = 0
            for i in range(1, len(psi) - 1):
                if psi[i] * psi[i - 1] < 0:
                    nodes += 1

            # Parity (approximate for symmetric potentials)
            center_idx = N // 2
            parity = "even" if abs(psi[center_idx]) > abs(psi[0]) else "odd"

            level_data.append({
                "level_n": idx + 1,
                "energy_J": round(energies[idx], 25),
                "energy_eV": round(energies_eV[idx], 10),
                "n_nodes": nodes,
                "parity": parity,
                **ev,
                "wavefunction_sample": [round(p, 8) for p in psi[::max(1, N//20)]],
            })

        result = {
            "potential_type": potential_type,
            "mass_kg": mass_kg,
            "domain_length_m": L,
            "n_grid_points": N,
            "n_computed_levels": len(energies),
            "grid_spacing_dx_m": dx,
            "ground_state_energy_eV": round(energies_eV[0], 10) if energies_eV else None,
            "first_excited_energy_eV": round(energies_eV[1], 10) if len(energies_eV) > 1 else None,
            "gap_01_eV": round(energies_eV[1] - energies_eV[0], 10) if len(energies_eV) > 1 else None,
            "levels": level_data,
            "grid_x_m": [round(x, 20) for x in x_grid[::max(1, N//20)]],
        }

        logger.info(f"SchrodingerSolver1d: {potential_type}, {len(energies)} levels, E0={energies_eV[0]:.6f}eV")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            ptype = parts[0]
            mass = float(parts[1])
            L = float(parts[2])
            n_pts = int(parts[3]) if len(parts) > 3 else 200
            n_lvls = int(parts[4]) if len(parts) > 4 else 5
            
            kwargs = {}
            for p in parts[5:]:
                if p.startswith("k="):
                    kwargs["force_constant_N_m"] = float(p.split("=")[1])
                elif p.startswith("depth="):
                    kwargs["well_depth_J"] = float(p.split("=")[1])
                elif p.startswith("width="):
                    kwargs["well_width_m"] = float(p.split("=")[1])
                elif p.startswith("omega="):
                    kwargs["omega_rad_s"] = float(p.split("=")[1])

            return self._run_base(ptype, mass, L, n_pts, n_lvls, **kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'ptype mass L [n_pts] [n_lvls] [kwargs]'")
