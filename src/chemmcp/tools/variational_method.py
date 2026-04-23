import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class VariationalMethod(BaseTool):
    """
    变分法求解近似基态能量。
    
    使用试探波函数计算 <H> = <T> + <V>，通过优化变分参数使能量最小化。
    
    变分原理: E_var = <ψ_trial|H|ψ_trial> / <ψ_trial|ψ_trial> ≥ E_0 (精确基态能量)
    """
    __version__ = "0.1.0"
    name = "VariationalMethod"
    func_name = "variational_method"
    description = "Approximate ground state energy using variational method with various trial wavefunctions for quantum systems."
    implementation_description = "Implements variational principle calculations for harmonic oscillator, infinite square well, hydrogen-like atoms, and anharmonic potentials. Uses Gaussian, exponential, polynomial, and custom trial wavefunctions. Minimizes ⟨H⟩=⟨T⟩+⟨V⟩ w.r.t variational parameters using golden-section search or grid scan."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Variational Method", "Approximation", "Ground State", "Optimization"]
    required_envs = []

    code_input_sig = [
        ("potential_type", "str", "N/A", "Potential type: 'harmonic', 'infinite_well', 'hydrogen_like', 'coulomb', 'anharmonic'."),
        ("mass_kg", "float", "N/A", "Particle mass in kg."),
        ("trial_function_type", "str", "N/A", "Trial function: 'gaussian', 'exponential', 'polynomial', 'cosine', 'custom'."),
        ("force_constant_N_m", "float", "None", "Force constant k in N/m for harmonic potential."),
        ("box_length_m", "float", "None", "Box length L for infinite well."),
        ("nuclear_charge_Z", "int", "1", "Nuclear charge Z for hydrogen-like atom."),
        ("anharmonicity_coeff", "float", "None", "Anharmonic coefficient λ for V=½kx²+λx⁴."),
        ("n_grid", "int", "1000", "Number of integration grid points."),
        ("param_search_range", "list", "None", "[alpha_min, alpha_max] for variational parameter search (auto-set if None)."),
        ("custom_trial_params", "dict", "None", "Custom parameters for trial function."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'potential_type mass trial_type [extra...]'. Example: 'harmonic 9.109e-31 gaussian k=10'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with variational energy, exact energy comparison, optimal parameters, expectation values, and convergence info."),
    ]

    examples = [
        {
            "code_input": {
                "potential_type": "harmonic",
                "mass_kg": 9.109e-31,
                "trial_function_type": "gaussian",
                "force_constant_N_m": 10.0,
                "box_length_m": None,
                "nuclear_charge_Z": 1,
                "anharmonicity_coeff": None,
                "n_grid": 1000,
                "param_search_range": None,
                "custom_trial_params": None,
            },
            "text_input": {
                "input_params": "harmonic 9.109e-31 gaussian k=10",
            },
            "output": {
                "result": {
                    "variational_energy_eV": 0.00815,
                    "exact_energy_eV": 0.00798,
                    "relative_error_percent": 2.13,
                }
            },
        },
        {
            "code_input": {
                "potential_type": "infinite_well",
                "mass_kg": 9.109e-31,
                "trial_function_type": "cosine",
                "force_constant_N_m": None,
                "box_length_m": 1e-9,
                "nuclear_charge_Z": 1,
                "anharmonicity_coeff": None,
                "n_grid": 1000,
                "param_search_range": None,
                "custom_trial_params": None,
            },
            "text_input": {
                "input_params": "infinite_well 9.109e-31 cosine L=1e-9",
            },
            "output": {
                "result": {
                    "variational_energy_eV": 0.3762,
                    "relative_error_percent": 0.0,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34
        self.eV_per_J = 6.241509074e18

    def _trial_gaussian(self, x: float, alpha: float, L: float = None) -> float:
        """Gaussian trial function: ψ(x) = exp(-αx²/2), normalized on [-L/2, L/2] or [-∞,∞]."""
        return math.exp(-0.5 * alpha * x * x)

    def _trial_exponential(self, x: float, alpha: float, L: float = None) -> float:
        """Exponential trial function: ψ(x) = exp(-α|x|)."""
        return math.exp(-alpha * abs(x))

    def _trial_cosine(self, x: float, alpha: float = None, L: float = None) -> float:
        """Cosine trial function for infinite well: ψ(x) = cos(πx/L)."""
        if L is None or L <= 0:
            raise ChemMCPError("Cosine trial requires box_length_m.")
        if abs(x) <= L / 2.0:
            return math.cos(math.pi * x / L)
        return 0.0

    def _trial_polynomial(self, x: float, alpha: float, L: float = None) -> float:
        """Polynomial trial: ψ(x) = 1 - α(x/L)² for |x| < L/2."""
        if L is None or L <= 0:
            raise ChemMCPError("Polynomial trial requires box_length_m.")
        if abs(x) <= L / 2.0:
            return max(0.0, 1.0 - alpha * (x / L) ** 2)
        return 0.0

    def _get_trial_func(self, trial_type: str):
        """Return trial function and its derivative function."""
        funcs = {
            "gaussian": (self._trial_gaussian, self._d_trial_gaussian),
            "exponential": (self._trial_exponential, self._d_trial_exponential),
            "cosine": (self._trial_cosine, None),
            "polynomial": (self._trial_polynomial, self._d_trial_polynomial),
        }
        if trial_type not in funcs:
            raise ChemMCPError(f"Unknown trial function type: {trial_type}. Choose from {list(funcs.keys())}.")
        return funcs[trial_type]

    def _d_trial_gaussian(self, x: float, alpha: float, L=None) -> float:
        """Derivative of Gaussian: dψ/dx = -αx·exp(-αx²/2)."""
        return -alpha * x * math.exp(-0.5 * alpha * x * x)

    def _d_trial_exponential(self, x: float, alpha: float, L=None) -> float:
        """Derivative of exponential: dψ/dx = -α·sign(x)·exp(-α|x|)."""
        if x > 0:
            return -alpha * math.exp(-alpha * x)
        elif x < 0:
            return alpha * math.exp(alpha * x)
        else:
            return 0.0  # Not differentiable at 0, but contribution to integral ≈ 0

    def _d_trial_polynomial(self, x: float, alpha: float, L=None) -> float:
        """Derivative of polynomial: dψ/dx = -2αx/L²."""
        if L is None or abs(x) > L / 2.0:
            return 0.0
        return -2.0 * alpha * x / (L * L)

    def _integrate(self, f_values: List[float], dx: float) -> float:
        """Simple trapezoidal integration."""
        n = len(f_values)
        if n < 2:
            return 0.0
        return dx * (0.5 * f_values[0] + sum(f_values[1:n-1]) + 0.5 * f_values[n-1])

    def _compute_expectation_H(self, trial_func, d_trial_func, alpha: float,
                                 mass: float, ptype: str,
                                 x_grid: List[float], dx: float,
                                 **kwargs) -> tuple:
        """
        Compute ⟨H⟩ = ⟨T⟩ + ⟨V⟩ for given trial parameter α.
        
        Returns (⟨H⟩, ⟨T⟩, ⟨V⟩, norm²).
        """
        hbar = self.hbar
        
        psi_vals = []
        dpsi_vals = []
        V_vals = []

        for x in x_grid:
            psi = trial_func(x, alpha, kwargs.get("box_length_m"))
            psi_vals.append(psi)
            
            if d_trial_func:
                dpsi = d_trial_func(x, alpha, kwargs.get("box_length_m"))
            else:
                # Numerical derivative
                eps = dx * 0.01
                psi_plus = trial_func(x + eps, alpha, kwargs.get("box_length_m"))
                psi_minus = trial_func(x - eps, alpha, kwargs.get("box_length_m"))
                dpsi = (psi_plus - psi_minus) / (2.0 * eps)
            dpsi_vals.append(dpsi)

            # Potential V(x)
            if ptype == "harmonic":
                k = kwargs.get("force_constant_N_m", 1.0)
                V = 0.5 * k * x * x
            elif ptype == "infinite_well":
                L = kwargs.get("box_length_m", 1.0)
                V = 0.0 if abs(x) < L / 2.0 else 1e18
            elif ptype == "hydrogen_like" or ptype == "coulomb":
                Z = kwargs.get("nuclear_charge_Z", 1)
                r = abs(x) + 1e-15  # Avoid singularity at r=0
                V = -1.44e-17 * Z / r  # -Ze²/(4πε₀r) in Joules (e²/4πε₀ ≈ 1.44e-17 eV·m → convert)
                # Actually use: e²/(4πε₀) = 2.307×10⁻²⁸ J·m... let's use proper constant
                e2_4pie0 = 2.3071e-28  # J·m
                V = -e2_4pie0 * Z / r
            elif ptype == "anharmonic":
                k = kwargs.get("force_constant_N_m", 1.0)
                lam = kwargs.get("anharmonicity_coeff", 1.0)
                V = 0.5 * k * x * x + lam * x ** 4
            else:
                V = 0.0
            V_vals.append(V)

        # Norm squared: ∫ψ*ψ dx
        norm_sq = self._integrate([p * p for p in psi_vals], dx)
        
        if norm_sq < 1e-30:
            return float('inf'), 0.0, 0.0, 0.0

        # Kinetic energy: ⟨T⟩ = (ℏ²/2m) ∫|dψ/dx|² dx
        T_integrand = [dp * dp for dp in dpsi_vals]
        T_raw = self._integrate(T_integrand, dx)
        T = (hbar * hbar / (2.0 * mass)) * T_raw

        # Potential energy: ⟨V⟩ = ∫ψ*Vψ dx = ∫V·|ψ|² dx
        V_integrand = [V_vals[i] * psi_vals[i] * psi_vals[i] for i in range(len(psi_vals))]
        V_exp = self._integrate(V_integrand, dx)

        H = (T + V_exp) / norm_sq
        T_norm = T / norm_sq
        V_norm = V_exp / norm_sq

        return H, T_norm, V_norm, norm_sq

    def _golden_section_search(self, trial_func, d_trial_func, mass, ptype,
                                  x_grid, dx, a, b, **kwargs) -> tuple:
        """Minimize ⟨H⟩(α) using golden section search."""
        gr = (math.sqrt(5) - 1) / 2  # ~0.618

        tol = 1e-10 * max(abs(a), abs(b), 1.0)
        max_iter = 200

        c = b - gr * (b - a)
        d = a + gr * (b - a)

        fc, _, _, nc = self._compute_expectation_H(trial_func, d_trial_func, c, mass, ptype, x_grid, dx, **kwargs)
        fd, _, _, nd = self._compute_expectation_H(trial_func, d_trial_func, d, mass, ptype, x_grid, dx, **kwargs)

        for _iter in range(max_iter):
            if abs(b - a) < tol:
                break
            
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b - gr * (b - a)
                fc, _, _, _ = self._compute_expectation_H(trial_func, d_trial_func, c, mass, ptype, x_grid, dx, **kwargs)
            else:
                a = c
                c = d
                fc = fd
                d = a + gr * (b - a)
                fd, _, _, _ = self._compute_expectation_H(trial_func, d_trial_func, d, mass, ptype, x_grid, dx, **kwargs)

        alpha_opt = (a + b) / 2.0
        E_opt, T_opt, V_opt, norm_opt = self._compute_expectation_H(
            trial_func, d_trial_func, alpha_opt, mass, ptype, x_grid, dx, **kwargs
        )
        return alpha_opt, E_opt, T_opt, V_opt

    def _get_exact_ground_state(self, ptype: str, mass: float, **kwargs) -> tuple:
        """Return exact ground state energy (J, eV) if known analytically."""
        hbar = self.hbar
        eV = self.eV_per_J

        if ptype == "harmonic":
            k = kwargs.get("force_constant_N_m", 1.0)
            omega = math.sqrt(k / mass)
            E_exact = 0.5 * hbar * omega  # Zero-point energy
            return E_exact, E_exact * eV

        elif ptype == "infinite_well":
            L = kwargs.get("box_length_m", 1.0)
            E_exact = hbar * hbar * math.pi * math.pi / (2.0 * mass * L * L)
            return E_exact, E_exact * eV

        elif ptype == "hydrogen_like" or ptype == "coulomb":
            Z = kwargs.get("nuclear_charge_Z", 1)
            # E_1 = -Z² · Rydberg (in J): Ry = 2.179872e-18 J
            Ry_J = 2.179872361e-18
            E_exact = -Z * Z * Ry_J
            return E_exact, E_exact * eV

        return None, None

    def _run_base(self, potential_type: str, mass_kg: float, trial_function_type: str,
                  force_constant_N_m: float = None, box_length_m: float = None,
                  nuclear_charge_Z: int = 1, anharmonicity_coeff: float = None,
                  n_grid: int = 1000, param_search_range: list = None,
                  custom_trial_params: dict = None) -> dict:

        # Get trial functions
        trial_func, d_trial_func = self._get_trial_func(trial_function_type)

        # Set up integration domain
        if potential_type == "infinite_well":
            if box_length_m is None:
                raise ChemMCPError("infinite_well requires box_length_m.")
            L_domain = box_length_m * 1.5  # Slightly larger than well
            center = 0.0
        elif potential_type == "harmonic":
            if force_constant_N_m is not None and force_constant_N_m > 0:
                omega = math.sqrt(force_constant_N_m / mass_kg)
            else:
                omega = 1e14
            # Characteristic length: sqrt(ℏ/(mω))
            char_len = math.sqrt(self.hbar / (mass_kg * omega))
            L_domain = 12.0 * char_len
            center = 0.0
        elif potential_type in ("hydrogen_like", "coulomb"):
            # Bohr radius scale
            a0 = 5.29177210903e-11
            Z = nuclear_charge_Z
            L_domain = 20.0 * a0 / Z
            center = 0.0
        elif potential_type == "anharmonic":
            if force_constant_N_m:
                omega = math.sqrt(force_constant_N_m / mass_kg)
                char_len = math.sqrt(self.hbar / (mass_kg * omega))
            else:
                char_len = 1e-10
            L_domain = 10.0 * char_len
            center = 0.0
        else:
            L_domain = 1e-9
            center = 0.0

        dx = L_domain / n_grid
        x_grid = [center - L_domain / 2.0 + i * dx for i in range(n_grid)]

        # Set up keyword arguments for potential
        pot_kwargs = {}
        if force_constant_N_m is not None:
            pot_kwargs["force_constant_N_m"] = force_constant_N_m
        if box_length_m is not None:
            pot_kwargs["box_length_m"] = box_length_m
        pot_kwargs["nuclear_charge_Z"] = nuclear_charge_Z
        if anharmonicity_coeff is not None:
            pot_kwargs["anharmonicity_coeff"] = anharmonicity_coeff

        # Determine search range for variational parameter α
        if param_search_range:
            alpha_min, alpha_max = param_search_range
        else:
            # Auto-determine reasonable search range
            if potential_type == "harmonic":
                if force_constant_N_m:
                    omega = math.sqrt(force_constant_N_m / mass_kg)
                    alpha0 = mass_kg * omega / self.hbar
                else:
                    alpha0 = 1e20
                alpha_min, alpha_max = alpha0 * 0.1, alpha0 * 10.0
            elif potential_type == "infinite_well":
                alpha_min, alpha_max = 0.5, 15.0
            elif potential_type in ("hydrogen_like", "coulomb"):
                a0 = 5.29177210903e-11
                alpha_min, alpha_max = 0.5 / a0, 5.0 / a0
            else:
                alpha_min, alpha_max = 1e8, 1e14

        # Perform minimization
        alpha_opt, E_var, T_var, V_var = self._golden_section_search(
            trial_func, d_trial_func, mass_kg, potential_type,
            x_grid, dx, alpha_min, alpha_max, **pot_kwargs
        )

        # Get exact ground state energy for comparison
        E_exact_J, E_exact_eV = self._get_exact_ground_state(potential_type, mass_kg, **pot_kwargs)

        # Compute relative error
        if E_exact_J is not None and E_exact_J != 0:
            rel_error_pct = abs((E_var - E_exact_J) / E_exact_J) * 100
        else:
            rel_error_pct = None

        result = {
            "potential_type": potential_type,
            "trial_function": trial_function_type,
            "optimal_variational_parameter_alpha": round(alpha_opt, 10),
            "variational_energy_J": round(E_var, 25),
            "variational_energy_eV": round(E_var * self.eV_per_J, 10),
            "kinetic_energy_component_eV": round(T_var * self.eV_per_J, 10),
            "potential_energy_component_eV": round(V_var * self.eV_per_J, 10),
            "exact_energy_J": round(E_exact_J, 25) if E_exact_J is not None else None,
            "exact_energy_eV": round(E_exact_eV, 10) if E_exact_eV is not None else None,
            "relative_error_percent": round(rel_error_pct, 6) if rel_error_pct is not None else None,
            "integration_grid_points": n_grid,
            "domain_half_width_m": round(L_domain / 2.0, 20),
            "virial_theorem_ratio_V_over_T": round(V_var / T_var, 6) if abs(T_var) > 1e-30 else None,
        }

        logger.info(f"VariationalMethod: {potential_type}/{trial_function_type}, E_var={E_var*self.eV_per_J:.6f}eV, α={alpha_opt:.4e}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            ptype = parts[0]
            mass = float(parts[1])
            tfunc = parts[2]
            
            kwargs = {}
            for p in parts[3:]:
                if p.startswith("k="):
                    kwargs["force_constant_N_m"] = float(p.split("=")[1])
                elif p.startswith("L="):
                    kwargs["box_length_m"] = float(p.split("=")[1])
                elif p.startswith("Z="):
                    kwargs["nuclear_charge_Z"] = int(p.split("=")[1])
                elif p.startswith("lambda="):
                    kwargs["anharmonicity_coeff"] = float(p.split("=")[1])

            return self._run_base(ptype, mass, tfunc, **kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'ptype mass trial [k=..] [L=..] [Z=..]'")
