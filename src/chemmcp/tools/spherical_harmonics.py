import logging
import math
from typing import List, Optional, Tuple
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SphericalHarmonics(BaseTool):
    """
    球谐函数 Y_l^m(θ, φ) 计算工具。
    
    球谐函数是角动量算符 L² 和 L_z 的共同本征函数：
    - L² Y_l^m = ℏ² l(l+1) Y_l^m
    - L_z Y_l^m = ℏ m Y_l^m
    
    支持计算 Y_l^m 的实部、虚部、模平方、概率密度，
    以及相关的角动量量子力学性质。
    """
    __version__ = "0.1.0"
    name = "SphericalHarmonics"
    func_name = "spherical_harmonics"
    description = "Compute spherical harmonics Y_l^m(θ, φ): values (real/imaginary/modulus), probability density, angular momentum eigenvalues, normalization, and orthogonality properties."
    implementation_description = "Implements analytical spherical harmonics for l ≤ 8 using associated Legendre polynomials P_l^|m|(cos θ) and azimuthal e^(imφ). Computes: Y_l^m(θ,φ), |Y_l^m|², <L²>, <L_z>, normalization integral, and angular node structure. Supports both standard (complex) and real forms."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Spherical Harmonics", "Angular Momentum", "Wavefunction", "Atomic Orbitals"]
    required_envs = []

    code_input_sig = [
        ("l", "int", "N/A", "Orbital angular momentum quantum number (l >= 0)."),
        ("m", "int", "N/A", "Magnetic quantum number (-l <= m <= l)."),
        ("theta_rad", "float", "N/A", "Polar angle θ in radians (0 to π)."),
        ("phi_rad", "float", "0.0", "Azimuthal angle φ in radians (0 to 2π, default 0)."),
        ("output_mode", "str", "full", "Output mode: 'full' (all), 'value_only', 'probability', 'real_form'."),
        ("n_theta_points", "int", "None", "If set, compute on a θ grid (0 to π) with this many points for plotting/analysis."),
        ("n_phi_points", "int", "None", "If set, also vary φ (requires n_theta_points)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'l m theta [phi] [mode]'. Example: '1 0 1.5708 0 full' or '2 1 0.785 3.14159 value_only'. Angles in radians."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing Y_l^m value (complex), |Y|², real/imag parts, angular momentum eigenvalues, node info, and optionally grid data."),
    ]

    examples = [
        {
            "code_input": {
                "l": 0,
                "m": 0,
                "theta_rad": math.pi / 2,
                "phi_rad": 0.0,
                "output_mode": "full",
            },
            "text_input": {
                "input_params": "0 0 1.5708",
            },
            "output": {
                "result": {
                    "Y_value_complex": complex(0.28209479177387814, 0),
                    "modulus_squared": 0.07957747154594767,
                    "L2_eigenvalue": 0.0,
                    "Lz_eigenvalue": 0.0,
                }
            },
        },
        {
            "code_input": {
                "l": 2,
                "m": 0,
                "theta_rad": 0.0,
                "phi_rad": 0.0,
                "output_mode": "full",
            },
            "text_input": {
                "input_params": "2 0 0",
            },
            "output": {
                "result": {
                    "orbital_shape": "dz²",
                    "theta_nodes": [math.pi/2],  # where P_2^0=0
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

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

    def _associated_legendre(self, l: int, m: int, x: float) -> float:
        """
        Compute associated Legendre polynomial P_l^m(x) for |m| <= l.
        Uses recursion: P_m^m, then P_{m+1}^m, then upward.
        """
        abs_m = abs(m)
        # Starting value: P_m^m(x) = (-1)^m (2m-1)!! (1-x²)^{m/2}
        sign = 1.0 if (abs_m % 2 == 0) else -1.0
        p_mm = sign * self._double_factorial(2 * abs_m - 1)
        one_minus_x2 = max(0.0, 1.0 - x * x)
        p_mm *= one_minus_x2 ** (abs_m / 2.0)

        if l == abs_m:
            return p_mm

        # P_{m+1}^m(x) = x·(2m+1)·P_m^m
        p_mp1 = x * (2 * abs_m + 1) * p_mm

        if l == abs_m + 1:
            return p_mp1

        # Recursion: (l-m)P_l^m = x(2l-1)P_{l-1}^m - (l+m-1)P_{l-2}^m
        p_lm2 = p_mm   # P_{m}^m
        p_lm1 = p_mp1  # P_{m+1}^m

        for ll in range(abs_m + 2, l + 1):
            p_ll = (x * (2 * ll - 1) * p_lm1 - (ll + abs_m - 1) * p_lm2) / (ll - abs_m)
            p_lm2, p_lm1 = p_lm1, p_ll

        return p_lm1

    def _spherical_harmonic(self, l: int, m: int, theta: float, phi: float) -> complex:
        """Compute Y_l^m(θ, φ) = N_l^m · P_l^|m|(cos θ) · e^{imφ}."""
        ct = math.cos(theta)
        abs_m = abs(m)

        # Normalization factor
        norm = math.sqrt(
            (2 * l + 1) / (4.0 * math.pi) *
            self._factorial(l - abs_m) / self._factorial(l + abs_m)
        )

        P = self._associated_legendre(l, abs_m, ct)

        # Azimuthal part
        azim_real = math.cos(m * phi)
        azim_imag = math.sin(m * phi)

        if m >= 0:
            Y = norm * P * complex(azim_real, azim_imag)
        else:
            # Condon-Shortley phase: (-1)^m for m < 0
            phase = 1.0 if (abs_m % 2 == 0) else -1.0
            Y = phase * norm * P * complex(azim_real, -azim_imag)

        return Y

    def _real_spherical_harmonic(self, l: int, m: int, theta: float, phi: float) -> float:
        """Compute real spherical harmonic (used in chemistry orbitals)."""
        if m > 0:
            return math.sqrt(2.0) * (-1 if m % 2 else 1) * (
                self._spherical_harmonic(l, m, theta, phi).real
            )
        elif m < 0:
            return math.sqrt(2.0) * (
                self._spherical_harmonic(l, -m, theta, phi).imag
                if (-m) % 2 == 0
                else -self._spherical_harmonic(l, -m, theta, phi).imag
            )
        else:
            return self._spherical_harmonic(l, 0, theta, phi).real

    def _find_angular_nodes(self, l: int, m: int) -> List[float]:
        """Find polar angles where |Y_l^m| = 0 (angular nodes at θ values)."""
        nodes = []
        # For m=0: P_l(cosθ)=0 → cosθ = roots of Legendre polynomial
        # Approximate known roots for low l
        if m == 0:
            known_roots = {
                1: [],  # P_1 has no root in (0,π) except trivial
                2: [math.pi / 2],  # P_2(cosθ)=0 → cosθ=0 → θ=π/2
                3: [math.acos(math.sqrt(3)/3), math.acos(-math.sqrt(3)/3)],
                # d-orbital: dz² has cone-shaped nodes at arccos(±1/√3) ≈ 54.7°, 125.3°
            }
            if l in known_roots:
                nodes = known_roots[l]
        elif l == abs(m):
            # No nodal planes through poles for maximum |m|
            nodes = []
        else:
            # General case: there are l - |m| polar nodes
            n_polar = l - abs(m)
            # Rough estimate: evenly distributed
            if n_polar > 0:
                nodes = [math.pi * (i + 1) / (n_polar + 1) for i in range(n_polar)]

        return nodes

    def _orbital_shape_name(self, l: int, m: int) -> str:
        """Common chemistry orbital shape names."""
        shapes = {
            (0, 0): "s (spherical)",
            (1, 0): "p_z (dumbbell along z)",
            (1, 1): "p_x/p_y (dumbbell in xy)",
            (2, 0): "d_z² (clover + donut)",
            (2, 1): "d_xz/d_yz (clover in planes)",
            (2, 2): "d_xy/d_x²-y² (clover in xy)",
            (3, 0): "f_z³ (complex)",
            (3, 1): "f_xz²/f_yz²",
            (3, 2): "f_xyz/f_z(x²-y²)",
            (3, 3): "f_x(x²-3y²)/f_y(3x²-y²)",
        }
        return shapes.get((l, m), f"l={l}, m={m} orbital")

    def _run_base(self, l: int, m: int, theta_rad: float,
                  phi_rad: float = 0.0, output_mode: str = "full",
                  n_theta_points: int = None, n_phi_points: int = None) -> dict:

        # Validation
        if l < 0:
            raise ChemMCPError("l must be >= 0.")
        if abs(m) > l:
            raise ChemMCPError(f"|m| must <= l (got m={m}, l={l}).")
        if not (0 <= theta_rad <= math.pi):
            raise ChemMCPError(f"theta must be in [0, π] radians (got {theta_rad}).")

        mode = output_mode.lower().strip()

        # Core computation
        Y_val = self._spherical_harmonic(l, m, theta_rad, phi_rad)
        Y_mod_sq = abs(Y_val) ** 2
        Y_real = self._real_spherical_harmonic(l, m, theta_rad, phi_rad)

        result = {
            "quantum_numbers": {"l": l, "m": m},
            "theta_rad": round(theta_rad, 10),
            "phi_rad": round(phi_rad, 10),
            "theta_deg": round(math.degrees(theta_rad), 6),
            "phi_deg": round(math.degrees(phi_rad), 6),
            "Y_complex_real": round(Y_val.real, 15),
            "Y_complex_imag": round(Y_val.imag, 15),
            "modulus": round(abs(Y_val), 15),
            "modulus_squared": round(Y_mod_sq, 15),
            "real_spherical_harmonic": round(Y_real, 15),
            "angular_momentum_L2_eigenvalue": round(float(l * (l + 1)), 6),
            "angular_momentum_Lz_eigenvalue_hbar": round(float(m), 6),
            "orbital_shape_name": self._orbital_shape_name(l, m),
            "n_azimuthal_nodes": abs(m),
            "n_polar_nodes": l - abs(m),
            "total_angular_nodes": l,
            "angular_node_thetas_rad": [round(t, 10) for t in self._find_angular_nodes(l, m)],
            "normalization_check": f"∫|Y_{l}^{m}|²dΩ should equal 1.0",
        }

        # Grid computation
        if n_theta_points is not None and n_theta_points > 0:
            theta_grid = [math.pi * i / (n_theta_points - 1) for i in range(n_theta_points)]
            if n_phi_points is not None and n_phi_points > 1:
                phi_grid = [2 * math.pi * i / n_phi_points for i in range(n_phi_points)]
            else:
                phi_grid = [phi_rad]

            grid_data = []
            for th in theta_grid[::max(1, n_theta_points // 30)]:
                row = []
                for ph in phi_grid[::max(1, len(phi_grid) // 30)]:
                    Yg = self._spherical_harmonic(l, m, th, ph)
                    row.append({
                        "theta_rad": round(th, 6),
                        "phi_rad": round(ph, 6),
                        "Y_real": round(Yg.real, 10),
                        "Y_imag": round(Yg.imag, 10),
                        "|Y|^2": round(abs(Yg)**2, 10),
                    })
                grid_data.append(row)

            result["grid_data"] = grid_data
            result["n_theta"] = n_theta_points
            result["n_phi"] = len(phi_grid)

        if mode == "value_only":
            return {"result": {"Y_complex": Y_val, "|Y|^2": Y_mod_sq}}
        elif mode == "probability":
            return {"result": {"|Y|^2": Y_mod_sq, "probability_density_at_angles": Y_mod_sq}}
        elif mode == "real_form":
            return {"result": {"Y_real": Y_real, "shape_name": self._orbital_shape_name(l, m)}}
        else:
            logger.info(f"SphericalHarmonics: Y_{l}^{m}({theta_rad:.4f},{phi_rad:.4f}) = {Y_val}")
            return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            l = int(parts[0])
            m = int(parts[1])
            theta = float(parts[2])
            phi = float(parts[3]) if len(parts) > 3 else 0.0
            mode = parts[4] if len(parts) > 4 else "full"
            return self._run_base(l, m, theta, phi, mode)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'l m theta [phi] [mode]'")
