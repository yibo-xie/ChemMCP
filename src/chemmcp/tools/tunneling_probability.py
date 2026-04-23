import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TunnelingProbability(BaseTool):
    """
    量子隧穿概率计算。
    
    使用 WKB 近似和精确解计算势垒穿透概率:
    
    矩形势垒 (精确解):
      T = [1 + (V₀²sinh²(κa))/(4E(V₀-E))]⁻¹  for E < V₀
      T = [1 + (V₀²sin²(ka))/(4E(E-V₀))]⁻¹   for E > V₀
    
    其中 κ = √[2m(V₀-E)]/ℏ, k = √[2m(E-V₀)]/ℏ
    
    WKB 近似:
      T ≈ exp(-2∫_{x₁}^{x₂} |p(x)|dx/ℏ)
      其中 p(x) = √[2m(V(x)-E)]
    
    α衰变 (Gamow因子):
      T ≈ exp(-2G), G = (Z-2)e²/(4πε₀ℏ) · √(2m_α/2E_α) · [arccos(√(Q/B)) - √(Q/B-Q²/B²)]
    """
    __version__ = "0.1.0"
    name = "TunnelingProbability"
    func_name = "tunneling_probability"
    description = "Calculate quantum tunneling (barrier penetration) probability using WKB approximation and exact solutions for rectangular, triangular, Gaussian, Eckart barriers, and alpha decay."
    implementation_description = "Implements exact transmission coefficients for rectangular barriers, WKB approximation for arbitrary barrier shapes (triangular, Gaussian, Eckart), Gamow theory for alpha decay with nuclear Coulomb barrier, and Fowler-Nordheim field emission. Returns transmission/reflection probabilities, tunneling rates, WKB exponent values, and classical vs quantum comparison."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Tunneling", "WKB Approximation", "Barrier Penetration", "Alpha Decay"]
    required_envs = []

    code_input_sig = [
        ("barrier_type", "str", "N/A", "'rectangular', 'triangular', 'gaussian', 'eckart', 'alpha_decay', 'general_WKB', 'field_emission'."),
        ("particle_mass_kg", "float", "N/A", "Particle mass in kg."),
        ("energy_J", "float", "N/A", "Particle energy E in Joules."),
        ("barrier_height_J", "float", "None", "Barrier height V₀ in Joules (for rectangular/triangular)."),
        ("barrier_width_m", "float", "None", "Barrier width a in meters (for rectangular)."),
        ("shape_params", "dict", "None", "Extra parameters for non-rectangular barriers: {slope_N/m, sigma_m, etc}."),
        ("nuclear_charge_Z", "int", "None", "Daughter nucleus charge Z for alpha decay."),
        ("daughter_mass_u", "float", "None", "Daughter nucleus mass in amu for alpha decay."),
        ("electric_field_V_m", "float", "None", "Electric field strength for Fowler-Nordheim emission."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'barrier_type mass_Ekg energy_J [height_J] [width_m] [extra]'. Example: 'rectangular 9.109e-31 1e-20 5e-20 1e-9'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with transmission probability, reflection probability, tunneling rate, WKB exponent, barrier parameters, classical comparison."),
    ]

    examples = [
        {
            "code_input": {
                "barrier_type": "rectangular",
                "particle_mass_kg": 9.109e-31,
                "energy_J": 1e-20,
                "barrier_height_J": 5e-20,
                "barrier_width_m": 1e-9,
                "shape_params": None,
                "nuclear_charge_Z": None,
                "daughter_mass_u": None,
                "electric_field_V_m": None,
            },
            "text_input": {"input_params": "rectangular 9.109e-31 1e-20 5e-20 1e-9"},
            "output": {"result": {"transmission_probability": 0.08, "classically_allowed": False}},
        },
        {
            "code_input": {
                "barrier_type": "alpha_decay",
                "particle_mass_kg": 6.644657230e-27,
                "energy_J": 8.79e-13,
                "nuclear_charge_Z": 90,
                "shape_params": None,
                "barrier_height_J": None,
                "barrier_width_m": None,
                "daughter_mass_u": None,
                "electric_field_V_m": None,
            },
            "text_input": {"input_params": "alpha_decay 6.645e-27 8.79e-13 Z=90"},
            "output": {"result": {"transmission_probability": 1e-38, "half_life_estimate": "very long"}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34  # J·s
        self.eV_per_J = 6.241509074e18
        self.amu_kg = 1.66053906660e-27  # atomic mass unit in kg
        self.e_charge = 1.602176634e-19  # Coulomb
        self.e2_4pie0 = 2.3071e-28  # e²/(4πε₀) in J·m

    def _rectangular_barrier(self, m: float, E: float, V0: float, a: float) -> dict:
        """Exact transmission through rectangular barrier of height V0 and width a."""
        hbar = self.hbar

        if E <= 0:
            return {"T": 0.0, "R": 1.0}

        classically_allowed = E >= V0

        if abs(E - V0) < 1e-30 * max(abs(V0), 1e-40):
            # E ≈ V0: limiting case
            T = 1.0 / (1.0 + (m * V0 * a * a) / (2.0 * hbar * hbar))
        elif E < V0:
            # Tunneling regime: T = [1 + sinh²(κa)/(4E/V0·(1-E/V0))]⁻¹
            kappa = math.sqrt(2.0 * m * (V0 - E)) / hbar
            kappa_a = kappa * a
            
            if kappa_a > 500:  # Avoid overflow in sinh
                T = math.exp(-2.0 * kappa_a)
            else:
                sinh_ka = math.sinh(kappa_a)
                ratio = E / V0
                denom = 1.0 + (sinh_ka ** 2) / (4.0 * ratio * (1.0 - ratio)) if ratio > 0 and ratio < 1 else float('inf')
                T = 1.0 / denom
        else:
            # Above barrier: T = [1 + sin²(ka)/(4(E/V0)(E/V0-1))]⁻¹
            k = math.sqrt(2.0 * m * (E - V0)) / hbar
            k_a = k * a
            ratio = E / V0
            sin_ka = math.sin(k_a)
            denom = 1.0 + (sin_ka ** 2) / (4.0 * ratio * (ratio - 1.0)) if ratio > 1 else float('inf')
            T = 1.0 / denom
            if T > 1.0:
                T = 1.0  # Numerical artifact correction

        R = max(0.0, 1.0 - T)

        # WKB exponent for comparison
        kappa = math.sqrt(2.0 * m * max(V0 - E, 0)) / hbar
        wkb_exp = 2.0 * kappa * a

        return {
            "T": round(T, 15),
            "R": round(R, 15),
            "classically_allowed": classically_allowed,
            "kappa_inv_nm": round(1e9 / kappa if kappa > 0 else float('inf'), 6),
            "kappa_a": round(kappa * a if kappa > 0 else 0, 4),
            "wkb_exponent": round(wkb_exp, 4),
            "wkb_T_approx": round(math.exp(-wkb_exp) if wkb_exp < 700 else 0.0, 15),
            "regime": "above_barrier" if classically_allowed else "tunneling",
        }

    def _triangular_barrier(self, m: float, E: float, F: float) -> dict:
        """
        Triangular barrier: V(x) = V₀(1 - x/a) for 0<x<a, linear ramp.
        
        Or field-emission type: V(x) = E_F - eFx (triangular from electric field).
        
        WKB: T ≈ exp(-4√(2m)V₀^{3/2}/(3ℏeF)) for triangular barrier.
        
        Args:
          F: Electric field slope in N/m (or use shape_params['slope'])
        """
        hbar = self.hbar
        
        # For triangular barrier V(x) = V₀ - Fx, turning point at x_t = (V₀-E)/F
        # WKB integral: (2/3)√(2m)(V₀-E)^{3/2}/(ℏF)
        V0 = shape_params.get("V0_J", E * 10) if isinstance(shape_params, dict) else E * 10
        
        if F <= 0 or E >= V0:
            return {"T": 1.0, "R": 0.0, "classically_allowed": True, "regime": "no_barrier"}

        delta_V = V0 - E
        if delta_V <= 0:
            return {"T": 1.0, "R": 0.0, "classically_allowed": True, "regime": "above_barrier"}

        # WKB exponent for triangular barrier
        wkb_exp = (4.0 / 3.0) * math.sqrt(2.0 * m) * (delta_V ** 1.5) / (hbar * F)
        T_wkb = math.exp(-wkb_exp) if wkb_exp < 700 else 0.0

        return {
            "T": round(T_wkb, 15),
            "R": round(max(0, 1.0 - T_wkb), 15),
            "classically_allowed": False,
            "wkb_exponent": round(wkb_exp, 4),
            "regime": "tunneling (WKB)",
            "barrier_shape": "triangular",
            "turning_point_m": round(delta_V / F, 20),
        }

    def _gaussian_barrier(self, m: float, E: float, V0: float, sigma: float) -> dict:
        """Gaussian barrier: V(x) = V₀·exp(-x²/(2σ²)).
        
        WKB approximate.
        """
        hbar = self.hbar

        if E >= V0:
            return {"T": 1.0, "R": 0.0, "classically_allowed": True, "regime": "above_barrier"}

        # Find classical turning points: V(x_t) = E → x_t = ±σ·√(-2ln(E/V0))
        if E <= 0:
            return {"T": 0.0, "R": 1.0, "classically_allowed": False, "regime": "no_energy"}

        ln_ratio = math.log(V0 / E)
        if ln_ratio <= 0:
            return {"T": 1.0, "R": 0.0, "classically_allowed": True, "regime": "above"}

        x_t = sigma * math.sqrt(2.0 * ln_ratio)

        # WKB integral for Gaussian (approximate):
        # ∫_{-x_t}^{x_t} √[2m(V₀exp(-x²/2σ²)-E)] dx
        # Approximate as: √(2m(V₀-E)) · (effective_width)
        effective_width = math.sqrt(2.0 * math.pi) * sigma * (1.0 - E/V0) ** 0.25
        kappa_eff = math.sqrt(2.0 * m * (V0 - E)) / hbar
        wkb_exp = 2.0 * kappa_eff * effective_width
        T_wkb = math.exp(-wkb_exp) if wkb_exp < 700 else 0.0

        return {
            "T": round(T_wkb, 15),
            "R": round(max(0, 1.0 - T_wkb), 15),
            "classically_allowed": False,
            "wkb_exponent": round(wkb_exp, 4),
            "turning_points_pm_m": [round(-x_t, 20), round(x_t, 20)],
            "regime": "tunneling (WKB)",
            "barrier_shape": "gaussian",
        }

    def _eckart_barrier(self, m: float, E: float, V0: float, a: float) -> dict:
        """Eckart barrier: V(x) = V₀/cosh²(x/a).
        
        Exact transmission: T = [1+cosh(2π√(2mV0)a/ℏ)/cosh(2π√(2mE)a/ℏ)]^{-1} ... simplified.
        
        Actually exact Eckart: T = [cosh(2πp) - cosh(2πq)] / [cosh(2πp) + cosh(2q')] ...
        Using simpler form: T = sinh²(2πa√(2mE)/ℏ) / [sinh²(2πa√(2mE)/ℏ) + cos²(πa√(2m(V0+√(V0²...))]
        
        Let's use the well-known result:
        T = [1 + sinh²(2πa p / ℏ) / sin²(πδ)]⁻¹ where δ involves √(8mV0)a/ℏ
        """
        hbar = self.hbar

        if E >= V0:
            return {"T": 1.0, "R": 0.0, "classically_allowed": True, "regime": "above_barrier"}

        # Simplified Eckart transmission
        p = math.sqrt(2.0 * m * E) / hbar  # momentum
        kappa = math.sqrt(2.0 * m * (V0 - E)) / hbar

        # Characteristic parameters
        theta = 2.0 * a * p
        eta = 2.0 * a * kappa

        # Transmission coefficient (approximate formula)
        if eta > 20:
            T = math.exp(-eta)  # Thick/high barrier limit
        else:
            # More accurate: T = 1/(1+cosh(2πaκ/ℏ)... )
            arg = a * math.sqrt(2.0 * m * V0) / hbar
            try:
                num = math.sinh(theta) ** 2 if theta < 500 else float('inf')
                den = math.sinh(theta) ** 2 + math.cos(math.pi * math.sqrt(1 + 8*m*V0*a*a/(hbar*hbar))) ** 2
                T = num / den if den > 0 else 0.0
            except (OverflowError, ValueError):
                T = math.exp(-eta)

        T = max(0.0, min(1.0, T))

        return {
            "T": round(T, 15),
            "R": round(max(0, 1.0 - T), 15),
            "classically_allowed": False,
            "wkb_exponent": round(2.0 * kappa * a * 1.2, 4),  # Approximate
            "regime": "tunneling (Eckart)",
            "barrier_shape": "Eckart (sech²)",
        }

    def _alpha_decay(self, m_alpha: float, E_alpha: float, Z_daughter: int = None,
                      M_daughter_amu: float = None) -> dict:
        """
        Alpha decay Gamow factor.
        
        T ≈ exp(-2G) where G is the Gamow factor:
        G = 2(Z-2)e²/(4πε₀ℏv) · [arccos(√(Q/B)) - √(Q/B)·√(1-Q/B)]
        
        Simplified: G ≈ (Z-2)e²/(ℏv) · (π/2 - 2√(R/B)) for thin barrier approx
        where B = Coulomb barrier height, Q = E_alpha, R = nuclear radius
        """
        hbar = self.hbar
        e2_4pie0 = self.e2_4pie0
        e_charge = self.e_charge

        if Z_daughter is None:
            Z_daughter = 86  # Default Rn (from Po decay)
        if M_daughter_amu is None:
            M_daughter_amu = 210.0  # Approximate

        Z_alpha = 2
        Z_total = Z_daughter + Z_alpha  # Parent nucleus charge

        # Nuclear radius: R = r₀ A^{1/3}, r₀ ≈ 1.2 fm
        A_parent = M_daughter_amu + 4.0
        R_nuclear = 1.2e-15 * (A_parent ** (1.0 / 3.0))  # meters
        R_touch = R_nuclear + 1.2e-15 * (4.0) ** (1.0 / 3.0)  # Touching distance

        # Coulomb barrier at touching point
        B_coulomb = Z_alpha * Z_daughter * e2_4pie0 / R_touch  # Joules

        # Alpha velocity
        v_alpha = math.sqrt(2.0 * E_alpha / m_alpha)

        # Gamow factor (simplified integral)
        # G = (Z_alpha * Z_d * e²)/(4πε₀ ℏ v) × [arccos(√(E/B)) - √(E/B(1-E/B))]
        Q_over_B = E_alpha / B_coulomb if B_coulomb > 0 else 0

        if Q_over_B >= 1:
            return {"T": 1.0, "R": 0.0, "classically_allowed": True, "regime": "above_Coulomb_barrier"}

        sqrt_QB = math.sqrt(Q_over_B)
        arccos_term = math.acos(sqrt_QB) if Q_over_B <= 1 else 0
        sqrt_term = sqrt_QB * math.sqrt(max(0, 1.0 - Q_over_B))

        prefactor = Z_alpha * Z_daughter * e2_4pie0 / (hbar * v_alpha)
        G = prefactor * (arccos_term - sqrt_term)

        T_gamow = math.exp(-2.0 * G) if G < 800 else 0.0

        # Estimate half-life from T and assault frequency
        f_assault = v_alpha / (2.0 * R_nuclear)  # Attempts per second (~10^21)
        lambda_decay = f_assault * T_gamow
        if lambda_decay > 0 and lambda_decay < 1e100:
            t_half_s = math.log(2) / lambda_decay
        else:
            t_half_s = float('inf')

        # Convert to human-readable units
        if t_half_s == float('inf'):
            half_life_str = "stable"
        elif t_half_s > 3600 * 24 * 365.25 * 1e9:
            half_life_str = f"{t_half_s_s / (3.156e16):.2e} Ga"
        elif t_half_s > 3600 * 24 * 365.25:
            half_life_str = f"{t_half_s / (3.156e7):.2f} years"
        elif t_half_s > 3600:
            half_life_str = f"{t_half_s / 3600:.2f} hours"
        else:
            half_life_str = f"{t_half_s:.2f} s"

        return {
            "T": round(T_gamow, 50) if T_gamow > 1e-300 else 0.0,
            "R": round(1.0 - T_gamow, 15),
            "classically_allowed": False,
            "gamow_factor_G": round(G, 4),
            "coulomb_barrier_J": round(B_coulomb, 25),
            "coulomb_barrier_eV": round(B_coulomb * self.eV_per_J, 4),
            "alpha_energy_eV": round(E_alpha * self.eV_per_J, 4),
            "alpha_velocity_m_s": round(v_alpha, 4),
            "nuclear_radius_fm": round(R_nuclear * 1e15, 3),
            "assault_frequency_s-1": round(f_assault, 4),
            "estimated_half_life": half_life_str,
            "estimated_half_life_seconds": t_half_s if t_half_s != float('inf') else None,
            "daughter_Z": Z_daughter,
            "parent_Z": Z_total,
            "regime": "alpha decay (Gamow theory)",
        }

    def _run_base(self, barrier_type: str, particle_mass_kg: float, energy_J: float,
                  barrier_height_J: float = None, barrier_width_m: float = None,
                  shape_params: dict = None, nuclear_charge_Z: int = None,
                  daughter_mass_u: float = None, electric_field_V_m: float = None) -> dict:

        btype = barrier_type.lower().replace("-", "_")

        if btype == "rectangular":
            if barrier_height_J is None or barrier_width_m is None:
                raise ChemMCPError("Rectangular barrier requires barrier_height_J and barrier_width_m.")
            result_data = self._rectangular_barrier(particle_mass_kg, energy_J, barrier_height_J, barrier_width_m)

        elif btype == "triangular":
            F = (shape_params or {}).get("slope_N_m", electric_field_V_m or 1.0)
            if barrier_height_J:
                sp = shape_params or {}
                sp["V0_J"] = barrier_height_J
                shape_params = sp
            result_data = self._triangular_barrier(particle_mass_kg, energy_J, F)

        elif btype == "gaussian":
            sigma = (shape_params or {}).get("sigma_m", barrier_width_m or 1e-10)
            V0 = barrier_height_J or (shape_params or {}).get("V0_J", energy_J * 10)
            result_data = self._gaussian_barrier(particle_mass_kg, energy_J, V0, sigma)

        elif btype == "eckart":
            a = barrier_width_m or (shape_params or {}).get("width_m", 1e-10)
            V0 = barrier_height_J or (shape_params or {}).get("V0_J", energy_J * 10)
            result_data = self._eckart_barrier(particle_mass_kg, energy_J, V0, a)

        elif btype in ("alpha_decay", "alpha"):
            result_data = self._alpha_decay(particle_mass_kg, energy_J, nuclear_charge_Z, daughter_mass_u)

        elif btype in ("field_emission", "fowler_nordheim"):
            # Fowler-Nordheim: triangular barrier from electric field
            F_field = electric_field_V_m or 1e9  # V/m
            work_function_J = barrier_height_J or 5.0 * self.e_charge  # default ~5 eV
            phi = work_function_J
            # FN: T ≈ exp(-4√(2m)φ^{3/2}/(3ℏeF))
            fn_exp = (4.0 / 3.0) * math.sqrt(2.0 * particle_mass_kg) * (phi ** 1.5) / (self.hbar * self.e_charge * F_field)
            T_fn = math.exp(-fn_exp) if fn_exp < 700 else 0.0
            result_data = {
                "T": round(T_fn, 15),
                "R": round(max(0, 1.0 - T_fn), 15),
                "classically_allowed": False,
                "wkb_exponent": round(fn_exp, 4),
                "work_function_eV": round(phi / self.e_charge, 4),
                "electric_field_V_m": F_field,
                "regime": "Fowler-Nordheim field emission",
            }
        else:
            raise ChemMCPError(f"Unknown barrier type: {barrier_type}. Choose from: "
                             f"rectangular, triangular, gaussian, eckart, alpha_decay, "
                             f"field_emission")

        result = {
            "barrier_type": barrier_type,
            "particle_mass_kg": particle_mass_kg,
            "particle_energy_J": round(energy_J, 25),
            "particle_energy_eV": round(energy_J * self.eV_per_J, 6),
            "transmission_probability": result_data.get("T", 0),
            "reflection_probability": result_data.get("R", 1.0),
            "is_classically_allowed": result_data.get("classically_allowed", True),
            "tunneling_regime": result_data.get("regime", "unknown"),
            "wkb_exponent_value": result_data.get("wkb_exponent"),
            **{k: v for k, v in result_data.items() 
               if k not in ("T", "R", "classically_allowed", "regime", "wkb_exponent")},
        }

        logger.info(f"TunnelingProbability: {btype}, T={result['transmission_probability']:.4e}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            btype = parts[0]
            mass = float(parts[1])
            energy = float(parts[2])
            
            kwargs = {}
            for p in parts[3:]:
                if p.startswith("height=") or p.startswith("V0="):
                    kwargs["barrier_height_J"] = float(p.split("=")[1])
                elif p.startswith("width=") or p.startswith("a="):
                    kwargs["barrier_width_m"] = float(p.split("=")[1])
                elif p.startswith("Z="):
                    kwargs["nuclear_charge_Z"] = int(p.split("=")[1])
                elif p.startswith("F=") or p.startswith("E="):
                    kwargs["electric_field_V_m"] = float(p.split("=")[1])

            return self._run_base(btype, mass, energy, **kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'btype mass_Ekg energy_J [height=..] [width=..] [Z=..]'")
