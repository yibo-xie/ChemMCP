import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SpectralLinewidth(BaseTool):
    """
    谱线宽度计算工具 — 自然展宽、碰撞（压力）展宽、多普勒展宽。
    
    三种机制独立贡献，总展宽可通过 Voigt 线型近似。
    """
    __version__                 = "0.1.0"
    name                        = "SpectralLinewidth"
    func_name                   = "calculate_linewidth"
    description                 = "Calculate spectral line width contributions: natural (lifetime), collisional (pressure), and Doppler (thermal) broadening."
    implementation_description  = "Natural: Δλ_nat = λ²/(2πcτ). Collisional: Δλ_coll ∝ P·σ·√(T/M). Doppler: Δλ_D = λ₀·√(8kT ln2 / (mc²)). Returns individual and total widths in wavelength & frequency units."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Spectroscopy", "Line Broadening", "Doppler", "Collisional", "Natural Width"]
    required_envs               = []

    code_input_sig   = [
        ("temperature_k",            "float",  "N/A",     "Temperature in Kelvin."),
        ("wavelength_m",             "float",  "N/A",     "Center wavelength in meters."),
        ("atomic_mass_kg",           "float",  "N/A",     "Atomic/molecular mass in kg."),
        ("natural_lifetime_s",       "float",  "None",    "Natural radiative lifetime in seconds (None → skip natural width)."),
        ("pressure_atm",             "float",  "None",    "Pressure in atmospheres (None → skip collisional width)."),
        ("collision_cross_section_m2","float", "1e-19",   "Collision cross-section in m²."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'T(K) lambda(m) mass_kg [tau_s] [P_atm] [sigma_m2]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with natural/collisional/doppler/total widths in Hz, m, nm."),
    ]

    examples         = [
        {
            "code_input": {
                "temperature_k":              300.0,
                "wavelength_m":               589.0e-9,
                "atomic_mass_kg":             23.0 * 1.66053906660e-27,  # Na atom
                "natural_lifetime_s":         16.3e-9,
                "pressure_atm":               1.0,
                "collision_cross_section_m2": 1e-19,
            },
            "text_input": {
                "input_params":               "300.0 589e-9 3.81924e-26 16.3e-9 1.0 1e-19",
            },
            "output": {
                "result": {
                    "temperature_K": 300.0,
                    "wavelength_nm": 589.0,
                    "doppler_width_Hz": 1506.7,
                    "doppler_width_m": 1.76e-12,
                    "natural_width_Hz": 9750000.0,
                    "total_width_approx_Hz": 9750015.07,
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """物理常数"""
        self.k_B = 1.380649e-23      # J/K
        self.h   = 6.62607015e-34     # J·s
        self.c   = 2.99792458e8       # m/s
        self.NA  = 6.02214076e23      # mol⁻¹
        self.kB_atm_L = 8.205736e-5   # atm·L/(mol·K)  (= R in these units)
        # R = k_B * N_A = 8.314 J/(mol·K)

    def _natural_width(self, wavelength_m: float, tau_s: float) -> dict:
        """自然展宽（洛伦兹线型）— 由海森堡测不准原理决定"""
        if tau_s is None or tau_s <= 0:
            return None
        
        # ΔE · τ ≈ ℏ/2  =>  Δν = 1/(2πτ),  Δλ = λ²/(2πcτ)
        delta_nu_Hz = 1.0 / (2.0 * math.pi * tau_s)
        delta_lambda_m = (wavelength_m ** 2) / (2.0 * math.pi * self.c * tau_s)
        
        return {
            "delta_nu_Hz":     delta_nu_Hz,
            "delta_lambda_m":  delta_lambda_m,
            "delta_lambda_nm": delta_lambda_m * 1e9,
            "FWHM_Hz":         delta_nu_Hz,          # Lorentzian FWHM = Δν
            "FWHM_m":          delta_lambda_m,
            "FWHM_nm":         delta_lambda_m * 1e9,
            "profile":         "Lorentzian",
        }

    def _doppler_width(self, temperature_k: float, wavelength_m: float, atomic_mass_kg: float) -> dict:
        """多普勒展宽（高斯线型）— 热运动导致"""
        T = temperature_k
        lam = wavelength_m
        m  = atomic_mass_kg
        
        if T <= 0 or m <= 0 or lam <= 0:
            raise ChemMCPError("T, λ, and mass must be positive.")
        
        # Doppler FWHM: Δν_D = ν₀ · √(8kT ln2 / (mc²))
        #            Δλ_D = λ₀ · √(8kT ln2 / (mc²))
        nu_0 = self.c / lam
        c_squared = self.c * self.c
        
        arg = (8.0 * self.k_B * T * math.log(2)) / (m * c_squared)
        if arg < 0:
            raise ChemMCPError("Invalid argument for sqrt in Doppler calculation.")
        
        sqrt_arg = math.sqrt(arg)
        
        delta_nu_Hz = nu_0 * sqrt_arg
        delta_lambda_m = lam * sqrt_arg
        
        return {
            "delta_nu_Hz":     delta_nu_Hz,
            "delta_lambda_m":  delta_lambda_m,
            "delta_lambda_nm": delta_lambda_m * 1e9,
            "FWHM_Hz":         delta_nu_Hz,
            "FWHM_m":          delta_lambda_m,
            "FWHM_nm":         delta_lambda_m * 1e9,
            "most_probable_speed_m_s": math.sqrt(2.0 * self.k_B * T / m),
            "profile":         "Gaussian",
        }

    def _collisional_width(
        self,
        temperature_k: float,
        pressure_atm: float,
        atomic_mass_kg: float,
        collision_cross_section_m2: float,
    ) -> dict:
        """碰撞（压力）展宽 — 洛伦兹线型"""
        if pressure_atm is None or pressure_atm <= 0:
            return None
        
        T = temperature_k
        P  = pressure_atm
        m  = atomic_mass_kg
        sigma = collision_cross_section_m2
        
        # Average relative speed: v_rel = √(8kT / (πμ)) ≈ √(8kT/(πm)) for same species
        v_rel = math.sqrt(8.0 * self.k_B * T / (math.pi * m))
        
        # Number density from ideal gas: n = P/(k_B·T)  [using SI: convert atm to Pa]
        # 1 atm = 101325 Pa
        P_Pa = pressure_atm * 101325.0
        n_density = P_Pa / (self.k_B * T)
        
        # Collision rate: γ_coll = n · σ · v_rel  (in s⁻¹, as FWHM in angular freq would be this)
        # For FWHM in Hz (not angular): γ_Hz = n · σ · v_rel / (2π)
        gamma_Hz = n_density * sigma * v_rel / (2.0 * math.pi)
        
        # Convert to wavelength width at a reference (approximate, needs central wavelength)
        # We'll return frequency-based; caller can convert with context
        delta_lambda_m_approx = None  # needs central wavelength
        
        return {
            "gamma_collision_rate_per_s": n_density * sigma * v_rel,
            "FWHM_Hz":         gamma_Hz,
            "FWHM_m":          None,  # needs wavelength context
            "pressure_atm":    P,
            "cross_section_m2": sigma,
            "mean_relative_speed_m_s": round(v_rel, 4),
            "number_density_m-3": round(n_density, 4),
            "profile":         "Lorentzian",
        }

    def _run_base(
        self,
        temperature_k: float,
        wavelength_m: float,
        atomic_mass_kg: float,
        natural_lifetime_s: Optional[float] = None,
        pressure_atm: Optional[float] = None,
        collision_cross_section_m2: Optional[float] = 1e-19,
    ) -> dict:
        """Core logic: calculate all three broadening mechanisms."""
        T  = float(temperature_k)
        lam = float(wavelength_m)
        m  = float(atomic_mass_kg)
        
        if T <= 0:
            raise ChemMCPError("Temperature must be > 0 K.")
        if lam <= 0:
            raise ChemMCPError("Wavelength must be > 0.")
        if m <= 0:
            raise ChemMCPError("Mass must be > 0.")

        sigma = float(collision_cross_section_m2) if collision_cross_section_m2 is not None else 1e-19

        # --- Calculate each contribution ---
        doppler = self._doppler_width(T, lam, m)
        natural = self._natural_width(lam, natural_lifetime_s) if natural_lifetime_s is not None else None
        collisional = self._collisional_width(T, pressure_atm, m, sigma) if pressure_atm is not None else None

        # --- Total width approximation ---
        # For Voigt profile, approximate total FWHM (Olivero algorithm approximation):
        # f_V ≈ f_L/2 + √(f_L²/4 + f_G²)
        fG = doppler["FWHM_Hz"] if doppler else 0.0
        fL = 0.0
        if natural:
            fL += natural["FWHM_Hz"]
        if collisional:
            fL += collisional["FWHM_Hz"]
        
        if fL > 0:
            total_fwhm_Hz = 0.5 * fL + math.sqrt(0.25 * fL * fL + fG * fG)
        else:
            total_fwhm_Hz = fG

        result = {
            "temperature_K": T,
            "wavelength_m": lam,
            "wavelength_nm": lam * 1e9,
            "atomic_mass_kg": m,
            "doppler": doppler,
            "natural": natural,
            "collisional": collisional,
            "total_FWHM_Hz": round(total_fwhm_Hz, 4),
            "dominant_broadening": self._identify_dominant(natural, collisional, doppler),
            "notes": "Total width uses Voigt approximation: fV ≈ fL/2 + √(fL²/4 + fG²)",
        }
        
        logger.info(f"Linewidth: T={T}K, λ={lam*1e9}nm, total_FWHM={total_fwhm_Hz:.2f}Hz")
        return result

    def _identify_dominant(self, natural, collisional, doppler) -> str:
        widths = []
        if natural and "FWMHz" not in str(natural):
            widths.append(("natural", natural["FWMHz"] if "FWMHz" in natural else natural.get("FWHM_Hz", 0)))
        if collisional:
            widths.append(("collisional", collisional.get("FWHM_Hz", 0)))
        if doppler:
            widths.append(("doppler", doppler.get("FWHM_Hz", 0)))
        
        if not widths:
            return "N/A"
        
        dominant = max(widths, key=lambda x: x[1])
        return f"{dominant[0]} ({dominant[1]:.2e} Hz)"

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need at least 'T lambda mass' params.")
            
            def pv(s, i):
                if i >= len(s): return None
                if s[i].lower() == "none": return None
                return float(s[i])
            
            T = float(parts[0])
            lam = float(parts[1])
            mass = float(parts[2])
            tau = pv(parts, 3)
            P = pv(parts, 4)
            sigma = pv(parts, 5) if pv(parts, 5) is not None else 1e-19
            
            return self._run_base(T, lam, mass, tau, P, sigma)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T(K) lambda(m) mass_kg [tau_s] [P_atm] [sigma_m2]'")
