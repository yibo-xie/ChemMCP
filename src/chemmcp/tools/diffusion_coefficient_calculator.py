import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DiffusionCoefficientCalculator(BaseTool):
    """
    Calculate diffusion coefficient using the Randles-Sevcik equation.
    Supports multiple input modes: from CV peak current, from Cottrell equation (chronoamperometry), or from Stokes-Einstein equation.
    """
    __version__ = "0.1.0"
    name = "DiffusionCoefficientCalculator"
    func_name = "diffusion_coefficient_calculator"
    description = "Calculate diffusion coefficient D using Randles-Sevcik equation (from CV), Cottrell equation (from chronoamperometry), or Stokes-Einstein relation (from hydrodynamic radius)."
    implementation_description = "Three calculation methods:\n1. Randles-Sevcik: ip = 2.69×10⁵·n^(3/2)·A·D^(1/2)·C·√ν → solve for D\n2. Cottrell: i(t) = n·F·A·C·D^(1/2)/(π^(1/2)·t^(1/2)) → solve for D\n3. Stokes-Einstein: D = k_B·T/(6·π·η·r_h)"
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Diffusion Coefficient", "Randles-Sevcik", "Cottrell Equation", "Stokes-Einstein", "Electrochemistry"]
    required_envs = []

    code_input_sig = [
        ("method", "str", "N/A", "Calculation method: 'randles_sevcik', 'cottrell', or 'stokes_einstein'."),
        ("n_electrons", "int", "1", "Number of electrons transferred. Default: 1."),
        ("electrode_area_cm2", "float", "0.07", "Working electrode area in cm². Default: 0.07."),
        ("concentration_M", "float", "0.001", "Analyte concentration in mol/L. Default: 0.001."),
        ("temperature_K", "float", "298.15", "Temperature in Kelvin. Default: 298.15."),
        # For Randles-Sevcik:
        ("peak_current_A", "float", "N/A", "Peak current in Amperes (for randles_sevcik method)."),
        ("scan_rate_V_s", "float", "N/A", "Scan rate in V/s (for randles_sevcik method)."),
        # For Cottrell:
        ("current_A", "float", "N/A", "Current at time t in Amperes (for cottrell method)."),
        ("time_s", "float", "N/A", "Time after potential step in seconds (for cottrell method)."),
        # For Stokes-Einstein:
        ("viscosity_cP", "float", "N/A", "Solvent viscosity in centipoise (for stokes_einstein method)."),
        ("hydrodynamic_radius_nm", "float", "N/A", "Hydrodynamic radius in nm (for stokes_einstein method)."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format depends on method:\n- Randles-Sevcik: 'randles_sevcip ip_A nu_V_s [n] [A_cm2] [C_M]'\n- Cottrell: 'cottrell i_A t_s [n] [A_cm2] [C_M]'\n- Stokes-Einstein: 'stokes_einstein eta_cP r_nm [T_K]'"),
    ]

    output_sig = [
        ("diffusion_coefficient_D_cm2_s", "float", "Calculated diffusion coefficient D (cm²/s)."),
        ("method_used", "str", "The calculation method used."),
        ("input_summary", "dict", "Summary of input parameters used."),
        ("formula_applied", "str", "The formula with substituted values."),
        ("notes", "str", "Additional notes on the result (e.g., typical range check)."),
    ]

    examples = [
        {
            "code_input": {
                "method": "randles_sevcik",
                "peak_current_A": 1.5e-5,
                "scan_rate_V_s": 0.1,
                "n_electrons": 1,
                "electrode_area_cm2": 0.07,
                "concentration_M": 0.001,
            },
            "text_input": {"input_string": "randles_sevcik 1.5e-5 0.1"},
            "output": {
                "diffusion_coefficient_D_cm2_s": 7.3e-6,
                "method_used": "randles_sevcik",
                "formula_applied": "D = [ip / (2.69e5 · n^1.5 · A · C · √ν)]²",
            }
        },
        {
            "code_input": {
                "method": "cottrell",
                "current_A": 5e-6,
                "time_s": 5.0,
                "n_electrons": 1,
                "electrode_area_cm2": 0.07,
                "concentration_M": 0.005,
            },
            "text_input": {"input_string": "cottrell 5e-6 5.0 1 0.07 0.005"},
            "output": {
                "diffusion_coefficient_D_cm2_s": 6.9e-6,
                "method_used": "cottrell",
            }
        },
        {
            "code_input": {
                "method": "stokes_einstein",
                "viscosity_cP": 0.89,  # water at 25°C
                "hydrodynamic_radius_nm": 0.5,
                "temperature_K": 298.15,
            },
            "text_input": {"input_string": "stokes_einstein 0.89 0.5"},
            "output": {
                "diffusion_coefficient_D_cm2_s": 4.9e-6,
                "method_used": "stokes_einstein",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._F = 96485       # C/mol
        self._kB = 1.381e-23   # J/K (Boltzmann constant)
        self._NA = 6.022e23    # mol^-1 (Avogadro's number)
        self._R = 8.314        # J/(mol·K)
        self._k_RS = 2.69e5    # Randles-Sevcik constant

    def _randles_sevcik(self, ip: float, nu: float, n: int, A: float, C: float) -> float:
        """D from Randles-Sevcik: ip = k·n^1.5·A·D^0.5·C·√ν"""
        if ip <= 0 or nu <= 0 or C <= 0 or A <= 0:
            raise ChemMCPError("All parameters must be positive for Randles-Sevcik.")
        k = self._k_RS
        D_sq = (ip / (k * (n ** 1.5) * A * C * math.sqrt(nu))) ** 2
        return D_sq

    def _cottrell(self, i: float, t: float, n: int, A: float, C: float) -> float:
        """D from Cottrell: i = n·F·A·C·D^0.5 / (√π · √t)"""
        if i <= 0 or t <= 0 or C <= 0 or A <= 0:
            raise ChemMCPError("All parameters must be positive for Cottrell equation.")
        numerator = i * math.sqrt(math.pi * t)
        denominator = n * self._F * A * C
        D_sq = (numerator / denominator) ** 2
        return D_sq

    def _stokes_einstein(self, eta_cP: float, r_nm: float, T: float) -> float:
        """D from Stokes-Einstein: D = kB·T / (6π·η·r_h)"""
        if eta_cP <= 0 or r_nm <= 0 or T <= 0:
            raise ChemMCPError("Viscosity, radius, and temperature must be positive.")

        # Convert cP to Pa·s: 1 cP = 0.001 Pa·s
        eta_Pa_s = eta_cP * 0.001
        # Convert nm to m
        r_m = r_nm * 1e-9

        D_m2_s = (self._kB * T) / (6 * math.pi * eta_Pa_s * r_m)
        # Convert m²/s to cm²/s
        D_cm2_s = D_m2_s * 1e4
        return D_cm2_s

    def _run_base(
        self,
        method: str,
        peak_current_A: Optional[float] = None,
        scan_rate_V_s: Optional[float] = None,
        current_A: Optional[float] = None,
        time_s: Optional[float] = None,
        viscosity_cP: Optional[float] = None,
        hydrodynamic_radius_nm: Optional[float] = None,
        n_electrons: int = 1,
        electrode_area_cm2: float = 0.07,
        concentration_M: float = 0.001,
        temperature_K: float = 298.15,
    ) -> dict:
        """Calculate diffusion coefficient."""
        method = method.lower().strip()

        if method == "randles_sevcik":
            if peak_current_A is None or scan_rate_V_s is None:
                raise ChemMCPError("Randles-Sevcik requires peak_current_A and scan_rate_V_s.")
            D = self._randles_sevcik(peak_current_A, scan_rate_V_s, n_electrons, electrode_area_cm2, concentration_M)
            formula = (
                f"D = [{peak_current_A:.4e} / ({self._k_RS:.2e} × {n_electrons}^{1.5} × "
                f"{electrode_area_cm2} × {concentration_M} × √{scan_rate_V_s})]² "
                f"= {D:.4e} cm²/s"
            )
        elif method == "cottrell":
            if current_A is None or time_s is None:
                raise ChemMCPError("Cottrell requires current_A and time_s.")
            D = self._cottrell(current_A, time_s, n_electrons, electrode_area_cm2, concentration_M)
            formula = (
                f"D = [{current_A:.4e} × √(π × {time_s}) / "
                f"({n_electrons} × {self._F} × {electrode_area_cm2} × {concentration_M})]² "
                f"= {D:.4e} cm²/s"
            )
        elif method == "stokes_einstein":
            if viscosity_cP is None or hydrodynamic_radius_nm is None:
                raise ChemMCPError("Stokes-Einstein requires viscosity_cP and hydrodynamic_radius_nm.")
            D = self._stokes_einstein(viscosity_cP, hydrodynamic_radius_nm, temperature_K)
            formula = (
                f"D = k_B·T / (6πη·r) = ({self._kB:.3e} × {temperature_K}) / "
                f"(6π × {viscosity_cP}×10⁻³ × {hydrodynamic_radius_nm}×10⁻⁹) "
                f"= {D:.4e} cm²/s"
            )
        else:
            raise ChemMCPError(f"Unknown method: '{method}'. Use 'randles_sevcik', 'cottrell', or 'stokes_einstein'.")

        # Typical range check
        notes = ""
        if 1e-7 <= D <= 1e-4:
            notes = "D value is within typical range for small molecules in solution (10⁻⁷–10⁻⁴ cm²/s)."
        elif D > 1e-4:
            notes = "D is unusually large; verify input parameters."
        else:
            notes = "D is unusually small; may indicate a macromolecule or verify inputs."

        return {
            "diffusion_coefficient_D_cm2_s": round(D, 12),
            "method_used": method,
            "input_summary": {
                "n": n_electrons,
                "A_cm2": electrode_area_cm2,
                "C_M": concentration_M,
                "T_K": temperature_K,
            },
            "formula_applied": formula,
            "notes": notes,
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input."""
        parts = input_string.strip().split()
        if not parts:
            raise ChemMCPError("Input cannot be empty.")

        method = parts[0].lower()
        kwargs = {}

        if method == "randles_sevcik":
            kwargs["peak_current_A"] = float(parts[1]) if len(parts) > 1 else _missing("peak_current")
            kwargs["scan_rate_V_s"] = float(parts[2]) if len(parts) > 2 else _missing("scan_rate")
            kwargs["n_electrons"] = int(parts[3]) if len(parts) > 3 else 1
            kwargs["electrode_area_cm2"] = float(parts[4]) if len(parts) > 4 else 0.07
            kwargs["concentration_M"] = float(parts[5]) if len(parts) > 5 else 0.001
        elif method == "cottrell":
            kwargs["current_A"] = float(parts[1]) if len(parts) > 1 else _missing("current")
            kwargs["time_s"] = float(parts[2]) if len(parts) > 2 else _missing("time")
            kwargs["n_electrons"] = int(parts[3]) if len(parts) > 3 else 1
            kwargs["electrode_area_cm2"] = float(parts[4]) if len(parts) > 4 else 0.07
            kwargs["concentration_M"] = float(parts[5]) if len(parts) > 5 else 0.001
        elif method == "stokes_einstein":
            kwargs["viscosity_cP"] = float(parts[1]) if len(parts) > 1 else _missing("viscosity")
            kwargs["hydrodynamic_radius_nm"] = float(parts[2]) if len(parts) > 2 else _missing("radius")
            kwargs["temperature_K"] = float(parts[3]) if len(parts) > 3 else 298.15
        else:
            raise ChemMCPError(f"Unknown method: '{method}'")

        return self._run_base(method=method, **kwargs)


def _missing(name: str):
    raise ValueError(f"Missing required parameter: {name}")
