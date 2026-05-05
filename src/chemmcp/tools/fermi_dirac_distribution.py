import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FermiDiracDistribution(BaseTool):
    """
    Fermi-Dirac 分布计算工具。
    
    费米子（电子等）的量子统计分布：
      f(E) = 1 / (exp((E - μ) / kT) + 1)
    
    当 T → 0 时，f(E) 变为阶跃函数：E < μ 时 f=1, E > μ 时 f=0。
    """
    __version__                 = "0.1.0"
    name                        = "FermiDiracDistribution"
    func_name                   = "calculate_fermi_dirac"
    description                 = "Calculate Fermi-Dirac distribution: occupation probability for fermions as a function of energy."
    implementation_description  = "Uses f(E)=1/(exp((E-μ)/kT)+1). Computes distribution at given energies or over a range around Fermi level. Handles T≈0 limit correctly."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Fermi-Dirac", "Quantum Statistics", "Fermions", "Solid State Physics"]
    required_envs               = []

    code_input_sig   = [
        ("temperature_k",            "float",  "N/A",     "Temperature in Kelvin."),
        ("fermi_energy_j",           "float",  "N/A",     "Fermi energy (chemical potential at T=0) in Joules."),
        ("energy_values_j",          "list",   "None",    "List of specific energies in J to evaluate f(E) at."),
        ("energy_range_scale_kt",    "float",  "5.0",     "Range around E_F in units of kT for auto-generated curve data."),
        ("num_curve_points",         "int",    "200",     "Number of points in auto-generated curve."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'T(K) E_F(J) [E1 E2 ...]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with fermi_energy, distribution values, curve data, and properties."),
    ]

    examples         = [
        {
            "code_input": {
                "temperature_k":        300.0,
                "fermi_energy_j":       9.0e-20,   # ~0.56 eV (typical metal)
                "energy_values_j":      None,
                "energy_range_scale_kt":"5.0",
                "num_curve_points":     200,
            },
            "text_input": {
                "input_params":         "300.0 9e-20",
            },
            "output": {
                "result": {
                    "temperature_K": 300.0,
                    "fermi_energy_J": 9e-20,
                    "fermi_energy_eV": 0.562,
                    "at_fermi_level": {"E_J": 9e-20, "f_E": 0.5},
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
        self.k_B = 1.380649e-23   # J/K
        self.eV  = 1.602176634e-19 # J/eV

    def _f_FD(self, E: float, mu: float, T: float) -> float:
        """Evaluate Fermi-Dirac function."""
        if T <= 0:
            # T → 0 limit: step function at E = μ
            return 1.0 if E < mu else (0.5 if abs(E - mu) < 1e-30 else 0.0)
        
        x = (E - mu) / (self.k_B * T)
        
        # Handle overflow/underflow
        if x > 700:
            return 0.0
        elif x < -700:
            return 1.0
        
        return 1.0 / (math.exp(x) + 1.0)

    def _run_base(
        self,
        temperature_k: float,
        fermi_energy_j: float,
        energy_values_j: Optional[List[float]] = None,
        energy_range_scale_kt: float = 5.0,
        num_curve_points: int = 200,
    ) -> dict:
        """Core logic."""
        T = float(temperature_k)
        E_F = float(fermi_energy_j)
        
        if T < 0:
            raise ChemMCPError("Temperature cannot be negative.")

        # --- Evaluate at specific energies ---
        specific_values = None
        if energy_values_j is not None:
            specific_values = []
            for E in energy_values_j:
                Ef = float(E)
                f_val = self._f_FD(Ef, E_F, T)
                specific_values.append({
                    "E_J":          Ef,
                    "E_eV":         round(Ef / self.eV, 12),
                    "f_E":          round(f_val, 12),
                    "delta_E_eV":   round((Ef - E_F) / self.eV, 12),
                    "delta_E_kT":   round((Ef - E_F) / (self.k_B * T), 6) if T > 0 else float('inf'),
                })

        # --- Auto-generate curve data around Fermi level ---
        kT = self.k_B * T
        curve_data = None
        
        if num_curve_points > 0:
            # Range: E_F ± scale * kT
            half_range = energy_range_scale_kt * kT
            E_min = E_F - half_range
            E_max = E_F + half_range
            dE = (E_max - E_min) / num_curve_points if num_curve_points > 0 else 0
            
            points = []
            for i in range(num_curve_points + 1):
                E = E_min + i * dE
                f_val = self._f_FD(E, E_F, T)
                points.append({
                    "E_J":     round(E, 20),
                    "E_eV":    round(E / self.eV, 14),
                    "f_E":     round(f_val, 12),
                    "delta_E_kT": round((E - E_F) / kT, 8) if kT > 0 else (1.0 if E >= E_F else -1.0),
                })
            
            curve_data = {
                "E_min_J":       E_min,
                "E_max_J":       E_max,
                "E_min_eV":      round(E_min / self.eV, 12),
                "E_max_eV":      round(E_max / self.eV, 12),
                "range_in_kT":   energy_range_scale_kt,
                "num_points":    len(points),
                "data":          points[:25],
                "total_points_computed": len(points),
            }

        # Value at Fermi level is always 0.5 (for any finite T, by definition when μ=E_F)
        result = {
            "temperature_K":           T,
            "kT_J":                    round(kT, 25),
            "kT_eV":                   round(kT / self.eV, 10),
            "fermi_energy_J":          E_F,
            "fermi_energy_eV":         round(E_F / self.eV, 12),
            "at_fermi_level":          {
                "E_J": E_F,
                "f_E": 0.5,
            },
            "specific_evaluations":    specific_values,
            "curve_data":              curve_data,
            "formula":                 "f(E) = 1/(exp((E-μ)/kT)+1)",
        }

        logger.info(f"Fermi-Dirac: T={T}K, E_F={E_F}J ({E_F/self.eV:.4f}eV)")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 'T E_F' params.")
            
            T = float(parts[0])
            E_F = float(parts[1])
            energies = [float(p) for p in parts[2:]] if len(parts) > 2 else None
            
            return self._run_base(T, E_F, energies)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T(K) E_F(J) [E1 E2 ...]'")
