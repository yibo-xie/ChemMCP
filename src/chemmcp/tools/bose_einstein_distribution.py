import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BoseEinsteinDistribution(BaseTool):
    """
    Bose-Einstein 分布计算工具。
    
    玻色子（光子、声子等）的量子统计分布：
      f(E) = 1 / (exp((E - μ) / kT) - 1)
    
    注意：对于玻色子，化学势 μ ≤ 0（或 μ < ε_0），否则 f(E) 可以为负（无物理意义）。
    光子气体：μ = 0。
    """
    __version__                 = "0.1.0"
    name                        = "BoseEinsteinDistribution"
    func_name                   = "calculate_bose_einstein"
    description                 = "Calculate Bose-Einstein distribution: occupation probability for bosons as a function of energy."
    implementation_description  = "Uses f(E)=1/(exp((E-μ)/kT)-1). Validates μ ≤ 0 to ensure physical results. Handles T→0 and E→μ singularities. Supports photon gas (μ=0) mode."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Bose-Einstein", "Quantum Statistics", "Bosons", "Photon Gas", "BE Condensation"]
    required_envs               = []

    code_input_sig   = [
        ("temperature_k",            "float",  "N/A",     "Temperature in Kelvin."),
        ("chemical_potential_mu_j",   "float",  "0.0",     "Chemical potential in Joules (must be ≤ 0; use 0 for photon gas)."),
        ("energy_values_j",          "list",   "None",    "List of specific energies in J to evaluate f(E) at."),
        ("energy_min_j",             "float",  "1e-23",   "Minimum energy for auto-generated curve data (J)."),
        ("energy_max_j",             "float",  "1e-20",   "Maximum energy for auto-generated curve data (J)."),
        ("num_curve_points",         "int",    "200",     "Number of points in curve."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'T(K) mu_J [E1 E2 ...]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with chemical_potential, distribution values, curve data, and BEC condition check."),
    ]

    examples         = [
        {
            "code_input": {
                "temperature_k":          5000.0,
                "chemical_potential_mu_j": 0.0,   # photon gas
                "energy_values_j":        None,
                "energy_min_j":           1e-23,
                "energy_max_j":           1e-20,
                "num_curve_points":       200,
            },
            "text_input": {
                "input_params":           "5000.0 0.0",
            },
            "output": {
                "result": {
                    "temperature_K": 5000.0,
                    "mu_J": 0.0,
                    "mode": "photon_gas",
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

    def _f_BE(self, E: float, mu: float, T: float) -> float:
        """Evaluate Bose-Einstein function."""
        if T <= 0:
            # At T=0: all bosons condense into ground state
            if abs(E - mu) < 1e-30:
                return float('inf')  # divergence at E=μ
            return 0.0 if E > mu else float('inf')
        
        x = (E - mu) / (self.k_B * T)
        
        if x < 1e-10:
            # E ≈ μ: singularity — return large number
            return 1e15
        
        exp_x = math.exp(x)
        
        if exp_x <= 1.0 + 1e-15:
            # Would give negative or infinite occupation
            return 1e15
        
        result = 1.0 / (exp_x - 1.0)
        
        if result < 0 or result > 1e15:
            return min(max(result, 0), 1e15)
        
        return result

    def _run_base(
        self,
        temperature_k: float,
        chemical_potential_mu_j: float = 0.0,
        energy_values_j: Optional[List[float]] = None,
        energy_min_j: Optional[float] = 1e-23,
        energy_max_j: Optional[float] = 1e-20,
        num_curve_points: int = 200,
    ) -> dict:
        """Core logic."""
        T = float(temperature_k)
        mu = float(chemical_potential_mu_j)
        
        if T < 0:
            raise ChemMCPError("Temperature cannot be negative.")
        if mu > 1e-15:
            raise ChemMCPError("Chemical potential μ must be ≤ 0 for Bose-Einstein statistics.")

        mode = "photon_gas" if abs(mu) < 1e-22 else "general_boson"

        # --- Evaluate at specific energies ---
        specific_values = None
        if energy_values_j is not None:
            specific_values = []
            for E in energy_values_j:
                Ef = float(E)
                f_val = self._f_BE(Ef, mu, T)
                specific_values.append({
                    "E_J":          Ef,
                    "E_eV":         round(Ef / self.eV, 12),
                    "f_E":          round(f_val, 8) if f_val < 1e14 else float('inf'),
                    "delta_E_eV":   round((Ef - mu) / self.eV, 12),
                })

        # --- Auto-generate curve data ---
        curve_data = None
        if num_curve_points > 0:
            E_min = float(energy_min_j)
            E_max = float(energy_max_j)
            
            if E_min <= mu:
                E_min = mu + 1e-25  # avoid exact singularity
            
            dE = (E_max - E_min) / num_curve_points
            
            points = []
            for i in range(num_curve_points + 1):
                E = E_min + i * dE
                f_val = self._f_BE(E, mu, T)
                points.append({
                    "E_J":     round(E, 22),
                    "E_eV":    round(E / self.eV, 14),
                    "f_E":     round(f_val, 8) if f_val < 1e14 else float('inf'),
                })
            
            curve_data = {
                "E_min_J":       E_min,
                "E_max_J":       E_max,
                "num_points":    len(points),
                "data":          points[:25],
                "total_points_computed": len(points),
            }

        # --- BEC condition check ---
        # For an ideal Bose gas in 3D box: BEC occurs when T ≤ T_c
        # where kT_c ∝ (n/ζ(3/2))^(2/3) · ℏ²/m
        # Here we just note the condition
        bec_note = (
            "For photon gas (μ=0): no BEC, but Planck distribution applies. "
            "For massive bosons with μ<0: BEC possible when n·λ_dB³ ≥ ζ(3/2) ≈ 2.612."
        )

        result = {
            "temperature_K":         T,
            "kT_J":                  round(self.k_B * T, 25),
            "kT_eV":                 round(self.k_B * T / self.eV, 10),
            "mu_J":                  mu,
            "mu_eV":                 round(mu / self.eV, 12),
            "mode":                  mode,
            "specific_evaluations":  specific_values,
            "curve_data":            curve_data,
            "bec_condition_note":    bec_note,
            "formula":               "f(E) = 1/(exp((E-μ)/kT)-1)",
        }

        logger.info(f"Bose-Einstein: T={T}K, μ={mu}J ({mu/self.eV:.4f}eV), mode={mode}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 'T mu' params.")
            
            T = float(parts[0])
            mu = float(parts[1])
            energies = [float(p) for p in parts[2:]] if len(parts) > 2 else None
            
            return self._run_base(T, mu, energies)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T(K) mu(J) [E1 E2 ...]'")
