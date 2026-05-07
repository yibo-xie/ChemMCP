import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R = 8.314462618   # J/(mol·K)
_NA = 6.02214076e23 # mol^-1
_KB = 1.380649e-23  # J/K

@ChemMCPManager.register_tool
class HeatCapacityCalculator(BaseTool):
    """
    热容计算工具 — 计算Cv、Cp，支持Debye模型和Einstein模型。
    
    用于固体热容、气体热容及Dulong-Petit极限分析。
    """
    __version__ = "0.1.0"
    name = "HeatCapacityCalculator"
    func_name = "calculate_heat_capacity"
    description = "Calculate heat capacity at constant volume (Cv) and constant pressure (Cp) using Debye model, Einstein model, ideal gas formulas, and Dulong-Petit limit."
    implementation_description = "Debye: Cv = 9·R·(T/θ_D)³ · ∫₀^(θ_D/T) x⁴·eˣ/(eˣ-1)² dx (numerical integration). Einstein: Cv = 3·R·(θ_E/T)² · e^(θ_E/T)/(e^(θ_E/T)-1)². Ideal gas: Cv_m = f/2·R, Cp_m = Cv_m + R."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Heat Capacity", "Debye Model", "Einstein Model", "Thermodynamics", "Solid State Physics"]
    required_envs = []

    code_input_sig = [
        ("model", "str", "N/A", "'debye', 'einstein', 'ideal_gas', or 'dulong_petit'"),
        ("temperature_k", "float", "N/A", "Temperature in Kelvin."),
        ("debye_temperature_k", "float", "N/A", "Debye temperature θ_D in Kelvin (for debye model)."),
        ("einstein_temperature_k", "float", "N/A", "Einstein temperature θ_E in K (for einstein model)."),
        ("degrees_of_freedom", "float", "N/A", "Number of degrees of freedom for ideal_gas (e.g., 3 for monatomic, 5 for diatomic)."),
        ("n_moles", "float", "1.0", "Amount of substance in moles."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: model|T|theta_D_or_theta_E|dof|n"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with Cv, Cp, Cv/R ratio, comparison to Dulong-Petit limit, and model-specific parameters."),
    ]

    examples = [
        {
            "code_input": {
                "model": "debye",
                "temperature_k": 200.0,
                "debye_temperature_k": 400.0,
            },
            "text_input": {
                "input_str": "debye|200|400||1"
            },
            "output": {
                "result": {
                    "model": "debye",
                    "T_K": 200.0,
                    "theta_D_K": 400.0,
                    "Cv_J_mol_K": "<value>",
                    "Cp_J_mol_K": "<value>",
                    "Cv_over_R": "<value>",
                    "Dulong_Petit_limit_J_mol_K": 24.94,
                    "approach_ratio": "<value>",
                }
            },
        },
        {
            "code_input": {
                "model": "ideal_gas",
                "temperature_k": 298.15,
                "degrees_of_freedom": 5,
            },
            "text_input": {
                "input_str": "ideal_gas|298.15||5|1"
            },
            "output": {
                "result": {
                    "model": "ideal_gas",
                    "Cv_J_mol_K": 20.79,
                    "Cp_J_mol_K": 29.10,
                    "gamma": 1.4,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _debye_integrand(self, x: float) -> float:
        """Integrand for Debye function: x^4 * e^x / (e^x - 1)^2."""
        if x < 1e-10:
            return x ** 2  # limit as x→0
        ex = math.exp(x)
        return (x ** 4 * ex) / ((ex - 1.0) ** 2)

    def _numerical_integration(self, upper_limit: int, n_steps: int = 10000) -> float:
        """Simpson's rule integration from 0 to upper_limit."""
        if upper_limit <= 0:
            return 0.0
        h = upper_limit / n_steps
        result = self._debye_integrand(0) + self._debye_integrand(upper_limit)
        for i in range(1, n_steps):
            x = i * h
            coeff = 4 if i % 2 == 1 else 2
            result += coeff * self._debye_integrand(x)
        return result * h / 3.0

    def _calc_debye(self, T: float, theta_D: float, n: float) -> dict:
        """Debye model heat capacity."""
        if theta_D <= 0:
            raise ChemMCPError("Debye temperature must be positive.")
        
        x_D = theta_D / T
        
        if x_D < 0.01:  # T >> θ_D → Dulong-Petit limit
            Cv = 3 * _R
        elif x_D > 50:  # T << θ_D → Cv → 0
            Cv = 0.0
        else:
            n_steps = max(10000, int(x_D * 25000))
            integral = self._numerical_integration(x_D, n_steps)
            Cv = 9 * _R * (1.0 / x_D ** 3) * integral
        
        Cv_mol = Cv * n
        Cp_mol = Cv_mol + n * _R  # approximate: Cp ≈ Cv + R for solids at moderate P
        dp_limit = 3 * _R
        
        return {
            "model": "debye",
            "T_K": T,
            "theta_D_K": theta_D,
            "T_over_theta_D": round(T / theta_D, 4),
            "Cv_J_mol_K": round(Cv_mol, 4),
            "Cp_J_mol_K": round(Cp_mol, 4),
            "Cv_over_R": round(Cv / _R, 4),
            "Dulong_Petit_limit_J_mol_K": round(dp_limit, 4),
            "approach_to_DP": round(Cv / dp_limit, 4),
            "n_moles": n,
        }

    def _calc_einstein(self, T: float, theta_E: float, n: float) -> dict:
        """Einstein solid heat capacity."""
        if theta_E <= 0:
            raise ChemMCPError("Einstein temperature must be positive.")
        
        x = theta_E / T
        if x > 500:
            Cv = 0.0
        elif x < 0.001:
            Cv = 3 * _R
        else:
            ex = math.exp(x)
            Cv = 3 * _R * (x ** 2) * ex / ((ex - 1.0) ** 2)
        
        Cv_mol = Cv * n
        Cp_mol = Cv_mol + n * _R
        
        return {
            "model": "einstein",
            "T_K": T,
            "theta_E_K": theta_E,
            "T_over_theta_E": round(T / theta_E, 4),
            "Cv_J_mol_K": round(Cv_mol, 4),
            "Cp_J_mol_K": round(Cp_mol, 4),
            "Cv_over_R": round(Cv / _R, 4),
            "Dulong_Petit_limit_J_mol_K": round(3 * _R, 4),
            "n_moles": n,
        }

    def _calc_ideal_gas(self, T: float, dof: float, n: float) -> dict:
        """Ideal gas heat capacity: Cv = (f/2)*R, Cp = Cv + R."""
        if dof <= 0:
            raise ChemMCPError("Degrees of freedom must be positive.")
        
        Cv = (dof / 2.0) * _R
        Cp = Cv + _R
        gamma = Cp / Cv if Cv > 0 else float('inf')
        
        return {
            "model": "ideal_gas",
            "T_K": T,
            "degrees_of_freedom": dof,
            "Cv_J_mol_K": round(Cv * n, 4),
            "Cp_J_mol_K": round(Cp * n, 4),
            "Cv_per_mole_J_mol_K": round(Cv, 4),
            "Cp_per_mole_J_mol_K": round(Cp, 4),
            "gamma": round(gamma, 4),
            "n_moles": n,
        }

    def _calc_dulong_petit(self, n: float) -> dict:
        """Dulong-Petit law: Cv ≈ 3R per mole of atoms."""
        Cv = 3 * _R * n
        Cp = Cv + n * _R
        
        return {
            "model": "dulong_petit",
            "Cv_J_mol_K": round(Cv, 4),
            "Cp_J_mol_K": round(Cp, 4),
            "validity": "High-temperature limit (T >> θ_D)",
            "n_moles": n,
        }

    def _run_base(self, model: str, temperature_k: float, debye_temperature_k: float = None,
                  einstein_temperature_k: float = None, degrees_of_freedom: float = None,
                  n_moles: float = 1.0) -> dict:
        m = model.lower().strip()
        
        if m == "debye":
            if debye_temperature_k is None:
                raise ChemMCPError("'debye' model requires debye_temperature_k.")
            return self._calc_debye(temperature_k, debye_temperature_k, n_moles)
        elif m == "einstein":
            if einstein_temperature_k is None:
                raise ChemMCPError("'einstein' model requires einstein_temperature_k.")
            return self._calc_einstein(temperature_k, einstein_temperature_k, n_moles)
        elif m == "ideal_gas":
            if degrees_of_freedom is None:
                raise ChemMCPError("'ideal_gas' model requires degrees_of_freedom.")
            return self._calc_ideal_gas(temperature_k, degrees_of_freedom, n_moles)
        elif m == "dulong_petit" or m == "dulong-petit":
            return self._calc_dulong_petit(n_moles)
        else:
            raise ChemMCPError(
                f"Unknown model: '{model}'. Options: 'debye', 'einstein', 'ideal_gas', 'dulong_petit'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            m = parts[0].strip()
            T = float(parts[1])
            tD = float(parts[2]) if len(parts) > 2 and parts[2].strip() else None
            tE = float(parts[3]) if len(parts) > 3 and parts[3].strip() else None
            dof = float(parts[4]) if len(parts) > 4 and parts[4].strip() else None
            n = float(parts[5]) if len(parts) > 5 else 1.0
            return self._run_base(m, T, tD, tE, dof, n)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'model|T|theta_D|theta_E|dof|n'")
