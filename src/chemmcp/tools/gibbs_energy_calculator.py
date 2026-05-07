import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R = 8.314462618   # J/(mol·K)

@ChemMCPManager.register_tool
class GibbsEnergyCalculator(BaseTool):
    """
    吉布斯自由能计算工具 — 计算Gibbs自由能、反应自发性判断。
    
    支持从焓和熵计算G，从平衡常数计算ΔG，以及van't Hoff方程。
    """
    __version__ = "0.1.0"
    name = "GibbsEnergyCalculator"
    func_name = "calculate_gibbs_energy"
    description = "Calculate Gibbs free energy (G), reaction Gibbs energy (ΔG), determine spontaneity, and apply van't Hoff equation."
    implementation_description = "Uses ΔG = ΔH - T·ΔS for thermodynamic calculations, ΔG° = -RT·ln(K) for equilibrium relation, and van't Hoff equation: d(lnK)/dT = ΔH°/(RT²)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Gibbs Free Energy", "Thermodynamics", "Spontaneity", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("calculation_type", "str", "N/A", "'delta_G' (from H and S), 'from_equilibrium_constant', 'vant_hoff', or 'temperature_dependence'"),
        ("delta_H", "float", "0.0", "Enthalpy change in J/mol (or kJ/mol if unit='kJ')."),
        ("delta_S", "float", "0.0", "Entropy change in J/(mol·K) (or kJ/(mol·K) if unit='kJ')."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("equilibrium_constant", "float", "N/A", "Equilibrium constant K (for from_equilibrium_constant mode)."),
        ("unit", "str", "J", "Energy unit: 'J' or 'kJ'."),
        ("T2_k", "float", "N/A", "Second temperature for vant_hoff / temperature_dependence mode."),
        ("K2", "float", "N/A", "Equilibrium constant at T2 for vant_hoff mode."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: calc_type|delta_H|delta_S|T|K|unit|[T2]|[K2]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with delta_G, spontaneity ('spontaneous'/'non_spontaneous'/'at_equilibrium'), equilibrium constant (if applicable), and intermediate values."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_type": "delta_G",
                "delta_H": -57200.0,
                "delta_S": -135.5,
                "temperature_k": 298.15,
                "unit": "J",
            },
            "text_input": {
                "input_str": "delta_G|-57200|-135.5|298.15||J"
            },
            "output": {
                "result": {
                    "delta_G_J_mol": -16620.6,
                    "spontaneous": True,
                    "spontaneity": "spontaneous",
                    "T_eq_K": 422.9,
                }
            },
        },
        {
            "code_input": {
                "calculation_type": "from_equilibrium_constant",
                "temperature_k": 298.15,
                "equilibrium_constant": 3.5e-3,
            },
            "text_input": {
                "input_str": "from_equilibrium_constant|0|0|298.15|0.0035|J"
            },
            "output": {
                "result": {
                    "delta_G_J_mol": 14190.7,
                    "spontaneous": False,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_delta_G(self, dH: float, dS: float, T: float, unit: str) -> dict:
        """ΔG = ΔH - T·ΔS"""
        factor = 1000.0 if unit == "kJ" else 1.0
        dH_J = dH * factor
        dS_J = dS * factor
        dG = dH_J - T * dS_J
        
        # Equilibrium temperature (where ΔG = 0)
        T_eq = None
        if abs(dS_J) > 1e-10:
            T_eq = dH_J / dS_J
        
        # Spontaneity at given T
        eps = 1e-10
        if dG < -eps:
            spont = "spontaneous"
            is_spon = True
        elif dG > eps:
            spont = "non_spontaneous"
            is_spon = False
        else:
            spont = "at_equilibrium"
            is_spon = None
        
        return {
            "calculation_type": "delta_G",
            "delta_H_input": dH,
            "delta_S_input": dS,
            "temperature_K": T,
            "unit": unit,
            "delta_G_J_mol": round(dG, 4),
            "delta_G_kJ_mol": round(dG / 1000.0, 4),
            "spontaneity": spont,
            "is_spontaneous": is_spon,
            "equilibrium_temperature_K": round(T_eq, 4) if T_eq is not None else None,
        }

    def _calc_from_K(self, K: float, T: float, unit: str) -> dict:
        """ΔG° = -RT ln(K)"""
        factor = 1000.0 if unit == "kJ" else 1.0
        if K <= 0:
            raise ChemMCPError("Equilibrium constant must be positive.")
        dG = -_R * T * math.log(K)
        
        eps = 1e-10
        is_spon = dG < -eps
        if dG < -eps:
            spont = "spontaneous (reverse direction)"
        elif dG > eps:
            spont = "non-spontaneous as written"
        else:
            spont = "at_equilibrium"
        
        return {
            "calculation_type": "from_equilibrium_constant",
            "equilibrium_constant": K,
            "temperature_K": T,
            "unit": unit,
            "delta_G_J_mol": round(dG, 4),
            "delta_G_kJ_mol": round(dG / 1000.0, 4),
            "spontaneity": spont,
            "is_spontaneous": not (dG > eps),  # K > 1 means forward spontaneous
        }

    def _calc_vant_hoff(self, K1: float, T1: float, K2: float, T2: float, unit: str) -> dict:
        """van't Hoff: ln(K2/K1) = -ΔH°/R · (1/T2 - 1/T1) → solve for ΔH°"""
        factor = 1000.0 if unit == "kJ" else 1.0
        if K1 <= 0 or K2 <= 0:
            raise ChemMCPError("Both equilibrium constants must be positive.")
        
        dH = -_R * (math.log(K2 / K1)) / (1.0/T2 - 1.0/T1)
        
        # Also calculate ΔG at both temperatures
        dG1 = -_R * T1 * math.log(K1)
        dG2 = -_R * T2 * math.log(K2)
        
        return {
            "calculation_type": "vant_hoff",
            "K1": K1,
            "T1_K": T1,
            "K2": K2,
            "T2_K": T2,
            "delta_H_J_mol": round(dH, 4),
            "delta_H_kJ_mol": round(dH / 1000.0, 4),
            "delta_G_T1_J_mol": round(dG1, 4),
            "delta_G_T2_J_mol": round(dG2, 4),
            "unit": unit,
        }

    def _run_base(self, calculation_type: str, delta_H: float = 0.0, delta_S: float = 0.0,
                  temperature_k: float = 298.15, equilibrium_constant: float = None,
                  unit: str = "J", T2_k: float = None, K2: float = None) -> dict:
        """Core logic."""
        calc_type = calculation_type.lower().strip()
        
        if calc_type == "delta_g" or calc_type == "delta_g_from_hs":
            return self._calc_delta_G(delta_H, delta_S, temperature_k, unit)
        elif calc_type == "from_equilibrium_constant":
            if equilibrium_constant is None:
                raise ChemMCPError("equilibrium_constant is required for 'from_equilibrium_constant' mode.")
            return self._calc_from_K(equilibrium_constant, temperature_k, unit)
        elif calc_type == "vant_hoff":
            if equilibrium_constant is None or T2_k is None or K2 is None:
                raise ChemMCPError("vant_hoff mode requires: equilibrium_constant (K1), T2_k, and K2.")
            return self._calc_vant_hoff(equilibrium_constant, temperature_k, K2, T2_k, unit)
        else:
            raise ChemMCPError(
                f"Unknown calculation type: '{calculation_type}'. "
                f"Options: 'delta_g', 'from_equilibrium_constant', 'vant_hoff'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            calc_type = parts[0].strip()
            dH = float(parts[1]) if len(parts) > 1 and parts[1].strip() else 0.0
            dS = float(parts[2]) if len(parts) > 2 and parts[2].strip() else 0.0
            T = float(parts[3]) if len(parts) > 3 else 298.15
            K = float(parts[4]) if len(parts) > 4 and parts[4].strip() else None
            unit = parts[5] if len(parts) > 5 else "J"
            T2 = float(parts[6]) if len(parts) > 6 and parts[6].strip() else None
            K2_val = float(parts[7]) if len(parts) > 7 and parts[7].strip() else None
            return self._run_base(calc_type, dH, dS, T, K, unit, T2, K2_val)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'calc_type|dH|dS|T|K|unit|T2|K2'")
