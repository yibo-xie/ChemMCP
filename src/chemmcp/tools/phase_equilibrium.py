import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R = 8.314462618     # J/(mol·K)
_R_Latm = 0.082057366 # L·atm/(mol·K)

@ChemMCPManager.register_tool
class PhaseEquilibrium(BaseTool):
    """
    相平衡计算工具 — 克劳修斯-克拉珀龙方程、相变温度/压力计算。
    
    支持气-液、气-固、固-液平衡计算。
    """
    __version__ = "0.1.0"
    name = "PhaseEquilibrium"
    func_name = "calculate_phase_equilibrium"
    description = "Calculate phase equilibrium using Clausius-Clapeyron equation: boiling/melting points at different pressures, vapor pressure curves, and phase boundary data."
    implementation_description = "Clausius-Clapeyron (integrated): ln(P2/P1) = -ΔH_vap/R · (1/T2 - 1/T1). Antoine equation: log10(P) = A - B/(T+C). Also supports solid-liquid (fusion) equilibria."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Phase Equilibrium", "Clausius-Clapeyron", "Boiling Point", "Vapor Pressure", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("calculation_type", "str", "N/A", "'clausius_clapeyron', 'boiling_point_at_pressure', 'pressure_at_temperature', 'triple_point_estimate', or 'antoine'"),
        ("delta_H_phase", "float", "N/A", "Enthalpy of phase change in J/mol (vaporization or fusion)."),
        ("T1_k", "float", "N/A", "Known temperature in Kelvin."),
        ("P1_atm", "float", "N/A", "Known pressure in atm at T1."),
        ("T2_k", "float", "N/A", "Target temperature K (for pressure calculation) or result."),
        ("P2_atm", "float", "N/A", "Target pressure atm (for temperature calculation) or result."),
        ("phase_type", "str", "vaporization", "'vaporization' or 'fusion'"),
        ("antoine_A", "float", "N/A", "Antoine A coefficient."),
        ("antoine_B", "float", "N/A", "Antoine B coefficient (in K)."),
        ("antoine_C", "float", "N/A", "Antoine C coefficient (in K)."),
        ("antoine_T_K", "float", "N/A", "Temperature for Antoine vapor pressure calc."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: calc_type|dH_phase|T1|P1|[T2]|[P2]|phase_type|[A,B,C,T_antoine]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with calculated temperature/pressure, phase boundary data, and intermediate values."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_type": "clausius_clapeyron",
                "delta_H_phase": 40660.0,
                "T1_k": 373.15,
                "P1_atm": 1.0,
                "T2_k": 368.15,
            },
            "text_input": {
                "input_str": "clausius_clapeyron|40660|373.15|1|368.15||vaporization"
            },
            "output": {
                "result": {
                    "P2_atm": 0.823,
                    "delta_H_J_mol": 40660,
                    "T1_K": 373.15,
                    "T2_K": 368.15,
                    "phase_type": "vaporization",
                }
            },
        },
        {
            "code_input": {
                "calculation_type": "boiling_point_at_pressure",
                "delta_H_phase": 38560.0,
                "T1_k": 353.25,
                "P1_atm": 1.0,
                "P2_atm": 0.5,
            },
            "text_input": {
                "input_str": "boiling_point_at_pressure|38560|353.25|1||0.5|vaporization"
            },
            "output": {
                "result": {
                    "T_boiling_K": "<value>",
                    "T_boiling_C": "<value>",
                }
            },
        },
        {
            "code_input": {
                "calculation_type": "antoine",
                "antoine_A": 8.07131,
                "antoine_B": 1730.63,
                "antoine_C": 233.426,
                "antoine_T_K": 300.0,
            },
            "text_input": {
                "input_str": "antoine|||||||vaporization|8.07131|1730.63|233.426|300"
            },
            "output": {
                "result": {
                    "vapor_pressure_atm": "<value>",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _clausius_clapeyron(self, dH: float, T1: float, P1: float, T2: float = None, P2: float = None) -> dict:
        """ln(P2/P1) = -ΔH/R * (1/T2 - 1/T1)."""
        if T2 is not None:
            # Solve for P2
            if T2 <= 0:
                raise ChemMCPError("Temperature must be positive.")
            ln_ratio = -dH / _R * (1.0/T2 - 1.0/T1)
            P2_calc = P1 * math.exp(ln_ratio)
            
            return {
                "method": "Clausius-Clapeyron",
                "given": {"T1_K": T1, "P1_atm": P1},
                "solved_for": "P2",
                "T2_K": T2,
                "P2_atm": round(P2_calc, 6),
                "ln_P2_over_P1": round(ln_ratio, 6),
                "delta_H_J_mol": dH,
            }
        
        elif P2 is not None:
            # Solve for T2
            if P2 <= 0:
                raise ChemMCPError("Pressure must be positive.")
            if P1 <= 0:
                raise ChemMCPError("Reference pressure must be positive.")
            ln_ratio = math.log(P2 / P1)
            inv_T2 = 1.0/T1 + _R/dH * ln_ratio
            
            if inv_T2 <= 0:
                raise ChemMCPError("Calculated temperature is non-physical (check input values).")
            
            T2_calc = 1.0 / inv_T2
            
            return {
                "method": "Clausius-Clapeyron",
                "given": {"T1_K": T1, "P1_atm": P1},
                "solved_for": "T2",
                "P2_atm": P2,
                "T2_K": round(T2_calc, 4),
                "T2_C": round(T2_calc - 273.15, 4),
                "ln_P2_over_P1": round(ln_ratio, 6),
                "delta_H_J_mol": dH,
            }
        
        else:
            raise ChemMCPError("Must provide either T2 or P2.")

    def _boiling_point_at_pressure(self, dH: float, T_b_normal: float, P_normal: float, P_target: float) -> dict:
        """Find boiling point at a different pressure using C-C equation."""
        return self._clausius_clapeyron(dH, T_b_normal, P_normal, P2=P_target)

    def _pressure_at_temperature(self, dH: float, T_known: float, P_known: float, T_target: float) -> dict:
        """Find vapor pressure at a given temperature."""
        return self._clausius_clapeyron(dH, T_known, P_known, T2=T_target)

    def _triple_point_estimate(self, dH_vap: float, dH_fus: float, T_b: float, T_m: float) -> dict:
        """
        Rough triple point estimate using intersection of sublimation and vaporization curves.
        Simplified: uses extrapolation.
        """
        # This is approximate — real triple point requires more data
        # Just return an estimate based on typical behavior
        T_tp_est = T_m * (1.0 - 0.05 * (dH_vap - dH_fus) / dH_vap)  # very rough
        
        return {
            "method": "Triple point estimate (approximate)",
            "note": "This is a rough estimate; accurate triple point determination requires experimental data.",
            "T_normal_boiling_K": T_b,
            "T_normal_melting_K": T_m,
            "delta_H_vap_J_mol": dH_vap,
            "delta_H_fus_J_mol": dH_fus,
            "T_triple_estimate_K": round(max(T_tp_est, T_m * 0.85), 2),
        }

    def _antoine(self, A: float, B: float, C: float, T_K: float) -> dict:
        """Antoine equation: log10(P/mmHg) = A - B/(T+C), convert to atm."""
        if T_K + C == 0:
            raise ChemMCPError("T + C cannot be zero.")
        log10_P_mmHg = A - B / (T_K + C)
        P_mmHg = 10 ** log10_P_mmHg
        P_atm = P_mmHg / 760.0
        
        return {
            "method": "Antoine equation",
            "A": A, "B": B, "C": C,
            "temperature_K": T_K,
            "log10_P_mmHg": round(log10_P_mmHg, 6),
            "vapor_pressure_mmHg": round(P_mmHg, 4),
            "vapor_pressure_atm": round(P_atm, 6),
        }

    def _run_base(self, calculation_type: str, delta_H_phase: float = None,
                  T1_k: float = None, P1_atm: float = None, T2_k: float = None,
                  P2_atm: float = None, phase_type: str = "vaporization",
                  antoine_A: float = None, antoine_B: float = None,
                  antoine_C: float = None, antoine_T_K: float = None) -> dict:
        ct = calculation_type.lower().strip()
        
        if ct in ("clausius_clapeyron", "cc"):
            if delta_H_phase is None or T1_k is None or P1_atm is None:
                raise ChemMCPError("Clausius-Clapeyron needs delta_H_phase, T1_k, P1_atm, plus T2_k or P2_atm.")
            return self._clausius_clapeyron(delta_H_phase, T1_k, P1_atm, T2_k, P2_atm)
        
        elif ct in ("boiling_point_at_pressure", "bp_at_p"):
            if any(x is None for x in [delta_H_phase, T1_k, P1_atm, P2_atm]):
                raise ChemMCPError("Need delta_H_phase, T1_k (normal boiling T), P1_atm (usually 1), P2_atm (target pressure).")
            return self._boiling_point_at_pressure(delta_H_phase, T1_k, P1_atm, P2_atm)
        
        elif ct in ("pressure_at_temperature", "p_at_t"):
            if any(x is None for x in [delta_H_phase, T1_k, P1_atm, T2_k]):
                raise ChemMCPError("Need delta_H_phase, T1_k, P1_atm, T2_k.")
            return self._pressure_at_temperature(delta_H_phase, T1_k, P1_atm, T2_k)
        
        elif ct in ("triple_point", "tp_estimate"):
            if any(x is None for x in [delta_H_phase, T1_k, T2_k]):
                raise ChemMCPError("Need delta_H_phase (as dH_vap), T1_k (as T_boiling), T2_k (as T_melting).")
            return self._triple_point_estimate(delta_H_phase, P1_atm or 6000, T1_k, T2_k)
        
        elif ct in ("antoine",):
            if any(x is None for x in [antoine_A, antoine_B, antoine_C, antoine_T_K]):
                raise ChemMCPError("Antoine needs coefficients A, B, C and temperature T.")
            return self._antoine(antoine_A, antoine_B, antoine_C, antoine_T_K)
        
        else:
            raise ChemMCPError(
                f"Unknown type: '{calculation_type}'. "
                f"Options: 'clausius_clapeyron', 'boiling_point_at_pressure', "
                f"'pressure_at_temperature', 'triple_point_estimate', 'antoine'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            ct = parts[0].strip()
            dH = float(parts[1]) if len(parts) > 1 and parts[1].strip() else None
            T1 = float(parts[2]) if len(parts) > 2 and parts[2].strip() else None
            P1 = float(parts[3]) if len(parts) > 3 and parts[3].strip() else None
            T2 = float(parts[4]) if len(parts) > 4 and parts[4].strip() else None
            P2 = float(parts[5]) if len(parts) > 5 and parts[5].strip() else None
            pt = parts[6] if len(parts) > 6 else "vaporization"
            AA = float(parts[7]) if len(parts) > 7 and parts[7].strip() else None
            AB = float(parts[8]) if len(parts) > 8 and parts[8].strip() else None
            AC = float(parts[9]) if len(parts) > 9 and parts[9].strip() else None
            AT = float(parts[10]) if len(parts) > 10 and parts[10].strip() else None
            return self._run_base(ct, dH, T1, P1, T2, P2, pt, AA, AB, AC, AT)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
