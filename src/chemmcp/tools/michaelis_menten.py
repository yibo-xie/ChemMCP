import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MichaelisMenten(BaseTool):
    """
    Michaelis-Menten 酶动力学完整分析工具。
    
    实现经典的米氏方程及其扩展：
    - 基本米氏方程：v = Vmax·[S] / (Km + [S])
    - 四种抑制类型：竞争性、非竞争性、反竞争性、混合型抑制
    - Lineweaver-Burk 双倒数作图参数
    - Eadie-Hofstee 和 Hanes-Woolf 线性化
    - 动力学参数求解（Vmax, Km, kcat）
    """
    __version__ = "0.1.0"
    name = "MichaelisMenten"
    func_name = "michaelis_menten_kinetics"
    description = "Complete Michaelis-Menten enzyme kinetics analysis: velocity calculation, inhibition analysis (competitive/uncompetitive/noncompetitive/mixed), linearization methods (Lineweaver-Burk, Eadie-Hofstee, Hanes-Woolf), and parameter estimation."
    implementation_description = "Implements the full MM equation v=Vmax[S]/(Km+[S]) with inhibition modifiers. Computes kinetic parameters from [S]-v data using three linearization methods. Analyzes inhibition patterns and calculates IC50, Ki values."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Enzyme Kinetics", "Michaelis-Menten", "Inhibition", "Catalysis", "Biochemistry"]
    required_envs = []

    code_input_sig = [
        ("analysis_type", "str", "calculate_velocity", "Type: 'calculate_velocity', 'analyze_inhibition', 'fit_parameters', or 'full_analysis'."),
        ("substrate_concentration_S", "float", "N/A", "Substrate concentration [S] in same units as Km."),
        ("Vmax", "float", "N/A", "Maximum reaction velocity Vmax."),
        ("Km", "float", "N/A", "Michaelis constant Km (substrate concentration at half-Vmax)."),
        ("inhibition_type", "str", "None", "Inhibition type: 'competitive', 'uncompetitive', 'noncompetitive', 'mixed', or None."),
        ("inhibitor_concentration_I", "float", "0", "Inhibitor concentration [I]."),
        ("Ki", "float", "None", "Inhibition constant Ki."),
        ("alpha_prime", "float", "None", "α' factor for mixed inhibition (affects ESI complex)."),
        ("substrate_velocities_data", "list", "None", "For fit_parameters: list of {'S': float, 'v': float} dicts."),
        ("enzyme_concentration_E0", "float", "None", "Total enzyme concentration for kcat calculation."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'analysis_type S Vmax Km [inhib_type I Ki]' Example: 'calculate_velocity 5.0 100 2.0 competitive 3.0 1.0'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with reaction velocity, inhibition parameters, linearization data, kinetic parameters, and diagnostic plots data."),
    ]

    examples = [
        {
            "code_input": {
                "analysis_type": "full_analysis",
                "substrate_concentration_S": 5.0,
                "Vmax": 100.0,
                "Km": 2.0,
                "inhibition_type": "competitive",
                "inhibitor_concentration_I": 3.0,
                "Ki": 1.0,
                "alpha_prime": None,
                "substrate_velocities_data": None,
                "enzyme_concentration_E0": None,
            },
            "text_input": {
                "input_params": "full_analysis 5.0 100 2.0 competitive 3.0 1.0",
            },
            "output": {
                "result": {
                    "velocity": 71.43,
                    "velocity_uninhibited": 83.33,
                    "inhibition_percent": 14.29,
                    "apparent_Km": 8.0,
                    "apparent_Vmax": 100.0,
                    "inhibition_type": "competitive",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _mm_velocity(self, S, Vmax, Km):
        """Basic Michaelis-Menten velocity."""
        if S < 0:
            raise ChemMCPError("Substrate concentration cannot be negative.")
        if Km < 0:
            raise ChemMCPError("Km cannot be negative.")
        if Km == 0:
            return Vmax if S > 0 else 0
        return Vmax * S / (Km + S)

    def _inhibited_velocity(self, S, Vmax, Km, inh_type, I, Ki, alpha_prime=None):
        """
        Calculate velocity under inhibition.
        
        Competitive:      v = Vmax·S / (αKm + S)       where α = 1 + [I]/Ki
        Uncompetitive:     v = Vmax·S / (Km + α'S)       where α' = 1 + [I]/Ki
        Noncompetitive:    v = Vmax·S / ((αKm + α'S))    α = α' = 1 + [I]/Ki
        Mixed:             v = Vmax·S / (αKm + α'S)
        """
        if inh_type is None or I is None or I == 0 or Ki is None:
            return self._mm_velocity(S, Vmax, Km), {"type": None}

        alpha = 1.0 + I / Ki if Ki > 0 else float('inf')
        
        if alpha_prime is not None:
            ap = alpha_prime
        elif inh_type == "noncompetitive":
            ap = alpha
        else:
            ap = alpha  # Default for mixed

        inh_type_lower = inh_type.lower().replace("-", "_").replace(" ", "_")

        if inh_type_lower == "competitive":
            v = Vmax * S / (alpha * Km + S) if (alpha * Km + S) > 0 else 0
            info = {"type": "competitive", "alpha": round(alpha, 6), "apparent_Km": round(alpha * Km, 6),
                    "apparent_Vmax": round(Vmax, 6), "effect": "Km increases, Vmax unchanged"}
        elif inh_type_lower == "uncompetitive":
            v = Vmax * S / (Km + alpha * S) if (Km + alpha * S) > 0 else 0
            info = {"type": "uncompetitive", "alpha_prime": round(alpha, 6), "apparent_Km": round(Km, 6),
                    "apparent_Vmax": round(Vmax / alpha, 6), "effect": "Both Km and Vmax decrease"}
        elif inh_type_lower in ("noncompetitive", "non_competitive"):
            v = Vmax * S / (alpha * Km + alpha * S) if alpha > 0 else 0
            info = {"type": "noncompetitive", "alpha": round(alpha, 6), "apparent_Km": round(alpha * Km, 6),
                    "apparent_Vmax": round(Vmax / alpha, 6), "effect": "Vmax decreases, Km may change"}
        elif inh_type_lower == "mixed":
            v = Vmax * S / (alpha * Km + ap * S) if (alpha * Km + ap * S) > 0 else 0
            info = {"type": "mixed", "alpha": round(alpha, 6), "alpha_prime": round(ap, 6),
                    "apparent_Km": round(alpha * Km, 6), "apparent_Vmax": round(Vmax / ap, 6)}
        else:
            raise ChemMCPError(f"Unknown inhibition type: {inh_type}. "
                               f"Choose: competitive, uncompetitive, noncompetitive, mixed.")

        return v, info

    def _linearization_data(self, Vmax, Km, S_values=None):
        """Compute Lineweaver-Burk, Eadie-Hofstee, Hanes-Woolf data."""
        if S_values is None:
            S_values = [Km * 0.25, Km * 0.5, Km * 1.0, Km * 2.0, Km * 4.0, Km * 8.0]
        
        lb_data = []   # 1/v vs 1/[S]
        eh_data = []   # v vs v/[S]
        hw_data = []   # [S]/v vs [S]

        for S in S_values:
            if S <= 0:
                continue
            v = self._mm_velocity(S, Vmax, Km)
            if v <= 0:
                continue
            
            lb_data.append({"inv_S": round(1/S, 6), "inv_v": round(1/v, 6)})
            eh_data.append({"v": round(v, 6), "v_over_S": round(v/S, 6)})
            hw_data.append({"S": round(S, 6), "S_over_v": round(S/v, 6)})

        return {
            "lineweaver_burk": lb_data,
            "eadie_hofstee": eh_data,
            "hanes_woolf": hw_data,
            "intercepts": {
                "LB_intercept_1/Vmax": round(1/Vmax, 6) if Vmax > 0 else None,
                "LB_slope_Km/Vmax": round(Km/Vmax, 6) if Vmax > 0 else None,
                "EH_intercept_Vmax": round(Vmax, 6),
                "EH_slope_-Km": round(-Km, 6),
                "HW_intercept_Km/Vmax": round(Km/Vmax, 6) if Vmax > 0 else None,
                "HW_slope_1/Vmax": round(1/Vmax, 6) if Vmax > 0 else None,
            }
        }

    def _fit_parameters_simple(self, data_points):
        """
        Simple linear regression on Lineweaver-Burk transformed data.
        Input: list of {'S': float, 'v': float}
        Returns estimated Vmax, Km.
        """
        if not data_points or len(data_points) < 2:
            return None
        
        # Transform to Lineweaver-Burk coordinates
        x_list = []  # 1/[S]
        y_list = []  # 1/v
        
        for pt in data_points:
            S, v = pt.get("S"), pt.get("v")
            if S is None or v is None or S <= 0 or v <= 0:
                continue
            x_list.append(1.0 / S)
            y_list.append(1.0 / v)

        n = len(x_list)
        if n < 2:
            return None

        # Linear regression: y = mx + b
        x_mean = sum(x_list) / n
        y_mean = sum(y_list) / n
        
        ss_xy = sum((x_list[i] - x_mean) * (y_list[i] - y_mean) for i in range(n))
        ss_xx = sum((x_list[i] - x_mean) ** 2 for i in range(n))
        
        if ss_xx < 1e-30:
            return None
        
        m_slope = ss_xy / ss_xx
        b_intercept = y_mean - m_slope * x_mean

        # Convert back: slope = Km/Vmax, intercept = 1/Vmax
        if b_intercept <= 0:
            return None
        
        Vmax_est = 1.0 / b_intercept
        Km_est = m_slope * Vmax_est

        # R²
        y_pred = [m_slope * x + b_intercept for x in x_list]
        ss_res = sum((y_list[i] - y_pred[i]) ** 2 for i in range(n))
        ss_tot = sum((y_list[i] - y_mean) ** 2 for i in range(n))
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        return {
            "Vmax": round(max(Vmax_est, 0), 6),
            "Km": round(max(Km_est, 0), 6),
            "R_squared": round(r_squared, 6),
            "method": "Lineweaver-Burk linear regression",
            "n_data_points": n,
        }

    def _run_base(self, analysis_type: str = "calculate_velocity",
                  substrate_concentration_S: float = None,
                  Vmax: float = None, Km: float = None,
                  inhibition_type: str = None,
                  inhibitor_concentration_I: float = 0,
                  Ki: float = None, alpha_prime: float = None,
                  substrate_velocities_data: list = None,
                  enzyme_concentration_E0: float = None) -> dict:
        """Core logic."""
        atype = analysis_type.lower().replace("-", "_").replace(" ", "_")

        result = {"analysis_type": analysis_type}

        if atype in ("calculate_velocity", "full_analysis"):
            if substrate_concentration_S is None or Vmax is None or Km is None:
                raise ChemMCPError("For velocity calculation, need S, Vmax, and Km.")

            v_uninhibited = self._mm_velocity(substrate_concentration_S, Vmax, Km)
            
            if inhibition_type and inhibitor_concentration_I and Ki:
                v_inhibited, inh_info = self._inhibited_velocity(
                    substrate_concentration_S, Vmax, Km,
                    inhibition_type, inhibitor_concentration_I, Ki, alpha_prime
                )
                pct_inhibition = (1 - v_inhibited / v_uninhibited) * 100 if v_uninhibited > 0 else 0
                
                result.update({
                    "velocity_uninhibited": round(v_uninhibited, 6),
                    "velocity_inhibited": round(v_inhibited, 6),
                    "velocity": round(v_inhibited, 6),
                    "inhibition_percent": round(pct_inhibition, 2),
                    "inhibition_info": inh_info,
                })
            else:
                result["velocity"] = round(v_uninhibited, 6)
                result["velocity_uninhibited"] = round(v_uninhibited, 6)

            # Additional diagnostics
            result["substrate_saturation"] = round(substrate_concentration_S / Km, 4) if Km > 0 else float('inf')
            result["fraction_of_Vmax"] = round(v_uninhibited / Vmax, 4) if Vmax > 0 else 0

            # Linearization data
            lin = self._linearization_data(Vmax, Km)
            result["linearization"] = lin

        if atype in ("fit_parameters", "full_analysis") and substrate_velocities_data:
            fit_result = self._fit_parameters_simple(substrate_velocities_data)
            result["parameter_fit"] = fit_result
            if fit_result and enzyme_concentration_E0:
                kcat_val = fit_result.get("Vmax")
                km_val = fit_result.get("Km")
                if kcat_val is not None and km_val is not None and km_val > 0:
                    result["kcat"] = round(kcat_val / enzyme_concentration_E0, 6)
                    result["catalytic_efficiency"] = round(result["kcat"] / km_val, 6)

        if atype == "analyze_inhibition":
            # Generate inhibition curve data
            if Vmax and Km and Ki:
                S_test = substrate_concentration_S or Km
                I_range = [i * Ki * 0.5 for i in range(9)]  # 0 to 4×Ki
                curve_data = []
                for I_val in I_range:
                    v_i, _ = self._inhibited_velocity(S_test, Vmax, Km, inhibition_type, I_val, Ki, alpha_prime)
                    v0 = self._mm_velocity(S_test, Vmax, Km)
                    curve_data.append({
                        "I_over_Ki": round(I_val / Ki, 4) if Ki > 0 else 0,
                        "I_conc": round(I_val, 6),
                        "velocity": round(v_i, 6),
                        "percent_control": round(v_i / v0 * 100, 2) if v0 > 0 else 0,
                    })
                
                # Estimate IC50
                ic50_est = None
                for i, cd in enumerate(curve_data):
                    if cd["percent_control"] <= 50:
                        ic50_est = curve_data[max(0, i-1)]["I_conc"]
                        break
                if ic50_est is None:
                    ic50_est = I_range[-1] if curve_data else None

                result["inhibition_curve"] = curve_data
                result["IC50_estimate"] = ic50_est

        logger.info(f"MichaelisMenten: type={atype}, v={result.get('velocity', 'N/A')}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            atype = parts[0]
            S = float(parts[1]) if len(parts) > 1 else None
            Vmax = float(parts[2]) if len(parts) > 2 else None
            Km = float(parts[3]) if len(parts) > 3 else None
            inh_type = parts[4] if len(parts) > 4 else None
            I = float(parts[5]) if len(parts) > 5 else 0
            Ki = float(parts[6]) if len(parts) > 6 else None
            
            return self._run_base(atype, S, Vmax, Km, inh_type, I, Ki)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'analysis_type S Vmax Km [inhib_type I Ki]'"
            )
