import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

R_ATM = 0.08206   # L·atm/(K·mol)
R_J = 8.314       # J/(mol·K)


@ChemMCPManager.register_tool
class EquilibriumConstant(BaseTool):
    """
    综合平衡常数计算工具：Kp、Kc、Ka、Kb、Kw 及相互转换。
    
    支持功能：
    - Kc ↔ Kp 转换（气相反应）
    - Ka ↔ Kb 转换（共轭酸碱对，Kw = Ka × Kb）
    - 从 ΔG° 计算 K
    - 从 ΔH° 和 ΔS° 计算 K
    - 水的离子积 Kw 随温度变化
    """
    __version__ = "0.1.0"
    name = "EquilibriumConstant"
    func_name = "calculate_equilibrium_constant"
    description = "Comprehensive equilibrium constant calculator: Kp/Kc conversion, Ka/Kb conjugate pairs, K from ΔG°, Kw vs temperature."
    implementation_description = "Supports: (1) Kc↔Kp via Kp=Kc(RT)^Δn; (2) Ka↔Kb via Kw=Ka·Kb; (3) K=exp(-ΔG°/RT); (4) K from ΔH°&ΔS° via ΔG°=ΔH°-TΔS°; (5) Kw temperature dependence."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Equilibrium", "Physical Chemistry", "Acid-Base", "Kp", "Kc", "Ka", "Kb"]
    required_envs = []

    code_input_sig = [
        ("calc_type", "str", "N/A", "Calculation type: 'kc_to_kp', 'kp_to_kc', 'ka_to_kb', 'kb_to_ka', 'deltaG_to_K', 'dHS_to_K', 'kw_at_T'."),
        ("value", "float", "N/A", "Input constant value (Kc, Kp, Ka, Kb, ΔG° in kJ/mol, or ΔH° in kJ/mol)."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (default 298.15)."),
        ("delta_n", "float", "0.0", "Change in gas moles Σν_gas(products) - Σν_gas(reactants) for Kc/Kp conversion."),
        ("delta_s", "float", "None", "ΔS° in J/(mol·K) (needed for dHS_to_K)."),
        ("kw_at_t", "float", "1e-14", "Kw value at given T (for ka/kb conversions, default 1e-14 at 25°C)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'calc_type value [temperature_k] [delta_n] [delta_s] [kw]'. Example: 'kc_to_kp 0.067 298 -1'"),
    ]

    output_sig = [
        ("result", "float", "Calculated equilibrium constant value."),
        ("explanation", "str", "Step-by-step calculation details with formula used."),
    ]

    examples = [
        {
            "code_input": {
                "calc_type": "kc_to_kp",
                "value": 0.067,
                "temperature_k": 298.0,
                "delta_n": -1.0,
                "delta_s": None,
                "kw_at_t": 1e-14,
            },
            "text_input": {
                "input_params": "kc_to_kp 0.067 298 -1",
            },
            "output": {
                "result": round(0.067 * (R_ATM * 298.0) ** (-1), 4),
                "explanation": f"2NO₂⇌N₂O₄: Kp = Kc(RT)^Δn = 0.067×(0.0821×298)⁻¹ = {round(0.067*(R_ATM*298)**(-1),4)}",
            },
        },
        {
            "code_input": {
                "calc_type": "ka_to_kb",
                "value": 1.8e-5,
                "temperature_k": 298.15,
                "delta_n": 0.0,
                "delta_s": None,
                "kw_at_t": 1e-14,
            },
            "text_input": {
                "input_params": "ka_to_kb 1.8e-5",
            },
            "output": {
                "result": 5.56e-10,
                "explanation": "Kb = Kw/Ka = 1.0×10⁻¹⁴ / 1.8×10⁻⁵ = 5.56×10⁻¹⁰ (NH₃/NH₄⁺ conjugate pair)",
            },
        },
        {
            "code_input": {
                "calc_type": "deltaG_to_K",
                "value": -33.0,
                "temperature_k": 298.15,
                "delta_n": 0.0,
                "delta_s": None,
                "kw_at_t": 1e-14,
            },
            "text_input": {
                "input_params": "deltaG_to_K -33.0 298.15",
            },
            "output": {
                "result": round(math.exp(33000 / (R_J * 298.15)), 4),
                "explanation": "K = exp(-ΔG°/RT) = exp(33000/(8.314×298.15)) ≈ large value",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, calc_type: str, value: float, temperature_k: float = 298.15,
                  delta_n: float = 0.0, delta_s: float = None,
                  kw_at_t: float = 1e-14) -> dict:
        """Core logic: compute equilibrium constants."""
        calc_type = calc_type.lower().strip()

        if calc_type == "kc_to_kp":
            # Kp = Kc * (R*T)^Δn
            rt = R_ATM * temperature_k
            if delta_n == 0:
                result = value
                explanation = f"Δn=0 → Kp = Kc = {result}"
            else:
                factor = rt ** delta_n
                result = value * factor
                explanation = (
                    f"Kp = Kc · (RT)^Δn\n"
                    f"Kc = {value}, R = {R_ATM} L·atm/(K·mol), T = {temperature_k} K, Δn = {delta_n}\n"
                    f"RT = {rt:.4f}\n"
                    f"(RT)^Δn = ({rt:.4f})^{delta_n} = {factor:.6g}\n"
                    f"Kp = {value} × {factor:.6g} = {result:.6g}"
                )

        elif calc_type == "kp_to_kc":
            # Kc = Kp / (R*T)^Δn = Kp * (R*T)^(-Δn)
            rt = R_ATM * temperature_k
            if delta_n == 0:
                result = value
                explanation = f"Δn=0 → Kc = Kp = {result}"
            else:
                factor = rt ** (-delta_n)
                result = value * factor
                explanation = (
                    f"Kc = Kp · (RT)^(-Δn)\n"
                    f"Kp = {value}, RT = {rt:.4f}, Δn = {delta_n}\n"
                    f"Kc = {value} × ({rt:.4f})^{(-delta_n)} = {result:.6g}"
                )

        elif calc_type == "ka_to_kb":
            # Kb = Kw / Ka
            if value <= 0:
                raise ChemMCPError("Ka must be positive.")
            result = kw_at_t / value
            pKa = -math.log10(value)
            pKb = -math.log10(result) if result > 0 else float('inf')
            explanation = (
                f"Kb = Kw / Ka\n"
                f"Kw = {kw_at_t:.0e}, Ka = {value:.4e}\n"
                f"Kb = {kw_at_t:.0e} / {value:.4e} = {result:.4e}\n"
                f"pKa = {pKa:.2f}, pKb = {pKb:.2f}"
            )

        elif calc_type == "kb_to_ka":
            # Ka = Kw / Kb
            if value <= 0:
                raise ChemMCPError("Kb must be positive.")
            result = kw_at_t / value
            pKb = -math.log10(value)
            pKa = -math.log10(result) if result > 0 else float('inf')
            explanation = (
                f"Ka = Kw / Kb\n"
                f"Kw = {kw_at_t:.0e}, Kb = {value:.4e}\n"
                f"Ka = {kw_at_t:.0e} / {value:.4e} = {result:.4e}\n"
                f"pKa = {pKa:.2f}, pKb = {pKb:.2f}"
            )

        elif calc_type == "deltag_to_k" or calc_type == "deltaG_to_K":
            # K = exp(-ΔG°/RT), ΔG° in kJ/mol → convert to J/mol
            dg_j_mol = abs(value) * 1000  # Convert kJ → J
            if value >= 0:
                # Non-spontaneous: K < 1
                result = math.exp(-dg_j_mol / (R_J * temperature_k))
                sign_str = "-"
            else:
                # Spontaneous: K > 1
                result = math.exp(dg_j_mol / (R_J * temperature_k))
                sign_str = "+"
            explanation = (
                f"K = exp(-ΔG°/RT)\n"
                f"ΔG° = {value} kJ/mol = {value*1000} J/mol\n"
                f"R = {R_J} J/(mol·K), T = {temperature_k} K\n"
                f"-ΔG°/RT = {-value*1000}/({R_J}×{temperature_k}) = {-value*1000/(R_J*temperature_k):.4f}\n"
                f"K = exp({-value*1000/(R_J*temperature_k):.4f}) = {result:.6g}"
            )

        elif calc_type == "dhs_to_k" or calc_type == "dHS_to_K":
            # ΔG° = ΔH° - TΔS°, then K = exp(-ΔG°/RT)
            if delta_s is None:
                raise ChemMCPError("delta_s (ΔS° in J/(mol·K)) is required for dHS_to_K mode.")
            dh_j = value * 1000  # kJ → J
            dg_j = dh_j - temperature_k * delta_s  # ΔG° in J/mol
            result = math.exp(-dg_j / (R_J * temperature_k))
            explanation = (
                f"ΔG° = ΔH° - TΔS°\n"
                f"ΔH° = {value} kJ/mol = {dh_j} J/mol\n"
                f"ΔS° = {delta_s} J/(mol·K), T = {temperature_k} K\n"
                f"ΔG° = {dh_j} - {temperature_k}×{delta_s} = {dg_j:.2f} J/mol = {dg_j/1000:.2f} kJ/mol\n"
                f"K = exp(-ΔG°/RT) = exp({-dg_j/(R_J*temperature_k):.4f}) = {result:.6g}"
            )

        elif calc_type == "kw_at_t":
            # Approximate Kw temperature dependence using van't Hoff
            # Kw at 25°C = 1.0e-14, ΔH°_ion ≈ 55.8 kJ/mol
            dH_ion = 55830  # J/mol, enthalpy of water autoionization
            T_ref = 298.15
            kw_ref = 1.0e-14
            ln_ratio = (-dH_ion / R_J) * (1.0 / temperature_k - 1.0 / T_ref)
            result = kw_ref * math.exp(ln_ratio)
            explanation = (
                f"van't Hoff equation for Kw:\n"
                f"ln(Kw_T/Kw_ref) = -(ΔH°/R)(1/T - 1/T_ref)\n"
                f"ΔH°_ion = 55.83 kJ/mol, T_ref = 298.15 K, Kw_ref = 1.0×10⁻¹⁴\n"
                f"ln(Kw_{temperature_k}/10⁻¹⁴) = -(55830/{R_J})(1/{temperature_k} - 1/298.15)\n"
                f"= {ln_ratio:.4f}\n"
                f"Kw({temperature_k} K) = {result:.4e}"
            )
            return {"result": result, "explanation": explanation}

        else:
            valid = [
                "kc_to_kp", "kp_to_kc", "ka_to_kb", "kb_to_ka",
                "deltag_to_k", "deltaG_to_K", "dhs_to_k", "dHS_to_K", "kw_at_t"
            ]
            raise ChemMCPError(f"Unknown calc_type '{calc_type}'. Valid types: {valid}")

        # Preserve scientific notation for very small/large values
        if abs(result) < 1e-4 or abs(result) > 1e6:
            result_out = result  # keep full precision for extreme values
        else:
            result_out = round(result, 6)

        logger.info(f"EquilibriumConstant: {calc_type} → {result:.6g}")
        return {"result": result_out, "explanation": explanation}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            ctype = parts[0]
            val = float(parts[1])
            T = float(parts[2]) if len(parts) > 2 else 298.15
            dn = float(parts[3]) if len(parts) > 3 else 0.0
            ds = float(parts[4]) if len(parts) > 4 else None
            kw = float(parts[5]) if len(parts) > 5 else 1e-14
            return self._run_base(ctype, val, T, dn, ds, kw)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'calc_type value [T] [dn] [ds] [kw]'")
