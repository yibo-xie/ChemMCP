import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TafelEquation(BaseTool):
    """
    Tafel方程分析电极动力学工具。
    η = a + b·log₁₀(j)  或  η = b·log(j/j₀)
    其中: η = 过电势, j = 电流密度, j₀ = 交换电流密度, b = Tafel斜率
    b = (2.303·R·T) / (α·n·F)
    支持计算：过电势、Tafel斜率、交换电流密度、传递系数等。
    """
    __version__      = "0.1.0"
    name             = "TafelEquation"
    func_name        = "tafel_equation"
    description      = "Analyze electrode kinetics using the Tafel equation: overpotential vs current density relationship."
    implementation_description = "Implements the Tafel equation: eta = a + b*log10(j) or eta = b*log10(j/j0). Computes Tafel slope, exchange current density, transfer coefficient, and overpotential from kinetic parameters."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Tafel Equation", "Electrode Kinetics", "Electrochemistry", "Overpotential", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("calc_mode", "str", "N/A", "What to calculate: 'overpotential', 'tafel_slope', 'exchange_current', 'transfer_coeff', or 'analyze'."),
        ("current_density_a_m2", "float", "N/A", "Current density j in A/m² (for most modes)."),
        ("j0_a_m2", "float", "N/A", "Exchange current density j₀ in A/m²."),
        ("alpha_transfer", "float", "0.5", "Charge transfer coefficient α (0 < α ≤ 1). Default: 0.5."),
        ("n_electrons", "int", "1", "Number of electrons in rate-determining step. Default: 1."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin. Default: 298.15 K (25°C)."),
        ("overpotential_v", "float", "0", "Overpotential η in V (for 'exchange_current' mode only)."),
        ("tafel_slope_v_dec", "float", "0", "Tafel slope b in V/decade (for 'analyze' mode only)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated string: 'mode;j;j0[;alpha][;n][;T][;eta_or_b]'. Example: 'overpotential;100;0.1;0.5;1;298.15'."),
    ]

    output_sig       = [
        ("calc_mode", "str", "Calculation mode used."),
        ("temperature_K", "float", "Temperature used (K)."),
        ("n_electrons", "int", "Number of electrons."),
        ("alpha_transfer", "float", "Charge transfer coefficient α."),
        ("j0_A_m2", "float", "Exchange current density j₀ (A/m²)."),
        ("j_A_m2", "float", "Current density j (A/m²)."),
        ("tafel_slope_V_per_dec", "float", "Tafel slope b (V/decade)."),
        ("overpotential_V", "float", "Overpotential η (V)."),
        ("tafel_constant_a_V", "float", "Tafel intercept a = -b·log₁₀(j₀) (V)."),
        ("kinetic_regime", "str", "'tafel' (high η), 'linear' (low η), or 'transition'."),
        ("summary", "str", "Human-readable summary of the Tafel analysis."),
    ]

    examples         = [
        {
            "code_input": {
                "calc_mode": "overpotential",
                "current_density_a_m2": 100.0,
                "j0_a_m2": 0.1,
                "alpha_transfer": 0.5,
                "n_electrons": 1,
                "temperature_k": 298.15,
                "overpotential_v": 0.0,
                "tafel_slope_v_dec": 0.0,
            },
            "text_input": {
                "input_params": "overpotential;100;0.1;0.5;1"
            },
            "output": {
                "calc_mode": "overpotential",
                "temperature_K": 298.15,
                "n_electrons": 1,
                "alpha_transfer": 0.5,
                "j0_A_m2": 0.1,
                "j_A_m2": 100.0,
                "tafel_slope_V_per_dec": 0.0592,
                "overpotential_V": 0.354,
                "tafel_constant_a_V": 0.1183,
                "kinetic_regime": "tafel",
                "summary": "Tafel: eta=0.354 V at j=100 A/m2, b=59.2 mV/dec.",
            }
        },
        {
            "code_input": {
                "calc_mode": "tafel_slope",
                "current_density_a_m2": 0.0,
                "j0_a_m2": 0.0,
                "alpha_transfer": 0.5,
                "n_electrons": 2,
                "temperature_k": 298.15,
                "overpotential_v": 0.0,
                "tafel_slope_v_dec": 0.0,
            },
            "text_input": {
                "input_params": "tafel_slope;;0.5;2;298.15"
            },
            "output": {
                "calc_mode": "tafel_slope",
                "temperature_K": 298.15,
                "n_electrons": 2,
                "alpha_transfer": 0.5,
                "j0_A_m2": 0.0,
                "j_A_m2": 0.0,
                "tafel_slope_V_per_dec": 0.0296,
                "overpotential_V": 0.0,
                "tafel_constant_a_V": 0.0,
                "kinetic_regime": "linear",
                "summary": "Tafel slope for n=2, alpha=0.5: b = 29.6 mV/dec.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618     # J/(mol·K)
        self.F = 96485.33212     # C/mol

    def _run_base(self, calc_mode: str, current_density_a_m2: float, j0_a_m2: float,
                  alpha_transfer: float = 0.5, n_electrons: int = 1,
                  temperature_k: float = 298.15,
                  overpotential_v: float = 0.0, tafel_slope_v_dec: float = 0.0) -> dict:
        """Analyze electrode kinetics using Tafel equation."""
        R = self.R
        F = self.F
        T = temperature_k
        alpha = alpha_transfer
        n = n_electrons

        if alpha <= 0 or alpha > 1:
            raise ChemMCPError("Transfer coefficient α must be in range (0, 1].")
        if n <= 0:
            raise ChemMCPError("Number of electrons must be positive.")
        if T <= 0:
            raise ChemMCPError("Temperature must be positive.")

        mode = calc_mode.lower().replace("-", "_")
        valid_modes = {"overpotential", "tafel_slope", "exchange_current", "transfer_coeff", "analyze"}
        if mode not in valid_modes:
            raise ChemMCPError(f"Unknown calc_mode '{calc_mode}'. Use: {valid_modes}")

        # Tafel slope: b = (2.303 * R * T) / (α * n * F)  [V/decade]
        b = (2.302585 * R * T) / (alpha * n * F)

        result = {
            "calc_mode": mode,
            "temperature_K": T,
            "n_electrons": n,
            "alpha_transfer": alpha,
            "tafel_slope_V_per_dec": round(b, 6),
        }

        if mode == "tafel_slope":
            result["summary"] = (
                f"Tafel slope calculation:\n"
                f"b = 2.303RT/(αnF) = 2.303×{R:.3f}×{T:.1f}/({alpha}×{n}×{F:.1f})\n"
                f"b = {b*1000:.1f} mV/decade"
            )
            return result

        j0 = j0_a_m2
        j = current_density_a_m2

        # Tafel intercept: a = -b * log10(j0)
        a = -b * math.log10(j0) if j0 > 0 else 0

        if mode == "overpotential":
            if j <= 0:
                raise ChemMCPError("Current density must be positive.")
            if j0 <= 0:
                raise ChemMCPError("Exchange current density must be positive.")
            # η = b * log10(j/j0)
            eta = b * math.log10(j / j0)
            regime = "tafel" if abs(eta) > 0.1 else ("transition" if abs(eta) > 0.01 else "linear")

            result.update({
                "j0_A_m2": j0,
                "j_A_m2": j,
                "overpotential_V": round(eta, 6),
                "tafel_constant_a_V": round(a, 6),
                "kinetic_regime": regime,
                "summary": (
                    f"Tafel equation: η = b·log₁₀(j/j₀)\n"
                    f"j = {j} A/m², j₀ = {j0} A/m²\n"
                    f"b = {b*1000:.1f} mV/dec\n"
                    f"η = {b:.4f} × log₁₀({j}/{j0}) = {eta:.4f} V\n"
                    f"Regime: {regime}"
                ),
            })

        elif mode == "exchange_current":
            eta = overpotential_v
            if abs(b) < 1e-12:
                raise ChemMCPError("Tafel slope is too small.")
            # j = j0 * 10^(η/b)
            log_ratio = eta / b
            ratio = 10.0 ** log_ratio
            j_calc = j0 * ratio

            result.update({
                "j0_A_m2": j0,
                "j_A_m2": round(j_calc, 6),
                "overpotential_V": round(eta, 6),
                "tafel_constant_a_V": round(a, 6),
                "kinetic_regime": "tafel" if abs(eta) > 0.1 else "linear",
                "summary": (
                    f"From overpotential η = {eta:.4f} V:\n"
                    f"log₁₀(j/j₀) = η/b = {eta:.4f}/{b:.4f} = {log_ratio:.4f}\n"
                    f"j = j₀ × 10^{log_ratio:.4f} = {j_calc:.4f} A/m²"
                ),
            })

        elif mode == "analyze":
            # Given experimental b and j0, analyze kinetics
            result.update({
                "j0_A_m2": j0,
                "tafel_constant_a_V": round(a, 6),
                "summary": (
                    f"Tafel kinetic analysis:\n"
                    f"j₀ = {j0} A/m²\n"
                    f"b = {b*1000:.1f} mV/dec (experimental: {tafel_slope_v_dec*1000:.1f} mV/dec)\n"
                    f"a = -b·log₁₀(j₀) = {a:.4f} V\n"
                    f"Equation: η = {a:.4f} + {b:.4f}·log₁₀(j)"
                ),
            })

        elif mode == "transfer_coeff":
            # Solve for α from measured b and known n
            if tafel_slope_v_dec > 0:
                alpha_calc = (2.302585 * R * T) / (n * F * tafel_slope_v_dec)
                result["alpha_transfer"] = round(alpha_calc, 4)
                result["summary"] = (
                    f"Solving for transfer coefficient:\n"
                    f"α = 2.303RT/(nFb) = {alpha_calc:.4f}"
                )

        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse semicolon-separated text input."""
        parts = input_params.strip().split(";")
        if len(parts) < 1:
            raise ChemMCPError("Text input requires at least calc_mode.")

        try:
            mode = parts[0].strip()
            j = float(parts[1].strip()) if len(parts) > 1 and parts[1].strip() else 0.0
            j0 = float(parts[2].strip()) if len(parts) > 2 and parts[2].strip() else 0.0
            alpha = float(parts[3].strip()) if len(parts) > 3 and parts[3].strip() else 0.5
            n = int(parts[4].strip()) if len(parts) > 4 and parts[4].strip() else 1
            T = float(parts[5].strip()) if len(parts) > 5 and parts[5].strip() else 298.15
            extra = float(parts[6].strip()) if len(parts) > 6 and parts[6].strip() else 0.0
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse values from '{input_params}': {e}")

        kwargs = {
            "calc_mode": mode, "current_density_a_m2": j, "j0_a_m2": j0,
            "alpha_transfer": alpha, "n_electrons": n, "temperature_k": T,
        }
        if mode == "exchange_current":
            kwargs["overpotential_v"] = extra
        elif mode in ("analyze", "transfer_coeff"):
            kwargs["tafel_slope_v_dec"] = extra

        return self._run_base(**kwargs)
