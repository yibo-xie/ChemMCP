import logging
import math
from typing import List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ArrheniusAnalyzer(BaseTool):
    """
    Arrhenius方程分析工具。
    根据不同温度下的速率常数数据，求活化能 Ea 和指前因子 A。
    支持线性拟合和作图数据输出。
    """
    __version__ = "0.1.0"
    name = "ArrheniusAnalyzer"
    func_name = "analyze_arrhenius"
    description = "Analyze Arrhenius equation: determine activation energy (Ea) and pre-exponential factor (A) from k vs T data."
    implementation_description = "Uses linearized form: ln(k) = ln(A) - Ea/(R·T). Performs linear regression of ln(k) vs 1/T to extract Ea and A."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Arrhenius", "Activation Energy", "Chemical Kinetics"]
    required_envs    = []

    code_input_sig = [
        ("temperatures_k", "str", "N/A", "Temperatures in Kelvin, comma-separated."),
        ("rate_constants", "str", "N/A", "Rate constants at each temperature, comma-separated (same order and length as temperatures_k)."),
        ("r_gas", "float", "8.314", "Gas constant R in J/(mol·K) (default 8.314)."),
        ("ea_unit", "str", "kJ/mol", "Unit for activation energy output: 'J/mol' or 'kJ/mol'."),
        # Optional: if you already have linearized data
        ("inverse_temperatures", "str", "", "Pre-computed 1/T values (K⁻¹), comma-separated (optional, overrides temperatures_k)."),
        ("ln_rate_constants", "str", "", "Pre-computed ln(k) values, comma-separated (optional, overrides rate_constants)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: temperatures_k rate_constants [r_gas ea_unit]."),
    ]

    output_sig = [
        ("activation_energy", "float", "Activation energy Ea in specified unit."),
        ("pre_exponential_factor", "float", "Pre-exponential factor A (same unit as input k)."),
        ("ea_unit", "str", "Unit of Ea."),
        ("a_unit", "str", "Unit of A."),
        ("r_squared", "float", "R² goodness of fit for the Arrhenius plot."),
        ("slope", "float", "Fitted slope (= -Ea/R)."),
        ("intercept", "float", "Fitted intercept (= ln(A))."),
        ("arrhenius_equation", "str", "The fitted Arrhenius equation string."),
        ("analysis", "str", "Detailed analysis including temperature dependence interpretation."),
        ("plot_data", "list", "Data points for Arrhenius plot: [{1/T, ln(k), ...}]."),
    ]

    examples         = [
        {
            "code_input": {
                "temperatures_k": '298,308,318,328,338',
                "rate_constants": '0.00012,0.00025,0.00049,0.00094,0.00178',
                "r_gas": 8.314,
                "ea_unit": 'kJ/mol',
                "inverse_temperatures": '',
                "ln_rate_constants": ''
            },
            "text_input": {
                "input_params": '298,308,318,328,338 0.00012,0.00025,0.00049,0.00094,0.00178'
            },
            "output": {
                "activation_energy": 53.4,
                "pre_exponential_factor": 250000000000.0,
                "ea_unit": 'kJ/mol',
                "a_unit": 's^-1',
                "r_squared": 0.9999,
                "slope": -6422.6,
                "intercept": 26.55,
                "arrhenius_equation": 'k = 2.5e11 * exp(-53400/(RT)) s^-1',
                "analysis": 'Strong temperature dependence.',
                "plot_data": []
            }
        },
        {
            "code_input": {
                "temperatures_k": '600,650,700,750',
                "rate_constants": '0.035,0.12,0.35,0.92',
                "r_gas": 8.314,
                "ea_unit": 'kJ/mol',
                "inverse_temperatures": '',
                "ln_rate_constants": ''
            },
            "text_input": {
                "input_params": '600,650,700,750 0.035,0.12,0.35,0.92 kJ/mol'
            },
            "output": {
                "activation_energy": 82.7,
                "pre_exponential_factor": 32000000000.0,
                "ea_unit": 'kJ/mol',
                "a_unit": 's^-1',
                "r_squared": 0.998,
                "slope": -9950,
                "intercept": 24.0,
                "arrhenius_equation": 'k = 3.2e10 * exp(-82700/(RT)) s^-1',
                "analysis": 'High activation energy.',
                "plot_data": []
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _linear_regression(self, x_vals: List[float], y_vals: List[float]) -> tuple:
        n = len(x_vals)
        sx = sum(x_vals)
        sy = sum(y_vals)
        sxy = sum(xi * yi for xi, yi in zip(x_vals, y_vals))
        sx2 = sum(xi * xi for xi in x_vals)
        sy2 = sum(yi * yi for yi in y_vals)
        denom = n * sx2 - sx * sx
        if abs(denom) < 1e-15:
            return 0.0, sy / n, 1.0
        slope = (n * sxy - sx * sy) / denom
        intercept = (sy - slope * sx) / n
        y_mean = sy / n
        ss_tot = sum((yi - y_mean) ** 2 for yi in y_vals)
        ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x_vals, y_vals))
        r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        return slope, intercept, max(0, r_sq)

    def _run_base(
        self,
        temperatures_k: str,
        rate_constants: str,
        r_gas: float = 8.314,
        ea_unit: str = "kJ/mol",
        inverse_temperatures: str = "",
        ln_rate_constants: str = "",
    ) -> dict:
        # Prepare data
        if inverse_temperatures and ln_rate_constants:
            inv_t = [float(x.strip()) for x in inverse_temperatures.split(",")]
            ln_k = [float(x.strip()) for x in ln_rate_constants.split(",")]
        else:
            temps = [float(t.strip()) for t in temperatures_k.split(",")]
            ks = [float(k.strip()) for k in rate_constants.split(",")]
            if len(temps) != len(ks):
                raise ChemMCPError("temperatures_k and rate_constants must have same length.")
            if any(t <= 0 for t in temps):
                raise ChemMCPError("All temperatures must be positive.")
            if any(k <= 0 for k in ks):
                raise ChemMCPError("All rate constants must be positive.")
            inv_t = [1.0 / t for t in temps]
            ln_k = [math.log(k) for k in ks]

        if len(inv_t) < 2:
            raise ChemMCPError("Need at least 2 data points.")

        # Linear regression: ln(k) vs 1/T
        slope, intercept, r_sq = self._linear_regression(inv_t, ln_k)

        # Extract parameters
        # slope = -Ea/R → Ea = -slope * R
        ea_j_mol = -slope * r_gas
        if ea_j_mol < 0:
            logger.warning(f"Negative Ea ({ea_j_mol}) — check data or consider non-Arrhenius behavior.")

        use_kj = ea_unit.lower().startswith("kj")
        ea_out = ea_j_mol / 1000.0 if use_kj else ea_j_mol

        # intercept = ln(A) → A = exp(intercept)
        a_val = math.exp(intercept)

        # Determine units from first k value
        if ln_rate_constants or rate_constants:
            try:
                k_first = float(rate_constants.split(",")[0]) if not ln_rate_constants else math.exp(float(ln_rate_constants.split(",")[0]))
            except (ValueError, IndexError):
                k_first = 1.0
        else:
            k_first = 1.0

        a_unit_str = "(same unit as k)"
        ea_unit_str = ea_unit

        # Generate plot data
        plot_data = []
        for it, lk in zip(inv_t, ln_k):
            plot_data.append({"inv_t": round(it, 8), "ln_k": round(lk, 6)})

        # Formatted equation
        if use_kj:
            eq_str = f"k = {a_val:.3e} · exp(-{ea_out:.1f}/(R·T))"
        else:
            eq_str = f"k = {a_val:.3e} · exp(-{ea_out:.1f}/(R·T))"

        # Interpretation
        if ea_out < 40:
            interp = "Low activation energy: reaction proceeds readily at room temperature."
        elif ea_out < 80:
            interp = "Moderate activation energy: moderate temperature sensitivity."
        elif ea_out < 160:
            interp = "High activation energy: strong temperature dependence; elevated temperature needed."
        else:
            interp = "Very high activation energy: reaction requires significant thermal energy or catalysis."

        analysis = (
            f"Arrhenius Analysis ({len(inv_t)} data points):\n"
            f"Linear fit: ln(k) = {intercept:.4f} + ({slope:.2f})·(1/T)\n"
            f"Ea = {ea_out:.2f} {ea_unit_str}\n"
            f"A = {a_val:.4e}\n"
            f"R² = {r_sq:.6f}\n\n"
            f"Equation: {eq_str}\n\n"
            f"Interpretation: {interp}"
        )

        return {
            "activation_energy": round(ea_out, 2),
            "pre_exponential_factor": round(a_val, 4),
            "ea_unit": ea_unit_str,
            "a_unit": a_unit_str,
            "r_squared": round(r_sq, 6),
            "slope": round(slope, 4),
            "intercept": round(intercept, 4),
            "arrhenius_equation": eq_str,
            "analysis": analysis,
            "plot_data": plot_data,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            temps = parts[0]
            ks = parts[1]
            kwargs = {"temperatures_k": temps, "rate_constants": ks}
            idx = 2
            if idx < len(parts):
                try:
                    kwargs["r_gas"] = float(parts[idx]); idx += 1
                except ValueError:
                    kwargs["ea_unit"] = parts[idx]; idx += 1
            if idx < len(parts):
                kwargs["ea_unit"] = parts[idx]; idx += 1
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
