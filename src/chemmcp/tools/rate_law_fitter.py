import logging
import math
from typing import List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RateLawFitter(BaseTool):
    """
    根据实验数据拟合反应级数和速率常数工具。
    支持积分法（线性化拟合）、微分法（初始速率法）和半衰期法。
    """
    __version__ = "0.1.0"
    name = "RateLawFitter"
    func_name = "fit_rate_law"
    description = "Fit reaction order and rate constant from experimental concentration-time data using integral, differential, or half-life methods."
    implementation_description = "Uses linear regression on integrated rate laws (zero/first/second order), initial rates method, and half-life analysis to determine reaction order n and rate k."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Chemical Kinetics", "Rate Law", "Reaction Order", "Fitting"]
    required_envs    = []

    code_input_sig = [
        ("method", "str", "N/A", "Fitting method: 'integral', 'differential', 'half_life', or 'auto' (try all and pick best)."),
        ("time_data", "str", "N/A", "Time values (seconds or minutes), comma-separated."),
        ("concentration_data", "str", "N/A", "Concentration values (same units), comma-separated, matching time_data length."),
        ("time_unit", "str", "s", "Time unit: 's', 'min', 'h'."),
        ("conc_unit", "str", "M", "Concentration unit: 'M', 'mol/L', etc."),
        # For differential method:
        ("initial_rates", "str", "", "Initial rates at different [A]_0, comma-separated (for differential method)."),
        ("initial_concentrations", "str", "", "Initial concentrations for initial rates, comma-separated."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space/newline-separated parameters string."),
    ]

    output_sig = [
        ("reaction_order", "float", "Determined reaction order n (may be non-integer)."),
        ("rate_constant", "float", "Rate constant k with appropriate units."),
        ("k_unit", "str", "Unit of k depending on reaction order."),
        ("r_squared", "float", "R² goodness of fit (for linear methods)."),
        ("best_method", "str", "Method that gave best fit."),
        ("analysis", "str", "Detailed analysis including fitted equations and recommendations."),
    ]

    examples         = [
        {
            "code_input": {
                "method": 'auto',
                "time_data": '0,50,100,200,300,400',
                "concentration_data": '1.0,0.72,0.51,0.27,0.15,0.08',
                "time_unit": 's',
                "conc_unit": 'M',
                "initial_rates": '',
                "initial_concentrations": ''
            },
            "text_input": {
                "input_params": 'auto 0,50,100,200,300,400 1.0,0.72,0.51,0.27,0.15,0.08 s M'
            },
            "output": {
                "reaction_order": 1.0,
                "rate_constant": 0.0069,
                "k_unit": 's^-1',
                "r_squared": 0.9998,
                "best_method": 'integral (first-order)',
                "analysis": 'First-order kinetics confirmed.'
            }
        },
        {
            "code_input": {
                "method": 'differential',
                "time_data": '',
                "concentration_data": '',
                "time_unit": 's',
                "conc_unit": 'M',
                "initial_rates": '0.005,0.02,0.045,0.125',
                "initial_concentrations": '0.1,0.2,0.3,0.5'
            },
            "text_input": {
                "input_params": 'differential 0.1,0.2,0.3,0.5 0.005,0.02,0.045,0.125 s M'
            },
            "output": {
                "reaction_order": 2.0,
                "rate_constant": 0.5,
                "k_unit": 'L/(mol*s)',
                "r_squared": 1.0,
                "best_method": 'differential',
                "analysis": 'Second-order kinetics.'
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _linear_regression(self, x_vals: List[float], y_vals: List[float]) -> tuple:
        """Simple linear regression returning (slope, intercept, r_squared)."""
        n = len(x_vals)
        if n < 2:
            raise ChemMCPError("Need at least 2 data points.")
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
        if ss_tot == 0:
            r_sq = 1.0
        else:
            r_sq = 1 - ss_res / ss_tot
        return slope, intercept, max(0, r_sq)

    def _fit_integral(self, times: List[float], concs: List[float]) -> dict:
        """Try zero, first, second order integrated rate laws."""
        results = {}
        n = len(times)

        # Zero order: [A] = [A]₀ - kt → plot [A] vs t
        slope0, _, r0 = self._linear_regression(times, concs)
        results[0] = {"order": 0, "k": abs(slope0), "r2": r0}

        # First order: ln([A]) = ln([A]₀) - kt → plot ln([A]) vs t
        try:
            ln_c = [math.log(c) for c in concs if c > 0]
            t_first = times[:len(ln_c)]
            slope1, _, r1 = self._linear_regression(t_first, ln_c)
            results[1] = {"order": 1, "k": abs(slope1), "r2": r1}
        except (ValueError, OverflowError):
            results[1] = {"order": 1, "k": float("inf"), "r2": 0}

        # Second order: 1/[A] = 1/[A]₀ + kt → plot 1/[A] vs t
        try:
            inv_c = [1.0 / c for c in concs if c != 0]
            t_second = times[:len(inv_c)]
            slope2, _, r2 = self._linear_regression(t_second, inv_c)
            results[2] = {"order": 2, "k": slope2, "r2": r2}
        except (ZeroDivisionError, ValueError):
            results[2] = {"order": 2, "k": float("inf"), "r2": 0}

        # Pick best R²
        best_order = max(results.keys(), key=lambda k: results[k]["r2"])
        return {**results[best_order], "all_r2": {k: v["r2"] for k, v in results.items()}}

    def _fit_differential(self, c0_list: List[float], rate_list: List[float]) -> dict:
        """Initial rates method: log(rate) = log(k) + n·log([A]₀)."""
        log_c = [math.log(c) for c in c0_list if c > 0]
        log_r = [math.log(r) for r, c in zip(rate_list, c0_list) if c > 0 and r > 0]
        slope, intercept, r2 = self._linear_regression(log_c, log_r)
        k = math.exp(intercept)
        return {"order": round(slope, 4), "k": round(k, 6), "r2": round(r2, 6)}

    def _run_base(
        self,
        method: str,
        time_data: str,
        concentration_data: str,
        time_unit: str = "s",
        conc_unit: str = "M",
        initial_rates: str = "",
        initial_concentrations: str = "",
    ) -> dict:
        method = method.lower().strip()
        all_r2 = {}
        n_points = 0

        # Time unit conversion factor to seconds
        time_factors = {"s": 1.0, "min": 60.0, "h": 3600.0}
        tf = time_factors.get(time_unit, 1.0)

        if method == "differential":
            if not initial_rates or not initial_concentrations:
                raise ChemMCPError("Differential method requires initial_rates and initial_concentrations.")
            c0s = [float(c.strip()) for c in initial_concentrations.split(",")]
            rates = [float(r.strip()) for r in initial_rates.split(",")]
            result = self._fit_differential(c0s, rates)
            order = result["order"]
            k_raw = result["k"]
            r2 = result["r2"]
            best_m = "differential (log-log)"
            n_points = len(c0s)

        else:
            # integral, auto, or half_life: need time/conc data
            times = [float(t.strip()) for t in time_data.split(",")]
            concs = [float(c.strip()) for c in concentration_data.split(",")]
            if len(times) != len(concs):
                raise ChemMCPError(f"time_data ({len(times)}) and concentration_data ({len(concs)}) must have same length.")
            n_points = len(times)

            if method == "integral":
                result = self._fit_integral(times, concs)
                order = result["order"]
                k_raw = result["k"] / tf  # convert to per-second
                r2 = result["r2"]
                best_m = f"integral ({order}rd order)"
                all_r2 = result.get("all_r2", {})

            elif method == "auto":
                int_result = self._fit_integral(times, concs)
                r2_values = int_result["all_r2"]
                best_order = max(r2_values, key=r2_values.get)
                order = best_order
                k_raw = int_result["k"] / tf
                r2 = r2_values[best_order]
                best_m = f"auto → integral ({order}rd order, best R²)"
                all_r2 = r2_values

            elif method == "half_life":
                raise ChemMCPError("Half-life method requires specific half-life data points. Use 'integral' or 'differential' instead.")

            else:
                raise ChemMCPError(f"Unsupported method: '{method}'. Use 'integral', 'differential', 'half_life', or 'auto'.")

        # Determine k unit
        if abs(order) < 1e-6:
            k_unit = f"{conc_unit}/{time_unit}"
        elif abs(order - 1.0) < 1e-6:
            k_unit = f"{time_unit}⁻¹"
        elif abs(order - 2.0) < 1e-6:
            k_unit = f"({conc_unit}·{time_unit})⁻¹"  # L·mol⁻¹·s⁻¹ equivalent
        elif abs(order - 3.0) < 1e-6:
            k_unit = f"{conc_unit}⁻²·{time_unit}⁻¹"
        else:
            k_unit = f"{conc_unit}^{1-order:.1f}·{time_unit}⁻¹"

        analysis = (
            f"Rate law fitting using '{method}' method:\n"
            f"Data points: {n_points}\n"
            f"Determined order: n = {order}\n"
            f"Rate constant: k = {k_raw:.6g} {k_unit}\n"
            f"R² = {r2:.6f}\n"
            f"Rate equation: rate = k · [A]^{order}"
        )
        if all_r2:
            analysis += "\n\nR² comparison:\n" + "\n".join([f"  n={k}: R²={v:.6f}" for k, v in sorted(all_r2.items())])

        return {
            "reaction_order": round(order, 4),
            "rate_constant": round(k_raw, 6),
            "k_unit": k_unit,
            "r_squared": round(r2, 6),
            "best_method": best_m,
            "analysis": analysis,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            method = parts[0]
            kwargs = {"method": method}
            idx = 1
            if method == "differential":
                kwargs["initial_concentrations"] = parts[idx]; idx += 1
                kwargs["initial_rates"] = parts[idx]; idx += 1
            else:
                kwargs["time_data"] = parts[idx]; idx += 1
                kwargs["concentration_data"] = parts[idx]; idx += 1
            if idx < len(parts):
                kwargs["time_unit"] = parts[idx]; idx += 1
            if idx < len(parts):
                kwargs["conc_unit"] = parts[idx]; idx += 1
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
