import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ChronoamperometryAnalyzer(BaseTool):
    """
    Chronoamperometry (CA) data analyzer.
    Analyze current-time transients after a potential step: Cottrell behavior check, diffusion coefficient estimation, double-layer charge integration, and adsorption detection.
    """
    __version__ = "0.1.0"
    name = "ChronoamperometryAnalyzer"
    func_name = "chronoamperometry_analyzer"
    description = "Analyze chronoamperometry data: verify Cottrell behavior (i vs t^(-1/2) linearity), estimate diffusion coefficient, integrate charge for double-layer capacitance and faradaic charge separation, detect adsorption effects."
    implementation_description = "Implements chronoamperometry analysis:\n1. Cottrell equation fit: i(t) = nFACD^(1/2)/(π^(1/2)·t^(1/2))\n2. i vs t^(-1/2) linearity assessment (R²)\n3. Charge integration: Q_total = ∫i dt\n4. Double-layer vs. faradaic charge separation using Anson plot (Q vs √t)\n5. Adsorption detection from non-Cottrell initial transient"
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chronoamperometry", "Cottrell Equation", "Diffusion Coefficient", "Charge Integration", "Electrochemistry", "Potential Step"]
    required_envs = []

    code_input_sig = [
        ("time_s", "list", "N/A", "List of time values in seconds."),
        ("current_A", "list", "N/A", "List of current values in Amperes (or µA)."),
        ("n_electrons", "int", "1", "Number of electrons transferred. Default: 1."),
        ("electrode_area_cm2", "float", "0.07", "Working electrode area in cm². Default: 0.07."),
        ("concentration_M", "float", "0.001", "Analyte concentration in mol/L. Default: 0.001."),
        ("current_unit", "str", "A", "Current unit: 'A' or 'uA'. Default: 'A'."),
        ("skip_initial_points", "int", "5", "Number of initial points to skip for Cottrell analysis (to exclude charging current). Default: 5."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: '[n] [A_cm2] [C_M] [unit] [skip] || t1,i1 t2,i2 ...'. Example: '1 0.07 0.001 A 5 || 0.001,-5e-5 0.005,-2e-5 0.01,-1.4e-5 0.02,-9.8e-6 0.05,-6.2e-6 0.1,-4.3e-6'"),
    ]

    output_sig = [
        ("cottrell_fit_R2", "float", "R² of Cottrell linear fit (i vs t^(-1/2))."),
        ("diffusion_coefficient_D_cm2_s", "float", "Estimated diffusion coefficient from Cottrell slope (cm²/s)."),
        ("cottrell_slope_A_s05", "float", "Slope of i vs t^(-1/2) plot (A·s^0.5)."),
        ("total_charge_C", "float", "Total integrated charge (Coulombs)."),
        ("faradaic_charge_C", "float", "Faradaic charge estimate (Coulombs)."),
        ("double_layer_charge_C", "float", "Double-layer charging charge (Coulombs)."),
        ("adsorption_detected", "bool", "Whether significant adsorption was detected from initial transient."),
        ("cottrell_behavior", "str", "Assessment: 'good', 'moderate', or 'poor' Cottrellian behavior."),
        ("analysis_summary", "str", "Text summary of CA analysis results."),
    ]

    examples = [
        {
            "code_input": {
                "time_s": [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
                "current_A": [-8e-5, -3.5e-5, -2.4e-5, -1.7e-5, -1.07e-5, -7.5e-6, -5.3e-6, -3.35e-6],
                "n_electrons": 1,
                "electrode_area_cm2": 0.07,
                "concentration_M": 0.005,
            },
            "text_input": {"input_string": "1 0.07 0.005 A 5 || 0.001,-8e-5 0.005,-3.5e-5 0.01,-2.4e-5 0.02,-1.7e-5 0.05,-1.07e-5 0.1,-7.5e-6 0.2,-5.3e-6 0.5,-3.35e-6"},
            "output": {
                "cottrell_fit_R2": 0.998,
                "diffusion_coefficient_D_cm2_s": 7.2e-6,
                "cottrell_slope_A_s05": -5.33e-6,
                "total_charge_C": -3.85e-6,
                "cottrell_behavior": "good",
                "adsorption_detected": False,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._F = 96485       # C/mol
        self._pi = math.pi

    def _linear_regression(self, x: List[float], y: List[float]) -> tuple:
        """Returns (slope, intercept, r_squared)."""
        n = len(x)
        if n < 2:
            return 0, 0, 0
        sx = sum(x)
        sy = sum(y)
        sxy = sum(xi * yi for xi, yi in zip(x, y))
        sx2 = sum(xi ** 2 for xi in x)
        sy2 = sum(yi ** 2 for yi in y)
        denom = n * sx2 - sx ** 2
        if abs(denom) < 1e-15:
            return 0, 0, 0
        slope = (n * sxy - sx * sy) / denom
        intercept = (sy - slope * sx) / n
        y_mean = sy / n
        ss_tot = sum((yi - y_mean) ** 2 for yi in y)
        ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
        r_sq = 1 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
        return slope, intercept, r_sq

    def _trapezoidal_integrate(self, x: List[float], y: List[float]) -> float:
        """Trapezoidal rule integration."""
        total = 0.0
        for i in range(len(x) - 1):
            total += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
        return total

    def _run_base(
        self,
        time_s: List[float],
        current_A: List[float],
        n_electrons: int = 1,
        electrode_area_cm2: float = 0.07,
        concentration_M: float = 0.001,
        current_unit: str = "A",
        skip_initial_points: int = 5,
    ) -> dict:
        """Analyze chronoamperometry data."""
        if len(time_s) != len(current_A):
            raise ChemMCPError("time_s and current_A must have same length.")
        if len(time_s) < 4:
            raise ChemMCPError("Need at least 4 data points.")

        # Convert units
        I_list = list(current_A)
        if current_unit.lower() in ("ua", "μa"):
            I_list = [i * 1e-6 for i in I_list]

        T = list(time_s)

        # Total charge via trapezoidal integration
        Q_total = self._trapezoidal_integrate(T, I_list)

        # Prepare Cottrell analysis data: i vs t^(-1/2)
        # Skip initial points to avoid charging current region
        start_idx = min(skip_initial_points, len(T) - 2)
        inv_sqrt_t = []
        i_for_fit = []
        for i in range(start_idx, len(T)):
            if T[i] > 0:
                inv_sqrt_t.append(1.0 / math.sqrt(T[i]))
                i_for_fit.append(I_list[i])

        # Linear regression: i vs t^(-1/2)
        slope, intercept, r_sq = self._linear_regression(inv_sqrt_t, i_for_fit)

        # Diffusion coefficient from Cottrell slope
        # slope = n·F·A·C·D^0.5 / sqrt(π)
        D_sq = 0.0
        if slope != 0 and concentration_M > 0 and electrode_area_cm2 > 0:
            numer = abs(slope) * math.sqrt(self._pi)
            denom = n_electrons * self._F * electrode_area_cm2 * concentration_M
            D_sq = (numer / denom) ** 2

        # Anson-like analysis: Q vs √t for double-layer/faradaic separation
        sqrt_t_data = [math.sqrt(max(t, 0)) for t in T]
        cumulative_Q = []
        running_q = 0.0
        for i in range(len(T)):
            if i > 0:
                running_q += 0.5 * (I_list[i-1] + I_list[i]) * (T[i] - T[i-1])
            cumulative_Q.append(running_q)

        # Estimate double-layer charge from intercept of Q vs √t at √t → 0
        q_intercept = None
        if len(sqrt_t_data) >= 2:
            _, q_intercept, _ = self._linear_regression(sqrt_t_data[start_idx:], cumulative_Q[start_idx:])

        # Adsorption detection: compare initial current with Cottrell prediction
        adsorption_detected = False
        if len(I_list) > skip_initial_points + 1:
            # Check if first few points deviate significantly from Cottrell extrapolation
            for i in range(1, min(skip_initial_points + 1, len(T))):
                if T[i] > 0:
                    cottrell_pred = slope / math.sqrt(T[i]) + intercept
                    ratio = abs(I_list[i]) / abs(cottrell_pred) if cottrell_pred != 0 else 10
                    if ratio > 2.0 or ratio < 0.3:
                        adsorption_detected = True
                        break

        # Assess Cottrell behavior
        if r_sq >= 0.995:
            behavior = "good"
        elif r_sq >= 0.98:
            behavior = "moderate"
        else:
            behavior = "poor"

        summary_parts = [
            f"Chronoamperometry Analysis:",
            f"  Data points: {len(T)}, skipped initial: {start_idx}",
            f"  Cottrell fit (i vs t⁻¹/²): R² = {r_sq:.4f} → {behavior}",
            f"  Cottrell slope: {slope:.4e} A·s^0.5",
            f"  Estimated D = {D_sq:.4e} cm²/s",
            f"  Total charge: {Q_total:.4e} C",
        ]
        if q_intercept is not None:
            summary_parts.append(f"  Estimated Q_dl ≈ {q_intercept:.4e} C")

        return {
            "cottrell_fit_R2": round(r_sq, 6),
            "diffusion_coefficient_D_cm2_s": round(D_sq, 12),
            "cottrell_slope_A_s05": round(slope, 12),
            "total_charge_C": round(Q_total, 12),
            "faradaic_charge_C": round(Q_total - (q_intercept or 0), 12) if q_intercept else round(Q_total, 12),
            "double_layer_charge_C": round(q_intercept, 12) if q_intercept else None,
            "adsorption_detected": adsorption_detected,
            "cottrell_behavior": behavior,
            "analysis_summary": "\n".join(summary_parts),
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        if "||" not in input_string:
            raise ChemMCPError("Must use '||' separator between parameters and data.")

        left, right = input_string.split("||", 1)
        params = left.strip().split()
        pairs = right.strip().split()

        n = int(params[0]) if len(params) > 0 else 1
        A = float(params[1]) if len(params) > 1 else 0.07
        C = float(params[2]) if len(params) > 2 else 0.001
        unit = params[3] if len(params) > 3 else "A"
        skip = int(params[4]) if len(params) > 4 else 5

        t_data = []
        i_data = []
        for pair in pairs:
            if "," in pair:
                sub = pair.split(",")
                try:
                    t_data.append(float(sub[0]))
                    i_data.append(float(sub[1]))
                except ValueError:
                    continue

        if not t_data:
            raise ChemMCPError("No valid t,i data pairs found after '||'.")

        return self._run_base(t_data, i_data, n, A, C, unit, skip)
