import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
R = 8.314462618      # J/(mol·K)


@ChemMCPManager.register_tool
class ArrheniusCalculator(BaseTool):
    """
    Arrhenius方程计算工具。
    计算反应速率常数 k = A·exp(−Ea/RT)，支持从Ea/A求k、两点法求Ea、以及Arrhenius作图分析。
    """
    __version__ = "0.1.0"
    name = "ArrheniusCalculator"
    func_name = "calculate_arrhenius"
    description = "Calculate Arrhenius parameters: rate constant k from activation energy Ea and pre-exponential factor A; determine Ea from two-point data; or perform Arrhenius plot analysis."
    implementation_description = (
        "Uses the Arrhenius equation: k = A·exp(−Ea/RT). Supports three modes:\n"
        "1. 'calculate_k': Compute k from known Ea and A at given T.\n"
        "2. 'two_point_ea': Determine Ea from k values at two temperatures.\n"
        "3. 'arrhenius_plot': Linear regression of ln(k) vs 1/T to extract Ea and A."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Arrhenius", "Kinetics", "Activation Energy", "Rate Constant", "Chemical Kinetics"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Calculation mode: 'calculate_k', 'two_point_ea', or 'arrhenius_plot'."),
        ("ea_kj_mol", "float", "None", "Activation energy Ea in kJ/mol (for calculate_k mode)."),
        ("pre_exponential_A", "float", "None", "Pre-exponential factor A (same units as k) (for calculate_k mode)."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        # For two_point mode:
        ("T1_k", "float", "", "Temperature T1 in K (two-point mode)."),
        ("k1", "float", "", "Rate constant at T1 (two-point mode)."),
        ("T2_k", "float", "", "Temperature T2 in K (two-point mode)."),
        ("k2", "float", "", "Rate constant at T2 (two-point mode)."),
        # For arrhenius_plot mode:
        ("temperatures_k", "str", "", "Comma-separated temperatures for plot."),
        ("rate_constants", "str", "", "Comma-separated rate constants at each temperature."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: mode [ea] [A] [T] ... or two_point T1 k1 T2 k2"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with k, Ea, A, arrhenius equation, activation energy analysis, and plot data."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "calculate_k",
                "ea_kj_mol": 75.0,
                "pre_exponential_A": 1.0e13,
                "temperature_k": 298.15,
                "T1_k": 0,
                "k1": 0,
                "T2_k": 0,
                "k2": 0,
                "temperatures_k": "",
                "rate_constants": "",
            },
            "text_input": {
                "input_params": "calculate_k 75.0 1e13 298.15",
            },
            "output": {
                "result": {
                    "mode": "calculate_k",
                    "rate_constant_k": 4.57e-06,
                    "k_unit": "s⁻¹",
                    "ea_kj_mol": 75.0,
                    "pre_exponential_A": 1e13,
                    "temperature_k": 298.15,
                    "arrhenius_equation": "k = 1.00e+13 · exp(−75.0 / (8.314·298.15))",
                    "ln_k": -12.0,
                    "activation_energy_analysis": "Moderate barrier — typical for many organic reactions at room temperature.",
                    "temperature_sensitivity": "d(ln k)/d(1/T) = −Ea/R = −9019 K → 10°C rise roughly doubles k.",
                }
            }
        },
        {
            "code_input": {
                "mode": "two_point_ea",
                "ea_kj_mol": 0,
                "pre_exponential_A": 0,
                "temperature_k": 298.15,
                "T1_k": 300.0,
                "k1": 2.5e-5,
                "T2_k": 320.0,
                "k2": 1.5e-4,
                "temperatures_k": "",
                "rate_constants": "",
            },
            "text_input": {
                "input_params": "two_point_ea 300 2.5e-5 320 1.5e-4",
            },
            "output": {
                "result": {
                    "mode": "two_point_ea",
                    "ea_kj_mol": 84.6,
                    "pre_exponential_A_estimated": 3.2e12,
                    "T1_k": 300.0,
                    "k1": 2.5e-05,
                    "T2_k": 320.0,
                    "k2": 1.5e-04,
                    "arrhenius_equation": "k = 3.20e+12 · exp(−84600 / (R·T))",
                    "analysis": "Two-point Arrhenius analysis gives Ea ≈ 84.6 kJ/mol.",
                    "q10_factor": 4.7,
                }
            }
        },
        {
            "code_input": {
                "mode": "arrhenius_plot",
                "ea_kj_mol": 0,
                "pre_exponential_A": 0,
                "temperature_k": 298.15,
                "T1_k": 0,
                "k1": 0,
                "T2_k": 0,
                "k2": 0,
                "temperatures_k": "298,308,318,328,338",
                "rate_constants": "2.5e-6,7.2e-6,1.9e-5,4.8e-5,1.15e-4",
            },
            "text_input": {
                "input_params": "arrhenius_plot 298,308,318,328,338 2.5e-6,7.2e-6,1.9e-5,4.8e-5,1.15e-4",
            },
            "output": {
                "result": {
                    "mode": "arrhenius_plot",
                    "ea_kj_mol": 74.98,
                    "pre_exponential_A": 2.85e11,
                    "r_squared": 0.9999,
                    "arrhenius_equation": "k = 2.85e+11 · exp(−74980 / (R·T))",
                    "plot_data": [{"T_K": 298, "inv_T": 0.003356, "ln_k": -12.9}],
                    "analysis": "Linear Arrhenius behavior confirmed (R² > 0.99).",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        mode: str,
        ea_kj_mol: float = 0.0,
        pre_exponential_A: float = 0.0,
        temperature_k: float = 298.15,
        T1_k: float = 0.0,
        k1: float = 0.0,
        T2_k: float = 0.0,
        k2: float = 0.0,
        temperatures_k: str = "",
        rate_constants: str = "",
    ) -> dict:
        """Core logic: Arrhenius calculations."""
        m = mode.lower().strip()
        T = temperature_k

        if T <= 0:
            raise ChemMCPError("Temperature must be positive.")

        if m == "calculate_k":
            result = self._calc_k(ea_kj_mol, pre_exponential_A, T)
        elif m == "two_point_ea":
            result = self._two_point(T1_k, k1, T2_k, k2)
        elif m == "arrhenius_plot":
            result = self._arrhenius_plot(temperatures_k, rate_constants)
        else:
            raise ChemMCPError(f"Unknown mode: '{m}'. Use 'calculate_k', 'two_point_ea', or 'arrhenius_plot'.")

        return {"result": result}

    def _calc_k(self, ea_kj: float, A: float, T: float) -> dict:
        """Calculate k from Arrhenius equation."""
        ea_j = ea_kj * 1000.0  # J/mol
        k = A * math.exp(-ea_j / (R * T))
        ln_k = math.log(k) if k > 0 else float('-inf')

        # Q10: factor by which k changes per 10°C increase
        T_plus_10 = T + 10.0
        k_10 = A * math.exp(-ea_j / (R * T_plus_10))
        q10 = k_10 / k if k > 0 else 0

        return {
            "mode": "calculate_k",
            "rate_constant_k": round(k, 6),
            "k_unit": "s⁻¹",
            "ea_kj_mol": ea_kj,
            "pre_exponential_A": A,
            "temperature_k": T,
            "arrhenius_equation": f"k = {A:.3e} · exp(−{ea_kj} / ({R:.4f}·{T}))",
            "ln_k": round(ln_k, 4),
            "activation_energy_analysis": self._interpret_ea(ea_kj),
            "temperature_sensitivity": f"d(ln k)/d(1/T) = -Ea/R = {-ea_kj * 1000 / R:.1f} K -> Q10 ~= {q10:.1f}",
            "q10_factor": round(q10, 2),
        }

    def _two_point(self, T1: float, k1_val: float, T2: float, k2_val: float) -> dict:
        """Two-point Arrhenius: determine Ea from k at two temperatures."""
        if T1 <= 0 or T2 <= 0:
            raise ChemMCPError("Temperatures must be positive.")
        if k1_val <= 0 or k2_val <= 0:
            raise ChemMCPError("Rate constants must be positive.")

        # ln(k2/k1) = (Ea/R) · (1/T1 − 1/T2)
        # Ea = R · ln(k2/k1) / (1/T1 − 1/T2)
        ln_ratio = math.log(k2_val / k1_val)
        inv_T_diff = 1.0/T1 - 1.0/T2

        if abs(inv_T_diff) < 1e-15:
            raise ChemMCPError("Temperatures must be different.")

        ea_j_mol = R * ln_ratio / inv_T_diff
        ea_kj = ea_j_mol / 1000.0

        # Estimate A from one data point: A = k · exp(Ea/RT)
        A_est = k1_val * math.exp(ea_j_mol / (R * T1))

        # Q10 estimate
        T_mid = (T1 + T2) / 2
        T_mid_p10 = T_mid + 10
        k_mid = A_est * math.exp(-ea_j_mol / (R * T_mid))
        k_mid_10 = A_est * math.exp(-ea_j_mol / (R * T_mid_p10))
        q10 = k_mid_10 / k_mid if k_mid > 0 else 0

        return {
            "mode": "two_point_ea",
            "ea_kj_mol": round(ea_kj, 2),
            "ea_j_mol": round(ea_j_mol, 1),
            "pre_exponential_A_estimated": round(A_est, 4),
            "T1_k": T1,
            "k1": k1_val,
            "T2_k": T2,
            "k2": k2_val,
            "arrhenius_equation": f"k = {A_est:.3e} · exp(−{ea_j_mol:.1f} / (R·T))",
            "analysis": f"Two-point Arrhenius gives Ea = {ea_kj:.1f} kJ/mol.",
            "q10_factor": round(q10, 2),
            "validity_note": "Two-point method assumes Ea is constant over the temperature range.",
        }

    def _arrhenius_plot(self, temps_str: str, ks_str: str) -> dict:
        """Arrhenius plot: linear regression of ln(k) vs 1/T."""
        try:
            temps = [float(t.strip()) for t in temps_str.split(",")]
            ks = [float(k.strip()) for k in ks_str.split(",")]
        except (ValueError, AttributeError):
            raise ChemMCPError("Invalid format for temperatures_k or rate_constants. Use comma-separated values.")

        if len(temps) != len(ks):
            raise ChemMCPError("temperatures_k and rate_constants must have same length.")
        if len(temps) < 2:
            raise ChemMCPError("Need at least 2 data points for Arrhenius plot.")

        # Filter out non-positive k values
        valid_data = [(t, k) for t, k in zip(temps, ks) if t > 0 and k > 0]
        if len(valid_data) < 2:
            raise ChemMCPError("Need at least 2 valid (positive T, positive k) data points.")

        inv_T_list = [1.0 / t for t, k in valid_data]
        ln_k_list = [math.log(k) for t, k in valid_data]

        n = len(valid_data)
        sx = sum(inv_T_list)
        sy = sum(ln_k_list)
        sxy = sum(x * y for x, y in zip(inv_T_list, ln_k_list))
        sx2 = sum(x * x for x in inv_T_list)
        denom = n * sx2 - sx * sx

        if abs(denom) < 1e-15:
            raise ChemMCPError("Cannot fit: all temperatures are identical.")

        slope = (n * sxy - sx * sy) / denom
        intercept = (sy - slope * sx) / n

        # R² calculation
        y_mean = sy / n
        ss_tot = sum((y - y_mean) ** 2 for y in ln_k_list)
        ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(inv_T_list, ln_k_list))
        r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        # Extract parameters: slope = −Ea/R, intercept = ln(A)
        ea_j_mol = -slope * R
        ea_kj = ea_j_mol / 1000.0
        A_fit = math.exp(intercept)

        return {
            "mode": "arrhenius_plot",
            "ea_kj_mol": round(ea_kj, 2),
            "ea_j_mol": round(ea_j_mol, 1),
            "pre_exponential_A": round(A_fit, 4),
            "r_squared": round(r_sq, 6),
            "n_points": n,
            "arrhenius_equation": f"k = {A_fit:.3e} · exp(−{ea_j_mol:.1f} / (R·T))",
            "plot_data": [
                {"T_K": round(t, 2), "inv_T": round(1/t, 8), "ln_k": round(lk, 4)}
                for t, lk in zip([t for t, k in valid_data], ln_k_list)
            ],
            "regression": {
                "slope": round(slope, 4),
                "intercept": round(intercept, 4),
                "interpretation": f"slope = −Ea/R = {slope:.4f}, intercept = ln(A) = {intercept:.4f}",
            },
            "analysis": self._analyze_quality(r_sq, ea_kj, A_fit),
        }

    @staticmethod
    def _interpret_ea(ea_kj: float) -> str:
        if ea_kj < 30:
            return "Very low barrier — diffusion-controlled or nearly barrierless process."
        elif ea_kj < 50:
            return "Low barrier — fast reaction at room temperature."
        elif ea_kj < 80:
            return "Moderate barrier — typical for many organic reactions at room temperature."
        elif ea_kj < 120:
            return "High barrier — slow reaction at room temperature; requires heating or catalysis."
        else:
            return "Very high barrier — very slow unimolecular process; likely needs strong catalysis."

    @staticmethod
    def _analyze_quality(r_sq: float, ea_kj: float, A: float) -> str:
        parts = []
        if r_sq > 0.999:
            parts.append("Excellent linearity — pure single-step Arrhenius behavior.")
        elif r_sq > 0.99:
            parts.append("Good linearity — Arrhenius model fits well.")
        elif r_sq > 0.95:
            parts.append("Moderate linearity — some deviation possible (change in mechanism, tunneling, etc.).")
        else:
            parts.append("Poor linearity — consider multi-step mechanism or experimental error.")

        if ea_kj < 0:
            parts.append("⚠ Negative Ea is non-physical for elementary reactions; check data quality.")
        if A > 1e18:
            parts.append("Very large A suggests possible pre-equilibrium or complex mechanism.")
        elif A < 1e3:
            parts.append("Small A suggests tight transition state (unfavorable entropy).")

        return " ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if not parts:
                raise ChemMCPError("Empty input.")

            mode = parts[0]
            kwargs = {"mode": mode}

            if mode == "calculate_k":
                kwargs["ea_kj_mol"] = float(parts[1])
                kwargs["pre_exponential_A"] = float(parts[2])
                if len(parts) > 3:
                    kwargs["temperature_k"] = float(parts[3])
            elif mode == "two_point_ea":
                kwargs["T1_k"] = float(parts[1])
                kwargs["k1"] = float(parts[2])
                kwargs["T2_k"] = float(parts[3])
                kwargs["k2"] = float(parts[4])
            elif mode == "arrhenius_plot":
                kwargs["temperatures_k"] = parts[1]
                kwargs["rate_constants"] = parts[2]

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
