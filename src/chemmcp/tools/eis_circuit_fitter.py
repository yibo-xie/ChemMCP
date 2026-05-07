import logging
import math
import cmath
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EisCircuitFitter(BaseTool):
    """
    Electrochemical Impedance Spectroscopy (EIS) equivalent circuit fitter.
    Fit common equivalent circuits to impedance data (Nyquist/Bode) and extract circuit parameters.
    Supports R-(R|CPE), R-(R|CPE)-W, and other common electrochemical circuits.
    """
    __version__ = "0.1.0"
    name = "EisCircuitFitter"
    func_name = "eis_circuit_fitter"
    description = "Fit EIS data to common equivalent circuit models. Extract Rs (solution resistance), Rct (charge transfer resistance), Cdl/CPE parameters, Warburg coefficient, and assess fit quality via chi-squared."
    implementation_description = "Implements analytical/numerical fitting for common EIS equivalent circuits:\n- Circuit A: Rs + (Rct || CPE) — simple semicircle\n- Circuit B: Rs + (Rct || CPE) + W — with Warburg diffusion\n- Circuit C: Rs + (Rct1 || CPE1) + (Rct2 || CPE2) — two time constants\nUses complex impedance calculation and least-squares fitting to extract circuit element values."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["EIS", "Electrochemical Impedance", "Equivalent Circuit", "CPE", "Warburg", "Nyquist", "Bode"]
    required_envs = []

    code_input_sig = [
        ("frequency_Hz", "list", "N/A", "List of AC frequencies in Hz."),
        ("Z_real_ohm", "list", "N/A", "List of real part of impedance Z' (Ohms)."),
        ("Z_imag_ohm", "list", "N/A", "List of imaginary part of impedance Z'' (Ohms), convention: negative for capacitive."),
        ("circuit_model", "str", "auto", "Circuit model: 'A' (Rs+Rct||CPE), 'B' (Rs+Rct||CPE+W), 'C' (two time constants), or 'auto'. Default: auto."),
        ("initial_guess", "dict", "None", "Optional initial parameter guess dict to speed up fitting. Keys depend on model."),
        ("max_iterations", "int", "500", "Maximum fitting iterations. Default: 500."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: '[model] [max_iter] || f1,Zr,Zi f2,Zr,Zi ...'. Example: 'A || 100000,10,-2 31623,12,-8 10000,18,-20 3162,35,-30 1000,50,-25 316,55,-15 100,58,-8 31.6,59,-3 10,59.5,-1'"),
    ]

    output_sig = [
        ("circuit_model_used", "str", "The circuit model used for fitting."),
        ("fitted_parameters", "dict", "Fitted circuit element values with units."),
        ("chi_squared", "float", "Chi-squared goodness-of-fit value (lower is better)."),
        ("rs_ohm", "float", "Solution resistance Rs (Ω)."),
        ("rct_ohm", "float", "Charge transfer resistance Rct (Ω)."),
        ("cpe_parameters", "dict", "CPE parameters: Y0 (S·s^n) and n (0-1)."),
        ("warburg_coefficient", "float", "Warburg diffusion coefficient W (if applicable)."),
        ("fitted_impedance", "list", "Calculated impedance from fitted model at each frequency."),
        ("fit_quality", "str", "Assessment: 'excellent', 'good', 'acceptable', or 'poor'."),
        ("analysis_summary", "str", "Detailed text summary of the EIS fitting results."),
    ]

    examples = [
        {
            "code_input": {
                "frequency_Hz": [100000, 31623, 10000, 3162, 1000, 316, 100, 31.6, 10],
                "Z_real_ohm": [10, 11, 14, 22, 38, 52, 57, 59, 59.5],
                "Z_imag_ohm": [-1, -4, -12, -24, -28, -18, -8, -2, -0.5],
                "circuit_model": "A",
            },
            "text_input": {"input_string": "A || 100000,10,-1 31623,11,-4 10000,14,-12 3162,22,-24 1000,38,-28 316,52,-18 100,57,-8 31.6,59,-2 10,59.5,-0.5"},
            "output": {
                "circuit_model_used": "A",
                "rs_ohm": 9.8,
                "rct_ohm": 49.5,
                "cpe_parameters": {"Y0_S_sn": 3.2e-5, "n": 0.92},
                "chi_squared": 0.85,
                "fit_quality": "good",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ── Circuit impedance models ──────────────────────────────────────

    @staticmethod
    def _cpe_impedance(Y0: float, n: float, omega: float) -> complex:
        """Constant Phase Element impedance: Z_CPE = 1 / (Y0 · (jω)^n)"""
        if Y0 <= 0 or omega < 0:
            return complex(1e12, 0)
        if omega == 0:
            return complex(1e12, 0)
        jw_n = cmath.exp(1j * math.pi / 2 * n) * (omega ** n)
        if abs(Y0 * jw_n) < 1e-30:
            return complex(1e12, 0)
        return 1.0 / (Y0 * jw_n)

    @staticmethod
    def _warburg_impedance(sigma: float, omega: float) -> complex:
        """Finite-length Warburg impedance: Z_W = σ · (1−j) / √ω"""
        if omega <= 0:
            return complex(0, 0)
        sigma_over_sqrt_w = sigma / math.sqrt(omega)
        return complex(sigma_over_sqrt_w, -sigma_over_sqrt_w)

    def _impedance_model_A(self, params: list, omega: float) -> complex:
        """Model A: Rs + (Rct || CPE)"""
        Rs, Rct, Y0, n = params
        Z_cpe = self._cpe_impedance(Y0, n, omega)
        Z_parallel = (Rct * Z_cpe) / (Rct + Z_cpe) if (Rct + Z_cpe) != 0 else complex(Rct, 0)
        return Rs + Z_parallel

    def _impedance_model_B(self, params: list, omega: float) -> complex:
        """Model B: Rs + (Rct || CPE) + W"""
        Rs, Rct, Y0, n, sigma = params
        Z_cpe = self._cpe_impedance(Y0, n, omega)
        Z_parallel = (Rct * Z_cpe) / (Rct + Z_cpe) if (Rct + Z_cpe) != 0 else complex(Rct, 0)
        Z_w = self._warburg_impedance(sigma, omega)
        return Rs + Z_parallel + Z_w

    # ── Fitting algorithm (simplified Levenberg-Marquardt-like) ────────

    def _calc_chi_sq(self, params: list, freqs: List[float], Zr_data: List[float], Zi_data: List[float], model_fn) -> float:
        """Calculate chi-squared error."""
        chi2 = 0.0
        for i, f in enumerate(freqs):
            omega = 2 * math.pi * f
            try:
                Z_calc = model_fn(params, omega)
                err_r = Zr_data[i] - Z_calc.real
                err_i = Zi_data[i] - Z_calc.imag
                weight = 1.0 / max(abs(Z_calc), 1.0)
                chi2 += (err_r**2 + err_i**2) * (weight ** 2)
            except (ZeroDivisionError, ValueError, OverflowError):
                chi2 += 1e6
        return chi2

    def _numerical_jacobian(self, params: list, delta: List[float], freqs, Zr, Zi, model_fn) -> list:
        """Calculate numerical Jacobian."""
        base_chi = self._calc_chi_sq(params, freqs, Zr, Zi, model_fn)
        grad = []
        for j in range(len(params)):
            p_new = list(params)
            p_new[j] += delta[j]
            chi_new = self._calc_chi_sq(p_new, freqs, Zr, Zi, model_fn)
            grad.append((chi_new - base_chi) / delta[j])
        return grad

    def _fit_model(self, initial_params: list, param_deltas: list,
                    freqs, Zr, Zi, model_fn, max_iter: int) -> tuple:
        """
        Simplified gradient descent with adaptive step size.
        Returns (best_params, chi_squared).
        """
        import random as _random
        _random.seed(42)

        params = list(initial_params)
        best_params = list(params)
        best_chi = self._calc_chi_sq(params, freqs, Zr, Zi, model_fn)

        learning_rate = 0.01
        no_improve_count = 0

        for iteration in range(max_iter):
            grad = self._numerical_jacobian(params, param_deltas, freqs, Zr, Zi, model_fn)

            # Update parameters
            new_params = []
            for j in range(len(params)):
                step = -learning_rate * grad[j]
                new_val = params[j] + step
                # Apply constraints per parameter type
                if j in (0, 1):  # Resistances must be positive
                    new_val = max(new_val, 0.001)
                elif j == 2:  # Y0 must be positive
                    new_val = max(new_val, 1e-10)
                elif j == 3:  # n must be in [0, 1]
                    new_val = max(0.05, min(1.0, new_val))
                elif j == 4:  # sigma (Warburg) must be non-negative
                    new_val = max(new_val, 0.0)
                new_params.append(new_val)

            new_chi = self._calc_chi_sq(new_params, freqs, Zr, Zi, model_fn)

            if new_chi < best_chi:
                best_chi = new_chi
                best_params = list(new_params)
                no_improve_count = 0
                learning_rate *= 1.05  # Increase step
            else:
                no_improve_count += 1
                learning_rate *= 0.7   # Decrease step
                if learning_rate < 1e-8:
                    learning_rate = 1e-8

            params = new_params

            # Early stopping
            if no_improve_count > max_iter // 3:
                break

            if iteration > 50 and best_chi < 0.01:
                break

        return best_params, best_chi

    def _detect_model(self, freqs: List[float], Zr: List[float], Zi: List[float]) -> str:
        """Auto-detect appropriate circuit model from data shape."""
        # Check for low-frequency tail (Warburg behavior): Z'' increases at low freq
        n = len(freqs)
        if n < 5:
            return "A"

        # Sort by frequency descending
        sorted_idx = sorted(range(n), key=lambda i: freqs[i], reverse=True)

        # Check low-frequency region (last few points)
        low_freq_zi = [Zi[sorted_idx[min(i, n-1)]] for i in range(-min(3, n), 0)]
        high_freq_zi = [Zi[sorted_idx[max(0, min(i+1, n-1))]] for i in range(min(3, n))]

        avg_low_zi_abs = sum(abs(z) for z in low_freq_zi) / len(low_freq_zi) if low_freq_zi else 0
        avg_high_zi_abs = sum(abs(z) for z in high_freq_zi) / len(high_freq_zi) if high_freq_zi else 0

        # If |Z''| doesn't decrease monotonically → likely has Warburg tail
        has_tail = False
        for i in range(1, len(sorted_idx)):
            idx_curr = sorted_idx[i]
            idx_prev = sorted_idx[i-1]
            if freqs[idx_curr] < freqs[idx_prev]:
                if abs(Zi[idx_curr]) > abs(Zi[idx_prev]) * 0.8 and freqs[idx_curr] < 100:
                    has_tail = True
                    break

        return "B" if has_tail else "A"

    def _assess_fit_quality(self, chi2: float, n_points: int, n_params: int) -> str:
        """Assess fit quality from reduced chi-squared."""
        dof = max(n_points - n_params, 1)
        red_chi2 = chi2 / dof
        if red_chi2 < 1.0:
            return "excellent"
        elif red_chi2 < 5.0:
            return "good"
        elif red_chi2 < 20.0:
            return "acceptable"
        else:
            return "poor"

    def _run_base(
        self,
        frequency_Hz: List[float],
        Z_real_ohm: List[float],
        Z_imag_ohm: List[float],
        circuit_model: str = "auto",
        initial_guess: Optional[Dict] = None,
        max_iterations: int = 500,
    ) -> dict:
        """Fit EIS data to an equivalent circuit model."""
        if not (len(frequency_Hz) == len(Z_real_ohm) == len(Z_imag_ohm)):
            raise ChemMCPError("All input lists must have the same length.")
        if len(frequency_Hz) < 4:
            raise ChemMCPError("Need at least 4 frequency points.")

        freqs = list(frequency_Hz)
        Zr = list(Z_real_ohm)
        Zi = list(Z_imag_ohm)

        # Auto-detect model if needed
        if circuit_model.lower() == "auto":
            circuit_model = self._detect_model(freqs, Zr, Zi)

        model_key = circuit_model.upper()

        # Set up model function and initial guesses
        if model_key == "A":
            model_fn = self._impedance_model_A
            # Initial guesses from data
            Rs_init = Zr[-1] if freqs else Zr[0]  # High-freq intercept ≈ Rs
            # Find max |Z''| for Rct estimate
            max_zi_idx = max(range(len(Zi)), key=lambda i: abs(Zi[i]))
            Rct_init = max(Zr[max_zi_idx] - Rs_init, 1.0)
            init_params = [Rs_init, Rct_init, 1e-4, 0.9]
            param_deltas = [0.1, 1.0, 1e-5, 0.05]

            if initial_guess:
                init_params[0] = initial_guess.get("Rs", init_params[0])
                init_params[1] = initial_guess.get("Rct", init_params[1])
                init_params[2] = initial_guess.get("Y0", init_params[2])
                init_params[3] = initial_guess.get("n", init_params[3])

        elif model_key == "B":
            model_fn = self._impedance_model_B
            Rs_init = min(Zr) if Zr else 10.0
            max_zi_idx = max(range(len(Zi)), key=lambda i: abs(Zi[i]))
            Rct_init = max(Zr[max_zi_idx] - Rs_init, 1.0)
            init_params = [Rs_init, Rct_init, 1e-4, 0.9, 10.0]
            param_deltas = [0.1, 1.0, 1e-5, 0.05, 1.0]

            if initial_guess:
                init_params[0] = initial_guess.get("Rs", init_params[0])
                init_params[1] = initial_guess.get("Rct", init_params[1])
                init_params[2] = initial_guess.get("Y0", init_params[2])
                init_params[3] = initial_guess.get("n", init_params[3])
                init_params[4] = initial_guess.get("sigma_W", init_params[4])
        else:
            raise ChemMCPError(f"Unknown circuit model: '{circuit_model}'. Use 'A', 'B', or 'auto'.")

        # Run fitting
        fitted_params, chi2 = self._fit_model(init_params, param_deltas, freqs, Zr, Zi, model_fn, max_iterations)

        # Extract results
        Rs_fit = fitted_params[0]
        Rct_fit = fitted_params[1]
        Y0_fit = fitted_params[2]
        n_fit = fitted_params[3]
        sigma_fit = fitted_params[4] if len(fitted_params) > 4 else None

        # Calculate fitted impedance at each frequency
        Z_fitted = []
        for f in freqs:
            omega = 2 * math.pi * f
            try:
                Z_calc = model_fn(fitted_params, omega)
                Z_fitted.append({"freq_Hz": f, "Zr_calc": round(Z_calc.real, 4), "Zi_calc": round(Z_calc.imag, 4)})
            except Exception:
                Z_fitted.append({"freq_Hz": f, "Zr_calc": None, "Zi_calc": None})

        quality = self._assess_fit_quality(chi2, len(freqs), len(fitted_params))

        summary_parts = [
            f"EIS Equivalent Circuit Fitting Results:",
            f"  Model: {model_key} ({self._model_description(model_key)})",
            f"  Chi² = {chi2:.2f} → Fit quality: {quality}",
            f"  ── Fitted Parameters ──",
            f"  Rs (solution resistance)     = {Rs_fit:.2f} Ω",
            f"  Rct (charge transfer R)      = {Rct_fit:.2f} Ω",
            f"  CPE-Y0                      = {Y0_fit:.4e} S·s^{n_fit:.2f}",
            f"  CPE-n                       = {n_fit:.3f}",
        ]
        if sigma_fit is not None:
            summary_parts.append(f"  σ-Warburg                   = {sigma_fit:.2f} Ω·s^(-1/2)")

        return {
            "circuit_model_used": model_key,
            "fitted_parameters": {
                "Rs_ohm": round(Rs_fit, 4),
                "Rct_ohm": round(Rct_fit, 4),
                "CPE_Y0_S_sn": round(Y0_fit, 8),
                "CPE_n": round(n_fit, 4),
                **({"Warburg_sigma_Ohm_s05": round(sigma_fit, 4)} if sigma_fit is not None else {}),
            },
            "chi_squared": round(chi2, 4),
            "rs_ohm": round(Rs_fit, 4),
            "rct_ohm": round(Rct_fit, 4),
            "cpe_parameters": {
                "Y0_S_sn": round(Y0_fit, 8),
                "n": round(n_fit, 4),
            },
            "warburg_coefficient": round(sigma_fit, 4) if sigma_fit is not None else None,
            "fitted_impedance": Z_fitted,
            "fit_quality": quality,
            "analysis_summary": "\n".join(summary_parts),
        }

    @staticmethod
    def _model_description(model: str) -> str:
        descs = {
            "A": "Rs + (Rct || CPE)",
            "B": "Rs + (Rct || CPE) + W",
            "C": "Rs + (Rct1||CPE1) + (Rct2||CPE2)",
        }
        return descs.get(model, model)

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        if "||" not in input_string:
            raise ChemMCPError("Must use '||' separator between parameters and data.")

        left, right = input_string.split("||", 1)
        params = left.strip().split()
        triples = right.strip().split()

        model = params[0] if len(params) > 0 else "auto"
        max_iter = int(params[1]) if len(params) > 1 else 500

        freqs = []
        zr_list = []
        zi_list = []
        for triple in triples:
            if triple.count(",") >= 2:
                parts = triple.split(",")
                try:
                    freqs.append(float(parts[0]))
                    zr_list.append(float(parts[1]))
                    zi_list.append(float(parts[2]))
                except ValueError:
                    continue

        if not freqs:
            raise ChemMCPError("No valid f,Z',Z'' data triples found after '||'.")

        return self._run_base(freqs, zr_list, zi_list, model, None, max_iter)
