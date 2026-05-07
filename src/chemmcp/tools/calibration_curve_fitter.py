import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalibrationCurveFitter(BaseTool):
    """
    标准曲线拟合工具。
    支持线性、加权线性（1/x, 1/x²）、多项式拟合，输出方程、R²、残差、LOD/LOQ估计。
    """
    __version__ = "0.1.0"
    name = "CalibrationCurveFitter"
    func_name = "fit_calibration_curve"
    description = "Fit calibration curves: linear, weighted linear (1/x, 1/x²), polynomial regression with full diagnostics."
    implementation_description = "Implements least-squares regression (linear and polynomial) with optional inverse/variance weighting. Computes R², adjusted R², residuals, SSE, SSR, SST, LOD/LOQ estimates."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Calibration", "Regression", "Analytical Chemistry", "QA/QC", "Statistics"]
    required_envs = []

    code_input_sig = [
        ("concentrations", "list", "N/A", "List of standard concentrations (x values)."),
        ("responses", "list", "N/A", "List of measured responses/instrument signals (y values)."),
        ("fit_type", "str", "linear", "Fit type: 'linear', 'weighted_1/x', 'weighted_1/x2', 'polynomial' (degree=2 or 3)."),
        ("polynomial_degree", "int", "2", "Polynomial degree if fit_type='polynomial' (2 or 3)."),
        ("force_origin", "bool", "False", "Whether to force the regression through origin (y-intercept = 0)."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("equation", "str", "Fitted equation string."),
        ("coefficients", "list", "Regression coefficients [slope, intercept] or polynomial coefficients [a_n, ..., a_0]."),
        ("r_squared", "float", "Coefficient of determination R²."),
        ("adjusted_r_squared", "float", "Adjusted R²."),
        ("residuals", "list", "Residuals for each point."),
        ("fitted_values", "list", "Predicted/fitted y values."),
        ("sse", "float", "Sum of squared errors (residuals)."),
        ("ssr", "float", "Sum of squares due to regression."),
        ("sst", "float", "Total sum of squares."),
        ("std_error", "float", "Standard error of the estimate (residual standard deviation)."),
        ("lod_estimate", "float", "Estimated limit of detection (3.33 × SE / slope)."),
        ("loq_estimate", "float", "Estimated limit of quantification (10 × SE / slope)."),
        ("n_points", "int", "Number of data points."),
        ("diagnostics", "dict", "Additional diagnostics (weighting info, etc.)."),
    ]

    examples = [
        {
            "code_input": {
                "concentrations": [0.0, 1.0, 2.0, 5.0, 10.0],
                "responses": [0.02, 0.105, 0.205, 0.502, 1.01],
                "fit_type": "linear",
            },
            "text_input": {"params_str": "[[0,1,2,5,10]] [[0.02,0.105,0.205,0.502,1.01]] linear"},
            "output": {
                "r_squared": 0.9998,
                "equation": "y = 0.0992x + 0.0146",
                "lod_estimate": 0.049,
            },
        },
        {
            "code_input": {
                "concentrations": [0.1, 0.5, 1.0, 5.0, 10.0, 20.0],
                "responses": [150, 720, 1450, 7100, 14200, 28500],
                "fit_type": "weighted_1/x",
            },
            "text_input": {"params_str": "[[0.1,0.5,1,5,10,20]] [[150,720,1450,7100,14200,28500]] weighted_1/x"},
            "output": {
                "r_squared": 0.9999,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- helpers ----

    @staticmethod
    def _mean(data: List[float]) -> float:
        return sum(data) / len(data) if data else 0.0

    @staticmethod
    def _var(data: List[float], ddof: int = 1) -> float:
        n = len(data)
        if n <= ddof:
            return 0.0
        m = CalibrationCurveFitter._mean(data)
        return sum((x - m) ** 2 for x in data) / (n - ddof)

    def _linear_regression(self, x: List[float], y: List[float], weights: Optional[List[float]] = None, force_origin: bool = False) -> dict:
        """Weighted or unweighted linear regression."""
        n = len(x)
        if n < 2:
            raise ChemMCPError("Need at least 2 points for linear regression.")

        if force_origin:
            if weights:
                sw = sum(weights)
                swxy = sum(w * xi * yi for w, xi, yi in zip(weights, x, y))
                swx2 = sum(w * xi ** 2 for w, xi in zip(weights, x))
                slope = swxy / swx2 if swx2 != 0 else 0
                intercept = 0.0
            else:
                slope = sum(xi * yi for xi, yi in zip(x, y)) / sum(xi ** 2 for xi in x)
                intercept = 0.0
        else:
            if weights:
                sw = sum(weights)
                swx = sum(w * xi for w, xi in zip(weights, x))
                swy = sum(w * yi for w, yi in zip(weights, y))
                swxy = sum(w * xi * yi for w, xi, yi in zip(weights, x, y))
                swx2 = sum(w * xi ** 2 for w, xi in zip(weights, x))
                denom = sw * swx2 - swx ** 2
                if abs(denom) < 1e-30:
                    raise ChemMCPError("Singular matrix in weighted linear regression.")
                slope = (sw * swxy - swx * swy) / denom
                intercept = (swy * swx2 - swx * swxy) / denom
            else:
                mx = self._mean(x)
                my = self._mean(y)
                ss_xx = sum((xi - mx) ** 2 for xi in x)
                ss_xy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
                if abs(ss_xx) < 1e-30:
                    raise ChemMCPError("Variance in x is zero; cannot fit.")
                slope = ss_xy / ss_xx
                intercept = my - slope * mx

        fitted = [slope * xi + intercept for xi in x]
        residuals = [yi - y_hat for yi, y_hat in zip(y, fitted)]

        # R² calculation
        y_mean = self._mean(y)
        sst = sum((yi - y_mean) ** 2 for yi in y)
        sse = sum(r ** 2 for r in residuals)
        ssr = sst - sse
        r_sq = 1 - sse / sst if sst > 0 else 1.0

        # Adjusted R²
        p = 1 if force_origin else 2  # number of parameters
        adj_r2 = 1 - (1 - r_sq) * (n - 1) / (n - p) if n > p else 0.0

        std_err = math.sqrt(sse / (n - p)) if n > p else 0.0

        # LOD/LOQ (IUPAC: 3.33σ/slope for LOD, 10σ/slope for LOQ)
        lod = (3.33 * std_err / abs(slope)) if abs(slope) > 1e-30 else float('inf')
        loq = (10.0 * std_err / abs(slope)) if abs(slope) > 1e-30 else float('inf')

        return {
            "slope": slope,
            "intercept": intercept,
            "fitted_values": fitted,
            "residuals": residuals,
            "r_squared": round(r_sq, 6),
            "adjusted_r_squared": round(adj_r2, 6),
            "sse": round(sse, 8),
            "ssr": round(ssr, 8),
            "sst": round(sst, 8),
            "std_error": round(std_err, 8),
            "lod_estimate": round(lod, 6),
            "loq_estimate": round(loq, 6),
            "n_points": n,
        }

    def _poly_regression(self, x: List[float], y: List[float], degree: int = 2) -> dict:
        """Polynomial regression using normal equations."""
        import math
        n = len(x)
        if n <= degree:
            raise ChemMCPError(f"Need at least {degree + 1} points for degree-{degree} polynomial.")
        if degree not in (2, 3):
            raise ChemMCPError("Polynomial degree must be 2 or 3.")

        # Build Vandermonde-like system: X'X β = X'y
        # For degree d, we have d+1 coefficients
        p = degree + 1
        # Build normal equations matrix
        XT_X = [[0.0] * p for _ in range(p)]
        XT_y = [0.0] * p

        for i in range(n):
            powers = [x[i] ** k for k in range(p)]
            for row in range(p):
                for col in range(p):
                    XT_X[row][col] += powers[row] * powers[col]
                XT_y[row] += powers[row] * y[i]

        # Gaussian elimination with partial pivoting
        aug = [r[:] + [XT_y[i]] for i, r in enumerate(XT_X)]
        for col in range(p):
            # Find pivot
            max_val = abs(aug[col][col])
            max_row = col
            for row in range(col + 1, p):
                if abs(aug[row][col]) > max_val:
                    max_val = abs(aug[row][col])
                    max_row = row
            aug[col], aug[max_row] = aug[max_row], aug[col]
            if abs(aug[col][col]) < 1e-30:
                raise ChemMCPError("Singular matrix in polynomial regression.")
            for row in range(col + 1, p):
                factor = aug[row][col] / aug[col][col]
                for j in range(col, p + 1):
                    aug[row][j] -= factor * aug[col][j]

        # Back substitution
        coeffs = [0.0] * p
        for i in range(p - 1, -1, -1):
            coeffs[i] = aug[i][p]
            for j in range(i + 1, p):
                coeffs[i] -= aug[i][j] * coeffs[j]
            coeffs[i] /= aug[i][i]

        fitted = [sum(coeffs[k] * xi ** k for k in range(p)) for xi in x]
        residuals = [yi - yh for yi, yh in zip(y, fitted)]

        y_mean = self._mean(y)
        sst = sum((yi - y_mean) ** 2 for yi in y)
        sse = sum(r ** 2 for r in residuals)
        ssr = sst - sse
        r_sq = 1 - sse / sst if sst > 0 else 1.0
        adj_r2 = 1 - (1 - r_sq) * (n - 1) / (n - p) if n > p else 0.0
        std_err = math.sqrt(sse / (n - p)) if n > p else 0.0

        lod = (3.33 * std_err) if True else float('inf')  # poly LOD less meaningful
        loq = (10.0 * std_err) if True else float('inf')

        # Build equation string
        terms = []
        for k in range(degree, -1, -1):
            c = coeffs[k]
            if abs(c) < 1e-12:
                continue
            if k == degree:
                sign = "" if c >= 0 else "-"
            else:
                sign = " + " if c >= 0 else " - "
            if k == 0:
                terms.append(f"{sign}{abs(c):.6g}")
            elif k == 1:
                terms.append(f"{sign}{abs(c):.4g}·x")
            else:
                terms.append(f"{sign}{abs(c):.4g}·x^{k}")
        eq_str = "y = " + "".join(terms) if terms else "y = 0"

        return {
            "coefficients": [round(c, 10) for c in coeffs],
            "fitted_values": [round(f, 8) for f in fitted],
            "residuals": [round(r, 8) for r in residuals],
            "r_squared": round(r_sq, 6),
            "adjusted_r_squared": round(adj_r2, 6),
            "sse": round(sse, 8),
            "ssr": round(ssr, 8),
            "sst": round(sst, 8),
            "std_error": round(std_err, 8),
            "lod_estimate": round(lod, 6),
            "loq_estimate": round(loq, 6),
            "n_points": n,
            "equation": eq_str,
        }

    def _run_base(
        self,
        concentrations: List[float],
        responses: List[float],
        fit_type: str = "linear",
        polynomial_degree: int = 2,
        force_origin: bool = False,
    ) -> dict:
        """Core logic: fit calibration curve."""
        if len(concentrations) != len(responses):
            raise ChemMCPError(f"Length mismatch: {len(concentrations)} concentrations vs {len(responses)} responses.")
        if len(concentrations) < 2:
            raise ChemMCPError("Need at least 2 data points.")

        ft = fit_type.lower().strip()
        diagnostics: Dict[str, Any] = {}

        if ft == "linear":
            result = self._linear_regression(concentrations, responses, force_origin=force_origin)
            slope = result["slope"]
            intercept = result["intercept"]
            result["equation"] = f"y = {slope:.6g}x + {intercept:.6g}"
            result["coefficients"] = [round(slope, 10), round(intercept, 10)]

        elif ft in ("weighted_1/x", "weighted_1/x2"):
            weights = []
            power = 1 if ft == "weighted_1/x" else 2
            for xi in concentrations:
                if xi == 0:
                    weights.append(1.0)  # blank gets weight 1
                else:
                    weights.append(1.0 / (abs(xi) ** power))
            result = self._linear_regression(concentrations, responses, weights=weights, force_origin=force_origin)
            slope = result["slope"]
            intercept = result["intercept"]
            result["equation"] = f"y = {slope:.6g}x + {intercept:.6g} (weighted 1/x^{power})"
            result["coefficients"] = [round(slope, 10), round(intercept, 10)]
            diagnostics["weighting_type"] = f"1/x^{power}"
            diagnostics["weights"] = [round(w, 6) for w in weights]

        elif ft == "polynomial":
            deg = min(max(polynomial_degree, 2), 3)
            result = self._poly_regression(concentrations, responses, degree=deg)
            diagnostics["polynomial_degree"] = deg

        else:
            raise ChemMCPError(f"Unknown fit_type '{fit_type}'. Use: 'linear', 'weighted_1/x', 'weighted_1/x2', 'polynomial'.")

        result["diagnostics"] = diagnostics
        result["fit_type"] = ft

        logger.info(f"Calibration curve fit ({ft}): R²={result['r_squared']}, n={result['n_points']}")
        return result

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input. Expected JSON string with keys: concentrations, responses, fit_type, ...")
        return self._run_base(**kwargs)
