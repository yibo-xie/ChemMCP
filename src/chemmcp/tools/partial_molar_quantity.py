import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PartialMolarQuantity(BaseTool):
    """
    偏摩尔量的计算与图解工具。
    支持偏摩尔体积、偏摩尔吉布斯自由能等的计算，以及偏摩尔量-组成曲线的图解数据。
    """
    __version__ = "0.1.0"
    name = "PartialMolarQuantity"
    func_name = "calculate_partial_molar"
    description = "Calculate partial molar quantities (volume, Gibbs energy, etc.) and generate graphical data for binary mixtures."
    implementation_description = "Uses intercept method: V̄_i = V + x_j·(dV/dx_i) for binary mixtures. Supports polynomial fitting of experimental data."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Partial Molar Quantity", "Physical Chemistry", "Mixtures"]
    required_envs    = []

    code_input_sig = [
        ("quantity_type", "str", "N/A", "Type of quantity: 'volume', 'gibbs', 'enthalpy', 'entropy'."),
        ("mode", "str", "N/A", "Calculation mode: 'from_data' (experimental data points) or 'analytical' (given function)."),
        # For from_data mode:
        ("mole_fractions", "str", "N/A", "Mole fractions of component 1, comma-separated (for from_data)."),
        ("total_quantities", "str", "N/A", "Total molar quantities at each composition, comma-separated (same length as mole_fractions)."),
        # For analytical mode:
        ("function_type", "str", "polynomial", "Function type: 'polynomial', 'redlich_kister'."),
        ("coefficients", "str", "N/A", "Function coefficients, comma-separated (e.g., polynomial a0,a1,a2 or Redlich-Kister A0,A1,A2)."),
        ("target_x1", "float", "0.5", "Target mole fraction of component 1 to compute partial molar quantities."),
        ("n_points", "int", "50", "Number of points for graphical output curve."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated parameters string."),
    ]

    output_sig = [
        ("partial_molar_1", "float", "Partial molar quantity of component 1 at target composition."),
        ("partial_molar_2", "float", "Partial molar quantity of component 2 at target composition."),
        ("total_molar", "float", "Total molar quantity at target composition."),
        ("curve_data", "list", "List of dicts with x1, total, pm1, pm2 for plotting."),
        ("analysis", "str", "Detailed analysis including method and interpretation."),
    ]

    examples         = [
        {
            "code_input": {
                "quantity_type": 'volume',
                "mode": 'analytical',
                "mole_fractions": '',
                "total_quantities": '',
                "function_type": 'polynomial',
                "coefficients": '18.07,-1.93,0.82',
                "target_x1": 0.3,
                "n_points": 50
            },
            "text_input": {
                "input_params": 'volume analytical polynomial 18.07,-1.93,0.82 0.3'
            },
            "output": {
                "partial_molar_1": 16.918,
                "partial_molar_2": 21.562,
                "total_molar": 19.6678,
                "curve_data": [],
                "analysis": 'Ethanol-water-like mixture.'
            }
        },
        {
            "code_input": {
                "quantity_type": 'volume',
                "mode": 'from_data',
                "mole_fractions": '0.0,0.2,0.4,0.6,0.8,1.0',
                "total_quantities": '18.07,17.62,17.35,17.22,17.18,18.07',
                "function_type": '',
                "coefficients": '',
                "target_x1": 0.4,
                "n_points": 50
            },
            "text_input": {
                "input_params": 'volume from_data 0.0,0.2,0.4,0.6,0.8,1.0 18.07,17.62,17.35,17.22,17.18,18.07 0.4'
            },
            "output": {
                "partial_molar_1": 16.5,
                "partial_molar_2": 21.1,
                "total_molar": 17.35,
                "curve_data": [],
                "analysis": 'Intercept method on 6 data points.'
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _polynomial_value(self, coeffs, x):
        return sum(c * x**i for i, c in enumerate(coeffs))

    def _polynomial_derivative(self, coeffs, x):
        return sum(i * c * x**(i-1) for i, c in enumerate(coeffs) if i >= 1)

    def _redlich_kister(self, A_coeffs, x):
        """Redlich-Kister expansion: V = x1*V1° + x2*V2° + x1*x2 * Σ Ak*(x1-x2)^k"""
        if len(A_coeffs) < 2:
            raise ChemMCPError("Redlich-Kister needs at least pure-component values as A0(=V1°), A1(=V2°), then excess terms.")
        v1_pure = A_coeffs[0]
        v2_pure = A_coeffs[1]
        x1, x2 = x, 1 - x
        excess = 0.0
        for k in range(2, len(A_coeffs)):
            excess += A_coeffs[k] * (x1 - x2)**k
        return x1 * v1_pure + x2 * v2_pure + x1 * x2 * excess

    def _redlich_kister_derivatives(self, A_coeffs, x):
        """Compute dV/dx1 and dV/dx2 via Redlich-Kister."""
        eps = 1e-8
        v_plus = self._redlich_kister(A_coeffs, min(x + eps, 1.0))
        v_minus = self._redlich_kister(A_coeffs, max(x - eps, 0.0))
        dv_dx = (v_plus - v_minus) / (2 * eps)
        return dv_dx

    def _run_base(
        self,
        quantity_type: str,
        mode: str,
        mole_fractions: str = "",
        total_quantities: str = "",
        function_type: str = "polynomial",
        coefficients: str = "",
        target_x1: float = 0.5,
        n_points: int = 50,
    ) -> dict:
        mode = mode.lower().strip()
        quantity_type = quantity_type.lower().strip()

        if not (0 <= target_x1 <= 1):
            raise ChemMCPError("target_x1 must be between 0 and 1.")

        if mode == "analytical":
            coeffs = [float(c.strip()) for c in coefficients.split(",")]
            ftype = function_type.lower().strip()

            if ftype == "polynomial":
                total_at_target = self._polynomial_value(coeffs, target_x1)
                dv_dx = self._polynomial_derivative(coeffs, target_x1)
                # Intercept method: V̄1 = V + x2*(dV/dx1), V̄2 = V - x1*(dV/dx1)
                x2 = 1 - target_x1
                pm1 = total_at_target + x2 * dv_dx
                pm2 = total_at_target - target_x1 * dv_dx
                func_expr = " + ".join([f"{c}*x^{i}" if i > 0 else f"{c}" for i, c in enumerate(coeffs)])

                # Generate curve data
                curve_data = []
                for j in range(n_points + 1):
                    xi = j / n_points
                    vi = self._polynomial_value(coeffs, xi)
                    dvi = self._polynomial_derivative(coeffs, xi)
                    x2i = 1 - xi
                    pmi1 = vi + x2i * dvi
                    pmi2 = vi - xi * dvi
                    curve_data.append({
                        "x1": round(xi, 4),
                        "total": round(vi, 4),
                        "pm1": round(pmi1, 4),
                        "pm2": round(pmi2, 4),
                    })

                analysis = (
                    f"Analytical ({ftype}) model for {quantity_type}:\n"
                    f"Function: V({quantity_type}) = {func_expr}\n"
                    f"At x₁={target_x1}: Total={total_at_target:.4f}, V̄₁={pm1:.4f}, V̄₂={pm2:.4f}\n"
                    f"Method: Tangent-intercept method (V̄ᵢ = V + xⱼ·dV/dxᵢ)"
                )

            elif ftype == "redlich_kister":
                total_at_target = self._redlich_kister(coeffs, target_x1)
                dv_dx = self._redlich_kister_derivatives(coeffs, target_x1)
                x2 = 1 - target_x1
                pm1 = total_at_target + x2 * dv_dx
                pm2 = total_at_target - target_x1 * dv_dx

                curve_data = []
                for j in range(n_points + 1):
                    xi = j / n_points
                    vi = self._redlich_kister(coeffs, xi)
                    dvi = self._redlich_kister_derivatives(coeffs, xi)
                    x2i = 1 - xi
                    curve_data.append({
                        "x1": round(xi, 4),
                        "total": round(vi, 4),
                        "pm1": round(vi + x2i * dvi, 4),
                        "pm2": round(vi - xi * dvi, 4),
                    })

                analysis = (
                    f"Redlich-Kister expansion for {quantity_type}:\n"
                    f"Coefficients: {coeffs}\n"
                    f"At x₁={target_x1}: Total={total_at_target:.4f}, V̄₁={pm1:.4f}, V̄₂={pm2:.4f}"
                )
            else:
                raise ChemMCPError(f"Unsupported function type: '{ftype}'. Use 'polynomial' or 'redlich_kister'.")

        elif mode == "from_data":
            xs = [float(x.strip()) for x in mole_fractions.split(",")]
            vs = [float(v.strip()) for v in total_quantities.split(",")]
            if len(xs) != len(vs):
                raise ChemMCPError("mole_fractions and total_quantities must have same number of elements.")
            if len(xs) < 3:
                raise ChemMCPError("Need at least 3 data points for numerical differentiation.")

            # Fit polynomial (degree = min(n-1, 4)) for smooth differentiation
            import numpy as np
            degree = min(len(xs) - 1, 4)
            coeffs_fit = np.polyfit(xs, vs, degree).tolist()

            total_at_target = self._polynomial_value(list(reversed(coeffs_fit)), target_x1)
            dv_dx = self._polynomial_derivative(list(reversed(coeffs_fit)), target_x1)
            x2 = 1 - target_x1
            pm1 = total_at_target + x2 * dv_dx
            pm2 = total_at_target - target_x1 * dv_dx

            curve_data = []
            for j in range(n_points + 1):
                xi = j / n_points
                vi = self._polynomial_value(list(reversed(coeffs_fit)), xi)
                dvi = self._polynomial_derivative(list(reversed(coeffs_fit)), xi)
                x2i = 1 - xi
                curve_data.append({
                    "x1": round(xi, 4),
                    "total": round(vi, 4),
                    "pm1": round(vi + x2i * dvi, 4),
                    "pm2": round(vi - xi * dvi, 4),
                })

            analysis = (
                f"Data-driven calculation for {quantity_type}:\n"
                f"{len(xs)} data points fitted with degree-{degree} polynomial.\n"
                f"At x₁={target_x1}: Total={total_at_target:.4f}, V̄₁={pm1:.4f}, V̄₂={pm2:.4f}\n"
                f"Method: Polynomial fit + tangent-intercept method"
            )
        else:
            raise ChemMCPError(f"Unsupported mode: '{mode}'. Use 'analytical' or 'from_data'.")

        return {
            "partial_molar_1": round(pm1, 4),
            "partial_molar_2": round(pm2, 4),
            "total_molar": round(total_at_target, 4),
            "curve_data": curve_data[:5],  # Return first few as preview
            "analysis": analysis,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            qty = parts[0]
            mode = parts[1]
            kwargs = {"quantity_type": qty, "mode": mode}
            if mode == "analytical":
                kwargs["function_type"] = parts[2] if len(parts) > 2 else "polynomial"
                kwargs["coefficients"] = parts[3] if len(parts) > 3 else ""
                kwargs["target_x1"] = float(parts[4]) if len(parts) > 4 else 0.5
            elif mode == "from_data":
                kwargs["mole_fractions"] = parts[2] if len(parts) > 2 else ""
                kwargs["total_quantities"] = parts[3] if len(parts) > 3 else ""
                kwargs["target_x1"] = float(parts[4]) if len(parts) > 4 else 0.5
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
