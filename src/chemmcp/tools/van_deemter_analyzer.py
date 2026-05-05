import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class VanDeemterAnalyzer(BaseTool):
    """
    Van Deemter 方程分析工具。
    拟合Van Deemter曲线（H = A + B/u + C·u），计算最佳线速度、最小塔板高度，优化色谱流速。
    """
    __version__ = "0.1.0"
    name = "VanDeemterAnalyzer"
    func_name = "analyze_van_deemter"
    description = "Analyze Van Deemter equation (H = A + B/u + C·u) to determine optimal linear velocity, minimum plate height, and flow rate optimization for chromatographic efficiency."
    implementation_description = (
        "Implements Van Deemter equation analysis: calculates A (eddy diffusion), B (longitudinal diffusion), "
        "C (mass transfer) coefficients from experimental data or user input. "
        "Computes optimal linear velocity (u_opt = √(B/C)), minimum plate height (H_min = A + 2√(B·C)), "
        "and generates Van Deemter curve data points. Supports Knox equation alternative form. "
        "Provides practical flow rate recommendations for specific column dimensions."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "Van Deemter", "Flow Optimization", "HPLC", "Plate Height", "Kinetics"]
    required_envs = []

    code_input_sig = [
        ("data_points", "None", "None", "Experimental data: [{'u_mm_s': float, 'H_um': float}, ...]. If None, use coefficient input."),
        ("A_coefficient", "float", "0.01", "Eddy diffusion term A (mm). Typical: 0.001-0.05 mm."),
        ("B_coefficient", "float", "0.0006", "Longitudinal diffusion B (mm²/s). Typical: 0.0001-0.005."),
        ("C_coefficient", "float", "0.001", "Mass transfer C (s). Typical: 0.0001-0.01 s."),
        ("particle_size_um", "float", "5.0", "Stationary phase particle size in μm (for reduced parameters)."),
        ("column_internal_diameter_mm", "float", "4.6", "Column ID for volumetric flow conversion."),
        ("diffusion_coefficient_cm2_s", "float", "1e-5", "Solute diffusion coefficient in mobile phase (cm²/s)."),
        ("velocity_range", "list", "None", "Custom velocity range [min, max] mm/s for curve generation."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'A=0.02 B=0.0008 C=0.002 dp=3 ID=4.6' or data-driven mode."),
    ]

    output_sig = [
        ("van_deemter_result", "dict", "Complete Van Deemter analysis: coefficients, optimal conditions, curve data, reduced parameters, and flow recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "A_coefficient": 0.015,
                "B_coefficient": 0.0008,
                "C_coefficient": 0.002,
                "particle_size_um": 3.0,
                "column_internal_diameter_mm": 4.6,
            },
            "text_input": {"input_params": "A=0.015 B=0.0008 C=0.002 dp=3 ID=4.6"},
            "output": {"van_deemter_result": {"u_optimal_mm_s": ..., "H_min_um": ..., "F_optimal_mL_min": ...}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _van_deemter_H(self, u: float, A: float, B: float, C: float) -> float:
        """H(u) = A + B/u + C·u."""
        if u <= 0:
            return float('inf')
        return A + B / u + C * u

    def _fit_coefficients(self, data_points: list) -> dict:
        """Least-squares fit of H = A + B/u + C·u from experimental data."""
        n = len(data_points)
        if n < 3:
            raise ChemMCPError("Need at least 3 data points for fitting.")

        # Transform to linear system: H = A*(1) + B*(1/u) + C*(u)
        H_vals = [d["H_um"] / 1000.0 for d in data_points]  # μm → mm
        u_vals = [d["u_mm_s"] for d in data_points]
        inv_u = [1.0 / u for u in u_vals]

        sum_H = sum(H_vals)
        sum_inv_u = sum(inv_u)
        sum_u = sum(u_vals)
        sum_H_inv_u = sum(h * iu for h, iu in zip(H_vals, inv_u))
        sum_H_u = sum(h * u for h, u in zip(H_vals, u_vals))
        sum_inv_u2 = sum(iu ** 2 for iu in inv_u)
        sum_u2 = sum(u ** 2 for u in u_vals)
        sum_inv_u_u = sum(iu * u for iu, u in zip(inv_u, u_vals))

        M = [[n, sum_inv_u, sum_u],
             [sum_inv_u, sum_inv_u2, sum_inv_u_u],
             [sum_u, sum_inv_u_u, sum_u2]]
        Y = [sum_H, sum_H_inv_u, sum_H_u]

        det = self._det3(M)
        if abs(det) < 1e-20:
            raise ChemMCPError("Singular matrix — check data (possibly duplicate velocities).")

        inv = self._inv3x3(M)
        A_val = sum(inv[0][i] * Y[i] for i in range(3))
        B_val = sum(inv[1][i] * Y[i] for i in range(3))
        C_val = sum(inv[2][i] * Y[i] for i in range(3))

        # R²
        H_mean = sum_H / n
        ss_tot = sum((h - H_mean) ** 2 for h in H_vals)
        ss_res = sum((h - self._van_deemter_H(u_vals[i], A_val, B_val, C_val)) ** 2
                     for i, h in enumerate(H_vals))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        return {"A_mm": round(A_val, 8), "B_mm2_s": round(B_val, 10),
                "C_s": round(C_val, 8), "R_squared": round(r2, 6), "n_data_points": n}

    @staticmethod
    def _det3(M):
        (a,b,c),(d,e,f),(g,h,i) = M
        return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)

    def _inv3x3(self, M):
        (a,b,c),(d,e,f),(g,h,i) = M
        det = self._det3(M)
        if abs(det) < 1e-30:
            raise ChemMCPError("Non-invertible matrix.")
        inv_det = 1.0 / det
        return [[(e*i-f*h)*inv_det, (c*h-b*i)*inv_det, (b*f-c*e)*inv_det],
                [(f*g-d*i)*inv_det, (a*i-c*g)*inv_det, (c*d-a*f)*inv_det],
                [(d*h-e*g)*inv_det, (b*g-a*h)*inv_det, (a*e-b*d)*inv_det]]

    def _generate_curve(self, A, B, C, u_min, u_max, n_points=50):
        u_opt = math.sqrt(B / C) if C > 0 else 10.0
        if u_min is None:
            u_min = max(0.1, u_opt * 0.05)
        if u_max is None:
            u_max = u_opt * 5.0
        return [{"u_mm_s": round(u_min + (u_max-u_min)*j/(n_points-1), 3),
                 "H_um": round(self._van_deemter_H(u_min+(u_max-u_min)*j/(n_points-1), A, B, C) * 1000, 2)}
                for j in range(n_points)]

    def _calc_reduced_parameters(self, A, B, C, dp_um, D_cm2_s):
        dp_mm = dp_um / 1000.0
        D_mm2_s = D_cm2_s * 100  # cm²/s → mm²/s
        h_min_raw = A + 2 * math.sqrt(B * C) if C > 0 else A
        u_opt = math.sqrt(B / C) if C > 0 else 1.0
        return {
            "reduced_A_h": round(A / dp_mm, 3) if dp_mm > 0 else 0,
            "reduced_B_2hnu": round(B / (D_mm2_s * dp_mm), 3) if (D_mm2_s > 0 and dp_mm > 0) else 0,
            "reduced_C_hnu": round(C * D_mm2_s / dp_mm, 3) if dp_mm > 0 else 0,
            "optimal_reduced_velocity_nu": round(u_opt * dp_mm / D_mm2_s, 2) if D_mm2_s > 0 else 0,
            "minimum_reduced_plate_height_h": round(h_min_raw / dp_mm, 2) if dp_mm > 0 else 0,
            "interpretation": self._interpret_h(h_min_raw / dp_mm if dp_mm > 0 else 99),
        }

    @staticmethod
    def _interpret_h(h):
        if h < 2: return "Excellent packing"
        elif h < 3: return "Good packing"
        elif h < 4: return "Acceptable"
        return "Suboptimal — check column quality"

    @staticmethod
    def _flow_rate_conversion(u_mm_s, col_id_mm):
        area_mm2 = math.pi * (col_id_mm / 2.0) ** 2
        return area_mm2 * u_mm_s * 60 / 1000  # mL/min

    def _get_recommendations(self, u_opt, H_min, reduced, dp, col_id):
        recs = []
        F_opt = self._flow_rate_conversion(u_opt, col_id)
        nu = reduced.get("optimal_reduced_velocity_nu", 0)
        h = reduced.get("minimum_reduced_plate_height_h", 99)

        recs.append(f"🎯 Optimal flow rate: {F_opt:.2f} mL/min (u_opt = {u_opt:.1f} mm/s)")
        recs.append({"h<2": "✓ Excellent column — world-class reduced plate height.",
                     "h<3": "✓ Good column — well within expected performance.",
                     "h<4": "⚠ Acceptable but not optimal."}.get(
                         f"h<{int(h)}" if h < 5 else "h>=5",
                         "🔶 Poor efficiency — consider replacing column."))

        if nu < 2:
            recs.append("B-term (longitudinal diffusion) dominates at optimum — normal for small particles.")
        elif nu > 10:
            recs.append("C-term (mass transfer) dominates — reduce flow rate or use smaller particles.")

        if dp <= 2:
            recs.append(f"UHPLC particles ({dp}μm) — ensure low extra-column volume (<15μL).")
        elif dp >= 5:
            recs.append("Consider sub-2μm or core-shell particles for higher efficiency.")

        return recs[:6]

    def _run_base(self, data_points=None, A_coefficient=0.01, B_coefficient=0.0006,
                  C_coefficient=0.001, particle_size_um=5.0,
                  column_internal_diameter_mm=4.6, diffusion_coefficient_cm2_s=1e-5,
                  velocity_range=None):
        """Core logic."""

        # Coefficients: from data fitting or user input
        if data_points:
            fit = self._fit_coefficients(data_points)
            A, B, C = fit["A_mm"], fit["B_mm2_s"], fit["C_s"]
            source, fit_info = "fitted_from_data", fit
        else:
            A, B, C = A_coefficient, B_coefficient, C_coefficient
            source = "user_provided"
            fit_info = {"A_mm": A, "B_mm2_s": B, "C_s": C}

        # Optimal conditions
        u_opt = math.sqrt(B / C) if C > 0 else float('inf')
        H_min = A + 2 * math.sqrt(B * C) if C > 0 else A
        F_opt = self._flow_rate_conversion(u_opt, column_internal_diameter_mm)

        # Reduced parameters (Knox-style)
        reduced = self._calc_reduced_parameters(A, B, C, particle_size_um, diffusion_coefficient_cm2_s)

        # Generate curve
        u_min = velocity_range[0] if velocity_range else None
        u_max = velocity_range[1] if velocity_range else None
        curve = self._generate_curve(A, B, C, u_min, u_max)

        # Term breakdown at u_opt
        terms = {
            "A_eddy_diffusion": round(A, 6),
            "B_longitudinal_at_uopt": round(B / u_opt, 6) if u_opt > 0 else None,
            "C_mass_transfer_at_uopt": round(C * u_opt, 6),
            "fraction_A_percent": round(A / H_min * 100, 1) if H_min > 0 else None,
            "fraction_B_percent": round((B / u_opt) / H_min * 100, 1) if (u_opt > 0 and H_min > 0) else None,
            "fraction_C_percent": round((C * u_opt) / H_min * 100, 1) if H_min > 0 else None,
        }

        # Estimate N for typical 150mm column
        N_est = int(round(150000 / (H_min * 1000))) if H_min > 0 else 0

        return {"van_deemter_result": {
            "equation": f"H = {A:.6f} + {B:.8f}/u + {C:.6f}·u",
            "coefficient_source": source,
            "coefficients_mm": fit_info,
            "optimal_conditions": {
                "optimal_linear_velocity_u_opt_mm_s": round(u_opt, 3),
                "minimum_plate_height_H_min_um": round(H_min * 1000, 1),
                "optimal_flow_rate_F_opt_mL_min": round(F_opt, 3),
                "estimated_N_for_150mm_column": N_est,
            },
            "term_breakdown_at_optimum": terms,
            "reduced_parameters": reduced,
            "van_deemter_curve": curve,
            "practical_recommendations": self._get_recommendations(
                u_opt, H_min, reduced, particle_size_um, column_internal_diameter_mm),
        }}

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        for part in input_params.strip().split():
            if "=" in part:
                k, v = part.split("=", 1)
                km = {"A": "A_coefficient", "B": "B_coefficient", "C": "C_coefficient",
                      "dp": "particle_size_um", "ID": "column_internal_diameter_mm"}.get(k, k)
                try: kwargs[km] = float(v)
                except ValueError: kwargs[km] = v
        return self._run_base(**kwargs)
