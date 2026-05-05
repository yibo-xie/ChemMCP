import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PlateNumberCalculator(BaseTool):
    """
    理论塔板数和柱效计算工具。
    计算色谱柱的理论塔板数(N)、塔板高度(H)、峰容量等柱效参数，评估和优化色谱系统性能。
    """
    __version__ = "0.1.0"
    name = "PlateNumberCalculator"
    func_name = "calculate_plate_number"
    description = "Calculate chromatographic column efficiency: theoretical plates (N), plate height (H), peak capacity, and related performance metrics."
    implementation_description = (
        "Implements standard chromatographic efficiency calculations including "
        "theoretical plate number (N) from peak width at half-height or baseline, "
        "van Deemter equation parameters (A, B, C terms), reduced plate height (h), "
        "and peak capacity for gradient separations."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "Column Efficiency", "Plate Number", "Van Deemter", "HPLC"]
    required_envs = []

    code_input_sig = [
        ("retention_time_min", "float", "N/A", "Peak retention time in minutes."),
        ("peak_width_baseline_min", "None", "None", "Peak width at baseline (between tangents) in minutes."),
        ("peak_width_half_height_min", "None", "None", "Peak width at half height in minutes."),
        ("dead_time_min", "None", "None", "Column dead time t₀ in minutes (for k calculation)."),
        ("column_length_mm", "None", "None", "Column length in mm (for H calculation)."),
        ("particle_size_um", "None", "None", "Particle size in μm (for reduced parameters)."),
        ("linear_velocity_mm_min", "None", "None", "Linear velocity in mm/min (for van Deemter analysis)."),
        ("van_deemter_params", "dict", "{}", "Known van Deemter coefficients: {'A': 1.0, 'B': 2.0, 'C': 0.05}."),
        ("gradient_time_min", "None", "None", "Gradient time in min (for peak capacity)."),
        ("k_range", "None", "None", "(k_first, k_last) for gradient peak capacity."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters as 'tR=5.0 Wb=0.3 L=100' or 'tR=5 Wh=0.15 t0=0.5 dp=1.7'."),
    ]

    output_sig = [
        ("efficiency", "dict", "Complete column efficiency analysis including N, H, h, k, resolution potential, and optimization suggestions."),
    ]

    examples = [
        {
            "code_input": {
                "retention_time_min": 5.0,
                "peak_width_half_height_min": 0.15,
                "dead_time_min": 0.8,
                "column_length_mm": 100,
                "particle_size_um": 1.7,
            },
            "text_input": {"input_params": "tR=5.0 Wh=0.15 t0=0.8 L=100 dp=1.7"},
            "output": {
                "efficiency": {"N": 17778, "H_um": 5.6, "k": 5.25, "reduced_h": 3.3}
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_N_wh(self, tR: float, Wh: float) -> float:
        """N from width at half-height: N = 5.54 * (tR/Wh)²"""
        return 5.54 * (tR / Wh) ** 2

    def _calc_N_wb(self, tR: float, Wb: float) -> float:
        """N from baseline width: N = 16 * (tR/Wb)²"""
        return 16.0 * (tR / Wb) ** 2

    def _calc_H(self, N: float, L_mm: float) -> float:
        """Plate height in μm."""
        if L_mm is None or L_mm <= 0:
            return None
        return (L_mm * 1000) / N  # mm → μm

    def _calc_reduced_h(self, H_um: float, dp_um: float) -> Optional[float]:
        """Reduced plate height h = H/dp."""
        if H_um is None or dp_um is None:
            return None
        return H_um / dp_um

    def _van_deemter(self, u: float, A: float = 1.0, B: float = 2.0, C: float = 0.05) -> dict:
        """H = A + B/u + C*u (u in mm/s)."""
        # Convert u from mm/min to mm/s
        u_s = u / 60.0
        H = A + B / max(u_s, 0.001) + C * u_s * 60  # Keep consistent units
        # Actually let's work in mm/min throughout
        H_val = A + B * 60.0 / max(u, 0.001) + C * u / 60.0
        u_opt = math.sqrt(B * 3600 / C) if C > 0 else u  # optimal velocity in mm/min
        H_min = A + 2 * math.sqrt(A * C * 3600 + B * C / 10000) if C > 0 else A

        return {
            "H_at_u": round(H_val, 3),
            "optimal_velocity_mm_per_min": round(u_opt, 1),
            "H_minimum": round(H_min, 3),
            "A_eddy_diffusion": A,
            "B_longitudinal": B,
            "C_mass_transfer": C,
            "equation": f"H = {A} + {B}/u + {C}·u",
        }

    def _peak_capacity_gradient(self, tg: float, t0: float, k_last: float = 20.0) -> float:
        """Gradient peak capacity nc = 1 + (√N/4)·ln(1+k_last) or simplified."""
        # Simplified: nc ≈ 1 + tg/t0 · (average slope factor)
        # More accurate: nc = (tg/t0) / (1 + average_k) approximation
        # Standard formula: nc = 1 + (√N/4.4) · ln((1+k_last)/(1+k_first))
        # Using practical formula: nc ≈ (tg / t0) * (gradient steepness correction)
        nc = 1 + (tg / t0) * 0.87 * math.log10(1 + k_last)
        return nc

    def _interpret_efficiency(self, N: float, h_red: Optional[float], k: Optional[float]) -> List[str]:
        notes = []
        if N < 2000:
            notes.append("⚠ Very low efficiency — check column condition, extra-column volume, or injection issues")
        elif N < 5000:
            notes.append("Below-average efficiency — consider replacing column or optimizing system")
        elif N < 15000:
            notes.append("Moderate efficiency — acceptable for routine analysis")
        elif N < 30000:
            notes.append("Good efficiency — well-performing system")
        else:
            notes.append("Excellent efficiency — high-performance UHPLC-grade system")

        if h_red is not None:
            if h_red < 2.0:
                notes.append("✓ Excellent reduced plate height (h < 2) — world-class packing")
            elif h_red < 3.0:
                notes.append("✓ Good reduced plate height (h < 3) — well-packed column")
            elif h_red < 4.0:
                notes.append("Acceptable reduced plate height (h < 4) — typical for conventional columns")
            else:
                notes.append("⚠ High reduced plate height (h > 4) — check for extra-band broadening")

        if k is not None:
            if k < 1:
                notes.append("⚠ Low retention factor (k < 1) — poor separation from void")
            elif k < 2:
                notes.append("Retention factor near lower limit of ideal range (2-10)")
            elif k <= 10:
                notes.append("✓ Retention factor in optimal range (2-10)")
            elif k <= 20:
                notes.append("Retention acceptable but on higher side (k > 10)")
            else:
                notes.append("Very high retention — long analysis time; consider stronger eluent")

        return notes

    def _run_base(self, retention_time_min: float,
                  peak_width_baseline_min: Optional[float] = None,
                  peak_width_half_height_min: Optional[float] = None,
                  dead_time_min: Optional[float] = None,
                  column_length_mm: Optional[float] = None,
                  particle_size_um: Optional[float] = None,
                  linear_velocity_mm_min: Optional[float] = None,
                  van_deemter_params: Optional[Dict] = None,
                  gradient_time_min: Optional[float] = None,
                  k_range: Optional[tuple] = None) -> dict:
        """Core calculation logic."""

        if peak_width_half_height_min is None and peak_width_baseline_min is None:
            raise ChemMCPError("Either peak_width_half_height_min or peak_width_baseline_min must be provided.")

        # Calculate N
        if peak_width_half_height_min is not None:
            N = self._calc_N_wh(retention_time_min, peak_width_half_height_min)
            method_N = "half-height (Wh)"
        else:
            N = self._calc_N_wb(retention_time_min, peak_width_baseline_min)
            method_N = "baseline (Wb)"

        # Derived quantities
        H_um = self._calc_H(N, column_length_mm)
        h_red = self._calc_reduced_h(H_um, particle_size_um)
        k_factor = None
        if dead_time_min:
            k_factor = (retention_time_min - dead_time_min) / dead_time_min

        # Asymmetry estimate (if both widths available)
        asymmetry = None
        if peak_width_baseline_min and peak_width_half_height_min:
            asymmetry = peak_width_baseline_min / (2 * peak_width_half_height_min)

        # Van Deemter analysis
        vd_result = None
        if linear_velocity_mm_min:
            params = van_deemter_params or {}
            vd_result = self._van_deemter(
                linear_velocity_mm_min,
                params.get("A", 1.0),
                params.get("B", 2.0),
                params.get("C", 0.05),
            )

        # Peak capacity (gradient)
        pc_result = None
        if gradient_time_min and dead_time_min:
            k_last = k_range[1] if k_range and len(k_range) > 1 else 20.0
            pc = self._peak_capacity_gradient(gradient_time_min, dead_time_min, k_last)
            pc_result = {
                "peak_capacity": round(pc, 1),
                "gradient_time_min": gradient_time_min,
                "dead_time_min": dead_time_min,
                "note": f"Can resolve approximately {int(pc)} peaks in this gradient window",
            }

        interpretation = self._interpret_efficiency(N, h_red, k_factor)

        result = {
            "efficiency": {
                "input_peak": {
                    "retention_time_min": retention_time_min,
                    "width_half_height_min": peak_width_half_height_min,
                    "width_baseline_min": peak_width_baseline_min,
                },
                "primary_metrics": {
                    "theoretical_plates_N": int(round(N)),
                    "method_used": method_N,
                    "plate_height_H_um": round(H_um, 1) if H_um else None,
                    "reduced_plate_height_h": round(h_red, 2) if h_red else None,
                    "retention_factor_k": round(k_factor, 3) if k_factor else None,
                    "asymmetry_factor_As": round(asymmetry, 2) if asymmetry else None,
                },
                "column_parameters": {
                    "length_mm": column_length_mm,
                    "particle_size_um": particle_size_um,
                    "dead_time_t0_min": dead_time_min,
                },
                "van_deemter_analysis": vd_result,
                "gradient_peak_capacity": pc_result,
                "quality_assessment": interpretation,
                "optimization_suggestions": self._get_optimization_suggestions(N, h_red, k_factor, particle_size_um),
            }
        }
        return result

    def _get_optimization_suggestions(self, N: float, h_red: Optional[float],
                                       k: Optional[float], dp: Optional[float]) -> List[str]:
        suggestions = []
        if N < 5000:
            suggestions.append("Replace column if >1000 injections or backpressure increased significantly")
            suggestions.append("Check for extra-column band broadening (tubing volume, detector cell)")
            suggestions.append("Reduce injection volume (<1% of peak volume)")
        if h_red and h_red > 4:
            suggestions.append("Use smaller ID tubing (0.12mm ID) to reduce extra-column effects")
            if dp and dp > 3:
                suggestions.append("Consider switching to sub-2μm or core-shell particles for better efficiency")
        if k is not None and k < 1:
            suggestions.append("Weaken initial mobile phase to increase k into 2-10 range")
        if k is not None and k > 20:
            suggestions.append("Strengthen mobile phase or use steeper gradient to reduce analysis time")
        return suggestions[:5]

    def _run_text(self, input_params: str) -> dict:
        """Parse key=value format input."""
        kwargs = {}
        for part in input_params.strip().split():
            if "=" in part:
                key, val = part.split("=", 1)
                try:
                    if key == "k_range":
                        vals = [float(x) for x in val.strip("()").split(",")]
                        kwargs[key] = tuple(vals)
                    elif key == "van_deemter_params":
                        pass  # Skip complex param in text mode
                    else:
                        kwargs[key] = float(val)
                except ValueError:
                    kwargs[key] = val
        if "retention_time_min" not in kwargs and "tR" in kwargs:
            kwargs["retention_time_min"] = kwargs.pop("tR")
        if "peak_width_half_height_min" not in kwargs and "Wh" in kwargs:
            kwargs["peak_width_half_height_min"] = kwargs.pop("Wh")
        if "peak_width_baseline_min" not in kwargs and "Wb" in kwargs:
            kwargs["peak_width_baseline_min"] = kwargs.pop("Wb")
        if "dead_time_min" not in kwargs and "t0" in kwargs:
            kwargs["dead_time_min"] = kwargs.pop("t0")
        if "column_length_mm" not in kwargs and "L" in kwargs:
            kwargs["column_length_mm"] = kwargs.pop("L")
        if "particle_size_um" not in kwargs and "dp" in kwargs:
            kwargs["particle_size_um"] = kwargs.pop("dp")

        return self._run_base(**kwargs)
