import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Typical extra-column volumes for common HPLC system components (μL)
COMPONENT_VOLUMES = {
    "injector_loop": 20,
    "injection_valve_to_column_tubing": 13,
    "column_to_detector_tubing": 13,
    "detector_flow_cell": 10,
    "guard_column": 50,
    "pre_column_filter": 3,
    "frits_inlet_outlet": 5,
}


@ChemMCPManager.register_tool
class DeadVolumeCalculator(BaseTool):
    """
    死体积和死时间计算工具。
    计算色谱系统总死体积、各组件贡献、死时间(t₀)，评估额外柱外体积对峰展宽的影响。
    """
    __version__ = "0.1.0"
    name = "DeadVolumeCalculator"
    func_name = "calculate_dead_volume"
    description = "Calculate chromatographic dead volume (V₀), dead time (t₀), and extra-column volume contributions from system components."
    implementation_description = (
        "Calculates total system dead volume by summing contributions from injector, tubing, "
        "detector cell, guard column, frits, and fittings. Computes dead time from flow rate. "
        "Assesses impact of extra-column volume on peak broadening using variance-based analysis. "
        "Provides optimization recommendations to minimize band broadening."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "Dead Volume", "Dead Time", "HPLC", "Band Broadening", "System Optimization"]
    required_envs = []

    code_input_sig = [
        ("flow_rate_mL_min", "float", "N/A", "Mobile phase flow rate in mL/min."),
        ("column_internal_diameter_mm", "float", "4.6", "Column internal diameter in mm."),
        ("column_length_mm", "float", "150.0", "Column length in mm."),
        ("tubing_id_mm", "float", "0.12", "Tubing internal diameter in mm (standard: 0.12-0.17)."),
        ("tubing_length_cm", "dict", "{}", "Tubing lengths: {'inlet_cm': 30, 'outlet_cm': 30}."),
        ("detector_cell_volume_uL", "float", "10.0", "Detector flow cell volume in μL."),
        ("injection_volume_uL", "float", "10.0", "Injection volume in μL."),
        ("guard_column_volume_uL", "None", "None", "Guard column volume in μL (if used)."),
        ("use_guard_column", "bool", "False", "Whether a guard column is installed."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'F=1.0 ID=4.6 L=150 tubing=0.12 detector=10'."),
    ]

    output_sig = [
        ("dead_volume_analysis", "dict", "Complete dead volume analysis including component breakdown, dead time, variance contribution, and optimization suggestions."),
    ]

    examples = [
        {
            "code_input": {
                "flow_rate_mL_min": 1.0,
                "column_internal_diameter_mm": 4.6,
                "column_length_mm": 150,
            },
            "text_input": {"input_params": "F=1.0 ID=4.6 L=150"},
            "output": {
                "dead_volume_analysis": {"total_dead_volume_uL": ..., "dead_time_min": ...}
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_tubing_volume(self, id_mm: float, length_cm: float) -> float:
        """Cylindrical tube volume in μL."""
        radius_mm = id_mm / 2.0
        volume_mm3 = math.pi * (radius_mm ** 2) * length_cm * 10  # cm → mm
        return volume_mm3 / 1000  # mm³ → μL

    def _calc_column_dead_volume(self, id_mm: float, length_mm: float) -> float:
        """Geometric column void volume (≈ 60% of total geometric volume for packed bed)."""
        radius_mm = id_mm / 2.0
        total_geom_vol_mm3 = math.pi * (radius_mm ** 2) * length_mm
        total_geom_vol_uL = total_geom_vol_mm3 / 1000
        # Packed column porosity ≈ 0.65-0.70; interstitial fraction ≈ 0.55-0.65
        void_fraction = 0.62  # typical for fully porous particles
        return total_geom_vol_uL * void_fraction

    def _calc_extra_column_variance(self, volumes: dict, peak_width_sec: float = 5.0) -> dict:
        """
        Estimate variance contribution of each extra-column component.
        σ²_total = Σσ²_i where σ_i ≈ V_i / √12 for rectangular profile.
        """
        variances = {}
        for name, vol_uL in volumes.items():
            sigma_sq = (vol_uL ** 2) / 12.0  # Variance of uniform distribution
            variances[name] = {"volume_uL": round(vol_uL, 1), "variance": round(sigma_sq, 2)}

        total_var = sum(v["variance"] for v in variances.values())
        total_sigma = math.sqrt(total_var)

        # Compare to typical peak variance
        # Peak variance ≈ (W/4)² for Gaussian peak
        # W in μL = flow_rate (μL/min) × width (min)
        peak_var_approx = (peak_width_sec / 60.0 * 1000 / 4) ** 2 if peak_width_sec > 0 else 1  # rough estimate
        extra_pct = (total_var / peak_var_approx * 100) if peak_var_approx > 0 else 0

        return {
            "component_variances": variances,
            "total_extra_column_variance": round(total_var, 2),
            "total_extra_column_sigma_uL": round(total_sigma, 2),
            "percent_of_peak_variance": round(extra_pct, 1),
            "acceptable": extra_pct < 10,  # Rule of thumb: extra-column < 10% of peak variance
        }

    def _run_base(self, flow_rate_mL_min: float,
                  column_internal_diameter_mm: float = 4.6,
                  column_length_mm: float = 150.0,
                  tubing_id_mm: float = 0.12,
                  tubing_length_cm: Optional[Dict] = None,
                  detector_cell_volume_uL: float = 10.0,
                  injection_volume_uL: float = 10.0,
                  guard_column_volume_uL: Optional[float] = None,
                  use_guard_column: bool = False) -> dict:
        """Core logic."""

        if tubing_length_cm is None:
            tubing_length_cm = {"inlet_cm": 30, "outlet_cm": 30}

        # Component-by-component breakdown
        inlet_tubing_vol = self._calc_tubing_volume(
            tubing_id_mm, tubing_length_cm.get("inlet_cm", 30))
        outlet_tubing_vol = self._calc_tubing_volume(
            tubing_id_mm, tubing_length_cm.get("outlet_cm", 30))

        components = {
            "Injection volume (estimated contribution)": injection_volume_uL / 3.0,  # ~1/3 contributes to band spread
            "Inlet tubing": round(inlet_tubing_vol, 1),
            "Outlet tubing": round(outlet_tubing_vol, 1),
            "Detector flow cell": detector_cell_volume_uL,
            "Frits & fittings": COMPONENT_VOLUMES["frits_inlet_outlet"],
        }

        if use_guard_column:
            gv = guard_column_volume_uL or COMPONENT_VOLUMES["guard_column"]
            components["Guard column"] = gv

        total_extra = sum(components.values())
        column_v0 = self._calc_column_dead_volume(column_internal_diameter_mm, column_length_mm)
        total_system_v0 = total_extra + column_v0

        # Dead times
        t0_system_min = total_system_v0 / flow_rate_mL_min
        t0_column_min = column_v0 / flow_rate_mL_min
        t0_extra_min = total_extra / flow_rate_mL_min

        # Variance analysis
        variance = self._calc_extra_column_variance(components)

        # Assessment
        extra_ratio = total_extra / column_v0 * 100 if column_v0 > 0 else 999
        if extra_ratio < 5:
            assessment = "✅ Excellent — minimal extra-column effects"
        elif extra_ratio < 15:
            assessment = "✅ Acceptable — standard HPLC configuration"
        elif extra_ratio < 30:
            assessment = "⚠️ Moderate — may affect early-eluting peaks"
        else:
            assessment = "🔶 High — significant band broadening expected"

        result = {
            "dead_volume_analysis": {
                "system_parameters": {
                    "flow_rate_mL_min": flow_rate_mL_min,
                    "column_ID_mm": column_internal_diameter_mm,
                    "column_length_mm": column_length_mm,
                    "tubing_ID_mm": tubing_id_mm,
                },
                "volume_breakdown": {
                    "column_void_volume_V0_uL": round(column_v0, 1),
                    "extra_column_volume_uL": round(total_extra, 1),
                    "total_system_dead_volume_uL": round(total_system_v0, 1),
                    "extra_column_percentage_of_column_V0": round(extra_ratio, 1),
                },
                "component_details_uL": {k: v for k, v in components.items()},
                "dead_time": {
                    "t0_column_min": round(t0_column_min, 3),
                    "t0_extra_column_min": round(t0_extra_min, 3),
                    "t0_total_system_min": round(t0_system_min, 3),
                },
                "variance_analysis": variance,
                "assessment": assessment,
                "recommendations": self._get_recommendations(
                    extra_ratio, column_internal_diameter_mm, tubing_id_mm,
                    detector_cell_volume_uL, injection_volume_uL,
                ),
            }
        }
        return result

    def _get_recommendations(self, extra_pct: float, col_id: float,
                              tub_id: float, det_vol: float, inj_vol: float) -> List[str]:
        recs = []

        if extra_pct > 20:
            recs.append("Reduce tubing ID — use 0.12mm or smaller ID tubing throughout.")
            recs.append("Shorten all connecting tubing to minimum necessary length.")
            if det_vol > 8:
                recs.append(f"Consider low-volume detector cell (current: {det_vol}μL; ≤5μL available).")
            if inj_vol > 20:
                recs.append("Reduce injection volume to minimize extra-column contribution.")
        elif extra_pct > 10:
            recs.append("Consider optimizing tubing lengths for UHPLC-grade performance.")

        if col_id <= 2.1:
            recs.append("Narrow-bore column detected — ensure all connections use 0.09-0.12mm ID tubing.")
            recs.append("Use micro-flow-rate compatible detector cell (≤2μL).")

        if tub_id > 0.17:
            recs.append(f"Tubing ID ({tub_id}mm) is large — consider switching to 0.12mm ID.")

        if not recs:
            recs.append("✓ Extra-column volume is well-controlled for this system configuration.")

        return recs[:6]

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                key_map = {
                    "F": "flow_rate_mL_min", "ID": "column_internal_diameter_mm",
                    "L": "column_length_mm", "tubing": "tubing_id_mm",
                    "detector": "detector_cell_volume_uL", "inj": "injection_volume_uL",
                }
                mapped = key_map.get(k, k)
                try:
                    kwargs[mapped] = float(v)
                except ValueError:
                    kwargs[mapped] = v
        return self._run_base(**kwargs)
