import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ─── Carrier Gas Properties ───
_GAS_PROPERTIES = {
    "helium": {
        "name": "Helium (He)",
        "molecular_weight_g_mol": 4.00,
        "viscosity_25c_uPa_s": 196.0,
        "diffusivity_cm2_s": 0.3-0.6,  # Approximate in carrier gas at 100°C
        "optimal_linear_velocity_cm_s": 20-40,  # Van Deemter optimum range
        "optimal_avg_velocity_cm_s": 25,
        "golay_optimum_ratio": 1.5,  # u_opt(golay) / u_opt(van deemter)
        "density_g_L_25c": 0.164,
        "thermal_conductivity_w_mK": 152.0,
        "ms_compatibility": "Excellent (inert)",
        "ionization_eV": 24.6,
        "safety": "Non-flammable; asphyxiation hazard in confined spaces",
        "cost": "High (limited global supply); consider recycling",
        "availability": "Supply-constrained since ~2022",
        "detector_performance": {
            "TID": "Good sensitivity",
            "FID": "Standard choice",
            "ECD": "Excellent baseline stability",
            "MS": "Best overall for MS detection",
            "NPD": "Good",
            "SCD/NCD": "Good",
        },
    },
    "hydrogen": {
        "name": "Hydrogen (H₂)",
        "molecular_weight_g_mol": 2.02,
        "viscosity_25c_uPa_s": 89.2,
        "diffusivity_cm2_s": 0.4-0.8,
        "optimal_linear_velocity_cm_s": 40-80,
        "optimal_avg_velocity_cm_s": 55,
        "golay_optimum_ratio": 1.3,
        "density_g_L_25c": 0.082,
        "thermal_conductivity_w_mK": 181.0,
        "ms_compatibility": "Good (reductive); may affect active compounds",
        "ionization_eV": 15.4,
        "safety": "⚠ Flammable (LEL 4% in air); requires safety measures",
        "cost": "Very low (generator or cylinder)",
        "availability": "Unlimited (can be generated on-site)",
        "detector_performance": {
            "TID": "Highest sensitivity (best κ difference vs organics)",
            "FID": "Excellent; slightly different response factors than He",
            "ECD": "Acceptable but baseline noise can be higher",
            "MS": "Good; watch for reduction of analytes",
            "NPD": "Good",
            "SCD/NCD": "Good",
        },
    },
    "nitrogen": {
        "name": "Nitrogen (N₂)",
        "molecular_weight_g_mol": 28.01,
        "viscosity_25c_uPa_s": 178.0,
        "diffusivity_cm2_s": 0.08-0.12,
        "optimal_linear_velocity_cm_s": 8-14,
        "optimal_avg_velocity_cm_s": 11,
        "golay_optimum_ratio": 2.0,
        "density_g_L_25c": 1.145,
        "thermal_conductivity_w_mK": 26.0,
        "ms_compatibility": "Good (inert)",
        "ionization_eV": 15.6,
        "safety": "Non-flammable; asphyxiation hazard",
        "cost": "Low",
        "availability": "Readily available",
        "detector_performance": {
            "TID": "Poor (κ too close to many organics)",
            "FID": "Acceptable but slow analysis",
            "ECD": "Excellent for ECD applications",
            "MS": "Good",
            "NPD": "Good",
            "SCD/NCD": "Good",
        },
    },
}

# Van Deemter equation parameters (approximate) for capillary GC columns
# H = A + B/u + C·u
# A = eddy diffusion (very small for capillary), B = longitudinal diffusion, C = mass transfer
_VAN_DEEMTER_GC = {
    "helium": {"A_um": 0.02, "B_um2_s": 200, "C_s_um": 0.00005},
    "hydrogen": {"A_um": 0.02, "B_um2_s": 80, "C_s_um": 0.00003},
    "nitrogen": {"A_um": 0.02, "B_um2_s": 600, "C_s_um": 0.00010},
}


@ChemMCPManager.register_tool
class GcCarrierGasSelector(BaseTool):
    """
    载气选择与流速优化工具。
    根据GC分析需求（检测器类型、分析速度、分离效率、安全性等）推荐最佳载气，
    并计算最优流速和线速度。
    """
    __version__ = "0.1.0"
    name = "GcCarrierGasSelector"
    func_name = "select_carrier_gas"
    description = "Select optimal GC carrier gas (He, H₂, N₂) and calculate optimal flow rate/linear velocity based on detector type, column dimensions, and analysis goals."
    implementation_description = (
        "Uses van Deemter and Golay equation parameters for each carrier gas to determine "
        "optimal linear velocity and flow rate. Considers detector compatibility, "
        "analysis speed requirements, safety constraints, and cost/availability factors."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["GC", "Carrier Gas", "Van Deemter", "Gas Chromatography", "Method Development"]
    required_envs = []

    code_input_sig = [
        ("detector_type", "str", "FID", "Detector: 'FID', 'TID', 'MS', 'ECD', 'NPD', 'SCD', 'NCD'."),
        ("column_length_m", "float", "30", "Column length in meters."),
        ("column_inner_diameter_mm", "float", "0.25", "Column inner diameter in mm."),
        ("film_thickness_um", "float", "0.25", "Stationary phase film thickness in μm."),
        ("analysis_goal", "str", "balanced", "Goal: 'speed' (fastest analysis), 'efficiency' (best resolution), 'balanced' (default)."),
        ("oven_temperature_c", "float", "150", "Average oven temperature for viscosity correction."),
        ("max_pressure_kpa", "None", "None", "System pressure limit (optional)."),
        ("sensitivity_requirement", "str", "standard", "Sensitivity: 'trace', 'standard', 'high_throughput'."),
        ("safety_constraint", "str", "none", "Safety: 'none', 'no_flammable', 'no_high_pressure'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Description: e.g., 'FID detection 30m x 0.25mm column need fast analysis'"),
    ]

    output_sig = [
        ("recommendation", "dict", "Complete carrier gas recommendation including gas selection rationale, optimal flow settings, van Deemter analysis, and safety notes."),
    ]

    examples = [
        {
            "code_input": {
                "detector_type": "FID",
                "column_length_m": 30,
                "column_inner_diameter_mm": 0.25,
                "analysis_goal": "balanced",
            },
            "text_input": {"input_params": "FID 30m x 0.25mm balanced"},
            "output": {
                "recommendation": {
                    "primary_gas": "helium",
                    "optimal_flow_ml_min": 1.1,
                    "optimal_linear_velocity_cm_s": 27,
                    "note": "See full output for details",
                    },
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_viscosity_at_T(self, gas_key: str, T_celsius: float) -> float:
        """Estimate gas viscosity at temperature using Sutherland's formula approximation."""
        base = _GAS_PROPERTIES[gas_key]["viscosity_25c_uPa_s"]
        T_kelvin = T_celsius + 273.15
        T_ref = 298.15
        # Simple power-law approximation for gases
        return base * (T_kelvin / T_ref) ** 0.7

    def _van_deemter_gc(self, gas_key: str, u_cm_s: float, T_celsius: float = 150) -> dict:
        """Calculate H from van Deemter/Golay equation for open tubular column."""
        p = _VAN_DEEMTER_GC[gas_key]
        A = p["A_um"]
        B = p["B_um2_s"]
        C = p["C_s_um"]

        # Temperature correction for B (longitudinal diffusion increases with T)
        T_factor = (T_celsius + 273.15) / 298.15
        B_T = B * T_factor

        u = max(u_cm_s, 0.1)
        H = A + B_T / u + C * u

        # Find optimum
        u_opt = math.sqrt(B_T / C) if C > 0 else 30
        H_min = A + 2 * math.sqrt(B_T * C)

        N_per_meter = 10000 / H_min if H_min > 0 else 0

        return {
            "H_um": round(H, 3),
            "u_cm_s": u,
            "H_minimum_um": round(H_min, 3),
            "optimal_u_cm_s": round(u_opt, 1),
            "plates_per_meter": int(N_per_meter),
            "parameters_used": f"A={A}, B={round(B_T,1)}, C={C}",
        }

    def _calc_flow_rate(self, u_cm_s: float, d_c_mm: float, T_col_c: float,
                         P_out_kpa: float = 101.325, T_room_c: float = 25) -> dict:
        """Convert linear velocity to volumetric flow rate at outlet."""
        r_cm = d_c_mm / 20  # mm → cm
        area_cm2 = math.pi * r_cm ** 2
        F_out = area_cm2 * u_cm_s  # cm³/s at column outlet conditions

        # Correct to room temperature (what flow meter reads)
        F_room = F_out * ((T_room_c + 273.15) / (T_col_c + 273.15))

        # Average pressure correction factor (compressibility)
        # For typical GC pressure drops, use simplified approach
        F_out_ml_min = F_out * 60
        F_room_ml_min = F_room * 60

        # Estimate inlet pressure needed (approximate)
        eta = self._calc_viscosity_at_T("helium", T_col_c)  # default He for estimate
        L_cm = 3000  # will be corrected by caller
        delta_P_est = 128 * eta * L_cm * u_cm_s / (d_c_mm / 10) ** 2 * 0.001  # rough kPa

        return {
            "flow_rate_outlet_ml_min": round(F_out_ml_min, 2),
            "flow_rate_room_temp_ml_min": round(F_room_ml_min, 2),
            "linear_velocity_cm_s": round(u_cm_s, 1),
            "column_cross_section_mm2": round(area_cm2 * 100, 3),
        }

    def _score_gas(self, gas_key: str, detector: str, goal: str,
                   safety: str, sensitivity: str) -> tuple:
        """Score a carrier gas for suitability (0-100)."""
        props = _GAS_PROPERTIES[gas_key]
        score = 50.0

        # Detector compatibility
        det_perf = props["detector_performance"].get(detector, "")
        if "excellent" in det_perf.lower() or "best" in det_perf.lower():
            score += 20
        elif "good" in det_perf.lower():
            score += 10
        elif "poor" in det_perf.lower():
            score -= 25

        # Analysis speed (H₂ is fastest due to high optimal velocity)
        if goal == "speed":
            if gas_key == "hydrogen":
                score += 20
            elif gas_key == "helium":
                score += 5
            elif gas_key == "nitrogen":
                score -= 15  # Very slow
        elif goal == "efficiency":
            if gas_key == "helium":
                score += 15
            elif gas_key == "hydrogen":
                score += 10
            elif gas_key == "nitrogen":
                score += 5  # N₂ has narrow optimum (good if you hit it exactly)

        # Safety
        if safety == "no_flammable":
            if gas_key == "hydrogen":
                score -= 40
            else:
                score += 5

        # Sensitivity
        if sensitivity == "trace":
            if gas_key == "helium":
                score += 10  # Best MS compatibility
            elif gas_key == "hydrogen":
                score -= 5  # May reduce some compounds

        # Availability/cost
        if gas_key == "helium":
            score -= 5  # Supply issues
        elif gas_key == "hydrogen":
            score += 10  # Unlimited availability
        elif gas_key == "nitrogen":
            score += 5  # Cheap and available

        return min(100, max(0, round(score, 1)))

    def _run_base(self, detector_type: str = "FID",
                  column_length_m: float = 30,
                  column_inner_diameter_mm: float = 0.25,
                  film_thickness_um: float = 0.25,
                  analysis_goal: str = "balanced",
                  oven_temperature_c: float = 150,
                  max_pressure_kpa: Optional[float] = None,
                  sensitivity_requirement: str = "standard",
                  safety_constraint: str = "none") -> dict:

        # Score each gas
        scores = {}
        for gas_key in _GAS_PROPERTIES:
            scores[gas_key] = self._score_gas(
                gas_key, detector_type, analysis_goal,
                safety_constraint, sensitivity_requirement
            )

        ranked = sorted(scores.items(), key=lambda x: -x[1])
        primary_gas = ranked[0][0]

        # Calculate optimal conditions for primary gas
        vd = self._van_deemter_gc(primary_gas, 30, oven_temperature_c)
        u_opt = vd["optimal_u_cm_s"]

        # Adjust for goal
        if analysis_goal == "speed":
            u_selected = u_opt * 1.8  # Above optimum for speed (moderate efficiency loss)
        elif analysis_goal == "efficiency":
            u_selected = u_opt * 0.9  # Near optimum
        else:
            u_selected = u_opt

        flow = self._calc_flow_rate(u_selected, column_inner_diameter_mm, oven_temperature_c)

        # Recalculate with actual column length
        L_cm = column_length_m * 100
        total_N = vd["plates_per_meter"] * column_length_m
        actual_vd = self._van_deemter_gc(primary_gas, u_selected, oven_temperature_c)
        tM_min = L_cm / u_selected / 60  # Dead time in minutes

        # Pressure drop estimate (Hagen-Poiseuille approximation)
        eta = self._calc_viscosity_at_T(primary_gas, oven_temperature_c) * 1e-7  # μPa·s → Pa·s
        d_cm = column_inner_diameter_mm / 10
        delta_P_pa = 128 * eta * L_cm * (u_selected * 1e-2) / (d_cm ** 2) / 1e5  # bar
        delta_P_kpa = delta_P_pa * 10

        result = {
            "recommendation": {
                "primary_recommendation": {
                    "gas": primary_gas,
                    "full_name": _GAS_PROPERTIES[primary_gas]["name"],
                    "suitability_score": scores[primary_gas],
                    "rationale": self._get_rationale(primary_gas, detector_type, analysis_goal),
                },
                "ranking_all_gases": [
                    {"gas": g, "score": s, "name": _GAS_PROPERTIES[g]["name"]}
                    for g, s in ranked
                ],
                "optimal_conditions": {
                    "linear_velocity_cm_s": round(u_selected, 1),
                    "flow_rate_ml_min": flow["flow_rate_outlet_ml_min"],
                    "flow_rate_at_detector_ml_min": flow["flow_rate_room_temp_ml_min"],
                    "estimated_inlet_pressure_kpa": round(101.325 + delta_P_kpa, 1),
                    "estimated_pressure_drop_kpa": round(delta_P_kpa, 1),
                    "pressure_ok": (
                        True if max_pressure_kpa is None else
                        (101.325 + delta_P_kpa) <= max_pressure_kpa
                    ),
                },
                "performance_prediction": {
                    "dead_time_tM_min": round(tM_min, 2),
                    "theoretical_plates_N": int(total_N),
                    "plate_height_H_um": actual_vd["H_um"],
                    "analysis_speed_category": (
                        "Fast (<5 min tM)" if tM_min < 5 else
                        "Moderate (5-15 min tM)" if tM_min < 15 else
                        "Slow (>15 min tM)"
                    ),
                },
                "van_deemter_analysis": {
                    "optimal_u_cm_s": vd["optimal_u_cm_s"],
                    "minimum_H_um": vd["H_minimum_um"],
                    "current_H_um": actual_vd["H_um"],
                    "efficiency_at_selected_vs_optimal": f"{actual_vd['H_um']}/{vd['H_minimum_um']} "
                                                       f"({vd['H_minimum_um']/max(actual_vd['H_um'],0.001)*100:.0f}% of optimum)"
                    if actual_vd["H_um"] > 0 else "N/A",
                },
                "gas_properties_reference": _GAS_PROPERTIES[primary_gas],
                "safety_notes": _GAS_PROPERTIES[primary_gas]["safety"],
                "alternative_considerations": {
                    "if_helium_unavailable": (
                        "Hydrogen is best alternative: faster analysis, unlimited supply. "
                        "Ensure proper safety (leak detector, ventilation)." if primary_gas != "hydrogen"
                        else "Consider helium for MS compatibility"
                    ),
                    "if_speed_critical": "Hydrogen recommended regardless of other factors (3-5× faster than N₂)",
                    "if_ms_detection": "Helium preferred; hydrogen acceptable with caution for reducible compounds",
                    "if_ECD_detection": "Nitrogen or helium (avoid hydrogen — higher baseline noise)",
                },
                "column_specifications_used": {
                    "length_m": column_length_m,
                    "inner_diameter_mm": column_inner_diameter_mm,
                    "film_thickness_um": film_thickness_um,
                    "phase_ratio_beta": round(column_inner_diameter_mm / (4 * film_thickness_um), 0),
                },
            }
        }
        return result

    def _get_rationale(self, gas: str, detector: str, goal: str) -> str:
        reasons = {
            "helium": f"Best overall performance for {detector}. Good balance of speed, efficiency, and detector compatibility.",
            "hydrogen": f"Fastest analysis speed ({goal} mode). Excellent TID/FID sensitivity. Ensure flammable gas safety protocols.",
            "nitrogen": "Safest option, low cost. Slower analysis but very flat van Deemter optimum near u_opt. Good for routine QA/QC.",
        }
        return reasons.get(gas, "Selected based on scoring criteria.")

    def _run_text(self, input_params: str) -> dict:
        text = input_params.lower()

        detector = "FID"
        for d in ["FID", "TID", "MS", "ECD", "NPD", "SCD", "NCD"]:
            if d.lower() in text:
                detector = d
                break

        goal = "balanced"
        if "fast" in text or "speed" in text:
            goal = "speed"
        elif "resolution" in text or "efficiency" in text:
            goal = "efficiency"

        safety = "no_flammable" if "safe" in text or "no_flam" in text else "none"

        # Extract column dimensions
        import re
        numbers = [float(x) for x in re.findall(r'[\d.]+', text)]
        L = numbers[0] if len(numbers) >= 1 else 30
        ID = numbers[1] if len(numbers) >= 2 else 0.25

        return self._run_base(detector, L, ID, 0.25, goal, 150, None, "standard", safety)
