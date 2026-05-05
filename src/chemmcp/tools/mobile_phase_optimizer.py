import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ─── Solvent Properties for Mobile Phase Optimization ───
_SOLVENT_PROPERTIES = {
    "water": {"name": "Water", "polarity_index": 10.2, "viscosity_cp": 0.89, "uv_cutoff_nm": 190, "strength_rp": 0.0},
    "acetonitrile": {"name": "Acetonitrile (ACN)", "polarity_index": 5.8, "viscosity_cp": 0.37, "uv_cutoff_nm": 190, "strength_rp": 0.652},
    "methanol": {"name": "Methanol (MeOH)", "polarity_index": 5.1, "viscosity_cp": 0.55, "uv_cutoff_nm": 205, "strength_rp": 0.729},
    "ethanol": {"name": "Ethanol (EtOH)", "polarity_index": 4.3, "viscosity_cp": 1.08, "uv_cutoff_nm": 210, "strength_rp": 0.688},
    "isopropanol": {"name": "Isopropanol (IPA)", "polarity_index": 3.9, "viscosity_cp": 2.05, "uv_cutoff_nm": 210, "strength_rp": 0.606},
    "tetrahydrofuran": {"name": "Tetrahydrofuran (THF)", "polarity_index": 4.0, "viscosity_cp": 0.46, "uv_cutoff_nm": 212, "strength_rp": 0.450},
    "dichloromethane": {"name": "Dichloromethane (DCM)", "polarity_index": 3.1, "viscosity_cp": 0.41, "uv_cutoff_nm": 233, "strength_rp": 0.309},
    "hexane": {"name": "n-Hexane", "polarity_index": 0.1, "viscosity_cp": 0.30, "uv_cutoff_nm": 201, "strength_rp": 0.01},
    "acetone": {"name": "Acetone", "polarity_index": 5.1, "viscosity_cp": 0.31, "uv_cutoff_nm": 330, "strength_rp": 0.564},
    "dmso": {"name": "DMSO", "polarity_index": 6.5, "viscosity_cp": 2.00, "uv_cutoff_nm": 268, "strength_rp": 0.612},
}

# Common buffer systems
_BUFFER_SYSTEMS = {
    "formic_acid": {"name": "Formic acid", "pKa": 3.75, "range": (2.5, 4.5), "conc_mM": [5, 10, 20], "ms_compatible": True, "notes": "Most common for LC-MS"},
    "acetic_acid": {"name": "Acetic acid / Ammonium acetate", "pKa": 4.76, "range": (3.8, 5.8), "conc_mM": [5, 10, 20], "ms_compatible": True, "notes": "Good for ESI-MS"},
    "ammonium_formate": {"name": "Ammonium formate", "pKa_formic": 3.75, "range": (3.0, 5.0), "conc_mM": [5, 10, 20, 50], "ms_compatible": True, "notes": "Volatile; excellent for MS"},
    "ammonium_bicarbonate": {"name": "Ammonium bicarbonate", "pKa": 10.3, "range": (8.0, 11.0), "conc_mM": [5, 10, 20], "ms_compatible": True, "notes": "Basic pH volatile buffer"},
    "phosphate": {"name": "Phosphate buffer", "pKa2": 7.21, "range": (5.8, 8.0), "conc_mM": [10, 25, 50], "ms_compatible": False, "notes": "NOT MS-compatible; good UV detection"},
    "trifluoroacetic_acid": {"name": "Trifluoroacetic acid (TFA)", "pKa": 0.3, "range": (0.5, 2.0), "conc_mM": [0.1, 0.5, 1], "ms_compatible": True, "notes": "Ion-pairing agent; suppresses MS signal"},
    "tfa_ion_pair": {"name": "TFA (ion-pairing mode)", "pKa": 0.3, "range": (0.5, 2.0), "conc_mM": [1, 5, 10], "ms_compatible": False, "notes": "Strong ion-pairing for basic compounds"},
}


@ChemMCPManager.register_tool
class MobilePhaseOptimizer(BaseTool):
    """
    流动相组成和梯度优化工具。
    推荐HPLC/GC流动相组成、缓冲体系、pH、梯度程序，并预测分离效果。
    """
    __version__ = "0.1.0"
    name = "MobilePhaseOptimizer"
    func_name = "optimize_mobile_phase"
    description = "Optimize HPLC mobile phase composition, buffer system, pH, and gradient program for optimal chromatographic separation."
    implementation_description = (
        "Uses solvent strength parameters (Snyder's P'), viscosity calculations, "
        "buffer selection rules, and gradient optimization heuristics to recommend "
        "optimal mobile phase conditions."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["HPLC", "Mobile Phase", "Gradient Optimization", "Method Development", "Chromatography"]
    required_envs = []

    code_input_sig = [
        ("stationary_phase", "str", "C18", "Stationary phase: 'C18', 'C8', 'Phenyl', 'HILIC', etc."),
        ("analyte_polarity", "str", "moderate", "Analyte polarity: 'non_polar', 'moderate', 'polar', 'ionic'."),
        ("target_compounds", "list", "N/A", "List of compound classes or specific properties."),
        ("detection_method", "str", "uv", "Detection: 'uv', 'ms', 'rid', 'cad', 'els'."),
        ("current_conditions", "dict", "{}", "Current method to improve upon (optional)."),
        ("optimization_goal", "str", "resolution", "Goal: 'resolution', 'speed', 'peak_shape', 'loadability'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Description: e.g., 'C18 column separate polar drugs need MS detection optimize resolution'"),
    ]

    output_sig = [
        ("optimization", "dict", "Complete mobile phase recommendation including organic modifier, buffer, pH, gradient program, and expected outcomes."),
    ]

    examples = [
        {
            "code_input": {
                "stationary_phase": "C18",
                "analyte_polarity": "moderate",
                "target_compounds": ["basic_drugs"],
                "detection_method": "ms",
                "optimization_goal": "resolution",
            },
            "text_input": {"input_params": "C18 basic drugs MS detection need resolution"},
            "output": {
                "optimization": {
                    "organic_modifier": "Acetonitrile (ACN)",
                    "buffer": "10 mM ammonium formate + 0.1% formic acid pH 3.5",
                    "gradient_program": [...],
                    "expected_outcome": "See details",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_elution_strength(self, org_percent: float, solvent_key: str) -> float:
        """Calculate effective elution strength."""
        s = _SOLVENT_PROPERTIES.get(solvent_key, _SOLVENT_PROPERTIES["acetonitrile"])
        return org_percent / 100.0 * s["strength_rp"]

    def _calc_viscosity(self, aq_frac: float, org_key: str) -> float:
        """Estimate binary mixture viscosity (log rule)."""
        eta_aq = _SOLVENT_PROPERTIES["water"]["viscosity_cp"]
        eta_org = _SOLVENT_PROPERTIES.get(org_key, _SOLVENT_PROPERTIES["acetonitrile"])["viscosity_cp"]
        # Logarithmic mixing rule
        if aq_frac <= 0 or aq_frac >= 1:
            return eta_org if aq_frac == 0 else eta_aq
        vis = aq_frac * math.log(eta_aq) + (1 - aq_frac) * math.log(eta_org)
        return round(math.exp(vis), 3)

    def _select_buffer(self, target_ph: float, detection: str, analyte_type: str = "neutral") -> dict:
        """Select appropriate buffer system."""
        candidates = []
        for buf_key, buf_data in _BUFFER_SYSTEMS.items():
            p_range = buf_data["range"]
            if p_range[0] <= target_ph <= p_range[1]:
                score = 50
                if detection == "ms" and not buf_data["ms_compatible"]:
                    score -= 40
                elif detection == "ms" and buf_data["ms_compatible"]:
                    score += 20
                if analyte_type == "basic" and "acid" in buf_key.lower():
                    score += 15
                if analyte_type == "acidic" and "base" in buf_key.lower() or "bicarb" in buf_key.lower():
                    score += 15
                candidates.append((buf_key, buf_data, score))

        candidates.sort(key=lambda x: -x[2])
        return candidates[0] if candidates else list(_BUFFER_SYSTEMS.items())[0]

    def _run_base(self, stationary_phase: str, analyte_polarity: str = "moderate",
                  target_compounds: Optional[List[str]] = None,
                  detection_method: str = "uv",
                  current_conditions: Optional[Dict] = None,
                  optimization_goal: str = "resolution") -> dict:
        """Core optimization logic."""
        target_compounds = target_compounds or []
        current_conditions = current_conditions or {}

        sp = stationary_phase.upper()
        is_hilic = sp in ("HILIC",)
        is_np = sp in ("CYANO", "CN", "NORMAL")

        # ── Organic Modifier Selection ──
        if is_hilic:
            org_key = "acetonitrile"  # ACN preferred for HILIC
            org_name = "Acetonitrile (ACN)"
            starting_b = 90  # High organic start for HILIC
            ending_b = 50
        elif is_np:
            org_key = "hexane"  # NP uses non-polar solvents
            org_name = "n-Hexane/IPA"
            starting_b = 95
            ending_b = 50
        else:
            # RP mode
            if analyte_polarity == "non_polar":
                org_key = "methanol"  # MeOH stronger for non-polar
                org_name = "Methanol (MeOH)"
            elif analyte_polarity == "polar":
                org_key = "acetonitrile"  # ACN better selectivity for polar
                org_name = "Acetonitrile (ACN)"
            else:
                org_key = "acetonitrile"  # Default ACN
                org_name = "Acetonitrile (ACN)"
            starting_b = 5
            ending_b = 95

        # ── Buffer Selection ──
        comp_text = " ".join(target_compounds).lower()
        analyte_type = "neutral"
        if any(w in comp_text for w in ["basic", "amine", "base"]):
            analyte_type = "basic"
        elif any(w in comp_text for w in ["acidic", "acid", "carboxyl"]):
            analyte_type = "acidic"

        if analyte_type == "basic":
            target_ph = 3.0  # Suppress silanol activity
        elif analyte_type == "acidic":
            target_ph = 5.5  # Keep acids protonated
        else:
            target_ph = 3.5  # Default slightly acidic

        buf_key, buf_data, buf_score = self._select_buffer(target_ph, detection_method, analyte_type)

        # ── Gradient Program ──
        if optimization_goal == "speed":
            grad_time = 5.0
            grad_steps = [(0, starting_b), (grad_time, ending_b)]
            hold_time = 1.0
            reeq_time = 1.0
        elif optimization_goal == "resolution":
            grad_time = 20.0
            grad_steps = [(0, starting_b), (2, starting_b), (grad_time - 2, ending_B := ending_b), (grad_time, ending_B)]
            hold_time = 3.0
            reeq_time = 5.0
        elif optimization_goal == "peak_shape":
            grad_time = 15.0
            grad_steps = [(0, starting_b), (1, starting_b), (grad_time, ending_b)]
            hold_time = 5.0
            reeq_time = 5.0
        else:  # general / loadability
            grad_time = 10.0
            grad_steps = [(0, starting_b), (grad_time, ending_b)]
            hold_time = 2.0
            reeq_time = 3.0

        total_run = grad_time + hold_time + reeq_time

        # Viscosity at midpoint
        mid_b = (starting_b + ending_b) / 2
        if is_hilic:
            visc = self._calc_viscosity(mid_b / 100, org_key)
        else:
            visc = self._calc_viscosity(1 - mid_b / 100, org_key)

        result = {
            "optimization": {
                "stationary_phase": stationary_phase,
                "mode": "HILIC" if is_hilic else ("Normal Phase" if is_np else "Reversed-Phase"),
                "organic_modifier": {
                    "solvent": org_name,
                    "key": org_key,
                    "elution_strength_P": _SOLVENT_PROPERTIES[org_key]["strength_rp"],
                    "viscosity_cp": _SOLVENT_PROPERTIES[org_key]["viscosity_cp"],
                    "uv_cutoff_nm": _SOLVENT_PROPERTIES[org_key]["uv_cutoff_nm"],
                    "alternative": "Methanol (MeOH)" if org_key == "acetonitrile" else "Acetonitrile (ACN)",
                },
                "aqueous_phase": {
                    "buffer_system": buf_data["name"],
                    "buffer_key": buf_key,
                    "recommended_ph": target_ph,
                    "concentration_mM": buf_data["conc_mM"][1],
                    "ms_compatible": buf_data["ms_compatible"],
                    "preparation_note": f"Dissolve {buf_data['conc_mM'][1]} mM buffer in water, adjust to pH {target_ph}",
                },
                "gradient_program": {
                    "total_time_min": round(total_run, 1),
                    "gradient_time_min": round(grad_time, 1),
                    "hold_time_min": hold_time,
                    "re_equilibration_time_min": reeq_time,
                    "steps": [{"time_min": t, "percent_b": b} for t, b in grad_steps],
                    "flow_rate_suggestion_ml_min": "0.3-0.5 mL/min (2.1mm ID) or 0.8-1.2 mL/min (4.6mm ID)",
                },
                "physical_properties": {
                    "estimated_viscosity_at_midpoint_cp": visc,
                    "backpressure_note": (
                        "Low backpressure (ACN/water)" if visc < 0.6 else
                        "Moderate backpressure (MeOH/water)" if visc < 0.8 else
                        "Higher backpressure — consider temperature increase"
                    ),
                    "temperature_recommendation_c": 40 if visc > 0.7 else 30,
                },
                "detection_considerations": {
                    "method": detection_method.upper(),
                    "compatible": "✓ Yes" if not (
                        (detection_method == "ms" and not buf_data["ms_compatible"]) or
                        (detection_method == "uv" and _SOLVENT_PROPERTIES[org_key]["uv_cutoff_nm"] > 220)
                    ) else "⚠ Check compatibility",
                    "wavelength_note": (
                        f"UV cutoff of {org_name} is {_SOLVENT_PROPERTIES[org_key]['uv_cutoff_nm']} nm"
                        if detection_method == "uv" else "N/A"
                    ),
                },
                "troubleshooting_tips": self._get_tips(analyte_polarity, analyte_type, optimization_goal),
                "improvement_over_current": self._compare_with_current(current_conditions, org_key, target_ph) if current_conditions else None,
            }
        }
        return result

    def _get_tips(self, polarity: str, atype: str, goal: str) -> List[str]:
        tips = []
        if atype == "basic":
            tips.extend([
                "For basic compounds: consider using charged surface hybrid (CSH) columns to reduce tailing",
                "Add 0.1% formic acid or use ammonium formate buffer pH ~3.0 to suppress silanol effects",
                "If tailing persists: try 10-20 mM ammonium bicarbonate pH 10 for very basic analytes on stable columns",
            ])
        if polarity == "polar":
            tips.append("Very polar compounds may show poor retention on RP — consider HILIC or ion-pairing")
        if goal == "resolution":
            tips.extend([
                "Try shallower gradients around the region of co-elution",
                "Consider adjusting column temperature ±10°C to fine-tune selectivity",
                "Alternative organic modifiers (THF, MeOH) can change selectivity significantly",
            ])
        if goal == "speed":
            tips.append("Shorter columns (50 mm) with sub-2μm particles can cut run time by 3-5×")
        return tips[:5]

    def _compare_with_current(self, current: Dict, org_key: str, ph: float) -> Dict[str, str]:
        improvements = []
        if "organic" in current:
            cur_org = current["organic"].lower()
            if "methanol" in cur_org and org_key == "acetonitrile":
                improvements.append("Switch MeOH → ACN: lower viscosity → lower pressure or higher efficiency")
            elif "acn" in cur_org and org_key == "methanol":
                improvements.append("Switch ACN → MeOH: different selectivity may resolve co-elution")
        if "ph" in current:
            try:
                cur_ph = float(current["ph"])
                if abs(cur_ph - ph) > 1.5:
                    improvements.append(f"Adjust pH {cur_ph} → {ph}: changes ionization → better peak shape/retention")
            except (ValueError, TypeError):
                pass
        return improvements if improvements else ["No major changes suggested from current method"]

    def _run_text(self, input_params: str) -> dict:
        """Parse natural language input."""
        text = input_params.lower()

        sp = "C18"
        for phase in ["C18", "C8", "Phenyl", "PFP", "HILIC", "Cyano"]:
            if phase.lower() in text:
                sp = phase
                break

        polarity = "moderate"
        if "non.polar" in text or "hydrophob" in text:
            polarity = "non_polar"
        elif "polar" in text:
            polarity = "polar"

        det = "uv"
        if "ms" in text or "mass.spec" in text:
            det = "ms"

        goal = "resolution"
        if "fast" in text or "speed" in text:
            goal = "speed"
        elif "shape" in text or "tailing" in text:
            goal = "peak_shape"

        targets = []
        if any(w in text for w in ["drug", "pharmaceutical"]):
            targets.append("drugs")
        if any(w in text for w in ["basic", "amine"]):
            targets.append("basic_compounds")
        if any(w in text for w in ["acidic", "carboxylic"]):
            targets.append("acidic_compounds")

        return self._run_base(sp, polarity, targets, det, {}, goal)
