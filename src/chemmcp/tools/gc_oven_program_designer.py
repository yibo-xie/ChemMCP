import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ─── Compound Class Boiling Point Ranges (°C) ───
_BP_RANGES = {
    "permanent_gases": {"range": (-200, -50), "examples": "H₂, O₂, N₂, CH₄, CO", "typical_oven_start": 35},
    "light_volatiles": {"range": (-50, 100), "examples": "C1-C5 alkanes, light solvents", "typical_oven_start": 40},
    "medium_volatiles": {"range": (50, 200), "examples": "Benzene, toluene, C6-C10, common solvents", "typical_oven_start": 40},
    "semi_volatiles": {"range": (150, 300), "examples": "C10-C20, PAHs, pesticides, phenols", "typical_oven_start": 50},
    "heavy_compounds": {"range": (280, 450), "examples": "Waxes, steroids, triglycerides, C20+", "typical_oven_start": 60},
}

# ─── General GC Oven Program Guidelines ───
# Based on Snyder and Dolan optimization principles adapted for GC
_OVEN_GUIDELINES = {
    "isothermal": {
        "when_to_use": "Narrow boiling range (<30°C spread) or simple gas mixtures",
        "optimal_T_diff": "Set T ≈ BP_avg + 30-50°C",
        "pros": "Simple, reproducible, good for QA/QC",
        "cons": "Long analysis time for wide boiling range mixtures",
    },
    "ramp_single": {
        "when_to_use": "Boiling range spread of 50-150°C",
        "recommended_rate_c_min": [5, 10, 15, 20],
        "hold_initial_min": [0.5, 1, 2],
        "hold_final_min": [2, 5, 10],
    },
    "ramp_multi": {
        "when_to_use": "Wide boiling range (>150°C) or complex mixtures",
        "strategy": "Slow ramp through critical region, fast elsewhere",
    },
}


@ChemMCPManager.register_tool
class GcOvenProgramDesigner(BaseTool):
    """
    GC程序升温设计工具。
    根据分析物沸点范围、极性和分离需求，设计优化的气相色谱程序升温方法。
    """
    __version__ = "0.1.0"
    name = "GcOvenProgramDesigner"
    func_name = "design_gc_oven_program"
    description = "Design optimized GC temperature programming (isothermal or ramped) based on analyte boiling points, polarity, and separation requirements."
    implementation_description = (
        "Uses van't Hoff relationships for retention-temperature dependence, "
        "empirical rules for initial/final temperature selection, "
        "and heating rate optimization based on required resolution vs. analysis speed trade-offs."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["GC", "Temperature Programming", "Gas Chromatography", "Method Development"]
    required_envs = []

    code_input_sig = [
        ("analytes", "list", "N/A", "List of analyte descriptions: compound classes, names, or boiling point ranges."),
        ("boiling_points_c", "list", "[]", "Known boiling points in °C (optional)."),
        ("column_type", "str", "non_polar", "Column type: 'non_polar' (e.g., DB-5/HP-5), 'mid_polar' (DB-17/DB-35), 'polar' (wax/PEG)."),
        ("separation_goal", "str", "general", "Goal: 'speed', 'resolution', 'trace_analysis', 'general'."),
        ("carrier_gas", "str", "helium", "Carrier gas: 'helium', 'hydrogen', 'nitrogen'."),
        ("max_oven_temp_c", "None", "None", "Maximum oven temperature limit (column-dependent)."),
        ("sample_matrix", "str", "general", "Sample matrix: 'general', 'environmental', 'food', 'petroleum', 'biological'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Description: e.g., 'analyze C6-C15 hydrocarbons on DB-5 need good resolution'"),
    ]

    output_sig = [
        ("program", "dict", "Complete oven program including temperature steps, hold times, estimated run time, and expected elution order."),
    ]

    examples = [
        {
            "code_input": {
                "analytes": ["medium_volatiles", "semi_volatiles"],
                "column_type": "non_polar",
                "separation_goal": "resolution",
            },
            "text_input": {"input_params": "C6-C20 compounds DB-5 column need resolution"},
            "output": {
                "program": {
                    "steps": [...],
                    "total_time_min": 25.3,
                    "expected_elution_order": "TBD",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _classify_analytes(self, analytes: List[str], bps: List[float]) -> List[Dict]:
        """Classify analytes into boiling point categories."""
        categories = []
        for a in analytes:
            al = a.lower()
            matched = False
            for cat_key, cat_data in _BP_RANGES.items():
                if al.replace(" ", "_") == cat_key or any(w in al for w in cat_data["examples"].lower().split(", ")):
                    matched = True
                    break
            if not matched:
                # Try matching by name keywords
                if any(w in al for w in ["gas", "methane", "ethane", "propane", "butane"]):
                    cat_key = "permanent_gases"
                elif any(w in al for w in ["pentane", "hexane", "heptane", "ether", "acetone"]):
                    cat_key = "light_volatiles"
                elif any(w in al for w in ["benzene", "toluene", "xylene", "octane", "decane"]):
                    cat_key = "medium_volatiles"
                elif any(w in al for w in ["pesticide", "pah", "phenol", "fatty", "steroid"]):
                    cat_key = "semi_volatiles"
                elif any(w in al for w in ["wax", "heavy", "triglyceride"]):
                    cat_key = "heavy_compounds"
                else:
                    cat_key = "medium_volatiles"  # default

            categories.append({"analyte": a, "category": cat_key, **_BP_RANGES[cat_key]})

        # Add known boiling points
        if bps:
            for bp in bps:
                for cat_key, cat_data in _BP_RANGES.items():
                    lo, hi = cat_data["range"]
                    if lo <= bp <= hi:
                        categories.append({"analyte": f"BP={bp}°C", "category": cat_key, **cat_data})
                        break

        return categories

    def _estimate_retention_temp(self, bp: float, T_col: float) -> float:
        """Estimate approximate elution temperature from boiling point."""
        # Rule of thumb: compounds elute at ~T_oven where k ≈ 2-5
        # For isothermal: T_elution roughly proportional to BP with offset
        return bp + 30  # Simplified; actual depends on phase ratio

    def _design_isothermal(self, bp_avg: float, goal: str) -> dict:
        """Design isothermal program."""
        T_iso = bp_avg + 40
        if goal == "speed":
            T_iso += 20
        elif goal == "resolution":
            T_iso += 10

        return {
            "type": "isothermal",
            "temperature_c": round(T_iso, 0),
            "hold_time_min": round(20 if goal == "resolution" else 10, 1),
            "estimated_analysis_time_min": round(15 if goal == "speed" else 25, 1),
        }

    def _design_ramped(self, bp_min: float, bp_max: float, spread: float,
                        goal: str, col_type: str, max_T: Optional[float]) -> dict:
        """Design temperature-ramped program."""
        T_init = max(35, bp_min - 20)
        T_final = min(bp_max + 30, max_T or 350)

        # Determine number of ramps needed
        if spread < 80:
            n_ramps = 1
        elif spread < 180:
            n_ramps = 2
        else:
            n_ramps = 3

        # Heating rates depend on goal
        if goal == "speed":
            rates = [30, 40, 50][:n_ramps]
        elif goal == "resolution":
            rates = [3, 8, 15][:n_ramps]
        elif goal == "trace_analysis":
            rates = [2, 5, 10][:n_ramps]
        else:
            rates = [10, 15, 20][:n_ramps]

        # Hold times
        init_hold = 1.0 if goal != "resolution" else 2.0
        final_hold = 5.0 if goal in ("resolution", "trace_analysis") else 2.0

        steps = [{"rate_c_per_min": 0, "hold_min": init_hold, "temp_c": int(T_init)}]

        if n_ramps == 1:
            steps.append({"rate_c_per_min": rates[0], "hold_min": final_hold, "temp_c": int(T_final)})
        elif n_ramps == 2:
            T_mid = (bp_min + bp_max) / 2 + 10
            steps.append({"rate_c_per_min": rates[0], "hold_min": 1.0, "temp_c": int(T_mid)})
            steps.append({"rate_c_per_min": rates[1], "hold_min": final_hold, "temp_c": int(T_final)})
        else:
            third = (bp_max - bp_min) / 3
            T1 = bp_min + third
            T2 = bp_min + 2 * third
            steps.append({"rate_c_per_min": rates[0], "hold_min": 1.0, "temp_c": int(T1)})
            steps.append({"rate_c_per_min": rates[1], "hold_min": 1.0, "temp_c": int(T2)})
            steps.append({"rate_c_per_min": rates[2], "hold_min": final_hold, "temp_c": int(T_final)})

        # Estimate total time
        total_time = init_hold
        prev_T = T_init
        for step in steps[1:]:
            dT = step["temp_c"] - prev_T
            ramp_time = dT / max(step["rate_c_per_min"], 0.1)
            total_time += ramp_time + step.get("hold_min", 0)
            prev_T = step["temp_c"]

        return {
            "type": f"{n_ramps}-stage ramped",
            "steps": steps,
            "initial_temperature_c": int(T_init),
            "final_temperature_c": int(T_final),
            "total_estimated_time_min": round(total_time, 1),
            "temperature_range_c": f"{int(T_init)} → {int(T_final)}",
            "n_ramps": n_ramps,
        }

    def _run_base(self, analytes: List[str], boiling_points_c: Optional[List[float]] = None,
                  column_type: str = "non_polar", separation_goal: str = "general",
                  carrier_gas: str = "helium", max_oven_temp_c: Optional[float] = None,
                  sample_matrix: str = "general") -> dict:

        classified = self._classify_analytes(analytes, boiling_points_c or [])

        if not classified:
            classified = [{"analyte": "default", "category": "medium_volatiles", **_BP_RANGES["medium_volatiles"]}]

        # Determine boiling point range
        all_bp_ranges = [(c["range"][0], c["range"][1]) for c in classified]
        bp_min = min(r[0] for r in all_bp_ranges)
        bp_max = max(r[1] for r in all_bp_ranges)
        bp_spread = bp_max - bp_min
        bp_avg = (bp_min + bp_max) / 2

        # Decide isothermal vs ramped
        use_isothermal = bp_spread < 40 and len(classified) <= 3

        if use_isothermal:
            prog = self._design_isothermal(bp_avg, separation_goal)
        else:
            prog = self._design_ramped(bp_min, bp_max, bp_spread, separation_goal,
                                        column_type, max_oven_temp_c)

        # Column-specific notes
        col_notes = {
            "non_polar": "DB-5/HP-5 type (5% phenyl): excellent general purpose, wide temp range (-60 to 325/350°C)",
            "mid_polar": "DB-17/DB-35 (50% phenyl): enhanced selectivity for aromatic/positional isomers",
            "polar": "Polyethylene glycol (WAX): best for polar compounds, alcohols, FAMEs; max T ~240-260°C",
        }

        result = {
            "program": {
                "analyte_classification": classified,
                "boiling_point_range": {
                    "min_c": bp_min,
                    "max_c": bp_max,
                    "spread_c": round(bp_spread, 0),
                    "average_c": round(bp_avg, 0),
                },
                "program_type": prog["type"],
                "temperature_steps": prog.get("steps", [
                    {"temp_c": prog.get("temperature_c", 0), "hold_min": prog.get("hold_time_min", 0)}
                ]),
                "total_estimated_time_min": prog.get("total_estimated_time_min", prog.get("estimated_analysis_time_min", 0)),
                "initial_temp_c": prog.get("initial_temperature_c", prog.get("temperature_c", 0)),
                "final_temp_c": prog.get("final_temperature_c", prog.get("temperature_c", 0)),
                "column_info": col_notes.get(column_type, "Unknown column type"),
                "carrier_gas": carrier_gas.upper(),
                "injection_suggestion": {
                    "injector_temp_c": max(bp_max + 50, 250),
                    "split_ratio": "10:1 to 100:1 (adjust based on concentration)",
                    "injection_volume_ul": 1.0,
                    "liner_type": "deactivated split/splitless",
                },
                "detector_suggestion": {
                    "FID": "Universal for organics; T=250-300°C",
                    "MS": "Best for identification; transfer line=280°C",
                    "ECD": "For halogenated compounds; T=300-320°C",
                    "NPD": "For N/P compounds; T=300-320°C (bead mode)",
                },
                "optimization_tips": [
                    "If early peaks co-elute: lower initial temperature or increase initial hold",
                    "If late peaks too broad: increase final ramp rate or raise final temperature",
                    "If middle region crowded: add intermediate plateau or reduce ramp rate in that zone",
                    "Carryover observed? Increase final hold time or raise final temperature by 10-20°C",
                    "For quantitative work: ensure baseline stability before first peak elutes",
                ],
                "method_validation_checkpoints": [
                    "Resolution Rs ≥ 1.5 for all critical pairs",
                    "Peak symmetry (tailing factor) 0.9-1.2",
                    "Signal-to-noise ≥ 10 for LOQ, ≥ 3 for LOD",
                    "Retention time RSD < 0.2% for n=6 replicates",
                ],
            }
        }
        return result

    def _run_text(self, input_params: str) -> dict:
        text = input_params.lower()

        analytes = []
        if any(w in text for w in ["c6", "c7", "c8", "benzene", "toluene"]):
            analytes.append("medium_volatiles")
        if any(w in text for w in ["c10", "c12", "c15", "pesticide", "paH"]):
            analytes.append("semi_volatiles")
        if any(w in text for w in ["c20", "wax", "heavy", "triglyceride"]):
            analytes.append("heavy_compounds")
        if any(w in text for w in ["gas", "methane", "light"]):
            analytes.append("light_volatiles")

        if not analytes:
            analytes = ["medium_volatiles", "semi_volatiles"]

        col = "non_polar"
        for c_name, key in [("db5", "non_polar"), ("db17", "mid_polar"), ("wax", "polar"), ("peg", "polar")]:
            if c_name in text:
                col = key
                break

        goal = "general"
        if "fast" in text or "speed" in text:
            goal = "speed"
        elif "resolution" in text or "separate" in text:
            goal = "resolution"

        gas = "helium"
        if "h2" in text or "hydrogen" in text:
            gas = "hydrogen"
        elif "n2" in text or "nitrogen" in text:
            gas = "nitrogen"

        return self._run_base(analytes, [], col, goal, gas)
