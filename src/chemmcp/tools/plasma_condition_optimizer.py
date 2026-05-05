"""
Plasma Condition Optimizer — 等离子体操作条件优化工具 (#328)

功能：
  根据待测元素和样品类型，优化ICP（ICP-OES/ICP-MS）操作参数，
  包括RF功率、雾化器气流量、辅助气流量、冷却气流量和采样深度等。
"""

import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 元素最佳等离子体条件参考数据库 ────────────────────────────
# 数据来源: 各仪器厂商应用手册 (Agilent, Thermo Fisher, PerkinElmer)
# 以及文献推荐值
ELEMENT_PLASMA_DB: Dict[str, Dict[str, Any]] = {
    # ── 易电离元素 (需要较低RF功率) ──
    "li": {"rf_optimal": 1200, "neb_opt": 0.90, "aux_opt": 0.9, "notes": "Low IP; reduce power to suppress ionization."},
    "na": {"rf_optimal": 1150, "neb_opt": 0.85, "aux_opt": 0.8, "notes": "Very low IP (5.14 eV); use cool plasma for Na/K."},
    "k":  {"rf_optimal": 1150, "neb_opt": 0.85, "aux_opt": 0.8, "notes": "Low IP; cool plasma recommended."},
    "rb": {"rf_optimal": 1150, "neb_opt": 0.85, "aux_opt": 0.8, "notes": "Very low IP (4.18 eV); cool plasma."},
    "cs": {"rf_optimal": 1100, "neb_opt": 0.80, "aux_opt": 0.8, "notes": "Extremely low IP (3.89 eV); requires cool plasma."},
    "ca": {"rf_optimal": 1300, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Moderate conditions; watch for oxide formation."},
    "sr": {"rf_optimal": 1300, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Similar to Ca; ionic lines benefit from higher power."},
    "ba": {"rf_optimal": 1350, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "Higher power needed for ionization."},

    # ── 过渡金属 ──
    "cr": {"rf_optimal": 1400, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Standard conditions work well."},
    "mn": {"rf_optimal": 1350, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Excellent sensitivity across wide power range."},
    "fe": {"rf_optimal": 1350, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Robust signal; standard conditions optimal."},
    "co": {"rf_optimal": 1400, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Standard conditions."},
    "ni": {"rf_optimal": 1400, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Standard conditions."},
    "cu": {"rf_optimal": 1350, "neb_opt": 0.90, "aux_opt": 0.9, "notes": "High sensitivity; lower nebulizer flow improves transport."},
    "zn": {"rf_optimal": 1350, "neb_opt": 0.90, "aux_opt": 0.9, "notes": "Lower nebulizer flow improves Zn sensitivity significantly."},

    # ── 难熔/高电离能元素 (需要较高功率) ──
    "ti": {"rf_optimal": 1500, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "Refractory; needs high power for complete atomization/ionization."},
    "v":  {"rf_optimal": 1500, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "Refractory oxides; high power essential."},
    "mo": {"rf_optimal": 1500, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "High first IP; benefits from high RF power."},
    "w":  {"rf_optimal": 1550, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "Very refractory; highest power recommended."},
    # ── 类金属 / 其他 ──
    "as": {"rf_optimal": 1450, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Hydride generation alternative for better DL."},
    "se": {"rf_optimal": 1450, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Hydride generation alternative; Ar2 interference at m/z 80."},
    "cd": {"rf_optimal": 1350, "neb_opt": 0.90, "aux_opt": 0.9, "notes": "Standard trace element conditions."},
    "sb": {"rf_optimal": 1400, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Standard conditions."},
    "pb": {"rf_optimal": 1400, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Standard conditions; widely used environmental analysis."},
    "hg": {"rf_optimal": 1300, "neb_opt": 0.85, "aux_opt": 0.8, "notes": "Cold vapor AAS or ICP-MS with gold trap preferred."},
    "be": {"rf_optimal": 1400, "neb_opt": 0.95, "aux_opt": 0.9, "notes": "Toxic! Standard ICP conditions."},
    "al": {"rf_optimal": 1450, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "Refractory Al2O3; higher power needed."},
    "si": {"rf_optimal": 1500, "neb_opt": 1.00, "aux_opt": 1.0, "notes": "Very refractory SiO2; high power required."},
}


# ── 样品类型默认条件调整 ──────────────────────────────────────
SAMPLE_TYPE_ADJUSTMENTS: Dict[str, Dict[str, Any]] = {
    "aqueous": {
        "rf_power_offset": 0,
        "neb_flow_offset": 0.0,
        "aux_flow_offset": 0.0,
        "coolant_flow": 15.0,   # L/min (Ar)
        "sampling_depth_mm": 8.0,
        "notes": "Standard aqueous solution conditions.",
    },
    "organic": {
        "rf_power_offset": +100,  # Higher power to burn organic solvent
        "neb_flow_offset": -0.10,  # Lower nebulizer flow for organic solvents
        "aux_flow_offset": +0.2,   # More auxiliary gas to prevent carbon buildup
        "coolant_flow": 14.0,      # Slightly reduced coolant
        "sampling_depth_mm": 10.0,  # Deeper sampling position
        "notes": "Use oxygen ashing (1–3% O2 in auxiliary gas) if available. "
                 "Consider chilled spray chamber for volatile organics.",
    },
    "high_matrix": {
        "rf_power_offset": +50,
        "neb_flow_offset": +0.05,
        "aux_flow_offset": 0.0,
        "coolant_flow": 15.0,
        "sampling_depth_mm": 7.0,
        "notes": "High TDS (>0.2%): use wider bore injector, consider internal standards, "
                 "dilute sample if possible.",
    },
    "environmental": {
        "rf_power_offset": 0,
        "neb_flow_offset": 0.0,
        "aux_flow_offset": 0.0,
        "coolant_flow": 15.0,
        "sampling_depth_mm": 8.0,
        "notes": "Environmental waters: check for dissolved solids, use appropriate internal standards.",
    },
}


@ChemMCPManager.register_tool
class PlasmaConditionOptimizer(BaseTool):
    """
    等离子体操作条件优化工具。
    基于待测元素和样品类型，推荐最优ICP操作参数。
    """
    __version__                = "0.1.0"
    name                       = "PlasmaConditionOptimizer"
    func_name                  = "optimize_plasma_conditions"
    description                = ("Optimize ICP operating conditions (RF power, gas flows, sampling depth) "
                                 "based on analyte elements and sample type.")
    implementation_description = ("Uses built-in database of element-specific optimal plasma conditions and "
                                 "sample-type adjustments to recommend a full set of ICP operating parameters.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["ICP", "Plasma Optimization", "ICP-OES", "ICP-MS",
                                   "Instrument Parameters", "Analytical Chemistry"]
    required_envs              = []

    code_input_sig = [
        ("analyte_elements",            "list",  "N/A",       "List of element symbols to analyze (e.g., ['Fe', 'Cu', 'Zn'])."),
        ("sample_type",                 "str",   "aqueous",   "Sample type: 'aqueous', 'organic', 'high_matrix', or 'environmental'."),
        ("optimization_goal",           "str",   "sensitivity", "Goal: 'sensitivity', 'stability', 'oxide_ratio', or 'doubly_charged'."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated: elem1,elem2,... [sample_type] [optimization_goal]"),
    ]

    output_sig = [
        ("rf_power_W",                  "float", "Recommended RF power (Watts)."),
        ("nebulizer_flow_Lmin",         "float", "Nebulizer gas flow rate (L/min Ar)."),
        ("auxiliary_flow_Lmin",         "float", "Auxiliary/plasma gas flow (L/min Ar)."),
        ("coolant_flow_Lmin",           "float", "Coolant gas flow (L/min Ar)."),
        ("sampling_depth_mm",           "float", "Sampling depth / torch position (mm)."),
        ("per_element_notes",          "dict",  "Element-specific notes and considerations."),
        ("tips",                        "list",  "Practical optimization tips."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_elements": ["Fe", "Cu", "Zn"],
                "sample_type": "aqueous",
                "optimization_goal": "sensitivity",
            },
            "text_input": {"input_params": "Fe,Cu,Zn aqueous sensitivity"},
            "output": {
                "rf_power_W": 1383.33,
                "nebulizer_flow_Lmin": 0.92,
            },
        },
        {
            "code_input": {
                "analyte_elements": ["Na", "K"],
                "sample_type": "aqueous",
                "optimization_goal": "stability",
            },
            "text_input": {"input_params": "Na,K aqueous stability"},
            "output": {
                "rf_power_W": 1150.0,
                "nebulizer_flow_Lmin": 0.85,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load databases."""
        self.elem_db = ELEMENT_PLASMA_DB
        self.sample_db = SAMPLE_TYPE_ADJUSTMENTS

    def _run_base(
        self,
        analyte_elements: List[str],
        sample_type: str = "aqueous",
        optimization_goal: str = "sensitivity",
    ) -> Dict[str, Any]:
        """Core optimization logic."""
        elem_keys = [e.strip().lower() for e in analyte_elements]
        stype = sample_type.strip().lower()
        goal = optimization_goal.strip().lower()

        # Validate sample type
        if stype not in self.sample_db:
            available = ", ".join(sorted(self.sample_db.keys()))
            raise ChemMCPError(f"Unknown sample type '{sample_type}'. Available: {available}")

        # Get sample type base adjustments
        sample_adj = self.sample_db[stype]

        # Collect per-element data
        elem_data_list = []
        rf_values = []
        neb_values = []

        for ek in elem_keys:
            if ek not in self.elem_db:
                available = ", ".join(sorted(self.elem_db.keys()))
                raise ChemMCPError(f"Element '{ek}' not found in plasma database. Available: {available}")

            edata = self.elem_db[ek]
            elem_data_list.append({"element": ek.upper(), **edata})
            rf_values.append(edata["rf_optimal"])
            neb_values.append(edata["neb_opt"])

        # Compute averaged optimal parameters with sample-type adjustment
        avg_rf = sum(rf_values) / len(rf_values) + sample_adj.get("rf_power_offset", 0)
        avg_neb = sum(neb_values) / len(neb_values) + sample_adj.get("neb_flow_offset", 0)

        # Get aux values similarly
        aux_values = [self.elem_db[e]["aux_opt"] for e in elem_keys]
        avg_aux = sum(aux_values) / len(aux_values) + sample_adj.get("aux_flow_offset", 0)

        # Goal-based adjustments
        goal_tips = []
        if goal == "stability":
            avg_rf -= 50  # Slightly lower power for robustness
            goal_tips.append("Reduced RF power by ~50 W for improved plasma stability and longer cone life.")
        elif goal == "oxide_ratio":
            avg_rf += 100  # Higher power reduces CeO/Ce ratio
            avg_aux += 0.1  # More auxiliary gas
            goal_tips.append("Increased RF power and auxiliary flow to minimize oxide formation (CeO+/Ce+ < 1.5%).")
        elif goal == "doubly_charged":
            avg_rf += 50  # Higher power increases doubly charged ratio
            avg_neb += 0.05
            goal_tips.append("Adjusted conditions to monitor/control doubly charged ion ratio (Ba++/Ba+).")

        # Clamp values to reasonable ranges
        avg_rf = max(800, min(1600, avg_rf))
        avg_neb = max(0.50, min(1.30, avg_neb))
        avg_aux = max(0.4, min(1.5, avg_aux))

        # Build output
        per_elem = {}
        for ed in elem_data_list:
            per_elem[ed["element"]] = {
                "recommended_rf_W": ed["rf_optimal"],
                "recommended_neb_Lmin": ed["neb_opt"],
                "notes": ed.get("notes", ""),
            }

        # General tips
        tips = list(goal_tips)
        tips.append(sample_adj.get("notes", ""))
        tips.append("Always perform tuning/matching after changing plasma conditions.")
        tips.append("Use internal standards (e.g., Sc, Ge, In, Bi, Lu) covering mass range for ICP-MS.")

        # Check for conflicting requirements
        rfs = [e["rf_optimal"] for e in elem_data_list]
        if max(rfs) - min(rfs) > 300:
            tips.append(
                f"⚠️ Wide RF power range ({min(rfs)}–{max(rfs)} W) requested by different elements. "
                f"Compromise value {avg_rf:.0f} W selected. Consider splitting into two methods."
            )

        logger.info(f"Plasma optimization for {[e.upper() for e in elem_keys]}: RF={avg_rf:.0f}W, "
                     f"Neb={avg_neb:.2f} L/min, Aux={avg_aux:.2f} L/min, sample={stype}")
        return {
            "rf_power_W": round(avg_rf, 1),
            "nebulizer_flow_Lmin": round(avg_neb, 2),
            "auxiliary_flow_Lmin": round(avg_aux, 2),
            "coolant_flow_Lmin": float(sample_adj.get("coolant_flow", 15.0)),
            "sampling_depth_mm": float(sample_adj.get("sampling_depth_mm", 8.0)),
            "per_element_notes": per_elem,
            "tips": tips,
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            elems = [x.strip() for x in parts[0].split(",") if x.strip()]
            stype = parts[1] if len(parts) > 1 else "aqueous"
            goal = parts[2] if len(parts) > 2 else "sensitivity"

            return self._run_base(elems, stype, goal)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
