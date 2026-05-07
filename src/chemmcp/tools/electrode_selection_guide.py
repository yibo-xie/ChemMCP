import logging
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Electrode database: (name, type, potential_range_V, pH_range, features, typical_uses)
ELECTRODE_DB = {
    # Working electrodes
    "glassy_carbon": {
        "type": "working",
        "potential_range_V": (-1.2, 1.5),
        "pH_range": (0, 14),
        "features": ["Chemically inert", "Wide potential window", "Good conductivity", "Low background current", "Renewable surface"],
        "typical_uses": ["Redox studies", "Organic electrochemistry", "Heavy metal detection", "HPLC-EC detection"],
        "cost": "moderate",
    },
    "platinum": {
        "type": "working",
        "potential_range_V": (-0.6, 1.3),
        "pH_range": (0, 14),
        "features": ["Excellent conductivity", "Catalytic for H2/O2 evolution", "Chemically stable", "Easy to clean"],
        "typical_uses": ["Hydrogen evolution", "Oxygen reduction", "Oxidation of organics", "Electrocatalysis"],
        "cost": "high",
    },
    "gold": {
        "type": "working",
        "potential_range_V": (-0.8, 1.4),
        "pH_range": (0, 14),
        "features": ["Surface can be modified with SAMs", "Good for thiol chemistry", "Catalytic for oxidation"],
        "typical_uses": ["Self-assembled monolayers", "Thiol-based sensors", "Cyanide detection", "Biosensors"],
        "cost": "high",
    },
    "bismuth_film": {
        "type": "working",
        "potential_range_V": (-1.4, -0.2),
        "pH_range": (2, 10),
        "features": ["Low background", "High sensitivity for metals", "Environmentally friendly alternative to Hg"],
        "typical_uses": ["Heavy metal stripping analysis (Pb, Cd, Zn)", "Anodic stripping voltammetry"],
        "cost": "low",
    },
    "screen_printed_carbon": {
        "type": "working",
        "potential_range_V": (-1.0, 1.0),
        "pH_range": (2, 12),
        "features": ["Disposable", "Portable", "Low cost", "Mass production"],
        "typical_uses": ["Point-of-care testing", "Field analysis", "Environmental monitoring", "Glucose sensing"],
        "cost": "very_low",
    },
    "boron_doped_diamond": {
        "type": "working",
        "potential_range_V": (-2.0, 2.5),
        "pH_range": (0, 14),
        "features": ["Extremely wide potential window", "Very low background current", "Chemically inert", "Long lifetime"],
        "typical_uses": ["Extreme potential applications", "Wastewater treatment", "Advanced oxidation", "Harsh media"],
        "cost": "high",
    },
    "mercury_film": {
        "type": "working",
        "potential_range_V": (-1.6, 0.1),
        "pH_range": (1, 10),
        "features": ["High hydrogen overpotential", "Excellent reproducibility", "Renewable surface"],
        "typical_uses": ["Stripping voltammetry", "Trace metal analysis", "Polarography"],
        "cost": "low",
    },
    # Reference electrodes
    "SCE": {
        "type": "reference",
        "potential_vs_SHE_V": 0.241,
        "pH_range": (0, 14),
        "features": ["Stable potential", "Common in aqueous solutions", "Contains KCl saturated solution"],
        "typical_uses": ["General aqueous electrochemistry", "CV, LSV, amperometry"],
        "notes": "Saturated Calomel Electrode; +0.241 V vs SHE at 25°C",
        "electrolyte": "KCl (sat'd)",
    },
    "Ag_AgCl_3M": {
        "type": "reference",
        "potential_vs_SHE_V": 0.210,
        "pH_range": (0, 13),
        "features": ["Compact", "Easy to use", "Stable", "No liquid junction issues (gel)"],
        "typical_uses": ["General purpose reference", "Biosensors", "Portable devices"],
        "notes": "Ag/AgCl (3M KCl); +0.210 V vs SHE at 25°C",
        "electrolyte": "KCl (3M)",
    },
    "Ag_AgCl_satd": {
        "type": "reference",
        "potential_vs_SHE_V": 0.199,
        "pH_range": (0, 13),
        "features": ["Most common reference electrode", "Inexpensive", "Robust"],
        "typical_uses": ["Standard lab work", "Teaching labs", "Routine measurements"],
        "notes": "Ag/AgCl (sat'd KCl); +0.199 V vs SHE at 25°C",
        "electrolyte": "KCl (sat'd)",
    },
    "SHE": {
        "type": "reference",
        "potential_vs_SHE_V": 0.000,
        "pH_range": (0, 14),
        "features": ["Primary standard", "Zero by definition", "Theoretical reference"],
        "typical_uses": ["Fundamental electrochemistry", "Potential reporting standard"],
        "notes": "Standard Hydrogen Electrode; not practically used as a lab electrode",
        "electrolyte": "H+ (a=1)",
    },
    "Cu_CuSO4": {
        "type": "reference",
        "potential_vs_SHE_V": 0.300,
        "pH_range": (3, 10),
        "features": ["Rugged", "Used in soil/corrosion", "Maintenance-free"],
        "typical_uses": ["Soil corrosion measurements", "Pipeline cathodic protection", "Field corrosion testing"],
        "notes": "Cu/CuSO4 (sat'd); +0.300 V vs SHE",
        "electrolyte": "CuSO4 (sat'd)",
    },
}


@ChemMCPManager.register_tool
class ElectrodeSelectionGuide(BaseTool):
    """
    Working electrode and reference electrode selection guide for electrochemical experiments.
    """
    __version__ = "0.1.0"
    name = "ElectrodeSelectionGuide"
    func_name = "electrode_selection_guide"
    description = "Guide the selection of working and reference electrodes based on experimental requirements including potential range, pH, analyte type, application scenario, and budget."
    implementation_description = "Uses a built-in database of common working electrodes (GC, Pt, Au, BDD, Bi-film, etc.) and reference electrodes (Ag/AgCl, SCE, SHE, etc.) to recommend optimal electrode combinations based on user-specified criteria."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Electrode Selection", "Working Electrode", "Reference Electrode", "Experimental Design"]
    required_envs = []

    code_input_sig = [
        ("application", "str", "N/A", "Application description or category (e.g., 'heavy metal detection', 'organic redox', 'biosensor', 'corrosion', 'catalysis')."),
        ("min_potential_V", "float", "-2.0", "Minimum required potential (V). Default: -2.0."),
        ("max_potential_V", "float", "2.5", "Maximum required potential (V). Default: 2.5."),
        ("pH", "float", "7.0", "Solution pH. Default: 7.0."),
        ("budget", "str", "any", "Budget level: 'low', 'moderate', 'high', 'any'. Default: any."),
        ("electrode_type", "str", "both", "Which electrode to select: 'working', 'reference', or 'both'. Default: both."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: 'application [min_V max_V] [pH] [budget] [type]'. Example: 'heavy metal stripping -1.4 0.1 5 low both'"),
    ]

    output_sig = [
        ("recommended_working", "list", "List of recommended working electrodes with details."),
        ("recommended_reference", "list", "List of recommended reference electrodes with details."),
        ("selection_rationale", "str", "Explanation of why these electrodes were recommended."),
        ("all_candidates_working", "list", "All working electrodes that meet basic criteria."),
        ("all_candidates_reference", "list", "All compatible reference electrodes."),
    ]

    examples = [
        {
            "code_input": {
                "application": "heavy metal stripping analysis",
                "min_potential_V": -1.4,
                "max_potential_V": 0.1,
                "pH": 5.0,
                "budget": "low",
                "electrode_type": "both",
            },
            "text_input": {"input_string": "heavy metal stripping -1.4 0.1 5 low both"},
            "output": {
                "recommended_working": [{"name": "bismuth_film", "match_reason": "Optimal for heavy metal stripping, wide negative range"}],
                "recommended_reference": [{"name": "Ag_AgCl_satd", "potential_vs_SHE_V": 0.199}],
                "selection_rationale": "For heavy metal stripping at acidic pH with low budget...",
            }
        },
        {
            "code_input": {
                "application": "organic redox study",
                "min_potential_V": -1.0,
                "max_potential_V": 1.5,
                "pH": 7.0,
                "budget": "moderate",
                "electrode_type": "working",
            },
            "text_input": {"input_string": "organic redox study -1.0 1.5 7 moderate working"},
            "output": {
                "recommended_working": [{"name": "glassy_carbon", "match_reason": "Wide potential window, inert surface ideal for organic redox"}],
                "selection_rationale": "Glassy carbon offers the best balance of potential window and cost for organic redox studies.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _score_electrode(self, name: str, info: dict, app: str, min_v: float, max_v: float, ph: float, budget: str) -> tuple:
        """Score an electrode against requirements. Returns (score, reasons)."""
        score = 0
        reasons = []

        # Check potential range compatibility
        if "potential_range_V" in info:
            e_min, e_max = info["potential_range_V"]
            if e_min <= min_v and e_max >= max_v:
                score += 30
                reasons.append("Full potential coverage")
            elif e_min <= min_v * 0.9 or e_max >= max_v * 1.1:
                score += 15
                reasons.append("Partial potential coverage")
            else:
                reasons.append(f"Limited range ({e_min}~{e_max} V)")

        # Check pH compatibility
        if "pH_range" in info:
            ph_min, ph_max = info["pH_range"]
            if ph_min <= ph <= ph_max:
                score += 20
                reasons.append("pH compatible")

        # Budget check
        if "cost" in info and budget != "any":
            cost_map = {"very_low": 1, "low": 2, "moderate": 3, "high": 4}
            req_cost = cost_map.get(budget, 99)
            elec_cost = cost_map.get(info["cost"], 3)
            if elec_cost <= req_cost:
                score += 15
                reasons.append("Within budget")
            else:
                reasons.append(f"Exceeds budget ({info['cost']})")

        # Application keyword matching
        app_lower = app.lower()
        match_score = 0
        for feature in info.get("features", []):
            for kw in app_lower.split():
                if kw.lower() in feature.lower() or feature.lower() in app_lower:
                    match_score += 5
        for use in info.get("typical_uses", []):
            for kw in app_lower.split():
                if kw.lower() in use.lower() or use.lower() in app_lower:
                    match_score += 8
        score += min(match_score, 30)
        if match_score > 0:
            reasons.append(f"Application relevance (+{min(match_score, 30)})")

        return score, reasons

    def _run_base(
        self,
        application: str,
        min_potential_V: float = -2.0,
        max_potential_V: float = 2.5,
        pH: float = 7.0,
        budget: str = "any",
        electrode_type: str = "both",
    ) -> dict:
        """Recommend electrodes based on experimental requirements."""
        rec_working = []
        rec_reference = []
        candidates_w = []
        candidates_ref = []
        rationale_parts = []

        for name, info in ELECTRODE_DB.items():
            etype = info.get("type", "")
            score, reasons = self._score_electrode(name, info, application, min_potential_V, max_potential_V, pH, budget)
            entry = {"name": name, "score": score, "reasons": reasons, **{k: v for k, v in info.items() if k != "features" and k != "typical_uses"}}
            entry["key_features"] = info.get("features", [])
            entry["typical_uses"] = info.get("typical_uses", [])

            if etype == "working":
                candidates_w.append(entry)
            elif etype == "reference":
                candidates_ref.append(entry)

        # Sort by score descending
        candidates_w.sort(key=lambda x: x["score"], reverse=True)
        candidates_ref.sort(key=lambda x: x["score"], reverse=True)

        if electrode_type in ("working", "both"):
            rec_working = candidates_w[:3]
        if electrode_type in ("reference", "both"):
            rec_reference = candidates_ref[:3]

        # Build rationale
        top_w = candidates_w[0] if candidates_w else None
        top_r = candidates_ref[0] if candidates_ref else None

        rationale = f"For '{application}' at pH {pH}, "
        if top_w:
            rationale += f"best working electrode: {top_w['name']} (score={top_w['score']}, {'; '.join(top_w['reasons'])}). "
        if top_r:
            rationale += f"Best reference: {top_r['name']} (score={top_r['score']}, {'; '.join(top_r['reasons'])})."

        return {
            "recommended_working": rec_working,
            "recommended_reference": rec_reference,
            "selection_rationale": rationale,
            "all_candidates_working": candidates_w,
            "all_candidates_reference": candidates_ref,
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input."""
        parts = input_string.strip().split()
        if not parts:
            raise ChemMCPError("Input cannot be empty.")

        application = parts[0]
        idx = 1

        min_v = float(parts[idx]) if idx < len(parts) else -2.0
        idx += 1
        max_v = float(parts[idx]) if idx < len(parts) else 2.5
        idx += 1
        ph = float(parts[idx]) if idx < len(parts) else 7.0
        idx += 1
        budget = parts[idx] if idx < len(parts) else "any"
        idx += 1
        etype = parts[idx] if idx < len(parts) else "both"

        return self._run_base(application, min_v, max_v, ph, budget, etype)
