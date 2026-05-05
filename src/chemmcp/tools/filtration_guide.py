import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Membrane filter database
FILTER_DATABASE = [
    {
        "material": "Cellulose Acetate (CA)",
        "pore_sizes_um": [0.1, 0.2, 0.22, 0.45, 0.8],
        "typical_pore": 0.45,
        "ph_range": "3-10",
        "max_temp_c": 75,
        "binding": "Low protein binding",
        "compatibility": ["aqueous", "buffer", "biological"],
        "incompatible": ["organic_solvent", "strong_acid", "strong_base"],
        "applications": ["particle_analysis", "sterile_filtration", "water_analysis", "HPLC_sample_prep"],
        "color": "white",
        "hydrophilic": True,
    },
    {
        "material": "Mixed Cellulose Ester (MCE)",
        "pore_sizes_um": [0.1, 0.22, 0.3, 0.45, 0.65, 0.8, 1.0, 1.2, 2.0, 3.0, 5.0, 8.0],
        "typical_pore": 0.45,
        "ph_range": "2-12",
        "max_temp_c": 75,
        "binding": "Moderate protein binding",
        "compatibility": ["aqueous", "buffer", "oil_base"],
        "incompatible": ["chlorinated_solvent", "strong_organic"],
        "applications": ["air_monitoring", "particle_counting", "microscopy", "gravimetric_analysis", "general_lab"],
        "color": "white",
        "hydrophilic": True,
    },
    {
        "material": "Nylon (Polyamide, NY)",
        "pore_sizes_um": [0.1, 0.2, 0.22, 0.45, 0.65, 0.8, 1.0, 2.0, 3.0],
        "typical_pore": 0.45,
        "ph_range": "3-14",
        "max_temp_c": 120,
        "binding": "High protein binding",
        "compatibility": ["aqueous", "alcohol", "hydrocarbon", "dilute_base", "most_organic"],
        "incompatible": ["concentrated_acid", "DMF", "DMSO"],
        "applications": ["HPLC_mobile_phase_filtration", "organic_solution_filtration", "sample_prep"],
        "color": "white",
        "hydrophilic": True,
    },
    {
        "material": "Polyethersulfone (PES)",
        "pore_sizes_um": [0.1, 0.2, 0.22, 0.45, 0.65],
        "typical_pore": 0.22,
        "ph_range": "1-14",
        "max_temp_c": 150,
        "binding": "Very low protein binding",
        "compatibility": ["aqueous", "buffer", "biological", "dilute_acid_base"],
        "incompatible": ["some_aromatic_solvent"],
        "applications": ["sterile_filtration", "cell_culture", "protein_solution", "pharmaceutical"],
        "color": "white/transparent",
        "hydrophilic": True,
    },
    {
        "material": "Polytetrafluoroethylene (PTFE)",
        "pore_sizes_um": [0.1, 0.2, 0.22, 0.45, 0.8, 1.0, 2.0, 3.0, 5.0, 10.0],
        "typical_pore": 0.45,
        "ph_range": "0-14",
        "max_temp_c": 260,
        "binding": "Very low binding (hydrophobic surface)",
        "compatibility": ["aggressive_chemicals", "strong_acid", "strong_base", "all_organic_solvents", "air_gas"],
        "incompatible": [],
        "applications": ["aggressive_chemical_filtration", "solvent_filtration", "air_monitoring", "corrosive_samples"],
        "color": "white/opaque",
        "hydrophilic": False,
    },
    {
        "material": "Glass Fiber (GF/F or GF/A)",
        "pore_sizes_um": [0.3, 0.7, 1.0, 1.6, 2.7],
        "typical_pore": 0.7,
        "ph_range": "0-14",
        "max_temp_c": 500,
        "binding": "Low binding",
        "compatibility": ["almost_all", "high_temp", "acid digestion"],
        "incompatible": ["HF (etches glass)"],
        "applications": ["pre-filtration", "heavy_load_samples", "air_sampling", "combustion_analysis", "TOM_analysis"],
        "color": "white",
        "hydrophilic": True,
    },
    {
        "material": "Polyvinylidene Fluoride (PVDF)",
        "pore_sizes_um": [0.1, 0.15, 0.22, 0.30, 0.45, 0.65, 0.80, 1.0, 2.0, 5.0],
        "typical_pore": 0.45,
        "ph_range": "5-10 (hydrophilic) / 1-14 (hydrophobic)",
        "max_temp_c": 135,
        "binding": "Low protein binding",
        "compatibility": ["aqueous", "alcohol", "dilute_acid_base"],
        "incompatible": ["strong_base", "PETRAethylene_glycol", "DMAC"],
        "applications": ["HPLC_prep", "sterile_filtration", "protein_binding_assays", "Western_blot"],
        "color": "white",
        "hydrophilic": True,
    },
    {
        "material": "Polycarbonate Track-Etched (PCTE)",
        "pore_sizes_um": [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 8.0, 10.0, 12.0],
        "typical_pore": 0.4,
        "ph_range": "1-14",
        "max_temp_c": 140,
        "binding": "Minimal non-specific binding",
        "compatibility": ["aqueous", "most_organic", "oil"],
        "incompatible": ["strong_oxidizing_agents"],
        "applications": ["particle_size_analysis", "microscopy", "electron_microscopy", "epa_methods", "critical_point_drying"],
        "color": "transparent/gray",
        "hydrophilic": False,
    },
]


@ChemMCPManager.register_tool
class FiltrationGuide(BaseTool):
    """
    滤膜选择指导：根据孔径、材质和目标分析物推荐合适的滤膜。
    """
    __version__ = "0.1.0"
    name = "FiltrationGuide"
    func_name = "recommend_filter"
    description = "Recommend the optimal filter membrane based on pore size requirements, sample chemistry, and target analytes."
    implementation_description = "Uses a rule-based filter database with material properties, compatibility, and application matching."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Filtration", "Sample Preparation", "Filter Selection", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("application_type", "str", "N/A", "Type of filtration application (e.g., 'HPLC_prep', 'sterile', 'particle_analysis', 'air_monitoring', 'general')."),
        ("pore_size_um", "float", "N/A", "Desired pore size in micrometers (μm). Use 0.0 if unsure."),
        ("sample_solvent", "str", "aqueous", "Sample solvent type: 'aqueous', 'organic', 'strong_acid', 'strong_base', 'mixed'."),
        ("target_analyte", "str", "general", "Target analyte class: 'protein', 'metal', 'organic_compound', 'particle', 'general'."),
        ("minimize_binding", "bool", "False", "Whether to prioritize low protein/analyte binding."),
        ("max_temperature_c", "float", "100.0", "Maximum operating temperature (°C)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'application_type pore_size_um [sample_solvent] [target_analyte] [minimize_binding] [max_temp]'"),
    ]

    output_sig = [
        ("primary_recommendation", "dict", "Best filter recommendation with details."),
        ("alternatives", "list", "Alternative filter options ranked by suitability."),
        ("selection_rationale", "str", "Explanation of selection criteria."),
        ("usage_notes", "list", "Important usage notes and precautions."),
    ]

    examples = [
        {
            "code_input": {
                "application_type": "HPLC_prep",
                "pore_size_um": 0.45,
                "sample_solvent": "organic",
            },
            "text_input": {
                "input_params": "HPLC_prep 0.45 organic",
            },
            "output": {
                "primary_recommendation": {"material": "Nylon (Polyamide, NY)"},
                "selection_rationale": "Nylon is compatible with organic solvents and widely used for HPLC mobile phase filtration.",
            },
        },
        {
            "code_input": {
                "application_type": "sterile",
                "pore_size_um": 0.22,
                "target_analyte": "protein",
                "minimize_binding": True,
            },
            "text_input": {
                "input_params": "sterile 0.22 aqueous protein True",
            },
            "output": {
                "primary_recommendation": {"material": "Polyethersulfone (PES)"},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        application_type: str,
        pore_size_um: float,
        sample_solvent: str = "aqueous",
        target_analyte: str = "general",
        minimize_binding: bool = False,
        max_temperature_c: float = 100.0,
    ) -> dict:
        """Core logic: recommend filter membrane."""
        app_key = application_type.lower().strip().replace(" ", "_")
        solvent_key = sample_solvent.lower().strip()
        analyte_key = target_analyte.lower().strip()

        # Score each filter material
        scored_filters = []
        for f in FILTER_DATABASE:
            score = self._score_filter(f, app_key, pore_size_um, solvent_key, analyte_key, minimize_binding, max_temperature_c)
            scored_filters.append((score, f))

        # Sort by score descending
        scored_filters.sort(key=lambda x: x[0], reverse=True)

        best_score, best = scored_filters[0]
        alternatives = []
        for score, f in scored_filters[1:4]:
            if score > 0.3:
                alternatives.append({"material": f["material"], "suitability_score": round(score, 2), "key_reasons": self._get_match_reasons(f, app_key, solvent_key)})

        rationale = self._build_rationale(best, app_key, pore_size_um, solvent_key, minimize_binding)
        notes = self._build_notes(best, solvent_key)

        logger.info(f"Filter recommended: {best['material']} for {application_type}, {pore_size_um}μm")
        return {
            "primary_recommendation": {
                "material": best["material"],
                "available_pore_sizes_um": best["pore_sizes_um"],
                "recommended_pore_um": best["typical_pore"] if pore_size_um == 0.0 else self._find_closest_pore(best, pore_size_um),
                "ph_range": best["ph_range"],
                "max_temperature_c": best["max_temp_c"],
                "binding_characteristics": best["binding"],
                "hydrophilic": best["hydrophilic"],
                "typical_applications": best["applications"],
                "compatibility": best["compatibility"],
                "incompatible_with": best["incompatible"],
            },
            "alternatives": alternatives,
            "selection_rationale": rationale,
            "usage_notes": notes,
        }

    def _score_filter(self, f, app, pore, solvent, analyte, min_bind, max_temp) -> float:
        score = 0.0

        # Application match (+3 max)
        if app in f["applications"]:
            score += 3.0
        elif any(app.replace("_", "") in a.replace("_", "") for a in f["applications"]):
            score += 1.5

        # Pore size availability (+2 max)
        if pore > 0:
            closest = min(f["pore_sizes_um"], key=lambda p: abs(p - pore))
            if abs(closest - pore) < 0.05:
                score += 2.0
            elif abs(closest - pore) < 0.2:
                score += 1.0
            else:
                score += 0.3
        else:
            score += 1.0  # Unknown pore size — neutral

        # Solvent compatibility (+3 max)
        if solvent in f["compatibility"]:
            score += 3.0
        elif solvent == "mixed":
            if len(f["incompatible"]) <= 1:
                score += 2.0
            else:
                score += 1.0
        elif solvent in f.get("incompatible", []):
            score -= 5.0  # Hard incompatible
        elif solvent in ("strong_acid", "strong_base"):
            try:
                ph_raw = f["ph_range"].split("/")[0].strip() if "/" in f["ph_range"] else f["ph_range"]
                ph_parts = ph_raw.split("-")
                ph_low = int(ph_parts[0].strip())
                ph_high = int(ph_parts[1].strip().split()[0])
            except (ValueError, IndexError):
                score += 0.5
                ph_low = 0
                ph_high = 14
            if solvent == "strong_acid" and ph_low <= 1:
                score += 2.0
            elif solvent == "strong_base" and ph_high >= 13:
                score += 2.0
            else:
                score -= 1.0

        # Binding priority (+1 max)
        if min_bind:
            if "low" in f["binding"].lower() or "very low" in f["binding"].lower():
                score += 1.5
            elif "minimal" in f["binding"].lower():
                score += 2.0
            elif "high" in f["binding"].lower():
                score -= 1.0

        # Temperature check
        if f["max_temp_c"] >= max_temp:
            score += 0.5
        else:
            score -= 2.0

        # Analyte-specific
        if analyte == "protein" and "low" in f["binding"].lower():
            score += 1.0
        if analyte == "particle" and "track" in f["material"].lower():
            score += 1.0

        return max(0.0, score)

    def _find_closest_pore(self, f, target):
        if not f["pore_sizes_um"]:
            return f["typical_pore"]
        return min(f["pore_sizes_um"], key=lambda p: abs(p - target))

    def _get_match_reasons(self, f, app, solvent):
        reasons = []
        if app in f["applications"]:
            reasons.append(f"Suitable for {app}")
        if solvent in f["compatibility"]:
            reasons.append(f"Compatible with {solvent}")
        return reasons

    def _build_rationale(self, f, app, pore, solvent, min_bind):
        parts = [f"{f['material']} selected because:"]
        if app in f["applications"]:
            parts.append(f"It is designed for {app.replace('_', ' ')} applications.")
        if pore > 0:
            closest = self._find_closest_pore(f, pore)
            parts.append(f"Available pore size {closest} μm matches your requirement ({pore} μm).")
        if solvent in f["compatibility"]:
            parts.append(f"Compatible with {solvent} samples.")
        if min_bind and ("low" in f["binding"].lower() or "minimal" in f["binding"].lower()):
            parts.append(f"{f['binding']} — suitable when minimizing analyte loss is important.")
        return " ".join(parts)

    def _build_notes(self, f, solvent):
        notes = [
            f"pH range: {f['ph_range']}",
            f"Max temperature: {f['max_temp_c']}°C",
            f"Binding characteristics: {f['binding']}",
        ]
        if f["incompatible"]:
            notes.append(f"Avoid contact with: {', '.join(f['incompatible'])}")
        if not f["hydrophilic"] and solvent == "aqueous":
            notes.append("This membrane is hydrophobic — pre-wetting with alcohol may be needed for aqueous solutions.")
        if f["material"].startswith("PTFE") and solvent == "aqueous":
            notes.append("For PTFE with aqueous samples, consider using hydrophilic-treated PTFE or pre-wet with methanol.")
        return notes

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            app = parts[0]
            pore = float(parts[1])
            solvent = parts[2] if len(parts) > 2 else "aqueous"
            analyte = parts[3] if len(parts) > 3 else "general"
            bind = parts[4].lower() == "true" if len(parts) > 4 else False
            temp = float(parts[5]) if len(parts) > 5 else 100.0
            return self._run_base(app, pore, solvent, analyte, bind, temp)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'app pore [solvent] [analyte] [bind_bool] [max_temp]'")
