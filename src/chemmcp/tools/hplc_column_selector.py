import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ─── HPLC Stationary Phase Database ───
# Based on manufacturer data (Agilent, Waters, Thermo, Phenomenex)
_STATIONARY_PHASES = {
    # --- Reversed-Phase (RP) Columns ---
    "C18": {
        "type": "Reversed-Phase",
        "mechanism": "Hydrophobic interaction",
        "polarity": "Non-polar",
        "typical_particle_size_um": [1.5, 1.7, 1.8, 2.5, 3, 5],
        "pore_size_nm": [80, 100, 120, 130, 300],
        "ph_range": (2.0, 8.0),
        "max_temp_c": 60,
        "max_pressure_bar": [400, 600, 800, 1000, 1200],
        "efficiency": "Very high",
        "retention": "Strong for non-polar compounds",
        "best_for": ["Small molecules", "Drug substances", "Natural products", "Environmental analysis"],
        "not_recommended_for": ["Very polar compounds without ion pairing", "Large biomolecules (>10 kDa)"],
        "endcapping": "Usually yes",
        "notes": "Most widely used HPLC column; ~70% of all HPLC applications.",
    },
    "C8": {
        "type": "Reversed-Phase",
        "mechanism": "Hydrophobic interaction (weaker than C18)",
        "polarity": "Non-polar (less than C18)",
        "typical_particle_size_um": [1.7, 2.5, 3, 5],
        "pore_size_nm": [80, 100, 120, 300],
        "ph_range": (2.0, 8.0),
        "max_temp_c": 60,
        "max_pressure_bar": [400, 600],
        "efficiency": "High",
        "retention": "Moderate; shorter retention than C18",
        "best_for": ["Medium polarity compounds", "Peptides", "Compounds too retained on C18"],
        "not_recommended_for": ["Very non-polar analytes requiring strong retention"],
        "endcapping": "Yes",
        "notes": "Good compromise between retention and analysis time.",
    },
    "C4": {
        "type": "Reversed-Phase",
        "mechanism": "Hydrophobic interaction (weak)",
        "polarity": "Weakly non-polar",
        "typical_particle_size_um": [1.7, 2.5, 3.5, 5],
        "pore_size_nm": [300],
        "ph_range": (2.0, 8.0),
        "max_temp_c": 80,
        "max_pressure_bar": [400, 600],
        "efficiency": "Moderate-High",
        "retention": "Weak; fast elution of large molecules",
        "best_for": ["Proteins", "Large peptides", "Biologics"],
        "not_recommended_for": ["Small molecule separations"],
        "endcapping": "Yes",
        "notes": "Wide pore (300 nm) essential for protein separations.",
    },
    "Phenyl": {
        "type": "Reversed-Phase",
        "mechanism": "Hydrophobic + π-π interaction / dipole-induced dipole",
        "polarity": "Non-polar with π-system",
        "typical_particle_size_um": [1.7, 2.5, 3, 5],
        "pore_size_nm": [80, 100, 120],
        "ph_range": (2.0, 8.0),
        "max_temp_c": 60,
        "max_pressure_bar": [400, 600],
        "efficiency": "High",
        "retention": "Selective for aromatic/unsaturated compounds",
        "best_for": ["Aromatic compounds", "Isomers", "PAHs", "Phenol derivatives", "Chiral recognition辅助"],
        "not_recommended_for": ["Aliphatic-only compounds"],
        "endcapping": "Yes",
        "notes": "π-π interactions provide different selectivity vs C18.",
    },
    "Cyano (CN)": {
        "type": "Mixed-Mode / Normal-Phase capable",
        "mechanism": "Dipole + weak hydrophobic",
        "polarity": "Moderately polar",
        "typical_particle_size_um": [3, 5],
        "pore_size_nm": [80, 100, 120],
        "ph_range": (1.5, 9.0),
        "max_temp_c": 60,
        "max_pressure_bar": [300, 400],
        "efficiency": "Moderate",
        "retention": "Variable depending on mobile phase mode",
        "best_for": ["RP-NP mode switching", "Mildly polar compounds", "Tannins", "Steroids"],
        "not_recommended_for": ["High-efficiency UHPLC applications"],
        "endcapping": "Partial",
        "notes": "Can operate in both RP and NP modes; unique selectivity.",
    },
    "PFP (Pentafluorophenyl)": {
        "type": "Reversed-Phase (Fluorinated)",
        "mechanism": "π-π, dipole-dipole, H-bonding, steric, ion-exchange",
        "polarity": "Fluorinated non-polar",
        "typical_particle_size_um": [1.7, 2.6, 3, 5],
        "pore_size_nm": [80, 100],
        "ph_range": (1.5, 9.0),
        "max_temp_c": 70,
        "max_pressure_bar": [600, 800, 1000],
        "efficiency": "Very high",
        "retention": "Highly selective; shape-sensitive",
        "best_for": ["Isomer separation", "Halogenated compounds", "Aromatic positional isomers", "Nucleosides/nucleotides"],
        "not_recommended_for": [],
        "endcapping": "Yes (fluorinated)",
        "notes": "Excellent alternative selectivity; growing popularity for challenging isomer separations.",
    },
    "HILIC": {
        "type": "Hydrophilic Interaction Liquid Chromatography",
        "mechanism": "Partitioning into water layer on polar surface",
        "polarity": "Polar",
        "typical_particle_size_um": [1.7, 1.8, 2.5, 3, 5],
        "pore_size_nm": [100, 120],
        "ph_range": (2.0, 8.0),
        "max_temp_c": 60,
        "max_pressure_bar": [400, 600, 800],
        "efficiency": "High",
        "retention": "Strong for polar compounds (opposite of RP)",
        "best_for": ["Very polar compounds", "Metabolites", "Glycans", "Small polar drugs", "Nucleotides"],
        "not_recommended_for": ["Non-polar compounds"],
        "endcapping": "N/A (polar surface)",
        "notes": "Mobile phase: ACN-rich (≥60%) with small amount of aqueous buffer.",
    },

    # --- Size Exclusion ---
    "SEC/GPC": {
        "type": "Size Exclusion / Gel Permeation",
        "mechanism": "Hydrodynamic volume / size exclusion",
        "polarity": "Various (depends on phase)",
        "typical_particle_size_um": [3, 5, 10],
        "pore_size_nm": [50, 100, 500, 1000, "wide range"],
        "ph_range": (2.0, 12.0),
        "max_temp_c": 80,
        "max_pressure_bar": [200, 300],
        "efficiency": "Moderate",
        "retention": "Larger molecules elute first",
        "best_for": ["Polymers", "Proteins MW distribution", "Aggregation studies"],
        "not_recommended_for": ["Small molecule separation"],
        "endcapping": "N/A",
        "notes": "Selection depends on MW range of interest.",
    },

    # --- Ion Exchange ---
    "IEX (SCX/SAX)": {
        "type": "Ion Exchange Chromatography",
        "mechanism": "Electrostatic attraction to charged groups",
        "polarity": "Ionic",
        "typical_particle_size_um": [3, 5, 10],
        "pore_size_nm": [80, 100, 300],
        "ph_range": (2.0, 12.0),
        "max_temp_c": 60,
        "max_pressure_bar": [200, 400],
        "efficiency": "Moderate-High",
        "retention": "Based on charge density",
        "best_for": ["Peptides", "Proteins", "Inorganic ions", "Charged metabolites"],
        "not_recommended_for": ["Neutral compounds"],
        "endcapping": "N/A",
        "notes": "SCX=strong cation exchange; SAX=strong anion exchange.",
    },
}

# Column dimension recommendations
_COLUMN_DIMENSIONS = {
    "analytical_standard": {"length_mm": 150, "id_mm": 4.6, "application": "Routine analysis"},
    "analytical_fast": {"length_mm": 50, "id_mm": 4.6, "application": "Fast screening"},
    "uhplc_standard": {"length_mm": 100, "id_mm": 2.1, "application": "UHPLC high efficiency"},
    "uhplc_fast": {"length_mm": 50, "id_mm": 2.1, "application": "UHPLC very fast"},
    "preparative_small": {"length_mm": 250, "id_mm": 10, "application": "Semi-preparative"},
    "preparative_large": {"length_mm": 250, "id_mm": 21.2, "application": "Preparative purification"},
}

# Particle size guidance
_PARTICLE_SIZE_GUIDANCE = {
    1.5: {"type": "Sub-2μm UHPLC", "pressure_category": "very_high", "max_efficiency": "~250000 plates/m", "note": "Requires UHPLC system (>1000 bar)"},
    1.7: {"type": "Core-shell / Sub-2μm UHPLC", "pressure_category": "high", "max_efficiency": "~220000 plates/m", "note": "Waters Acquity BEH standard"},
    1.8: {"type": "Sub-2μm UHPLC", "pressure_category": "high", "max_efficiency": "~210000 plates/m", "note": "Common sub-2μm size"},
    2.5: {"type": "Small particle HPLC/UHPLC", "pressure_category": "medium-high", "max_efficiency": "~180000 plates/m", "note": "Good balance of efficiency and pressure"},
    2.6: {"type": "Core-shell (Superficially Porous)", "pressure_category": "medium", "max_efficiency": "~200000 plates/m", "note": "Kinetex/Halo style; near-UHPLC efficiency at lower pressure"},
    2.7: {"type": "Core-shell (Superficially Porous)", "pressure_category": "medium", "max_efficiency": "~190000 plates/m", "note": "Similar to 2.6 μm core-shell"},
    3.0: {"type": "Standard HPLC", "pressure_category": "medium", "max_efficiency": "~150000 plates/m", "note": "Traditional HPLC standard"},
    3.5: {"type": "Standard HPLC", "pressure_category": "medium", "max_efficiency": "~140000 plates/m", "note": "Waters XBridge common size"},
    5.0: {"type": "Conventional HPLC", "pressure_category": "low", "max_efficiency": "~100000 plates/m", "note": "Most common; compatible with most systems"},
}


@ChemMCPManager.register_tool
class HplcColumnSelector(BaseTool):
    """
    HPLC色谱柱选择工具。
    根据分析物性质、分离需求、仪器条件推荐最合适的色谱柱（固定相、粒径、柱长等）。
    """
    __version__ = "0.1.0"
    name = "HplcColumnSelector"
    func_name = "select_hplc_column"
    description = "Select the optimal HPLC column including stationary phase, particle size, column dimensions, and pore size based on analyte properties and separation requirements."
    implementation_description = (
        "Uses a rule-based expert system combining stationary phase chemistry knowledge, "
        "particle size efficiency-pressure trade-offs, and column dimension optimization rules "
        "to recommend the best HPLC column configuration."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["HPLC", "Chromatography", "Column Selection", "Method Development", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("analytes", "list", "N/A", "List of analyte descriptions or compound classes, e.g., ['small_molecule_drug', 'peptide', 'aromatic_isomers']."),
        ("molecular_weight_range", "str", "unknown", "MW range as 'min-max' Da, e.g., '100-500'."),
        ("polarity", "str", "unknown", "Overall polarity: 'non_polar', 'moderate', 'polar', 'ionic'."),
        ("separation_goal", "str", "general", "Goal: 'speed', 'resolution', 'throughput', 'preparation', 'isomer_separation'."),
        ("instrument_type", "str", "hplc", "Instrument: 'hplc' (<400 bar), 'uhplc' (>600 bar), or 'any'."),
        ("special_requirements", "list", "[]", "Special needs: ['high_ph_stable', 'high_temp', 'large_molecule', etc.]"),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Description of separation problem. Example: 'separate aromatic drug isomers MW 200-400 need high resolution UHPLC'"),
    ]

    output_sig = [
        ("recommendation", "dict", "Complete column recommendation including primary choice, alternatives, rationale, and method suggestions."),
    ]

    examples = [
        {
            "code_input": {
                "analytes": ["small_molecule_drug", "aromatic_compounds"],
                "molecular_weight_range": "200-400",
                "polarity": "moderate",
                "separation_goal": "resolution",
                "instrument_type": "uhplc",
            },
            "text_input": {"input_params": "drug-like compounds MW 200-400 moderate polarity need resolution on UHPLC"},
            "output": {
                "recommendation": {
                    "primary_choice": {"stationary_phase": "C18", "particle_size_um": 1.7, "column_length_mm": 100, "inner_diameter_mm": 2.1},
                    "rationale": "See full output",
                    "alternatives": [...],
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _score_phase(self, phase_data: dict, analytes: List[str], polarity: str,
                     mw_range_str: str, goal: str) -> float:
        """Score a stationary phase for suitability (0-100)."""
        score = 50.0

        # Polarity matching
        pol_map = {
            "non_polar": ["C18", "C8", "C4", "Phenyl", "PFP (Pentafluorophenyl)"],
            "moderate": ["C18", "C8", "Phenyl", "PFP (Pentafluorophenyl)", "Cyano (CN)"],
            "polar": ["HILIC", "Cyano (CN)", "PFP (Pentafluorophenyl)"],
            "ionic": ["IEX (SCX/SAX)", "HILIC", "PFP (Pentafluorophenyl)"],
        }
        if polarity in pol_map:
            for key in pol_map[polarity]:
                if key.lower() in phase_data.get("type", "").lower() or key == list(_STATIONARY_PHASES.keys())[list(_STATIONARY_PHASES.values()).index(phase_data)]:
                    score += 20

        # Analyte matching
        for a in analytes:
            al = a.lower()
            best_for_str = " ".join(phase_data.get("best_for", []))
            not_rec_str = " ".join(phase_data.get("not_recommended", []))
            if any(w in best_for_str for w in al.replace("_", " ").split()):
                score += 15
            if any(w in not_rec_str for w in al.replace("_", " ").split()):
                score -= 25

        # Goal adjustment
        if goal == "speed":
            if phase_data.get("efficiency") in ("Very high", "High"):
                score += 10
        elif goal == "resolution":
            if phase_data.get("efficiency") in ("Very high", "High"):
                score += 15
        elif goal == "isomer_separation":
            if "PFP" in phase_data.get("type", "") or "Phenyl" in phase_data.get("type", ""):
                score += 20

        return min(100, max(0, score))

    def _run_base(self, analytes: List[str], molecular_weight_range: str = "unknown",
                  polarity: str = "unknown", separation_goal: str = "general",
                  instrument_type: str = "hplc",
                  special_requirements: Optional[List[str]] = None) -> dict:
        """Core selection logic."""
        special_requirements = special_requirements or []

        # Parse MW range
        mw_min, mw_max = 0, float('inf')
        if molecular_weight_range != "unknown" and "-" in molecular_weight_range:
            try:
                parts = molecular_weight_range.split("-")
                mw_min, mw_max = int(parts[0]), int(parts[1])
            except ValueError:
                pass

        # Score each stationary phase
        scored = []
        for phase_name, phase_data in _STATIONARY_PHASES.items():
            sc = self._score_phase(phase_data, analytes, polarity, molecular_weight_range, separation_goal)
            scored.append((phase_name, phase_data, sc))
        scored.sort(key=lambda x: -x[2])

        # Determine optimal dimensions
        pressure_limit = {"hplc": 400, "uhplc": 1000, "any": 1000}.get(instrument_type, 400)

        # Particle size selection based on instrument and goal
        if instrument_type == "uhplc" and separation_goal in ("resolution", "isomer_separation"):
            rec_ps = 1.7
        elif instrument_type == "uhplc" and separation_goal == "speed":
            rec_ps = 1.8
        elif separation_goal == "speed":
            rec_ps = 2.7  # Core-shell for speed on HPLC
        elif separation_goal == "resolution":
            rec_ps = 1.7 if instrument_type == "uhplc" else 2.6
        else:
            rec_ps = 3.0 if instrument_type == "hplc" else 1.7

        ps_info = _PARTICLE_SIZE_GUIDANCE.get(rec_ps, _PARTICLE_SIZE_GUIDANCE[5.0])

        # Column length
        if separation_goal == "speed":
            rec_len = 50
            rec_id = 2.1 if instrument_type == "uhplc" else 4.6
        elif separation_goal in ("resolution", "isomer_separation"):
            rec_len = 100 if instrument_type == "uhplc" else 150
            rec_id = 2.1 if instrument_type == "uhplc" else 4.6
        elif separation_goal == "preparation":
            rec_len = 250
            rec_id = 21.2 if "prep" in " ".join(special_requirements).lower() else 10
        else:
            rec_len = 100 if instrument_type == "uhplc" else 150
            rec_id = 2.1 if instrument_type == "uhplc" else 4.6

        # Pore size from MW
        if mw_max > 10000:
            rec_pore = 300
        elif mw_max > 3000:
            rec_pore = 130
        else:
            rec_pore = 100

        # Build recommendation
        primary_name, primary_data, primary_score = scored[0]
        result = {
            "recommendation": {
                "primary_choice": {
                    "stationary_phase": primary_name,
                    "full_type": primary_data["type"],
                    "mechanism": primary_data["mechanism"],
                    "particle_size_um": rec_ps,
                    "particle_info": ps_info,
                    "column_length_mm": rec_len,
                    "inner_diameter_mm": rec_id,
                    "pore_size_nm": rec_pore,
                    "ph_range": primary_data["ph_range"],
                    "max_temperature_c": primary_data["max_temp_c"],
                    "suitability_score": round(primary_score, 1),
                },
                "rationale": (
                    f"Selected {primary_name} ({primary_data['type']}) as best match for {', '.join(analytes)}. "
                    f"{primary_data['mechanism']} mechanism provides appropriate selectivity. "
                    f"Particle size {rec_ps} μm ({ps_info['type']}) balances efficiency and backpressure."
                ),
                "top_alternatives": [
                    {"stationary_phase": n, "type": d["type"], "score": round(s, 1)}
                    for n, d, s in scored[1:4]
                ],
                "method_suggestions": {
                    "starting_mobile_phase": self._suggest_mobile_phase(primary_name, polarity),
                    "flow_rate_ml_min": round(self._estimate_flow_rate(rec_id, rec_ps), 1),
                    "injection_volume_ul": round(rec_id ** 2 * 0.01, 1),
                    "column_dead_volume_ul": round(3.14 * (rec_id / 2) ** 2 * rec_len / 1000, 2),
                    "estimated_backpressure_bar": self._estimate_pressure(rec_ps, rec_len, rec_id, primary_name),
                    "gradient_recommendation": self._suggest_gradient(separation_goal, polarity),
                },
                "special_considerations": [
                    primary_data.get("notes", ""),
                    f"Instrument compatibility: {'✓ Compatible' if ps_info['pressure_category'] in ('low', 'medium') or instrument_type == 'uhplc' else '⚠ Verify pressure limit'}",
                ],
                "all_phases_scored": [{"phase": n, "score": round(s, 1)} for n, _, s in scored],
            }
        }
        return result

    def _suggest_mobile_phase(self, phase: str, polarity: str) -> str:
        """Suggest starting mobile phase."""
        if phase == "HILIC":
            return "ACN/Water (70:30) with 10 mM ammonium formate pH 3.0"
        elif phase in ("IEX (SCX/SAX)"):
            return "Salt gradient in buffer (e.g., 0-500 mM NaCl in 20 mM phosphate)"
        elif phase == "Cyano (CN)" and polarity == "polar":
            return "Hexane/IPA (90:10) for NP mode"
        elif polarity in ("polar", "ionic"):
            return "Water/ACN (+0.1% formic acid), start 95:5"
        else:
            return "Water/ACN (+0.1% formic acid), start 60:40"

    def _estimate_flow_rate(self, id_mm: float, ps_um: float) -> float:
        """Estimate optimal flow rate."""
        return 0.2 * id_mm ** 2 / ps_um * 3.0

    def _estimate_pressure(self, ps_um: float, len_mm: float, id_mm: float, phase: str) -> str:
        """Rough estimate of backpressure."""
        eta = 0.89  # ~cP for water/ACN
        u = 300  # linear velocity mm/min approx
        p_est = 2500 * eta * u * len_mm / (ps_um ** 2)
        if p_est < 100:
            return f"<100 (very low)"
        elif p_est < 400:
            return f"~{int(p_est)} (standard HPLC)"
        elif p_est < 800:
            return f"~{int(p_est)} (requires UHPLC)"
        else:
            return f"~{int(p_est)} (requires high-pressure UHPLC)"

    def _suggest_gradient(self, goal: str, polarity: str) -> str:
        if goal == "speed":
            return "Steep gradient: 5% B → 95% B in 5 min"
        elif goal == "resolution":
            return "Shallow gradient: 5% B → 95% B in 20-30 min"
        elif goal == "isomer_separation":
            return "Isocratic or very shallow gradient recommended; fine-tune %B ±2%"
        else:
            return "Standard gradient: 5% B → 95% B in 10-15 min"

    def _run_text(self, input_params: str) -> dict:
        """Parse natural language input."""
        text = input_params.lower()

        # Simple keyword extraction
        analytes = []
        if any(w in text for w in ["drug", "pharmaceutical", "api"]):
            analytes.append("small_molecule_drug")
        if any(w in text for w in ["peptide", "protein", "biologic"]):
            analytes.append("peptide")
        if any(w in text for w in ["aromatic", "phenol", "benzene", "pah"]):
            analytes.append("aromatic_compounds")
        if any(w in text for w in ["isomer", "stereoisomer", "positional"]):
            analytes.append("positional_isomers")
        if any(w in text for w in ["polar", "metabolite", "glycan"]):
            analytes.append("polar_compound")

        if not analytes:
            analytes = ["general_organic"]

        polarity = "unknown"
        if "non-polar" in text or "nonpolar" in text or "hydrophob" in text:
            polarity = "non_polar"
        elif "polar" in text or "hydrophilic" in text:
            polarity = "polar"
        elif "ionic" in text or "charged" in text:
            polarity = "ionic"
        elif "moderate" in text:
            polarity = "moderate"

        goal = "general"
        if "fast" in text or "speed" in text or "quick" in text:
            goal = "speed"
        elif "resolution" in text or "separate" in text:
            goal = "resolution"
        elif "isomer" in text:
            goal = "isomer_separation"
        elif "prep" in text or "purif" in text:
            goal = "preparation"

        inst = "uhplc" if "uhplc" in text else ("hplc" if "hplc" in text else "any")

        return self._run_base(analytes, "unknown", polarity, goal, inst)
