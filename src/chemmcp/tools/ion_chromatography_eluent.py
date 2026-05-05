import logging
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common IC eluent recipes database
# Concentrations in mM unless noted otherwise
ELUENT_RECIPES = {
    # --- Isocratic Eluents ---
    "carbonate_bicarbonate": {
        "name": "Carbonate/Bicarbonate (Standard Anion)",
        "application": "General anion analysis (F⁻, Cl⁻, NO₂⁻, Br⁻, NO₃⁻, PO₄³⁻, SO₄²⁻)",
        "composition": {"Na2CO3": 4.0, "NaHCO3": 1.0},
        "ph": "~10.4",
        "flow_rate_mL_min": 1.0,
        "suppressor": "Chemical or membrane (ASRS)",
        "compatible_columns": ["IonPac AS11-HC", "IonPac AS22", "IonPac AS14"],
    },
    "hydroxide_isocratic": {
        "name": "Sodium Hydroxide (KOH Generator)",
        "application": "Full-range anion analysis with gradient capability",
        "composition": {"NaOH": 10.0},
        "ph": ">12",
        "flow_rate_mL_min": 1.0,
        "suppressor": "Anion self-regenerating suppressor (ASRS)",
        "compatible_columns": ["IonPac AS19", "IonPac AS20", "IonPac AS11-HC"],
    },
    "methanesulfonic_acid": {
        "name": "Methanesulfonic Acid (MSA)",
        "application": "General cation analysis (Li⁺, Na⁺, NH₄⁺, K⁺, Ca²⁺, Mg²⁺)",
        "composition": {"MSA": 20.0},
        "ph": "~2.5",
        "flow_rate_mL_min": 0.8,
        "suppressor": "Cation self-regenerating suppressor (CSRS)",
        "compatible_columns": ["IonPac CS12A", "IonPac CS16", "IonPac CS19"],
    },
    "nitric_acid_dilute": {
        "name": "Dilute Nitric Acid",
        "application": "Alkali/alkaline earth metal cations",
        "composition": {"HNO3": 30.0},
        "ph": "~1.5",
        "flow_rate_mL_min": 1.0,
        "suppressor": "CSRS",
        "compatible_columns": ["IonPac CS12A", "IonPac CS14A"],
    },
    # --- Gradient Eluents ---
    "hydroxide_gradient_anion": {
        "name": "NaOH Gradient (Anion)",
        "application": "Complex sample matrices with wide polarity range",
        "composition_start": {"NaOH": 1.0},
        "composition_end": {"NaOH": 60.0},
        "gradient_time_min": 25,
        "ph": ">12",
        "flow_rate_mL_min": 1.0,
        "suppressor": "ASRS-300 (4mm) or ASRS-500 (2mm)",
        "compatible_columns": ["IonPac AS19", "IonPac AS24"],
    },
    # --- Specialty Eluents ---
    "tetraborate": {
        "name": "Sodium Tetraborate (Borate)",
        "application": "Bromide/Iodide separation, oxyhalides",
        "composition": {"Na2B4O7": 10.0},
        "ph": "~9.2",
        "flow_rate_mL_min": 1.0,
        "suppressor": "ASRS",
        "compatible_columns": ["IonPac AS9-HC", "IonPac AS16"],
    },
    "p_hydroxybenzoic_acid": {
        "name": "p-Hydroxybenzoic Acid (PHBA)",
        "application": "Alkaline earth and transition metals (non-suppressed)",
        "composition": {"PHBA": 4.0, "glycine": 0.18, "H2SO4": 0.03},
        "ph": "~3.8",
        "flow_rate_mL_min": 0.9,
        "suppressor": "None (direct conductivity)",
        "compatible_columns": ["IonPac CS5A"],
    },
}

# Analyte-specific eluent recommendations
ANALYTE_ELUENT_MAP = {
    "fluoride": ("carbonate_bicarbonate", "F⁻ elutes early (~3 min), good peak shape on AS14/AS22"),
    "chloride": ("carbonate_bicarbonate", "Cl⁻ is a standard analyte in most anion methods"),
    "nitrite": ("carbonate_bicarbonate", "NO₂⁻ close to Cl⁻; use high-efficiency column for resolution"),
    "bromide": ("hydroxide_gradient_anion", "Br⁻ needs stronger eluent; NaOH gradient recommended"),
    "nitrate": ("hydroxide_gradient_anion", "NO₃⁻ strongly retained; gradient improves speed"),
    "sulfate": ("hydroxide_gradient_anion", "SO₄²⁻ highly retained; may need 30-50mM NaOH"),
    "phosphate": ("hydroxide_gradient_anion", "PO₄³⁻ divalent/trivalent; requires strong elution"),
    "lithium": ("methanesulfonic_acid", "Li⁺ elutes first in cation method; MSA isocratic"),
    "sodium": ("methanesulfonic_acid", "Na⁺ standard cation analyte"),
    "ammonium": ("methanesulfonic_acid", "NH₄⁺ well-resolved from K⁺ on CS12A/CS16"),
    "potassium": ("methanesulfonic_acid", "K⁺ standard cation analyte"),
    "calcium": ("methanesulfonic_acid", "Ca²⁺ divalent; longer retention than monovalents"),
    "magnesium": ("methanesulfonic_acid", "Mg²⁺ similar to Ca²⁺; use isocratic MSA"),
}


# Molecular weights (g/mol)
MW = {
    "Na2CO3": 105.99, "NaHCO3": 84.01, "NaOH": 40.00,
    "MSA": 96.11, "HNO3": 63.01, "Na2B4O7": 201.22,
    "PHBA": 138.12, "glycine": 75.07, "H2SO4": 98.08,
}


@ChemMCPManager.register_tool
class IonChromatographyEluent(BaseTool):
    """
    离子色谱淋洗液配制工具。
    根据目标分析物、分离模式、检测器类型，推荐淋洗液配方、浓度和配制步骤。
    """
    __version__ = "0.1.0"
    name = "IonChromatographyEluent"
    func_name = "prepare_ic_eluent"
    description = "Recommend and calculate ion chromatography eluent preparation: composition, concentration, pH, preparation steps, and compatibility."
    implementation_description = (
        "Contains a built-in recipe library of common IC eluents (carbonate/bicarbonate, hydroxide, "
        "MSA, nitric acid, borate, PHBA). Calculates exact masses/volumes needed to prepare a given "
        "volume of eluent at specified concentration. Provides step-by-step preparation instructions "
        "and column/detector compatibility guidance."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Ion Chromatography", "Eluent Preparation", "IC", "Mobile Phase", "Anion", "Cation"]
    required_envs = []

    code_input_sig = [
        ("target_analytes", "list", "N/A", "List of target analyte names (e.g., ['Cl-', 'NO3-', 'SO4(2-)'])."),
        ("eluent_type", "None", "None", "Type of eluent: 'anion', 'cation', 'auto' (auto-select based on analytes), or specific recipe name."),
        ("final_volume_L", "float", "1.0", "Final volume of eluent to prepare in liters."),
        ("concentration_factor", "float", "1.0", "Concentration multiplier (1.0 = standard concentration)."),
        ("mode", "str", "'isocratic'", "Chromatographic mode: 'isocratic' or 'gradient'."),
        ("suppress", "bool", "True", "Whether chemical suppression is used."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'analytes=Cl-,NO3-,SO4(2-) volume=1 mode=isocratic'."),
    ]

    output_sig = [
        ("eluent_recipe", "dict", "Complete eluent preparation guide including composition, masses/volumes, steps, and notes."),
    ]

    examples = [
        {
            "code_input": {
                "target_analytes": ["Cl-", "NO3-", "SO4(2-)"],
                "final_volume_L": 0.5,
            },
            "text_input": {"input_params": "analytes=Cl-,NO3-,SO4(2-) volume=0.5"},
            "output": {
                "eluent_recipe": {"recipe_name": "...", "components": [...], "preparation_steps": [...]}
            },
        },
        {
            "code_input": {
                "target_analytes": ["Na+", "K+", "Ca(2+)", "Mg(2+)"],
                "eluent_type": "cation",
                "final_volume_L": 1.0,
            },
            "text_input": {"input_params": "analytes=Na+,K+,Ca(2+),Mg(2+) type=cation"},
            "output": {
                "eluent_recipe": {"recipe_name": "Methanesulfonic Acid (MSA)", "components": [...], "preparation_steps": [...]}
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _detect_mode(self, analytes: List[str]) -> str:
        """Auto-detect if analytes are anions or cations."""
        anion_indicators = ["-", "fluoride", "chloride", "bromide", "iodide", "nitrate", "nitrite",
                           "sulfate", "phosphate", "borate", "acetate", "formate", "oxalate",
                           "carbonate", "chlorate", "perchlorate", "thiocyanate"]
        cation_indicators = ["+", "lithium", "sodium", "ammonium", "potassium", "calcium",
                            "magnesium", "barium", "strontium"]

        anion_count = sum(1 for a in analytes if any(ind in a.lower() for ind in anion_indicators))
        cation_count = sum(1 for a in analytes if any(ind in a.lower() for ind in cation_indicators))

        if cation_count > anion_count:
            return "cation"
        return "anion"

    def _select_recipe(self, analytes: List[str], eluent_type: Optional[str],
                       mode: str) -> dict:
        """Select best eluent recipe based on analytes and type."""
        if eluent_type and eluent_type in ELUENT_RECIPES:
            return ELUENT_RECIPES[eluent_type]

        detected_mode = self._detect_mode(analytes)

        if eluent_type == "anion" or (not eluent_type and detected_mode == "anion"):
            if mode == "gradient":
                return ELUENT_RECIPES["hydroxide_gradient_anion"]
            # Check if any analyte needs gradient
            needs_gradient = any(a.lower() in ["sulfate", "so4", "phosphate", "po4", "iodide", "i-"] for a in analytes)
            if needs_gradient:
                return ELUENT_RECIPES["hydroxide_gradient_anion"]
            return ELUENT_RECIPES["carbonate_bicarbonate"]

        elif eluent_type == "cation" or (not eluent_type and detected_mode == "cation"):
            return ELUENT_RECIPES["methanesulfonic_acid"]

        return ELUENT_RECIPES["carbonate_bicarbonate"]

    def _calc_mass(self, compound: str, mmol_per_L: float, volume_L: float) -> dict:
        """Calculate mass needed for given molar concentration."""
        mw = MW.get(compound)
        if mw is None:
            raise ChemMCPError(f"Molecular weight unknown for '{compound}'.")
        mol_L = mmol_per_L / 1000.0  # mM → M
        total_mol = mol_L * volume_L
        mass_g = total_mol * mw
        return {
            "compound": compound,
            "mw_g_mol": round(mw, 2),
            "concentration_mmol_L": mmol_per_L,
            "mass_grams": round(mass_g, 4),
            "mass_mg": round(mass_g * 1000, 1),
        }

    def _generate_steps(self, recipe: dict, components_calc: list, volume_L: float,
                        suppress: bool) -> List[str]:
        """Generate step-by-step preparation instructions."""
        steps = [
            f"📋 Prepare approximately {volume_L * 0.8:.1f} L of deionized water (resistivity ≥18.2 MΩ·cm).",
        ]

        has_gradient = "composition_end" in recipe

        if not has_gradient:
            for comp in components_calc:
                steps.append(
                    f"⚖️ Weigh {comp['mass_grams']} g ({comp['mass_mg']:.1f} mg) of {comp['compound']}."
                )
            steps.append(
                f"💧 Transfer weighed compounds to a {volume_L * 1000:.0f} mL volumetric flask."
            )
            steps.append("🔢 Dissolve completely in ~80% final volume of DI water.")
            steps.append(f"📏 Dilute to mark ({volume_L * 1000:.0f} mL) with DI water.")
        else:
            steps.append("This is a gradient eluent — prepare two separate solutions:")
            steps.append(f"  Solution A (weak): see component calculations below")
            steps.append(f"  Solution B (strong): see component calculations below")
            steps.append("  Program pump to blend A→B according to gradient profile.")

        if suppress:
            steps.append(
                f"✅ Install/regenerate suppressor ({recipe.get('suppressor', 'ASRS')}). "
                "Ensure regenerant flow is active."
            )

        steps.append("🌀 Degas eluent by helium sparging (5 min) or vacuum degassing.")
        steps.append("🔬 Filter through 0.45 μm (or 0.2 μm) nylon membrane filter.")
        steps.append("💾 Place in reservoir; allow temperature equilibration before use.")

        return steps

    def _get_analyte_notes(self, analytes: List[str]) -> List[str]:
        """Get analyte-specific retention/separation notes."""
        notes = []
        for a in analytes:
            key_lower = a.lower().replace("(2-)", "").replace("+", "").replace("-", "")
            if key_lower in ANALYTE_ELUENT_MAP:
                _, note = ANALYTE_ELUENT_MAP[key_lower]
                notes.append(f"  • {a}: {note}")
        return notes

    def _run_base(self, target_analytes: list,
                  eluent_type: Optional[str] = None,
                  final_volume_L: float = 1.0,
                  concentration_factor: float = 1.0,
                  mode: str = "isocratic",
                  suppress: bool = True) -> dict:
        """Core logic."""

        recipe = self._select_recipe(target_analytes, eluent_type, mode)

        # Calculate component masses
        components_calc = []
        base_comp = recipe.get("composition", {})
        for comp_name, conc_mM in base_comp.items():
            adjusted_conc = conc_mM * concentration_factor
            calc = self._calc_mass(comp_name, adjusted_conc, final_volume_L)
            components_calc.append(calc)

        # Gradient end-point calculation
        end_components = None
        if "composition_end" in recipe:
            end_components = []
            for comp_name, conc_mM in recipe["composition_end"].items():
                adjusted_conc = conc_mM * concentration_factor
                end_components.append(self._calc_mass(comp_name, adjusted_conc, final_volume_L))

        # Steps
        steps = self._generate_steps(recipe, components_calc, final_volume_L, suppress)

        # Analyte-specific notes
        analyte_notes = self._get_analyte_notes(target_analytes)

        result = {
            "eluent_recipe": {
                "recipe_name": recipe["name"],
                "application": recipe["application"],
                "mode": "gradient" if "composition_end" in recipe else "isocratic",
                "ph": recipe.get("ph", "N/A"),
                "flow_rate_mL_min": recipe.get("flow_rate_mL_min"),
                "suppression": {
                    "required": suppress,
                    "suppressor_type": recipe.get("suppressor", "Unknown"),
                },
                "recommended_columns": recipe.get("compatible_columns", []),
                "volume_prepared_L": final_volume_L,
                "concentration_multiplier": concentration_factor,
                "components": components_calc,
                "gradient_end_components": end_components,
                "preparation_steps": steps,
                "analyte_specific_notes": analyte_notes,
                "safety_notes": [
                    "Wear PPE (gloves, safety glasses) when handling acids/bases.",
                    "Always add acid/base to water slowly while stirring.",
                    "Prepare fresh daily for hydroxide eluents (absorbs CO₂).",
                    "Degas thoroughly to avoid baseline noise.",
                ],
            }
        }
        return result

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                if k == "analytes":
                    kwargs["target_analytes"] = v.split(",")
                elif k == "volume":
                    kwargs["final_volume_L"] = float(v)
                elif k == "type":
                    kwargs["eluent_type"] = v
                elif k == "mode":
                    kwargs["mode"] = v
                elif k == "factor":
                    kwargs["concentration_factor"] = float(v)
                elif k == "suppress":
                    kwargs["suppress"] = v.lower() != "false"
        return self._run_base(**kwargs)
