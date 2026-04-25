import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AlkalineEarthProperties(BaseTool):
    """
    碱土金属（第2族）性质查询工具。
    覆盖 Be, Mg, Ca, Sr, Ba, Ra 的物理/化学性质、特征反应、化合物类型等。
    """
    __version__ = "0.1.0"
    name = "AlkalineEarthProperties"
    func_name = "get_alkaline_earth_properties"
    description = "Query properties of alkaline earth metals (Group 2: Be, Mg, Ca, Sr, Ba, Ra), including physical data, chemical reactivity trends, characteristic reactions, and compound types."
    implementation_description = "Built-in database of Group 2 element properties covering atomic structure, physical constants, chemical behavior (increasing reactivity down group), hard/soft acid character, and applications."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Alkaline Earth", "Group 2", "Periodic Trends", "Element Properties"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'Ca', 'calcium', or 'all' for all)."),
        ("property_type", "str", "all", "Property category: 'physical', 'chemical', 'reactions', 'compounds', 'trends', 'applications', or 'all'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element [property_type]'. Example: 'Ca reactions' or 'all trends'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested property data."),
    ]

    examples = [
        {
            "code_input": {"element": "Ca", "property_type": "reactions"},
            "text_input": {"input_params": "Ca reactions"},
            "output": {
                "result": {
                    "element": "Ca",
                    "name": "Calcium",
                    "reactions": {
                        "water": "Ca + 2H₂O → Ca(OH)₂ + H₂↑ (moderate, warm water)",
                        "oxygen": "2Ca + O₂ → 2CaO (burns with brick-red flame)",
                        "halogen": "Ca + F₂ → CaF₂ (vigorous); Ca + Cl₂ → CaCl₂",
                        "acid": "Ca + 2HCl → CaCl₂ + H₂↑",
                        "nitrogen": "3Ca + N₂ → Ca₃N₂ (heated)",
                    }
                }
            }
        },
    ]

    DATABASE = {
        "Be": {
            "name": "Beryllium", "number": 4, "atomic_mass": 9.012,
            "electron_config": "[He] 2s²", "metallic_radius_pm": 112,
            "density_g_cm3": 1.85, "melting_point_c": 1287, "boiling_point_c": 2469,
            "first_ionization_kjmol": 899.5, "second_ionization_kjmol": 1757.1,
            "electronegativity_pauling": 1.57, "E_standard_V": -1.85,
            "ionic_radius_pm": 45,
            "flame_color": "No color (white)",
            "characteristics": [
                "Lightest alkaline earth metal; very hard and brittle",
                "Amphoteric oxide (unique in G2 — others have basic oxides)",
                "Covalent character in compounds (high charge density, Fajans' rules)",
                "Does not react with water; protective BeO layer",
                "High charge/size ratio → strong polarizing power (hard acid)",
                "Toxic! Berylliosis from inhalation of dust/fumes",
                "Used in X-ray windows (Be transparent to X-rays), aerospace alloys (Cu-Be), nuclear reactors (moderator)",
            ],
            "reactions": {
                "water": "No reaction at room temp (BeO protection)",
                "acid": "Be + 2HCl → BeCl₂ + H₂↑ (slow)",
                "base": "Be + 2NaOH + 2H₂O → Na₂[Be(OH)₄] + H₂↑ (amphoteric! unique in G2)",
                "oxygen": "2Be + O₂ → 2BeO (protective layer)",
            },
            "anomaly": "Shows diagonal relationship with Al (amphoteric, covalent compounds)",
        },
        "Mg": {
            "name": "Magnesium", "number": 12, "atomic_mass": 24.305,
            "electron_config": "[Ne] 3s²", "metallic_radius_pm": 160,
            "density_g_cm3": 1.738, "melting_point_c": 650, "boiling_point_c": 1090,
            "first_ionization_kjmol": 737.7, "second_ionization_kjmol": 1450.7,
            "electronegativity_pauling": 1.31, "E_standard_V": -2.37,
            "ionic_radius_pm": 72,
            "flame_color": "Brilliant white",
            "characteristics": [
                "Light structural metal (alloys: duralumin, magnalium, electron)",
                "Reacts slowly with cold water, readily with steam/hot water",
                "Burns with intense white light (photography flashbulbs historically)",
                "Essential for life (chlorophyll contains Mg²⁺ at center of porphyrin ring)",
                "Grignard reagents (RMgX) — cornerstone of organic synthesis",
                "Reducing agent in Ti/K production (Kroll process, Pidgeon process)",
            ],
            "reactions": {
                "water": "Mg + 2H₂O(steam) → Mg(OH)₂ + H₂↑ (cold water very slow)",
                "oxygen": "2Mg + O₂ → 2MgO (burns with brilliant white light)",
                "nitrogen": "3Mg + N₂ → Mg₃N₂ (burning in air)",
                "carbon_dioxide": "2Mg + CO₂ → 2MgO + C (burns in CO₂! important fire safety note)",
                "acid": "Mg + 2HCl → MgCl₂ + H₂↑ (vigorous)",
            },
            "applications": ["Structural alloys (aircraft, cars)", "Fireworks (white sparks)", "Grignard reagents", "Medicine (antacids, laxatives)", "Desulfurization in steel"],
        },
        "Ca": {
            "name": "Calcium", "number": 20, "atomic_mass": 40.078,
            "electron_config": "[Ar] 4s²", "metallic_radius_pm": 197,
            "density_g_cm3": 1.55, "melting_point_c": 842, "boiling_point_c": 1484,
            "first_ionization_kjmol": 589.8, "second_ionization_kjmol": 1145.4,
            "electronegativity_pauling": 1.00, "E_standard_V": -2.87,
            "ionic_radius_pm": 100,
            "flame_color": "Brick-red / orange-red",
            "characteristics": [
                "5th most abundant element in Earth's crust (~3.4%)",
                "Essential for living organisms (bones, teeth, cell signaling, muscle function)",
                "Reacts steadily with water (less vigorously than Na but more than Mg)",
                "Softer than Mg; can be cut with knife",
                "Never found free in nature (always as compounds: CaCO₃ limestone, CaSO₄ gypsum, CaF₂ fluorite)",
            ],
            "reactions": {
                "water": "Ca + 2H₂O → Ca(OH)₂ + H₂↑ (moderate rate, warm water faster)",
                "oxygen": "2Ca + O₂ → 2CaO (burns with brick-red flame)",
                "halogen": "Ca + F₂ → CaF₂ (very vigorous); Ca + Cl₂ → CaCl₂",
                "nitrogen": "3Ca + N₂ → Ca₃N₂ (at high temperature)",
                "acid": "Ca + 2HCl → CaCl₂ + H₂↑ (vigorous)",
            },
            "applications": ["Cement/concrete (CaO from CaCO₃)", "Metallurgy (deoxidizer)", "Biology (bones: hydroxyapatite Ca₁₀(PO₄)₆(OH)₂)", "Food supplements", "Reducing agent"],
        },
        "Sr": {
            "name": "Strontium", "number": 38, "atomic_mass": 87.62,
            "electron_config": "[Kr] 5s²", "metallic_radius_pm": 215,
            "density_g_cm3": 2.64, "melting_point_c": 777, "boiling_point_c": 1382,
            "first_ionization_kjmol": 549.5, "second_ionization_kjmol": 1064.3,
            "electronegativity_pauling": 0.95, "E_standard_V": -2.89,
            "ionic_radius_pm": 118,
            "flame_color": "Crimson red",
            "characteristics": [
                "Softer than Ca; more reactive with water",
                "Similar to calcium chemically (can substitute in biological systems — ⁹⁰Sr is dangerous bone-seeking radioisotope)",
                "Strontianite (SrCO₃) and celestite (SrSO₄) are main ores",
                "Compounds give crimson red flame (fireworks)",
            ],
            "reactions": {
                "water": "Sr + 2H₂O → Sr(OH)₂ + H₂↑ (more vigorous than Ca)",
                "oxygen": "2Sr + O₂ → 2SrO (burns with crimson flame)",
                "halogen": "Sr + F₂ → SrF₂ (extremely vigorous)",
            },
            "applications": ["Fireworks (crimson red: SrCO₃, Sr(NO₃)₂)", "CRT displays (historical)", "Toothpaste for sensitive teeth (SrCl₂)"],
        },
        "Ba": {
            "name": "Barium", "number": 56, "atomic_mass": 137.33,
            "electron_config": "[Xe] 6s²", "metallic_radius_pm": 217,
            "density_g_cm3": 3.51, "melting_point_c": 727, "boiling_point_c": 1845,
            "first_ionization_kjmol": 502.9, "second_ionization_kjmol": 965.2,
            "electronegativity_pauling": 0.89, "E_standard_V": -2.91,
            "ionic_radius_pm": 135,
            "flame_color": "Apple green",
            "characteristics": [
                "Most reactive non-radioactive alkaline earth metal",
                "Reacts vigorously with water (comparable to Group 1 reactivity)",
                "All barium compounds soluble except BaSO₄ (used in X-ray imaging — 'barium meal')",
                "Ba²⁺ is highly toxic (blocks K⁺ channels, affects muscles including heart)",
                "Stored under oil to prevent oxidation",
            ],
            "reactions": {
                "water": "Ba + 2H₂O → Ba(OH)₂ + H₂↑ (vigorous, exothermic)",
                "oxygen": "2Ba + O₂ → 2BaO (burns with apple-green flame)",
                "halogen": "Ba + F₂ → BaF₂ (spontaneous combustion)",
                "air": "Tarnishes rapidly; forms BaO and Ba₃N₂",
                "nitrogen": "3Ba + N₂ → Ba₃N₂ (readily at room temperature)",
            },
            "applications": ["X-ray contrast agent (BaSO₄ — insoluble, non-toxic)", "Fireworks (green: BaCl₂)", "Vacuum tubes (getter)", "Drilling fluids (BaSO₄ weighting agent)"],
        },
        "Ra": {
            "name": "Radium", "number": 88, "atomic_mass": 226.03,
            "electron_config": "[Rn] 7s²", "metallic_radius_pm": 223,
            "density_g_cm3": 5.5, "melting_point_c": 700, "boiling_point_c": 1140,
            "first_ionization_kjmol": 509.3, "second_ionization_kjmol": 979,
            "electronegativity_pauling": 0.90, "E_standard_V": -2.92,
            "ionic_radius_pm": 148,
            "flame_color": "Brillient red (carmines)",
            "characteristics": [
                "Only radioactive Group 2 element (all isotopes radioactive)",
                "Most stable isotope: ²²⁶Ra (t½ = 1600 years)",
                "Intensely radioactive — glows blue due to ionized air around it",
                "Decay product of ²³⁸U (present in uranium ores)",
                "Historically used in luminous paints (watch dials) — caused radiation poisoning in workers",
                "Marie Curie discovered it (1898) and named it after radius (ray)",
            ],
            "reactions": {
                "water": "Ra + 2H₂O → Ra(OH)₂ + H₂↑ (very vigorous, solution strongly basic)",
                "air": "Tarnishes rapidly, reacts with N₂ forming Ra₃N₂",
            },
            "applications": ["Historical: luminous paints, cancer treatment (radiotherapy)", "Modern: limited due to extreme radioactivity hazards; mostly replaced by safer sources"],
            "safety": "EXTREMELY DANGEROUS — alpha emitter, bone-seeker like Sr-90, causes cancer. Requires specialized handling.",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "all") -> dict:
        """Core logic."""
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()

        if element == "All":
            result = {}
            for sym in self.DATABASE:
                result[sym] = self._filter_data(self.DATABASE[sym], prop_type)
            return {"result": result}

        if element not in self.DATABASE:
            found = None
            for sym, data in self.DATABASE.items():
                if sym.lower() == element.lower() or data["name"].lower() == element.lower():
                    found = sym
                    break
            if not found:
                valid = ", ".join(list(self.DATABASE.keys()) + ["All"])
                raise ChemMCPError(f"Element '{element}' not found. Options: {valid}")
            element = found

        data = self.DATABASE[element]
        return {"result": {**{"element": element, "name": data["name"]}, **self._filter_data(data, prop_type)}}

    def _filter_data(self, data: dict, prop_type: str) -> dict:
        if prop_type == "all":
            return {k: v for k, v in data.items() if k not in ("number", "atomic_mass")}
        elif prop_type == "physical":
            keys = ("density_g_cm3", "melting_point_c", "boiling_point_c", "metallic_radius_pm",
                     "first_ionization_kjmol", "second_ionization_kjmol", "electronegativity_pauling",
                     "ionic_radius_pm", "E_standard_v")
            return {k: data.get(k) for k in keys if k in data}
        elif prop_type == "chemical":
            return {"reactions": data.get("reactions"), "characteristics": data.get("characteristics"),
                    "anomaly": data.get("anomaly")}
        elif prop_type == "reactions":
            return {"reactions": data.get("reactions", {})}
        elif prop_type == "trends":
            return self._get_trends()
        elif prop_type == "applications":
            return {"applications": data.get("applications", [])}
        else:
            return {k: v for k, v in data.items() if k != "number"}

    def _get_trends(self) -> dict:
        return {
            "group_trends": {
                "atomic_radius": "Increases down group: Be(112) < Mg(160) < Ca(197) < Sr(215) < Ba(217) pm",
                "ionic_radius_M2p": "Be²⁺(45) < Mg²⁺(72) < Ca²⁺(100) < Sr²⁺(118) < Ba²⁺(135) < Ra²⁺(148) pm",
                "ionization_energy": "Decreases down group (both IE1 and IE2) — outer electrons farther from nucleus",
                "electronegativity": "Decreases: Be(1.57) > Mg(1.31) > Ca(1.00) > Sr(0.95) > Ba(0.89) ≈ Ra(0.90)",
                "density": "Generally increases: Be(1.85) < Mg(1.74) < Ca(1.55) < Sr(2.64) < Ba(3.51) < Ra(~5.5)",
                "reactivity_with_water": "Increases dramatically: Be(none) << Mg(slow/steam) < Ca(moderate) < Sr(vigorous) < Ba(very vigorous)",
                "oxide_basicity": "Increases: BeO(amphoteric) < MgO(weakly basic) < CaO < SrO < BaO(strongly basic)",
                "solubility_of_sulfates": "DECREASES down group: BeSO₄(soluble) > MgSO₄ > CaSO₄(slightly sol.) > SrSO₄(insol.) > BaSO₄(very insol.) > RaSO₄",
                "solubility_of_hydroxides": "INCREASES down group: Be(OH)₂(amphoteric) < Mg(OH)₂(insol.) < Ca(OH)₂(slightly sol.) < Sr(OH)₂ < Ba(OH)₂(soluble)",
                "thermal_stability_of_carbonates_nitrates": "INCREASES down group (higher decomposition T needed)",
                "character_of_compounds": "More ionic down group (lower polarizing power); Be compounds most covalent",
            },
            "key_anomalies": [
                "Be shows diagonal relationship with Al (amphoteric oxide, covalent compounds, Be₂C like Al₄C₃)",
                "Be does not react with water (unlike rest of group)",
                "Be forms beryllates with base (amphoteric behavior unique in G2)",
                "Mg burns in CO₂ (important: cannot use CO₂ extinguisher on Mg fires!)",
                "Density anomaly: Ca(1.55) < Mg(1.738) — Ca has different metallic structure",
                "Ba density > Sr despite trend (close-packed structure difference)",
                "Sulfate solubility reverses the normal trend (usually increases down groups)",
            ]
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            elem = parts[0] if parts else "All"
            prop = parts[1] if len(parts) > 1 else "all"
            return self._run_base(elem, prop)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'element [property_type]'")
