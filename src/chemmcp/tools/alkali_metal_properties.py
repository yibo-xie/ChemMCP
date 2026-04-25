import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AlkaliMetalProperties(BaseTool):
    """
    碱金属（第1族）典型性质与反应查询工具。
    覆盖 Li, Na, K, Rb, Cs, Fr 的物理/化学性质、特征反应、焰色反应等。
    """
    __version__ = "0.1.0"
    name = "AlkaliMetalProperties"
    func_name = "get_alkali_metal_properties"
    description = "Query typical properties and reactions of alkali metals (Group 1: Li, Na, K, Rb, Cs, Fr), including physical properties, chemical reactivity, flame tests, and characteristic reactions."
    implementation_description = "Uses a built-in database of alkali metal properties including atomic data, physical constants, chemical reactivity trends (increasing down the group), common reactions with water/oxygen/halogens, flame colors, and biological roles."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Alkali Metals", "Group 1", "Periodic Trends", "Element Properties", "Reactivity"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'Na', 'sodium', 'potassium', or 'all' for all)."),
        ("property_type", "str", "all", "Property category: 'physical', 'chemical', 'reactions', 'flame_test', 'trends', 'biological', or 'all'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'element [property_type]'. Example: 'Na reactions' or 'all trends'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested property data for the specified element(s)."),
    ]

    examples = [
        {
            "code_input": {"element": "Na", "property_type": "reactions"},
            "text_input": {"input_params": "Na reactions"},
            "output": {
                "result": {
                    "element": "Na",
                    "name": "Sodium",
                    "reactions": {
                        "with_water": "2Na + 2H2O → 2NaOH + H2↑ (vigorous, melts into sphere)",
                        "with_oxygen": "4Na + O2 → 2Na2O (main); 2Na + O2 → Na2O2 (heated)",
                        "with_chlorine": "2Na + Cl2 → 2NaCl (burns with yellow-orange flame)",
                        "with_acid": "2Na + 2HCl → 2NaCl + H2↑ (explosive)",
                    },
                }
            }
        },
        {
            "code_input": {"element": "all", "property_type": "flame_test"},
            "text_input": {"input_params": "all flame_test"},
            "output": {
                "result": {
                    "flame_colors": {
                        "Li": "Crimson red",
                        "Na": "Intense yellow (persistent)",
                        "K": "Lilac (through cobalt glass)",
                        "Rb": "Red-violet",
                        "Cs": "Blue-violet",
                        "Fr": "Unknown (predicted: similar to Cs)",
                    }
                }
            }
        },
    ]

    # Comprehensive alkali metal database
    DATABASE = {
        "Li": {
            "name": "Lithium", "number": 3, "atomic_mass": 6.94,
            "electron_config": "[He] 2s¹", "metallic_radius_pm": 152,
            "density_g_cm3": 0.534, "melting_point_c": 180.5, "boiling_point_c": 1342,
            "first_ionization_kjmol": 520.2, "electronegativity_pauling": 0.98,
            "E_standard_V": -3.04,
            "flame_color": "Crimson red",
            "reactions": {
                "water": "2Li + 2H₂O → 2LiOH + H₂↑ (steady, less vigorous than Na)",
                "oxygen": "4Li + O₂ → 2LiO₂ (mainly forms Li₂O on burning; also Li₂O₂)",
                "nitrogen": "6Li + N₂ → 2Li₃N (directly combines with N₂, unique among alkali metals)",
                "hydrogen": "2Li + H₂ → 2LiH (forms stable hydride)",
                "halogen_general": "2Li + X₂ → 2LiX (X = F, Cl, Br, I)",
            },
            "characteristics": [
                "Lightest solid element (density 0.534 g/cm³)",
                "Hardest alkali metal",
                "Highest melting/boiling point of Group 1",
                "Only alkali metal that reacts directly with nitrogen",
                "Smallest ionic radius (Li⁺: 76 pm) → high charge density → polarizing power",
                "Forms covalent character in some compounds (anomalous behavior)",
                "Used in lithium-ion batteries, psychiatric medication (Li₂CO₃), greases",
            ],
            "biological_role": "Essential trace element; mood stabilization (Li₂CO₃ for bipolar disorder)",
            "oxide_type": "Li₂O (normal oxide)",
        },
        "Na": {
            "name": "Sodium", "number": 11, "atomic_mass": 22.99,
            "electron_config": "[Ne] 3s¹", "metallic_radius_pm": 186,
            "density_g_cm3": 0.968, "melting_point_c": 97.8, "boiling_point_c": 883,
            "first_ionization_kjmol": 495.8, "electronegativity_pauling": 0.93,
            "E_standard_V": -2.71,
            "flame_color": "Intense yellow (persistent, masks other colors)",
            "reactions": {
                "water": "2Na + 2H₂O → 2NaOH + H₂↑ (vigorous, melts into shiny sphere, may ignite H₂)",
                "oxygen": "4Na + O₂ → 2Na₂O (room temp); 2Na + O₂ → Na₂O₂ (heated/burning, yellow)",
                "halogen_general": "2Na + X₂ → 2NaX (vigorous, often explosive with F₂)",
                "acid": "2Na + 2HCl → 2NaCl + H₂↑ (explosive)",
                "ammonia": "2Na + 2NH₃(l) → 2NaNH₂ + H₂↑ (deep blue solution, solvated electrons)",
            },
            "characteristics": [
                "Most abundant alkali metal in Earth's crust (2.36%)",
                "Soft enough to cut with knife, stored under oil",
                "Essential for life (Na⁺ electrolyte balance, nerve function)",
                "Sodium lamp: 589 nm D-line (same as flame color)",
                "Reacts so vigorously with water it can ignite hydrogen gas",
            ],
            "biological_role": "Major extracellular cation; nerve impulse transmission, osmotic balance",
            "oxide_type": "Na₂O / Na₂O₂ (peroxide when heated)",
        },
        "K": {
            "name": "Potassium", "number": 19, "atomic_mass": 39.10,
            "electron_config": "[Ar] 4s¹", "metallic_radius_pm": 227,
            "density_g_cm3": 0.862, "melting_point_c": 63.5, "boiling_point_c": 759,
            "first_ionization_kjmol": 418.8, "electronegativity_pauling": 0.82,
            "E_standard_V": -2.93,
            "flame_color": "Lilac (viewed through cobalt glass to filter Na contamination)",
            "reactions": {
                "water": "2K + 2H₂O → 2KOH + H₂↑ (very violent, H₂ usually ignites immediately)",
                "oxygen": "K + O₂ → KO₂ (forms superoxide predominantly!)",
                "halogen_general": "2K + X₂ → 2KX (more violent than Na)",
                "acid": "2K + 2HCl → 2KCl + H₂↑ (extremely violent)",
            },
            "characteristics": [
                "More reactive than Na; lower density than water (floats!)",
                "Forms superoxide KO₂ as main product with O₂ (unlike Na which needs heating)",
                "Essential nutrient for plants (K⁺) and animals",
                "Found in minerals: sylvite (KCl), carnallite (KMgCl₃·6H₂O)",
                "Important fertilizer component (N-P-K)",
            ],
            "biological_role": "Major intracellular cation; enzyme activation, osmotic regulation, heart function",
            "oxide_type": "KO₂ (superoxide, predominant product)",
        },
        "Rb": {
            "name": "Rubidium", "number": 37, "atomic_mass": 85.47,
            "electron_config": "[Kr] 5s¹", "metallic_radius_pm": 248,
            "density_g_cm3": 1.532, "melting_point_c": 39.3, "boiling_point_c": 688,
            "first_ionization_kjmol": 403.0, "electronegativity_pauling": 0.82,
            "E_standard_V": -2.98,
            "flame_color": "Red-violet",
            "reactions": {
                "water": "2Rb + 2H₂O → 2RbOH + H₂↑ (explosive)",
                "oxygen": "Rb + O₂ → RbO₂ (superoxide)",
                "air": "Spontaneously ignites in air (pyrophoric)",
            },
            "characteristics": [
                "More reactive than K; ignites spontaneously in air",
                "Soft, waxy metal, golden color",
                "Denser than water (sinks, unlike Na/K)",
                "Used in atomic clocks (⁸⁷Rb frequency standard)",
                "Very rare (~90 ppm in crust)",
            ],
            "biological_role": "Can partially substitute for K in biological systems (not essential)",
            "oxide_type": "RbO₂ (superoxide)",
        },
        "Cs": {
            "name": "Cesium", "number": 55, "atomic_mass": 132.91,
            "electron_config": "[Xe] 6s¹", "metallic_radius_pm": 265,
            "density_g_cm3": 1.873, "melting_point_c": 28.5, "boiling_point_c": 671,
            "first_ionization_kjmol": 375.7, "electronegativity_pauling": 0.79,
            "E_standard_V": -3.03,
            "flame_color": "Blue-violet",
            "reactions": {
                "water": "2Cs + 2H₂O → 2CsOH + H₂↑ (explosive at -116°C!)",
                "oxygen": "Cs + O₂ → CsO₂ (superoxide)",
                "ice": "Explodes on contact with ice (even at very low temperatures)",
                "halogens": "2Cs + F₂ → 2CsF (spontaneous combustion in cold F₂ gas)",
            },
            "characteristics": [
                "Most reactive naturally occurring metal (excluding Fr)",
                "Lowest melting point of all metals (28.5°C, near room temperature!)",
                "Largest metallic radius of non-radioactive elements",
                "Liquid just above body temperature",
                "Used in atomic clocks (most accurate, ¹³³Cs defines the second)",
                "Photoelectric cells (lowest work function of any metal)",
                "Drilling fluids (cesium formate brine)",
            ],
            "biological_role": "No known biological role; chemically similar to K but too reactive",
            "oxide_type": "CsO₂ (superoxide)",
        },
        "Fr": {
            "name": "Francium", "number": 87, "atomic_mass": 223.02,
            "electron_config": "[Rn] 7s¹", "metallic_radius_pm": 270,
            "density_g_cm3": 1.87, "melting_point_c": 27, "boiling_point_c": 677,
            "first_ionization_kjmol": 380, "electronegativity_pauling": 0.70,
            "E_standard_V": -3.1,
            "flame_color": "Unknown (predicted: similar to Cs, blue-violet)",
            "reactions": {
                "water": "Predicted: extremely violent, more explosive than Cs",
                "air": "Decays rapidly (longest isotope ²²³Fr: t½ = 22 min)",
            },
            "characteristics": [
                "Second rarest naturally occurring element (estimated <30 g in Earth's crust)",
                "All isotopes radioactive (most stable: ²²³Fr, t½ = 22.3 min)",
                "Never been seen in bulk (too little produced, decays too fast)",
                "Properties predicted from periodic trends",
                "Produced from actinium decay (alpha decay of ²²⁷Ac)",
            ],
            "biological_role": "None (radioactive, no stable isotopes)",
            "oxide_type": "Unknown (predicted: FrO₂ superoxide)",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "all") -> dict:
        """Core logic: query alkali metal properties."""
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()

        if element == "All":
            result = {}
            for sym, data in self.DATABASE.items():
                result[sym] = self._filter_properties(data, prop_type)
            return {"result": result}

        if element not in self.DATABASE:
            # Try case-insensitive search
            found = None
            for sym, data in self.DATABASE.items():
                if sym.lower() == element.lower() or data["name"].lower() == element.lower():
                    found = sym
                    break
            if not found:
                valid = ", ".join(list(self.DATABASE.keys()) + ["All"])
                raise ChemMCPError(f"Element '{element}' not found. Valid options: {valid}")
            element = found

        data = self.DATABASE[element]
        filtered = self._filter_properties(data, prop_type)

        if prop_type == "trends":
            return {"result": self._get_trends()}

        return {"result": {**{"element": element, "name": data["name"]}, **filtered}}

    def _filter_properties(self, data: dict, prop_type: str) -> dict:
        """Filter data by property type."""
        if prop_type == "all":
            return {k: v for k, v in data.items()
                    if k not in ("number", "atomic_mass", "electron_config")}
        elif prop_type == "physical":
            return {k: v for k, v in data.items()
                    if k in ("density_g_cm3", "melting_point_c", "boiling_point_c",
                              "metallic_radius_pm", "first_ionization_kjmol",
                              "electronegativity_pauling", "E_standard_v")}
        elif prop_type == "chemical":
            return {k: v for k, v in data.items()
                    if k in ("electron_config", "oxide_type", "reactions", "characteristics")}
        elif prop_type == "reactions":
            return {"reactions": data.get("reactions", {})}
        elif prop_type == "flame_test":
            return {"flame_color": data.get("flame_color"), "element": data["name"]}
        elif prop_type == "biological":
            return {"biological_role": data.get("biological_role")}
        elif prop_type == "trends":
            return self._get_trends()
        else:
            raise ChemMCPError(f"Unknown property type: {prop_type}")

    def _get_trends(self) -> dict:
        """Return Group 1 periodic trends."""
        return {
            "group_trends": {
                "atomic_radius": "Increases down group (Li: 152 pm → Cs: 265 pm) due to additional electron shells",
                "ionic_radius_Li_to_Cs": "Li⁺(76) < Na⁺(102) < K⁺(138) < Rb⁺(152) < Cs⁺(167) pm",
                "ionization_energy": "Decreases down group (Li: 520 → Cs: 376 kJ/mol) — outer electron farther from nucleus, better shielded",
                "electronegativity": "Decreases (Li: 0.98 → Cs: 0.79) — weaker attraction for bonding electrons",
                "density": "Generally increases (Li: 0.534 → Cs: 1.873 g/cm³); exception: K (0.862) < Na (0.968)",
                "melting_boiling_point": "Decreases down group — metallic bonding weakens as atomic size increases",
                "reactivity": "Increases dramatically down group — easier to lose the single valence s-electron",
                "character_of_oxides": "Basicity increases: Li₂O (less basic) → CsOH (strongest base)",
                "hydration_energy": "Decreases (Li⁺ most strongly hydrated due to small size/high charge density)",
                "flame_color_wavelength": "Increases in wavelength (red → violet) as ionization energy decreases",
            },
            "key_anomalies": [
                "Li shows diagonal relationship with Mg (similar ionic radius, charge density)",
                "Li is the only Group 1 element that forms a nitride (Li₃N) directly",
                "Li forms covalent organolithium compounds unlike other ionic alkali compounds",
                "Li salts are often anhydrous while others form hydrates (high hydration energy of Li⁺)",
                "K has lower density than Na (anomaly in density trend)",
                "Flame test: Na's intense yellow persists and can mask other colors (use cobalt blue glass)",
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
