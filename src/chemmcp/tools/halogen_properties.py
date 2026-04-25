import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HalogenProperties(BaseTool):
    """
    卤素（第17族）性质与相互反应查询工具。
    覆盖 F, Cl, Br, I, At 的物理化学性质、氧化能力递变、卤素间置换反应、
    卤化氢性质、含氧酸等。
    """
    __version__ = "0.1.0"
    name = "HalogenProperties"
    func_name = "get_halogen_properties"
    description = "Query properties of halogens (Group 17: F, Cl, Br, I, At), including physical data, oxidizing power trends, inter-halogen displacement reactions, hydrogen halide properties, oxyacids, and characteristic reactions."
    implementation_description = "Built-in database covering all halogen properties: atomic structure, physical states, electronegativity, electron affinity, redox potentials, inter-halogen reactions (displacement series), hydrogen halide properties (acidity, stability, reducing strength), oxyacid naming/strength trends, and industrial applications."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Halogens", "Group 17", "Redox", "Periodic Trends", "Hydrogen Halides", "Oxyacids"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'Cl', 'chlorine', or 'all' for all)."),
        ("property_type", "str", "all", "'physical', 'reactions', 'displacement', 'hydrogen_halide', 'oxyacid', 'trends', 'applications', or 'all'."),
        ("query_element2", "str", "", "Optional second element for inter-halogen reaction queries (e.g., 'Cl' with query_element2='Br' for displacement)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element [property_type] [element2]'. Example: 'Cl displacement Br' or 'all trends'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested property data."),
    ]

    examples = [
        {
            "code_input": {"element": "Cl", "property_type": "reactions", "query_element2": ""},
            "text_input": {"input_params": "Cl reactions"},
            "output": {
                "result": {
                    "element": "Cl", "name": "Chlorine",
                    "reactions": {
                        "with_H2": "H2 + Cl2 → 2HCl (explosive in light; chain reaction)",
                        "with_water": "Cl2 + H2O ⇌ HCl + HOCl (disproportionation)",
                        "with_NaOH_cold": "Cl2 + 2NaOH → NaCl + NaClO + H2O",
                        "with_NaOH_hot": "3Cl2 + 6NaOH → 5NaCl + NaClO3 + 3H2O",
                    }
                }
            }
        },
        {
            "code_input": {"element": "all", "property_type": "trends", "query_element2": ""},
            "text_input": {"input_params": "all trends"},
            "output": {"result": {"group_trends": {}, "displacement_series": []}},
        },
    ]

    DATABASE = {
        "F": {
            "name": "Fluorine", "number": 9, "atomic_mass": 19.00,
            "electron_config": "[He] 2s² 2p⁵", "atomic_radius_pm": 42,
            "density_g_cm3_L": 0.001696, "melting_point_c": -219.6, "boiling_point_c": -188.1,
            "electronegativity_pauling": 3.98,  # highest of all elements
            "electron_affinity_kjmol": -328.0,
            "E_standard_X2/X_V": 2.87,  # strongest oxidizing agent
            "state_stp": "Pale yellow gas",
            "color": "Pale yellow",
            "characteristics": [
                "MOST electronegative element (Pauling 3.98) — attracts electrons more strongly than any other element",
                "MOST powerful oxidizing agent among halogens (E° = +2.87 V)",
                "Smallest atomic radius among halogens — very high charge density",
                "Reacts with ALL elements except He, Ne, Ar under appropriate conditions (even noble gases! XeF₂, KrF₂)",
                "F-F bond is unusually weak (158 kJ/mol) due to lone-pair repulsion — makes F₂ highly reactive despite strong bond to other elements",
                "Forms hydrogen bonds like water (HF has very high boiling point for its MW)",
                "Extremely corrosive and toxic; attacks glass (reacts with SiO₂)",
                "No higher oxyacids (no HFO₄ equivalent) — cannot expand octet effectively",
            ],
            "reactions": {
                "with_H2": "H₂ + F₂ → 2HF (explosive even in dark/cold; chain reaction with ΔH = -537 kJ/mol)",
                "with_water": "2F₂ + 2H₂O → 4HF + O₂ (unlike other halogens — fluorine oxidizes water!)",
                "with_metal": "2Fe + 3F₂ → 2FeF₃ (forms highest oxidation state; unlike Cl which gives FeCl₂ also)",
                "with_glass": "SiO₂ + 2F₂ → SiF₄ + O₂ (attacks glass; stored in Ni/ Cu/Monel containers)",
                "with_noble_gas": "Xe + F₂ → XeF₂ (and XeF₄, XeF₆); Kr + F₂ → KrF₂",
                "disproportionation_in_base": "Does NOT disproportionate (already highest oxidation state possible for halogen)",
            },
            "hydrogen_halide": {
                "formula": "HF", "name": "Hydrogen fluoride",
                "state": "Gas (bp 19.5°C — liquid near room temp due to H-bonding!)",
                "acidity": "WEAK acid in dilute solution (strong H-bonds hold HF together); forms stable H₂F⁺ and HF₂⁻ complexes",
                "stability": "MOST stable hydrogen halide (H-F bond: 565 kJ/mol, strongest single bond known)",
                "reducing_power": "NONE — F⁻ is the weakest reducing agent among halide ions",
                "toxicity": "HIGHLY toxic and corrosive; causes severe burns and systemic fluoride poisoning",
                "special": "Etches glass (used for glass marking); forms azeotrope with water (35.35% HF, bp 120°C)",
            },
            "applications": ["Teflon/PTFE production (polymerization of C₂F₄)", "Uranium enrichment (UF₆ gaseous diffusion)", "Toothpaste (NaF for cavity prevention)", "Pharmaceuticals (fluorinated drugs often more lipophilic/bioavailable)", "Refrigerants (CFCs/HFCs — being phased out)"],
        },
        "Cl": {
            "name": "Chlorine", "number": 17, "atomic_mass": 35.45,
            "electron_config": "[Ne] 3s² 3p⁵", "atomic_radius_pm": 79,
            "density_g_cm3_L": 0.003214, "melting_point_c": -101.5, "boiling_point_c": -34.04,
            "electronegativity_pauling": 3.16,
            "electron_affinity_kjmol": -349.0,  # highest EA actually
            "E_standard_X2/X_V": 1.36,
            "state_stp": "Greenish-yellow gas",
            "color": "Greenish-yellow",
            "characteristics": [
                "Second most electronegative/halogen; very strong oxidizing agent",
                "First element isolated in pure form (Scheele 1774)",
                "Discovered as a component of marine salt (Greek chloros = green-yellow)",
                "Used in disinfection (kills bacteria by oxidizing cellular components)",
                "Supports combustion of many elements (e.g., burning Fe wool in Cl₂ gives FeCl₃)",
                "Intermediate reactivity between F₂ (extreme) and Br₂ (moderate)",
            ],
            "reactions": {
                "with_H2": "H₂ + Cl₂ → 2HCl (explosive in UV light or above 240°C; photochemical chain reaction)",
                "with_water": "Cl₂ + H₂O ⇌ HCl + HOCl (disproportionation; equilibrium lies left, K ≈ 4 × 10⁻⁴)",
                "with_cold_dilute_NaOH": "Cl₂ + 2NaOH(cold,dil) → NaCl + NaClO + H₂O (bleach formation)",
                "with_hot_conc_NaOH": "3Cl₂ + 6NaOH(hot,conc) → 5NaCl + NaClO₃ + 3H₂O (chlorate formation)",
                "with_metal": "2Fe + 3Cl₂ → 2FeCl₃ (directly to +3 state; unlike O₂ which needs conditions)",
                "displacement": "Cl₂ can displace Br⁻ and I⁻ from their salts (but not F⁻)",
            },
            "hydrogen_halide": {
                "formula": "HCl", "name": "Hydrogen chloride",
                "state": "Gas (bp -85°C); aqueous: hydrochloric acid",
                "acidity": "STRONG acid (completely dissociates in water); pKa ≈ -7",
                "stability": "Very stable (H-Cl bond: 431 kJ/mol)",
                "reducing_power_of_Cl_minus": "Weak reducing agent (Cl⁻ can be oxidized by MnO₂, KMnO₄, etc.)",
                "special": "Stomach acid (~0.5% HCl); important industrial chemical",
            },
            "applications": ["Water treatment/disinfection", "PVC plastic (polyvinyl chloride)", "Bleaching (paper, textiles)", "Production of organic chemicals (solvents, pesticides)", "Swimming pool sanitation"],
        },
        "Br": {
            "name": "Bromine", "number": 35, "atomic_mass": 79.90,
            "electron_config": "[Ar] 3d¹⁰ 4s² 4p⁵", "atomic_radius_pm": 94,
            "density_g_cm3_L": 3.1028, "melting_point_c": -7.2, "boiling_point_c": 58.8,
            "electronegativity_pauling": 2.96,
            "electron_affinity_kjmol": -324.6,
            "E_standard_X2/X_V": 1.07,
            "state_stp": "Reddish-brown liquid",
            "color": "Reddish-brown",
            "characteristics": [
                "Only nonmetallic element that is liquid at room temperature",
                "Name from Greek bromos (stench) — extremely pungent odor",
                "Volatile liquid (significant vapor pressure at RT) — fumes are highly irritating/toxic",
                "Less reactive than chlorine but still a strong oxidizing agent",
                "Found mainly in seawater (65 ppm as Br⁻) and salt brines",
            ],
            "reactions": {
                "with_H2": "H₂ + Br₂ ⇌ 2HBr (equilibrium; requires T > 200°C and Pt catalyst; less vigorous than Cl₂)",
                "with_water": "Br₂ + H₂O ⇌ HBr + HOBr (similar to Cl but K even smaller, ~7×10⁻⁹)",
                "with_NaOH": "Br₂ + 2NaOH(cold) → NaBr + NaBrO + H₂O",
                "displacement": "Br₂ can displace I⁻ from iodide salts (but not Cl⁻ or F⁻)",
            },
            "hydrogen_halide": {
                "formula": "HBr", "name": "Hydrogen bromide",
                "state": "Gas (bp -66°C); colorless fuming gas",
                "acidity": "STRONG acid (stronger than HCl because H-Br bond weaker; pKa ≈ -9)",
                "stability": "Moderately stable (H-Br bond: 366 kJ/mol)",
                "reducing_power_of_Br_minus": "Moderate reducing agent (stronger than Cl⁻, weaker than I⁻)",
            },
            "applications": ["Flame retardants (brominated compounds)", "Water treatment (sanitizing pools/spas)", "Pharmaceuticals (AgBr in photography historically)", "Agriculture (fumigants, pesticides)", "Oil drilling fluids"],
        },
        "I": {
            "name": "Iodine", "number": 53, "atomic_mass": 126.90,
            "electron_config": "[Kr] 4d¹⁰ 5s² 5p⁵", "atomic_radius_pm": 114,
            "density_g_cm3_S": 4.93, "melting_point_c": 113.7, "boiling_point_c": 184.3,
            "electronegativity_pauling": 2.66,
            "electron_affinity_kjmol": -295.2,
            "E_standard_X2/X_V": 0.54,
            "state_stp": "Shiny gray-black solid (purple vapor)",
            "color": "Purple (vapor), dark gray-black (solid)",
            "characteristics": [
                "Least reactive halogen (weakest oxidizing agent among stable halogens)",
                "Sublimes readily (solid → purple vapor without melting at atmospheric pressure)",
                "Only slightly soluble in water (0.03 g/100mL); much more soluble in KI solution (I₃⁻ complex) or organic solvents",
                "Essential trace element (thyroid hormone synthesis — thyroxine T₄ contains iodine)",
                "Deficiency causes goiter; excess causes hyperthyroidism",
                "Iodine-starch test: deep blue-black complex (qualitative test for starch AND iodine)",
            ],
            "reactions": {
                "with_H2": "H₂ + I₂ ⇌ 2HI (equilibrium strongly favors reactants; requires high T, catalyst; reversible)",
                "with_water": "Very slight disproportionation (negligible compared to Cl, Br)",
                "displacement": "I₂ CANNOT displace any other halide ion (weakest oxidizer in group)",
                "with_thiosulfate": "I₂ + 2S₂O₃²⁻ → 2I⁻ + S₄O₆²⁻ (iodometric titration basis)",
                "with_starch": "I₂ + starch → deep blue complex (analytical detection method)",
            },
            "hydrogen_halide": {
                "formula": "HI", "name": "Hydrogen iodide",
                "state": "Gas (bp -35°C); colorless",
                "acidity": "STRONG acid (strongest hydrogen halide acid; pKa ≈ -10; H-I bond weakest at 298 kJ/mol)",
                "stability": "LEAST stable HI (decomposes to H₂ + I₂ on standing; light-sensitive)",
                "reducing_power_of_I_minus": "Strongest reducing agent among halide ions (I⁻ easily oxidized by air, Fe³⁺, etc.)",
            },
            "applications": ["Antiseptic/iodine tincture (wound disinfection)", "Thyroid treatment/prevention (iodized salt, KI tablets)", "Contrast media (X-ray imaging)", "Catalyst (organic synthesis, e.g., carbonyl protection)", "LCD displays (polarizing films use iodine-doped PVA)"],
        },
        "At": {
            "name": "Astatine", "number": 85, "atomic_mass": 210.0,
            "electron_config": "[Xe] 4f¹⁴ 5d¹⁰ 6s² 6p⁵", "atomic_radius_pm": 127,
            "density_g_cm3_S": 7.0, "melting_point_c": 302, "boiling_point_c": 337,
            "electronegativity_pauling": 2.20,
            "E_standard_X2/X_V": 0.3,
            "state_stp": "Metallic solid (predicted)",
            "color": "Black/dark (predicted, probably darker than iodine)",
            "characteristics": [
                "Rarest naturally occurring element (<1 g total in Earth's crust at any time)",
                "All isotopes radioactive (longest-lived: ²¹⁰At, t½ = 8.1 hours)",
                "Properties predicted from periodic trends; never observed in macroscopic quantities",
                "May show some metallic character (closer to Po than other halogens)",
                "Produced in nuclear reactors / particle accelerators; decays too fast for practical use",
                "Name from Greek astatos (unstable) — reflects its radioactivity",
            ],
            "reactions": {
                "predicted_behavior": "Weakest halogen oxidizing agent; would be displaced by ALL other halogens",
                "expected_compounds": "At₂, HAt (should be strongest HI-type acid), metal astatides",
            },
            "applications": ["Medical research (alpha-emitter for targeted alpha therapy — ²¹¹At)", "Theoretical interest only — no practical bulk applications"],
        },
    }

    # Inter-halogen displacement series (oxidizing power)
    DISPLACEMENT_SERIES = "F₂ > Cl₂ > Br₂ > I₂"  # Each can displace those to its right

    # Displacement reactions
    DISPLACEMENT_REACTIONS = {
        ("F2", "Cl"): "F₂ + 2NaCl → 2NaF + Cl₂ (F displaces Cl⁻)",
        ("F2", "Br"): "F₂ + 2NaBr → 2NaF + Br₂",
        ("F2", "I"): "F₂ + 2NaI → 2NaF + I₂",
        ("Cl2", "Br"): "Cl₂ + 2KBr → 2KCl + Br₂ (yellow-orange solution turns brown-red)",
        ("Cl2", "I"): "Cl₂ + 2KI → 2KCl + I₂ (brown solution; starch turns blue-black if I₂ present)",
        ("Br2", "I"): "Br₂ + 2NaI → 2NaBr + I₂",
    }

    # Oxyacid trends
    OXYACID_TRENDS = {
        "naming_pattern": "perhalic(+VII) > halic(V) > halous(III) > hypohalous(I)",
        "strength_trend_acidic_H": "HClO₄ > HClO₃ > HClO₂ > HClO (same halogen: more O = stronger acid)",
        "strength_trend_same_oxidation": "HClO₄ > HBrO₄ > HIO₄ (different halogens: higher EN = stronger oxoacid)",
        "stability_trend": "Increases with oxidation number and down the group",
        "oxidizing_power_trend": "Opposite to stability trend: lower oxidation state = stronger oxidizer",
        "known_oxyacids": {
            "F": "None (no stable oxyacids; OF, OF₂ exist but are fluorides not true oxyacids)",
            "Cl": "HClO(hypochlorous), HClO₂(chlorous), HClO₃(chloric), HClO₄(perchloric)",
            "Br": "HBrO(hypobromous), HBrO₂(bromous — unstable), HBrO₃(bromic), HBrO₄(perbromic)",
            "I": "HIO(hypoiodous), HIO₃(iodic), HIO₅(periodic / metaperiodic), H₇IO₇(orthoperiodic)",
        }
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "all", query_element2: str = "") -> dict:
        """Core logic."""
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()
        elem2 = query_element2.strip().capitalize() if query_element2 else ""

        if prop_type == "displacement":
            return self._get_displacement(element, elem2)

        if prop_type == "trends":
            return {"result": self._get_all_trends()}

        if element == "All":
            result = {}
            for sym in self.DATABASE:
                result[sym] = self._filter_data(self.DATABASE[sym], prop_type)
            return {"result": result}

        if element not in self.DATABASE:
            found = None
            for sym, d in self.DATABASE.items():
                if sym.lower() == element.lower() or d["name"].lower() == element.lower():
                    found = sym
                    break
            if not found:
                raise ChemMCPError(f"Element '{element}' not found. Options: {list(self.DATABASE.keys()) + ['All']}")
            element = found

        data = self.DATABASE[element]
        return {"result": {**{"element": element, "name": data["name"]}, **self._filter_data(data, prop_type)}}

    def _filter_data(self, data: dict, prop_type: str) -> dict:
        if prop_type == "all":
            return {k: v for k, v in data.items() if k not in ("number", "atomic_mass")}
        elif prop_type == "physical":
            keys = ("density_g_cm3_L", "density_g_cm3_S", "melting_point_c", "boiling_point_c",
                     "atomic_radius_pm", "electronegativity_pauling", "electron_affinity_kjmol",
                     "E_standard_x2_x_v", "state_stp", "color")
            return {k: data.get(k) for k in keys if k in data}
        elif prop_type == "reactions":
            return {"reactions": data.get("reactions", {})}
        elif prop_type == "hydrogen_halide":
            return {"hydrogen_halide": data.get("hydrogen_halide", {})}
        elif prop_type == "applications":
            return {"applications": data.get("applications", [])}
        elif prop_type == "trends":
            return self._get_all_trends()
        else:
            return {k: v for k, v in data.items() if k not in ("number", "atomic_mass")}

    def _get_displacement(self, el1: str, el2: str) -> dict:
        """Get inter-halogen displacement info."""
        key = (el1.upper(), el2.capitalize()) if el2 else None
        if key and key in self.DISPLACEMENT_REACTIONS:
            return {"result": {
                "reaction": self.DISPLACEMENT_REACTIONS[key],
                "series": self.DISPLACEMENT_SERIES,
                "explanation": f"{el1} is a stronger oxidizing agent than {el2}, so it can displace {el2}⁻ from its salts.",
            }}
        elif el2:
            # Check reverse
            rev_key = (el2.upper(), el1.capitalize())
            if rev_key in self.DISPLACEMENT_REACTIONS:
                return {"result": {
                    "reaction": f"No displacement: {el2} CANNOT displace {el1}. Reverse reaction possible: {self.DISPLACEMENT_REACTIONS[rev_key]}",
                    "series": self.DISPLACEMENT_SERIES,
                    "explanation": f"{el2} is a WEAKER oxidizing agent than {el1}; cannot displace {el1}⁻.",
                }}
        return {"result": {
            "displacement_series": self.DISPLACEMENT_SERIES,
            "all_reactions": self.DISPLACEMENT_REACTIONS,
            "note": "Each halogen can displace halide ions of any halogen below it in the series.",
        }}

    def _get_all_trends(self) -> dict:
        return {
            "group_trends": {
                "physical_state_at_STP": "Gas(F, Cl) → Liquid(Br) → Solid(I, At) — increasing London dispersion forces",
                "color": "Pale yellow(F) → Greenish-yellow(Cl) → Reddish-brown(Br) → Purple(I) → Dark(At)",
                "atomic_radius": "Increases down group: F(42) < Cl(79) < Br(94) < I(114) < At(~127) pm",
                "electronegativity": "Decreases: F(3.98) > Cl(3.16) > Br(2.96) > I(2.66) > At(2.20)",
                "electron_affinity": "F(328) < Cl(349) > Br(325) > I(295) — Cl has highest EA (anomaly: F's small size causes repulsion)",
                "bond_energy_X-X": "F-F(158) < I-I(151) < Br(193) < Cl(242) kJ/mol — F-F anomalously weak (lone pair repulsion)",
                "oxidizing_strength_X2": "F₂ >> Cl₂ > Br₂ > I₂ >> At₂ (E° decreases down group)",
                "reducing_strength_X-minus": "F⁻ << Cl⁻ < Br⁻ < I⁻ (reverse trend — larger ions easier to oxidize)",
                "HX_bond_strength": "H-F(565) > H-Cl(431) > H-Br(366) > H-I(298) kJ/mol",
                "HX_acidity": "HF(weak) < HCl(strong) < HBr(stronger) < HI(strongest) — bond strength dominates",
                "HX_stability": "HF > HCl > HBr > HI (thermal decomposition increases down group)",
                "HX_reducing_power": "HF(none) < HCl < HBr < HI (HI strongest reducer)",
                "solubility_in_water": "All react (F₂ vigorously oxidizes water; others undergo disproportionation to varying degrees)",
            },
            "key_anomalies": [
                "F-F bond is weakest halogen bond (lone-pair repulsion in small F atom) — makes F₂ paradoxically reactive",
                "Cl has higher electron affinity than F (EA anomaly — small size of F causes electron-electron repulsion)",
                "HF is a weak acid while other HX are strong acids (strong H-bonding in HF inhibits dissociation)",
                "F has no positive oxidation states in compounds (always -1); others have +1, +3, +5, +7",
                "I₂ is soluble in nonpolar solvents (violet) but only slightly in water (gives brown color)",
                "Br is the only liquid nonmetal at STP",
            ]
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            elem = parts[0] if parts else "All"
            prop = parts[1] if len(parts) > 1 else "all"
            elem2 = parts[2] if len(parts) > 2 else ""
            return self._run_base(elem, prop, elem2)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'element [property_type] [element2]'")
