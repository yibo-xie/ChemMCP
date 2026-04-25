import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ActinideProperties(BaseTool):
    """
    锕系元素性质查询工具（第89-103号：Ac-Lr）。
    覆盖放射性、氧化态多样性（尤其U的+3到+6）、核化学/燃料循环、5f轨道特性。
    """
    __version__ = "0.1.0"
    name = "ActinideProperties"
    func_name = "get_actinide_properties"
    description = "Query actinide (5f-block, elements 89-103) properties: radioactivity, diverse oxidation states, nuclear chemistry, 5f orbital characteristics, key isotopes, and applications."
    implementation_description = "Built-in database of actinides covering radioisotope data, oxidation state diversity (wider than lanthanides), nuclear fuel cycle, coordination chemistry, and safety information."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Actinides", "5f-Block", "Radioactivity", "Nuclear Chemistry", "Uranium", "Plutonium"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'U', 'uranium', 'Pu', or 'all')."),
        ("property_type", "str", "all", "'oxidation', 'isotopes', 'applications', 'safety', 'trends', or 'all'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element [property_type]'. Example: 'U oxidation' or 'all trends'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested data."),
    ]

    examples = [
        {
            "code_input": {"element": "U", "property_type": "oxidation"},
            "text_input": {"input_params": "U oxidation"},
            "output": {"result": {"element": "Uranium", "common_ox_states": [+3, +4, +5, +6] }}
        },
    ]

    DATABASE = {
        "Ac": {"name": "Actinium", "Z": 89, "config": "[Rn] 6d¹ 7s²",
               "most_stable_isotope": ("^227Ac", "21.77 y", "β⁻"), "common_ox": [+3],
               "notes": "No 5f electrons; all isotopes radioactive; used in targeted alpha therapy"},
        "Th": {"name": "Thorium", "Z": 90, "config": "[Rn] 6d² 7s²",
               "most_stable_isotope": ("^232Th", "1.405×10¹⁰ y", "α → ^208Pb"),
               "common_ox": [+4], "notes": "Potential breeder fuel (^232Th+n→^233U); ThO₂ in ceramics; lung carcinogen if inhaled"},
        "Pa": {"name": "Protactinium", "Z": 91, "config": "[Rn] 5f² 6d¹ 7s²",
               "most_stable_isotope": ("^231Pa", "32,760 y", "α"),
               "common_ox": [+4, +5], "notes": "Rare; found in trace in U ores; part of ^235U decay chain"},
        "U": {"name": "Uranium", "Z": 92, "config": "[Rn] 5f³ 6d¹ 7s²",
             "most_stable_isotope": ("^238U", "4.468×10⁹ y", "α"),
             "common_ox": [+3, +4, +5, +6], "key_species": {
                 "U3+": "Reducing; green/purple solutions",
                 "U4+": "Green; UO2 is nuclear fuel form",
                 "UO2+": "Unstable, disproportionates to U4+ + UO2^2+",
                 "UO2^2+": "URANYL — linear O=U=O! Yellow-green; very stable; basis of yellowcake (U3O8)",
             },
             "notes": "^235U (0.72%) is fissile — powers reactors & weapons; enrichment via centrifugation; depleted U for armor/shielding"},
        "Np": {"name": "Neptunium", "Z": 93, "config": "[Rn] 5f⁴ 6d¹ 7s²",
              "most_stable_isotope": ("^237Np", "2.144×10⁶ y", "α"),
              "common_ox": [+3,+4,+5,+6,+7], "notes": "First transuranium element (1940); all ox states +3 to +7 accessible"},
        "Pu": {"name": "Plutonium", "Z": 94, "config": "[Rn] 5f⁶ 7s²",
              "most_stable_isotope": ("^244Pu", "8.08×10⁷ y", "α"),
              "weapon_grade": ("^239Pu", "24,110 y", "α/SF", "FISSILE — nuclear weapons material"),
              "common_ox": [+3,+4,+5,+6,+7], "colors": {"Pu(III)": "blue/violet", "Pu(IV)": "brown", "Pu(V)": "pink", "Pu(VI)": "yellow-tan"},
              "notes": "Extremely toxic (chemical + alpha radiation); pyrophoric; 6 allotropes; critical mass ~10 kg; RTG power source for spacecraft"},
        "Am": {"name": "Americium", "Z": 95, "config": "[Rn] 5f⁷ 7s²",
              "most_stable_isotope": ("^243Am", "7370 y", "α"),
              "common_ox": [+2,+3,+4,+5,+6], "dominant": "+3",
              "notes": "^241Am (t½=432y) used in SMOKE DETECTORS (ionization chamber); also industrial gauges"},
        "Cm": {"name": "Curium", "Z": 96, "config": "[Rn] 5f⁷ 6d¹ 7s²",
              "most_stable_isotope": ("^247Cm", "1.56×10⁷ y", "α"),
              "common_ox": [+3,+4], "dominant": "+3",
              "notes": "Alpha source for space missions (Cm-244); strong heat source from alpha decay; glows purplish-pink in darkness"},
        "Bk": {"name": "Berkelium", "Z": 97, "config": "[Rn] 5f⁹ 7s²",
              "most_stable_isotope": ("^247Bk", "1380 y", "α"),
              "common_ox": [+3,+4], "dominant": "+3",
              "notes": "First actinide where +3 is more stable than +4 (unlike lighter actinides)"},
        "Cf": {"name": "Californium", "Z": 98, "config": "[Rn] 5f¹⁰ 7s²",
              "most_stable_isotope": ("^251Cf", "898 y", "α"),
              "common_ox": [+2,+3,+4], "dominant": "+3",
              "notes": "^252Cf (t½=2.6y) spontaneous fission neutron source — portable for well-logging, cancer therapy (neutron brachytherapy); extremely rare/expensive"},
        "Es": {"name": "Einsteinium", "Z": 99, "config": "[Rn] 5f¹¹ 7s²",
              "most_stable_isotope": ("^252Es", "1.29 y", "α/EC/β⁺/SF"),
              "common_ox": [+2,+3], "dominant": "+3",
              "notes": "Produced only in trace amounts in high-flux nuclear reactors; first synthesized from 'Mike' thermonuclear test debris (1952)"},
        "Fm": {"name": "Fermium", "Z": 100, "config": "[Rn] 5f¹² 7s²",
              "most_stable_isotope": ("^257Fm", "100.5 d", "SF/α"),
              "common_ox": [+2,+3], "dominant": "+3",
              "notes": "No practical applications; studied for nuclear structure physics"},
        "Md": {"name": "Mendelevium", "Z": 101, "config": "[Rn] 5f¹³ 7s²",
              "most_stable_isotope": ("^258Md", "51.5 d", "α/EC"),
              "common_ox": [+2,+3], "dominant": "+3",
              "notes": "Only a few atoms ever produced at a time; no bulk properties known"},
        "No": {"name": "Nobelium", "Z": 102, "config": "[Rn] 5f¹⁴ 7s²",
              "most_stable_isotope": ("^259No", "58 min", "α/Ec/SF"),
              "common_ox": [+2,+3], "notes": "Chemistry difficult due to short half-lives; gas-phase studies suggest +2 more stable than expected"},
        "Lr": {"name": "Lawrencium", "Z": 103, "config": "[Rn] 5f¹⁴ 7s² 7p¹",
              "most_stable_isotope": ("^262Lr", "3.6 h", "α/EC/SF"),
              "common_ox": [+3], "notes": "Last actinide; may have p-electron in valence (controversial); chemistry barely characterized"},
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "all") -> dict:
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()

        if prop_type == "trends":
            return {"result": self._get_trends()}

        if element == "All":
            result = {}
            for sym in self.DATABASE:
                d = self.DATABASE[sym]
                result[sym] = {k: v for k, v in d.items() if k != "Z"}
            return {"result": result}

        if element not in self.DATABASE:
            raise ChemMCPError(f"Element '{element}' not found. Options: {list(self.DATABASE.keys()) + ['All']}")

        data = self.DATABASE[element]
        return {"result": {**{"element": element}, **data}}

    def _get_trends(self) -> dict:
        return {
            "actinide_vs_lanthanide": {
                "5f_vs_4f": "5f orbitals are MORE EXTENDED (less shielded) than 4f → participate more in bonding → wider range of oxidation states and more covalent character in actinides",
                "oxidation_state_range": "Actinides (early): up to +7 (Np, Pu) vs Lanthanides: mostly +3 (some +2/+4)",
                "covalency": "Early actinides (Th-U) show significant covalent bonding; late actinides (Am onwards) behave more like lanthanides (+3 dominant)",
                "ionic_radii": "Actinide contraction similar to lanthanide contraction but less regular (due to variable oxidation states)",
            },
            "oxidation_state_trend": {
                "early_actinides_Ac_U": "Wide range of accessible ox states (+3 to +6 for U; Np/Pu reach +7)",
                "mid_actinides_Pu_Am": "+3 becomes increasingly dominant; higher ox states become stronger oxidizers",
                "late_actinides_Cf_Lr": "+3 dominates (like lanthanides); +2 appears for some (Md, No)",
                "pattern": "Range of stable oxidation states DECREASES across the series (similar to d-block but starting from much higher maximum)",
            },
            "metallic_properties": {
                "density": "Very high (Th 11.7, U 19.1, Pu 19.8 g/cm³) — among densest metals",
                "structures": "Multiple allotropes (Pu has 6!); complex phase diagrams",
                "pyrophoricity": "Many finely divided actinides ignite spontaneously in air (especially Pu)",
            },
            "radioactivity_trends": {
                "half_lives": "Decrease dramatically after U (^238U: 4.5 By → ^244Pu: 80 My → ^251Cf: 898 y → ^262Lr: 3.6 h)",
                "decay_modes": "Alpha dominant for heavy isotopes; spontaneous fission important for Cf/Fm onwards; beta decay for neutron-deficient isotopes; electron capture common for proton-rich isotopes",
                "criticality": "Fissile isotopes (^233U, ^235U, ^239U, ^239Pu) can sustain chain reactions; minimum critical mass decreases with better reflectors/tamper",
            },
            "separation_challenges": [
                "PUREX process (Plutonium URanium EXtraction): TBP/kerosene solvent extraction — separates U and Pu from fission products and each other",
                "Actinide-lanthanide separation: Very difficult (similar ionic radii/chemistry); uses soft-donor ligands that prefer slightly more covalent actinides (e.g., BTP, CMPO, TODGA)",
                "Higher oxidation state separations: U(VI)/Pu(VI) extract into organic phase while trivalent lanthanides/fission products stay aqueous",
            ],
            "environmental_health": [
                "All actinides are RADIOACTIVE and most are chemically TOXIC (heavy metal poisoning + radiation damage)",
                "Alpha emitters (most actinides) are especially dangerous if inhaled/ingested (internal exposure damages tissue locally)",
                "Bone-seekers: Ra, Pu, Am (analogous to Ca²⁺) — accumulate in skeleton with decades-long biological half-life",
                "Kidney seekers: U (uranyl ion forms complexes in renal tubules) — causes nephrotoxicity",
                "Criticality accidents: Accumulation of fissile material (esp. ^239Pu, ^235U) in solution can accidentally achieve critical mass → fatal radiation burst",
            ]
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            elem = parts[0] if parts else "U"
            prop = parts[1] if len(parts) > 1 else "all"
            return self._run_base(elem, prop)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'element [property_type]'")
