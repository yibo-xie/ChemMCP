import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LanthanideProperties(BaseTool):
    """
    镧系元素性质查询工具（第57-71号：La-Lu）。
    重点覆盖镧系收缩效应、氧化态（+3为主）、光谱性质、磁性、分离方法及应用。
    """
    __version__ = "0.1.0"
    name = "LanthanideProperties"
    func_name = "get_lanthanide_properties"
    description = "Query lanthanide (4f-block, elements 57-71) properties: lanthanide contraction effect, oxidation states (+3 dominant), spectral properties (f-f transitions), magnetic moments, separation methods (ion exchange), and applications."
    implementation_description = "Built-in database of all 15 lanthanides (La-Lu) plus Sc and Y (often grouped as 'rare earths'). Covers the lanthanide contraction phenomenon, ionic radii trends, f-orbital shielding effects, characteristic sharp-line spectra, magnetic data, and industrial/technological uses."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Lanthanides", "Rare Earth", "f-Block", "Lanthanide Contraction", "Magnetism", "Spectroscopy"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'Eu', 'europium', or 'all' for all)."),
        ("property_type", "str", "all", "'physical', 'oxidation', 'spectral', 'magnetic', 'applications', 'contraction', 'trends', or 'all'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element [property_type]'. Example: 'Eu oxidation' or 'all contraction'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested data."),
    ]

    examples = [
        {
            "code_input": {"element": "Eu", "property_type": "oxidation"},
            "text_input": {"input_params": "Eu oxidation"},
            "output": {"result": {"element": "Europium", "common_ox_states": [+2, +3] }}
        },
    ]

    # Lanthanide data: symbol → properties
    DATABASE = {
        "La": {"name": "Lanthanum", "Z": 57, "atomic_mass": 138.91, "config": "[Xe] 5d¹ 6s²",
               "common_ox": [+3], "notable_ox": [], "M3_radius_pm": 103.2,
               "color_aq_La3p": "Colorless", "ground_term": "^1D₀",
               "notes": "No 4f electrons; often grouped with lanthanides despite being d¹"},
        "Ce": {"name": "Cerium", "Z": 58, "atomic_mass": 140.12, "config": "[Xe] 4f¹ 5d¹ 6s²",
               "common_ox": [+3, +4], "notable_ox": "+4 stable (unique! — 4f⁰ is stable config)",
               "M3_radius_pm": 101, "color_aq_Ce3p": "Colorless", "color_solid_Ce4p": "Yellow-orange",
               "ground_term": "^2F_{5/2}",
               "notes": "Only lanthanide with accessible +4 in aqueous chemistry; Ce(IV) is strong oxidizer; used in self-cleaning ovens (catalytic oxidizer); mischmetal component"},
        "Pr": {"name": "Praseodymium", "Z": 59, "atomic_mass": 140.91, "config": "[Xe] 4f³ 6s²",
               "common_ox": [+3, +4], "notable_ox": "+4 only in solids (PrO₂)",
               "M3_radius_pm": 99, "color_aq_Pr3p": "Yellow-green", "ground_term": "^4I_{9/2}",
               "notes": "Pr(III) salts yellow-green; Pr₆O₁₁ mixed oxide; used in high-strength alloys, magnets (PrFeB)"},
        "Nd": {"name": "Neodymium", "Z": 60, "atomic_mass": 144.24, "config": "[Xe] 4f⁴ 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 98.3, "color_aq_Nd3p": "Purple/violet-red", "ground_term": "^5I_{4}",
               "notes": "MOST IMPORTANT commercial lanthanide; Nd₂Fe₁₄B magnets (strongest permanent magnets); Nd:YAG lasers (1064 nm); purple glass"},
        "Pm": {"name": "Promethium", "Z": 61, "atomic_mass": 145.0, "config": "[Xe] 4f⁵ 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 97, "color_aq_Pm3p": "Unknown (radioactive)", "ground_term": "^6H_{5/2}",
               "notes": "ONLY radioactive lanthanide (no stable isotopes); longest: ^145Pm t½=17.7y; found in trace amounts in uranium ores; used in nuclear batteries (beta source → light via phosphor)"},
        "Sm": {"name": "Samarium", "Z": 62, "atomic_mass": 150.36, "config": "[Xe] 4f⁶ 6s²",
               "common_ox": [+2, +3], "notable_ox": "+2 relatively stable (half-filled f-shell approaching)",
               "M3_radius_pm": 95.8, "color_aq_Sm3p": "Pale yellow", "ground_term": "^7F₀",
               "notes": "SmCo₅ permanent magnets (high-temperature performance); Sm(II) can exist (unusual); neutron absorber (control rods)"},
        "Eu": {"name": "Europium", "Z": 63, "atomic_mass": 151.96, "config": "[Xe] 4f⁷ 6s²",
               "common_ox": [+2, +3], "notable_ox": "+2 very stable (half-filled 4f⁷ shell!)",
               "M3_radius_pm": 94.7, "color_aq_Eu3p": "Pale pink", "color_aq_Eu2p": "Colorless",
               "ground_term": "^7F₀",
               "notes": "MOST notable +2 state among lanthanides (4f⁷ half-fill stability); Eu²⁺ similar to Sr²⁺/Ba²⁺ (similar ionic radius); red/blue phosphors in LEDs/screens; euro banknote anti-counterfeiting ink"},
        "Gd": {"name": "Gadolinium", "Z": 64, "atomic_mass": 157.25, "config": "[Xe] 4f⁷ 5d¹ 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 93.8, "color_aq_Gd3p": "Colorless", "ground_term": "^8S_{7/2}",
               "notes": "Half-filled 4f⁷ shell → S-state (no orbital contribution to μeff); highest thermal neutron capture cross-section of any element (49,000 barns) — used in nuclear reactor control rods and MRI contrast agents (Gd-DTPA); Curie temperature 292K (ferromagnetic below RT!)"},
        "Tb": {"name": "Terbium", "Z": 65, "atomic_mass": 158.93, "config": "[Xe] 4f⁹ 6s²",
               "common_ox": [+3, +4], "notable_ox": "+4 in TbO₂/TbF₄",
               "M3_radius_pm": 92.3, "color_aq_Tb3p": "Colorless/pale pink", "ground_term": "^7F₆",
               "notes": "Tb(IV) exists (4f⁸ → one hole from half-fill); green phosphor in fluorescent lights (most efficient green phosphor); magnetostrictive material (Terfenol-D)"},
        "Dy": {"name": "Dysprosium", "Z": 66, "atomic_mass": 162.50, "config": "[Xe] 4f¹⁰ 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 91.2, "color_aq_Dy3p": "Yellow", "ground_term": "^{15/2}H_{9/2}",
               "notes": "Highest magnetic moment of any element at room temp; Dy has highest magnetic susceptibility of all elements; used in NdFeB magnets (improves coercivity); Néel temperature 179K (antiferromagnetic → paramagnetic)"},
        "Ho": {"name": "Holmium", "Z": 67, "atomic_mass": 164.93, "config": "[Xe] 4f¹¹ 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 90.1, "color_aq_Ho3p": "Yellow-brown/pink", "ground_term": "^{16/2}I_{8}",
               "notes": "Has highest magnetic moment of ANY naturally occurring element (μeff ≈ 10.6 BM); Ho laser at 2.1 μm (eye-safe wavelength for LIDAR)"},
        "Er": {"name": "Erbium", "Z": 68, "atomic_mass": 167.26, "config": "[Xe] 4f¹² 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 89.0, "color_aq_Er3p": "Pink", "ground_term": "^{15/2}I_{15/2}",
               "notes": "Er:doped fiber amplifiers (EDFA — backbone of internet! optical signal amplification at 1550 nm); Er:YAG laser (2940 nm — absorbed by water, used in dentistry/dermatology); pink color in glasses"},
        "Tm": {"name": "Thulium", "Z": 69, "atomic_mass": 168.93, "config": "[Xe] 4f¹³ 6s²",
               "common_ox": [+3], "notable_ox": ["+2 (very rare)"],
               "M3_radius_pm": 88.0, "color_aq_Tm3p": "Green", "ground_term": "^3H₆",
               "notes": "Rarest lanthanide after promethium; portable X-ray sources (^170Tm emits γ-rays); Tm-doped solid-state lasers"},
        "Yb": {"name": "Ytterbium", "Z": 70, "atomic_mass": 173.05, "config": "[Xe] 4f¹⁴ 6s²",
               "common_ox": [+2, +3], "notable_ox": "+2 quite stable (filled 4f¹⁴ shell)",
               "M3_radius_pm": 86.8, "color_aq_Yb3p": "Colorless", "color_aq_Yb2p": "Colorless",
               "ground_term": "^1S₀",
               "notes": "Yb(II) stable (filled f-shell like Eu(II)); Yb:YAG lasers (1030 nm); pressure ionization studies; used in metallurgy as deoxidizing agent"},
        "Lu": {"name": "Lutetium", "Z": 71, "atomic_mass": 174.97, "config": "[Xe] 4f¹⁴ 5d¹ 6s²",
               "common_ox": [+3], "notable_ox": [],
               "M3_radius_pm": 86.1, "color_aq_Lu3p": "Colorless", "ground_term": "^1D₂",
               "notes": "Last lanthanide; filled 4f subshell + 5d¹; heaviest and hardest lanthanide; actually a d-block element by electron config but grouped with lanthanides chemically; used in PET scanners (^176Lu) and as catalyst in petroleum refining"},
    }

    CONTRACTION_DATA = {
        "definition": "The **lanthanide contraction** is the greater-than-expected decrease in ionic (and atomic) radii across the lanthanide series (La³⁺ 103.2 pm → Lu³⁺ 86.1 pm, total decrease ~17 pm across 14 elements).",
        "cause": "Poor shielding of nuclear charge by 4f electrons. Although each successive element adds one proton (+1 nuclear charge) and one 4f electron, the 4f orbital is very diffuse and does not shield effectively. The effective nuclear charge (Z_eff) increases steadily, pulling all electrons (including outer 5s/5p/6s) closer.",
        "consequences": [
            "Mo/Hf and Tc/Ta have nearly identical radii (hard to separate — always found together in minerals)",
            "Zr (160 pm) and Hf (159 pm) are almost identical in size — extremely difficult chemical separation",
            "Second/third row transition metals are very similar in size to their second-row congeners",
            "Au is more noble than Ag (relativistic + contraction effects combined)",
            "Heavy post-transition metals (Tl, Pb, Bi) show 'inert pair effect' (6s electrons reluctant to participate in bonding)",
            "All Ln³⁺ ions have similar sizes → similar chemistry → hard to separate (requires ion exchange/solvent extraction)",
        ],
        "radii_data": {
            "La3+": 103.2, "Ce3+": 101.0, "Pr3+": 99.0, "Nd3+": 98.3, "Pm3+": 97.0,
            "Sm3+": 95.8, "Eu3+": 94.7, "Gd3+": 93.8, "Tb3+": 92.3, "Dy3+": 91.2,
            "Ho3+": 90.1, "Er3+": 89.0, "Tm3+": 88.0, "Yb3+": 86.8, "Lu3+": 86.1,
        }
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "all") -> dict:
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()

        if prop_type == "contraction":
            return {"result": self.CONTRACTION_DATA}

        if prop_type == "trends":
            return {"result": self._get_trends()}

        if element == "All":
            result = {}
            for sym in self.DATABASE:
                result[sym] = {k: v for k, v in self.DATABASE[sym].items() if k != "Z"}
            return {"result": result}

        if element not in self.DATABASE:
            raise ChemMCPError(f"Element '{element}' not found. Options: {list(self.DATABASE.keys()) + ['All']}")

        data = self.DATABASE[element]
        return {"result": {**{"element": element}, **data}}

    def _get_trends(self) -> dict:
        return {
            "oxidation_states": {
                "dominant": "+3 for ALL lanthanides (the hallmark of lanthanide chemistry)",
                "exceptions": {
                    "+4": "Ce(IV) most common/stable (4f⁰); also Pr(IV), Tb(IV) in oxides/fluorides",
                    "+2": "Eu(II), Yb(II) most stable (f⁷, f¹⁴ half-filled/filled shells); Sm(II), Tm(II) less so",
                },
                "stability_pattern": "+4 favored when it gives empty/half-filled/filled 4f shell; +2 similarly",
            },
            "ionic_radii": "Steadily decrease La³⁺(103.2) → Lu³⁺(86.1) pm (lanthanide contraction)",
            "basicity_of_oxides_hydroxides": "Decreases left→right (larger ions = more basic): La(OH)₃ > Lu(OH)₃",
            "solubility": "Generally similar; slight decrease in solubility of hydroxycarboxylates across series",
            "complex_stability": "Increases slightly across series (smaller ions form stronger complexes — chelate effect)",
            "spectral_characteristics": {
                "type": "f-f transitions (Laporte-forbidden → weak but VERY SHARP lines)",
                "why_sharp": "4f orbitals are inner (shielded by 5s²5p⁶) → barely affected by ligand field → narrow energy levels → sharp absorption/emission lines",
                "contrast_with_d_block": "d-d transitions are broad (d orbitals exposed to ligands); f-f transitions are atom-like (sharp)",
                "application": "Sharp emission lines make lanthanides ideal for lasers (Nd:YAG), phosphors (Eu red, Tb green), and fiber optics (Er 1550nm)",
            },
            "magnetic_properties": {
                "origin": "Unpaired 4f electrons contribute spin AND orbital angular momentum (unlike 3d where ligands quench orbital contribution)",
                "formula": "μeff = g_J √[J(J+1)] μB (not spin-only!)",
                "exceptions": "La³⁺(f⁰) and Lu³⁺(f¹⁴) diamagnetic; Gd³⁺(f⁷, ^8S_{7/2}) has only spin contribution (L=0)",
                "temperature_dependence": "Most follow Curie-Weiss law closely",
            },
            "separation_methods": [
                "Ion-exchange chromatography (using cation exchange resin + complexing eluent like citrate/ammonium α-HIBA — different formation constants give differential retention)",
                "Solvent extraction (organophosphorus extractants like TBP, HDEHP — slightly different distribution coefficients)",
                "Fractional crystallization (historical — exploits small solubility differences)",
                "Redox separation (for Ce/Eu which have accessible +4/+2 states)",
            ],
            "abundance_in_earths_crust_ppm": "Ce(66) > Nd(38) > La(32) > Y(31) > Sm(7) > Gd(6) > Pr(6.5) > Dy(4.5) > Er(3) > Yb(2.7) > Ho(1.3) > Tb(1.1) > Tm(0.45) > Lu(0.8) > Eu(2) > Pm(trace)",
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            elem = parts[0] if parts else "All"
            prop = parts[1] if len(parts) > 1 else "all"
            return self._run_base(elem, prop)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'element [property_type]'")
