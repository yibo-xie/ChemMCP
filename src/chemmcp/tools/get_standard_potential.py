import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetStandardPotential(BaseTool):
    """
    Query standard electrode potentials (E°) for common redox half-reactions.
    Database covers standard reduction potentials at 25°C (298 K), 1 M, 1 atm.
    """
    __version__ = "0.1.0"
    name = "GetStandardPotential"
    func_name = "get_standard_potential"
    description = "Query standard electrode potential (E°) for a given redox half-reaction or couple."
    implementation_description = "Uses a built-in database of standard reduction potentials (vs SHE, Standard Hydrogen Electrode) at 25°C. Covers ~100 common half-reactions including metals, nonmetals, and complex ions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Standard Potential", "Redox", "Electrode", "SHE"]
    required_envs = []

    code_input_sig = [
        ("query", "str", "N/A", "Half-reaction query: species name, formula, or redox couple (e.g., 'Cu2+/Cu', 'Fe3+/Fe2+', 'MnO4-/Mn2+', 'Zn2+/Zn', 'O2/H2O')."),
        ("as_oxidation", "bool", "False", "If True, return oxidation potential (E°ox = -E°red). Default returns reduction potential."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Half-reaction query string."),
    ]

    output_sig = [
        ("half_reaction", "str", "The balanced half-reaction equation with E° value."),
        ("E0_V", "float", "Standard electrode potential in Volts (vs SHE)."),
        ("description", "str", "Brief description of the redox couple."),
    ]

    examples = [
        {
            "code_input": {"query": "Cu2+/Cu", "as_oxidation": False},
            "text_input": {"query": "Cu2+/Cu"},
            "output": {
                "half_reaction": "Cu²⁺ + 2e⁻ ⇌ Cu(s)    E° = +0.337 V",
                "E0_V": 0.337,
                "description": "Copper(II)/Copper couple. Common reference electrode.",
            }
        },
        {
            "code_input": {"query": "MnO4- Mn2+", "as_oxidation": False},
            "text_input": {"query": "MnO4- Mn2+"},
            "output": {
                "half_reaction": "MnO₄⁻ + 8H⁺ + 5e⁻ ⇌ Mn²⁺ + 4H₂O    E° = +1.507 V",
                "E0_V": 1.507,
                "description": "Permanganate/Manganese(II) couple in acidic medium. Strong oxidizing agent.",
            }
        },
        {
            "code_input": {"query": "Zn2+/Zn", "as_oxidation": False},
            "text_input": {"query": "Zn2+/Zn"},
            "output": {
                "half_reaction": "Zn²⁺ + 2e⁻ ⇌ Zn(s)    E° = -0.763 V",
                "E0_V": -0.763,
                "description": "Zinc/Zinc ion couple. Strong reducing agent.",
            }
        },
        {
            "code_input": {"query": "Fe3+/Fe2+", "as_oxidation": False},
            "text_input": {"query": "Fe3+/Fe2+"},
            "output": {
                "half_reaction": "Fe³⁺ + e⁻ ⇌ Fe²⁺    E° = +0.771 V",
                "E0_V": 0.771,
                "description": "Iron(III)/Iron(II) couple. Important in biological redox systems.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build standard electrode potential database.
        All values are standard REDUCTION potentials vs SHE at 25°C.
        Source: CRC Handbook / IUPAC standard values.
        """
        self._potentials = [
            # ── Very strong oxidizing agents (E° > +1.5 V) ──
            ("F2(g) + 2e⁻ ⇌ 2F⁻", 2.87, "Fluorine/Fluoride"),
            ("H4XeO6 + 2H⁺ + 2e⁻ ⇌ XeO3 + 3H2O", 3.0, "Perxenate/Xenate"),
            ("O3 + 2H⁺ + 2e⁻ ⇌ O2 + H2O", 2.07, "Ozone/Oxygen"),
            ("S2O8^2- + 2e⁻ ⇌ 2SO4^2-", 2.01, "Peroxodisulfate/Sulfate"),
            ("Co³⁺ + e⁻ ⇌ Co²⁺", 1.92, "Cobalt(III)/Cobalt(II) (acidic)"),
            ("H2O2 + 2H⁺ + 2e⁻ ⇌ 2H2O", 1.776, "Hydrogen Peroxide/Water"),
            ("Au⁺ + e⁻ ⇌ Au(s)", 1.692, "Gold(I)/Gold"),
            ("MnO4⁻ + 4H⁺ + 3e⁻ ⇌ MnO2(s) + 2H2O", 1.70, "Permanganate/Manganese Dioxide"),
            ("Ce⁴⁺ + e⁻ ⇌ Ce³⁺", 1.61, "Cerium(IV)/Cerium(III)"),
            ("2HClO + 2H⁺ + 2e⁻ ⇌ Cl2(g) + 2H2O", 1.63, "Hypochlorous Acid/Chlorine"),
            ("MnO4⁻ + 8H⁺ + 5e⁻ ⇌ Mn²⁺ + 4H2O", 1.507, "Permanganate/Manganese(II)"),
            ("Cl2(g) + 2e⁻ ⇌ 2Cl⁻", 1.358, "Chlorine/Chloride"),
            ("Cr2O7²⁻ + 14H⁺ + 6e⁻ ⇌ 2Cr³⁺ + 7H2O", 1.33, "Dichromate/Chromium(III)"),
            ("O2(g) + 4H⁺ + 4e⁻ ⇌ 2H2O", 1.229, "Oxygen/Water (acidic)"),
            ("MnO2(s) + 4H⁺ + 2e⁻ ⇌ Mn²⁺ + 2H2O", 1.23, "Manganese Dioxide/Manganese(II)"),
            ("Br2(l) + 2e⁻ ⇌ 2Br⁻", 1.065, "Bromine/Bromide"),
            ("NO3⁻ + 4H⁺ + 3e⁻ ⇌ NO(g) + 2H2O", 0.96, "Nitrate/Nitric Oxide"),
            ("NO3⁻ + 3H⁺ + 2e⁻ ⇌ HNO2 + H2O", 0.94, "Nitrite/Nitrous Acid"),
            ("2Hg²⁺ + 2e⁻ ⇌ Hg2²⁺", 0.92, "Mercury(II)/Mercury(I)"),

            # ── Moderate oxidizing agents (+0.5 to +1.0 V) ──
            ("Ag⁺ + e⁻ ⇌ Ag(s)", 0.7996, "Silver(I)/Silver"),
            ("Fe³⁺ + e⁻ ⇌ Fe²⁺", 0.771, "Iron(III)/Iron(II)"),
            ("O2(g) + 2H⁺ + 2e⁻ ⇌ H2O2(aq)", 0.695, "Oxygen/Hydrogen Peroxide"),
            ("I2(s) + 2e⁻ ⇌ 2I⁻", 0.535, "Iodine/Iodide"),
            ("O2(g) + 2H2O + 4e⁻ ⇌ 4OH⁻", 0.401, "Oxygen/Hydroxide (basic)"),
            ("Cu²⁺ + 2e⁻ ⇌ Cu(s)", 0.337, "Copper(II)/Copper"),
            ("S4O6²⁻ + 2e⁻ ⇌ 2S2O3²⁻", 0.08, "Tetrathionate/Thiosulfate"),
            ("2H⁺ + 2e⁻ ⇌ H2(g)", 0.0000, "Standard Hydrogen Electrode (SHE) — REFERENCE"),
            ("Sn⁴⁺ + 2e⁻ ⇌ Sn²⁺", 0.15, "Tin(IV)/Tin(II)"),
            ("S(s) + 2H⁺ + 2e⁻ ⇌ H2S(g)", 0.14, "Sulfur/Hydrogen Sulfide"),
            ("AgCl(s) + e⁻ ⇌ Ag(s) + Cl⁻", 0.222, "Silver Chloride/Silver"),
            ("Hg2Cl2(s) + 2e⁻ ⇌ 2Hg(l) + 2Cl⁻", 0.268, "Calomel/Mercury (SCE ref)"),
            ("Cu²⁺ + e⁻ ⇌ Cu⁺", 0.153, "Copper(II)/Copper(I)"),

            # ── Weak oxidizing / reducing boundary (-0.5 to +0.15 V) ──
            ("AgI(s) + e⁻ ⇌ Ag(s) + I⁻", -0.152, "Silver Iodide/Silver"),
            ("Sn²⁺ + 2e⁻ ⇌ Sn(s)", -0.14, "Tin(II)/Tin"),
            ("Pb²⁺ + 2e⁻ ⇌ Pb(s)", -0.126, "Lead(II)/Lead"),
            ("Fe³⁺ + 3e⁻ ⇌ Fe(s)", -0.036, "Iron(III)/Iron"),
            ("2H⁺ + 2e⁻ ⇌ H2(g)", 0.000, "SHE Reference (defined as 0 V)"),

            # ── Reducing agents (E° < 0 V) ──
            ("Ni²⁺ + 2e⁻ ⇌ Ni(s)", -0.25, "Nickel/Nickel ion"),
            ("Co²⁺ + 2e⁻ ⇌ Co(s)", -0.28, "Cobalt/Cobalt ion"),
            ("Cd²⁺ + 2e⁻ ⇌ Cd(s)", -0.403, "Cadmium/Cadmium ion"),
            ("Fe²⁺ + 2e⁻ ⇌ Fe(s)", -0.44, "Iron(II)/Iron"),
            ("Cr³⁺ + 3e⁻ ⇌ Cr(s)", -0.74, "Chromium(III)/Chromium"),
            ("Zn²⁺ + 2e⁻ ⇌ Zn(s)", -0.763, "Zinc/Zinc ion"),
            ("2H2O + 2e⁻ ⇌ H2(g) + 2OH⁻", -0.828, "Water/Hydroxide (basic)"),
            ("Mn²⁺ + 2e⁻ ⇌ Mn(s)", -1.18, "Manganese(II)/Manganese"),
            ("Al³⁺ + 3e⁻ ⇌ Al(s)", -1.66, "Aluminum/Aluminum ion"),
            ("Mg²⁺ + 2e⁻ ⇌ Mg(s)", -2.37, "Magnesium/Magnesium ion"),
            ("Na⁺ + e⁻ ⇌ Na(s)", -2.714, "Sodium/Sodium ion"),
            ("Ca²⁺ + 2e⁻ ⇌ Ca(s)", -2.87, "Calcium/Calcium ion"),
            ("K⁺ + e⁻ ⇌ K(s)", -2.931, "Potassium/Potassium ion"),
            ("Li⁺ + e⁻ ⇌ Li(s)", -3.04, "Lithium/Lithium ion"),

            # ── Organic / biochem relevant ──
            ("NAD⁺ + 2H⁺ + 2e⁻ ⇌ NADH", -0.32, "NAD⁺/NADH (biochemical, pH 7)"),
            ("Quinone + 2H⁺ + 2e⁻ ⇌ Hydroquinone", 0.699, "Quinone/Hydroquinone"),
            ("Cytochrome c (Fe³⁺) + e⁻ ⇌ Cytochrome c (Fe²⁺)", 0.254, "Cytochrome c (pH 7)"),
            ("FAD + 2H⁺ + 2e⁻ ⇌ FADH2", -0.22, "FAD/FADH2 (approximate)"),
            ("Acetaldehyde + 2H⁺ + 2e⁻ ⇌ Ethanol", -0.197, "Acetaldehyde/Ethanol"),
            ("Pyruvate + 2H⁺ + 2e⁻ ⇌ Lactate", -0.19, "Pyruvate/Lactate"),
            ("Oxaloacetate + 2H⁺ + 2e⁻ ⇌ Malate", -0.166, "Oxaloacetate/Malate"),
            ("Fumarate + 2H⁺ + 2e⁻ ⇌ Succinate", 0.031, "Fumarate/Succinate"),

            # ── Additional useful couples ──
            ("AgBr(s) + e⁻ ⇌ Ag(s) + Br⁻", 0.071, "Silver Bromide/Silver"),
            ("Ag2S(s) + 2e⁻ ⇌ 2Ag(s) + S²⁻", -0.69, "Silver Sulfide/Silver"),
            ("Au³⁺ + 3e⁻ ⇌ Au(s)", 1.498, "Gold(III)/Gold"),
            ("AuCl4⁻ + 3e⁻ ⇌ Au(s) + 4Cl⁻", 1.002, "Tetrachloroaurate/Gold"),
            ("Ba²⁺ + 2e⁻ ⇌ Ba(s)", -2.91, "Barium/Barium ion"),
            ("Cs⁺ + e⁻ ⇌ Cs(s)", -3.026, "Cesium/Cesium ion"),
            ("Rb⁺ + e⁻ ⇌ Rb(s)", -2.98, "Rubidium/Rubidium ion"),
            ("Sr²⁺ + 2e⁻ ⇌ Sr(s)", -2.89, "Strontium/Strontium ion"),
            ("Ca²⁺ + 2e⁻ ⇌ Ca(s)", -2.87, "Calcium/Calcium ion"),
            ("La³⁺ + 3e⁻ ⇌ La(s)", -2.38, "Lanthanum/Lanthanum ion"),
            ("Mg²⁺ + 2e⁻ ⇌ Mg(s)", -2.37, "Magnesium/Magnesium ion"),
            ("Be²⁺ + 2e⁻ ⇌ Be(s)", -1.85, "Beryllium/Beryllium ion"),
            ("Al³⁺ + 3e⁻ ⇌ Al(s)", -1.66, "Aluminum/Aluminum ion"),
            ("Ti²⁺ + 2e⁻ ⇌ Ti(s)", -1.63, "Titanium/Titanium ion"),
            ("V²⁺ + 2e⁻ ⇌ V(s)", -1.13, "Vanadium/Vanadium ion"),
            ("Cr²⁺ + 2e⁻ ⇌ Cr(s)", -0.91, "Chromium(II)/Chromium"),
            ("TiO²⁺ + 2H⁺ + e⁻ ⇌ Ti³⁺ + H2O", 0.10, "Titanium(IV)/Titanium(III)"),
            ("VO2⁺ + 2H⁺ + e⁻ ⇌ VO²⁺ + H2O", 1.00, "Vanadium(V)/Vanadium(IV)"),
            ("VO²⁺ + 2H⁺ + e⁻ ⇌ V³⁺ + H2O", 0.34, "Vanadium(IV)/Vanadium(III)"),
            ("V³⁺ + e⁻ ⇌ V²⁺", -0.26, "Vanadium(III)/Vanadium(II)"),
            ("V(eV/V)2+ + 2e⁻ ⇌ V(s)", -1.13, "Vanadium(II)/Vanadium metal"),
            ("S(s) + 2e⁻ ⇌ S²⁻", -0.48, "Sulfur/Sulfide"),
            ("SO4²⁻ + 4H⁺ + 2e⁻ ⇌ H2SO3 + H2O", 0.17, "Sulfite/Sulfurous Acid"),
            ("H2SO3 + 4H⁺ + 4e⁻ ⇌ S(s) + 3H2O", 0.45, "Sulfur/Sulfurous Acid"),
            ("S2O8²⁻ + 2e⁻ ⇌ 2SO4²⁻", 2.01, "Peroxodisulfate/Sulfate"),
            ("HSO5⁻ + 2H⁺ + 2e⁻ ⇌ HSO4⁻ + H2O", 1.81, "Peroxymonosulfate/Sulfate"),
            ("HOCl + H⁺ + e⁻ ⇌ ½Cl2 + H2O", 1.63, "Hypochlorous Acid/Chlorine"),
            ("ClO⁻ + H2O + 2e⁻ ⇌ Cl⁻ + 2OH⁻", 0.81, "Hypochlorite/Chloride (basic)"),
            ("ClO2 + e⁻ ⇌ ClO2⁻", 0.95, "Chlorine Dioxide/Chlorite"),
            ("ClO2⁻ + 2H2O + 4e⁻ ⇌ Cl⁻ + 4OH⁻", 0.78, "Chlorite/Chloride (basic)"),
            ("ClO3⁻ + 6H⁺ + 6e⁻ ⇌ Cl⁻ + 3H2O", 1.45, "Chlorate/Chloride"),
            ("ClO4⁻ + 2H⁺ + 2e⁻ ⇌ ClO3⁻ + H2O", 1.20, "Perchlorate/Chlorate"),
            ("BrO3⁻ + 6H⁺ + 6e⁻ ⇌ Br⁻ + 3H2O", 1.44, "Bromate/Bromide"),
            ("IO3⁻ + 6H⁺ + 5e⁻ ⇌ ½I2 + 3H2O", 1.20, "Iodate/Iodine"),
            ("IO3⁻ + 3H2O + 6e⁻ ⇌ I⁻ + 6OH⁻", 0.26, "Iodate/Iodide (basic)"),
            ("H5IO6 + H⁺ + 2e⁻ ⇌ IO3⁻ + 3H2O", 1.60, "Periodic Acid/Iodate"),
            ("MnO4²⁻ + 2H2O + 2e⁻ ⇌ MnO2(s) + 4OH⁻", 0.60, "Manganate/Manganese Dioxide (basic)"),
            ("MnO4⁻ + e⁻ ⇌ MnO4²⁻", 0.56, "Permanganate/Manganate"),
            ("Co(NH3)6³⁺ + e⁻ ⇌ Co(NH3)6²⁺", 0.11, "Cobalt(III)/(II) ammine"),
            ("Co(CN)6³⁻ + e⁻ ⇌ Co(CN)6⁴⁻", -0.83, "Cobalt(III)/(II) cyanide"),
            ("Fe(CN)6³⁻ + e⁻ ⇌ Fe(CN)6⁴⁻", 0.36, "Ferricyanide/Ferrocyanide"),
            ("Fe(edta)⁻ + e⁻ ⇌ Fe(edta)²⁻", 0.12, "Iron(III)/(II)-EDTA"),
            ("PbO2(s) + 4H⁺ + 2e⁻ ⇌ Pb²⁺ + 2H2O", 1.455, "Lead Dioxide/Lead(II)"),
            ("Pb²⁺ + 2e⁻ ⇌ Pb(s)", -0.126, "Lead(II)/Lead"),
            ("PbSO4(s) + 2e⁻ ⇌ Pb(s) + SO4²⁻", -0.356, "Lead Sulfate/Lead"),
            ("PbCl2(s) + 2e⁻ ⇌ Pb(s) + 2Cl⁻", -0.267, "Lead Chloride/Lead"),
            ("Hg2²⁺ + 2e⁻ ⇌ 2Hg(l)", 0.792, "Mercury(I)/Mercury"),
            ("Hg²⁺ + 2e⁻ ⇌ Hg(l)", 0.851, "Mercury(II)/Mercury"),
            ("Sn⁴⁺ + 2e⁻ ⇌ Sn²⁺", 0.151, "Tin(IV)/Tin(II)"),
            ("Sn²⁺ + 2e⁻ ⇌ Sn(s)", -0.138, "Tin(II)/Tin"),
            ("NiOOH + H2O + e⁻ ⇌ Ni(OH)2 + OH⁻", 0.49, "Nickel oxyhydroxide/Ni(OH)2 (battery)"),
            ("O2(g) + H2O + 2e⁻ ⇌ HO2⁻ + OH⁻", -0.065, "Oxygen/Hydroperoxide (basic)"),
            ("HO2⁻ + H2O + 2e⁻ ⇌ 3OH⁻", 0.88, "Hydroperoxide/Hydroxide (basic)"),
        ]

    # Common aliases for redox couples (user-friendly notation → database key)
    _aliases = {
        "fe3+/fe2+": ("Fe³⁺ + e⁻ ⇌ Fe²⁺", 0.771, "Iron(III)/Iron(II)"),
        "cu2+/cu": ("Cu²⁺ + 2e⁻ ⇌ Cu(s)", 0.337, "Copper(II)/Copper"),
        "zn2+/zn": ("Zn²⁺ + 2e⁻ ⇌ Zn(s)", -0.763, "Zinc/Zinc ion"),
        "ag+/ag": ("Ag⁺ + e⁻ ⇌ Ag(s)", 0.7996, "Silver(I)/Silver"),
        "ni2+/ni": ("Ni²⁺ + 2e⁻ ⇌ Ni(s)", -0.25, "Nickel/Nickel ion"),
        "co2+/co": ("Co²⁺ + 2e⁻ ⇌ Co(s)", -0.28, "Cobalt/Cobalt ion"),
        "cd2+/cd": ("Cd²⁺ + 2e⁻ ⇌ Cd(s)", -0.403, "Cadmium/Cadmium ion"),
        "i2/i-": ("I₂(s) + 2e⁻ ⇌ 2I⁻", 0.535, "Iodine/Iodide"),
        "br2/br-": ("Br₂(l) + 2e⁻ ⇌ 2Br⁻", 1.065, "Bromine/Bromide"),
        "cl2/cl-": ("Cl₂(g) + 2e⁻ ⇌ 2Cl⁻", 1.358, "Chlorine/Chloride"),
        "f2/f-": ("F₂(g) + 2e⁻ ⇌ 2F⁻", 2.87, "Fluorine/Fluoride"),
        "h+/h2": ("2H⁺ + 2e⁻ ⇌ H₂(g)", 0.000, "Standard Hydrogen Electrode (SHE) — REFERENCE"),
        "o2/h2o": ("O₂(g) + 4H⁺ + 4e⁻ ⇌ 2H₂O", 1.229, "Oxygen/Water (acidic)"),
        "o2/oh-": ("O₂(g) + 2H₂O + 4e⁻ ⇌ 4OH⁻", 0.401, "Oxygen/Hydroxide (basic)"),
        "mno4-/mn2+": ("MnO₄⁻ + 8H⁺ + 5e⁻ ⇌ Mn²⁺ + 4H₂O", 1.507, "Permanganate/Manganese(II)"),
        "cr2o72-/cr3+": ("Cr₂O₇²⁻ + 14H⁺ + 6e⁻ ⇌ 2Cr³⁺ + 7H₂O", 1.33, "Dichromate/Chromium(III)"),
        "no3-/no": ("NO₃⁻ + 4H⁺ + 3e⁻ ⇌ NO(g) + 2H₂O", 0.96, "Nitrate/Nitric Oxide"),
        "s4o62-/s2o32-": ("S₄O₆²⁻ + 2e⁻ ⇌ 2S₂O₃²⁻", 0.08, "Tetrathionate/Thiosulfate"),
        "sn4+/sn2+": ("Sn⁴⁺ + 2e⁻ ⇌ Sn²⁺", 0.151, "Tin(IV)/Tin(II)"),
        "li+/li": ("Li⁺ + e⁻ ⇌ Li(s)", -3.04, "Lithium/Lithium ion"),
        "k+/k": ("K⁺ + e⁻ ⇌ K(s)", -2.931, "Potassium/Potassium ion"),
        "na+/na": ("Na⁺ + e⁻ ⇌ Na(s)", -2.714, "Sodium/Sodium ion"),
        "ca2+/ca": ("Ca²⁺ + 2e⁻ ⇌ Ca(s)", -2.87, "Calcium/Calcium ion"),
        "mg2+/mg": ("Mg²⁺ + 2e⁻ ⇌ Mg(s)", -2.37, "Magnesium/Magnesium ion"),
        "al3+/al": ("Al³⁺ + 3e⁻ ⇌ Al(s)", -1.66, "Aluminum/Aluminum ion"),
        "au3+/au": ("Au³⁺ + 3e⁻ ⇌ Au(s)", 1.498, "Gold(III)/Gold"),
        "ce4+/ce3+": ("Ce⁴⁺ + e⁻ ⇌ Ce³⁺", 1.61, "Cerium(IV)/Cerium(III)"),
        "pb2+/pb": ("Pb²⁺ + 2e⁻ ⇌ Pb(s)", -0.126, "Lead(II)/Lead"),
    }

    def _run_base(self, query: str, as_oxidation: bool = False) -> dict:
        """Query standard electrode potential."""
        q = query.strip()
        q_lower = q.lower().replace(' ', '/').replace('_', '/')  # Normalize separators

        # Try exact alias match first
        if q_lower in self._aliases:
            hr_str, e0_val, desc = self._aliases[q_lower]
            result_e0 = -e0_val if as_oxidation else e0_val
            return {
                "half_reaction": f"{hr_str}    E° = {result_e0:+.3f} V ({'oxidation' if as_oxidation else 'reduction'})",
                "E0_V": round(result_e0, 4),
                "description": desc,
            }

        q_lower_for_fuzzy = q

        best_match = None
        best_score = 0

        for hr_str, e0_val, desc in self._potentials:
            # Build searchable text from the half-reaction
            search_text = f"{hr_str} {desc}".lower()

            # Exact match check
            if q_lower == search_text.strip() or q_lower in hr_str.lower():
                result_e0 = -e0_val if as_oxidation else e0_val
                return {
                    "half_reaction": f"{hr_str}    E° = {result_e0:+.3f} V ({'oxidation' if as_oxidation else 'reduction'})",
                    "E0_V": round(result_e0, 4),
                    "description": desc,
                }

            # Score partial matches
            score = self._match_score(q_lower, search_text)
            if score > best_score:
                best_score = score
                best_match = (hr_str, e0_val, desc)

        if best_match and best_score >= 0.3:
            hr_str, e0_val, desc = best_match
            result_e0 = -e0_val if as_oxidation else e0_val
            logger.info(f"Fuzzy match (score={best_score:.2f}): '{q}' → '{hr_str}'")
            return {
                "half_reaction": f"{hr_str}    E° = {result_e0:+.3f} V ({'oxidation' if as_oxidation else 'reduction'})",
                "E0_V": round(result_e0, 4),
                "description": desc + " [fuzzy matched]",
            }

        raise ChemMCPError(
            f"Cannot find standard potential for '{query}'. "
            f"Try searching by species name or redox couple (e.g., 'Cu2+/Cu', 'Fe3+/Fe2+', "
            f"'Zn2+/Zn', 'MnO4-/Mn2+', 'Ag+/Ag', 'O2/H2O', 'F2/F-', 'Li+/Li'). "
            f"Database contains ~100 common half-reactions."
        )

    def _run_text(self, query: str) -> dict:
        return self._run_base(query)

    @staticmethod
    def _match_score(query: str, text: str) -> float:
        """Calculate match score between query and entry text."""
        import re
        # Split on common delimiters but preserve meaningful tokens
        q_tokens = set(re.findall(r'[A-Za-z0-9]+', query))
        t_tokens = set(re.findall(r'[A-Za-z0-9]+', text))
        if not q_tokens:
            return 0
        overlap = len(q_tokens & t_tokens)
        return overlap / len(q_tokens)
