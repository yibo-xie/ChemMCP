import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DissolvePrecipitate(BaseTool):
    """
    分析沉淀溶解条件。
    支持酸溶解、碱溶解（两性氢氧化物）、配位溶解、氧化还原溶解等多种机制。
    """
    __version__ = "0.1.0"
    name = "DissolvePrecipitate"
    func_name = "dissolve_precipitate"
    description = "Analyze conditions to dissolve a precipitate, including acid dissolution (for carbonates, sulfides, hydroxides), base dissolution for amphoteric hydroxides, complex ion formation (e.g., AgCl in NH3), and redox dissolution."
    implementation_description = "Built-in database of common precipitates with their viable dissolution methods. For each method, calculates minimum reagent concentration using Ksp/Kf equilibrium where applicable. Covers: acid dissolution, amphoteric dissolution in strong base, complex ion formation with ligands (NH3, CN-, OH-, Cl-, NH3), and oxidative dissolution."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Dissolution", "Complex Ion", "pH", "Precipitation", "Le Chatelier", "Kf"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "N/A", "Chemical formula of the precipitate to dissolve, e.g., 'AgCl', 'Al(OH)3', 'CuS', 'CaCO3'."),
        ("method", "str", "auto", "Dissolution method: 'auto' (recommend all), 'acid', 'base', 'complex', 'redox', or specific reagent like 'NH3', 'HNO3'."),
        ("ligand_conc", "float", "None", "Optional: ligand concentration (mol/L) for complex dissolution. If None, calculates minimum required."),
        ("acid_conc", "float", "None", "Optional: acid concentration (mol/L) for acid dissolution. If None, calculates minimum required."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Format: 'compound [method] [reagent_conc]'. E.g., 'AgCl', 'Al(OH)3 base', 'CuS redox 2.0'."),
    ]

    output_sig = [
        ("compound", "str", "Compound formula analyzed."),
        ("viable_methods", "list", "List of all viable dissolution methods with details."),
        ("recommended_method", "str", "The recommended best method."),
        ("equilibrium_analysis", "str", "Detailed equilibrium calculations and reasoning."),
        ("minimum_reagent_concentration", "dict", "Minimum concentration of each reagent needed for complete dissolution (mol/L)."),
    ]

    examples = [
        {
            "code_input": {"compound": "AgCl", "method": "auto", "ligand_conc": None, "acid_conc": None},
            "text_input": {"input_str": "AgCl"},
            "output": {
                "compound": "AgCl",
                "viable_methods": [
                    {"method": "Complexation with ammonia", "reaction": "AgCl(s) + 2NH3 → [Ag(NH3)2]+ + Cl-", "kf": 1.7e7, "min_NH3": "~0.5 M"},
                    {"method": "Complexation with cyanide", "reaction": "AgCl(s) + 2CN- → [Ag(CN)2]- + Cl-", "kf": 5.6e18, "min_CN": "~0.003 M"},
                ],
                "recommended_method": "Complexation with ammonia (safe, effective; Kf=1.7×10⁷)",
                "equilibrium_analysis": "K = Ksp × Kf = (1.77×10⁻¹⁰) × (1.7×10⁷) = 3.0×10⁻³. With [NH3]=1M: S ≈ √(K×[NH3]²) ≈ 0.055 mol/L.",
                "minimum_reagent_concentration": {"NH3": "≥0.5 M for practical dissolution", "CN-": "≥0.003 M (toxic)"},
            },
        },
        {
            "code_input": {"compound": "Al(OH)3", "method": "auto", "ligand_conc": None, "acid_conc": None},
            "text_input": {"input_str": "Al(OH)3"},
            "output": {
                "compound": "Al(OH)3",
                "viable_methods": [
                    {"method": "Strong base (amphoteric)", "reaction": "Al(OH)3(s) + OH- → [Al(OH)4]-", "note": "Amphoteric hydroxide"},
                    {"method": "Acid dissolution", "reaction": "Al(OH)3(s) + 3H+ → Al3+ + 3H2O", "note": "Standard acid-base neutralization"},
                ],
                "recommended_method": "Strong base (NaOH) — Al(OH)3 is distinctly amphoteric",
                "equilibrium_analysis": "As amphoteric: K(formation of aluminate) is large. pH > ~10 effectively dissolves Al(OH)3.",
                "minimum_reagent_concentration": {"NaOH": "≥0.1 M", "HCl/HNO3": "≥0.01 M"},
            },
        },
        {
            "code_input": {"compound": "CuS", "method": "redox", "ligand_conc": None, "acid_conc": None},
            "text_input": {"input_str": "CuS redox"},
            "output": {
                "compound": "CuS",
                "viable_methods": [
                    {"method": "Oxidative dissolution (hot HNO3)", "reaction": "3CuS + 8HNO3 → 3Cu(NO3)2 + 3S + 2NO + 4H2O", "note": "Ksp too small for simple acid/complex dissolution"},
                ],
                "recommended_method": "Hot concentrated HNO3 (oxidative)",
                "equilibrium_analysis": "CuS Ksp=6×10⁻³⁶, extremely insoluble. Non-oxidizing acids cannot provide sufficient [S²⁻] reduction. HNO3 oxidizes S²⁻ to S°, driving dissolution.",
                "minimum_reagent_concentration": {"HNO3": "conc. hot (~6-12 M)"},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize precipitate dissolution database."""
        # Format: {compound: {methods: [...], ksp, molar_mass, notes}}
        self._db = {
            # --- Halides ---
            "AgCl": {
                "ksp": 1.77e-10, "mw": 143.32,
                "methods": [
                    {
                        "type": "complex",
                        "name": "Ammonia complexation",
                        "reaction": "AgCl(s) + 2NH₃(aq) ⇌ [Ag(NH₃)₂]⁺(aq) + Cl⁻(aq)",
                        "ligand": "NH3", "coord": 2,
                        "kf": 1.7e7,
                        "k_overall": 1.77e-10 * 1.7e7,  # Ksp × Kf
                        "min_conc": 0.5,  # approximate min [NH3] in M
                        "notes": "Forms colorless diammine silver(I) complex. Most practical lab method.",
                    },
                    {
                        "type": "complex",
                        "name": "Cyanide complexation",
                        "reaction": "AgCl(s) + 2CN⁻(aq) ⇌ [Ag(CN)₂]⁻(aq) + Cl⁻(aq)",
                        "ligand": "CN-", "coord": 2,
                        "kf": 5.6e18,
                        "k_overall": 1.77e-10 * 5.6e18,
                        "min_conc": 0.003,
                        "notes": "Very effective but highly toxic. Used in electroplating baths.",
                    },
                    {
                        "type": "complex",
                        "name": "Thiosulfate complexation",
                        "reaction": "AgCl(s) + 2S₂O₃²⁻(aq) ⇌ [Ag(S₂O₃)₂]³⁻(aq) + Cl⁻(aq)",
                        "ligand": "S2O3^2-", "coord": 2,
                        "kf": 6.9e13,  # approx for Ag(S2O3)2^3-
                        "k_overall": 1.77e-10 * 6.9e13,
                        "min_conc": 0.01,
                        "notes": "Used in photography fixer. Forms stable thiosilver complex.",
                    },
                ],
            },
            "AgBr": {
                "ksp": 5.35e-13, "mw": 187.77,
                "methods": [
                    {
                        "type": "complex", "name": "Ammonia complexation",
                        "reaction": "AgBr(s) + 2NH₃ ⇌ [Ag(NH₃)₂]⁺ + Br⁻",
                        "ligand": "NH3", "coord": 2, "kf": 1.7e7,
                        "k_overall": 5.35e-13 * 1.7e7, "min_conc": 2.0,
                        "notes": "Requires higher [NH3] than AgCl due to lower Ksp.",
                    },
                    {
                        "type": "complex", "name": "Thiosulfate (photography fixer)",
                        "reaction": "AgBr(s) + 2S₂O₃²⁻ ⇌ [Ag(S₂O₃)₂]³⁻ + Br⁻",
                        "ligand": "S2O3^2-", "coord": 2, "kf": 6.9e13,
                        "k_overall": 5.35e-13 * 6.9e13, "min_conc": 0.05,
                        "notes": "Standard photographic fixing agent for AgBr emulsions.",
                    },
                    {
                        "type": "complex", "name": "Cyanide complexation",
                        "reaction": "AgBr(s) + 2CN⁻ ⇌ [Ag(CN)₂]⁻ + Br⁻",
                        "ligand": "CN-", "coord": 2, "kf": 5.6e18,
                        "k_overall": 5.35e-13 * 5.6e18, "min_conc": 0.001,
                        "notes": "Very effective but extremely toxic.",
                    },
                ],
            },
            "AgI": {
                "ksp": 8.52e-17, "mw": 234.77,
                "methods": [
                    {
                        "type": "complex", "name": "Cyanide complexation",
                        "reaction": "AgI(s) + 2CN⁻ ⇌ [Ag(CN)₂]⁻ + I⁻",
                        "ligand": "CN-", "coord": 2, "kf": 5.6e18,
                        "k_overall": 8.52e-17 * 5.6e18, "min_conc": 0.02,
                        "notes": "NH3 insufficient due to very low Ksp. CN- or Na2S2O3 needed.",
                    },
                    {
                        "type": "complex", "name": "Thiosulfate complexation",
                        "reaction": "AgI(s) + 2S₂O₃²⁻ ⇌ [Ag(S₂O₃)₂]³⁻ + I⁻",
                        "ligand": "S2O3^2-", "coord": 2, "kf": 6.9e13,
                        "k_overall": 8.52e-17 * 6.9e13, "min_conc": 0.3,
                        "notes": "Possible but requires concentrated thiosulfate.",
                    },
                ],
            },

            # --- Hydroxides (including amphoteric) ---
            "Al(OH)3": {
                "ksp": 3.28e-34, "mw": 78.00,
                "methods": [
                    {
                        "type": "base", "name": "Strong base (amphoteric dissolution)",
                        "reaction": "Al(OH)₃(s) + OH⁻(aq) ⇌ [Al(OH)₄]⁻(aq)",
                        "ligand": "OH-", "coord": 1, "kf": 1.1e33,
                        "k_overall": 3.28e-34 * 1.1e33, "min_conc": 0.05,
                        "notes": "Classic amphoteric behavior. Dissolves in excess NaOH/KOH to form aluminate.",
                    },
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Al(OH)₃(s) + 3H⁺(aq) ⇌ Al³⁺(aq) + 3H₂O(l)",
                        "ligand": "H+", "coord": 3, "kf": None,  # essentially irreversible
                        "k_overall": None, "min_conc": 0.01,
                        "notes": "Dissolves readily in both strong acids and strong bases.",
                    },
                ],
            },
            "Zn(OH)2": {
                "ksp": 3.00e-17, "mw": 99.40,
                "methods": [
                    {
                        "type": "base", "name": "Strong base (amphoteric)",
                        "reaction": "Zn(OH)₂(s) + 2OH⁻ ⇌ [Zn(OH)₄]²⁻(aq)",
                        "ligand": "OH-", "coord": 2, "kf": 2.9e15,
                        "k_overall": 3.00e-17 * 2.9e15, "min_conc": 0.03,
                        "notes": "Amphoteric. Dissolves in excess NaOH to form zincate.",
                    },
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Zn(OH)₂(s) + 2H⁺ ⇌ Zn²⁺(aq) + 2H₂O(l)",
                        "ligand": "H+", "coord": 2, "kf": None, "min_conc": 0.01,
                        "notes": "Also dissolves in acids (e.g., HCl).",
                    },
                    {
                        "type": "complex", "name": "Ammonia complexation",
                        "reaction": "Zn(OH)₂(s) + 4NH₃ ⇌ [Zn(NH₃)₄]²⁺(aq) + 2OH⁻(aq)",
                        "ligand": "NH3", "coord": 4, "kf": 2.9e9,
                        "k_overall": None, "min_conc": 1.0,
                        "notes": "Dissolves in excess ammonia solution.",
                    },
                ],
            },
            "Cr(OH)3": {
                "ksp": 6.30e-31, "mw": 103.02,
                "methods": [
                    {
                        "type": "base", "name": "Strong base (amphoteric)",
                        "reaction": "Cr(OH)₃(s) + OH⁻ ⇌ [Cr(OH)₄]⁻(aq)",
                        "ligand": "OH-", "coord": 1, "kf": 1e28,  # approx
                        "k_overall": 6.30e-31 * 1e28, "min_conc": 0.1,
                        "notes": "Amphoteric. Dissolves slowly in concentrated NaOH.",
                    },
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Cr(OH)₃(s) + 3H⁺ ⇌ Cr³⁺(aq) + 3H₂O(l)",
                        "ligand": "H+", "coord": 3, "min_conc": 0.01,
                    },
                ],
            },
            "Pb(OH)2": {
                "ksp": 1.43e-20, "mw": 241.21,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution (or nitric acid)",
                        "reaction": "Pb(OH)₂(s) + 2H⁺ ⇌ Pb²⁺(aq) + 2H₂O(l)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Weakly amphoteric; dissolves in hot concentrated NaOH slowly.",
                    },
                    {
                        "type": "base", "name": "Concentrated strong base (slow)",
                        "reaction": "Pb(OH)₂(s) + OH⁻ ⇌ [Pb(OH)₃]⁻(aq)",
                        "ligand": "OH-", "coord": 1, "kf": 1e4,  # weakly amphoteric
                        "k_overall": 1.43e-20 * 1e4, "min_conc": 2.0,
                        "notes": "Only slightly amphoteric; requires very concentrated base.",
                    },
                ],
            },
            "Cu(OH)2": {
                "ksp": 2.20e-20, "mw": 97.56,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Cu(OH)₂(s) + 2H⁺ ⇌ Cu²⁺(aq) + 2H₂O(l)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Dissolves readily in dilute acids.",
                    },
                    {
                        "type": "complex", "name": "Ammonia complexation",
                        "reaction": "Cu(OH)₂(s) + 4NH₃ ⇌ [Cu(NH₃)₄]²⁺(aq) + 2OH⁻(aq)",
                        "ligand": "NH3", "coord": 4, "kf": 2.1e13,
                        "k_overall": None, "min_conc": 0.5,
                        "notes": "Dissolves in excess NH3 to form deep blue [Cu(NH3)4]2+ complex.",
                    },
                ],
            },
            "Fe(OH)3": {
                "ksp": 2.79e-39, "mw": 106.87,
                "methods": [
                    {
                        "type": "acid", "name": "Strong acid dissolution",
                        "reaction": "Fe(OH)₃(s) + 3H⁺ ⇌ Fe³⁺(aq) + 3H₂O(l)",
                        "ligand": "H+", "coord": 3, "min_conc": 0.001,
                        "notes": "Not amphoteric. Only dissolves in strong acid.",
                    },
                ],
            },
            "Mg(OH)2": {
                "ksp": 5.61e-12, "mw": 58.33,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Mg(OH)₂(s) + 2H⁺ ⇌ Mg²⁺(aq) + 2H₂O(l)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Not amphoteric. Dissolves in acids (even weak acids like acetic acid).",
                    },
                ],
            },

            # --- Carbonates ---
            "CaCO3": {
                "ksp": 3.36e-9, "mw": 100.09,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution (CO2 evolution)",
                        "reaction": "CaCO₃(s) + 2H⁺(aq) ⇌ Ca²⁺(aq) + CO₂(g) + H₂O(l)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Effervesces (bubbles CO2) with any strong acid. Also dissolves in water saturated with CO2 (forms Ca(HCO3)2).",
                    },
                ],
            },
            "BaCO3": {
                "ksp": 2.58e-9, "mw": 197.34,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution (CO2 evolution)",
                        "reaction": "BaCO₃(s) + 2H⁺ ⇌ Ba²⁺(aq) + CO₂(g) + H₂O(l)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Effervesces with strong acids like HCl.",
                    },
                ],
            },
            "Ag2CO3": {
                "ksp": 8.46e-12, "mw": 275.75,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Ag₂CO₃(s) + 2H⁺ ⇌ 2Ag⁺(aq) + CO₂(g) + H₂O(l)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Also forms Ag2O on heating. Dissolves in nitric acid.",
                    },
                ],
            },

            # --- Sulfides ---
            "CuS": {
                "ksp": 6.00e-36, "mw": 95.61,
                "methods": [
                    {
                        "type": "redox", "name": "Oxidative dissolution (hot conc. HNO3)",
                        "reaction": "3CuS(s) + 8HNO₃(conc.) → 3Cu(NO₃)₂(aq) + 3S(s) + 2NO(g) + 4H₂O(l)",
                        "ligand": "HNO3", "min_conc": 6.0,
                        "notes": "Ksp far too small for non-oxidizing acid. Hot concentrated HNO3 oxidizes S²⁻ to S°, shifting equilibrium.",
                    },
                    {
                        "type": "redox", "name": "Aqua regia / royal water",
                        "reaction": "CuS(s) + oxidant → Cu²⁺ + S/SO₄²⁻",
                        "ligand": "aqua_regia", "min_conc": None,
                        "notes": "Aqua regia (HCl:HNO3 = 3:1) can dissolve most stubborn sulfides.",
                    },
                ],
            },
            "HgS": {
                "ksp": 4.00e-53, "mw": 232.66,
                "methods": [
                    {
                        "type": "redox", "name": "Aqua regia (oxidation + complexation)",
                        "reaction": "3HgS(s) + 2HNO₃(conc.) + 12HCl(conc.) → 3H₂[HgCl₄](aq) + 3S(s) + 2NO(g) + 4H₂O(l)",
                        "ligand": "aqua_regia", "min_conc": None,
                        "notes": "Most insoluble common sulfide. Only aqua regia (or bromine water + HCl) can dissolve it. Oxidizes S²⁻ while Cl- complexes Hg2+.",
                    },
                    {
                        "type": "complex", "name": "Sulfide leaching (Na2S + NaOH)",
                        "reaction": "HgS(s) + S²⁻(aq) ⇌ [HgS₂]²⁻(aq)",
                        "ligand": "S^2-", "coord": 1, "kf": 1e52,  # polysulfide
                        "k_overall": 4.00e-53 * 1e52, "min_conc": 0.5,
                        "notes": "Polysulfide dissolves HgS via formation of thiomercurate(II) complex.",
                    },
                ],
            },
            "ZnS": {
                "ksp": 2.50e-22, "mw": 97.46,
                "methods": [
                    {
                        "type": "acid", "name": "Dilute strong acid (pH < 0.5)",
                        "reaction": "ZnS(s) + 2H⁺ ⇌ Zn²⁺(aq) + H₂S(g/aq)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.3,
                        "notes": "Moderately insoluble. Dissolves in dilute HCl (unlike CuS, HgS). H2S gas released (rotten egg smell).",
                    },
                ],
            },
            "MnS": {
                "ksp": 3.00e-14, "mw": 87.00,
                "methods": [
                    {
                        "type": "acid", "name": "Even weak acid (acetic acid works)",
                        "reaction": "MnS(s) + 2H⁺ ⇌ Mn²⁺(aq) + H₂S(g/aq)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.01,
                        "notes": "Relatively soluble sulfide. Even CH3COOH can dissolve it. Distinguishes from more insoluble sulfides.",
                    },
                ],
            },
            "CdS": {
                "ksp": 8.00e-27, "mw": 144.47,
                "methods": [
                    {
                        "type": "acid", "name": "Concentrated strong acid (>6M HCl, hot)",
                        "reaction": "CdS(s) + 2H⁺ ⇌ Cd²⁺(aq) + H₂S(g/aq)",
                        "ligand": "H+", "coord": 2, "min_conc": 2.0,
                        "notes": "Needs concentrated hot HCl. Between ZnS and CuS in acid solubility hierarchy.",
                    },
                ],
            },
            "PbS": {
                "ksp": 9.04e-29, "mw": 239.27,
                "methods": [
                    {
                        "type": "redox", "name": "Oxidative dissolution (HNO3)",
                        "reaction": "3PbS(s) + 8HNO₃ → 3Pb(NO₃)₂ + 3S + 2NO + 4H₂O",
                        "ligand": "HNO3", "min_conc": 4.0,
                        "notes": "Very insoluble. Requires hot concentrated HNO3 for oxidation.",
                    },
                ],
            },
            "FeS": {
                "ksp": 6.30e-19, "mw": 87.91,
                "methods": [
                    {
                        "type": "acid", "name": "Dilute strong acid",
                        "reaction": "FeS(s) + 2H⁺ ⇌ Fe²⁺(aq) + H₂S(g/aq)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.05,
                        "notes": "Releases H2S gas. Dissolves in moderately dilute acids.",
                    },
                ],
            },
            "Ag2S": {
                "ksp": 6.30e-50, "mw": 247.80,
                "methods": [
                    {
                        "type": "redox", "name": "Hot concentrated HNO3 (oxidation)",
                        "reaction": "3Ag₂S(s) + 8HNO₃(conc., hot) → 6AgNO₃ + 3S + 2NO + 4H₂O",
                        "ligand": "HNO3", "min_conc": 8.0,
                        "notes": "Extremely insoluble. Requires hot concentrated HNO3 or fusion with Na2CO3.",
                    },
                ],
            },

            # --- Sulfates ---
            "BaSO4": {
                "ksp": 1.08e-10, "mw": 233.39,
                "methods": [
                    {
                        "type": "complex", "name": "Convert to carbonate then acid-dissolve",
                        "reaction": "BaSO₄(s) + Na₂CO₃(sat., hot) ⇌ BaCO₃(s) + SO₄²⁻(aq); then BaCO₃ + 2H⁺ → Ba²⁺ + CO₂ + H₂O",
                        "ligand": "CO3^2-", "min_conc": None,
                        "notes": "Cannot be dissolved by simple acid/base. Must first convert to carbonate via prolonged boiling with saturated Na2CO3 (Ksp conversion), then dissolve BaCO3 with acid.",
                    },
                ],
            },
            "PbSO4": {
                "ksp": 2.53e-8, "mw": 303.26,
                "methods": [
                    {
                        "type": "complex", "name": "Acetate buffer / conversion to basic salt",
                        "reaction": "PbSO₄(s) + 2CH₃COO⁻(aq) ⇌ Pb(CH₃COO)₂(aq) + SO₄²⁻(aq)",
                        "ligand": "CH3COO-", "min_conc": 1.0,
                        "notes": "Dissolves in concentrated ammonium acetate or hot sodium acetate solution.",
                    },
                    {
                        "type": "complex", "name": "Conversion to carbonate",
                        "reaction": "PbSO₄(s) + CO₃²⁻(hot) → PbCO₃(s) + SO₄²⁻(aq); then acid",
                        "ligand": "CO3^2-", "min_conc": None,
                        "notes": "Can also convert to carbonate then dissolve.",
                    },
                ],
            },
            "CaSO4": {
                "ksp": 4.93e-5, "mw": 136.14,
                "methods": [
                    {
                        "type": "complex", "name": "Convert to less soluble form or use excess water",
                        "reaction": "CaSO₄(s) ⇌ Ca²⁺(aq) + SO₄²⁻(aq)  (moderately soluble)",
                        "ligand": None, "min_conc": None,
                        "notes": "Moderately soluble (Ksp~5×10⁻⁵). Dissolves in large volume of water or with (NH4)2SO4 as common ion effect reversal via complex formation.",
                    },
                ],
            },

            # --- Chromates ---
            "BaCrO4": {
                "ksp": 1.17e-10, "mw": 253.33,
                "methods": [
                    {
                        "type": "redox", "name": "Acid reduction (reduce Cr(VI) to Cr(III))",
                        "reaction": "2BaCrO₄(s) + 6H⁺ + 3SO₃²⁻ → 2Ba²⁺ + 2Cr³⁺ + 3SO₄²⁻ + 3H₂O",
                        "ligand": "acid + reducing agent", "min_conc": None,
                        "notes": "Reduce Cr(VI) to Cr(III) in acidic medium to destroy chromate ion, then Ba2+ goes into solution.",
                    },
                    {
                        "type": "acid", "name": "Strong acid (converts chromate to dichromate, shifts equilibrium)",
                        "reaction": "2BaCrO₄(s) + 2H⁺ ⇌ 2Ba²⁺ + Cr₂O₇²⁻ + H₂O",
                        "ligand": "H+", "coord": 2, "min_conc": 1.0,
                        "notes": "In strong acid, CrO4^2- converts to Cr2O7^2-, reducing [CrO4^2-] and driving dissolution. Limited effectiveness.",
                    },
                ],
            },
            "PbCrO4": {
                "ksp": 2.8e-13, "mw": 323.20,
                "methods": [
                    {
                        "type": "redox", "name": "Basic reduction (NaOH + reducing agent)",
                        "reaction": "PbCrO₄(s) + OH⁻ → PbO(s) + CrO₄²⁻(aq); then reduce Cr(VI)",
                        "ligand": "NaOH", "min_conc": 2.0,
                        "notes": "Chrome yellow pigment. Dissolves in NaOH to release chromate, then reduce Cr(VI) to soluble Cr(III).",
                    },
                ],
            },

            # --- Phosphates ---
            "Ca3(PO4)2": {
                "ksp": 2.07e-33, "mw": 310.18,
                "methods": [
                    {
                        "type": "acid", "name": "Strong acid dissolution",
                        "reaction": "Ca₃(PO₄)₂(s) + 6H⁺ ⇌ 3Ca²⁺(aq) + 2H₃PO₄(aq)",
                        "ligand": "H+", "coord": 6, "min_conc": 0.1,
                        "notes": "Dissolves in strong acids to form phosphoric acid.",
                    },
                ],
            },
            "Ag3PO4": {
                "ksp": 8.89e-17, "mw": 418.58,
                "methods": [
                    {
                        "type": "acid", "name": "Acid dissolution",
                        "reaction": "Ag₃PO₄(s) + 3H⁺ ⇌ 3Ag⁺(aq) + H₃PO₄(aq)",
                        "ligand": "H+", "coord": 3, "min_conc": 0.1,
                        "notes": "Yellow precipitate. Dissolves in HNO3.",
                    },
                ],
            },

            # --- Fluorides ---
            "CaF2": {
                "ksp": 3.45e-11, "mw": 78.08,
                "methods": [
                    {
                        "type": "complex", "name": "Complexation with borate/aluminum salts in hot conc. acid",
                        "reaction": "CaF₂(s) + 2H⁺ ⇌ Ca²⁺(aq) + 2HF(aq); HF complexes with H₃BO₃ or Al³⁺",
                        "ligand": "H+ + complexing agent", "min_conc": 1.0,
                        "notes": "Very difficult to dissolve. Use conc. H2SO4 to form HF gas, or conc. HClO4. In analysis, often converted to sulfate by fusion with SiO2.",
                    },
                    {
                        "type": "complex", "name": "Aluminum salt complexation (in acid)",
                        "reaction": "6F⁻(from CaF₂) + Al³⁺ → [AlF₆]³⁻(aq)",
                        "ligand": "Al3+", "coord": 6, "kf": 6.9e19,
                        "k_overall": 3.45e-11 * 6.9e19, "min_conc": 0.1,
                        "notes": "Al3+ complexes F- strongly, pulling dissolution forward.",
                    },
                ],
            },

            # --- Oxalates ---
            "CaC2O4": {
                "ksp": 2.32e-9, "mw": 128.10,
                "methods": [
                    {
                        "type": "redox", "name": "Oxidative dissolution (permanganate/sulfuric acid)",
                        "reaction": "5CaC₂O₄(s) + 2MnO₄⁻ + 16H⁺ → 5Ca²⁺ + 10CO₂ + 2Mn²⁺ + 8H₂O",
                        "ligand": "KMnO4/H2SO4", "min_conc": None,
                        "notes": "Used in redox titrations. Permanganate oxidizes oxalate to CO2 in hot acidic solution.",
                    },
                    {
                        "type": "acid", "name": "Strong acid dissolution",
                        "reaction": "CaC₂O₄(s) + 2H⁺ ⇌ Ca²⁺(aq) + H₂C₂O₄(aq)",
                        "ligand": "H+", "coord": 2, "min_conc": 0.5,
                        "notes": "Dissolves in strong mineral acids (oxalic acid is moderately weak).",
                    },
                ],
            },
        }

        self._aliases = {
            "silver chloride": "AgCl", "silver bromide": "AgBr", "silver iodide": "AgI",
            "aluminum hydroxide": "Al(OH)3", "zinc hydroxide": "Zn(OH)2",
            "chromium(iii) hydroxide": "Cr(OH)3", "lead hydroxide": "Pb(OH)2",
            "copper(ii) hydroxide": "Cu(OH)2", "iron(iii) hydroxide": "Fe(OH)3",
            "magnesium hydroxide": "Mg(OH)2",
            "calcium carbonate": "CaCO3", "barium carbonate": "BaCO3", "silver carbonate": "Ag2CO3",
            "copper(ii) sulfide": "CuS", "mercury sulfide": "HgS", "zinc sulfide": "ZnS",
            "manganese(ii) sulfide": "MnS", "cadmium sulfide": "CdS", "lead sulfide": "PbS",
            "ferrous sulfide": "FeS", "silver sulfide": "Ag2S",
            "barium sulfate": "BaSO4", "lead sulfate": "PbSO4", "calcium sulfate": "CaSO4",
            "barium chromate": "BaCrO4", "lead chromate": "PbCrO4",
            "calcium phosphate": "Ca3(PO4)2", "silver phosphate": "Ag3PO4",
            "calcium fluoride": "CaF2", "calcium oxalate": "CaC2O4",
        }

    def _run_base(self, compound: str, method: str = "auto",
                  ligand_conc: Optional[float] = None,
                  acid_conc: Optional[float] = None) -> dict:
        """Core logic: analyze dissolution methods for a precipitate."""
        key = self._resolve(compound)
        if key not in self._db:
            raise ChemMCPError(
                f"Compound '{compound}' not found in dissolution database. "
                f"Available ({len(self._db)}): {sorted(self._db.keys())[:25]}{'...' if len(self._db)>25 else ''}"
            )

        info = self._db[key]
        methods = info["methods"]

        # Filter by requested method type
        if method and method.lower() != "auto":
            mt = method.lower()
            filtered = [m for m in methods if mt in m.get("type", "").lower() or mt in m.get("name", "").lower()]
            if not filtered:
                filtered = [m for m in methods if mt in str(m.get("ligand", "")).lower()]
            if filtered:
                methods = filtered

        # Build detailed method list
        viable = []
        recommended = None
        best_score = -1
        min_reagents = {}

        for m in methods:
            entry = {
                "method": m.get("name", "Unknown"),
                "reaction": m.get("reaction", ""),
                "type": m.get("type", ""),
            }
            if "kf" in m:
                entry["formation_constant_Kf"] = m["kf"]
            if "k_overall" in m and m["k_overall"] is not None:
                entry["overall_equilibrium_K"] = f"{m['k_overall']:.3e}"
            if "min_conc" in m:
                entry["approximate_min_concentration_M"] = m["min_conc"]
                lig = m.get("ligand", "")
                if lig:
                    min_reagents[lig] = f"≥{m['min_conc']} M"
            if "notes" in m:
                entry["notes"] = m["notes"]
            viable.append(entry)

            # Score for recommendation (prefer safe, simple methods)
            score = 0
            if m.get("type") == "acid":
                score += 3  # simplest
            elif m.get("type") == "base":
                score += 2
            elif m.get("type") == "complex":
                lig_name = m.get("ligand", "")
                if "CN" in lig_name:
                    score -= 2  # toxic penalty
                elif "NH3" in lig_name:
                    score += 2  # relatively safe
                else:
                    score += 1
            elif m.get("type") == "redox":
                score += 0  # last resort

            if score > best_score:
                best_score = score
                recommended = m.get("name", "Unknown")

        # Equilibrium analysis
        ksp = info["ksp"]
        eq_parts = [f"Compound: {key}, Ksp = {ksp:.3e}"]
        for m in info["methods"]:
            ko = m.get("k_overall")
            kf_val = m.get("kf")
            if ko is not None and kf_val:
                eq_parts.append(
                    f"\n• {m['name']}: K = Ksp × Kf = ({ksp:.2e}) × ({kf_val:.2e}) = {ko:.3e}"
                )
                if ligand_conc is not None and m.get("coord"):
                    coord = m["coord"]
                    # Rough solubility estimate: S ≈ (K × [L]^coord)^(1/n) simplified
                    try:
                        s_est = (ko * (ligand_conc ** coord)) ** 0.5 if coord else ko ** 0.5
                        eq_parts.append(f"  At [L]={ligand_conc} M, estimated S ≈ {s_est:.3e} mol/L")
                    except (OverflowError, ValueError):
                        pass

        eq_analysis = "\n".join(eq_parts)

        logger.info(f"DissolvePrecipitate: {key} → {len(viable)} methods, recommended={recommended}")
        return {
            "compound": key,
            "viable_methods": viable,
            "recommended_method": recommended or viable[0]["method"] if viable else "No method found",
            "equilibrium_analysis": eq_analysis,
            "minimum_reagent_concentration": min_reagents,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        parts = input_str.strip().split()
        compound = parts[0]
        method = parts[1] if len(parts) > 1 else "auto"
        conc = float(parts[2]) if len(parts) > 2 else None
        return self._run_base(compound, method, conc)

    def _resolve(self, name: str) -> str:
        n = name.strip()
        if n in self._db:
            return n
        nl = n.lower()
        if nl in self._aliases:
            return self._aliases[nl]
        for k in self._db:
            if k.lower() == nl:
                return k
        return n
