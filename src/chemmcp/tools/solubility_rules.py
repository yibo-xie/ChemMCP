import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SolubilityRules(BaseTool):
    """
    查询溶解性规则（定性判断）。
    基于标准通用化学溶解性规则表，预测离子化合物在水中的溶解性。
    覆盖所有常见离子类型：硝酸盐、乙酸盐、卤化物、硫酸盐、碳酸盐、磷酸盐、氢氧化物、硫化物、铬酸盐等。
    """
    __version__ = "0.1.0"
    name = "SolubilityRules"
    func_name = "solubility_rules"
    description = "Query solubility rules for qualitative prediction of whether an ionic compound dissolves in water at room temperature. Covers all common ion types with exceptions and detailed rule explanations."
    implementation_description = "Comprehensive built-in solubility rule set based on standard general chemistry guidelines. Rule priority system: always-soluble ions first, then conditional rules with exceptions. Returns solubility verdict, applicable rule, and all relevant exception information."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Solubility Rules", "Qualitative Analysis", "Ionic Compounds", "General Chemistry", "Solubility Guidelines"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "None", "Chemical formula to check (e.g., 'NaCl', 'BaSO4', 'AgCl', 'CaCO3'). If empty or None, returns all rules."),
        ("cation", "str", "None", "Cation to look up rules for (e.g., 'Na+', 'Ba2+', 'Ag+'). Optional if compound is given."),
        ("anion", "str", "None", "Anion to look up rules for (e.g., 'Cl-', 'SO4^2-', 'CO3^2-'). Optional if compound is given."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Compound formula or ion name. E.g., 'NaCl', 'BaSO4', 'PbI2', 'NH4NO3'. Leave empty for all rules."),
    ]

    output_sig = [
        ("solubility", "str", "'soluble', 'slightly soluble', or 'insoluble' in water at 25°C."),
        ("rule_applied", "str", "The primary solubility rule that determines the result."),
        ("exceptions", "list", "List of applicable exceptions (if any)."),
        ("compound_analyzed", "str", "The compound that was analyzed."),
        ("cation", "str", "Cation identified from compound."),
        ("anion", "str", "Anion identified from compound."),
        ("detailed_rules", "list", "All relevant solubility rules for reference."),
        ("notes", "str", "Additional notes about this specific compound."),
    ]

    examples = [
        {
            "code_input": {"compound": "NaCl", "cation": None, "anion": None},
            "text_input": {"input_str": "NaCl"},
            "output": {
                "solubility": "soluble",
                "rule_applied": "Group 1 cations (Na+) form soluble compounds with all common anions.",
                "exceptions": [],
                "compound_analyzed": "NaCl",
                "cation": "Na+",
                "anion": "Cl-",
                "notes": "Common table salt. Solubility: ~36 g/100mL water at 25°C.",
                "detailed_rules": [],
            },
        },
        {
            "code_input": {"compound": "BaSO4", "cation": None, "anion": None},
            "text_input": {"input_str": "BaSO4"},
            "output": {
                "solubility": "insoluble",
                "rule_applied": "Most sulfate salts are soluble EXCEPT those of Sr2+, Ba2+, Pb2+, Ca2+ (slightly).",
                "exceptions": ["Ba2+ is a listed exception for sulfate solubility."],
                "compound_analyzed": "BaSO4",
                "cation": "Ba2+",
                "anion": "SO4^2-",
                "notes": "Barium sulfate is used in X-ray imaging (barium meal) due to its insolubility.",
                "detailed_rules": [],
            },
        },
        {
            "code_input": {"compound": "AgCl", "cation": None, "anion": None},
            "text_input": {"input_str": "AgCl"},
            "output": {
                "solubility": "insoluble",
                "rule_applied": "Most chloride, bromide, and iodide salts are soluble EXCEPT those of Ag+, Pb2+, Hg2^2+.",
                "exceptions": ["Ag+ forms insoluble halides (except AgF)."],
                "compound_analyzed": "AgCl",
                "cation": "Ag+",
                "anion": "Cl-",
                "notes": "White precipitate; darkens on exposure to light. Ksp=1.77×10⁻¹⁰.",
                "detailed_rules": [],
            },
        },
        {
            "code_input": {"compound": None, "cation": None, "anion": None},
            "text_input": {"input_str": ""},
            "output": {
                "solubility": "N/A (reference mode)",
                "rule_applied": "Full rule reference requested",
                "exceptions": [],
                "compound_analyzed": "ALL COMPOUNDS (reference)",
                "cation": "N/A",
                "anion": "N/A",
                "detailed_rules": "[full list of all solubility rules]",
                "notes": "Complete solubility rule reference returned.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize comprehensive solubility rules database."""
        # --- ALWAYS SOLUBLE (no exceptions) ---
        self._always_soluble_cations = [
            "Group 1 alkali metal ions: Li+, Na+, K+, Rb+, Cs+",
            "Ammonium ion: NH4+",
        ]
        self._always_soluble_anions = [
            "Nitrate: NO3-",
            "Acetate: CH3COO- (or C2H3O2-)",
            "Chlorate: ClO3-",
            "Perchlorate: ClO4-",
        ]

        # --- CONDITIONAL RULES ---
        # Each rule: {anion_group, general_rule, exceptions, slightly_soluble_exceptions, notes}
        self._rules = [
            {
                "group": "Halides (Cl-, Br-, I-)",
                "anions": ["Cl-", "Br-", "I-"],
                "general_rule": "Most chloride, bromide, and iodide salts are SOLUBLE.",
                "solubility": "soluble",
                "exceptions": {
                    "Ag+": "insoluble (AgCl white, AgBr pale yellow, AgI yellow)",
                    "Pb2+": "insoluble (PbCl2 moderately cold-soluble, PbBr2/PbI2 insol.)",
                    "Hg2^2+": "insoluble (Hg2Cl2 white / calomel)",
                    "Cu+": "insoluble (CuCl, CuBr, CuI)",
                    "Tl+": "insoluble (TlCl slightly)",
                    "Hg2+": "insoluble (mercuric halides)",
                },
                "slightly_soluble": {
                    "Pb2+": "PbCl2 has moderate solubility in cold water (~0.01 M); more soluble hot",
                },
                "notes": "All fluorides (F-) follow different rules — see Fluorides section. AgF is soluble!",
            },
            {
                "group": "Fluorides (F-)",
                "anions": ["F-"],
                "general_rule": "Most fluoride salts are INSOLUBLE or slightly soluble.",
                "solubility": "generally insoluble",
                "exceptions": {
                    "Group 1 (Li+, Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE",
                    "Ag+": "SOLUBLE (unlike other silver halides!)",
                    "Sn2+": "soluble",
                    "Bi3+": "soluble",
                },
                "slightly_soluble": {
                    "Ca2+": "slightly soluble",
                    "Sr2+": "slightly soluble",
                    "Ba2+": "slightly soluble",
                    "Pb2+": "slightly soluble",
                    "Mg2+": "slightly soluble to insoluble",
                    "Fe2+": "slightly soluble to insoluble",
                    "Ni2+": "slightly soluble to insoluble",
                    "Zn2+": "slightly soluble to insoluble",
                    "Mn2+": "slightly soluble to insoluble",
                },
                "notes": "F- is small and highly charged → high lattice energy → many fluorides insoluble. CaF2 (fluorite) Ksp=3.45×10⁻¹¹.",
            },
            {
                "group": "Sulfates (SO4^2-)",
                "anions": ["SO4^2-", "SO4(2-)"],
                "general_rule": "Most sulfate salts are SOLUBLE.",
                "solubility": "soluble",
                "exceptions": {
                    "Sr2+": "INSOLUBLE (SrSO4)",
                    "Ba2+": "INSOLUBLE (BaSO4) — very important in qualitative analysis",
                    "Pb2+": "INSOLUBLE (PbSO4)",
                    "Ca2+": "SLIGHTLY SOLUBLE (CaSO4, Ksp≈5×10⁻⁵)",
                    "Ra2+": "INSOLUBLE (RaSO4)",
                    "Ag2+": "SLIGHTLY SOLUBLE (Ag2SO4)",
                    "Hg2^2+": "SLIGHTLY/INSOLUBLE",
                    "Hg2+": "SLIGHTLY/INSOLUBLE",
                },
                "slightly_soluble": {
                    "Ca2+": "CaSO4 · 2H2O (gypsum) — slightly soluble (~0.02 M)",
                    "Ag2+": "Ag2SO4 moderately soluble",
                },
                "notes": "Al2(SO4)3, Cr2(SO4)3, Fe2(SO4)3 hydrolyze in water (acidic solutions).",
            },
            {
                "group": "Hydroxides (OH-)",
                "anions": ["OH-"],
                "general_rule": "Most hydroxide salts are INSOLUBLE.",
                "solubility": "insoluble",
                "exceptions": {
                    "Group 1 (Li+, Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE (decomposes on heating)",
                    "Ba2+": "MODERATELY SOLUBLE (Ba(OH)2)",
                    "Sr2+": "SLIGHTLY SOLUBLE (Sr(OH)2)",
                    "Ca2+": "SLIGHTLY SOLUBLE (Ca(OH)2)",
                    "Tl+": "soluble",
                },
                "slightly_soluble": {
                    "Ca(OH)2, Sr(OH)2, Ba(OH)2": "slightly to moderately soluble (strong bases from limited solubility)",
                },
                "notes": "Al(OH)3, Zn(OH)2, Cr(OH)3, Pb(OH)2, Sn(OH)2 are AMPHOTERIC (dissolve in excess strong base). Group 2 hydroxides become more soluble down the group.",
            },
            {
                "group": "Sulfides (S^2-)",
                "anions": ["S^2-", "S2-", "S(2-)"],
                "general_rule": "Most sulfide salts are INSOLUBLE.",
                "solubility": "insoluble",
                "exceptions": {
                    "Group 1 (Li+, Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE",
                    "Mg2+": "SLIGHTLY SOLUBLE (MgS hydrolyzes)",
                    "Ca2+": "SOLUBLE (alkaline earth sulfides react with water)",
                    "Sr2+": "SOLUBLE (alkaline earth sulfides react with water)",
                    "Ba2+": "SOLUBLE (alkaline earth sulfides react with water)",
                    "Mn2+": "SLIGHTLY SOLUBLE (easily dissolved by dilute acid)",
                },
                "slightly_soluble": {
                    "MnS": "Dissolves even in acetic acid (very slightly soluble)",
                },
                "notes": "Acid solubility hierarchy: MnS > FeS/ZnS > CdS > PbS/CuS > HgS (most insoluble). Used extensively in qualitative analysis group separation.",
            },
            {
                "group": "Carbonates (CO3^2-)",
                "anions": ["CO3^2-", "CO3(2-)"],
                "general_rule": "Most carbonate salts are INSOLUBLE.",
                "solubility": "insoluble",
                "exceptions": {
                    "Group 1 (Li+, Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE (decomposes on heating)",
                },
                "slightly_soluble": {},
                "notes": "All carbonates react with acid: CO3^2- + 2H+ → CO2↑ + H2O. Used as gas-evolving test in qualitative analysis.",
            },
            {
                "group": "Phosphates (PO4^3-)",
                "anions": ["PO4^3-", "PO4(3-)"],
                "general_rule": "Most phosphate salts are INSOLUBLE.",
                "solubility": "insoluble",
                "exceptions": {
                    "Group 1 (Li+, Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE",
                },
                "slightly_soluble": {},
                "notes": "Most phosphates are insoluble. Important biologically (bones: Ca3(PO4)2, DNA/RNA backbone).",
            },
            {
                "group": "Chromates (CrO4^2-)",
                "anions": ["CrO4^2-", "CrO4(2-)"],
                "general_rule": "Most chromate salts are INSOLUBLE.",
                "solubility": "insoluble",
                "exceptions": {
                    "Group 1 (Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE",
                    "Mg2+": "SOLUBLE (MgCrO4)",
                    "Ca2+": "SLIGHTLY SOLUBLE",
                },
                "slightly_soluble": {
                    "CaCrO4": "slightly soluble",
                },
                "notes": "BaCrO4 (yellow), PbCrO4 (chrome yellow) are important insoluble chromates used as pigments and indicators.",
            },
            {
                "group": "Oxalates (C2O4^2-)",
                "anions": ["C2O4^2-", "C2O4(2-)"],
                "general_rule": "Most oxalate salts are INSOLUBLE.",
                "solubility": "insoluble",
                "exceptions": {
                    "Group 1 (Na+, K+, etc.)": "SOLUBLE",
                    "NH4+": "SOLUBLE",
                    "Mg2+": "SOLUBLE",
                    "Fe3+": "forms soluble complex [Fe(C2O4)3]^3-",
                },
                "slightly_soluble": {},
                "notes": "CaC2O4 is the major component of kidney stones. Used in redox titrations with KMnO4.",
            },
        ]

        # Compound-specific database for quick lookup
        self._compound_db = {
            # Nitrates (all soluble)
            "NaNO3": ("soluble", "All nitrates are soluble."), "KNO3": ("soluble", "All nitrates are soluble."),
            "Ca(NO3)2": ("soluble", "All nitrates are soluble."), "Ba(NO3)2": ("soluble", "All nitrates are soluble."),
            "AgNO3": ("soluble", "All nitrates are soluble."), "Pb(NO3)2": ("soluble", "All nitrates are soluble."),
            "Al(NO3)3": ("soluble", "All nitrates are soluble."), "Fe(NO3)3": ("soluble", "All nitrates are soluble."),
            "NH4NO3": ("soluble", "All ammonium compounds are soluble."),
            # Acetates (all soluble)
            "CH3COONa": ("soluble", "All acetates are soluble."), "CH3COOK": ("soluble", "All acetates are soluble."),
            "CH3COONH4": ("soluble", "All acetates & ammonium compounds are soluble."),
            "Pb(CH3COO)2": ("soluble", "All acetates are soluble."),
            # Halides
            "NaCl": ("soluble", "Group 1 compounds are soluble."), "NaBr": ("soluble", "Group 1 compounds are soluble."),
            "NaI": ("soluble", "Group 1 compounds are soluble."), "KCl": ("soluble", "Group 1 compounds are soluble."),
            "KBr": ("soluble", "Group 1 compounds are soluble."), "KI": ("soluble", "Group 1 compounds are soluble."),
            "NH4Cl": ("soluble", "Ammonium compounds are soluble."), "NH4Br": ("soluble", "Ammonium compounds are soluble."),
            "CaCl2": ("soluble", "Ca2+ not an exception for chlorides."), "CaBr2": ("soluble", "Ca2+ not an exception for bromides."),
            "BaCl2": ("soluble", "Ba2+ not an exception for chlorides."), "AlCl3": ("soluble", "Al3+ not an exception for halides."),
            "FeCl3": ("soluble", "Fe3+ not an exception for halides."), "CuCl2": ("soluble", "Cu2+ not an exception for halides."),
            "ZnCl2": ("soluble", "Zn2+ not an exception for halides."), "MgCl2": ("soluble", "Mg2+ not an exception for halides."),
            "AgCl": ("insoluble", "Ag+ is an exception: AgCl, AgBr, AgI are insoluble."),
            "AgBr": ("insoluble", "Ag+ is an exception: AgCl, AgBr, AgI are insoluble."),
            "AgI": ("insoluble", "Ag+ is an exception: AgCl, AgBr, AgI are insoluble."),
            "AgF": ("soluble", "Exception to the exception: AgF IS soluble (fluorides differ)."),
            "PbCl2": ("slightly soluble", "Pb2+ is an exception: PbCl2 is moderately cold-soluble."),
            "PbI2": ("insoluble", "Pb2+ is an exception: PbI2 is insoluble (bright yellow)."),
            "Hg2Cl2": ("insoluble", "Hg2^2+ is an exception: calomel is insoluble."),
            "CuCl": ("insoluble", "Cu(I) halides are insoluble."),
            # Sulfates
            "Na2SO4": ("soluble", "Group 1 compounds are soluble."), "K2SO4": ("soluble", "Group 1 compounds are soluble."),
            "(NH4)2SO4": ("soluble", "Ammonium compounds are soluble."),
            "MgSO4": ("soluble", "Mg2+ not a sulfate exception."), "Al2(SO4)3": ("soluble", "Al3+ not a sulfate exception."),
            "FeSO4": ("soluble", "Fe2+ not a sulfate exception."), "Fe2(SO4)3": ("soluble", "Fe3+ not a sulfate exception."),
            "CuSO4": ("soluble", "Cu2+ not a sulfate exception."), "ZnSO4": ("soluble", "Zn2+ not a sulfate exception."),
            "Ag2SO4": ("slightly soluble", "Ag2+ is a slight exception for sulfates."),
            "CaSO4": ("slightly soluble", "Ca2+ is a slight exception for sulfates (gypsum)."),
            "BaSO4": ("insoluble", "Ba2+ is a key exception for sulfates."),
            "PbSO4": ("insoluble", "Pb2+ is an exception for sulfates."),
            "SrSO4": ("insoluble", "Sr2+ is an exception for sulfates."),
            # Hydroxides
            "NaOH": ("soluble", "Group 1 hydroxides are soluble."), "KOH": ("soluble", "Group 1 hydroxides are soluble."),
            "Ba(OH)2": ("moderately soluble", "Ba(OH)2 is sufficiently soluble to be a strong base."),
            "Ca(OH)2": ("slightly soluble", "Ca(OH)2 is slightly soluble but enough for limewater."),
            "Sr(OH)2": ("slightly soluble", "Sr(OH)2 is slightly soluble."),
            "Mg(OH)2": ("insoluble", "Most hydroxides are insoluble; Mg(OH)2 is milk of magnesia."),
            "Al(OH)3": ("insoluble", "Most hydroxides are insoluble; Al(OH)3 is amphoteric."),
            "Zn(OH)2": ("insoluble", "Most hydroxides are insoluble; Zn(OH)2 is amphoteric."),
            "Fe(OH)2": ("insoluble", "Most hydroxides are insoluble."), "Fe(OH)3": ("insoluble", "Very insoluble gelatinous ppt."),
            "Cu(OH)2": ("insoluble", "Most hydroxides are insoluble."), "Cr(OH)3": ("insoluble", "Amphoteric hydroxide."),
            "Pb(OH)2": ("insoluble", "Weakly amphoteric."), "Ni(OH)2": ("insoluble",), "Co(OH)2": ("insoluble",),
            "Mn(OH)2": ("insoluble",), "NH3OH (or NH4OH)": ("N/A", "Not a standard compound; NH3 + H2O ⇌ NH4+ + OH-"),
            # Carbonates
            "Na2CO3": ("soluble", "Group 1 carbonates are soluble."), "K2CO3": ("soluble", "Group 1 carbonates are soluble."),
            "(NH4)2CO3": ("soluble", "Ammonium carbonate is soluble."),
            "CaCO3": ("insoluble", "Most carbonates are insoluble (limestone, marble, shells)."),
            "BaCO3": ("insoluble", "Most carbonates are insoluble."), "MgCO3": ("insoluble", "Most carbonates are insoluble."),
            "PbCO3": ("insoluble", "Most carbonates are insoluble."), "Ag2CO3": ("insoluble", "Most carbonates are insoluble."),
            "FeCO3": ("insoluble", "Most carbonates are insoluble."), "ZnCO3": ("insoluble", "Most carbonates are insoluble."),
            "Al2(CO3)3": ("N/A", "Does not exist in water (hydrolyzes to Al(OH)3 + CO2)."),
            # Phosphates
            "Na3PO4": ("soluble", "Group 1 phosphates are soluble."), "K3PO4": ("soluble", "Group 1 phosphates are soluble."),
            "(NH4)3PO4": ("soluble", "Ammonium phosphate is soluble."),
            "Ca3(PO4)2": ("insoluble", "Most phosphates are insoluble (bone mineral)."),
            "Ag3PO4": ("insoluble", "Most phosphates are insoluble (yellow ppt)."),
            # Sulfides
            "Na2S": ("soluble", "Group 1 sulfides are soluble."), "K2S": ("soluble", "Group 1 sulfides are soluble."),
            "(NH4)2S": ("soluble", "Ammonium sulfide is soluble."),
            "ZnS": ("insoluble", "Most sulfides are insoluble."), "CuS": ("insoluble", "Very insoluble sulfide."),
            "CdS": ("insoluble", "Most sulfides are insoluble (yellow)."), "PbS": ("insoluble", "Very insoluble (galena)."),
            "FeS": ("insoluble", "Most sulfides are insoluble."), "HgS": ("insoluble", "Most insoluble common sulfide."),
            "Ag2S": ("insoluble", "Extremely insoluble."), "MnS": ("slightly soluble", "Dissolves in dilute acids."),
            "CaS": ("soluble", "Alkaline earth sulfides react with water (hydrolysis)."),
            # Chromates
            "Na2CrO4": ("soluble", "Group 1 chromates are soluble."), "K2CrO4": ("soluble", "Group 1 chromates are soluble."),
            "BaCrO4": ("insoluble", "Most chromates are insoluble (yellow ppt)."),
            "PbCrO4": ("insoluble", "Most chromates are insoluble (chrome yellow pigment)."),
            "Ag2CrO4": ("insoluble brick-red", "Used as indicator in Mohr titration."),
            "SrCrO4": ("insoluble", "Most chromates are insoluble."),
            # Oxalates
            "Na2C2O4": ("soluble", "Group 1 oxalates are soluble."), "K2C2O4": ("soluble", "Group 1 oxalates are soluble."),
            "CaC2O4": ("insoluble", "Most oxalates are insoluble (kidney stone component)."),
            # Fluorides
            "NaF": ("soluble", "Group 1 fluorides are soluble."), "KF": ("soluble", "Group 1 fluorides are soluble."),
            "NH4F": ("soluble", "Ammonium fluoride is soluble."),
            "CaF2": ("insoluble", "Most fluorides are insoluble (fluorite, Ksp=3.45×10⁻¹¹)."),
            "MgF2": ("insoluble", "Most fluorides are insoluble."), "PbF2": ("insoluble", "Most fluorides are insoluble."),
            "AgF": ("soluble", "AgF is an exception — it IS soluble unlike AgCl/AgBr/AgI."),
        }

        # Cation/anion parser helpers
        self._common_cations = {
            "na+": "Na+", "k+": "K+", "li+": "Li+", "rb+": "Rb+", "cs+": "Cs+",
            "nh4+": "NH4+", "ammonium": "NH4+",
            "mg2+": "Mg2+", "ca2+": "Ca2+", "sr2+": "Sr2+", "ba2+": "Ba2+",
            "al3+": "Al3+", "zn2+": "Zn2+", "fe2+": "Fe2+", "fe3+": "Fe3+",
            "cu2+": "Cu2+", "cu+": "Cu+", "ag+": "Ag+", "pb2+": "Pb2+",
            "hg2^2+": "Hg2^2+", "hg22+": "Hg2^2+", "hg2+": "Hg2+",
            "ni2+": "Ni2+", "co2+": "Co2+", "mn2+": "Mn2+", "cr3+": "Cr3+",
            "cd2+": "Cd2+", "sn2+": "Sn2+", "sn4+": "Sn4+", "bi3+": "Bi3+",
            "tl+": "Tl+", "ra2+": "Ra2+",
        }
        self._common_anions = {
            "no3-": "NO3-", "nitrate": "NO3-",
            "ch3coo-": "CH3COO-", "acetate": "CH3COO-", "c2h3o2-": "CH3COO-",
            "clo3-": "ClO3-", "clo4-": "ClO4-", "perchlorate": "ClO4-",
            "cl-": "Cl-", "br-": "Br-", "i-": "I-", "f-": "F-",
            "so4^2-": "SO4^2-", "so4(2-)": "SO4^2-", "sulfate": "SO4^2-",
            "oh-": "OH-", "hydroxide": "OH-",
            "s^2-": "S^2-", "s2-": "S^2-", "sulfide": "S^2-",
            "co3^2-": "CO3^2-", "co3(2-)": "CO3^2-", "carbonate": "CO3^2-",
            "po4^3-": "PO4^3-", "po4(3-)": "PO4^3-", "phosphate": "PO4^3-",
            "cro4^2-": "CrO4^2-", "cro4(2-)": "CrO4^2-", "chromate": "CrO4^2-",
            "c2o4^2-": "C2O4^2-", "c2o4(2-)": "C2O4^2-", "oxalate": "C2O4^2-",
        }

    def _run_base(self, compound: Optional[str] = None,
                  cation: Optional[str] = None,
                  anion: Optional[str] = None) -> dict:
        """Core logic: query solubility rules."""
        # If no compound provided, return all rules
        if not compound and not cation and not anion:
            return self._return_all_rules()

        # Parse compound into cation + anion if needed
        if compound:
            cat, ani = self._parse_compound(compound)
            if cat and ani:
                cation = cat
                anion = ani

        if not cation or not anion:
            raise ChemMCPError(
                f"Cannot determine cation/anion from '{compound}'. "
                f"Provide both cation and anion explicitly."
            )

        # Resolve names
        cat_key = self._resolve_cation(cation)
        ani_key = self._resolve_anion(anion)

        # Check direct compound lookup first
        comp_key = compound.strip() if compound else f"{cat_key.replace('+','')}{ani_key}"
        if comp_key in self._compound_db:
            solub, rule = self._compound_db[comp_key]
            exc_list = self._get_exceptions_for(cat_key, ani_key)
            notes = self._get_notes(comp_key)
            return self._build_result(
                solub, rule, exc_list, comp_key, cat_key, ani_key, notes
            )

        # Apply rule-based determination
        solub, rule_text, exc_list, notes = self._apply_rules(cat_key, ani_key)

        return self._build_result(
            solub, rule_text, exc_list,
            compound or f"{cat_key}({ani_key})",
            cat_key, ani_key, notes
        )

    def _build_result(self, solub, rule, exc_list, comp, cat, ani, notes):
        """Build standardized output."""
        detailed_rules = []
        for r in self._rules:
            entry = {
                "group": r["group"],
                "rule": r["general_rule"],
                "exceptions": list(r["exceptions"].keys()) if r.get("exceptions") else [],
            }
            detailed_rules.append(entry)

        return {
            "solubility": solub,
            "rule_applied": rule,
            "exceptions": exc_list,
            "compound_analyzed": comp,
            "cation": cat,
            "anion": ani,
            "detailed_rules": detailed_rules,
            "notes": notes,
        }

    def _return_all_rules(self) -> dict:
        """Return complete solubility rule reference."""
        detailed = []
        for r in self._rules:
            entry = {
                "group": r["group"],
                "general_rule": r["general_rule"],
                "typical_solubility": r["solubility"],
                "exceptions": r.get("exceptions", {}),
                "slightly_soluble": r.get("slightly_soluble", {}),
                "notes": r.get("notes", ""),
            }
            detailed.append(entry)

        summary_lines = [
            "=== COMPLETE SOLUBILITY RULES REFERENCE ===",
            "",
            "ALWAYS SOLUBLE:",
            "  • All Group 1 (alkali metal) compounds: Li+, Na+, K+, Rb+, Cs+",
            "  • All ammonium (NH4+) compounds",
            "  • All nitrates (NO3-), acetates (CH3COO-), chlorates (ClO3-), perchlorates (ClO4-)",
            "",
            "CONDITIONALLY SOLUBULE — see detailed_rules for each group:",
        ]
        for r in self._rules:
            summary_lines.append(f"  • {r['group']}: {r['general_rule']}")

        return {
            "solubility": "N/A (reference mode)",
            "rule_applied": "Full rule reference requested",
            "exceptions": [],
            "compound_analyzed": "ALL COMPOUNDS (reference)",
            "cation": "N/A",
            "anion": "N/A",
            "detailed_rules": detailed,
            "notes": "\n".join(summary_lines),
        }

    def _parse_compound(self, compound: str) -> tuple:
        """Parse compound formula into (cation, anion) using regex and lookup tables."""
        import re
        c = compound.strip()

        # Try direct lookup first
        if c in self._compound_db:
            pass  # will try to extract from name

        # Case-insensitive search for cation patterns
        best_cat, best_ani = None, None
        c_lower = c.lower()

        # Sort by length descending to match longer names first (e.g., 'NH4+' before 'Na+')
        cat_items = sorted(self._common_cations.items(), key=lambda x: len(x[0]), reverse=True)
        ani_items = sorted(self._common_anions.items(), key=lambda x: len(x[0]), reverse=True)

        for alias, key in cat_items:
            clean_alias = alias.replace("+", "").replace("^", "").replace("2", "")
            if clean_alias in c_lower or alias.lower() in c_lower:
                best_cat = key
                break

        for alias, key in ani_items:
            clean_alias = alias.replace("-", "").replace("^", "").replace("(", "").replace(")", "").replace("2", "").replace("3", "").replace("4", "")
            if clean_alias in c_lower or alias.lower() in c_lower:
                best_ani = key
                break

        return best_cat, best_ani

    def _apply_rules(self, cat_key: str, ani_key: str) -> tuple:
        """Apply solubility rules to determine solubility."""
        # Rule 1: Always-soluble cations
        if cat_key in ("Li+", "Na+", "K+", "Rb+", "Cs+", "NH4+"):
            exc = []
            note = f"{cat_key} compounds are always soluble (with any common anion)."
            return "soluble", f"Group 1 / NH4+ rule: {cat_key} compounds are always soluble.", exc, note

        # Rule 2: Always-soluble anions
        if ani_key in ("NO3-", "CH3COO-", "ClO3-", "ClO4-"):
            aname = {"NO3-": "nitrate", "CH3COO-": "acetate", "ClO3-": "chlorate", "ClO4-": "perchlorate"}
            return "soluble", f"All {aname.get(ani_key, ani_key)} salts are soluble.", [], ""

        # Rule 3: Apply conditional rules
        for rule in self._rules:
            if ani_key in rule["anions"]:
                exceptions = rule.get("exceptions", {})
                slightly = rule.get("slightly_soluble", {})

                # Check if cation is an exception
                is_exception = False
                exc_desc = None
                for exc_key, exc_desc in exceptions.items():
                    if cat_key in exc_key or exc_key == cat_key:
                        is_exception = True
                        break

                if is_exception:
                    # Check if it's in slightly soluble instead
                    for sl_key in (slightly or {}):
                        if cat_key in sl_key:
                            return (
                                "slightly soluble",
                                f"{rule['general_rule']}. {cat_key} is a SLIGHTLY SOLUBILITY exception.",
                                [f"{cat_key}: {exc_desc}"],
                                rule.get("notes", ""),
                            )
                    return (
                        "insoluble",
                        f"{rule['general_rule']}. {cat_key} is a listed EXCEPTION.",
                        [f"{cat_key}: {exc_desc}"],
                        rule.get("notes", ""),
                    )
                else:
                    return (
                        rule["solubility"],
                        rule["general_rule"],
                        [],
                        rule.get("notes", ""),
                    )

        # Fallback: unknown combination
        return (
            "unknown",
            f"No specific rule found for {cat_key} + {ani_key}. Consult specialized references.",
            [],
            "This combination may require experimental verification or specialized knowledge.",
        )

    def _get_exceptions_for(self, cat_key: str, ani_key: str) -> list:
        """Get exception descriptions for a cation-anion pair."""
        for rule in self._rules:
            if ani_key in rule["anions"]:
                excs = []
                for exc_key, desc in rule.get("exceptions", {}).items():
                    if cat_key in exc_key or exc_key == cat_key:
                        excs.append(desc)
                return excs
        return []

    def _get_notes(self, compound: str) -> str:
        """Get compound-specific notes."""
        notes_map = {
            "BaSO4": "Used in 'barium meal' X-ray imaging due to extreme insolubility. Ksp=1.08×10⁻¹⁰.",
            "AgCl": "Light-sensitive white precipitate (used in early photography). Curdles/darkens on exposure. Ksp=1.77×10⁻¹⁰.",
            "CaCO3": "Limestone, marble, seashells. Dissolves in acid with CO₂ effervescence. Ksp=3.36×10⁻⁹.",
            "PbI2": "Bright yellow precipitate ('gold rain' demo). Ksp=9.8×10⁻⁹.",
            "HgS": "Most insoluble common salt (Ksp=4×10⁻⁵³). Only aqua regia can dissolve it.",
            "CaF2": "Fluorite mineral. Very difficult to dissolve (requires conc. H2SO4 or boric acid fusion).",
            "BaCrO4": "Yellow precipitate. Insoluble in acetic acid (distinguishes from PbCrO4).",
            "PbCrO4": "Chrome yellow pigment. Used in pigments and as a pH indicator.",
            "Ag2CrO4": "Brick-red precipitate. Used as indicator in Mohr method (argentometric) titration.",
            "CaC2O4": "Main component of kidney stones. Redox-titratable with permanganate.",
            "Ca(OH)2": "Limewater. Slightly soluble but enough to be strongly basic. Used in CO2 test.",
            "Al(OH)3": "Gelatinous white precipitate. AMPHOTERIC — dissolves in both acid AND excess base.",
            "Mg(OH)2": "Milk of magnesia (antacid). White gelatinous precipitate.",
            "CuS": "Black precipitate. Extremely insoluble — requires hot concentrated HNO3 to dissolve.",
            "MnS": "Flesh-colored precipitate. Uniquely dissolves in even CH3COOH (acetic acid).",
        }
        return notes_map.get(compound, "")

    def _resolve_cation(self, name: str) -> str:
        n = name.strip()
        if n in self._common_cations.values():
            return n
        nl = n.lower().replace(" ", "")
        if nl in self._common_cations:
            return self._common_cations[nl]
        for k, v in self._common_cations.items():
            if k.replace("+","").replace("^","") == nl:
                return v
        return n

    def _resolve_anion(self, name: str) -> str:
        n = name.strip()
        if n in self._common_anions.values():
            return n
        nl = n.lower().replace(" ", "")
        if nl in self._common_anions:
            return self._common_anions[nl]
        for k, v in self._common_anions.items():
            if k.replace("-","").replace("^","").replace("(","").replace(")","") == nl:
                return v
        return n

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        s = input_str.strip()
        if not s:
            return self._run_base()
        return self._run_base(compound=s)
