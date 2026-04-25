import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PrecipitationPrediction(BaseTool):
    """
    预测沉淀生成（基于溶度积 Ksp）。
    通过比较离子积 Qsp 与 Ksp 判断是否会产生沉淀。
    """
    __version__ = "0.1.0"
    name = "PrecipitationPrediction"
    func_name = "predict_precipitation"
    description = "Predict whether a precipitate will form based on solubility product constant (Ksp). Compares ion product Qsp with Ksp to determine precipitation."
    implementation_description = "Uses built-in Ksp database of common salts at 25°C. Calculates Qsp from given ion concentrations and compares with Ksp: if Qsp > Ksp → precipitate forms; Qsp < Ksp → no precipitate; Qsp ≈ Ksp → saturated."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Precipitation", "Ksp", "Solubility Equilibrium", "Ionic Product"]
    required_envs = []

    code_input_sig = [
        ("cations", "list", "N/A", "List of cation species, e.g., ['Ag+', 'Ba2+', 'Pb2+']."),
        ("anions", "list", "N/A", "List of anion species, e.g., ['Cl-', 'SO4^2-', 'I-']."),
        ("concentrations", "dict", "N/A", "Dictionary of ion concentrations in mol/L, e.g., {'Ag+': 0.01, 'Cl-': 0.01}."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Semicolon-separated: cations;anions;concentrations. E.g., 'Ag+,Ba2+;Cl-,SO4^2-;Ag+=0.01,Cl-=0.001,Ba2+=0.1,SO4^2-=0.01'."),
    ]

    output_sig = [
        ("result", "dict", "Detailed precipitation analysis including will_precipitate, Qsp, Ksp, precipitate_formula, explanation."),
    ]

    examples = [
        {
            "code_input": {
                "cations": ["Ag+"],
                "anions": ["Cl-"],
                "concentrations": {"Ag+": 0.01, "Cl-": 0.01},
            },
            "text_input": {
                "input_str": "Ag+;Cl-;Ag+=0.01,Cl-=0.01",
            },
            "output": {
                "result": {
                    "will_precipitate": True,
                    "precipitates": [{"formula": "AgCl", "Qsp": 1e-4, "Ksp": 1.77e-10, "verdict": "Qsp >> Ksp, precipitate WILL form"}],
                    "explanation": "Mixing Ag+ (0.01 M) and Cl- (0.01 M) gives Qsp = [Ag+][Cl-] = 1×10⁻⁴, which is far greater than Ksp(AgCl) = 1.77×10⁻¹⁰. A white AgCl precipitate will form.",
                }
            },
        },
        {
            "code_input": {
                "cations": ["Ba2+"],
                "anions": ["SO4^2-"],
                "concentrations": {"Ba2+": 0.001, "SO4^2-": 0.001},
            },
            "text_input": {
                "input_str": "Ba2+;SO4^2-;Ba2+=0.001,SO4^2-=0.001",
            },
            "output": {
                "result": {
                    "will_precipitate": True,
                    "precipitates": [{"formula": "BaSO4", "Qsp": 1e-6, "Ksp": 1.08e-10, "verdict": "Qsp >> Ksp, precipitate WILL form"}],
                    "explanation": "Qsp = [Ba²⁺][SO₄²⁻] = 1×10⁻⁶ > Ksp(BaSO₄) = 1.08×10⁻¹⁰. White BaSO₄ precipitate will form.",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Ksp database for common salts at 25°C."""
        # Ksp values at 25°C (298 K)
        # Format: {precipitate_formula: {"Ksp": value, "cation": X+, "anion": Y-, "stoich": (m, n)}}
        # For salt: XmYn ⇌ mX^n+ + nY^m-, Ksp = [X^n+]^m * [Y^m-]^n
        self._ksp_db = {
            # Halides
            "AgCl":     {"Ksp": 1.77e-10,   "cation": "Ag+",  "anion": "Cl-",   "stoich": (1, 1)},
            "AgBr":     {"Ksp": 5.35e-13,   "cation": "Ag+",  "anion": "Br-",   "stoich": (1, 1)},
            "AgI":      {"Ksp": 8.52e-17,   "cation": "Ag+",  "anion": "I-",    "stoich": (1, 1)},
            "PbCl2":    {"Ksp": 1.7e-5,     "cation": "Pb2+", "anion": "Cl-",   "stoich": (1, 2)},
            "PbI2":     {"Ksp": 9.8e-9,     "cation": "Pb2+", "anion": "I-",    "stoich": (1, 2)},
            "Hg2Cl2":   {"Ksp": 1.43e-18,   "cation": "Hg2^2+", "anion": "Cl-", "stoich": (1, 2)},
            "CuCl":     {"Ksp": 1.72e-7,    "cation": "Cu+",  "anion": "Cl-",   "stoich": (1, 1)},
            "CuI":      {"Ksp": 1.27e-12,   "cation": "Cu+",  "anion": "I-",    "stoich": (1, 1)},
            # Sulfates
            "BaSO4":    {"Ksp": 1.08e-10,   "cation": "Ba2+", "anion": "SO4^2-", "stoich": (1, 1)},
            "PbSO4":    {"Ksp": 2.53e-8,    "cation": "Pb2+", "anion": "SO4^2-", "stoich": (1, 1)},
            "CaSO4":    {"Ksp": 4.93e-5,    "cation": "Ca2+", "anion": "SO4^2-", "stoich": (1, 1)},
            "SrSO4":    {"Ksp": 3.44e-7,    "cation": "Sr2+", "anion": "SO4^2-", "stoich": (1, 1)},
            "Ag2SO4":   {"Ksp": 1.20e-5,    "cation": "Ag+",  "anion": "SO4^2-", "stoich": (2, 1)},
            # Carbonates
            "CaCO3":    {"Ksp": 3.36e-9,    "cation": "Ca2+", "anion": "CO3^2-", "stoich": (1, 1)},
            "BaCO3":    {"Ksp": 2.58e-9,    "cation": "Ba2+", "anion": "CO3^2-", "stoich": (1, 1)},
            "MgCO3":    {"Ksp": 6.82e-6,    "cation": "Mg2+", "anion": "CO3^2-", "stoich": (1, 1)},
            "PbCO3":    {"Ksp": 7.40e-14,   "cation": "Pb2+", "anion": "CO3^2-", "stoich": (1, 1)},
            "Ag2CO3":   {"Ksp": 8.46e-12,   "cation": "Ag+",  "anion": "CO3^2-", "stoich": (2, 1)},
            "SrCO3":    {"Ksp": 5.60e-10,   "cation": "Sr2+", "anion": "CO3^2-", "stoich": (1, 1)},
            "ZnCO3":    {"Ksp": 1.46e-10,   "cation": "Zn2+", "anion": "CO3^2-", "stoich": (1, 1)},
            # Hydroxides
            "Fe(OH)2":  {"Ksp": 4.87e-17,   "cation": "Fe2+", "anion": "OH-",   "stoich": (1, 2)},
            "Fe(OH)3":  {"Ksp": 2.79e-39,   "cation": "Fe3+", "anion": "OH-",   "stoich": (1, 3)},
            "Cu(OH)2":  {"Ksp": 2.20e-20,   "cation": "Cu2+", "anion": "OH-",   "stoich": (1, 2)},
            "Mg(OH)2":  {"Ksp": 5.61e-12,   "cation": "Mg2+", "anion": "OH-",   "stoich": (1, 2)},
            "Ca(OH)2":  {"Ksp": 5.02e-6,    "cation": "Ca2+", "anion": "OH-",   "stoich": (1, 2)},
            "Al(OH)3":  {"Ksp": 3.28e-34,   "cation": "Al3+", "anion": "OH-",   "stoich": (1, 3)},
            "Zn(OH)2":  {"Ksp": 3.00e-17,   "cation": "Zn2+", "anion": "OH-",   "stoich": (1, 2)},
            "Pb(OH)2":  {"Ksp": 1.43e-20,   "cation": "Pb2+", "anion": "OH-",   "stoich": (1, 2)},
            "Cr(OH)3":  {"Ksp": 6.30e-31,   "cation": "Cr3+", "anion": "OH-",   "stoich": (1, 3)},
            "Mn(OH)2":  {"Ksp": 2.06e-13,   "cation": "Mn2+", "anion": "OH-",   "stoich": (1, 2)},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "cation": "Ni2+", "anion": "OH-",   "stoich": (1, 2)},
            "Co(OH)2":  {"Ksp": 5.92e-15,   "cation": "Co2+", "anion": "OH-",   "stoich": (1, 2)},
            # Phosphates
            "Ca3(PO4)2":{"Ksp": 2.07e-33,   "cation": "Ca2+", "anion": "PO4^3-", "stoich": (3, 2)},
            "Ag3PO4":  {"Ksp": 8.89e-17,   "cation": "Ag+",  "anion": "PO4^3-", "stoich": (3, 1)},
            # Sulfides
            "FeS":      {"Ksp": 6.30e-19,   "cation": "Fe2+", "anion": "S^2-",   "stoich": (1, 1)},
            "MnS":      {"Ksp": 3.00e-14,   "cation": "Mn2+", "anion": "S^2-",   "stoich": (1, 1)},
            "ZnS":      {"Ksp": 2.50e-22,   "cation": "Zn2+", "anion": "S^2-",   "stoich": (1, 1)},
            "CdS":      {"Ksp": 8.00e-27,   "cation": "Cd2+", "anion": "S^2-",   "stoich": (1, 1)},
            "PbS":      {"Ksp": 9.04e-29,   "cation": "Pb2+", "anion": "S^2-",   "stoich": (1, 1)},
            "CuS":      {"Ksp": 6.00e-36,   "cation": "Cu2+", "anion": "S^2-",   "stoich": (1, 1)},
            "HgS":      {"Ksp": 4.00e-53,   "cation": "Hg2+", "anion": "S^2-",   "stoich": (1, 1)},
            "Ag2S":     {"Ksp": 6.30e-50,   "cation": "Ag+",  "anion": "S^2-",   "stoich": (2, 1)},
            # Chromates
            "BaCrO4":  {"Ksp": 1.17e-10,   "cation": "Ba2+", "anion": "CrO4^2-", "stoich": (1, 1)},
            "PbCrO4":  {"Ksp": 2.8e-13,    "cation": "Pb2+", "anion": "CrO4^2-", "stoich": (1, 1)},
            "Ag2CrO4": {"Ksp": 1.12e-12,   "cation": "Ag+",  "anion": "CrO4^2-", "stoich": (2, 1)},
            # Others
            "AgSCN":   {"Ksp": 1.03e-12,   "cation": "Ag+",  "anion": "SCN-",  "stoich": (1, 1)},
            "CaF2":    {"Ksp": 3.45e-11,   "cation": "Ca2+", "anion": "F-",    "stoich": (1, 2)},
            "PbF2":    {"Ksp": 3.3e-8,     "cation": "Pb2+", "anion": "F-",    "stoich": (1, 2)},
            "CaC2O4":  {"Ksp": 2.32e-9,    "cation": "Ca2+", "anion": "C2O4^2-", "stoich": (1, 1)},
        }

        # Ion name normalization map
        self._ion_aliases = {
            "so4": "so4^2-", "sulfate": "so4^2-",
            "co3": "co3^2-", "carbonate": "co3^2-",
            "oh": "oh-", "hydroxide": "oh-",
            "po4": "po4^3-", "phosphate": "po4^3-",
            "cro4": "cro4^2-", "chromate": "cro4^2-",
            "no3": "no3-", "nitrate": "no3-",
            "s2": "s^2-", "sulfide": "s^2-",
            "scn": "scn-", "thiocyanate": "scn-",
            "c2o4": "c2o4^2-", "oxalate": "c2o4^2-",
        }

    def _run_base(self, cations: List[str], anions: List[str], concentrations: Dict[str, float]) -> dict:
        """Core logic: predict precipitation for all possible cation-anion pairs."""
        if not cations or not anions:
            raise ChemMCPError("Both cations and anions must be provided.")
        if not concentrations:
            raise ChemMCPError("Concentrations dictionary cannot be empty.")

        results = []
        any_precipitate = False

        for cat in cations:
            cat_norm = self._normalize_ion(cat)
            for an in anions:
                an_norm = self._normalize_ion(an)

                # Find matching precipitate in Ksp database
                match = self._find_precipitate(cat_norm, an_norm)
                if match is None:
                    continue

                formula, ksp_data = match

                # Get concentrations
                cat_conc = self._get_concentration(cat_norm, concentrations)
                an_conc = self._get_concentration(an_norm, concentrations)

                if cat_conc is None or an_conc is None:
                    results.append({
                        "formula": formula,
                        "status": "missing_concentration",
                        "message": f"Missing concentration for {cat_norm} or {an_norm}",
                    })
                    continue

                # Calculate Qsp
                m, n = ksp_data["stoich"]
                qsp = (cat_conc ** m) * (an_conc ** n)
                ksp = ksp_data["Ksp"]

                # Determine verdict
                tolerance = 1.5  # factor for "approximately equal"
                if qsp > ksp * tolerance:
                    verdict = "Qsp > Ksp, precipitate WILL form"
                    will_ppt = True
                    any_precipitate = True
                elif qsp >= ksp / tolerance:
                    verdict = "Qsp ≈ Ksp, solution is SATURATED (at equilibrium)"
                    will_ppt = True
                    any_precipitate = True
                else:
                    verdict = "Qsp < Ksp, NO precipitate (unsaturated)"
                    will_ppt = False

                results.append({
                    "formula": formula,
                    "cation": cat_norm,
                    "anion": an_norm,
                    "Qsp": qsp,
                    "Ksp": ksp,
                    "Qsp_Ksp_ratio": round(qsp / ksp, 4),
                    "will_precipitate": will_ppt,
                    "verdict": verdict,
                })

        # Build explanation
        if any_precipitate:
            ppt_list = [r["formula"] for r in results if r.get("will_precipitate")]
            explanation = (
                f"Precipitation predicted: {', '.join(ppt_list)}. "
                f"Out of {len(results)} possible salt(s), {sum(1 for r in results if r.get('will_precipitate'))} will precipitate."
            )
        elif results:
            explanation = (
                f"No precipitation expected for the given ion concentrations. "
                f"All {len(results)} checked salts have Qsp < Ksp."
            )
        else:
            explanation = "No matching salts found in Ksp database for the given ion combinations."

        return {
            "will_precipitate": any_precipitate,
            "precipitates": results,
            "explanation": explanation,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse semicolon-separated text input."""
        try:
            parts = input_str.strip().split(";")
            if len(parts) < 3:
                raise ValueError("Expected format: cations;anions;concentrations")

            cations = [c.strip() for c in parts[0].split(",") if c.strip()]
            anions = [a.strip() for a in parts[1].split(",") if a.strip()]
            conc_parts = parts[2].split(",")
            concentrations = {}
            for cp in conc_parts:
                cp = cp.strip()
                if "=" not in cp:
                    continue
                ion, val = cp.split("=", 1)
                concentrations[ion.strip()] = float(val.strip())

            return self._run_base(cations, anions, concentrations)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Expected format: 'cations;anions;conc_dict'")

    def _normalize_ion(self, ion: str) -> str:
        """Normalize ion name."""
        s = ion.strip().lower()
        return self._ion_aliases.get(s, ion.strip())

    def _find_precipitate(self, cation: str, anion: str) -> Optional[tuple]:
        """Find precipitate formula from Ksp database."""
        for formula, data in self._ksp_db.items():
            cat_match = self._ion_matches(data["cation"], cation)
            an_match = self._ion_matches(data["anion"], anion)
            if cat_match and an_match:
                return (formula, data)
        return None

    def _ion_matches(self, db_ion: str, query_ion: str) -> bool:
        """Check if two ion names refer to the same ion."""
        d = db_ion.lower().replace("^", "").replace("_", "")
        q = query_ion.lower().replace("^", "").replace("_", "")
        return d == q

    def _get_concentration(self, ion: str, conc_dict: Dict[str, float]) -> Optional[float]:
        """Get concentration for an ion, trying various name formats."""
        # Exact match
        if ion in conc_dict:
            return conc_dict[ion]
        # Case-insensitive
        for key, val in conc_dict.items():
            if key.lower() == ion.lower():
                return val
        # Normalized match
        norm_key = self._normalize_ion(key)
        if norm_key.lower() == ion.lower().replace("^", "").replace("_", ""):
            return val
        return None
