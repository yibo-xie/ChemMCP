import logging
import math
from typing import Dict, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class WillPrecipitate(BaseTool):
    """
    判断是否生成沉淀（Qsp vs Ksp）。
    针对单一沉淀物，精确计算离子积Qsp并与Ksp比较。
    """
    __version__ = "0.1.0"
    name = "WillPrecipitate"
    func_name = "will_precipitate"
    description = "Determine whether a precipitate will form by comparing ion product Qsp with solubility product Ksp for a single compound. Returns Qsp, Ksp, ratio, verdict, and detailed explanation."
    implementation_description = "Calculates Qsp from given ion concentrations, compares with Ksp from built-in database. Qsp > Ksp → precipitate; Qsp ≈ Ksp → saturated; Qsp < Ksp → no precipitate."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Precipitation", "Qsp", "Ksp", "Ion Product", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "N/A", "Chemical formula of potential precipitate, e.g., 'AgCl', 'BaSO4', 'CaCO3'."),
        ("ion_concentrations", "dict", "N/A", "Dictionary of ion concentrations in mol/L, e.g., {'Ag+': 0.01, 'Cl-': 0.01} or {'Ba2+': 0.001, 'SO4^2-': 0.001}."),
        ("ksp_override", "float", "None", "Optional: override Ksp value. If not provided, uses database."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Format: 'compound;ion1=conc1,ion2=conc2'. E.g., 'AgCl;Ag+=0.01,Cl-=0.01' or 'BaSO4;Ba2+=0.001,SO4^2-=0.0005'."),
    ]

    output_sig = [
        ("compound", "str", "Compound formula."),
        ("qsp", "float", "Calculated ion product Qsp."),
        ("ksp", "float", "Solubility product Ksp used."),
        ("qsp_ksp_ratio", "float", "Qsp/Ksp ratio."),
        ("verdict", "str", "'precipitate', 'no_precipitate', or 'saturated'."),
        ("explanation", "str", "Detailed explanation of the result."),
    ]

    examples = [
        {
            "code_input": {
                "compound": "AgCl",
                "ion_concentrations": {"Ag+": 0.01, "Cl-": 0.01},
                "ksp_override": None,
            },
            "text_input": {
                "input_str": "AgCl;Ag+=0.01,Cl-=0.01",
            },
            "output": {
                "compound": "AgCl",
                "qsp": 1e-4,
                "ksp": 1.77e-10,
                "qsp_ksp_ratio": 564971.75,
                "verdict": "precipitate",
                "explanation": "Qsp=[Ag+][Cl-]=1×10⁻⁴ >> Ksp=1.77×10⁻¹⁰ (ratio > 5.6×10⁵). White AgCl precipitate WILL form.",
            },
        },
        {
            "code_input": {
                "compound": "BaSO4",
                "ion_concentrations": {"Ba2+": 1e-6, "SO4^2-": 1e-6},
                "ksp_override": None,
            },
            "text_input": {
                "input_str": "BaSO4;Ba2+=1e-6,SO4^2-=1e-6",
            },
            "output": {
                "compound": "BaSO4",
                "qsp": 1e-12,
                "ksp": 1.08e-10,
                "qsp_ksp_ratio": 0.0093,
                "verdict": "no_precipitate",
                "explanation": "Qsp=[Ba²⁺][SO₄²⁻]=1×10⁻¹² < Ksp=1.08×10⁻¹⁰ (ratio=0.009). Solution is unsaturated, NO precipitate.",
            },
        },
        {
            "code_input": {
                "compound": "PbI2",
                "ion_concentrations": {"Pb2+": 0.01, "I-": 0.02},
                "ksp_override": None,
            },
            "text_input": {
                "input_str": "PbI2;Pb2+=0.01,I-=0.02",
            },
            "output": {
                "compound": "PbI2",
                "qsp": 4e-6,
                "ksp": 9.8e-9,
                "qsp_ksp_ratio": 408.16,
                "verdict": "precipitate",
                "explanation": "Qsp=[Pb²⁺][I-]²=(0.01)(0.02)²=4×10⁻⁶ >> Ksp=9.8×10⁻⁹. Yellow PbI₂ precipitate WILL form.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Ksp database."""
        self._ksp_db = {
            # Halides
            "AgCl":     {"Ksp": 1.77e-10,   "cat": "Ag+",  "ani": "Cl-",   "m": 1, "n": 1},
            "AgBr":     {"Ksp": 5.35e-13,   "cat": "Ag+",  "ani": "Br-",   "m": 1, "n": 1},
            "AgI":      {"Ksp": 8.52e-17,   "cat": "Ag+",  "ani": "I-",    "m": 1, "n": 1},
            "PbCl2":    {"Ksp": 1.7e-5,     "cat": "Pb2+", "ani": "Cl-",   "m": 1, "n": 2},
            "PbI2":     {"Ksp": 9.8e-9,     "cat": "Pb2+", "ani": "I-",    "m": 1, "n": 2},
            "Hg2Cl2":   {"Ksp": 1.43e-18,   "cat": "Hg2^2+","ani": "Cl-",  "m": 1, "n": 2},
            "CuCl":     {"Ksp": 1.72e-7,    "cat": "Cu+",  "ani": "Cl-",   "m": 1, "n": 1},
            "CuI":      {"Ksp": 1.27e-12,   "cat": "Cu+",  "ani": "I-",    "m": 1, "n": 1},
            # Sulfates
            "BaSO4":    {"Ksp": 1.08e-10,   "cat": "Ba2+", "ani": "SO4^2-", "m": 1, "n": 1},
            "PbSO4":    {"Ksp": 2.53e-8,    "cat": "Pb2+", "ani": "SO4^2-", "m": 1, "n": 1},
            "CaSO4":    {"Ksp": 4.93e-5,    "cat": "Ca2+", "ani": "SO4^2-", "m": 1, "n": 1},
            "SrSO4":    {"Ksp": 3.44e-7,    "cat": "Sr2+", "ani": "SO4^2-", "m": 1, "n": 1},
            "Ag2SO4":   {"Ksp": 1.20e-5,    "cat": "Ag+",  "ani": "SO4^2-", "m": 2, "n": 1},
            # Carbonates
            "CaCO3":    {"Ksp": 3.36e-9,    "cat": "Ca2+", "ani": "CO3^2-", "m": 1, "n": 1},
            "BaCO3":    {"Ksp": 2.58e-9,    "cat": "Ba2+", "ani": "CO3^2-", "m": 1, "n": 1},
            "MgCO3":    {"Ksp": 6.82e-6,    "cat": "Mg2+", "ani": "CO3^2-", "m": 1, "n": 1},
            "PbCO3":    {"Ksp": 7.40e-14,   "cat": "Pb2+", "ani": "CO3^2-", "m": 1, "n": 1},
            "Ag2CO3":   {"Ksp": 8.46e-12,   "cat": "Ag+",  "ani": "CO3^2-", "m": 2, "n": 1},
            "SrCO3":    {"Ksp": 5.60e-10,   "cat": "Sr2+", "ani": "CO3^2-", "m": 1, "n": 1},
            "ZnCO3":    {"Ksp": 1.46e-10,   "cat": "Zn2+", "ani": "CO3^2-", "m": 1, "n": 1},
            # Hydroxides
            "Fe(OH)2":  {"Ksp": 4.87e-17,   "cat": "Fe2+", "ani": "OH-",   "m": 1, "n": 2},
            "Fe(OH)3":  {"Ksp": 2.79e-39,   "cat": "Fe3+", "ani": "OH-",   "m": 1, "n": 3},
            "Cu(OH)2":  {"Ksp": 2.20e-20,   "cat": "Cu2+", "ani": "OH-",   "m": 1, "n": 2},
            "Mg(OH)2":  {"Ksp": 5.61e-12,   "cat": "Mg2+", "ani": "OH-",   "m": 1, "n": 2},
            "Ca(OH)2":  {"Ksp": 5.02e-6,    "cat": "Ca2+", "ani": "OH-",   "m": 1, "n": 2},
            "Al(OH)3":  {"Ksp": 3.28e-34,   "cat": "Al3+", "ani": "OH-",   "m": 1, "n": 3},
            "Zn(OH)2":  {"Ksp": 3.00e-17,   "cat": "Zn2+", "ani": "OH-",   "m": 1, "n": 2},
            "Pb(OH)2":  {"Ksp": 1.43e-20,   "cat": "Pb2+", "ani": "OH-",   "m": 1, "n": 2},
            "Cr(OH)3":  {"Ksp": 6.30e-31,   "cat": "Cr3+", "ani": "OH-",   "m": 1, "n": 3},
            "Mn(OH)2":  {"Ksp": 2.06e-13,   "cat": "Mn2+", "ani": "OH-",   "m": 1, "n": 2},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "cat": "Ni2+", "ani": "OH-",   "m": 1, "n": 2},
            "Co(OH)2":  {"Ksp": 5.92e-15,   "cat": "Co2+", "ani": "OH-",   "m": 1, "n": 2},
            # Phosphates
            "Ca3(PO4)2":{"Ksp": 2.07e-33,   "cat": "Ca2+", "ani": "PO4^3-", "m": 3, "n": 2},
            "Ag3PO4":  {"Ksp": 8.89e-17,   "cat": "Ag+",  "ani": "PO4^3-", "m": 3, "n": 1},
            # Sulfides
            "FeS":      {"Ksp": 6.30e-19,   "cat": "Fe2+", "ani": "S^2-",   "m": 1, "n": 1},
            "MnS":      {"Ksp": 3.00e-14,   "cat": "Mn2+", "ani": "S^2-",   "m": 1, "n": 1},
            "ZnS":      {"Ksp": 2.50e-22,   "cat": "Zn2+", "ani": "S^2-",   "m": 1, "n": 1},
            "CdS":      {"Ksp": 8.00e-27,   "cat": "Cd2+", "ani": "S^2-",   "m": 1, "n": 1},
            "PbS":      {"Ksp": 9.04e-29,   "cat": "Pb2+", "ani": "S^2-",   "m": 1, "n": 1},
            "CuS":      {"Ksp": 6.00e-36,   "cat": "Cu2+", "ani": "S^2-",   "m": 1, "n": 1},
            "HgS":      {"Ksp": 4.00e-53,   "cat": "Hg2+", "ani": "S^2-",   "m": 1, "n": 1},
            "Ag2S":     {"Ksp": 6.30e-50,   "cat": "Ag+",  "ani": "S^2-",   "m": 2, "n": 1},
            # Chromates
            "BaCrO4":  {"Ksp": 1.17e-10,   "cat": "Ba2+", "ani": "CrO4^2-", "m": 1, "n": 1},
            "PbCrO4":  {"Ksp": 2.8e-13,    "cat": "Pb2+", "ani": "CrO4^2-", "m": 1, "n": 1},
            "Ag2CrO4": {"Ksp": 1.12e-12,   "cat": "Ag+",  "ani": "CrO4^2-", "m": 2, "n": 1},
            # Fluorides/Others
            "CaF2":    {"Ksp": 3.45e-11,   "cat": "Ca2+", "ani": "F-",    "m": 1, "n": 2},
            "PbF2":    {"Ksp": 3.3e-8,     "cat": "Pb2+", "ani": "F-",    "m": 1, "n": 2},
            "CaC2O4":  {"Ksp": 2.32e-9,    "cat": "Ca2+", "ani": "C2O4^2-", "m": 1, "n": 1},
        }

    def _run_base(self, compound: str, ion_concentrations: Dict[str, float],
                  ksp_override: Optional[float] = None) -> dict:
        """Core logic: compare Qsp vs Ksp."""
        key = self._resolve(compound)
        if key not in self._ksp_db:
            raise ChemMCPError(
                f"Unknown compound '{compound}'. Available: {sorted(self._ksp_db.keys())[:25]}"
            )

        d = self._ksp_db[key]
        ksp = ksp_override if ksp_override is not None else d["Ksp"]
        m, n = d["m"], d["n"]

        # Get ion concentrations
        cat_conc = self._get_ion_conc(d["cat"], ion_concentrations)
        ani_conc = self._get_ion_conc(d["ani"], ion_concentrations)
        if cat_conc is None or ani_conc is None:
            missing = d["cat"] if cat_conc is None else d["ani"]
            raise ChemMCPError(f"Missing concentration for ion: {missing}")

        # Calculate Qsp = [cation]^m * [anion]^n
        qsp = (cat_conc ** m) * (ani_conc ** n)
        ratio = qsp / ksp if ksp > 0 else float('inf')

        # Determine verdict
        tolerance = 2.0
        if qsp > ksp * tolerance:
            verdict = "precipitate"
        elif qsp >= ksp / tolerance:
            verdict = "saturated"
        else:
            verdict = "no_precipitate"

        # Build Qsp expression
        qsp_expr = f"[{d['cat']}]^{m}[{d['ani']}]^{n} = ({cat_conc:.3e})^{m} × ({ani_conc:.3e})^{n} = {qsp:.3e}"

        explanation = (
            f"For {key}: {qsp_expr}\n"
            f"Ksp({key}) = {ksp:.3e}\n"
            f"Qsp/Ksp = {ratio:.4f} → {verdict.upper().replace('_', ' ')}\n"
        )

        if verdict == "precipitate":
            explanation += (
                f"Since Qsp ({qsp:.3e}) > Ksp ({ksp:.3e}), "
                f"the solution is supersaturated and {key} precipitate will form."
            )
        elif verdict == "saturated":
            explanation += (
                f"Qsp ≈ Ksp, the solution is at saturation equilibrium."
            )
        else:
            explanation += (
                f"Since Qsp ({qsp:.3e}) < Ksp ({ksp:.3e}), "
                f"the solution is unsaturated and no {key} precipitate will form."
            )

        logger.info(f"WillPrecipitate: {key} Qsp={qsp:.3e} vs Ksp={ksp:.3e} → {verdict}")
        return {
            "compound": key,
            "qsp": round(qsp, int(max(0, 3 - math.log10(max(abs(qsp), 1e-30))))),
            "ksp": ksp,
            "qsp_ksp_ratio": round(ratio, 4),
            "verdict": verdict,
            "explanation": explanation,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split(";")
            compound = parts[0].strip()
            conc_parts = parts[1].split(",") if len(parts) > 1 else []
            concentrations = {}
            for cp in conc_parts:
                cp = cp.strip()
                if "=" not in cp:
                    continue
                ion, val = cp.split("=", 1)
                concentrations[ion.strip()] = float(val.strip())
            return self._run_base(compound, concentrations)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'compound;ion1=conc1,ion2=conc2'")

    def _resolve(self, name: str) -> str:
        n = name.strip()
        if n in self._ksp_db:
            return n
        for k in self._ksp_db:
            if k.lower() == n.lower():
                return k
        return n

    def _get_ion_conc(self, ion: str, conc_dict: dict):
        if ion in conc_dict:
            return conc_dict[ion]
        for key, val in conc_dict.items():
            if key.replace("^", "").replace("_", "").lower() == ion.replace("^", "").replace("_", "").lower():
                return val
        return None
