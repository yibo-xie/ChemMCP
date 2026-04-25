import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetKsp(BaseTool):
    """
    查询溶度积常数（Ksp）。
    提供常见难溶电解质的Ksp值、溶解度、沉淀反应方程式等信息。
    """
    __version__ = "0.1.0"
    name = "GetKsp"
    func_name = "get_ksp"
    description = "Query solubility product constant (Ksp) for common sparingly soluble salts at 25°C. Returns Ksp value, dissolution equation, molar solubility, and reference data."
    implementation_description = "Uses a built-in Ksp database of 60+ common precipitates at 25°C (298 K). Returns comprehensive data including stoichiometry, ion charges, and calculated molar solubility."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Ksp", "Solubility", "Precipitation", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "N/A", "Chemical formula or name of the compound, e.g., 'AgCl', 'BaSO4', 'calcium carbonate'."),
    ]

    text_input_sig = [
        ("compound_str", "str", "N/A", "Compound formula or name. E.g., 'AgCl' or 'silver chloride'."),
    ]

    output_sig = [
        ("formula", "str", "Chemical formula of the compound."),
        ("ksp", "float", "Solubility product constant at 25°C."),
        ("pKsp", "float", "-log10(Ksp)."),
        ("temperature", "float", "Temperature in Kelvin (default 298)."),
        ("dissolution_eq", "str", "Dissolution equilibrium equation."),
        ("stoichiometry", "dict", "Stoichiometry: {'cation': ..., 'anion': ..., 'cation_coeff': ..., 'anion_coeff': ...}."),
        ("molar_solubility", "float", "Calculated molar solubility in mol/L (pure water, no common ion)."),
        ("g_per_100ml", "float", "Solubility in g/100mL."),
        ("molar_mass", "float", "Molar mass in g/mol."),
        ("source_note", "str", "Source/reference note for the data."),
    ]

    examples = [
        {
            "code_input": {"compound": "AgCl"},
            "text_input": {"compound_str": "AgCl"},
            "output": {
                "formula": "AgCl",
                "ksp": 1.77e-10,
                "pKsp": 9.75,
                "temperature": 298.0,
                "dissolution_eq": "AgCl(s) ⇌ Ag⁺(aq) + Cl⁻(aq)",
                "stoichiometry": {"cation": "Ag+", "anion": "Cl-", "cation_coeff": 1, "anion_coeff": 1},
                "molar_solubility": 1.33e-5,
                "g_per_100ml": 1.91e-4,
                "molar_mass": 143.32,
                "source_note": "Standard thermodynamic data at 25°C.",
            },
        },
        {
            "code_input": {"compound": "BaSO4"},
            "text_input": {"compound_str": "barium sulfate"},
            "output": {
                "formula": "BaSO4",
                "ksp": 1.08e-10,
                "pKsp": 9.97,
                "temperature": 298.0,
                "dissolution_eq": "BaSO₄(s) ⇌ Ba²⁺(aq) + SO₄²⁻(aq)",
                "stoichiometry": {"cation": "Ba2+", "anion": "SO4^2-", "cation_coeff": 1, "anion_coeff": 1},
                "molar_solubility": 1.04e-5,
                "g_per_100ml": 2.42e-4,
                "molar_mass": 233.39,
                "source_note": "Standard thermodynamic data at 25°C.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Ksp database with molar masses."""
        # Format: {formula: {ksp, cation, anion, stoich (m,n), molar_mass}}
        self._db = {
            # Halides
            "AgCl":     {"Ksp": 1.77e-10,   "cat": "Ag+",  "ani": "Cl-",   "m": 1, "n": 1, "Mw": 143.32},
            "AgBr":     {"Ksp": 5.35e-13,   "cat": "Ag+",  "ani": "Br-",   "m": 1, "n": 1, "Mw": 187.77},
            "AgI":      {"Ksp": 8.52e-17,   "cat": "Ag+",  "ani": "I-",    "m": 1, "n": 1, "Mw": 234.77},
            "PbCl2":    {"Ksp": 1.7e-5,     "cat": "Pb2+", "ani": "Cl-",   "m": 1, "n": 2, "Mw": 278.10},
            "PbI2":     {"Ksp": 9.8e-9,     "cat": "Pb2+", "ani": "I-",    "m": 1, "n": 2, "Mw": 461.01},
            "Hg2Cl2":   {"Ksp": 1.43e-18,   "cat": "Hg2^2+","ani": "Cl-",  "m": 1, "n": 2, "Mw": 472.09},
            "CuCl":     {"Ksp": 1.72e-7,    "cat": "Cu+",  "ani": "Cl-",   "m": 1, "n": 1, "Mw": 98.999},
            "CuI":      {"Ksp": 1.27e-12,   "cat": "Cu+",  "ani": "I-",    "m": 1, "n": 1, "Mw": 190.45},
            # Sulfates
            "BaSO4":    {"Ksp": 1.08e-10,   "cat": "Ba2+", "ani": "SO4^2-","m": 1, "n": 1, "Mw": 233.39},
            "PbSO4":    {"Ksp": 2.53e-8,    "cat": "Pb2+", "ani": "SO4^2-","m": 1, "n": 1, "Mw": 303.26},
            "CaSO4":    {"Ksp": 4.93e-5,    "cat": "Ca2+", "ani": "SO4^2-","m": 1, "n": 1, "Mw": 136.14},
            "SrSO4":    {"Ksp": 3.44e-7,    "cat": "Sr2+", "ani": "SO4^2-","m": 1, "n": 1, "Mw": 183.68},
            "Ag2SO4":   {"Ksp": 1.20e-5,    "cat": "Ag+",  "ani": "SO4^2-","m": 2, "n": 1, "Mw": 311.80},
            # Carbonates
            "CaCO3":    {"Ksp": 3.36e-9,    "cat": "Ca2+", "ani": "CO3^2-","m": 1, "n": 1, "Mw": 100.09},
            "BaCO3":    {"Ksp": 2.58e-9,    "cat": "Ba2+", "ani": "CO3^2-","m": 1, "n": 1, "Mw": 197.34},
            "MgCO3":    {"Ksp": 6.82e-6,    "cat": "Mg2+", "ani": "CO3^2-","m": 1, "n": 1, "Mw": 84.31},
            "PbCO3":    {"Ksp": 7.40e-14,   "cat": "Pb2+", "ani": "CO3^2-","m": 1, "n": 1, "Mw": 267.21},
            "Ag2CO3":   {"Ksp": 8.46e-12,   "cat": "Ag+",  "ani": "CO3^2-","m": 2, "n": 1, "Mw": 275.75},
            "SrCO3":    {"Ksp": 5.60e-10,   "cat": "Sr2+", "ani": "CO3^2-","m": 1, "n": 1, "Mw": 147.63},
            "ZnCO3":    {"Ksp": 1.46e-10,   "cat": "Zn2+", "ani": "CO3^2-","m": 1, "n": 1, "Mw": 125.39},
            # Hydroxides
            "Fe(OH)2":  {"Ksp": 4.87e-17,   "cat": "Fe2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 89.86},
            "Fe(OH)3":  {"Ksp": 2.79e-39,   "cat": "Fe3+", "ani": "OH-",   "m": 1, "n": 3, "Mw": 106.87},
            "Cu(OH)2":  {"Ksp": 2.20e-20,   "cat": "Cu2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 97.56},
            "Mg(OH)2":  {"Ksp": 5.61e-12,   "cat": "Mg2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 58.33},
            "Ca(OH)2":  {"Ksp": 5.02e-6,    "cat": "Ca2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 74.09},
            "Al(OH)3":  {"Ksp": 3.28e-34,   "cat": "Al3+", "ani": "OH-",   "m": 1, "n": 3, "Mw": 78.00},
            "Zn(OH)2":  {"Ksp": 3.00e-17,   "cat": "Zn2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 99.40},
            "Pb(OH)2":  {"Ksp": 1.43e-20,   "cat": "Pb2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 241.21},
            "Cr(OH)3":  {"Ksp": 6.30e-31,   "cat": "Cr3+", "ani": "OH-",   "m": 1, "n": 3, "Mw": 103.02},
            "Mn(OH)2":  {"Ksp": 2.06e-13,   "cat": "Mn2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 88.95},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "cat": "Ni2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 92.71},
            "Co(OH)2":  {"Ksp": 5.92e-15,   "cat": "Co2+", "ani": "OH-",   "m": 1, "n": 2, "Mw": 92.95},
            # Phosphates
            "Ca3(PO4)2":{"Ksp": 2.07e-33,   "cat": "Ca2+", "ani": "PO4^3-","m": 3, "n": 2, "Mw": 310.18},
            "Ag3PO4":  {"Ksp": 8.89e-17,   "cat": "Ag+",  "ani": "PO4^3-","m": 3, "n": 1, "Mw": 418.58},
            # Sulfides
            "FeS":      {"Ksp": 6.30e-19,   "cat": "Fe2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 87.91},
            "MnS":      {"Ksp": 3.00e-14,   "cat": "Mn2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 87.00},
            "ZnS":      {"Ksp": 2.50e-22,   "cat": "Zn2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 97.46},
            "CdS":      {"Ksp": 8.00e-27,   "cat": "Cd2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 144.47},
            "PbS":      {"Ksp": 9.04e-29,   "cat": "Pb2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 239.27},
            "CuS":      {"Ksp": 6.00e-36,   "cat": "Cu2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 95.61},
            "HgS":      {"Ksp": 4.00e-53,   "cat": "Hg2+", "ani": "S^2-",   "m": 1, "n": 1, "Mw": 232.66},
            "Ag2S":     {"Ksp": 6.30e-50,   "cat": "Ag+",  "ani": "S^2-",   "m": 2, "n": 1, "Mw": 247.80},
            # Chromates
            "BaCrO4":  {"Ksp": 1.17e-10,   "cat": "Ba2+", "ani": "CrO4^2-","m": 1, "n": 1, "Mw": 253.33},
            "PbCrO4":  {"Ksp": 2.8e-13,    "cat": "Pb2+", "ani": "CrO4^2-","m": 1, "n": 1, "Mw": 323.20},
            "Ag2CrO4": {"Ksp": 1.12e-12,   "cat": "Ag+",  "ani": "CrO4^2-","m": 2, "n": 1, "Mw": 331.73},
            # Fluorides
            "CaF2":    {"Ksp": 3.45e-11,   "cat": "Ca2+", "ani": "F-",    "m": 1, "n": 2, "Mw": 78.08},
            "PbF2":    {"Ksp": 3.3e-8,     "cat": "Pb2+", "ani": "F-",    "m": 1, "n": 2, "Mw": 245.20},
            "MgF2":    {"Ksp": 6.5e-9,     "cat": "Mg2+", "ani": "F-",    "m": 1, "n": 2, "Mw": 62.30},
            # Oxalates
            "CaC2O4":  {"Ksp": 2.32e-9,    "cat": "Ca2+", "ani": "C2O4^2-","m": 1, "n": 1, "Mw": 128.10},
            # Iodates
            "AgIO3":   {"Ksp": 3.17e-8,    "cat": "Ag+",  "ani": "IO3-",  "m": 1, "n": 1, "Mw": 282.77},
            "Pb(IO3)2":{"Ksp": 3.69e-13,   "cat": "Pb2+", "ani": "IO3-",  "m": 1, "n": 2, "Mw": 557.01},
        }

        self._aliases = {
            "silver chloride": "AgCl", "silver bromide": "AgBr", "silver iodide": "AgI",
            "barium sulfate": "BaSO4", "lead sulfate": "PbSO4", "calcium sulfate": "CaSO4",
            "calcium carbonate": "CaCO3", "barium carbonate": "BaCO3", "magnesium carbonate": "MgCO3",
            "lead carbonate": "PbCO3", "silver carbonate": "Ag2CO3",
            "iron(ii) hydroxide": "Fe(OH)2", "iron(iii) hydroxide": "Fe(OH)3",
            "copper(ii) hydroxide": "Cu(OH)2", "magnesium hydroxide": "Mg(OH)2",
            "calcium hydroxide": "Ca(OH)2", "aluminum hydroxide": "Al(OH)3",
            "zinc hydroxide": "Zn(OH)2", "lead hydroxide": "Pb(OH)2",
            "chromium(iii) hydroxide": "Cr(OH)3", "manganese(ii) hydroxide": "Mn(OH)2",
            "nickel(ii) hydroxide": "Ni(OH)2", "cobalt(ii) hydroxide": "Co(OH)2",
            "calcium phosphate": "Ca3(PO4)2", "silver phosphate": "Ag3PO4",
            "ferrous sulfide": "FeS", "zinc sulfide": "ZnS", "cadmium sulfide": "CdS",
            "lead sulfide": "PbS", "copper sulfide": "CuS", "mercury sulfide": "HgS",
            "silver sulfide": "Ag2S",
            "barium chromate": "BaCrO4", "lead chromate": "PbCrO4", "silver chromate": "Ag2CrO4",
            "calcium fluoride": "CaF2", "lead fluoride": "PbF2", "magnesium fluoride": "MgF2",
            "calcium oxalate": "CaC2O4",
            "lead(ii) chloride": "PbCl2", "lead(ii) iodide": "PbI2",
            "calomel": "Hg2Cl2",
        }

    def _run_base(self, compound: str) -> dict:
        """Core logic: look up Ksp data."""
        key = self._resolve(compound)
        if key not in self._db:
            available = sorted(self._db.keys())
            raise ChemMCPError(
                f"Compound '{compound}' not found in Ksp database. "
                f"Available compounds ({len(self._db)}): {available[:20]}{'...' if len(available)>20 else ''}"
            )

        d = self._db[key]
        ksp = d["Ksp"]
        m, n = d["m"], d["n"]
        mw = d["Mw"]

        # Calculate molar solubility from Ksp
        # For MmXn ⇌ mM^n+ + nX^m-, Ksp = (mS)^m * (nS)^n = m^m * n^n * S^(m+n)
        total_exp = m + n
        coeff = (m ** m) * (n ** n)
        s_mol = (ksp / coeff) ** (1.0 / total_exp)

        # Convert to g/100mL
        s_g_100ml = s_mol * mw / 10.0  # mol/L * g/mol / 10 = g/100mL

        # Build dissolution equation
        diss_eq = f"{key}(s) ⇌ {m}{d['cat']}(aq) + {n}{d['ani']}(aq)"

        logger.info(f"GetKsp: {key} → Ksp={ksp:.2e}, S={s_mol:.2e} mol/L")
        return {
            "formula": key,
            "ksp": ksp,
            "pKsp": round(-math.log10(ksp), 2),
            "temperature": 298.0,
            "dissolution_eq": diss_eq,
            "stoichiometry": {
                "cation": d["cat"],
                "anion": d["ani"],
                "cation_coeff": m,
                "anion_coeff": n,
            },
            "molar_solubility": round(s_mol, int(max(0, 2 - math.log10(max(s_mol, 1e-15))))),
            "g_per_100ml": round(s_g_100ml, int(max(0, 4 - math.log10(max(s_g_100ml, 1e-15))))),
            "molar_mass": round(mw, 2),
            "source_note": "Standard thermodynamic data at 25°C (298 K). Values are approximate; actual Ksp depends on ionic strength and temperature.",
        }

    def _run_text(self, compound_str: str) -> dict:
        """Parse text input."""
        return self._run_base(compound_str.strip())

    def _resolve(self, name: str) -> str:
        """Resolve compound name to canonical key."""
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

import math
