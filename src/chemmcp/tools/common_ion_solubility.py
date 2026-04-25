import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CommonIonSolubility(BaseTool):
    """
    同离子效应对溶解度的影响。
    计算在含有共同离子的溶液中，难溶盐的溶解度变化（与纯水中的溶解度比较）。
    """
    __version__ = "0.1.0"
    name = "CommonIonSolubility"
    func_name = "common_ion_solubility"
    description = "Calculate the effect of a common ion on the solubility of a sparingly soluble salt. Compares solubility in pure water vs. in a solution containing a common ion, with quantitative reduction factor."
    implementation_description = "For salt MmXn with Ksp, with common ion X at concentration C: Ksp = (mS)^m × (nS + C)^n. Solves for S numerically/analytically. When C >> nS, approximates as S ≈ [Ksp / (m^m × C^n)]^(1/m). Also calculates pure water solubility for comparison."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Common Ion Effect", "Solubility", "Ksp", "Equilibrium", "Le Chatelier"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "N/A", "Compound formula, e.g., 'AgCl', 'BaSO4', 'PbI2', 'Ag3PO4'."),
        ("common_ion", "str", "N/A", "The common ion present in solution, e.g., 'Cl-', 'SO4^2-', 'I-', 'Ag+'."),
        ("common_ion_conc", "float", "N/A", "Concentration of the common ion in mol/L."),
        ("temperature", "float", "298.0", "Temperature in Kelvin (default 25°C)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Format: 'compound common_ion concentration'. E.g., 'AgCl Cl- 0.1' or 'PbI2 I- 0.01'."),
    ]

    output_sig = [
        ("compound", "str", "Compound formula."),
        ("ksp", "float", "Solubility product constant used."),
        ("stoichiometry", "str", "Stoichiometry type (e.g., 'MX 1:1', 'MX2 1:2')."),
        ("solubility_pure_water", "float", "Molar solubility in pure water (mol/L)."),
        ("solubility_with_common_ion", "float", "Molar solubility in presence of common ion (mol/L)."),
        ("reduction_factor", "float", "Ratio: S_pure / S_common (>1 means reduction)."),
        ("percent_decrease", "float", "Percentage decrease in solubility due to common ion effect."),
        ("ion_concentrations", "dict", "Saturated concentrations of each ion with common ion present."),
        ("explanation", "str", "Step-by-step calculation and Le Chatelier explanation."),
    ]

    examples = [
        {
            "code_input": {"compound": "AgCl", "common_ion": "Cl-", "common_ion_conc": 0.1, "temperature": 298.0},
            "text_input": {"input_str": "AgCl Cl- 0.1"},
            "output": {
                "compound": "AgCl",
                "ksp": 1.77e-10,
                "stoichiometry": "MX (1:1)",
                "solubility_pure_water": 1.33e-5,
                "solubility_with_common_ion": 1.77e-9,
                "reduction_factor": 7514.1,
                "percent_decrease": 98.67,
                "explanation": "AgCl ⇌ Ag+ + Cl-. Pure water: S=√Ksp=1.33×10⁻⁵ M. With 0.1M Cl-: Ksp=S×(S+0.1)≈S×0.1, S≈1.77×10⁻⁹ M. Solubility reduced by ~7500×.",
                "ion_concentrations": {"Ag+": 1.77e-9, "Cl-": 0.1},
            },
        },
        {
            "code_input": {"compound": "BaSO4", "common_ion": "SO4^2-", "common_ion_conc": 0.01, "temperature": 298.0},
            "text_input": {"input_str": "BaSO4 SO4^2- 0.01"},
            "output": {
                "compound": "BaSO4",
                "ksp": 1.08e-10,
                "stoichiometry": "MX (1:1)",
                "solubility_pure_water": 1.04e-5,
                "solubility_with_common_ion": 1.08e-8,
                "reduction_factor": 962.9,
                "percent_decrease": 99.90,
                "explanation": "BaSO4 ⇌ Ba²⁺ + SO₄²⁻. Pure water: S=√Ksp=1.04×10⁻⁵ M. With 0.01M SO₄²⁻: S≈Ksp/[SO₄²⁻]=1.08×10⁻⁸ M. ~963-fold reduction.",
                "ion_concentrations": {"Ba2+": 1.08e-8, "SO4^2-": 0.01},
            },
        },
        {
            "code_input": {"compound": "PbI2", "common_ion": "I-", "common_ion_conc": 0.01, "temperature": 298.0},
            "text_input": {"input_str": "PbI2 I- 0.01"},
            "output": {
                "compound": "PbI2",
                "ksp": 9.8e-9,
                "stoichiometry": "MX2 (1:2)",
                "solubility_pure_water": 1.35e-3,
                "solubility_with_common_ion": 9.8e-5,
                "reduction_factor": 13.78,
                "percent_decrease": 92.74,
                "explanation": "PbI2 ⇌ Pb²⁺ + 2I⁻; Ksp=[Pb²⁺][I⁻]². Pure water: S=∛(Ksp/4)=1.35×10⁻³ M. With 0.01M I⁻: Ksp=S×(0.01)², S=9.8×10⁻⁵ M. ~13.8-fold reduction.",
                "ion_concentrations": {"Pb2+": 9.8e-5, "I-": 0.01},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Ksp database."""
        self._db = {
            # Format: {formula: {Ksp, m, n, Mw}}
            # m = cation coefficient, n = anion coefficient
            "AgCl":     {"Ksp": 1.77e-10,   "m": 1, "n": 1, "cat": "Ag+",  "ani": "Cl-",   "Mw": 143.32},
            "AgBr":     {"Ksp": 5.35e-13,   "m": 1, "n": 1, "cat": "Ag+",  "ani": "Br-",   "Mw": 187.77},
            "AgI":      {"Ksp": 8.52e-17,   "m": 1, "n": 1, "cat": "Ag+",  "ani": "I-",    "Mw": 234.77},
            "PbCl2":    {"Ksp": 1.7e-5,     "m": 1, "n": 2, "cat": "Pb2+", "ani": "Cl-",   "Mw": 278.10},
            "PbI2":     {"Ksp": 9.8e-9,     "m": 1, "n": 2, "cat": "Pb2+", "ani": "I-",    "Mw": 461.01},
            "Hg2Cl2":   {"Ksp": 1.43e-18,   "m": 1, "n": 2, "cat": "Hg2^2+","ani": "Cl-",  "Mw": 472.09},
            "CuCl":     {"Ksp": 1.72e-7,    "m": 1, "n": 1, "cat": "Cu+",  "ani": "Cl-",   "Mw": 98.999},
            "CuI":      {"Ksp": 1.27e-12,   "m": 1, "n": 1, "cat": "Cu+",  "ani": "I-",    "Mw": 190.45},
            "BaSO4":    {"Ksp": 1.08e-10,   "m": 1, "n": 1, "cat": "Ba2+", "ani": "SO4^2-","Mw": 233.39},
            "PbSO4":    {"Ksp": 2.53e-8,    "m": 1, "n": 1, "cat": "Pb2+", "ani": "SO4^2-","Mw": 303.26},
            "CaSO4":    {"Ksp": 4.93e-5,    "m": 1, "n": 1, "cat": "Ca2+", "ani": "SO4^2-","Mw": 136.14},
            "SrSO4":    {"Ksp": 3.44e-7,    "m": 1, "n": 1, "cat": "Sr2+", "ani": "SO4^2-","Mw": 183.68},
            "Ag2SO4":   {"Ksp": 1.20e-5,    "m": 2, "n": 1, "cat": "Ag+",  "ani": "SO4^2-","Mw": 311.80},
            "CaCO3":    {"Ksp": 3.36e-9,    "m": 1, "n": 1, "cat": "Ca2+", "ani": "CO3^2-","Mw": 100.09},
            "BaCO3":    {"Ksp": 2.58e-9,    "m": 1, "n": 1, "cat": "Ba2+", "ani": "CO3^2-","Mw": 197.34},
            "MgCO3":    {"Ksp": 6.82e-6,    "m": 1, "n": 1, "cat": "Mg2+", "ani": "CO3^2-","Mw": 84.31},
            "PbCO3":    {"Ksp": 7.40e-14,   "m": 1, "n": 1, "cat": "Pb2+", "ani": "CO3^2-","Mw": 267.21},
            "Ag2CO3":   {"Ksp": 8.46e-12,   "m": 2, "n": 1, "cat": "Ag+",  "ani": "CO3^2-","Mw": 275.75},
            "SrCO3":    {"Ksp": 5.60e-10,   "m": 1, "n": 1, "cat": "Sr2+", "ani": "CO3^2-","Mw": 147.63},
            "ZnCO3":    {"Ksp": 1.46e-10,   "m": 1, "n": 1, "cat": "Zn2+", "ani": "CO3^2-","Mw": 125.39},
            "Fe(OH)2":  {"Ksp": 4.87e-17,   "m": 1, "n": 2, "cat": "Fe2+", "ani": "OH-",   "Mw": 89.86},
            "Fe(OH)3":  {"Ksp": 2.79e-39,   "m": 1, "n": 3, "cat": "Fe3+", "ani": "OH-",   "Mw": 106.87},
            "Cu(OH)2":  {"Ksp": 2.20e-20,   "m": 1, "n": 2, "cat": "Cu2+", "ani": "OH-",   "Mw": 97.56},
            "Mg(OH)2":  {"Ksp": 5.61e-12,   "m": 1, "n": 2, "cat": "Mg2+", "ani": "OH-",   "Mw": 58.33},
            "Ca(OH)2":  {"Ksp": 5.02e-6,    "m": 1, "n": 2, "cat": "Ca2+", "ani": "OH-",   "Mw": 74.09},
            "Al(OH)3":  {"Ksp": 3.28e-34,   "m": 1, "n": 3, "cat": "Al3+", "ani": "OH-",   "Mw": 78.00},
            "Zn(OH)2":  {"Ksp": 3.00e-17,   "m": 1, "n": 2, "cat": "Zn2+", "ani": "OH-",   "Mw": 99.40},
            "Pb(OH)2":  {"Ksp": 1.43e-20,   "m": 1, "n": 2, "cat": "Pb2+", "ani": "OH-",   "Mw": 241.21},
            "Cr(OH)3":  {"Ksp": 6.30e-31,   "m": 1, "n": 3, "cat": "Cr3+", "ani": "OH-",   "Mw": 103.02},
            "Mn(OH)2":  {"Ksp": 2.06e-13,   "m": 1, "n": 2, "cat": "Mn2+", "ani": "OH-",   "Mw": 88.95},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "m": 1, "n": 2, "cat": "Ni2+", "ani": "OH-",   "Mw": 92.71},
            "Co(OH)2":  {"Ksp": 5.92e-15,   "m": 1, "n": 2, "cat": "Co2+", "ani": "OH-",   "Mw": 92.95},
            "Ca3(PO4)2":{"Ksp": 2.07e-33,   "m": 3, "n": 2, "cat": "Ca2+", "ani": "PO4^3-","Mw": 310.18},
            "Ag3PO4":  {"Ksp": 8.89e-17,   "m": 3, "n": 1, "cat": "Ag+",  "ani": "PO4^3-","Mw": 418.58},
            "FeS":      {"Ksp": 6.30e-19,   "m": 1, "n": 1, "cat": "Fe2+", "ani": "S^2-",   "Mw": 87.91},
            "MnS":      {"Ksp": 3.00e-14,   "m": 1, "n": 1, "cat": "Mn2+", "ani": "S^2-",   "Mw": 87.00},
            "ZnS":      {"Ksp": 2.50e-22,   "m": 1, "n": 1, "cat": "Zn2+", "ani": "S^2-",   "Mw": 97.46},
            "CdS":      {"Ksp": 8.00e-27,   "m": 1, "n": 1, "cat": "Cd2+", "ani": "S^2-",   "Mw": 144.47},
            "PbS":      {"Ksp": 9.04e-29,   "m": 1, "n": 1, "cat": "Pb2+", "ani": "S^2-",   "Mw": 239.27},
            "CuS":      {"Ksp": 6.00e-36,   "m": 1, "n": 1, "cat": "Cu2+", "ani": "S^2-",   "Mw": 95.61},
            "HgS":      {"Ksp": 4.00e-53,   "m": 1, "n": 1, "cat": "Hg2+", "ani": "S^2-",   "Mw": 232.66},
            "Ag2S":     {"Ksp": 6.30e-50,   "m": 2, "n": 1, "cat": "Ag+",  "ani": "S^2-",   "Mw": 247.80},
            "BaCrO4":  {"Ksp": 1.17e-10,   "m": 1, "n": 1, "cat": "Ba2+", "ani": "CrO4^2-","Mw": 253.33},
            "PbCrO4":  {"Ksp": 2.8e-13,    "m": 1, "n": 1, "cat": "Pb2+", "ani": "CrO4^2-","Mw": 323.20},
            "Ag2CrO4": {"Ksp": 1.12e-12,   "m": 2, "n": 1, "cat": "Ag+",  "ani": "CrO4^2-","Mw": 331.73},
            "CaF2":    {"Ksp": 3.45e-11,   "m": 1, "n": 2, "cat": "Ca2+", "ani": "F-",    "Mw": 78.08},
            "PbF2":    {"Ksp": 3.3e-8,     "m": 1, "n": 2, "cat": "Pb2+", "ani": "F-",    "Mw": 245.20},
            "MgF2":    {"Ksp": 6.5e-9,     "m": 1, "n": 2, "cat": "Mg2+", "ani": "F-",    "Mw": 62.30},
            "CaC2O4":  {"Ksp": 2.32e-9,    "m": 1, "n": 1, "cat": "Ca2+", "ani": "C2O4^2-","Mw": 128.10},
        }

    def _run_base(self, compound: str, common_ion: str, common_ion_conc: float,
                  temperature: float = 298.0) -> dict:
        """Core logic: calculate common ion effect on solubility."""
        key = compound.strip()
        if key not in self._db:
            raise ChemMCPError(
                f"Unknown compound '{compound}'. Available ({len(self._db)}): "
                f"{sorted(self._db.keys())[:25]}{'...' if len(self._db)>25 else ''}"
            )

        d = self._db[key]
        ksp = d["Ksp"]
        m, n = d["m"], d["n"]
        cat, ani = d["cat"], d["ani"]

        if common_ion_conc <= 0:
            raise ChemMCPError("Common ion concentration must be positive.")

        # Determine which ion is the common ion
        is_cation_common = self._ion_match(common_ion, cat)
        is_anion_common = self._ion_match(common_ion, ani)

        if not is_cation_common and not is_anion_common:
            raise ChemMCPError(
                f"Common ion '{common_ion}' does not match either ion of {key} ({cat}, {ani})."
            )

        # Pure water solubility
        total_exp = m + n
        coeff = (m ** m) * (n ** n)
        s_pure = (ksp / coeff) ** (1.0 / total_exp)

        # Solubility with common ion
        if is_anion_common:
            # Common anion: Ksp = (m*S)^m * (n*S + C)^n ≈ (m*S)^m * C^n when C >> nS
            C = common_ion_conc
            # Approximate: S ≈ [Ksp / (m^m * C^n)]^(1/m)
            s_common = (ksp / ((m ** m) * (C ** n))) ** (1.0 / m)
            cat_conc = m * s_common
            ani_eff = n * s_common + C
        else:
            # Common cation: Ksp = (m*S + C)^m * (n*S)^n ≈ C^m * (n*S)^n when C >> mS
            C = common_ion_conc
            s_common = (ksp / ((C ** m) * (n ** n))) ** (1.0 / n)
            cat_eff = m * s_common + C
            ani_conc = n * s_common

        # Reduction metrics
        if s_common > 0:
            reduction = s_pure / s_common
            pct_decrease = (1 - s_common / s_pure) * 100
        else:
            reduction = float('inf')
            pct_decrease = 100.0

        # Stoichiometry label
        stoich_label = f"M{m}X{n} ({m}:{n})" if (m > 1 or n > 1) else f"MX ({m}:{n})"

        # Build explanation
        common_label = common_ion
        if is_anion_common:
            ion_concs = {cat: round(cat_conc, int(max(0, 2 - math.log10(max(abs(cat_conc), 1e-15))))),
                         ani: round(ani_eff, int(max(0, 2 - math.log10(max(abs(ani_eff), 1e-15)))))}
            explanation = (
                f"{key} ⇌ {m}{cat}(aq) + {n}{ani}(aq); Ksp = {ksp:.3e}\n"
                f"Pure water: S = ∛(Ksp/{coeff}) = {s_pure:.3e} mol/L\n"
                f"With [{common_label}] = {common_ion_conc} M (common ion):\n"
                f"  Ksp = ({m}S)^{m} × ({n}S + {common_ion_conc})^{n} ≈ ({m}S)^{m} × ({common_ion_conc})^{n}\n"
                f"  S ≈ [{ksp:.3e} / ({m**m} × {common_ion_conc}**{n})]^(1/{m}) = {s_common:.3e} mol/L\n"
                f"Reduction factor: {reduction:.1f}× | Solubility decreased by {pct_decrease:.1f}%\n"
                f"[Le Chatelier: Adding {common_label} shifts equilibrium LEFT, reducing dissolution]"
            )
        else:
            ion_concs = {cat: round(cat_eff, int(max(0, 2 - math.log10(max(abs(cat_eff), 1e-15))))),
                         ani: round(ani_conc, int(max(0, 2 - math.log10(max(abs(ani_conc), 1e-15)))))}
            explanation = (
                f"{key} ⇌ {m}{cat}(aq) + {n}{ani}(aq); Ksp = {ksp:.3e}\n"
                f"Pure water: S = ∛(Ksp/{coeff}) = {s_pure:.3e} mol/L\n"
                f"With [{common_label}] = {common_ion_conc} M (common ion):\n"
                f"  Ksp = ({m}S + {common_ion_conc})^{m} × ({n}S)^{n} ≈ ({common_ion_conc})^{m} × ({n}S)^{n}\n"
                f"  S ≈ [{ksp:.3e} / ({common_ion_conc}**{m} × {n**n})]^(1/{n}) = {s_common:.3e} mol/L\n"
                f"Reduction factor: {reduction:.1f}× | Solubility decreased by {pct_decrease:.1f}%\n"
                f"[Le Chatelier: Adding {common_label} shifts equilibrium LEFT, reducing dissolution]"
            )

        logger.info(f"CommonIonSolubility: {key} + {common_ion}@{common_ion_conc} → S={s_common:.3e} (pure={s_pure:.3e})")
        return {
            "compound": key,
            "ksp": ksp,
            "stoichiometry": stoich_label,
            "solubility_pure_water": round(s_pure, int(max(0, 2 - math.log10(max(s_pure, 1e-15))))),
            "solubility_with_common_ion": round(s_common, int(max(0, 2 - math.log10(max(s_common, 1e-15))))),
            "reduction_factor": round(reduction, 1),
            "percent_decrease": round(pct_decrease, 2),
            "ion_concentrations": ion_concs,
            "explanation": explanation,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        parts = input_str.strip().split()
        if len(parts) < 3:
            raise ChemMCPError(f"Need compound, common ion, and concentration. Got: '{input_str}'. Format: 'compound ion conc'")
        compound = parts[0]
        ion = parts[1]
        try:
            conc = float(parts[2])
        except ValueError:
            raise ChemMCPError(f"Concentration must be a number. Got: '{parts[2]}'")
        temp = float(parts[3]) if len(parts) > 3 else 298.0
        return self._run_base(compound, ion, conc, temp)

    @staticmethod
    def _ion_match(ion1: str, ion2: str) -> bool:
        """Check if two ion strings refer to the same ion."""
        i1 = ion1.replace("^", "").replace("_", "").replace(" ", "").lower()
        i2 = ion2.replace("^", "").replace("_", "").replace(" ", "").lower()
        return i1 == i2
