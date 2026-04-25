import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalculateSolubility(BaseTool):
    """
    从Ksp计算溶解度。
    支持不同化学计量比的难溶盐（MX, MX2, M2X, MX3, M3X2等）。
    """
    __version__ = "0.1.0"
    name = "CalculateSolubility"
    func_name = "calculate_solubility_from_ksp"
    description = "Calculate molar solubility (mol/L) and mass solubility (g/100mL) from Ksp value for sparingly soluble salts with various stoichiometries."
    implementation_description = "For salt MmXn: Ksp = (mS)^m × (nS)^n, solve for S. Supports direct Ksp input or compound name lookup. Handles all common stoichiometry types."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Solubility", "Ksp", "Equilibrium", "Calculation"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "N/A", "Compound formula (e.g., 'AgCl', 'PbI2', 'Ag3PO4'). Used to look up Ksp and stoichiometry."),
        ("ksp_value", "float", "None", "Optional: directly provide Ksp value. If given, compound is only used for molar mass lookup."),
        ("temperature", "float", "298.0", "Temperature in Kelvin."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Format: 'compound [ksp_value]'. E.g., 'AgCl' or 'PbI2 9.8e-9'."),
    ]

    output_sig = [
        ("compound", "str", "Compound formula."),
        ("ksp_used", "float", "Ksp value used in calculation."),
        ("stoichiometry", "str", "Stoichiometry type description (e.g., 'MX (1:1)', 'MX2 (1:2)')."),
        ("molar_solubility", "float", "Molar solubility S in mol/L."),
        ("g_per_100ml", "float", "Mass solubility in g/100mL."),
        ("ion_concentrations", "dict", "Saturated concentrations of each ion: {ion: conc mol/L}."),
        ("formula_detail", "str", "Step-by-step formula derivation and calculation."),
    ]

    examples = [
        {
            "code_input": {"compound": "AgCl", "ksp_value": None, "temperature": 298.15},
            "text_input": {"input_str": "AgCl"},
            "output": {
                "compound": "AgCl",
                "ksp_used": 1.77e-10,
                "stoichiometry": "MX (1:1)",
                "molar_solubility": 1.33e-5,
                "g_per_100ml": 1.91e-4,
                "ion_concentrations": {"Ag+": 1.33e-5, "Cl-": 1.33e-5},
                "formula_detail": "AgCl ⇌ Ag+ + Cl-; Ksp = [Ag+][Cl-] = S×S = S²; S = √Ksp = √(1.77×10⁻¹⁰) = 1.33×10⁻⁵ mol/L",
            },
        },
        {
            "code_input": {"compound": "PbI2", "ksp_value": None, "temperature": 298.15},
            "text_input": {"input_str": "PbI2"},
            "output": {
                "compound": "PbI2",
                "ksp_used": 9.8e-9,
                "stoichiometry": "MX2 (1:2)",
                "molar_solubility": 1.35e-3,
                "g_per_100ml": 6.22e-4,
                "ion_concentrations": {"Pb2+": 1.35e-3, "I-": 2.69e-3},
                "formula_detail": "PbI2 ⇌ Pb²⁺ + 2I⁻; Ksp = [Pb²⁺][I⁻]² = S×(2S)² = 4S³; S = ∛(Ksp/4) = ∛(9.8×10⁻⁹/4) = 1.35×10⁻³ mol/L",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Ksp + molar mass database."""
        self._db = {
            # Same format as GetKsp
            "AgCl":     {"Ksp": 1.77e-10,   "m": 1, "n": 1, "Mw": 143.32},
            "AgBr":     {"Ksp": 5.35e-13,   "m": 1, "n": 1, "Mw": 187.77},
            "AgI":      {"Ksp": 8.52e-17,   "m": 1, "n": 1, "Mw": 234.77},
            "PbCl2":    {"Ksp": 1.7e-5,     "m": 1, "n": 2, "Mw": 278.10},
            "PbI2":     {"Ksp": 9.8e-9,     "m": 1, "n": 2, "Mw": 461.01},
            "BaSO4":    {"Ksp": 1.08e-10,   "m": 1, "n": 1, "Mw": 233.39},
            "CaCO3":    {"Ksp": 3.36e-9,    "m": 1, "n": 1, "Mw": 100.09},
            "Mg(OH)2":  {"Ksp": 5.61e-12,   "m": 1, "n": 2, "Mw": 58.33},
            "Fe(OH)3":  {"Ksp": 2.79e-39,   "m": 1, "n": 3, "Mw": 106.87},
            "Ag3PO4":  {"Ksp": 8.89e-17,   "m": 3, "n": 1, "Mw": 418.58},
            "Ag2CrO4": {"Ksp": 1.12e-12,   "m": 2, "n": 1, "Mw": 331.73},
            "Ca3(PO4)2":{"Ksp": 2.07e-33,   "m": 3, "n": 2, "Mw": 310.18},
            "Al(OH)3":  {"Ksp": 3.28e-34,   "m": 1, "n": 3, "Mw": 78.00},
            "CuS":      {"Ksp": 6.00e-36,   "m": 1, "n": 1, "Mw": 95.61},
            "ZnS":      {"Ksp": 2.50e-22,   "m": 1, "n": 1, "Mw": 97.46},
            "BaCrO4":  {"Ksp": 1.17e-10,   "m": 1, "n": 1, "Mw": 253.33},
            "CaF2":    {"Ksp": 3.45e-11,   "m": 1, "n": 2, "Mw": 78.08},
            "PbSO4":    {"Ksp": 2.53e-8,    "m": 1, "n": 1, "Mw": 303.26},
            "Cu(OH)2":  {"Ksp": 2.20e-20,   "m": 1, "n": 2, "Mw": 97.56},
            "HgS":      {"Ksp": 4.00e-53,   "m": 1, "n": 1, "Mw": 232.66},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "m": 1, "n": 2, "Mw": 92.71},
            "CaC2O4":  {"Ksp": 2.32e-9,    "m": 1, "n": 1, "Mw": 128.10},
            "Ag2CO3":   {"Ksp": 8.46e-12,   "m": 2, "n": 1, "Mw": 275.75},
            "Fe(OH)2":  {"Ksp": 4.87e-17,   "m": 1, "n": 2, "Mw": 89.86},
            "Mn(OH)2":  {"Ksp": 2.06e-13,   "m": 1, "n": 2, "Mw": 88.95},
            "Co(OH)2":  {"Ksp": 5.92e-15,   "m": 1, "n": 2, "Mw": 92.95},
            "Zn(OH)2":  {"Ksp": 3.00e-17,   "m": 1, "n": 2, "Mw": 99.40},
            "Pb(OH)2":  {"Ksp": 1.43e-20,   "m": 1, "n": 2, "Mw": 241.21},
            "Cr(OH)3":  {"Ksp": 6.30e-31,   "m": 1, "n": 3, "Mw": 103.02},
            "SrCO3":    {"Ksp": 5.60e-10,   "m": 1, "n": 1, "Mw": 147.63},
            "SrSO4":    {"Ksp": 3.44e-7,    "m": 1, "n": 1, "Mw": 183.68},
            "CdS":      {"Ksp": 8.00e-27,   "m": 1, "n": 1, "Mw": 144.47},
            "PbS":      {"Ksp": 9.04e-29,   "m": 1, "n": 1, "Mw": 239.27},
            "Ag2S":     {"Ksp": 6.30e-50,   "m": 2, "n": 1, "Mw": 247.80},
            "PbCrO4":  {"Ksp": 2.8e-13,    "m": 1, "n": 1, "Mw": 323.20},
            "MgF2":    {"Ksp": 6.5e-9,     "m": 1, "n": 2, "Mw": 62.30},
            "PbF2":    {"Ksp": 3.3e-8,     "m": 1, "n": 2, "Mw": 245.20},
            "Ag2SO4":   {"Ksp": 1.20e-5,    "m": 2, "n": 1, "Mw": 311.80},
            "BaCO3":    {"Ksp": 2.58e-9,    "m": 1, "n": 1, "Mw": 197.34},
            "MgCO3":    {"Ksp": 6.82e-6,    "m": 1, "n": 1, "Mw": 84.31},
            "PbCO3":    {"Ksp": 7.40e-14,   "m": 1, "n": 1, "Mw": 267.21},
            "ZnCO3":    {"Ksp": 1.46e-10,   "m": 1, "n": 1, "Mw": 125.39},
            "CaSO4":    {"Ksp": 4.93e-5,    "m": 1, "n": 1, "Mw": 136.14},
            "FeS":      {"Ksp": 6.30e-19,   "m": 1, "n": 1, "Mw": 87.91},
            "MnS":      {"Ksp": 3.00e-14,   "m": 1, "n": 1, "Mw": 87.00},
            "Hg2Cl2":   {"Ksp": 1.43e-18,   "m": 1, "n": 2, "Mw": 472.09},
            "CuCl":     {"Ksp": 1.72e-7,    "m": 1, "n": 1, "Mw": 98.999},
            "CuI":      {"Ksp": 1.27e-12,   "m": 1, "n": 1, "Mw": 190.45},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "m": 1, "n": 2, "Mw": 92.71},
            "Ag3PO4":  {"Ksp": 8.89e-17,   "m": 3, "n": 1, "Mw": 418.58},
        }

    def _run_base(self, compound: str, ksp_value: Optional[float] = None,
                  temperature: float = 298.0) -> dict:
        """Core logic: calculate solubility from Ksp."""
        key = compound.strip()
        if key not in self._db:
            raise ChemMCPError(
                f"Unknown compound '{compound}'. Available: {sorted(self._db.keys())[:30]}"
            )

        d = self._db[key]
        ksp = ksp_value if ksp_value is not None else d["Ksp"]
        m, n = d["m"], d["n"]
        mw = d["Mw"]

        if ksp <= 0:
            raise ChemMCPError("Ksp must be positive.")

        # General formula: Ksp = (m*S)^m * (n*S)^n = m^m * n^n * S^(m+n)
        total_exp = m + n
        coeff = (m ** m) * (n ** n)
        s_mol = (ksp / coeff) ** (1.0 / total_exp)

        # Ion concentrations at saturation
        cat_conc = m * s_mol
        ani_conc = n * s_mol

        # g/100mL
        s_g_100ml = s_mol * mw / 10.0

        # Stoichiometry label
        stoich_label = f"M{m}X{n} ({m}:{n})" if (m > 1 or n > 1) else f"MX ({m}:{n})"

        # Formula detail
        cat_sym = "M^m+" if m == 1 else f"M^{m}+"
        ani_sym = "X^n-" if n == 1 else f"X^{n}-"
        formula_detail = (
            f"{key} ⇌ {m}M + {n}X; "
            f"Ksp = [M]^m·[X]^n = ({m}S)^{m} · ({n}S)^{n} = "
            f"{coeff}·S^{total_exp}; "
            f"S = ({ksp:.2e}/{coeff})^(1/{total_exp}) = {s_mol:.3e} mol/L"
        )

        logger.info(f"CalculateSolubility: {compound} → S={s_mol:.3e} mol/L")
        return {
            "compound": key,
            "ksp_used": ksp,
            "stoichiometry": stoich_label,
            "molar_solubility": round(s_mol, int(max(0, 2 - math.log10(max(s_mol, 1e-15))))),
            "g_per_100ml": round(s_g_100ml, int(max(0, 4 - math.log10(max(s_g_100ml, 1e-15))))),
            "ion_concentrations": {
                "cation": round(cat_conc, int(max(0, 2 - math.log10(max(cat_conc, 1e-15))))),
                "anion": round(ani_conc, int(max(0, 2 - math.log10(max(ani_conc, 1e-15))))),
            },
            "formula_detail": formula_detail,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        parts = input_str.strip().split()
        compound = parts[0]
        ksp_val = float(parts[1]) if len(parts) > 1 else None
        temp = float(parts[2]) if len(parts) > 2 else 298.0
        return self._run_base(compound, ksp_val, temp)
