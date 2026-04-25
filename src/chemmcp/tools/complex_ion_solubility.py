import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ComplexIonSolubility(BaseTool):
    """
    配离子形成对溶解度的影响。
    当存在能与阳离子形成配合物的配体时，难溶盐的溶解度可能显著增加（如AgCl在氨水中）。
    """
    __version__ = "0.1.0"
    name = "ComplexIonSolubility"
    func_name = "complex_ion_solubility"
    description = "Calculate how complex ion (coordination complex) formation affects solubility of sparingly soluble salts. When a ligand that complexes the cation is present, solubility can increase dramatically."
    implementation_description = "For MX(s) ⇌ M+ + X- (Ksp) and M+ + nL ⇌ [MLn]+ (Kf, formation constant). Overall: MX(s) + nL ⇌ [MLn]+ + X-, K = Ksp × Kf. Solves for S at given [L]. Built-in database of common metal-ligand formation constants."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Complex Ion", "Solubility", "Formation Constant", "Kf", "Coordination Chemistry", "Ligand"]
    required_envs = []

    code_input_sig = [
        ("compound", "str", "N/A", "Compound formula, e.g., 'AgCl', 'AgBr', 'Cu(OH)2', 'ZnS'."),
        ("ligand", "str", "N/A", "Ligand that forms a complex with the cation, e.g., 'NH3', 'CN-', 'OH-', 'NH3' for ammine complexes."),
        ("ligand_concentration", "float", "N/A", "Free ligand concentration in mol/L."),
        ("formation_constant", "float", "None", "Optional: override Kf value. If None, uses built-in database."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Format: 'compound ligand concentration'. E.g., 'AgCl NH3 1.0' or 'Cu(OH)2 NH3 2.0'."),
    ]

    output_sig = [
        ("compound", "str", "Compound formula analyzed."),
        ("ksp", "float", "Solubility product constant."),
        ("ligand", "str", "Ligand used."),
        ("kf", "float", "Formation constant (stability constant) used."),
        ("overall_k", "float", "Overall equilibrium constant K = Ksp × Kf."),
        ("solubility_pure_water", "float", "Molar solubility in pure water (mol/L)."),
        ("solubility_with_ligand", "float", "Molar solubility with ligand present (mol/L)."),
        ("enhancement_factor", "float", "Ratio: S_ligand / S_pure (>1 means enhancement)."),
        ("dominant_species", "str", "The dominant dissolved species at equilibrium."),
        ("stepwise_equilibria", "list", "Step-by-step equilibrium equations and calculations."),
        ("explanation", "str", "Detailed explanation of the complexation-enhanced dissolution mechanism."),
    ]

    examples = [
        {
            "code_input": {"compound": "AgCl", "ligand": "NH3", "ligand_concentration": 1.0, "formation_constant": None},
            "text_input": {"input_str": "AgCl NH3 1.0"},
            "output": {
                "compound": "AgCl",
                "ksp": 1.77e-10,
                "ligand": "NH3",
                "kf": 1.7e7,
                "overall_k": 3.01e-3,
                "solubility_pure_water": 1.33e-5,
                "solubility_with_ligand": 0.0548,
                "enhancement_factor": 4120,
                "dominant_species": "[Ag(NH3)2]+",
                "explanation": "AgCl dissolves ~4000× more in 1M NH3 due to Ag(NH3)2+ complex formation (K=Ksp×Kf=3.0×10⁻³).",
                "stepwise_equilibria": ["① AgCl(s) ⇌ Ag+ + Cl-; Ksp=1.77e-10", "② Ag+ + 2NH3 ⇌ [Ag(NH3)2]+; Kf=1.7e7", "③ Overall: K=3.01e-3"]
            },
        },
        {
            "code_input": {"compound": "AgBr", "ligand": "NH3", "ligand_concentration": 1.0, "formation_constant": None},
            "text_input": {"input_str": "AgBr NH3 1.0"},
            "output": {
                "compound": "AgBr",
                "ksp": 5.35e-13,
                "ligand": "NH3",
                "kf": 1.7e7,
                "overall_k": 9.10e-6,
                "solubility_pure_water": 7.32e-7,
                "solubility_with_ligand": 0.00302,
                "enhancement_factor": 4126,
                "dominant_species": "[Ag(NH3)2]+",
                "explanation": "AgBr is less soluble than AgCl in pure water, but NH3 complexation enhances both similarly.",
                "stepwise_equilibria": ["① AgBr(s) ⇌ Ag+ + Br-; Ksp=5.35e-13", "② Ag+ + 2NH3 ⇌ [Ag(NH3)2]+; Kf=1.7e7", "③ Overall: K=9.10e-6"]
            },
        },
        {
            "code_input": {"compound": "Cu(OH)2", "ligand": "NH3", "ligand_concentration": 2.0, "formation_constant": None},
            "text_input": {"input_str": "Cu(OH)2 NH3 2.0"},
            "output": {
                "compound": "Cu(OH)2",
                "ksp": 2.20e-20,
                "ligand": "NH3",
                "kf": 2.1e13,
                "overall_k": 462.0,
                "solubility_pure_water": 1.76e-7,
                "solubility_with_ligand": 0.0429,
                "enhancement_factor": 243750,
                "dominant_species": "[Cu(NH3)4]2+ (deep blue)",
                "explanation": "Cu(OH)2 is very insoluble but dissolves readily in excess NH3 to form deep blue [Cu(NH3)4]2+.",
                "stepwise_equilibria": ["① Cu(OH)2(s) ⇌ Cu2+ + 2OH-; Ksp=2.20e-20", "② Cu2+ + 4NH3 ⇌ [Cu(NH3)4]2+; Kf=2.1e13", "③ Overall: K=462.0"]
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Ksp and formation constant databases."""
        # Ksp database (same as other tools)
        self._ksp_db = {
            "AgCl":     {"Ksp": 1.77e-10,   "m": 1, "n": 1, "cat": "Ag+",  "ani": "Cl-",   "Mw": 143.32},
            "AgBr":     {"Ksp": 5.35e-13,   "m": 1, "n": 1, "cat": "Ag+",  "ani": "Br-",   "Mw": 187.77},
            "AgI":      {"Ksp": 8.52e-17,   "m": 1, "n": 1, "cat": "Ag+",  "ani": "I-",    "Mw": 234.77},
            "PbI2":     {"Ksp": 9.8e-9,     "m": 1, "n": 2, "cat": "Pb2+", "ani": "I-",    "Mw": 461.01},
            "Cu(OH)2":  {"Ksp": 2.20e-20,   "m": 1, "n": 2, "cat": "Cu2+", "ani": "OH-",   "Mw": 97.56},
            "Zn(OH)2":  {"Ksp": 3.00e-17,   "m": 1, "n": 2, "cat": "Zn2+", "ani": "OH-",   "Mw": 99.40},
            "Al(OH)3":  {"Ksp": 3.28e-34,   "m": 1, "n": 3, "cat": "Al3+", "ani": "OH-",   "Mw": 78.00},
            "Cr(OH)3":  {"Ksp": 6.30e-31,   "m": 1, "n": 3, "cat": "Cr3+", "ani": "OH-",   "Mw": 103.02},
            "Fe(OH)3":  {"Ksp": 2.79e-39,   "m": 1, "n": 3, "cat": "Fe3+", "ani": "OH-",   "Mw": 106.87},
            "HgS":      {"Ksp": 4.00e-53,   "m": 1, "n": 1, "cat": "Hg2+", "ani": "S^2-",   "Mw": 232.66},
            "ZnS":      {"Ksp": 2.50e-22,   "m": 1, "n": 1, "cat": "Zn2+", "ani": "S^2-",   "Mw": 97.46},
            "CdS":      {"Ksp": 8.00e-27,   "m": 1, "n": 1, "cat": "Cd2+", "ani": "S^2-",   "Mw": 144.47},
            "BaSO4":    {"Ksp": 1.08e-10,   "m": 1, "n": 1, "cat": "Ba2+", "ani": "SO4^2-","Mw": 233.39},
            "CaF2":    {"Ksp": 3.45e-11,   "m": 1, "n": 2, "cat": "Ca2+", "ani": "F-",    "Mw": 78.08},
            "PbCl2":    {"Ksp": 1.7e-5,     "m": 1, "n": 2, "cat": "Pb2+", "ani": "Cl-",   "Mw": 278.10},
            "Ni(OH)2":  {"Ksp": 5.48e-16,   "m": 1, "n": 2, "cat": "Ni2+", "ani": "OH-",   "Mw": 92.71},
            "Co(OH)2":  {"Ksp": 5.92e-15,   "m": 1, "n": 2, "cat": "Co2+", "ani": "OH-",   "Mw": 92.95},
            "Mn(OH)2":  {"Ksp": 2.06e-13,   "m": 1, "n": 2, "cat": "Mn2+", "ani": "OH-",   "Mw": 88.95},
            "Mg(OH)2":  {"Ksp": 5.61e-12,   "m": 1, "n": 2, "cat": "Mg2+", "ani": "OH-",   "Mw": 58.33},
            "CaCO3":    {"Ksp": 3.36e-9,    "m": 1, "n": 1, "cat": "Ca2+", "ani": "CO3^2-","Mw": 100.09},
            "Ag2S":     {"Ksp": 6.30e-50,   "m": 2, "n": 1, "cat": "Ag+",  "ani": "S^2-",   "Mw": 247.80},
        }

        # Formation constants (Kf / stability constants) database
        # Format: {cation: {ligand: {kf, coord, complex_formula}}}
        self._kf_db = {
            "Ag+": {
                "NH3":   {"kf": 1.7e7,    "coord": 2, "formula": "[Ag(NH₃)₂]⁺",   "name": "diamminesilver(I)"},
                "CN-":   {"kf": 5.6e18,   "coord": 2, "formula": "[Ag(CN)₂]⁻",    "name": "dicyanoargentate(I)"},
                "S2O3^2-": {"kf": 6.9e13,  "coord": 2, "formula": "[Ag(S₂O₃)₂]³⁻", "name": "dithiosilver(I)"},
                "Cl-":   {"kf": 2.0e5,    "coord": 3, "formula": "[AgCl₃]²⁻",     "name": "trichloroargentate(II) (weak)"},
                "NH3_stepwise": {"kf": [2.0e3, 6.6e3, 1.0e4], "coord": 2, "formula": "[Ag(NH₃)₂]⁺", "name": "stepwise β2=1.7e7"},
            },
            "Cu2+": {
                "NH3":   {"kf": 2.1e13,   "coord": 4, "formula": "[Cu(NH₃)₄]²⁺",  "name": "tetraamminecopper(II) deep blue"},
                "CN-":   {"kf": 1.0e24,   "coord": 2, "formula": "[Cu(CN)₄]³⁻",   "name": "tetracyanocuprate(I) (redox)"},
                "OH-":   {"kf": 5.0e18,   "coord": 4, "formula": "[Cu(OH)₄]²⁻",   "name": "tetrahydroxocuprate(II)"},
                "en":    {"kf": 1.0e20,   "coord": 3, "formula": "[Cu(en)₃]²⁺",   "name": "tris(ethylenediamine)copper(II)"},
                "Cl-":   {"kf": 1.0e5,    "coord": 4, "formula": "[CuCl₄]²⁻",     "name": "tetrachlorocuprate(II) green/yellow"},
            },
            "Zn2+": {
                "NH3":   {"kf": 2.9e9,    "coord": 4, "formula": "[Zn(NH₃)₄]²⁺",  "name": "tetraamminezinc(II)"},
                "OH-":   {"kf": 2.9e15,   "coord": 4, "formula": "[Zn(OH)₄]²⁻",   "name": "tetrahydroxozincate(II) / zincate"},
                "CN-":   {"kf": 1.0e16,   "coord": 4, "formula": "[Zn(CN)₄]²⁻",   "name": "tetracyanozincate(II)"},
            },
            "Al3+": {
                "OH-":   {"kf": 1.1e33,   "coord": 4, "formula": "[Al(OH)₄]⁻",    "name": "tetrahydroxoaluminate(III) / aluminate"},
                "F-":    {"kf": 1.0e13,   "coord": 6, "formula": "[AlF₆]³⁻",     "name": "hexafluoroaluminate(III)"},
            },
            "Cr3+": {
                "OH-":   {"kf": 1e28,     "coord": 4, "formula": "[Cr(OH)₄]⁻",    "name": "tetrahydroxochromate(III)"},
            },
            "Fe3+": {
                "CN-":   {"kf": 1.0e37,   "coord": 6, "formula": "[Fe(CN)₆]³⁻",   "name": "hexacyanoferrate(III) red"},
                "SCN-":  {"kf": 1.2e2,    "coord": 1, "formula": "[Fe(SCN)]²⁺",    "name": "thiocyanatoiron(III) blood-red"},
                "OH-":   {"kf": 2.0e20,   "coord": 2, "formula": "[Fe(OH)₄]⁻",    "name": "tetrahydroxoferrate(III)"},
                "C2O4^2-": {"kf": 1.0e15, "coord": 3, "formula": "[Fe(C₂O₄)₃]³⁻", "name": "tris(oxalato)ferrate(III) green"},
            },
            "Fe2+": {
                "CN-":   {"kf": 1.0e24,   "coord": 6, "formula": "[Fe(CN)₆]⁴⁻",   "name": "hexacyanoferrate(II) yellow"},
                "en":    {"kf": 4.0e9,    "coord": 3, "formula": "[Fe(en)₃]²⁺",   "name": "tris(ethylenediamine)iron(II)"},
            },
            "Hg2+": {
                "CN-":   {"kf": 2.5e41,   "coord": 4, "formula": "[Hg(CN)₄]²⁻",   "name": "tetracyanomercurate(II)"},
                "I-":    {"kf": 6.8e29,   "coord": 4, "formula": "[HgI₄]²⁻",     "name": "tetraiodomercurate(II) intense red"},
                "OH-":   {"kf": 1e25,     "coord": 4, "formula": "[Hg(OH)₄]²⁻",   "name": "tetrahydroxomercurate(II)"},
                "Cl-":   {"kf": 1.0e16,   "coord": 4, "formula": "[HgCl₄]²⁻",     "name": "tetrachloromercurate(II)"},
                "S^2-":  {"kf": 1e52,     "coord": 2, "formula": "[HgS₂]²⁻",     "name": "dithiomercurate(II)"},
                "NH3":   {"kf": 1.9e19,   "coord": 4, "formula": "[Hg(NH₃)₄]²⁺",  "name": "tetraamminemercury(II)"},
            },
            "Ni2+": {
                "NH3":   {"kf": 5.5e8,    "coord": 6, "formula": "[Ni(NH₃)₆]²⁺",  "name": "hexaamminenickel(II) blue-violet"},
                "CN-":   {"kf": 2.0e31,   "coord": 4, "formula": "[Ni(CN)₄]²⁻",   "name": "tetracyanonickelate(II)"},
                "en":    {"kf": 1.0e18,   "coord": 3, "formula": "[Ni(en)₃]²⁺",   "name": "tris(ethylenediamine)nickel(II)"},
            },
            "Co2+": {
                "NH3":   {"kf": 1.0e5,    "coord": 6, "formula": "[Co(NH₃)₆]²⁺",  "name": "hexaamminecobalt(II) pink"},
                "CN-":   {"kf": 1.0e19,   "coord": 6, "formula": "[Co(CN)₆]⁴⁻",   "name": "tetracyanocobaltate(II)"},
                "en":    {"kf": 1.0e13,   "coord": 3, "formula": "[Co(en)₃]²⁺",   "name": "tris(ethylenediamine)cobalt(II)"},
            },
            "Pb2+": {
                "OH-":   {"kf": 1e18,     "coord": 3, "formula": "[Pb(OH)₃]⁻",    "name": "trihydroxoplumbate(II)"},
                "CH3COO-": {"kf": 1e4,    "coord": 2, "formula": "Pb(CH₃COO)₂(aq)", "name": "lead acetate (soluble)"},
                "I-":    {"kf": 1e4,     "coord": 4, "formula": "[PbI₄]²⁻",     "name": "tetraiodoplumbate(II) yellow"},
            },
            "Cd2+": {
                "NH3":   {"kf": 1.0e7,    "coord": 4, "formula": "[Cd(NH₃)₄]²⁺",  "name": "tetraamminecadmium(II)"},
                "CN-":   {"kf": 1.0e18,   "coord": 4, "formula": "[Cd(CN)₄]²⁻",   "name": "tetracyanocadmate(II)"},
                "OH-":   {"kf": 1e11,     "coord": 4, "formula": "[Cd(OH)₄]²⁻",   "name": "tetrahydroxocadmate(II)"},
            },
            "Ca2+": {
                "C2O4^2-": {"kf": 1.7e7,  "coord": 1, "formula": "[Ca(C₂O₄)](aq)", "name": "calcium oxalate complex"},
                "EDTA":  {"kf": 1.0e11,   "coord": 1, "formula": "[Ca(EDTA)]²⁻",  "name": "calcium EDTA complex"},
            },
            "Ba2+": {
                "EDTA":  {"kf": 5.8e7,   "coord": 1, "formula": "[Ba(EDTA)]²⁻",  "name": "barium EDTA complex"},
            },
            "Mn2+": {
                "NH3":   {"kf": 1e1,      "coord": 1, "formula": "[Mn(NH₃)]²⁺",   "name": "very weak ammonia complex"},
                "C2O4^2-": {"kf": 1e4,    "coord": 1, "formula": "[Mn(C₂O₄)](aq)", "name": "manganese oxalate complex"},
            },
        }

        # Ligand name normalization
        self._ligand_aliases = {
            "nh3": "NH3", "ammonia": "NH3", "ammonia solution": "NH3",
            "cn-": "CN-", "cyanide": "CN-", "cn": "CN-",
            "oh-": "OH-", "hydroxide": "OH-", "naoh": "OH-",
            "s2o3^2-": "S2O3^2-", "thiosulfate": "S2O3^2-", "s2o3(2-)": "S2O3^2-",
            "cl-": "Cl-", "chloride": "Cl-", "cl": "Cl-",
            "i-": "I-", "iodide": "I-", "i": "I-",
            "scn-": "SCN-", "thiocyanate": "SCN-", "scn": "SCN-",
            "c2o4^2-": "C2O4^2-", "oxalate": "C2O4^2-", "oxalic acid": "C2O4^2-",
            "en": "en", "ethylenediamine": "en", "ethylene diamine": "en",
            "edta": "EDTA", "edta4-": "EDTA",
            "ch3coo-": "CH3COO-", "acetate": "CH3COO-", "acetate ion": "CH3COO-",
            "f-": "F-", "fluoride": "F-", "f": "F-",
            "s^2-": "S^2-", "sulfide": "S^2-", "s2-": "S^2-",
        }

    def _run_base(self, compound: str, ligand: str, ligand_concentration: float,
                  formation_constant: Optional[float] = None) -> dict:
        """Core logic: calculate solubility enhancement by complexation."""
        key = compound.strip()
        if key not in self._ksp_db:
            raise ChemMCPError(
                f"Unknown compound '{compound}'. Available ({len(self._ksp_db)}): "
                f"{sorted(self._ksp_db.keys())[:25]}{'...' if len(self._ksp_db)>25 else ''}"
            )

        d = self._ksp_db[key]
        ksp = d["Ksp"]
        m, n = d["m"], d["n"]
        cat = d["cat"]

        if ligand_concentration <= 0:
            raise ChemMCPError("Ligand concentration must be positive.")

        # Resolve ligand name
        lig_key = self._resolve_ligand(ligand)

        # Get Kf
        if formation_constant is not None:
            kf = formation_constant
            kf_info = {"kf": kf, "coord": 2, "formula": f"[ML_n]^+", "name": "user-provided"}
        else:
            # Look up Kf for this cation-ligand pair
            # Try exact match first, then with/without charge
            kf_info = None
            for cat_key_try in [cat, cat.replace('+', ''), cat.split('+')[0]]:
                if cat_key_try in self._kf_db:
                    cat_found = cat_key_try
                    break
            else:
                raise ChemMCPError(
                    f"No formation constant data for cation '{cat}'. "
                    f"Available cations: {sorted(self._kf_db.keys())}. "
                    f"Provide formation_constant explicitly."
                )
            
            if lig_key not in self._kf_db[cat_found]:
                available_ligs = list(self._kf_db[cat_found].keys())
                raise ChemMCPError(
                    f"No Kf data for {cat} + {ligand}. "
                    f"Available ligands for {cat}: {available_ligs}. "
                    f"Provide formation_constant explicitly."
                )
            kf_info = self._kf_db[cat_found][lig_key]
            kf = kf_info["kf"]

        coord = kf_info.get("coord", 2)
        complex_name = kf_info.get("formula", "[ML_n]+")

        # Pure water solubility
        total_exp = m + n
        coeff = (m ** m) * (n ** n)
        s_pure = (ksp / coeff) ** (1.0 / total_exp)

        # Overall reaction: MmXn(s) + m*coord*L ⇌ m[ML_coord]^(charge) + n*X^-
        # Simplified for MX type (m=n=1):
        # MX(s) + nL ⇌ [ML_n]^+ + X^-; K = Ksp * Kf
        # At equilibrium: S = [X-] = [[ML_n]^+]; [L]_free ≈ L_total (if L >> n*S)
        #
        # For general case, approximate:
        # K_overall = Ksp * Kf
        # For MX (1:1): S ≈ sqrt(Ksp * Kf * [L]^coord) when coordination consumes ligand

        L = ligand_concentration
        k_overall = ksp * kf

        # Calculate solubility with ligand
        # For MX (most common case, m=n=1):
        # Ksp = [M+][X-], Kf = [ML_n]/([M+][L]^n)
        # Mass balance: [M+]_total = S = [M+] + [ML_n]
        # Charge/material balance gives us: S = [X-]
        # From Ksp*Kf = [ML_n][X-]/[L]^n → if [L]≈L_initial:
        #   S^2 = Ksp*Kf*L^n → S = sqrt(Ksp*Kf*L^n)
        if m == 1 and n == 1:
            s_ligand = math.sqrt(k_overall * (L ** coord))
        elif m == 1:
            # MXn type: Ksp = [M][X]^n; mass: S = [M]_total; [X] = nS
            # With complexation: [M]_total = [M_free] + [ML]
            # Approximate: S ≈ (Ksp * Kf * L^coord)^(1/(n+1)) ... simplified
            try:
                s_ligand = (k_overall * (L ** coord)) ** (1.0 / (n + 1))
            except (OverflowError, ValueError, ZeroDivisionError):
                s_ligand = k_overall ** 0.5
        elif n == 1:
            # MmX type
            try:
                s_ligand = (k_overall * (L ** coord)) ** (1.0 / (m + 1))
            except (OverflowError, ValueError, ZeroDivisionError):
                s_ligand = k_overall ** 0.5
        else:
            # General case - rough estimate
            try:
                s_ligand = (k_overall * (L ** coord)) ** (1.0 / max(m, n) + 1)
            except (OverflowError, ValueError, ZeroDivisionError):
                s_ligand = k_overall ** 0.5

        # Sanity cap: solubility can't exceed total possible (e.g., pure solid amount)
        # But mathematically it can be > 1 M for very favorable complexation
        if s_ligand < 0:
            s_ligand = abs(s_ligand)

        # Enhancement factor
        if s_pure > 0:
            enhancement = s_ligand / s_pure
        else:
            enhancement = float('inf')

        # Build stepwise equilibria
        steps = []
        steps.append(f"① Dissolution: {key}(s) ⇌ {m}{cat}(aq) + {n}{d['ani']}(aq); Ksp = {ksp:.3e}")
        steps.append(f"② Complexation: {cat} + {coord}{lig_key} ⇌ {complex_name}; Kf = {kf:.3e}")
        steps.append(f"③ Overall: {key}(s) + {coord}{lig_key} ⇌ {complex_name} + {n}{d['ani']}(aq); K = Ksp×Kf = {k_overall:.3e}")
        steps.append(
            f"④ At [{lig_key}] = {L} M: S ≈ √(K×[{L}]^{coord}) "
            f"= √({k_overall:.3e} × {L**coord:.3e}) = {s_ligand:.4e} mol/L"
        )

        explanation = (
            f"The presence of {lig_key} dramatically increases {key}'s solubility through complex ion formation.\n\n"
            f"Mechanism: As {key} dissolves slightly, free {cat} ions are immediately "
            f"sequestered by {lig_key} to form stable {complex_name} (Kf={kf:.2e}).\n"
            f"This removes {cat} from solution, shifting the dissolution equilibrium rightward "
            f"(Le Chatelier's principle), causing more {key} to dissolve.\n\n"
            f"In pure water: S = {s_pure:.3e} mol/L\n"
            f"In {L} M {lig_key}: S = {s_ligand:.4e} mol/L ({enhancement:.1f}× enhancement)\n\n"
            f"Dominant species in solution: {complex_name}"
        )

        logger.info(f"ComplexIonSolubility: {key} + {lig_key}@{L}M → S={s_ligand:.4e} (pure={s_pure:.3e}, {enhancement:.0f}×)")
        return {
            "compound": key,
            "ksp": ksp,
            "ligand": lig_key,
            "kf": kf,
            "overall_k": round(k_overall, int(max(0, 3 - math.log10(max(abs(k_overall), 1e-100))))),
            "solubility_pure_water": round(s_pure, int(max(0, 2 - math.log10(max(s_pure, 1e-15))))),
            "solubility_with_ligand": round(s_ligand, int(max(0, 3 - math.log10(max(s_ligand, 1e-15))))),
            "enhancement_factor": round(enhancement, 1),
            "dominant_species": complex_name,
            "stepwise_equilibria": steps,
            "explanation": explanation,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        parts = input_str.strip().split()
        if len(parts) < 3:
            raise ChemMCPError(
                f"Need compound, ligand, and concentration. Got: '{input_str}'. "
                f"Format: 'compound ligand concentration'"
            )
        compound = parts[0]
        ligand = parts[1]
        try:
            conc = float(parts[2])
        except ValueError:
            raise ChemMCPError(f"Ligand concentration must be numeric. Got: '{parts[2]}'")
        kf_override = float(parts[3]) if len(parts) > 3 else None
        return self._run_base(compound, ligand, conc, kf_override)

    def _resolve_ligand(self, name: str) -> str:
        """Resolve ligand name/alias to canonical key."""
        n = name.strip()
        if n in self._kf_db.get("Ag+", {}):
            return n
        nl = n.lower()
        if nl in self._ligand_aliases:
            return self._ligand_aliases[nl]
        # Try direct match across all cations
        for cat_data in self._kf_db.values():
            if n in cat_data:
                return n
            for key in cat_data:
                if key.lower() == nl:
                    return key
        return n  # return as-is, will fail gracefully in _run_base
