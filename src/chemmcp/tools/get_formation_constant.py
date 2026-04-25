"""
配合物稳定常数（逐级/累积）查询工具
Query formation/stability constants (stepwise Kn and cumulative βn) for coordination complexes.
"""
import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetFormationConstant(BaseTool):
    """
    查询配合物的逐级形成常数(Kn)和累积稳定常数(βn)。
    内置常见配合物体系的文献数据。
    """
    __version__ = "0.1.0"
    name = "GetFormationConstant"
    func_name = "get_formation_constant"
    description = "Query stepwise formation constants (Kn) and cumulative stability constants (βn) for coordination complexes."
    implementation_description = "Uses a built-in database of literature log K values for common complex systems. Returns both stepwise Kn and cumulative βn values."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Equilibrium", "Stability Constants", "Complexes"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal ion symbol, e.g., 'Cu2+', 'Ag+', 'Fe3+'"),
        ("ligand", "str", "N/A", "Ligand formula or name, e.g., 'NH3', 'CN-', 'en', 'OH-'."),
        ("step", "int", "0", "Specific step number to query (1-indexed). 0 or None returns all steps."),
        ("constant_type", "str", "all", "Type of constant: 'stepwise' (Kn), 'cumulative' (βn), or 'all' (both)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion ligand [step] [type]', e.g., 'Cu2+ NH3' or 'Ag+ CN- 2 cumulative'."),
    ]

    output_sig = [
        ("metal_ion", "str", "The metal ion queried."),
        ("ligand", "str", "The ligand queried."),
        ("max_coordination", "int", "Maximum coordination number for this system."),
        ("stepwise_constants", "dict", "Stepwise formation constants K1, K2, ... Kn as {n: log_Kn}."),
        ("cumulative_constants", "dict", "Cumulative stability constants β1, β2, ... βn as {n: log_βn}."),
        ("overall_formation", "str", "Overall formation reaction and total βn."),
        ("source_note", "str", "Data source / conditions (typically 25°C, I≈0)."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Cu2+",
                "ligand": "NH3",
                "step": 0,
                "constant_type": "all",
            },
            "text_input": {
                "query": "Cu2+ NH3"
            },
            "output": {
                "metal_ion": "Cu2+",
                "ligand": "NH3",
                "max_coordination": 4,
                "stepwise_constants": {"1": "4.15", "2": "3.50", "3": "2.89", "4": "2.13"},
                "cumulative_constants": {"1": "4.15", "2": "7.65", "3": "10.54", "4": "12.67"},
                "overall_formation": "[Cu(NH3)4]2+: log β4 = 12.67",
                "source_note": "Data at 25°C, I ≈ 0; typical textbook values.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Ag+",
                "ligand": "NH3",
                "step": 2,
                "constant_type": "stepwise",
            },
            "text_input": {
                "query": "Ag+ NH3 2 stepwise"
            },
            "output": {
                "metal_ion": "Ag+",
                "ligand": "NH3",
                "max_coordination": 2,
                "stepwise_constants": {"2": "3.87"},
                "cumulative_constants": {},
                "overall_formation": "[Ag(NH3)2]+: log β2 = 7.23",
                "source_note": "Data at 25°C, I ≈ 0.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize the formation constant database with literature values."""
        # Format: metal_ligand -> { "max_n": int, "log_K": [K1, K2, ...], "notes": str }
        # All values are log10(K) or log10(β) at ~25°C, ionic strength ≈ 0
        self._db = {
            # --- Cu(II) systems ---
            "cu2+_nh3": {
                "metal": "Cu2+", "ligand": "NH3", "max_n": 4,
                "log_K": [4.15, 3.50, 2.89, 2.13],
                "notes": "[Cu(NH3)4]2+ deep blue; data from Greenwood & Earnshaw"
            },
            "cu2+_cn": {
                "metal": "Cu2+", "ligand": "CN-", "max_n": 4,
                "log_K": [None, 5.48, None, None],  # stepwise often incomplete
                "log_beta": [None, 19.6, None, 30.3],
                "notes": "[Cu(CN)4]3- very stable"
            },
            "cu2+_en": {
                "metal": "Cu2+", "ligand": "en", "max_n": 2,
                "log_K": [10.71, 9.25],
                "notes": "ethylenediamine; strong chelate effect"
            },
            "cu2+_oh": {
                "metal": "Cu2+", "ligand": "OH-", "max_n": 4,
                "log_K": [6.5, 5.4, 4.4, 3.5],
                "notes": "hydroxo complexes; pH dependent"
            },

            # --- Ag(I) systems ---
            "ag+_nh3": {
                "metal": "Ag+", "ligand": "NH3", "max_n": 2,
                "log_K": [3.36, 3.87],
                "notes": "[Ag(NH3)2]+; Tollens' reagent"
            },
            "ag+_cn": {
                "metal": "Ag+", "ligand": "CN-", "max_n": 2,
                "log_K": [None, None],  # usually given as beta
                "log_beta": [None, 21.1],
                "notes": "[Ag(CN)2]- extremely stable; used in silver extraction"
            },
            "ag+_s2o3": {
                "metal": "Ag+", "ligand": "S2O3^2-", "max_n": 2,
                "log_K": [8.82, 4.64],
                "notes": "thiosulfate; photographic fixer"
            },
            "ag+_i": {
                "metal": "Ag+", "ligand": "I-", "max_n": 3,
                "log_K": [6.58, 5.11, 0.66],
                "notes": "[AgI3]2-"
            },
            "ag+_nh3_alt": {
                "metal": "Ag+", "ligand": "NH3", "max_n": 2,
                "log_K": [3.31, 3.91],
                "notes": "alternative dataset"
            },

            # --- Fe(III) systems ---
            "fe3+_cn": {
                "metal": "Fe3+", "ligand": "CN-", "max_n": 6,
                "log_K": [None, None, None, None, None, None],
                "log_beta": [None, None, None, None, None, 43.9],
                "notes": "[Fe(CN)6]3-, ferricyanide; extremely stable"
            },
            "fe3+_scn": {
                "metal": "Fe3+", "ligand": "SCN-", "max_n": 5,
                "log_K": [2.36, 1.41, 0.62, 0.18, -0.34],
                "notes": "blood-red [Fe(SCN)(H2O)5]2+; stepwise decreasing"
            },
            "fe3+_f": {
                "metal": "Fe3+", "ligand": "F-", "max_n": 6,
                "log_K": [5.28, 4.02, 2.99, 1.50, -0.04, -1.14],
                "notes": "[FeF6]3- colorless; weak field"
            },

            # --- Fe(II) systems ---
            "fe2+_cn": {
                "metal": "Fe2+", "ligand": "CN-", "max_n": 6,
                "log_K": [None, None, None, None, None, None],
                "log_beta": [None, None, None, None, None, 35.4],
                "notes": "[Fe(CN)6]4-, ferrocyanide"
            },
            "fe2+_phen": {
                "metal": "Fe2+", "ligand": "phen", "max_n": 3,
                "log_K": [5.85, 5.45, 9.55],
                "notes": "1,10-phenanthroline; ferroin indicator (red)"
            },
            "fe2+_bpy": {
                "metal": "Fe2+", "ligand": "bpy", "max_n": 3,
                "log_K": [4.20, 3.65, 8.80],
                "notes": "2,2'-bipyridine"
            },

            # --- Co(III) systems ---
            "co3+_nh3": {
                "metal": "Co3+", "ligand": "NH3", "max_n": 6,
                "log_K": [None, None, None, None, None, None],
                "log_beta": [None, None, None, None, None, 33.7],
                "notes": "[Co(NH3)6]3+ kinetically inert, very stable"
            },
            "co3+_cn": {
                "metal": "Co3+", "ligand": "CN-", "max_n": 6,
                "log_K": [None, None, None, None, None, None],
                "log_beta": [None, None, None, None, None, 54.0],
                "notes": "[Co(CN)6]3- extraordinarily stable"
            },
            "co2+_nh3": {
                "metal": "Co2+", "ligand": "NH3", "max_n": 6,
                "log_K": [2.11, 1.74, 1.09, 0.76, 0.12, -0.51],
                "notes": "[Co(NH3)6]2+ pink; labile"
            },

            # --- Ni(II) systems ---
            "ni2+_nh3": {
                "metal": "Ni2+", "ligand": "NH3", "max_n": 6,
                "log_K": [2.80, 2.24, 1.73, 1.19, 0.81, 0.34],
                "notes": "[Ni(NH3)6]2+ violet-blue"
            },
            "ni2+_cn": {
                "metal": "Ni2+", "ligand": "CN-", "max_n": 4,
                "log_K": [None, None, None, None],
                "log_beta": [None, None, None, 30.3],
                "notes": "[Ni(CN)4]2- square planar"
            },
            "ni2+_en": {
                "metal": "Ni2+", "ligand": "en", "max_n": 3,
                "log_K": [7.96, 6.32, 4.49],
                "notes": "[Ni(en)3]2+ chelate; violet"
            },

            # --- Zn(II) systems ---
            "zn2+_nh3": {
                "metal": "Zn2+", "ligand": "NH3", "max_n": 4,
                "log_K": [2.37, 2.44, 2.50, 2.15],
                "notes": "[Zn(NH3)4]2+"
            },
            "zn2+_cn": {
                "metal": "Zn2+", "ligand": "CN-", "max_n": 4,
                "log_K": [None, None, None, None],
                "log_beta": [None, None, None, 16.7],
                "notes": "[Zn(CN)4]2-"
            },
            "zn2+_oh": {
                "metal": "Zn2+", "ligand": "OH-", "max_n": 4,
                "log_K": [4.4, 4.9, 4.1, 1.0],
                "notes": "zincate formation at high pH"
            },
            "zn2+_en": {
                "metal": "Zn2+", "ligand": "en", "max_n": 3,
                "log_K": [5.92, 5.15, 1.86],
                "notes": "chelate effect"
            },

            # --- Al(III) systems ---
            "al3+_f": {
                "metal": "Al3+", "ligand": "F-", "max_n": 6,
                "log_K": [6.97, 5.47, 3.95, 2.67, 1.07, 0.21],
                "notes": "[AlF6]3- colorless; used in cryolite"
            },
            "al3+_oh": {
                "metal": "Al3+", "ligand": "OH-", "max_n": 4,
                "log_K": [9.0, 8.4, 7.5, 6.3],
                "notes": "aluminate [Al(OH)4]- at high pH"
            },

            # --- Cd(II) systems ---
            "cd2+_cn": {
                "metal": "Cd2+", "ligand": "CN-", "max_n": 4,
                "log_K": [5.54, 5.08, 4.57, 2.78],
                "notes": "[Cd(CN)4]2-"
            },
            "cd2+_nh3": {
                "metal": "Cd2+", "ligand": "NH3", "max_n": 4,
                "log_K": [2.65, 2.10, 1.44, 0.93],
                "notes": "[Cd(NH3)4]2+"
            },
            "cd2+_cl": {
                "metal": "Cd2+", "ligand": "Cl-", "max_n": 4,
                "log_K": [1.95, 0.55, 0.10, 0.31],
                "notes": "[CdCl4]2-"
            },
            "cd2+_i": {
                "metal": "Cd2+", "ligand": "I-", "max_n": 4,
                "log_K": [2.38, 0.94, 0.53, 0.14],
                "notes": "[CdI4]2-"
            },

            # --- Hg(II) systems ---
            "hg2+_cl": {
                "metal": "Hg2+", "ligand": "Cl-", "max_n": 4,
                "log_K": [6.74, 6.48, 0.85, 1.00],
                "notes": "[HgCl4]2-"
            },
            "hg2+_i": {
                "metal": "Hg2+", "ligand": "I-", "max_n": 4,
                "log_K": [12.87, 10.95, 3.17, 2.23],
                "notes": "[HgI4]2-; Nessler's reagent base"
            },
            "hg2+_cn": {
                "metal": "Hg2+", "ligand": "CN-", "max_n": 4,
                "log_K": [None, None, None, None],
                "log_beta": [None, None, None, 41.5],
                "notes": "[Hg(CN)4]2- extremely stable"
            },
            "hg2+_s2o3": {
                "metal": "Hg2+", "ligand": "S2O3^2-", "max_n": 4,
                "log_K": [None, None, None, None],
                "log_beta": [None, None, None, 29.8],
                "notes": "[Hg(S2O3)4]6-"
            },
            "hg2+_nh3": {
                "metal": "Hg2+", "ligand": "NH3", "max_n": 4,
                "log_K": [8.8, 8.70, 1.00, 0.78],
                "notes": "[Hg(NH3)4]2+"
            },

            # --- Cr(III) system ---
            "cr3+_nh3": {
                "metal": "Cr3+", "ligand": "NH3", "max_n": 6,
                "log_K": [None, None, None, None, None, None],
                "log_beta": [None, None, None, None, None, 26.0],
                "notes": "[Cr(NH3)6]3+ kinetically inert"
            },
            "cr3+_scn": {
                "metal": "Cr3+", "ligand": "SCN-", "max_n": 6,
                "log_K": [3.52, 1.67, 0.83, 0.37, 0.12, -0.14],
                "notes": "[Cr(SCN)6]3-"
            },

            # --- Pb(II) system ---
            "pb2+_i": {
                "metal": "Pb2+", "ligand": "I-", "max_n": 4,
                "log_K": [2.16, 1.27, 0.86, 0.56],
                "notes": "[PbI4]2- yellow precipitate (bright yellow)"
            },
            "pb2+_oh": {
                "metal": "Pb2+", "ligand": "OH-", "max_n": 3,
                "log_K": [6.2, 3.6, 2.5],
                "notes": "plumbite [Pb(OH)3]-"
            },

            # --- EDTA systems (hexadentate) ---
            "ca2+_edta": {
                "metal": "Ca2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [10.69],
                "notes": "CaY2-; important in water hardness"
            },
            "mg2+_edta": {
                "metal": "Mg2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [8.64],
                "notes": "MgY2-"
            },
            "fe3+_edta": {
                "metal": "Fe3+", "ligand": "EDTA", "max_n": 1,
                "log_K": [25.1],
                "notes": "FeY-; yellow"
            },
            "cu2+_edta": {
                "metal": "Cu2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [18.80],
                "notes": "CuY2-; deep blue"
            },
            "ni2+_edta": {
                "metal": "Ni2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [18.56],
                "notes": "NiY2-; blue"
            },
            "zn2+_edta": {
                "metal": "Zn2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [16.50],
                "notes": "ZnY2-"
            },
            "co2+_edta": {
                "metal": "Co2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [16.26],
                "notes": "CoY2-; rose/pink"
            },
            "cd2+_edta": {
                "metal": "Cd2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [16.46],
                "notes": "CdY2-"
            },
            "pb2+_edta": {
                "metal": "Pb2+", "ligand": "EDTA", "max_n": 1,
                "log_K": [18.3],
                "notes": "PbY2-"
            },
            "al3+_edta": {
                "metal": "Al3+", "ligand": "EDTA", "max_n": 1,
                "log_K": [16.5],
                "notes": "AlY-; slow formation"
            },

            # --- Oxalate systems ---
            "fe3+_c2o4": {
                "metal": "Fe3+", "ligand": "C2O4^2-", "max_n": 3,
                "log_K": [9.4, 6.8, 4.0],
                "notes": "[Fe(C2O4)3]3- green"
            },
            "fe2+_c2o4": {
                "metal": "Fe2+", "ligand": "C2O4^2-", "max_n": 3,
                "log_K": [4.2, 3.1, 1.0],
                "notes": "[Fe(C2O4)3]4-"
            },
            "co3+_c2o4": {
                "metal": "Co3+", "ligand": "C2O4^2-", "max_n": 3,
                "log_K": [None, None, None],
                "log_beta": [None, None, ~20],
                "notes": "[Co(C2O4)3]3- green"
            },

            # --- Acetylacetonate (acac) ---
            "co3+_acac": {
                "metal": "Co3+", "ligand": "acac", "max_n": 3,
                "log_K": [None, None, None],
                "log_beta": [None, None, ~15],
                "notes": "[Co(acac)3]"
            },
            "cr3+_acac": {
                "metal": "Cr3+", "ligand": "acac", "max_n": 3,
                "log_K": [None, None, None],
                "log_beta": [None, None, ~14],
                "notes": "[Cr(acac)3]"
            },
        }

        # Ligand name normalization map
        self._ligand_aliases = {
            "ammonia": "nh3", "ammine": "nh3", "ammonia": "nh3",
            "cyanide": "cn", "cyano": "cn", "cn-": "cn", "cn": "cn",
            "hydroxide": "oh", "hydroxo": "oh", "oh-": "oh", "oh": "oh",
            "fluoride": "f", "fluoro": "f", "f-": "f", "f": "f",
            "chloride": "cl", "chloro": "cl", "cl-": "cl", "cl": "cl",
            "iodide": "i", "iodo": "i", "i-": "i", "i": "i",
            "thiosulfate": "s2o3", "s2o3^2-": "s2o3", "s2o3": "s2o3",
            "thiocyanate": "scn", "scn-": "scn", "thiocyanato": "scn", "scn": "scn",
            "oxalate": "c2o4", "ox": "c2o4", "c2o4^2-": "c2o4", "c2o4": "c2o4",
            "ethylenediamine": "en", "en": "en",
            "phenanthroline": "phen", "phen": "phen",
            "bipyridine": "bpy", "bpy": "bpy",
            "edta": "edta", "edta": "edta",
            "acetylacetonate": "acac", "acac": "acac",
            "water": "h2o", "aqua": "h2o", "h2o": "h2o",
            "carbonyl": "co", "co": "co",
            "nitro": "no2", "nitro": "no2", "no2-": "no2", "no2": "no2",
        }

        # Metal ion normalization
        self._metal_aliases = {
            "cu(ii)": "cu2+", "copper(ii)": "cu2+", "copper(2+)": "cu2+",
            "cu(i)": "cu+", "copper(i)": "cu+", "copper(1+)": "cu+",
            "fe(ii)": "fe2+", "iron(ii)": "fe2+", "iron(2+)": "fe2+",
            "fe(iii)": "fe3+", "iron(iii)": "fe3+", "iron(3+)": "fe3+",
            "co(ii)": "co2+", "cobalt(ii)": "co2+", "cobalt(2+)": "co2+",
            "co(iii)": "co3+", "cobalt(iii)": "co3+", "cobalt(3+)": "co3+",
            "ni(ii)": "ni2+", "nickel(ii)": "ni2+", "nickel(2+)": "ni2+",
            "zn(ii)": "zn2+", "zinc(ii)": "zn2+", "zinc(2+)": "zn2+",
            "ag(i)": "ag+", "silver(i)": "ag+", "silver(1+)": "ag+",
            "cd(ii)": "cd2+", "cadmium(ii)": "cd2+", "cadmium(2+)": "cd2+",
            "hg(ii)": "hg2+", "mercury(ii)": "hg2+", "mercury(2+)": "hg2+",
            "al(iii)": "al3+", "aluminum(iii)": "al3+", "aluminum(3+)": "al3+",
            "cr(iii)": "cr3+", "chromium(iii)": "cr3+", "chromium(3+)": "cr3+",
            "pb(ii)": "pb2+", "lead(ii)": "pb2+", "lead(2+)": "pb2+",
            "ca(ii)": "ca2+", "calcium(ii)": "ca2+", "calcium(2+)": "ca2+",
            "mg(ii)": "mg2+", "magnesium(ii)": "mg2+", "magnesium(2+)": "mg2+",
        }

    def _normalize_key(self, metal_ion: str, ligand: str) -> str:
        """Normalize metal + ligand into a lookup key."""
        m = metal_ion.lower().replace(" ", "")
        l = ligand.lower().replace(" ", "").replace("^", "")

        # Apply metal aliases
        m = self._metal_aliases.get(m, m)

        # Apply ligand aliases
        l = self._ligand_aliases.get(l, l)

        return f"{m}_{l}"

    def _run_base(self, metal_ion: str, ligand: str, step: int = 0,
                  constant_type: str = "all") -> dict:
        """Query formation constants for a metal-ligand system."""
        key = self._normalize_key(metal_ion, ligand)

        if key not in self._db:
            # Try case-insensitive search
            matches = [k for k in self._db.keys()
                       if metal_ion.lower() in k.lower() or ligand.lower() in k.lower()]
            if matches:
                raise ChemMCPError(
                    f"Exact match not found for ({metal_ion}, {ligand}). "
                    f"Did you mean: {matches}? "
                    f"Available metals: {set(k.split('_')[0] for k in self._db.keys())}"
                )
            raise ChemMCPError(
                f"No data found for metal='{metal_ion}', ligand='{ligand}'. "
                f"Available systems include: Cu2+/NH3, Ag+/CN-, Fe3+/SCN-, Co3+/NH3, "
                f"Ni2+/NH3, Zn2+/NH3, Al3+/F-, Cd2+/CN-, Hg2+/I-, various EDTA/oxalate systems."
            )

        entry = self._db[key]
        max_n = entry["max_n"]
        log_K_list = entry.get("log_K", [])
        log_beta_list = entry.get("log_beta", [])

        # Compute cumulative betas from stepwise if not provided
        computed_betas = {}
        cumsum = 0.0
        for i, k_val in enumerate(log_K_list):
            if k_val is not None:
                cumsum += k_val
            computed_betas[i + 1] = round(cumsum, 2) if k_val is not None else None

        # Use provided betas if available (override computed)
        final_betas = dict(computed_betas)
        if log_beta_list:
            for i, b_val in enumerate(log_beta_list):
                if b_val is not None:
                    final_betas[i + 1] = round(b_val, 2) if isinstance(b_val, (int, float)) else b_val

        # Build stepwise result
        stepwise_result = {}
        for i, k_val in enumerate(log_K_list):
            if k_val is not None:
                stepwise_result[str(i + 1)] = f"{k_val:.2f}" if isinstance(k_val, float) else str(k_val)
            else:
                stepwise_result[str(i + 1)] = "N/A"

        # Build cumulative result
        cumulative_result = {}
        for n in range(1, max_n + 1):
            val = final_betas.get(n)
            cumulative_result[str(n)] = f"{val:.2f}" if isinstance(val, (int, float)) else (val or "N/A")

        # Filter by step if requested
        if step > 0:
            if constant_type == "stepwise":
                s = str(step)
                if s in stepwise_result:
                    stepwise_result = {s: stepwise_result[s]}
                else:
                    stepwise_result = {s: "N/A"}
                cumulative_result = {}
            elif constant_type == "cumulative":
                s = str(step)
                if s in cumulative_result:
                    cumulative_result = {s: cumulative_result[s]}
                else:
                    cumulative_result = {s: "N/A"}
                stepwise_result = {}

        # Overall formation string
        overall = f"[{entry['metal']}({entry['ligand']}){max_n}]{self._get_complex_charge(entry['metal'], entry['ligand'])}"
        last_beta = final_betas.get(max_n, "N/A")
        if isinstance(last_beta, (int, float)):
            overall += f": log β{max_n} = {last_beta:.2f}"
        else:
            overall += f": log β{max_n} = {last_beta}"

        logger.info(f"Queried formation constants: {key} → {overall}")

        return {
            "metal_ion": entry["metal"],
            "ligand": entry["ligand"],
            "max_coordination": max_n,
            "stepwise_constants": stepwise_result if constant_type in ("all", "stepwise") else {},
            "cumulative_constants": cumulative_result if constant_type in ("all", "cumulative") else {},
            "overall_formation": overall,
            "source_note": entry.get("notes", "Data at 25°C, I ≈ 0."),
        }

    def _get_complex_charge(self, metal: str, ligand: str) -> str:
        """Estimate complex charge (simplified)."""
        metal_charge = 0
        m = metal.replace("+", "").replace("-", "")
        if "+" in metal:
            metal_charge = metal.count("+")
        if "-" in metal:
            metal_charge = -metal.count("-")

        # Simple ligand charge estimation
        ligand_charges = {
            "nh3": 0, "h2o": 0, "en": 0, "phen": 0, "bpy": 0, "acac": -1,
            "cn": -1, "oh": -1, "f": -1, "cl": -1, "i": -1, "scn": -1,
            "s2o3": -2, "c2o4": -2, "edta": -4, "no2": -1, "co": 0,
        }
        lg = ligand.lower().replace("^", "").replace("2-", "-").replace("3-", "-").replace("4-", "-")
        lg = self._ligand_aliases.get(lg, lg)
        lg_charge = ligand_charges.get(lg, 0)

        total = metal_charge + lg_charge  # simplified (single ligand)
        if total > 0:
            return "+" + str(total) if total > 1 else "+"
        elif total < 0:
            return str(total)
        return ""

    def _run_text(self, query: str) -> dict:
        """Parse text query: 'metal_ion ligand [step] [type]'"""
        parts = query.strip().split()
        if len(parts) < 2:
            raise ChemMCPError(
                "Format: 'metal_ion ligand [step] [constant_type]'. "
                "Example: 'Cu2+ NH3', 'Ag+ CN- 2 cumulative'"
            )

        metal_ion = parts[0]
        ligand = parts[1]
        step = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 0
        const_type = parts[3] if len(parts) > 3 else "all"

        return self._run_base(metal_ion, ligand, step, const_type)
