"""
AAS Flame Selector — 原子吸收火焰类型选择工具 (#324)

功能：
  根据待测元素和样品基体特征，推荐最优的原子吸收光谱(AAS)火焰类型。
  支持空气-乙炔火焰和氧化亚氮-乙炔火焰两种主要类型。
"""

import logging
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 元素火焰推荐数据库 ────────────────────────────────────────
# 数据来源: Skoog & Holler "Principles of Instrumental Analysis"
# 以及各仪器厂商应用手册
ELEMENT_FLAME_DB: Dict[str, Dict[str, Any]] = {
    # ═══ 空气-乙炔火焰 (Air-C₂H₂) 元素 ═══
    # 温度 ~2300–2450 K
    "ag": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 328.1,
           "fuel_ratio": "lean", "notes": "Standard air-C2H2; good sensitivity."},
    "al": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 309.3,
           "fuel_ratio": "rich", "notes": "Forms refractory Al2O3; needs N2O-C2H2."},
    "au": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 242.8,
           "fuel_ratio": "lean", "notes": "Excellent sensitivity in air-C2H2."},
    "b":  {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 249.8,
           "fuel_ratio": "rich", "notes": "Refractory oxide; requires N2O-C2H2."},
    "ba": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 553.6,
           "fuel_ratio": "rich", "notes": "Forms stable BaO; N2O-C2H2 required."},
    "be": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 234.9,
           "fuel_ratio": "rich", "notes": "Toxic! Refractory BeO; N2O-C2H2 mandatory."},
    "bi": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 223.1,
           "fuel_ratio": "stoichiometric", "notes": "Good sensitivity in air-C2H2."},
    "ca": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 422.7,
           "fuel_ratio": "lean", "notes": "Add LaCl3 or EDTA to suppress PO4³⁻ interference."},
    "cd": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 228.8,
           "fuel_ratio": "lean", "notes": "Very high sensitivity; DLP < 0.001 ppm."},
    "co": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 240.7,
           "fuel_ratio": "lean", "notes": "Standard air-C2H2 measurement."},
    "cr": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 357.9,
           "fuel_ratio": "rich", "notes": "Rich flame improves Cr atom formation."},
    "cs": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 852.1,
           "fuel_ratio": "lean", "notes": "Ionization suppressor (KCl) recommended."},
    "cu": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 324.8,
           "fuel_ratio": "lean", "notes": "Excellent sensitivity; one of the best AAS elements."},
    "fe": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 248.3,
           "fuel_ratio": "lean", "notes": "Standard air-C2H2; add NH4Cl for Si interference."},
    "ga": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 287.4,
           "fuel_ratio": "lean", "notes": "N2O-C2H2 gives slightly better sensitivity."},
    "hg": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 253.7,
           "fuel_ratio": "lean", "notes": "Cold vapor AAS recommended for trace analysis."},
    "in": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 303.9,
           "fuel_ratio": "lean", "notes": "Good sensitivity in air-C2H2."},
    "k":  {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 766.5,
           "fuel_ratio": "lean", "notes": "Ionization easily suppressed with CsCl or Na as ionizer."},
    "li": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 670.8,
           "fuel_ratio": "lean", "notes": "Ionization suppressor needed (KCl)."},
    "mg": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 285.2,
           "fuel_ratio": "lean", "notes": "Excellent sensitivity; add Sr or La for Al/Si/P interference."},
    "mn": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 279.5,
           "fuel_ratio": "lean", "notes": "Standard air-C2H2; very stable signal."},
    "mo": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 313.3,
           "fuel_ratio": "rich", "notes": "Refractory MoOx; N2O-C2H2 recommended."},
    "na": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 589.0,
           "fuel_ratio": "lean", "notes": "Use ionization buffer (KCl). Emission mode also common."},
    "ni": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 232.0,
           "fuel_ratio": "lean", "notes": "Good sensitivity; multiple lines available."},
    "pb": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 283.3,
           "fuel_ratio": "stoichiometric", "notes": "Platform atomization improves precision."},
    "pd": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 244.8,
           "fuel_ratio": "lean", "notes": "Standard air-C2H2 analysis."},
    "pt": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 266.0,
           "fuel_ratio": "lean", "notes": "Moderate sensitivity in air-C2H2."},
    "rb": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 780.0,
           "fuel_ratio": "lean", "notes": "Ionization suppressor essential (K)."},
    "sb": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 217.6,
           "fuel_ratio": "lean", "notes": "N2O-C2H2 can improve sensitivity slightly."},
    "se": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 196.1,
           "fuel_ratio": "lean", "notes": "Hydride generation AAS preferred for better sensitivity."},
    "si": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 251.6,
           "fuel_ratio": "rich", "notes": "Forms refractory SiO2; N2O-C2H2 mandatory."},
    "sn": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 224.6,
           "fuel_ratio": "stoichiometric", "notes": "Hydride generation AAS recommended for trace levels."},
    "sr": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 460.7,
           "fuel_ratio": "rich", "notes": "Forms refractory SrO; N2O-C2H2 required."},
    "te": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 214.3,
           "fuel_ratio": "lean", "notes": "Hydride generation AAS preferred."},
    "ti": {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 364.3,
           "fuel_ratio": "rich", "notes": "Highly refractory TiO2; only N2O-C2H2 works."},
    "tl": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 276.8,
           "fuel_ratio": "lean", "notes": "Ionization suppressor recommended."},
    "v":  {"flame": "nitrous_oxide-acetylene", "temp_k": 2950, "wavelength_nm": 318.5,
           "fuel_ratio": "rich", "notes": "Refractory V-oxides; N2O-C2H2 required."},
    "zn": {"flame": "air-acetylene", "temp_k": 2300, "wavelength_nm": 213.9,
           "fuel_ratio": "lean", "notes": "Excellent sensitivity; DLP < 0.001 ppm."},
}


@ChemMCPManager.register_tool
class AasFlameSelector(BaseTool):
    """
    原子吸收火焰类型选择工具。
    根据元素符号和样品基体，推荐最优火焰类型、燃/氧化剂比及注意事项。
    """
    __version__                = "0.1.0"
    name                       = "AasFlameSelector"
    func_name                  = "select_flame_type"
    description                = ("Select optimal flame type (air-acetylene or nitrous oxide-acetylene) "
                                 "for Atomic Absorption Spectroscopy (AAS) based on element and sample matrix.")
    implementation_description = ("Uses built-in database of ~50 elements with recommended flame types, "
                                 "temperatures, fuel/oxidant ratios, and analytical notes.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["AAS", "Atomic Absorption", "Flame Spectroscopy",
                                   "Analytical Chemistry", "Elemental Analysis"]
    required_envs              = []

    code_input_sig = [
        ("element_symbol",              "str",   "N/A",         "Element symbol (e.g., 'Fe', 'Ca', 'Al'). Case-insensitive."),
        ("sample_matrix",               "str",   "general",     "Sample matrix type: 'general', 'aqueous', 'organic', 'high_salt'."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",         "Space-separated: element_symbol [sample_matrix]"),
    ]

    output_sig = [
        ("recommended_flame",           "str",   "Recommended flame: 'air-acetylene' or 'nitrous_oxide-acetylene'."),
        ("flame_temperature_K",         "float", "Approximate flame temperature in Kelvin."),
        ("analytical_wavelength_nm",    "float", "Most common AAS wavelength for this element (nm)."),
        ("fuel_to_oxidant_ratio",       "str",   "Recommended fuel-to-oxidant ratio: 'lean', 'stoichiometric', or 'rich'."),
        ("reason",                      "str",   "Explanation for the recommendation."),
        ("considerations",              "str",   "Special considerations, interferences, and tips."),
        ("element_info",                "dict",  "Full database entry for this element."),
    ]

    examples = [
        {
            "code_input": {
                "element_symbol": "Al",
                "sample_matrix": "aqueous",
            },
            "text_input": {"input_params": "Al aqueous"},
            "output": {
                "recommended_flame": "nitrous_oxide-acetylene",
                "flame_temperature_K": 2950.0,
            },
        },
        {
            "code_input": {
                "element_symbol": "Cu",
                "sample_matrix": "general",
            },
            "text_input": {"input_params": "Cu"},
            "output": {
                "recommended_flame": "air-acetylene",
                "flame_temperature_K": 2300.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load element-flame database."""
        self.db = ELEMENT_FLAME_DB

    def _run_base(
        self,
        element_symbol: str,
        sample_matrix: str = "general",
    ) -> Dict[str, Any]:
        """Core logic."""
        key = element_symbol.strip().lower()
        if key not in self.db:
            available = ", ".join(sorted(self.db.keys()))
            raise ChemMCPError(f"Element '{element_symbol}' not found in database. Available elements: {available}")

        info = self.db[key]
        flame = info["flame"]
        temp = info["temp_k"]
        wavelength = info["wavelength_nm"]
        ratio = info["fuel_ratio"]
        notes = info.get("notes", "")

        # Build reason based on chemistry
        if flame == "nitrous_oxide-acetylene":
            reason = (
                f"Element {key.upper()} forms thermally stable/refractory oxides that require "
                f"the higher temperature of N₂O-C₂H₂ (~{temp} K) for complete atomization. "
                f"Air-C₂H₂ (~2300 K) is insufficient."
            )
        else:
            reason = (
                f"Element {key.upper()} can be efficiently atomized in air-C₂H₂ flame "
                f"(~{temp} K). No refractory oxide formation issues."
            )

        # Matrix-specific considerations
        matrix_notes = ""
        matrix = sample_matrix.strip().lower()
        if matrix == "organic":
            matrix_notes = (
                "For organic samples: consider using a leaner flame to reduce carbon buildup, "
                "or perform acid digestion prior to analysis. Background correction (D₂ or Zeeman) is recommended."
            )
        elif matrix == "high_salt":
            matrix_notes = (
                "High salt matrix: use matrix-matching standards, consider standard addition method, "
                "and increase burner slot height to reduce clogging risk."
            )
        elif matrix == "aqueous":
            matrix_notes = "Aqueous samples are generally compatible with the recommended flame conditions."

        # Combine notes
        all_considerations = "; ".join(filter(None, [notes, matrix_notes]))

        logger.info(f"AAS flame selection for {key.upper()}: {flame} @ {temp}K, matrix={matrix}")
        return {
            "recommended_flame": flame,
            "flame_temperature_K": float(temp),
            "analytical_wavelength_nm": float(wavelength),
            "fuel_to_oxidant_ratio": ratio,
            "reason": reason,
            "considerations": all_considerations,
            "element_info": info,
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")
            elem = parts[0]
            matrix = parts[1] if len(parts) > 1 else "general"
            return self._run_base(elem, matrix)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'element_symbol [sample_matrix]'")
