"""
ICP-MS Isotope Selector — ICP-MS同位素选择与多原子干扰规避工具 (#327)

功能：
  根据待测元素、基体组成和仪器类型，选择最优同位素用于ICP-MS分析，
  规避多原子离子干扰（如 ArO⁺, ArN⁺, ArAr⁺, 氧化物/双电荷干扰等）。
"""

import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── ICP-MS 同位素数据库 ──────────────────────────────────────
# 数据来源: IUPAC Isotopic Compositions, NIST, 各仪器厂商应用手册
# 包含: 同位素质量, 丰度(%), 主要多原子干扰, 推荐程度
ICP_MS_ISOTOPE_DB: Dict[str, List[Dict[str, Any]]] = {
    "li": [
        {"mass": 6, "abundance": 7.59, "interferences": [], "notes": "Low abundance; use for Li isotope ratio work."},
        {"mass": 7, "abundance": 92.41, "interferences": [], "notes": "Recommended for quantification. No significant interferences."},
    ],
    "b": [
        {"mass": 10, "abundance": 19.9, "interferences": [], "notes": "Clean; lower abundance."},
        {"mass": 11, "abundance": 80.1, "interferences": [], "notes": "Recommended — no polyatomic interferences."},
    ],
    "mg": [
        {"mass": 24, "abundance": 78.99, "interferences": ["C2 (24)", "Na2 (48→charge)"],
         "notes": "Na matrix causes Na2++ interference at m/z 12 (not direct). Generally OK."},
        {"mass": 25, "abundance": 10.00, "interferences": ["C2H (25)"],
         "notes": "C2H from organic solvent/CO2; minor in clean samples."},
        {"mass": 26, "abundance": 11.01, "interferences": ["CN (26)", "C2H2 (26)"],
         "notes": "CN band common in plasma; use collision cell or choose Mg-24."},
    ],
    "cr": [
        {"mass": 50, "abundance": 4.345, "interferences": ["TiV (50)", "ClO+? no", "Ar3C? rare"],
         "notes": "Very low abundance; Ti/V overlap possible."},
        {"mass": 52, "abundance": 83.79, "interferences": ["ArC (52)", "ClOH (52)", "S2 (64?no)"],
         "notes": "Most abundant but ArC (40+12) interferes! Use Cr-53 or CRC mode."},
        {"mass": 53, "abundance": 9.501, "interferences": ["ArC1H (53)", "ClO (53)"],
         "notes": "Recommended if Cl present; moderate abundance."},
        {"mass": 54, "abundance": 2.365, "interferences": ["ArN (54)", "Fe (54 isobaric!)"],
         "notes": "Fe-54 isobaric overlap — avoid if Fe >> Cr."},
    ],
    "mn": [
        {"mass": 55, "abundance": 100, "interferences": ["ArNH (55)", "ArO1H (55)"],
         "notes": "Only stable isotope. Use He/H2 collision mode if needed."},
    ],
    "co": [
        {"mass": 59, "abundance": 100, "interferences": ["ArNa (59)", "CaO (59)", "NO3H?"],
         "notes": "Only stable isotope. CaO and ArNa are main concerns. Use KED/CRC."},
    ],
    "ni": [
        {"mass": 58, "abundance": 68.08, "interferences": ["CaO (58)", "Fe (58 isobaric!)", "KCl (58)"],
         "notes": "Fe-58 isobaric! Avoid if Fe-rich matrix."},
        {"mass": 60, "abundance": 26.22, "interferences": ["CaO (60)", "Ni (60)"], "notes": "CaO interference common."},
        {"mass": 61, "abundance": 1.140, "interferences": [], "notes": "Clean but very low abundance."},
        {"mass": 62, "abundance": 3.634, "interferences": ["CaO (62)"], "notes": "CaO interference."},
    ],
    "cu": [
        {"mass": 63, "abundance": 69.17, "interferences": ["ArNa (63)", "PO2 (63)", "TiO (64?no)"],
         "notes": "ArNa from Na matrix; PO2 from P-containing samples."},
        {"mass": 65, "abundance": 30.83, "interferences": ["ArMg (65)", "S2 (64?no), CaOH (65)"],
         "notes": "If Na-free matrix, Cu-65 may be cleaner. Use pair ratio check."},
    ],
    "zn": [
        {"mass": 64, "abundance": 48.63, "interferences": ["Ni (64 isobaric!)", "MgAr (64)", "S2 (64)"],
         "notes": "Ni-64 isobaric! Major issue in Ni alloys/samples."},
        {"mass": 66, "abundance": 27.90, "interferences": ["TiO (66)", "Ca2O? no", "Ba2++(66)"],
         "notes": "TiO interference in Ti-containing samples."},
        {"mass": 67, "abundance": 4.102, "interferences": ["TiO1H (67)", "SO2H (67)"],
         "notes": "Lower abundance but often cleaner."},
        {"mass": 68, "abundance": 18.75, "interferences": ["Zn (68)", "MgAr (68)", "Ca2O? no", "Sr2++(68)"],
         "notes": "MgAr, Sr doubly charged."}
    ],
    "as": [
        {"mass": 75, "abundance": 100, "interferences": ["ArCl (75)!!!", "CaAr (80?no)"],
         "notes": "CRITICAL: ArCl overlaps exactly! Must use CRC (He/H2) or reaction cell (O2)."},
    ],
    "se": [
        {"mass": 77, "abundance": 7.63, "interferences": ["ArCl1H (77)!!!", "Se (77)"],
         "notes": "ArClH severe if Cl present. Use Se-82 or CRC."},
        {"mass": 78, "abundance": 23.77, "interferences": ["Kr (78 isobaric!!!)", "Se (78)"],
         "notes": "Kr from Ar gas supply! High-purity Ar required or use Se-82."},
        {"mass": 80, "abundance": 49.61, "interferences": ["Ar2 (80!!!)", "Se (80)"],
         "notes": "Ar2 dimer background! Use cold plasma or CRC."},
        {"mass": 82, "abundance": 8.73, "interferences": ["BrH (82)", "Rb (82 isobaric)", "Kr (82 weak)", "Se (82)"],
         "notes": "RECOMMENDED: Cleanest Se isotope despite low abundance."},
    ],
    "sr": [
        {"mass": 84, "abundance": 0.56, "interferences": ["Kr (84 isobaric!!!)", "Sr (84)"],
         "notes": "Kr interference. Low abundance anyway."},
        {"mass": 86, "abundance": 9.86, "interferences": ["Kr (86 weak)", "Sr (86)"],
         "notes": "Moderate abundance; Kr from Ar gas."},
        {"mass": 87, "abundance": 7.00, "interferences": ["Rb (87 isobaric!!!)", "Sr (87)"],
         "notes": "Rb isobaric overlap."},
        {"mass": 88, "abundance": 82.58, "interferences": ["Sr (88)", "Kr (88 trace)"],
         "notes": "Most abundant; generally clean with high-purity Ar."},
    ],
    "cd": [
        {"mass": 106, "abundance": 27.64, "interferences": ["Zr (106 isobaric!!!)", "Cd (106)"],
         "notes": "Zr isobaric! Avoid in Zr-containing samples."},
        {"mass": 108, "abundance": 36.04, "interferences": ["Mo (108?) no Mo=98", "Zr (108?) no"],
         "notes": "Generally clean; recommended for many applications."},
        {"mass": 110, "abundance": 12.49, "interferences": ["Pd (110 isobaric!!!)", "Cd (110)"],
         "notes": "Pd isobaric! Avoid if Pd catalysts used."},
        {"mass": 111, "abundance": 12.80, "interferences": ["MoO (111)? Mo=96", "Cd (111)"],
         "notes": "Possible MoO if Mo present."},
        {"mass": 112, "abundance": 24.13, "interferences": ["Sn (112 isobaric!!!)", "Cd (112)"],
         "notes": "Sn isobaric! Common in environmental samples."},
        {"mass": 113, "abundance": 12.22, "interferences": ["In (113 isobaric!!!)", "Cd (113)"],
         "notes": "In isobaric overlap."},
        {"mass": 114, "abundance": 28.73, "interferences": ["Sn (114 isobaric!!!)", "Cd (114)"],
         "notes": "Sn isobaric again. Cd-108 or Cd-111 often best choices."},
    ],
    "sb": [
        {"mass": 121, "abundance": 57.21, "interferences": ["Sb (121)"], "notes": "Generally clean."},
        {"mass": 123, "abundance": 42.79, "interferences": ["Sb (123)"], "notes": "Generally clean. Use 121/123 ratio for confirmation."},
    ],
    "ba": [
        {"mass": 135, "abundance": 6.592, "interferences": ["Ba (135)"], "notes": "Low abundance."},
        {"mass": 137, "abundance": 11.232, "interferences": ["Ba (137)"], "notes": "Lanthanide oxide interferences possible."},
        {"mass": 138, "abundance": 71.70, "interferences": ["La (138 isobaric!!!)", "Xe (138!!!)", "Ba (138)"],
         "notes": "Xe from Ar gas! La isobaric. Need high-res or CRC."},
    ],
    "pb": [
        {"mass": 204, "abundance": 1.4, "interferences": ["Hg (204 isobaric!!!)", "Pb (204)"],
         "notes": "Hg isobaric! Very low abundance."},
        {"mass": 206, "abundance": 24.1, "interferences": ["Pb (206)"], "notes": "RECOMMENDED — clean, good abundance."},
        {"mass": 207, "abundance": 22.1, "interferences": ["Pb (207)"], "notes": "Clean; single isotope spike for IDMS."},
        {"mass": 208, "abundance": 52.4, "interferences": ["Pb (208)"], "notes": "Highest abundance; Rn decay product. Overlaps with Hg-204? no."},
    ],
    "u": [
        {"mass": 238, "abundance": 99.27, "interferences": ["U (238)", "Pt (238? no Pt=195)", "Hg2 (238??)", "PbO2 (238?)"],
         "notes": "Primary U isotope. 238U/235U ratio used for enrichment assessment."},
        {"mass": 235, "abundance": 0.72, "interferences": ["U (235)", "238U tail (Hg? no)"],
         "notes": "For nuclear forensics / enrichment. 233U as internal standard."},
    ],
}


# ── 常见多原子干扰数据库 (matrix → [polyatomic ion, mass]) ────
POLYATOMIC_DB: Dict[str, List[Dict[str, Any]]] = {
    "cl": [{"ion": "ArCl+", "masses": [75, 77], "affects": ["As", "Se"]},
           {"ion": "ClO+", "masses": [51, 53], "affects": ["Cr", "V"]},
           {"ion": "Cl2+", "masses": [70, 72, 74], "affects": ["Ge", "Se"]}],
    "s":  [{"ion": "SO+", "masses": [48], "affects": ["Ti"]},
           {"ion": "SO2+", "masses": [64, 66], "affects": ["Zn"]},
           {"ion": "ArS+", "masses": [72, 74], "affects": ["Ge", "Se"]}],
    "ca":[{"ion": "CaO+", "masses": [56, 57, 58, 59, 60], "affects": ["Fe", "Co", "Ni"]},
          {"ion": "CaOH+", "masses": [57, 58, 59, 61], "affects": ["Fe", "Ni"]},
          {"ion": "ArCa+", "masses": [80, 82, 83, 84, 86], "affects": ["Se", "Sr", "Kr"]}],
    "fe":[{"ion": "ArO+", "masses": [56], "affects": ["Fe itself - background"]},
          {"ion": "FeO+", "masses": [72, 74], "affects": ["Ge"]},
          {"ion": "FeH+", "masses": [57], "affects": ["Fe"]}],
    "na":[{"ion": "ArNa+", "masses": [63, 65], "affects": ["Cu", "Cu"]},
          {"ion": "Na2+", "masses": [46], "affects": []},  # charge state issues
          {"ion": "NaO+", "masses": [39], "affects": ["K"]}],
    "k": [{"ion": "ArK+", "masses": [79, 81], "affects": ["Br"]},
          {"ion": "KO+", "masses": [55, 59], "affects": ["Mn", "Co"]}],
    "ar":[{"ion": "Ar2+", "masses": [76, 78, 80, 82, 84, 86], "affects": ["Se", "Kr", "Sr"]},
          {"ion": "ArN+", "masses": [54], "affects": ["Cr"]},
          {"ion": "ArC+", "masses": [52], "affects": ["Cr"]},
          {"ion": "ArC1H+", "masses": [53], "affects": ["Cr"]}],
}


@ChemMCPManager.register_tool
class IcpMsIsotopeSelector(BaseTool):
    """
    ICP-MS同位素选择与多原子干扰规避工具。
    选择最优分析同位素，识别和规避多原子离子干扰。
    """
    __version__                = "0.1.0"
    name                       = "IcpMsIsotopeSelector"
    func_name                  = "select_isotope"
    description                = ("Select optimal isotope for ICP-MS analysis avoiding "
                                 "polyatomic interferences based on sample matrix and instrument type.")
    implementation_description = ("Uses built-in database of isotopic abundances, known polyatomic interferences "
                                 "(ArCl+, ArO+, oxides, doubly-charged ions), and recommends the best isotope "
                                 "with mitigation strategies for quadrupole, sector field, or ICP-MS/MS instruments.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["ICP-MS", "Mass Spectrometry", "Isotope Selection",
                                   "Polyatomic Interference", "Elemental Analysis"]
    required_envs              = []

    code_input_sig = [
        ("analyte_element",             "str",   "N/A",       "Element symbol (e.g., 'As', 'Se', 'Fe')."),
        ("known_matrix",                "list",  "[]",        "List of matrix elements/ions that may cause polyatomic interferences."),
        ("instrument_type",             "str",   "quad",      "Instrument type: 'quad' (Q-ICP-MS), 'sector' (SF-ICP-MS), or 'icpmsms'."),
        ("resolution_mode",             "str",   "low",       "Resolution mode: 'low' (~300), 'medium' (~4000), or 'high' (~10000)."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated: element [matrix_elem1,matrix_elem2,...] [instrument_type] [resolution]"),
    ]

    output_sig = [
        ("recommended_isotope",         "int",   "Recommended isotope mass number."),
        ("abundance",                   "float", "Natural abundance (%) of the recommended isotope."),
        ("polyatomic_interferences",     "list",  "List of potential polyatomic interferences at this isotope."),
        ("alternative_isotopes",        "list",  "List of alternative isotopes with pros/cons."),
        ("interference_reduction_tips",  "list",  "Specific recommendations to reduce/remove interferences."),
        ("risk_assessment",            "str",    "Overall risk assessment summary."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_element": "As",
                "known_matrix": ["Cl"],
                "instrument_type": "quad",
                "resolution_mode": "low",
            },
            "text_input": {"input_params": "As Cl quad low"},
            "output": {
                "recommended_isotope": 75,
                "risk_assessment": "",
            },
        },
        {
            "code_input": {
                "analyte_element": "Se",
                "known_matrix": [],
                "instrument_type": "quad",
            },
            "text_input": {"input_params": "Se [] quad"},
            "output": {
                "recommended_isotope": 82,
                "abundance": 8.73,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load isotope and polyatomic databases."""
        self.isotope_db = ICP_MS_ISOTOPE_DB
        self.poly_db = POLYATOMIC_DB

    def _run_base(
        self,
        analyte_element: str,
        known_matrix: Optional[List[str]] = None,
        instrument_type: str = "quad",
        resolution_mode: str = "low",
    ) -> Dict[str, Any]:
        """Core selection logic."""
        key = analyte_element.strip().lower()
        matrix_lower = [m.strip().lower() for m in (known_matrix or [])]
        inst = instrument_type.strip().lower()
        res = resolution_mode.strip().lower()

        if key not in self.isotope_db:
            available = ", ".join(sorted(self.isotope_db.keys()))
            raise ChemMCPError(f"Element '{analyte_element}' not found. Available: {available}")

        isotopes = self.isotope_db[key]

        # Score each isotope
        scored = []
        for iso in isotopes:
            mass = iso["mass"]
            ab = iso.get("abundance", 0)
            intfs = iso.get("interferences", [])
            notes = iso.get("notes", "")

            score = ab * 1.0  # Base score = abundance

            # Check which interferences are relevant to this matrix
            active_intfs = []
            for intf in intfs:
                # Normalize interference string
                intf_lower = intf.lower()
                for m in matrix_lower:
                    if m == "ar" or m in intf_lower or any(c.isalpha() and m in c.lower() for c in intf.split() if c.isalpha()):
                        active_intfs.append(intf)
                        break

            # Penalty for active interferences
            penalty = len(active_intfs) * 30

            # Critical interferences get extra penalty
            critical_keywords = ["!!!", "isobaric", "critical", "exact overlap"]
            for ai in active_intfs:
                if any(kw in ai.lower() for kw in critical_keywords):
                    penalty += 20

            # Bonus for clean isotopes
            if not active_intfs:
                score += 20

            # Instrument type adjustment
            if inst == "icpmsms":
                # MS/MS can eliminate most polyatomic interferences using reaction gases
                penalty = int(penalty * 0.2)  # Reduce penalty significantly
            elif inst == "sector" and res in ("medium", "high"):
                # High-resolution SF can resolve many overlaps
                penalty = int(penalty * 0.3)

            # Resolution adjustment for quad
            if inst == "quad" and res == "low":
                penalty += 5  # Quad LR can't resolve anything

            final_score = max(score - penalty, 0)
            scored.append((final_score, iso, active_intfs))

        # Sort by score descending
        scored.sort(key=lambda x: x[0], reverse=True)

        # Best isotope
        best_score, best_iso, best_intfs = scored[0]

        # Build alternatives
        alts = []
        for _, iso, intfs in scored[1:]:
            alts.append({
                "mass": iso["mass"],
                "abundance": iso.get("abundance", 0),
                "potential_interferences": iso.get("interferences", []),
                "active_matrix_interferences": intfs,
                "notes": iso.get("notes", ""),
            })

        # Generate reduction tips
        tips = []
        if best_intfs:
            tips.append(f"⚠️ Active interferences detected at m/z {best_iso['mass']}: {', '.join(best_intfs)}")
            if inst == "quad":
                tips.append("Use Collision/Reaction Cell (CRC) with He or H2 to attenuate polyatomic ions.")
                tips.append("Consider Kinetic Energy Discrimination (KED) mode.")
            elif inst == "sector":
                tips.append(f"Switch to medium/high resolution (R≥4000) to resolve spectral overlaps.")
            elif inst == "icpmsms":
                tips.append("Use appropriate reaction gas (e.g., O2 for As→AsO, H2 for Se removal).")

        # Matrix-specific tips
        if "cl" in matrix_lower:
            tips.append("High Cl matrix: consider dilution, matrix separation, or use O2 reaction gas to convert As→AsO+.")
        if "ca" in matrix_lower:
            tips.append("Ca matrix: expect CaO+/CaOH+ interferences on transition metals (m/z 56–61).")
        if "s" in matrix_lower:
            tips.append("Sulfur matrix: SO+/SO2+ can interfere with Zn, Ge isotopes.")

        # General recommendations
        if not best_intfs:
            tips.append("✅ No significant matrix-induced polyatomic interferences expected at selected isotope.")
        tips.append("Always verify with blank and matrix-matched calibration standards.")
        tips.append("Monitor secondary isotope(s) for ratio confirmation if available.")

        # Risk assessment
        if best_score < 20:
            risk = "HIGH RISK: Selected isotope has significant unresolved interferences. Strongly recommend alternative approach."
        elif best_score < 50:
            risk = "MODERATE RISK: Some interferences expected. Mitigation strategies recommended."
        elif best_intfs:
            risk = "LOW-MODERATE RISK: Minor interferences manageable with proper technique."
        else:
            risk = "LOW RISK: Clean isotope with minimal interference concerns."

        logger.info(f"ICP-MS isotope selection for {key.upper()}: m/z={best_iso['mass']}, "
                     f"abundance={best_iso.get('abundance','?')}%, score={best_score:.1f}")
        return {
            "recommended_isotope": int(best_iso["mass"]),
            "abundance": float(best_iso.get("abundance", 0)),
            "polyatomic_interferences": best_intfs,
            "alternative_isotopes": alts,
            "interference_reduction_tips": tips,
            "risk_assessment": risk,
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            elem = parts[0]
            matrix = []
            inst = "quad"
            res = "low"

            idx = 1
            while idx < len(parts):
                p = parts[idx]
                if p in ("quad", "sector", "icpmsms"):
                    inst = p
                elif p in ("low", "medium", "high"):
                    res = p
                elif "," in p or p.isalpha():
                    matrix.extend([x.strip() for x in p.split(",") if x.strip()])
                idx += 1

            return self._run_base(elem, matrix or None, inst, res)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
