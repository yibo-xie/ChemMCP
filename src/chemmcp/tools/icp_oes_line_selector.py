"""
ICP-OES Line Selector — ICP-OES分析线选择工具 (#326)

功能：
  根据待测元素、基体组成和检测限要求，推荐最优ICP-OES分析线，
  在灵敏度和抗干扰能力之间权衡。
"""

import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── ICP-OES 元素分析线数据库 ────────────────────────────────
# 数据来源: NIST Atomic Spectra Database, 各仪器厂商应用手册
# 每条线包含: 波长(nm), 灵敏度等级, 主要干扰, 推荐观测模式
ICP_OES_LINE_DB: Dict[str, List[Dict[str, Any]]] = {
    "ag": [
        {"wavelength_nm": 328.068, "sensitivity": "very_high", "interferences": ["VH (OH band)", "UH"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Most sensitive; check for V and U in matrix."},
        {"wavelength_nm": 338.289, "sensitivity": "high", "interferences": ["Ti", "W"],
         "viewing": "axial", "dl_ppb": 2, "notes": "Alternative line; less interference than 328 nm."},
        ],
    "al": [
        {"wavelength_nm": 396.152, "sensitivity": "high", "interferences": ["Ca", "Mo", "V"],
         "viewing": "radial", "dl_ppb": 5, "notes": "Most sensitive; Ca is major interferent."},
        {"wavelength_nm": 309.271, "sensitivity": "medium", "interferences": ["Mg", "V", "Fe"],
         "viewing": "radial", "dl_ppb": 10, "notes": "Less sensitive but fewer interferences."},
        {"wavelength_nm": 167.079, "sensitivity": "medium", "interferences": ["Cu"],
         "viewing": "radial", "dl_ppb": 8, "notes": "UV region; good for high-Ca samples."},
        ],
    "as": [
        {"wavelength_nm": 189.042, "sensitivity": "high", "interferences": ["Mn (distant)"],
         "viewing": "axial", "dl_ppb": 5, "notes": "Primary As line; vacuum/ purge required."},
        {"wavelength_nm": 193.696, "sensitivity": "medium", "interferences": ["Al", "C (molecular)"],
         "viewing": "axial", "dl_ppb": 15, "notes": "Al interference possible."},
        {"wavelength_nm": 197.197, "sensitivity": "low", "interferences": ["Pb (weak)", "Ar (background)"],
         "viewing": "radial", "dl_ppb": 40, "notes": "Cleaner spectral region."},
        ],
    "au": [
        {"wavelength_nm": 242.795, "sensitivity": "very_high", "interferences": ["Mn (weak)"],
         "viewing": "axial", "dl_ppb": 2, "notes": "Most sensitive Au line."},
        {"wavelength_nm": 267.595, "sensitivity": "high", "interferences": ["Co", "W"],
         "viewing": "axial", "dl_ppb": 4, "notes": "Alternative with moderate sensitivity."},
        ],
    "ba": [
        {"wavelength_nm": 455.403, "sensitivity": "very_high", "interferences": ["Ni (weak)"],
         "viewing": "axial", "dl_ppb": 0.05, "notes": "Very high sensitivity; best for trace Ba."},
        {"wavelength_nm": 493.409, "sensitivity": "high", "interferences": ["Fe (weak)"],
         "viewing": "axial", "dl_ppb": 0.1, "notes": "Alternative high-sensitivity line."},
        {"wavelength_nm": 233.527, "sensitivity": "medium", "interferences": ["Fe", "Co"],
         "viewing": "radial", "dl_ppb": 1, "notes": "Fewer interferences for complex matrices."},
        ],
    "ca": [
        {"wavelength_nm": 393.366, "sensitivity": "very_high", "interferences": ["Al (weak)"],
         "viewing": "axial", "dl_ppb": 0.01, "notes": "Ionic line; extremely sensitive."},
        {"wavelength_nm": 396.847, "sensitivity": "very_high", "interferences": ["Mo (weak)"],
         "viewing": "axial", "dl_ppb": 0.02, "notes": "Paired with 393.366 for confirmation."},
        {"wavelength_nm": 422.673, "sensitivity": "high", "interferences": ["CN band"],
         "viewing": "radial", "dl_ppb": 0.5, "notes": "Atomic line; useful for high-Ca samples."},
        ],
    "cd": [
        {"wavelength_nm": 226.502, "sensitivity": "very_high", "interferences": ["Ru (distant)", "Fe (weak)"],
         "viewing": "axial", "dl_ppb": 0.3, "notes": "Most sensitive Cd line."},
        {"wavelength_nm": 214.438, "sensitivity": "high", "interferences": ["Pt", "Ir (weak)"],
         "viewing": "axial", "dl_ppb": 0.5, "notes": "Alternative; Pt/Ir interference rare."},
        {"wavelength_nm": 228.802, "sensitivity": "high", "interferences": ["As (overlap!)", "Ni (close)"],
         "viewing": "axial", "dl_ppb": 0.4, "notes": "As 228.812 overlaps — avoid if As present."},
        ],
    "co": [
        {"wavelength_nm": 228.616, "sensitivity": "very_high", "interferences": ["Ti", "Ce"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Most sensitive Co line."},
        {"wavelength_nm": 238.892, "sensitivity": "high", "interferences": ["Fe (distant)", "W"],
         "viewing": "axial", "dl_ppb": 2, "notes": "Commonly used alternative."},
        ],
    "cr": [
        {"wavelength_nm": 205.552, "sensitivity": "very_high", "interferences": ["B (weak)", "Pt"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Most sensitive Cr line."},
        {"wavelength_nm": 267.716, "sensitivity": "high", "interferences": ["Pt", "Mn (weak)"],
         "viewing": "axial", "dl_ppb": 2, "notes": "UV-Vis region; widely used."},
        {"wavelength_nm": 357.869, "sensitivity": "high", "interferences": ["V", "Ga", "Mn"],
         "viewing": "radial", "dl_ppb": 3, "notes": "AAS-familiar wavelength; V interference common."},
        ],
    "cu": [
        {"wavelength_nm": 324.754, "sensitivity": "very_high", "interferences": ["Nb (minor)"],
         "viewing": "axial", "dl_ppb": 0.5, "notes": "Most sensitive Cu line; excellent S/N."},
        {"wavelength_nm": 327.396, "sensitivity": "very_high", "interferences": ["Nb (minor)"],
         "viewing": "axial", "dl_ppb": 0.6, "notes": "Paired line for ratio confirmation."},
        {"wavelength_nm": 213.598, "sensitivity": "medium", "interferences": ["Ni", "Pd (weak)"],
         "viewing": "radial", "dl_ppb": 5, "notes": "Fewer interferences for complex matrices."},
        ],
    "fe": [
        {"wavelength_nm": 238.204, "sensitivity": "very_high", "interferences": ["Pt", "Sn (weak)"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Most sensitive Fe line."},
        {"wavelength_nm": 259.940, "sensitivity": "very_high", "interferences": ["Mn (distant)", "Sb"],
         "viewing": "axial", "dl_ppb": 1.5, "notes": "Very commonly used Fe line."},
        {"wavelength_nm": 240.489, "sensitivity": "high", "interferences": ["Pt", "Rh"],
         "viewing": "radial", "dl_ppb": 3, "notes": "Moderate sensitivity with fewer issues."},
        ],
    "k": [
        {"wavelength_nm": 766.491, "sensitivity": "very_high", "interferences": ["CN band tail"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Doublet with 769.896 nm; very intense."},
        {"wavelength_nm": 769.896, "sensitivity": "very_high", "interferences": ["Ar (weak)"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Second K doublet component."},
        ],
    "mg": [
        {"wavelength_nm": 279.553, "sensitivity": "very_high", "interferences": ["Ti (weak)", "Zr"],
         "viewing": "axial", "dl_ppb": 0.01, "notes": "Triplet member; most sensitive Mg line."},
        {"wavelength_nm": 280.270, "sensitivity": "very_high", "interferences": ["Zr (weak)"],
         "viewing": "axial", "dl_ppb": 0.02, "notes": "Triplet member."},
        {"wavelength_nm": 285.213, "sensitivity": "high", "interferences": ["Na (ion) weak", "Fe (distant)"],
         "viewing": "radial", "dl_ppb": 0.1, "notes": "Atomic line; AAS-familiar wavelength."},
        ],
    "mn": [
        {"wavelength_nm": 257.610, "sensitivity": "very_high", "interferences": ["Fe (distant)", "Bi (weak)"],
         "viewing": "axial", "dl_ppb": 0.2, "notes": "Most sensitive Mn line; widely used."},
        {"wavelength_nm": 259.373, "sensitivity": "high", "interferences": ["Fe (distant)", "Sb"],
         "viewing": "axial", "dl_ppb": 0.4, "notes": "Alternative high-sensitivity line."},
        {"wavelength_nm": 293.306, "sensitivity": "medium", "interferences": ["V (weak)", "Cr"],
         "viewing": "radial", "dl_ppb": 1, "notes": "Fewer interferences."},
        ],
    "mo": [
        {"wavelength_nm": 202.030, "sensitivity": "very_high", "interferences": ["W (minor)"],
         "viewing": "axial", "dl_ppb": 2, "notes": "Most sensitive Mo line."},
        {"wavelength_nm": 204.593, "sensitivity": "high", "interferences": ["Al (weak)", "Rb"],
         "viewing": "axial", "dl_ppb": 4, "notes": "Alternative line."},
        {"wavelength_nm": 281.615, "sensitivity": "medium", "interferences": ["Mn (weak)", "Co"],
         "viewing": "radial", "dl_ppb": 10, "notes": "Lower sensitivity but cleaner."},
        ],
    "na": [
        {"wavelength_nm": 589.592, "sensitivity": "very_high", "interferences": ["OH band (weak)"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Na D1 line; very intense."},
        {"wavelength_nm": 588.995, "sensitivity": "very_high", "interferences": ["OH band (weak)"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Na D2 line; paired with D1."},
        ],
    "ni": [
        {"wavelength_nm": 231.604, "sensitivity": "very_high", "interferences": ["Co (close)", "Zr (weak)"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Most sensitive Ni line."},
        {"wavelength_nm": 221.647, "sensitivity": "high", "interferences": ["Co", "Si (distant)"],
         "viewing": "axial", "dl_ppb": 2, "notes": "Alternative; Co is main concern."},
        {"wavelength_nm": 341.476, "sensitivity": "medium", "interferences": ["Pt (weak)", "Zr"],
         "viewing": "radial", "dl_ppb": 5, "notes": "Fewer interferences."},
        ],
    "pb": [
        {"wavelength_nm": 220.353, "sensitivity": "very_high", "interferences": ["Pd (close)", "Sn (distant)"],
         "viewing": "axial", "dl_ppb": 2, "notes": "Most sensitive Pb line."},
        {"wavelength_nm": 217.000, "sensitivity": "high", "interferences": ["Co (weak)", "Ni (distant)"],
         "viewing": "axial", "dl_ppb": 4, "notes": "Alternative; widely used."},
        {"wavelength_nm": 283.306, "sensitivity": "high", "interferences": ["Fe (overlap!)", "Mg (weak)"],
         "viewing": "radial", "dl_ppb": 5, "notes": "Fe 283.305 nearly overlaps — avoid if Fe >> Pb."},
        ],
    "se": [
        {"wavelength_nm": 196.089, "sensitivity": "high", "interferences": ["Fe (close)", "Al (band)"],
         "viewing": "axial", "dl_ppb": 10, "notes": "Primary Se line; vacuum/purge needed."},
        {"wavelength_nm": 203.985, "sensitivity": "medium", "interferences": ["Cr (weak)"],
         "viewing": "axial", "dl_ppb": 25, "notes": "Alternative Se line."},
        ],
    "si": [
        {"wavelength_nm": 251.611, "sensitivity": "high", "interferences": ["V (close)", "Fe trace"],
         "viewing": "radial", "dl_ppb": 5, "notes": "Most sensitive Si line."},
        {"wavelength_nm": 212.412, "sensitivity": "medium", "interferences": ["Sb (weak)"],
         "viewing": "radial", "dl_ppb": 10, "notes": "UV region; fewer interferences."},
        ],
    "sr": [
        {"wavelength_nm": 407.771, "sensitivity": "very_high", "interferences": ["Cr (weak)", "CN band"],
         "viewing": "axial", "dl_ppb": 0.03, "notes": "Ionic line; most sensitive Sr line."},
        {"wavelength_nm": 421.552, "sensitivity": "high", "interferences": ["Rb (weak)", "Cr"],
         "viewing": "axial", "dl_ppb": 0.1, "notes": "Alternative ionic line."},
        ],
    "ti": [
        {"wavelength_nm": 334.941, "sensitivity": "very_high", "interferences": ["Zr (close)", "Ca (weak)"],
         "viewing": "axial", "dl_ppb": 0.3, "notes": "Most sensitive Ti line."},
        {"wavelength_nm": 336.121, "sensitivity": "high", "interferences": ["Nb (weak)", "Zr"],
         "viewing": "axial", "dl_ppb": 0.5, "notes": "Alternative high-sensitivity line."},
        {"wavelength_nm": 349.854, "sensitivity": "medium", "interferences": ["Cr (weak)", "W"],
         "viewing": "radial", "dl_ppb": 2, "notes": "Fewer interferences."},
        ],
    "v": [
        {"wavelength_nm": 292.402, "sensitivity": "very_high", "interferences": ["Fe (distant)", "Mo (weak)"],
         "viewing": "axial", "dl_ppb": 0.5, "notes": "Most sensitive V line."},
        {"wavelength_nm": 309.311, "sensitivity": "high", "interferences": ["OH band", "Mg (ion)"],
         "viewing": "axial", "dl_ppb": 1, "notes": "Alternative; OH band background possible."},
        {"wavelength_nm": 310.230, "sensitivity": "medium", "interferences": ["Ni (weak)", "Al"],
         "viewing": "radial", "dl_ppb": 3, "notes": "Moderate sensitivity."},
        ],
    "zn": [
        {"wavelength_nm": 206.200, "sensitivity": "very_high", "interferences": ["Cr (close)", "Bi (weak)"],
         "viewing": "axial", "dl_ppb": 0.3, "notes": "Most sensitive Zn line."},
        {"wavelength_nm": 213.856, "sensitivity": "high", "interferences": ["Ni (close)", "Cu (weak)", "Fe (3pm!)"],
         "viewing": "axial", "dl_ppb": 0.5, "notes": "Fe 213.859 only 3 pm away — use HR or alternate line."},
        {"wavelength_nm": 202.548, "sensitivity": "medium", "interferences": ["Mg (weak)"],
         "viewing": "radial", "dl_ppb": 2, "notes": "Fewer interferences."},
        ],
}


@ChemMCPManager.register_tool
class IcpOesLineSelector(BaseTool):
    """
    ICP-OES分析线选择工具。
    基于元素、基体和检测要求，推荐最优分析波长及备选线。
    """
    __version__                = "0.1.0"
    name                       = "IcpOesLineSelector"
    func_name                  = "select_analytical_line"
    description                = ("Select optimal analytical wavelength for ICP-OES analysis "
                                 "considering sensitivity vs. interference trade-offs.")
    implementation_description = ("Uses built-in database of analytical lines for ~30 elements with sensitivity ratings, "
                                 "known spectral interferences, detection limits, and recommended viewing modes.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["ICP-OES", "Emission Spectroscopy", "Analytical Chemistry",
                                   "Line Selection", "Elemental Analysis"]
    required_envs              = []

    code_input_sig = [
        ("analyte_element",             "str",   "N/A",       "Element symbol (e.g., 'Fe', 'Cd', 'Pb')."),
        ("matrix_composition",          "list",  "[]",        "List of matrix elements that may cause interference."),
        ("concentration_range",         "str",   "trace",     "Concentration range: 'ultra_trace', 'trace', 'minor', or 'major'."),
        ("prefer_viewing_mode",         "str",   "auto",      "Preferred viewing mode: 'axial', 'radial', or 'auto' (let tool decide)."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated: element [matrix_elem1,matrix_elem2,...] [conc_range] [view_mode]"),
    ]

    output_sig = [
        ("recommended_wavelength_nm",   "float", "Recommended primary analytical wavelength (nm)."),
        ("alternative_lines",           "list",  "List of alternative wavelengths with details."),
        ("sensitivity_rating",          "str",   "Sensitivity of the recommended line."),
        ("potential_interferences",     "list",  "Known interferences from the matrix at this line."),
        ("recommended_viewing_mode",    "str",   "Recommended viewing orientation: 'axial' or 'radial'."),
        ("selection_rationale",         "str",   "Explanation for why this line was selected."),
        ("detection_limit_estimate_ppb","float", "Estimated detection limit (ppb) for the selected line."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_element": "Cd",
                "matrix_composition": ["Fe"],
                "concentration_range": "trace",
            },
            "text_input": {"input_params": "Cd Fe trace"},
            "output": {
                "recommended_wavelength_nm": 214.438,
                "sensitivity_rating": "high",
            },
        },
        {
            "code_input": {
                "analyte_element": "Fe",
                "matrix_composition": [],
                "concentration_range": "ultra_trace",
            },
            "text_input": {"input_params": "Fe [] ultra_trace"},
            "output": {
                "recommended_wavelength_nm": 238.204,
                "sensitivity_rating": "very_high",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load ICP-OES line database."""
        self.db = ICP_OES_LINE_DB

    def _run_base(
        self,
        analyte_element: str,
        matrix_composition: Optional[List[str]] = None,
        concentration_range: str = "trace",
        prefer_viewing_mode: str = "auto",
    ) -> Dict[str, Any]:
        """Core selection logic."""
        key = analyte_element.strip().lower()
        matrix_lower = [m.strip().lower() for m in (matrix_composition or [])]

        if key not in self.db:
            available = ", ".join(sorted(self.db.keys()))
            raise ChemMCPError(f"Element '{analyte_element}' not found in ICP-OES database. Available: {available}")

        lines = self.db[key]
        if not lines:
            raise ChemMCPError(f"No analytical lines available for element '{key}'.")

        # Sensitivity priority based on concentration range
        sens_priority = {
            "ultra_trace": ["very_high", "high", "medium", "low"],
            "trace":       ["very_high", "high", "medium", "low"],
            "minor":       ["high", "medium", "very_high", "low"],  # avoid detector saturation
            "major":       ["medium", "low", "high", "very_high"],  # use less sensitive lines
        }
        pref_order = sens_priority.get(concentration_range.strip().lower(), sens_priority["trace"])

        # Score each line
        scored_lines = []
        for line in lines:
            score = 0
            sens = line["sensitivity"]

            # Base score from sensitivity ranking
            try:
                idx = pref_order.index(sens)
                score += (len(pref_order) - idx) * 10
            except ValueError:
                score += 1

            # Penalty for matrix interferences
            has_matrix_interference = False
            interfs_in_matrix = []
            for intf in line.get("interferences", []):
                intf_key = intf.lower().replace(" ", "").replace("(", "").replace(")", "").replace("+", "").replace("band", "").strip()
                # Check if any matrix element matches this interference
                for m in matrix_lower:
                    if m == intf_key or m in intf_key or intf_key in m.replace("^", "").replace("-", ""):
                        has_matrix_interference = True
                        interfs_in_matrix.append(intf)
                        break

            if has_matrix_interference:
                score -= 15  # Significant penalty for matrix interference
            else:
                score += 5  # Bonus for clean line

            # Viewing mode preference
            viewing = line.get("viewing", "axial")
            if prefer_viewing_mode != "auto" and viewing != prefer_viewing_mode:
                score -= 2  # Small penalty for non-preferred viewing mode

            scored_lines.append((score, line, interfs_in_matrix))

        # Sort by score descending
        scored_lines.sort(key=lambda x: x[0], reverse=True)

        # Select best line
        best_score, best_line, best_intfs = scored_lines[0]

        # Build alternatives list
        alt_lines = []
        for _, line, intfs in scored_lines[1:]:
            alt_lines.append({
                "wavelength_nm": line["wavelength_nm"],
                "sensitivity": line["sensitivity"],
                "interferences": line.get("interferences", []),
                "matrix_interferences_found": intfs,
                "viewing": line.get("viewing", "radial"),
                "detection_limit_ppb": line.get("dl_ppb", None),
            })

        # Determine viewing mode
        rec_viewing = best_line.get("viewing", "axial")
        if prefer_viewing_mode != "auto":
            rec_viewing = prefer_viewing_mode

        # Build rationale
        rationale_parts = []
        rationale_parts.append(
            f"Selected λ={best_line['wavelength_nm']:.3f} nm ({best_line['sensitivity']} sensitivity)."
        )
        if best_intfs:
            rationale_parts.append(
                f"⚠️ Potential matrix interferences detected: {', '.join(best_intfs)}. "
                f"Consider alternative lines if accuracy is compromised."
            )
        else:
            rationale_parts.append("No significant matrix interferences expected.")
        if concentration_range.lower() == "ultra_trace":
            rationale_parts.append("Prioritized maximum sensitivity for ultra-trace analysis.")
        elif concentration_range.lower() == "major":
            rationale_parts.append("Used lower-sensitivity line to avoid detector saturation.")

        logger.info(f"ICP-OES line selection for {key.upper()}: λ={best_line['wavelength_nm']} nm, "
                     f"score={best_score}, matrix_intfs={best_intfs}")
        return {
            "recommended_wavelength_nm": float(best_line["wavelength_nm"]),
            "alternative_lines": alt_lines,
            "sensitivity_rating": best_line["sensitivity"],
            "potential_interferences": best_intfs,
            "recommended_viewing_mode": rec_viewing,
            "selection_rationale": " ".join(rationale_parts),
            "detection_limit_estimate_ppb": float(best_line.get("dl_ppb", 0)),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            elem = parts[0]
            matrix = []
            conc = "trace"
            viewing = "auto"

            idx = 1
            while idx < len(parts):
                p = parts[idx]
                if p in ("ultra_trace", "trace", "minor", "major"):
                    conc = p
                elif p in ("axial", "radial"):
                    viewing = p
                elif "," in p or p.isalpha():
                    matrix.extend([x.strip() for x in p.split(",") if x.strip()])
                idx += 1

            return self._run_base(elem, matrix or None, conc, viewing)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
