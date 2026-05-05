"""
Internal Standard Selector — 内标物选择建议
基于结构相似性、响应因子、色谱行为等推荐合适的内标物
"""
import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见分析类型及其推荐内标数据库 ────────────────────────────────
# 格式: (内标名称, 适用范围, 结构特点, CAS号(可选), 分子量, 备注)
IS_DATABASE: Dict[str, List[dict]] = {
    "gc": [
        {"name": "n-Tetradecane (C14)", "mw": 198.39, "bp_C": 254, "polarity": "non-polar",
         "best_for": "Non-polar to medium-polar compounds; hydrocarbons, fatty acid methyl esters.",
         "notes": "Very stable; well-resolved from most analytes on non-polar columns."},
        {"name": "n-Dodecane (C12)", "mw": 170.33, "bp_C": 216, "polarity": "non-polar",
         "best_for": "Volatile to semi-volatile non-polar compounds.",
         "notes": "Good for GC-MS; minimal fragmentation interference."},
        {"name": "5α-Cholestane", "mw": 372.67, "bp_C": "decomposes", "polarity": "non-polar",
         "best_for": "Steroids, lipids, high-MW non-polar compounds.",
         "notes": "Excellent for lipid analysis by GC-FID/MS."},
        {"name": "Methyl undecanoate (C11:0 FAME IS)", "mw": 200.32, "bp_C": 243, "polarity": "medium",
         "best_for": "Fatty acid methyl ester (FAME) analysis.",
         "notes": "Odd-chain FAME not found naturally in most samples."},
        {"name": "Biphenyl", "mw": 154.21, "bp_C": 255, "polarity": "low-polarity",
         "best_for": "Aromatic compounds; pesticide residue analysis.",
         "notes": "Commonly used in EPA methods."},
        {"name": "1-Bromonaphthalene", "mw": 207.07, "bp_C": 281, "polarity": "medium",
         "best_for": "Halogenated compounds; ECD detection.",
         "notes": "Strong ECD response; good for organochlorine analysis."},
        {"name": "Isotopically labeled analog", "mw": "N/A", "bp_C": "≈ analyte", "polarity": "identical",
         "best_for": "ALL analytes — gold standard when available.",
         "notes": "Best choice: d₃-, d₅-, ¹³C-labeled version of analyte. Compensates for matrix effects and extraction loss."},
        {"name": "Camphor", "mw": 152.23, "bp_C": 204, "polarity": "medium",
         "best_for": "Oxygen-containing terpenes; chiral separations.",
         "notes": "Distinctive MS pattern; well-separated from many volatiles."},
    ],
    "lc_uv": [
        {"name": "Caffeine", "mw": 194.19, "lambda_max_nm": 273, "polarity": "medium",
         "best_for": "Medium-polar compounds at UV 220-280 nm.",
         "notes": "Good UV chromophore; stable; inexpensive."},
        {"name": "p-Hydroxybenzoic acid n-butyl ester (paraben)", "mw": 194.23, "lambda_max_nm": 257,
         "polarity": "medium", "best_for": "Phenolic compounds; acidic analytes.",
         "notes": "Well-retained on reversed-phase columns."},
        {"name": "4-Nitroaniline", "mw": 138.12, "lambda_max_nm": 380, "polarity": "medium",
         "best_for": "Visible-range detection; aromatic amines.",
         "notes": "Strong visible absorbance; distinct from most analytes."},
        {"name": "Phenacetin", "mw": 179.22, "lambda_max_nm": 244, "polarity": "medium-low",
         "best_for": "Drug/pharmaceutical analysis by HPLC-UV.",
         "notes": "Classic pharmaceutical IS; well-characterized."},
        {"name": "Benzoic acid", "mw": 122.12, "lambda_max_nm": 230, "polarity": "medium",
         "best_for": "Acidic compounds; food additive analysis.",
         "notes": "Simple structure; sharp peak; cheap."},
        {"name": "Isotopically labeled analog", "mw": "N/A", "lambda_max_nm": "≈ analyte",
         "polarity": "identical", "best_for": "ALL analytes — gold standard.",
         "notes": "LC-MS/MS ideal IS; compensates for ionization suppression/enhancement."},
    ],
    "lc_ms": [
        {"name": "Chloramphenicol-d₅", "mw": 326.56, "polarity": "medium",
         "best_for": "Small molecule drugs; metabolites; broad utility LC-MS IS.",
         "notes": "Stable isotope labeled; good ionization in both ESI+/-."},
        {"name": "Tolbutamide", "mw": 270.35, "polarity": "medium",
         "best_for": "ESI+ mode pharmaceutical analysis.",
         "notes": "Good proton affinity; commonly used in bioanalysis."},
        {"name": "Diclofenac-d₄", "mw": 300.11, "polarity": "acidic",
         "best_for": "Acidic drugs; NSAID-type compounds; ESI- mode.",
         "notes": "Stable label; excellent for acidic analytes."},
        {"name": "Propranolol-d₇", "mw": 299.41, "polarity": "basic",
         "best_for": "Basic drugs; amines; β-blockers; ESI+ mode.",
         "notes": "Labeled basic IS; matches ionization of basic analytes."},
        {"name": "Reserpine", "mw": 608.68, "polarity": "basic",
         "best_for": "LC-MS system suitability test; tuning standard.",
         "notes": "Often used as system performance IS; high MW reference."},
        {"name": "Isotopically labeled analog of analyte", "mw": "N/A", "polarity": "identical",
         "best_for": "ALL analytes — ALWAYS first choice for LC-MS/MS.",
         "notes": "Compensates for extraction recovery AND matrix effects. Expensive but worth it."},
    ],
    "gc_ms": [
        {"name": "Deuterated PAH standards (e.g., Phenanthrene-d₁₀)", "mw": 188.29, "polarity": "non-polar",
         "best_for": "PAH analysis; environmental samples.",
         "notes": "EPA Method 8270 style IS mixture available commercially."},
        {"name": "Terphenyl-d₁₄", "mw": 242.34, "polarity": "non-polar",
         "best_for": "High-temperature GC-MS; semivolatile organics.",
         "notes": "Thermal stability; late eluting — avoids coelution."},
        {"name": "PCB 209 (4,4'-Dichlorobiphenyl)", "mw": 360.87, "polarity": "non-polar",
         "best_for": "PCB analysis; organochlorine pesticides.",
         "notes": "Not found environmentally; excellent GC-ECD/MS IS."},
        {"name": "Isotopically labeled analog", "mw": "N/A", "polarity": "identical",
         "best_for": "ALL analytes — best possible choice.",
         "notes": "GC-MS/MS quantification benefits greatly from labeled IS."},
    ],
    "ic": [
        {"name": "Lithium (Li⁺)", "mw": 6.94, "charge": "+1",
         "best_for": "Cation analysis (Na⁺, K⁺, NH₄⁺, etc.) by IC.",
         "notes": "Well-separated from common cations; low background."},
        {"name": "2,2-Bipyridine", "mw": 156.18, "charge": "neutral→complex",
         "best_for": "Transition metal IC analysis.",
         "notes": "Forms complexes with metals; useful for post-column derivatization."},
        {"name": "Nitrate (NO₃⁻) or Bromide (Br⁻)", "mw": 62.0 / 79.9, "charge": "-1",
         "best_for": "Anion analysis (when not target analytes).",
         "notes": "Choose an anion not present in sample that elutes near region of interest."},
    ],
}


@ChemMCPManager.register_tool
class InternalStandardSelector(BaseTool):
    """
    内标物选择工具：根据待测物信息、分析方法、基质类型，
    推荐最合适的内标物，并给出选择理由和注意事项。
    """
    __version__ = "0.1.0"
    name = "InternalStandardSelector"
    func_name = "select_internal_standard"
    description = "Recommend suitable internal standards based on analyte properties, analytical method (GC/LC/IC), and structural similarity considerations."
    implementation_description = "Uses a built-in database of common internal standards organized by analytical technique (GC, LC-UV, LC-MS, GC-MS, IC). Ranks candidates by structural similarity, chromatographic behavior, response factor compatibility, and practical availability."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Internal Standard", "Chromatography", "GC", "LC", "Method Development", "Quantitative Analysis"]
    required_envs = []

    code_input_sig = [
        ("analyte", "str", "", "Name or description of the analyte (e.g., 'caffeine', 'PAH', 'pesticide')."),
        ("analysis_type", "str", "lc_ms", "Analysis type: 'gc', 'gc_ms', 'lc_uv', 'lc_ms', 'ic'."),
        ("analyte_mw", "float", "0", "Analyte molecular weight (0 if unknown)."),
        ("analyte_polarity", "str", "unknown", "Analyte polarity: 'non-polar', 'low', 'medium', 'high', 'acidic', 'basic', 'unknown'."),
        ("detection_method", "str", "", "Specific detector: 'fid', 'ecd', 'ms', 'uv', 'fluorescence', 'cad', ''."),
        ("matrix_type", "str", "general", "Sample matrix: 'biological', 'environmental', 'food', 'pharmaceutical', 'water', 'soil', 'general'."),
        ("has_isotope_labeled", "bool", "False", "Whether isotopically labeled IS of analyte is commercially available."),
        ("chromatography_details", "str", "", "Extra details about column, mobile phase, etc."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., 'acetophenone gc_ms' or 'dopamine lc_ms basic'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with recommended internal standards, ranking scores, selection rationale, and practical notes."),
    ]

    examples = [
        {
            "code_input": {
                "analyte": "acetophenone",
                "analysis_type": "gc",
                "analyte_mw": 120.15,
                "analyte_polarity": "medium",
                "detection_method": "fid",
                "matrix_type": "environmental",
                "has_isotope_labeled": False,
                "chromatography_details": "",
            },
            "text_input": {
                "input_params": "acetophenone gc",
            },
            "output": {
                "result": {
                    "mode": "recommendation",
                    "analyte": "acetophenone",
                    "note": "IS recommendations with scoring.",
                }
            }
        },
        {
            "code_input": {
                "analyte": "dopamine",
                "analysis_type": "lc_ms",
                "analyte_mw": 153.18,
                "analyte_polarity": "basic",
                "detection_method": "ms",
                "matrix_type": "biological",
                "has_isotope_labeled": True,
                "chromatography_details": "HILIC column, ESI+ mode",
            },
            "text_input": {
                "input_params": "dopamine lc_ms basic biological",
            },
            "output": {
                "result": {
                    "mode": "recommendation",
                    "analyte": "dopamine",
                    "note": "Dopamine-specific IS recommendations.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _get_candidates(self, analysis_type: str) -> List[dict]:
        key = analysis_type.lower().strip()
        mapping = {
            "gc": "gc", "gc_fid": "gc", "gc_ecd": "gc",
            "gc_ms": "gc_ms",
            "lc_uv": "lc_uv", "lc_uvis": "lc_uv", "hplc_uv": "lc_uv",
            "lc_ms": "lc_ms", "lc_msms": "lc_ms", "uPLC_ms": "lc_ms",
            "ic": "ic", "ion_chromatography": "ic",
        }
        db_key = mapping.get(key, "lc_ms")
        if db_key not in IS_DATABASE:
            raise ChemMCPError(
                f"Unknown analysis type '{analysis_type}'. "
                f"Supported: gc, gc_ms, lc_uv, lc_ms, ic."
            )
        return IS_DATABASE[db_key]

    def _score_candidate(self, candidate: dict, analyte: str, mw: float,
                         polarity: str, detection: str, matrix: str,
                         has_labeled: bool) -> dict:
        """Score an IS candidate from 0-100."""
        score = 50  # base score

        name_lower = candidate["name"].lower()

        # Gold standard: isotopically labeled
        if "isotope" in name_lower or "labeled" in name_lower:
            if has_labeled:
                score += 40
            else:
                score += 25  # still recommend but note availability

        # Polarity match
        cand_pol = candidate.get("polarity", "").lower()
        pol = polarity.lower()
        if pol != "unknown" and cand_pol != "identical":
            if pol == cand_pol or (pol in cand_pol or cand_pol in pol):
                score += 15
            elif (pol in ("acidic", "basic") and cand_pol == "medium") or \
                 (cand_pol in ("acidic", "basic") and pol == "medium"):
                score += 8
            elif pol == "non-polar" and cand_pol in ("non-polar", "low-polarity"):
                score += 12
            elif pol == "non-polar" and cand_pol == "medium":
                score += 3

        # MW proximity (within 20-200% range preferred)
        if mw > 0:
            cand_mw = candidate.get("mw", 0)
            if isinstance(cand_mw, (int, float)) and cand_mw > 0:
                ratio = mw / cand_mw
                if 0.5 <= ratio <= 2.0:
                    score += 10
                elif 0.2 <= ratio <= 5.0:
                    score += 5

        # Detection compatibility
        det = detection.lower() if detection else ""
        if "ecd" in det and "bromo" in name_lower:
            score += 15  # brominated compounds great for ECD
        if "uv" in det and "lambda_max" in candidate:
            score += 5  # has UV data

        # Matrix considerations
        if matrix == "biological" and "isotope" in name_lower:
            score += 10  # especially important for complex matrices
        if matrix == "environmental" and ("cholestane" in name_lower or "terphenyl" in name_lower):
            score += 5

        return min(score, 100)

    def _run_base(self, analyte: str = "", analysis_type: str = "lc_ms",
                  analyte_mw: float = 0.0, analyte_polarity: str = "unknown",
                  detection_method: str = "", matrix_type: str = "general",
                  has_isotope_labeled: bool = False,
                  chromatography_details: str = "") -> dict:

        if not analyte.strip():
            raise ChemMCPError("Analyte name must be provided.")

        candidates = self._get_candidates(analysis_type)

        # Score all candidates
        scored = []
        for cand in candidates:
            s = self._score_candidate(
                cand, analyte, analyte_mw, analyte_polarity,
                detection_method, matrix_type, has_isotope_labeled
            )
            scored.append({**cand, "score": s})

        scored.sort(key=lambda x: x["score"], reverse=True)

        top = scored[0]
        return {"result": {
            "mode": "recommendation",
            "analyte": analyte,
            "analysis_type": analysis_type,
            "analyte_mw": analyte_mw,
            "analyte_polarity": analyte_polarity,
            "detection_method": detection_method or "not specified",
            "matrix_type": matrix_type,
            "top_recommendation": {
                "name": top["name"],
                "score": top["score"],
                "mw": top.get("mw", "N/A"),
                "polarity": top.get("polarity", "N/A"),
                "best_for": top.get("best_for", ""),
                "notes": top.get("notes", ""),
            },
            "all_candidates_ranked": [
                {"name": c["name"], "score": c["score"], "reason": c.get("best_for", "")}
                for c in scored[:6]
            ],
            "selection_criteria_applied": [
                "Structural & chemical similarity to analyte",
                "Chromatographic separation (no co-elution)",
                "Similar response factor / detection sensitivity",
                "Stability under analysis conditions",
                "Absence from sample matrix",
                "Commercial availability & cost",
            ],
            "important_notes": self._generate_notes(analyte, analysis_type, has_isotope_labeled),
            "validation_checklist": [
                "✓ Confirm IS does NOT exist in unspiked blank samples",
                "✓ Verify baseline resolution from nearest analyte peak (Rs > 1.5)",
                "✓ Check IS stability in sample solvent and autosampler (≥24h)",
                "✓ Verify linear response over expected concentration range",
                "✓ Test precision with IS (RSD < 3% at mid-level QC)",
                "✓ Document retention time and response factor relative to analyte",
            ],
        }}

    @staticmethod
    def _generate_notes(analyte: str, analysis_type: str, has_labeled: bool) -> List[str]:
        notes = []
        if has_labeled:
            notes.append("🏆 Isotopically labeled IS is available — this should be your FIRST choice.")
            notes.append("   Labeled IS compensates for: extraction loss, matrix effects, injection variability.")
        else:
            notes.append("⚠ No isotopically labeled IS specified. Structural analog is second-best option.")
            notes.append("   Consider custom synthesis or commercial sourcing of labeled compound.")

        if analysis_type in ("lc_ms", "lc_msms"):
            notes.append("For LC-MS: monitor for in-source fragmentation or adduct formation overlap.")
        if analysis_type in ("gc", "gc_ms"):
            notes.append("For GC: ensure IS thermal stability at max oven temperature.")

        notes.append(f"Always run matrix-matched calibration with {analyte} + selected IS.")
        return notes

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            analyte = parts[0] if parts else ""
            atype = parts[1] if len(parts) > 1 else "lc_ms"
            polarity = parts[2] if len(parts) > 2 else "unknown"
            matrix = parts[3] if len(parts) > 3 else "general"
            return self._run_base(analyte=analyte, analysis_type=atype,
                                   analyte_polarity=polarity, matrix_type=matrix)
        except IndexError:
            raise ChemMCPError(f"Failed to parse text input '{input_params}'. Need at least: analyte analysis_type.")
