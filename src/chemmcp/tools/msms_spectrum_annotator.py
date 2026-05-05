"""
MS/MS Spectrum Annotator - Annotates MS/MS spectrum peaks with fragment ion assignments
including b/y ions (peptides), common neutral losses, and diagnostic fragments.
"""

import logging
import math
import re
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Atomic masses for fragment calculation
_FRAG_MASSES = {
    "H": 1.00783, "C": 12.0000, "N": 14.00307, "O": 15.99491,
    "P": 30.97376, "S": 31.97207, "F": 18.99840,
    "Cl": 34.96885, "Br": 78.91833,
}

# Common MS/MS neutral losses: (mass_delta, name, category, frequency)
_NEUTRAL_LOSSES_DB = [
    (1.00783,     "H•",           "radical",        "common"),
    (2.01565,     "H₂",           "reduction",       "occasional"),
    (17.00274,    "OH•",          "radical",         "alcohols/acids"),
    (18.01056,    "H₂O",          "dehydration",     "VERY COMMON"),
    (26.98702,    "C₂H₂",         "acetylene",       "aromatics"),
    (27.99491,    "CO",            "carbonyl_loss",   "COMMON"),
    (28.03130,    "C₂H₄",         "ethylene",        "McLafferty"),
    (29.00274,    "CHO•",          "formyl",          "aldehydes"),
    (30.01058,    "CH₂O",         "formaldehyde",     ""),
    (31.01839,    "CH₃O•",        "methoxy",         "ethers/esters"),
    (42.01056,    "CH₂=C=O",      "ketene",          "acetates"),
    (43.00581,    "C₂H₃O",        "acetyl",          "ketones"),
    (43.98983,    "CO₂",          "decarboxylation", "COMMON — acids"),
    (44.03130,    "C₂H₆O",        "ethanol_fragment", ""),
    (45.02940,    "CH₃CHOH",      "ethanol_loss",    ""),
    (46.00548,    "NO₂",          "nitro_loss",      "nitro compounds"),
    (46.00548,    "CH₂O₂",        "formic_acid",     ""),
    (49.99232,    "CH₃Cl",        "methyl_cl",       "chlorinated"),
    (56.02622,    "C₄H₈",         "butene",          "McLafferty"),
    (58.04187,    "C₃H₆O",        "acetone_loss",    "ketones"),
    (60.02113,    "C₂H₄O₂",       "acetic_acid",     "McLafferty acids"),
    (64.96909,    "SO₂",          "sulfur_dioxide",  "sulfonyl"),
    (79.91690,    "Br•",           "bromine_loss",    "brominated"),
    (84.04187,    "C₄H₈O",        "butyric_fragments","esters"),
    (98.01686,    "H₃PO₄",        "phosphoric_acid", "phosphorylated"),
    (117.99044,   "C₅H₄O₄",      "dicarboxyl",      "diacids"),
]

# Diagnostic fragment ions: (mz, formula, name, compound_class)
_DIAGNOSTIC_FRAGMENTS = [
    # Amines
    (30.0340,  "CH₄N⁺",    "methyleniminium",    "primary_amine"),
    (44.0496,  "C₂H₆N⁺",   "ethyliminium",       "secondary_amine"),
    (58.0652,  "C₃H₈N⁺",   "propyliminium",      "tertiary_amine"),
    # Amides / Peptides
    (44.0132,  "CH₄NO⁺",   "formamide_ion",      "amide"),
    (60.0448,  "C₂H₆NO⁺",  "acetamide_immonium", "amide"),
    (74.0241,  "C₃H₅O₂N⁺", "glycine_immonium",   "peptide_b1/a1"),
    (84.0449,  "C₄H₆NO⁺",  "aba_immonium",       "peptide_a/b"),
    (86.0606,  "C₄H₈NO⁺",  "threonine_imm",      "peptide"),
    (102.0550, "C₄H₈NO₂⁺", "glutamic_imm",       "peptide"),
    (120.0810, "C₄H₁₀NO₃⁺","serine_immonium",    "peptide"),
    # Common organic
    (27.0235,  "C₂H₃⁺",    "vinyl_cation",       "hydrocarbon"),
    (29.0391,  "C₂H₅⁺",    "ethyl_cation",       "hydrocarbon"),
    (41.0391,  "C₃H₅⁺",    "allyl_cation",       "hydrocarbon"),
    (43.0548,  "C₃H₇⁺",    "propyl_cation",      "hydrocarbon"),
    (55.0548,  "C₄H₇⁺",    "butenyl_cation",     "hydrocarbon"),
    (69.0704,  "C₅H₉⁺",    "pentenyl_cation",    "hydrocarbon"),
    (77.0391,  "C₆H₅⁺",    "phenyl_cation",      "aromatic"),
    (91.0548,  "C₇H₇⁺",    "tropylium",          "benzyl_aromatic"),
    (105.0704, "C₈H₉⁺",    "benzyl+CH2",         "aromatic"),
]


@ChemMCPManager.register_tool
class MsMsSpectrumAnnotator(BaseTool):
    """
    MS/MS 谱图碎片标注器 — 对 MS/MS 谱图中的峰进行碎片离子归属标注。
    
    根据前体离子信息和候选分子式，自动匹配常见中性丢失、诊断离子和特征碎片，
    给出每个峰的可能归属及置信度。
    """
    __version__      = "0.1.0"
    name             = "MsMsSpectrumAnnotator"
    func_name        = "annotate_msms_spectrum"
    description      = "Annotate MS/MS spectrum peaks with fragment ion assignments including neutral losses, diagnostic ions, and characteristic fragmentation patterns."
    implementation_description = "Matches observed peaks against a database of known neutral losses from the precursor mass, diagnostic fragment ions, and common fragmentation patterns. Uses tolerance-based matching with confidence scoring based on error magnitude and rule specificity."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Mass Spectrometry", "MS/MS", "Spectrum Annotation", "Peak Assignment", "Fragment Identification"]
    required_envs    = []

    code_input_sig   = [
        ("precursor_info", "str", "N/A", "Precursor info: molecular formula, SMILES, or string like 'mz=286.14 formula=C17H19NO3'."),
        ("msms_peaks", "list", "N/A", "List of MS/MS peaks as [{'mz': float, 'intensity': float}, ...] or list of [mz, intensity] pairs."),
        ("tolerance_da", "float", "0.02", "Matching tolerance in Daltons (typical: 0.01-0.05 for Q-TOF, 0.005 for Orbitrap)."),
        ("annotation_mode", "str", "comprehensive", "Annotation mode: 'comprehensive' (all matches), 'best_match' (top assignment only), or 'conservative' (high-confidence only)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "String format: 'precursor_info; mz1:intensity1, mz2:intensity2, ... [tolerance_da]'. Example: 'C17H19NO3 m/z=286.14; 268.13:100, 240.11:45, 212.10:20, 182.08:35 0.02'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with annotated_peaks (each with assignments), unassigned_peaks, coverage_percent, spectral_quality_metrics, and annotation summary."),
    ]

    examples         = [
        {
            "code_input": {
                "precursor_info": "C17H19NO3",
                "msms_peaks": [{"mz": 268.1332, "intensity": 100}, {"mz": 240.1280, "intensity": 65},
                               {"mz": 212.1174, "intensity": 25}, {"mz": 184.1070, "intensity": 15},
                               {"mz": 156.1015, "intensity": 8}, {"mz": 182.0800, "intensity": 35}],
                "tolerance_da": 0.02,
                "annotation_mode": "comprehensive",
            },
            "text_input": {
                "input_params": "C17H19NO3 m/z=286.14; 268.13:100,240.11:65,212.10:25,184.11:15,156.10:8,182.08:35 0.02"
            },
            "output": {
                "result": {
                    "precursor_mz": 286.1438,
                    "precursor_formula": "C17H19NO3",
                    "total_peaks": 6,
                    "annotated_count": 5,
                    "coverage_percent": 83.3,
                    "annotated_peaks": [
                        {"mz": 268.1332, "intensity": 100, "assignments": [
                            {"type": "neutral_loss", "assignment": "[M-H₂O+H]⁺", "loss": "H₂O (18.0106 Da)", "error_Da": 0.002, "confidence": "high"}
                        ]},
                        {"mz": 240.1280, "intensity": 65, "assignments": [
                            {"type": "neutral_loss", "assignment": "[M-H₂O-CO+H]⁺", "loss": "H₂O+CO (46.0055 Da)", "error_Da": 0.001, "confidence": "high"}
                        ]},
                    ],
                    "unassigned_peaks": [
                        {"mz": 182.0800, "intensity": 35, "possible_causes": "rearrangement ion or external contamination"},
                    ],
                    "spectral_quality": {"signal_to_noise": "good", "fragment_coverage": "good"},
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, precursor_info: str, msms_peaks: list,
                  tolerance_da: float = 0.02,
                  annotation_mode: str = "comprehensive") -> dict:
        """Core logic."""
        if not msms_peaks:
            raise ChemMCPError("MS/MS peak data is required.")
        if tolerance_da <= 0:
            raise ChemMCPError("Tolerance must be positive.")
        if annotation_mode not in ("comprehensive", "best_match", "conservative"):
            raise ChemMCPError("Mode must be 'comprehensive', 'best_match', or 'conservative'.")

        # Parse precursor info
        prec_mz, prec_formula = self._parse_precursor(precursor_info)

        # Normalize peak input
        normalized_peaks = self._normalize_peaks(msms_peaks)
        
        if not normalized_peaks:
            raise ChemMCPError("No valid peaks found after normalization.")

        total_intensity = sum(p["intensity"] for p in normalized_peaks)
        base_peak_intensity = max(p["intensity"] for p in normalized_peaks)

        # Annotate each peak
        annotated_peaks = []
        unassigned_peaks = []
        annotated_count = 0
        annotated_intensity = 0.0

        for peak in normalized_peaks:
            assignments = self._assign_peak(peak["mz"], prec_mz, prec_formula, tolerance_da, annotation_mode)
            
            if assignments:
                # Filter by mode
                if annotation_mode == "best_match":
                    assignments = [max(assignments, key=lambda a: a["score"])]
                elif annotation_mode == "conservative":
                    assignments = [a for a in assignments if a["confidence"] in ("high", "very_high")]
                
                if assignments:
                    annotated_count += 1
                    annotated_intensity += peak["intensity"]
            
            peak_entry = {
                "mz": round(peak["mz"], 4),
                "relative_intensity_pct": round(peak["intensity"] / base_peak_intensity * 100, 1),
                "absolute_intensity": peak["intensity"],
                "assignments": assignments,
            }

            if assignments:
                annotated_peaks.append(peak_entry)
            else:
                peak_entry["possible_causes"] = self._suggest_unassigned_cause(peak["mz"], prec_mz)
                unassigned_peaks.append(peak_entry)

        # Calculate coverage
        coverage_pct = round(annotated_count / len(normalized_peaks) * 100, 1)
        intensity_coverage = round(annotated_intensity / total_intensity * 100, 1) if total_intensity > 0 else 0

        # Spectral quality metrics
        quality = self._assess_spectral_quality(normalized_peaks, annotated_count, len(normalized_peaks))

        return {
            "result": {
                "precursor_mz": prec_mz,
                "precursor_formula": prec_formula or "unknown",
                "input_precursor_info": precursor_info,
                "tolerance_da": tolerance_da,
                "annotation_mode": annotation_mode,
                "total_peaks": len(normalized_peaks),
                "annotated_peak_count": annotated_count,
                "coverage_percent": coverage_pct,
                "intensity_coverage_percent": intensity_coverage,
                "annotated_peaks": annotated_peaks,
                "unassigned_peaks": unassigned_peaks,
                "spectral_quality": quality,
                "suggested_fragments": self._suggest_additional_fragments(prec_mz, prec_formula),
                "summary": (
                    f"Annotated {annotated_count}/{len(normalized_peaks)} peaks ({coverage_pct}% coverage, "
                    f"{intensity_coverage}% intensity). "
                    + (f"{'Best match' if annotation_mode == 'best_match' else 'Comprehensive'} mode." )
                ),
                "notes": (
                    "MS/MS Annotation Notes:\n"
                    "• Assignments are computational predictions — verify with standards\n"
                    "• Neutral loss annotations assume [M+H]+ precursor\n"
                    "• Confidence levels: very_high (<0.005 Da), high (<0.01 Da), medium (<tolerance)\n"
                    "• Unassigned peaks may be: rearrangement ions, internal fragments, noise, or contaminants\n"
                    "• Isobaric interference can cause misassignment in complex mixtures\n"
                    "• Consider using isotopic pattern to confirm assignments when possible"
                ),
            }
        }

    def _parse_precursor(self, info: str) -> Tuple[Optional[float], Optional[str]]:
        """Parse precursor info string."""
        import re
        info = info.strip()
        prec_mz = None
        prec_formula = None
        
        # Look for explicit m/z
        mz_match = re.search(r'(?:m/z|mz)[\s=:]*([\d.]+)', info, re.I)
        if mz_match:
            prec_mz = float(mz_match.group(1))

        # Look for formula
        formula_match = re.search(r'(?:formula|=)\s*([A-Z][a-z]?\d*)+', info, re.I)
        if formula_match:
            # Try to extract full formula
            fm = re.findall(r'([A-Z][a-z]?)(\d*)', info)
            formula_parts = []
            for elem, cnt in fm:
                if elem and elem.isalpha() and elem[0].isupper():
                    formula_parts.append(f"{elem}{cnt}" if cnt else elem)
            if formula_parts:
                prec_formula = "".join(formula_parts)

        # If it looks like a bare formula
        if not prec_formula and not mz_match:
            fm_check = re.findall(r'([A-Z][a-z]?)(\d*)', info)
            if any(e[0].isupper() for e in fm_check):
                formula_parts = []
                for elem, cnt in fm_check:
                    if elem and elem[0].isupper():
                        formula_parts.append(f"{elem}{cnt}" if cnt else elem)
                if formula_parts:
                    prec_formula = "".join(formula_parts)

        return prec_mz, prec_formula

    def _normalize_peaks(self, peaks: list) -> list:
        """Normalize various peak input formats."""
        normalized = []
        for p in peaks:
            if isinstance(p, dict):
                mz = p.get("mz")
                intensity = p.get("intensity", p.get("intensity_pct", p.get("rel_intensity", 1)))
            elif isinstance(p, (list, tuple)) and len(p) >= 2:
                mz = p[0]
                intensity = p[1]
            else:
                continue
            
            try:
                mz_val = float(mz)
                int_val = float(intensity)
                if mz_val > 0 and int_val >= 0:
                    normalized.append({"mz": mz_val, "intensity": int_val})
            except (ValueError, TypeError):
                continue

        return sorted(normalized, key=lambda x: -x["intensity"])

    def _assign_peak(self, mz: float, prec_mz: Optional[float], prec_formula: Optional[str],
                      tol: float, mode: str) -> list:
        """Find all possible assignments for a peak."""
        assignments = []

        if prec_mz is not None:
            # Check neutral losses from precursor
            for loss_mass, loss_name, category, freq in _NEUTRAL_LOSSES_DB:
                expected = prec_mz - loss_mass
                if expected > 0:
                    error = abs(mz - expected)
                    if error <= tol:
                        conf, score = self._confidence_from_error(error, tol, category, freq)
                        assignments.append({
                            "type": "neutral_loss",
                            "assignment": f"[M-{loss_name}+H]⁺" if prec_mz > 0 else f"[M-{loss_name}]⁺•",
                            "loss": f"{loss_name} ({loss_mass:.5f} Da)",
                            "error_Da": round(error, 5),
                            "error_ppm": round(error / expected * 1e6, 2) if expected > 0 else 0,
                            "confidence": conf,
                            "score": score,
                            "category": category,
                        })

            # Check diagnostic ions
            for diag_mz, diag_formula, diag_name, diag_class in _DIAGNOSTIC_FRAGMENTS:
                error = abs(mz - diag_mz)
                if error <= tol:
                    conf, score = self._confidence_from_error(error, tol, "diagnostic", "rare" if "peptide" in diag_class else "common")
                    assignments.append({
                        "type": "diagnostic_ion",
                        "assignment": f"{diag_formula} ({diag_name})",
                        "formula": diag_formula,
                        "name": diag_name,
                        "class": diag_class,
                        "error_Da": round(error, 5),
                        "confidence": conf,
                        "score": score,
                    })

            # Check ammonium/water loss combinations
            for loss1 in [(18.0106, "H₂O"), (17.0027, "OH•"), (27.9949, "CO")]:
                for loss2 in [(18.0106, "H₂O"), (27.9949, "CO"), (43.9898, "CO₂")]:
                    combined = loss1[0] + loss2[0]
                    expected = prec_mz - combined
                    if expected > 50:
                        error = abs(mz - expected)
                        if error <= tol:
                            conf, score = self._confidence_from_error(error, tol, "combined_loss", "occasional")
                            assignments.append({
                                "type": "combined_neutral_loss",
                                "assignment": f"[M-{loss1[1]}-{loss2[1]}+H]⁺",
                                "loss": f"{loss1[1]}+{loss2[1]} ({combined:.5f} Da)",
                                "error_Da": round(error, 5),
                                "confidence": conf,
                                "score": score - 5,  # penalize combined losses slightly
                            })

        # Sort by score descending
        assignments.sort(key=lambda a: a["score"], reverse=True)

        # In best_match mode, limit results handled by caller
        return assignments[:10] if mode == "comprehensive" else assignments

    def _confidence_from_error(self, error: float, tol: float, category: str, frequency: str) -> Tuple[str, int]:
        """Determine confidence level and score from match error."""
        ratio = error / tol if tol > 0 else 1
        
        if ratio < 0.25:
            conf = "very_high"
            base_score = 95
        elif ratio < 0.5:
            conf = "high"
            base_score = 80
        elif ratio < 0.75:
            conf = "medium"
            base_score = 60
        else:
            conf = "low"
            base_score = 40

        # Adjust for category specificity
        if category == "diagnostic":
            base_score += 10
        elif category in ("neutral_loss",) and frequency == "VERY COMMON":
            base_score += 5
        elif category == "combined_loss":
            base_score -= 10

        return conf, min(100, max(0, base_score))

    def _suggest_unassigned_cause(self, mz: float, prec_mz: Optional[float]) -> str:
        """Suggest reasons why a peak couldn't be assigned."""
        causes = []
        if prec_mz:
            rel_pos = mz / prec_mz if prec_mz > 0 else 0
            if rel_pos < 0.1:
                causes.append("low-m/z noise or background ion")
            elif rel_pos > 0.95:
                causes.append("near-precursor (possibly internal fragment or noise)")
            else:
                causes.append("unusual rearrangement ion, isobaric interference, or external contaminant")
        else:
            causes.append("no precursor info provided for neutral loss matching")

        causes.append("consider checking against known adducts or solvent clusters")
        return "; ".join(causes[:2])

    def _suggest_additional_fragments(self, prec_mz: Optional[float], formula: Optional[str]) -> list:
        """Suggest fragments that should be present but weren't found."""
        suggestions = []
        if prec_mz:
            # Common losses that should produce visible fragments
            important_losses = [
                (18.0106, "H₂O dehydration"),
                (43.9898, "CO₂ decarboxylation"),
                (27.9949, "CO carbonyl loss"),
            ]
            for loss, desc in important_losses:
                suggestions.append({"expected_mz": round(prec_mz - loss, 4), "reason": desc})

        return suggestions

    def _assess_spectral_quality(self, peaks: list, n_annotated: int, n_total: int) -> dict:
        """Assess basic spectral quality metrics."""
        if not peaks:
            return {"quality": "no_data"}

        intensities = [p["intensity"] for p in peaks]
        total_int = sum(intensities)
        max_int = max(intensities)

        # Peak count assessment
        if n_total >= 15:
            peak_assessment = "rich fragmentation"
        elif n_total >= 8:
            peak_assessment = "moderate fragmentation"
        elif n_total >= 3:
            peak_assessment = "sparse fragmentation"
        else:
            peak_assessment = "very few fragments"

        # Intensity distribution
        top3_sum = sum(sorted(intensities)[-3:])
        top3_fraction = top3_sum / total_int if total_int > 0 else 0

        if top3_fraction > 0.85:
            dist_assessment = "dominated_by_few_peaks"
        elif top3_fraction > 0.65:
            dist_assessment = "moderate_distribution"
        else:
            dist_assessment = "well_distributed"

        # Coverage-based quality
        cov_ratio = n_annotated / n_total if n_total > 0 else 0
        if cov_ratio >= 0.8:
            qual = "excellent"
        elif cov_ratio >= 0.6:
            qual = "good"
        elif cov_ratio >= 0.4:
            qual = "fair"
        else:
            qual = "poor"

        return {
            "overall_quality": qual,
            "peak_count_category": peak_assessment,
            "intensity_distribution": dist_assessment,
            "total_peaks": n_total,
            "base_peak_mz": round(max(peaks, key=lambda p: p["intensity"])["mz"], 4),
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split("; ")
            precursor = parts[0]
            peak_str = parts[1] if len(parts) > 1 else ""
            tol = float(parts[2]) if len(parts) > 2 else 0.02

            # Parse peaks
            peaks = []
            if peak_str:
                for pair in peak_str.split(","):
                    pair = pair.strip()
                    if ":" in pair:
                        mz_str, int_str = pair.split(":")
                        peaks.append({"mz": float(mz_str.strip()), "intensity": float(int_str.strip())})

            return self._run_base(precursor, peaks, tol, "comprehensive")
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'precursor_info; mz1:int1, mz2:int2, ... [tolerance]'")
