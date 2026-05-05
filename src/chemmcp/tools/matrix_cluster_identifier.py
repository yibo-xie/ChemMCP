"""
Matrix Cluster Identifier - Identifies matrix-related cluster ions in MALDI and ESI spectra
including common matrix peaks, solvent clusters, and background ions.
"""

import logging
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# MALDI matrix database: matrix_name → {clusters, adducts, fragments, notes}
_MALDI_MATRICES: Dict[str, Dict] = {
    "CHCA": {
        "full_name": "α-Cyano-4-hydroxycinnamic acid",
        "molecular_weight": 189.0426,
        "common_clusters": [
            # (mz, formula, description)
            (190.0504,  "[M+H]⁺",    "Matrix protonated molecule"),
            (212.0323,  "[M+Na]⁺",   "Matrix sodium adduct"),
            (228.0267,  "[M+K]⁺",    "Matrix potassium adduct"),
            (172.0400,  "[M-H₂O+H]⁺","Matrix dehydration"),
            (144.0451,  "[M-H₂O-CO+H]⁺", "Matrix -H₂O-CO"),
            (117.0345,  "C₈H₅O⁺",   "Characteristic fragment"),
            (89.0244,   "C₆H₅O⁺",   "Fragment"),
            (379.0856,  "[2M+H]⁺",   "Matrix dimer"),
            (401.0675,  "[2M+Na]⁺",  "Dimer sodium adduct"),
            (188.0422,  "[M]⁺•",     "Matrix radical"),
            (146.0294,  "[M-HCOOH+H]⁺", "Matrix - formic acid"),
            (164.0355,  "[M-H₂O-H+]? or [M-CHCN+H]?", "Minor fragment"),
        ],
        "background_range": (50, 500),
        "notes": "Most common MALDI matrix for peptides/proteins; strong UV absorption at 355 nm",
    },
    "DHB": {
        "full_name": "2,5-Dihydroxybenzoic acid",
        "molecular_weight": 154.0266,
        "common_clusters": [
            (155.0339,  "[M+H]⁺",    "Matrix protonated molecule"),
            (177.0158,  "[M+Na]⁺",   "Sodium adduct"),
            (193.0098,  "[M+K]⁺",    "Potassium adduct"),
            (137.0234,  "[M-H₂O+H]⁺","Dehydration product"),
            (109.0285,  "[M-H₂O-CO+H]⁺", "-H₂O-CO loss"),
            (137.0234,  "[M-H₂O+H]⁺", "Dehydrated"),
            (91.0184,   "C₇H₇O⁺? or fragment", "Aromatic fragment"),
            (309.0541,  "[2M+H]⁺",   "Dimer"),
            (331.0360,  "[2M+Na]⁺",  "Dimer + Na"),
            (123.0179,  "[M-CH₂O₂+H]⁺", "-32 Da (CH₂O₂) loss"),
            (97.0239,   "C₅H₅O₂⁺ or C₆H₅O⁺?", "Small fragment"),
        ],
        "background_range": (50, 600),
        "notes": "Good for carbohydrates and glycans; forms fine crystals",
    },
    "SA": {
        "full_name": "Sinapinic acid",
        "molecular_weight": 224.0683,
        "common_clusters": [
            (225.0756,  "[M+H]⁺",    "Protonated matrix"),
            (247.0575,  "[M+Na]⁺",   "Sodium adduct"),
            (263.0515,  "[M+K]⁺",    "Potassium adduct"),
            (207.0650,  "[M-H₂O+H]⁺","Dehydration"),
            (179.0701,  "[M-H₂O-CO+H]⁺", "-H₂O-CO"),
            (449.1435,  "[2M+H]⁺",   "Dimer"),
            (214.0549,  "[M-CH₃+H]⁺? or demethylated", "Demethylated"),
            (189.0550,  "fragment?",  "Possible fragment"),
        ],
        "background_range": (100, 700),
        "notes": "Preferred for intact proteins; good for high MW compounds",
    },
    "DCTB": {
        "full_name": "trans-2-[3-(4-tert-butylphenyl)-2-methyl-2-propenylidene]malononitrile",
        "molecular_weight": 278.1464,
        "common_clusters": [
            (279.1537,  "[M+H]⁺",    "Protonated"),
            (301.1356,  "[M+Na]⁺",   "Sodium adduct"),
            (557.2981,  "[2M+H]⁺",   "Dimer"),
            (264.1272,  "[M-CH₃+H]⁺? or fragment", "Fragment"),
        ],
        "background_range": (100, 800),
        "notes": "Specialized for polymers and synthetic materials analysis",
    },
    "9AA": {
        "full_name": "9-Aminoacridine",
        "molecular_weight": 194.0920,
        "common_clusters": [
            (195.0993,  "[M+H]⁺",    "Protonated"),
            (217.0812,  "[M+Na]⁺",   "Sodium adduct"),
            (176.0787,  "[M-NH₃+H]⁺? or fragment", "Ammonia loss"),
            (389.1910,  "[2M+H]⁺",   "Dimer"),
        ],
        "background_range": (80, 500),
        "notes": "Commonly used for negative mode metabolite imaging",
    },
}

# ESI background / solvent cluster ions
_ESI_BACKGROUND_IONS: Dict[str, List[Dict]] = {
    "positive": [
        {"mz": 84.0419,  "formula": "C₅H₆N⁺",      "name": "Acetonitrile-derived cluster"},
        {"mz": 102.0550,  "formula": "C₅H₈NO₂⁺",    "name": "ACN/water cluster"},
        {"mz": 112.0867,  "formula": "C₇H₁₁N⁺",     "name": "Unknown ESI background"},
        {"mz": 130.1000,  "formula": "C₈H₁₂NO⁺",    "name": "Solvent cluster"},
        {"mz": 132.0809,  "formula": "C₆H₁₁NO₂⁺",   "name": "Plasticizer/phthalate background"},
        {"mz": 149.0237,  "formula": "C₈H₅O₃⁺",     "name": "Phthalate fragment (common contaminant)"},
        {"mz": 183.0552,  "formula": "C₉H₉O₄⁺",     "name": "PEG-related background"},
        {"mz": 195.0877,  "formula": "C₁₀H₁₂NO₃⁺",  "name": "Triton/PEG background"},
        {"mz": 207.0827,  "formula": "C₁₀H₁₂NO₄⁺",  "name": "Background ion"},
        {"mz": 215.1633,  "formula": "C₁₂H₂₀O₃⁺",   "name": "PEG fragment (m=4)"},
        {"mz": 223.1489,  "formula": "C₁₂H₁₉NO₃⁺",  "name": "Triton-related"},
        {"mz": 255.1913,  "formula": "C₁₄H₂₆O₃⁺",   "name": "PEG fragment (m=5)"},
        {"mz": 279.1591,  "formula": "C₁₅H₂₂NO₃⁺",  "name": "Triton/PEG background"},
        {"mz": 391.2843,  "formula": "C₂₀H₄₀O₆⁺",   "name": "PEG oligomer (n≈8)"},
        {"mz": 132.0222,  "formula": "C₆H₅NO₂⁺",    "name": "Benzamide contamination"},
    ],
    "negative": [
        {"mz": 91.0031,   "formula": "C₆H₅O₂⁻",     "name": "Benzoate / formate cluster"},
        {"mz": 96.9601,   "formula": "H₂PO₄⁻",       "name": "Phosphate (ubiquitous in negative mode)"},
        {"mz": 117.0185,  "formula": "C₅HO₄⁻",      "name": "Phthalate fragment"},
        {"mz": 119.0341,  "formula": "C₅H₄O₃⁻",     "name": "Terephthalate"},
        {"mz": 159.1441,  "formula": "C₈H₂₀O₃⁻",    "name": "PEG fragment"},
        {"mz": 175.1240,  "formula": "C₈H₁₈O₄⁻",    "name": "PEG fragment"},
        {"mz": 191.1189,  "formula": "C₈H₁₆O₅⁻",    "name": "PEG fragment"},
        {"mz": 207.1138,  "formula": "C₈H₁₄O₆⁻",    "name": "PEG fragment"},
        {"mz": 267.1820,  "formula": "C₁₄H₂₆O₄⁻",   "name": "PEG fragment"},
        {"mz": 283.1770,  "formula": "C₁₄H₂₄O₅⁻",   "name": "PEG fragment"},
        {"mz": 353.2220,  "formula": "C₁₈H₃₂O₅⁻",   "name": "PEG oligomer"},
    ],
}

# Common contaminant database
_CONTAMINANTS: List[Dict] = [
    {"mz": 149.0237, "formula": "C₈H₅O₃⁺", "name": "Phthalate (plasticizer) — VERY COMMON", "source": "lab plasticware, tubing"},
    {"mz": 195.0877, "formula": "C₁₀H₁₂NO₃⁺", "name": "Triton X-100 residue", "source": "lab detergent"},
    {"mz": 207.0827, "formula": "C₁₀H₁₂NO₄⁺", "name": "Triton-related", "source": "lab detergent"},
    {"mz": 221.0612, "formula": "C₉H₁₀O₅⁺", "name": "Unknown lab contaminant", "source": "unknown"},
    {"mz": 332.2392, "formula": "C₁₆H₃₃NO₄⁺", "name": "Triton/Polyethylene glycol", "source": "detergent"},
    {"mz": 391.2843, "formula": "C₂₀H₄₀O₆⁺", "name": "PEG (n=8)", "source": "cosmetics, solvents"},
    {"mz": 435.3160, "formula": "C₂₂H₄₅O₈⁺", "name": "PEG (n=9)", "source": "cosmetics, solvents"},
]


@ChemMCPManager.register_tool
class MatrixClusterIdentifier(BaseTool):
    """
    基质簇离子识别器 — 识别 MALDI 或 ESI 质谱中的基质相关簇离子和背景峰。
    
    包含常见 MALDI 基质（CHCA、DHB、SA 等）的特征峰、ESI 溶剂簇离子、
    以及常见实验室污染物的数据库，帮助区分分析物信号与背景干扰。
    """
    __version__      = "0.1.0"
    name             = "MatrixClusterIdentifier"
    func_name        = "identify_matrix_clusters"
    description      = "Identify matrix-related cluster ions, background peaks, and common contaminants in MALDI and ESI mass spectra."
    implementation_description = "Matches observed m/z values against comprehensive databases of MALDI matrix clusters (CHCA, DHB, SA, DCTB, 9AA), ESI solvent/additive clusters, and ubiquitous laboratory contaminants (phthalates, PEG, Triton). Provides assignment confidence and cleaning recommendations."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Mass Spectrometry", "Matrix Clusters", "MALDI", "Background Ions", "Contaminants"]
    required_envs    = []

    code_input_sig   = [
        ("matrix_name", "str", "DHB", "MALDI matrix name ('CHCA', 'DHB', 'SA', 'DCTB', '9AA') or 'ESI' for electrospray background."),
        ("ionization_mode", "str", "positive", "Ionization mode: 'positive' or 'negative'."),
        ("mz_range", "tuple", "(50, 2000)", "m/z range of interest as (min_mz, max_mz)."),
        ("observed_peaks", "list", "None", "List of observed m/z values to check against databases. If None, returns expected clusters for the given matrix/mode."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'matrix_name [ionization_mode] [mz_min,mz_max] [peak1,peak2,...]'. Example: 'DHB positive 50,800 155.03,177.02,379.05'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with expected_cluster_ions list, matched_assignments (if peaks provided), background_peak_list, contaminant_warnings, and cleaning_recommendations."),
    ]

    examples         = [
        {
            "code_input": {
                "matrix_name": "DHB",
                "ionization_mode": "positive",
                "mz_range": (50, 800),
                "observed_peaks": [155.034, 177.016, 309.054, 286.144],
            },
            "text_input": {
                "input_params": "DHB positive 50,800 155.03,177.02,309.05,286.14"
            },
            "output": {
                "result": {
                    "matrix_name": "DHB",
                    "matrix_full_name": "2,5-Dihydroxybenzoic acid",
                    "expected_clusters": 11,
                    "matched_assignments": [
                        {"observed_mz": 155.034, "assignment": "[M+H]⁺ (DHB)", "type": "matrix_ion", "error_ppm": 0.6},
                        {"observed_mz": 177.016, "assignment": "[M+Na]⁺ (DHB)", "type": "matrix_adduct", "error_ppm": 1.1},
                        {"observed_mz": 309.054, "assignment": "[2M+H]⁺ (DHB dimer)", "type": "matrix_dimer", "error_ppm": 0.0},
                    ],
                    "unassigned_peaks": [{"mz": 286.144, "note": "Not a known DHB cluster — possible analyte signal"}],
                    "contaminant_warnings": [],
                    "cleaning_recommendations": ["Peak at m/z 286.14 is likely an analyte signal (not a matrix cluster)"],
                }
            },
        },
        {
            "code_input": {
                "matrix_name": "ESI",
                "ionization_mode": "positive",
                "observed_peaks": [149.024, 286.144, 391.285],
            },
            "text_input": {
                "input_params": "ESI positive 0,2000 149.02,286.14,391.29"
            },
            "output": {
                "result": {
                    "matched_assignments": [
                        {"observed_mz": 149.024, "assignment": "Phthalate (plasticizer contaminant)", "type": "contaminant"},
                        {"observed_mz": 391.285, "assignment": "PEG oligomer (n≈8)", "type": "contaminant"},
                    ],
                    "cleaning_recommendations": ["Check plasticware for phthalate source", "Avoid PEG-containing consumables"],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, matrix_name: str = "DHB", ionization_mode: str = "positive",
                  mz_range: Tuple[int, int] = (50, 2000),
                  observed_peaks: Optional[List[float]] = None) -> dict:
        """Core logic."""
        matrix_upper = matrix_name.upper()
        
        # Determine if MALDI or ESI
        is_maldi = matrix_upper in _MALDI_MATRICES
        is_esi = matrix_upper == "ESI"

        if not is_maldi and not is_esi:
            available = list(_MALDI_MATRICES.keys()) + ["ESI"]
            raise ChemMCPError(f"Unknown matrix '{matrix_name}'. Available: {available}")

        mode = ionization_mode.lower()
        if mode not in ("positive", "negative"):
            raise ChemMCPError("Ionization mode must be 'positive' or 'negative'.")

        mz_lo, mz_hi = mz_range

        # Build expected cluster list
        expected_clusters = []
        matrix_info = None

        if is_maldi:
            matrix_info = _MALDI_MATRICES[matrix_upper]
            for entry in matrix_info["common_clusters"]:
                emz, formula, desc = entry
                if mz_lo <= emz <= mz_hi:
                    expected_clusters.append({
                        "mz": round(emz, 4),
                        "formula": formula,
                        "description": desc,
                        "type": "matrix_ion" if "dimer" not in desc.lower() else "matrix_dimer",
                        "source": f"{matrix_info['full_name']}",
                    })

        # Add ESI background ions
        esi_bg = _ESI_BACKGROUND_IONS.get(mode, [])
        for ion in esi_bg:
            emz = ion["mz"]
            if mz_lo <= emz <= mz_hi:
                expected_clusters.append({
                    "mz": round(emz, 4),
                    "formula": ion["formula"],
                    "description": ion["name"],
                    "type": "esi_background",
                    "source": "ESI solvent/suppressant",
                })

        # Add contaminants (both modes always checked)
        for cont in _CONTAMINANTS:
            cmz = cont["mz"]
            if mz_lo <= cmz <= mz_hi:
                expected_clusters.append({
                    "mz": round(cmz, 4),
                    "formula": cont["formula"],
                    "description": cont["name"],
                    "type": "contaminant",
                    "source": cont.get("source", "laboratory environment"),
                })

        # Match observed peaks if provided
        matched = []
        unassigned = []

        if observed_peaks:
            # Tolerance for matching (ppm-based with floor)
            tolerance_da = 0.02  # 20 mDa default

            for obs_mz in observed_peaks:
                best_match = None
                best_error = float("inf")

                for cluster in expected_clusters:
                    error = abs(obs_mz - cluster["mz"])
                    if error < tolerance_da and error < best_error:
                        best_error = error
                        best_match = dict(cluster)
                        best_match["observed_mz"] = round(obs_mz, 4)
                        best_match["error_Da"] = round(error, 5)
                        best_match["error_ppm"] = round(error / obs_mz * 1e6, 2) if obs_mz > 0 else 0

                if best_match:
                    matched.append(best_match)
                else:
                    unassigned.append({
                        "mz": round(obs_mz, 4),
                        "note": "Not matched to any known matrix/contaminant peak — possibly analyte signal",
                    })

        # Generate cleaning recommendations
        recommendations = self._generate_recommendations(matched, unassigned, matrix_name, mode)

        # Contaminant warnings
        contaminant_matches = [m for m in matched if m.get("type") == "contaminant"]

        return {
            "result": {
                "matrix_name": matrix_name,
                "matrix_full_name": matrix_info["full_name"] if matrix_info else "ESI Background",
                "matrix_mw": matrix_info["molecular_weight"] if matrix_info else None,
                "ionization_mode": mode,
                "mz_range": (mz_lo, mz_hi),
                "total_expected_clusters": len(expected_clusters),
                "expected_cluster_ions": sorted(expected_clusters, key=lambda x: x["mz"]),
                "matched_assignments": matched if observed_peaks else [],
                "unassigned_peaks": unassigned if observed_peaks else [],
                "match_summary": (
                    f"{len(matched)}/{len(observed_peaks)} peaks matched" if observed_peaks
                    else "No observed peaks provided — showing expected cluster database"
                ),
                "contaminant_warnings": [
                    f"⚠ {m['description']} at m/z {m['observed_mz']} — source: {m.get('source', '?')}"
                    for m in contaminant_matches
                ] if contaminant_matches else [],
                "cleaning_recommendations": recommendations,
                "notes": (
                    "Matrix Cluster Identification Notes:\n"
                    "• Always run blank samples to identify background peaks\n"
                    "• MALDI matrix signals should be present in all spots/spectra\n"
                    "• Analyte peaks are those NOT matching matrix/background\n"
                    "• Phthalates (m/z 149.02) are extremely common — use glassware when possible\n"
                    "• PEG series appears as repeating 44.03 Da spacings\n"
                    "• Consider recrystallizing MALDI matrix to reduce background"
                    if is_maldi else
                    "ESI Background Notes:\n"
                    "• Run solvent blanks regularly to monitor background\n"
                    "• Phthalates from plastics are the most common contaminant\n"
                    "• PEG from cosmetics/hand cream produces characteristic series\n"
                    "• Sodium/potassium adducts indicate salt contamination\n"
                    "• Use fresh solvents and avoid plasticware near the ion source"
                ),
            }
        }

    def _generate_recommendations(self, matched: list, unassigned: list, matrix: str, mode: str) -> List[str]:
        """Generate cleaning/recommendation suggestions."""
        recs = []

        if any(m.get("type") == "contaminant" for m in matched):
            sources = set(m.get("source", "?") for m in matched if m.get("type") == "contaminant")
            recs.append(f"Contaminants detected from: {', '.join(sources)}. Review sample handling and consumables.")

        phthalate_found = any("phthalate" in str(m).lower() or m.get("mz") == 149.0237 for m in matched)
        if phthalate_found:
            recs.append("Replace plastic vials/tubing with glass or PTFE alternatives to eliminate phthalate source.")

        peg_found = any("peg" in str(m).lower() for m in matched)
        if peg_found:
            recs.append("PEG contamination detected — avoid hand creams, detergents, and PEG-based consumables near MS lab.")

        dimer_found = any("dimer" in str(m.get("type", "")).lower() for m in matched)
        if dimer_found:
            if matrix.upper() in _MALDI_MATRICES:
                recs.append(f"Matrix dimer detected — reduce matrix concentration or improve crystallization for {matrix}.")

        if unassigned:
            analyte_candidates = [u for u in unassigned if "analyte" in u.get("note", "").lower()]
            if analyte_candidates:
                recs.append(f"{len(analyte_candidates)} potential analyte signal(s) detected (not matching matrix/background).")
            else:
                recs.append(f"{len(unassigned)} unassigned peak(s) may warrant further investigation.")

        if not recs:
            recs.append("Spectrum appears clean relative to known matrix/contaminant database.")

        return recs

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            matrix = parts[0]
            mode = parts[1] if len(parts) > 1 else "positive"
            
            # Parse mz range
            range_part = parts[2] if len(parts) > 2 else "50,2000"
            rlo, rhi = map(float, range_part.split(","))
            mz_range = (int(rlo), int(rhi))

            # Parse observed peaks
            peaks = None
            if len(parts) > 3:
                peaks = [float(p) for p in parts[3].split(",")]

            return self._run_base(matrix, mode, mz_range, peaks)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'matrix_name [mode] [mz_min,mz_max] [peak1,peak2,...]'")
