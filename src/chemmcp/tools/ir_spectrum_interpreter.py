"""
IR Spectrum Interpreter — 红外光谱峰归属与官能团识别
输入峰位列表或光谱描述，识别官能团并归属特征吸收
"""
import logging
from typing import Optional, List, Dict, Any, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── IR 特征吸收数据库 (cm⁻¹) ───────────────────────────────────────
# 每条: (波数范围下限, 上限), {官能团信息}
IR_PEAK_DB: List[dict] = [
    # ── O-H 伸缩振动 (宽而强) ──
    {"range_cm-1": (3600, 3200), "group": "O-H stretch (alcohol/phenol)",
     "shape": "broad, strong", "intensity": "s", "notes": "Free OH ~3600; H-bonded broadens to 3300. Phenol slightly lower."},
    {"range_cm-1": (3400, 2400), "group": "O-H stretch (carboxylic acid)",
     "shape": "very broad (200-1000cm⁻¹ wide)", "intensity": "vs", "notes": "Characteristic very broad feature; often obscures C-H region. Carboxylic acid diagnostic."},

    # ── N-H 伸缩振动 ──
    {"range_cm-1": (3500, 3300), "group": "N-H stretch (primary amine)",
     "shape": "doublet (asym + sym)", "intensity": "m-s", "notes": "Primary: two bands (~3350, ~3460). Secondary: one band (~3300)."},
    {"range_cm-1": (3350, 3310), "group": "N-H stretch (secondary amine/amide)",
     "shape": "single sharp band", "intensity": "m", "notes": "Secondary amine or amide N-H stretch. Amide A band for proteins."},
    {"range_cm-1": (3300, 3030), "group": "N-H stretch (amide/protein)",
     "shape": "medium-strong", "intensity": "m", "notes": "Amide A band; protein secondary structure analysis region."},

    # ─═ C-H 伸缩振动 ═─
    {"range_cm-1": (3100, 3000), "group": "=C-H stretch (aromatic/alkene)",
     "shape": "sharp-medium", "intensity": "m", "notes": "Aromatic C-H just above 3000; alkene =C-H 3080-3020. Distinguishes sp² from sp³ C-H."},
    {"range_cm-1": (3000, 2840), "group": "-C-H stretch (alkane)",
     "shape": "sharp", "intensity": "m-s", "notes": "sp³ C-H below 3000. Symmetric/asymmetric CH₂ at ~2926/2853. CH₃ at ~2962/2872."},
    {"range_cm-1": (2830, 2695), "group": "C-H stretch (aldehyde)",
     "shape": "doublet (Fermi resonance)", "intensity": "w-m", "notes": "Aldehyde C-H Fermi doublet at ~2820 & 2720 cm⁻¹. Diagnostic for aldehydes."},
    {"range_cm-1": (2800, 2700), "group": "O=C-H stretch (formic acid derivative)",
     "shape": "weak", "intensity": "w", "notes": "Formate/formic acid characteristic."},

    # ─═ 三键区域 ═─
    {"range_cm-1": (2260, 2100), "group": "C≡C stretch (alkyne)",
     "shape": "sharp, medium", "intensity": "variable", "notes": "Terminal alkyne ~2140-2100 (stronger). Internal alkyne ~2260-2190 (weak/absent if symmetric)."},
    {"range_cm-1": (2260, 2220), "group": "C≡N stretch (nitrile)",
     "shape": "sharp, medium-strong", "intensity": "m", "notes": "Sharp medium-intensity band. Saturated nitrile ~2250; conjugated/aromatic nitrile ~2235-2220."},
    {"range_cm-1": (2150, 2110), "group": "Cumulene / isocyanide / diazo",
     "shape": "sharp", "intensity": "variable", "notes": "Isocyanides ~2120-2180 (very strong). Diazo compounds ~2100-2270."},

    # ─═ 羰基 C=O 区域 (最重要!) ═─
    {"range_cm-1": (1810, 1780), "group": "C=O stretch (acid chloride / anhydride asym)",
     "shape": "sharp, strong", "intensity": "vs", "notes": "Acid chloride highest frequency (~1800). Anhydride asymmetric stretch here."},
    {"range_cm-1": (1760, 1720), "group": "C=O stretch (ester / acid anhydride sym / γ-lactone)",
     "shape": "sharp, strong", "intensity": "vs", "notes": "Ester ~1740 (saturated) / 1718 (conjugated). Anhydride sym stretch ~1750-1730. γ-Lactone ~1760."},
    {"range_cm-1": (1740, 1720), "group": "C=O stretch (aldehyde / formate ester)",
     "shape": "sharp, strong", "intensity": "vs", "notes": "Aldehyde ~1735-1720 (saturated); lower if conjugated (~1705)."},
    {"range_cm-1": (1730, 1700), "group": "C=O stretch (carboxylic acid / α,β-unsaturated ketone)",
     "shape": "broad-ish, strong", "intensity": "s-vs", "notes": "Carboxylic acid ~1720-1680 (H-bonding lowers and broadens). α,β-Unsat ketone ~1680."},
    {"range_cm-1": (1715, 1690), "group": "C=O stretch (ketone)",
     "shape": "sharp, strong", "intensity": "vs", "notes": "Saturated ketone ~1715. Lowered by ring strain (cyclobutanone ~1780) or conjugation."},
    {"range_cm-1": (1690, 1640), "group": "C=O stretch (amide I / carboxylate / conjugated carbonyl)",
     "shape": "strong", "intensity": "s-vs", "notes": "Amide I (protein backbone) ~1650 (β-sheet) / 1655 (α-helix). Conjugation lowers to 1680-1660. Urea ~1645."},
    {"range_cm-1": (1670, 1640), "group": "C=O stretch (urea / strongly conjugated)",
     "shape": "strong", "intensity": "s", "notes": "Urea carbonyl. Strongly conjugated systems push lowest."},

    # ─═ C=C / C=N 双键区域 ═─
    {"range_cm-1": (1680, 1630), "group": "C=C stretch (alkene)",
     "shape": "variable intensity", "intensity": "variable", "notes": "Variable intensity (often weak unless conjugated). Terminal alkene stronger. cis weaker than trans."},
    {"range_cm-1": (1620, 1590), "group": "C=C stretch (aromatic ring quadrant stretch)",
     "shape": "medium-sharp", "intensity": "m", "notes": "Aromatic quadrant stretch. Often multiple bands in 1600-1450 fingerprint region."},
    {"range_cm-1": (1610, 1580), "group": "C=C aromatic + conjugation",
     "shape": "sharp", "intensity": "m", "notes": "Aromatic semicircle stretch. Substitution pattern affects number/intensity of bands."},
    {"range_cm-1": (1650, 1580), "group": "N-H bend (primary amine) + C=N stretch (imine/quinoxaline)",
     "shape": "medium", "intensity": "m", "notes": "Primary amine N-H bend (scissoring) ~1640. Imine C=N ~1660-1640."},

    # ─═ 指纹区重要峰 ═─
    {"range_cm-1": (1560, 1510), "group": "Aromatic semicircle stretch / N-O asymmetric (nitro)",
     "shape": "sharp", "intensity": "m-s", "notes": "Nitro compounds: NO₂ asym str ~1550-1510 (strong). Aromatic ~1500."},
    {"range_cm-1": (1485, 1445), "group": "CH₂ scissoring / CH₃ bend (umbrella) / aromatic",
     "shape": "medium", "intensity": "m", "notes": "CH₂ scissor ~1465. CH₃ umbrella bend ~1450. Aromatic ~1450."},
    {"range_cm-1": (1440, 1390), "group": "CH₃ bend (asymmetric) / C-H rock / gem-dimethyl",
     "shape": "medium", "intensity": "m", "notes": "gem-Dimethyl splitting ~1385/1370. Characteristic for t-butyl/isopropyl groups."},
    {"range_cm-1": (1380, 1350), "group": "CH₃ symmetric bend (umbrella) / NO₂ sym stretch / S=O sym",
     "shape": "medium", "intensity": "m", "notes": "NO₂ symmetric ~1350 (strong). Sulfoxide S=O sym ~1350-1310."},
    {"range_cm-1": (1330, 1290), "group": "S=O asymmetric stretch (sulfone) / nitro aromatics",
     "shape": "strong", "intensity": "s", "notes": "Sulfone asym S=O₂ ~1325-1300 (very strong)."},
    {"range_cm-1": (1300, 1150), "group": "C-O stretch (alcohol/ester/ether) / amide III",
     "shape": "strong", "intensity": "s-vs", "notes": "Alcohol C-O ~1200-1050 (broader). Ester C-O(-C) ~1250-1150. Ether ~1140-1020. Amide III ~1300-1200."},
    {"range_cm-1": (1120, 1030), "group": "C-O stretch (alcohol/primary-secondary) / Si-O-Si",
     "shape": "strong", "intensity": "s", "notes": "Primary alcohol ~1050. Secondary ~1100. Tertiary ~1150. Siloxane ~1080-1010 (very strong, broad)."},
    {"range_cm-1": (1000, 950), "group": "C=C bending (alkene) / =C-H out-of-plane (vinyl)",
     "shape": "variable", "intensity": "variable", "notes": "Vinyl =C-H oop ~995/910 (terminal alkene diagnostic). Ring breathing mode for mono-substituted benzene ~1000-990."},
    {"range_cm-1": (900, 650), "group": "=C-H out-of-plane (aromatic substitution pattern)",
     "shape": "sharp, strong", "intensity": "s-vs", "notes": "CRITICAL for determining aromatic substitution pattern! Mono-subst: ~750/690. ortho: ~750. meta: ~880/780. para: ~850."},
    {"range_cm-1": (770, 630), "group": "C-Cl stretch / C-Br stretch (halogen)",
     "shape": "strong", "intensity": "s", "notes": "C-Cl ~770-700 (strong). C-Br ~650-500. Multiple halogens give multiple bands."},
    {"range_cm-1": (600, 450), "group": "Heavy atom stretches / skeletal bends / metal-oxygen",
     "shape": "variable", "intensity": "variable", "notes": "Low-frequency region. Metal-oxygen bonds appear here. Useful for inorganic complexes."},
]


@ChemMCPManager.register_tool
class IrSpectrumInterpreter(BaseTool):
    """
    红外光谱解析工具：根据输入的峰位列表（波数 cm⁻¹），
    自动匹配官能团，给出归属建议和置信度评分。
    """
    __version__ = "0.1.0"
    name = "IrSpectrumInterpreter"
    func_name = "interpret_ir_spectrum"
    description = "Interpret IR spectrum peaks by identifying functional groups from wavenumber (cm⁻¹) data. Supports peak list input and provides assignments with confidence scores."
    implementation_description = "Matches input peak positions against a comprehensive database of IR absorption ranges covering O-H, N-H, C-H, triple bond, carbonyl, double bond, fingerprint, and out-of-plane bending regions. Ranks matches by proximity to expected range centers."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["IR Spectroscopy", "FTIR", "Functional Group Identification", "Structural Elucidation", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("peaks", "list", "required", "List of observed peak positions in cm⁻¹ (e.g., [3400, 2920, 2850, 1710, 1600, 1500, 1450])."),
        ("intensities", "list", "[]", "Optional intensity labels matching each peak: 'vs', 's', 'm', 'w', 'br' (broad)."),
        ("sample_type", "str", "organic", "Sample type hint: 'organic', 'polymer', 'protein', 'inorganic_complex', 'unknown'."),
        ("confidence_threshold", "float", "0.3", "Minimum match confidence (0-1) to include in results."),
        ("mode", "str", "identify", "Mode: 'identify' (peak→group), 'reverse_lookup' (group→expected peaks), 'full_report'."),
        ("target_groups", "list", "[]", "For reverse_lookup mode: list of functional group names to find expected peaks."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., '3400 2920 2850 1710 1600 1500 1450 1370 1250 1050 850 720' or 'identify protein'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with identified functional groups, peak assignments, confidence scores, and structural interpretation summary."),
    ]

    examples = [
        {
            "code_input": {
                "peaks": [3400, 3050, 2920, 2850, 1710, 1600, 1500, 1450, 1370,
                          1260, 1170, 1100, 830, 720],
                "intensities": [],
                "sample_type": "organic",
                "confidence_threshold": 0.3,
                "mode": "identify",
                "target_groups": [],
            },
            "text_input": {
                "input_params": "3400 3050 2920 2850 1710 1600 1500 1450 1370 1260 1170 1100 830 720",
            },
            "output": {
                "result": {
                    "mode": "ir_interpretation",
                    "n_peaks_input": 14,
                    "note": "Example output showing functional group identification.",
                }
            }
        },
        {
            "code_input": {
                "peaks": [3290, 3060, 2920, 2850, 1650, 1540, 1450, 1240, 1080, 720],
                "intensities": ["br", "", "", "", "s", "s", "", "s", "s", ""],
                "sample_type": "protein",
                "confidence_threshold": 0.25,
                "mode": "identify",
                "target_groups": [],
            },
            "text_input": {
                "input_params": "3290 br 3060 2920 2850 1650 s 1540 s 1450 1240 s 1080 s 720 protein",
            },
            "output": {
                "result": {
                    "mode": "ir_interpretation",
                    "sample_type": "protein",
                    "note": "Protein/amide-focused interpretation.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _match_peak(peak: float, db_entry: dict) -> float:
        """Calculate match confidence (0-1) based on how close peak is to range center."""
        lo, hi = db_entry["range_cm-1"]
        center = (lo + hi) / 2
        half_width = (hi - lo) / 2

        if lo <= peak <= hi:
            # Within range: score based on distance from center
            dist = abs(peak - center)
            return max(0.3, 1.0 - (dist / half_width) * 0.5)
        elif peak < lo:
            # Below range: penalize by distance
            gap = lo - peak
            return max(0, 0.5 - gap / 200)
        else:
            # Above range
            gap = peak - hi
            return max(0, 0.5 - gap / 200)

    def _identify_peaks(self, peaks: List[float], intensities: List[str] = None,
                        threshold: float = 0.3) -> List[dict]:
        """Match each peak against the database."""
        if intensities is None:
            intensities = []

        all_matches = []
        for i, peak in enumerate(peaks):
            peak_int = intensities[i].lower() if i < len(intensities) else ""
            best_match = None
            best_score = 0

            for entry in IR_PEAK_DB:
                score = self._match_peak(peak, entry)
                if score >= threshold and score > best_score:
                    best_score = score
                    best_match = {**entry, "match_score": round(score, 3),
                                 "observed_peak": peak}

            if best_match:
                all_matches.append(best_match)

        # Sort by match score descending
        all_matches.sort(key=lambda x: x["match_score"], reverse=True)
        return all_matches

    def _summarize_structure(self, matches: List[dict], sample_type: str) -> dict:
        """Generate a structural interpretation summary."""
        found_groups = set()
        has_oh = False
        has_nh = False
        has_carbonyl = False
        has_aromatic = False
        has_alkane = False
        has_alkene = False
        has_triple_bond = False
        has_nitro = False
        has_c_o = False
        has_halogen = False

        for m in matches:
            g = m["group"].lower()
            found_groups.add(m["group"])
            if "o-h" in g or "alcohol" in g or "carboxylic" in g or "phenol" in g:
                has_oh = True
            if "n-h" in g or "amine" in g or "amide" in g:
                has_nh = True
            if "c=o" in g or "carbonyl" in g or "carboxyl" in g:
                has_carbonyl = True
            if "aromatic" in g:
                has_aromatic = True
            if "alkane" in g or "-c-h" in g and "alkene" not in g:
                has_alkane = True
            if "alkene" in g or "c=c" in g and "aromatic" not in g:
                has_alkene = True
            if "c≡" in g or "triple" in g or "nitrile" in g:
                has_triple_bond = True
            if "nitro" in g:
                has_nitro = True
            if "c-o" in g and "carbonyl" not in g:
                has_c_o = True
            if "cl" in g or "br" in g or "halogen" in g:
                has_halogen = True

        # Build structural description
        features = []
        if sample_type == "protein":
            features.append("Protein/amide profile detected")
            features.append("Key bands: Amide I (~1650), Amide II (~1540), Amide III (~1240)")
        else:
            parts = []
            if has_oh:
                parts.append("hydroxyl/alcohol/carboxylic acid")
            if has_nh:
                parts.append("amine/amide")
            if has_carbonyl:
                parts.append("carbonyl compound")
            if has_aromatic:
                parts.append("aromatic ring(s)")
            if has_alkene:
                parts.append("alkene")
            if has_triple_bond:
                parts.append("triple bond (alkyne/nitrile)")
            if has_nitro:
                parts.append("nitro group")
            if has_c_o:
                parts.append("ether/ester C-O")
            if has_halogen:
                parts.append("halogen substituent")

            if parts:
                features.append(f"Spectrum suggests: {', '.join(parts)}")
            else:
                features.append("Limited characteristic features detected. Consider checking fingerprint region more carefully.")

        return {
            "functional_groups_found": sorted(found_groups),
            "key_features": {
                "has_OH": has_oh, "has_NH": has_nh, "has_carbonyl": has_carbonyl,
                "has_aromatic": has_aromatic, "has_alkene": has_alkene,
                "has_triple_bond": has_triple_bond, "has_nitro": has_nitro,
                "has_C_O_single": has_c_o, "has_halogen": has_halogen,
            },
            "structural_summary": features,
            "next_steps": [
                "Confirm assignments by comparing with reference spectrum (if available).",
                "Check for combination/overtone bands in 2000-1650 cm⁻¹ region.",
                "Consider running complementary analysis (NMR, MS) for full structural confirmation.",
                "Note that some peaks may be overtones or combinations rather than fundamental vibrations.",
            ],
        }

    def _run_base(self, peaks: list = None, intensities: list = None,
                  sample_type: str = "organic", confidence_threshold: float = 0.3,
                  mode: str = "identify", target_groups: list = None) -> dict:

        if target_groups is None:
            target_groups = []
        if peaks is None:
            raise ChemMCPError("'peaks' list is required for identify mode.")

        if len(peaks) == 0:
            raise ChemMCPError("'peaks' list cannot be empty.")

        if mode == "reverse_lookup":
            return self._reverse_lookup(target_groups)

        matches = self._identify_peaks(peaks, intensities, confidence_threshold)
        summary = self._summarize_structure(matches, sample_type)

        # Peak-by-peak assignment table
        assignment_table = []
        for i, peak in enumerate(peaks):
            peak_int = intensities[i] if intensities and i < len(intensities) else ""
            assigned = [m for m in matches if m.get("observed_peak") == peak]
            assignment_table.append({
                "peak_cm-1": peak,
                "intensity": peak_int or "not specified",
                "assignments": [{"group": a["group"], "score": a["match_score"],
                                  "shape_hint": a.get("shape", ""),
                                  "notes": a.get("notes", "")} for a in assigned[:3]],
            })

        return {"result": {
            "mode": "ir_interpretation",
            "input": {
                "n_peaks": len(peaks),
                "peaks_cm-1": peaks,
                "intensities": intensities or "none provided",
                "sample_type": sample_type,
                "threshold": confidence_threshold,
            },
            "assignments": assignment_table,
            "all_matches_sorted_by_confidence": matches,
            "structural_summary": summary,
            "spectral_regions_coverage": self._check_region_coverage(peaks),
            "database_info": {
                "total_entries_in_db": len(IR_PEAK_DB),
                "wavenumber_range_cm-1": [450, 3700],
            },
        }}

    def _reverse_lookup(self, target_groups: list) -> dict:
        """Given functional group names, find expected IR peaks."""
        results = []
        for tgt in target_groups:
            tgt_lower = tgt.lower()
            matches = [e for e in IR_PEAK_DB if any(kw in e["group"].lower() for kw in tgt_lower.split())]
            results.append({"target": tgt, "expected_peaks": matches})

        return {"result": {
            "mode": "reverse_lookup",
            "targets_searched": target_groups,
            "results": results,
        }}

    @staticmethod
    def _check_region_coverage(peaks: list) -> dict:
        """Check which spectral regions have peaks."""
        regions = {
            "O-H/N-H stretch (3700-3100)": any(p > 3100 for p in peaks),
            "C-H stretch (3100-2800)": any(2800 <= p <= 3100 for p in peaks),
            "Triple bond/C≡N (2300-2050)": any(2050 <= p <= 2300 for p in peaks),
            "Carbonyl C=O (1820-1640)": any(1640 <= p <= 1820 for p in peaks),
            "C=C/C=N (1680-1570)": any(1570 <= p <= 1680 for p in peaks),
            "Fingerprint (1500-400)": any(p < 1500 for p in peaks),
        }
        covered = sum(1 for v in regions.values() if v)
        return {
            "regions": regions,
            "regions_with_peaks": f"{covered}/{len(regions)}",
            "coverage_pct": round(covered / len(regions) * 100, 0),
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()

            # Check for special modes
            if parts[0].lower() == "protein":
                return self._run_base(
                    peaks=[3290, 3060, 2920, 2850, 1650, 1540, 1450, 1240, 1080, 720],
                    intensities=["br", "", "", "", "s", "s", "", "s", "s", ""],
                    sample_type="protein", confidence_threshold=0.25)

            # Parse peak list (numbers, optionally followed by intensity labels)
            peaks = []
            intensities = []
            for p in parts:
                p_lower = p.lower().strip()
                if p_lower in ("vs", "s", "m", "w", "br", "vb"):
                    intensities.append(p_lower)
                else:
                    try:
                        peaks.append(float(p))
                    except ValueError:
                        # Could be sample type at end
                        continue

            stype = "organic"
            for p in parts[-3:]:
                if p.lower() in ("organic", "polymer", "protein", "inorganic_complex", "unknown"):
                    stype = p.lower()

            return self._run_base(peaks=peaks, intensities=intensities,
                                   sample_type=stype)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
