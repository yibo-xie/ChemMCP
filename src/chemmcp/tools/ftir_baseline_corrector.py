"""
FTIR Baseline Corrector — FTIR 基线校正参数建议
识别基线类型并推荐校正方法和参数
"""
import logging
from typing import Optional, List, Dict, Any, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 基线类型数据库 ────────────────────────────────────────────────
BASELINE_TYPES: Dict[str, dict] = {
    "flat": {
        "name": "Flat / Horizontal Baseline",
        "description": "Ideal baseline: constant offset across spectrum. Usually indicates good sample preparation.",
        "common_causes": ["Clean ATR measurement", "Well-prepared KBr pellet", "Good background subtraction"],
        "correction_method": "Simple offset subtraction (linear, 0-order).",
        "recommended_params": {"method": "offset", "anchor_points": "any reference region"},
        "software_settings": {"oplab": "Baseline: Linear/0-order", "omnic": "Baseline: Flat correction",
                              "python": "scipy.signal.detrend or simple mean subtraction"},
        "difficulty": "trivial",
    },
    "tilted_linear": {
        "name": "Tilted Linear Baseline",
        "description": "Baseline has a constant slope across the spectral range. Common in transmission measurements with scattering.",
        "common_causes": ["Light scattering from particles", "Non-uniform film thickness",
                          "Source intensity drift", "Detector nonlinearity"],
        "correction_method": "Linear (1st-order polynomial) fit through anchor points in transparent regions.",
        "recommended_params": {"method": "linear_poly", "order": 1,
                               "anchor_regions": [[4000-3800], [2400-2200], [2000-1900], [1800-1770]]},
        "software_settings": {"oplab": "Baseline: Linear (1st order)", "omnic": "Baseline: Linear",
                              "python": "numpy.polyfit(order=1) on anchor points"},
        "difficulty": "easy",
    },
    "curved_concave_up": {
        "name": "Curved Baseline (Concave Up / U-shaped)",
        "description": "Baseline curves upward at edges. Very common in ATR and diffuse reflectance.",
        "common_causes": ["ATR penetration depth variation with wavelength",
                          "Diffuse reflectance scattering (Kubelka-Munk effect)",
                          "Refractive index dispersion (anomalous dispersion near strong bands)",
                          "Sample surface roughness"],
        "correction_method": "Polynomial or rubber-band baseline. Use higher-order polynomial (2-4) or concave hull method.",
        "recommended_params": {"method": "polynomial", "order": "3-4",
                               "algorithm": "rubber_band (concave hull) or modified polynomial",
                               "anchor_regions": "transparent regions between major absorption bands"},
        "software_settings": {"oplab": "Baseline: Polynomial (3rd order)", "omnic": "Baseline: Concave rubber band",
                              "python": "baselinelib.modpoly or scipy.optimize.curve_fit (3rd order poly)"},
        "difficulty": "moderate",
    },
    "curved_concave_down": {
        "name": "Curved Baseline (Concave Down / Inverted U)",
        "description": "Baseline peaks in the middle and drops at edges. Less common but occurs in specific conditions.",
        "common_causes": ["Thermal emission from hot sample", "Instrumental stray light pattern",
                          "Interference fringes (etalon effect)", "Fluorescence background (in Raman/IR combo)"],
        "correction_method": "Polynomial fit (2nd-4th order) or asymmetric least squares (ALS).",
        "recommended_params": {"method": "asymmetric_least_squares", "order": 2,
                               "als_params": {"lam": 1e5, "p": 0.001, "niter": 15}},
        "software_settings": {"oplab": "Baseline: Polynomial + manual adjustment",
                              "omnic": "Baseline: Concave rubber band (may need inversion)",
                              "python": "baselinelib.asls or custom ALS implementation"},
        "difficulty": "moderate",
    },
    "wavy_oscillating": {
        "name": "Wavy / Oscillating Baseline",
        "description": "Baseline shows sinusoidal-like oscillations superimposed on broader trend.",
        "common_causes": ["Interference fringes (channel spectra) from parallel surfaces",
                          "Etalon effect (multiple internal reflections)",
                          "Thin-film interference patterns"],
        "correction_method": "First correct broad trend (polynomial), then apply fringe removal (FFT filter, derivative, or apodization).",
        "recommended_params": {"method": "two_step",
                               "step1": "polynomial baseline (order 3-4)",
                               "step2": "fft_filter or Savitzky-Golay derivative to remove periodic component",
                               "fringe_period_estimate": "measure peak-to-peak distance in cm⁻¹"},
        "software_settings": {"oplab": "Apozidation → Baseline correction", "omnic": "Advanced ATR correction → Fringe removal",
                              "python": "scipy.fft for fringe period detection + filtering"},
        "difficulty": "advanced",
    },
    "stepped": {
        "name": "Stepped / Discontinuous Baseline",
        "description": "Baseline shows abrupt steps or discontinuities at certain wavenumbers.",
        "common_causes": ["Detector gain change / range switching",
                          "Source change (mid-IR vs near-IR source switch)",
                          "Beam splitter artifact", "Filter wheel transition point"],
        "correction_method": "Piecewise polynomial or spline fitting. Treat each segment independently, then smooth the junction.",
        "recommended_params": {"method": "piecewise_polynomial_or_spline",
                               "segments": "identify step locations; fit each segment separately",
                               "smoothing": "cubic spline through segment anchors"},
        "software_settings": {"oplab": "Manual multi-point baseline with segment handling",
                              "omnic": "Manual baseline with careful point selection at each segment",
                              "python": "scipy.interpolate.CubicSpline with knot placement at step boundaries"},
        "difficulty": "moderate-to-advanced",
    },
    "fluorescence_background": {
        "name": "Fluorescence Background (Rising Toward Low Wavenumber)",
        "description": "Baseline rises sharply toward low wavenumber (long wavelength). Common in NIR-excited systems.",
        "common_causes": ["Sample fluorescence (especially with 1064nm Nd:YAG excitation)",
                          "Thermal emission from dark/hot samples",
                          "Strong scattering tail extending into measured range"],
        "correction_method": "Asymmetric Least Squares (ALS) or penalized least squares optimized for one-sided backgrounds.",
        "recommended_params": {"method": "asymmetric_least_squares",
                               "als_params": {"lam": 1e6, "p": 0.001, "niter": 20},
                               "alternative": "Polynomial (order 5-6) if fluorescence is smooth"},
        "software_settings": {"oplab": "Baseline: High-order polynomial or automatic fluorescence removal",
                              "omnic": "Advanced baseline: Automatic fluorescence correction",
                              "python": "baselinelib.asls (optimized parameters) or airPLS algorithm"},
        "difficulty": "moderate",
    },
}


# ── 常见样品类型的典型基线特征 ────────────────────────────────────
SAMPLE_BASELINE_PROFILES: Dict[str, dict] = {
    "atr_crystal": {
        "name": "ATR (Single Reflection Crystal)",
        "typical_baseline": "curved_concave_up",
        "why": "ATR penetration depth ∝ λ; longer wavelengths penetrate deeper → more absorption offset.",
        "default_correction": "Rubber band / concave hull (most FTIR software default for ATR).",
        "tips": [
            "Always collect fresh background before sample measurement.",
            "Ensure good crystal contact (uniform pressure).",
            "Use ATR correction algorithm built into software as first step.",
            "For quantitative work, consider using internal standard band ratio method.",
        ],
    },
    "kbr_pellet": {
        "name": "KBr Pellet Transmission",
        "typical_baseline": "flat or tilted_linear",
        "why": "Transmission through homogeneous matrix; quality depends on grinding/mixing.",
        "default_correction": "Linear or low-order polynomial usually sufficient.",
        "tips": [
            "Grind sample:KBr thoroughly (~1:100 ratio) to <2μm particles.",
            "Ensure pellet is transparent (no visible scattering).",
            "Moisture in KBr gives broad O-H ~3400 — dry KBr at 110°C before use.",
            "Christiansen effects can cause anomalous baselines — reduce particle size.",
        ],
    },
    "thin_film": {
        "name": "Thin Film (NaCl/CaF₂/ZnSe windows)",
        "typical_baseline": "wavy_oscillating or tilted_linear",
        "why": "Thin films produce interference fringes (etalons) from internal reflections.",
        "default_correction": "Fringe removal + polynomial baseline.",
        "tips": [
            "Fringe spacing Δν̃ = 1/(2nd) where n=refractive index, d=film thickness.",
            "Wedge the slightly to suppress fringes (non-parallel windows).",
            "Use apodization function in software for fringe suppression.",
        ],
    },
    "diffuse_reflectance_drs": {
        "name": "Diffuse Reflectance (DRIFTS)",
        "typical_baseline": "curved_concave_up (strong curvature)",
        "why": "Kubelka-Munk transformation needed; strong scattering dependence on wavelength.",
        "default_correction": "Kubelka-Munk conversion + polynomial/rubber-band baseline.",
        "tips": [
            "Mix with KBr powder (~5-10% sample).",
            "Use Kubelka-Munk units (F(R)) not raw reflectance for quantitative analysis.",
            "Particle size affects both scattering and apparent band shape.",
        ],
    },
    "gas_cell": {
        "name": "Gas Cell Transmission",
        "typical_baseline": "flat or stepped",
        "why": "Clean gas cell with well-aligned windows gives flat baseline.",
        "default_correction": "Usually minimal correction needed; check for etalon fringes.",
        "tips": [
            "Evacuate/purge cell properly before background scan.",
            "Check for atmospheric CO₂/H₂O vapor contamination in single-beam spectrum.",
            "Path length affects band intensities proportionally.",
        ],
    },
    "seir_attenuated_total": {
        "name": "SEIR (Surface-Enhanced IR)",
        "typical_baseline": "curved_concave_up or complex",
        "why": "Metal substrate contribution + enhanced near-field effects.",
        "default_correction": "High-order polynomial or reference substrate subtraction.",
        "tips": [
            "Always measure bare substrate under identical conditions for reference.",
            "Enhancement is distance-dependent (<10nm from surface).",
            "Baseline may differ significantly from normal IR of same material.",
        ],
    },
}


@ChemMCPManager.register_tool
class FtirBaselineCorrector(BaseTool):
    """
    FTIR 基线校正工具：根据光谱类型和基线形状，推荐最佳校正方法、
    算法参数和软件操作步骤。
    """
    __version__ = "0.1.0"
    name = "FtirBaselineCorrector"
    func_name = "correct_ftir_baseline"
    description = "Recommend FTIR baseline correction methods, algorithms, and parameters based on spectrum type, baseline shape, and sample preparation method."
    implementation_description = "Contains a knowledge base of 7 baseline types with recommended correction methods (polynomial, rubber-band, ALS, FFT), algorithm-specific parameters, software settings for common FTIR packages (OP Lab, OMNIC, Python), and sample-type-specific guidance."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["FTIR", "Baseline Correction", "Spectroscopy", "Data Processing", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("spectrum_type", "str", "", "Sample/measurement type: 'atr_crystal', 'kbr_pellet', 'thin_film', 'diffuse_reflectance_drs', 'gas_cell', 'seir'."),
        ("baseline_type", "str", "", "Observed baseline shape: 'flat', 'tilted_linear', 'curved', 'curved_concave_up', 'curved_concave_down', 'wavy_oscillating', 'stepped', 'fluorescence_background', 'auto-detect'."),
        ("wavenumber_range", "list", "[4000, 400]", "Measured spectral range in cm⁻¹ [max, min]."),
        ("major_absorption_bands", "list", "[]", "Approximate positions of major absorption bands (to avoid as anchor points)."),
        ("software", "str", "", "FTIR software: 'oplab', 'omnic', 'python', 'matlab', 'origin', 'general'."),
        ("correction_goal", "str", "quantitative", "Goal: 'quantitative' (accurate area/height), 'qualitative' (identification), 'publishing' (figure-quality)."),
        ("quality_issues", "list", "[]", "Additional issues: ['noise_high', 'saturation', 'water_vapor', 'co2_atmospheric']."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., 'atr_crystal curved python' or 'kbr_pellet tilted oplab quantitative'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with recommended correction method, parameters, software settings, step-by-step protocol, and validation checklist."),
    ]

    examples = [
        {
            "code_input": {
                "spectrum_type": "atr_crystal",
                "baseline_type": "curved",
                "wavenumber_range": [4000, 400],
                "major_absorption_bands": [3400, 2920, 2850, 1710, 1600, 1500, 1450],
                "software": "python",
                "correction_goal": "quantitative",
                "quality_issues": [],
            },
            "text_input": {
                "input_params": "atr_crystal curved python",
            },
            "output": {
                "result": {
                    "mode": "baseline_correction_recommendation",
                    "spectrum_type": "ATR crystal",
                    "note": "Recommended baseline correction protocol.",
                }
            }
        },
        {
            "code_input": {
                "spectrum_type": "kbr_pellet",
                "baseline_type": "tilted_linear",
                "wavenumber_range": [4000, 400],
                "major_absorption_bands": [],
                "software": "omnic",
                "correction_goal": "qualitative",
                "quality_issues": ["noise_high"],
            },
            "text_input": {
                "input_params": "kbr_pellet tilted omnic qualitative noise_high",
            },
            "output": {
                "result": {
                    "mode": "baseline_correction_recommendation",
                    "spectrum_type": "KBr pellet",
                    "note": "KBr pellet baseline protocol.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _resolve_baseline_type(self, specified: str, spectrum_type: str) -> str:
        """If user says 'curved' or 'auto-detect', resolve to specific type."""
        s = specified.lower().strip()
        if s in ("auto", "auto-detect", ""):
            # Infer from spectrum type
            profile = SAMPLE_BASELINE_PROFILES.get(spectrum_type, {})
            return profile.get("typical_baseline", "curved_concave_up")
        if s == "curved":
            return "curved_concave_up"  # most common curved type
        mapping = {
            "flat": "flat", "horizontal": "flat",
            "tilted": "tilted_linear", "linear": "tilted_linear", "slope": "tilted_linear",
            "u_shaped": "curved_concave_up", "atypical": "curved_concave_up",
            "inverted_u": "curved_concave_down", "dome": "curved_concave_down",
            "wavy": "wavy_oscillating", "fringe": "wavy_oscillating", "interference": "wavy_oscillating",
            "step": "stepped", "discontinuous": "stepped",
            "fluorescence": "fluorescence_background", "rising": "fluorescence_background",
        }
        return mapping.get(s, s)

    def _run_base(self, spectrum_type: str = "", baseline_type: str = "",
                  wavenumber_range: list = None, major_absorption_bands: list = None,
                  software: str = "", correction_goal: str = "quantitative",
                  quality_issues: list = None) -> dict:

        if wavenumber_range is None:
            wavenumber_range = [4000, 400]
        if major_absorption_bands is None:
            major_absorption_bands = []
        if quality_issues is None:
            quality_issues = []

        resolved_type = self._resolve_baseline_type(baseline_type, spectrum_type)

        # Get baseline type info
        bl_info = BASELINE_TYPES.get(resolved_type, BASELINE_TYPES["curved_concave_up"])

        # Get sample type info
        sample_info = SAMPLE_BASELINE_PROFILES.get(spectrum_type, {})

        # Determine optimal anchor regions (avoiding absorption bands)
        full_range = wavenumber_range
        anchor_suggestions = self._suggest_anchors(full_range, major_absorption_bands)

        # Build recommendation
        rec_method = bl_info["correction_method"]
        params = dict(bl_info["recommended_params"])
        params["anchor_points"] = anchor_suggestions

        sw_settings = bl_info["software_settings"].get(software.lower(), bl_info["software_settings"].get("general", "Manual baseline correction"))

        # Quality-issue-specific advice
        quality_advice = self._get_quality_advice(quality_issues)

        return {"result": {
            "mode": "baseline_correction_recommendation",
            "input_summary": {
                "spectrum_type": spectrum_type or "not specified",
                "baseline_type_detected": resolved_type,
                "wavenumber_range_cm-1": full_range,
                "n_major_bands_avoided": len(major_absorption_bands),
                "software": software or "not specified",
                "goal": correction_goal,
                "quality_issues": quality_issues or ["none reported"],
            },
            "baseline_classification": {
                "type_name": bl_info["name"],
                "description": bl_info["description"],
                "common_causes": bl_info["common_causes"],
                "difficulty": bl_info["difficulty"],
            },
            "recommended_correction": {
                "method": rec_method,
                "parameters": params,
                "software_specific_settings": sw_settings,
            },
            "sample_specific_guidance": sample_info.get("tips", []) if sample_info else [
                "No specific guidance for this sample type. Use general recommendations above.",
            ],
            "anchor_point_selection": {
                "strategy": "Select points in transparent (non-absorbing) regions",
                "suggested_anchor_regions_cm-1": anchor_suggestions,
                "regions_to_avoid": major_absorption_bands,
                "n_anchors_recommended": max(5, len(anchor_suggestions)),
            },
            "quality_issue_advice": quality_advice,
            "validation_checklist": [
                "✓ After correction, verify that baseline regions are at or near zero intensity.",
                "✓ Check that real absorption bands retain their true shape (no distortion from over-correction).",
                "✓ Compare corrected spectrum against a known reference if available.",
                f"✓ For '{correction_goal}' goal: verify band areas/heights are reproducible (±3%).",
                "✓ Document all correction parameters for reproducibility.",
                "✓ Save both raw and corrected spectra (never overwrite original data).",
            ],
            "step_by_step_protocol": self._generate_protocol(
                resolved_type, software, correction_goal, anchor_suggestions),
        }}

    @staticmethod
    def _suggest_anchors(wavenumber_range: list, avoid_peaks: list) -> List[str]:
        """Suggest anchor point regions avoiding major absorption bands."""
        hi, lo = wavenumber_range

        # Standard transparent regions in mid-IR
        candidate_regions = [
            (4000, 3800),  # above O-H/N-H
            (2500, 2400),  # gap between C-H and triple bond
        ]

        # Add dynamic regions between absorption bands
        sorted_peaks = sorted([p for p in avoid_peaks if lo <= p <= hi])
        prev_boundary = hi
        for pk in sorted_peaks:
            if prev_boundary - pk > 80:  # gap > 80 cm⁻¹
                candidate_regions.append((prev_boundary, pk + 40))
            prev_boundary = pk - 40
        if prev_boundary - lo > 80:
            candidate_regions.append((prev_boundary, lo))

        # Format as readable strings
        return [f"[{r[0]:.0f}-{r[1]:.0f}]" for r in candidate_regions[:8]]

    @staticmethod
    def _get_quality_advice(issues: list) -> dict:
        advice = {}
        for issue in issues:
            il = issue.lower()
            if "noise" in il:
                advice["noise"] = (
                    "Apply smoothing BEFORE baseline correction "
                    "(Savitzky-Golay, window=9-17, poly order 2-3). "
                    "Do not over-smooth — preserve genuine spectral features."
                )
            elif "saturate" in il:
                advice["saturation"] = (
                    "Saturated bands cannot be recovered. Exclude saturated regions "
                    "from baseline fitting. Reduce sample amount or path length for re-measurement."
                )
            elif "water" in il:
                advice["water_vapor"] = (
                    "Atmospheric H₂O vapor: subtract reference water spectrum or "
                    "purge instrument with dry air/N₂. Key bands: 3756/3652 (bend+rot), "
                    "1590 (bending), 3756-3580 (stretch manifold)."
                )
            elif "co2" in il:
                advice["co2"] = (
                    "CO₂ atmospheric band at ~2350 and 667 cm⁻¹. Purge with dry N₂ "
                    "or subtract atmosphere reference. Usually narrow and easily handled."
                )
        return advice

    @staticmethod
    def _generate_protocol(bl_type: str, software: str, goal: str,
                           anchors: list) -> List[str]:
        sw = software.lower() if software else "general"

        protocol = [
            f"1. Load raw spectrum into {sw} software.",
            "2. Inspect full-range spectrum visually; note baseline shape and problem regions.",
            f"3. Select baseline correction type appropriate for '{bl_type}'.",
        ]

        if goal == "quantitative":
            protocol.append("4a. For QUANTITATIVE analysis:")
            protocol.append("   - Choose consistent anchor points across all samples.")
            protocol.append("   - Use same correction parameters for sample and reference.")
            protocol.append("   - Validate with standard of known concentration.")
        else:
            protocol.append("4b. For QUALITATIVE/PUBLISHING analysis:")
            protocol.append("   - Prioritize visual appearance while preserving band shapes.")
            protocol.append("   - May use higher-order polynomial for cleaner look.")

        protocol.extend([
            f"5. Place anchor points in transparent regions: {', '.join(anchors[:5])}",
            "6. Apply baseline correction; inspect result.",
            "7. If over-corrected (bands go negative): reduce polynomial order or adjust anchors.",
            "8. If under-corrected (residual slope/curve): add more anchors or increase order.",
            "9. Export corrected spectrum (new file, never overwrite raw data).",
            "10. Record all parameters in lab notebook/LIMS.",
        ])
        return protocol

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            stype = parts[0] if parts else ""
            btype = parts[1] if len(parts) > 1 else ""
            sw = parts[2] if len(parts) > 2 else ""
            goal = parts[3] if len(parts) > 3 else "quantitative"
            issues = parts[4:] if len(parts) > 4 else []
            return self._run_base(spectrum_type=stype, baseline_type=btype,
                                   software=sw, correction_goal=goal,
                                   quality_issues=issues)
        except IndexError:
            raise ChemMCPError(f"Failed to parse text input '{input_params}'. "
                                f"Need: spectrum_type [baseline_type] [software] [goal] [issue1 ...]")
