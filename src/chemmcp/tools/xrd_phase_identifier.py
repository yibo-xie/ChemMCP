"""
XRD Phase Identifier — XRD物相鉴定辅助工具 (#330)
"""

import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


XRD_PHASE_DB: Dict[str, Dict[str, Any]] = {
    "quartz": {
        "pdf_id": "00-046-1045", "name": "Quartz (α-SiO₂)", "formula": "SiO₂",
        "crystal_system": "Trigonal", "space_group": "P3₂21",
        "lattice": {"a": 4.913, "b": 4.913, "c": 5.405, "alpha": 90, "beta": 90, "gamma": 120},
        "color": "colorless/white", "density_g_cm3": 2.65,
        "peaks": [
            (20.85, 22, "100", 4.255), (26.64, 100, "101", 3.343),
            (36.54, 17, "110", 2.458), (39.46, 8, "102", 2.282),
            (42.44, 6, "111", 2.128), (50.14, 13, "112", 1.818),
            (54.87, 3, "103", 1.672), (59.96, 8, "211", 1.541),
            (67.74, 8, "212", 1.381), (68.14, 10, "203", 1.375),
        ],
    },
    "cristobalite": {
        "pdf_id": "00-039-1425", "name": "Cristobalite (β-SiO₂)", "formula": "SiO₂",
        "crystal_system": "Tetragonal", "space_group": "P4₁2₁2",
        "lattice": {"a": 4.971, "b": 4.971, "c": 6.928, "alpha": 90, "beta": 90, "gamma": 90},
        "color": "white", "density_g_cm3": 2.33,
        "peaks": [
            (21.98, 30, "101", 4.040), (28.42, 50, "111", 3.138),
            (31.46, 100, "110/101", 2.841), (36.14, 35, "102", 2.483),
            (41.25, 12, "200", 2.188), (56.53, 15, "212", 1.627),
        ],
    },
    "rutile": {
        "pdf_id": "00-021-1276", "name": "Rutile (TiO₂)", "formula": "TiO₂",
        "crystal_system": "Tetragonal", "space_group": "P4₂/mnm",
        "lattice": {"a": 4.593, "b": 4.593, "c": 2.959, "alpha": 90, "beta": 90, "gamma": 90},
        "color": "dark brown/black", "density_g_cm3": 4.25,
        "peaks": [
            (27.45, 100, "110", 3.247), (36.09, 55, "101", 2.487),
            (39.19, 20, "200", 2.297), (41.23, 25, "111", 2.189),
            (54.32, 63, "211", 1.687), (56.64, 30, "220", 1.624),
            (62.74, 16, "002", 1.480), (69.01, 15, "301", 1.360),
        ],
    },
    "anatase": {
        "pdf_id": "00-021-1272", "name": "Anatase (TiO₂)", "formula": "TiO₂",
        "crystal_system": "Tetragonal", "space_group": "I4₁/amd",
        "lattice": {"a": 3.784, "b": 3.784, "c": 9.515, "alpha": 90, "beta": 90, "gamma": 90},
        "color": "white/brownish", "density_g_cm3": 3.89,
        "peaks": [
            (25.31, 100, "101", 3.516), (37.80, 40, "004", 2.378),
            (38.57, 20, "200", 2.332), (48.05, 20, "105", 1.893),
            (53.89, 20, "211", 1.699), (55.06, 30, "204", 1.666),
            (70.31, 10, "215", 1.338), (75.03, 14, "303", 1.264),
        ],
    },
    "corundum": {
        "pdf_id": "00-046-1212", "name": "Corundum (α-Al₂O₃)", "formula": "Al₂O₃",
        "crystal_system": "Trigonal", "space_group": "R-3c",
        "lattice": {"a": 4.759, "b": 4.759, "c": 12.991, "alpha": 90, "beta": 90, "gamma": 120},
        "color": "white/red/blue", "density_g_cm3": 3.98,
        "peaks": [
            (25.58, 70, "012", 3.479), (35.15, 100, "104", 2.552),
            (37.77, 80, "110", 2.379), (43.36, 50, "113", 2.085),
            (52.55, 90, "024", 1.740), (57.50, 60, "116", 1.601),
            (61.31, 40, "214", 1.510), (66.52, 40, "300", 1.404),
        ],
    },
    "calcite": {
        "pdf_id": "00-047-1743", "name": "Calcite (CaCO₃)", "formula": "CaCO₃",
        "crystal_system": "Trigonal", "space_group": "R-3c",
        "lattice": {"a": 4.990, "b": 4.990, "c": 17.061, "alpha": 90, "beta": 90, "gamma": 120},
        "color": "colorless/white", "density_g_cm3": 2.71,
        "peaks": [
            (23.05, 18, "012", 3.855), (29.40, 100, "104", 3.035),
            (35.96, 14, "113", 2.495), (39.40, 18, "202", 2.285),
            (43.15, 18, "018", 2.095), (47.12, 17, "116", 1.928),
            (48.50, 18, "211", 1.876), (57.42, 6, "122", 1.604),
        ],
    },
    "hematite": {
        "pdf_id": "00-033-0664", "name": "Hematite (α-Fe₂O₃)", "formula": "Fe₂O₃",
        "crystal_system": "Trigonal", "space_group": "R-3c",
        "lattice": {"a": 5.035, "b": 5.035, "c": 13.747, "alpha": 90, "beta": 90, "gamma": 120},
        "color": "red/dark red", "density_g_cm3": 5.26,
        "peaks": [
            (24.14, 20, "012", 3.683), (33.15, 100, "104", 2.700),
            (35.61, 60, "110", 2.519), (40.85, 30, "113", 2.207),
            (49.48, 25, "024", 1.841), (54.09, 40, "116", 1.693),
            (62.45, 15, "214", 1.485),
        ],
    },
    "magnetite": {
        "pdf_id": "00-019-0629", "name": "Magnetite (Fe₃O₄)", "formula": "Fe₃O₄",
        "crystal_system": "Cubic", "space_group": "Fd-3m",
        "lattice": {"a": 8.396, "b": 8.396, "c": 8.396, "alpha": 90, "beta": 90, "gamma": 90},
        "color": "black, magnetic", "density_g_cm3": 5.18,
        "peaks": [
            (30.10, 32, "220", 2.967), (35.43, 100, "311", 2.532),
            (43.07, 20, "400", 2.101), (53.73, 28, "422", 1.705),
            (56.96, 18, "511/333", 1.615), (62.57, 14, "440", 1.482),
        ],
    },
    "nacl": {
        "pdf_id": "00-005-0628", "name": "Halite (NaCl)", "formula": "NaCl",
        "crystal_system": "Cubic", "space_group": "Fm-3m",
        "lattice": {"a": 5.640, "b": 5.640, "c": 5.640, "alpha": 90, "beta": 90, "gamma": 90},
        "color": "colorless/white", "density_g_cm3": 2.16,
        "peaks": [
            (27.36, 100, "111", 3.258), (31.70, 55, "200", 2.821),
            (45.45, 14, "220", 1.994), (53.87, 6, "311", 1.701),
            (56.47, 6, "222", 1.628), (66.26, 3, "400", 1.410),
        ],
    },
    "zno_wurtzite": {
        "pdf_id": "00-036-1451", "name": "Zincite (ZnO, wurtzite)", "formula": "ZnO",
        "crystal_system": "Hexagonal", "space_group": "P6₃mc",
        "lattice": {"a": 3.250, "b": 3.250, "c": 5.207, "alpha": 90, "beta": 90, "gamma": 120},
        "color": "white/pale yellow", "density_g_cm3": 5.61,
        "peaks": [
            (31.77, 100, "100", 2.815), (34.42, 44, "002", 2.603),
            (36.25, 21, "101", 2.476), (47.54, 5, "102", 1.912),
            (56.59, 11, "110", 1.626), (62.86, 11, "103", 1.477),
            (69.09, 6, "112", 1.360),
        ],
    },
}

RADIATION_WAVELENGTHS: Dict[str, float] = {
    "CuKa": 1.5406, "CuKa1": 1.5406, "CuKa2": 1.5444,
    "CoKa": 1.7890, "MoKa": 0.7107, "FeKa": 1.9360, "CrKa": 2.2909,
}


def _two_theta_to_d(two_theta_deg: float, wavelength_A: float) -> float:
    theta_rad = math.radians(two_theta_deg / 2.0)
    if math.sin(theta_rad) == 0:
        return float("inf")
    return wavelength_A / (2.0 * math.sin(theta_rad))


@ChemMCPManager.register_tool
class XrdPhaseIdentifier(BaseTool):
    __version__                = "0.1.0"
    name                       = "XrdPhaseIdentifier"
    func_name                  = "identify_phases"
    description                = ("Assist XRD phase identification by matching observed "
                                 "d-spacing/intensity patterns against a built-in reference database.")
    implementation_description = ("Uses Bragg's law for d-spacing conversion and pattern matching against "
                                 "an internal database of common phases. Computes Figure of Merit (FOM).")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["XRD", "X-ray Diffraction", "Phase Identification",
                                   "Crystallography", "Materials Science"]
    required_envs              = []

    code_input_sig = [
        ("observed_peaks",               "list",  "N/A",       "List of [{'two_theta_deg': float, 'intensity': float}, ...]."),
        ("radiation_type",               "str",   "CuKa",      "Radiation: 'CuKa', 'CoKa', 'MoKa', 'FeKa', 'CrKa'."),
        ("candidate_phases",             "list",  "None",      "List of phase names to search (None = all)."),
        ("tolerance_two_theta_deg",      "float", "0.05",      "Matching tolerance in degrees 2θ."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated 2theta,intensity pairs [radiation] [tolerance]. E.g.: '26.64,100 36.54,17 CuKa'"),
    ]

    output_sig = [
        ("matched_phases",               "list",  "List of matched phases with scores."),
        ("figure_of_merit",             "float", "Figure of Merit (lower is better; < 15 is good)."),
        ("unidentified_peaks",           "list",  "Unmatched observed peaks."),
        ("phase_details",                "dict",  "Detailed match info per phase."),
    ]

    examples = [
        {
            "code_input": {
                "observed_peaks": [
                    {"two_theta_deg": 26.64, "intensity": 100},
                    {"two_theta_deg": 20.85, "intensity": 22},
                    {"two_theta_deg": 36.54, "intensity": 17},
                    {"two_theta_deg": 50.14, "intensity": 13},
                ],
                "radiation_type": "CuKa",
            },
            "text_input": {"input_params": "26.64,100 20.85,22 36.54,17 50.14,13 CuKa"},
            "output": {
                "matched_phases": [{"phase_name": "Quartz (α-SiO₂)", "match_score": 0.95}],
                "figure_of_merit": 2.5,
            },
        },
        {
            "code_input": {
                "observed_peaks": [
                    {"two_theta_deg": 27.45, "intensity": 100},
                    {"two_theta_deg": 36.09, "intensity": 55},
                    {"two_theta_deg": 54.32, "intensity": 63},
                ],
                "radiation_type": "CuKa",
                "candidate_phases": ["rutile"],
            },
            "text_input": {"input_params": "27.45,100 36.09,55 54.32,63 CuKa rutile"},
            "output": {
                "matched_phases": [{"phase_name": "Rutile (TiO₂)", "match_score": 0.98}],
                "figure_of_merit": 1.2,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.phase_db = XRD_PHASE_DB
        self.wavelengths = RADIATION_WAVELENGTHS

    def _run_base(
        self,
        observed_peaks: List[Dict[str, float]],
        radiation_type: str = "CuKa",
        candidate_phases: Optional[List[str]] = None,
        tolerance_two_theta_deg: float = 0.05,
    ) -> Dict[str, Any]:
        """Core phase identification logic."""
        rad_key = radiation_type.replace("-", "").replace("_", "").upper()
        # Normalize radiation key
        norm_rad = None
        for k in self.wavelengths:
            if k.replace("-", "").replace("_", "").upper() == rad_key or \
               k.upper() == rad_key or radiation_type.upper() == k.upper():
                norm_rad = k
                break
        if not norm_rad:
            available = ", ".join(sorted(self.wavelengths.keys()))
            raise ChemMCPError(f"Unknown radiation type '{radiation_type}'. Available: {available}")
        wavelength_A = self.wavelengths[norm_rad]

        # Determine search set
        if candidate_phases:
            search_keys = [p.strip().lower().replace(" ", "_").replace("-", "_")
                          for p in candidate_phases]
        else:
            search_keys = list(self.phase_db.keys())

        # Convert observed peaks to d-spacing for robustness across radiation types
        obs_d_list = []
        for op in observed_peaks:
            tt = op.get("two_theta_deg", op.get("two_theta", 0))
            intensity = op.get("intensity", 1)
            d = _two_theta_to_d(tt, wavelength_A)
            obs_d_list.append({"two_theta": tt, "d_spacing": d, "intensity": intensity})

        # Match each phase
        phase_results = []
        for pkey in search_keys:
            if pkey not in self.phase_db:
                continue
            phase = self.phase_db[pkey]
            ref_peaks = phase["peaks"]

            matched_ref_peaks = []
            unmatched_obs = list(range(len(obs_d_list)))
            total_delta_2t = 0.0
            match_count = 0

            for ri, ref_peak in enumerate(ref_peaks):
                ref_tt = ref_peak[0]
                ref_int = ref_peak[1]
                ref_d = ref_peak[3] if len(ref_peak) > 3 else _two_theta_to_d(ref_tt, wavelength_A)

                best_match_idx = -1
                best_delta = tolerance_two_theta_deg * 2

                for oi, obs in enumerate(obs_d_list):
                    if oi not in unmatched_obs:
                        continue
                    # Match by d-spacing (more robust than 2θ when radiation differs)
                    delta_d = abs(obs["d_spacing"] - ref_d)
                    delta_2t = abs(obs["two_theta"] - ref_tt)
                    # Use d-spacing as primary criterion
                    d_tolerance = ref_d * 0.01  # 1% d-spacing tolerance
                    if delta_d < d_tolerance or delta_2t < tolerance_two_theta_deg:
                        if delta_2t < best_delta:
                            best_delta = delta_2t
                            best_match_idx = oi

                if best_match_idx >= 0:
                    matched_ref_peaks.append({
                        "ref_two_theta": ref_tt,
                        "ref_intensity": ref_int,
                        "ref_hkl": ref_peak[2],
                        "obs_two_theta": obs_d_list[best_match_idx]["two_theta"],
                        "obs_intensity": obs_d_list[best_match_idx]["intensity"],
                        "delta_2theta": round(best_delta, 4),
                        "d_spacing": round(ref_d, 4),
                    })
                    total_delta_2t += best_delta ** 2
                    unmatched_obs.remove(best_match_idx)
                    match_count += 1

            if match_count >= 2:  # Need at least 2 peak matches
                # Calculate match score
                rms = math.sqrt(total_delta_2t / max(match_count, 1))
                intensity_correlation = min(1.0, match_count / len(ref_peaks))
                match_score = round(intensity_correlation * (1.0 - min(rms, 1.0)), 3)

                phase_results.append((match_score, pkey, phase, matched_ref_peaks, match_count))

        # Sort by score descending
        phase_results.sort(key=lambda x: x[0], reverse=True)

        # Build output
        matched_phases = []
        details = {}
        for score, pkey, phase, matches, mcount in phase_results:
            entry = {
                "phase_name": phase["name"],
                "pdf_id": phase["pdf_id"],
                "formula": phase["formula"],
                "crystal_system": phase["crystal_system"],
                "lattice_params": phase["lattice"],
                "match_score": score,
                "matched_peaks_count": mcount,
                "total_ref_peaks": len(phase["peaks"]),
                "peak_matches": matches,
            }
            matched_phases.append(entry)
            details[phase["name"]] = entry

        # Figure of Merit (FOM) calculation — simplified Smith-Snyder FOM
        # FOM = (N_matched / N_observed) × (1 / avg_delta_2theta)
        n_obs = len(obs_d_list)
        if matched_phases:
            best = phase_results[0]
            avg_delta = math.sqrt(best[3] and sum(m["delta_2theta"]**2 for m in best[3]) / max(len(best[3]), 1) or 0)
            fom = (n_obs / max(n_obs, 1)) * (1.0 / max(avg_delta, 0.001)) * (1 - best[0])
            fom = round(fom, 2)
        else:
            fom = 999.0

        # Unidentified peaks
        unidentified = []
        for i, obs in enumerate(obs_d_list):
            identified = False
            for _, _, _, matches, _ in phase_results:
                for m in matches:
                    if abs(m["obs_two_theta"] - obs["two_theta"]) < tolerance_two_theta_deg * 2:
                        identified = True
                        break
                if identified:
                    break
            if not identified:
                unidentified.append({"two_theta_deg": obs["two_theta"], "d_spacing_A": round(obs["d_spacing"], 4),
                                   "intensity": obs["intensity"]})

        logger.info(f"XRD phase ID: {len(matched_phases)} phases matched, FOM={fom}, {len(unidentified)} unidentified peaks")
        return {
            "matched_phases": matched_phases,
            "figure_of_merit": fom,
            "unidentified_peaks": unidentified,
            "phase_details": details,
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            peaks = []
            radiation = "CuKa"
            candidates = None
            tolerance = 0.05

            idx = 0
            while idx < len(parts):
                p = parts[idx]
                if "," in p:
                    sub = p.split(",")
                    try:
                        peaks.append({"two_theta_deg": float(sub[0]), "intensity": float(sub[1])})
                    except (ValueError, IndexError):
                        pass
                elif p.upper() in ("CUKA", "CUKA1", "CUKA2", "COKA", "MOKA", "FEKA", "CRKA"):
                    radiation = p
                elif p in self.phase_db or p.replace("_", " ") in [v.get("name","") for v in self.phase_db.values()]:
                    candidates = [p]
                else:
                    try:
                        tolerance = float(p)
                    except ValueError:
                        candidates = [p]
                idx += 1

            if not peaks:
                raise ValueError("No peak data found. Format: '2theta,intensity 2theta,intensity ... [radiation]'")

            return self._run_base(peaks, radiation, candidates, tolerance)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. "
                               f"Format: '2theta,intensity 2theta,intensity ... [radiation] [tolerance] [phase_name]'")
