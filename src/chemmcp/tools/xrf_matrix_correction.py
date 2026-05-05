"""
XRF Matrix Correction — XRF基体校正系数计算工具 (#329)

功能：
  计算X射线荧光光谱(XRF)分析中的基体效应校正系数，
  支持基本参数法(FP)、经验α系数法和康普顿散射比率法。
"""

import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 元素X射线数据 ──────────────────────────────────────────────
XRF_ELEMENT_DATA: Dict[str, Dict[str, Any]] = {
    "na": {"Z": 11, "density": 0.97, "ka_keV": 1.041, "kb_keV": 1.072},
    "mg": {"Z": 12, "density": 1.74, "ka_keV": 1.254, "kb_keV": 1.297},
    "al": {"Z": 13, "density": 2.70, "ka_keV": 1.487, "kb_keV": 1.557},
    "si": {"Z": 14, "density": 2.33, "ka_keV": 1.740, "kb_keV": 1.836},
    "p":  {"Z": 15, "density": 1.82, "ka_keV": 2.013, "kb_keV": 2.139},
    "s":  {"Z": 16, "density": 2.07, "ka_keV": 2.308, "kb_keV": 2.464},
    "k":  {"Z": 19, "density": 0.86, "ka_keV": 3.314, "kb_keV": 3.590},
    "ca": {"Z": 20, "density": 1.55, "ka_keV": 3.692, "kb_keV": 4.013},
    "ti": {"Z": 22, "density": 4.51, "ka_keV": 4.511, "kb_keV": 4.932},
    "v":  {"Z": 23, "density": 6.11, "ka_keV": 4.952, "kb_keV": 5.427},
    "cr": {"Z": 24, "density": 7.19, "ka_keV": 5.415, "kb_keV": 5.947},
    "mn": {"Z": 25, "density": 7.43, "ka_keV": 5.899, "kb_keV": 6.491},
    "fe": {"Z": 26, "density": 7.87, "ka_keV": 6.404, "kb_keV": 7.058},
    "co": {"Z": 27, "density": 8.90, "ka_keV": 6.930, "kb_keV": 7.649},
    "ni": {"Z": 28, "density": 8.91, "ka_keV": 7.478, "kb_keV": 8.265},
    "cu": {"Z": 29, "density": 8.96, "ka_keV": 8.048, "kb_keV": 8.905},
    "zn": {"Z": 30, "density": 7.14, "ka_keV": 8.639, "kb_keV": 9.572},
    "as": {"Z": 33, "density": 5.73, "ka_keV": 10.544, "kb_keV": 11.727},
    "sr": {"Z": 38, "density": 2.64, "ka_keV": 14.165, "kb_keV": 15.836},
    "zr": {"Z": 40, "density": 6.52, "ka_keV": 15.777, "kb_keV": 17.666},
    "mo": {"Z": 42, "density": 10.22, "ka_keV": 17.479, "kb_keV": 19.608},
    "cd": {"Z": 48, "density": 8.65, "ka_keV": 23.174, "kb_keV": 26.093},
    "sb": {"Z": 51, "density": 6.68, "ka_keV": 26.359, "kb_keV": 29.802},
    "ba": {"Z": 56, "density": 3.59, "ka_keV": 32.194, "kb_keV": 36.378},
    "pb": {"Z": 82, "density": 11.34, "ka_keV": 74.970, "kb_keV": 85.000},
}


def _approximate_mac(element_z: int, photon_energy_keV: float) -> float:
    """Approximate mass absorption coefficient (cm²/g). μ/ρ ∝ Z⁴/(A·E³)."""
    if photon_energy_keV <= 0:
        return 0.0
    a = 2.0 * element_z * 1.02
    c = 9.0e-3
    mac = c * (element_z ** 4) / (a * (photon_energy_keV ** 3))
    return mac


@ChemMCPManager.register_tool
class XrfMatrixCorrection(BaseTool):
    """
    XRF基体校正系数计算工具。
    计算吸收-增强效应的校正系数，支持多种校正算法。
    """
    __version__                = "0.1.0"
    name                       = "XrfMatrixCorrection"
    func_name                  = "calculate_matrix_correction"
    description                = ("Calculate X-ray fluorescence (XRF) matrix correction coefficients "
                                 "using fundamental parameters, empirical alpha coefficients, or Compton ratio methods.")
    implementation_description = ("Implements three matrix correction approaches: (1) Fundamental Parameters (FP) "
                                 "using mass absorption coefficients, (2) Lachance-Traill (alpha coefficient) method, "
                                 "(3) Compton scatter ratio method for absorption compensation.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["XRF", "X-ray Fluorescence", "Matrix Correction",
                                   "Analytical Chemistry", "Quantitative Analysis"]
    required_envs              = []

    code_input_sig = [
        ("analyte_element",             "str",   "N/A",       "Analyte element symbol (e.g., 'Fe', 'Pb', 'Zn')."),
        ("matrix_elements",             "list",  "[]",        "List of matrix element symbols."),
        ("concentrations",              "dict",  "{}",        "Concentrations of all elements including analyte (weight fraction, 0-1)."),
        ("correction_method",           "str",   "fundamental", "Method: 'fundamental', 'empirical_alpha', or 'compton_ratio'."),
        ("incident_energy_keV",         "float", "None",      "Incident X-ray energy in keV (optional; auto-selected if None)."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated: analyte matrix_elem1,matrix_elem2,... [method] [energy_keV]"),
    ]

    output_sig = [
        ("correction_coefficients",     "dict",  "Alpha coefficients or FP correction factors per element pair."),
        ("corrected_concentration",     "float", "Matrix-corrected concentration of the analyte."),
        ("absorption_enhancement_effects","list","List of identified absorption and enhancement effects."),
        ("matrix_effect_summary",       "str",   "Human-readable summary of matrix effects and their magnitude."),
        ("total_correction_factor",     "float", "Overall correction factor applied to raw intensity."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_element": "Fe",
                "matrix_elements": ["Si", "Al", "Ca"],
                "concentrations": {"Fe": 0.10, "Si": 0.40, "Al": 0.15, "Ca": 0.05},
                "correction_method": "fundamental",
            },
            "text_input": {"input_params": "Fe Si,Al,Ca fundamental"},
            "output": {
                "total_correction_factor": 1.0,
                "matrix_effect_summary": "",
            },
        },
        {
            "code_input": {
                "analyte_element": "Cr",
                "matrix_elements": ["Fe", "Ni"],
                "concentrations": {"Cr": 0.18, "Fe": 0.70, "Ni": 0.10},
                "correction_method": "empirical_alpha",
            },
            "text_input": {"input_params": "Cr Fe,Ni empirical_alpha"},
            "output": {
                "total_correction_factor": 1.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load element database."""
        self.elem_data = XRF_ELEMENT_DATA

    def _run_fundamental(
        self,
        analyte_key: str,
        matrix_keys: List[str],
        concs: Dict[str, float],
        energy_keV: Optional[float],
    ) -> Dict[str, Any]:
        """Fundamental Parameters (FP) method for matrix correction."""
        # Get analyte data
        if analyte_key not in self.elem_data:
            raise ChemMCPError(f"Element '{analyte_key}' not in XRF database.")

        analyte = self.elem_data[analyte_key]
        z_a = analyte["Z"]

        # Auto-select incident energy if not provided (use Rh Kα = 20.216 keV as common tube line)
        if energy_keV is None:
            energy_keV = 20.216

        # Analyte fluorescence energy
        e_fluor = analyte.get("ka_keV", 6.4)

        # Calculate total mass absorption coefficient of the matrix at:
        # (a) incident energy E0 (for primary beam absorption)
        # (b) fluorescence energy Ef (for outgoing beam absorption)
        mu_total_incident = 0.0
        mu_total_fluor = 0.0

        all_elements = set([analyte_key] + matrix_keys)
        for elem_key in all_elements:
            c = concs.get(elem_key, 0.0)
            if elem_key in self.elem_data:
                z_e = self.elem_data[elem_key]["Z"]
                mu_i = _approximate_mac(z_e, energy_keV)
                mu_f = _approximate_mac(z_e, e_fluor)
                mu_total_incident += c * mu_i
                mu_total_fluor += c * mu_f

        # Absorption correction factor A = 1 / (mu_rho at E0 + mu_rho at Ef)
        # Simplified: use total absorption effect
        absorption_factor = mu_total_incident + mu_total_fluor

        # Enhancement check: does any matrix element have its edge just above the analyte edge?
        enhancement_effects = []
        for mk in matrix_keys:
            if mk in self.elem_data:
                me = self.elem_data[mk]
                mk_ka = me.get("ka_keV", 0)
                # If matrix element can emit X-rays that excite the analyte (enhancement)
                if mk_ka > analyte.get("ka_edge", analyte.get("ka_keV", 0) * 0.9):
                    # Check if matrix Ka energy > analyte K-edge
                    if mk_ka > e_fluor * 0.8:  # rough threshold
                        enhancement_effects.append(
                            f"Secondary enhancement possible from {mk.upper()} Kα ({mk_ka:.2f} keV) "
                            f"exciting {analyte_key.upper()}"
                        )

        # Absorption effects list
        absorption_effects = []
        for mk in matrix_keys:
            if mk in self.elem_data:
                z_m = self.elem_data[mk]["Z"]
                c_mk = concs.get(mk, 0)
                if c_mk > 0.01:  # Only significant concentrations
                    mu_m_i = _approximate_mac(z_m, energy_keV)
                    mu_m_f = _approximate_mac(z_m, e_fluor)
                    absorption_effects.append(
                        f"{mk.upper()} (c={c_mk:.1%}): μ/ρ@E₀≈{mu_m_i:.1f} cm²/g, "
                        f"μ/ρ@E_f≈{mu_m_f:.1f} cm²/g"
                    )

        # Total correction factor (simplified model)
        # In real FP: C_i = (I_measured / sensitivity_i) × absorption_factor × enhancement_factor
        # Here we compute a relative correction factor
        raw_conc = concs.get(analyte_key, 0.1)
        # Simplified: correction ≈ 1 / exp(-μρ × thickness) → approximate as linear factor
        # For thin samples: C_corrected ≈ C_raw × (1 + Σα_j×C_j)
        alpha_coeffs = {}
        total_alpha_correction = 0.0
        for mk in matrix_keys:
            if mk in self.elem_data:
                z_m = self.elem_data[mk]["Z"]
                c_mk = concs.get(mk, 0)
                # Alpha coefficient approximation: α_ij ∝ (μ_j(E0) - μ_j(Ef)) / μ_i(Ef)
                mu_j_i = _approximate_mac(z_m, energy_keV)
                mu_j_f = _approximate_mac(z_m, e_fluor)
                mu_a_f = _approximate_mac(z_a, e_fluor)
                if mu_a_f > 0:
                    alpha = (mu_j_i - mu_j_f) / max(mu_a_f, 0.1)
                else:
                    alpha = 0.0
                alpha_coeffs[f"{analyte_key.upper()}-{mk.upper()}"] = round(alpha, 4)
                total_alpha_correction += alpha * c_mk

        correction_factor = 1.0 + total_alpha_correction
        corrected_conc = raw_conc * correction_factor

        # Summary
        summary_parts = [
            f"XRF Matrix Correction (FP method) for **{analyte_key.upper()}**:",
            f"  Incident energy: {energy_keV:.2f} keV | Fluorescence energy: {e_fluor:.2f} keV",
            f"  Matrix absorption (μ/ρ): incident={mu_total_incident:.2f}, fluorescent={mu_total_fluor:.2f}",
        ]
        if abs(total_alpha_correction) > 0.01:
            summary_parts.append(f"  Net α-correction: {total_alpha_correction:+.4f} (factor={correction_factor:.4f})")
        else:
            summary_parts.append("  Matrix effect is negligible (< 1%).")
        if enhancement_effects:
            summary_parts.append(f"  ⚠️ Enhancement effects: {'; '.join(enhancement_effects)}")

        logger.info(f"XRF-FP correction for {analyte_key.upper()}: factor={correction_factor:.4f}")
        return {
            "correction_coefficients": alpha_coeffs,
            "corrected_concentration": round(corrected_conc, 6),
            "absorption_enhancement_effects": absorption_effects + enhancement_effects,
            "matrix_effect_summary": "\n".join(summary_parts),
            "total_correction_factor": round(correction_factor, 6),
        }

    def _run_empirical_alpha(
        self,
        analyte_key: str,
        matrix_keys: List[str],
        concs: Dict[str, float],
    ) -> Dict[str, Any]:
        """Lachance-Traill empirical alpha coefficient method."""
        # Use literature/typical alpha values for common pairs
        # α_ij represents the influence of element j on analyte i
        # Positive α = absorption effect (signal depression)
        # Negative α = enhancement effect (signal increase)

        # Simplified influence coefficients based on atomic number differences
        if analyte_key not in self.elem_data:
            raise ChemMCPError(f"Element '{analyte_key}' not in XRF database.")
        z_a = self.elem_data[analyte_key]["Z"]

        alpha_coeffs = {}
        total_correction = 0.0
        effects_list = []

        for mk in matrix_keys:
            if mk in self.elem_data:
                z_m = self.elem_data[mk]["Z"]
                c_mk = concs.get(mk, 0)

                # Empirical alpha estimation based on Z difference
                dz = z_m - z_a
                if abs(dz) <= 2:
                    alpha = 0.0  # Similar Z: minimal interaction
                elif dz > 2:
                    # Heavier matrix element: absorbs analyte fluorescence (positive α)
                    alpha = 0.5 * math.log10(1 + abs(dz)) * (1 + c_mk)
                else:
                    # Lighter matrix element: may cause enhancement (negative α)
                    alpha = -0.3 * math.log10(1 + abs(dz)) * (c_mk ** 0.5)

                pair_key = f"{analyte_key.upper()}-{mk.upper()}"
                alpha_coeffs[pair_key] = round(alpha, 4)
                total_correction += alpha * c_mk

                if abs(alpha) > 0.1:
                    direction = "absorption (signal ↓)" if alpha > 0 else "enhancement (signal ↑)"
                    effects_list.append(f"{mk.upper()}: α={alpha:+.3f} [{direction}]")

        raw_conc = concs.get(analyte_key, 0.1)
        correction_factor = 1.0 / (1.0 + total_correction) if (1.0 + total_correction) != 0 else 1.0
        corrected_conc = raw_conc * correction_factor

        summary = (
            f"Empirical α-coefficient correction for {analyte_key.upper()}: "
            f"Σ(α_j×C_j)={total_correction:+.4f}, factor={correction_factor:.4f}"
        )

        logger.info(f"XRF-alpha correction for {analyte_key.upper()}: factor={correction_factor:.4f}")
        return {
            "correction_coefficients": alpha_coeffs,
            "corrected_concentration": round(corrected_conc, 6),
            "absorption_enhancement_effects": effects_list,
            "matrix_effect_summary": summary,
            "total_correction_factor": round(correction_factor, 6),
        }

    def _run_compton_ratio(
        self,
        analyte_key: str,
        matrix_keys: List[str],
        concs: Dict[str, float],
    ) -> Dict[str, Any]:
        """Compton scatter ratio method for absorption compensation."""
        if analyte_key not in self.elem_data:
            raise ChemMCPError(f"Element '{analyte_key}' not in XRF database.")

        # Compton ratio method: uses ratio of Compton peak to Rayleigh peak
        # to estimate average mass absorption coefficient
        # R_C/R_R ∝ <Z> of the sample

        # Compute average Z of matrix
        all_elems = set([analyte_key] + matrix_keys)
        avg_z = 0.0
        total_c = 0.0
        for ek in all_elems:
            c = concs.get(ek, 0)
            if ek in self.elem_data:
                avg_z += self.elem_data[ek]["Z"] * c
                total_c += c

        avg_z = avg_z / max(total_c, 0.01)

        # Mass absorption coefficient scales roughly with Z^3–4
        # Higher average Z → more absorption → lower signal
        # Correction: multiply by (avg_z / reference_z)^power
        reference_z = 14.0  # Si as typical reference (e.g., glass standard)
        power = 3.5  # Empirical exponent
        compton_correction = (reference_z / max(avg_z, 1.0)) ** power

        raw_conc = concs.get(analyte_key, 0.1)
        corrected_conc = raw_conc / max(compton_correction, 0.1)

        effects = [f"Average matrix Z = {avg_z:.1f}", f"Compton-based absorption correction: {compton_correction:.4f}"]
        summary = (
            f"Compton ratio correction for {analyte_key.upper()}: "
            f"<Z>={avg_z:.1f}, correction={compton_correction:.4f}"
        )

        return {
            "correction_coefficients": {"compton_ratio": round(compton_correction, 4)},
            "corrected_concentration": round(corrected_conc, 6),
            "absorption_enhancement_effects": effects,
            "matrix_effect_summary": summary,
            "total_correction_factor": round(1.0 / max(compton_correction, 0.001), 6),
        }

    def _run_base(
        self,
        analyte_element: str,
        matrix_elements: Optional[List[str]] = None,
        concentrations: Optional[Dict[str, float]] = None,
        correction_method: str = "fundamental",
        incident_energy_keV: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Core logic dispatcher."""
        key = analyte_element.strip().lower()
        matrix = [m.strip().lower() for m in (matrix_elements or [])]
        concs = dict(concentrations) if concentrations else {}

        # Ensure analyte has a default concentration
        if key not in concs:
            concs[key] = 0.1  # Assume 10% for calculation

        method = correction_method.strip().lower()

        if method == "fundamental":
            return self._run_fundamental(key, matrix, concs, incident_energy_keV)
        elif method == "empirical_alpha" or method == "empirical":
            return self._run_empirical_alpha(key, matrix, concs)
        elif method == "compton_ratio" or method == "compton":
            return self._run_compton_ratio(key, matrix, concs)
        else:
            raise ChemMCPError(
                f"Unknown correction method '{correction_method}'. "
                "Use 'fundamental', 'empirical_alpha', or 'compton_ratio'."
            )

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            analyte = parts[0]
            matrix = []
            method = "fundamental"
            energy = None

            idx = 1
            while idx < len(parts):
                p = parts[idx]
                if p in ("fundamental", "empirical_alpha", "empirical", "compton_ratio", "compton"):
                    method = p
                elif "," in p:
                    matrix.extend([x.strip() for x in p.split(",") if x.strip()])
                else:
                    try:
                        energy = float(p)
                    except ValueError:
                        matrix.append(p.strip())
                idx += 1

            return self._run_base(analyte, matrix or None, None, method, energy)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
