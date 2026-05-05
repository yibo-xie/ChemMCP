"""
Fluorescence Quantum Yield Calculator — 荧光量子产率计算工具 (#322)

功能：
  1. 比较法 (Comparative method): Φ_s = Φ_r × (I_s/I_r) × (A_r/A_s) × (η_s²/η_r²)
  2. 绝对法 (Absolute method): Φ = 发射光子数 / 吸收光子数
"""

import logging
import math
from typing import Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见标准参考物质的量子产率 ────────────────────────────────
REFERENCE_STANDARDS: Dict[str, Dict[str, float]] = {
    "quinine_sulfate_0.5M_H2SO4": {
        "quantum_yield": 0.54,
        "refractive_index": 1.34,
        "excitation_nm": 348,
        "emission_nm": 450,
    },
    "fluorescein_0.1M_NaOH": {
        "quantum_yield": 0.93,
        "refractive_index": 1.34,
        "excitation_nm": 494,
        "emission_nm": 514,
    },
    "rhodamine_6G_ethanol": {
        "quantum_yield": 0.95,
        "refractive_index": 1.36,
        "excitation_nm": 488,
        "emission_nm": 550,
    },
    "rhodamine_B_ethanol": {
        "quantum_yield": 0.49,
        "refractive_index": 1.36,
        "excitation_nm": 554,
        "emission_nm": 577,
    },
    "910_anthracene_ethanol": {
        "quantum_yield": 0.27,
        "refractive_index": 1.36,
        "excitation_nm": 356,
        "emission_nm": 402,
    },
    "tsp_water": {   # tris(2,2'-bipyridyl)ruthenium(II) chloride
        "quantum_yield": 0.04,
        "refractive_index": 1.33,
        "excitation_nm": 452,
        "emission_nm": 610,
    },
}


@ChemMCPManager.register_tool
class FluorescenceQuantumYield(BaseTool):
    """
    荧光量子产率计算工具。
    支持比较法和绝对法两种计算模式。
    """
    __version__                = "0.1.0"
    name                       = "FluorescenceQuantumYield"
    func_name                  = "calculate_quantum_yield"
    description                = ("Calculate fluorescence quantum yield using comparative method "
                                 "(with reference standard) or absolute method.")
    implementation_description = (
        "Comparative method formula: Φ_s = Φ_r × (I_s/I_r) × (A_r/A_s) × (η_s²/η_r²). "
        "Absolute method formula: Φ = N_emitted / N_absorbed."
    )
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Fluorescence", "Photoluminescence", "Spectroscopy",
                                   "Quantum Yield", "Analytical Chemistry"]
    required_envs              = []

    code_input_sig = [
        ("method",                      "str",   "comparative", "Calculation method: 'comparative' or 'absolute'."),
        # --- comparative method params ---
        ("sample_integrated_area",      "float", "N/A",         "Integrated fluorescence emission area of sample (a.u.)."),
        ("reference_integrated_area",   "float", "N/A",         "Integrated fluorescence emission area of reference standard (a.u.)."),
        ("sample_absorbance",           "float", "N/A",         "Absorbance of sample at excitation wavelength."),
        ("reference_absorbance",        "float", "N/A",         "Absorbance of reference at excitation wavelength."),
        ("sample_refractive_index",     "float", "1.33",        "Refractive index of sample solvent."),
        ("reference_refractive_index",  "float", "1.33",        "Refractive index of reference solvent."),
        ("reference_standard_name",     "str",   "quinine_sulfate_0.5M_H2SO4", "Name of the reference standard."),
        # --- absolute method params ---
        ("photons_emitted",             "float", "N/A",         "Number of emitted photons (absolute method)."),
        ("photons_absorbed",            "float", "N/A",         "Number of absorbed photons (absolute method)."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated string. For comparative: 'comparative I_s I_r A_s A_r [n_s n_r std_name]'. "
         "For absolute: 'absolute N_emit N_absorb'."),
    ]

    output_sig = [
        ("quantum_yield",               "float", "Calculated fluorescence quantum yield (Φ), dimensionless, 0–1 range."),
        ("method_used",                 "str",   "Which method was used for calculation."),
        ("details",                     "dict",  "Intermediate values and correction factors applied."),
    ]

    examples = [
        {
            "code_input": {
                "method": "comparative",
                "sample_integrated_area": 850000.0,
                "reference_integrated_area": 920000.0,
                "sample_absorbance": 0.052,
                "reference_absorbance": 0.050,
                "sample_refractive_index": 1.33,
                "reference_refractive_index": 1.34,
                "reference_standard_name": "quinine_sulfate_0.5M_H2SO4",
            },
            "text_input": {"input_params": "comparative 850000 920000 0.052 0.050 1.33 1.34 quinine_sulfate_0.5M_H2SO4"},
            "output": {
                "quantum_yield": 0.496,
                "method_used": "comparative",
                "details": {},
            },
        },
        {
            "code_input": {
                "method": "absolute",
                "photons_emitted": 42000.0,
                "photons_absorbed": 100000.0,
            },
            "text_input": {"input_params": "absolute 42000 100000"},
            "output": {
                "quantum_yield": 0.42,
                "method_used": "absolute",
                "details": {},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load reference standards database."""
        self.standards = REFERENCE_STANDARDS

    def _run_comparative(
        self,
        sample_area: float,
        ref_area: float,
        sample_abs: float,
        ref_abs: float,
        sample_ri: float = 1.33,
        ref_ri: float = 1.33,
        ref_std_name: str = "quinine_sulfate_0.5M_H2SO4",
    ) -> Dict[str, Any]:
        """Comparative (relative) quantum yield calculation."""
        # Validate inputs
        if ref_area <= 0:
            raise ChemMCPError("Reference integrated area must be positive.")
        if sample_abs <= 0:
            raise ChemMCPError("Sample absorbance must be positive.")
        if ref_abs <= 0:
            raise ChemMCPError("Reference absorbance must be positive.")

        # Get reference quantum yield
        if ref_std_name not in self.standards:
            available = ", ".join(sorted(self.standards.keys()))
            raise ChemMCPError(f"Unknown reference standard '{ref_std_name}'. Available: {available}")
        phi_r = self.standards[ref_std_name]["quantum_yield"]

        # Refractive index correction factor
        ri_correction = (sample_ri ** 2) / (ref_ri ** 2)

        # Main formula: Φ_s = Φ_r × (I_s/I_r) × (A_r/A_s) × (η_s²/η_r²)
        phi_s = phi_r * (sample_area / ref_area) * (ref_abs / sample_abs) * ri_correction

        logger.info(
            f"Comparative QY: Φ={phi_r} × ({sample_area}/{ref_area}) × "
            f"({ref_abs}/{sample_abs}) × ({ri_correction:.4f}) = {phi_s:.6f}"
        )
        return {
            "quantum_yield": round(phi_s, 6),
            "method_used": "comparative",
            "details": {
                "phi_reference": phi_r,
                "intensity_ratio": round(sample_area / ref_area, 6),
                "absorbance_ratio": round(ref_abs / sample_abs, 6),
                "ri_correction_factor": round(ri_correction, 6),
                "reference_standard": ref_std_name,
            },
        }

    def _run_absolute(self, photons_emitted: float, photons_absorbed: float) -> Dict[str, Any]:
        """Absolute quantum yield calculation."""
        if photons_absorbed <= 0:
            raise ChemMCPError("Photons absorbed must be positive.")
        phi = photons_emitted / photons_absorbed
        if phi > 1.0:
            logger.warning(f"Quantum yield {phi:.4f} > 1.0; check input data or calibration.")
        logger.info(f"Absolute QY: {photons_emitted}/{photons_absorbed} = {phi:.6f}")
        return {
            "quantum_yield": round(min(phi, 1.0), 6),
            "method_used": "absolute",
            "details": {
                "photons_emitted": photons_emitted,
                "photons_absorbed": photons_absorbed,
            },
        }

    def _run_base(
        self,
        method: str = "comparative",
        # comparative
        sample_integrated_area: Optional[float] = None,
        reference_integrated_area: Optional[float] = None,
        sample_absorbance: Optional[float] = None,
        reference_absorbance: Optional[float] = None,
        sample_refractive_index: float = 1.33,
        reference_refractive_index: float = 1.33,
        reference_standard_name: str = "quinine_sulfate_0.5M_H2SO4",
        # absolute
        photons_emitted: Optional[float] = None,
        photons_absorbed: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Core logic dispatcher."""
        method = method.strip().lower()
        if method == "comparative":
            for name, val in [
                ("sample_integrated_area", sample_integrated_area),
                ("reference_integrated_area", reference_integrated_area),
                ("sample_absorbance", sample_absorbance),
                ("reference_absorbance", reference_absorbance),
            ]:
                if val is None:
                    raise ChemMCPError(f"Parameter '{name}' is required for comparative method.")
            return self._run_comparative(
                sample_area=sample_integrated_area,
                ref_area=reference_integrated_area,
                sample_abs=sample_absorbance,
                ref_abs=reference_absorbance,
                sample_ri=sample_refractive_index,
                ref_ri=reference_refractive_index,
                ref_std_name=reference_standard_name,
            )
        elif method == "absolute":
            if photons_emitted is None or photons_absorbed is None:
                raise ChemMCPError("Parameters 'photons_emitted' and 'photons_absorbed' are required for absolute method.")
            return self._run_absolute(photons_emitted, photons_absorbed)
        else:
            raise ChemMCPError(f"Unknown method '{method}'. Use 'comparative' or 'absolute'.")

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            method = parts[0].strip().lower()

            if method == "comparative":
                if len(parts) < 5:
                    raise ValueError("Comparative needs at least: method I_s I_r A_s A_r")
                return self._run_base(
                    method="comparative",
                    sample_integrated_area=float(parts[1]),
                    reference_integrated_area=float(parts[2]),
                    sample_absorbance=float(parts[3]),
                    reference_absorbance=float(parts[4]),
                    sample_refractive_index=float(parts[5]) if len(parts) > 5 else 1.33,
                    reference_refractive_index=float(parts[6]) if len(parts) > 6 else 1.33,
                    reference_standard_name=parts[7] if len(parts) > 7 else "quinine_sulfate_0.5M_H2SO4",
                )
            elif method == "absolute":
                if len(parts) < 3:
                    raise ValueError("Absolute needs: method N_emit N_absorb")
                return self._run_base(
                    method="absolute",
                    photons_emitted=float(parts[1]),
                    photons_absorbed=float(parts[2]),
                )
            else:
                raise ValueError(f"Unknown method: {method}")
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
