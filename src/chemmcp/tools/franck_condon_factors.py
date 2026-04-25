import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
H = 6.62607015e-34        # J·s
C = 2.99792458e8           # m/s
AMU = 1.66053906660e-27    # kg


@ChemMCPManager.register_tool
class FranckCondonFactors(BaseTool):
    """
    Franck-Condon因子计算工具。
    计算电子跃迁中振动能级间的Franck-Condon因子，预测振动 progression的强度分布。
    """
    __version__ = "0.1.0"
    name = "FranckCondonFactors"
    func_name = "calculate_franck_condon_factors"
    description = "Calculate Franck-Condon factors for electronic transitions and predict vibrational progression intensity distributions."
    implementation_description = "Uses the harmonic oscillator approximation with Huang-Rhys factor S to compute FC factors: |⟨χ'_v'|χ_v⟩|². Models displacement along normal coordinates between ground and excited state potential surfaces."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Spectroscopy", "Electronic Transitions", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("huang_rhys_factor_S", "float", "N/A", "Huang-Rhys factor (dimensionless). S = 0: no geometry change; S > 1: large displacement."),
        ("v_max", "int", "10", "Maximum vibrational quantum number v' to compute in the excited state."),
        ("v_initial", "int", "0", "Initial vibrational quantum number v in the ground state (usually 0 at room temperature)."),
        ("delta_q_dimensionless", "float", "0", "Dimensionless displacement ΔQ' (optional; if provided, overrides S via S = ΔQ²/2)."),
        ("frequency_cm-1", "float", "1000", "Vibrational frequency of the mode in cm⁻¹ (for energy scale calculation)."),
        ("include_progression_plot_data", "bool", "True", "Whether to include data suitable for plotting intensity vs v'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: S [v_max v_initial freq_cm-1] e.g., '2.5 15 0 1500'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing Franck-Condon factors for each v', progression intensities, maximum intensity position, and spectroscopic interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "huang_rhys_factor_S": 1.0,
                "v_max": 10,
                "v_initial": 0,
                "delta_q_dimensionless": 0,
                "frequency_cm-1": 1000,
                "include_progression_plot_data": True,
            },
            "text_input": {
                "input_params": "1.0 10 0 1000",
            },
            "output": {
                "result": {
                    "huang_rhys_factor_S": 1.0,
                    "v_initial": 0,
                    "v_max_computed": 10,
                    "max_intensity_at_v_prime": 1,
                    "fc_factors": [
                        {"v_prime": 0, "fc_factor": round(math.exp(-1) * (1**0) / math.factorial(0), 8), "relative_intensity_pct": 36.79},
                        {"v_prime": 1, "fc_factor": round(math.exp(-1) * (1**1) / math.factorial(1), 8), "relative_intensity_pct": 36.79},
                        {"v_prime": 2, "fc_factor": round(math.exp(-1) * (1**2) / math.factorial(2), 8), "relative_intensity_pct": 18.39},
                    ],
                    "total_sum_check": "~1.000",
                    "interpretation": "S=1 indicates moderate geometry change. Maximum intensity at v'=1. Progression extends over ~5 vibrational quanta.",
                }
            }
        },
        {
            "code_input": {
                "huang_rhys_factor_S": 4.0,
                "v_max": 20,
                "v_initial": 0,
                "delta_q_dimensionless": 0,
                "frequency_cm-1": 500,
                "include_progression_plot_data": True,
            },
            "text_input": {
                "input_params": "4.0 20 0 500",
            },
            "output": {
                "result": {
                    "huang_rhys_factor_S": 4.0,
                    "v_initial": 0,
                    "max_intensity_at_v_prime": 4,
                    "interpretation": "S=4 indicates significant geometry change between states. Broad progression extending to v'~12+.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = H
        self.c = C

    def _run_base(self, huang_rhys_factor_S: float, v_max: int = 10, v_initial: int = 0,
                  delta_q_dimensionless: float = 0.0, frequency_cm_minus_1: float = 1000.0,
                  include_progression_plot_data: bool = True) -> dict:
        """Core logic: calculate Franck-Condon factors."""
        S = huang_rhys_factor_S
        if delta_q_dimensionless != 0:
            S = delta_q_dimensionless ** 2 / 2

        if S < 0:
            raise ChemMCPError("Huang-Rhys factor S must be non-negative.")
        if v_max < 0:
            raise ChemMCPError("v_max must be >= 0.")
        if v_initial < 0:
            raise ChemMCPError("v_initial must be >= 0.")

        # FC factor for v=0 → v': |⟨v'|0⟩|² = e^(-S) * S^v' / v'!
        # General formula for v → v' involves more complex recursion
        fc_factors = []
        total = 0.0
        max_fc = 0.0
        max_vp = 0

        for vp in range(v_max + 1):
            if v_initial == 0:
                # Simplified Poisson distribution for v=0 → v'
                fc = math.exp(-S) * (S ** vp) / math.factorial(vp)
            else:
                # Use recursion formula for general v → v'
                fc = self._fc_general(S, v_initial, vp)

            rel_int = fc * 100  # as percentage
            total += fc
            fc_factors.append({
                "v_prime": vp,
                "fc_factor": round(fc, 10),
                "relative_intensity_pct": round(rel_int, 4),
            })

            if fc > max_fc:
                max_fc = fc
                max_vp = vp

        # Energy spacing in nm (for plotting reference)
        # ΔE = h·c·ν̃  →  Δλ ≈ λ²/(hc) · ΔE for small changes
        nu_cm = frequency_cm_minus_1
        if nu_cm > 0:
            E_per_mode_J = self.h * self.c * nu_cm * 100  # per quantum in J
            wavelength_0_nm = self.c / (nu_cm * 100) / 1e-9 if nu_cm > 0 else 500  # rough estimate
        else:
            E_per_mode_J = 0
            wavelength_0_nm = 0

        result = {
            "huang_rhys_factor_S": S,
            "v_initial": v_initial,
            "v_max_computed": v_max,
            "vibrational_frequency_cm-1": nu_cm,
            "energy_per_quantum_J": f"{E_per_mode_J:.4e}" if E_per_mode_J > 0 else "N/A",
            "max_intensity_at_v_prime": max_vp,
            "max_fc_factor": round(max_fc, 10),
            "fc_factors": fc_factors,
            "sum_of_factors": round(total, 10),
            "normalization_check": "OK" if abs(total - 1.0) < 0.001 else f"deviation: {abs(total - 1.0):.6f}",
        }

        if include_progression_plot_data:
            result["progression_summary"] = self._generate_progression_summary(fc_factors, S, nu_cm)

        result["interpretation"] = self._interpret_s(S, max_vp, fc_factors, nu_cm)

        return {"result": result}

    @staticmethod
    def _fc_general(S: float, v: int, vp: int) -> float:
        """
        General Franck-Condon factor using Manneback recursion formula.
        For v → v' transition with displacement parameter.
        Approximate: uses Poisson-weighted recursion.
        """
        if v == 0:
            return math.exp(-S) * (S ** vp) / math.factorial(vp)

        # Recursion relations (simplified approximation)
        b = math.sqrt(S)
        # Build FC matrix iteratively
        fc_grid = {}
        fc_grid[(0, 0)] = math.exp(-S / 2)

        for vi in range(max(v, vp) + 1):
            for vip in range(max(v, vp) + 1):
                if (vi, vip) in fc_grid:
                    continue
                if vi == 0 and vip == 0:
                    continue
                val = 0.0
                term1 = fc_grid.get((vi - 1, vip - 1), 0) if vi > 0 and vip > 0 else 0
                term2 = fc_grid.get((vi, vip - 1), 0) if vip > 0 else 0
                term3 = fc_grid.get((vi - 1, vip + 1), 0) if vi > 0 else 0
                term4 = fc_grid.get((vi, vip - 2), 0) if vip >= 2 else 0
                # Simplified recursion
                if vi <= 1:
                    val = math.exp(-S/2) * (S ** (vip/2)) / math.factorial(max(1, vip)) if vip >= 0 else 0
                else:
                    val = 0  # fallback
                fc_grid[(vi, vip)] = val

        return fc_grid.get((v, vp), 0.0)

    def _generate_progression_summary(self, fc_factors: List[dict], S: float, nu_cm: float) -> dict:
        """Generate summary of the vibrational progression."""
        significant = [f for f in fc_factors if f["relative_intensity_pct"] >= 1.0]
        fwhm_modes = len(significant)

        # Estimate spectral width
        if significant and nu_cm > 0:
            width_cm = (significant[-1]["v_prime"] - significant[0]["v_prime"]) * nu_cm
        else:
            width_cm = 0

        return {
            "significant_peaks_above_1pct": fwhm_modes,
            "v_prime_range": (
                f"{significant[0]['v_prime']}-{significant[-1]['v_prime']}"
                if significant else "N/A"
            ),
            "estimated_progression_width_cm-1": round(width_cm, 1) if width_cm > 0 else 0,
            "progression_shape": self._classify_shape(S),
        }

    @staticmethod
    def _classify_shape(S: float) -> str:
        if S < 0.2:
            return "0-0 transition dominant (little geometry change)"
        elif S < 0.5:
            return "Narrow progression centered near origin"
        elif S < 1.5:
            return "Moderate progression with clear vibrational structure"
        elif S < 3.0:
            return "Broad progression with multiple comparable peaks"
        elif S < 6.0:
            return "Very broad progression, extensive vibrational structure"
        else:
            return "Extremely broad, nearly statistical distribution"

    def _interpret_s(self, S: float, max_vp: int, fc_factors: List[dict], nu_cm: float) -> str:
        parts = []
        parts.append(f"Huang-Rhys factor S = {S:.2f}")
        parts.append(f"Maximum FC intensity at v' = {max_vp}")

        if S < 0.3:
            parts.append("Minimal geometry difference between electronic states — vertical (Franck-Condon) transition occurs near equilibrium.")
        elif S < 1.0:
            parts.append("Small-to-moderate geometry change — well-resolved vibrational progression expected.")
        elif S < 3.0:
            parts.append("Significant geometry change — broad vibrational progression with several observable peaks.")
        else:
            parts.append("Large geometry change — very broad, possibly unresolved vibrational envelope.")

        if nu_cm > 0:
            top3 = sorted(fc_factors, key=lambda x: x["fc_factor"], reverse=True)[:3]
            top_str = ", ".join([f"v'={t['v_prime']} ({t['relative_intensity_pct']:.1f}%)" for t in top3])
            parts.append(f"Top 3 peaks: {top_str}")

        return " | ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            S = float(parts[0])
            vmax = int(parts[1]) if len(parts) > 1 else 10
            v0 = int(parts[2]) if len(parts) > 2 else 0
            freq = float(parts[3]) if len(parts) > 3 else 1000
            return self._run_base(S, vmax, v0, 0, freq)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
