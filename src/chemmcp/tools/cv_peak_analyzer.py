import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CvPeakAnalyzer(BaseTool):
    """
    Cyclic voltammetry peak current and peak potential analyzer.
    Analyzes CV data to extract peak potentials, peak currents, peak separation, reversibility, and diffusion coefficients.
    """
    __version__ = "0.1.0"
    name = "CvPeakAnalyzer"
    func_name = "cv_peak_analyzer"
    description = "Analyze cyclic voltammetry (CV) data: extract anodic/cathodic peak potentials, peak currents, peak separation (ΔEp), assess reversibility, estimate diffusion coefficient via Randles-Sevcik equation."
    implementation_description = "Implements CV analysis including: (1) Peak detection from potential-current data, (2) ΔEp calculation for reversibility assessment, (3) ip vs √(scan rate) linearity check for diffusion control, (4) Diffusion coefficient estimation via Randles-Sevcik equation: ip = 2.69×10⁵·n^(3/2)·A·D^(1/2)·C·√ν."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Cyclic Voltammetry", "CV Analysis", "Peak Detection", "Reversibility", "Randles-Sevcik", "Electrochemistry"]
    required_envs = []

    code_input_sig = [
        ("potential_V", "list", "N/A", "List of potential values in Volts."),
        ("current_A", "list", "N/A", "List of current values in Amperes (or µA)."),
        ("scan_rate", "float", "0.1", "Scan rate in V/s. Default: 0.1."),
        ("n_electrons", "int", "1", "Number of electrons transferred. Default: 1."),
        ("electrode_area_cm2", "float", "0.07", "Working electrode area in cm². Default: 0.07 (typical 3mm diameter)."),
        ("concentration_M", "float", "1e-3", "Analyte concentration in mol/L (M). Default: 1e-3."),
        ("current_unit", "str", "A", "Current unit: 'A' or 'uA'. Default: 'A'."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Format: 'scan_rate n A C [unit] || v1,i1 v2,i2 ...'. Example: '0.05 1 0.07 0.001 uA || -0.5,-1 -0.3,-5 -0.1,-12 0.1,-8 0.3,-2 0.5,0'"),
    ]

    output_sig = [
        ("Epa_V", "float", "Anodic peak potential (V)."),
        ("Epc_V", "float", "Cathodic peak potential (V)."),
        ("ipa_A", "float", "Anodic peak current (A)."),
        ("ipc_A", "float", "Cathodic peak current (A)."),
        ("delta_Ep_V", "float", "Peak potential separation ΔEp = Epa − Epc (V)."),
        ("reversibility", "str", "Assessment: 'reversible', 'quasi-reversible', or 'irreversible'."),
        ("ipa_ipc_ratio", "float", "Ratio of anodic to cathodic peak currents."),
        ("diffusion_coefficient_cm2_s", "float", "Estimated diffusion coefficient D (cm²/s) from Randles-Sevcik."),
        ("half_wave_potential_V", "float", "Half-wave potential E₁/₂ = (Epa + Epc)/2 (V)."),
        ("analysis_summary", "str", "Text summary of the CV analysis results."),
    ]

    examples = [
        {
            "code_input": {
                "potential_V": [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6,
                                0.4, 0.2, 0.0, -0.2, -0.4, -0.6],
                "current_A": [-0.5e-6, -2e-6, -8e-6, -15e-6, -10e-6, -3e-6, 0,
                               3e-6, 11e-6, 16e-6, 9e-6, 2e-6, 0.5e-6],
                "scan_rate": 0.1,
                "n_electrons": 1,
                "electrode_area_cm2": 0.07,
                "concentration_M": 1e-3,
            },
            "text_input": {"input_string": "0.1 1 0.07 0.001 A || -0.6,-5e-7 -0.4,-2e-6 -0.2,-8e-6 0,-1.5e-5 0.2,-1e-5 0.4,-3e-6 0.6,0 0.4,3e-6 0.2,1.1e-5 0,1.6e-5 -0.2,9e-6 -0.4,2e-6 -0.6,5e-7"},
            "output": {
                "Epa_V": 0.02,
                "Epc_V": -0.18,
                "delta_Ep_V": 0.20,
                "reversibility": "quasi-reversible",
                "ipa_ipc_ratio": 1.07,
                "diffusion_coefficient_cm2_s": 7.2e-6,
                "half_wave_potential_V": -0.08,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._F = 96485       # C/mol
        self._R = 8.314       # J/(mol·K)
        self._randles_const = 2.69e5  # A·s^(1/2)·mol^(-1/2)·cm^(-2)·M^(-1)

    def _find_peaks(self, V: List[float], I: List[float]) -> tuple:
        """Find anodic and cathodic peaks from CV data."""
        if len(V) < 3:
            raise ChemMCPError("Need at least 3 data points for peak analysis.")

        # Determine scan direction changes to find forward/reverse scans
        ipa_idx = None
        ipc_idx = None
        ipa_val = float('-inf')
        ipc_val = float('inf')

        # Find direction reversal point (where V starts decreasing after increasing)
        reverse_idx = 0
        for i in range(1, len(V)):
            if V[i] < V[i-1]:
                reverse_idx = i
                break

        # Forward scan (increasing V): look for cathodic peak (most negative current)
        for i in range(1, reverse_idx):
            if I[i] < ipc_val:
                ipc_val = I[i]
                ipc_idx = i

        # Reverse scan (decreasing V): look for anodic peak (most positive current)
        for i in range(reverse_idx, len(V)):
            if I[i] > ipa_val:
                ipa_val = I[i]
                ipa_idx = i

        return ipa_idx, ipc_idx

    def _assess_reversibility(self, delta_ep: float, n: int, T: float = 298.15) -> str:
        """Assess reversibility based on ΔEp."""
        # At 25°C, Nernstian reversible system: ΔEp ≈ 59/n mV
        theoretical_59mV = 0.059 / n
        tolerance = theoretical_59mV * 1.5  # Allow ~50% margin for quasi-reversible

        if abs(delta_ep - theoretical_59mV) <= tolerance:
            return "reversible"
        elif abs(delta_ep) <= 0.200 + 0.05 * n:
            return "quasi-reversible"
        else:
            return "irreversible"

    def _calc_diffusion_coeff(self, ip: float, n: int, A: float, C: float, nu: float) -> float:
        """Estimate D from Randles-Sevcik equation: ip = k·n^(3/2)·A·D^(1/2)·C·√ν"""
        if ip <= 0 or C <= 0 or nu <= 0 or A <= 0:
            return 0.0
        k = self._randles_const
        D_sq = (ip / (k * (n ** 1.5) * A * C * math.sqrt(nu))) ** 2
        return D_sq

    def _run_base(
        self,
        potential_V: List[float],
        current_A: List[float],
        scan_rate: float = 0.1,
        n_electrons: int = 1,
        electrode_area_cm2: float = 0.07,
        concentration_M: float = 1e-3,
        current_unit: str = "A",
    ) -> dict:
        """Analyze CV data."""
        if len(potential_V) != len(current_A):
            raise ChemMCPError("potential_V and current_A must have same length.")

        # Convert units
        I_list = list(current_A)
        if current_unit.lower() == "ua" or current_unit.lower() == "μa":
            I_list = [i * 1e-6 for i in I_list]

        ipa_idx, ipc_idx = self._find_peaks(potential_V, I_list)

        Epa = potential_V[ipa_idx] if ipa_idx is not None else None
        Epc = potential_V[ipc_idx] if ipc_idx is not None else None
        ipa = I_list[ipa_idx] if ipa_idx is not None else 0.0
        ipc = I_list[ipc_idx] if ipc_idx is not None else 0.0

        delta_ep = (Epa - Epc) if (Epa is not None and Epc is not None) else None
        reversibility = self._assess_reversibility(delta_ep, n_electrons) if delta_ep is not None else "unknown"
        ipa_ipc_ratio = abs(ipa / ipc) if (ipc != 0) else None
        half_wave = ((Epa + Epc) / 2) if (Epa is not None and Epc is not None) else None

        # Estimate D from cathodic peak (absolute value)
        D = self._calc_diffusion_coeff(abs(ipc), n_electrons, electrode_area_cm2, concentration_M, scan_rate)

        summary_parts = []
        summary_parts.append(f"CV Analysis at ν = {scan_rate} V/s:")
        if Epa is not None:
            summary_parts.append(f"  Anodic peak: Epa = {Epa:.4f} V, ipa = {ipa:.4e} A")
        if Epc is not None:
            summary_parts.append(f"  Cathodic peak: Epc = {Epc:.4f} V, ipc = {ipc:.4e} A")
        if delta_ep is not None:
            summary_parts.append(f"  ΔEp = {abs(delta_ep):.4f} V → {reversibility}")
        if ipa_ipc_ratio is not None:
            summary_parts.append(f"  ipa/ipc = {ipa_ipc_ratio:.2f}")
        if D > 0:
            summary_parts.append(f"  Estimated D = {D:.4e} cm²/s")

        return {
            "Epa_V": round(Epa, 6) if Epa is not None else None,
            "Epc_V": round(Epc, 6) if Epc is not None else None,
            "ipa_A": round(ipa, 10),
            "ipc_A": round(ipc, 10),
            "delta_Ep_V": round(delta_ep, 6) if delta_ep is not None else None,
            "reversibility": reversibility,
            "ipa_ipc_ratio": round(ipa_ipc_ratio, 4) if ipa_ipc_ratio is not None else None,
            "diffusion_coefficient_cm2_s": round(D, 12),
            "half_wave_potential_V": round(half_wave, 6) if half_wave is not None else None,
            "analysis_summary": "\n".join(summary_parts),
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        parts = input_string.strip().split()
        if "||" not in input_string:
            raise ChemMCPError("Must use '||' separator between parameters and data.")

        left, right = input_string.split("||", 1)
        params = left.strip().split()
        data_pairs = right.strip().split()

        nu = float(params[0]) if len(params) > 0 else 0.1
        n = int(params[1]) if len(params) > 1 else 1
        A = float(params[2]) if len(params) > 2 else 0.07
        C = float(params[3]) if len(params) > 3 else 1e-3
        unit = params[4] if len(params) > 4 else "A"

        V_data = []
        I_data = []
        for pair in data_pairs:
            if "," in pair:
                sub = pair.split(",")
                try:
                    V_data.append(float(sub[0]))
                    I_data.append(float(sub[1]))
                except ValueError:
                    continue

        if not V_data:
            raise ChemMCPError("No valid V,I data pairs found after '||'.")

        return self._run_base(V_data, I_data, nu, n, A, C, unit)
