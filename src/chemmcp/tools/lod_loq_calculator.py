import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LodLoqCalculator(BaseTool):
    """
    检出限和定量限计算工具。
    支持 3σ/10σ 法（空白标准差法）和信噪比（S/N）法。
    """
    __version__ = "0.1.0"
    name = "LodLoqCalculator"
    func_name = "calculate_lod_loq"
    description = "Calculate Limit of Detection (LOD) and Limit of Quantification (LOQ) using 3σ/10σ method or S/N method."
    implementation_description = "Supports two methods: (1) Blank standard deviation method: LOD=3σ/slope, LOQ=10σ/slope; (2) Signal-to-noise method: LOD=3×noise/SN. Also supports calibration curve-based approach per ICH Q2(R1)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["LOD", "LOQ", "Detection Limit", "Analytical Chemistry", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("method", "str", "blank_std", "Method: 'blank_std' (3σ/10σ from blank measurements), 'calibration' (from calibration curve residual SD), 'sn' (signal-to-noise ratio)."),
        ("blank_data", "list", "", "Blank measurement values (for blank_std method)."),
        ("blank_std", "float", "", "Standard deviation of blank (alternative to blank_data)."),
        ("slope", "float", "", "Calibration curve slope (required for blank_std and calibration methods)."),
        ("signal_noise_ratio", "float", "", "S/N ratio of low-concentration sample (for sn method)."),
        ("low_sample_signal", "float", "", "Signal of a low-concentration sample near LOD (for sn method)."),
        ("n_blank", "int", "10", "Number of blank measurements (used for reporting)."),
        ("confidence_factor_lod", "float", "3.0", "Factor for LOD (default 3.0, IUPAC recommends 3.3 for calibration-based)."),
        ("confidence_factor_loq", "float", "10.0", "Factor for LOQ (default 10.0)."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("lod", "float", "Limit of Detection in concentration units."),
        ("loq", "float", "Limit of Quantification in concentration units."),
        ("method_used", "str", "Method used for calculation."),
        ("intermediate_values", "dict", "Intermediate calculation values (σ, slope, etc.)."),
        ("explanation", "str", "Step-by-step explanation."),
    ]

    examples = [
        {
            "code_input": {
                "method": "blank_std",
                "blank_data": [0.002, 0.003, 0.001, 0.004, 0.002, 0.003, 0.001, 0.002, 0.003, 0.002],
                "slope": 0.0992,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"lod": 0.091, "loq": 0.302},
        },
        {
            "code_input": {
                "method": "sn",
                "low_sample_signal": 0.15,
                "signal_noise_ratio": 3.5,
                "slope": 0.0992,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"lod": 0.129},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _std(data: List[float], ddof: int = 1) -> float:
        n = len(data)
        if n <= ddof:
            return 0.0
        m = sum(data) / n
        return math.sqrt(sum((x - m) ** 2 for x in data) / (n - ddof))

    @staticmethod
    def _mean(data: List[float]) -> float:
        return sum(data) / len(data) if data else 0.0

    def _run_base(
        self,
        method: str = "blank_std",
        blank_data: Optional[List[float]] = None,
        blank_std: Optional[float] = None,
        slope: Optional[float] = None,
        signal_noise_ratio: Optional[float] = None,
        low_sample_signal: Optional[float] = None,
        n_blank: int = 10,
        confidence_factor_lod: float = 3.0,
        confidence_factor_loq: float = 10.0,
    ) -> dict:
        """Core logic: calculate LOD and LOQ."""
        meth = method.lower().strip()
        intermediate: Dict[str, Any] = {}

        if meth == "blank_std":
            # Determine σ from data or provided std
            if blank_data is not None and len(blank_data) > 0:
                sigma = self._std(blank_data)
                n_blank = len(blank_data)
                intermediate["blank_mean"] = round(self._mean(blank_data), 8)
                intermediate["blank_n"] = n_blank
            elif blank_std is not None:
                sigma = blank_std
                intermediate["provided_blank_std"] = sigma
            else:
                raise ChemMCPError("For 'blank_std' method, provide either blank_data or blank_std.")

            if slope is None or slope == 0:
                raise ChemMCPError("Slope is required and must be non-zero for blank_std method.")

            intermediate["sigma"] = round(sigma, 8)
            intermediate["slope"] = slope

            lod = (confidence_factor_lod * sigma) / abs(slope)
            loq = (confidence_factor_loq * sigma) / abs(slope)

            explanation = (
                f"Method: Blank standard deviation (n={n_blank})\n"
                f"  σ_blank = {sigma:.6g}\n"
                f"  Slope = {slope:.6g}\n"
                f"  LOD = {confidence_factor_lod} × σ / |slope| = {confidence_factor_lod} × {sigma:.6g} / {abs(slope):.6g} = {lod:.4g}\n"
                f"  LOQ = {confidence_factor_loq} × σ / |slope| = {confidence_factor_loq} × {sigma:.6g} / {abs(slope):.6g} = {loq:.4g}"
            )

        elif meth == "calibration":
            # Calibration curve residual standard deviation approach (ICH Q2(R1))
            if slope is None or slope == 0:
                raise ChemMCPError("Slope is required for calibration method.")
            if blank_data is not None and len(blank_data) > 0:
                sigma = self._std(blank_data, ddof=1)
            elif blank_std is not None:
                sigma = blank_std
            else:
                raise ChemMCPError("For 'calibration' method, provide blank_data or blank_std as residual SD.")

            intermediate["residual_sd"] = round(sigma, 8)
            intermediate["slope"] = slope

            # ICH uses 3.3 factor for LOD from calibration
            lod = (confidence_factor_lod * sigma) / abs(slope)
            loq = (confidence_factor_loq * sigma) / abs(slope)

            explanation = (
                f"Method: Calibration curve residual SD\n"
                f"  Sy/x (residual SD) = {sigma:.6g}\n"
                f"  Slope = {slope:.6g}\n"
                f"  LOD = {confidence_factor_lod} × Sy/x / |slope| = {lod:.4g}\n"
                f"  LOQ = {confidence_factor_loq} × Sy/x / |slope| = {loq:.4g}"
            )

        elif meth == "sn":
            # Signal-to-Noise method
            if signal_noise_ratio is None or signal_noise_ratio <= 0:
                raise ChemMCPError("signal_noise_ratio must be positive for S/N method.")
            if low_sample_signal is None:
                raise ChemMCPError("low_sample_signal is required for S/N method.")
            if slope is None or slope == 0:
                raise ChemMCPError("Slope is required for S/N method.")

            sn = signal_noise_ratio
            sig = low_sample_signal
            noise = sig / sn  # estimated noise level

            intermediate["signal"] = sig
            intermediate["sn_ratio"] = sn
            intermediate["estimated_noise"] = round(noise, 8)

            # LOD concentration: signal at LOD = 3 × noise
            lod_conc_signal = 3.0 * noise
            lod = lod_conc_signal / abs(slope)
            loq_conc_signal = 10.0 * noise
            loq = loq_conc_signal / abs(slope)

            explanation = (
                f"Method: Signal-to-Noise\n"
                f"  Sample signal = {sig:.6g}, S/N = {sn:.2f}\n"
                f"  Estimated noise = signal/SN = {sig:.6g}/{sn:.2f} = {noise:.6g}\n"
                f"  LOD signal threshold = 3 × noise = {3*noise:.6g} → LOD = {lod:.4g}\n"
                f"  LOQ signal threshold = 10 × noise = {10*noise:.6g} → LOQ = {loq:.4g}"
            )

        else:
            raise ChemMCPError(f"Unknown method '{method}'. Use: 'blank_std', 'calibration', or 'sn'.")

        logger.info(f"LOD/LOQ ({meth}): LOD={lod:.4g}, LOQ={loq:.4g}")
        return {
            "lod": round(lod, 6),
            "loq": round(loq, 6),
            "method_used": meth,
            "intermediate_values": intermediate,
            "explanation": explanation,
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
