import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EquilibriumConstantThermo(BaseTool):
    """
    从热力学数据计算平衡常数。
    使用公式: ΔG° = -RT ln(K)  =>  K = exp(-ΔG° / RT)
    也支持从 ΔH° 和 ΔS° 计算: ΔG° = ΔH° - TΔS°, 然后 K = exp(-ΔG° / RT)
    """
    __version__                = "0.1.0"
    name                       = "EquilibriumConstantThermo"
    func_name                  = "equilibrium_constant_thermo"
    description                = "Calculate equilibrium constant K from thermodynamic data (ΔG°, or ΔH° & ΔS°)."
    implementation_description = "Uses ΔG° = -RT ln(K) to compute K. If ΔH° and ΔS° are provided instead, first computes ΔG° = ΔH° - TΔS°."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Thermodynamics", "Equilibrium", "Physical Chemistry"]
    required_envs              = []

    code_input_sig             = [
        ("temperature_k",     "float", "N/A",   "Temperature in Kelvin."),
        ("delta_g",           "float", "N/A",   "Standard Gibbs free energy change in kJ/mol. Provide this OR delta_h + delta_s."),
        ("delta_h",           "float", "None",  "Standard enthalpy change in kJ/mol (optional, used with delta_s)."),
        ("delta_s",           "float", "None",  "Standard entropy change in J/(mol·K) (optional, used with delta_h)."),
        ("r_gas",             "float", "8.314", "Gas constant in J/(mol·K). Default: 8.314."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Space-separated: 'temperature_k delta_g [delta_h delta_s r_gas]'."),
    ]

    output_sig                 = [
        ("K",                 "float", "Equilibrium constant K (dimensionless)."),
        ("delta_g_calculated","float", "Calculated ΔG° in kJ/mol (if derived from H and S)."),
        ("method",            "str",   "Method used: 'direct' (from ΔG°) or 'derived' (from ΔH° & ΔS°)."),
    ]

    examples                   = [
        {
            "code_input": {
                "temperature_k": 298.15,
                "delta_g": -23.7,
                "delta_h": None,
                "delta_s": None,
                "r_gas": 8.314,
            },
            "text_input": {
                "input_params": "298.15 -23.7",
            },
            "output": {
                "K": 13600.0,
                "delta_g_calculated": -23.7,
                "method": "direct",
            }
        },
        {
            "code_input": {
                "temperature_k": 298.15,
                "delta_g": None,
                "delta_h": -57.2,
                "delta_s": -112.5,
                "r_gas": 8.314,
            },
            "text_input": {
                "input_params": "298.15 None -57.2 -112.5",
            },
            "output": {
                "K": 6800.0,
                "delta_g_calculated": -23.65,
                "method": "derived",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, temperature_k: float, delta_g: float = None,
                  delta_h: float = None, delta_s: float = None,
                  r_gas: float = 8.314) -> dict:
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")

        method = "direct"

        if delta_g is not None:
            dg = delta_g  # kJ/mol
        elif delta_h is not None and delta_s is not None:
            # ΔG° = ΔH° - TΔS°; note: ΔH in kJ/mol, ΔS in J/(mol·K)
            dg = delta_h - temperature_k * delta_s / 1000.0
            method = "derived"
        else:
            raise ChemMCPError("Either delta_g, or both delta_h and delta_s must be provided.")

        # K = exp(-ΔG° / RT); ΔG° in kJ/mol → convert to J/mol
        exponent = -dg * 1000.0 / (r_gas * temperature_k)

        if exponent > 700:
            K = float('inf')
        elif exponent < -700:
            K = 0.0
        else:
            K = math.exp(exponent)

        logger.info(f"T={temperature_k}K, ΔG°={dg} kJ/mol -> K={K:.4e} ({method})")
        return {
            "K": round(K, 4),
            "delta_g_calculated": round(dg, 4),
            "method": method,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least temperature and delta_g.")
            t = float(parts[0])
            dg = None if parts[1].lower() == "none" else float(parts[1])
            dh = None
            ds = None
            rg = 8.314
            if len(parts) >= 4:
                dh = None if parts[2].lower() == "none" else float(parts[2])
                ds = None if parts[3].lower() == "none" else float(parts[3])
            if len(parts) >= 5:
                rg = float(parts[4])
            return self._run_base(t, dg, dh, ds, rg)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
