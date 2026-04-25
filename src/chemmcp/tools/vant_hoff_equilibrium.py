import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

R = 8.314  # J/(mol·K)


@ChemMCPManager.register_tool
class VantHoffEquilibrium(BaseTool):
    """
    范特霍夫方程：计算平衡常数的温度依赖性。
    
    ln(K2/K1) = -ΔH°/R × (1/T2 - 1/T1)
    """
    __version__ = "0.1.0"
    name = "VantHoffEquilibrium"
    func_name = "vant_hoff_calc"
    description = "Calculate equilibrium constant at a different temperature using the Van't Hoff equation (temperature dependence of K)."
    implementation_description = "Uses ln(K2/K1) = -ΔH°/R × (1/T2 - 1/T1). Assumes ΔH° is constant over the temperature range. R = 8.314 J/(mol·K)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Equilibrium", "Van't Hoff", "Temperature Dependence"]
    required_envs = []

    code_input_sig = [
        ("k1", "float", "N/A", "Equilibrium constant at temperature T1."),
        ("t1", "float", "N/A", "Initial temperature in Kelvin."),
        ("t2", "float", "N/A", "Target temperature in Kelvin."),
        ("delta_h", "float", "N/A", "Standard enthalpy change ΔH° in J/mol (assumed constant over T range)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'k1 t1 t2 delta_h'. Example: '0.01 300 350 -50000'"),
    ]

    output_sig = [
        ("k2", "float", "Equilibrium constant at temperature T2."),
        ("ln_ratio", "float", "Natural log of K2/K1."),
        ("explanation", "str", "Step-by-step calculation with formula."),
    ]

    examples = [
        {
            "code_input": {
                "k1": 6.7e-2,
                "t1": 298.15,
                "t2": 350.0,
                "delta_h": -57000.0,
            },
            "text_input": {
                "input_params": "0.067 298.15 350 -57000",
            },
            "output": {
                "k2": 0.5982,
                "ln_ratio": 2.1894,
                "explanation": "ln(K2/0.067) = -(-57000)/8.314×(1/350-1/298.15) → K2 = 0.067×exp(2.189) = 0.598",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, k1: float, t1: float, t2: float, delta_h: float) -> dict:
        """Core logic: Van't Hoff equation."""
        if k1 <= 0:
            raise ChemMCPError("K1 must be positive.")
        if t1 <= 0 or t2 <= 0:
            raise ChemMCPError("Temperatures must be positive in Kelvin.")
        if t1 == t2:
            raise ChemMCPError("T1 and T2 must be different.")

        ln_ratio = -(delta_h / R) * (1.0 / t2 - 1.0 / t1)
        k2 = k1 * math.exp(ln_ratio)

        explanation = (
            f"Van't Hoff Equation:\n"
            f"ln(K2/K1) = -ΔH°/R × (1/T2 - 1/T1)\n"
            f"= -({delta_h}/{R}) × (1/{t2} - 1/{t1})\n"
            f"= {ln_ratio:.4f}\n"
            f"K2 = K1 × exp(ln_ratio) = {k1} × {math.exp(ln_ratio):.4f}\n"
            f"= {k2:.6f}"
        )

        logger.info(f"Van't Hoff: K1={k1} @ {t1}K → K2={k2:.6f} @ {t2}K (ΔH={delta_h} J/mol)")
        return {"k2": round(k2, 6), "ln_ratio": round(ln_ratio, 6), "explanation": explanation}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 4:
                raise ValueError("Need 4 parameters: k1 t1 t2 delta_h")
            k1 = float(parts[0])
            t1 = float(parts[1])
            t2 = float(parts[2])
            delta_h = float(parts[3])
            return self._run_base(k1, t1, t2, delta_h)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'k1 t1 t2 delta_h'")
