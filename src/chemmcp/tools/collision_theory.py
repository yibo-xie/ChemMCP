import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CollisionTheory(BaseTool):
    """
    碰撞理论计算反应速率常数。
    基于 Arrhenius 方程结合碰撞理论：k = A * f * exp(-Ea/RT)
    其中 A 为指前因子（碰撞频率），f 为方位因子，Ea 为活化能。
    """
    __version__ = "0.1.0"
    name = "CollisionTheory"
    func_name = "collision_theory_rate"
    description = "Calculate reaction rate constant using collision theory (Arrhenius equation with steric factor)."
    implementation_description = "Uses the Arrhenius equation k = A * f * exp(-Ea/RT) where A is pre-exponential factor, f is steric factor, Ea is activation energy, R is gas constant, T is temperature."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Collision Theory", "Rate Constant", "Arrhenius"]
    required_envs = []

    code_input_sig = [
        ("pre_exponential_factor_A", "float", "N/A", "Pre-exponential factor (frequency factor), units depend on reaction order."),
        ("activation_energy_Ea", "float", "N/A", "Activation energy in J/mol."),
        ("temperature_K", "float", "N/A", "Absolute temperature in Kelvin."),
        ("steric_factor_f", "float", "1.0", "Steric (probability) factor, between 0 and 1. Default is 1.0."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'A Ea T [f]' where A=pre-exponential factor, Ea=activation energy (J/mol), T=temperature (K), f=steric factor (optional, default 1.0)."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing rate_constant_k, exponential_term, collision_frequency_factor, and input parameters summary."),
    ]

    examples = [
        {
            "code_input": {
                "pre_exponential_factor_A": 1e10,
                "activation_energy_Ea": 50000.0,
                "temperature_K": 300.0,
                "steric_factor_f": 1.0,
            },
            "text_input": {
                "input_params": "1e10 50000.0 300.0 1.0",
            },
            "output": {
                "result": {
                    "rate_constant_k": 5759863694.5,
                    "exponential_term": 0.57598636945,
                    "collision_frequency_factor": 1e10,
                    "steric_factor": 1.0,
                }
            },
        },
        {
            "code_input": {
                "pre_exponential_factor_A": 5e9,
                "activation_energy_Ea": 75000.0,
                "temperature_K": 298.0,
                "steric_factor_f": 0.01,
            },
            "text_input": {
                "input_params": "5e9 75000.0 298.0 0.01",
            },
            "output": {
                "result": {
                    "rate_constant_k": 3262.8,
                    "exponential_term": 6.5257e-05,
                    "collision_frequency_factor": 5e9,
                    "steric_factor": 0.01,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314  # J/(mol·K)

    def _run_base(self, pre_exponential_factor_A: float, activation_energy_Ea: float,
                  temperature_K: float, steric_factor_f: float = 1.0) -> dict:
        if temperature_K <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")
        if activation_energy_Ea < 0:
            raise ChemMCPError("Activation energy cannot be negative.")
        if steric_factor_f < 0 or steric_factor_f > 1:
            raise ChemMCPError("Steric factor must be between 0 and 1.")

        exp_term = math.exp(-activation_energy_Ea / (self.R * temperature_K))
        k = pre_exponential_factor_A * steric_factor_f * exp_term

        logger.info(f"CollisionTheory: A={pre_exponential_factor_A}, Ea={activation_energy_Ea}J/mol, "
                     f"T={temperature_K}K, f={steric_factor_f} -> k={k:.6g}")

        return {
            "rate_constant_k": round(k, 6),
            "exponential_term": exp_term,
            "collision_frequency_factor": pre_exponential_factor_A,
            "steric_factor": steric_factor_f,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need at least 3 parameters: A Ea T")
            A = float(parts[0])
            Ea = float(parts[1])
            T = float(parts[2])
            f = float(parts[3]) if len(parts) > 3 else 1.0
            return self._run_base(A, Ea, T, f)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'A Ea T [f]'")
