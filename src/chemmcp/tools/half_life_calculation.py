import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HalfLifeCalculation(BaseTool):
    __version__ = "0.1.0"
    name = "HalfLifeCalculation"
    func_name = "calculate_half_life"
    description = "Calculate radioactive decay parameters: remaining amount, decay constant, elapsed time, or initial amount using half-life equations."
    implementation_description = "Uses the radioactive decay law N(t) = N0*e^(-lambda*t) where lambda = ln(2)/t_half."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Nuclear Chemistry", "Radioactive Decay", "Half-Life", "Kinetics"]
    required_envs = []

    code_input_sig = [
        ("calc_type", "str", "N/A", "Type of calculation."),
        ("half_life", "float", "N/A", "Half-life of the nuclide."),
        ("initial_amount", "float", "N/A", "Initial quantity."),
        ("time_elapsed", "float", "N/A", "Time elapsed."),
        ("remaining_amount", "float", "0", "Remaining quantity after decay."),
        ("decay_constant", "float", "0", "Decay constant lambda."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: calc_type half_life initial_amount time_elapsed [remaining_amount] [decay_constant]."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing all calculated parameters."),
    ]

    examples = [
        {
            "code_input": {
                "calc_type": "remaining_amount",
                "half_life": 5730.0,
                "initial_amount": 100.0,
                "time_elapsed": 17190.0,
                "remaining_amount": 0.0,
                "decay_constant": 0.0,
            },
            "text_input": {"input_params": "remaining_amount 5730 100 17190"},
            "output": {
                "result": {
                    "calc_type": "remaining_amount",
                    "remaining_amount": 12.5,
                    "half_lives_elapsed": 3.0,
                }
            }
        },
        {
            "code_input": {
                "calc_type": "decay_constant",
                "half_life": 4.468e9,
                "initial_amount": 1.0,
                "time_elapsed": 0.0,
                "remaining_amount": 0.0,
                "decay_constant": 0.0,
            },
            "text_input": {"input_params": "decay_constant 4.468e9 1 0"},
            "output": {
                "result": {
                    "calc_type": "decay_constant",
                    "decay_constant": 1.551e-10,
                }
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.LN2 = math.log(2)

    def _run_base(self, calc_type, half_life, initial_amount, time_elapsed,
                  remaining_amount=0.0, decay_constant=0.0):
        if half_life <= 0:
            raise ChemMCPError("Half-life must be positive.")
        lam = self.LN2 / half_life
        n_halves = time_elapsed / half_life if time_elapsed > 0 else None

        result = {
            "calc_type": calc_type,
            "half_life": half_life,
            "initial_amount": initial_amount,
            "time_elapsed": time_elapsed,
            "decay_constant": round(lam, 10),
        }

        if calc_type == "remaining_amount":
            if initial_amount <= 0:
                raise ChemMCPError("initial_amount must be positive.")
            rem = initial_amount * math.exp(-lam * time_elapsed)
            result["remaining_amount"] = round(rem, 6)
            result["half_lives_elapsed"] = round(n_halves, 4) if n_halves is not None else None
            result["formula_used"] = "N(t) = N0 * e^(-lambda*t)"
        elif calc_type == "decay_constant":
            result["formula_used"] = "lambda = ln(2) / t_1/2"
        elif calc_type == "elapsed_time":
            if remaining_amount <= 0:
                raise ChemMCPError("remaining_amount must be positive.")
            if remaining_amount >= initial_amount:
                raise ChemMCPError("remaining_amount must be less than initial_amount.")
            t = -math.log(remaining_amount / initial_amount) / lam
            result["elapsed_time"] = round(t, 6)
            result["half_lives_elapsed"] = round(t / half_life, 4)
            result["formula_used"] = "t = -ln(N/N0) / lambda"
        elif calc_type == "initial_amount":
            if remaining_amount <= 0:
                raise ChemMCPError("remaining_amount must be positive.")
            n0 = remaining_amount / math.exp(-lam * time_elapsed)
            result["initial_amount"] = round(n0, 6)
            result["formula_used"] = "N0 = N(t) / e^(-lambda*t)"
        elif calc_type == "half_life_from_decay":
            if decay_constant <= 0:
                raise ChemMCPError("decay_constant must be positive.")
            t_half = self.LN2 / decay_constant
            result["calculated_half_life"] = round(t_half, 6)
            result["formula_used"] = "t_1/2 = ln(2) / lambda"
        else:
            raise ChemMCPError(f"Unknown calc_type '{calc_type}'.")
        logger.info(f"Half-life calculation ({calc_type}): {result}")
        return result

    def _run_text(self, input_params):
        parts = input_params.split()
        if len(parts) < 4:
            raise ValueError("Need at least 4 params.")
        kw = {
            "calc_type": parts[0],
            "half_life": float(parts[1]),
            "initial_amount": float(parts[2]),
            "time_elapsed": float(parts[3]),
        }
        if len(parts) > 4:
            kw["remaining_amount"] = float(parts[4])
        if len(parts) > 5:
            kw["decay_constant"] = float(parts[5])
        return self._run_base(**kw)
