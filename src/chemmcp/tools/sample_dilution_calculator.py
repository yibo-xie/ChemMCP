import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SampleDilutionCalculator(BaseTool):
    """
    样品稀释计算器：计算稀释比例、终浓度和所需溶剂体积。
    基于 C1V1 = C2V2 稀释公式。
    """
    __version__ = "0.1.0"
    name = "SampleDilutionCalculator"
    func_name = "calculate_dilution"
    description = "Calculate dilution ratio, final concentration, and required solvent volume for sample preparation."
    implementation_description = "Uses C1*V1 = C2*V2 dilution formula and serial dilution calculations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Sample Preparation", "Dilution", "Analytical Chemistry", "Laboratory"]
    required_envs = []

    code_input_sig = [
        ("initial_concentration", "float", "N/A", "Initial concentration of the stock solution (same units as desired final concentration)."),
        ("initial_volume", "float", "N/A", "Initial volume of the stock solution (mL)."),
        ("final_volume", "float", "N/A", "Desired final total volume after dilution (mL)."),
        ("dilution_steps", "int", "1", "Number of serial dilution steps (1 for simple dilution)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'initial_concentration initial_volume final_volume [dilution_steps]'."),
    ]

    output_sig = [
        ("dilution_factor", "float", "Total dilution factor (V_final / V_initial)."),
        ("solvent_volume_needed", "float", "Volume of solvent to add (mL)."),
        ("final_concentration", "float", "Final concentration after dilution."),
        ("dilution_ratio", "str", "Dilution ratio in format '1:N'."),
        ("step_details", "list", "Details for each dilution step (if serial dilution)."),
    ]

    examples = [
        {
            "code_input": {
                "initial_concentration": 1000.0,
                "initial_volume": 1.0,
                "final_volume": 100.0,
                "dilution_steps": 1,
            },
            "text_input": {
                "input_params": "1000.0 1.0 100.0 1",
            },
            "output": {
                "dilution_factor": 100.0,
                "solvent_volume_needed": 99.0,
                "final_concentration": 10.0,
                "dilution_ratio": "1:100",
                "step_details": [{"step": 1, "dilution_factor": 100.0, "solvent_added_ml": 99.0, "conc": 10.0}],
            },
        },
        {
            "code_input": {
                "initial_concentration": 500.0,
                "initial_volume": 2.0,
                "final_volume": 50.0,
                "dilution_steps": 1,
            },
            "text_input": {
                "input_params": "500.0 2.0 50.0",
            },
            "output": {
                "dilution_factor": 25.0,
                "solvent_volume_needed": 48.0,
                "final_concentration": 20.0,
                "dilution_ratio": "1:25",
                "step_details": [{"step": 1, "dilution_factor": 25.0, "solvent_added_ml": 48.0, "conc": 20.0}],
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        initial_concentration: float,
        initial_volume: float,
        final_volume: float,
        dilution_steps: int = 1,
    ) -> dict:
        """Core logic: calculate dilution parameters."""
        if initial_concentration <= 0:
            raise ChemMCPError("Initial concentration must be positive.")
        if initial_volume <= 0:
            raise ChemMCPError("Initial volume must be positive.")
        if final_volume <= initial_volume:
            raise ChemMCPError("Final volume must be greater than initial volume.")
        if dilution_steps < 1:
            raise ChemMCPError("Dilution steps must be >= 1.")

        total_dilution_factor = final_volume / initial_volume
        solvent_needed = final_volume - initial_volume
        final_conc = (initial_concentration * initial_volume) / final_volume

        # Simplify dilution ratio
        ratio = self._simplify_ratio(1, int(total_dilution_factor))

        # Serial dilution step details
        step_details = []
        if dilution_steps == 1:
            step_details.append({
                "step": 1,
                "dilution_factor": round(total_dilution_factor, 4),
                "solvent_added_ml": round(solvent_needed, 4),
                "conc": round(final_conc, 6),
            })
        else:
            # Equal dilution factor per step
            step_factor = math.pow(total_dilution_factor, 1.0 / dilution_steps)
            current_conc = initial_concentration
            current_vol = initial_volume
            for i in range(dilution_steps):
                step_final_vol = current_vol * step_factor
                solvent_step = step_final_vol - current_vol
                current_conc = (current_conc * current_vol) / step_final_vol
                step_details.append({
                    "step": i + 1,
                    "dilution_factor": round(step_factor, 4),
                    "solvent_added_ml": round(solvent_step, 4),
                    "conc": round(current_conc, 6),
                })
                current_vol = step_final_vol

        logger.info(f"Dilution calculated: factor={total_dilution_factor}, final_conc={final_conc}")
        return {
            "dilution_factor": round(total_dilution_factor, 4),
            "solvent_volume_needed": round(solvent_needed, 4),
            "final_concentration": round(final_conc, 6),
            "dilution_ratio": ratio,
            "step_details": step_details,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            c1 = float(parts[0])
            v1 = float(parts[1])
            v2 = float(parts[2])
            steps = int(parts[3]) if len(parts) > 3 else 1
            return self._run_base(c1, v1, v2, steps)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'C1 V1 V2 [steps]'")

    @staticmethod
    def _simplify_ratio(numerator: int, denominator: int) -> str:
        from math import gcd
        if denominator <= 0:
            return f"1:{numerator}"
        g = gcd(numerator, denominator)
        simplified_denom = denominator // g
        return f"1:{simplified_denom}"
