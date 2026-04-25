import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LeChatelierPrediction(BaseTool):
    """
    预测 Le Chatelier 原理下的平衡移动方向。
    当系统受到扰动（浓度、压力、温度变化）时，预测平衡如何移动以抵消该扰动。
    """
    __version__                = "0.1.0"
    name                       = "LeChatelierPrediction"
    func_name                  = "le_chatelier_prediction"
    description                = "Predict equilibrium shift direction using Le Chatelier's principle under disturbances."
    implementation_description = "Analyzes concentration, pressure, temperature, and inert gas disturbances; predicts shift direction and reasoning based on Le Chatelier's principle."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Le Chatelier", "Equilibrium", "Physical Chemistry", "Prediction"]
    required_envs              = []

    code_input_sig             = [
        ("reaction_type",     "str",   "N/A",   "Reaction type: 'exothermic' or 'endothermic'."),
        ("delta_n_gas",       "float", "N/A",   "Change in moles of gas (Σproducts - Σreactants) for gaseous species."),
        ("disturbance",       "str",   "N/A",   "Type of disturbance: 'increase_conc', 'decrease_conc', 'increase_pressure', 'decrease_pressure', 'increase_temp', 'decrease_temp', 'add_inert_gas_constant_v', 'add_inert_gas_constant_p', 'catalyst'."),
        ("species_affected",  "str",   "None",  "Name of the species affected (required for concentration changes)."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Format: 'reaction_type delta_n_gas disturbance [species_affected]'."),
    ]

    output_sig                 = [
        ("shift_direction",   "str",   "'forward' (→ products), 'backward' (→ reactants), or 'no_shift'."),
        ("reasoning",         "str",   "Detailed explanation of the prediction."),
        ("K_effect",          "str",   "'K increases', 'K decreases', or 'K unchanged'."),
    ]

    examples                   = [
        {
            "code_input": {
                "reaction_type": "exothermic",
                "delta_n_gas": -1.0,
                "disturbance": "increase_temp",
                "species_affected": None,
            },
            "text_input": {
                "input_params": "exothermic -1 increase_temp",
            },
            "output": {
                "shift_direction": "backward",
                "reasoning": "For an exothermic reaction, heat is treated as a product. Increasing temperature adds heat → system shifts to consume it (toward reactants).",
                "K_effect": "K decreases",
            }
        },
        {
            "code_input": {
                "reaction_type": "endothermic",
                "delta_n_gas": 2.0,
                "disturbance": "increase_pressure",
                "species_affected": None,
            },
            "text_input": {
                "input_params": "endothermic 2 increase_pressure",
            },
            "output": {
                "shift_direction": "backward",
                "reasoning": "Δn_gas = +2 (more moles of gas on product side). Increasing pressure favors the side with fewer moles of gas (reactants).",
                "K_effect": "K unchanged",
            }
        },
        {
            "code_input": {
                "reaction_type": "exothermic",
                "delta_n_gas": 0.0,
                "disturbance": "increase_conc",
                "species_affected": "N2",
            },
            "text_input": {
                "input_params": "exothermic 0 increase_conc N2",
            },
            "output": {
                "shift_direction": "forward",
                "reasoning": "Increasing concentration of reactant N2 shifts equilibrium toward products to consume the added N2.",
                "K_effect": "K unchanged",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, delta_n_gas: float, disturbance: str,
                  species_affected: str = None) -> dict:
        reaction_type = reaction_type.lower()
        disturbance = disturbance.lower()

        if reaction_type not in ("exothermic", "endothermic"):
            raise ChemMCPError("reaction_type must be 'exothermic' or 'endothermic'.")

        shift = "no_shift"
        reasoning = ""
        K_effect = "K unchanged"

        d = disturbance
        if d in ("increase_conc", "decrease_conc"):
            if not species_affected:
                raise ChemMCPError("species_affected is required for concentration disturbances.")
            K_effect = "K unchanged"
            if d == "increase_conc":
                shift = "forward"
                reasoning = (
                    f"Increasing concentration of {species_affected}. "
                    f"The system shifts to consume the added {species_affected}, "
                    f"moving equilibrium toward products."
                )
            else:
                shift = "backward"
                reasoning = (
                    f"Decreasing concentration of {species_affected}. "
                    f"The system shifts to produce more {species_affected}, "
                    f"moving equilibrium toward reactants."
                )

        elif d == "increase_pressure" or d == "decrease_pressure":
            K_effect = "K unchanged"
            if abs(delta_n_gas) < 1e-6:
                shift = "no_shift"
                reasoning = (
                    f"Δn_gas = 0: equal moles of gas on both sides. "
                    f"Pressure change has no effect on equilibrium position."
                )
            elif d == "increase_pressure":
                if delta_n_gas > 0:
                    shift = "backward"
                    side_fewer = "reactants"
                else:
                    shift = "forward"
                    side_fewer = "products"
                reasoning = (
                    f"Δn_gas = {delta_n_gas:+.1f} ({'more' if delta_n_gas > 0 else 'fewer'} gas moles on product side). "
                    f"Increasing pressure favors the side with fewer gas moles ({side_fewer})."
                )
            else:
                if delta_n_gas > 0:
                    shift = "forward"
                    side_more = "products"
                else:
                    shift = "backward"
                    side_more = "reactants"
                reasoning = (
                    f"Δn_gas = {delta_n_gas:+.1f} ({'more' if delta_n_gas > 0 else 'fewer'} gas moles on product side). "
                    f"Decreasing pressure favors the side with more gas moles ({side_more})."
                )

        elif d in ("increase_temp", "decrease_temp"):
            if reaction_type == "exothermic":
                # Heat is a product
                if d == "increase_temp":
                    shift = "backward"
                    K_effect = "K decreases"
                    reasoning = (
                        "Exothermic reaction: heat is a product. "
                        "Adding heat (↑T) shifts equilibrium toward reactants to consume it. K decreases."
                    )
                else:
                    shift = "forward"
                    K_effect = "K increases"
                    reasoning = (
                        "Exothermic reaction: heat is a product. "
                        "Removing heat (↓T) shifts equilibrium toward products to produce more heat. K increases."
                    )
            else:  # endothermic
                # Heat is a reactant
                if d == "increase_temp":
                    shift = "forward"
                    K_effect = "K increases"
                    reasoning = (
                        "Endothermic reaction: heat is a reactant. "
                        "Adding heat (↑T) shifts equilibrium toward products to consume it. K increases."
                    )
                else:
                    shift = "backward"
                    K_effect = "K decreases"
                    reasoning = (
                        "Endothermic reaction: heat is a reactant. "
                        "Removing heat (↓T) shifts equilibrium toward reactants to produce more heat. K decreases."
                    )

        elif d in ("add_inert_gas_constant_v", "add_inert_gas_constant_p"):
            K_effect = "K unchanged"
            if d == "add_inert_gas_constant_v":
                shift = "no_shift"
                reasoning = (
                    "Adding inert gas at constant volume does not change partial pressures of reacting species. "
                    "No shift occurs."
                )
            else:
                shift = "no_shift"
                reasoning = (
                    "Adding inert gas at constant pressure dilutes the system but Q=K remains valid at new partial pressures. "
                    "No shift occurs (equivalent to reducing all partial pressures proportionally)."
                )

        elif d == "catalyst":
            K_effect = "K unchanged"
            shift = "no_shift"
            reasoning = (
                "A catalyst speeds up both forward and reverse reactions equally. "
                "It reduces time to reach equilibrium but does not change the equilibrium position or K value."
            )

        else:
            raise ChemMCPError(f"Unknown disturbance type: {disturbance}")

        logger.info(f"Le Chatelier: disturbance={d}, shift={shift}")
        return {
            "shift_direction": shift,
            "reasoning": reasoning,
            "K_effect": K_effect,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need: reaction_type delta_n_gas disturbance [species]")
            rtype = parts[0]
            dn = float(parts[1])
            dist = parts[2]
            species = parts[3] if len(parts) > 3 else None
            return self._run_base(rtype, dn, dist, species)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
