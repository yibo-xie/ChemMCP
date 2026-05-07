import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LeChatelierAnalyzer(BaseTool):
    """
    勒夏特列原理综合分析工具。
    
    当平衡系统受到扰动时，系统将调节以抵消该扰动。
    全面分析：浓度、压力/体积、温度、惰性气体、催化剂的影响，
    并给出平衡移动方向、K值变化及定量分析。
    """
    __version__ = "0.1.0"
    name = "LeChatelierAnalyzer"
    func_name = "analyze_le_chatelier"
    description = "Comprehensive Le Chatelier's principle analyzer: predict equilibrium shift direction, K value change, and quantitative analysis for concentration, pressure, temperature, inert gas, and catalyst disturbances."
    implementation_description = "Implements full Le Chatelier analysis: (1) Concentration changes shift Q vs K; (2) Pressure/volume effects via Δn_gas; (3) Temperature effects on exo/endothermic reactions; (4) Inert gas at constant V vs constant P; (5) Catalyst effect on rates only. Outputs shift direction, reasoning, K effect, and reaction quotient comparison."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Le Chatelier", "Equilibrium", "Physical Chemistry", "Chemical Equilibrium", "Prediction"]
    required_envs = []

    code_input_sig = [
        ("reaction_type", "str", "N/A", "Reaction thermicity: 'exothermic' (ΔH<0) or 'endothermic' (ΔH>0)."),
        ("delta_n_gas", "float", "N/A", "Change in moles of gas: Σν(gas, products) - Σν(gas, reactants)."),
        ("disturbance", "str", "N/A", "Disturbance type: 'increase_conc', 'decrease_conc', 'increase_pressure', 'decrease_pressure', 'increase_volume', 'decrease_volume', 'increase_temp', 'decrease_temp', 'add_inert_constant_V', 'add_inert_constant_P', 'add_catalyst'."),
        ("species", "str", "None", "Species name for concentration disturbances (e.g., 'NO2', 'N2O4')."),
        ("species_role", "str", "reactant", "Role of species: 'reactant' or 'product' (for conc disturbances)."),
        ("current_Q", "float", "None", "Current reaction quotient Q (optional, for quantitative analysis)."),
        ("K_value", "float", "None", "Equilibrium constant K value (optional, for Q vs K comparison)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'reaction_type delta_n disturbance [species] [role] [Q] [K]'. Example: 'exothermic -1 increase_temp' or 'endothermic 1 decrease_pressure N2 0.5 2.0'"),
    ]

    output_sig = [
        ("shift_direction", "str", "'forward'(→products), 'backward'(→reactants), or 'no_shift'."),
        ("reasoning", "str", "Detailed step-by-step explanation of the prediction."),
        ("K_effect", "str", "'K increases', 'K decreases', or 'K unchanged'."),
        ("Q_vs_K", "str", "Comparison of Q and K after disturbance (if Q/K provided)."),
        ("quantitative_note", "str", "Additional quantitative insights if applicable."),
    ]

    examples = [
        {
            "code_input": {
                "reaction_type": "exothermic",
                "delta_n_gas": -1.0,
                "disturbance": "increase_temp",
                "species": None,
                "species_role": "reactant",
                "current_Q": None,
                "K_value": None,
            },
            "text_input": {
                "input_params": "exothermic -1 increase_temp",
            },
            "output": {
                "shift_direction": "backward",
                "reasoning": "Exothermic reaction treats heat as a product. Increasing T → system shifts to absorb heat (toward reactants).",
                "K_effect": "K decreases",
            },
        },
        {
            "code_input": {
                "reaction_type": "endothermic",
                "delta_n_gas": 1.0,
                "disturbance": "increase_pressure",
                "species": None,
                "species_role": "reactant",
                "current_Q": 0.5,
                "K_value": 2.0,
            },
            "text_input": {
                "input_params": "endothermic 1 increase_pressure",
            },
            "output": {
                "shift_direction": "backward",
                "reasoning": "Δn_gas > 0: increasing pressure favors side with fewer gas moles (reactants).",
                "K_effect": "K unchanged",
                "Q_vs_K": "Q increases but remains < K initially; system shifts to reduce total gas moles.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, delta_n_gas: float, disturbance: str,
                  species: str = None, species_role: str = "reactant",
                  current_Q: float = None, K_value: float = None) -> dict:
        """Core logic: analyze Le Chatelier response."""
        rtype = reaction_type.lower().strip()
        dist = disturbance.lower().strip()

        # Validate inputs
        if rtype not in ("exothermic", "endothermic"):
            raise ChemMCPError("reaction_type must be 'exothermic' or 'endothermic'.")

        valid_disturbances = [
            "increase_conc", "decrease_conc",
            "increase_pressure", "decrease_pressure",
            "increase_volume", "decrease_volume",
            "increase_temp", "decrease_temp",
            "add_inert_constant_v", "add_inert_constant_p",
            "add_catalyst",
        ]
        if dist not in valid_disturbances:
            raise ChemMCPError(f"Unknown disturbance '{dist}'. Valid: {valid_disturbances}")

        shift = "no_shift"
        k_effect = "K unchanged"
        reasoning_parts = []
        q_vs_k = "N/A"
        quant_note = ""

        # --- Concentration Disturbances ---
        if dist in ("increase_conc", "decrease_conc"):
            if not species:
                raise ChemMCPError("Species name is required for concentration disturbances.")

            is_increase = dist == "increase_conc"
            is_reactant = species_role.lower().strip() == "reactant"

            if is_increase and is_reactant:
                shift = "forward"
                reasoning_parts.append(
                    f"增加反应物 [{species}] 的浓度。"
                    f"根据勒夏特列原理，系统将消耗增加的反应物，向生成物方向移动（正向移动）。"
                    f"Q < K，反应正向进行以重新建立平衡。"
                )
            elif is_increase and not is_reactant:
                shift = "backward"
                reasoning_parts.append(
                    f"增加生成物 [{species}] 的浓度。"
                    f"系统将消耗增加的生成物，向反应物方向移动（逆向移动）。"
                    f"Q > K，反应逆向进行。"
                )
            elif not is_increase and is_reactant:
                shift = "backward"
                reasoning_parts.append(
                    f"减少反应物 [{species}] 的浓度。"
                    f"系统将补充减少的反应物，逆向移动以产生更多反应物。"
                    f"Q > K，反应逆向进行。"
                )
            else:
                shift = "forward"
                reasoning_parts.append(
                    f"减少生成物 [{species}] 的浓度。"
                    f"系统将补充减少的生成物，正向移动。"
                    f"Q < K，反应正向进行。"
                )

            # Q/K analysis
            if current_Q is not None and K_value is not None:
                if shift == "forward":
                    q_vs_k = f"Q ({current_Q}) < K ({K_value}) → 正向移动"
                else:
                    q_vs_k = f"Q ({current_Q}) > K ({K_value}) → 逆向移动"

        # --- Pressure / Volume Disturbances ---
        elif dist in ("increase_pressure", "decrease_volume"):
            # Increasing pressure = decreasing volume → same effect
            if delta_n_gas > 0:
                shift = "backward"
                reasoning_parts.append(
                    f"增压（或减小体积）。Δn_gas = {delta_n_gas} > 0（生成物气体分子数更多）。"
                    f"系统向气体分子数少的一侧（反应物方向）移动以降低压力。"
                )
            elif delta_n_gas < 0:
                shift = "forward"
                reasoning_parts.append(
                    f"增压（或减小体积）。Δn_gas = {delta_n_gas} < 0（反应物气体分子数更多）。"
                    f"系统向气体分子数少的一侧（生成物方向）移动。"
                )
            else:
                shift = "no_shift"
                reasoning_parts.append(
                    f"增压（或减小体积）。Δn_gas = {delta_n_gas} = 0（两侧气体分子数相等）。"
                    f"压力变化对平衡位置无影响。"
                )
            k_effect = "K unchanged"

        elif dist in ("decrease_pressure", "increase_volume"):
            if delta_n_gas > 0:
                shift = "forward"
                reasoning_parts.append(
                    f"减压（或增大体积）。Δn_gas = {delta_n_gas} > 0。"
                    f"系统向气体分子数多的一侧（生成物方向）移动以增大压力。"
                )
            elif delta_n_gas < 0:
                shift = "backward"
                reasoning_parts.append(
                    f"减压（或增大体积）。Δn_gas = {delta_n_gas} < 0。"
                    f"系统向气体分子数多的一侧（反应物方向）移动。"
                )
            else:
                shift = "no_shift"
                reasoning_parts.append(
                    f"减压（或增大体积）。Δn_gas = 0，压力变化无影响。"
                )
            k_effect = "K unchanged"

        # --- Temperature Disturbances ---
        elif dist == "increase_temp":
            if rtype == "exothermic":
                shift = "backward"
                k_effect = "K decreases"
                reasoning_parts.append(
                    f"升温 + 放热反应。热视为「生成物」。"
                    f"升温→系统吸热→向逆反应方向（吸热方向）移动。"
                    f"K 减小（平衡常数随温度变化：放热反应升温K减小）。"
                )
            else:
                shift = "forward"
                k_effect = "K increases"
                reasoning_parts.append(
                    f"升温 + 吸热反应。热视为「反应物」。"
                    f"升温→系统吸热→向正反应方向（吸热方向）移动。"
                    f"K 增大（吸热反应升温K增大）。"
                )

        elif dist == "decrease_temp":
            if rtype == "exothermic":
                shift = "forward"
                k_effect = "K increases"
                reasoning_parts.append(
                    f"降温 + 放热反应。系统放热→正向移动（放热方向）。"
                    f"K 增大（放热反应降温K增大）。"
                )
            else:
                shift = "backward"
                k_effect = "K decreases"
                reasoning_parts.append(
                    f"降温 + 吸热反应。系统放热→逆向移动（放热方向=逆向）。"
                    f"K 减小（吸热反应降温K减小）。"
                )

        # --- Inert Gas ---
        elif dist == "add_inert_constant_v":
            shift = "no_shift"
            k_effect = "K unchanged"
            reasoning_parts.append(
                "恒容条件下加入惰性气体。总压增大但各组分分压不变。"
                "平衡不移动（Q值不变，仍等于K）。"
            )

        elif dist == "add_inert_constant_p":
            if abs(delta_n_gas) < 1e-10:
                shift = "no_shift"
                reasoning_parts.append(
                    "恒压条件下加入惰性气体，但 Δn_gas = 0。体积膨胀对各组分分压影响相同，平衡不移动。"
                )
            elif delta_n_gas > 0:
                shift = "forward"
                reasoning_parts.append(
                    f"恒压加入惰性气体→体积膨胀→各组分分压降低→等效于减压。"
                    f"Δn_gas = {delta_n_gas} > 0，向气体分子数多的方向（正向）移动。"
                )
            else:
                shift = "backward"
                reasoning_parts.append(
                    f"恒压加入惰性气体→体积膨胀→等效于减压。"
                    f"Δn_gas = {delta_n_gas} < 0，向气体分子数多的方向（逆向）移动。"
                )
            k_effect = "K unchanged"

        # --- Catalyst ---
        elif dist == "add_catalyst":
            shift = "no_shift"
            k_effect = "K unchanged"
            reasoning_parts.append(
                "加入催化剂。催化剂同等加快正逆反应速率，缩短达到平衡的时间。"
                "不改变平衡位置和K值。"
            )
            quant_note = "催化剂降低活化能 Ea，使正逆反应速率常数 k₁ 和 k₋₁ 同等倍数增加，故 K = k₁/k₋₁ 不变。"

        reasoning = "\n".join(reasoning_parts)

        logger.info(f"LeChatelierAnalyzer: {dist} on {rtype}(Δn={delta_n_gas}) → {shift}")
        return {
            "shift_direction": shift,
            "reasoning": reasoning,
            "K_effect": k_effect,
            "Q_vs_K": q_vs_k,
            "quantitative_note": quant_note,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            rtype = parts[0]
            dn = float(parts[1])
            dist = parts[2]
            species = parts[3] if len(parts) > 3 else None
            role = parts[4] if len(parts) > 4 else "reactant"
            Q = float(parts[5]) if len(parts) > 5 else None
            K = float(parts[6]) if len(parts) > 6 else None
            return self._run_base(rtype, dn, dist, species, role, Q, K)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'reaction_type delta_n disturbance [species] [role] [Q] [K]'")
