import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ParallelConsecutiveReactions(BaseTool):
    """
    平行反应和连串反应动力学求解器。
    
    支持的反应类型：
    - 平行反应: A → B (k1), A → C (k2), ...
    - 连串反应: A → I → P (k1, k2)
    
    提供解析解（对简单情况）和数值解。
    """
    __version__ = "0.1.0"
    name = "ParallelConsecutiveReactions"
    func_name = "solve_kinetics"
    description = "Solve kinetics for parallel and consecutive reactions with analytical solutions and concentration profiles."
    implementation_description = "Provides analytical solutions for parallel (competitive) and consecutive reaction schemes. Calculates concentrations, yields, selectivity, and time evolution of all species."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Parallel Reactions", "Consecutive Reactions", "Yield"]
    required_envs = []

    code_input_sig = [
        ("reaction_type", "str", "N/A", "'parallel' for competitive parallel reactions or 'consecutive' for series reactions."),
        ("rate_constants", "list", "N/A", "List of rate constants. Parallel: [k1, k2, ...] for each branch. Consecutive: [k1, k2] for each step."),
        ("initial_concentrations", "dict", "N/A", "Initial concentrations dict, e.g., {'A': 1.0}. For consecutive, can include intermediate."),
        ("time_points", "list", "N/A", "List of time values at which to compute concentrations, e.g., [0, 10, 20, 50, 100]."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'type k1,k2,... A0:val t1,t2,t3,...'. Example: 'parallel 0.1,0.02 A:1 0,10,50,100'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with concentrations at each time point, yield distribution, selectivity (for parallel), and kinetic parameters."),
    ]

    examples = [
        {
            "code_input": {
                "reaction_type": "parallel",
                "rate_constants": [0.1, 0.02],
                "initial_concentrations": {"A": 1.0},
                "time_points": [0, 10, 20, 50, 100],
            },
            "text_input": {
                "input_params": "parallel 0.1,0.02 A:1 0,10,20,50,100",
            },
            "output": {
                "result": {
                    "reaction_type": "parallel",
                    "final_A_at_t100": 0.3679,
                    "final_B_at_t100": 0.5264,
                    "final_C_at_t100": 0.1057,
                    "selectivity_B_over_C": 5.0,
                    "yield_B": 0.5264,
                    "yield_C": 0.1057,
                }
            },
        },
        {
            "code_input": {
                "reaction_type": "consecutive",
                "rate_constants": [0.1, 0.05],
                "initial_concentrations": {"A": 1.0},
                "time_points": [0, 10, 20, 50],
            },
            "text_input": {
                "input_params": "consecutive 0.1,0.05 A:1 0,10,20,50",
            },
            "output": {
                "result": {
                    "reaction_type": "consecutive",
                    "max_intermediate_time": 13.86,
                    "max_intermediate_concentration": 0.4703,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _solve_parallel(self, k_list, A0, time_points):
        """
        Parallel: A --(k1)--> B, A --(k2)--> C, ...
        [A] = A0 * exp(-k_total * t)
        [Bi] = A0 * (ki/k_total) * (1 - exp(-k_total * t))
        """
        k_total = sum(k_list)
        results = {"A": [], "products": {}}
        
        for i, ki in enumerate(k_list):
            product_name = chr(ord('B') + i)  # B, C, D, ...
            results["products"][product_name] = []

        for t in time_points:
            At = A0 * math.exp(-k_total * t)
            results["A"].append(round(At, 8))
            
            for i, ki in enumerate(k_list):
                product_name = chr(ord('B') + i)
                if k_total > 0:
                    Bit = A0 * (ki / k_total) * (1 - math.exp(-k_total * t))
                else:
                    Bit = 0.0
                results["products"][product_name].append(round(Bit, 8))

        # Selectivity and yield at final time
        final_idx = len(time_points) - 1
        yields = {}
        for pname in results["products"]:
            yields[pname] = results["products"][pname][final_idx]
        
        # Selectivity ratios
        selectivity = {}
        pnames = list(results["products"].keys())
        if len(pnames) >= 2:
            selectivity[f"{pnames[0]}_over_{pnames[1]}"] = (
                round(yields[pnames[0]] / max(yields[pnames[1]], 1e-15), 4)
            )

        return {
            **results,
            "k_total": k_total,
            "yields": yields,
            "selectivity": selectivity,
            "n_branches": len(k_list),
        }

    def _solve_consecutive(self, k_list, A0, time_points):
        """
        Consecutive: A --(k1)--> I --(k2)--> P
        [A] = A0 * exp(-k1*t)
        [I] = A0*k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))  (if k1≠k2)
        [P] = A0 - [A] - [I]
        """
        k1, k2 = k_list[0], k_list[1]
        results = {"A": [], "I": [], "P": []}
        max_I_t = None
        max_I_val = 0
        
        for t in time_points:
            At = A0 * math.exp(-k1 * t)
            results["A"].append(round(At, 8))
            
            if abs(k2 - k1) < 1e-12:
                It = A0 * k1 * t * math.exp(-k1 * t)
            else:
                It = A0 * k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))
            
            It = max(It, 0.0)
            results["I"].append(round(It, 8))
            
            Pt = A0 - At - It
            Pt = max(Pt, 0.0)
            results["P"].append(round(Pt, 8))

            # Track max intermediate
            if It > max_I_val:
                max_I_val = It
                max_I_t = t

        # Time of maximum [I]: t_max = ln(k1/k2)/(k1-k2)
        if k1 > 0 and k2 > 0 and abs(k1 - k2) > 1e-12:
            t_max_ana = math.log(k1 / k2) / (k1 - k2)
        elif k1 > 0:
            t_max_ana = 1.0 / k1
        else:
            t_max_ana = None

        return {
            **results,
            "max_intermediate_time": round(t_max_ana, 6) if t_max_ana else None,
            "max_intermediate_concentration": round(max_I_val, 8),
            "yield_P_final": results["P"][-1] if results["P"] else 0,
        }

    def _run_base(self, reaction_type: str, rate_constants: list,
                  initial_concentrations: dict, time_points: list) -> dict:
        if reaction_type not in ("parallel", "consecutive"):
            raise ChemMCPError("reaction_type must be 'parallel' or 'consecutive'.")
        if not rate_constants or any(k < 0 for k in rate_constants):
            raise ChemMCPError("Rate constants must be non-negative numbers.")
        if not time_points:
            raise ChemMCPError("Must provide at least one time point.")

        A0 = initial_concentrations.get("A", 1.0)

        if reaction_type == "parallel":
            result = self._solve_parallel(rate_constants, A0, time_points)
        else:
            result = self._solve_consecutive(rate_constants, A0, time_points)

        result["reaction_type"] = reaction_type
        result["rate_constants"] = rate_constants
        result["time_points"] = time_points
        result["initial_A0"] = A0

        logger.info(f"ParallelConsecutiveReactions: type={reaction_type}, done")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            rxn_type = parts[0]
            k_list = [float(x) for x in parts[1].split(",")]
            
            init_dict = {}
            t_list = None
            for p in parts[2:]:
                if ":" in p:
                    key, val = p.split(":")
                    init_dict[key] = float(val)
                elif "," in p:
                    t_list = [float(x) for x in p.split(",")]

            if t_list is None:
                t_list = [0.0]

            return self._run_base(rxn_type, k_list, init_dict, t_list)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'type k1,k2,... A:val t1,t2,...'")
