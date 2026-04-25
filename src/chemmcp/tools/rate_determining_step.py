import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RateDeterminingStep(BaseTool):
    """
    分析并识别多步反应中的速控步骤（Rate-Determining Step, RDS）。
    通过比较各步骤的速率常数，识别最慢的步骤，并推导总反应速率方程。
    """
    __version__ = "0.1.0"
    name = "RateDeterminingStep"
    func_name = "analyze_rds"
    description = "Identify and analyze the rate-determining step (RDS) in a multi-step reaction mechanism."
    implementation_description = "Compares rate constants of each step to identify the slowest (rate-determining) step. Derives the overall rate law based on which step is RDS. Supports analysis with or without pre-equilibrium assumptions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "RDS", "Mechanism", "Rate Law"]
    required_envs = []

    code_input_sig = [
        ("mechanism_steps", "list", "N/A", "List of mechanism steps, each as a dict: {'reactants': str, 'products': str, 'k': float, 'reversible': bool}. Example: [{'reactants':'A->B','k':0.1,'reversible':False}, {'reactants':'B->C','k':10,'reversible':False}]."),
        ("has_pre_equilibrium", "bool", "False", "Whether fast pre-equilibrium steps exist before RDS."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Semicolon-separated steps: 'step1_desc;k1;rev?;step2_desc;k2;rev?;...' where rev?=true/false. Example: 'A->B;0.1;false;B->C;10;false'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with RDS index, RDS description, overall rate law, rate constant ratios, and justification."),
    ]

    examples = [
        {
            "code_input": {
                "mechanism_steps": [
                    {"reactants": "A -> B", "products": "", "k": 0.001, "reversible": False},
                    {"reactants": "B + C -> D", "products": "", "k": 5.0, "reversible": False},
                    {"reactants": "D -> E", "products": "", "k": 100.0, "reversible": False},
                ],
                "has_pre_equilibrium": False,
            },
            "text_input": {
                "input_params": "A->B;0.001;false;B+C->D;5;false;D->E;100;false",
            },
            "output": {
                "result": {
                    "rds_step_index": 1,
                    "rds_step_description": "A -> B",
                    "rds_rate_constant": 0.001,
                    "overall_rate_law": "rate = k1[A]",
                    "slowest_to_fastest_ratio": 100000,
                    "justification": "Step 1 has the smallest rate constant (k=0.001), making it the bottleneck.",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, mechanism_steps: list, has_pre_equilibrium: bool = False) -> dict:
        if not mechanism_steps:
            raise ChemMCPError("Mechanism steps list cannot be empty.")

        n_steps = len(mechanism_steps)
        k_values = []
        for i, step in enumerate(mechanism_steps):
            k = step.get("k")
            if k is None or k < 0:
                raise ChemMCPError(f"Step {i+1}: invalid rate constant.")
            k_values.append(k)

        # Find RDS (smallest k for forward irreversible steps)
        min_k = min(k_values)
        rds_idx = k_values.index(min_k)
        rds_step = mechanism_steps[rds_idx]

        # Sort to get ranking
        sorted_ks = sorted(enumerate(k_values), key=lambda x: x[1])
        
        # Rate constant ratios relative to RDS
        ratios = [k / min_k if min_k > 0 else float('inf') for k in k_values]

        # Derive approximate rate law based on RDS
        step_desc = rds_step.get("reactants", f"Step {rds_idx+1}")
        
        # Simple rate law derivation from RDS reactants
        if has_pre_equilibrium and rds_idx > 0:
            rate_law = f"rate = k{rds_idx+1} * K_eq * [pre-equilibrium reactants]"
        else:
            # Extract species from step description
            import re
            species = re.findall(r'[A-Z][a-z0-9]*', step_desc)
            if species:
                rate_law = f"rate = k{rds_idx+1}" + "".join(f"[{s}]" for s in species)
            else:
                rate_law = f"rate = k{rds_idx+1} * [reactants]"

        # Build justification
        justification = (
            f"Step {rds_idx+1} ('{step_desc}') has the smallest rate constant "
            f"(k={min_k}), making it the bottleneck. "
        )
        if max(ratios) > 1000:
            justification += f"It is {max(ratios):.0f}x slower than the fastest step — clearly the RDS."
        elif max(ratios) > 10:
            justification += f"It is {max(ratios):.1f}x slower than the fastest step — likely the RDS."
        else:
            justification += "Rate constants are within an order of magnitude — no clear single RDS."

        result = {
            "rds_step_index": rds_idx + 1,  # 1-indexed for user display
            "rds_step_description": step_desc,
            "rds_rate_constant": min_k,
            "overall_rate_law": rate_law,
            "rate_constants": k_values,
            "rate_constant_ratios": [round(r, 2) for r in ratios],
            "step_ranking": [(idx+1, k) for idx, k in sorted_ks],
            "slowest_to_fastest_ratio": round(max(ratios) if ratios else 0, 2),
            "has_pre_equilibrium": has_pre_equilibrium,
            "justification": justification,
            "n_total_steps": n_steps,
        }

        logger.info(f"RateDeterminingStep: RDS = Step {rds_idx+1}, k={min_k}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split(";")
            steps = []
            i = 0
            while i < len(parts) - 2:
                desc = parts[i].strip()
                k_val = float(parts[i+1].strip())
                rev = parts[i+2].strip().lower() == "true" if i+2 < len(parts) else False
                steps.append({"reactants": desc, "products": "", "k": k_val, "reversible": rev})
                i += 3
            return self._run_base(steps)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'desc;k;rev;desc;k;rev;...'")
