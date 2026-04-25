import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ReactionQuotient(BaseTool):
    """
    计算反应商 Q 并判断反应方向。
    Q 与 K 的比较决定反应进行的方向：
      - Q < K → 正向进行（→ 产物）
      - Q = K → 平衡状态
      - Q > K → 逆向进行（→ 反应物）
    """
    __version__                = "0.1.0"
    name                       = "ReactionQuotient"
    func_name                  = "reaction_quotient"
    description                = "Calculate reaction quotient Q and predict reaction direction by comparing with K."
    implementation_description = "Computes Q from current concentrations/pressures using mass action expression, compares with equilibrium constant K to determine reaction direction."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Reaction Quotient", "Equilibrium", "Physical Chemistry", "Direction"]
    required_envs              = []

    code_input_sig             = [
        ("current_state",     "dict",  "N/A",   "Dict of current concentrations or partial pressures, e.g. {'A': 0.1, 'B': 0.2, 'C': 0.01, 'D': 0.01}."),
        ("stoichiometry",     "dict",  "N/A",   "Stoichiometric coefficients (negative for reactants, positive for products), e.g. {'A': -1, 'B': -3, 'C': 2, 'D': 1}."),
        ("K_eq",              "float", "N/A",   "Equilibrium constant Kc or Kp."),
        ("mode",              "str",   "kc",     "'kc' for concentrations, 'kp' for partial pressures."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Format: 'A=0.1 B=0.2 | C=0.01 D=0.01 | A=-1 B=-3 C=2 D=1 | K=10 kc'"),
    ]

    output_sig                 = [
        ("Q",                 "float", "Reaction quotient value."),
        ("K",                 "float", "Equilibrium constant value."),
        ("direction",         "str",   "'forward', 'backward', or 'equilibrium'."),
        ("comparison",        "str",   "String showing Q vs K relationship (e.g., 'Q < K')."),
        ("expression",        "str",   "The mass action expression used."),
        ("explanation",       "str",   "Plain-language explanation of what this means for the reaction."),
    ]

    examples                   = [
        {
            "code_input": {
                "current_state": {"N2": 0.01, "H2": 0.01, "NH3": 0.5},
                "stoichiometry": {"N2": -1, "H2": -3, "NH3": 2},
                "K_eq": 1500.0,
                "mode": "kc",
            },
            "text_input": {
                "input_params": "N2=0.01 H2=0.01 NH3=0.5 | N2=-1 H2=-3 NH3=2 | K=1500 kc",
            },
            "output": {
                "Q": 500000.0,
                "K": 1500.0,
                "direction": "backward",
                "comparison": "Q >> K",
                "expression": "Qc = [NH3]^2 / ([N2] × [H2]^3)",
                "explanation": "Q (500000) > K (1500): too much product relative to equilibrium. Reaction will proceed in reverse to form more reactants.",
            }
        },
        {
            "code_input": {
                "current_state": {"NO2": 0.010, "N2O4": 0.100},
                "stoichiometry": {"NO2": 2, "N2O4": -1},
                "K_eq": 6.7,
                "mode": "kc",
            },
            "text_input": {
                "input_params": "NO2=0.01 N2O4=0.1 | NO2=2 N2O4=-1 | K=6.7 kc",
            },
            "output": {
                "Q": 0.001,
                "K": 6.7,
                "direction": "forward",
                "comparison": "Q < K",
                "expression": "Qc = [NO2]^2 / [N2O4]",
                "explanation": "Q (0.001) < K (6.7): more reactant than at equilibrium. Reaction will proceed forward to form more product.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, current_state: dict, stoichiometry: dict, K_eq: float,
                  mode: str = "kc") -> dict:
        mode = mode.lower()
        if mode not in ("kc", "kp"):
            raise ChemMCPError("Mode must be 'kc' or 'kp'.")
        if K_eq < 0:
            raise ChemMCPError("K must be non-negative.")

        species = list(stoichiometry.keys())
        bracket = "[" if mode == "kc" else "P_"

        # Build mass action expression
        num_parts = []
        denom_parts = []
        numerator = 1.0
        denominator = 1.0

        for s in species:
            coef = stoichiometry[s]
            val = current_state.get(s, 0.0)
            if val < 0:
                raise ChemMCPError(f"Value for {s} cannot be negative.")

            abs_coef = abs(coef)
            if coef > 0:
                numerator *= math.pow(val, coef)
                if abs_coef == 1:
                    num_parts.append(f"{bracket}{s}{']' if mode == 'kc' else ''}")
                else:
                    num_parts.append(f"{bracket}{s}{']' if mode == 'kc' else ''}^{abs_coef}")
            elif coef < 0:
                denominator *= math.pow(val, -coef)
                if abs_coef == 1:
                    denom_parts.append(f"{bracket}{s}{']' if mode == 'kc' else ''}")
                else:
                    denom_parts.append(f"{bracket}{s}{']' if mode == 'kc' else ''}^{abs_coef}")

        if denominator == 0:
            raise ChemMCPError("Denominator is zero — a reactant value is zero.")

        Q = numerator / denominator

        # Determine direction
        tol = 1e-6 * max(K_eq, 1.0)
        if abs(Q - K_eq) <= tol:
            direction = "equilibrium"
            comparison = "Q ≈ K"
            explanation = (
                f"Q ({Q:.4g}) ≈ K ({K_eq:.4g}): system is at or very near equilibrium. "
                f"No net reaction direction."
            )
        elif Q < K_eq:
            direction = "forward"
            comparison = "Q < K"
            explanation = (
                f"Q ({Q:.4g}) < K ({K_eq:.4g}): there are relatively more reactants than at equilibrium. "
                f"The reaction will proceed in the forward direction (toward products) to reach equilibrium."
            )
        else:
            direction = "backward"
            comparison = "Q > K"
            explanation = (
                f"Q ({Q:.4g}) > K ({K_eq:.4g}): there are relatively more products than at equilibrium. "
                f"The reaction will proceed in the reverse direction (toward reactants) to reach equilibrium."
            )

        label = "Qc" if mode == "kc" else "Qp"
        K_label = "Kc" if mode == "kc" else "Kp"

        num_str = " × ".join(num_parts) if num_parts else "1"
        denom_str = " × ".join(denom_parts) if denom_parts else "1"
        expr = f"{label} = {num_str}" + (f" / {denom_str}" if denom_str != "1" else "")

        logger.info(f"Reaction Quotient: Q={Q:.4e}, K={K_eq}, direction={direction}")
        return {
            "Q": round(Q, 6),
            "K": round(K_eq, 6),
            "direction": direction,
            "comparison": comparison,
            "expression": expr,
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            # Format: "A=0.1 B=0.2 | A=-1 B=-3 C=2 D=1 | K=10 kc" (3 or 4 parts)
            parts = input_params.split("|")
            if len(parts) < 3:
                raise ValueError("Need format: 'state | stoich | K [mode]'")

            state = {}
            for item in parts[0].strip().split():
                name, val = item.split("=")
                state[name.strip()] = float(val)

            stoich = {}
            for item in parts[1].strip().split():
                name, val = item.split("=")
                stoich[name.strip()] = float(val)

            # Part 2: "K=val kc" or "K=val"
            k_mode_part = parts[2].strip().split()
            K_eq = float(k_mode_part[0].split("=")[1])
            mode = k_mode_part[1] if len(k_mode_part) > 1 else "kc"
            return self._run_base(state, stoich, K_eq, mode)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'A=v B=v | A=c B=c | K=val kc/kp'")
