import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalculateEquilibriumConstant(BaseTool):
    """
    根据平衡时的浓度/分压计算平衡常数 Kc 或 Kp。
    对于反应 aA + bB ⇌ cC + dD:
      Kc = [C]^c × [D]^d / ([A]^a × [B]^b)
      Kp = (P_C)^c × (P_D)^d / ((P_A)^a × (P_B)^b)
    """
    __version__                = "0.1.0"
    name                       = "CalculateEquilibriumConstant"
    func_name                  = "calculate_equilibrium_constant"
    description                = "Calculate equilibrium constant Kc or Kp from equilibrium concentrations or partial pressures."
    implementation_description = "Computes Kc from molar concentrations or Kp from partial pressures using the mass action expression."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Equilibrium", "Kc", "Kp", "Physical Chemistry"]
    required_envs              = []

    code_input_sig             = [
        ("products",          "list",  "N/A",   "List of dicts: [{'name': 'C', 'coefficient': 2, 'value': 0.3}, ...] for products at equilibrium."),
        ("reactants",         "list",  "N/A",   "List of dicts: [{'name': 'A', 'coefficient': 1, 'value': 0.1}, ...] for reactants at equilibrium."),
        ("mode",              "str",   "kc",     "'kc' for concentration, 'kp' for partial pressure."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Format: 'reactant1:coef:value reactant2:coef:value | product1:coef:value | kc/kp'"),
    ]

    output_sig                 = [
        ("K",                 "float", "Equilibrium constant value."),
        ("expression",        "str",   "The mass action expression used."),
        ("mode",              "str",   "'Kc' or 'Kp'."),
    ]

    examples                   = [
        {
            "code_input": {
                "products": [{"name": "NO2", "coefficient": 2, "value": 0.056}],
                "reactants": [{"name": "N2O4", "coefficient": 1, "value": 0.032}],
                "mode": "kc",
            },
            "text_input": {
                "input_params": "N2O4:1:0.032 | NO2:2:0.056 | kc",
            },
            "output": {
                "K": 0.098,
                "expression": "Kc = [NO2]^2 / [N2O4]",
                "mode": "Kc",
            }
        },
        {
            "code_input": {
                "products": [{"name": "H2", "coefficient": 3, "value": 1.2}, {"name": "N2", "coefficient": 1, "value": 0.4}],
                "reactants": [{"name": "NH3", "coefficient": 2, "value": 0.8}],
                "mode": "kp",
            },
            "text_input": {
                "input_params": "NH3:2:0.8 | H2:3:1.2 N2:1:0.4 | kp",
            },
            "output": {
                "K": 1.296,
                "expression": "Kp = (P_H2)^3 × (P_N2) / (P_NH3)^2",
                "mode": "Kp",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, products: list, reactants: list, mode: str = "kc") -> dict:
        mode = mode.lower()
        if mode not in ("kc", "kp"):
            raise ChemMCPError("Mode must be 'kc' or 'kp'.")

        # Build numerator (products) and denominator (reactants)
        num_expr_parts = []
        denom_expr_parts = []

        numerator = 1.0
        for p in products:
            coef = p["coefficient"]
            val = p["value"]
            if val < 0:
                raise ChemMCPError(f"Concentration/pressure cannot be negative for {p['name']}.")
            numerator *= math.pow(val, coef)
            if coef == 1:
                num_expr_parts.append(f"[{p['name']}]" if mode == "kc" else f"P_{p['name']}")
            else:
                num_expr_parts.append(f"[{p['name']}]^{coef}" if mode == "kc" else f"P_{p['name']}^{coef}")

        denominator = 1.0
        for r in reactants:
            coef = r["coefficient"]
            val = r["value"]
            if val < 0:
                raise ChemMCPError(f"Concentration/pressure cannot be negative for {r['name']}.")
            denominator *= math.pow(val, coef)
            if coef == 1:
                denom_expr_parts.append(f"[{r['name']}]" if mode == "kc" else f"P_{r['name']}")
            else:
                denom_expr_parts.append(f"[{r['name']}^{coef}" if mode == "kc" else f"P_{r['name']}^{coef}")

        if denominator == 0:
            raise ChemMCPError("Denominator is zero — a reactant concentration/pressure is zero.")

        K = numerator / denominator

        label = "Kc" if mode == "kc" else "Kp"
        bracket = "[" if mode == "kc" else "P_"
        num_str = " × ".join(num_expr_parts) if len(num_expr_parts) > 1 else (num_expr_parts[0] if num_expr_parts else "1")
        denom_str = " × ".join(denom_expr_parts) if len(denom_expr_parts) > 1 else (denom_expr_parts[0] if denom_expr_parts else "1")
        expression = f"{label} = {num_str} / {denom_str}"

        logger.info(f"Calculated {label} = {K:.4e}")
        return {
            "K": round(K, 6),
            "expression": expression,
            "mode": label,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            # Format: "reactant1:coef:value reactant2:coef:value | product1:coef:value product2:coef:value | kc/kp"
            parts = input_params.split("|")
            if len(parts) < 3:
                raise ValueError("Need format: 'reactants | products | mode'")

            reactants = []
            for item in parts[0].strip().split():
                name, coef, val = item.split(":")
                reactants.append({"name": name.strip(), "coefficient": float(coef), "value": float(val)})

            products = []
            for item in parts[1].strip().split():
                name, coef, val = item.split(":")
                products.append({"name": name.strip(), "coefficient": float(coef), "value": float(val)})

            mode = parts[2].strip()
            return self._run_base(products, reactants, mode)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'R1:c:v R2:c:v | P1:c:v P2:c:v | kc/kp'")
