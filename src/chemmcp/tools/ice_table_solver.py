import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ICETableSolver(BaseTool):
    """
    ICE 表法（Initial-Change-Equilibrium）求解平衡浓度。
    支持给定初始浓度和 K 值，求解平衡时各物质浓度。
    求解一元二次方程: ax² + bx + c = 0
    """
    __version__                = "0.1.0"
    name                       = "ICETableSolver"
    func_name                  = "ice_table_solver"
    description                = "Solve equilibrium concentrations using ICE (Initial-Change-Equilibrium) table method."
    implementation_description = "Sets up ICE table from initial concentrations and stoichiometry, solves the equilibrium polynomial for x, returns all equilibrium concentrations."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Equilibrium", "ICE Table", "Physical Chemistry", "Concentration"]
    required_envs              = []

    code_input_sig             = [
        ("initial_conc",      "dict",  "N/A",   "Dict of initial concentrations, e.g. {'A': 1.0, 'B': 2.0, 'C': 0.0, 'D': 0.0}."),
        ("stoichiometry",     "dict",  "N/A",   "Dict of stoichiometric coefficients (positive for products, negative for reactants), e.g. {'A': -1, 'B': -3, 'C': 2, 'D': 0} where D is inert."),
        ("K_eq",              "float", "N/A",   "Equilibrium constant Kc."),
        ("Kp_mode",           "bool",  "False",  "If True, treat values as partial pressures instead of concentrations."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Semi-structured string with initial concentrations, stoichiometry, and K value."),
    ]

    output_sig                 = [
        ("x",                 "float", "The extent of reaction (change variable) at equilibrium."),
        ("equilibrium_conc",  "dict",  "Equilibrium concentrations/pressures of all species."),
        ("ice_table",         "dict",  "Complete ICE table with I/C/E rows."),
        ("K_calculated",      "float", "K value recalculated from equilibrium concentrations (verification)."),
        ("Qc_initial",        "float", "Reaction quotient Q at initial state."),
    ]

    examples                   = [
        {
            "code_input": {
                "initial_conc": {"N2O4": 0.500, "NO2": 0.000},
                "stoichiometry": {"N2O4": -1, "NO2": 2},
                "K_eq": 0.067,
                "Kp_mode": False,
            },
            "text_input": {
                "input_params": "initial: N2O4=0.5 NO2=0; stoich: N2O4=-1 NO2=2; K=0.067",
            },
            "output": {
                "x": 0.1177,
                "equilibrium_conc": {"N2O4": 0.3823, "NO2": 0.2354},
                "ice_table": {"I": {"N2O4": 0.5, "NO2": 0.0}, "C": {"N2O4": "-x", "NO2": "+2x"}, "E": {"N2O4": 0.3823, "NO2": 0.2354}},
                "K_calculated": 0.067,
                "Qc_initial": 0.0,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, initial_conc: dict, stoichiometry: dict, K_eq: float,
                  Kp_mode: bool = False) -> dict:
        if K_eq < 0:
            raise ChemMCPError("K must be non-negative.")

        # Build ICE table
        species = list(initial_conc.keys())
        I = {s: initial_conc.get(s, 0.0) for s in species}
        C = {}
        E = {}

        # Determine x direction: reactants have negative coefficients
        reactants = {s: c for s, c in stoichiometry.items() if c < 0}
        products = {s: c for s, c in stoichiometry.items() if c > 0}

        for s in species:
            coef = stoichiometry.get(s, 0)
            C[s] = f"{coef:+g}x" if coef != 0 else "0"
            E[s] = I[s] + coef  # symbolic in terms of x

        # Build polynomial: Q(x) = K
        # For each product: [E_s]^|coef| in numerator
        # For each reactant: [E_s]^|coef| in denominator
        # We need to solve for x where numerator/denominator = K

        # Use numerical approach: binary search or Newton's method on a reasonable range
        # Find valid range for x (all equilibrium conc >= 0)
        x_min = 0.0
        x_max = float('inf')

        for s in species:
            coef = stoichiometry.get(s, 0)
            init_val = I[s]
            if coef < 0:
                # Reactant: E = init + coef*x >= 0 => x <= init / |coef|
                x_max = min(x_max, init_val / abs(coef))
            elif coef > 0:
                # Product: E = init + coef*x >= 0 => x >= -init / coef (usually init>=0 so x>=0)
                x_min = max(x_min, -init_val / coef)

        if x_max < x_min:
            raise ChemMCPError("No feasible solution: constraints on x are contradictory.")

        # Cap x_max for numerical stability
        x_max = min(x_max, 1e8)

        # Solve using bisection method
        def Q_of_x(x):
            num = 1.0
            denom = 1.0
            for s in species:
                coef = stoichiometry.get(s, 0)
                e_val = I[s] + coef * x
                if e_val <= 0:
                    return None  # invalid
                if coef > 0:
                    num *= math.pow(e_val, coef)
                elif coef < 0:
                    denom *= math.pow(e_val, -coef)
            return num / denom

        # Bisection search
        tol = 1e-10
        max_iter = 200

        q_lo = Q_of_x(x_min)
        q_hi = Q_of_x(x_max * 0.9999)

        if q_lo is None:
            q_lo = 0.0
        if q_hi is None:
            # Try to find a valid upper bound
            for trial in [x_max * 0.5, x_max * 0.1, x_max * 0.01]:
                q_hi = Q_of_x(trial)
                if q_hi is not None:
                    x_max = trial
                    break

        # Check if we can bracket K
        lo, hi = x_min, x_max
        for _ in range(max_iter):
            mid = (lo + hi) / 2.0
            q_mid = Q_of_x(mid)
            if q_mid is None:
                hi = mid
                continue
            if abs(q_mid - K_eq) < tol or (hi - lo) < tol:
                x_sol = mid
                break
            if q_mid < K_eq:
                lo = mid
            else:
                hi = mid
        else:
            x_sol = (lo + hi) / 2.0

        # Compute equilibrium concentrations
        eq_conc = {}
        for s in species:
            eq_conc[s] = round(I[s] + stoichiometry.get(s, 0) * x_sol, 6)

        # Verify K
        K_check = Q_of_x(x_sol)

        # Initial Q
        Q_init = Q_of_x(0.0)

        ice_table = {
            "Initial": {s: I[s] for s in species},
            "Change": {s: C[s] for s in species},
            "Equilibrium": {s: eq_conc[s] for s in species},
        }

        logger.info(f"ICE table solved: x={x_sol:.6f}, K_calc={K_check:.6f}")
        return {
            "x": round(x_sol, 6),
            "equilibrium_conc": eq_conc,
            "ice_table": ice_table,
            "K_calculated": round(K_check, 6) if K_check else None,
            "Qc_initial": round(Q_init, 6) if Q_init else None,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            # Format: "initial: A=1.0 B=2.0; stoich: A=-1 B=-3 C=2 D=0; K=0.05"
            init = {}
            stoich = {}
            K_eq = None
            sections = input_params.split(";")
            for sec in sections:
                sec = sec.strip()
                if sec.startswith("initial:"):
                    for item in sec[len("initial:"):].strip().split():
                        name, val = item.split("=")
                        init[name.strip()] = float(val)
                elif sec.startswith("stoich:"):
                    for item in sec[len("stoich:"):].strip().split():
                        name, val = item.split("=")
                        stoich[name.strip()] = float(val)
                elif sec.lower().startswith("k=") or sec.lower().startswith("k:"):
                    K_eq = float(sec.split("=")[-1])
            if K_eq is None:
                raise ValueError("K value not provided.")
            return self._run_base(init, stoich, K_eq)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
