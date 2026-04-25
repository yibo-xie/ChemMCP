import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PressureEffectEquilibrium(BaseTool):
    """
    分析压力变化对气相平衡的影响。
    计算压力变化后的新平衡组成，以及各组分分压和摩尔分数的变化。
    """
    __version__                = "0.1.0"
    name                       = "PressureEffectEquilibrium"
    func_name                  = "pressure_effect_equilibrium"
    description                = "Analyze the effect of pressure change on gas-phase equilibrium composition."
    implementation_description = "Given initial total pressure and equilibrium constant Kp, calculates new equilibrium partial pressures and mole fractions after a pressure change, using the ideal gas law and mass action expression."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Pressure", "Equilibrium", "Gas Phase", "Physical Chemistry"]
    required_envs              = []

    code_input_sig             = [
        ("initial_total_p",   "float", "N/A",   "Initial total pressure in atm."),
        ("new_total_p",       "float", "N/A",   "New total pressure after change in atm."),
        ("Kp",                "float", "N/A",   "Pressure equilibrium constant Kp."),
        ("stoichiometry",     "dict",  "N/A",   "Stoichiometric coefficients for gaseous species, e.g. {'A': -1, 'B': -1, 'C': 2}. Negative for reactants, positive for products."),
        ("initial_mole_frac", "dict",  "None",  "Initial mole fractions at old equilibrium (optional). If None, will solve from Kp and P_initial."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Semi-structured string with pressures, Kp, and stoichiometry."),
    ]

    output_sig                 = [
        ("old_equilibrium",   "dict",  "Old equilibrium: partial pressures and mole fractions."),
        ("new_equilibrium",   "dict",  "New equilibrium: partial pressures and mole fractions."),
        ("shift_direction",   "str",   "'forward', 'backward', or 'no_shift'."),
        ("degree_of_dissociation_change", "str", "Description of how degree of dissociation/reaction changed."),
    ]

    examples                   = [
        {
            "code_input": {
                "initial_total_p": 1.0,
                "new_total_p": 5.0,
                "Kp": 0.067,
                "stoichiometry": {"N2O4": -1, "NO2": 2},
                "initial_mole_frac": None,
            },
            "text_input": {
                "input_params": "P_old=1.0 P_new=5.0 Kp=0.067; N2O4=-1 NO2=2",
            },
            "output": {
                "old_equilibrium": {"partial_pressures": {"N2O4": 0.688, "NO2": 0.312}, "mole_fractions": {"N2O4": 0.688, "NO2": 0.312}},
                "new_equilibrium": {"partial_pressures": {"N2O4": 4.38, "NO2": 0.62}, "mole_fractions": {"N2O4": 0.876, "NO2": 0.124}},
                "shift_direction": "backward",
                "degree_of_dissociation_change": "Dissociation decreased as pressure increased (shift toward fewer moles).",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, initial_total_p: float, new_total_p: float, Kp: float,
                  stoichiometry: dict, initial_mole_frac: dict = None) -> dict:
        if initial_total_p <= 0 or new_total_p <= 0:
            raise ChemMCPError("Pressures must be positive.")
        if Kp < 0:
            raise ChemMCPError("Kp must be non-negative.")

        species = list(stoichiometry.keys())
        reactants = [s for s in species if stoichiometry[s] < 0]
        products = [s for s in species if stoichiometry[s] > 0]
        delta_n = sum(stoichiometry[s] for s in species if stoichiometry[s] > 0) + \
                  sum(stoichiometry[s] for s in species if stoichiometry[s] < 0)

        # --- Solve old equilibrium ---
        # Use degree of dissociation approach for simple cases
        # For general case: use numerical solver

        def solve_equilibrium(P_total, K_val):
            """Solve for equilibrium partial pressures at given total pressure."""
            # For a reaction like aA ⇌ bB (or similar), use numerical approach
            # Let's parameterize by extent of reaction ξ
            # Simple approach: assume one dominant reactant, find its conversion

            # Try range of extent values
            best_x = None
            best_diff = float('inf')

            # Determine reasonable bounds for x
            # Assume we start with mostly reactants
            x_lo = 0.0
            x_hi = 1.0  # max fraction of reactant that can react

            for _ in range(100):
                mid = (x_lo + x_hi) / 2.0

                # Calculate partial pressures based on extent
                # Simplified model: start with 1 mol total, all reactant(s)
                # At extent x: moles of each species
                n = {}
                # Initial: normalize to 1 mol of first reactant
                n_init = {}
                for s in species:
                    if stoichiometry[s] < 0:
                        n_init[s] = abs(stoichiometry[s])  # stoich amounts
                    else:
                        n_init[s] = 0.0

                total_init = sum(n_init.values())
                for s in species:
                    n_init[s] /= total_init  # normalize to 1 mol total

                # At extent x (fraction of limiting reactant reacted)
                n_eq = {}
                for s in species:
                    n_eq[s] = n_init[s] + stoichiometry[s] * x * min(
                        n_init[s] / abs(stoichiometry[s]) for s in reactants
                        if stoichiometry[s] < 0
                    ) if any(stoichiometry[s] < 0 for s in species) else n_init[s]

                # Recalculate more robustly
                pass

            # Use simpler direct approach: iterate on one variable
            # For reaction with Δn_gas ≠ 0, we can express everything in terms of alpha (dissociation fraction)
            return self._solve_by_iteration(P_total, K_val, stoichiometry)

        def _solve_by_iteration(self, P_tot, K_val, stoich):
            """Iterative solution for equilibrium at pressure P_tot."""
            sp = list(stoich.keys())
            rcts = [s for s in sp if stoich[s] < 0]
            prods = [s for s in sp if stoich[s] > 0]

            # Parameterize by extent xi (moles of reference reactant consumed)
            # Start with 1 mol of each reactant in stoich ratio
            ref_rct = rcts[0] if rcts else sp[0]
            ref_coef = abs(stoich[ref_rct])

            # Initial moles (in stoich proportions)
            n0 = {}
            for s in sp:
                if stoich[s] < 0:
                    n0[s] = abs(stoich[s]) / ref_coef  # normalized
                else:
                    n0[s] = 0.0

            # Search over xi
            lo, hi = 0.0, n0[ref_rct]
            best_xi = 0.0
            best_err = float('inf')

            for _ in range(200):
                mid = (lo + hi) / 2.0
                n = {}
                valid = True
                for s in sp:
                    n[s] = n0[s] + stoich[s] * mid
                    if n[s] < -1e-10:
                        valid = False
                        break
                    n[s] = max(n[s], 0.0)
                if not valid:
                    hi = mid
                    continue

                n_total = sum(n.values())
                if n_total < 1e-15:
                    lo = mid
                    continue

                # Partial pressures
                P = {s: n[s] / n_total * P_tot for s in sp}

                # Calculate Q
                num = 1.0
                den = 1.0
                for s in sp:
                    c = stoich[s]
                    if c > 0:
                        num *= P[s] ** c
                    elif c < 0:
                        den *= P[s] ** (-c)

                Q = num / den if den > 0 else float('inf')
                err = abs(Q - K_val)

                if err < best_err:
                    best_err = err
                    best_xi = mid

                if err < 1e-8:
                    break
                if Q < K_val:
                    lo = mid  # need more products → increase extent
                else:
                    hi = mid

            # Final calculation at best_xi
            n_final = {}
            for s in sp:
                n_final[s] = max(n0[s] + stoich[s] * best_xi, 0.0)
            n_tot = sum(n_final.values())
            P_final = {s: n_final[s] / n_tot * P_tot for s in sp} if n_tot > 0 else {s: 0.0 for s in sp}
            mf = {s: n_final[s] / n_tot if n_tot > 0 else 0.0 for s in sp}

            return {
                "partial_pressures": {s: round(P_final[s], 6) for s in sp},
                "mole_fractions": {s: round(mf[s], 6) for s in sp},
                "extent": round(best_xi, 6),
            }

        old_eq = _solve_by_iteration(self, initial_total_p, Kp, stoichiometry)
        new_eq = _solve_by_iteration(self, new_total_p, Kp, stoichiometry)

        # Determine shift direction
        if delta_n > 0:
            if new_total_p > initial_total_p:
                shift = "backward"
            elif new_total_p < initial_total_p:
                shift = "forward"
            else:
                shift = "no_shift"
        elif delta_n < 0:
            if new_total_p > initial_total_p:
                shift = "forward"
            elif new_total_p < initial_total_p:
                shift = "backward"
            else:
                shift = "no_shift"
        else:
            shift = "no_shift"

        desc_map = {
            "forward": "Reaction shifted toward products (more moles of gas favored at lower pressure).",
            "backward": "Reaction shifted toward reactants (fewer moles of gas favored at higher pressure).",
            "no_shift": "No shift: equal moles of gas on both sides.",
        }

        logger.info(f"Pressure effect: {initial_total_p}→{new_total_p} atm, shift={shift}")
        return {
            "old_equilibrium": old_eq,
            "new_equilibrium": new_eq,
            "shift_direction": shift,
            "degree_of_dissociation_change": desc_map.get(shift, ""),
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            # Format: "P_old=1.0 P_new=5.0 Kp=0.067; A=-1 B=2"
            p_old = p_new = Kp = None
            stoich = {}
            sections = input_params.split(";")
            for sec in sections:
                sec = sec.strip()
                if sec.startswith("P_old=") or sec.startswith("p_old="):
                    p_old = float(sec.split("=")[1])
                elif sec.startswith("P_new=") or sec.startswith("p_new="):
                    p_new = float(sec.split("=")[1])
                elif sec.upper().startswith("KP=") or sec.upper().startswith("KP:"):
                    Kp = float(sec.split("=")[1])
                else:
                    for item in sec.split():
                        name, val = item.split("=")
                        stoich[name.strip()] = float(val)
            if None in (p_old, p_new, Kp):
                raise ValueError("Need P_old, P_new, and Kp.")
            return self._run_base(p_old, p_new, Kp, stoich)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
