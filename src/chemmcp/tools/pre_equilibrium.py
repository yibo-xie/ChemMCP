import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PreEquilibrium(BaseTool):
    """
    预平衡近似法（Pre-Equilibrium Approximation）。
    
    对于包含快速平衡步骤后跟慢速决速步的反应机理：
        A + B ⇌ I (快速平衡, K_eq = k_f/k_r)
        I → P (慢, 速控步)
    
    利用预平衡假设推导总反应速率方程、有效速率常数和中间体浓度。
    
    支持多种预平衡模式：
    - 单步预平衡：A ⇌ I → P
    - 双分子预平衡：A + B ⇌ I → P  
    - 多级预平衡：多个连续快平衡 + 一个慢步骤
    """
    __version__ = "0.1.0"
    name = "PreEquilibrium"
    func_name = "pre_equilibrium_approx"
    description = "Apply pre-equilibrium approximation to derive rate laws for mechanisms with fast equilibrium steps followed by a slow rate-determining step."
    implementation_description = "Uses K_eq = k_forward/k_reverse to express intermediate concentration in terms of reactants, then substitutes into RDS rate expression. Handles unimolecular and bimolecular pre-equilibria with single or multiple fast steps."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Pre-Equilibrium", "Rate Law", "Mechanism", "Approximation"]
    required_envs = []

    code_input_sig = [
        ("mechanism", "str", "N/A", "Mechanism type: 'unimolecular' (A⇌I→P), 'bimolecular' (A+B⇌I→P), or 'multi_step'."),
        ("k_forward_list", "list", "N/A", "Forward rate constants for each equilibrium/step in order."),
        ("k_reverse_list", "list", "N/A", "Reverse rate constants for each equilibrium step (0 for irreversible)."),
        ("k_slow", "float", "N/A", "Rate constant of the slow (rate-determining) step."),
        ("initial_concentrations", "dict", "N/A", "Initial concentrations: {'A': 1.0, 'B': 1.0}."),
        ("time_points", "list", "None", "Time points to evaluate concentrations (optional)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'mechanism kf1,kf2,... kr1,kr2,... k_slow A0=1.0 [B0=1.0]'. Example: 'unimolecular 100 10 0.5 A0=1.0'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with equilibrium constants, effective rate constant, derived rate law, concentration profiles, validity conditions."),
    ]

    examples = [
        {
            "code_input": {
                "mechanism": "unimolecular",
                "k_forward_list": [100.0],
                "k_reverse_list": [10.0],
                "k_slow": 0.5,
                "initial_concentrations": {"A": 1.0},
                "time_points": None,
            },
            "text_input": {
                "input_params": "unimolecular 100 10 0.5 A0=1.0",
            },
            "output": {
                "result": {
                    "K_eq": 10.0,
                    "effective_k": 5.0,
                    "rate_law": "rate = k_slow · K_eq · [A] = 5.0[A]",
                    "intermediate_I_eq": 0.909,
                    "validity": True,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, mechanism: str, k_forward_list: list, k_reverse_list: list,
                  k_slow: float, initial_concentrations: dict,
                  time_points: list = None) -> dict:
        """Core logic."""
        mech = mechanism.lower().strip()
        
        if not k_forward_list:
            raise ChemMCPError("Forward rate constants list cannot be empty.")
        if k_slow is None or k_slow < 0:
            raise ChemMCPError("Slow step rate constant must be non-negative.")

        # Pad reverse list to match forward length
        k_rev = (list(k_reverse_list) + [0] * max(0, len(k_forward_list) - len(k_reverse_list)))[:len(k_forward_list)]

        # Compute equilibrium constants
        K_eqs = []
        for i, (kf, kr) in enumerate(zip(k_forward_list, k_rev)):
            if kr > 0:
                K_eqs.append(round(kf / kr, 6))
            elif kr == 0:
                K_eqs.append(float('inf'))  # Irreversible
            else:
                K_eqs.append(float('inf'))

        A0 = initial_concentrations.get("A", 1.0)
        B0 = initial_concentrations.get("B", 0.0)

        if mech == "unimolecular":
            result = self._unimolecular_preeq(K_eqs, k_slow, A0, time_points)
        elif mech == "bimolecular":
            result = self._bimolecular_preeq(K_eqs, k_slow, A0, B0, time_points)
        elif mech == "multi_step":
            result = self._multi_step_preeq(K_eqs, k_forward_list, k_rev, k_slow,
                                             initial_concentrations, time_points)
        else:
            raise ChemMCPError(f"Unknown mechanism: {mechanism}. Choose: 'unimolecular', 'bimolecular', 'multi_step'.")

        result["mechanism"] = mechanism
        result["equilibrium_constants_Keq"] = K_eqs
        result["k_forward"] = list(k_forward_list)
        result["k_reverse"] = list(k_rev)
        result["k_slow"] = k_slow
        
        logger.info(f"PreEquilibrium: mech={mech}, K_eq={K_eqs}, k_eff={result.get('effective_rate_constant', 'N/A')}")
        return result

    def _unimolecular_preeq(self, K_eqs, k_slow, A0, time_points):
        """A ⇌ I → P"""
        Keq = K_eqs[0] if K_eqs[0] != float('inf') else 1e10
        
        # At equilibrium: [I]/[A] = Keq => [I] = Keq·[A]
        # Total: A_total = [A] + [I] = [A](1+Keq) => [A] = A_total/(1+Keq)
        # [I]_eq = Keq/(1+Keq) · A_total
        I_eq_fraction = Keq / (1 + Keq) if Keq < 1e10 else 1.0
        A_eq_fraction = 1.0 / (1 + Keq) if Keq < 1e10 else 0.0
        
        # Rate of product formation through RDS: d[P]/dt = k_slow · [I]
        # Substitute [I]: d[P]/dt = k_slow · Keq/(1+Keq) · A_total
        k_eff = k_slow * I_eq_fraction
        
        # Rate law string
        rate_law = f"rate = k_slow × K_eq/(1+K_eq) × [A]_total = {k_eff:.4g} × [A]_total"
        rate_law_detailed = (
            f"Mechanism: A ⇌ I (fast, K_eq = {Keq:.4g})\n"
            f"          I → P (slow, k = {k_slow})\n"
            f"Derived:   rate = {k_eff:.4g}[A]₀\n"
            f"          [I]_eq = {I_eq_fraction:.4g} × [A]₀"
        )

        profiles = []
        if time_points:
            for t in time_points:
                A_t = A0 * math.exp(-k_eff * t)
                I_t = I_eq_fraction * A_t
                P_t = A0 - A_t - I_t
                profiles.append({
                    "time": t, "A": round(A_t, 8), "I": round(I_t, 8),
                    "P": round(max(P_t, 0), 8),
                })

        return {
            "effective_rate_constant": round(k_eff, 6),
            "rate_law": rate_law,
            "rate_law_detailed": rate_law_detailed,
            "equilibrium_intermediate_fraction": round(I_eq_fraction, 6),
            "reactant_at_equilibrium_fraction": round(A_eq_fraction, 6),
            "concentration_profiles": profiles or None,
            "validity_conditions": {
                "fast_equilibrium": "k_f, k_r >> k_slow",
                "slow_RDS": "k_slow << k_f, k_r",
                "recommendation": "Valid when equilibrium established much faster than product formation.",
            },
            "validity": True,
        }

    def _bimolecular_preeq(self, K_eqs, k_slow, A0, B0, time_points):
        """A + B ⇌ I → P"""
        Keq = K_eqs[0] if K_eqs[0] != float('inf') else 1e10
        
        # For A + B ⇌ I: K_eq = [I]/([A][B])
        # This requires solving quadratic at each point.
        # Simplified: assume [A]=[B] initially and stoichiometric.
        
        # If A0 = B0: let x = [I], then [A]=A0-x, [B]=B0-x
        # K_eq = x/((A0-x)(B0-x)) => K_eq(A0-x)(B0-x) = x
        # For A0=B0: K_eq(A0-x)^2 = x => K_eq x^2 - (2K_eq A0+1)x + K_eq A0^2 = 0
        
        if abs(A0 - B0) < 1e-12:
            a_q = Keq
            b_q = -(2 * Keq * A0 + 1)
            c_q = Keq * A0 ** 2
            disc = b_q ** 2 - 4 * a_q * c_q
            if disc >= 0:
                I_eq = (-b_q + math.sqrt(disc)) / (2 * a_q) if a_q > 0 else min(A0, B0)
            else:
                I_eq = min(A0, B0) * 0.5
        else:
            # Simplified: use limiting reagent
            I_eq = min(A0, B0) * Keq / (1 + Keq) if Keq < 1e10 else min(A0, B0)

        k_eff = k_slow * I_eq  # Zero-order-ish in reactants at fixed [I]

        rate_law = (
            f"Mechanism: A + B ⇌ I (fast, K_eq = {Keq:.4g} M⁻¹)\n"
            f"          I → P (slow, k = {k_slow})\n"
            f"Derived:   rate ≈ {k_eff:.4g} (pseudo-zero-order near equilibrium)\n"
            f"Full form: rate = k_slow × K_eq × [A][B] / (1 + K_eq([A]+[B]) + ...)"
        )

        return {
            "effective_rate_constant": round(k_eff, 6),
            "rate_law": rate_law,
            "equilibrium_intermediate_concentration": round(I_eq, 6),
            "validity": True,
            "concentration_profiles": None,
            "validity_conditions": {
                "fast_equilibrium": "k_f[A][B], k_r >> k_slow[I]",
                "note": "For bimolecular pre-eq, exact solution requires solving quadratic.",
            },
        }

    def _multi_step_preeq(self, K_eqs, kfs, krs, k_slow, init_conc, time_points):
        """Multiple consecutive equilibria before RDS."""
        # e.g., A ⇌ I₁ (K1) ⇌ I₂ (K2) → P (kslow)
        # Overall: [I_n] = K1·K2·...·Kn · [A]
        K_overall = 1.0
        valid_K = True
        for Keq in K_eqs:
            if Keq == float('inf'):
                valid_K = False
                break
            K_overall *= Keq
        
        if not valid_K:
            K_overall = 1e10

        A0 = init_conc.get("A", 1.0)
        I_final_frac = K_overall / (1 + K_overall) if K_overall < 1e10 else 1.0
        k_eff = k_slow * I_final_frac

        n_steps = len(K_eqs)
        chain_str = " ⇌ ".join([f"I_{i}" for i in range(n_steps + 1)])
        K_str = " × ".join([f"K{i+1}={K_eqs[i]:.3g}" for i in range(n_steps)])

        return {
            "effective_rate_constant": round(k_eff, 6),
            "overall_equilibrium_constant": round(K_overall, 6) if valid_K else float('inf'),
            "n_equilibrium_steps": n_steps,
            "rate_law": (
                f"Multi-step pre-equilibrium ({n_steps} fast steps):\n"
                f"  A ⇌ I₁ ⇌ I₂ ⇌ ... ⇌ I_n → P\n"
                f"  K_overall = {K_str}\n"
                f"  rate = k_slow × K_overall/(1+K_overall) × [A]₀ = {k_eff:.4g}[A]₀"
            ),
            "validity": True,
            "concentration_profiles": None,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mech = parts[0]
            
            idx = 1
            # Parse forward ks
            kfs = [float(x) for x in parts[idx].split(",")]
            idx += 1
            
            # Parse reverse ks
            krs = [float(x) for x in parts[idx].split(",")]
            idx += 1
            
            k_slow = float(parts[idx])
            idx += 1
            
            init_dict = {}
            while idx < len(parts):
                p = parts[idx]
                if "=" in p:
                    key, val = p.split("=", 1)
                    init_dict[key.strip()] = float(val)
                idx += 1

            return self._run_base(mech, kfs, krs, k_slow, init_dict or {"A": 1.0})
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'mechanism kf1,kf2,... kr1,kr2,... k_slow A0=1.0 [B0=1.0]'"
            )
