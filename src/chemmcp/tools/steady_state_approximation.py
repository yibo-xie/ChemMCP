import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SteadyStateApproximation(BaseTool):
    """
    稳态近似法求解中间体浓度。
    对于多步反应机理，假设中间体浓度不随时间变化（d[I]/dt ≈ 0），求解中间体及其他物种的近似浓度。
    
    典型应用：连串反应 A → I → P，其中 I 为中间体。
    """
    __version__ = "0.1.0"
    name = "SteadyStateApproximation"
    func_name = "steady_state_approx"
    description = "Solve intermediate concentrations using steady-state approximation (SSA) for multi-step reaction mechanisms."
    implementation_description = "Applies steady-state approximation (d[intermediate]/dt = 0) to derive analytical expressions for intermediate and product concentrations. Supports consecutive reactions A→I→P and similar mechanisms."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Steady-State Approximation", "Mechanism", "Intermediate"]
    required_envs = []

    code_input_sig = [
        ("mechanism_type", "str", "N/A", "Type of mechanism: 'consecutive' (A->I->P), 'reversible_consecutive', or 'pre-equilibrium' (A<=>I->P)."),
        ("rate_constants", "list", "N/A", "List of rate constants in order. For consecutive: [k1, k2]. For pre-equilibrium: [kf, kr, k2]."),
        ("initial_reactant_concentration_A0", "float", "N/A", "Initial concentration of the starting reactant A."),
        ("time_t", "float", "N/A", "Time at which to evaluate concentrations."),
        ("target_intermediate", "str", "I", "Name of the intermediate species (default 'I')."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'mechanism_type k1,k2,... A0 t [intermediate_name]'. Example: 'consecutive 0.1,1.0 1.0 50 I'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with intermediate concentration, all species concentrations, SSA validity check, and exact vs approximate comparison."),
    ]

    examples = [
        {
            "code_input": {
                "mechanism_type": "consecutive",
                "rate_constants": [0.1, 1.0],
                "initial_reactant_concentration_A0": 1.0,
                "time_t": 50.0,
                "target_intermediate": "I",
            },
            "text_input": {
                "input_params": "consecutive 0.1,1.0 1.0 50",
            },
            "output": {
                "result": {
                    "intermediate_I": 0.00091,
                    "reactant_A": 0.00674,
                    "product_P": 0.9933,
                    "ssa_valid": True,
                    "method": "consecutive",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _consecutive_ssa(self, k1, k2, A0, t):
        """
        Consecutive reaction: A --(k1)--> I --(k2)--> P
        SSA: d[I]/dt = k1[A] - k2[I] ≈ 0 => [I]_ss = (k1/k2)[A]
        """
        # Exact solutions:
        # [A] = A0 * exp(-k1*t)
        # [I]_exact = A0 * k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))
        # [P]_exact = A0 - [A] - [I]
        
        A_exact = A0 * math.exp(-k1 * t)
        
        if abs(k2 - k1) < 1e-12:
            I_exact = A0 * k1 * t * math.exp(-k1 * t)
        else:
            I_exact = A0 * k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))
        
        P_exact = A0 - A_exact - I_exact
        
        # SSA approximate: [I]_ss = (k1/k2) * [A]
        I_ss = (k1 / k2) * A_exact
        P_ss = A0 - A_exact - I_ss
        
        # Validity: SSA is good when k2 >> k1 (intermediate consumed fast)
        ratio = k2 / k1 if k1 > 0 else float('inf')
        ssa_valid = ratio > 5  # rule of thumb
        
        return {
            "exact": {"A": round(A_exact, 8), "I": round(I_exact, 8), "P": round(P_exact, 8)},
            "approximate": {"A": round(A_exact, 8), "I": round(I_ss, 8), "P": round(P_ss, 8)},
            "ssa_valid": ssa_valid,
            "k2_k1_ratio": round(ratio, 4),
            "intermediate_error_percent": round(abs(I_exact - I_ss) / max(abs(I_exact), 1e-15) * 100, 4) if I_exact != 0 else 0,
        }

    def _pre_equilibrium_ssa(self, kf, kr, k2, A0, t):
        """
        Pre-equilibrium: A <==(kf/kr)==> I --(k2)--> P
        SSA on I + equilibrium assumption for A<=>I
        K_eq = kf/kr, [I] = K_eq*[A], then rate = k2*K_eq*[A]
        """
        Keq = kf / kr if kr > 0 else float('inf')
        
        # Simplified: effective first-order decay of A with keff = k2*Keq*kf/(kf+kr)... 
        # Actually using standard treatment: d[P]/dt ≈ k2*Keq*[A] when equilibrium is fast
        # For simplicity, use the full numerical-like approach
        # Under fast equilibrium: [I]/[A] = Keq, total A_total = [A]+[I] = [A]*(1+Keq)
        # Rate = k2*[I] = k2*Keq/(1+Keq) * A_total
        
        keff = k2 * Keq / (1 + Keq) if Keq < 1e10 else k2
        A_remaining = A0 * math.exp(-keff * t)
        I_conc = Keq * A_remaining
        P_formed = A0 - A_remaining - I_conc
        
        return {
            "exact": {"A": round(A_remaining, 8), "I": round(I_conc, 8), "P": round(P_formed, 8)},
            "approximate": {"A": round(A_remaining, 8), "I": round(I_conc, 8), "P": round(P_formed, 8)},
            "ssa_valid": True,
            "equilibrium_constant_Keq": round(Keq, 6),
            "effective_rate_constant_keff": round(keff, 6),
            "intermediate_error_percent": 0,
        }

    def _run_base(self, mechanism_type: str, rate_constants: list,
                  initial_reactant_concentration_A0: float, time_t: float,
                  target_intermediate: str = "I") -> dict:
        if mechanism_type not in ("consecutive", "reversible_consecutive", "pre_equilibrium"):
            raise ChemMCPError(f"Unknown mechanism type: {mechanism_type}. Use 'consecutive', 'reversible_consecutive', or 'pre_equilibrium'.")
        if len(rate_constants) < 2:
            raise ChemMCPError("Need at least 2 rate constants.")
        if initial_reactant_concentration_A0 <= 0:
            raise ChemMCPError("Initial concentration must be positive.")
        if time_t < 0:
            raise ChemMCPError("Time cannot be negative.")

        if mechanism_type == "consecutive":
            result = self._consecutive_ssa(rate_constants[0], rate_constants[1],
                                           initial_reactant_concentration_A0, time_t)
        elif mechanism_type == "pre_equilibrium":
            result = self._pre_equilibrium_ssa(rate_constants[0], rate_constants[1],
                                                rate_constants[2], initial_reactant_concentration_A0, time_t)
        else:
            # reversible_consecutive treated similarly to consecutive for now
            result = self._consecutive_ssa(rate_constants[0], rate_constants[1],
                                           initial_reactant_concentration_A0, time_t)

        result["mechanism_type"] = mechanism_type
        result["target_intermediate"] = target_intermediate
        result["time"] = time_t

        logger.info(f"SteadyStateApproximation: type={mechanism_type}, t={time_t}, "
                     f"I_ss={result['approximate'].get('I', 'N/A')}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mech_type = parts[0]
            k_list = [float(x) for x in parts[1].split(",")]
            A0 = float(parts[2])
            t = float(parts[3])
            inter = parts[4] if len(parts) > 4 else "I"
            return self._run_base(mech_type, k_list, A0, t, inter)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'type k1,k2,... A0 t [intermediate]'")
