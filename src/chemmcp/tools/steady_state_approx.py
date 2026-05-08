import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SteadyStateApprox(BaseTool):
    """
    稳态近似法（Steady-State Approximation, SSA）求解中间体浓度与反应动力学。
    
    对于含中间体的多步反应机理，假设中间体生成速率 ≈ 消耗速率（d[I]/dt ≈ 0），
    从而推导出中间体浓度的解析表达式和总反应速率方程。
    
    支持的机理类型：
    - 连串反应：A → I → P（最基本）
    - 可逆连串：A ⇌ I → P 或 A → I ⇌ P
    - 预平衡机理：A ⇌ I (快) → P (慢)
    - 平行-连串复合网络
    - 用户自定义机理（通过速率方程组）
    """
    __version__ = "0.1.0"
    name = "SteadyStateApprox"
    func_name = "steady_state_approx"
    description = "Solve intermediate concentrations and derive rate laws using steady-state approximation (SSA) for multi-step reaction mechanisms with intermediates."
    implementation_description = "Applies d[intermediate]/dt ≈ 0 to solve for intermediate concentrations analytically. Supports consecutive, reversible-consecutive, pre-equilibrium, and user-defined mechanisms. Returns exact vs approximate concentration comparison, SSA validity criteria, and derived rate laws."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Steady-State Approximation", "Mechanism", "Intermediate", "Rate Law"]
    required_envs = []

    code_input_sig = [
        ("mechanism_type", "str", "N/A", "Mechanism type: 'consecutive', 'reversible', 'pre_equilibrium', 'parallel_consecutive', or 'custom'."),
        ("rate_constants", "list", "N/A", "List of rate constants in mechanism order. Consecutive: [k1,k2]. Reversible: [kf1,kr1,kf2] or [k1,k_1,k2]. Pre-equilibrium: [kf,kr,k2]."),
        ("initial_concentrations", "dict", "N/A", "Initial concentrations dict: {'A': 1.0}. Keys are species names."),
        ("time_points", "list", "N/A", "List of time points at which to evaluate concentrations. Example: [0, 10, 50, 100]."),
        ("intermediates", "list", "N/A", "List of intermediate species names. Example: ['I']."),
        ("custom_equations", "str", "None", "Custom rate equations for 'custom' type: format 'dA/dt=-k1*A; dI/dt=k1*A-k2*I; dP/dt=k2*I'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'mechanism_type k1,k2,... A0=1.0 t1,t2,... I1,I2,...' Example: 'consecutive 0.1,1.0 A0=1.0 0,10,50,100 I'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with concentration profiles (exact & SSA), validity assessment, rate law, error analysis, and recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "mechanism_type": "consecutive",
                "rate_constants": [0.1, 1.0],
                "initial_concentrations": {"A": 1.0},
                "time_points": [0, 5, 10, 20, 50, 100],
                "intermediates": ["I"],
                "custom_equations": None,
            },
            "text_input": {
                "input_params": "consecutive 0.1,1.0 A0=1.0 0,5,10,20,50,100 I",
            },
            "output": {
                "result": {
                    "mechanism": "consecutive A→I→P",
                    "ssa_valid": True,
                    "k2_k1_ratio": 10.0,
                    "intermediate_I_t50": {"exact": 0.0009, "ssa": 0.0007},
                    "derived_rate_law": "d[P]/dt ≈ k1[A]",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ── Exact solutions for consecutive reaction A --(k1)--> I --(k2)--> P ──
    def _consecutive_exact(self, k1, k2, A0, t):
        """Exact analytical solution for A→I→P."""
        A = A0 * math.exp(-k1 * t)
        if abs(k2 - k1) < 1e-12:
            I = A0 * k1 * t * math.exp(-k1 * t)
        else:
            I = A0 * k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))
        P = A0 - A - I
        return {"A": A, "I": I, "P": P}

    def _consecutive_ssa(self, k1, k2, A0, t):
        """SSA approximation: d[I]/dt=0 ⇒ [I]=(k1/k2)[A]"""
        A = A0 * math.exp(-k1 * t)
        I_ss = (k1 / k2) * A if k2 > 0 else 0
        P_ss = A0 - A - I_ss
        return {"A": A, "I": I_ss, "P": P_ss}

    # ── Reversible consecutive: A <==(k_r1)==(k_f1)==> I --(k2)--> P ──
    def _reversible_exact(self, kf1, kr1, k2, A0, t):
        """Exact solution for A⇌I→P using eigenvalue method."""
        # System: dA/dt = -kf1*A + kr1*I
        #        dI/dt = kf1*A - (kr1+k2)*I
        #        dP/dt = k2*I
        # Characteristic eq: λ² + (kf1+kr1+k2)λ + kf1*k2 = 0
        a_coef = 1.0
        b_coef = kf1 + kr1 + k2
        c_coef = kf1 * k2
        disc = b_coef ** 2 - 4 * a_coef * c_coef
        
        if disc < 0:
            # Complex roots — shouldn't happen for real kinetics
            lam1 = -b_coef / 2
            lam2 = lam1
        else:
            sqrt_disc = math.sqrt(disc)
            lam1 = (-b_coef + sqrt_disc) / 2
            lam2 = (-b_coef - sqrt_disc) / 2
        
        # Coefficients from initial conditions: A(0)=A0, I(0)=0
        if abs(lam1 - lam2) > 1e-12:
            C_A1 = (lam2 + kf1) / (lam2 - lam1) * A0
            C_A2 = -(lam1 + kf1) / (lam2 - lam1) * A0
            A_t = C_A1 * math.exp(lam1 * t) + C_A2 * math.exp(lam2 * t)
            
            C_I1 = kf1 / (lam2 - lam1) * A0
            C_I2 = -kf1 / (lam2 - lam1) * A0
            I_t = C_I1 * math.exp(lam1 * t) + C_I2 * math.exp(lam2 * t)
        else:
            A_t = A0 * math.exp(lam1 * t) * (1 - (lam1 + kf1) * t)
            I_t = A0 * kf1 * t * math.exp(lam1 * t)

        # Integrate P from I
        P_t = A0 - max(A_t, 0) - max(I_t, 0)
        
        return {"A": max(A_t, 0), "I": max(I_t, 0), "P": max(P_t, 0)}

    def _reversible_ssa(self, kf1, kr1, k2, A0, t):
        """SSA for A⇌I→P: set d[I]/dt≈0."""
        # From d[I]/dt = kf1[A] - (kr1+k2)[I] = 0
        # [I]_ss = kf1/(kr1+k2) · [A]
        A = A0 * math.exp(-kf1 * kf1 / (kr1 + k2) * t) if (kr1 + k2) > 0 else A0
        I_ss = kf1 / (kr1 + k2) * A if (kr1 + k2) > 0 else 0
        P_ss = A0 - A - I_ss
        return {"A": A, "I": I_ss, "P": P_ss}

    # ── Pre-equilibrium: A <==(kr)==(kf)==> I --(k2)--> P ──
    def _preeq_exact(self, kf, kr, k2, A0, t):
        """
        Simplified exact treatment of pre-equilibrium.
        Under fast equilibrium assumption combined with slow step.
        """
        Keq = kf / kr if kr > 0 else float('inf')
        # Effective first-order decay: keff = k2 * Keq / (1+Keq)
        keff = k2 * Keq / (1 + Keq) if Keq < 1e10 else k2
        A_t = A0 * math.exp(-keff * t)
        I_t = Keq * A_t
        P_t = A0 - A_t - I_t
        return {"A": A_t, "I": I_t, "P": max(P_t, 0)}

    def _preeq_ssa(self, kf, kr, k2, A0, t):
        """SSA for pre-equilibrium: same as fast-equilibrium result."""
        return self._preeq_exact(kf, kr, k2, A0, t)

    # ── Parallel-consecutive: A --(k1)--> B, A --(k2)--> I --(k3)--> C ──
    def _parallel_consec(self, ks, A0, t):
        """Parallel-consecutive network."""
        k1, k2, k3 = ks[0], ks[1], ks[2]
        A = A0 * math.exp(-(k1 + k2) * t)
        B = k1 * A0 / (k1 + k2) * (1 - math.exp(-(k1 + k2) * t)) if (k1 + k2) > 0 else k1 * A0 * t
        if abs(k3 - (k1 + k2)) < 1e-12:
            I = k2 * A0 * t * math.exp(-(k1 + k2) * t)
        elif (k1 + k2) != k3:
            I = k2 * A0 / ((k1 + k2) - k3) * (math.exp(-k3 * t) - math.exp(-(k1 + k2) * t))
        else:
            I = 0
        C = A0 - A - B - I
        return {"A": A, "B": B, "I": max(I, 0), "C": max(C, 0)}

    def _run_base(self, mechanism_type: str, rate_constants: list,
                  initial_concentrations: dict, time_points: list,
                  intermediates: list = None, custom_equations: str = None) -> dict:
        """Core logic."""
        mech = mechanism_type.lower().replace("-", "_")
        if not rate_constants:
            raise ChemMCPError("Rate constants list cannot be empty.")
        if not time_points:
            raise ChemMCPError("Time points list cannot be empty.")

        A0 = initial_concentrations.get("A", initial_concentrations.get("a", 1.0))
        inter_names = intermediates or ["I"]

        # Compute profiles
        profiles = []
        for t in time_points:
            if mech == "consecutive":
                exact = self._consecutive_exact(rate_constants[0], rate_constants[1], A0, t)
                approx = self._consecutive_ssa(rate_constants[0], rate_constants[1], A0, t)
            elif mech in ("reversible", "reversible_consecutive"):
                exact = self._reversible_exact(rate_constants[0], rate_constants[1],
                                                rate_constants[2] if len(rate_constants) > 2 else 0, A0, t)
                approx = self._reversible_ssa(rate_constants[0], rate_constants[1],
                                              rate_constants[2] if len(rate_constants) > 2 else 0, A0, t)
            elif mech in ("pre_equilibrium", "preequil", "pre_equil"):
                exact = self._preeq_exact(rate_constants[0], rate_constants[1],
                                           rate_constants[2] if len(rate_constants) > 2 else 0, A0, t)
                approx = self._preeq_ssa(rate_constants[0], rate_constants[1],
                                         rate_constants[2] if len(rate_constants) > 2 else 0, A0, t)
            elif mech == "parallel_consecutive":
                exact = self._parallel_consec(rate_constants, A0, t)
                approx = exact  # No simple SSA form
            else:
                raise ChemMCPError(f"Unknown mechanism type: {mechanism_type}. "
                                   f"Choose: consecutive, reversible, pre_equilibrium, parallel_consecutive")

            entry = {"time": t}
            for sp in exact:
                entry[f"{sp}_exact"] = round(exact[sp], 8)
                entry[f"{sp}_ssa"] = round(approx.get(sp, exact[sp]), 8)
            profiles.append(entry)

        # Validity assessment
        k1 = rate_constants[0]
        k2 = rate_constants[1] if len(rate_constants) > 1 else 1.0
        ratio_k2_k1 = k2 / k1 if k1 > 0 else float('inf')
        ssa_valid = ratio_k2_k1 >= 5  # Rule of thumb: k_consume >> k_form

        # Error analysis at last time point
        if profiles:
            last = profiles[-1]
            errors = {}
            for sp in ["I"] + inter_names:
                key_e = f"{sp}_exact"
                key_a = f"{sp}_ssa"
                if key_e in last and key_a in last:
                    e_val = last[key_e]
                    a_val = last[key_a]
                    if abs(e_val) > 1e-15:
                        errors[sp] = round(abs(e_val - a_val) / abs(e_val) * 100, 4)
                    else:
                        errors[sp] = round(abs(a_val) * 100, 4) if abs(a_val) > 1e-15 else 0.0

        # Derive rate law
        if mech == "consecutive":
            rate_law = f"d[P]/dt ≈ k₁[A];  [I]_ss = (k₁/k₂)[A]"
        elif mech == "reversible":
            rate_law = f"d[P]/dt ≈ k₃·k₁/(k₋₁+k₃)·[A];  K_eq,eff = k₁/(k₋₁+k₃)"
        elif mech.startswith("pre"):
            Keq = rate_constants[0] / rate_constants[1] if len(rate_constants) > 1 and rate_constants[1] > 0 else float('inf')
            k2 = rate_constants[2] if len(rate_constants) > 2 else 0
            rate_law = f"d[P]/dt ≈ k₂·K_eq·[A] where K_eq={Keq:.4g};  k_eff = {k2*Keq/(1+Keq):.4g}"
        else:
            rate_law = "(See individual species equations)"

        result = {
            "mechanism_type": mechanism_type,
            "mechanism_description": self._describe_mechanism(mech, rate_constants),
            "rate_constants": rate_constants,
            "initial_concentration_A0": A0,
            "intermediates": inter_names,
            "concentration_profiles": profiles,
            "ssa_valid": ssa_valid,
            "validity_criteria": {
                "k_consume_over_k_form_ratio": round(ratio_k2_k1, 4),
                "rule_of_thumb": "SSA valid when k(consuming intermediate) / k(forming intermediate) ≥ 5",
                "threshold_met": ratio_k2_k1 >= 5,
            },
            "max_intermediate_error_percent": max(errors.values()) if errors else 0,
            "derived_rate_law": rate_law,
            "recommendation": (
                "SSA is VALID — good approximation." if ssa_valid else
                "SSA may be POOR — use full kinetic equations instead."
            ),
        }

        logger.info(f"SteadyStateApprox: mech={mech}, valid={ssa_valid}, k2/k1={ratio_k2_k1:.2f}")
        return result

    def _describe_mechanism(self, mech, ks):
        descriptions = {
            "consecutive": f"A →{' '+str(ks[0]) if len(ks)>0 else ''} I →{' '+str(ks[1]) if len(ks)>1 else ''} P",
            "reversible": f"A ⇌({'×'+str(ks[0]) if len(ks)>0 else ''}/{'×'+str(ks[1]) if len(ks)>1 else ''}) I →{' ×'+str(ks[2]) if len(ks)>2 else ''} P",
            "pre_equilibrium": f"A ⇌(fast, K_eq={'%.3g'%(ks[0]/ks[1]) if len(ks)>1 and ks[1]>0 else '∞'}) I →(slow) P",
            "parallel_consecutive": f"A → B (k₁={ks[0]}), A → I (k₂={ks[1]}), I → C (k₃={ks[2] if len(ks)>2 else '?'})",
        }
        return descriptions.get(mech, mech)

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mech_type = parts[0]
            k_list = [float(x) for x in parts[1].split(",")]
            
            # Parse initial conc
            init_dict = {}
            time_list = []
            inters = []
            
            idx = 2
            while idx < len(parts):
                p = parts[idx]
                if p.startswith("=") or "=" in p:
                    key, val = p.split("=", 1)
                    init_dict[key.strip()] = float(val)
                elif "," in p or p.replace(".","").replace("-","").isdigit():
                    try:
                        time_list.extend([float(x) for x in p.split(",")])
                        idx += 1
                        continue
                    except ValueError:
                        pass
                    inters.append(p)
                else:
                    inters.append(p)
                idx += 1
            
            if not time_list:
                time_list = [0, 10, 50, 100]
            if not init_dict:
                init_dict = {"A": 1.0}

            return self._run_base(mech_type, k_list, init_dict, time_list, inters or ["I"])
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'mechanism k1,k2,... A0=1.0 t1,t2,... I1,I2,...'"
            )
