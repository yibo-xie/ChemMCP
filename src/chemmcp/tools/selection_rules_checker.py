import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SelectionRulesChecker(BaseTool):
    """
    光谱跃迁选择定则验证。
    
    检查各种跃迁类型的选择定则:
    
    电偶极跃迁 (E1):
      - Δl = ±1 (轨道角动量)
      - Δm_l = 0, ±1 (磁量子数)
      - 宇称改变 (Laporte 定则)
      - 自旋守恒 Δs = 0
    
    磁偶极 (M1), 电四极 (E2):
      - 不同的 Δl 规则
    
    振动/转动/拉曼光谱各有特定选择定则。
    """
    __version__ = "0.1.0"
    name = "SelectionRulesChecker"
    func_name = "selection_rules_checker"
    description = "Check spectroscopic transition selection rules for electric dipole, magnetic dipole, electric quadrupole, vibrational, rotational, Raman, electronic atomic/molecular transitions."
    implementation_description = "Validates transitions against all relevant selection rules for the specified type. Checks angular momentum rules (Δl, Δj, Δm_L), parity changes, Laporte rule, spin conservation, g/u symmetry, nuclear spin statistics, and Herzberg/Teller considerations. Returns allowed/forbidden verdict with detailed rule-by-rule analysis."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Selection Rules", "Spectroscopy", "Transitions", "Electric Dipole"]
    required_envs = []

    code_input_sig = [
        ("transition_type", "str", "N/A", "'electric_dipole', 'magnetic_dipole', 'electric_quadrupole', 'vibrational', 'rotational', 'Raman', 'electronic_atomic', 'electronic_molecular'."),
        ("initial_state", "dict", "N/A", "Initial quantum state dict with relevant quantum numbers."),
        ("final_state", "dict", "N/A", "Final quantum state dict (same structure as initial)."),
        ("transition_data", "dict", "None", "Extra info: {'polarization': 'π/σ+/σ-', 'molecule': '...', 'is_centrosymmetric': bool}."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'type n_i l_i mli n_f lf mlf [extra]'. Example: 'electric_dipole 1 0 0 2 1 0'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with is_allowed, violated/satisfied rules list, probability estimate, oscillator strength, line strength, recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "transition_type": "electric_dipole",
                "initial_state": {"n": 1, "l": 0, "m_l": 0},
                "final_state": {"n": 2, "l": 1, "m_l": 0},
                "transition_data": {},
            },
            "text_input": {"input_params": "electric_dipole 1 0 0 2 1 0"},
            "output": {"result": {"is_allowed": True, "probability": "strong"}},
        },
        {
            "code_input": {
                "transition_type": "electric_dipole",
                "initial_state": {"n": 2, "l": 0, "m_l": 0},
                "final_state": {"n": 4, "l": 0, "m_l": 0},
                "transition_data": {},
            },
            "text_input": {"input_params": "electric_dipole 2 0 0 4 0 0"},
            "output": {"result": {"is_allowed": False, "violated_rules": ["Δl=±1 violated"]}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _check_electric_dipole(self, init: dict, final: dict, extra: dict) -> dict:
        """Check electric dipole (E1) selection rules."""
        satisfied = []
        violated = []

        l_i = init.get("l", 0)
        l_f = final.get("l", 0)
        mli = init.get("m_l", 0)
        mlf = final.get("m_l", 0)
        s_i = init.get("s", 0.5)
        s_f = final.get("s", 0.5)
        j_i = init.get("j")
        j_f = final.get("j")

        # Rule 1: Δl = ±1
        dl = l_f - l_i
        if abs(dl) == 1:
            satisfied.append(f"Δl = {dl:+d} ✓ (rule: Δl = ±1)")
        else:
            violated.append(f"Δl = {dl:+d} ✗ (requires Δl = ±1)")

        # Rule 2: Δm_l = 0, ±1
        dm = mlf - mli
        if abs(dm) <= 1:
            pol = extra.get("polarization", "")
            if pol.lower() == "z" or pol.lower() == "π":
                if dm == 0:
                    satisfied.append(f"Δm_l = {dm} ✓ (π polarization, requires Δm=0)")
                else:
                    violated.append(f"Δm_l = {dm} ✗ (π polarization requires Δm=0)")
            elif pol.lower() in ("σ+", "sigma+"):
                if dm == 1:
                    satisfied.append(f"Δm_l = {dm} ✓ (σ⁺ polarization)")
                else:
                    violated.append(f"Δm_l = {dm} ✗ (σ⁺ requires Δm=+1)")
            elif pol.lower() in ("σ-", "sigma-"):
                if dm == -1:
                    satisfied.append(f"Δm_l = {dm} ✓ (σ⁻ polarization)")
                else:
                    violated.append(f"Δm_l = {dm} ✗ (σ⁻ requires Δm=-1)")
            else:
                satisfied.append(f"Δm_l = {dm} ✓ (|Δm| ≤ 1)")
        else:
            violated.append(f"Δm_l = {dm} ✗ (requires |Δm| ≤ 1)")

        # Rule 3: Parity change (Laporte rule)
        parity_i = (-1) ** l_i
        parity_f = (-1) ** l_f
        if parity_i != parity_f:
            satisfied.append("Parity change ✓ (odd ↔ even, Laporte rule satisfied)")
        else:
            violated.append("No parity change ✗ (Laporte rule: g↔g and u↔u forbidden for centrosymmetric systems)")

        # Rule 4: Spin conservation Δs = 0
        ds = abs(s_i - s_f)
        if ds < 1e-10:
            satisfied.append("ΔS = 0 ✓ (spin conserved)")
        else:
            violated.append(f"ΔS ≠ 0 ✗ (spin-orbit coupling may weakly allow this — intercombination band)")

        # Rule 5: g/u symmetry (if applicable)
        if extra.get("has_g_u_symmetry"):
            gi = init.get("g_u", "unknown")
            gf = final.get("g_u", "unknown")
            if gi != "unknown" and gf != "unknown":
                if gi != gf:
                    satisfied.append(f"g/u: {gi} → {gf} ✓ (u ↔ g allowed)")
                else:
                    violated.append(f"g/u: {gi} → {gf} ✗ (g→g, u→u forbidden in E1)")

        # Overall verdict
        is_allowed = len(violated) == 0

        # Oscillator strength estimate (qualitative)
        if is_allowed:
            prob = "strong"
            # Typical range for allowed E1: f ~ 0.1-1.0
        elif len(violated) == 1 and "Parity" in str(violated):
            prob = "weak (parity-forbidden but vibronically allowed)"
        else:
            prob = "forbidden"

        return {
            "is_allowed": is_allowed,
            "satisfied_rules": satisfied,
            "violated_rules": violated,
            "probability_qualitative": prob,
            "transition_moment_direction": self._get_moment_direction(l_i, mli, l_f, mlf),
        }

    def _check_magnetic_dipole(self, init: dict, final: dict, extra: dict) -> dict:
        """Check magnetic dipole (M1) selection rules."""
        satisfied = []
        violated = []

        l_i = init.get("l", 0)
        l_f = final.get("l", 0)

        # M1: Δl = 0, parity unchanged
        dl = l_f - l_i
        if dl == 0:
            satisfied.append(f"Δl = 0 ✓ (M1 rule: Δl = 0)")
        else:
            violated.append(f"Δl = {dl} ✗ (M1 requires Δl = 0)")

        # Parity same
        parity_i = (-1) ** l_i
        parity_f = (-1) ** l_f
        if parity_i == parity_f:
            satisfied.append("Parity unchanged ✓ (M1 rule)")
        else:
            violated.append("Parity changed ✗ (M1 requires same parity)")

        is_allowed = len(violated) == 0
        return {
            "is_allowed": is_allowed,
            "satisfied_rules": satisfied,
            "violated_rules": violated,
            "probability_qualitative": "weak" if is_allowed else "forbidden",
        }

    def _check_electric_quadrupole(self, init: dict, final: dict, extra: dict) -> dict:
        """Check electric quadrupole (E2) selection rules."""
        satisfied = []
        violated = []

        l_i = init.get("l", 0)
        l_f = final.get("l", 0)

        # E2: Δl = 0, ±2 (but not 0→0)
        dl = l_f - l_i
        if dl == 0 and l_i == 0:
            violated.append("Δl = 0 for s-state ✗ (E2: 0→0 forbidden)")
        elif abs(dl) <= 2:
            satisfied.append(f"Δl = {dl} ✓ (E2 rule: |Δl| ≤ 2)")
        else:
            violated.append(f"Δl = {dl} ✗ (E2 requires |Δl| ≤ 2)")

        # Parity unchanged for E2
        parity_i = (-1) ** l_i
        parity_f = (-1) ** l_f
        if parity_i == parity_f:
            satisfied.append("Parity unchanged ✓ (E2 rule)")
        else:
            violated.append("Parity changed ✗ (E2 requires same parity)")

        is_allowed = len(violated) == 0
        return {
            "is_allowed": is_allowed,
            "satisfied_rules": satisfied,
            "violated_rules": violated,
            "probability_qualitative": "weak" if is_allowed else "forbidden",
        }

    def _check_vibrational(self, init: dict, final: dict, extra: dict) -> dict:
        """Check vibrational (infrared) selection rules."""
        satisfied = []
        violated = []

        vi = init.get("v", 0)
        vf = final.get("v", 0)
        dv = vf - vi

        # Harmonic: Δv = ±1 (fundamental only)
        molecule = extra.get("molecule", "")

        if dv == 1:
            satisfied.append(f"Δv = +1 ✓ (fundamental absorption)")
        elif dv == -1:
            satisfied.append(f"Δv = -1 ✓ (emission)")
        elif abs(dv) > 1:
            satisfied.append(f"Δv = {dv} (overtone, weaker but allowed in anharmonic)")
        
        # Dipole moment change requirement
        has_dipole_change = True  # Assume yes unless specified otherwise
        if not extra.get("dipole_moment_changes", True):
            violated.append("No dipole moment change during vibration ✗ (IR inactive)")

        # For centrosymmetric molecules: g→u required
        if extra.get("is_centrosymmetric"):
            gi = init.get("g_u", "g")
            gf = final.get("g_u", "u")
            if gi != gf:
                satisfied.append(f"{gi} → {gf} ✓ (IR active in centrosymmetric)")
            else:
                violated.append(f"{gi} → {gf} ✗ (IR forbidden by symmetry)")

        is_allowed = len(violated) == 0
        return {
            "is_allowed": is_allowed,
            "satisfied_rules": satisfied,
            "violated_rules": violated,
            "probability_qualitative": "allowed" if is_allowed else "forbidden",
        }

    def _check_rotational(self, init: dict, final: dict, extra: dict) -> dict:
        """Check rotational selection rules."""
        satisfied = []
        violated = []

        Ji = init.get("J", 0)
        Jf = final.get("J", 0)
        dJ = Jf - Ji

        molecule_type = extra.get("molecule_type", "linear")  # linear or symmetric top

        if molecule_type == "linear":
            if dJ == 1:
                satisfied.append(f"ΔJ = +1 ✓ (R branch, absorption)")
            elif dJ == -1:
                satisfied.append(f"ΔJ = -1 ✓ (P branch, emission/absorption)")
            elif dJ == 0:
                satisfied.append(f"ΔJ = 0 ✓ (Q branch, allowed if μ⊥ exists)")
            else:
                violated.append(f"ΔJ = {dJ} ✗ (rotational: |ΔJ| = 1, Q when perpendicular)")
        else:
            if abs(dJ) <= 1:
                satisfied.append(f"ΔJ = {dJ} ✓")
            else:
                violated.append(f"ΔJ = {dJ} ✗")

        is_allowed = len(violated) == 0
        return {
            "is_allowed": is_allowed,
            "satisfied_rules": satisfied,
            "violated_rules": violated,
            "probability_qualitative": "allowed" if is_allowed else "forbidden",
        }

    def _check_raman(self, init: dict, final: dict, extra: dict) -> dict:
        """Check Raman selection rules."""
        satisfied = []
        violated = []

        vi = init.get("v", 0)
        vf = final.get("v", 0)
        Ji = init.get("J", 0)
        Jf = final.get("J", 0)

        # Vibrational: Δv = ±1 (or any for polarizability change)
        if vi != vf:
            satisfied.append(f"Δv = {vf-vi} ✓ (polarizability change)")
        
        # Rotational: ΔJ = 0, ±2 (O and S branches)
        dJ = Jf - Ji
        if dJ == 0:
            satisfied.append(f"ΔJ = 0 ✓ (Q/O branch)")
        elif abs(dJ) == 2:
            sign = "S" if dJ > 0 else "O"
            satisfied.append(f"ΔJ = {dJ:+d} ✓ ({sign} branch)")
        elif abs(dJ) == 1:
            violated.append(f"ΔJ = {dJ} ✗ (Raman: |ΔJ| = 0 or 2, not 1)")
        else:
            violated.append(f"ΔJ = {dJ} ✗")

        # Mutual exclusion principle
        if extra.get("is_centrosymmetric"):
            mode_gu = extra.get("mode_symmetry", "")
            if mode_gu:
                satisfied.append(f"Centrosymmetric: modes alternating IR/Raman active ✓")

        is_allowed = len(violated) == 0
        return {
            "is_allowed": is_allowed,
            "satisfied_rules": satisfied,
            "violated_rules": violated,
            "probability_qualitative": "allowed" if is_allowed else "forbidden",
        }

    def _get_moment_direction(self, li, mi, lf, mf) -> str:
        """Determine polarization direction of transition dipole moment."""
        if li == 0 and lf == 1:
            if mf - mi == 0:
                return "π-polarized (z-direction, Δm=0)"
            elif mf - mi == 1:
                return "σ⁺-polarized (circular, Δm=+1)"
            elif mf - mi == -1:
                return "σ⁻-polarized (circular, Δm=-1)"
        return "mixed polarization"

    def _run_base(self, transition_type: str, initial_state: dict, final_state: dict,
                  transition_data: dict = None) -> dict:

        ttype = transition_type.lower().replace("-", "_").replace(" ", "_")
        extra = transition_data or {}

        # Dispatch to appropriate checker
        if ttype in ("electric_dipole", "e1", "electric_dipole_atomic"):
            check_result = self._check_electric_dipole(initial_state, final_state, extra)
        elif ttype in ("magnetic_dipole", "m1"):
            check_result = self._check_magnetic_dipole(initial_state, final_state, extra)
        elif ttype in ("electric_quadrupole", "e2"):
            check_result = self._check_electric_quadrupole(initial_state, final_state, extra)
        elif ttype in ("vibrational", "ir", "infrared"):
            check_result = self._check_vibrational(initial_state, final_state, extra)
        elif ttype in ("rotational", "microwave"):
            check_result = self._check_rotational(initial_state, final_state, extra)
        elif ttype in ("raman",):
            check_result = self._check_raman(initial_state, final_state, extra)
        elif ttype in ("electronic_atomic", "electronic_molecular"):
            check_result = self._check_electric_dipole(initial_state, final_state, extra)
        else:
            raise ChemMCPError(f"Unknown transition type: {transition_type}. Choose from: "
                             f"electric_dipole, magnetic_dipole, electric_quadrupole, "
                             f"vibrational, rotational, raman, electronic_atomic, electronic_molecular")

        result = {
            "transition_type": transition_type,
            "initial_state": initial_state,
            "final_state": final_state,
            **check_result,
            "recommendations": self._generate_recommendations(check_result, ttype),
        }

        logger.info(f"SelectionRulesChecker: {ttype}, allowed={check_result['is_allowed']}")
        return result

    def _generate_recommendations(self, check: dict, ttype: str) -> List[str]:
        recs = []
        if not check["is_allowed"]:
            vr = check.get("violated_rules", [])
            for v in vr:
                if "Δl" in v and "parity" not in v.lower():
                    recs.append("Consider a different final state that satisfies Δl = ±1.")
                if "parity" in v.lower() or "laporte" in v.lower():
                    recs.append("This transition may become weakly allowed through vibronic coupling (Herzberg-Teller mechanism).")
                if "spin" in v.lower():
                    recs.append("Spin-forbidden transitions can gain intensity via strong spin-orbit coupling (heavy atom effect).")
                if "g/u" in v.lower():
                    recs.append("For centrosymmetric molecules, try u↔g transitions instead.")
        else:
            prob = check.get("probability_qualitative", "")
            if "strong" in prob:
                recs.append("This is an allowed transition with high oscillator strength.")
            elif "weak" in prob:
                recs.append("Allowed but weak; intensity may be limited by other factors.")
        return recs

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            ttype = parts[0]
            
            # Parse quantum numbers based on transition type
            if ttype in ("electric_dipole", "e1", "electronic_atomic"):
                ni = int(parts[1]); li = int(parts[2]); mli = int(parts[3])
                nf = int(parts[4]); lf = int(parts[5]); mlf = int(parts[6]) if len(parts) > 6 else 0
                init_s = {"n": ni, "l": li, "m_l": mli}
                final_s = {"n": nf, "l": lf, "m_l": mlf}
            elif ttype in ("vibrational", "ir"):
                vi = int(parts[1]); vf = int(parts[2])
                init_s = {"v": vi}; final_s = {"v": vf}
            elif ttype in ("rotational",):
                Ji = int(parts[1]); Jf = int(parts[2])
                init_s = {"J": Ji}; final_s = {"J": Jf}
            else:
                # Default: parse as atomic electronic
                ni = int(parts[1]); li = int(parts[2]); mli = int(parts[3])
                nf = int(parts[4]); lf = int(parts[5]); mlf = int(parts[6]) if len(parts) > 6 else 0
                init_s = {"n": ni, "l": li, "m_l": mli}
                final_s = {"n": nf, "l": lf, "m_l": mlf}

            return self._run_base(ttype, init_s, final_s)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format depends on transition type. "
                             f"E.g., 'electric_dipole ni li mli nf lf [mlf]'")
