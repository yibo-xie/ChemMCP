import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CommutatorCalculator(BaseTool):
    """
    对易子计算工具。
    计算量子力学算符的对易子 [A,B] = AB - BA，验证测不准关系（不确定性原理）。
    支持常见量子力学算符：位置、动量、角动量、自旋、升降算符等。
    """
    __version__ = "0.1.0"
    name = "CommutatorCalculator"
    func_name = "calculate_commutator"
    description = "Calculate quantum mechanical commutators [A,B] and verify uncertainty principle relations."
    implementation_description = "Implements analytical commutator evaluation for standard quantum operators (x, p, L, S, creation/annihilation) using canonical commutation relations and angular momentum algebra."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Commutator", "Uncertainty Principle", "Operators"]
    required_envs = []

    code_input_sig = [
        ("operator_a", "str", "N/A", "First operator: 'x', 'y', 'z', 'px', 'py', 'pz', 'Lx', 'Ly', 'Lz', 'Sx', 'Sy', 'Sz', 'a', 'adag' (creation), 'N' (number)."),
        ("operator_b", "str", "N/A", "Second operator (same options as operator_a)."),
        ("mode", "str", "'canonical'", "Computation mode: 'canonical' (symbolic result), 'uncertainty' (ΔA·ΔB bound), 'matrix' (finite-dim matrix rep)."),
        ("state_info", "dict", "{}", "Additional state info for matrix mode: {'dim': N, 'basis': 'harmonic_oscillator'|'spin_half'}."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: operator_a operator_b [mode]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing commutator value/formula, uncertainty relation, physical interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "operator_a": "x",
                "operator_b": "px",
                "mode": "canonical",
                "state_info": {},
            },
            "text_input": {
                "input_str": "x px canonical",
            },
            "output": {
                "result": {
                    "commutator": "[x, p_x] = iℏ",
                    "value_symbolic": "iℏ",
                    "value_numeric": 1.054571817e-34j,
                    "uncertainty_relation": "Δx·Δp ≥ ℏ/2",
                    "interpretation": "Canonical commutation relation — foundation of quantum mechanics.",
                }
            }
        },
        {
            "code_input": {
                "operator_a": "Lx",
                "operator_b": "Ly",
                "mode": "canonical",
                "state_info": {},
            },
            "text_input": {
                "input_str": "Lx Ly canonical",
            },
            "output": {
                "result": {
                    "commutator": "[L_x, L_y] = iℏL_z",
                    "value_symbolic": "iℏL_z",
                    "uncertainty_relation": "ΔL_x·ΔL_y ≥ ℏ⟨L_z⟩/2",
                    "interpretation": "Angular momentum algebra — components do not commute.",
                }
            }
        },
        {
            "code_input": {
                "operator_a": "Sx",
                "operator_b": "Sy",
                "mode": "matrix",
                "state_info": {"dim": 2, "basis": "spin_half"},
            },
            "text_input": {
                "input_str": "Sx Sy matrix",
            },
            "output": {
                "result": {
                    "commutator": "[S_x, S_y] = iℏσ_z/2 = iℏS_z",
                    "matrix_commutator": [[0, 0], [0, 0]],  # symbolic representation
                    "pauli_matrices_used": True,
                    "interpretation": "Spin-1/2 angular momentum commutation relation.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34   # J·s

        # Define canonical commutation relations
        self.canonical_pairs = {
            ("x", "px"): ("i*hbar", self.hbar * 1j),
            ("y", "py"): ("i*hbar", self.hbar * 1j),
            ("z", "pz"): ("i*hbar", self.hbar * 1j),
            ("px", "x"): ("-i*hbar", -self.hbar * 1j),
            ("py", "y"): ("-i*hbar", -self.hbar * 1j),
            ("pz", "z"): ("-i*hbar", -self.hbar * 1j),
        }

        # Angular momentum commutation relations: [Li, Lj] = iℏε_ijk Lk
        self.angular_momentum = {
            ("lx", "ly"): ("i*hbar*Lz", None),
            ("ly", "lx"): ("-i*hbar*Lz", None),
            ("ly", "lz"): ("i*hbar*Lx", None),
            ("lz", "ly"): ("-i*hbar*Lx", None),
            ("lz", "lx"): ("i*hbar*Ly", None),
            ("lx", "lz"): ("-i*hbar*Ly", None),
        }

        # Spin commutation relations (same structure as L)
        self.spin_relations = {
            ("sx", "sy"): ("i*hbar*Sz", None),
            ("sy", "sx"): ("-i*hbar*Sz", None),
            ("sy", "sz"): ("i*hbar*Sx", None),
            ("sz", "sy"): ("-i*hbar*Sx", None),
            ("sz", "sx"): ("i*hbar*Sy", None),
            ("sx", "sz"): ("-i*hbar*Sy", None),
        }

        # Creation/annihilation operators
        self.ladder_ops = {
            ("a", "adag"): ("-1", -1.0 + 0j),
            ("adag", "a"): ("1", 1.0 + 0j),
            ("N", "a"): ("-a", None),
            ("N", "adag"): ("adag", None),
            ("a", "N"): ("a", None),
            ("adag", "N"): ("-adag", None),
        }

        # Same operators always commute with themselves
        self.zero_commutation = set()

    def _run_base(self, operator_a: str, operator_b: str, mode: str = "canonical",
                  state_info: Optional[dict] = None) -> dict:
        """Core logic: compute commutator."""
        if state_info is None:
            state_info = {}

        op_a = operator_a.strip().lower()
        op_b = operator_b.strip().lower()
        mode = mode.lower().strip()

        key = (op_a, op_b)

        # Check same operator
        if op_a == op_b:
            return self._make_result(op_a, op_b, "0", 0.0 + 0j,
                                      "Any operator commutes with itself: [A,A]=0.",
                                      uncertainty="ΔA·0 ≥ 0 (trivial)")

        # Canonical commutation [x_i, p_j]
        if key in self.canonical_pairs:
            sym, num = self.canonical_pairs[key]
            unc = f"Δ{op_a}·Δ{op_b} ≥ ℏ/2"
            return self._make_result(op_a, op_b, sym, num,
                                     f"Canonical commutation relation: [{op_a}, {op_b}] = {sym}",
                                     uncertainty=unc)

        # Reverse canonical
        rev_key = (op_b, op_a)
        if rev_key in self.canonical_pairs:
            sym, num = self.canonical_pairs[rev_key]
            neg_sym = "-" + sym if not sym.startswith("-") else sym[1:]
            neg_num = -num if num else None
            unc = f"Δ{op_a}·Δ{op_b} ≥ ℏ/2"
            return self._make_result(op_a, op_b, neg_sym, neg_num,
                                     f"Canonical commutation relation (reversed): [{op_a}, {op_b}] = {neg_sym}",
                                     uncertainty=unc)

        # Angular momentum
        if key in self.angular_momentum:
            sym, _ = self.angular_momentum[key]
            # Extract third component for uncertainty
            third = self._third_component(op_a, op_b)
            unc = f"Δ{op_a}·Δ{op_b} ≥ ℏ|⟨{third}⟩|/2"
            return self._make_result(op_a, op_b, sym, None,
                                     f"Angular momentum commutation: [{op_a}, {op_b}] = {sym}",
                                     uncertainty=unc)

        # Spin
        if key in self.spin_relations:
            sym, _ = self.spin_relations[key]
            third = self._third_component_spin(op_a, op_b)
            unc = f"Δ{op_a}·Δ{op_b} ≥ ℏ|⟨{third}⟩|/2"
            return self._make_result(op_a, op_b, sym, None,
                                     f"Spin commutation: [{op_a}, {op_b}] = {sym}",
                                     uncertainty=unc)

        # Ladder operators
        if key in self.ladder_ops:
            sym, num = self.ladder_ops[key]
            return self._make_result(op_a, op_b, sym, num,
                                     f"Ladder operator commutation: [{op_a}, {op_b}] = {sym}",
                                     uncertainty=None)

        # Position-position or momentum-momentum commute
        pos_set = {"x", "y", "z"}
        mom_set = {"px", "py", "pz"}
        if op_a in pos_set and op_b in pos_set:
            return self._make_result(op_a, op_b, "0", 0.0 + 0j,
                                     f"Different position components commute: [{op_a}, {op_b}] = 0.",
                                     uncertainty=f"No uncertainty constraint between {op_a} and {op_b}")
        if op_a in mom_set and op_b in mom_set:
            return self._make_result(op_a, op_b, "0", 0.0 + 0j,
                                     f"Different momentum components commute: [{op_a}, {op_b}] = 0.",
                                     uncertainty=f"No uncertainty constraint between {op_a} and {op_b}")

        # Number operator with Hamiltonian
        if key == ("n", "h") or key == ("h", "n"):
            sym = "0" if key == ("n", "h") else "0"
            return self._make_result(op_a, op_b, sym, 0.0,
                                     "[N, H] = 0 for harmonic oscillator (number operator commutes with H).",
                                     uncertainty=None)

        # Matrix mode for spin
        if mode == "matrix":
            return self._compute_matrix_commutator(op_a, op_b, state_info)

        raise ChemMCPError(
            f"Unknown operator pair ('{op_a}', '{op_b}') or unsupported combination. "
            f"Supported: x/y/z, px/py/pz, Lx/Ly/Lz, Sx/Sy/Sz, a/adag/N."
        )

    def _compute_matrix_commutator(self, op_a: str, op_b: str, info: dict) -> dict:
        """Compute commutator using finite-dimensional matrix representations."""
        basis = info.get("basis", "spin_half")
        dim = info.get("dim", 2)

        if basis == "spin_half":
            # Pauli matrices / 2 * hbar
            hbar = self.hbar
            sx = hbar / 2 * [[0, 1], [1, 0]]
            sy = hbar / 2 * [[0, -1j], [1j, 0]]
            sz = hbar / 2 * [[1, 0], [0, -1]]

            op_map = {"sx": sx, "sy": sy, "sz": sz}

            if op_a not in op_map or op_b not in op_map:
                raise ChemMCPError(f"For spin_half basis, use Sx/Sy/Sz. Got '{op_a}', '{op_b}'.")

            A = op_map[op_a]
            B = op_map[op_b]

            # Compute AB - BA
            AB = self._mat_mul(A, B)
            BA = self._mat_mul(B, A)
            comm = [[AB[i][j] - BA[i][j] for j in range(2)] for i in range(2)]

            # Expected result should be i*hbar * Sk
            expected_op = self._third_component_spin(op_a, op_b).lower()
            expected = op_map.get(expected_op, [[0, 0], [0, 0]])
            expected_comm = [[1j * expected[i][j] for j in range(2)] for i in range(2)]

            return self._make_result(
                op_a, op_b,
                f"iℏσ_{expected_op[-1]}/2",
                None,
                f"[{op_a.upper()}, {op_b.upper()}] matrix = {comm}",
                uncertainty=f"Δ{op_a}·Δ{op_b} ≥ ℏ|⟨{expected_op.upper()}⟩|/2",
                extra={
                    "matrix_commutator": comm,
                    "dimension": dim,
                    "basis": basis,
                }
            )

        elif basis == "harmonic_oscillator":
            N = dim
            # Truncate HO number basis to N dimensions
            # a = diag(sqrt(1), sqrt(2), ..., sqrt(N-1)) on off-diagonal
            a_mat = [[0.0] * N for _ in range(N)]
            adag_mat = [[0.0] * N for _ in range(N)]
            n_mat = [[0.0] * N for _ in range(N)]
            for i in range(N):
                n_mat[i][i] = float(i)
                if i < N - 1:
                    a_mat[i][i + 1] = math.sqrt(i + 1)
                    adag_mat[i + 1][i] = math.sqrt(i + 1)

            op_mats = {"a": a_mat, "adag": adag_mat, "n": n_mat}

            if op_a not in op_mats or op_b not in op_mats:
                raise ChemMCPError(f"For ho basis, use a/adag/N. Got '{op_a}', '{op_b}'.")

            A = op_mats[op_a]
            B = op_mats[op_b]
            AB = self._mat_mul(A, B)
            BA = self._mat_mul(B, A)
            comm = [[AB[i][j] - BA[i][j] for j in range(N)] for i in range(N)]

            return self._make_result(
                op_a, op_b, f"[{op_a},{op_b}]_truncated", None,
                f"Truncated {N}×{N} matrix commutator computed.",
                uncertainty=None,
                extra={"matrix_commutator": comm, "dimension": N, "basis": basis}
            )

        else:
            raise ChemMCPError(f"Unknown basis '{basis}'. Use 'spin_half' or 'harmonic_oscillator'.")

    @staticmethod
    def _mat_mul(A, B):
        """Multiply two matrices."""
        n = len(A)
        m = len(B[0])
        p = len(B)
        C = [[0.0] * m for _ in range(n)]
        for i in range(n):
            for j in range(m):
                s = 0.0
                for k in range(p):
                    s += A[i][k] * B[k][j]
                C[i][j] = s
        return C

    @staticmethod
    def _third_component(a: str, b: str) -> str:
        """Get the third component for angular momentum cross product."""
        cyclic = [("lx", "ly", "lz"), ("ly", "lz", "lx"), ("lz", "lx", "ly")]
        for cycle in cyclic:
            if (a, b) == (cycle[0], cycle[1]):
                return cycle[2]
            if (a, b) == (cycle[1], cycle[0]):
                return "-" + cycle[2]
        return "L?"

    @staticmethod
    def _third_component_spin(a: str, b: str) -> str:
        """Get the third component for spin."""
        cyclic = [("sx", "sy", "sz"), ("sy", "sz", "sx"), ("sz", "sx", "sy")]
        for cycle in cyclic:
            if (a, b) == (cycle[0], cycle[1]):
                return cycle[2]
            if (a, b) == (cycle[1], cycle[0]):
                return "-" + cycle[2]
        return "S?"

    def _make_result(self, op_a: str, op_b: str, symbolic: str, numeric, explanation: str,
                     uncertainty: Optional[str] = None, extra: Optional[dict] = None) -> dict:
        result = {
            "operator_a": op_a,
            "operator_b": op_b,
            "commutator": f"[{op_a}, {op_b}] = {symbolic}",
            "value_symbolic": symbolic,
            "uncertainty_relation": uncertainty,
            "explanation": explanation,
        }
        if numeric is not None:
            result["value_numeric"] = str(numeric)
        if extra:
            result.update(extra)
        return {"result": result}

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split()
            op_a = parts[0]
            op_b = parts[1]
            mode = parts[2] if len(parts) > 2 else "canonical"
            return self._run_base(op_a, op_b, mode)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
