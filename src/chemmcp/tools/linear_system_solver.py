import logging
import json
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LinearSystemSolver(BaseTool):
    """
    线性方程组求解器，用于化学平衡计算、质量守恒、最小二乘拟合等。
    支持 Ax=b 的精确求解和超定/欠定系统的最小二乘解。
    """
    __version__ = "0.1.0"
    name = "LinearSystemSolver"
    func_name = "solve_linear_system"
    description = "Solve linear systems Ax=b for chemical equilibrium calculations, mass conservation, stoichiometry, and least-squares fitting."
    implementation_description = "Uses numpy.linalg.solve for square systems, numpy.linalg.lstsq for overdetermined/underdetermined systems. Supports multiple right-hand sides."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Linear Algebra", "Equation Solver", "Chemical Equilibrium", "Mass Conservation"]
    required_envs = []

    code_input_sig = [
        ("matrix_a", "list", "N/A", "Coefficient matrix A (m x n) as list of lists."),
        ("vector_b", "list", "N/A", "Right-hand side vector b (length m) or matrix."),
        ("method", "str", "auto", "Solution method: 'auto', 'solve' (square), 'lstsq' (least squares)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string: '{\"matrix_a\": [[2,1],[1,3]], \"vector_b\": [5,6]}'"),
    ]

    output_sig = [
        ("solution", "list", "Solution vector x."),
        ("residual_norm", "float", "L2 norm of residual ||Ax - b||."),
        ("system_type", "str", "Type of system: 'square', 'overdetermined', or 'underdetermined'."),
        ("rank", "int", "Rank of coefficient matrix A."),
    ]

    examples = [
        {
            "code_input": {
                "matrix_a": [[3, 1], [1, 2], [4, 3]],
                "vector_b": [9, 8, 17],
                "method": "auto",
            },
            "text_input": {
                "input_str": '{"matrix_a": [[3,1],[1,2],[4,3]], \"vector_b\": [9,8,17]}',
            },
            "output": {
                "solution": [2.0, 3.0],
                "residual_norm": 0.0,
                "system_type": "overdetermined",
                "rank": 2,
            },
        },
        {
            "code_input": {
                "matrix_a": [[2, 1], [1, 3]],
                "vector_b": [5, 6],
                "method": "auto",
            },
            "text_input": {
                "input_str": '{"matrix_a": [[2,1],[1,3]], \"vector_b\": [5,6]}',
            },
            "output": {
                "solution": [1.8, 1.4],
                "residual_norm": 0.0,
                "system_type": "square",
                "rank": 2,
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, matrix_a: list, vector_b: list, method: str = "auto") -> dict:
        """Core logic: solve linear system."""
        A = np.array(matrix_a, dtype=float)
        b = np.array(vector_b, dtype=float)

        if A.ndim != 2:
            raise ChemMCPError("Coefficient matrix must be 2D.")
        if A.shape[0] == 0:
            raise ChemMCPError("Matrix cannot be empty.")
        if len(b.shape) == 1:
            b = b.reshape(-1, 1)
        if A.shape[0] != b.shape[0]:
            raise ChemMCPError(f"Dimension mismatch: A has {A.shape[0]} rows but b has {b.shape[0]} elements.")

        m, n = A.shape

        # Determine system type
        if m == n:
            system_type = "square"
        elif m > n:
            system_type = "overdetermined"
        else:
            system_type = "underdetermined"

        rank_A = int(np.linalg.matrix_rank(A))

        # Choose method
        if method == "auto":
            use_lstsq = (m != n)
        else:
            use_lstsq = (method == "lstsq")

        try:
            if not use_lstsq and m == n:
                x = np.linalg.solve(A, b)
                residual = float(np.linalg.norm(A @ x - b))
            else:
                x, residuals, rank_val, sv = np.linalg.lstsq(A, b, rcond=None)
                residual = float(np.sqrt(residuals[0])) if len(residuals) > 0 and residuals[0] > 0 else 0.0
        except np.linalg.LinAlgError as e:
            # Fallback to least squares
            x, residuals, rank_val, sv = np.linalg.lstsq(A, b, rcond=None)
            residual = float(np.sqrt(residuals[0])) if len(residuals) > 0 and residuals[0] > 0 else 0.0

        solution = x.flatten().tolist()

        logger.info(f"Solved {system_type} linear system ({m}x{n}), rank={rank_A}, residual={residual:.2e}")
        return {
            "solution": [round(v, 6) for v in solution],
            "residual_norm": round(residual, 6),
            "system_type": system_type,
            "rank": rank_A,
        }

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            matrix_a = params.get("matrix_a")
            vector_b = params.get("vector_b")
            method = params.get("method", "auto")
            return self._run_base(matrix_a, vector_b, method=method)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
