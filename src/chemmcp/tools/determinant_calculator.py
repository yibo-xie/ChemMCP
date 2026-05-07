import logging
import json
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DeterminantCalculator(BaseTool):
    """
    行列式计算工具，用于久期方程求解、判断矩阵奇异性等。
    支持小矩阵的余子式展开和大矩阵的LU分解法。
    """
    __version__ = "0.1.0"
    name = "DeterminantCalculator"
    func_name = "compute_determinant"
    description = "Compute matrix determinant for secular equation solving, singularity checking, and quantum chemistry applications."
    implementation_description = "Uses numpy.linalg.det for general matrices. Also provides cofactor expansion details for small matrices (up to 3x3)."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Linear Algebra", "Determinant", "Secular Equation", "Matrix"]
    required_envs = []

    code_input_sig = [
        ("matrix", "list", "N/A", "Square matrix as list of lists."),
        ("show_details", "bool", "False", "Show cofactor expansion steps for matrices up to 3x3."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string: '{\"matrix\": [[1,2],[3,4]], \"show_details\": false}'."),
    ]

    output_sig = [
        ("determinant", "float", "The determinant value."),
        ("matrix_size", "int", "Dimension of the matrix."),
        ("is_singular", "bool", "Whether the matrix is singular (det ≈ 0)."),
        ("details", "str or null", "Cofactor expansion steps if requested and matrix ≤ 3x3."),
    ]

    examples = [
        {
            "code_input": {
                "matrix": [[1, 2], [3, 4]],
                "show_details": False,
            },
            "text_input": {
                "input_str": '{"matrix": [[1,2],[3,4]], "show_details": false}',
            },
            "output": {
                "determinant": -2.0,
                "matrix_size": 2,
                "is_singular": False,
                "details": None,
            },
        },
        {
            "code_input": {
                "matrix": [[6, 1, 1], [4, -2, 5], [2, 8, 7]],
                "show_details": True,
            },
            "text_input": {
                "input_str": '{"matrix": [[6,1,1],[4,-2,5],[2,8,7]], "show_details": true}',
            },
            "output": {
                "determinant": -306.0,
                "matrix_size": 3,
                "is_singular": False,
                "details": "det = 6*((-2)*7-5*8) - 1*(4*7-5*2) + 1*(4*8-(-2)*2) = 6*(-54) - 18 + 36 = -306",
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

    def _cofactor_2x2(self, A: np.ndarray) -> str:
        """Generate cofactor expansion string for 2x2 matrix."""
        det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
        return f"det = ({A[0,0]})*({A[1,1]}) - ({A[0,1]})*({A[1,0]}) = {det}"

    def _cofactor_3x3(self, A: np.ndarray) -> str:
        """Sarrus rule / cofactor expansion for 3x3."""
        a, b, c = A[0, 0], A[0, 1], A[0, 2]
        d, e, f = A[1, 0], A[1, 1], A[1, 2]
        g, h, i = A[2, 0], A[2, 1], A[2, 2]
        terms = [
            f"{a}*(({e})*({i})-({f})*({h}))",
            f"-{b}*(({d})*({i})-({f})*({g}))",
            f"+{c}*(({d})*({h})-({e})*({g}))",
        ]
        det_val = round(float(np.linalg.det(A)), 6)
        return f"det = {' '.join(terms)} = {det_val}"

    def _run_base(self, matrix: list, show_details: bool = False) -> dict:
        """Core logic: compute determinant."""
        A = np.array(matrix, dtype=float)

        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ChemMCPError("Input must be a square matrix.")
        if A.shape[0] == 0:
            raise ChemMCPError("Matrix cannot be empty.")

        n = A.shape[0]
        det = float(np.linalg.det(A))
        det_rounded = round(det, 6)
        is_singular = abs(det) < 1e-10

        details = None
        if show_details and n <= 3:
            if n == 1:
                details = f"det = {A[0, 0]}"
            elif n == 2:
                details = self._cofactor_2x2(A)
            elif n == 3:
                details = self._cofactor_3x3(A)

        logger.info(f"Computed determinant of {n}x{n} matrix: {det_rounded}")
        return {
            "determinant": det_rounded,
            "matrix_size": n,
            "is_singular": is_singular,
            "details": details,
        }

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            matrix = params.get("matrix")
            show_details = params.get("show_details", False)
            return self._run_base(matrix, show_details=show_details)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
