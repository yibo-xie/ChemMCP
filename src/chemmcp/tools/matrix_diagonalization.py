import logging
import json
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MatrixDiagonalization(BaseTool):
    """
    矩阵对角化工具，用于Hückel方法、Hartree-Fock计算核心等量子化学应用。
    返回对角矩阵 D 和变换矩阵 P，使得 A = P * D * P^(-1)。
    """
    __version__ = "0.1.0"
    name = "MatrixDiagonalization"
    func_name = "diagonalize_matrix"
    description = "Diagonalize a matrix for Hückel method, Hartree-Fock calculations, and quantum chemistry applications."
    implementation_description = "Uses numpy.linalg.eigh (symmetric) or numpy.linalg.eig (general) to compute P, D such that A = P @ D @ P_inv. Returns diagonal entries of D and transformation matrices."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Linear Algebra", "Diagonalization", "Quantum Chemistry", "Huckel Method", "Hartree-Fock"]
    required_envs = []

    code_input_sig = [
        ("matrix", "list", "N/A", "Square matrix as list of lists."),
        ("symmetric", "bool", "True", "Whether the matrix is symmetric/Hermitian."),
        ("return_pinv", "bool", "False", "Whether to also return the inverse of the transformation matrix."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string: '{\"matrix\": [[0,1,0],[1,0,1],[0,1,0]], \"symmetric\": true}'."),
    ]

    output_sig = [
        ("diagonal_values", "list", "Diagonal elements (eigenvalues) sorted descending by absolute value."),
        ("transformation_matrix", "list", "Matrix P whose columns are eigenvectors."),
        ("inverse_transformation", "list or null", "P^(-1) if requested, else null."),
        ("is_diagonalizable", "bool", "Whether the matrix is diagonalizable over reals."),
    ]

    examples = [
        {
            "code_input": {
                "matrix": [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
                "symmetric": True,
                "return_pinv": False,
            },
            "text_input": {
                "input_str": '{"matrix": [[0,1,0],[1,0,1],[0,1,0]], "symmetric": true}',
            },
            "output": {
                "diagonal_values": [1.414214, 0.0, -1.414214],
                "transformation_matrix": [[0.5, 0.707107, 0.5], [0.707107, 0.0, -0.707107], [0.5, -0.707107, 0.5]],
                "inverse_transformation": None,
                "is_diagonalizable": True,
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

    def _run_base(self, matrix: list, symmetric: bool = True, return_pinv: bool = False) -> dict:
        """Core logic: diagonalize matrix A → P D P^{-1}."""
        A = np.array(matrix, dtype=float)

        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ChemMCPError("Input must be a square matrix.")
        if A.shape[0] == 0:
            raise ChemMCPError("Matrix cannot be empty.")

        n = A.shape[0]

        try:
            if symmetric:
                eigenvalues, P = np.linalg.eigh(A)
            else:
                eigenvalues, P = np.linalg.eig(A)
        except np.linalg.LinAlgError as e:
            raise ChemMCPError(f"Diagonalization failed: {str(e)}")

        # Sort by absolute value descending
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        eigenvalues = eigenvalues[idx].real.astype(float)
        P = P[:, idx].real.astype(float)

        # Check diagonalizability: n linearly independent eigenvectors
        rank = int(np.linalg.matrix_rank(P))
        is_diag = (rank == n)

        result = {
            "diagonal_values": [round(ev, 6) for ev in eigenvalues.tolist()],
            "transformation_matrix": P.tolist(),
            "inverse_transformation": None,
            "is_diagonalizable": is_diag,
        }

        if return_pinv:
            try:
                P_inv = np.linalg.inv(P)
                result["inverse_transformation"] = P_inv.tolist()
            except np.linalg.LinAlgError:
                result["inverse_transformation"] = None

        logger.info(f"Diagonalized {n}x{n} matrix, diagonalizable={is_diag}")
        return result

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            matrix = params.get("matrix")
            symmetric = params.get("symmetric", True)
            return_pinv = params.get("return_pinv", False)
            return self._run_base(matrix, symmetric=symmetric, return_pinv=return_pinv)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
