import logging
import json
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MatrixEigenvalueSolver(BaseTool):
    """
    矩阵特征值求解器，用于分子轨道能级、振动模式分析等化学应用。
    支持实对称/非对称矩阵的特征值和特征向量计算。
    """
    __version__ = "0.1.0"
    name = "MatrixEigenvalueSolver"
    func_name = "solve_eigenvalues"
    description = "Solve matrix eigenvalues and eigenvectors for molecular orbital energy levels, vibrational mode analysis, and quantum chemistry applications."
    implementation_description = "Uses numpy.linalg.eig for general matrices and numpy.linalg.eigh for symmetric (Hermitian) matrices. Returns eigenvalues (sorted) and corresponding eigenvectors."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Linear Algebra", "Eigenvalues", "Quantum Chemistry", "Molecular Orbitals", "Vibrational Modes"]
    required_envs = []

    code_input_sig = [
        ("matrix", "list", "N/A", "Square matrix as list of lists, e.g., [[1,2],[3,4]]."),
        ("symmetric", "bool", "True", "Whether the matrix is symmetric (use eigh for better numerical stability)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string of the matrix and options: '{\"matrix\": [[1,2],[3,4]], \"symmetric\": true}'."),
    ]

    output_sig = [
        ("eigenvalues", "list", "Sorted eigenvalues (descending by absolute value)."),
        ("eigenvectors", "list", "Corresponding eigenvectors (each as a list)."),
        ("matrix_size", "int", "Dimension of the input matrix (n x n)."),
    ]

    examples = [
        {
            "code_input": {
                "matrix": [[2, -1, 0], [-1, 2, -1], [0, -1, 2]],
                "symmetric": True,
            },
            "text_input": {
                "input_str": '{"matrix": [[2,-1,0],[-1,2,-1],[0,-1,2]], "symmetric": true}',
            },
            "output": {
                "eigenvalues": [3.414214, 2.0, 0.585786],
                "eigenvectors": [[-0.5, 0.707107, 0.5], [0.707107, 0.0, -0.707107], [-0.5, -0.707107, 0.5]],
                "matrix_size": 3,
            },
        },
        {
            "code_input": {
                "matrix": [[4, 1], [2, 3]],
                "symmetric": False,
            },
            "text_input": {
                "input_str": '{"matrix": [[4,1],[2,3]], "symmetric": false}',
            },
            "output": {
                "eigenvalues": [5.0, 2.0],
                "eigenvectors": [[0.707107, 0.447214], [0.707107, -0.894427]],
                "matrix_size": 2,
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

    def _run_base(self, matrix: list, symmetric: bool = True) -> dict:
        """Core logic: compute eigenvalues and eigenvectors."""
        A = np.array(matrix, dtype=float)

        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ChemMCPError("Input must be a square matrix (n x n).")
        if A.shape[0] == 0:
            raise ChemMCPError("Matrix cannot be empty.")

        n = A.shape[0]

        try:
            if symmetric:
                eigenvalues, eigenvectors = np.linalg.eigh(A)
            else:
                eigenvalues, eigenvectors = np.linalg.eig(A)
        except np.linalg.LinAlgError as e:
            raise ChemMCPError(f"Eigenvalue computation failed: {str(e)}")

        # Sort by absolute value descending
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        eigenvalues = eigenvalues[idx].real.astype(float)
        eigenvectors = eigenvectors[:, idx].real.astype(float)

        # Normalize eigenvectors to unit length
        for i in range(n):
            norm = np.linalg.norm(eigenvectors[:, i])
            if norm > 1e-10:
                eigenvectors[:, i] /= norm

        logger.info(f"Computed {n}x{n} matrix eigenvalues: {eigenvalues.tolist()}")
        return {
            "eigenvalues": [round(ev, 6) for ev in eigenvalues.tolist()],
            "eigenvectors": [row.tolist() for row in eigenvectors.T],
            "matrix_size": n,
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse JSON text input."""
        try:
            params = json.loads(input_str)
            matrix = params.get("matrix")
            symmetric = params.get("symmetric", True)
            return self._run_base(matrix, symmetric=symmetric)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
