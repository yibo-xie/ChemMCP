import logging
import json
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SvdDecomposition(BaseTool):
    """
    奇异值分解工具，用于化学计量学、光谱解析、主成分分析等。
    分解 A = U * Σ * V^T，返回奇异值和左右奇异向量。
    """
    __version__ = "0.1.0"
    name = "SvdDecomposition"
    func_name = "svd_decompose"
    description = "Singular Value Decomposition (SVD) for chemometrics, spectral analysis, principal component analysis, and data compression."
    implementation_description = "Uses numpy.linalg.svd to compute full SVD: A = U @ diag(S) @ Vt. Returns singular values (sorted descending), left/right singular vectors, and optional low-rank approximation."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Linear Algebra", "SVD", "Chemometrics", "PCA", "Spectral Analysis"]
    required_envs = []

    code_input_sig = [
        ("matrix", "list", "N/A", "Matrix (m x n) as list of lists."),
        ("full_matrices", "bool", "False", "Return full U/V matrices or economy-sized."),
        ("rank_approximation", "int", "0", "Low-rank approximation rank (0 = full rank)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string: '{\"matrix\": [[1,2,3],[4,5,6]], \"full_matrices\": false}'."),
    ]

    output_sig = [
        ("singular_values", "list", "Singular values sorted in descending order."),
        ("u_matrix", "list", "Left singular vectors (U)."),
        ("vt_matrix", "list", "Right singular vectors transposed (V^T)."),
        ("effective_rank", "int", "Numerical rank (number of singular values > threshold)."),
        ("reconstruction_error", "float or null", "Frobenius norm error if rank_approximation > 0."),
    ]

    examples = [
        {
            "code_input": {
                "matrix": [[1, 2], [3, 4], [5, 6]],
                "full_matrices": False,
                "rank_approximation": 0,
            },
            "text_input": {
                "input_str": '{"matrix": [[1,2],[3,4],[5,6]], "full_matrices": false}',
            },
            "output": {
                "singular_values": [9.525518, 0.514301],
                "u_matrix": [[-0.229847, 0.883452], [-0.524744, 0.240782], [-0.819642, -0.401889]],
                "vt_matrix": [[-0.619629, -0.78454], [-0.78454, 0.619629]],
                "effective_rank": 2,
                "reconstruction_error": None,
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

    def _run_base(self, matrix: list, full_matrices: bool = False, rank_approximation: int = 0) -> dict:
        """Core logic: compute SVD."""
        A = np.array(matrix, dtype=float)

        if A.ndim != 2:
            raise ChemMCPError("Input must be a 2D matrix.")
        if A.shape[0] == 0 or A.shape[1] == 0:
            raise ChemMCPError("Matrix cannot be empty.")

        try:
            U, s, Vt = np.linalg.svd(A, full_matrices=full_matrices)
        except np.linalg.LinAlgError as e:
            raise ChemMCPError(f"SVD computation failed: {str(e)}")

        # Effective rank: count singular values above threshold
        effective_rank = int(np.sum(s > max(s.max() * 1e-10, 1e-10)))

        result = {
            "singular_values": [round(sv, 6) for sv in s.tolist()],
            "u_matrix": [[round(v, 6) for v in row] for row in U.tolist()],
            "vt_matrix": [[round(v, 6) for v in row] for row in Vt.tolist()],
            "effective_rank": effective_rank,
            "reconstruction_error": None,
        }

        # Low-rank approximation
        if 0 < rank_approximation < len(s):
            k = min(rank_approximation, len(s))
            U_k = U[:, :k]
            s_k = s[:k]
            Vt_k = Vt[:k, :]
            A_approx = U_k @ np.diag(s_k) @ Vt_k
            error = float(np.linalg.norm(A - A_approx, 'fro'))
            result["reconstruction_error"] = round(error, 6)

        logger.info(f"SVD computed for {A.shape} matrix, effective_rank={effective_rank}")
        return result

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            matrix = params.get("matrix")
            full_mat = params.get("full_matrices", False)
            rank_appr = params.get("rank_approximation", 0)
            return self._run_base(matrix, full_matrices=full_mat, rank_approximation=rank_appr)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
