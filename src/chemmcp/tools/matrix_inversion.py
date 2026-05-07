import logging
import json
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MatrixInversion(BaseTool):
    """
    矩阵求逆工具，用于最小二乘拟合、正规方程求解等化学计算。
    支持普通逆矩阵和伪逆（Moore-Penrose）。
    """
    __version__ = "0.1.0"
    name = "MatrixInversion"
    func_name = "invert_matrix"
    description = "Compute matrix inverse or pseudo-inverse for least squares fitting, normal equation solving, and regression analysis."
    implementation_description = "Uses numpy.linalg.inv for regular inverse and numpy.linalg.pinv for Moore-Penrose pseudo-inverse (handles singular/rectangular matrices)."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["Linear Algebra", "Matrix Inverse", "Least Squares", "Regression"]
    required_envs = []

    code_input_sig = [
        ("matrix", "list", "N/A", "Square matrix as list of lists."),
        ("use_pseudo_inverse", "bool", "False", "Use Moore-Penrose pseudo-inverse (for singular/ill-conditioned matrices)."),
        ("condition_number", "bool", "False", "Also compute and return condition number."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string: '{\"matrix\": [[4,7],[2,6]], \"use_pseudo_inverse\": false}'."),
    ]

    output_sig = [
        ("inverse", "list", "The inverted matrix."),
        ("is_singular", "bool", "Whether the original matrix is singular/ill-conditioned."),
        ("condition_number", "float or null", "Condition number if requested, else null."),
        ("matrix_size", "tuple", "Shape of the original matrix."),
    ]

    examples = [
        {
            "code_input": {
                "matrix": [[4, 7], [2, 6]],
                "use_pseudo_inverse": False,
            },
            "text_input": {
                "input_str": '{"matrix": [[4,7],[2,6]], "use_pseudo_inverse": false}',
            },
            "output": {
                "inverse": [[0.6, -0.7], [-0.2, 0.4]],
                "is_singular": False,
                "condition_number": None,
                "matrix_size": [2, 2],
            },
        },
        {
            "code_input": {
                "matrix": [[1, 2], [2, 4]],
                "use_pseudo_inverse": True,
            },
            "text_input": {
                "input_str": '{"matrix": [[1,2],[2,4]], "use_pseudo_inverse": true}',
            },
            "output": {
                "inverse": [[0.1, 0.2], [0.2, 0.4]],
                "is_singular": True,
                "condition_number": None,
                "matrix_size": [2, 2],
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

    def _run_base(self, matrix: list, use_pseudo_inverse: bool = False, condition_number: bool = False) -> dict:
        """Core logic: compute matrix inverse."""
        A = np.array(matrix, dtype=float)

        if A.ndim != 2:
            raise ChemMCPError("Input must be a 2D matrix.")
        if A.shape[0] == 0 or A.shape[1] == 0:
            raise ChemMCPError("Matrix cannot be empty.")

        shape = list(A.shape)
        is_singular = False

        try:
            if use_pseudo_inverse:
                A_inv = np.linalg.pinv(A)
                # Check if pseudo-inverse was needed due to singularity
                if A.shape[0] == A.shape[1]:
                    try:
                        np.linalg.det(A)
                        cond = np.linalg.cond(A)
                        if cond > 1e12:
                            is_singular = True
                    except np.linalg.LinAlgError:
                        is_singular = True
            else:
                det_val = float(np.linalg.det(A))
                if abs(det_val) < 1e-10:
                    raise ChemMCPError(
                        f"Matrix is singular (det ≈ {det_val:.2e}). "
                        "Use use_pseudo_inverse=True to compute pseudo-inverse."
                    )
                A_inv = np.linalg.inv(A)
        except np.linalg.LinAlgError as e:
            if not use_pseudo_inverse:
                raise ChemMCPError(f"Matrix inversion failed: {str(e)}. Try use_pseudo_inverse=True.")
            A_inv = np.linalg.pinv(A)
            is_singular = True

        result = {
            "inverse": [[round(v, 6) for v in row] for row in A_inv.tolist()],
            "is_singular": is_singular,
            "condition_number": None,
            "matrix_size": shape,
        }

        if condition_number and A.shape[0] == A.shape[1]:
            try:
                result["condition_number"] = round(float(np.linalg.cond(A)), 6)
            except Exception:
                pass

        logger.info(f"Inverted {shape} matrix, singular={is_singular}")
        return result

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            matrix = params.get("matrix")
            use_pinv = params.get("use_pseudo_inverse", False)
            cond_num = params.get("condition_number", False)
            return self._run_base(matrix, use_pseudo_inverse=use_pinv, condition_number=cond_num)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
