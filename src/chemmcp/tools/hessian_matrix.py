import logging
import math
from typing import Dict, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

from .partial_derivative import _safe_eval

logger = logging.getLogger(__name__)


def _compute_eigenvalues_2x2(matrix: List[List[float]]) -> List[float]:
    """Compute eigenvalues of a 2x2 matrix analytically."""
    a, b = matrix[0][0], matrix[0][1]
    c, d = matrix[1][0], matrix[1][1]
    trace = a + d
    det_val = a * d - b * c
    discriminant = trace ** 2 - 4 * det_val
    if discriminant < 0:
        discriminant = 0.0
    sqrt_disc = math.sqrt(discriminant)
    return [round((trace + sqrt_disc) / 2, 8), round((trace - sqrt_disc) / 2, 8)]


def _compute_eigenvalues_general(matrix: List[List[float]]) -> List[float]:
    """Compute eigenvalues using power iteration for the dominant eigenvalue and deflation."""
    n = len(matrix)
    if n == 1:
        return [round(matrix[0][0], 8)]
    if n == 2:
        return _compute_eigenvalues_2x2(matrix)

    # For larger matrices, use numpy if available
    try:
        import numpy as np
        arr = np.array(matrix, dtype=float)
        eigenvalues = np.linalg.eigvalsh(arr) if _is_symmetric(arr) else np.linalg.eigvals(arr)
        return sorted([round(float(ev), 8) for ev in eigenvalues])
    except ImportError:
        # Fallback: use QR algorithm (simplified)
        return _qr_eigenvalues(matrix)


def _is_symmetric(mat) -> bool:
    try:
        import numpy as np
        return np.allclose(mat, mat.T, atol=1e-10)
    except ImportError:
        n = len(mat)
        for i in range(n):
            for j in range(i + 1, n):
                if abs(mat[i][j] - mat[j][i]) > 1e-10:
                    return False
        return True


def _qr_eigenvalues(matrix: List[List[float]], max_iter: int = 1000, tol: float = 1e-10) -> List[float]:
    """Simple QR iteration for eigenvalue estimation."""
    try:
        import numpy as np
        A = np.array(matrix, dtype=float)
        n = A.shape[0]
        for _ in range(max_iter):
            Q, R = np.linalg.qr(A)
            A = R @ Q
            # Check for convergence (lower triangular)
            off_diag = sum(abs(A[i, j]) for i in range(n) for j in range(i))
            if off_diag < tol:
                break
        return sorted([round(float(A[i, i]), 8) for i in range(n)])
    except ImportError:
        # Pure Python fallback — just return diagonal as approximation
        n = len(matrix)
        return [round(matrix[i][i], 8) for i in range(n)]


def _determinant(matrix: List[List[float]]) -> float:
    """Compute determinant recursively."""
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    try:
        import numpy as np
        return float(np.linalg.det(np.array(matrix)))
    except ImportError:
        det = 0.0
        for j in range(n):
            minor = [row[:j] + row[j+1:] for row in matrix[1:]]
            det += ((-1) ** j) * matrix[0][j] * _determinant(minor)
        return det


@ChemMCPManager.register_tool
class HessianMatrix(BaseTool):
    """
    Hessian矩阵计算工具 —— 振动频率、过渡态确认。
    计算多元函数的Hessian矩阵（二阶偏导数矩阵）。
    """
    __version__ = "0.1.0"
    name = "HessianMatrix"
    func_name = "hessian_matrix"
    description = (
        "Compute the Hessian matrix (matrix of second partial derivatives). "
        "Used for vibrational frequency analysis and transition state confirmation."
    )
    implementation_description = (
        "Computes ∂²f/∂xᵢ∂xⱼ using central finite differences. "
        "Also returns eigenvalues and determinant for saddle point / minimum analysis."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Numerical Methods", "Vibrational Analysis", "Transition State", "Hessian"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Mathematical expression string."),
        ("variables", "list", "N/A", "List of variable names."),
        ("eval_point", "dict", "N/A", "Dictionary mapping variable names to float values."),
        ("step_size", "float", "1e-5", "Step size for numerical differentiation (default: 1e-5)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A",
         "Space-separated: 'func_expr var1,var2,... val1,val2,... [step_size]'. "
         "Example: 'x**2+x*y+y**2 x,y 1,1'"),
    ]

    output_sig = [
        ("hessian", "list", "Hessian matrix as list of lists (2D array)."),
        ("eigenvalues", "list", "Eigenvalues of the Hessian matrix."),
        ("determinant", "float", "Determinant of the Hessian matrix."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "x**2 + x*y + y**2",
                "variables": ["x", "y"],
                "eval_point": {"x": 1.0, "y": 1.0},
                "step_size": 1e-5,
            },
            "text_input": {
                "input_str": "x**2+x*y+y**2 x,y 1,1",
            },
            "output": {
                "hessian": [[2.0, 1.0], [1.0, 2.0]],
                "eigenvalues": [3.0, 1.0],
                "determinant": 3.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        func_expr: str,
        variables: List[str],
        eval_point: Dict[str, float],
        step_size: float = 1e-5,
    ) -> Dict:
        n = len(variables)
        h = step_size

        # Build Hessian matrix
        hessian = []
        for i in range(n):
            row = []
            vi = variables[i]
            for j in range(n):
                vj = variables[j]
                if i == j:
                    # Second pure derivative: ∂²f/∂vi²
                    p_pp = dict(eval_point); p_pp[vi] += h
                    p_pm = dict(eval_point); p_pm[vi] += h; p_pm[vj] += h  # not used for diagonal
                    p_mp = dict(eval_point); p_mp[vi] -= h
                    p_mm = dict(eval_point); p_mm[vi] -= h

                    f_pp = _safe_eval(func_expr, p_pp)
                    f_p = _safe_eval(func_expr, eval_point)
                    f_mp = _safe_eval(func_expr, p_mp)
                    d2 = (f_pp - 2 * f_p + f_mp) / (h * h)
                else:
                    # Mixed partial: ∂²f/∂vi∂vj
                    p_pp = dict(eval_point); p_pp[vi] += h; p_pp[vj] += h
                    p_pm = dict(eval_point); p_pm[vi] += h; p_pm[vj] -= h
                    p_mp = dict(eval_point); p_mp[vi] -= h; p_mp[vj] += h
                    p_mm = dict(eval_point); p_mm[vi] -= h; p_mm[vj] -= h

                    f_pp = _safe_eval(func_expr, p_pp)
                    f_pm = _safe_eval(func_expr, p_pm)
                    f_mp = _safe_eval(func_expr, p_mp)
                    f_mm = _safe_eval(func_expr, p_mm)
                    d2 = (f_pp - f_pm - f_mp + f_mm) / (4 * h * h)
                row.append(round(d2, 8))
            hessian.append(row)

        eigenvalues = _compute_eigenvalues_general(hessian)
        det_val = round(_determinant(hessian), 8)

        logger.info(f"Hessian at {eval_point}: eigenvalues={eigenvalues}, det={det_val}")
        return {
            "hessian": hessian,
            "eigenvalues": eigenvalues,
            "determinant": det_val,
        }

    def _run_text(self, input_str: str) -> Dict:
        try:
            parts = input_str.strip().split()
            if len(parts) < 3:
                raise ValueError("Need at least 3 parts: func_expr variables values")

            func_expr = parts[0]
            variables = parts[1].split(",")
            values = [float(v) for v in parts[2].split(",")]
            step_size = float(parts[3]) if len(parts) > 3 else 1e-5

            if len(variables) != len(values):
                raise ValueError(f"Mismatch: {len(variables)} variables but {len(values)} values")

            eval_point = dict(zip(variables, values))
            return self._run_base(func_expr, variables, eval_point, step_size)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
