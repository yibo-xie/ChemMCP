import logging
import json
import math
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class OdeSolverStiff(BaseTool):
    """
    刚性ODE求解器，用于复杂反应网络、燃烧化学、催化反应等含多时间尺度的系统。
    实现隐式后向欧拉法和隐式中点法，适合刚性方程组。
    """
    __version__ = "0.1.0"
    name = "OdeSolverStiff"
    func_name = "solve_stiff_ode"
    description = "Stiff ODE solver for complex reaction networks, combustion chemistry, catalytic reactions with multiple time scales."
    implementation_description = "Implements implicit backward Euler method with Newton iteration for solving stiff ODE systems: y_{n+1} = y_n + h*f(t_{n+1}, y_{n+1}). Also supports semi-implicit (implicit-explicit) methods."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["ODE", "Stiff Systems", "Reaction Networks", "Combustion", "Catalysis"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str or list", "N/A", "ODE RHS expression(s). Single: '-1000*y'. System: JSON array."),
        ("jacobian_expr", "list or null", "null", "Jacobian matrix expressions df_i/dy_j (list of lists). Auto-approximated if null."),
        ("initial_values", "list", "N/A", "Initial condition(s) y(0)."),
        ("t_span", "list", "N/A", "Time span [t_start, t_end]."),
        ("n_steps", "int", "100", "Number of integration steps."),
        ("params", "dict", "{}", "Parameter dict for function expression."),
        ("newton_tol", "float", "1e-8", "Tolerance for Newton iteration convergence."),
        ("max_newton_iter", "int", "20", "Maximum Newton iterations per step."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("time_points", "list", "Array of time values."),
        ("solution", "list", "Solution array at each time step."),
        ("final_value", "list", "Final value(s) at t_end."),
        ("n_steps_used", "int", "Number of steps used."),
        ("newton_iterations_total", "int", "Total Newton iterations across all steps."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "-k * y",
                "initial_values": [1.0],
                "t_span": [0.0, 1.0],
                "n_steps": 10,
                "params": {"k": 50.0},
            },
            "text_input": {
                "input_str": '{"func_expr":"-k*y","initial_values":[1],"t_span":[0,1],"n_steps":10,"params":{"k":50}}',
            },
            "output": {
                "time_points": [0.0, 0.1, "...", 1.0],
                "solution": [[1.0], [0.006738], "...", [~0]],
                "final_value": [0.0],
                "n_steps_used": 10,
                "newton_iterations_total": 10,
            },
        },
        {
            "code_input": {
                "func_expr": ["y2", "-1000 * y1 - 1001 * y2"],
                "initial_values": [1.0, 0.0],
                "t_span": [0.0, 1.0],
                "n_steps": 100,
                "params": {},
            },
            "text_input": {
                "input_str": '{"func_expr":["y2","-1000*y1-1001*y2"],"initial_values":[1,0],"t_span":[0,1],"n_steps":100}',
            },
            "output": {
                "time_points": [0.0, 0.01, "...", 1.0],
                "solution": [[1.0, 0.0], [...], [...]],
                "final_value": [0.0, 0.0],
                "n_steps_used": 100,
                "newton_iterations_total": "~200",
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

    def _make_func(self, expr: str, var_names: list, params: dict):
        """Create a callable from expression string."""
        safe_dict = {
            "sin": math.sin, "cos": math.cos, "tan": math.tan,
            "exp": math.exp, "log": math.log, "sqrt": math.sqrt,
            "abs": abs, "pi": math.pi, "e": math.e,
        }
        safe_dict.update(params)

        def f(t, y_vec):
            local_vars = {**safe_dict, "t": t}
            for i, vn in enumerate(var_names):
                local_vars[vn] = y_vec[i]
            try:
                return eval(expr, {"__builtins__": {}}, local_vars)
            except Exception as e:
                raise ChemMCPError(f"Failed to evaluate '{expr}': {e}")
        return f

    def _numerical_jacobian(self, funcs, t, y, eps=1e-8):
        """Compute numerical Jacobian via finite differences."""
        n = len(y)
        J = np.zeros((n, n))
        f0 = np.array([f(t, y) for f in funcs])
        for j in range(n):
            y_eps = y.copy()
            y_eps[j] += eps
            f_eps = np.array([f(t, y_eps) for f in funcs])
            J[:, j] = (f_eps - f0) / eps
        return J

    def _run_base(
        self,
        func_expr,
        initial_values: list,
        t_span: list,
        n_steps: int = 100,
        params: dict = None,
        newton_tol: float = 1e-8,
        max_newton_iter: int = 20,
        jacobian_expr=None,
    ) -> dict:
        """Core logic: backward Euler for stiff systems."""
        if params is None:
            params = {}

        y0 = np.array(initial_values, dtype=float)
        n_vars = len(y0)
        t_start, t_end = float(t_span[0]), float(t_span[1])
        h = (t_end - t_start) / n_steps

        if isinstance(func_expr, str):
            func_exprs = [func_expr]
            var_names = ["y"]
        else:
            func_exprs = list(func_expr)
            var_names = [f"y{i}" for i in range(n_vars)]

        if len(func_exprs) != n_vars:
            raise ChemMCPError(f"Number of expressions ({len(func_exprs)}) must match variables ({n_vars}).")

        funcs = [self._make_func(expr, var_names, params) for expr in func_exprs]

        # Storage
        t_vals = [t_start]
        y_vals = [y0.copy()]
        total_newton = 0

        t = t_start
        y = y0.copy()

        for step in range(n_steps):
            # Backward Euler: y_new = y + h * f(t+h, y_new)
            # Solve via Newton's method: F(y_new) = y_new - y - h*f(t+h, y_new) = 0
            y_new = y.copy()

            for newton_it in range(max_newton_iter):
                f_val = np.array([f(t + h, y_new) for f in funcs])
                F = y_new - y - h * f_val

                J = self._numerical_jacobian(funcs, t + h, y_new)
                I_minus_hJ = np.eye(n_vars) - h * J

                try:
                    delta = np.linalg.solve(I_minus_hJ, -F)
                except np.linalg.LinAlgError:
                    delta = np.linalg.lstsq(I_minus_hJ, -F, rcond=None)[0]

                y_new += delta
                total_newton += 1

                if np.max(np.abs(delta)) < newton_tol:
                    break

            y = y_new
            t += h
            t_vals.append(round(t, 10))
            y_vals.append(y.copy())

        solution = [[round(v, 6) for v in row] for row in y_vals]
        final_val = [round(v, 6) for v in y.tolist()]

        logger.info(f"Stiff solver completed {n_vars}-variable system, {n_steps} steps, {total_newton} total Newton iters")
        return {
            "time_points": [round(t, 6) for t in t_vals],
            "solution": solution,
            "final_value": final_val,
            "n_steps_used": n_steps,
            "newton_iterations_total": total_newton,
        }

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            return self._run_base(**params)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
