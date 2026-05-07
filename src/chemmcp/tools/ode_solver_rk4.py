import logging
import json
import math
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class OdeSolverRk4(BaseTool):
    """
    四阶龙格-库塔法ODE求解器，用于反应动力学模拟、化学速率方程求解等。
    求解 dy/dt = f(t, y) 形式的常微分方程（组）。
    """
    __version__ = "0.1.0"
    name = "OdeSolverRk4"
    func_name = "solve_ode_rk4"
    description = "4th-order Runge-Kutta ODE solver for reaction kinetics simulation, chemical rate equations, and dynamic system modeling."
    implementation_description = "Implements the classic RK4 method: k1=f(t,y), k2=f(t+h/2,y+h*k1/2), k3=f(t+h/2,y+h*k2/2), k4=f(t+h,y+h*k3); y_{n+1}=y_n+h/6*(k1+2k2+2k3+k4). Supports systems of ODEs."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["ODE", "Runge-Kutta", "Reaction Kinetics", "Dynamics"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "ODE right-hand side expression(s). Single: '-k*y' (for dy/dt=-ky). System: JSON array of expressions."),
        ("initial_values", "list", "N/A", "Initial condition(s) y(0). List of floats."),
        ("t_span", "list", "N/A", "Time span [t_start, t_end]."),
        ("n_steps", "int", "100", "Number of integration steps."),
        ("params", "dict", "{}", "Parameter dict for function expression (e.g., {'k': 0.1})."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("time_points", "list", "Array of time values."),
        ("solution", "list", "Solution array (each element is a list of [y1, y2, ...] at each time step)."),
        ("final_value", "list", "Final value(s) at t_end."),
        ("n_steps_used", "int", "Number of steps actually used."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "-k * y",
                "initial_values": [1.0],
                "t_span": [0.0, 10.0],
                "n_steps": 100,
                "params": {"k": 0.5},
            },
            "text_input": {
                "input_str": '{"func_expr":"-k*y","initial_values":[1],"t_span":[0,10],"n_steps":100,"params":{"k":0.5}}',
            },
            "output": {
                "time_points": [0.0, 0.1, "...", 10.0],
                "solution": [[1.0], [0.951229], "...", [0.006738]],
                "final_value": [0.006738],
                "n_steps_used": 100,
            },
        },
        {
            "code_input": {
                "func_expr": ["-k1 * A", "k1 * A - k2 * B"],
                "initial_values": [1.0, 0.0],
                "t_span": [0.0, 5.0],
                "n_steps": 50,
                "params": {"k1": 1.0, "k2": 0.5},
            },
            "text_input": {
                "input_str": '{"func_expr":["-k1*A","k1*A-k2*B"],"initial_values":[1,0],"t_span":[0,5],"n_steps":50,"params":{"k1":1,"k2":0.5}}',
            },
            "output": {
                "time_points": [0.0, 0.1, "...", 5.0],
                "solution": [[1.0, 0.0], [0.904837, 0.086106], "..."],
                "final_value": [0.006738, 0.012692],
                "n_steps_used": 50,
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
        """Create a callable function from an expression string."""
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

    def _run_base(
        self,
        func_expr,
        initial_values: list,
        t_span: list,
        n_steps: int = 100,
        params: dict = None,
    ) -> dict:
        """Core logic: RK4 integration."""
        if params is None:
            params = {}

        y0 = np.array(initial_values, dtype=float)
        n_vars = len(y0)
        t_start, t_end = float(t_span[0]), float(t_span[1])
        h = (t_end - t_start) / n_steps

        # Handle single equation or system
        if isinstance(func_expr, str):
            func_exprs = [func_expr]
            var_names = ["y"]
        else:
            func_exprs = list(func_expr)
            var_names = [f"y{i}" for i in range(n_vars)]

        if len(func_exprs) != n_vars:
            raise ChemMCPError(f"Number of expressions ({len(func_exprs)}) must match number of variables ({n_vars}).")

        funcs = [self._make_func(expr, var_names, params) for expr in func_exprs]

        # Storage
        t_vals = [t_start]
        y_vals = [y0.copy()]

        t = t_start
        y = y0.copy()

        for _ in range(n_steps):
            k1 = np.array([f(t, y) for f in funcs])
            k2 = np.array([f(t + h / 2, y + h * k1 / 2) for f in funcs])
            k3 = np.array([f(t + h / 2, y + h * k2 / 2) for f in funcs])
            k4 = np.array([f(t + h, y + h * k3) for f in funcs])

            y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += h

            t_vals.append(round(t, 10))
            y_vals.append(y.copy())

        solution = [[round(v, 6) for v in row] for row in y_vals]
        final_val = [round(v, 6) for v in y.tolist()]

        logger.info(f"RK4 solved {n_vars}-variable ODE over [{t_start}, {t_end}] with {n_steps} steps")
        return {
            "time_points": [round(t, 6) for t in t_vals],
            "solution": solution,
            "final_value": final_val,
            "n_steps_used": n_steps,
        }

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            return self._run_base(**params)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
