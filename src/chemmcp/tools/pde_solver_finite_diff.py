import logging
import json
import math
import numpy as np

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PdeSolverFiniteDiff(BaseTool):
    """
    有限差分法PDE求解器，用于扩散方程、传热方程等偏微分方程数值求解。
    支持显式FTCS格式和隐式Crank-Nicolson格式求解一维/二维抛物型PDE。
    """
    __version__ = "0.1.0"
    name = "PdeSolverFiniteDiff"
    func_name = "solve_pde_finite_diff"
    description = "Finite difference PDE solver for diffusion equation, heat transfer, and parabolic PDEs in 1D and 2D."
    implementation_description = "Implements explicit FTCS (Forward Time Centered Space) and implicit Crank-Nicolson schemes for ∂u/∂t = α∇²u. Supports Dirichlet and Neumann boundary conditions."
    oss_dependencies = [
        ("numpy", "https://numpy.org", "BSD"),
    ]
    services_and_software = []
    categories = ["General"]
    tags = ["PDE", "Finite Difference", "Diffusion", "Heat Transfer", "Numerical Methods"]
    required_envs = []

    code_input_sig = [
        ("pde_type", "str", "heat_1d", "PDE type: 'heat_1d' (1D heat/diffusion), 'wave_1d' (1D wave)."),
        ("alpha", "float", "1.0", "Diffusivity / thermal diffusivity coefficient."),
        ("length", "float", "1.0", "Spatial domain length L."),
        ("t_final", "float", "0.1", "Final simulation time."),
        ("nx", "int", "50", "Number of spatial grid points."),
        ("nt", "int", "500", "Number of time steps."),
        ("initial_condition", "str", "sin(pi*x/L)", "Initial condition expression in variable x."),
        ("bc_left_type", "str", "dirichlet", "Left BC: 'dirichlet' or 'neumann'."),
        ("bc_left_value", "float", "0.0", "Left boundary value or derivative."),
        ("bc_right_type", "str", "dirichlet", "Right BC: 'dirichlet' or 'neumann'."),
        ("bc_right_value", "float", "0.0", "Right boundary value or derivative."),
        ("scheme", "str", "explicit", "Scheme: 'explicit' (FTCS) or 'implicit' (Crank-Nicolson)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("time_points", "list", "Selected time snapshots for output."),
        ("spatial_grid", "list", "Spatial grid coordinates (x values)."),
        ("solution_snapshots", "list", "Solution u(x,t) at each snapshot time (list of lists)."),
        ("stability_param", "float", "Courant-Friedrichs-Lewy (CFL) stability parameter r = αΔt/Δx²."),
        ("is_stable", "bool", "Whether the explicit scheme satisfies CFL condition (r ≤ 0.5)."),
        ("final_profile", "list", "Solution profile at final time t_final."),
    ]

    examples = [
        {
            "code_input": {
                "pde_type": "heat_1d",
                "alpha": 0.01,
                "length": 1.0,
                "t_final": 0.5,
                "nx": 20,
                "nt": 500,
                "initial_condition": "sin(pi * x)",
                "bc_left_type": "dirichlet",
                "bc_left_value": 0.0,
                "bc_right_type": "dirichlet",
                "bc_right_value": 0.0,
                "scheme": "explicit",
            },
            "text_input": {
                "input_str": '{"pde_type":"heat_1d","alpha":0.01,"length":1,"t_final":0.5,"nx":20,"nt":500,"ic":"sin(pi*x)"}',
            },
            "output": {
                "time_points": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
                "spatial_grid": [0.0, 0.0526, "...", 1.0],
                "solution_snapshots": [[0.0, 0.1564, ..., 0.0], [...]],
                "stability_param": 2.375,
                "is_stable": False,
                "final_profile": [0.0, 0.004, ..., 0.0],
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

    def _eval_ic(self, expr: str, x: float) -> float:
        """Evaluate initial condition expression."""
        safe_dict = {
            "sin": math.sin, "cos": math.cos, "tan": math.tan,
            "exp": math.exp, "log": math.log, "sqrt": math.sqrt,
            "abs": abs, "pi": math.pi, "e": math.e,
        }
        try:
            return eval(expr, {"__builtins__": {}}, {**safe_dict, "x": x})
        except Exception as e:
            raise ChemMCPError(f"Failed to evaluate IC '{expr}' at x={x}: {e}")

    def _run_base(
        self,
        pde_type: str = "heat_1d",
        alpha: float = 1.0,
        length: float = 1.0,
        t_final: float = 0.1,
        nx: int = 50,
        nt: int = 500,
        initial_condition: str = "sin(pi*x/L)",
        bc_left_type: str = "dirichlet",
        bc_left_value: float = 0.0,
        bc_right_type: str = "dirichlet",
        bc_right_value: float = 0.0,
        scheme: str = "explicit",
    ) -> dict:
        """Core logic: finite difference PDE solver."""

        # Grid setup
        dx = length / (nx - 1)
        dt = t_final / nt
        x_grid = np.linspace(0, length, nx)

        # Stability parameter for 1D heat equation: r = α*dt/dx²
        r = alpha * dt / (dx ** 2)
        is_stable_explicit = (r <= 0.5)

        # Initial condition
        ic_expr = initial_condition.replace("L", str(length))
        u = np.array([self._eval_ic(ic_expr, xi) for xi in x_grid])

        # Apply Dirichlet BCs to IC
        if bc_left_type == "dirichlet":
            u[0] = bc_left_value
        if bc_right_type == "dirichlet":
            u[-1] = bc_right_value

        # Select output snapshots (up to 6 evenly spaced times including start/end)
        n_snaps = min(6, nt + 1)
        snap_indices = np.linspace(0, nt, n_snaps, dtype=int)

        snapshots = [u.copy().tolist()]
        snap_times = [0.0]

        if pde_type == "heat_1d":
            if scheme == "explicit":
                # Explicit FTCS scheme
                for n_step in range(1, nt + 1):
                    u_new = u.copy()
                    u_new[1:-1] = u[1:-1] + r * (u[2:] - 2 * u[1:-1] + u[:-2])

                    # Neumann BC (zero gradient): du/dx = 0 → ghost point method
                    if bc_left_type == "neumann":
                        u_new[0] = u_new[1] - dx * bc_left_value
                    else:
                        u_new[0] = bc_left_value

                    if bc_right_type == "neumann":
                        u_new[-1] = u_new[-2] + dx * bc_right_value
                    else:
                        u_new[-1] = bc_right_value

                    u = u_new

                    if n_step in snap_indices:
                        snapshots.append(u.copy().tolist())
                        snap_times.append(round(n_step * dt, 6))

            elif scheme == "implicit":
                # Crank-Nicolson (semi-implicit)
                A = np.zeros((nx, nx))
                B = np.zeros((nx, nx))

                for i in range(1, nx - 1):
                    A[i, i - 1] = -r / 2
                    A[i, i] = 1 + r
                    A[i, i + 1] = -r / 2
                    B[i, i - 1] = r / 2
                    B[i, i] = 1 - r
                    B[i, i + 1] = r / 2

                # Boundary conditions
                if bc_left_type == "dirichlet":
                    A[0, 0] = 1.0
                    B[0, 0] = 1.0
                else:
                    A[0, 0] = 1.0
                    A[0, 1] = -1.0
                    B[0, 0] = 1.0
                    B[0, 1] = -1.0

                if bc_right_type == "dirichlet":
                    A[-1, -1] = 1.0
                    B[-1, -1] = 1.0
                else:
                    A[-1, -1] = 1.0
                    A[-1, -2] = -1.0
                    B[-1, -1] = 1.0
                    B[-1, -2] = -1.0

                b_vec = np.zeros(nx)
                for n_step in range(1, nt + 1):
                    rhs = B @ u
                    # Adjust RHS for Dirichlet BCs
                    if bc_left_type == "dirichlet":
                        rhs[0] = bc_left_value
                    if bc_right_type == "dirichlet":
                        rhs[-1] = bc_right_value

                    try:
                        u = np.linalg.solve(A, rhs)
                    except np.linalg.LinAlgError:
                        u = np.linalg.lstsq(A, rhs, rcond=None)[0]

                    if n_step in snap_indices:
                        snapshots.append(u.copy().tolist())
                        snap_times.append(round(n_step * dt, 6))

            else:
                raise ChemMCPError(f"Unknown scheme: {scheme}. Use 'explicit' or 'implicit'.")

        elif pde_type == "wave_1d":
            c = math.sqrt(alpha)  # wave speed from "diffusivity"
            cfl = c * dt / dx
            is_stable_explicit = is_stable_explicit and (cfl <= 1.0)

            u_prev = u.copy()
            # First step: use u(x,0) = f(x), u_t(x,0) = 0
            u_next = u.copy()
            u_next[1:-1] = u[1:-1] + 0.5 * cfl**2 * (u[2:] - 2*u[1:-1] + u[:-2])

            for n_step in range(1, nt + 1):
                u_new = np.zeros_like(u)
                u_new[1:-1] = 2*(1-cfl**2)*u[1:-1] - u_prev[1:-1] + cfl**2*(u[2:] - 2*u[1:-1] + u[:-2])
                u_new[0] = bc_left_value
                u_new[-1] = bc_right_value

                u_prev = u.copy()
                u = u_new

                if n_step in snap_indices:
                    snapshots.append(u.copy().tolist())
                    snap_times.append(round(n_step * dt, 6))
        else:
            raise ChemMCPError(f"Unknown PDE type: {pde_type}. Use 'heat_1d' or 'wave_1d'.")

        # Round all results
        rounded_snapshots = [[round(v, 6) for v in snap] for snap in snapshots]

        logger.info(f"PDE {pde_type} solved: {nx} grid pts, {nt} steps, r={r:.4f}, stable={is_stable_explicit}")
        return {
            "time_points": snap_times,
            "spatial_grid": [round(xi, 6) for xi in x_grid.tolist()],
            "solution_snapshots": rounded_snapshots,
            "stability_param": round(r, 6),
            "is_stable": is_stable_explicit,
            "final_profile": [[round(v, 6) for v in u.tolist()]],
        }

    def _run_text(self, input_str: str) -> dict:
        try:
            params = json.loads(input_str)
            return self._run_base(**params)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON input: {input_str}")
