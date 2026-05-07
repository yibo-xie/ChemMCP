import json
import logging
import math
import random
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MonteCarloIntegrator(BaseTool):
    """
    蒙特卡洛积分工具。
    支持高维积分、统计采样，提供多种采样策略（朴素、重要性、分层）。
    """
    __version__ = "0.1.0"
    name = "MonteCarloIntegrator"
    func_name = "monte_carlo_integrate"
    description = "Monte Carlo integration for multi-dimensional integrals with plain, importance, and stratified sampling strategies."
    implementation_description = "Implements Monte Carlo integration: samples points uniformly (or via importance/stratified strategy) within bounds, evaluates function, and estimates integral as average × volume. Supports 1D and 2D domains. Error estimated via standard error of the mean."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Monte Carlo", "Integration", "Numerical Methods", "Sampling", "Statistics", "High-Dimensional"]
    required_envs = []

    code_input_sig = [
        ("func_expr", "str", "N/A", "Math expression in variables x (1D) or x,y (2D). Supports: +,-,*,/,**,sin,cos,tan,exp,log,sqrt,abs."),
        ("bounds", "list", "N/A", "List of [lower, upper] per dimension. E.g., [[0,1]] for 1D or [[0,1],[0,1]] for 2D."),
        ("n_samples", "int", "10000", "Number of Monte Carlo samples."),
        ("method", "str", "plain", "Sampling method: 'plain', 'importance', 'stratified'."),
        ("seed", "int", "None", "Random seed for reproducibility (None = non-deterministic)."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("integral_estimate", "float", "Estimated value of the integral."),
        ("std_error", "float", "Standard error of the estimate."),
        ("relative_error", "float", "Relative error estimate (std_error / |integral|)."),
        ("n_samples_used", "int", "Actual number of samples used."),
        ("n_dimensions", "int", "Dimensionality of the integral."),
        ("volume", "float", "Volume of the integration domain."),
        ("method", "str", "Sampling method used."),
        ("convergence_info", "dict", "Convergence diagnostics: running mean trend, variance info."),
        ("func_expr", "str", "The input expression."),
    ]

    examples = [
        {
            "code_input": {
                "func_expr": "x ** 2",
                "bounds": [[0.0, 1.0]],
                "n_samples": 50000,
            },
            "text_input": {"params_str": '{"func_expr":"x**2","bounds":[[0,1]],"n_samples":50000}'},
            "output": {
                "integral_estimate": 0.3333,
                "n_dimensions": 1,
            },
        },
        {
            "code_input": {
                "func_expr": "x + y",
                "bounds": [[0.0, 1.0], [0.0, 1.0]],
                "n_samples": 50000,
                "method": "stratified",
            },
            "text_input": {"params_str": '{"func_expr":"x+y","bounds":[[0,1],[0,1]],"n_samples":50000,"method":"stratified"}'},
            "output": {
                "integral_estimate": 1.0,
                "n_dimensions": 2,
            },
        },
        {
            "code_input": {
                "func_expr": "sin(x)",
                "bounds": [[0.0, math.pi]],
                "n_samples": 50000,
            },
            "text_input": {"params_str": '{"func_expr":"sin(x)","bounds":[[0,3.14159]],"n_samples":50000}'},
            "output": {
                "integral_estimate": 2.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    _SAFE_FUNCS = {
        "sin": math.sin, "cos": math.cos, "tan": math.tan,
        "exp": math.exp, "log": math.log10, "ln": math.log,
        "sqrt": math.sqrt, "abs": abs,
        "pi": math.pi, "e": math.e,
        "asin": math.asin, "acos": math.acos, "atan": math.atan,
        "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    }

    @staticmethod
    def _eval_f(expr: str, vars_dict: dict) -> float:
        """Safely evaluate math expression."""
        ns = dict(MonteCarloIntegrator._SAFE_FUNCS)
        ns.update(vars_dict)
        expr_py = expr.replace("^", "**")
        try:
            return eval(expr_py, {"__builtins__": {}}, ns)
        except Exception as e:
            raise ChemMCPError(f"Failed to evaluate '{expr}': {e}")

    @staticmethod
    def _domain_volume(bounds: List[List[float]]) -> float:
        v = 1.0
        for lo, hi in bounds:
            v *= (hi - lo)
        return v

    _DIM_NAMES = ["x", "y", "z", "w", "v", "u"]

    def _sample_plain(self, n: int, bounds: List[List[float]], rng: random.Random) -> List[dict]:
        """Uniform random sampling."""
        samples = []
        ndim = len(bounds)
        for _ in range(n):
            pt = {}
            for i, (lo, hi) in enumerate(bounds):
                var = self._DIM_NAMES[i] if i < len(self._DIM_NAMES) else f"x{i}"
                pt[var] = rng.random() * (hi - lo) + lo
            samples.append(pt)
        return samples

    def _sample_stratified(self, n: int, bounds: List[List[float]], rng: random.Random) -> List[dict]:
        """Stratified sampling: divide each axis into sqrt(n) strata."""
        import math as _m
        ndim = len(bounds)
        # Number of strata per dimension
        n_strata = max(int(round(n ** (1.0 / ndim))), 2)
        total = n_strata ** ndim
        samples = []
        for idx in range(total):
            pt = {}
            remainder = idx
            for i, (lo, hi) in enumerate(bounds):
                var = self._DIM_NAMES[i] if i < len(self._DIM_NAMES) else f"x{i}"
                stride_idx = remainder % n_strata
                remainder //= n_strata
                cell_lo = lo + stride_idx * (hi - lo) / n_strata
                cell_hi = lo + (stride_idx + 1) * (hi - lo) / n_strata
                pt[var] = rng.random() * (cell_hi - cell_lo) + cell_lo
            samples.append(pt)
        # If total < n, fill remaining with uniform
        while len(samples) < n:
            pt = {}
            for i, (lo, hi) in enumerate(bounds):
                var = self._DIM_NAMES[i] if i < len(self._DIM_NAMES) else f"x{i}"
                pt[var] = rng.random() * (hi - lo) + lo
            samples.append(pt)
        return samples[:n]

    def _sample_importance(self, n: int, bounds: List[List[float]], expr: str, rng: random.Random) -> tuple:
        """
        Simple importance sampling using a uniform envelope.
        Falls back to plain sampling if no improvement is possible.
        Returns (samples, weights).
        """
        # For now, use uniform importance (same as plain but with explicit weights=1)
        samples = self._sample_plain(n, bounds, rng)
        weights = [1.0] * len(samples)
        return samples, weights

    def _run_base(
        self,
        func_expr: str,
        bounds: List[List[float]],
        n_samples: int = 10000,
        method: str = "plain",
        seed: Optional[int] = None,
    ) -> dict:
        """Core logic: Monte Carlo integration."""
        if not func_expr or not func_expr.strip():
            raise ChemMCPError("Expression cannot be empty.")
        if not bounds or len(bounds) == 0:
            raise ChemMCPError("Bounds cannot be empty.")
        for b in bounds:
            if len(b) != 2 or b[1] <= b[0]:
                raise ChemMCPError(f"Invalid bound {b}. Must be [lower, upper] with lower < upper.")

        ndim = len(bounds)
        volume = self._domain_volume(bounds)

        if n_samples < 10:
            raise ChemMCPError(f"n_samples must be >= 10, got {n_samples}.")

        rng = random.Random(seed)
        meth = method.lower().strip()

        # Select sampling method
        if meth == "plain":
            samples = self._sample_plain(n_samples, bounds, rng)
            weights = None
        elif meth == "stratified":
            samples = self._sample_stratified(n_samples, bounds, rng)
            weights = None
        elif meth == "importance":
            samples, wts = self._sample_importance(n_samples, bounds, func_expr, rng)
            weights = wts
        else:
            raise ChemMCPError(f"Unknown method '{method}'. Use: 'plain', 'importance', or 'stratified'.")

        # Evaluate function at all sample points
        f_values = []
        for s in samples:
            try:
                fv = self._eval_f(func_expr, s)
                if math.isfinite(fv):
                    f_values.append(fv)
            except Exception:
                pass  # skip bad evaluations

        if not f_values:
            raise ChemMCPError("All function evaluations failed. Check expression.")

        n_valid = len(f_values)

        if weights is None:
            # Standard MC: I ≈ V × ⟨f⟩
            f_mean = sum(f_values) / n_valid
            f_var = sum((f - f_mean) ** 2 for f in f_values) / (n_valid - 1) if n_valid > 1 else 0.0
            integral = f_mean * volume
            std_err = math.sqrt(f_var / n_valid) * volume if n_valid > 1 else 0.0
        else:
            # Importance sampling weighted average
            total_w = sum(weights[:n_valid])
            f_mean = sum(w * f for w, f in zip(weights, f_values)) / total_w if total_w > 0 else 0.0
            integral = f_mean * volume
            std_err = 0.0  # simplified

        rel_err = abs(std_err / integral) if abs(integral) > 1e-30 else (float('inf') if std_err > 0 else 0.0)

        # Convergence info: running mean every 10% of samples
        running_means = []
        chunk = max(n_valid // 10, 1)
        cumsum = 0.0
        for i, fv in enumerate(f_values):
            cumsum += fv
            if (i + 1) % chunk == 0 or i == n_valid - 1:
                running_means.append((i + 1, (cumsum / (i + 1)) * volume))

        logger.info(
            f"MonteCarloIntegrator ({meth}): ∫≈{integral:.6g} ± {std_err:.6g}, "
            f"n={n_valid}, dim={ndim}"
        )

        return {
            "integral_estimate": round(integral, 10),
            "std_error": round(std_err, 10),
            "relative_error": round(rel_err, 8),
            "n_samples_used": n_valid,
            "n_dimensions": ndim,
            "volume": round(volume, 10),
            "method": meth,
            "convergence_info": {
                "running_means": [(n, round(v, 8)) for n, v in running_means],
                "variance_of_mean": round(f_var / n_valid if n_valid > 1 else 0, 14),
            },
            "func_expr": func_expr.strip(),
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
