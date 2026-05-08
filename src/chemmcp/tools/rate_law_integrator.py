import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Integrated rate law equations for different reaction orders
RATE_LAW_DATA = {
    0: {
        "name": "Zero-order kinetics",
        "differential": "−d[A]/dt = k",
        "integrated": "[A] = [A]₀ − kt",
        "half_life": "t₁/₂ = [A]₀ / (2k)",
        "units_of_k": "mol·L⁻¹·s⁻¹ (or M·s⁻¹)",
        "linear_plot": "[A] vs t → slope = −k",
    },
    1: {
        "name": "First-order kinetics",
        "differential": "−d[A]/dt = k[A]",
        "integrated": "ln[A] = ln[A]₀ − kt",
        "half_life": "t₁/₂ = ln(2) / k ≈ 0.693 / k",
        "units_of_k": "s⁻¹",
        "linear_plot": "ln[A] vs t → slope = −k",
    },
    2: {
        "name": "Second-order kinetics (single reactant)",
        "differential": "−d[A]/dt = k[A]²",
        "integrated": "1/[A] = 1/[A]₀ + kt",
        "half_life": "t₁/₂ = 1 / (k·[A]₀)",
        "units_of_k": "L·mol⁻¹·s⁻¹ (or M⁻¹·s⁻¹)",
        "linear_plot": "1/[A] vs t → slope = +k",
    },
}


@ChemMCPManager.register_tool
class RateLawIntegrator(BaseTool):
    """
    速率方程积分工具。
    对零级、一级、二级动力学进行积分，计算任意时刻浓度、半衰期、以及达到指定浓度所需时间。
    """
    __version__ = "0.1.0"
    name = "RateLawIntegrator"
    func_name = "integrate_rate_law"
    description = "Integrate rate laws for zero, first, and second order reactions. Calculate concentration at time t, half-life, and time to reach a given fraction or concentration."
    implementation_description = (
        "Uses analytical integrated rate laws: [A]=[A]₀−kt (zero-order), ln[A]=ln[A]₀−kt (first-order), "
        "1/[A]=1/[A]₀+kt (second-order). Calculates half-life (t₁/₂), fractional lifetime, "
        "and supports inverse calculation (time to reach target concentration)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Rate Law", "Integrated Rate Law", "Half-Life", "Chemical Kinetics"]
    required_envs = []

    code_input_sig = [
        ("order", "int", "N/A", "Reaction order: 0 (zero), 1 (first), or 2 (second)."),
        ("k", "float", "N/A", "Rate constant with proper units matching the reaction order."),
        ("initial_concentration", "float", "N/A", "Initial concentration [A]₀ in mol/L (M)."),
        ("time_s", "float", "0.0", "Reaction time in seconds."),
        ("target_fraction", "float", "None", "Optional target fraction remaining (e.g., 0.5 for half-life, 0.1 for 90% consumed)."),
        ("target_concentration", "float", "None", "Optional target concentration to calculate the required time."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: order k initial_concentration [time] [target_fraction]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing concentration, half_life, integrated_equation, and kinetic analysis."),
    ]

    examples = [
        {
            "code_input": {
                "order": 1,
                "k": 3.33e-4,
                "initial_concentration": 0.100,
                "time_s": 3600,
                "target_fraction": None,
                "target_concentration": None,
            },
            "text_input": {
                "input_params": "1 3.33e-4 0.100 3600",
            },
            "output": {
                "result": {
                    "order": 1,
                    "rate_constant_k": 3.33e-4,
                    "initial_concentration_M": 0.1,
                    "time_s": 3600,
                    "concentration_M": 0.301,
                    "fraction_remaining": 0.301,
                    "percent_consumed": 69.9,
                    "half_life_s": 2081.6,
                    "integrated_equation": "ln[A] = ln(0.100) − (3.33e-4)·t",
                    "kinetic_analysis": "First-order decay; ~70% consumed in 1 hour.",
                }
            }
        },
        {
            "code_input": {
                "order": 2,
                "k": 2.5e-3,
                "initial_concentration": 0.200,
                "time_s": 1800,
                "target_fraction": 0.25,
                "target_concentration": None,
            },
            "text_input": {
                "input_params": "2 2.5e-3 0.200 1800 0.25",
            },
            "output": {
                "result": {
                    "order": 2,
                    "rate_constant_k": 2.5e-3,
                    "initial_concentration_M": 0.2,
                    "time_s": 1800,
                    "concentration_M": 0.182,
                    "fraction_remaining": 0.911,
                    "half_life_s": 2000.0,
                    "time_to_fraction_0.25_s": 12000.0,
                    "integrated_equation": "1/[A] = 1/0.200 + (2.5e-3)·t",
                    "kinetic_analysis": "Second-order; slower than first-order equivalent due to [A]² dependence.",
                }
            }
        },
        {
            "code_input": {
                "order": 0,
                "k": 5.0e-6,
                "initial_concentration": 0.500,
                "time_s": 72000,
                "target_fraction": None,
                "target_concentration": 0.10,
            },
            "text_input": {
                "input_params": "0 5.0e-6 0.500 72000 None 0.10",
            },
            "output": {
                "result": {
                    "order": 0,
                    "rate_constant_k": 5e-06,
                    "initial_concentration_M": 0.5,
                    "time_s": 72000,
                    "concentration_M": 0.14,
                    "fraction_remaining": 0.28,
                    "half_life_s": 50000.0,
                    "time_to_target_conc_s": 80000.0,
                    "integrated_equation": "[A] = 0.500 − (5e-6)·t",
                    "kinetic_analysis": "Zero-order: constant rate until reactant depleted.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.rate_laws = dict(RATE_LAW_DATA)

    def _run_base(
        self,
        order: int,
        k: float,
        initial_concentration: float,
        time_s: float = 0.0,
        target_fraction: Optional[float] = None,
        target_concentration: Optional[float] = None,
    ) -> dict:
        """Core logic: integrate rate law and compute all quantities."""
        if order not in RATE_LAW_DATA:
            raise ChemMCPError(f"Unsupported order: {order}. Supported: {list(RATE_LAW_DATA.keys())}")
        if k < 0:
            raise ChemMCPError("Rate constant k must be non-negative.")
        if initial_concentration < 0:
            raise ChemMCPError("Initial concentration must be non-negative.")

        A0 = initial_concentration
        t = max(time_s, 0.0)
        law = RATE_LAW_DATA[order]

        # Calculate concentration at time t
        At = self._calc_concentration(order, k, A0, t)

        # Ensure non-negative
        At = max(At, 0.0)

        # Fraction remaining
        frac = At / A0 if A0 > 0 else 0.0

        # Half-life
        t_half = self._calc_half_life(order, k, A0)

        result = {
            "order": order,
            "rate_constant_k": k,
            "k_units": law["units_of_k"],
            "initial_concentration_M": round(A0, 6),
            "time_s": t,
            "concentration_M": round(At, 6),
            "fraction_remaining": round(frac, 6),
            "percent_consumed": round((1 - frac) * 100, 2),
            "half_life_s": round(t_half, 2),
            "integrated_equation": self._format_equation(law["integrated"], A0, k),
            "differential_form": law["differential"],
            "linear_plot_method": law["linear_plot"],
        }

        # Target fraction → time to reach it
        if target_fraction is not None:
            if not (0 < target_fraction <= 1):
                raise ChemMCPError("target_fraction must be between 0 and 1.")
            t_frac = self._time_for_fraction(order, k, A0, target_fraction)
            result[f"time_to_fraction_{target_fraction}_s"] = round(t_frac, 2)

        # Target concentration → time to reach it
        if target_concentration is not None:
            if target_concentration < 0:
                raise ChemMCPError("target_concentration must be non-negative.")
            if target_concentration > A0:
                raise ChemMCPError("target_concentration cannot exceed initial concentration for irreversible reaction.")
            t_target = self._time_for_concentration(order, k, A0, target_concentration)
            result["time_to_target_conc_s"] = round(t_target, 2)

        # Kinetic analysis summary
        result["kinetic_analysis"] = self._generate_analysis(order, k, A0, t, At, frac)

        return {"result": result}

    @staticmethod
    def _calc_concentration(order: int, k: float, A0: float, t: float) -> float:
        """Calculate [A] at time t."""
        if order == 0:
            return A0 - k * t
        elif order == 1:
            return A0 * math.exp(-k * t)
        elif order == 2:
            denom = 1 + k * A0 * t
            return A0 / denom if denom > 0 else 0.0
        else:
            raise ValueError(f"Unsupported order: {order}")

    @staticmethod
    def _calc_half_life(order: int, k: float, A0: float) -> float:
        """Calculate half-life t₁/₂."""
        if order == 0:
            return A0 / (2 * k) if k > 0 else float('inf')
        elif order == 1:
            return math.log(2) / k if k > 0 else float('inf')  # 0.693/k
        elif order == 2:
            return 1 / (k * A0) if k > 0 and A0 > 0 else float('inf')

    @staticmethod
    def _time_for_fraction(order: int, k: float, A0: float, f: float) -> float:
        """Calculate time to reach fraction f of [A]₀ remaining."""
        if order == 0:
            return A0 * (1 - f) / k if k > 0 else float('inf')
        elif order == 1:
            return -math.log(f) / k if k > 0 and f > 0 else float('inf')
        elif order == 2:
            return (1/f - 1) / (k * A0) if k > 0 and A0 > 0 and f > 0 else float('inf')

    @staticmethod
    def _time_for_concentration(order: int, k: float, A0: float, At: float) -> float:
        """Calculate time to reach target concentration [A]ₜ."""
        if order == 0:
            return (A0 - At) / k if k > 0 else float('inf')
        elif order == 1:
            if At <= 0:
                return float('inf')
            return math.log(A0 / At) / k if k > 0 else float('inf')
        elif order == 2:
            if At <= 0:
                return float('inf')
            return (1/At - 1/A0) / (k) if k > 0 and A0 > 0 else float('inf')

    @staticmethod
    def _format_equation(eq_template: str, A0: float, k: float) -> str:
        """Format integrated equation with actual values."""
        eq = eq_template.replace("[A]₀", f"{A0:.4g}")
        eq = eq.replace("k", f"{k:.4g}")
        return eq

    @staticmethod
    def _generate_analysis(order: int, k: float, A0: float, t: float, At: float, frac: float) -> str:
        parts = []
        name = RATE_LAW_DATA[order]["name"]

        if order == 0:
            parts.append(f"Zero-order: constant consumption rate of {k:.4e} M/s.")
            if At <= 0:
                parts.append("Reactant fully consumed — zero-order stops when [A]=0.")
            else:
                total_time = A0 / k
                parts.append(f"Total reaction duration: {total_time:.1f} s ({total_time/3600:.2f} h).")
        elif order == 1:
            parts.append(f"First-order: half-life is constant ({math.log(2)/k:.1f} s), independent of [A].")
            tau = 1 / k  # mean lifetime
            parts.append(f"Mean lifetime τ = 1/k = {tau:.1f} s.")
            if frac < 0.01:
                parts.append(">99% consumed — effectively complete.")
            elif frac < 0.37:
                parts.append(f"More than one half-life elapsed (~{-math.log(frac):.1f} half-lives).")
        elif order == 2:
            parts.append(f"Second-order: half-life depends on [A]₀ (t½ = 1/(k·[A]₀)).")
            parts.append("Higher initial concentration → faster consumption (rate ∝ [A]²).")
            if frac < 0.5:
                parts.append("More than one half-life elapsed.")

        # Time scale characterization
        if t > 0:
            n_halves = 0
            if order == 1 and k > 0:
                n_halves = k * t / math.log(2)
            if n_halves > 5:
                parts.append(f"~{n_halves:.1f} half-lives elapsed — reaction well advanced.")
            elif n_halves > 1:
                parts.append(f"~{n_halves:.1f} half-lives elapsed.")

        return " ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if len(parts) < 3:
                raise ChemMCPError("Need at least: order k initial_concentration")

            kwargs = {
                "order": int(parts[0]),
                "k": float(parts[1]),
                "initial_concentration": float(parts[2]),
            }
            if len(parts) > 3 and parts[3].lower() != "none":
                kwargs["time_s"] = float(parts[3])
            if len(parts) > 4 and parts[4].lower() != "none":
                kwargs["target_fraction"] = float(parts[4])
            if len(parts) > 5 and parts[5].lower() != "none":
                kwargs["target_concentration"] = float(parts[5])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
