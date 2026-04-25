import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HalfLifeCalculator(BaseTool):
    """
    计算各级反应的半衰期工具。
    支持零级、一级、二级和n级反应（n≠1）的半衰期计算。
    """
    __version__ = "0.1.0"
    name = "HalfLifeCalculator"
    func_name = "calculate_half_life"
    description = "Calculate half-life (t₁/₂) for zero, first, second, and nth order reactions."
    implementation_description = "Uses standard formulas: t₁/₂ = [A]₀/(2k) (zero), ln(2)/k (first), 1/(k[A]₀) (second), and general formula for nth order."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Half-Life", "Chemical Kinetics", "Reaction Order"]
    required_envs    = []

    code_input_sig = [
        ("reaction_order", "float", "N/A", "Reaction order n (0, 1, 2, or any positive number)."),
        ("rate_constant", "float", "N/A", "Rate constant k (unit depends on reaction order)."),
        ("initial_concentration", "float", "1.0", "Initial concentration [A]₀ (required for n ≠ 1)."),
        ("time_unit", "str", "s", "Time unit for output: 's', 'min', 'h'."),
        ("conc_unit", "str", "M", "Concentration unit for input."),
        # For fractional life calculation
        ("fraction_remaining", "float", "0.5", "Fraction remaining (default 0.5 for half-life; use e.g., 0.25 for quarter-life)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: reaction_order rate_constant [initial_concentration time_unit conc_unit fraction_remaining]."),
    ]

    output_sig = [
        ("half_life", "float", "Half-life t₁/₂ in specified time unit."),
        ("half_life_unit", "str", "Unit of half-life."),
        ("formula_used", "str", "The formula used for calculation."),
        ("analysis", "str", "Detailed analysis including interpretation of result."),
        # Additional useful outputs
        ("fractional_life", "float", "The fractional life calculated (for non-0.5 fractions)."),
        ("percent_decomposed", "float", "Percentage decomposed after this time."),
    ]

    examples         = [
        {
            "code_input": {
                "reaction_order": 1,
                "rate_constant": 0.00693,
                "initial_concentration": 1.0,
                "time_unit": 's',
                "conc_unit": 'M',
                "fraction_remaining": 0.5
            },
            "text_input": {
                "input_params": '1 0.00693'
            },
            "output": {
                "half_life": 100.04,
                "half_life_unit": 's',
                "formula_used": 't_1/2 = ln(2) / k',
                "analysis": 'First-order t_1/2 ~ 100 s.',
                "fractional_life": 100.04,
                "percent_decomposed": 50.0
            }
        },
        {
            "code_input": {
                "reaction_order": 2,
                "rate_constant": 0.05,
                "initial_concentration": 0.1,
                "time_unit": 'min',
                "conc_unit": 'M',
                "fraction_remaining": 0.5
            },
            "text_input": {
                "input_params": '2 0.05 0.1 min'
            },
            "output": {
                "half_life": 200.0,
                "half_life_unit": 'min',
                "formula_used": 't_1/2 = 1 / (k*[A]_0)',
                "analysis": 'Second-order t_1/2 depends on [A]_0.',
                "fractional_life": 200.0,
                "percent_decomposed": 50.0
            }
        },
        {
            "code_input": {
                "reaction_order": 0,
                "rate_constant": 0.01,
                "initial_concentration": 1.0,
                "time_unit": 's',
                "conc_unit": 'M',
                "fraction_remaining": 0.5
            },
            "text_input": {
                "input_params": '0 0.01 1.0 s'
            },
            "output": {
                "half_life": 50.0,
                "half_life_unit": 's',
                "formula_used": 't_1/2 = [A]_0 / (2k)',
                "analysis": 'Zero-order constant rate.',
                "fractional_life": 50.0,
                "percent_decomposed": 50.0
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        reaction_order: float,
        rate_constant: float,
        initial_concentration: float = 1.0,
        time_unit: str = "s",
        conc_unit: str = "M",
        fraction_remaining: float = 0.5,
    ) -> dict:
        if rate_constant <= 0:
            raise ChemMCPError("Rate constant must be positive.")
        if fraction_remaining <= 0 or fraction_remaining >= 1:
            raise ChemMCPError("Fraction remaining must be between 0 and 1 (exclusive).")
        if initial_concentration <= 0:
            raise ChemMCPError("Initial concentration must be positive.")

        n = reaction_order
        k = rate_constant
        a0 = initial_concentration
        f = fraction_remaining

        # Time unit conversion factor to seconds
        time_factors = {"s": 1.0, "min": 60.0, "h": 3600.0}
        tf = time_factors.get(time_unit.lower(), 1.0)
        # Convert k to per-time_unit
        # Note: user provides k in their desired time unit, so we use as-is

        abs_n = abs(n)

        if abs_n < 1e-6:
            # Zero order: [A] = [A]₀ - kt → when [A] = f·[A]₀: t = (1-f)·[A]₀/k
            t_half = (1 - f) * a0 / k
            formula = f"t = (1-{f})·[A]₀ / k"
            interp = (
                f"Zero-order reaction: rate is independent of concentration.\n"
                f"Half-life is proportional to initial concentration.\n"
                f"Total reaction time (complete): {a0 / k:.2f} {time_unit}"
            )

        elif abs(abs_n - 1.0) < 1e-6:
            # First order: ln([A]₀/[A]) = kt → t = ln(1/f)/k
            t_half = math.log(1.0 / f) / k
            formula = f"t = ln(1/{f}) / k"
            if abs(f - 0.5) < 1e-6:
                formula = "t₁/₂ = ln(2) / k ≈ 0.693/k"
            interp = (
                f"First-order reaction: half-life is independent of initial concentration.\n"
                f"This is characteristic of radioactive decay and many unimolecular reactions.\n"
                f"After 3 half-lives: {(1-0.5**3)*100:.1f}% decomposed\n"
                f"After 7 half-lives: {(1-0.5**7)*100:.2f}% decomposed (>99%)"
            )

        elif abs(abs_n - 2.0) < 1e-6:
            # Second order: 1/[A] - 1/[A]₀ = kt → t = (1/f - 1)/(k·[A]₀)
            t_half = (1.0 / f - 1.0) / (k * a0)
            formula = f"t = (1/{f} - 1) / (k·[A]₀)"
            if abs(f - 0.5) < 1e-6:
                formula = "t₁/₂ = 1 / (k·[A]₀)"
            interp = (
                f"Second-order reaction: half-life inversely proportional to [A]₀.\n"
                f"As reaction proceeds, each successive half-life doubles.\n"
                f"If [A]₀ halves, t₁/₂ doubles."
            )
        else:
            # General nth order (n ≠ 1):
            # ∫ d[A]/[A]^n from a0 to f*a0 = k*t
            # (1/(n-1)) · (1/(f·a0)^{n-1} - 1/a0^{n-1}) = k*t
            # t = [1/(k·(n-1))] · [(f·a0)^{1-n} - a0^{1-n}]
            #     = [a0^{1-n} / (k·(n-1))] · [f^{1-n} - 1]
            exponent = 1 - n
            t_half = (a0 ** exponent) / (k * (n - 1)) * (f ** (-exponent) - 1) if abs(n - 1) > 1e-10 else float("inf")
            formula = f"t = [A]₀^{1-n} / (k·(n-1)) · (f^{1-n} - 1)"
            interp = (
                f"{n}th-order reaction (n ≠ 1).\n"
                f"Half-life depends on both k and [A]₀.\n"
                f"For n>1: t₁/₂ decreases as [A]₀ decreases.\n"
                f"For 0<n<1: t₁/₂ increases as [A]₀ decreases."
            )

        pct_decomposed = round((1 - f) * 100, 1)

        analysis = (
            f"Half-Life Calculation:\n"
            f"Reaction order: n = {n}\n"
            f"Rate constant: k = {k}\n"
            f"Initial concentration: [A]₀ = {a0} {conc_unit}\n"
            f"Formula: {formula}\n"
            f"Result: t{'₁/₂' if abs(f-0.5)<1e-6 else ''} = {t_half:.4g} {time_unit}\n"
            f"Percent decomposed: {pct_decomposed}%\n\n"
            + interp
        )

        return {
            "half_life": round(t_half, 6),
            "half_life_unit": time_unit,
            "formula_used": formula,
            "analysis": analysis,
            "fractional_life": round(t_half, 6),
            "percent_decomposed": pct_decomposed,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            n = float(parts[0])
            k = float(parts[1])
            kwargs = {"reaction_order": n, "rate_constant": k}
            idx = 2
            if idx < len(parts):
                kwargs["initial_concentration"] = float(parts[idx]); idx += 1
            if idx < len(parts):
                kwargs["time_unit"] = parts[idx]; idx += 1
            if idx < len(parts):
                kwargs["conc_unit"] = parts[idx]; idx += 1
            if idx < len(parts):
                kwargs["fraction_remaining"] = float(parts[idx]); idx += 1
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
