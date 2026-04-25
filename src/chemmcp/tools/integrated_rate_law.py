import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IntegratedRateLaw(BaseTool):
    """
    积分速率方程求解工具。
    支持零级、一级、二级反应的积分速率方程计算，包括浓度-时间关系、时间-浓度反解等。
    """
    __version__ = "0.1.0"
    name = "IntegratedRateLaw"
    func_name = "solve_integrated_rate_law"
    description = "Solve integrated rate laws for zero, first, and second order reactions: compute concentration at time t, or time to reach a given concentration."
    implementation_description = "Implements integrated rate equations: [A]=[A]₀-kt (zero), [A]=[A]₀·exp(-kt) (first), 1/[A]=1/[A]₀+kt (second). Supports forward and inverse calculations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Integrated Rate Law", "Chemical Kinetics", "Concentration"]
    required_envs    = []

    code_input_sig = [
        ("reaction_order", "int", "N/A", "Reaction order: 0, 1, or 2."),
        ("calculation_type", "str", "N/A", "What to calculate: 'concentration' ([A] at time t), 'time' (t to reach given [A]), 'fraction_remaining' (fraction at time t), or 'summary' (full table)."),
        ("rate_constant", "float", "N/A", "Rate constant k."),
        ("initial_concentration", "float", "N/A", "Initial concentration [A]₀."),
        # For concentration calculation:
        ("time_value", "float", "N/A", "Time value for calculation (for 'concentration' and 'fraction_remaining' types)."),
        # For time calculation:
        ("target_concentration", "float", "N/A", "Target concentration (for 'time' type)."),
        ("target_fraction", "float", "0.5", "Target fraction remaining (alternative to target_concentration)."),
        # Units
        ("time_unit", "str", "s", "Time unit: 's', 'min', 'h'."),
        ("conc_unit", "str", "M", "Concentration unit."),
        # For summary table:
        ("n_time_points", "int", "10", "Number of time points in summary table."),
        ("max_time", "float", "", "Maximum time for summary (default: ~5 half-lives)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated parameters string."),
    ]

    output_sig = [
        ("result", "float", "The calculated value (concentration, time, or fraction)."),
        ("result_unit", "str", "Unit of the result."),
        ("formula_used", "str", "The integrated rate law formula used."),
        ("analysis", "str", "Detailed step-by-step solution."),
        ("data_table", "list", "Summary table of concentration vs time (for 'summary' type)."),
    ]

    examples         = [
        {
            "code_input": {
                "reaction_order": 1,
                "calculation_type": 'concentration',
                "rate_constant": 0.00693,
                "initial_concentration": 1.0,
                "time_value": 200,
                "target_concentration": 0.0,
                "target_fraction": 0.5,
                "time_unit": 's',
                "conc_unit": 'M',
                "n_time_points": 10,
                "max_time": 500
            },
            "text_input": {
                "input_params": '1 concentration 0.00693 1.0 200 s'
            },
            "output": {
                "result": 0.25,
                "result_unit": 'M',
                "formula_used": '[A] = [A]_0 * exp(-kt)',
                "analysis": 'After 200s, 25% remains.',
                "data_table": []
            }
        },
        {
            "code_input": {
                "reaction_order": 2,
                "calculation_type": 'time',
                "rate_constant": 0.05,
                "initial_concentration": 0.2,
                "time_value": 0.0,
                "target_concentration": 0.05,
                "target_fraction": 0.25,
                "time_unit": 'min',
                "conc_unit": 'M',
                "n_time_points": 10,
                "max_time": 500
            },
            "text_input": {
                "input_params": '2 time 0.05 0.2 min 0.25'
            },
            "output": {
                "result": 300.0,
                "result_unit": 'min',
                "formula_used": 't = (1/[A] - 1/[A]_0) / k',
                "analysis": 'Second-order time to reach 25%.',
                "data_table": []
            }
        },
        {
            "code_input": {
                "reaction_order": 1,
                "calculation_type": 'summary',
                "rate_constant": 0.01,
                "initial_concentration": 1.0,
                "time_value": 0.0,
                "target_concentration": 0.0,
                "target_fraction": 0.5,
                "time_unit": 's',
                "conc_unit": 'M',
                "n_time_points": 5,
                "max_time": 500
            },
            "text_input": {
                "input_params": '1 summary 0.01 1.0 500 s 5'
            },
            "output": {
                "result": 0.0,
                "result_unit": '',
                "formula_used": '[A] = [A]_0 * exp(-kt)',
                "analysis": 'Summary table generated.',
                "data_table": [{'t': 0, 'conc': 1.0, 'fraction': 1.0}]
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _conc_at_t(self, n: int, k: float, a0: float, t: float) -> float:
        """Calculate [A] at time t."""
        if n == 0:
            return max(0, a0 - k * t)
        elif n == 1:
            return a0 * math.exp(-k * t)
        elif n == 2:
            denom = 1 + k * a0 * t
            if denom <= 0:
                return 0.0
            return a0 / denom
        else:
            raise ChemMCPError(f"Reaction order {n} not supported. Use 0, 1, or 2.")

    def _t_to_reach(self, n: int, k: float, a0: float, target: float) -> float:
        """Calculate time to reach target concentration."""
        if n == 0:
            if abs(k) < 1e-15:
                raise ChemMCPError("k cannot be zero for zero-order.")
            return (a0 - target) / k
        elif n == 1:
            if target <= 0:
                raise ChemMCPError("Target must be positive for first-order.")
            return -math.log(target / a0) / k
        elif n == 2:
            if target <= 0:
                raise ChemMCPError("Target must be positive for second-order.")
            return (1.0 / target - 1.0 / a0) / k
        else:
            raise ChemMCPError(f"Reaction order {n} not supported.")

    def _run_base(
        self,
        reaction_order: int,
        calculation_type: str,
        rate_constant: float,
        initial_concentration: float,
        time_value: float = 0.0,
        target_concentration: float = 0.0,
        target_fraction: float = 0.5,
        time_unit: str = "s",
        conc_unit: str = "M",
        n_time_points: int = 10,
        max_time: float = 0.0,
    ) -> dict:
        n = reaction_order
        calc_type = calculation_type.lower().strip()
        k = rate_constant
        a0 = initial_concentration

        if n not in (0, 1, 2):
            raise ChemMCPError(f"Reaction order must be 0, 1, or 2. Got {n}.")
        if k < 0:
            raise ChemMCPError("Rate constant must be non-negative.")

        # Formula strings
        formulas = {
            0: "[A] = [A]₀ − k·t",
            1: "[A] = [A]₀ · exp(−kt)",
            2: "[A] = [A]₀ / (1 + k·[A]₀·t)",
        }
        formula_str = formulas[n]

        if calc_type == "concentration":
            result_val = self._conc_at_t(n, k, a0, time_value)
            frac = result_val / a0 if a0 != 0 else 0
            result_unit = conc_unit
            analysis = (
                f"Integrated Rate Law (order {n}):\n"
                f"Formula: {formula_str}\n"
                f"[A]₀ = {a0} {conc_unit}, k = {k}, t = {time_value} {time_unit}\n"
                f"[A]({time_value} {time_unit}) = {result_val:.6g} {conc_unit}\n"
                f"Fraction remaining: {frac:.4f} ({frac*100:.2f}%)\n"
                f"Decomposed: {(1-frac)*100:.2f}%"
            )

        elif calc_type == "time":
            tgt = target_concentration if target_concentration > 0 else target_fraction * a0
            result_val = self._t_to_reach(n, k, a0, tgt)
            result_unit = time_unit
            frac = tgt / a0
            analysis = (
                f"Time to reach target concentration:\n"
                f"Formula: rearranged from {formula_str}\n"
                f"[A]₀ = {a0} {conc_unit}, k = {k}\n"
                f"Target: [A] = {tgt:.6g} {conc_unit} ({frac*100:.2f}% remaining)\n"
                f"t = {result_val:.4g} {time_unit}"
            )

        elif calc_type == "fraction_remaining":
            result_val = self._conc_at_t(n, k, a0, time_value) / a0
            result_unit = "dimensionless"
            analysis = (
                f"Fraction remaining at t = {time_value} {time_unit}:\n"
                f"[A]/[A]₀ = {result_val:.6g}"
            )

        elif calc_type == "summary":
            # Generate summary table
            if max_time <= 0:
                # Estimate ~5 half-lives
                if n == 0:
                    max_time = a0 / k if k > 0 else 500
                elif n == 1:
                    max_time = 5 * math.log(2) / k if k > 0 else 500
                elif n == 2:
                    max_time = 5.0 / (k * a0) if k > 0 and a0 > 0 else 500

            table = []
            dt = max_time / n_time_points
            for i in range(n_time_points + 1):
                t = i * dt
                c = self._conc_at_t(n, k, a0, t)
                fr = c / a0 if a0 != 0 else 0
                table.append({
                    "t": round(t, 4),
                    "conc": round(c, 6),
                    "fraction": round(fr, 6),
                    "decomposed_pct": round((1 - fr) * 100, 2),
                })

            result_val = 0.0
            result_unit = ""
            analysis = (
                f"Summary table ({n_time_points+1} points, order={n}, k={k}, [A]₀={a0}):\n"
                + "\n".join([f"  t={row['t']:.2f} {time_unit}: [A]={row['conc']:.4g} ({row['fraction']:.2%}, {row['decomposed_pct']:.1f}% decomposed)" for row in table])
            )
            data_table = table
            return {
                "result": result_val,
                "result_unit": result_unit,
                "formula_used": formula_str,
                "analysis": analysis,
                "data_table": data_table[:10],
            }

        else:
            raise ChemMCPError(f"Unsupported calculation_type: '{calc_type}'. Use 'concentration', 'time', 'fraction_remaining', or 'summary'.")

        return {
            "result": round(result_val, 6),
            "result_unit": result_unit,
            "formula_used": formula_str,
            "analysis": analysis,
            "data_table": [],
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            n = int(parts[0])
            calc_type = parts[1]
            kwargs = {"reaction_order": n, "calculation_type": calc_type}
            idx = 2
            if idx < len(parts):
                kwargs["rate_constant"] = float(parts[idx]); idx += 1
            if idx < len(parts):
                kwargs["initial_concentration"] = float(parts[idx]); idx += 1
            if calc_type == "concentration":
                if idx < len(parts):
                    kwargs["time_value"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["time_unit"] = parts[idx]; idx += 1
            elif calc_type == "time":
                if idx < len(parts):
                    try:
                        kwargs["target_concentration"] = float(parts[idx])
                    except ValueError:
                        kwargs["target_fraction"] = float(parts[idx])
                    idx += 1
                if idx < len(parts):
                    kwargs["time_unit"] = parts[idx]; idx += 1
            elif calc_type == "summary":
                if idx < len(parts):
                    kwargs["max_time"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["time_unit"] = parts[idx]; idx += 1
                if idx < len(parts):
                    kwargs["n_time_points"] = int(parts[idx]); idx += 1
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
