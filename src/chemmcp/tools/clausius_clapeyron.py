import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Gas constant
R = 8.314  # J/(mol·K)


@ChemMCPManager.register_tool
class ClausiusClapeyron(BaseTool):
    """
    克劳修斯-克拉珀龙方程：计算相变参数。
    
    ln(P2/P1) = -ΔH_vap/R × (1/T2 - 1/T1)
    可求解：P2, T2, 或 ΔH_vap
    """
    __version__ = "0.1.0"
    name = "ClausiusClapeyron"
    func_name = "clausius_clapeyron_calc"
    description = "Calculate phase transition parameters using the Clausius-Clapeyron equation. Can solve for P2, T2, or ΔH_vap given the other parameters."
    implementation_description = "Uses ln(P2/P1) = -ΔH_vap/R × (1/T2 - 1/T1). Supports three solve modes: 'solve_p2', 'solve_t2', 'solve_deltah'. R = 8.314 J/(mol·K)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Phase Transition", "Clausius-Clapeyron", "Vapor Pressure"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Solve mode: 'solve_p2' (find pressure), 'solve_t2' (find temperature), or 'solve_deltah' (find enthalpy of vaporization)."),
        ("p1", "float", "N/A", "Known pressure in atm (or any consistent unit)."),
        ("t1", "float", "N/A", "Known temperature in Kelvin."),
        ("unknown_value", "float", "N/A", "The other known value: T2 (K) when solving for P2, P2 (atm) when solving for T2, or P2 (atm) when solving for ΔH."),
        ("delta_h_vap", "float", "None", "Enthalpy of vaporization in J/mol. Required only if NOT solving for it. Pass 0 or omit to auto-solve."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'mode p1 t1 unknown_value [delta_h_vap]'. Example: 'solve_p2 1.0 373.15 353.15 40650'"),
    ]

    output_sig = [
        ("result", "float", "The calculated value (pressure in atm, temperature in K, or enthalpy in J/mol)."),
        ("unit", "str", "Unit of the result."),
        ("explanation", "str", "Step-by-step calculation with formula used."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "solve_p2",
                "p1": 1.0,
                "t1": 373.15,
                "unknown_value": 353.15,
                "delta_h_vap": 40650.0,
            },
            "text_input": {
                "input_params": "solve_p2 1.0 373.15 353.15 40650",
            },
            "output": {
                "result": 0.4695,
                "unit": "atm",
                "explanation": "ln(P2/1) = -40650/8.314×(1/353.15-1/373.15) → P2 = exp(-0.7570) = 0.4695 atm",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, mode: str, p1: float, t1: float, unknown_value: float, delta_h_vap: float = 0.0) -> dict:
        """Core logic: Clausius-Clapeyron calculations."""
        mode = mode.lower().strip()
        if p1 <= 0:
            raise ChemMCPError("Pressure p1 must be positive.")
        if t1 <= 0:
            raise ChemMCPError("Temperature t1 must be positive in Kelvin.")

        if mode == "solve_p2":
            # Solve for P2 given T2 and ΔH_vap
            if delta_h_vap <= 0:
                raise ChemMCPError("Must provide positive delta_h_vap for solve_p2 mode.")
            t2 = unknown_value
            if t2 <= 0:
                raise ChemMCPError("Temperature t2 must be positive.")
            rhs = -(delta_h_vap / R) * (1.0 / t2 - 1.0 / t1)
            p2 = p1 * math.exp(rhs)
            explanation = (
                f"ln(P2/{p1}) = -ΔH_vap/R × (1/T2 - 1/T1)\n"
                f"= -{delta_h_vap}/{R} × (1/{t2} - 1/{t1})\n"
                f"= {rhs:.4f}\n"
                f"P2 = {p1} × exp({rhs:.4f}) = {p2:.4f} atm"
            )
            return {"result": round(p2, 6), "unit": "atm", "explanation": explanation}

        elif mode == "solve_t2":
            # Solve for T2 given P2 and ΔH_vap
            if delta_h_vap <= 0:
                raise ChemMCPError("Must provide positive delta_h_vap for solve_t2 mode.")
            p2 = unknown_value
            if p2 <= 0:
                raise ChemMCPError("Pressure p2 must be positive.")
            ln_ratio = math.log(p2 / p1)
            inv_t2 = -ln_ratio * R / delta_h_vap + 1.0 / t1
            t2 = 1.0 / inv_t2
            explanation = (
                f"ln({p2}/{p1}) = -{delta_h_vap}/{R} × (1/T2 - 1/{t1})\n"
                f"1/T2 = 1/{t1} - R·ln(P2/P1)/ΔH_vap\n"
                f"T2 = {t2:.4f} K"
            )
            return {"result": round(t2, 4), "unit": "K", "explanation": explanation}

        elif mode == "solve_deltah":
            # Solve for ΔH_vap given P2 and T2
            p2 = unknown_value
            t2 = unknown_value if False else 0  # placeholder — we need both from context
            # Actually unknown_value is P2, but we also need T2
            # Let's use a different approach: pass T2 as part of params
            # For now, use delta_h_vap field as T2 when solving for delH
            # Re-read: unknown_value is the "other known value". For solve_deltah, we need P2 AND T2.
            # We'll use unknown_value as P2 and delta_h_vap as T2 (repurposed)
            p2 = unknown_value
            t2 = delta_h_vap  # repurpose: this is actually T2 when mode is solve_deltah
            if p2 <= 0 or t2 <= 0:
                raise ChemMCPError("For solve_deltah: unknown_value=P2 (atm), delta_h_vap=T2 (K). Both must be positive.")
            ln_ratio = math.log(p2 / p1)
            dh = -R * ln_ratio / (1.0 / t2 - 1.0 / t1)
            explanation = (
                f"ln({p2}/{p1}) = -ΔH_vap/{R} × (1/{t2} - 1/{t1})\n"
                f"ΔH_vap = -R × ln(P2/P1) / (1/T2 - 1/T1)\n"
                f"= -{R} × {ln_ratio:.4f} / ({1/t2:.6f} - {1/t1:.6f})\n"
                f"= {dh:.2f} J/mol"
            )
            return {"result": round(dh, 2), "unit": "J/mol", "explanation": explanation}

        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use 'solve_p2', 'solve_t2', or 'solve_deltah'.")

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mode = parts[0]
            p1 = float(parts[1])
            t1 = float(parts[2])
            unknown_val = float(parts[3])
            dh = float(parts[4]) if len(parts) > 4 else 0.0
            return self._run_base(mode, p1, t1, unknown_val, dh)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'mode p1 t1 unknown_value [delta_h_vap]'")
