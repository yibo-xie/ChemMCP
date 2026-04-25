import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EnzymeKinetics(BaseTool):
    """
    Michaelis-Menten 酶动力学参数求解工具。
    使用三种线性化方法（Lineweaver-Burk、Eadie-Hofstee、Hanes-Woolf）拟合 Vmax 和 Km。
    """
    __version__ = "0.1.0"
    name = "EnzymeKinetics"
    func_name = "enzyme_kinetics_fit"
    description = "Solve Michaelis-Menten enzyme kinetics parameters (Vmax, Km) using multiple linearization methods."
    implementation_description = "Implements Lineweaver-Burk (double reciprocal), Eadie-Hofstee, and Hanes-Woolf linear regression methods to estimate Vmax and Km from [S] vs v data. Also computes kcat if enzyme concentration is provided."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Enzyme Kinetics", "Michaelis-Menten", "Biochemistry", "Linear Regression"]
    required_envs = []

    code_input_sig = [
        ("substrate_concentrations_S", "list", "N/A", "List of substrate concentrations [S] (same units, e.g., mM)."),
        ("reaction_velocities_v", "list", "N/A", "List of reaction velocities v (same units, e.g., mM/s), corresponding to each [S]."),
        ("enzyme_concentration_E0", "float", "N/A", "Total enzyme concentration (optional, for kcat calculation). Provide as 0 or omit if unknown."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'S1,S2,S3,... v1,v2,v3,... [E0]' where S=substrate conc, v=velocity, E0=enzyme conc (optional)."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with Vmax, Km from each method, best fit, kcat (if E0 given), and R² values."),
    ]

    examples = [
        {
            "code_input": {
                "substrate_concentrations_S": [0.5, 1.0, 2.0, 4.0, 8.0],
                "reaction_velocities_v": [0.21, 0.35, 0.53, 0.71, 0.85],
                "enzyme_concentration_E0": None,
            },
            "text_input": {
                "input_params": "0.5,1,2,4,8 0.21,0.35,0.53,0.71,0.85",
            },
            "output": {
                "result": {
                    "Vmax_best": 1.02,
                    "Km_best": 2.05,
                    "best_method": "Hanes-Woolf",
                    "kcat": None,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _linear_regression(self, x_list, y_list):
        """Simple linear regression: y = slope * x + intercept."""
        n = len(x_list)
        if n < 2:
            raise ChemMCPError("Need at least 2 data points for regression.")
        sum_x = sum(x_list)
        sum_y = sum(y_list)
        sum_xy = sum(xi * yi for xi, yi in zip(x_list, y_list))
        sum_x2 = sum(xi ** 2 for xi in x_list)
        sum_y2 = sum(yi ** 2 for yi in y_list)

        denom = n * sum_x2 - sum_x ** 2
        if abs(denom) < 1e-15:
            raise ChemMCPError("Cannot fit: x values are nearly identical.")

        slope = (n * sum_xy - sum_x * sum_y) / denom
        intercept = (sum_y - slope * sum_x) / n

        # R² calculation
        y_mean = sum_y / n
        ss_tot = sum((yi - y_mean) ** 2 for yi in y_list)
        ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x_list, y_list))
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        return {"slope": slope, "intercept": intercept, "r_squared": round(r_squared, 6)}

    def _lineweaver_burk(self, S_data, v_data):
        """1/v = (Km/Vmax)*(1/S) + 1/Vmax"""
        inv_s = [1.0 / s if s != 0 else float('inf') for s in S_data]
        inv_v = [1.0 / v if v != 0 else float('inf') for v in v_data]
        # Filter out infinities
        filtered = [(x, y) for x, y in zip(inv_s, inv_v) if x != float('inf') and y != float('inf')]
        if len(filtered) < 2:
            return {"Vmax": None, "Km": None, "r_squared": 0}
        xs, ys = zip(*filtered)
        result = self._linear_regression(list(xs), list(ys))
        slope = result["slope"]
        intercept = result["intercept"]
        if intercept == 0:
            return {"Vmax": None, "Km": None, "r_squared": result["r_squared"]}
        Vmax = 1.0 / intercept
        Km = slope * Vmax
        return {"Vmax": round(Vmax, 6), "Km": round(Km, 6), "r_squared": result["r_squared"]}

    def _eadie_hofstee(self, S_data, v_data):
        """v = -Km*(v/S) + Vmax"""
        v_over_s = [v / s if s != 0 else 0 for v, s in zip(v_data, S_data)]
        result = self._linear_regression(v_over_s, list(v_data))
        slope = result["slope"]
        intercept = result["intercept"]
        Vmax = intercept
        Km = -slope
        return {"Vmax": round(Vmax, 6), "Km": round(max(Km, 0), 6), "r_squared": result["r_squared"]}

    def _hanes_woolf(self, S_data, v_data):
        """S/v = (1/Vmax)*S + Km/Vmax"""
        s_over_v = [s / v if v != 0 else 0 for s, v in zip(S_data, v_data)]
        result = self._linear_regression(list(S_data), s_over_v)
        slope = result["slope"]
        intercept = result["intercept"]
        if slope == 0 or slope == float('inf'):
            return {"Vmax": None, "Km": None, "r_squared": result["r_squared"]}
        Vmax = 1.0 / slope
        Km = intercept * Vmax
        return {"Vmax": round(Vmax, 6), "Km": round(Km, 6), "r_squared": result["r_squared"]}

    def _run_base(self, substrate_concentrations_S: list, reaction_velocities_v: list,
                  enzyme_concentration_E0: float = None) -> dict:
        if len(substrate_concentrations_S) != len(reaction_velocities_v):
            raise ChemMCPError("S and v lists must have the same length.")
        if len(substrate_concentrations_S) < 3:
            raise ChemMCPError("Need at least 3 data points for reliable fitting.")
        if any(s <= 0 for s in substrate_concentrations_S):
            raise ChemMCPError("All substrate concentrations must be positive.")
        if any(v <= 0 for v in reaction_velocities_v):
            raise ChemMCPError("All velocities must be positive.")

        lb = self._lineweaver_burk(substrate_concentrations_S, reaction_velocities_v)
        eh = self._eadie_hofstee(substrate_concentrations_S, reaction_velocities_v)
        hw = self._hanes_woolf(substrate_concentrations_S, reaction_velocities_v)

        methods = {
            "Lineweaver-Burk": lb,
            "Eadie-Hofstee": eh,
            "Hanes-Woolf": hw,
        }

        # Pick best by R² (prefer Hanes-Woolf for statistical robustness, then Eadie-Hofstee)
        valid_methods = {name: m for name, m in methods.items()
                         if m.get("Vmax") is not None and m.get("Vmax") > 0}
        if not valid_methods:
            raise ChemMCPError("Could not fit data with any method.")

        best_name = max(valid_methods.keys(), key=lambda n: valid_methods[n]["r_squared"])
        best = valid_methods[best_name]

        kcat_val = None
        if enzyme_concentration_E0 and enzyme_concentration_E0 > 0 and best["Vmax"]:
            kcat_val = round(best["Vmax"] / enzyme_concentration_E0, 6)

        logger.info(f"EnzymeKinetics: best={best_name}, Vmax={best['Vmax']}, Km={best['Km']}")

        return {
            "Vmax_best": best["Vmax"],
            "Km_best": best["Km"],
            "best_method": best_name,
            "kcat": kcat_val,
            "Lineweaver_Burk": lb,
            "Eadie_Hofstee": eh,
            "Hanes_Woolf": hw,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            S_list = [float(x) for x in parts[0].split(",")]
            v_list = [float(x) for x in parts[1].split(",")]
            E0 = float(parts[2]) if len(parts) > 2 else None
            return self._run_base(S_list, v_list, E0)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'S1,S2,... v1,v2,... [E0]'")
