import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
NA = 6.02214076e23         # mol⁻¹


@ChemMCPManager.register_tool
class BeerLambertCalculator(BaseTool):
    """
    Beer-Lambert定律计算工具。
    计算吸光度、浓度、透射率、摩尔消光系数等参数，支持多组分体系。
    A = ε · b · c, T = 10^(-A)
    """
    __version__ = "0.1.0"
    name = "BeerLambertCalculator"
    func_name = "calculate_beer_lambert"
    description = "Calculate absorbance, concentration, transmittance, and molar absorptivity using the Beer-Lambert law. Supports single-component and multi-component systems."
    implementation_description = "Implements Beer-Lambert law (A = εbc) with full variable solving: given any two of {A, ε, b, c}, compute the third. Also handles transmittance (T), percent transmittance (%T), and multi-component additive absorbance."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "UV-Vis", "Analytical Chemistry", "Beer-Lambert", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Calculation mode: 'absorbance', 'concentration', 'epsilon', 'pathlength', 'transmittance', 'multi_component', or 'dilution'."),
        ("absorbance_A", "float", "0", "Absorbance (unitless)."),
        ("epsilon_M_cm", "float", "0", "Molar extinction coefficient ε in M⁻¹·cm⁻¹."),
        ("pathlength_cm", "float", "1.0", "Path length b in cm (standard cuvette: 1 cm)."),
        ("concentration_M", "float", "0", "Concentration c in mol/L (M)."),
        ("transmittance_T", "float", "0", "Transmittance T (0 to 1)."),
        ("percent_transmittance", "float", "0", "Percent transmittance %T (0 to 100)."),
        ("components", "list", "[]", "List of component dicts for multi-component mode: [{'epsilon': x, 'c': y}, ...]."),
        ("initial_volume_mL", "float", "0", "Initial volume for dilution calculations."),
        ("final_volume_mL", "float", "0", "Final volume for dilution calculations."),
        ("wavelength_nm", "float", "0", "Optional wavelength context (nm)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Mode-specific parameters. E.g., 'absorbance 100000 1.0 0.001' or 'transmittance_from_A 0.5'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing calculated values, intermediate steps, and practical interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "absorbance",
                "absorbance_A": 0,
                "epsilon_M_cm": 100000,
                "pathlength_cm": 1.0,
                "concentration_M": 0.001,
                "transmittance_T": 0,
                "percent_transmittance": 0,
                "components": [],
                "initial_volume_mL": 0,
                "final_volume_mL": 0,
                "wavelength_nm": 280,
            },
            "text_input": {
                "input_params": "absorbance 100000 1.0 0.001 280",
            },
            "output": {
                "result": {
                    "mode": "absorbance",
                    "given": {"epsilon_M_cm": 100000, "pathlength_cm": 1.0, "concentration_M": 0.001},
                    "absorbance_A": 100.0,
                    "note": "A > 2 is outside linear range; sample needs dilution.",
                    "transmittance_T": 10**(-100),
                    "percent_transmittance": 10**(-98),
                    "formula": "A = ε × b × c = 100000 × 1.0 × 0.001 = 100.0",
                    "valid_range_check": "WARNING: A = 100 >> 2. Dilute sample and remeasure.",
                }
            }
        },
        {
            "code_input": {
                "mode": "concentration",
                "absorbance_A": 0.654,
                "epsilon_M_cm": 15000,
                "pathlength_cm": 1.0,
                "concentration_M": 0,
                "transmittance_T": 0,
                "percent_transmittance": 0,
                "components": [],
                "initial_volume_mL": 0,
                "final_volume_mL": 0,
                "wavelength_nm": 260,
            },
            "text_input": {
                "input_params": "concentration 0.654 15000 1.0 260",
            },
            "output": {
                "result": {
                    "mode": "concentration",
                    "given": {"absorbance_A": 0.654, "epsilon_M_cm": 15000, "pathlength_cm": 1.0},
                    "concentration_M": round(0.654 / (15000 * 1.0), 8),
                    "transmittance_T": round(10 ** (-0.654), 6),
                    "percent_transmittance": round(10 ** (-0.654) * 100, 3),
                    "formula": "c = A / (ε × b) = 0.654 / (15000 × 1.0)",
                    "range_check": "A = 0.654 is within optimal range (0.2-0.8). Good measurement.",
                }
            }
        },
        {
            "code_input": {
                "mode": "transmittance",
                "absorbance_A": 0,
                "epsilon_M_cm": 0,
                "pathlength_cm": 0,
                "concentration_M": 0,
                "transmittance_T": 0.35,
                "percent_transmittance": 0,
                "components": [],
                "initial_volume_mL": 0,
                "final_volume_mL": 0,
                "wavelength_nm": 0,
            },
            "text_input": {
                "input_params": "transmittance 0.35",
            },
            "output": {
                "result": {
                    "mode": "transmittance → absorbance",
                    "transmittance_T": 0.35,
                    "percent_transmittance": 35.0,
                    "absorbance_A": round(-math.log10(0.35), 4),
                    "formula": "A = -log₁₀(T) = -log₁₀(0.35)",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, mode: str, absorbance_A: float = 0.0, epsilon_M_cm: float = 0.0,
                  pathlength_cm: float = 1.0, concentration_M: float = 0.0,
                  transmittance_T: float = 0.0, percent_transmittance: float = 0.0,
                  components: List[dict] = None, initial_volume_mL: float = 0.0,
                  final_volume_mL: float = 0.0, wavelength_nm: float = 0.0) -> dict:
        """Core logic."""
        if components is None:
            components = []

        m = mode.lower().replace("-", "_").replace(" ", "_")

        if m == "absorbance":
            return self._calc_absorbance(epsilon_M_cm, pathlength_cm, concentration_M, wavelength_nm)
        elif m == "concentration":
            return self._calc_concentration(absorbance_A, epsilon_M_cm, pathlength_cm, wavelength_nm)
        elif m == "epsilon":
            return self._calc_epsilon(absorbance_A, pathlength_cm, concentration_M, wavelength_nm)
        elif m == "pathlength":
            return self._calc_pathlength(absorbance_A, epsilon_M_cm, concentration_M)
        elif m in ("transmittance", "transmittance_from_a"):
            return self._calc_from_transmittance(transmittance_T, percent_transmittance)
        elif m == "multi_component":
            return self._calc_multi_component(components, pathlength_cm)
        elif m == "dilution":
            return self._calc_dilution(absorbance_A, initial_volume_mL, final_volume_mL)
        else:
            raise ChemMCPError(
                f"Unknown mode '{mode}'. Use: 'absorbance', 'concentration', 'epsilon', "
                "'pathlength', 'transmittance', 'multi_component', or 'dilution'."
            )

    def _calc_absorbance(self, eps: float, b: float, c: float, wl: float) -> dict:
        """A = ε · b · c"""
        if eps <= 0:
            raise ChemMCPError("Molar extinction coefficient must be positive.")
        if b <= 0:
            raise ChemMCPError("Path length must be positive.")
        if c < 0:
            raise ChemMCPError("Concentration cannot be negative.")

        A = eps * b * c
        T = 10 ** (-A)
        pct_T = T * 100

        result = {
            "mode": "absorbance",
            "given": {"epsilon_M_cm": eps, "pathlength_cm": b, "concentration_M": c},
            "absorbance_A": round(A, 6),
            "transmittance_T": round(T, 10),
            "percent_transmittance": round(pct_T, 6),
            "formula": f"A = ε × b × c = {eps} × {b} × {c} = {round(A, 4)}",
        }

        if wl:
            result["wavelength_nm"] = wl

        result["valid_range_check"] = self._check_range(A)
        return {"result": result}

    def _calc_concentration(self, A: float, eps: float, b: float, wl: float) -> dict:
        """c = A / (ε · b)"""
        if A < 0:
            raise ChemMCPError("Absorbance cannot be negative.")
        if eps <= 0:
            raise ChemMCPError("Molar extinction coefficient must be positive.")
        if b <= 0:
            raise ChemMCPError("Path length must be positive.")

        c = A / (eps * b)
        T = 10 ** (-A)
        pct_T = T * 100

        result = {
            "mode": "concentration",
            "given": {"absorbance_A": A, "epsilon_M_cm": eps, "pathlength_cm": b},
            "concentration_M": round(c, 10),
            "concentration_uM": round(c * 1e6, 4),
            "concentration_nM": round(c * 1e9, 4),
            "transmittance_T": round(T, 10),
            "percent_transmittance": round(pct_T, 4),
            "formula": f"c = A / (ε × b) = {A} / ({eps} × {b}) = {round(c, 8)} M",
        }

        if wl:
            result["wavelength_nm"] = wl

        result["range_check"] = self._check_range(A)
        return {"result": result}

    def _calc_epsilon(self, A: float, b: float, c: float, wl: float) -> dict:
        """ε = A / (b · c)"""
        if A < 0:
            raise ChemMCPError("Absorbance cannot be negative.")
        if b <= 0:
            raise ChemMCPError("Path length must be positive.")
        if c <= 0:
            raise ChemMCPError("Concentration must be positive.")

        eps = A / (b * c)
        T = 10 ** (-A)

        result = {
            "mode": "epsilon",
            "given": {"absorbance_A": A, "pathlength_cm": b, "concentration_M": c},
            "epsilon_M_cm": round(eps, 4),
            "log10_epsilon": round(math.log10(eps), 3) if eps > 0 else None,
            "transmittance_T": round(T, 10),
            "formula": f"ε = A / (b × c) = {A} / ({b} × {c}) = {round(eps, 2)} M⁻¹cm⁻¹",
        }
        if wl:
            result["wavelength_nm"] = wl
        result["classification"] = self._classify_epsilon(eps)
        return {"result": result}

    def _calc_pathlength(self, A: float, eps: float, c: float) -> dict:
        """b = A / (ε · c)"""
        if A < 0:
            raise ChemMCPError("Absorbance cannot be negative.")
        if eps <= 0 or c <= 0:
            raise ChemMCPError("ε and c must be positive.")

        b = A / (eps * c)
        return {"result": {
            "mode": "pathlength",
            "given": {"absorbance_A": A, "epsilon_M_cm": eps, "concentration_M": c},
            "pathlength_cm": round(b, 6),
            "pathlength_mm": round(b * 10, 4),
            "formula": f"b = A / (ε × c) = {A} / ({eps} × {c}) = {round(b, 4)} cm",
        }}

    def _calc_from_transmittance(self, T: float, pct_T: float) -> dict:
        """Convert between T, %T, and A."""
        if pct_T > 0:
            T = pct_T / 100.0
        if T <= 0 or T > 1:
            raise ChemMCPError(f"Transmittance must be in range (0, 1] (or %T in (0, 100]). Got T={T}.")

        A = -math.log10(T)
        return {"result": {
            "mode": "transmittance ↔ absorbance",
            "transmittance_T": round(T, 10),
            "percent_transmittance": round(T * 100, 6),
            "absorbance_A": round(A, 6),
            "formula": f"A = -log₁₀(T) = -log₁₀({T:.6f}) = {round(A, 4)}",
            "light_absorbed_pct": round((1 - T) * 100, 4),
            "range_check": self._check_range(A),
        }}

    def _calc_multi_component(self, components: List[dict], b: float) -> dict:
        """A_total = Σ(ε_i · b · c_i)"""
        if not components:
            raise ChemMCPError("At least one component required for multi_component mode.")

        total_A = 0
        comp_results = []
        for i, comp in enumerate(components):
            eps = comp.get("epsilon", 0)
            conc = comp.get("c", 0) or comp.get("concentration", 0)
            name = comp.get("name", f"Component_{i+1}")
            Ai = eps * b * conc
            total_A += Ai
            comp_results.append({
                "name": name,
                "epsilon_M_cm": eps,
                "concentration_M": conc,
                "absorbance_i": round(Ai, 6),
                "fraction_of_total": round(Ai / total_A * 100, 2) if total_A > 0 else 0,
            })

        T_total = 10 ** (-total_A)

        return {"result": {
            "mode": "multi_component",
            "pathlength_cm": b,
            "components": comp_results,
            "total_absorbance_A": round(total_A, 6),
            "total_transmittance_T": round(T_total, 10),
            "total_percent_transmittance": round(T_total * 100, 4),
            "formula": "A_total = Σ(εᵢ × b × cᵢ)",
            "additivity_note": "Absorbances are additive at each wavelength.",
            "range_check": self._check_range(total_A),
        }}

    def _calc_dilution(self, A_initial: float, V_initial: float, V_final: float) -> dict:
        """Calculate absorbance after dilution: A_final = A_initial × (V_initial/V_final)"""
        if V_initial <= 0:
            raise ChemMCPError("Initial volume must be positive.")
        if V_final <= 0:
            raise ChemMCPError("Final volume must be positive.")
        if V_final < V_initial:
            raise ChemMCPError("Final volume must be >= initial volume for dilution.")

        dilution_factor = V_initial / V_final
        A_final = A_initial * dilution_factor
        df_inverse = V_final / V_initial

        return {"result": {
            "mode": "dilution",
            "initial_absorbance_A": A_initial,
            "initial_volume_mL": V_initial,
            "final_volume_mL": V_final,
            "dilution_factor": round(dilution_factor, 6),
            "fold_dilution": f"1:{round(df_inverse, 1)}",
            "final_absorbance_A": round(A_final, 6),
            "formula": f"A_final = A_initial × (V_i/V_f) = {A_initial} × ({V_initial}/{V_final}) = {round(A_final, 4)}",
            "range_check": self._check_range(A_final),
        }}

    @staticmethod
    def _check_range(A: float) -> str:
        """Check if absorbance is in optimal measurement range."""
        if A < 0:
            return "ERROR: Negative absorbance is physically impossible."
        elif A < 0.01:
            return "Below detection limit — too dilute."
        elif A < 0.2:
            return "Low absorbance — consider using more concentrated sample."
        elif A <= 0.8:
            return "✓ Optimal range (0.2-0.8) for accurate measurements."
        elif A <= 2.0:
            return "Acceptable but approaching upper limit of linearity; some error expected."
        elif A <= 3.0:
            return "⚠ High absorbance — significant deviation from linearity. Dilute sample."
        else:
            return "✗ WAY too high — detector saturation. Must dilute sample."

    @staticmethod
    def _classify_epsilon(eps: float) -> str:
        """Classify molar absorptivity strength."""
        if eps < 10:
            return "n→π* transition (forbidden/weak)"
        elif eps < 100:
            return "Weak absorption (forbidden transition)"
        elif eps < 1000:
            return "Moderate (partially allowed)"
        elif eps < 10000:
            return "Strong (allowed π→π* transition)"
        else:
            return f"Very strong (ε={eps:.0e}); typical of allowed π→π* or charge-transfer transitions"

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            mode = parts[0]

            if mode == "absorbance":
                eps = float(parts[1]) if len(parts) > 1 else 0
                b = float(parts[2]) if len(parts) > 2 else 1.0
                c = float(parts[3]) if len(parts) > 3 else 0
                wl = float(parts[4]) if len(parts) > 4 else 0
                return self._run_base(mode, 0, eps, b, c, 0, 0, [], 0, 0, wl)
            elif mode == "concentration":
                A = float(parts[1]) if len(parts) > 1 else 0
                eps = float(parts[2]) if len(parts) > 2 else 0
                b = float(parts[3]) if len(parts) > 3 else 1.0
                wl = float(parts[4]) if len(parts) > 4 else 0
                return self._run_base(mode, A, eps, b, 0, 0, 0, [], 0, 0, wl)
            elif mode in ("transmittance", "t"):
                T = float(parts[1]) if len(parts) > 1 else 0
                return self._run_base("transmittance", 0, 0, 0, 0, T, 0)
            elif mode == "epsilon":
                A = float(parts[1]) if len(parts) > 1 else 0
                b = float(parts[2]) if len(parts) > 2 else 1.0
                c = float(parts[3]) if len(parts) > 3 else 0
                return self._run_base(mode, A, 0, b, c)
            elif mode == "dilution":
                A = float(parts[1]) if len(parts) > 1 else 0
                Vi = float(parts[2]) if len(parts) > 2 else 0
                Vf = float(parts[3]) if len(parts) > 3 else 0
                return self._run_base(mode, A, 0, 0, 0, 0, 0, [], Vi, Vf)
            else:
                raise ValueError(f"Unknown mode: {mode}")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
