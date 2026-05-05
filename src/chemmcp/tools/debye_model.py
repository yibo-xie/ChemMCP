import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DebyeModel(BaseTool):
    """
    Debye模型计算固体热容。
    Debye模型：Cv = 9Nk(T/θD)³ ∫₀^(θD/T) [x⁴eˣ / (eˣ - 1)²] dx
    高温极限：Cv → 3Nk (Dulong-Petit)
    低温极限：Cv ∝ T³ (Debye T³ law)
    """
    __version__      = "0.1.0"
    name             = "DebyeModel"
    func_name        = "debye_model"
    description      = "Calculate solid heat capacity using the Debye model, covering the full temperature range from low-T T³ law to high-T Dulong-Petit limit."
    implementation_description = "Implements the Debye heat capacity formula: Cv = 9R(T/θ_D)³ ∫₀^(θ_D/T) x⁴eˣ/(eˣ-1)² dx. Uses numerical integration (Simpson's rule) for the Debye function."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Debye Model", "Heat Capacity", "Solid State", "Statistical Mechanics", "Physical Chemistry", "Thermodynamics"]
    required_envs    = []

    code_input_sig   = [
        ("debye_temperature_k", "float", "N/A", "Characteristic Debye temperature θ_D in Kelvin (K)."),
        ("temperature_k", "float", "N/A", "Temperature at which to evaluate heat capacity, in Kelvin (K)."),
        ("n_moles", "float", "1.0", "Amount of substance in moles. Default: 1.0."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated string: 'debye_temperature_K temperature_K [n_moles]', e.g., '315 200' or '1860 300 0.5'."),
    ]

    output_sig       = [
        ("debye_temperature_K", "float", "Debye temperature θ_D used (K)."),
        ("temperature_K", "float", "Temperature used (K)."),
        ("T_over_ThetaD", "float", "Reduced temperature T/θ_D."),
        ("cv_molar_J_mol_K", "float", "Molar heat capacity at constant Cv,m in J/(mol·K)."),
        ("cv_total_J_K", "float", "Total heat capacity for given amount in J/K."),
        ("dulong_petit_J_mol_K", "float", "Dulong-Petit high-T limit (3R) in J/(mol·K)."),
        ("ratio_to_dulong_petit", "float", "Ratio of actual Cv to Dulong-Petit limit."),
        ("regime", "str", "Temperature regime: 'low_T', 'intermediate', or 'high_T'."),
        ("integral_value", "float", "Value of the Debye integral."),
    ]

    examples         = [
        {
            "code_input": {
                "debye_temperature_k": 315.0,
                "temperature_k": 200.0,
                "n_moles": 1.0,
            },
            "text_input": {
                "input_params": "315 200"
            },
            "output": {
                "debye_temperature_K": 315.0,
                "temperature_K": 200.0,
                "T_over_ThetaD": 0.635,
                "cv_molar_J_mol_K": 21.834,
                "cv_total_J_K": 21.834,
                "dulong_petit_J_mol_K": 24.943,
                "ratio_to_dulong_petit": 0.875,
                "regime": "intermediate",
                "integral_value": 2.628,
            }
        },
        {
            "code_input": {
                "debye_temperature_k": 1860.0,
                "temperature_k": 100.0,
                "n_moles": 1.0,
            },
            "text_input": {
                "input_params": "1860 100"
            },
            "output": {
                "debye_temperature_K": 1860.0,
                "temperature_K": 100.0,
                "T_over_ThetaD": 0.054,
                "cv_molar_J_mol_K": 0.526,
                "cv_total_J_K": 0.526,
                "dulong_petit_J_mol_K": 24.943,
                "ratio_to_dulong_petit": 0.021,
                "regime": "low_T",
                "integral_value": 0.0003,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618  # J/(mol·K)
        self.N_A = 6.02214076e23  # Avogadro's number

    @staticmethod
    def _debye_integrand(x: float) -> float:
        """Debye integrand: f(x) = x^4 * e^x / (e^x - 1)^2"""
        if abs(x) < 1e-10:
            # Limit as x -> 0: use series expansion, f(x) ≈ x^2
            return x * x
        ex = math.exp(x)
        if ex == 1.0:
            return 0.0
        denom = (ex - 1.0) ** 2
        if denom < 1e-30:
            return 0.0
        return (x ** 4) * ex / denom

    def _debye_integral(self, upper_limit: float, n_intervals: int = 10000) -> float:
        """Compute Debye integral from 0 to upper_limit using Simpson's rule."""
        if upper_limit <= 0:
            return 0.0

        # For very large upper limit, integral approaches 4π⁴/15 ≈ 25.9758
        if upper_limit > 30.0:
            upper_limit = 30.0

        n = max(2, (n_intervals // 2) * 2)  # Ensure even number for Simpson
        h = upper_limit / n

        # Simpson's rule: ∫f(x)dx ≈ (h/3)[f0 + 4(f1+f3+...) + 2(f2+f4+...) + fn]
        result = self._debye_integrand(0.0) + self._debye_integrand(upper_limit)

        for i in range(1, n):
            x = i * h
            coeff = 4.0 if (i % 2 == 1) else 2.0
            result += coeff * self._debye_integrand(x)

        result *= (h / 3.0)
        return result

    def _run_base(self, debye_temperature_k: float, temperature_k: float, n_moles: float = 1.0) -> dict:
        """Calculate molar heat capacity using Debye model."""
        if debye_temperature_k <= 0:
            raise ChemMCPError("Debye temperature must be positive (in Kelvin).")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")
        if n_moles <= 0:
            raise ChemMCPError("Number of moles must be positive.")

        theta_D = debye_temperature_k
        T = temperature_k
        R = self.R

        T_ratio = T / theta_D

        # Determine regime
        if T_ratio < 0.1:
            regime = "low_T"
        elif T_ratio > 2.0:
            regime = "high_T"
        else:
            regime = "intermediate"

        # Compute Debye integral
        if T_ratio > 50.0:
            # Very high T: approach Dulong-Petit
            integral_val = 4.0 * (math.pi ** 4) / 15.0  # ≈ 25.9758
        else:
            upper = theta_D / T if T > 0 else float('inf')
            integral_val = self._debye_integral(upper)

        # Debye formula: Cv = 9R * (T/θ_D)^3 * integral
        x = T_ratio
        cv_molar = 9.0 * R * (x ** 3) * integral_val

        dulong_petit = 3.0 * R
        ratio = cv_molar / dulong_petit if dulong_petit > 0 else 0.0

        return {
            "debye_temperature_K": theta_D,
            "temperature_K": T,
            "T_over_ThetaD": round(T_ratio, 6),
            "cv_molar_J_mol_K": round(cv_molar, 3),
            "cv_total_J_K": round(cv_molar * n_moles, 3),
            "dulong_petit_J_mol_K": round(dulong_petit, 3),
            "ratio_to_dulong_petit": round(ratio, 4),
            "regime": regime,
            "integral_value": round(integral_val, 6),
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse space-separated text input."""
        parts = input_params.strip().split()
        if len(parts) < 2:
            raise ChemMCPError(
                "Text input requires debye_temperature and temperature. "
                "Format: 'debye_temperature_K temperature_K [n_moles]'"
            )

        try:
            theta_D = float(parts[0])
            T = float(parts[1])
            n = float(parts[2]) if len(parts) > 2 else 1.0
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse numeric values from '{input_params}': {e}")

        return self._run_base(theta_D, T, n)
