import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EinsteinModel(BaseTool):
    """
    Einstein模型计算固体热容。
    Einstein模型：将固体中N个原子视为3N个独立的量子谐振子。
    Cv = 3Nk(θE/T)² * e^(θE/T) / (e^(θE/T) - 1)²
    高温极限：Cv → 3Nk (Dulong-Petit)
    低温极限：Cv 指数衰减趋于0（比Debye模型下降更快）
    """
    __version__      = "0.1.0"
    name             = "EinsteinModel"
    func_name        = "einstein_model"
    description      = "Calculate solid heat capacity using the Einstein model of independent quantum harmonic oscillators."
    implementation_description = "Implements the Einstein solid heat capacity: Cv = 3R(θ_E/T)^2 * exp(θ_E/T) / (exp(θ_E/T)-1)^2 per mole. Each atom treated as a 3D isotropic harmonic oscillator with characteristic Einstein temperature θ_E."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Einstein Model", "Heat Capacity", "Solid State", "Statistical Mechanics", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("einstein_temperature_k", "float", "N/A", "Characteristic Einstein temperature θ_E in Kelvin (K)."),
        ("temperature_k", "float", "N/A", "Temperature at which to evaluate heat capacity, in Kelvin (K)."),
        ("n_moles", "float", "1.0", "Amount of substance in moles. Default: 1.0."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated string: 'einstein_temperature_K temperature_K [n_moles]', e.g., '1200 300' or '200 100 2'."),
    ]

    output_sig       = [
        ("einstein_temperature_K", "float", "Einstein temperature θ_E used (K)."),
        ("temperature_K", "float", "Temperature used (K)."),
        ("T_over_ThetaE", "float", "Reduced temperature T/θ_E."),
        ("cv_molar_J_mol_K", "float", "Molar heat capacity at constant volume Cv,m in J/(mol·K)."),
        ("cv_total_J_K", "float", "Total heat capacity for given amount in J/K."),
        ("dulong_petit_J_mol_K", "float", "Dulong-Petit high-T limit (3R) in J/(mol·K)."),
        ("ratio_to_dulong_petit", "float", "Ratio of actual Cv to Dulong-Petit limit."),
        ("regime", "str", "Temperature regime: 'low_T', 'intermediate', or 'high_T'."),
        ("explanation", "str", "Brief explanation of the result."),
    ]

    examples         = [
        {
            "code_input": {
                "einstein_temperature_k": 1200.0,
                "temperature_k": 300.0,
                "n_moles": 1.0,
            },
            "text_input": {
                "input_params": "1200 300"
            },
            "output": {
                "einstein_temperature_K": 1200.0,
                "temperature_K": 300.0,
                "T_over_ThetaE": 0.250,
                "cv_molar_J_mol_K": 16.315,
                "cv_total_J_K": 16.315,
                "dulong_petit_J_mol_K": 24.943,
                "ratio_to_dulong_petit": 0.654,
                "regime": "low_T",
                "explanation": "At T << θ_E, Einstein model predicts exponential drop-off in Cv.",
            }
        },
        {
            "code_input": {
                "einstein_temperature_k": 200.0,
                "temperature_k": 400.0,
                "n_moles": 1.0,
            },
            "text_input": {
                "input_params": "200 400"
            },
            "output": {
                "einstein_temperature_K": 200.0,
                "temperature_K": 400.0,
                "T_over_ThetaE": 2.000,
                "cv_molar_J_mol_K": 24.612,
                "cv_total_J_K": 24.612,
                "dulong_petit_J_mol_K": 24.943,
                "ratio_to_dulong_petit": 0.987,
                "regime": "high_T",
                "explanation": "At T >> θ_E, Cv approaches the Dulong-Petit limit of 3R.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618  # J/(mol·K)

    def _run_base(self, einstein_temperature_k: float, temperature_k: float, n_moles: float = 1.0) -> dict:
        """Calculate molar heat capacity using Einstein model."""
        if einstein_temperature_k <= 0:
            raise ChemMCPError("Einstein temperature must be positive (in Kelvin).")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")
        if n_moles <= 0:
            raise ChemMCPError("Number of moles must be positive.")

        theta_E = einstein_temperature_k
        T = temperature_k
        R = self.R

        x = theta_E / T  # θ_E / T

        # Determine regime
        if x > 5.0:
            regime = "low_T"
        elif x < 0.5:
            regime = "high_T"
        else:
            regime = "intermediate"

        # Einstein formula: Cv = 3R * x^2 * e^x / (e^x - 1)^2
        try:
            ex = math.exp(x)
            denom = (ex - 1.0) ** 2
            if denom < 1e-30:
                cv_molar = 0.0
            else:
                cv_molar = 3.0 * R * (x * x) * ex / denom
        except OverflowError:
            cv_molar = 0.0

        dulong_petit = 3.0 * R
        ratio = cv_molar / dulong_petit if dulong_petit > 0 else 0.0

        return {
            "einstein_temperature_K": theta_E,
            "temperature_K": T,
            "T_over_ThetaE": round(T / theta_E, 6),
            "cv_molar_J_mol_K": round(cv_molar, 3),
            "cv_total_J_K": round(cv_molar * n_moles, 3),
            "dulong_petit_J_mol_K": round(dulong_petit, 3),
            "ratio_to_dulong_petit": round(ratio, 4),
            "regime": regime,
            "explanation": f"At T/θ_E = {T/theta_E:.3f}, Cv = {cv_molar:.3f} J/(mol·K), approaching Dulong-Petit limit as T increases.",
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse space-separated text input."""
        parts = input_params.strip().split()
        if len(parts) < 2:
            raise ChemMCPError(
                "Text input requires einstein_temperature and temperature. "
                "Format: 'einstein_temperature_K temperature_K [n_moles]'"
            )

        try:
            theta_E = float(parts[0])
            T = float(parts[1])
            n = float(parts[2]) if len(parts) > 2 else 1.0
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse numeric values from '{input_params}': {e}")

        return self._run_base(theta_E, T, n)
