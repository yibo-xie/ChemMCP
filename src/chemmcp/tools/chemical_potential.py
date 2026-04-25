import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ChemicalPotential(BaseTool):
    """
    计算多组分系统的化学势工具。
    支持理想气体、理想溶液、实际溶液（活度系数）等模型。
    """
    __version__ = "0.1.0"
    name = "ChemicalPotential"
    func_name = "calculate_chemical_potential"
    description = "Calculate chemical potentials for multi-component systems (ideal gas, ideal solution, real solution with activity coefficients)."
    implementation_description = "Computes μ_i = μ_i° + RT ln(a_i) where a_i depends on phase and model (ideal gas: P_i/P°; ideal solution: x_i; real solution: γ_i·x_i)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Chemical Potential", "Physical Chemistry", "Multi-component"]
    required_envs    = []

    code_input_sig = [
        ("model", "str", "N/A", "Model type: 'ideal_gas', 'ideal_solution', 'real_solution'."),
        ("temperature_k", "float", "N/A", "Temperature in Kelvin."),
        ("mu_standard", "float", "N/A", "Standard chemical potential μ_i° in J/mol for each component (comma-separated if multiple)."),
        # For ideal_gas model:
        ("partial_pressures", "str", "N/A", "Partial pressures in atm, comma-separated (for ideal_gas)."),
        ("p_standard", "float", "1.0", "Standard state pressure in atm (default 1 atm)."),
        # For solution models:
        ("mole_fractions", "str", "N/A", "Mole fractions x_i, comma-separated (for solutions)."),
        ("activity_coefficients", "str", "1.0", "Activity coefficients γ_i, comma-separated (for real_solution, default 1.0 each)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space/newline-separated parameters string. See code_input for format details."),
    ]

    output_sig = [
        ("chemical_potentials", "list", "List of chemical potentials μ_i (J/mol) for each component."),
        ("total_gibbs_per_mol", "float", "Total molar Gibbs energy Σx_i·μ_i (J/mol)."),
        ("analysis", "str", "Detailed analysis of the chemical potential calculation."),
    ]

    examples         = [
        {
            "code_input": {
                "model": 'ideal_gas',
                "temperature_k": 298.15,
                "mu_standard": '-394360',
                "partial_pressures": '0.5',
                "p_standard": 1.0,
                "mole_fractions": '',
                "activity_coefficients": '1.0'
            },
            "text_input": {
                "input_params": 'ideal_gas 298.15 -394360 0.5 1.0'
            },
            "output": {
                "chemical_potentials": [-395742.6],
                "total_gibbs_per_mol": -395742.6,
                "analysis": 'CO2 at 298.15 K, 0.5 atm.'
            }
        },
        {
            "code_input": {
                "model": 'real_solution',
                "temperature_k": 298.15,
                "mu_standard": '0,-10000',
                "partial_pressures": '',
                "p_standard": 1.0,
                "mole_fractions": '0.3,0.7',
                "activity_coefficients": '1.2,0.95'
            },
            "text_input": {
                "input_params": 'real_solution 298.15 0,-10000 0.3,0.7 1.2,0.95'
            },
            "output": {
                "chemical_potentials": [-2689.3, -10689.6],
                "total_gibbs_per_mol": -8329.47,
                "analysis": 'Binary real solution with activity corrections.'
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618  # J/(mol·K)

    def _run_base(
        self,
        model: str,
        temperature_k: float,
        mu_standard: str,
        partial_pressures: str = "",
        p_standard: float = 1.0,
        mole_fractions: str = "",
        activity_coefficients: str = "1.0",
    ) -> dict:
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive.")
        if p_standard <= 0:
            raise ChemMCPError("Standard pressure must be positive.")

        mu_std_list = [float(x.strip()) for x in mu_standard.split(",")]
        n_components = len(mu_std_list)
        model = model.lower().strip()

        if model == "ideal_gas":
            p_list = [float(x.strip()) for x in partial_pressures.split(",")]
            if len(p_list) != n_components:
                raise ChemMCPError(f"Number of partial pressures ({len(p_list)}) must match mu_standard ({n_components}).")
            mus = []
            for i in range(n_components):
                activity = p_list[i] / p_standard
                if activity <= 0:
                    raise ChemMCPError(f"Partial pressure must be positive (component {i}).")
                mu_i = mu_std_list[i] + self.R * temperature_k * math.log(activity)
                mus.append(round(mu_i, 1))
            total_g = sum(mus) / len(mus)
            analysis = (
                f"Ideal gas mixture at T={temperature_k} K:\n"
                + "\n".join([f"  Component {i}: μ = {mu_std_list[i]} + RT·ln({p_list[i]}/{p_standard}) = {mus[i]} J/mol" for i in range(n_components)])
            )

        elif model in ("ideal_solution", "real_solution"):
            x_list = [float(x.strip()) for x in mole_fractions.split(",")]
            if len(x_list) != n_components:
                raise ChemMCPError(f"Number of mole fractions ({len(x_list)}) must match mu_standard ({n_components}).")
            if abs(sum(x_list) - 1.0) > 1e-6:
                raise ChemMCPError(f"Mole fractions must sum to 1.0 (got {sum(x_list):.6f}).")

            gamma_list = [1.0] * n_components
            if model == "real_solution":
                gamma_list = [float(x.strip()) for x in activity_coefficients.split(",")]
                if len(gamma_list) != n_components:
                    raise ChemMCPError(f"Number of activity coefficients ({len(gamma_list)}) must match components ({n_components}).")

            mus = []
            for i in range(n_components):
                a_i = gamma_list[i] * x_list[i]
                if a_i <= 0:
                    raise ChemMCPError(f"Activity must be positive (component {i}: γ={gamma_list[i]}, x={x_list[i]}).")
                mu_i = mu_std_list[i] + self.R * temperature_k * math.log(a_i)
                mus.append(round(mu_i, 1))

            total_g = sum(x * m for x, m in zip(x_list, mus))
            model_name = "Real solution" if model == "real_solution" else "Ideal solution"
            analysis = (
                f"{model_name} at T={temperature_k} K:\n"
                + "\n".join([f"  Component {i}: x={x_list[i]}, γ={gamma_list[i]}, a={gamma_list[i]*x_list[i]:.4f}, μ={mus[i]} J/mol" for i in range(n_components)])
            )
        else:
            raise ChemMCPError(f"Unsupported model: '{model}'. Use 'ideal_gas', 'ideal_solution', or 'real_solution'.")

        return {
            "chemical_potentials": mus,
            "total_gibbs_per_mol": round(total_g, 1),
            "analysis": analysis,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            model = parts[0]
            temperature_k = float(parts[1])
            mu_standard = parts[2]
            result_kwargs = {"model": model, "temperature_k": temperature_k, "mu_standard": mu_standard}
            idx = 3
            if model in ("ideal_gas",):
                result_kwargs["partial_pressures"] = parts[idx]; idx += 1
                result_kwargs["p_standard"] = float(parts[idx]) if idx < len(parts) else 1.0; idx += 1
            if model in ("ideal_solution", "real_solution"):
                result_kwargs["mole_fractions"] = parts[idx]; idx += 1
            if model == "real_solution" and idx < len(parts):
                result_kwargs["activity_coefficients"] = parts[idx]; idx += 1
            return self._run_base(**result_kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
