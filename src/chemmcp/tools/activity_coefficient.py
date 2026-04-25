import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ActivityCoefficient(BaseTool):
    """
    计算溶液中组分的活度系数。
    
    支持方法：
    - Debye-Hückel极限公式: log₁₀(γ±) = -A|z₊z₋|√I
    - 扩展Debye-Hückel: log₁₀(γ±) = -A|z₊z₋|√I / (1 + Ba√I)
    """
    __version__ = "0.1.0"
    name = "ActivityCoefficient"
    func_name = "calculate_activity_coefficient"
    description = "Calculate mean ionic activity coefficient in electrolyte solutions using Debye-Hückel theory."
    implementation_description = "Supports 'debye_huckel_limiting' and 'extended_debye_huckel' methods. Uses A=0.509 (aqueous, 25°C), B=0.328 (nm⁻¹) for extended form. Calculates γ± for given ionic strength and ion charges."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Activity Coefficient", "Solution Chemistry", "Electrolyte", "Debye-Hückel"]
    required_envs = []

    code_input_sig = [
        ("method", "str", "extended_debye_huckel", "Method: 'debye_huckel_limiting' or 'extended_debye_huckel'."),
        ("ionic_strength", "float", "N/A", "Ionic strength I in mol/L."),
        ("z_plus", "int", "N/A", "Charge number of cation (e.g., +1 for Na⁺, +2 for Ca²⁺)."),
        ("z_minus", "int", "N/A", "Charge number of anion (e.g., -1 for Cl⁻, -2 for SO₄²⁻)."),
        ("a_ion_size", "float", "0.3", "Ion size parameter a in nm (for extended Debye-Hückel). Typical values: 0.3-0.9 nm."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'method ionic_strength z_plus z_minus [a_ion_size]'. Example: 'extended_debye_huckel 0.01 1 -1 0.4'"),
    ]

    output_sig = [
        ("gamma_pm", "float", "Mean ionic activity coefficient γ± (dimensionless)."),
        ("log_gamma", "float", "Base-10 logarithm of γ±."),
        ("explanation", "str", "Calculation details with formula used."),
    ]

    examples = [
        {
            "code_input": {
                "method": "debye_huckel_limiting",
                "ionic_strength": 0.005,
                "z_plus": 1,
                "z_minus": -1,
                "a_ion_size": 0.3,
            },
            "text_input": {
                "input_params": "debye_huckel_limiting 0.005 1 -1",
            },
            "output": {
                "gamma_pm": 0.9207,
                "log_gamma": -0.03596,
                "explanation": "Debye-Hückel limiting: log10(γ±) = -0.509×|1×(-1)|×√0.005 = -0.0360 → γ± = 0.921",
            },
        },
        {
            "code_input": {
                "method": "extended_debye_huckel",
                "ionic_strength": 0.05,
                "z_plus": 2,
                "z_minus": -1,
                "a_ion_size": 0.4,
            },
            "text_input": {
                "input_params": "extended_debye_huckel 0.05 2 -1 0.4",
            },
            "output": {
                "gamma_pm": 0.4586,
                "log_gamma": -0.3385,
                "explanation": "Extended Debye-Hückel: log10(γ±) = -0.509×|2×(-1)|×√0.05/(1+0.328×0.4×√0.05) = -0.338 → γ± = 0.459",
            },
        },
    ]

    # Constants for aqueous solution at 25°C
    A = 0.509   # (mol/L)^(-1/2)
    B = 0.328   # nm⁻¹

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, method: str, ionic_strength: float, z_plus: int, z_minus: int, a_ion_size: float = 0.3) -> dict:
        """Core logic: calculate mean ionic activity coefficient."""
        if ionic_strength < 0:
            raise ChemMCPError("Ionic strength must be non-negative.")
        if z_plus == 0 or z_minus == 0:
            raise ChemMCPError("Ion charges cannot be zero.")
        if a_ion_size <= 0:
            raise ChemMCPError("Ion size parameter must be positive.")

        method = method.lower().strip()
        sqrt_I = math.sqrt(ionic_strength)
        abs_z_prod = abs(z_plus * z_minus)

        if method == "debye_huckel_limiting":
            log_gamma = -self.A * abs_z_prod * sqrt_I
            formula_str = f"log₁₀(γ±) = -A·|z₊·z₋|·√I = -{self.A} × {abs_z_prod} × √{ionic_strength}"
        elif method == "extended_debye_huckel":
            denominator = 1.0 + self.B * a_ion_size * sqrt_I
            log_gamma = -self.A * abs_z_prod * sqrt_I / denominator
            formula_str = (
                f"log₁₀(γ±) = -A·|z₊·z₋|·√I / (1 + B·a·√I)\n"
                f"= -{self.A} × {abs_z_prod} × √{ionic_strength} / (1 + {self.B} × {a_ion_size} × √{ionic_strength})"
            )
        else:
            raise ChemMCPError(f"Unknown method '{method}'. Use 'debye_huckel_limiting' or 'extended_debye_huckel'.")

        gamma_pm = 10.0 ** log_gamma

        explanation = (
            f"{method.replace('_', ' ').title()}:\n"
            f"{formula_str}\n"
            f"= {log_gamma:.6f}\n"
            f"γ± = 10^{log_gamma:.6f} = {gamma_pm:.6f}"
        )

        logger.info(f"ActivityCoeff: method={method}, I={ionic_strength}, z+={z_plus}, z-={z_minus}, γ±={gamma_pm:.6f}")
        return {
            "gamma_pm": round(gamma_pm, 6),
            "log_gamma": round(log_gamma, 6),
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            method = parts[0]
            I = float(parts[1])
            zp = int(parts[2])
            zm = int(parts[3])
            a_size = float(parts[4]) if len(parts) > 4 else 0.3
            return self._run_base(method, I, zp, zm, a_size)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'method ionic_strength z_plus z_minus [a_ion_size]'")
