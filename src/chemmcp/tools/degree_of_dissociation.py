import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DegreeOfDissociation(BaseTool):
    """
    计算弱电解质的解离度（α）。
    支持精确求解和近似公式。
    """
    __version__ = "0.1.0"
    name = "DegreeOfDissociation"
    func_name = "calculate_degree_of_dissociation"
    description = "Calculate the degree of dissociation (α) for weak electrolytes (weak acids or bases)."
    implementation_description = "Solves the equilibrium expression Ka = C0·α²/(1-α) exactly using quadratic formula, with sqrt(Ka/C0) approximation available."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Equilibrium", "Acid-Base", "Physical Chemistry", "Dissociation"]
    required_envs = []

    code_input_sig = [
        ("C0", "float", "N/A", "Initial concentration of the weak electrolyte (mol/L or M)."),
        ("K_eq", "float", "N/A", "Equilibrium constant (Ka for acids, Kb for bases)."),
        ("method", "str", "exact", "Calculation method: 'exact' (quadratic) or 'approximate' (sqrt(K/C))."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'C0 K_eq [method]'. Example: '0.1 1.8e-5 exact'"),
    ]

    output_sig = [
        ("alpha", "float", "Degree of dissociation (0 to 1)."),
        ("H_conc", "float", "Concentration of dissociated H+ (for acids) or OH- (for bases) in mol/L."),
        ("HA_eq_conc", "float", "Equilibrium concentration of undissociated species in mol/L."),
        ("ph", "float", "pH of the solution (calculated from [H+])."),
        ("method_used", "str", "The method actually used for calculation."),
    ]

    examples = [
        {
            "code_input": {
                "C0": 0.1,
                "K_eq": 1.8e-5,
                "method": "exact"
            },
            "text_input": {
                "input_params": "0.1 1.8e-5 exact"
            },
            "output": {
                "alpha": 0.0134,
                "H_conc": 0.00134,
                "HA_eq_conc": 0.09866,
                "ph": 2.873,
                "method_used": "exact"
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, C0: float, K_eq: float, method: str = "exact") -> dict:
        """
        核心逻辑：计算解离度 α
        对于弱酸 HA ⇌ H+ + A-:
          Ka = C0·α² / (1-α)
        精确解: α = [-Ka + sqrt(Ka² + 4·Ka·C0)] / (2·C0)
        近似解(当α很小时): α ≈ sqrt(Ka/C0)
        """
        if C0 <= 0:
            raise ChemMCPError("Initial concentration C0 must be positive.")
        if K_eq <= 0:
            raise ChemMCPError("Equilibrium constant K_eq must be positive.")

        method = method.lower() if method else "exact"

        if method == "approximate":
            # α ≈ sqrt(Ka/C0)
            alpha = math.sqrt(K_eq / C0)
            # 检查近似是否合理 (α < 0.05 即 5%)
            if alpha >= 0.05:
                logger.warning(f"Approximation may be inaccurate: α={alpha:.4f} ≥ 0.05. Use 'exact' method.")
        else:
            # 精确求解二次方程: C0·α² + Ka·α - Ka = 0
            # α = [-Ka + sqrt(Ka² + 4·Ka·C0)] / (2·C0)
            discriminant = K_eq ** 2 + 4 * K_eq * C0
            alpha = (-K_eq + math.sqrt(discriminant)) / (2 * C0)

        if alpha < 0 or alpha > 1:
            raise ChemMCPError(f"Calculated α={alpha:.6f} is out of valid range [0, 1]. Check input values.")

        H_conc = C0 * alpha
        HA_eq = C0 * (1 - alpha)

        # pH = -log10([H+])
        if H_conc > 0:
            ph = -math.log10(H_conc)
        else:
            ph = 7.0

        logger.info(f"Degree of dissociation: α={alpha:.6f}, pH={ph:.4f}, method={method}")

        return {
            "alpha": round(alpha, 6),
            "H_conc": round(H_conc, 8),
            "HA_eq_conc": round(HA_eq, 8),
            "ph": round(ph, 4),
            "method_used": method,
        }

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need at least C0 and K_eq. Format: 'C0 K_eq [method]'")
            C0 = float(parts[0])
            K_eq = float(parts[1])
            method = parts[2] if len(parts) > 2 else "exact"
            return self._run_base(C0, K_eq, method)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'C0 K_eq [method]'")
