import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CommonIonEffect(BaseTool):
    """
    分析同离子效应对弱酸/弱碱解离的影响。
    同离子效应会抑制弱电解质的解离，降低解离度。
    """
    __version__ = "0.1.0"
    name = "CommonIonEffect"
    func_name = "analyze_common_ion_effect"
    description = "Analyze the common ion effect on weak acid/base dissociation, comparing with and without the common ion present."
    implementation_description = "Calculates new equilibrium with common ion using Ka = [H+][A-]/[HA], where [A-]initial > 0 from the common ion salt."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Equilibrium", "Acid-Base", "Common Ion Effect", "Le Chatelier"]
    required_envs = []

    code_input_sig = [
        ("Ka", "float", "N/A", "Acid dissociation constant (Ka). For bases, provide Kb."),
        ("C0_weak_acid", "float", "N/A", "Initial concentration of the weak acid/base (mol/L)."),
        ("common_ion_conc", "float", "N/A", "Concentration of the common ion from added salt (mol/L)."),
        ("ion_type", "str", "anion", "Type of common ion: 'anion' (conjugate base) or 'cation' (for weak bases)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'Ka C0_weak_acid common_ion_conc ion_type'. Example: '1.8e-5 0.1 0.05 anion'"),
    ]

    output_sig = [
        ("alpha_original", "float", "Degree of dissociation WITHOUT common ion."),
        ("alpha_new", "float", "Degree of dissociation WITH common ion (should be lower)."),
        ("ph_original", "float", "pH without common ion."),
        ("ph_new", "float", "pH with common ion."),
        ("h_conc_original", "float", "[H+] without common ion (mol/L)."),
        ("h_conc_new", "float", "[H+] with common ion (mol/L)."),
        ("suppression_ratio", "float", "Ratio of new/original [H+] showing suppression factor."),
        ("explanation", "str", "Text explanation of the common ion effect."),
    ]

    examples = [
        {
            "code_input": {
                "Ka": 1.8e-5,
                "C0_weak_acid": 0.1,
                "common_ion_conc": 0.05,
                "ion_type": "anion"
            },
            "text_input": {
                "input_params": "1.8e-5 0.1 0.05 anion"
            },
            "output": {
                "alpha_original": 0.0134,
                "alpha_new": 0.000036,
                "ph_original": 2.87,
                "ph_new": 4.44,
                "h_conc_original": 0.00134,
                "h_conc_new": 3.6e-5,
                "suppression_ratio": 0.027,
                "explanation": "Common ion suppresses dissociation.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, Ka: float, C0_weak_acid: float, common_ion_conc: float, ion_type: str = "anion") -> dict:
        """
        核心逻辑：分析同离子效应
        
        无同离子时：HA ⇌ H+ + A-, Ka = C0·α²/(1-α)
        有同离子时（加入NaA提供[A-] = Cs）：
          HA ⇌ H+ + A-
          初始:  C0     0     Cs
          平衡:  C0-x   x     Cs+x
          Ka = x(Cs+x)/(C0-x)
          
          当 Cs >> x 时近似: x ≈ Ka·C0/Cs
        """
        if Ka <= 0:
            raise ChemMCPError("Ka must be positive.")
        if C0_weak_acid <= 0:
            raise ChemMCPError("Weak acid/base concentration must be positive.")
        if common_ion_conc < 0:
            raise ChemMCPError("Common ion concentration cannot be negative.")

        ion_type = ion_type.lower() if ion_type else "anion"

        # ---- 无同离子时的原始状态 ----
        disc_orig = Ka ** 2 + 4 * Ka * C0_weak_acid
        alpha_orig = (-Ka + math.sqrt(disc_orig)) / (2 * C0_weak_acid)
        h_orig = C0_weak_acid * alpha_orig
        ph_orig = -math.log10(h_orig) if h_orig > 0 else 7.0

        if common_ion_conc == 0:
            return {
                "alpha_original": round(alpha_orig, 6),
                "alpha_new": round(alpha_orig, 6),
                "ph_original": round(ph_orig, 4),
                "ph_new": round(ph_orig, 4),
                "h_conc_original": round(h_orig, 10),
                "h_conc_new": round(h_orig, 10),
                "suppression_ratio": 1.0,
                "explanation": "No common ion added. Dissociation is at normal level."
            }

        # ---- 有同离子时 ----
        # 精确求解: Ka = x(Cs+x)/(C0-x), 其中 x=[H+]
        # 展开: Ka(C0-x) = x(Cs+x)
        #       Ka·C0 - Ka·x = Cs·x + x²
        #       x² + (Cs+Ka)x - Ka·C0 = 0
        a_coef = 1.0
        b_coef = common_ion_conc + Ka
        c_coef = -Ka * C0_weak_acid

        discriminant = b_coef ** 2 - 4 * a_coef * c_coef
        if discriminant < 0:
            raise ChemMCPError("No real solution for equilibrium with given parameters.")

        h_new = (-b_coef + math.sqrt(discriminant)) / (2 * a_coef)
        if h_new < 0:
            h_new = 0.0

        alpha_new = h_new / C0_weak_acid if C0_weak_acid > 0 else 0.0
        ph_new = -math.log10(h_new) if h_new > 0 else 7.0

        suppression = h_new / h_orig if h_orig > 0 else 1.0

        explanation = (
            f"Adding {common_ion_conc} M of common ion suppresses dissociation. "
            f"[H+] drops from {h_orig:.2e} to {h_new:.2e} M "
            f"(suppressed by factor of {suppression:.4f}). "
            f"pH increases from {ph_orig:.2f} to {ph_new:.2f}."
        )

        logger.info(f"Common ion effect: α {alpha_orig:.6f}→{alpha_new:.6f}, pH {ph_orig:.2f}→{ph_new:.2f}")

        return {
            "alpha_original": round(alpha_orig, 6),
            "alpha_new": round(alpha_new, 6),
            "ph_original": round(ph_orig, 4),
            "ph_new": round(ph_new, 4),
            "h_conc_original": round(h_orig, 10),
            "h_conc_new": round(h_new, 10),
            "suppression_ratio": round(suppression, 6),
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 3:
                raise ValueError("Need Ka, C0, common_ion_conc. Format: 'Ka C0 common_ion [ion_type]'")
            Ka = float(parts[0])
            C0 = float(parts[1])
            c_common = float(parts[2])
            ion_type = parts[3] if len(parts) > 3 else "anion"
            return self._run_base(Ka, C0, c_common, ion_type)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'Ka C0 common_ion [ion_type]'")
