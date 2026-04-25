import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BufferCapacity(BaseTool):
    """
    计算缓冲容量（Buffer Capacity, β）。
    使用 van Slyke 公式计算缓冲溶液抵抗 pH 变化的能力。
    """
    __version__ = "0.1.0"
    name = "BufferCapacity"
    func_name = "calculate_buffer_capacity"
    description = "Calculate buffer capacity (β) using the van Slyke equation: β = 2.303 × ([H+] + [OH-] + C_total·Ka·[H+]/(Ka+[H+])²)."
    implementation_description = "Implements van Slyke buffer capacity formula to quantify buffer's resistance to pH change."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Buffer Capacity", "van Slyke", "Acid-Base", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("ha_conc", "float", "N/A", "Concentration of weak acid [HA] in mol/L."),
        ("a_conc", "float", "N/A", "Concentration of conjugate base [A-] in mol/L."),
        ("Ka", "float", "N/A", "Acid dissociation constant Ka."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'ha_conc a_conc Ka'. Example: '0.1 0.1 1.8e-5'"),
    ]

    output_sig = [
        ("beta", "float", "Buffer capacity β (mol·L⁻¹·pH⁻¹)."),
        ("ph", "float", "Current pH of the buffer."),
        ("ph_effective_low", "float", "Lower bound of effective buffer range (pKa - 1)."),
        ("ph_effective_high", "float", "Upper bound of effective buffer range (pKa + 1)."),
        ("total_buffer_conc", "float", "Total concentration C_total = [HA] + [A-]."),
        ("explanation", "str", "Interpretation of the buffer capacity value."),
    ]

    examples = [
        {
            "code_input": {
                "ha_conc": 0.1,
                "a_conc": 0.1,
                "Ka": 1.8e-5,
            },
            "text_input": {
                "input_params": "0.1 0.1 1.8e-5"
            },
            "output": {
                "beta": 0.115,
                "ph": 4.76,
                "ph_effective_low": 3.76,
                "ph_effective_high": 5.76,
                "total_buffer_conc": 0.2,
                "explanation": "Good buffer capacity.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Kw = 1.0e-14

    def _run_base(self, ha_conc: float, a_conc: float, Ka: float) -> dict:
        """
        核心逻辑：van Slyke 缓冲容量公式
        
        β = dCb/dpH = -dCa/dpH = 2.303 × ([H+] + [OH-] + C_total·Ka·[H+]/(Ka+[H+])²)
        
        其中:
          [H+] 由 Henderson-Hasselbalch 计算
          C_total = [HA] + [A-]
        """
        if ha_conc < 0 or a_conc < 0:
            raise ChemMCPError("Concentrations cannot be negative.")
        if Ka <= 0:
            raise ChemMCPError("Ka must be positive.")

        # 计算当前 pH 和 [H+]
        pKa = -math.log10(Ka)
        if a_conc > 0 and ha_conc > 0:
            ph = pKa + math.log10(a_conc / ha_conc)
        elif a_conc == 0:
            # 纯弱酸
            disc = Ka ** 2 + 4 * Ka * ha_conc
            h = (-Ka + math.sqrt(disc)) / 2
            ph = -math.log10(h) if h > 0 else 7.0
        else:
            ph = 14.0  # 纯共轭碱

        h_conc = 10 ** (-ph)
        oh_conc = self.Kw / h_conc

        C_total = ha_conc + a_conc

        # van Slyke 公式
        term_acid = h_conc           # 强酸贡献
        term_base = oh_conc          # 强碱贡献
        denom = (Ka + h_conc) ** 2
        if denom > 0:
            term_buffer = C_total * Ka * h_conc / denom
        else:
            term_buffer = 0.0

        beta = 2.303 * (term_acid + term_base + term_buffer)

        # 有效缓冲范围：pKa ± 1
        effective_low = max(0, pKa - 1)
        effective_high = min(14, pKa + 1)

        # 解释
        if beta >= 0.1:
            quality = "good buffer capacity"
        elif beta >= 0.01:
            quality = "moderate buffer capacity"
        else:
            quality = "low buffer capacity"

        explanation = (
            f"β = {beta:.4f} mol·L⁻¹·pH⁻¹ → {quality}. "
            f"Effective range: pH {effective_low:.2f}–{effective_high:.2f}. "
            f"Higher total concentration and ratio near 1:1 give maximum capacity."
        )

        logger.info(f"Buffer capacity: β={beta:.6f} at pH={ph:.4f}")

        return {
            "beta": round(beta, 6),
            "ph": round(ph, 4),
            "ph_effective_low": round(effective_low, 4),
            "ph_effective_high": round(effective_high, 4),
            "total_buffer_conc": round(C_total, 8),
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) < 3:
                raise ValueError("Need ha_conc, a_conc, Ka.")
            ha = float(parts[0])
            a_ = float(parts[1])
            ka = float(parts[2])
            return self._run_base(ha, a_, ka)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'ha_conc a_conc Ka'")
