import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BufferPreparation(BaseTool):
    """
    缓冲溶液配制计算（Henderson-Hasselbalch 方程）。
    根据目标 pH 计算所需的酸和共轭碱的量。
    """
    __version__ = "0.1.0"
    name = "BufferPreparation"
    func_name = "prepare_buffer"
    description = "Calculate buffer solution preparation using the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])."
    implementation_description = "Uses Henderson-Hasselbalch equation to compute required amounts of acid and conjugate base for a target pH."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Buffer", "Henderson-Hasselbalch", "Acid-Base", "Solution Preparation"]
    required_envs = []

    code_input_sig = [
        ("target_ph", "float", "N/A", "Target pH value for the buffer."),
        ("pKa", "float", "N/A", "pKa of the weak acid component."),
        ("total_volume_l", "float", "0.001", "Total desired volume of buffer in liters."),
        ("total_concentration_m", "float", "0.1", "Total concentration of buffer components (mol/L)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'target_ph pKa [total_volume_L] [total_conc_M]'. Example: '4.75 4.76 0.5 0.1'"),
    ]

    output_sig = [
        ("target_ph", "float", "Requested target pH."),
        ("pKa", "float", "pKa used in calculation."),
        ("ratio_base_acid", "float", "The ratio [A-]/[HA] needed to achieve target pH."),
        ("ha_concentration_m", "float", "Required concentration of acid (HA) in mol/L."),
        ("a_concentration_m", "float", "Required concentration of conjugate base (A-) in mol/L."),
        ("achieved_ph", "float", "Actual pH achieved with calculated ratio (should match target)."),
        ("volume_l", "float", "Total buffer volume in L."),
        ("explanation", "str", "Step-by-step preparation instructions."),
    ]

    examples = [
        {
            "code_input": {
                "target_ph": 4.75,
                "pKa": 4.76,
                "total_volume_l": 0.5,
                "total_concentration_m": 0.1,
            },
            "text_input": {
                "input_params": "4.75 4.76 0.5 0.1"
            },
            "output": {
                "target_ph": 4.75,
                "pKa": 4.76,
                "ratio_base_acid": 0.977,
                "ha_concentration_m": 0.0506,
                "a_concentration_m": 0.0494,
                "achieved_ph": 4.75,
                "volume_l": 0.5,
                "explanation": "Buffer preparation instructions.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, target_ph: float, pKa: float, total_volume_l: float = 0.001, total_concentration_m: float = 0.1) -> dict:
        """
        核心逻辑：Henderson-Hasselbalch 方程
        pH = pKa + log([A-]/[HA])
        => [A-]/[HA] = 10^(pH - pKa)
        
        设总浓度 Ct = [HA] + [A-]
        则 [HA] = Ct / (1 + r), [A-] = r·Ct / (1 + r), 其中 r = [A-]/[HA]
        """
        if target_ph < 0 or target_ph > 14:
            raise ChemMCPError("pH must be between 0 and 14.")
        if pKa <= 0:
            raise ChemMCPError("pKa must be positive.")
        if total_volume_l <= 0:
            raise ChemMCPError("Volume must be positive.")
        if total_concentration_m <= 0:
            raise ChemMCPError("Concentration must be positive.")

        # 计算比例 r = [A-]/[HA]
        ratio = 10 ** (target_ph - pKa)

        # 各组分浓度
        C_total = total_concentration_m
        ha_conc = C_total / (1 + ratio)
        a_conc = ratio * C_total / (1 + ratio)

        # 验证：反算 pH
        achieved_ph = pKa + math.log10(ratio) if ratio > 0 else pKa

        # 配制说明
        mol_ha = ha_conc * total_volume_l
        mol_a = a_conc * total_volume_l

        explanation = (
            f"To prepare {total_volume_l*1000:.1f} mL of {C_total} M buffer at pH {target_ph:.2f}:\n"
            f"1. Ratio [A-]/[HA] = 10^(pH - pKa) = 10^({target_ph:.2f} - {pKa:.2f}) = {ratio:.4f}\n"
            f"2. [HA] = {ha_conc:.6f} M ({mol_ha*1000:.3f} mmol)\n"
            f"3. [A-] = {a_conc:.6f} M ({mol_a*1000:.3f} mmol)\n"
            f"4. Mix components and dilute to {total_volume_l*1000:.1f} mL with water.\n"
            f"5. Effective buffer range: pKa ± 1 → {pKa-1:.2f} to {pKa+1:.2f}"
        )

        logger.info(f"Buffer prep: pH={target_ph}, pKa={pKa}, ratio={ratio:.4f}")

        return {
            "target_ph": round(target_ph, 4),
            "pKa": round(pKa, 4),
            "ratio_base_acid": round(ratio, 6),
            "ha_concentration_m": round(ha_conc, 8),
            "a_concentration_m": round(a_conc, 8),
            "achieved_ph": round(achieved_ph, 4),
            "volume_l": round(total_volume_l, 6),
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need at least target_ph and pKa.")
            target_ph = float(parts[0])
            pKa = float(parts[1])
            vol = float(parts[2]) if len(parts) > 2 else 0.001
            conc = float(parts[3]) if len(parts) > 3 else 0.1
            return self._run_base(target_ph, pKa, vol, conc)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'target_ph pKa [vol_L] [conc_M]'")
