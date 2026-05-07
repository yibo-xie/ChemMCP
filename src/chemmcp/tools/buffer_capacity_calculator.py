"""
缓冲容量计算工具 (Buffer Capacity Calculator)

基于 van Slyke 方程计算缓冲溶液的缓冲容量（β值）。
"""

import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BufferCapacityCalculator(BaseTool):
    """
    缓冲容量计算工具。基于 van Slyke 方程计算缓冲溶液的缓冲能力。
    缓冲容量 β = dCb/dpH，表示使1L溶液改变1个pH单位所需加入的强酸或强碱的物质的量。
    """
    __version__ = "0.1.0"
    name = "BufferCapacityCalculator"
    func_name = "calculate_buffer_capacity"
    description = "Calculate buffer capacity (β) using the van Slyke equation for acid-base buffer systems."
    implementation_description = (
        "Uses van Slyke equation: β = 2.303 × Ca × (Ka×[H+]/(Ka+[H+])²) for a weak acid/conjugate base buffer, "
        "where Ca is total buffer concentration. Also supports approximate calculation via ΔCb/ΔpH."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Buffer", "Physical Chemistry", "Acid-Base", "Solution Chemistry", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("total_buffer_conc", "float", "N/A", "Total concentration of the buffer system (mol/L), i.e., [HA] + [A-]."),
        ("ph_initial", "float", "N/A", "Initial pH of the buffer solution."),
        ("pka", "float", "N/A", "pKa of the weak acid in the buffer system."),
        ("delta_ph", "float", "0.1", "pH change for approximate calculation (default 0.1). Used for ΔCb/ΔpH method."),
        ("calculation_mode", "str", "exact", "Calculation mode: 'exact' (van Slyke differential) or 'approximate' (ΔCb/ΔpH)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'total_buffer_conc ph_initial pka [delta_ph] [calculation_mode]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary with keys: beta(buffer_capacity, eq/L/pH), buffer_range(optimal_pH_range), "
         "max_beta(maximum_buffer_capacity), ha_ratio([HA]/[A-] ratio), detailed_calculation(str)"),
    ]

    examples = [
        {
            "code_input": {
                "total_buffer_conc": 0.10,
                "ph_initial": 4.76,
                "pka": 4.76,
                "delta_ph": 0.1,
                "calculation_mode": "exact",
            },
            "text_input": {
                "input_params": "0.10 4.76 4.76 0.1 exact"
            },
            "output": {
                "result": {
                    "beta": 0.0576,
                    "buffer_range": (3.76, 5.76),
                    "max_beta": 0.0576,
                }
            },
        },
        {
            "code_input": {
                "total_buffer_conc": 0.20,
                "ph_initial": 7.40,
                "pka": 7.21,  # H2PO4-/HPO4^2- system
                "delta_ph": 0.05,
                "calculation_mode": "exact",
            },
            "text_input": {
                "input_params": "0.20 7.40 7.21 0.05 exact"
            },
            "output": {
                "result": {
                    "beta": 0.4618,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        total_buffer_conc: float,
        ph_initial: float,
        pka: float,
        delta_ph: float = 0.1,
        calculation_mode: str = "exact",
    ) -> dict:
        """
        核心逻辑：缓冲容量计算

        Parameters:
            total_buffer_conc: 缓冲体系总浓度 (mol/L) = [HA] + [A-]
            ph_initial: 初始 pH 值
            pka: 弱酸的 pKa 值
            delta_ph: pH变化量（用于近似计算）
            calculation_mode: 计算模式 ('exact' 或 'approximate')

        Returns:
            dict: 缓冲容量及相关参数
        """
        # 输入验证
        if total_buffer_conc <= 0:
            raise ChemMCPError("Total buffer concentration must be positive.")
        if ph_initial < 0 or ph_initial > 14:
            raise ChemMCPError("pH must be between 0 and 14.")
        if delta_ph <= 0:
            raise ChemMCPError("Delta pH must be positive.")

        # 计算 H+ 浓度
        h_conc = 10.0 ** (-ph_initial)
        # 计算 Ka
        ka = 10.0 ** (-pka)

        # van Slyke 方程：β = 2.303 × C_total × (Ka × [H+]) / (Ka + [H+])²
        denominator = (ka + h_conc) ** 2
        if denominator == 0:
            raise ChemMCPError("Division by zero in buffer capacity calculation.")

        beta = 2.303 * total_buffer_conc * (ka * h_conc) / denominator

        # 最大缓冲容量（在 pH = pKa 时）
        max_beta = 2.303 * total_buffer_conc * 0.25  # 当 [H+] = Ka 时，Ka[H+]/(Ka+[H+])² = 1/4

        # 有效缓冲范围（通常为 pKa ± 1）
        buffer_range_low = pka - 1.0
        buffer_range_high = pka + 1.0

        # 计算 [HA]/[A-] 比率（Henderson-Hasselbalch方程）
        # pH = pKa + log([A-]/[HA]) => [A-]/[HA] = 10^(pH-pKa)
        if ph_initial != pka:
            ratio_a_minus_ha = 10.0 ** (ph_initial - pka)
            ha_fraction = 1.0 / (1.0 + ratio_a_minus_ha)
            a_minus_fraction = ratio_a_minus_ha / (1.0 + ratio_a_minus_ha)
        else:
            ha_fraction = 0.5
            a_minus_fraction = 0.5

        # 近似计算模式（ΔCb/ΔpH）
        approx_beta = None
        if calculation_mode == "approximate":
            # 模拟加入少量强碱后的pH变化
            ph_new = ph_initial + delta_ph
            h_new = 10.0 ** (-ph_new)
            # 新的 [A-]/[HA] 比率
            if h_new > 0:
                ratio_new = ka / h_new  # from Ka = [H+][A-]/[HA]
                # [A-]new = C_total × ratio_new/(1+ratio_new)
                a_minus_new = total_buffer_conc * ratio_new / (1.0 + ratio_new)
                a_minus_old = total_buffer_conc * a_minus_fraction
                delta_cb = a_minus_new - a_minus_old  # 加入的强碱量
                approx_beta = abs(delta_cb / delta_ph)

        # 构建详细计算过程
        calc_detail = (
            f"Buffer System Parameters:\n"
            f"  Total buffer concentration (Ca): {total_buffer_conc} mol/L\n"
            f"  Initial pH: {ph_initial}\n"
            f"  pKa: {pka}\n"
            f"  [H+]: {h_conc:.4e} mol/L\n"
            f"  Ka: {ka:.4e}\n\n"
            f"van Slyke Equation:\n"
            f"  β = 2.303 × Ca × (Ka×[H+]) / (Ka+[H+])²\n"
            f"  β = 2.303 × {total_buffer_conc} × ({ka:.4e}×{h_conc:.4e}) / ({ka:.4e}+{h_conc:.4e})²\n"
            f"  β = {beta:.6f} eq/(L·pH)\n\n"
            f"Maximum buffer capacity (at pH=pKa):\n"
            f"  β_max = 2.303 × Ca × 0.25 = {max_beta:.6f} eq/(L·pH)\n\n"
            f"Effective buffer range: pH {buffer_range_low:.2f} ~ {buffer_range_high:.2f}\n"
            f"  [HA] fraction: {ha_fraction:.4f} ({ha_fraction*100:.1f}%)\n"
            f"  [A-] fraction: {a_minus_fraction:.4f} ({a_minus_fraction*100:.1f}%)"
        )

        result = {
            "beta": round(beta, 6),
            "max_beta": round(max_beta, 6),
            "buffer_range": (round(buffer_range_low, 2), round(buffer_range_high, 2)),
            "ha_fraction": round(ha_fraction, 6),
            "a_minus_fraction": round(a_minus_fraction, 6),
            "detailed_calculation": calc_detail,
        }

        if approx_beta is not None:
            result["approximate_beta"] = round(approx_beta, 6)

        logger.info(f"Buffer capacity at pH={ph_initial}: β={beta:.6f} eq/(L·pH)")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 3:
                raise ValueError(f"Need at least 3 parameters, got {len(parts)}.")

            ca = float(parts[0])
            ph = float(parts[1])
            pka = float(parts[2])
            dph = float(parts[3]) if len(parts) > 3 else 0.1
            mode = parts[4] if len(parts) > 4 else "exact"

            return self._run_base(ca, ph, pka, dph, mode)
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
