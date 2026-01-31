import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

# 设置日志
logger = logging.getLogger(__name__)

@ChemMCPManager.register_tool
class KcKpConverter(BaseTool):
    """
    化学平衡常数 Kc 与 Kp 相互转换工具。
    公式: Kp = Kc * (R * T)^delta_n
    R = 0.08206 L·atm/(K·mol)
    """
    __version__      = "0.1.0"
    name             = "KcKpConverter"
    func_name        = "convert_kc_kp"
    description      = "Convert between concentration equilibrium constant (Kc) and pressure equilibrium constant (Kp)."
    implementation_description = "Uses the thermodynamic relation Kp = Kc(RT)^delta_n where R = 0.08206 L·atm/(K·mol)."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Thermodynamics", "Physical Chemistry", "Equilibrium"]
    required_envs    = []  # 纯计算工具，不需要 API Key

    # 逻辑输入：数值、开氏温度、气体分子数变化量、转换方向
    code_input_sig   = [
        ("value", "float", "N/A", "The value of the equilibrium constant to be converted."),
        ("temperature_k", "float", "N/A", "Absolute temperature in Kelvin (K)."),
        ("delta_n", "float", "N/A", "Change in moles of gas (products - reactants)."),
        ("direction", "str", "kc_to_kp", "Conversion direction: 'kc_to_kp' or 'kp_to_kc'."),
    ]

    # 文本输入：空格分隔的字符串
    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated string: 'value temperature_k delta_n direction'."),
    ]

    # 输出
    output_sig       = [
        ("result", "float", "The converted equilibrium constant."),
    ]

    # 示例
    examples         = [
        {
            "code_input": {
                "value": 2.0,
                "temperature_k": 300.0,
                "delta_n": 1.0,
                "direction": "kc_to_kp"
            },
            "text_input": {
                "input_params": "2.0 300.0 1.0 kc_to_kp"
            },
            "output": {
                "result": 49.236
            }
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        # 此工具为纯数学计算，不需要外部客户端初始化
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """
        初始化常数
        """
        self.R = 0.08206  # L·atm/(K·mol)

    def _run_base(self, value: float, temperature_k: float, delta_n: float, direction: str = "kc_to_kp") -> float:
        """
        核心逻辑：实现 Kc 和 Kp 的转换
        """
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be a positive value in Kelvin.")

        # 计算 (R * T)^delta_n
        rt_factor = math.pow(self.R * temperature_k, delta_n)

        if direction.lower() == "kc_to_kp":
            # Kp = Kc * (RT)^dn
            result = value * rt_factor
        elif direction.lower() == "kp_to_kc":
            # Kc = Kp / (RT)^dn
            if rt_factor == 0:
                raise ChemMCPError("Invalid calculation: (RT)^delta_n resulted in zero.")
            result = value / rt_factor
        else:
            raise ChemMCPError("Invalid direction. Use 'kc_to_kp' or 'kp_to_kc'.")

        logger.info(f"Converted {value} with T={temperature_k}, dn={delta_n} ({direction}) -> {result}")
        return round(result, 6)

    def _run_text(self, input_params: str) -> float:
        """
        解析文本输入并调用核心逻辑
        """
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Insufficient parameters.")
            
            value = float(parts[0])
            temp = float(parts[1])
            dn = float(parts[2])
            direction = parts[3] if len(parts) > 3 else "kc_to_kp"
            
            return self._run_base(value, temp, dn, direction)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format should be 'value temp dn direction'")