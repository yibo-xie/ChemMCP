"""
亨利定律计算工具
计算气体在液体中的溶解度（C = k_H × P_gas）。
"""
import logging
from typing import Dict, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HenryLaw(BaseTool):
    """
    亨利定律计算工具。

    根据亨利定律 C = k_H × P_gas 计算气体溶解度，
    支持多种单位制转换。
    """
    __version__                 = "0.1.0"
    name                        = "HenryLaw"
    func_name                   = "calculate_henry_law"
    description                 = "Calculate gas solubility in liquids using Henry's law (C = k_H × P_gas), with support for multiple unit conventions."
    implementation_description  = "Applies Henry's law: C = k_H × P_gas. Supports common unit forms: M/atm, atm·L/mol, and g/(100mL·atm). Provides concentration in multiple units."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Henry's Law", "Gas Solubility", "Solution Chemistry", "Physical Chemistry"]
    required_envs               = []

    code_input_sig = [
        ("henry_constant",          "float", "N/A",     "Henry's law constant. Unit depends on 'units' parameter."),
        ("gas_partial_pressure_atm", "float", "N/A",     "Partial pressure of the gas in atmospheres (atm)."),
        ("units",                   "str",   "M/atm",   "Unit convention for Henry's constant: 'M/atm', 'atm·L/mol', or 'g/100mL/atm'."),
        ("solute_molar_mass_g_mol",  "float", "None",    "Molar mass of solute in g/mol (needed for g/100mL conversion; None to skip)."),
        ("temperature_k",           "float", "298.15",   "Temperature in Kelvin (for reference)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'henry_constant pressure_atm [units] [molar_mass_g_mol] [T(K)]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing dissolved_concentration, concentration_in_various_units, henry_constant_info, and explanation."),
    ]

    examples = [
        {
            "code_input": {
                "henry_constant": 0.0013,
                "gas_partial_pressure_atm": 2.0,
                "units": "M/atm",
                "solute_molar_mass_g_mol": 28.0,
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_params": "0.0013 2.0 M/atm 28.0 298.15",
            },
            "output": {
                "result": {
                    "dissolved_concentration_M": 0.0026,
                    "dissolved_concentration_g_100mL": "...",
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        henry_constant: float,
        gas_partial_pressure_atm: float,
        units: str = "M/atm",
        solute_molar_mass_g_mol: Optional[float] = None,
        temperature_k: float = 298.15,
    ) -> Dict[str, Any]:
        """核心逻辑：应用亨利定律计算气体溶解度。"""
        # ---- 输入验证 ----
        if gas_partial_pressure_atm < 0:
            raise ChemMCPError("Gas partial pressure cannot be negative.")
        if henry_constant < 0:
            raise ChemMCPError("Henry's constant cannot be negative.")
        valid_units = ["M/atm", "atm·L/mol", "g/100mL/atm", "atm*L/mol", "g/(100mL*atm)"]
        units_normalized = units.replace("·", "*").replace("（", "(").replace("）", ")")
        if units_normalized not in valid_units:
            # 宽容匹配
            for vu in valid_units:
                if units.lower().replace(" ", "").replace("·", "*") == vu.lower().replace(" ", ""):
                    units_normalized = vu
                    break

        # ---- 根据不同单位制计算浓度 ----
        if units_normalized in ("M/atm",):
            # C (M) = k_H (M/atm) * P (atm)
            conc_M = henry_constant * gas_partial_pressure_atm
        elif units_normalized in ("atm·L/mol", "atm*L/mol"):
            # k_H' (atm·L/mol): P = k_H' * x ≈ k_H' * C (for dilute)
            # C (mol/L) = P / k_H'
            if henry_constant > 0:
                conc_M = gas_partial_pressure_atm / henry_constant
            else:
                conc_M = float("inf")
        elif units_normalized in ("g/100mL/atm", "g/(100mL*atm)"):
            # C (g/100mL) = k_H * P
            conc_g_100ml = henry_constant * gas_partial_pressure_atm
            if solute_molar_mass_g_mol and solute_molar_mass_g_mol > 0:
                # 转换为 M: g/100mL → g/L → mol/L
                conc_g_L = conc_g_100ml * 10.0
                conc_M = conc_g_L / solute_molar_mass_g_mol
            else:
                conc_M = None
        else:
            raise ChemMCPError(f"Unsupported unit type: {units}. Use 'M/atm', 'atm·L/mol', or 'g/100mL/atm'.")

        # ---- 构建多单位结果 ----
        result = {
            "henry_constant_value": henry_constant,
            "henry_constant_units": units,
            "gas_partial_pressure_atm": gas_partial_pressure_atm,
            "temperature_K": temperature_k,
        }

        if conc_M is not None:
            result["dissolved_concentration_M"] = round(conc_M, 10)

            # 转换为其他单位
            if solute_molar_mass_g_mol and solute_molar_mass_g_mol > 0:
                conc_g_L = conc_M * solute_molar_mass_g_mol
                result["dissolved_concentration_g_L"] = round(conc_g_L, 6)
                result["dissolved_concentration_g_100mL"] = round(conc_g_L / 10.0, 6)
                result["molality_approx_mol_kg"] = round(conc_M, 6)  # 稀溶液近似

            # 计算 Bunsen 系数近似 (α = V_gas / (V_liquid * P))
            # 在 STP 下：α (mL gas/mL liquid/atm) ≈ C(M) * 22400 mL/mol / 1000 mL/L
            bunsen_coeff = conc_M * 22.4  # 近似值
            result["bunsen_coefficient_approx"] = round(bunsen_coeff, 6)

        if units_normalized in ("g/100mL/atm", "g/(100mL*atm)") and 'conc_g_100ml' in dir():
            result["dissolved_concentration_g_100mL_direct"] = round(conc_g_100ml, 6)

        result["explanation"] = (
            f"Henry's law: C = k_H × P = {henry_constant} × {gas_partial_pressure_atm} "
            f"= {conc_M if conc_M is not None else conc_g_100ml} ({units}). "
            f"At T={temperature_k}K."
        )

        logger.info(f"Henry's law: k_H={henry_constant} ({units}), P={gas_partial_pressure_atm} atm → "
                     f"C={conc_M} M")
        return result

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least: henry_constant pressure_atm")

            kwargs = {
                "henry_constant": float(parts[0]),
                "gas_partial_pressure_atm": float(parts[1]),
            }
            if len(parts) > 2: kwargs["units"] = parts[2]
            if len(parts) > 3: kwargs["solute_molar_mass_g_mol"] = float(parts[3])
            if len(parts) > 4: kwargs["temperature_k"] = float(parts[4])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'k_H P [units] [molar_mass] [T(K)]'"
            )
