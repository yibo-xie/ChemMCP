"""
杠杆定则相组成计算工具
在两相区中使用杠杆定则计算各相的质量分数和组成分布。
"""
import logging
from typing import Dict, Any, Optional, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LeverRuleCalculator(BaseTool):
    """
    杠杆定则相组成计算工具。

    在二元相图的两相区中，使用杠杆定则计算各相的质量分数：
    W_α = (C_β - C_0) / (C_β - C_α)
    W_β = (C_0 - C_α) / (C_β - C_α)
    """
    __version__                 = "0.1.0"
    name                        = "LeverRuleCalculator"
    func_name                   = "apply_lever_rule"
    description                 = "Apply the lever rule to calculate phase mass fractions and compositions in a two-phase region of a binary phase diagram."
    implementation_description  = "Uses the lever rule formula: W₁=(C₂-C₀)/(C₂-C₁), W₂=(C₀-C₁)/(C₂-C₁). Calculates mass/mole fractions of each phase and the distribution of components between phases."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Lever Rule", "Phase Diagram", "Materials Science", "Thermodynamics", "Phase Fractions"]
    required_envs               = []

    code_input_sig = [
        ("overall_composition",  "float", "N/A",     "Overall composition (mass or mole fraction of component A/B), range 0 to 1."),
        ("composition_phase1",   "float", "N/A",     "Composition of phase 1 (boundary value on one side), same basis as overall."),
        ("composition_phase2",   "float", "N/A",     "Composition of phase 2 (boundary value on other side), same basis as overall."),
        ("phase1_name",          "str",   "α",        "Name of phase 1 (e.g., 'α', 'Liquid', 'solid solution')."),
        ("phase2_name",          "str",   "β",        "Name of phase 2 (e.g., 'β', 'Solid', 'vapor')."),
        ("total_mass_g",         "float", "100.0",    "Total mass of the system in grams (optional, for absolute mass calculation)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'overall_comp comp_phase1 comp_phase2 [phase1_name] [phase2_name] [total_mass_g]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing fraction_phase1, fraction_phase2, mass_distribution, component_distribution, and explanation."),
    ]

    examples = [
        {
            "code_input": {
                "overall_composition": 0.35,
                "composition_phase1": 0.15,
                "composition_phase2": 0.65,
                "phase1_name": "α",
                "phase2_name": "β",
                "total_mass_g": 100.0,
            },
            "text_input": {
                "input_params": "0.35 0.15 0.65 α β 100.0",
            },
            "output": {
                "result": {
                    "fraction_phase1": 0.5,
                    "fraction_phase2": 0.5,
                }
            },
        },
        # Hypoeutectic steel example (ferrite + cementite)
        {
            "code_input": {
                "overall_composition": 0.008,
                "composition_phase1": 0.000,
                "composition_phase2": 0.067,
                "phase1_name": "α-ferrite",
                "phase2_name": "Fe₃C",
                "total_mass_g": 500.0,
            },
            "text_input": {
                "input_params": "0.008 0.000 0.067 alpha-ferrite Fe3C 500.0",
            },
            "output": {
                "result": {
                    "fraction_phase1": 0.881,
                    "fraction_phase2": 0.119,
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
        overall_composition: float,
        composition_phase1: float,
        composition_phase2: float,
        phase1_name: str = "α",
        phase2_name: str = "β",
        total_mass_g: Optional[float] = None,
    ) -> Dict[str, Any]:
        """核心逻辑：应用杠杆定则。"""
        # ---- 输入验证 ----
        if not (0 <= overall_composition <= 1):
            raise ChemMCPError("Overall composition must be between 0 and 1.")
        if not (0 <= composition_phase1 <= 1):
            raise ChemMCPError("Phase 1 composition must be between 0 and 1.")
        if not (0 <= composition_phase2 <= 1):
            raise ChemMCPError("Phase 2 composition must be between 0 and 1.")

        # 检查是否在两相区内
        c_min = min(composition_phase1, composition_phase2)
        c_max = max(composition_phase1, composition_phase2)
        if not (c_min <= overall_composition <= c_max):
            logger.warning(
                f"Overall composition ({overall_composition}) is outside the two-phase "
                f"range [{c_min}, {c_max}]. Result may indicate single-phase region."
            )

        # 确认哪个是低组成边界、哪个是高组成边界
        if composition_phase1 < composition_phase2:
            C_low, C_high = composition_phase1, composition_phase2
            name_low, name_high = phase1_name, phase2_name
        else:
            C_low, C_high = composition_phase2, composition_phase1
            name_low, name_high = phase2_name, phase1_name

        # ---- 杠杆定则 ----
        denom = C_high - C_low
        if abs(denom) < 1e-15:
            raise ChemMCPError("Phase compositions are identical; cannot apply lever rule.")

        # 低组成相（富B）的分数
        w_low = (C_high - overall_composition) / denom
        # 高组成相（富A）的分数
        w_high = (overall_composition - C_low) / denom

        # 归一化（消除数值误差）
        total_w = w_low + w_high
        w_low /= total_w
        w_high /= total_w

        # 映射回原始相名
        if composition_phase1 < composition_phase2:
            w_phase1, w_phase2 = w_low, w_high
        else:
            w_phase1, w_phase2 = w_high, w_low

        result = {
            "overall_composition": overall_composition,
            "phase1_name": phase1_name,
            "phase2_name": phase2_name,
            "composition_phase1": composition_phase1,
            "composition_phase2": composition_phase2,
            "fraction_phase1": round(w_phase1, 6),
            "fraction_phase2": round(w_phase2, 6),
            "explanation": (
                f"Lever rule applied:\n"
                f"  • W_{phase1_name} = ({composition_phase2} - {overall_composition}) / "
                f"({composition_phase2} - {composition_phase1}) = {w_phase1:.4f}\n"
                f"  • W_{phase2_name} = ({overall_composition} - {composition_phase1}) / "
                f"({composition_phase2} - {composition_phase1}) = {w_phase2:.4f}\n"
                f"  • The system contains {w_phase1*100:.1f}% {phase1_name} and "
                f"{w_phase2*100:.1f}% {phase2_name}."
            ),
        }

        # 如果提供了总质量，计算各相的绝对质量
        if total_mass_g is not None and total_mass_g > 0:
            result["total_mass_g"] = total_mass_g
            result["mass_phase1_g"] = round(w_phase1 * total_mass_g, 4)
            result["mass_phase2_g"] = round(w_phase2 * total_mass_g, 4)

            # 计算组分在各相中的质量分布
            # 假设组成是质量分数
            mass_A_in_p1 = w_phase1 * total_mass_g * composition_phase1
            mass_B_in_p1 = w_phase1 * total_mass_g * (1 - composition_phase1)
            mass_A_in_p2 = w_phase2 * total_mass_g * composition_phase2
            mass_B_in_p2 = w_phase2 * total_mass_g * (1 - composition_phase2)

            result["component_distribution"] = {
                f"component_A_in_{phase1_name}_g": round(mass_A_in_p1, 4),
                f"component_B_in_{phase1_name}_g": round(mass_B_in_p1, 4),
                f"component_A_in_{phase2_name}_g": round(mass_A_in_p2, 4),
                f"component_B_in_{phase2_name}_g": round(mass_B_in_p2, 4),
            }

        logger.info(f"Lever rule: {phase1_name}={w_phase1:.4f}, {phase2_name}={w_phase2:.4f}, "
                     f"C0={overall_composition}")

        return result

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need: overall_comp comp_phase1 comp_phase2")

            kwargs = {
                "overall_composition": float(parts[0]),
                "composition_phase1": float(parts[1]),
                "composition_phase2": float(parts[2]),
            }
            if len(parts) > 3: kwargs["phase1_name"] = parts[3]
            if len(parts) > 4: kwargs["phase2_name"] = parts[4]
            if len(parts) > 5: kwargs["total_mass_g"] = float(parts[5])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'overall_comp comp1 comp2 [name1] [name2] [mass_g]'"
            )
