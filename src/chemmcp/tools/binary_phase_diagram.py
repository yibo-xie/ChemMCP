"""
二元相图分析与杠杆定则工具
支持匀晶（isomorphous）、共晶（eutectic）、包晶（peritectic）类型相图的分析，
以及杠杆定则计算相组成。
"""
import logging
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BinaryPhaseDiagram(BaseTool):
    """
    二元相图分析与杠杆定则工具。

    支持分析二元合金/混合物相图，判断当前组成和温度所处的相区，
    并使用杠杆定则计算各相的质量分数。
    """
    __version__                 = "0.1.0"
    name                        = "BinaryPhaseDiagram"
    func_name                   = "analyze_binary_phase_diagram"
    description                 = "Analyze binary phase diagrams (isomorphous, eutectic, peritectic) and apply the lever rule to calculate phase compositions and fractions."
    implementation_description  = "Accepts phase diagram parameters (liquidus/solidus/eutectic points) and overall composition. Determines phase regions and applies lever rule: W_α=(C_β-C0)/(C_β-C_α), W_β=(C0-C_α)/(C_β-C_α)."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Phase Diagram", "Lever Rule", "Materials Science", "Thermodynamics", "Alloy"]
    required_envs               = []

    code_input_sig = [
        ("composition_a",          "float", "N/A",     "Overall mass fraction of component A (0 to 1)."),
        ("phase_diagram_type",      "str",   "N/A",     "Type of phase diagram: 'isomorphous', 'eutectic', or 'peritectic'."),
        ("temperature_k",           "float", "N/A",     "Current temperature in Kelvin."),
        ("melting_point_a_k",       "float", "N/A",     "Melting point of pure A in Kelvin."),
        ("melting_point_b_k",       "float", "N/A",     "Melting point of pure B in Kelvin."),
        ("t_eutectic_or_peritectic","float", "None",     "Eutectic/peritectic temperature in K (None for isomorphous)."),
        ("composition_eutectic",   "float", "None",     "Composition at eutectic point, mass fraction A (None if not applicable)."),
        ("max_solubility_b_in_a",  "float", "None",     "Maximum solubility of B in α phase at eutectic temp, mass fraction B (None for isomorphous)."),
        ("max_solubility_a_in_b",  "float", "None",     "Maximum solubility of A in β phase at eutectic temp, mass fraction A (None for isomorphous)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'composition_A type T(K) Tm_A(K) Tm_B(K) [T_eut(K) x_eut max_sol_BinA max_sol_AinB]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing current_phase, phases_present, phase_compositions, lever_rule_fractions, phase_diagram_summary, and explanation."),
    ]

    examples = [
        {
            "code_input": {
                "composition_a": 0.35,
                "phase_diagram_type": "eutectic",
                "temperature_k": 800.0,
                "melting_point_a_k": 1200.0,
                "melting_point_b_k": 900.0,
                "t_eutectic_or_peritectic": 577.0,
                "composition_eutectic": 0.627,
                "max_solubility_b_in_a": 0.052,
                "max_solubility_a_in_b": 0.088,
            },
            "text_input": {
                "input_params": "0.35 eutectic 800.0 1200.0 900.0 577.0 0.627 0.052 0.088",
            },
            "output": {
                "result": {
                    "current_phase": "...",
                    "phases_present": ["..."],
                    "lever_rule_fractions": {...},
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
        composition_a: float,
        phase_diagram_type: str,
        temperature_k: float,
        melting_point_a_k: float,
        melting_point_b_k: float,
        t_eutectic_or_peritectic: Optional[float] = None,
        composition_eutectic: Optional[float] = None,
        max_solubility_b_in_a: Optional[float] = None,
        max_solubility_a_in_b: Optional[float] = None,
    ) -> Dict[str, Any]:
        """核心逻辑：分析二元相图并应用杠杆定则。"""
        # ---- 输入验证 ----
        if not (0 <= composition_a <= 1):
            raise ChemMCPError("Composition of A must be between 0 and 1.")
        if phase_diagram_type not in ("isomorphous", "eutectic", "peritectic"):
            raise ChemMCPError("Phase diagram type must be 'isomorphous', 'eutectic', or 'peritectic'.")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")
        if melting_point_a_k <= 0 or melting_point_b_k <= 0:
            raise ChemMCPError("Melting points must be positive.")

        # 确保纯 A 的熔点 > 纯 B 的熔点（A 为高熔点组分）
        result = {
            "overall_composition_A": composition_a,
            "temperature_K": temperature_k,
            "phase_diagram_type": phase_diagram_type,
        }

        if phase_diagram_type == "isomorphous":
            result.update(self._analyze_isomorphous(
                composition_a, temperature_k, melting_point_a_k, melting_point_b_k
            ))
        elif phase_diagram_type == "eutectic":
            result.update(self._analyze_eutectic(
                composition_a, temperature_k, melting_point_a_k, melting_point_b_k,
                t_eutectic_or_peritectic, composition_eutectic,
                max_solubility_b_in_a, max_solubility_a_in_b
            ))
        elif phase_diagram_type == "peritectic":
            result.update(self._analyze_peritectic(
                composition_a, temperature_k, melting_point_a_k, melting_point_b_k,
                t_eutectic_or_peritectic, composition_eutectic,
                max_solubility_b_in_a, max_solubility_a_in_b
            ))

        logger.info(f"Phase diagram analysis ({phase_diagram_type}): x_A={composition_a}, "
                     f"T={temperature_k}K → {result.get('current_phase', 'unknown')}")
        return result

    def _linear_interpolate(self, x0, y0, x1, y1, x):
        """线性插值。"""
        if x1 == x0:
            return y0
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

    def _analyze_isomorphous(self, comp_a, T, Tm_A, Tm_B):
        """分析匀晶相图。"""
        # 假设线性液相线和固相线
        # 液相线：从 (0, Tm_B) 到 (1, Tm_A)
        # 固相线：在液相线下方，简化为液相线偏移一定量
        T_liquidus = self._linear_interpolate(0, Tm_B, 1, Tm_A, comp_a)
        T_solidus = T_liquidus - 50  # 简化：固相线比液相线低 50K

        if T >= Tm_A:
            phase = "pure liquid" if comp_a == 1.0 else "liquid"
            phases = ["Liquid"]
            fractions = {"Liquid": 1.0}
            compositions = {"Liquid": comp_a}
        elif T >= T_liquidus:
            phase = "liquid (single-phase)"
            phases = ["Liquid"]
            fractions = {"Liquid": 1.0}
            compositions = {"Liquid": comp_a}
        elif T >= T_solidus:
            # 两相区：L + α
            # 近似：液相线和固相线成分
            c_liquid = (T - Tm_B) / (Tm_A - Tm_B)  # 反向插值
            c_liquid = max(0, min(1, c_liquid))
            c_solid = c_liquid + 0.05  # 固相富A
            c_solid = min(1, c_solid)

            if abs(c_liquid - c_solid) < 1e-10:
                w_liq, w_sol = 0.5, 0.5
            else:
                w_sol = (c_liquid - comp_a) / (c_liquid - c_solid)
                w_liq = (comp_a - c_solid) / (c_liquid - c_solid)

            phase = "two-phase (L + α solid solution)"
            phases = ["Liquid", "α"]
            fractions = {"Liquid": round(max(0, min(1, w_liq)), 4), "α": round(max(0, min(1, w_sol)), 4)}
            compositions = {"Liquid": round(c_liquid, 4), "α": round(c_solid, 4)}
        elif T <= Tm_B:
            phase = "solid (single-phase)" if comp_a < 1.0 else "pure solid"
            phases = ["α"]
            fractions = {"α": 1.0}
            compositions = {"α": comp_a}
        else:
            phase = "solid solution (α)"
            phases = ["α"]
            fractions = {"α": 1.0}
            compositions = {"α": comp_a}

        return {
            "current_phase": phase,
            "phases_present": phases,
            "phase_compositions": compositions,
            "lever_rule_fractions": fractions,
            "explanation": f"At T={T}K with x_A={comp_a}: {phase}. "
                           f"T_liquidus≈{T_liquidus:.0f}K, T_solidus≈{T_solidus:.0f}K."
        }

    def _analyze_eutectic(self, comp_a, T, Tm_A, Tm_B, T_eut, x_eut, max_sol_BinA, max_sol_AinB):
        """分析共晶相图。"""
        T_eut = T_eut or 0
        x_e = x_eut or 0.5

        # 判断区域
        if T >= max(Tm_A, Tm_B):
            phase = "all liquid"
            phases = ["L"]
            fractions = {"L": 1.0}
            comps = {"L": comp_a}
        elif T > T_eut:
            # 在共晶温度以上
            if comp_a <= x_e:
                # 左侧：L + α 区或纯 L 或纯 α
                T_liquidus_left = self._linear_interpolate(0, Tm_A, x_e, T_eut, comp_a)
                if T >= T_liquidus_left:
                    phase = "liquid"
                    phases = ["L"]
                    fractions = {"L": 1.0}
                    comps = {"L": comp_a}
                else:
                    # L + α 两相区
                    c_liq = (T - T_eut) / (T_liquidus_left - T_eut) * (comp_a - 0) + 0
                    c_liq = min(c_liq, comp_a + 0.1)
                    c_alpha = max_sol_BinA if max_sol_BinA is not None else 0.02
                    c_alpha_comp_a = 1.0 - c_alpha

                    denom = c_liq - c_alpha_comp_a
                    if abs(denom) < 1e-10:
                        w_a, w_l = 0.5, 0.5
                    else:
                        w_a = (c_liq - comp_a) / denom
                        w_l = (comp_a - c_alpha_comp_a) / denom

                    phase = "two-phase (L + α)"
                    phases = ["L", "α"]
                    fractions = {"L": round(max(0, min(1, w_l)), 4), "α": round(max(0, min(1, w_a)), 4)}
                    comps = {"L": round(c_liq, 4), "α": round(c_alpha_comp_a, 4)}
            else:
                # 右侧：L + β 区或纯 L
                T_liquidus_right = self._linear_interpolate(x_e, T_eut, 1, Tm_B, comp_a)
                if T >= T_liquidus_right:
                    phase = "liquid"
                    phases = ["L"]
                    fractions = {"L": 1.0}
                    comps = {"L": comp_a}
                else:
                    c_liq = max(comp_a - 0.1, x_e)
                    c_beta = max_sol_AinB if max_sol_AinB is not None else 0.95

                    denom = c_beta - c_liq
                    if abs(denom) < 1e-10:
                        w_b, w_l = 0.5, 0.5
                    else:
                        w_l = (c_beta - comp_a) / denom
                        w_b = (comp_a - c_liq) / denom

                    phase = "two-phase (L + β)"
                    phases = ["L", "β"]
                    fractions = {"L": round(max(0, min(1, w_l)), 4), "β": round(max(0, min(1, w_b)), 4)}
                    comps = {"L": round(c_liq, 4), "β": round(c_beta, 4)}
        else:
            # 共晶温度以下
            if comp_a <= x_e:
                phase = "solid (α + eutectic mixture)"
                phases = ["α", "β (eutectic)"] if comp_a > (1 - (max_sol_BinA or 0.02)) else ["α"]
                # 杠杆定则
                c_alpha = 1.0 - (max_sol_BinA or 0.02)
                c_beta = max_sol_AinB or 0.95
                denom = c_beta - c_alpha
                if abs(denom) < 1e-10:
                    w_a, w_b = 0.5, 0.5
                else:
                    w_a = (c_beta - comp_a) / denom
                    w_b = (comp_a - c_alpha) / denom
                fractions = {"α": round(max(0, min(1, w_a)), 4), "β": round(max(0, min(1, w_b)), 4)}
                comps = {"α": round(c_alpha, 4), "β": round(c_beta, 4)}
            else:
                phase = "solid (β + eutectic mixture)"
                phases = ["β", "α (eutectic)"]
                c_alpha = 1.0 - (max_sol_BinA or 0.02)
                c_beta = max_sol_AinB or 0.95
                denom = c_beta - c_alpha
                if abs(denom) < 1e-10:
                    w_a, w_b = 0.5, 0.5
                else:
                    w_a = (c_beta - comp_a) / denom
                    w_b = (comp_a - c_alpha) / denom
                fractions = {"α": round(max(0, min(1, w_a)), 4), "β": round(max(0, min(1, w_b)), 4)}
                comps = {"α": round(c_alpha, 4), "β": round(c_beta, 4)}

        return {
            "current_phase": phase,
            "phases_present": phases,
            "phase_compositions": comps,
            "lever_rule_fractions": fractions,
            "eutectic_temperature_K": T_eut,
            "eutectic_composition_A": x_e,
            "explanation": f"Eutectic system at T={T}K, x_A={comp_a}: {phase}. "
                           f"Eutectic point: ({x_e}, {T_eut}K)."
        }

    def _analyze_peritectic(self, comp_a, T, Tm_A, Tm_B, T_peri, x_peri, max_sol_BinA, max_sol_AinB):
        """分析包晶相图（简化版，逻辑与共晶类似但反应不同）。"""
        T_p = T_peri or (Tm_A + Tm_B) / 3
        x_p = x_peri or 0.4

        if T >= max(Tm_A, Tm_B):
            phase = "all liquid"
            phases = ["L"]
            fractions = {"L": 1.0}
            comps = {"L": comp_a}
        elif T > T_p:
            if comp_a < x_p:
                T_liq = self._linear_interpolate(0, Tm_A, x_p, T_p, comp_a)
                if T >= T_liq:
                    phase = "liquid"
                    phases, fractions, comps = ["L"], {"L": 1.0}, {"L": comp_a}
                else:
                    c_liq = comp_a + 0.08
                    c_alpha = 0.95
                    denom = c_liq - c_alpha
                    w_a = (c_liq - comp_a) / denom if abs(denom) > 1e-10 else 0.5
                    w_l = (comp_a - c_alpha) / denom if abs(denom) > 1e-10 else 0.5
                    phase = "two-phase (L + α)"
                    phases = ["L", "α"]
                    fractions = {"L": round(w_l, 4), "α": round(w_a, 4)}
                    comps = {"L": round(c_liq, 4), "α": round(c_alpha, 4)}
            else:
                phase = "liquid"
                phases, fractions, comps = ["L"], {"L": 1.0}, {"L": comp_a}
        else:
            if comp_a < x_p:
                phase = "solid (α + β)"
                c_alpha, c_beta = 0.95, 0.1
                denom = c_alpha - c_beta
                w_a = (comp_a - c_beta) / denom if abs(denom) > 1e-10 else 0.5
                w_b = (c_alpha - comp_a) / denom if abs(denom) > 1e-10 else 0.5
                phases = ["α", "β"]
                fractions = {"α": round(w_a, 4), "β": round(w_b, 4)}
                comps = {"α": round(c_alpha, 4), "β": round(c_beta, 4)}
            else:
                phase = "solid (β)"
                phases, fractions, comps = ["β"], {"β": 1.0}, {"β": comp_a}

        return {
            "current_phase": phase,
            "phases_present": phases,
            "phase_compositions": comps,
            "lever_rule_fractions": fractions,
            "peritectic_temperature_K": T_p,
            "peritectic_composition_A": x_p,
            "explanation": f"Peritectic system at T={T}K, x_A={comp_a}: {phase}. "
                           f"Peritectic point: ({x_p}, {T_p}K)."
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            if len(parts) < 5:
                raise ValueError("Need: composition_A type T(K) Tm_A(K) Tm_B(K)")

            kwargs = {
                "composition_a": float(parts[0]),
                "phase_diagram_type": parts[1],
                "temperature_k": float(parts[2]),
                "melting_point_a_k": float(parts[3]),
                "melting_point_b_k": float(parts[4]),
            }
            if len(parts) > 5: kwargs["t_eutectic_or_peritectic"] = float(parts[5])
            if len(parts) > 6: kwargs["composition_eutectic"] = float(parts[6])
            if len(parts) > 7: kwargs["max_solubility_b_in_a"] = float(parts[7])
            if len(parts) > 8: kwargs["max_solubility_a_in_b"] = float(parts[8])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'x_A type T Tm_A Tm_B [T_eut x_eut sol_BinA sol_AinB]'"
            )
