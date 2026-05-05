"""
低共熔点（Eutectic Point）确定与分析工具
基于理想溶液模型计算二元体系的共晶温度和组成。
"""
import logging
import math
from typing import Dict, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EutecticPointFinder(BaseTool):
    """
    低共熔点确定工具。

    基于理想溶液的液相线方程，计算二元体系的共晶温度和共晶组成。
    使用公式: ln(x_i) = (ΔH_fus,i / R) × (1/T - 1/T_m,i)
    """
    __version__                 = "0.1.0"
    name                        = "EutecticPointFinder"
    func_name                   = "find_eutectic_point"
    description                 = "Determine the eutectic point (temperature and composition) of a binary system using ideal solution liquidus equations."
    implementation_description  = "Solves the intersection of two liquidus curves: ln(x_A) = ΔH_fus,A/R·(1/T - 1/Tm,A) and ln(x_B) = ΔH_fus,B/R·(1/T - 1/Tm,B), with x_A + x_B = 1. Uses numerical bisection to find T_eut and x_eut."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Eutectic Point", "Phase Diagram", "Materials Science", "Thermodynamics", "Binary System"]
    required_envs               = []

    code_input_sig = [
        ("component_a_name",            "str",   "N/A",     "Name of component A."),
        ("component_b_name",            "str",   "N/A",     "Name of component B."),
        ("melting_point_a_k",           "float", "N/A",     "Melting point of pure A in Kelvin."),
        ("melting_point_b_k",           "float", "N/A",     "Melting point of pure B in Kelvin."),
        ("enthalpy_of_fusion_a_jmol",   "float", "N/A",     "Enthalpy of fusion of A in J/mol."),
        ("enthalpy_of_fusion_b_jmol",   "float", "N/A",     "Enthalpy of fusion of B in J/mol."),
        ("ideal_eutectic_composition",  "float", "None",    "Optional known experimental eutectic composition (mass or mole fraction A). If provided, only calculates T_eut at this composition."),
        ("n_iterations",                "int",   "1000",    "Max iterations for numerical solver."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'A_name B_name Tm_A(K) Tm_B(K) dHfus_A(J/mol) dHfus_B(J/mol) [x_eut_known] [n_iter]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing eutectic_temperature_k, eutectic_composition_A, phase_diagram_summary, and details."),
    ]

    examples = [
        {
            "code_input": {
                "component_a_name": "Bi",
                "component_b_name": "Cd",
                "melting_point_a_k": 544.5,
                "melting_point_b_k": 594.2,
                "enthalpy_of_fusion_a_jmol": 10900.0,
                "enthalpy_of_fusion_b_jmol": 6100.0,
                "ideal_eutectic_composition": None,
                "n_iterations": 1000,
            },
            "text_input": {
                "input_params": "Bi Cd 544.5 594.2 10900.0 6100.0",
            },
            "output": {
                "result": {
                    "eutectic_temperature_k": "... (between melting points)",
                    "eutectic_composition_A": "... (0-1)",
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
        """初始化气体常数。"""
        self.R = 8.314462618  # J/(mol·K)

    def _liquidus_composition(self, T: float, Tm: float, dH_fus: float) -> float:
        """根据液相线方程计算给定温度下的组分摩尔分数。

        理想溶液液相线: ln(x) = (ΔH_fus/R) * (1/Tm - 1/T)
        当 T < Tm 时，ln(x) < 0，即 x < 1（液相中该组分的摩尔分数小于1）
        """
        if T <= 0 or T >= Tm:
            return None
        try:
            ln_x = (dH_fus / self.R) * (1.0 / Tm - 1.0 / T)
            x = math.exp(ln_x)
            return min(x, 1.0) if x > 0 else None
        except (OverflowError, ValueError):
            return None

    def _run_base(
        self,
        component_a_name: str,
        component_b_name: str,
        melting_point_a_k: float,
        melting_point_b_k: float,
        enthalpy_of_fusion_a_jmol: float,
        enthalpy_of_fusion_b_jmol: float,
        ideal_eutectic_composition: Optional[float] = None,
        n_iterations: int = 1000,
    ) -> Dict[str, Any]:
        """核心逻辑：求解共晶点。"""
        # ---- 输入验证 ----
        if melting_point_a_k <= 0 or melting_point_b_k <= 0:
            raise ChemMCPError("Melting points must be positive.")
        if enthalpy_of_fusion_a_jmol <= 0 or enthalpy_of_fusion_b_jmol <= 0:
            raise ChemMCPError("Enthalpies of fusion must be positive.")

        Tm_A = melting_point_a_k
        Tm_B = melting_point_b_k
        dH_A = enthalpy_of_fusion_a_jmol
        dH_B = enthalpy_of_fusion_b_jmol

        # 确保 Tm_A > Tm_B (A 为高熔点组分)
        if Tm_A < Tm_B:
            Tm_A, Tm_B = melting_point_b_k, melting_point_a_k
            dH_A, dH_B = enthalpy_of_fusion_b_jmol, enthalpy_of_fusion_a_jmol
            component_a_name, component_b_name = component_b_name, component_a_name
            logger.info("Swapped A/B so that A has higher melting point.")

        if ideal_eutectic_composition is not None:
            # 已知共晶组成，直接计算该组成对应的温度
            x_known = ideal_eutectic_composition
            if not (0 < x_known < 1):
                raise ChemMCPError("Eutectic composition must be between 0 and 1.")

            # 在两个液相线上分别求 T，取较低者（共晶温度）
            T_from_A = self._solve_T_for_x(x_known, Tm_A, dH_A)
            T_from_B = self._solve_T_for_x(1.0 - x_known, Tm_B, dH_B)

            T_eut = min(
                t for t in [T_from_A, T_from_B] if t is not None and t > 0
            ) if any(t for t in [T_from_A, T_from_B] if t is not None) else None

            return {
                "component_A": component_a_name,
                "component_B": component_b_name,
                "eutectic_temperature_K": round(T_eut, 2) if T_eut else None,
                "eutectic_composition_A": x_known,
                "method": "direct calculation at known composition",
                "phase_diagram_summary": (
                    f"Eutectic point for {component_a_name}-{component_b_name} system: "
                    f"x_A = {x_known:.3f}, T_eut ≈ {T_eut:.1f} K ({T_eut - 273.15:.1f} °C) "
                    f"(calculated at specified composition)."
                ),
            }

        # ---- 数值求解：找两条液相线的交点 ----
        # 搜索范围：T ∈ [min(Tm)*0.3, min(Tm_A, Tm_B)]
        T_low = min(Tm_A, Tm_B) * 0.3
        T_high = min(Tm_A, Tm_B) * 0.999

        best_T = None
        best_x = None
        best_residual = float("inf")

        # 使用更细的网格搜索
        for i in range(n_iterations):
            T = T_low + (T_high - T_low) * i / (n_iterations - 1)

            x_A_on_curve = self._liquidus_composition(T, Tm_A, dH_A)
            x_B_on_curve = self._liquidus_composition(T, Tm_B, dH_B)

            if x_A_on_curve is None or x_B_on_curve is None:
                continue

            residual = abs(x_A_on_curve + x_B_on_curve - 1.0)

            if residual < best_residual:
                best_residual = residual
                best_T = T
                best_x = x_A_on_curve

        # 用二分法精化（在最佳点附近搜索符号变化）
        if best_T is not None and best_residual > 1e-8:
            # 在 best_T 附近寻找 f(T) 符号变化的区间
            step = (T_high - T_low) / n_iterations * 5
            T_left = max(T_low, best_T - step * 50)
            T_right = min(T_high, best_T + step * 50)

            # 在更小区间内重新找符号变化
            f_left = self._f_total(T_left, Tm_A, dH_A, Tm_B, dH_B)
            f_right = self._f_total(T_right, Tm_A, dH_A, Tm_B, dH_B)

            if f_left is not None and f_right is not None:
                if f_left * f_right < 0:
                    # 存在根，二分法
                    refined = self._bisection_root(T_left, T_right, Tm_A, dH_A, Tm_B, dH_B)
                    if refined is not None:
                        x_A_final = self._liquidus_composition(refined, Tm_A, dH_A)
                        if x_A_final is not None:
                            best_T = refined
                            best_x = x_A_final
                            best_residual = abs(x_A_final + (self._liquidus_composition(refined, Tm_B, dH_B) or 0) - 1.0)

        logger.info(f"Eutectic point found: T={best_T:.1f}K, x_A={best_x:.4f}")

        return {
            "component_A": component_a_name,
            "component_B": component_b_name,
            "eutectic_temperature_K": round(best_T, 2) if best_T else None,
            "eutectic_composition_A": round(best_x, 6) if best_x else None,
            "eutectic_composition_B": round(1.0 - best_x, 6) if best_x else None,
            "residual_error": round(best_residual, 8),
            "method": "numerical intersection of liquidus curves",
            "input_data": {
                f"Tm_{component_a_name}_K": melting_point_a_k,
                f"Tm_{component_b_name}_K": melting_point_b_k,
                f"dHfus_{component_a_name}_J_mol": enthalpy_of_fusion_a_jmol,
                f"dHfus_{component_b_name}_J_mol": enthalpy_of_fusion_b_jmol,
            },
            "phase_diagram_summary": (
                f"Eutectic point for {component_a_name}-{component_b_name} system:\n"
                f"  • Eutectic temperature: {best_T:.1f} K ({best_T - 273.15:.1f} °C)\n"
                f"  • Eutectic composition: x_{component_a_name} = {best_x:.4f}, x_{component_b_name} = {1-best_x:.4f}\n"
                f"  • Method: Numerical solution of ideal liquidus curve intersection\n"
                f"  • Residual |x_A + x_B - 1|: {best_residual:.2e}"
            ) if best_T else "Failed to converge to a eutectic point.",
        }

    def _solve_T_for_x(self, x: float, Tm: float, dH_fus: float) -> Optional[float]:
        """反解液相线方程：已知 x 求 T。ln(x) = (ΔH_fus/R) * (1/Tm - 1/T)"""
        if x <= 0 or x >= 1:
            return None
        try:
            ln_x = math.log(x)
            # 1/T = 1/Tm - R·ln(x)/ΔH_fus
            inv_T = 1.0 / Tm - self.R * ln_x / dH_fus
            if inv_T > 0:
                return 1.0 / inv_T
        except (ValueError, ZeroDivisionError):
            pass
        return None

    def _f_total(self, T, Tm_A, dH_A, Tm_B, dH_B):
        """计算 f(T) = x_A(T) + x_B(T) - 1。"""
        x_A = self._liquidus_composition(T, Tm_A, dH_A)
        x_B = self._liquidus_composition(T, Tm_B, dH_B)
        if x_A is None or x_B is None:
            return None
        return x_A + x_B - 1.0

    def _bisection_root(self, T_a, T_b, Tm_A, dH_A, Tm_B, dH_B, max_iter=200, tol=1e-12):
        """二分法求 f(T)=0 的根。"""
        f_a = self._f_total(T_a, Tm_A, dH_A, Tm_B, dH_B)
        f_b = self._f_total(T_b, Tm_A, dH_A, Tm_B, dH_B)
        if f_a is None or f_b is None:
            return None

        for _ in range(max_iter):
            T_mid = (T_a + T_b) / 2.0
            f_mid = self._f_total(T_mid, Tm_A, dH_A, Tm_B, dH_B)
            if f_mid is None:
                break
            if abs(f_mid) < tol or (T_b - T_a) < 1e-12:
                return T_mid
            if f_a * f_mid < 0:
                T_b = T_mid
                f_b = f_mid
            else:
                T_a = T_mid
                f_a = f_mid
        return (T_a + T_b) / 2.0

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            if len(parts) < 6:
                raise ValueError("Need: A_name B_name Tm_A(K) Tm_B(K) dHfus_A(J/mol) dHfus_B(J/mol)")

            kwargs = {
                "component_a_name": parts[0],
                "component_b_name": parts[1],
                "melting_point_a_k": float(parts[2]),
                "melting_point_b_k": float(parts[3]),
                "enthalpy_of_fusion_a_jmol": float(parts[4]),
                "enthalpy_of_fusion_b_jmol": float(parts[5]),
            }
            if len(parts) > 6: kwargs["ideal_eutectic_composition"] = float(parts[6])
            if len(parts) > 7: kwargs["n_iterations"] = int(parts[7])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'A_name B_name Tm_A Tm_B dHfus_A dHfus_B [x_eut] [n_iter]'"
            )
