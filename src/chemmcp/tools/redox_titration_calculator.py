import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RedoxTitrationCalculator(BaseTool):
    """
    氧化还原滴定计算工具：电位曲线、当量点电位、Nernst方程。
    
    支持常见氧化还原体系（Ce⁴⁺/Fe²⁺、MnO₄⁻/Fe²⁺、Cr₂O₇²⁻/Fe²⁺等）。
    """
    __version__ = "0.1.0"
    name             = "RedoxTitrationCalculator"
    func_name        = "redox_titration"
    description      = "Calculate redox titration curves, equivalence point potential, and electrode potential using Nernst equation."
    implementation_description = "Applies Nernst equation E = E°' + (RT/nF)ln(ox/red) for each half-reaction. Computes system potential at each point of the titration curve. R=8.314 J/(mol·K), F=96485 C/mol, T=298.15K by default."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Redox Titration", "Nernst Equation", "Potentiometric Titration", "Analytical Chemistry"]
    required_envs    = []

    # 常见氧化还原电对的标准电极电势 (V, vs SHE)
    REDOX_COUPLES = {
        "Ce4_Ce3":     {"name": "Ce⁴⁺/Ce³⁺",      "E0": 1.61,   "n": 1, "ox": "Ce⁴⁺", "red": "Ce³⁺"},
        "Fe3_Fe2":     {"name": "Fe³⁺/Fe²⁺",       "E0": 0.771,  "n": 1, "ox": "Fe³⁺", "red": "Fe²⁺"},
        "MnO4_Mn2":    {"name": "MnO₄⁻/Mn²⁺",      "E0": 1.51,   "n": 5, "ox": "MnO₄⁻", "red": "Mn²⁺"},
        "Cr2O7_3Cr3":  {"name": "Cr₂O₇²⁻/Cr³⁺",    "E0": 1.33,   "n": 6, "ox": "Cr₂O₇²⁻", "red": "Cr³⁺"},
        "I2_2I":       {"name": "I₂/2I⁻",           "E0": 0.5355, "n": 2, "ox": "I₂", "red": "I⁻"},
        "Br2_2Br":     {"name": "Br₂/2Br⁻",         "E0": 1.066,  "n": 2, "ox": "Br₂", "red": "Br⁻"},
        "Sn4_Sn2":     {"name": "Sn⁴⁺/Sn²⁺",        "E0": 0.154,  "n": 2, "ox": "Sn⁴⁺", "red": "Sn²⁺"},
        "VO2_VO":      {"name": "VO₂⁺/VO²⁺",        "E0": 1.00,   "n": 1, "ox": "VO₂⁺", "red": "VO²⁺"},
        "NO3_NO":      {"name": "NO₃⁻/NO",          "E0": 0.96,   "n": 3, "ox": "NO₃⁻", "red": "NO"},
        "H2O2_H2O":    {"name": "H₂O₂/H₂O",         "E0": 1.776,  "n": 2, "ox": "H₂O₂", "red": "H₂O"},
        "Fe2_Fe":      {"name": "Fe²⁺/Fe",           "E0": -0.447, "n": 2, "ox": "Fe²⁺", "red": "Fe"},
        "Cu2_Cu":      {"name": "Cu²⁺/Cu",           "E0": 0.337,  "n": 2, "ox": "Cu²⁺", "red": "Cu"},
        "Ag_Ag":       {"name": "Ag⁺/Ag",            "E0": 0.7996, "n": 1, "ox": "Ag⁺", "red": "Ag"},
        "Sce":         {"name": "SCE (reference)",   "E0": 0.242,  "n": 0, "ox": "", "red": ""},
    }

    R = 8.314       # J/(mol·K)
    F = 96485       # C/mol
    T_DEFAULT = 298.15  # K (25°C)

    code_input_sig   = [
        ("oxidant_couple_key", "str", "N/A", "Key of oxidant couple from database (e.g., 'Ce4_Ce3', 'MnO4_Mn2', 'Cr2O7_3Cr3')."),
        ("reductant_couple_key", "str", "N/A", "Key of reductant couple from database (e.g., 'Fe3_Fe2')."),
        ("reductant_concentration_mol_L", "float", "N/A", "Initial concentration of reductant (analyte) in mol/L."),
        ("reductant_volume_ml", "float", "N/A", "Volume of reductant solution (mL)."),
        ("oxidant_concentration_mol_L", "float", "N/A", "Concentration of oxidant titrant in mol/L."),
        ("temperature_k_or_None", "float_or_None", "None", "Temperature in Kelvin (default: 298.15 K / 25°C)."),
        ("n_points", "int", "100", "Number of points on the titration curve."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'oxidant_key reductant_key C_red V_red C_ox [T_K] [n_points]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with equivalence data, potential at key points, full titration curve (volume vs E), and reaction stoichiometry."),
    ]

    examples         = [
        {
            "code_input": {
                "oxidant_couple_key": "Ce4_Ce3",
                "reductant_couple_key": "Fe3_Fe2",
                "reductant_concentration_mol_L": 0.10,
                "reductant_volume_ml": 25.0,
                "oxidant_concentration_mol_L": 0.10,
                "temperature_k_or_None": None,
                "n_points": 50,
            },
            "text_input": {
                "input_params": "Ce4_Ce3 Fe3_Fe2 0.10 25.0 0.10 50",
            },
            "output": {
                "result": {
                    "equivalence_point": {"volume_ml": "...", "potential_V": "..."},
                    "curve_data": [{"volume_ml": ..., "potential_V": ...}, ...],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _nernst(self, E0: float, n: int, Q: float, T: float) -> float:
        """Nernst方程: E = E° + (RT/nF) * ln(Q)"""
        if n == 0 or Q <= 0:
            return E0
        return E0 + (self.R * T) / (n * self.F) * math.log(Q)

    def _run_base(
        self,
        oxidant_couple_key: str,
        reductant_couple_key: str,
        reductant_concentration_mol_L: float,
        reductant_volume_ml: float,
        oxidant_concentration_mol_L: float,
        temperature_k_or_None: Optional[float] = None,
        n_points: int = 100,
    ) -> dict:
        """核心逻辑：氧化还原滴定曲线计算"""
        if oxidant_couple_key not in self.REDOX_COUPLES:
            raise ChemMCPError(f"Unknown oxidant couple '{oxidant_couple_key}'. Available: {list(self.REDOX_COUPLES.keys())}")
        if reductant_couple_key not in self.REDOX_COUPLES:
            raise ChemMCPError(f"Unknown reductant couple '{reductant_couple_key}'. Available: {list(self.REDOX_COUPLES.keys())}")

        ox_couple = self.REDOX_COUPLES[oxidant_couple_key]
        red_couple = self.REDOX_COUPLES[reductant_couple_key]

        T = temperature_k_or_None if temperature_k_or_None else self.T_DEFAULT

        Cr = float(reductant_concentration_mol_L)
        Vr = float(reductant_volume_ml)
        Co = float(oxidant_concentration_mol_L)

        n_ox = ox_couple["n"]  # 氧化剂半反应电子数
        n_red = red_couple["n"]  # 还原剂半反应电子数

        # 化学计量关系：n_red × mol_reductant = n_ox × mol_oxidant
        Ve = Cr * Vr * n_red / (Co * n_ox) if Co > 0 and n_ox > 0 else 0

        RT_over_F = self.R * T / self.F  # ≈ 0.02569 V at 298K

        # ---- 滴定曲线 ----
        curve_data = []
        V_max = Ve * 1.4

        for i in range(n_points + 1):
            Vox = V_max * i / n_points
            Vtotal = Vr + Vox
            if Vtotal < 1e-12:
                Vtotal = 1e-10

            if Vox < 1e-12:
                # 滴定前：纯还原态
                E = red_couple["E0"]  # 近似为标准电势（无氧化态时无法精确计算）
            elif Vox < Ve - 1e-10:
                # 当量点前：用还原剂电对的 Nernst 方程
                f = Vox / Ve  # 反应进度比例
                # 剩余还原剂 : 生成的氧化态产物 = (1-f) : f
                ratio_ox_red = f / (1 - f) if (1 - f) > 1e-20 else 1e20
                E = self._nernst(red_couple["E0"], n_red, ratio_ox_red, T)
            elif abs(Vox - Ve) < 1e-10:
                # 当量点电位
                E_eq = (n_ox * ox_couple["E0"] + n_red * red_couple["E0"]) / (n_ox + n_red)
                E = E_eq
            else:
                # 过量氧化剂：用氧化剂电对的 Nernst 方程
                excess_f = (Vox - Ve) / Ve
                ratio_ox_red = excess_f if excess_f > 1e-20 else 1e20
                E = self._nernst(ox_couple["E0"], n_ox, ratio_ox_red, T)

            curve_data.append({
                "volume_ml": round(Vox, 4),
                "potential_V_she": round(E, 6),
                "potential_V_vs_sce": round(E - 0.242, 6),  # 转换为 vs SCE
            })

        # 关键点电位
        E_eq_point = (n_ox * ox_couple["E0"] + n_red * red_couple["E0"]) / (n_ox + n_red)
        E_half_eq_reductant = red_couple["E0"]
        E_half_eq_oxidant = ox_couple["E0"]

        result = {
            "titration_system": f"{ox_couple['name']} titrating {red_couple['name']}",
            "temperature_K": round(T, 2),
            "RT_F_V": round(RT_over_F, 6),
            "half_reactions": {
                "oxidant": ox_couple,
                "reductant": red_couple,
            },
            "stoichiometry": {
                "n_oxidant_electrons": n_ox,
                "n_reductant_electrons": n_red,
                "equivalence_volume_ml": round(Ve, 6),
            },
            "equivalence_point": {
                "volume_ml": round(Ve, 6),
                "potential_V_SHE": round(E_eq_point, 6),
                "potential_V_SCE": round(E_eq_point - 0.242, 6),
                "formula": f"E_eq = ({n_ox}×{ox_couple['E0']} + {n_red}×{red_couple['E0']}) / ({n_ox}+{n_red})",
            },
            "key_potentials": {
                "initial": round(red_couple["E0"], 4),
                "half_equivalence_reductant_side": round(E_half_eq_reductant, 4),
                "equivalence": round(E_eq_point, 6),
                "half_equivalence_oxidant_side": round(E_half_eq_oxidant, 4),
                "excess_10pct": round(
                    self._nernst(ox_couple["E0"], n_ox, 0.10, T), 4
                ),
            },
            "titration_curve": curve_data,
        }

        logger.info(f"Redox titration: {ox_couple['name']} + {red_couple['name']}, Ve={Ve:.2f}mL, Eeq={E_eq_point:.3f}V")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            ox_key = parts[0]
            red_key = parts[1]
            Cr = float(parts[2])
            Vr = float(parts[3])
            Co = float(parts[4])
            T = float(parts[5]) if len(parts) > 5 and parts[5].lower() != "none" else None
            np_ = int(parts[6]) if len(parts) > 6 else 100
            return self._run_base(ox_key, red_key, Cr, Vr, Co, T, np_)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'oxidant_key reductant_key C_red V_red C_ox [T_K] [n_points]'")
