import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AcidBaseTitrationCalculator(BaseTool):
    """
    酸碱滴定计算工具：当量点、pH曲线、指示剂选择建议。
    
    支持强酸/弱酸与强碱/弱碱的滴定体系，计算完整滴定曲线。
    """
    __version__ = "0.1.0"
    name             = "AcidBaseTitrationCalculator"
    func_name        = "calculate_acid_base_titration"
    description      = "Calculate acid-base titration curves, equivalence point, and pH at key points for various titration systems."
    implementation_description = "Solves acid-base equilibrium for strong/strong, strong/weak, weak/strong, weak/weak titration systems. Computes pH curve data points using charge balance and mass balance equations."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Acid-Base Titration", "Equivalence Point", "pH Curve", "Analytical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("acid_type", "str", "N/A", "Type of analyte: 'strong' or 'weak'."),
        ("base_type", "str", "N/A", "Type of titrant: 'strong' or 'weak'."),
        ("acid_concentration", "float", "N/A", "Initial concentration of acid (mol/L)."),
        ("acid_volume_ml", "float", "N/A", "Initial volume of acid solution (mL)."),
        ("base_concentration", "float", "N/A", "Concentration of base titrant (mol/L)."),
        ("ka_or_None", "float_or_None", "None", "Ka value if acid is weak (e.g., 1.8e-5 for acetic acid)."),
        ("kb_or_None", "float_or_None", "None", "Kb value if base is weak (e.g., 1.8e-5 for ammonia)."),
        ("n_points", "int", "100", "Number of points on the titration curve."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'acid_type base_type Ca Va Cb [ka] [kb] [n_points]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with equivalence point, pH curve data, key pH values, and indicator suggestion."),
    ]

    examples         = [
        {
            "code_input": {
                "acid_type": "weak",
                "base_type": "strong",
                "acid_concentration": 0.10,
                "acid_volume_ml": 25.0,
                "base_concentration": 0.10,
                "ka_or_None": 1.8e-5,
                "kb_or_None": None,
                "n_points": 50,
            },
            "text_input": {
                "input_params": "weak strong 0.10 25.0 0.10 1.8e-5 50",
            },
            "output": {
                "result": {
                    "equivalence_point": {"volume_ml": "...", "ph": "..."},
                    "key_ph_values": {...},
                    "curve_data": [{"volume_ml": ..., "ph": ...}, ...],
                }
            },
        },
    ]

    # Kw at 25°C
    KW = 1.0e-14

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        acid_type: str,
        base_type: str,
        acid_concentration: float,
        acid_volume_ml: float,
        base_concentration: float,
        ka_or_None: Optional[float] = None,
        kb_or_None: Optional[float] = None,
        n_points: int = 100,
    ) -> dict:
        """核心逻辑：酸碱滴定曲线计算"""
        # 输入验证
        if acid_type not in ("strong", "weak"):
            raise ChemMCPError("acid_type must be 'strong' or 'weak'")
        if base_type not in ("strong", "weak"):
            raise ChemMCPError("base_type must be 'strong' or 'weak'")
        if acid_concentration <= 0 or base_concentration <= 0:
            raise ChemMCPError("Concentrations must be positive.")
        if acid_volume_ml <= 0:
            raise ChemMCPError("Volume must be positive.")

        Ca = float(acid_concentration)
        Va = float(acid_volume_ml)
        Cb = float(base_concentration)
        ka = ka_or_None
        kb = kb_or_None

        # 当量点体积
        Ve = Ca * Va / Cb

        # 计算滴定曲线
        curve_data = []
        V_max = Ve * 1.4  # 计算到140%当量点
        V_min = 0.0

        for i in range(n_points + 1):
            Vb = V_min + (V_max - V_min) * i / n_points
            ph = self._calculate_ph(Vb, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve)
            curve_data.append({
                "volume_ml": round(Vb, 4),
                "ph": round(ph, 4),
            })

        # 关键点的 pH 值
        key_ph = self._get_key_ph(Ca, Va, Cb, acid_type, base_type, ka, kb, Ve)

        # 当量点 pH
        eq_ph = self._calculate_ph(Ve, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve)

        # 指示剂建议
        indicator_suggestion = self._suggest_indicator(eq_ph)

        result = {
            "titration_system": f"{acid_type.capitalize()} acid + {base_type.capitalize()} base",
            "parameters": {
                "acid_conc_mol_L": Ca,
                "acid_volume_ml": Va,
                "base_conc_mol_L": Cb,
                "ka": ka,
                "kb": kb,
            },
            "equivalence_point": {
                "volume_ml": round(Ve, 6),
                "ph": round(eq_ph, 4),
            },
            "key_ph_values": key_ph,
            "indicator_suggestion": indicator_suggestion,
            "curve_data": curve_data,
        }

        logger.info(f"Acid-base titration: {acid_type}+{base_type}, Ve={Ve:.2f}mL, pH_eq={eq_ph:.2f}")
        return result

    def _calculate_ph(self, Vb: float, Ca: float, Va: float, Cb: float,
                      acid_type: str, base_type: str, ka: Optional[float],
                      kb: Optional[float], Ve: float) -> float:
        """计算加入 Vb mL 碱后的 pH"""
        total_vol = Va + Vb
        if total_vol == 0:
            total_vol = 1e-10

        if acid_type == "strong" and base_type == "strong":
            return _ph_strong_strong(Vb, Ca, Va, Cb, Ve, total_vol, self.KW)
        elif acid_type == "strong" and base_type == "weak":
            return _ph_strong_weak(Vb, Ca, Va, Cb, Ve, total_vol, self.KW, kb)
        elif acid_type == "weak" and base_type == "strong":
            return _ph_weak_strong(Vb, Ca, Va, Cb, Ve, total_vol, self.KW, ka)
        else:  # weak + weak
            return _ph_weak_weak(Vb, Ca, Va, Cb, Ve, total_vol, self.KW, ka, kb)

    def _get_key_ph(self, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve):
        """获取关键点的 pH"""
        keys = {}
        # 初始 pH (Vb=0)
        keys["initial"] = round(self._calculate_ph(0, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve), 4)
        # 半当量点
        keys["half_equivalence"] = round(self._calculate_ph(Ve / 2, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve), 4)
        # 当量点
        keys["equivalence"] = round(self._calculate_ph(Ve, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve), 4)
        # 过量 10%
        keys["excess_10pct"] = round(self._calculate_ph(Ve * 1.1, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve), 4)
        # 过量 20%
        keys["excess_20pct"] = round(self._calculate_ph(Ve * 1.2, Ca, Va, Cb, acid_type, base_type, ka, kb, Ve), 4)
        return keys

    def _suggest_indicator(self, eq_ph: float) -> dict:
        """根据当量点 pH 建议指示剂"""
        indicators = [
            {"name": "Thymol Blue (1st)", "range_low": 1.2, "range_high": 2.8, "color_change": "Red → Yellow"},
            {"name": "Methyl Orange", "range_low": 3.1, "range_high": 4.4, "color_change": "Red → Yellow"},
            {"name": "Bromocresol Green", "range_low": 3.8, "range_high": 5.4, "color_change": "Yellow → Blue"},
            {"name": "Methyl Red", "range_low": 4.4, "range_high": 6.2, "color_change": "Red → Yellow"},
            {"name": "Bromothymol Blue", "range_low": 6.0, "range_high": 7.6, "color_change": "Yellow → Blue"},
            {"name": "Phenol Red", "range_low": 6.4, "range_high": 8.0, "color_change": "Yellow → Red"},
            {"name": "Cresol Red", "range_low": 7.2, "range_high": 8.8, "color_change": "Yellow → Red"},
            {"name": "Phenolphthalein", "range_low": 8.2, "range_high": 10.0, "color_change": "Colorless → Pink"},
            {"name": "Thymolphthalein", "range_low": 9.3, "range_high": 10.5, "color_change": "Colorless → Blue"},
        ]
        suitable = []
        for ind in indicators:
            if ind["range_low"] <= eq_ph <= ind["range_high"]:
                suitable.append(ind)

        best = suitable[0] if suitable else min(indicators, key=lambda x: abs((x["range_low"]+x["range_high"])/2 - eq_ph))
        return {
            "equivalence_ph": round(eq_ph, 4),
            "best_indicator": best["name"],
            "transition_range": f"{best['range_low']}-{best['range_high']}",
            "color_change": best["color_change"],
            "all_suitable": [s["name"] for s in suitable] if suitable else [best["name"]],
        }

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            acid_t = parts[0]
            base_t = parts[1]
            Ca = float(parts[2])
            Va = float(parts[3])
            Cb = float(parts[4])
            ka = float(parts[5]) if len(parts) > 5 and parts[5].lower() != "none" else None
            kb = float(parts[6]) if len(parts) > 6 and parts[6].lower() != "none" else None
            np_ = int(parts[7]) if len(parts) > 7 else 100
            return self._run_base(acid_t, base_t, Ca, Va, Cb, ka, kb, np_)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'acid_type base_type Ca Va Cb [ka] [kb] [n_points]'")


# ============================================================
# 内部函数：各种滴定体系的 pH 计算
# ============================================================

def _ph_strong_strong(Vb, Ca, Va, Cb, Ve, Vtot, Kw):
    """强酸 + 强碱"""
    excess_H = Ca * Va / Vtot - Cb * Vb / Vtot
    if abs(excess_H) < 1e-12:
        return 7.0  # 中性
    elif excess_H > 0:
        return -math.log10(excess_H)
    else:
        return 14.0 + math.log10(-excess_H)


def _ph_strong_weak(Vb, Ca, Va, Cb, Ve, Vtot, Kw, kb):
    """强酸 + 弱碱（如 HCl 滴定 NH3）"""
    if Vb < 1e-12:
        return -math.log10(Ca) if Ca > 0 else 7.0

    mol_acid = Ca * Va  # 强酸的物质的量
    mol_base = Cb * Vb  # 弱碱的物质的量

    if mol_base <= mol_acid:
        # 剩余强酸
        excess_H = (mol_acid - mol_base) / Vtot
        if excess_H > 1e-10:
            return -math.log10(excess_H)
        return 7.0
    else:
        # 过量弱碱
        excess_base_mol = mol_base - mol_acid
        C_b = excess_base_mol / Vtot
        if kb is None or kb <= 0:
            return 14.0  # 默认
        OH = math.sqrt(kb * C_b) if C_b > 0 else 0
        if OH > 1e-10:
            return 14.0 + math.log10(OH)
        return 7.0


def _ph_weak_strong(Vb, Ca, Va, Cb, Ve, Vtot, Kw, ka):
    """弱酸 + 强碱（最常见）"""
    if Vb < 1e-12:
        # 纯弱酸溶液
        if ka is None or ka <= 0:
            return 7.0
        H = math.sqrt(ka * Ca) if Ca > 0 else 0
        return -math.log10(H) if H > 1e-10 else 7.0 - 0.5 * math.log10(Kw / Ca) if Ca > 0 else 7.0

    mol_acid = Ca * Va
    mol_base_added = Cb * Vb

    if mol_base_added < mol_acid:
        # 缓冲区：HA + A-
        fraction = mol_base_added / mol_acid  # 被中和的比例
        if ka is None or ka <= 0 or fraction >= 1.0:
            return 7.0
        # Henderson-Hasselbalch: pH = pKa + log([A-]/[HA])
        pKa = -math.log10(ka)
        if fraction > 0 and fraction < 1.0:
            ratio = fraction / (1.0 - fraction)
            return pKa + math.log10(ratio)
        elif fraction <= 0:
            H = math.sqrt(ka * Ca * (1 - fraction))
            return -math.log10(H) if H > 1e-10 else 7.0
        else:
            return 7.0
    elif abs(mol_base_added - mol_acid) < 1e-12:
        # 当量点：生成共轭碱 A-
        C_salt = mol_acid / Vtot
        if ka is None or ka <= 0:
            return 7.0
        Kb = Kw / ka
        OH = math.sqrt(Kb * C_salt) if C_salt > 0 else 0
        if OH > 1e-10:
            return 14.0 + math.log10(OH)
        return 7.0
    else:
        # 过量强碱
        excess_OH = (mol_base_added - mol_acid) / Vtot
        if excess_OH > 1e-10:
            return 14.0 + math.log10(excess_OH)
        return 7.0


def _ph_weak_weak(Vb, Ca, Va, Cb, Ve, Vtot, Kw, ka, kb):
    """弱酸 + 弱碱（较少见）"""
    if Vb < 1e-12:
        if ka and ka > 0:
            H = math.sqrt(ka * Ca)
            return -math.log10(H) if H > 1e-10 else 7.0
        return 7.0

    mol_acid = Ca * Va
    mol_base = Cb * Vb

    if mol_base < mol_acid:
        if ka and ka > 0:
            fraction = mol_base / mol_acid
            pKa = -math.log10(ka)
            if 0 < fraction < 1:
                return pKa + math.log10(fraction / (1 - fraction))
            H = math.sqrt(ka * Ca)
            return -math.log10(H) if H > 1e-10 else 7.0
        return 7.0
    elif abs(mol_base - mol_acid) < 1e-12:
        # 近似中性偏碱性（取决于 Ka, Kb）
        if ka and kb:
            return 7.0 + 0.5 * (math.log10(kb) + math.log10(ka) - math.log10(Kw))
        return 7.0
    else:
        if kb and kb > 0:
            excess_Cb = (mol_base - mol_acid) / Vtot
            OH = math.sqrt(kb * excess_Cb)
            return 14.0 + math.log10(OH) if OH > 1e-10 else 7.0
        return 7.0
