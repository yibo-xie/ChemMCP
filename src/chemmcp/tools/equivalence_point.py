import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EquivalencePoint(BaseTool):
    """
    计算滴定等当点（Equivalence Point）的 pH 值。
    包括等当点体积、pH、指示剂建议和 pH 突跃范围。
    """
    __version__ = "0.1.0"
    name = "EquivalencePoint"
    func_name = "calculate_equivalence_point"
    description = "Calculate the equivalence point of an acid-base titration: volume, pH, indicator suggestion, and pH jump range."
    implementation_description = "Determines equivalence point properties based on titration type (strong/weak acid-base combinations) using equilibrium calculations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Titration", "Equivalence Point", "Indicator", "Acid-Base"]
    required_envs = []

    code_input_sig = [
        ("analyte_type", "str", "N/A", "Analyte type: 'strong_acid', 'weak_acid', 'strong_base', 'weak_base'."),
        ("analyte_conc", "float", "N/A", "Analyte concentration (mol/L)."),
        ("analyte_volume_ml", "float", "N/A", "Analyte volume (mL)."),
        ("titrant_type", "str", "N/A", "Titrant type: 'strong_base' or 'strong_acid'."),
        ("titrant_conc", "float", "N/A", "Titrant concentration (mol/L)."),
        ("Ka", "float", "None", "Ka for weak acids."),
        ("Kb", "float", "None", "Kb for weak bases."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'atype conc vol_ml ttype tconc [Ka] [Kb]'. Example: 'weak_acid 0.1 25 strong_base 0.1 1.8e-5'"),
    ]

    output_sig = [
        ("eq_volume_ml", "float", "Titrant volume at equivalence point (mL)."),
        ("ph_at_eq", "float", "pH at the equivalence point."),
        ("indicator_suggestion", "str", "Recommended indicator(s) with color transition range."),
        ("ph_jump_low", "float", "Lower bound of the pH jump region (±0.1% from eq point)."),
        ("ph_jump_high", "float", "Upper bound of the pH jump region."),
        ("titration_summary", "str", "Summary of the titration system and key findings."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_type": "weak_acid",
                "analyte_conc": 0.1,
                "analyte_volume_ml": 25,
                "titrant_type": "strong_base",
                "titrant_conc": 0.1,
                "Ka": 1.8e-5,
                "Kb": None,
            },
            "text_input": {
                "input_params": "weak_acid 0.1 25 strong_base 0.1 1.8e-5"
            },
            "output": {
                "eq_volume_ml": 25.0,
                "ph_at_eq": 8.72,
                "indicator_suggestion": "Phenolphthalein (pH 8.2–10.0)",
                "ph_jump_low": 7.8,
                "ph_jump_high": 9.7,
                "titration_summary": "Weak acid + strong base titration.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Kw = 1.0e-14

        # 常用指示剂数据库: (名称, 变色范围低, 变色范围高, 酸色→碱色)
        self.indicators = [
            ("Methyl Orange", 3.1, 4.4, "Red → Yellow"),
            ("Methyl Red", 4.4, 6.2, "Red → Yellow"),
            ("Bromothymol Blue", 6.0, 7.6, "Yellow → Blue"),
            ("Phenol Red", 6.4, 8.0, "Yellow → Red"),
            ("Phenolphthalein", 8.2, 10.0, "Colorless → Pink"),
            ("Thymolphthalein", 9.3, 10.5, "Colorless → Blue"),
        ]

    def _run_base(self, analyte_type: str, analyte_conc: float, analyte_volume_ml: float,
                  titrant_type: str, titrant_conc: float, Ka: float = None, Kb: float = None) -> dict:
        """核心逻辑：计算等当点"""
        if analyte_conc <= 0 or titrant_conc <= 0:
            raise ChemMCPError("Concentrations must be positive.")
        if analyte_volume_ml <= 0:
            raise ChemMCPError("Volume must be positive.")

        Va = analyte_volume_ml / 1000.0  # L
        Ca = analyte_conc
        Ct = titrant_conc

        # 等当点体积
        V_eq_ml = (Ca * Va) / Ct * 1000
        V_total_L = (analyte_volume_ml + V_eq_ml) / 1000.0

        atype = analyte_type.lower().strip()
        ttype = titrant_type.lower().strip()

        # 计算 pH
        if atype in ("strong_acid", "sa") and ttype in ("strong_base", "sb"):
            ph_eq = 7.0
            summary = "Strong acid + strong base → neutral salt solution at eq point."
        elif atype in ("weak_acid", "wa") and ttype in ("strong_base", "sb"):
            if Ka is None:
                raise ChemMCPError("Ka required for weak_acid + strong_base.")
            ph_eq = self._calc_ph_wa_sb_eq(Ca, Va, V_total_L, Ka)
            summary = (
                f"Weak acid + strong base → conjugate base hydrolyzes at eq point. "
                f"[A-] = {Ca*Va/V_total_L:.4f} M, Kh = Kw/Ka = {self.Kw/Ka:.2e}, pH > 7."
            )
        elif atype in ("strong_base", "sb") and ttype in ("strong_acid", "sa"):
            ph_eq = 7.0
            summary = "Strong base + strong acid → neutral salt solution at eq point."
        elif atype in ("weak_base", "wb") and ttype in ("strong_acid", "sa"):
            if Kb is None:
                raise ChemMCPError("Kb required for weak_base + strong_acid.")
            ph_eq = self._calc_ph_wb_sa_eq(Ca, Va, V_total_L, Kb)
            summary = (
                f"Weak base + strong acid → conjugate acid dissociates at eq point. "
                f"[BH+] = {Ca*Va/V_total_L:.4f} M, Ka(conj) = Kw/Kb = {self.Kw/Kb:.2e}, pH < 7."
            )
        else:
            raise ChemMCPError(f"Unsupported combination: {atype} + {ttype}")

        # pH 突跃范围（等当点前后各 0.1%）
        ph_low, ph_high = self._calc_ph_jump_range(
            atype, ttype, Ca, Va, Ct, V_eq_ml, Ka, Kb
        )

        # 指示剂选择
        indicator = self._suggest_indicator(ph_low, ph_high)

        logger.info(f"Equivalence point: V_eq={V_eq_ml:.2f} mL, pH={ph_eq:.2f}")

        return {
            "eq_volume_ml": round(V_eq_ml, 4),
            "ph_at_eq": round(ph_eq, 4),
            "indicator_suggestion": indicator,
            "ph_jump_low": round(ph_low, 2),
            "ph_jump_high": round(ph_high, 2),
            "titration_summary": summary,
        }

    def _calc_ph_wa_sb_eq(self, Ca: float, Va: float, V_total: float, Ka: float) -> float:
        """弱酸+强碱 等当点：纯 A- 水解"""
        c_A = Ca * Va / V_total
        Kh = self.Kw / Ka
        disc = Kh ** 2 + 4 * Kh * c_A
        oh = (-Kh + math.sqrt(disc)) / 2
        poh = -math.log10(oh) if oh > 0 else 7.0
        return 14 - poh

    def _calc_ph_wb_sa_eq(self, Ca: float, Va: float, V_total: float, Kb: float) -> float:
        """弱碱+强酸 等当点：纯 BH+ 解离"""
        c_BH = Ca * Va / V_total
        Ka_conj = self.Kw / Kb
        disc = Ka_conj ** 2 + 4 * Ka_conj * c_BH
        h = (-Ka_conj + math.sqrt(disc)) / 2
        return -math.log10(h) if h > 0 else 7.0

    def _calc_ph_jump_range(self, atype: str, ttype: str, Ca: float, Va: float,
                              Ct: float, V_eq: float, Ka: float, Kb: float) -> tuple:
        """计算 pH 突跃范围（±0.1% 等当量）"""
        # 使用 TitrationCurve 的方法来计算
        try:
            # -0.1%
            V_99 = V_eq * 0.999
            V_99_L = V_99 / 1000.0
            V_tot_99 = (Va * 1000 + V_99) / 1000.0
            if atype == "weak_acid":
                ph_lo, _ = self._ph_wa_sb(Ca, Va/1000.0, Ct, V_99_L, V_tot_99, Ka, V_eq)
            elif atype == "weak_base":
                ph_lo, _ = self._ph_wb_sa(Ca, Va/1000.0, Ct, V_99_L, V_tot_99, Kb, V_eq)
            else:
                ph_lo, _ = self._ph_sa_sb(Ca, Va/1000.0, Ct, V_99_L, V_tot_99)

            # +0.1%
            V_101 = V_eq * 1.001
            V_101_L = V_101 / 1000.0
            V_tot_101 = (Va * 1000 + V_101) / 1000.0
            if atype == "weak_acid":
                ph_hi, _ = self._ph_wa_sb(Ca, Va/1000.0, Ct, V_101_L, V_tot_101, Ka, V_eq)
            elif atype == "weak_base":
                ph_hi, _ = self._ph_wb_sa(Ca, Va/1000.0, Ct, V_101_L, V_tot_101, Kb, V_eq)
            else:
                ph_hi, _ = self._ph_sa_sb(Ca, Va/1000.0, Ct, V_101_L, V_tot_101)

            return ph_lo, ph_hi
        except Exception:
            return 3.0, 11.0

    def _suggest_indicator(self, ph_low: float, ph_high: float) -> str:
        """选择合适的指示剂"""
        best = []
        for name, lo, hi, transition in self.indicators:
            # 指示剂变色范围应在突跃范围内或接近
            center = (lo + hi) / 2
            jump_center = (ph_low + ph_high) / 2
            if lo >= ph_low - 0.5 and hi <= ph_high + 0.5:
                best.append((name, lo, hi, transition, abs(center - jump_center)))

        if best:
            best.sort(key=lambda x: x[4])
            name, lo, hi, trans, _ = best[0]
            return f"{name} (pH {lo}–{hi}, {trans})"
        else:
            # 找最接近的
            jump_center = (ph_low + ph_high) / 2
            closest = min(self.indicators, key=lambda x: abs((x[1]+x[2])/2 - jump_center))
            return f"{closest[0]} (pH {closest[1]}–{closest[2]}, {closest[3]}) [approximate]"

    # 复用 TitrationCurve 的 pH 计算方法
    def _ph_sa_sb(self, Ca, Va, Ct, Vb, Vt): return self._ph_sa_sb_static(Ca, Va, Ct, Vb, Vt)
    def _ph_wa_sb(self, Ca, Va, Ct, Vb, Vt, Ka, Ve): return self._ph_wa_sb_static(Ca, Va, Ct, Vb, Vt, Ka, Ve)
    def _ph_wb_sa(self, Ca, Va, Ct, Vb, Vt, Kb, Ve): return self._ph_wb_sb_static(Ca, Va, Ct, Vb, Vt, Kb, Ve)

    @staticmethod
    def _ph_sa_sb_static(Ca, Va, Ct, Vb, Vt):
        mol_H, mol_OH = Ca*Va, Ct*Vb
        ex = mol_H - mol_OH
        if abs(ex) < 1e-12: return 7.0, "neutral"
        elif ex > 0: return -math.log10(ex/Vt), "acidic"
        else: return 14 + math.log10(abs(ex)/Vt), "basic"

    @staticmethod
    def _ph_wa_sb_static(Ca, Va, Ct, Vb, Vt, Ka, Ve):
        Kw = 1e-14; mol_HA, mol_OH = Ca*Va, Ct*Vb
        if mol_OH <= 0:
            C_HA = mol_HA/Vt; d = Ka**2+4*Ka*C_HA; h = (-Ka+math.sqrt(d))/2
            return -math.log10(h) if h>0 else 7.0, "pure HA"
        elif mol_OH < mol_HA:
            c_ha = (mol_HA-mol_OH)/Vt; c_a = mol_OH/Vt
            pKa = -math.log10(Ka); ph = pKa + math.log10(c_a/c_ha) if c_a>0 and c_ha>0 else 7.0
            return ph, "buffer"
        else:
            oh = (mol_OH-mol_HA)/Vt; poh = -math.log10(oh) if oh>0 else 7.0
            return 14-poh, "excess base"

    @staticmethod
    def _ph_wb_sb_static(Ca, Va, Ct, Vb, Vt, Kb, Ve):
        Kw = 1e-14; mol_B, mol_H = Ca*Va, Ct*Vb
        if mol_H <= 0:
            C_B = mol_B/Vt; d = Kb**2+4*Kb*C_B; oh = (-Kb+math.sqrt(d))/2
            return 14-(-math.log10(oh) if oh>0 else 7.0), "pure B"
        elif mol_H < mol_B:
            c_b = (mol_B-mol_H)/Vt; c_bh = mol_H/Vt
            pKb = -math.log10(Kb); poh = pKb + math.log10(c_b/c_bh) if c_b>0 and c_bh>0 else 7.0
            return 14-poh, "buffer"
        else:
            h = (mol_H-mol_B)/Vt
            return -math.log10(h) if h>0 else 7.0, "excess acid"

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) < 5:
                raise ValueError("Need atype, conc, vol_ml, ttype, tconc.")
            atype = parts[0]; conc = float(parts[1]); vol = float(parts[2])
            ttype = parts[3]; tconc = float(parts[4])
            ka = float(parts[5]) if len(parts) > 5 else None
            kb = float(parts[6]) if len(parts) > 6 else None
            return self._run_base(atype, conc, vol, ttype, tconc, ka, kb)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
