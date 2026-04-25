import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TitrationCurve(BaseTool):
    """
    生成滴定曲线数据点。
    支持强酸-强碱、弱酸-强碱、强酸-弱碱三种滴定类型。
    """
    __version__ = "0.1.0"
    name = "TitrationCurve"
    func_name = "generate_titration_curve"
    description = "Generate titration curve data points for acid-base titrations (strong/weak combinations)."
    implementation_description = "Calculates pH at each point of a titration curve using equilibrium equations. Supports strong acid-strong base, weak acid-strong base, and strong acid-weak base."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Titration", "pH Curve", "Acid-Base", "Equivalence Point"]
    required_envs = []

    code_input_sig = [
        ("analyte_type", "str", "N/A", "Analyte type: 'strong_acid', 'weak_acid', 'strong_base', 'weak_base'."),
        ("analyte_conc", "float", "N/A", "Analyte concentration (mol/L)."),
        ("analyte_volume_ml", "float", "N/A", "Initial analyte volume (mL)."),
        ("titrant_type", "str", "N/A", "Titrant type: 'strong_base' (for acids) or 'strong_acid' (for bases)."),
        ("titrant_conc", "float", "N/A", "Titrant concentration (mol/L)."),
        ("num_points", "int", "51", "Number of data points to generate along the curve."),
        ("Ka", "float", "None", "Ka of weak acid (needed if analyte is weak_acid)."),
        ("Kb", "float", "None", "Kb of weak base (needed if analyte is weak_base)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'analyte_type conc vol_ml titrant_type titrant_conc [num_points] [Ka/Kb]'. Example: 'strong_acid 0.1 50 strong_base 0.1 21'"),
    ]

    output_sig = [
        ("curve_points", "list", "List of dicts with volume_added_ml, ph, dominant_species at each point."),
        ("equivalence_volume_ml", "float", "Volume of titrant at equivalence point (mL)."),
        ("titration_type", "str", "Description of the titration system."),
        ("analyte_info", "dict", "Summary of analyte parameters."),
        ("titrant_info", "dict", "Summary of titrant parameters."),
    ]

    examples = [
        {
            "code_input": {
                "analyte_type": "strong_acid",
                "analyte_conc": 0.1,
                "analyte_volume_ml": 50,
                "titrant_type": "strong_base",
                "titrant_conc": 0.1,
                "num_points": 11,
                "Ka": None,
                "Kb": None,
            },
            "text_input": {
                "input_params": "strong_acid 0.1 50 strong_base 0.1 11"
            },
            "output": {
                "equivalence_volume_ml": 50.0,
                "titration_type": "Strong Acid + Strong Base",
                "curve_points": [
                    {"volume_added_ml": 0, "ph": 1.0, "dominant_species": "HCl"},
                    {"volume_added_ml": 25, "ph": 1.48, "dominant_species": "HCl + NaCl"},
                    {"volume_added_ml": 50, "ph": 7.0, "dominant_species": "NaCl"},
                    {"volume_added_ml": 75, "ph": 12.3, "dominant_species": "NaOH + NaCl"},
                ],
                "analyte_info": {"type": "strong_acid", "conc_mol_L": 0.1, "volume_ml": 50},
                "titrant_info": {"type": "strong_base", "conc_mol_L": 0.1},
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Kw = 1.0e-14

    def _run_base(self, analyte_type: str, analyte_conc: float, analyte_volume_ml: float,
                  titrant_type: str, titrant_conc: float, num_points: int = 51,
                  Ka: float = None, Kb: float = None) -> dict:
        """核心逻辑：生成滴定曲线"""
        # 参数校验
        if analyte_conc <= 0 or titrant_conc <= 0:
            raise ChemMCPError("Concentrations must be positive.")
        if analyte_volume_ml <= 0:
            raise ChemMCPError("Volume must be positive.")
        num_points = max(3, min(500, num_points))

        Va = analyte_volume_ml / 1000.0  # L
        Ca = analyte_conc
        Ct = titrant_conc

        # 等当点体积
        V_eq = (Ca * Va) / Ct * 1000  # mL

        # 确定滴定类型
        atype = analyte_type.lower().strip()
        ttype = titrant_type.lower().strip()

        # 生成体积点（从 0 到 V_eq*1.2）
        V_max = V_eq * 1.2
        curve_points = []

        for i in range(num_points):
            Vb = V_max * i / (num_points - 1) if num_points > 1 else V_eq / 2
            Vb_L = Vb / 1000.0

            ph, species = self._calc_ph_at_point(
                atype, ttype, Ca, Va, Ct, Vb_L, Ka, Kb, V_eq
            )
            curve_points.append({
                "volume_added_ml": round(Vb, 4),
                "ph": round(ph, 4),
                "dominant_species": species,
            })

        tit_type_str = f"{atype.replace('_', ' ').title()} + {ttype.replace('_', ' ').title()}"

        logger.info(f"Generated titration curve: {num_points} points, V_eq={V_eq:.2f} mL")

        return {
            "curve_points": curve_points,
            "equivalence_volume_ml": round(V_eq, 4),
            "titration_type": tit_type_str,
            "analyte_info": {"type": atype, "conc_mol_L": Ca, "volume_ml": analyte_volume_ml},
            "titrant_info": {"type": ttype, "conc_mol_L": Ct},
        }

    def _calc_ph_at_point(self, atype: str, ttype: str, Ca: float, Va: float,
                           Ct: float, Vb: float, Ka: float, Kb: float, V_eq: float) -> tuple:
        """计算某一体积点的 pH"""
        V_total = Va + Vb
        if V_total == 0:
            V_total = 1e-15

        if atype in ("strong_acid", "sa") and ttype in ("strong_base", "sb"):
            return self._ph_sa_sb(Ca, Va, Ct, Vb, V_total)
        elif atype in ("weak_acid", "wa") and ttype in ("strong_base", "sb"):
            if Ka is None:
                raise ChemMCPError("Ka required for weak_acid titration.")
            return self._ph_wa_sb(Ca, Va, Ct, Vb, V_total, Ka, V_eq)
        elif atype in ("strong_base", "sb") and ttype in ("strong_acid", "sa"):
            return self._ph_sb_sa(Ca, Va, Ct, Vb, V_total)
        elif atype in ("weak_base", "wb") and ttype in ("strong_acid", "sa"):
            if Kb is None:
                raise ChemMCPError("Kb required for weak_base titration.")
            return self._ph_wb_sa(Ca, Va, Ct, Vb, V_total, Kb, V_eq)
        else:
            raise ChemMCPError(f"Unsupported combination: {atype} + {ttype}")

    def _ph_sa_sb(self, Ca: float, Va: float, Ct: float, Vb: float, V_total: float) -> tuple:
        """强酸 + 强碱"""
        mol_H = Ca * Va
        mol_OH = Ct * Vb
        excess_H = mol_H - mol_OH

        if abs(excess_H) < 1e-12:
            ph = 7.0
            species = "Neutral (salt + water)"
        elif excess_H > 0:
            h = excess_H / V_total
            ph = -math.log10(h) if h > 0 else 7.0
            species = f"Excess H+ ({h:.3e} M)"
        else:
            oh = abs(excess_H) / V_total
            poh = -math.log10(oh) if oh > 0 else 7.0
            ph = 14 - poh
            species = f"Excess OH- ({oh:.3e} M)"

        return ph, species

    def _ph_wa_sb(self, Ca: float, Va: float, Ct: float, Vb: float,
                  V_total: float, Ka: float, V_eq: float) -> tuple:
        """弱酸 + 强碱"""
        mol_HA_initial = Ca * Va
        mol_OH_added = Ct * Vb
        V_eq_L = V_eq / 1000.0

        if mol_OH_added <= 0:
            # 纯弱酸
            C_HA = mol_HA_initial / V_total
            disc = Ka ** 2 + 4 * Ka * C_HA
            h = (-Ka + math.sqrt(disc)) / 2
            ph = -math.log10(h) if h > 0 else 7.0
            species = "HA (mostly undissociated)"
        elif mol_OH_added < mol_HA_initial:
            # 缓冲区：HA + A-
            mol_HA_left = mol_HA_initial - mol_OH_added
            mol_A_formed = mol_OH_added
            c_ha = mol_HA_left / V_total
            c_a = mol_A_formed / V_total
            pKa = -math.log10(Ka)
            if c_a > 0 and c_ha > 0:
                ph = pKa + math.log10(c_a / c_ha)
            elif c_a == 0:
                disc = Ka ** 2 + 4 * Ka * c_ha
                h = (-Ka + math.sqrt(disc)) / 2
                ph = -math.log10(h) if h > 0 else 7.0
            else:
                ph = 14.0
            species = f"Buffer region (HA + A-)"
        elif abs(mol_OH_added - mol_HA_initial) < 1e-12:
            # 等当点：纯共轭碱 A-
            c_A = mol_HA_initial / V_total
            Kh = self.Kw / Ka
            disc = Kh ** 2 + 4 * Kh * c_A
            oh = (-Kh + math.sqrt(disc)) / 2
            poh = -math.log10(oh) if oh > 0 else 7.0
            ph = 14 - poh
            species = "Equivalence point (A- hydrolysis)"
        else:
            # 过量强碱
            mol_OH_excess = mol_OH_added - mol_HA_initial
            oh = mol_OH_excess / V_total
            poh = -math.log10(oh) if oh > 0 else 7.0
            ph = 14 - poh
            species = f"Excess OH-"

        return ph, species

    def _ph_sb_sa(self, Ca: float, Va: float, Ct: float, Vb: float, V_total: float) -> tuple:
        """强碱 + 强酸（与 sa+sb 对称）"""
        mol_OH = Ca * Va
        mol_H = Ct * Vb
        excess_OH = mol_OH - mol_H

        if abs(excess_OH) < 1e-12:
            ph = 7.0
            species = "Neutral (salt + water)"
        elif excess_OH > 0:
            oh = excess_OH / V_total
            poh = -math.log10(oh) if oh > 0 else 7.0
            ph = 14 - poh
            species = f"Excess OH-"
        else:
            h = abs(excess_OH) / V_total
            ph = -math.log10(h) if h > 0 else 7.0
            species = f"Excess H+ ({h:.3e} M)"

        return ph, species

    def _ph_wb_sa(self, Ca: float, Va: float, Ct: float, Vb: float,
                  V_total: float, Kb: float, V_eq: float) -> tuple:
        """弱碱 + 强酸"""
        mol_B_initial = Ca * Va
        mol_H_added = Ct * Vb

        if mol_H_added <= 0:
            # 纯弱碱
            C_B = mol_B_initial / V_total
            disc = Kb ** 2 + 4 * Kb * C_B
            oh = (-Kb + math.sqrt(disc)) / 2
            poh = -math.log10(oh) if oh > 0 else 7.0
            ph = 14 - poh
            species = "B (mostly unprotonated)"
        elif mol_H_added < mol_B_initial:
            # 缓冲区：B + BH+
            mol_B_left = mol_B_initial - mol_H_added
            mol_BH_formed = mol_H_added
            c_b = mol_B_left / V_total
            c_bh = mol_BH_formed / V_total
            pKb = -math.log10(Kb)
            pKa_conj = 14 - pKb
            if c_b > 0 and c_bh > 0:
                poh = pKb + math.log10(c_b / c_bh)
                ph = 14 - poh
            else:
                ph = 1.0
            species = f"Buffer region (B + BH+)"
        elif abs(mol_H_added - mol_B_initial) < 1e-12:
            # 等当点：共轭酸 BH+
            c_BH = mol_B_initial / V_total
            Ka_conj = self.Kw / Kb
            disc = Ka_conj ** 2 + 4 * Ka_conj * c_BH
            h = (-Ka_conj + math.sqrt(disc)) / 2
            ph = -math.log10(h) if h > 0 else 7.0
            species = "Equivalence point (BH+ dissociation)"
        else:
            mol_H_excess = mol_H_added - mol_B_initial
            h = mol_H_excess / V_total
            ph = -math.log10(h) if h > 0 else 7.0
            species = f"Excess H+"

        return ph, species

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) < 5:
                raise ValueError("Need analyte_type, conc, vol_ml, titrant_type, titrant_conc.")
            atype = parts[0]
            conc = float(parts[1])
            vol = float(parts[2])
            ttype = parts[3]
            tconc = float(parts[4])
            npts = int(parts[5]) if len(parts) > 5 else 51
            ka_val = float(parts[6]) if len(parts) > 6 else None
            kb_val = float(parts[7]) if len(parts) > 7 else None
            return self._run_base(atype, conc, vol, ttype, tconc, npts, ka_val, kb_val)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
