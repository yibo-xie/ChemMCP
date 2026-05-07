import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PrecipitationTitrationCalculator(BaseTool):
    """
    沉淀滴定计算工具：Mohr法、Volhard法、Fajans法。
    
    支持AgNO₃滴定卤素离子，计算当量点、pAg/pX曲线和指示剂变色点。
    """
    __version__ = "0.1.0"
    name             = "PrecipitationTitrationCalculator"
    func_name        = "precipitation_titration"
    description      = "Calculate precipitation titration parameters for Mohr, Volhard, and Fajans methods with pX/pAg curves."
    implementation_description = "Solves solubility equilibrium for Ag⁺ titrating halides. Computes pAg/pX at each point using Ksp values. Supports Mohr (CrO₄²⁻ indicator), Volhard (back-titration with SCN⁻/Fe³⁺), and Fajans (adsorption indicator) methods."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Precipitation Titration", "Mohr Method", "Volhard Method", "Analytical Chemistry"]
    required_envs    = []

    # 常见沉淀的 Ksp 数据（25°C）
    KSP_DATA = {
        "Cl": {"ion": "Cl⁻", "formula": "AgCl",   "Ksp": 1.8e-10,  "color_white": True},
        "Br": {"ion": "Br⁻", "formula": "AgBr",   "Ksp": 5.0e-13,  "color_white": False, "pale_yellow": True},
        "I":  {"ion": "I⁻",  "formula": "AgI",    "Ksp": 8.3e-17,  "color_yellow": True},
        "SCN":{"ion": "SCN⁻","formula": "AgSCN",  "Ksp": 1.0e-12,  "color_white": True},
        # Mohr 法指示剂相关
        "Ag2CrO4": {"formula": "Ag₂CrO₄", "Ksp": 1.12e-12, "indicator": True},
    }

    # 沉淀滴定方法配置
    METHOD_CONFIGS = {
        "mohr": {
            "name": "Mohr Method",
            "description": "Direct titration of Cl⁻/Br⁻ with AgNO₃ using K₂CrO₄ indicator",
            "applicable_analytes": ["Cl", "Br"],
            "indicator": "K₂CrO₄ (chromate)",
            "indicator_reaction": "2Ag⁺ + CrO₄²⁻ → Ag₂CrO₄↓(brick red)",
            "ph_range": (6.5, 10.5),
            "note": "Must be neutral or weakly alkaline; acidic conditions dissolve CrO₄²⁻ precipitate.",
        },
        "volhard": {
            "name": "Volhard Method",
            "description": "Back-titration: excess AgNO₃ + NH₄SCN titrated with Fe³⁺ indicator",
            "applicable_analytes": ["Cl", "Br", "I", "SCN"],
            "indicator": "Fe³⁺ (ferric alum / iron(III) ammonium sulfate)",
            "indicator_reaction": "Fe³⁺ + SCN⁻ → [Fe(SCN)]²⁺ (red)",
            "ph_range": (0, 1),
            "note": "Acidic medium (HNO₃) to prevent Fe³⁺ hydrolysis and Ag⁺ interference.",
        },
        "fajans": {
            "name": "Fajans Method",
            "description": "Direct titration with adsorption indicator (fluorescein, dichlorofluorescein)",
            "applicable_analytes": ["Cl", "Br", "I"],
            "indicator": "Fluorescein / Dichlorofluorescein",
            "indicator_reaction": "Adsorption of indicator on precipitate surface causes color change",
            "ph_range": (7, 10),
            "note": "pH depends on specific indicator used.",
        },
    }

    code_input_sig   = [
        ("method_type", "str", "N/A", "Method type: 'mohr', 'volhard', or 'fajans'."),
        ("analyte_key", "str", "N/A", "Analyte ion key from KSP data (e.g., 'Cl', 'Br', 'I')."),
        ("analyte_concentration_mol_L", "float", "N/A", "Initial concentration of analyte (mol/L)."),
        ("analyte_volume_ml", "float", "N/A", "Volume of analyte solution (mL)."),
        ("titrant_concentration_mol_L", "float", "N/A", "Concentration of AgNO₃ titrant (mol/L)."),
        ("n_points", "int", "100", "Number of points on the titration curve."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'method analyte C_analyte V_analyte C_titrant [n_points]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with equivalence point, pAg/pX curve, method-specific details, and feasibility assessment."),
    ]

    examples         = [
        {
            "code_input": {
                "method_type": "mohr",
                "analyte_key": "Cl",
                "analyte_concentration_mol_L": 0.050,
                "analyte_volume_ml": 25.00,
                "titrant_concentration_mol_L": 0.050,
                "n_points": 50,
            },
            "text_input": {
                "input_params": "mohr Cl 0.050 25.00 0.050 50",
            },
            "output": {
                "result": {
                    "equivalence_point": {"volume_ml": "...", "pAg": "...", "pX": "..."},
                    "curve_data": [{"volume_ml": ..., "pAg": ..., "pX": ...}, ...],
                    "method_details": {...}
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
        method_type: str,
        analyte_key: str,
        analyte_concentration_mol_L: float,
        analyte_volume_ml: float,
        titrant_concentration_mol_L: float,
        n_points: int = 100,
    ) -> dict:
        """核心逻辑：沉淀滴定计算"""
        # 验证输入
        mtype = method_type.lower()
        if mtype not in self.METHOD_CONFIGS:
            raise ChemMCPError(f"Unknown method '{method_type}'. Available: {list(self.METHOD_CONFIGS.keys())}")

        # Case-insensitive lookup
        akey = analyte_key
        if akey not in self.KSP_DATA or akey == "Ag2CrO4":
            # Try capitalized
            akey = analyte_key.capitalize()
        if akey not in self.KSP_DATA or akey == "Ag2CrO4":
            akey = analyte_key.upper()
        if akey not in self.KSP_DATA or akey == "AG2CRO4":
            available = [k for k in self.KSP_DATA.keys() if k != "AG2CRO4"]
            raise ChemMCPError(f"Unknown analyte '{analyte_key}'. Available: {available}")

        method = self.METHOD_CONFIGS[mtype]
        analyte = self.KSP_DATA[akey]

        Ca = float(analyte_concentration_mol_L)
        Va = float(analyte_volume_ml)
        Ct = float(titrant_concentration_mol_L)
        Ksp = analyte["Ksp"]

        # 当量点体积
        Ve = Ca * Va / Ct if Ct > 0 else 0

        # ---- 滴定曲线：pAg 和 pX ----
        curve_data = []
        V_max = Ve * 1.4

        for i in range(n_points + 1):
            Vag = V_max * i / n_points
            Vtotal = Va + Vag
            if Vtotal < 1e-12:
                Vtotal = 1e-10

            if Vag < 1e-12:
                # 滴定前：纯分析物溶液
                X_conc = Ca  # 近似（忽略解离）
                pX = -math.log10(X_conc) if X_conc > 1e-14 else 14
                Ag_conc = Ksp / X_conc if X_conc > 1e-20 else math.sqrt(Ksp)
                pAg = -math.log10(Ag_conc) if Ag_conc > 1e-14 else 14
            elif Vag < Ve - 1e-10:
                # 当量点前：过量分析物
                mol_X_initial = Ca * Va
                mol_Ag_added = Ct * Vag
                mol_X_remaining = mol_X_initial - mol_Ag_added
                X_free = mol_X_remaining / Vtotal
                Ag_free = Ksp / X_free if X_free > 1e-30 else 0
                pX = -math.log10(X_free) if X_free > 1e-14 else 14
                pAg = -math.log10(Ag_free) if Ag_free > 1e-14 else 14
            elif abs(Vag - Ve) < 1e-10:
                # 当量点：[Ag⁺] = [X⁻] = √Ksp
                sqrt_Ksp = math.sqrt(Ksp)
                pAg_eq = -math.log10(sqrt_Ksp) if sqrt_Ksp > 1e-14 else 14
                pX_eq = pAg_eq
                pAg = pAg_eq
                pX = pX_eq
            else:
                # 过量 AgNO₃
                mol_Ag_excess = (Vag - Ve) * Ct
                Ag_free = mol_Ag_excess / Vtotal
                X_free = Ksp / Ag_free if Ag_free > 1e-30 else 0
                pAg = -math.log10(Ag_free) if Ag_free > 1e-14 else 14
                pX = -math.log10(X_free) if X_free > 1e-14 else 14

            curve_data.append({
                "volume_agno3_ml": round(Vag, 4),
                "pAg": round(pAg, 4),
                "pX": round(pX, 4),
            })

        # 当量点数据
        sqrt_Ksp = math.sqrt(Ksp)
        pAg_eq = round(-math.log10(sqrt_Ksp), 4) if sqrt_Ksp > 1e-14 else 14
        pX_eq = pAg_eq

        # ---- 方法特定信息 ----
        method_details = None
        if mtype == "mohr":
            method_details = self._mohr_details(analyte, pAg_eq)
        elif mtype == "volhard":
            method_details = self._volhard_details()
        elif mtype == "fajans":
            method_details = self._fajans_details()

        # 突跃范围估算（±0.1%当量点的pAg变化）
        # 在99%和101%当量点处
        jump_99 = None
        jump_101 = None
        for pt in curve_data:
            if abs(pt["volume_agno3_ml"] - Ve * 0.99) < V_max / n_points * 0.6:
                jump_99 = pt["pAg"]
            if abs(pt["volume_agno3_ml"] - Ve * 1.01) < V_max / n_points * 0.6:
                jump_101 = pt["pAg"]

        result = {
            "method": method["name"],
            "analyte": analyte["ion"],
            "precipitate_formula": analyte["formula"],
            "Ksp": f"{Ksp:.2e}",
            "pKsp": round(-math.log10(Ksp), 4),
            "conditions": {
                "C_analyte_mol_L": Ca,
                "V_analyte_ml": Va,
                "C_titrant_mol_L": Ct,
            },
            "equivalence_point": {
                "volume_ml": round(Ve, 6),
                "pAg": pAg_eq,
                "pX": pX_eq,
                "[Ag+]_eq_mol_L": f"{sqrt_Ksp:.2e}",
                "[X-]_eq_mol_L": f"{sqrt_Ksp:.2e}",
            },
            "jump_range": {
                "pAg_at_99%": jump_99,
                "pAg_at_101%": jump_101,
                "delta_pAg": round((jump_101 or 0) - (jump_99 or 0), 2) if jump_99 and jump_101 else None,
            } if jump_99 and jump_101 else None,
            "method_specific": method_details,
            "titration_curve": curve_data,
        }

        logger.info(f"Precipitation titration ({method['name']}): {analyte['ion']}, Ve={Ve:.2f}mL, pAg_eq={pAg_eq}")
        return result

    def _mohr_details(self, analyte: dict, pAg_eq: float) -> dict:
        """Mohr法详细信息"""
        Ksp_Ag2CrO4 = self.KSP_DATA["Ag2CrO4"]["Ksp"]
        # CrO₄²⁻开始沉淀的条件: [Ag⁺]² × [CrO₄²⁻] ≥ Ksp(Ag₂CrO₄)
        # 设 [CrO₄²⁻] ≈ 0.005~0.01 mol/L (常用浓度)
        c_chromate_typical = 0.005
        Ag_for_indicator = math.sqrt(Ksp_Ag2CrO4 / c_chromate_typical)
        pAg_indicator = round(-math.log10(Ag_for_indicator), 4)

        return {
            "indicator": "Potassium chromate (K₂CrO₄)",
            "typical_concentration_M": "0.005 ~ 0.01",
            "indicator_precipitation_condition": f"[Ag⁺] ≥ √(Ksp_Ag2CrO4/[CrO4²-]) = {Ag_for_indicator:.2e} M (pAg ≤ {pAg_indicator})",
            "pAg_at_indicator_color_change": pAg_indicator,
            "pAg_at_equivalence": pAg_eq,
            "indicator_error_possible": pAg_indicator > pAg_eq,
            "note": "Indicator should change color slightly after equivalence point for accurate results.",
            "ph_requirement": f"pH {self.METHOD_CONFIGS['mohr']['ph_range'][0]}-{self.METHOD_CONFIGS['mohr']['ph_range'][1]}",
        }

    def _volhard_details(self) -> dict:
        """Volhard法详细信息"""
        return {
            "type": "Back-titration",
            "procedure": (
                "1. Add known excess AgNO₃ to sample containing analyte\n"
                "2. Filter or add nitrobenzene (for Cl⁻)\n"
                "3. Titrate excess Ag⁺ with NH₄SCN using Fe³⁺ as indicator\n"
                "4. Red color of [Fe(SCN)]²⁺ indicates endpoint"
            ),
            "indicator": "Iron(III) ammonium sulfate / Ferric alum (Fe³⁺)",
            "indicator_concentration": "~0.015 M Fe³⁺",
            "medium": "HNO₃ (nitric acid), pH 0-1",
            "note_for_Cl": "For Cl⁻ analysis, must filter AgCl precipitate or add nitrobenzene to prevent AgCl reacting with SCN⁻.",
        }

    def _fajans_details(self) -> dict:
        """Fajans法详细信息"""
        return {
            "type": "Adsorption indicator method",
            "principle": "Indicator dye is adsorbed onto precipitate surface at endpoint, causing color change.",
            "common_indicators": [
                {"name": "Fluorescein", "ph_range": "7-10", "analytes": "Cl⁻, Br⁻, I⁻", "color_change": "Yellow-green → Pink"},
                {"name": "Dichlorofluorescein", "ph_range": "4-10", "analytes": "Cl⁻, Br⁻, I⁻", "color_change": "Yellow-green → Red"},
                {"name": "Eosin", "ph_range": "2-10", "analytes": "Br⁻, I⁻, SCN⁻", "color_change": "Orange → Red-violet"},
            ],
            "note": "Avoid strong illumination (some indicators are light-sensitive); colloidal precipitate needed for good endpoint detection.",
        }

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            mtype = parts[0]
            akey = parts[1]
            Ca = float(parts[2])
            Va = float(parts[3])
            Ct = float(parts[4])
            np_ = int(parts[5]) if len(parts) > 5 else 100
            return self._run_base(mtype, akey, Ca, Va, Ct, np_)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'method analyte C_analyte V_analyte C_titrant [n_points]'")
