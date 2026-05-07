import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ComplexometricTitrationHelper(BaseTool):
    """
    配位滴定计算工具：EDTA滴定、条件稳定常数、pM计算。
    
    支持常见金属离子的EDTA滴定，包括条件稳定常数K'f的计算和滴定曲线生成。
    """
    __version__ = "0.1.0"
    name             = "ComplexometricTitrationHelper"
    func_name        = "complexometric_titration"
    description      = "Calculate EDTA complexometric titration parameters: conditional formation constant, pM at equivalence, minimum pH requirement, and titration curve."
    implementation_description = "Computes conditional stability constant K'_f = α_Y(H⁻) × K_f for EDTA-metal complexes at given pH, determines minimum pH for accurate titration (ΔpM ≥ 0.2, TE ≤ 0.1%), and generates pM vs volume curve data."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Complexometry", "EDTA", "Titration", "Conditional Constant", "Analytical Chemistry"]
    required_envs    = []

    # 常见金属离子与EDTA配合物的 lgKf 数据（25°C, I=0.1）
    METAL_DATA = {
        "Mg": {"ion": "Mg²⁺", "charge": 2, "lg_kf": 8.7,  "color": "colorless"},
        "Ca": {"ion": "Ca²⁺", "charge": 2, "lg_kf": 10.69, "color": "colorless"},
        "Ba": {"ion": "Ba²⁺", "charge": 2, "lg_kf": 7.86, "color": "colorless"},
        "Zn": {"ion": "Zn²⁺", "charge": 2, "lg_kf": 16.50, "color": "colorless"},
        "Cu": {"ion": "Cu²⁺", "charge": 2, "lg_kf": 18.80, "color": "blue"},
        "Pb": {"ion": "Pb²⁺", "charge": 2, "lg_kf": 18.04, "color": "colorless"},
        "Mn": {"ion": "Mn²⁺", "charge": 2, "lg_kf": 13.87, "color": "pale_pink"},
        "Fe2": {"ion": "Fe²⁺", "charge": 2, "lg_kf": 14.32, "color": "pale_green"},
        "Fe3": {"ion": "Fe³⁺", "charge": 3, "lg_kf": 25.1,  "color": "yellow"},
        "Ni": {"ion": "Ni²⁺", "charge": 2, "lg_kf": 18.62, "color": "green"},
        "Co": {"ion": "Co²⁺", "charge": 2, "lg_kf": 16.31, "color": "pink"},
        "Cd": {"ion": "Cd²⁺", "charge": 2, "lg_kf": 16.46, "color": "colorless"},
        "Hg": {"ion": "Hg²⁺", "charge": 2, "lg_kf": 21.80, "color": "colorless"},
        "Al": {"ion": "Al³⁺", "charge": 3, "lg_kf": 16.3,  "color": "colorless"},
        "Cr": {"ion": "Cr³⁺", "charge": 3, "lg_kf": 23.0,  "color": "green_violet"},
        "Bi": {"ion": "Bi³⁺", "charge": 3, "lg_kf": 27.94, "color": "colorless"},
        "Sr": {"ion": "Sr²⁺", "charge": 2, "lg_kf": 8.63,  "color": "colorless"},
        "Ag": {"ion": "Ag⁺",  "charge": 1, "lg_kf": 7.32,  "color": "colorless"},
    }

    # EDTA 的 lgα_Y(H⁻) 与 pH 关系数据表（常用 pH 范围）
    # 来源：分析化学教材标准数据
    _ALPHA_Y_H_TABLE = [
        (0.0, 23.64), (0.5, 23.37), (1.0, 22.96), (1.5, 22.39), (2.0, 21.72),
        (2.5, 20.90), (3.0, 19.91), (3.5, 18.74), (4.0, 17.44), (4.5, 16.02),
        (5.0, 14.45), (5.5, 12.73), (6.0, 10.95), (6.5, 9.19),  (7.0, 7.50),
        (7.5, 6.00),  (8.0, 4.65),  (8.5, 3.55),  (9.0, 2.70),  (9.5, 2.00),
        (10.0, 1.47), (10.5, 1.09), (11.0, 0.81), (11.5, 0.62), (12.0, 0.46),
        (12.5, 0.35), (13.0, 0.26), (13.5, 0.20), (14.0, 15),
    ]

    code_input_sig   = [
        ("metal_symbol", "str", "N/A", "Metal ion symbol from database (e.g., 'Mg', 'Ca', 'Zn', 'Cu', 'Fe3')."),
        ("metal_concentration_mol_L", "float", "N/A", "Concentration of metal ion solution (mol/L)."),
        ("metal_volume_ml", "float", "N/A", "Volume of metal ion solution (mL)."),
        ("edta_concentration_mol_L", "float", "N/A", "Concentration of EDTA titrant (mol/L)."),
        ("ph", "float", "N/A", "pH of the solution during titration."),
        ("auxiliary_conc_or_None", "float_or_None", "None", "Concentration of auxiliary complexing agent (mol/L) if used."),
        ("auxiliary_lgf_or_None", "float_or_None", "None", "log(formation constant) of metal-auxiliary ligand complex if used."),
        ("n_points", "int", "100", "Number of points on the pM titration curve."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'metal C_metal V_mol C_edta pH [aux_conc] [aux_lgf] [n_points]'"),
    ]

    output_sig       = [
        ("result", "dict", "Dictionary with conditional K'f, pM at eq, min pH, titration curve, and feasibility assessment."),
    ]

    examples         = [
        {
            "code_input": {
                "metal_symbol": "Zn",
                "metal_concentration_mol_L": 0.010,
                "metal_volume_ml": 20.00,
                "edta_concentration_mol_L": 0.010,
                "ph": 10.0,
                "auxiliary_conc_or_None": None,
                "auxiliary_lgf_or_None": None,
                "n_points": 50,
            },
            "text_input": {
                "input_params": "Zn 0.010 20.00 0.010 10.0 50",
            },
            "output": {
                "result": {
                    "conditional_log_kf_prime": "...",
                    "pm_at_equivalence": "...",
                    "minimum_ph": "...",
                    "feasible": True,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _get_alpha_y_h(self, ph: float) -> float:
        """插值获取 lgα_Y(H⁻)"""
        table = self._ALPHA_Y_H_TABLE
        if ph <= table[0][0]:
            return table[0][1]
        if ph >= table[-1][0]:
            return table[-1][1]
        for i in range(len(table) - 1):
            if table[i][0] <= ph <= table[i+1][0]:
                t = (ph - table[i][0]) / (table[i+1][0] - table[i][0])
                return table[i][1] + t * (table[i+1][1] - table[i][1])
        return table[-1][1]

    def _run_base(
        self,
        metal_symbol: str,
        metal_concentration_mol_L: float,
        metal_volume_ml: float,
        edta_concentration_mol_L: float,
        ph: float,
        auxiliary_conc_or_None: Optional[float] = None,
        auxiliary_lgf_or_None: Optional[float] = None,
        n_points: int = 100,
    ) -> dict:
        """核心逻辑：配位滴定计算"""
        # Case-insensitive lookup: try exact, then capitalized
        metal_key = metal_symbol
        if metal_key not in self.METAL_DATA:
            # Try capitalizing first letter: "zn" -> "Zn"
            metal_key = metal_symbol.capitalize()
        if metal_key not in self.METAL_DATA:
            # Try upper as last resort
            metal_key = metal_symbol.upper()
        if metal_key not in self.METAL_DATA:
            available = ", ".join(sorted(self.METAL_DATA.keys()))
            raise ChemMCPError(f"Unknown metal '{metal_symbol}'. Available: {available}")

        metal = self.METAL_DATA[metal_key]
        Cm = float(metal_concentration_mol_L)
        Vm = float(metal_volume_ml)
        Ce = float(edta_concentration_mol_L)

        lg_Kf = metal["lg_kf"]
        Kf = 10 ** lg_Kf

        # ---- 条件稳定常数 ----
        lg_alpha_y_h = self._get_alpha_y_h(ph)
        alpha_y_h = 10 ** lg_alpha_y_h

        # 辅助配位效应
        lg_alpha_L = 0.0
        if auxiliary_conc_or_None is not None and auxiliary_lgf_or_None is not None:
            # α_M(L) = 1 + β[L], 简化处理
            beta = 10 ** auxiliary_lgf_or_None
            alpha_L = 1 + beta * auxiliary_conc_or_None
            lg_alpha_L = math.log10(alpha_L) if alpha_L > 1 else 0

        lg_alpha_total = max(lg_alpha_y_h, lg_alpha_L)  # 近似
        lg_Kf_prime = lg_Kf - lg_alpha_y_h - lg_alpha_L
        Kf_prime = 10 ** lg_Kf_prime if lg_Kf_prime > -30 else 0

        # 当量点体积
        Ve = Cm * Vm / Ce if Ce > 0 else 0

        # 当量点 pM'
        Cm_eq = Cm * Vm / (Vm + Ve) if (Vm + Ve) > 0 else 0
        if Kf_prime > 0 and Cm_eq > 0:
            pm_eq = 0.5 * (lg_Kf_prime + math.log10(Cm_eq))
        else:
            pm_eq = 0

        # 最小pH要求（要求 lgK'f ≥ 8 才能准确滴定）
        min_ph_for_titration = None
        for ph_val, lg_alpha in self._ALPHA_Y_H_TABLE:
            lg_kf_test = lg_Kf - lg_alpha - lg_alpha_L
            if lg_kf_test >= 8:
                min_ph_for_titration = ph_val
                break

        # 滴定曲线
        curve_data = []
        V_max = Ve * 1.4
        for i in range(n_points + 1):
            Ve_added = V_max * i / n_points
            Vtotal = Vm + Ve_added
            if Vtotal == 0:
                Vtotal = 1e-10

            if Ve_added < 1e-12:
                # 滴定前：游离金属浓度 ≈ Cm（忽略解离）
                pM_val = -math.log10(Cm) if Cm > 1e-14 else 14
            elif Ve_added < Ve - 1e-10:
                # 当量点前：剩余金属为主
                fraction = Ve_added / Ve
                Cm_free = Cm * Vm * (1 - fraction) / Vtotal
                pM_val = -math.log10(Cm_free) if Cm_free > 1e-14 else 14
            elif abs(Ve_added - Ve) < 1e-10:
                # 当量点
                pM_val = pm_eq
            else:
                # 过量EDTA
                excess_edta = (Ve_added - Ve) * Ce / Vtotal
                # [MY] ≈ Cm*Vm/(Vm+Ve_added), [Y'] = excess_edta
                c_my = Cm * Vm / Vtotal
                if Kf_prime > 0:
                    m_free = c_my / (Kf_prime * excess_edta) if excess_edta > 1e-20 else math.sqrt(c_my / Kf_prime)
                else:
                    m_free = c_my
                pM_val = -math.log10(m_free) if m_free > 1e-14 else 14

            curve_data.append({
                "volume_edta_ml": round(Ve_added, 4),
                "pm": round(pM_val, 4),
            })

        # 可行性判断
        feasible = lg_Kf_prime >= 8
        sharpness = "sharp" if lg_Kf_prime >= 10 else "moderate" if lg_Kf_prime >= 6 else "poor"

        result = {
            "metal_ion": metal["ion"],
            "symbol": metal_key,
            "formation_constant": {
                "log_Kf": lg_Kf,
                "Kf": f"{Kf:.2e}",
            },
            "conditions": {
                "ph": round(ph, 4),
                "C_metal_mol_L": Cm,
                "V_metal_ml": Vm,
                "C_edta_mol_L": Ce,
            },
            "acid_effect": {
                "lg_alpha_Y_H": round(lg_alpha_y_h, 4),
                "alpha_Y_H": f"{alpha_y_h:.2e}",
            },
            "auxiliary_effect": {
                "lg_alpha_M_L": round(lg_alpha_L, 4) if lg_alpha_L else 0,
            } if auxiliary_conc_or_None else None,
            "conditional_constant": {
                "log_Kf_prime": round(lg_Kf_prime, 4),
                "Kf_prime": f"{Kf_prime:.2e}" if Kf_prime > 0 else "~0",
            },
            "equivalence_point": {
                "volume_edta_ml": round(Ve, 6),
                "C_metal_at_eq_mol_L": round(Cm_eq, 8),
                "pM_at_equivalence": round(pm_eq, 4),
            },
            "feasibility": {
                "feasible": feasible,
                "min_ph_required": round(min_ph_for_titration, 2) if min_ph_for_titration else "N/A (even at pH 14)",
                "current_ph_sufficient": ph >= (min_ph_for_titration or 99),
                "sharpness": sharpness,
                "recommendation": (
                    f"Titration is {'feasible' if feasible else 'not feasible'} at pH={ph}. "
                    f"lgK'f = {lg_Kf_prime:.2f} ({sharpness} endpoint)."
                ),
            },
            "titration_curve": curve_data,
        }

        logger.info(f"Complexometric titration: {metal['ion']}, pH={ph}, lgK'f={lg_Kf_prime:.2f}, feasible={feasible}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            metal = parts[0]
            Cm = float(parts[1])
            Vm = float(parts[2])
            Ce = float(parts[3])
            ph = float(parts[4])
            aux_c = float(parts[5]) if len(parts) > 5 and parts[5].lower() != "none" else None
            aux_f = float(parts[6]) if len(parts) > 6 and parts[6].lower() != "none" else None
            np_ = int(parts[7]) if len(parts) > 7 else 100
            return self._run_base(metal, Cm, Vm, Ce, ph, aux_c, aux_f, np_)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'metal C_metal V_mol C_edta pH [aux_conc] [aux_lgf] [n_points]'")
