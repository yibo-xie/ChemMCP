"""
沸点升高计算工具（依数性性质）
计算非挥发性溶质引起的溶液沸点升高：ΔTb = i × Kb × m
"""
import logging
from typing import Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BoilingPointElevation(BaseTool):
    """
    沸点升高计算工具。

    基于依数性公式 ΔTb = i × Kb × m 计算溶液的沸点升高值和新沸点。
    """
    __version__                 = "0.1.0"
    name                        = "BoilingPointElevation"
    func_name                   = "calculate_boiling_point_elevation"
    description                 = "Calculate boiling point elevation of a solution using the colligative property formula ΔTb = i × Kb × m."
    implementation_description  = "Uses the boiling point elevation formula for non-volatile solutes: ΔTb = i·Kb·m, where i is the van't Hoff factor, Kb is the ebullioscopic constant, and m is molality."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Colligative Properties", "Boiling Point", "Solution Chemistry", "Physical Chemistry"]
    required_envs               = []

    code_input_sig = [
        ("molality",              "float", "N/A",     "Molality of solution in mol solute / kg solvent."),
        ("ebullioscopic_constant_kb", "float", "0.512", "Ebullioscopic constant (Kb) in K·kg/mol (water = 0.512)."),
        ("vanthoff_factor_i",     "float", "1.0",     "van't Hoff factor (i=1 for nonelectrolyte, i≈2 for NaCl, etc.)."),
        ("normal_boiling_point_k","float", "373.15",   "Normal boiling point of pure solvent in Kelvin."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'molality [Kb] [i] [normal_bp_K]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing delta_tb_k, new_boiling_point_k, and explanation."),
    ]

    examples = [
        {
            "code_input": {
                "molality": 0.5,
                "ebullioscopic_constant_kb": 0.512,
                "vanthoff_factor_i": 1.0,
                "normal_boiling_point_k": 373.15,
            },
            "text_input": {
                "input_params": "0.5 0.512 1.0 373.15",
            },
            "output": {
                "result": {
                    "delta_tb_k": 0.256,
                    "new_boiling_point_k": 373.406,
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
        molality: float,
        ebullioscopic_constant_kb: float = 0.512,
        vanthoff_factor_i: float = 1.0,
        normal_boiling_point_k: float = 373.15,
    ) -> Dict[str, Any]:
        """核心逻辑：ΔTb = i × Kb × m"""
        if molality < 0:
            raise ChemMCPError("Molality cannot be negative.")
        if ebullioscopic_constant_kb < 0:
            raise ChemMCPError("Ebullioscopic constant cannot be negative.")
        if vanthoff_factor_i < 0:
            raise ChemMCPError("van't Hoff factor cannot be negative.")
        if normal_boiling_point_k <= 0:
            raise ChemMCPError("Normal boiling point must be positive.")

        delta_tb = vanthoff_factor_i * ebullioscopic_constant_kb * molality
        new_bp = normal_boiling_point_k + delta_tb

        logger.info(f"Boiling point elevation: ΔTb={delta_tb:.4f}K, new T_b={new_bp:.4f}K "
                     f"(m={molality}, Kb={ebullioscopic_constant_kb}, i={vanthoff_factor_i})")

        return {
            "delta_tb_k": round(delta_tb, 6),
            "new_boiling_point_k": round(new_bp, 6),
            "new_boiling_point_c": round(new_bp - 273.15, 4),
            "parameters_used": {
                "molality_mol_kg": molality,
                "ebullioscopic_constant_Kb": ebullioscopic_constant_kb,
                "vanthoff_factor_i": vanthoff_factor_i,
                "normal_boiling_point_K": normal_boiling_point_k,
            },
            "explanation": (
                f"ΔTb = i × Kb × m = {vanthoff_factor_i} × {ebullioscopic_constant_kb} × {molality} "
                f"= {delta_tb:.4f} K. New boiling point = {new_bp:.4f} K ({new_bp - 273.15:.2f} °C)."
            ),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入。"""
        try:
            parts = input_params.split()
            kwargs = {"molality": float(parts[0])}
            if len(parts) > 1: kwargs["ebullioscopic_constant_kb"] = float(parts[1])
            if len(parts) > 2: kwargs["vanthoff_factor_i"] = float(parts[2])
            if len(parts) > 3: kwargs["normal_boiling_point_k"] = float(parts[3])
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'molality [Kb] [i] [bp_K]'")
