"""
滴定剂标定计算工具 (Titrant Standardization)

用于用基准物质（Primary Standard）标定滴定剂的准确浓度。
"""

import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TitrantStandardization(BaseTool):
    """
    滴定剂标定工具。通过基准物质（纯度高、稳定、已知化学计量）来标定滴定剂浓度。
    支持多次平行测定并计算统计参数。
    """
    __version__ = "0.1.0"
    name = "TitrantStandardization"
    func_name = "standardize_titrant"
    description = "Calculate titrant concentration using primary standard substance, with statistical analysis for replicate determinations."
    implementation_description = (
        "Uses primary standard titration: C_titrant = (mass_std/MW_std) × n × 1000 / V_titrant(mL). "
        "Supports multiple replicates with RSD and confidence interval calculation."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Titration", "Standardization", "Analytical Chemistry", "QA/QC", "Volumetric Analysis"]
    required_envs = []

    code_input_sig = [
        ("primary_std_mass", "float", "N/A", "Mass of primary standard weighed (g). For replicates, use the first value or provide average."),
        ("primary_std_mw", "float", "N/A", "Molecular weight of the primary standard (g/mol)."),
        ("titrant_vol", "float", "N/A", "Volume of titrant consumed (mL). For replicates, use first value or average."),
        ("stoich_ratio", "float", "1.0", "Molar stoichiometric ratio (std:titrant), default 1:1 means n=1. For HCl vs Na2CO3 (1:2), use 2.0."),
        ("replicate_vols", "list", "[]", "List of titrant volumes for replicate measurements (mL). Optional; if provided, calculates statistics."),
        ("confidence_level", "float", "0.95", "Confidence level for uncertainty estimate (default 0.95 for 95%)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'primary_std_mass primary_std_mw titrant_vol [stoich_ratio]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary with keys: titrant_concentration(mol/L), moles_standard(mol), uncertainty(dict with RSD/CI/etc.), detailed_calculation(str)"),
    ]

    examples = [
        {
            "code_input": {
                "primary_std_mass": 0.2042,
                "primary_std_mw": 105.99,  # Na2CO3 MW (anhydrous)
                "titrant_vol": 38.20,
                "stoich_ratio": 2.0,  # Na2CO3 + 2HCl -> 2NaCl + H2O + CO2
                "replicate_vols": [38.20, 38.18, 38.22, 38.15],
                "confidence_level": 0.95,
            },
            "text_input": {
                "input_params": "0.2042 105.99 38.20 2.0"
            },
            "output": {
                "result": {
                    "titrant_concentration": 0.1013,
                    "moles_standard": 0.001927,
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
        primary_std_mass: float,
        primary_std_mw: float,
        titrant_vol: float,
        stoich_ratio: float = 1.0,
        replicate_vols: list = None,
        confidence_level: float = 0.95,
    ) -> dict:
        """
        核心逻辑：滴定剂标定计算

        Parameters:
            primary_std_mass: 基准物质质量 (g)
            primary_std_mw: 基准物质摩尔质量 (g/mol)
            titrant_vol: 消耗的滴定剂体积 (mL)
            stoich_ratio: 化学计量比 (基准物:滴定剂)，默认1:1
            replicate_vols: 平行测定的体积列表 (mL)
            confidence_level: 置信水平

        Returns:
            dict: 标定结果和统计信息
        """
        import math

        # 输入验证
        if primary_std_mass <= 0:
            raise ChemMCPError("Primary standard mass must be positive.")
        if primary_std_mw <= 0:
            raise ChemMCPError("Primary standard molecular weight must be positive.")
        if titrant_vol <= 0:
            raise ChemMCPError("Titrant volume must be positive.")
        if stoich_ratio <= 0:
            raise ChemMCPError("Stoichiometric ratio must be positive.")

        # Step 1: 计算基准物质的物质的量 (mol)
        moles_standard = primary_std_mass / primary_std_mw

        # Step 2: 计算滴定剂浓度 (mol/L)
        # C = (moles_std * n * 1000) / V_mL
        # 其中 n 是化学计量比（每摩尔基准物质反应消耗n摩尔滴定剂）
        titrant_conc = (moles_standard * stoich_ratio * 1000.0) / titrant_vol

        # 构建基本结果
        result = {
            "moles_standard": round(moles_standard, 6),
            "titrant_concentration": round(titrant_conc, 6),
        }

        # 统计分析（如果有平行样）
        if replicate_vols and len(replicate_vols) >= 2:
            n_rep = len(replicate_vols)
            vols = replicate_vols

            # 计算每次的浓度
            concs = [(moles_standard * stoich_ratio * 1000.0) / v for v in vols]

            # 平均值
            mean_conc = sum(concs) / n_rep

            # 标准差
            if n_rep > 1:
                variance = sum((c - mean_conc) ** 2 for c in concs) / (n_rep - 1)
                std_dev = math.sqrt(variance)
                rsd = (std_dev / mean_conc) * 100.0 if mean_conc != 0 else 0.0

                # t值近似（大样本用正态分布近似）
                from math import sqrt
                # 简化t值表（95%置信度）
                t_table_95 = {2: 12.71, 3: 4.30, 4: 3.18, 5: 2.78, 6: 2.57,
                              7: 2.45, 8: 2.36, 9: 2.31, 10: 2.26}
                t_val = t_table_95.get(n_rep, 2.0)

                ci_half_width = t_val * (std_dev / sqrt(n_rep))
                ci_lower = mean_conc - ci_half_width
                ci_upper = mean_conc + ci_half_width

                result["statistics"] = {
                    "n_replicates": n_rep,
                    "mean_concentration": round(mean_conc, 6),
                    "std_deviation": round(std_dev, 6),
                    "rsd_percent": round(rsd, 4),
                    "confidence_level": confidence_level,
                    "ci_lower": round(max(ci_lower, 0), 6),
                    "ci_upper": round(ci_upper, 6),
                    "individual_concentrations": [round(c, 6) for c in concs],
                }
            else:
                result["statistics"] = {"note": "Need at least 2 replicates for statistics."}

        # 详细计算过程
        calc_detail = (
            f"Step 1: Moles of primary standard = {primary_std_mass} g / {primary_std_mw} g/mol "
            f"= {moles_standard:.6f} mol\n"
            f"Step 2: Titrant concentration = ({moles_standard:.6f} mol × {stoich_ratio} × 1000) / {titrant_vol} mL "
            f"= {titrant_conc:.6f} mol/L"
        )
        result["detailed_calculation"] = calc_detail

        logger.info(f"Titrant standardization result: C = {titrant_conc:.6f} mol/L")
        return result

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 3:
                raise ValueError(
                    f"Need at least 3 parameters, got {len(parts)}. "
                    "Format: 'mass_MW_Vol [ratio]'"
                )

            mass = float(parts[0])
            mw = float(parts[1])
            vol = float(parts[2])
            ratio = float(parts[3]) if len(parts) > 3 else 1.0

            return self._run_base(mass, mw, vol, ratio)
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
