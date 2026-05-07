"""
返滴定计算工具 (Back Titration Solver)

用于计算返滴定（回滴定）中待测物的含量或纯度。
适用场景：反应速度较慢、无合适指示剂、待测物为不溶性固体等。
"""

import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BackTitrationSolver(BaseTool):
    """
    返滴定计算工具。通过加入过量标准溶液与待测物反应，再用另一种标准溶液回滴剩余量，
    从而计算待测物的纯度或含量。
    """
    __version__ = "0.1.0"
    name = "BackTitrationSolver"
    func_name = "solve_back_titration"
    description = "Calculate analyte purity or content using back titration (return titration) method."
    implementation_description = (
        "Uses back titration principle: excess standard solution is added to react with analyte, "
        "then the remaining excess is back-titrated with another standard solution. "
        "Formula: moles_analyte = (C_excess*V_excess - C_back*V_back) / n; purity = (moles*MW/mass)*100%"
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Titration", "Analytical Chemistry", "Quantitative Analysis", "Volumetric Analysis"]
    required_envs = []

    code_input_sig = [
        ("analyte_mass", "float", "N/A", "Mass of the sample containing the analyte (g)."),
        ("analyte_mw", "float", "N/A", "Molecular weight of the analyte (g/mol)."),
        ("titrant_conc", "float", "N/A", "Concentration of the excess titrant added (mol/L)."),
        ("excess_titrant_vol", "float", "N/A", "Volume of the excess titrant added (mL)."),
        ("back_titrant_conc", "float", "N/A", "Concentration of the back titrant used for return titration (mol/L)."),
        ("back_titrant_vol", "float", "N/A", "Volume of the back titrant consumed in return titration (mL)."),
        ("stoich_ratio", "float", "1.0", "Molar stoichiometric ratio (analyte:titrant), default 1:1 means n=1."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'analyte_mass analyte_mw titrant_conc excess_titrant_vol back_titrant_conc back_titrant_vol [stoich_ratio]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary with keys: moles_excess_titrant(mol), moles_reacted_titrant(mol), moles_analyte(mol), mass_analyte(g), analyte_purity(%), detailed_calculation(str)"),
    ]

    examples = [
        {
            "code_input": {
                "analyte_mass": 0.5000,
                "analyte_mw": 100.09,  # CaCO3 MW
                "titrant_conc": 0.1000,
                "excess_titrant_vol": 50.00,
                "back_titrant_conc": 0.1050,
                "back_titrant_vol": 12.50,
                "stoich_ratio": 1.0,
            },
            "text_input": {
                "input_params": "0.5000 100.09 0.1000 50.00 0.1050 12.50 1.0"
            },
            "output": {
                "result": {
                    "moles_excess_titrant": 0.001313,
                    "moles_reacted_titrant": 0.003687,
                    "moles_analyte": 0.003687,
                    "mass_analyte": 0.3688,
                    "analyte_purity": 73.77,
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
        analyte_mass: float,
        analyte_mw: float,
        titrant_conc: float,
        excess_titrant_vol: float,
        back_titrant_conc: float,
        back_titrant_vol: float,
        stoich_ratio: float = 1.0,
    ) -> dict:
        """
        核心逻辑：返滴定计算

        Parameters:
            analyte_mass: 待测样品质量 (g)
            analyte_mw: 待测物摩尔质量 (g/mol)
            titrant_conc: 过量滴定剂浓度 (mol/L)
            excess_titrant_vol: 加入的过量滴定剂体积 (mL)
            back_titrant_conc: 回滴用滴定剂浓度 (mol/L)
            back_titrant_vol: 回滴消耗的体积 (mL)
            stoich_ratio: 化学计量比 (待测物:滴定剂)，默认1:1

        Returns:
            dict: 包含各步计算结果和最终纯度的字典
        """
        # 输入验证
        if analyte_mass <= 0:
            raise ChemMCPError("Analyte mass must be positive.")
        if analyte_mw <= 0:
            raise ChemMCPError("Analyte molecular weight must be positive.")
        if titrant_conc <= 0 or back_titrant_conc <= 0:
            raise ChemMCPError("Titrant concentrations must be positive.")
        if excess_titrant_vol < 0 or back_titrant_vol < 0:
            raise ChemMCPError("Volumes cannot be negative.")
        if stoich_ratio <= 0:
            raise ChemMCPError("Stoichiometric ratio must be positive.")

        # Step 1: 计算回滴消耗的过量滴定剂的物质的量 (mol)
        # C_back * V_back (L) = moles of back titrant = moles of remaining excess titrant
        moles_excess_remaining = back_titrant_conc * (back_titrant_vol / 1000.0)

        # Step 2: 计算与待测物反应掉的滴定剂的物质的量 (mol)
        total_moles_excess_added = titrant_conc * (excess_titrant_vol / 1000.0)
        moles_reacted_titrant = total_moles_excess_added - moles_excess_remaining

        if moles_reacted_titrant < 0:
            raise ChemMCPError(
                f"Calculated reacted moles ({moles_reacted_titrant:.6f} mol) is negative. "
                "Check input values: back titrant volume may exceed the excess added."
            )

        # Step 3: 根据化学计量比计算待测物的物质的量 (mol)
        # analyte : titrant = 1 : n (stoich_ratio is n here)
        # moles_analyte = moles_reacted_titrant / stoich_ratio
        moles_analyte = moles_reacted_titrant / stoich_ratio

        # Step 4: 计算待测物质量 (g)
        mass_analyte = moles_analyte * analyte_mw

        # Step 5: 计算纯度 (%)
        analyte_purity = (mass_analyte / analyte_mass) * 100.0

        # 构建详细计算过程字符串
        calc_detail = (
            f"Step 1: Moles of excess titrant remaining = {back_titrant_conc} × {back_titrant_vol}/1000 "
            f"= {moles_excess_remaining:.6f} mol\n"
            f"Step 2: Total moles of excess titrant added = {titrant_conc} × {excess_titrant_vol}/1000 "
            f"= {total_moles_excess_added:.6f} mol\n"
            f"Step 3: Moles of titrant reacted with analyte = {total_moles_excess_added:.6f} - {moles_excess_remaining:.6f} "
            f"= {moles_reacted_titrant:.6f} mol\n"
            f"Step 4: Moles of analyte = {moles_reacted_titrant:.6f} / {stoich_ratio} "
            f"= {moles_analyte:.6f} mol\n"
            f"Step 5: Mass of analyte = {moles_analyte:.6f} × {analyte_mw} "
            f"= {mass_analyte:.4f} g\n"
            f"Step 6: Analyte purity = ({mass_analyte:.4f} / {analyte_mass}) × 100% "
            f"= {analyte_purity:.2f}%"
        )

        logger.info(f"Back titration result: purity={analyte_purity:.2f}%, mass_analyte={mass_analyte:.4f}g")

        return {
            "moles_excess_titrant": round(moles_excess_remaining, 6),
            "moles_reacted_titrant": round(moles_reacted_titrant, 6),
            "moles_analyte": round(moles_analyte, 6),
            "mass_analyte": round(mass_analyte, 4),
            "analyte_purity": round(analyte_purity, 2),
            "detailed_calculation": calc_detail,
        }

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入并调用核心逻辑"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 6:
                raise ValueError(
                    f"Need at least 6 parameters, got {len(parts)}. "
                    "Format: 'mass MW C_excess V_excess C_back V_back [ratio]'"
                )

            analyte_mass = float(parts[0])
            analyte_mw = float(parts[1])
            titrant_conc = float(parts[2])
            excess_titrant_vol = float(parts[3])
            back_titrant_conc = float(parts[4])
            back_titrant_vol = float(parts[5])
            stoich_ratio = float(parts[6]) if len(parts) > 6 else 1.0

            return self._run_base(
                analyte_mass, analyte_mw, titrant_conc, excess_titrant_vol,
                back_titrant_conc, back_titrant_vol, stoich_ratio,
            )
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
