"""
稳健性实验设计工具 (Robustness DOE Designer)

生成 Plackett-Burman 筛选设计矩阵，用于分析方法稳健性评估。
"""

import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 标准 Plackett-Burman 设计矩阵（生成向量）
# 每个设计的首行（生成器），其余行通过循环移位得到
PB_DESIGNS = {
    # N=3: 2因子 + 1虚拟
    3: {
        "generators": [[1, 1, -1]],
        "max_factors": 2,
    },
    # N=4: 3因子（部分析因）
    4: {
        "generators": [[1, 1, 1, -1]],
        "max_factors": 3,
    },
    # N=7: 4-7因子（实际是8-1=7行）
    7: {
        "generators": [[1, 1, 1, -1, 1, -1, -1]],
        "max_factors": 7,
    },
    # N=11: 8-11因子（最常用）
    11: {
        "generators": [[1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1]],
        "max_factors": 11,
    },
    # N=15: 12-15因子
    15: {
        "generators": [[1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, -1]],
        "max_factors": 15,
    },
    # N=19: 16-19因子
    19: {
        "generators": [[1, 1, -1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1, 1, 1]],
        "max_factors": 19,
    },
    # N=23: 20-23因子
    23: {
        "generators": [[1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1]],
        "max_factors": 23,
    },
}


@ChemMCPManager.register_tool
class RobustnessDoeDesigner(BaseTool):
    """
    稳健性实验设计工具。基于 Plackett-Burman 筛选设计，
    为分析方法稳健性评估生成完整的实验设计方案。
    """
    __version__ = "0.1.0"
    name = "RobustnessDoeDesigner"
    func_name = "design_robustness_doe"
    description = "Generate Plackett-Burman screening design matrix for analytical method robustness evaluation."
    implementation_description = (
        "Generates Plackett-Burman (PB) design matrices for screening main effects in robustness studies. "
        "Supports standard PB designs with N=3/4/7/11/15/19/23 runs. "
        "Includes effect estimation formulas, aliasing structure analysis, and ICH-compliant factor recommendations."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["DOE", "Robustness", "Experimental Design", "ICH", "Analytical Chemistry", "Quality by Design"]
    required_envs = []

    code_input_sig = [
        ("factors", "list", "N/A",
         "List of factor dicts: [{'name': 'pH', 'low': 3.0, 'high': 5.0, 'unit': ''}, ...]. "
         "Each factor must have name, low, high values."),
        ("num_runs", "int", "0",
         "Number of PB runs (0=auto-select based on factor count). Options: 3,4,7,11,15,19,23."),
        ("response_name", "str", "Response",
         "Name of the response variable to evaluate (e.g., 'Resolution', 'Tailing_Factor', 'Plate_Count')."),
        ("include_center_point", "bool", "True",
         "Whether to include center point(s) in the design for curvature detection."),
        ("center_replicates", "int", "3",
         "Number of center point replicates if include_center_point is True."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "JSON-like string or semicolon-separated factors. "
         "Format: 'factor1_name,low,high;factor2_name,low,high;... [num_runs]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing: design_matrix(list of dicts with run/factor levels), "
         "effect_estimation_formula, aliasing_structure, analysis_instructions, report"),
    ]

    examples = [
        {
            "code_input": {
                "factors": [
                    {"name": "Organic_Ratio", "low": 38.0, "high": 42.0, "unit": "%"},
                    {"name": "Buffer_pH", "low": 2.8, "high": 3.2, "unit": ""},
                    {"name": "Flow_Rate", "low": 0.9, "high": 1.1, "unit": "mL/min"},
                    {"name": "Column_Temp", "low": 28, "high": 32, "unit": "°C"},
                    {"name": "Wavelength", "low": 252, "high": 258, "unit": "nm"},
                    {"name": "Injection_Vol", "low": 8, "high": 12, "unit": "μL"},
                    {"name": "Gradient_Start", "low": -2, "high": 2, "unit": "%"},
                ],
                "num_runs": 0,
                "response_name": "Resolution_RS1",
                "include_center_point": True,
                "center_replicates": 3,
            },
            "text_input": {
                "input_params": "Organic_Ratio,38,42%;Buffer_pH,2.8,3.2;Flow_Rate,0.9,1.1,mL/min;Column_Temp,28,32,C;Wavelength,252,258,nm;Injection_Vol,8,12,uL 11"
            },
            "output": {
                "result": {
                    "design_info": {"num_runs": 14, "num_factors": 7},
                    "design_matrix": [...],
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
        factors: list,
        num_runs: int = 0,
        response_name: str = "Response",
        include_center_point: bool = True,
        center_replicates: int = 3,
    ) -> dict:
        """
        核心逻辑：生成 Plackett-Burman 设计矩阵

        Parameters:
            factors: 因子列表，每个因子包含 name, low, high, unit
            num_runs: PB设计运行次数（0=自动选择）
            response_name: 响应变量名
            include_center_point: 是否包含中心点
            center_replicates: 中心点重复次数

        Returns:
            dict: 完整的实验设计方案
        """
        # 输入验证
        if not factors or len(factors) < 2:
            raise ChemMCPError("Need at least 2 factors for Plackett-Burman design.")
        num_factors = len(factors)

        for i, f in enumerate(factors):
            if "name" not in f or "low" not in f or "high" not in f:
                raise ChemMCPError(f"Factor {i} missing required fields (name, low, high).")
            if f["low"] == f["high"]:
                raise ChemMCPError(f"Factor '{f['name']}' has identical low and high values.")

        # 自动选择运行次数
        if num_runs == 0:
            num_runs = self._select_pb_runs(num_factors)
        elif num_runs not in PB_DESIGNS:
            available = sorted(PB_DESIGNS.keys())
            raise ChemMCPError(f"Invalid num_runs={num_runs}. Available: {available}")

        pb_design = PB_DESIGNS[num_runs]
        max_factors = pb_design["max_factors"]
        generator = pb_design["generators"][0]

        if num_factors > max_factors:
            raise ChemMCPError(
                f"Too many factors ({num_factors}) for {num_runs}-run PB design (max {max_factors}). "
                f"Use at least {self._select_pb_runs(num_factors)} runs."
            )

        # 生成 PB 矩阵
        matrix = self._generate_pb_matrix(generator, num_runs)

        # 构建设计表（将编码值映射到实际水平）
        design_table = []
        for run_idx, row in enumerate(matrix):
            run_data = {"Run": run_idx + 1}
            for col_idx, factor in enumerate(factors):
                level_code = row[col_idx] if col_idx < len(row) else row[-1]
                actual_level = factor["high"] if level_code == 1 else factor["low"]
                run_data[factor["name"]] = actual_level
                run_data[f"{factor['name']}_coded"] = level_code
            design_table.append(run_data)

        # 添加中心点
        center_points = []
        if include_center_point:
            for cp_idx in range(center_replicates):
                cp_run = {"Run": len(design_table) + cp_idx + 1, "_type": "Center_Point"}
                for factor in factors:
                    center_val = (factor["low"] + factor["high"]) / 2.0
                    cp_run[factor["name"]] = center_val
                    cp_run[f"{factor['name']}_coded"] = 0
                center_points.append(cp_run)

        total_runs = len(design_table) + len(center_points)

        # 效应估计公式
        effect_formula = self._generate_effect_formula(factors, num_runs, response_name)

        # 别名结构分析
        aliasing = self._analyze_aliasing(num_factors, num_runs)

        # 分析指导
        analysis_instructions = self._get_analysis_instructions(response_name, num_factors, total_runs, center_replicates)

        # 因子汇总
        factor_summary = []
        for f in factors:
            range_val = abs(f["high"] - f["low"])
            unit_str = f.get("unit", "")
            factor_summary.append({
                "name": f["name"],
                "low": f["low"],
                "high": f["high"],
                "center": round((f["low"] + f["high"]) / 2.0, 4),
                "range": round(range_val, 4),
                "unit": unit_str,
                "delta_for_effect": round(range_val / 2.0, 4),
            })

        result = {
            "design_info": {
                "design_type": "Plackett-Burman Screening Design",
                "num_factors": num_factors,
                "pb_runs": num_runs,
                "center_point_runs": center_replicates if include_center_point else 0,
                "total_runs": total_runs,
                "response_variable": response_name,
            },
            "factor_summary": factor_summary,
            "design_matrix": design_table,
            "center_points": center_points if include_center_point else [],
            "effect_estimation": effect_formula,
            "aliasing_structure": aliasing,
            "analysis_instructions": analysis_instructions,
            "report_text": self._generate_report(
                num_factors, num_runs, total_runs, factors, response_name, aliasing
            ),
        }

        logger.info(f"Generated PB design: {num_factors} factors, {total_runs} runs ({num_runs} PB + {len(center_points)} CP)")
        return result

    def _select_pb_runs(self, n_factors: int) -> int:
        """根据因子数自动选择合适的PB运行次数"""
        if n_factors <= 2:
            return 3
        elif n_factors <= 3:
            return 4
        elif n_factors <= 7:
            return 7
        elif n_factors <= 11:
            return 11
        elif n_factors <= 15:
            return 15
        elif n_factors <= 19:
            return 19
        else:
            return 23

    def _generate_pb_matrix(self, generator: list, n: int) -> list:
        """从生成向量循环移位生成完整PB矩阵"""
        matrix = [generator[:]]
        for i in range(1, n):
            row = generator[-i:] + generator[:-i]
            matrix.append(row)
        return matrix

    def _generate_effect_formula(self, factors: list, n_runs: int, response: str) -> dict:
        """生成效应估计公式"""
        factor_names = [f["name"] for f in factors]
        formula_parts = []
        for fn in factor_names:
            formula_parts.append(f"E({fn}) = (Σ(R_i × X_{{{fn}}})) / (N/2)")

        return {
            "general_formula": "Effect = (Σ(Response at High) - Σ(Response at Low)) / (N/2)",
            "per_factor_formulas": formula_parts,
            "interpretation": (
                "|Effect| > critical value → significant factor. "
                "For α=0.05 with N runs, approximate critical effect ≈ 2×SD_response/√N."
            ),
            "pareto_chart_guidance": "Sort |Effect| descending; plot Pareto chart to identify critical factors.",
        }

    def _analyze_aliasing(self, n_factors: int, n_runs: int) -> dict:
        """分析别名结构"""
        # 在PB设计中，主效应与二阶及更高阶交互作用混淆
        degrees_of_freedom = n_runs - 1
        main_effects_df = n_factors
        remaining_df = degrees_of_freedom - main_effects_df
        dummy_factors = max(0, remaining_df)

        return {
            "total_degrees_of_freedom": degrees_of_freedom,
            "main_effects_dof": min(main_effects_df, degrees_of_freedom),
            "dummy_or_error_dof": max(0, remaining_df),
            "aliasing_pattern": (
                "In Plackett-Burman designs, each main effect is aliased with "
                "all two-factor and higher-order interactions involving other factors. "
                "This is a Resolution III design."
            ),
            "note": (
                f"With {n_factors} factors in {n_runs}-run PB design, there are "
                f"{dummy_factors} dummy factor columns available for error estimation."
            ),
        }

    def _get_analysis_instructions(self, response: str, n_factors: int, total_runs: int, n_cp: int) -> list:
        """返回分析步骤指导"""
        steps = [
            f"1. Execute all {total_runs} experiments in randomized order.",
            f"2. Record the response value ('{response}') for each run.",
            "3. Calculate the main effect for each factor using the effect estimation formula.",
            "4. Estimate experimental error from dummy factor columns (if available) or center points.",
            "5. Compare |Effect| vs critical value (or use normal/half-normal probability plot).",
            "6. Identify significant factors (|Effect| > critical value).",
            "7. For significant factors, verify that system suitability criteria are still met.",
            "8. Document conclusions: method is robust if no factor significantly affects the response.",
        ]
        if n_cp > 0:
            steps.insert(3, f"3b. Check curvature: compare center point mean vs factorial point mean (t-test).")
        return steps

    def _generate_report(self, n_factors, n_runs, total_runs, factors, response, aliasing) -> str:
        """生成报告文本"""
        lines = [
            f"═══ PLACKETT-BURMAN ROBUSTNESS DESIGN ═══",
            f"",
            f"Design Type: Plackett-Burman Screening",
            f"Factors: {n_factors}",
            f"PB Runs: {n_runs} | Center Points: {total_runs - n_runs} | Total: {total_runs}",
            f"Response Variable: {response}",
            f"",
            f"─── Factors ───",
        ]
        for f in factors:
            unit = f.get("unit", "")
            lines.append(f"  • {f['name']}: {f['low']} ~ {f['high']} {unit}")

        lines.extend([
            f"",
            f"─── Aliasing Structure ───",
            f"  {aliasing['aliasing_pattern']}",
            f"",
            f"─── Notes ───",
            f"  • Randomize run order before execution.",
            f"  • This is a screening design for main effects only.",
            f"  • Significant factors should be studied further with response surface methodology.",
        ])
        return "\n".join(lines)

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            factors = []
            runs = 0

            for part in parts:
                if ',' in part and not part.replace(',', '').replace('.', '').replace('-', '').replace('%', '').isdigit():
                    sub_parts = part.split(',')
                    if len(sub_parts) >= 3:
                        unit = sub_parts[3] if len(sub_parts) > 3 else ""
                        # 清理单位字符串
                        unit = unit.replace('uL', 'μL')
                        factors.append({
                            "name": sub_parts[0],
                            "low": float(sub_parts[1]),
                            "high": float(sub_parts[2]),
                            "unit": unit,
                        })
                elif part.isdigit() and int(part) in PB_DESIGNS:
                    runs = int(part)

            if len(factors) < 2:
                raise ValueError(f"Need at least 2 factors, got {len(factors)}.")

            return self._run_base(factors, runs)
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
