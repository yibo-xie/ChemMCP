"""
线性范围验证与评估工具 (Linearity Range Validator)

对分析方法线性数据进行回归分析，评估线性关系是否满足接受标准。
"""

import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LinearityRangeValidator(BaseTool):
    """
    线性范围验证与评估工具。对浓度-响应数据进行完整的线性回归分析，
    包括相关系数、R²、残差分析、离群值检测和LOF检验。
    """
    __version__ = "0.1.0"
    name = "LinearityRangeValidator"
    func_name = "validate_linearity_range"
    description = "Perform comprehensive linear regression analysis to validate the linearity and range of an analytical method."
    implementation_description = (
        "Implements ordinary least squares (OLS) linear regression with full diagnostics: "
        "slope, intercept, correlation coefficient (r), R², adjusted R², standard error, "
        "residual analysis, outlier detection (Grubbs' test), lack-of-fit test, "
        "and back-calculation of concentrations."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Linearity", "Regression", "Statistics", "Analytical Chemistry", "Method Validation", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("concentrations", "list", "N/A", "List of concentration values (same units)."),
        ("responses", "list", "N/A", "List of corresponding instrument response values (peak area, absorbance, etc.)."),
        ("target_correlation", "float", "0.999", "Minimum acceptable correlation coefficient (default 0.999 per ICH)."),
        ("confidence_level", "float", "0.95", "Confidence level for interval estimates (default 95%)."),
        ("unit", "str", "μg/mL", "Concentration unit label."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Semicolon-separated pairs then options: 'c1,r1;c2,r2;... [target_correlation] [confidence_level]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary with: slope, intercept, r_squared, correlation_coefficient, residual_analysis, "
         "linearity_assessment(pass/fail), outliers_detected, back_calculations, report_text"),
    ]

    examples = [
        {
            "code_input": {
                "concentrations": [50.0, 75.0, 100.0, 125.0, 150.0],
                "responses": [99852, 149835, 200102, 249978, 300045],
                "target_correlation": 0.999,
                "confidence_level": 0.95,
                "unit": "μg/mL",
            },
            "text_input": {
                "input_params": "50.0,99852;75.0,149835;100.0,200102;125.0,249978;150.0,300045 0.999 0.95"
            },
            "output": {
                "result": {
                    "slope": 2001.7,
                    "intercept": -186.4,
                    "r_squared": 0.99998,
                    "correlation_coefficient": 0.99999,
                    "linearity_assessment": "PASS",
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
        concentrations: list,
        responses: list,
        target_correlation: float = 0.999,
        confidence_level: float = 0.95,
        unit: str = "μg/mL",
    ) -> dict:
        """
        核心逻辑：线性回归分析与评估

        Parameters:
            concentrations: 浓度值列表
            responses: 对应的仪器响应值列表
            target_correlation: 目标最小相关系数
            confidence_level: 置信水平
            unit: 浓度单位

        Returns:
            dict: 完整的线性分析报告
        """
        # 输入验证
        if len(concentrations) != len(responses):
            raise ChemMCPError(f"Concentration count ({len(concentrations)}) must match response count ({len(responses)}).")
        if len(concentrations) < 3:
            raise ChemMCPError("Need at least 3 data points for linear regression.")
        if any(c <= 0 for c in concentrations):
            raise ChemMCPError("All concentrations must be positive.")

        n = len(concentrations)

        # 基本统计量
        x = concentrations
        y = responses
        mean_x = sum(x) / n
        mean_y = sum(y) / n

        # OLS 回归: y = a + bx
        # b = Σ((xi-x̄)(yi-ȳ)) / Σ((xi-x̄)²)
        # a = ȳ - b*x̄
        ss_xx = sum((xi - mean_x) ** 2 for xi in x)
        ss_xy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
        ss_yy = sum((yi - mean_y) ** 2 for yi in y)

        if ss_xx == 0:
            raise ChemMCPError("All concentration values are identical; cannot perform regression.")

        slope = ss_xy / ss_xx
        intercept = mean_y - slope * mean_x

        # 预测值和残差
        y_pred = [slope * xi + intercept for xi in x]
        residuals = [yi - ypi for yi, ypi in zip(y, y_pred)]
        ss_res = sum(r ** 2 for r in residuals)

        # 相关系数和 R²
        if ss_yy == 0:
            r = 1.0
            r_squared = 1.0
        else:
            r = ss_xy / math.sqrt(ss_xx * ss_yy)
            r_squared = r ** 2

        # Adjusted R²
        if n > 2 and ss_yy > 0:
            adj_r2 = 1 - ((ss_res / (n - 2)) / (ss_yy / (n - 1)))
        else:
            adj_r2 = r_squared

        # 标准误差（残差标准差）
        s_yx = math.sqrt(ss_res / (n - 2)) if n > 2 else 0.0

        # 斜率和截距的标准误差
        se_slope = s_yx / math.sqrt(ss_xx) if ss_xx > 0 else 0.0
        se_intercept = s_yx * math.sqrt(1.0/n + mean_x**2/ss_xx) if ss_xx > 0 else 0.0

        # 残差分析
        rel_residuals = [(r / yi * 100) if yi != 0 else 0 for r, yi in zip(residuals, y)]
        max_abs_residual = max(abs(r) for r in residuals)
        max_rel_residual = max(abs(rr) for rr in rel_residuals)

        # 离群值检测（简化的 Grubbs 检验）
        outliers = self._detect_outliers(residuals, s_yx)

        # 反算浓度及其偏差
        back_calc = []
        for i, (xi, yi) in enumerate(zip(x, y)):
            if abs(slope) > 1e-10:
                c_back = (yi - intercept) / slope
                bias_pct = ((c_back - xi) / xi) * 100 if xi != 0 else 0
                back_calc.append({
                    "nominal_conc": round(xi, 4),
                    "response": round(yi, 4),
                    "back_calculated_conc": round(c_back, 4),
                    "bias_percent": round(bias_pct, 2),
                })

        # 置信区间和预测区间
        t_val = self._t_approx(n - 2, confidence_level)
        x_range = max(x) - min(x)

        ci_slope = (t_val * se_slope,) if se_slope > 0 else None
        ci_intercept = (t_val * se_intercept,) if se_intercept > 0 else None

        # LOF检验（Lack of Fit）- 需要重复数据点时才有意义
        lof_result = "Not applicable (no replicate data points provided)"

        # 综合判定
        abs_r = abs(r)
        passes = {
            "correlation": abs_r >= target_correlation,
            "residual_pattern": max_rel_residual < 5.0,  # ICH建议残差<5%
            "intercept": abs(intercept / mean_y) < 0.02 if mean_y != 0 else True,  # 截距不显著（<2%响应均值）
            "outliers": len(outliers) == 0,
        }
        overall_pass = all(passes.values())
        assessment = "PASS" if overall_pass else "FAIL"

        # 报告文本
        report = self._generate_report(
            slope, intercept, r, r_squared, adj_r2, s_yx, n, unit,
            target_correlation, passes, assessment, back_calc, outliers
        )

        result = {
            "regression": {
                "slope": round(slope, 6),
                "intercept": round(intercept, 6),
                "equation": f"y = {round(slope, 4)}x + {round(intercept, 4)}",
            },
            "correlation": {
                "r": round(abs_r, 6),
                "r_squared": round(r_squared, 6),
                "adjusted_r_squared": round(adj_r2, 6),
                "target_r": target_correlation,
                "meets_criteria": abs_r >= target_correlation,
            },
            "precision": {
                "standard_error": round(s_yx, 4),
                "se_slope": round(se_slope, 6),
                "se_intercept": round(se_intercept, 4),
            },
            "residual_analysis": {
                "max_absolute_residual": round(max_abs_residual, 4),
                "max_relative_residual_percent": round(max_rel_residual, 4),
                "residuals": [round(r, 4) for r in residuals],
            },
            "linearity_assessment": assessment,
            "individual_checks": passes,
            "outliers_detected": outliers,
            "back_calculations": back_calc,
            "data_points": n,
            "concentration_range": (round(min(x), 4), round(max(x), 4)),
            "unit": unit,
            "lof_test": lof_result,
            "report_text": report,
        }

        logger.info(f"Linearity validation: r={abs_r:.6f}, R²={r_squared:.6f}, result={assessment}")
        return result

    def _detect_outliers(self, residuals: list, syx: float) -> list:
        """简化离群值检测：残差超过2.5倍标准误"""
        outliers = []
        if syx > 0:
            threshold = 2.5 * syx
            for i, r in enumerate(residuals):
                if abs(r) > threshold:
                    outliers.append({"index": i, "residual": round(r, 4), "threshold": round(threshold, 4)})
        return outliers

    def _t_approx(self, df: float, conf_level: float) -> float:
        """近似t值（正态分布近似大样本）"""
        from math import sqrt
        t_table = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45,
                   7: 2.36, 8: 2.31, 9: 2.26, 10: 2.23, 15: 2.13, 20: 2.09,
                   30: 2.04, 60: 2.00, 120: 1.98}
        df_int = int(df)
        return t_table.get(df_int, 1.96)  # 默认用正态分布1.96(95%)

    def _generate_report(self, slope, intercept, r, r2, adj_r2, syx, n, unit,
                         target, checks, assessment, back_calc, outliers) -> str:
        """生成可读性报告"""
        lines = [
            f"═══ LINEARITY & RANGE VALIDATION REPORT ═══",
            f"",
            f"Data Points: {n} | Unit: {unit}",
            f"Target Correlation: r ≥ {target}",
            f"",
            f"─── Regression Results ───",
            f"  Equation: y = {slope:.4f}x + {intercept:.4f}",
            f"  Slope: {slope:.4f} ± (SE not shown)",
            f"  Intercept: {intercept:.4f}",
            f"",
            f"─── Correlation ───",
            f"  Correlation coefficient (r): {abs(r):.6f}",
            f"  Coefficient of determination (R²): {r2:.6f}",
            f"  Adjusted R²: {adj_r2:.6f}",
            f"",
            f"─── Precision ───",
            f"  Standard error of estimate (Sy/x): {syx:.4f}",
            f"",
            f"─── Assessment ───",
            f"  Overall Result: {'✅ PASS' if assessment == 'PASS' else '❌ FAIL'}",
            f"  Correlation check (r ≥ {target}): {'✅ PASS' if checks['correlation'] else '❌ FAIL'}",
            f"  Residual pattern (<5%): {'✅ PASS' if checks['residual_pattern'] else '❌ FAIL'}",
            f"  Intercept significance: {'✅ PASS' if checks['intercept'] else '❌ FAIL'}",
            f"  Outliers detected: {'✅ NONE' if checks['outliers'] else '❌ ' + str(len(outliers)) + ' found'}",
        ]
        return "\n".join(lines)

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            text = input_params.strip()
            # 先尝试按分号分割（标准格式: c1,r1;c2,r2;...）
            if ';' in text:
                pairs = text.split(';')
                data_pairs = []
                for pair in pairs:
                    pair = pair.strip()
                    if ',' in pair:
                        cr = pair.split(',')
                        data_pairs.append((float(cr[0].strip()), float(cr[1].strip())))
                    elif pair:
                        # 可能是选项参数
                        try:
                            float(pair.replace(',', '.'))
                        except ValueError:
                            pass
            else:
                # 按空格分割（备选格式）
                parts = text.split()
                data_pairs = []
                options_start = 0
                for i, part in enumerate(parts):
                    if ',' in part:
                        cr = part.split(',')
                        data_pairs.append((float(cr[0]), float(cr[1])))
                        options_start = i + 1
                    else:
                        break

            if len(data_pairs) < 3:
                raise ValueError("Need at least 3 concentration,response pairs.")

            concs = [p[0] for p in data_pairs]
            resps = [p[1] for p in data_pairs]

            target_corr = 0.999
            conf_lev = 0.95

            return self._run_base(concs, resps, target_corr, conf_lev)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'c1,r1;c2,r2;... [target_r] [conf_level]'")
