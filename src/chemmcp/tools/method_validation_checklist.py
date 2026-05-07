"""
方法验证参数清单生成工具 (Method Validation Checklist)

基于 ICH Q2(R1) / FDA 分析程序验证指南，生成分析方法验证参数清单。
"""

import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ICH Q2(R1) 验证参数定义库
ICH_VALIDATION_PARAMETERS = {
    "Specificity": {
        "description": "Ability to assess unequivocally the analyte in the presence of components expected to be present (impurities, degradants, matrix).",
        "acceptance_criteria": {
            "chromatographic": "Peak purity index ≥ 0.999; Resolution between analyte and closest peak ≥ 1.5; No interference at analyte retention time from blank/matrix.",
            "spectroscopic": "Recovery of analyte in presence of interferents 98-102%; No spectral overlap.",
        },
        "recommended_tests": [
            "Blank injection (diluent/solvent)",
            "Placebo/sample matrix injection",
            "Forced degradation study (acid, base, oxidation, thermal, photolytic, humidity)",
            "Spiking with known impurities/degradants",
            "Peak purity assessment (PDA/MS)",
            "Resolution evaluation for all critical pairs",
        ],
    },
    "Linearity": {
        "description": "Demonstrates a proportional relationship between analyte concentration and instrument response over a specified range.",
        "acceptance_criteria": {
            "default": "Correlation coefficient (r) ≥ 0.999; Y-intercept not significantly different from zero (p > 0.05); Residual plot shows random distribution.",
            "loose": "r ≥ 0.995 for trace analysis or biological samples.",
        },
        "recommended_tests": [
            "Prepare minimum 5 concentration levels (e.g., 50%, 75%, 100%, 125%, 150% of target)",
            "Each level in triplicate",
            "Linear regression analysis (slope, intercept, r²)",
            "Residual analysis (plot residuals vs concentration)",
            "Lack-of-fit test (if appropriate)",
            "Back-calculate concentrations and report %bias",
        ],
    },
    "Range": {
        "description": "The interval between upper and lower levels of analyte that has been demonstrated to be determined with suitable precision, accuracy, and linearity.",
        "acceptance_criteria": {
            "default": "The range is established by confirming that the method provides acceptable precision, accuracy, and linearity across the range.",
        },
        "recommended_tests": [
            "Define range based on intended use (e.g., 80-120% of assay concentration, LOQ-150% for impurities)",
            "Verify precision and accuracy at range limits",
        ],
    },
    "Accuracy": {
        "description": "Closeness of test results to the true value; often expressed as recovery (%).",
        "acceptance_criteria": {
            "assay_drug_substance": "Mean recovery 98.0-102.0% at each level (n=3).",
            "assay_drug_product": "Mean recovery 98.0-102.0% at each level (n=3, 9 determinations total at 3 levels).",
            "impurities": "Mean recovery 80-120% at LOQ, 90-110% at 50-120% of specification.",
            "content_uniformity": "Mean recovery 85-115% (per USP <905>).",
        },
        "recommended_tests": [
            "Minimum 3 levels (e.g., 80%, 100%, 120% of target), n=3 per level (9 determinations total)",
            "Spike known amounts into placebo/matrix",
            "Calculate %recovery at each level",
            "Report mean recovery and RSD",
        ],
    },
    "Precision": {
        "description": "Degree of scatter between individual measurements under prescribed conditions.",
        "sub_categories": {
            "Repeatability": {
                "description": "Precision under same operating conditions over short interval (intra-assay).",
                "criteria": "RSD ≤ 2.0% for assay; RSD ≤ 10.0% for impurity at specification level; RSD ≤ 20.0% at LOQ.",
                "tests": ["Minimum 6 replicates at 100% concentration OR 9 determinations across 3 levels (3×3)."],
            },
            "Intermediate Precision": {
                "description": "Within-laborations variation: different days, analysts, instruments, reagents.",
                "criteria": "RSD ≤ 3.0% for assay; ANOVA showing no significant analyst/day/equipment effect (p > 0.05).",
                "tests": [
                    "Two analysts, each performing independent validation on different days",
                    "Different instruments where applicable",
                    "Different reagent lots where applicable",
                    "ANOVA for statistical comparison",
                ],
            },
            "Reproducibility": {
                "description": "Precision between laboratories (usually assessed via transfer).",
                "criteria": "Agreed acceptance criteria between labs; typically RSD ≤ 5.0%.",
                "tests": ["Method transfer with side-by-side testing."],
            },
        },
    },
    "Detection Limit (LOD)": {
        "description": "Lowest amount of analyte detectable but not necessarily quantified.",
        "acceptance_criteria": {
            "default": "Signal-to-noise ratio (S/N) ≥ 3:1 for chromatographic methods.",
        },
        "recommended_tests": [
            "Visual evaluation: Compare low-concentration samples with blank",
            "S/N ratio: Measure peak height vs baseline noise",
            "Standard deviation of response: LOD = 3.3 × σ/S",
            "Standard deviation of blank: LOD = 3.3 × σ_blank/S",
            "Calibration line method: LOD = 3.3 × SD_intercept / slope",
        ],
    },
    "Quantitation Limit (LOQ)": {
        "description": "Lowest amount of analyte that can be quantitatively determined with suitable precision and accuracy.",
        "acceptance_criteria": {
            "default": "S/N ≥ 10:1; Accuracy 80-120%; Precision (RSD) ≤ 10-20% at LOQ level.",
        },
        "recommended_tests": [
            "S/N ratio: S/N ≥ 10:1",
            "Standard deviation: LOQ = 10 × σ/S",
            "Verify with actual samples: prepare and analyze LOQ level in at least 6 replicates",
            "Confirm precision (RSD ≤ 10-20%) and accuracy (80-120%) at LOQ",
        ],
    },
    "Robustness": {
        "description": "Capacity of a method to remain unaffected by small, deliberate variations in method parameters.",
        "acceptance_criteria": {
            "default": "System suitability criteria met under all varied conditions; No significant change in resolution, tailing factor, or retention time.",
        },
        "recommended_tests": [
            "Organic phase composition (±2% absolute)",
            "pH of aqueous buffer (±0.2 units)",
            "Column temperature (±5°C)",
            "Flow rate (±10%)",
            "Wavelength (±2 nm for UV detection)",
            "Different columns (same type, different lot or supplier)",
            "Gradient program variations",
            "Injection volume variation",
            "Statistical evaluation (e.g., Plackett-Burman design recommended)",
        ],
    },
}

# 方法类型到验证参数的映射
METHOD_TYPE_PARAMS = {
    "HPLC": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "UPLC": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "GC": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "UV_Vis": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ"],
    "Titration": ["Specificity", "Accuracy", "Precision", "Robustness"],
    "Karl_Fischer": ["Specificity", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "AA": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "ICP_MS": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "ICP_OES": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "LC_MS": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "GC_MS": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
    "Dissolution": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "Robustness"],
    "General": ["Specificity", "Linearity", "Range", "Accuracy", "Precision", "LOD", "LOQ", "Robustness"],
}


@ChemMCPManager.register_tool
class MethodValidationChecklist(BaseTool):
    """
    方法验证参数清单生成工具。根据 ICH Q2(R1)/FDA 指南，
    为指定的分析方法类型生成完整的验证参数清单和接受标准。
    """
    __version__ = "0.1.0"
    name = "MethodValidationChecklist"
    func_name = "generate_validation_checklist"
    description = "Generate analytical method validation parameter checklist based on ICH Q2(R1)/FDA guidelines."
    implementation_description = (
        "Maps method types to ICH Q2(R1) validation parameters (Specificity, Linearity, Range, "
        "Accuracy, Precision, LOD, LOQ, Robustness) with detailed acceptance criteria and recommended tests."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Validation", "ICH Q2(R1)", "FDA", "Analytical Chemistry", "QA/QC", "Regulatory"]
    required_envs = []

    code_input_sig = [
        ("method_type", "str", "N/A", "Type of analytical method (HPLC/GC/UV_Vis/Titration/Karl_Fisher/AA/ICP_MS/LC_MS/Dissolution/etc.)."),
        ("analyte_name", "str", "Unknown Analyte", "Name of the target analyte."),
        ("matrix_type", "str", "N/A", "Sample matrix description (e.g., 'Tablet', 'Capsule', 'API', 'Plasma', 'Water')."),
        ("regulatory_framework", "str", "ICH_Q2R1", "Regulatory framework: 'ICH_Q2R1', 'FDA', 'USP', or 'All'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'method_type [analyte_name] [matrix_type] [regulatory_framework]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing: method_info, validation_parameters(dict with each param's details), "
         "test_schedule(recommended order), regulatory_references, total_tests_count"),
    ]

    examples = [
        {
            "code_input": {
                "method_type": "HPLC",
                "analyte_name": "Paracetamol",
                "matrix_type": "Tablet",
                "regulatory_framework": "ICH_Q2R1",
            },
            "text_input": {
                "input_params": "HPLC Paracetamol Tablet ICH_Q2R1"
            },
            "output": {
                "result": {
                    "method_info": {"method_type": "HPLC", "analyte_name": "Paracetamol"},
                    "validation_parameters": {"Specificity": "see_full_output"},
                    "total_tests_count": 40,
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
        analyte_name: str = "Unknown Analyte",
        matrix_type: str = "N/A",
        regulatory_framework: str = "ICH_Q2R1",
    ) -> dict:
        """
        核心逻辑：生成方法验证参数清单

        Parameters:
            method_type: 分析方法类型
            analyte_name: 待测物名称
            matrix_type: 样品基质类型
            regulatory_framework: 法规框架

        Returns:
            dict: 完整的验证参数清单
        """
        # 标准化方法类型
        method_key = method_type.upper().replace("-", "_").replace(" ", "_")
        if method_key not in METHOD_TYPE_PARAMS:
            # 尝试模糊匹配
            for key in METHOD_TYPE_PARAMS:
                if key in method_key or method_key in key:
                    method_key = key
                    break
            else:
                method_key = "General"
                logger.warning(f"Unknown method type '{method_type}', using General template.")

        # 获取该方法的验证参数列表
        param_names = METHOD_TYPE_PARAMS[method_key]

        # 构建每个参数的详细信息
        validation_params = {}
        total_tests = 0

        for param_name in param_names:
            if param_name in ICH_VALIDATION_PARAMETERS:
                param_data = ICH_VALIDATION_PARAMETERS[param_name]
                param_detail = {
                    "description": param_data.get("description", ""),
                    "acceptance_criteria": param_data.get("acceptance_criteria", "See guidelines"),
                    "recommended_tests": param_data.get("recommended_tests", []),
                    "status": "Pending",
                }

                # 特殊处理 Precision（包含子类别）
                if param_name == "Precision":
                    param_detail["sub_categories"] = param_data.get("sub_categories", {})

                validation_params[param_name] = param_detail
                total_tests += len(param_data.get("recommended_tests", []))

        # 推荐测试顺序（按逻辑依赖关系排序）
        test_order = self._get_test_order(param_names)

        # 法规引用
        references = self._get_references(regulatory_framework)

        result = {
            "method_info": {
                "method_type": method_type,
                "method_key": method_key,
                "analyte_name": analyte_name,
                "matrix_type": matrix_type,
                "regulatory_framework": regulatory_framework,
            },
            "validation_parameters": validation_params,
            "test_schedule": test_order,
            "regulatory_references": references,
            "total_recommended_tests": total_tests,
            "summary": f"Method validation checklist for {analyte_name} by {method_type} ({matrix_type}). "
                       f"Total {len(param_names)} parameters with ~{total_tests} individual tests.",
        }

        logger.info(f"Generated validation checklist for {method_type}: {len(param_names)} parameters")
        return result

    def _get_test_order(self, param_names: list) -> list:
        """返回推荐的测试执行顺序"""
        priority_map = {
            "Specificity": 1,
            "Linearity": 2,
            "Range": 3,
            "LOD": 4,
            "LOQ": 5,
            "Accuracy": 6,
            "Precision": 7,
            "Robustness": 8,
        }
        ordered = sorted(
            [(name, priority_map.get(name, 99)) for name in param_names],
            key=lambda x: x[1]
        )
        return [name for name, _ in ordered]

    def _get_references(self, framework: str) -> list:
        """获取法规引用"""
        ref_map = {
            "ICH_Q2R1": [
                "ICH Q2(R1): Validation of Analytical Procedures: Text and Methodology (2005)",
                "ICH Q2(R1) Step 4 version, November 2005",
            ],
            "FDA": [
                "FDA Guidance for Industry: Analytical Procedures and Methods Validation (2000, Draft 2015)",
                "FDA Guidance: Bioanalytical Method Validation (2018)",
            ],
            "USP": [
                "USP <1225> Validation of Compendial Procedures",
                "USP <1033> Biological Assay Validation",
                "USP <621> Chromatography",
                "USP <1058> Analytical Instrument Qualification",
            ],
            "All": [
                "ICH Q2(R1): Validation of Analytical Procedures (2005)",
                "FDA Guidance: Analytical Procedures and Methods Validation (2015 Draft)",
                "USP <1225> Validation of Compendial Procedures",
            ],
        }
        return ref_map.get(framework, ref_map["ICH_Q2R1"])

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            if len(parts) < 1:
                raise ValueError("Need at least method_type.")

            method_type = parts[0]
            analyte = parts[1] if len(parts) > 1 else "Unknown Analyte"
            matrix = parts[2] if len(parts) > 2 else "N/A"
            framework = parts[3] if len(parts) > 3 else "ICH_Q2R1"

            return self._run_base(method_type, analyte, matrix, framework)
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
