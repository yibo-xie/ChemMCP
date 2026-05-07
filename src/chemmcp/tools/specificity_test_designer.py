"""
专属性测试方案设计工具 (Specificity Test Designer)

基于 ICH Q2(R1) 指南设计分析方法专属性测试方案，
包括强制降解条件、峰纯度要求、分离度标准等。
"""

import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# 强制降解条件参考库（基于ICH Q1B和行业实践）
FORCED_DEGRADATION_CONDITIONS = {
    "acidic": {
        "name": "Acid Degradation",
        "description": "Hydrolysis under acidic conditions",
        "typical_conditions": [
            {"condition": "0.1 M HCl, room temperature, 1 hour"},
            {"condition": "0.1 M HCl, 60°C, 30 minutes"},
            {"condition": "1 M HCl, room temperature, 1 hour (stress)"},
            {"condition": "0.01 M HCl, reflux, 30 minutes (mild)"},
        ],
        "target_degradation": "5-20% degradation (aim for ~10-15%)",
        "neutralization": "Neutralize with equimolar NaOH before analysis",
        "acceptance": "Peak purity of analyte ≥ 0.999; No co-eluting degradants at analyte RT",
    },
    "basic": {
        "name": "Base Degradation",
        "description": "Hydrolysis under basic conditions",
        "typical_conditions": [
            {"condition": "0.1 M NaOH, room temperature, 1 hour"},
            {"condition": "0.1 M NaOH, 60°C, 30 minutes"},
            {"condition": "0.01 M NaOH, room temperature, 24 hours (mild)"},
        ],
        "target_degradation": "5-20% degradation (aim for ~10-15%)",
        "neutralization": "Neutralize with equimolar HCl before analysis",
        "acceptance": "Peak purity of analyte ≥ 0.999; No co-eluting degradants at analyte RT",
    },
    "oxidative": {
        "name": "Oxidative Degradation",
        "description": "Oxidation with hydrogen peroxide",
        "typical_conditions": [
            {"condition": "3% H2O2, room temperature, 1 hour (standard)"},
            {"condition": "3% H2O2, 40°C, 30 minutes (accelerated)"},
            {"condition": "6% H2O2, room temperature, 30 minutes (stress)"},
            {"condition": "0.3% H2O2, room temperature, 24 hours (mild)"},
            {"condition": "15% H2O2, room temperature, 10 minutes (extreme stress)"},
        ],
        "target_degradation": "5-20% degradation (aim for ~10-15%)",
        "neutralization": "May dilute to quench excess peroxide; avoid if unstable",
        "acceptance": "Peak purity of analyte ≥ 0.999; Mass balance within 95-105%",
    },
    "thermal": {
        "name": "Thermal (Dry Heat) Degradation",
        "description": "Solid-state thermal stress",
        "typical_conditions": [
            {"condition": "Solid drug substance, 60°C, 3 days (dry heat)"},
            {"condition": "Solid drug substance, 80°C, 24 hours (stress)"},
            {"condition": "Solid drug substance, 105°C, 6 hours (extreme)"},
        ],
        "target_degradation": "5-20% degradation",
        "neutralization": "N/A (solid state)",
        "acceptance": "Peak purity of analyte ≥ 0.999; Identify degradation products",
    },
    "photolytic": {
        "name": "Photolytic Degradation",
        "description": "Light exposure per ICH Q1B",
        "typical_conditions": [
            {"condition": "ICH Q1B Option 1: D65/ID65 lamp, ≥1.2×10^6 lux hours + ≥200 W·h/m² UV"},
            {"condition": "ICH Q1B Option 2: Cool white fluorescent + near UV fluorescent lamp"},
            {"condition": "UVA/UVB chamber, 24-48 hours exposure"},
        ],
        "target_degradation": "Visible change or 5-20% degradation",
        "neutralization": "N/A",
        "acceptance": "Compare with dark control; Peak purity ≥ 0.999; Report % degradation",
    },
    "humidity": {
        "name": "Humidity Degradation",
        "description": "High humidity stress (solid state)",
        "typical_conditions": [
            {"condition": "75% RH ± 5%, 25°C, 7 days (solid in open dish)"},
            {"condition": "90% RH ± 5%, 25°C, 5 days (stress)"},
            {"condition": "40°C/75% RH, 7 days (accelerated humidity)"},
        ],
        "target_degradation": "5-20% degradation or visible change",
        "neutralization": "N/A",
        "acceptance": "Peak purity ≥ 0.999; Compare with dry control sample",
    },
}


# 方法类型特定的接受标准
METHOD_SPECIFICITY_CRITERIA = {
    "HPLC": {
        "resolution": "≥ 1.5 between all critical pairs (analyte vs nearest peak)",
        "tailing_factor": "≤ 2.0 (preferably ≤ 1.5)",
        "peak_purity_index": "≥ 0.999 (PDA detection)",
        "theoretical_plates": "≥ 2000 (for analyte peak)",
        "mass_balance": "95.0-105.0% (sum of analyte + degradants)",
    },
    "UPLC": {
        "resolution": "≥ 1.5 between all critical pairs",
        "tailing_factor": "≤ 2.0 (preferably ≤ 1.5)",
        "peak_purity_index": "≥ 0.999",
        "theoretical_plates": "≥ 5000 (for analyte peak)",
        "mass_balance": "95.0-105.0%",
    },
    "GC": {
        "resolution": "≥ 1.0 between all critical pairs",
        "tailing_factor": "≤ 2.0",
        "peak_purity": "MS confirmation or orthogonal detection recommended",
        "mass_balance": "95.0-105.0%",
    },
    "UV_Vis": {
        "specificity": "Recovery 98-102% in presence of interferents; No spectral overlap confirmed by derivative spectra or multi-wavelength analysis",
    },
    "Titration": {
        "specificity": "No interference from placebo/excipients; Recovery 99-101%",
    },
    "General": {
        "resolution": "Baseline separation from all known components",
        "peak_purity": "Confirmed by orthogonal method or PDA/MS",
        "recovery_in_matrix": "98-102% in presence of expected interferents",
    },
}


@ChemMCPManager.register_tool
class SpecificityTestDesigner(BaseTool):
    """
    专属性测试方案设计工具。基于 ICH Q2(R1) 指南，为分析方法设计完整的专属性测试方案。
    包括强制降解实验设计、峰纯度要求、分离度标准、加样协议等。
    """
    __version__ = "0.1.0"
    name = "SpecificityTestDesigner"
    func_name = "design_specificity_test"
    description = "Design a comprehensive specificity/selectivity test plan for analytical methods per ICH Q2(R1)."
    implementation_description = (
        "Generates forced degradation study design (acid/base/oxidative/thermal/photolytic/humidity), "
        "peak purity requirements, resolution criteria, spiking protocol, and acceptance criteria "
        "tailored to the analytical method type and known impurities."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Specificity", "Forced Degradation", "ICH Q2(R1)", "Analytical Chemistry", "Method Validation", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("analyte_name", "str", "N/A", "Name of the target analyte."),
        ("known_impurities", "list", "[]", "List of known impurity/degradant names to evaluate separation against."),
        ("method_type", "str", "HPLC", "Analytical method type (HPLC/UPLC/GC/UV_Vis/Titration/etc.)."),
        ("matrix_components", "list", "[]", "List of matrix components/placebo ingredients that could interfere."),
        ("degradation_types", "list", "[]",
         "Specific degradation types to include. Empty list = all 6 types (acidic/basic/oxidative/thermal/photolytic/humidity)."),
        ("include_mass_balance", "bool", "True", "Whether to include mass balance assessment in the test plan."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "String format: 'analyte_name [method_type] [impurity1;impurity2;...]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing: forced_degradation_plan(dict per type), acceptance_criteria, "
         "spiking_protocol, resolution_requirements, mass_balance_guidance, test_schedule, report"),
    ]

    examples = [
        {
            "code_input": {
                "analyte_name": "Aspirin",
                "known_impurities": ["Salicylic acid", "Acetic anhydride residue"],
                "method_type": "HPLC",
                "matrix_components": ["Starch", "Microcrystalline cellulose", "Magnesium stearate"],
                "degradation_types": [],
                "include_mass_balance": True,
            },
            "text_input": {
                "input_params": "Aspirin HPLC Salicylic_acid;Acetic_anhydride"
            },
            "output": {
                "result": {
                    "analyte_name": "Aspirin",
                    "forced_degradation_tests": {...},
                    "acceptance_criteria": {...},
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
        analyte_name: str,
        known_impurities: list = None,
        method_type: str = "HPLC",
        matrix_components: list = None,
        degradation_types: list = None,
        include_mass_balance: bool = True,
    ) -> dict:
        """
        核心逻辑：设计专属性测试方案

        Parameters:
            analyte_name: 待测物名称
            known_impurities: 已知杂质列表
            method_type: 分析方法类型
            matrix_components: 基质成分列表
            degradation_types: 指定的降解类型（空=全部）
            include_mass_balance: 是否包含质量平衡评估

        Returns:
            dict: 完整的专属性测试方案
        """
        # 处理默认值
        known_impurities = known_impurities or []
        matrix_components = matrix_components or []

        # 确定降解类型
        if not degradation_types:
            degradation_types = list(FORCED_DEGRADATION_CONDITIONS.keys())
        else:
            # 标准化输入
            deg_map = {
                "acid": "acidic", "acidic": "acidic", "base": "basic", "basic": "basic",
                "ox": "oxidative", "oxidative": "oxidative", "oxidation": "oxidative",
                "thermal": "thermal", "heat": "thermal", "thermo": "thermal",
                "photo": "photolytic", "photolytic": "photolytic", "light": "photolytic",
                "humidity": "humidity", "humid": "humidity", "rh": "humidity",
            }
            degradation_types = [deg_map.get(d.lower(), d.lower()) for d in degradation_types]
            degradation_types = [d for d in degradation_types if d in FORCED_DEGRADATION_CONDITIONS]

        # 获取方法类型的接受标准
        method_key = method_type.upper().replace("-", "_").replace(" ", "_")
        criteria = METHOD_SPECIFICITY_CRITERIA.get(method_key, METHOD_SPECIFICITY_CRITERIA["General"])

        # 构建强制降解实验计划
        fd_plan = {}
        for deg_type in degradation_types:
            if deg_type in FORCED_DEGRADATION_CONDITIONS:
                cond = FORCED_DEGRADATION_CONDITIONS[deg_type]
                fd_plan[deg_type] = {
                    "name": cond["name"],
                    "description": cond["description"],
                    "recommended_conditions": cond["typical_conditions"],
                    "target_degradation": cond["target_degradation"],
                    "neutralization": cond.get("neutralization", "N/A"),
                    "acceptance": cond.get("acceptance", criteria.get("resolution", "See general criteria")),
                    "status": "Pending",
                }

        # 加样协议（Spiking Protocol）
        spiking_protocol = self._design_spiking_protocol(analyte_name, known_impurities, matrix_components)

        # 分离度要求
        resolution_req = self._get_resolution_requirements(analyte_name, known_impurities, criteria)

        # 质量平衡指导
        mass_balance = None
        if include_mass_balance:
            mass_balance = {
                "description": (
                    "Mass balance assesses whether the sum of analyte and all degradation products "
                    "accounts for the initial amount loaded, confirming no undetected degradants."
                ),
                "formula": "%Mass Balance = (%Area_analyte + Σ%Area_degradants) × 100 / Initial_area_untreated",
                "acceptance_range": "95.0 - 105.0%",
                "procedure": [
                    "1. Prepare untreated reference solution (100% concentration).",
                    "2. Subject samples to each forced degradation condition.",
                    "3. Analyze all samples using the same method.",
                    "4. Calculate normalized % area for analyte and all degradant peaks.",
                    "5. Sum areas and compare to untreated reference.",
                    "6. If mass balance is outside 95-105%, investigate for non-eluting or poorly responding species.",
                ],
                "note": "For UV detection, ensure similar molar absorptivity or use correction factors.",
            }

        # 测试时间表
        schedule = self._generate_schedule(fd_plan, spiking_protocol, matrix_components, analyte_name)

        result = {
            "analyte_name": analyte_name,
            "method_type": method_type,
            "test_plan_summary": f"Specificity test plan for {analyte_name} by {method_type}. "
                                  f"Includes {len(fd_plan)} forced degradation tests, "
                                  f"{len(known_impurities)} known impurities, "
                                  f"{len(matrix_components)} matrix components.",
            "forced_degradation_plan": fd_plan,
            "acceptance_criteria": criteria,
            "spiking_protocol": spiking_protocol,
            "resolution_requirements": resolution_req,
            "mass_balance_assessment": mass_balance,
            "test_schedule": schedule,
            "report_text": self._generate_report(
                analyte_name, method_type, fd_plan, criteria, known_impurities, matrix_components
            ),
        }

        logger.info(f"Generated specificity test plan for {analyte_name}: {len(fd_plan)} degradation types")
        return result

    def _design_spiking_protocol(self, analyte: str, impurities: list, matrix: list) -> dict:
        """设计加样/干扰实验协议"""
        protocol = {
            "blank_interference_test": {
                "description": "Inject blank/diluent to confirm no interference at analyte retention time.",
                "samples": ["Diluent/Solvent blank", "Placebo blank (if applicable)"],
                "acceptance": "No peak > 0.05% of analyte response at analyte RT",
            },
        }

        if impurities:
            protocol["impurity_spiking"] = {
                "description": f"Spike {analyte} with known impurities at specification level or above.",
                "samples": [
                    f"{analyte} at target level + each individual impurity at LOQ/specification level",
                    f"Mixture of {analyte} + all known impurities combined",
                ],
                "acceptance": f"Resolution between {analyte} and each impurity ≥ 1.5; All impurities quantifiable",
            }

        if matrix:
            protocol["matrix_interference"] = {
                "description": "Analyze spiked matrix sample to check for interference from excipients/components.",
                "samples": [
                    "Matrix blank (placebo only)",
                    f"Matrix + {analyte} at 100% target level",
                    f"Matrix + {analyte} + all known impurities",
                ],
                "acceptance": "Recovery of analyte 98-102%; No co-elution with matrix components",
            }

        return protocol

    def _get_resolution_requirements(self, analyte: str, impurities: list, criteria: dict) -> dict:
        """生成分离度要求"""
        req = {
            "general_criteria": criteria,
            "critical_pairs": [],
        }
        for imp in impurities:
            req["critical_pairs"].append({
                "pair": f"{analyte} vs {imp}",
                "minimum_resolution": criteria.get("resolution", "≥ 1.5"),
                "priority": "Critical",
            })
        return req

    def _generate_schedule(self, fd_plan: dict, spiking: dict, matrix: list, analyte: str = "analyte") -> list:
        """生成推荐测试顺序"""
        steps = []
        step_num = 1

        # Step 1: System suitability
        steps.append({"step": step_num, "test": "System Suitability Test (SST)", "type": "Pre-test"})
        step_num += 1

        # Step 2: Blank tests
        steps.append({"step": step_num, "test": "Blank/Diluent injection", "type": "Interference Check"})
        step_num += 1

        if matrix:
            steps.append({"step": step_num, "test": "Placebo/Matrix blank injection", "type": "Interference Check"})
            step_num += 1

        # Step 3: Reference standard
        steps.append({"step": step_num, "test": f"Untreated {analyte or 'analyte'} reference solution", "type": "Reference"})
        step_num += 1

        # Step 4: Forced degradation (order matters)
        deg_order = ["acidic", "basic", "oxidative", "thermal", "photolytic", "humidity"]
        for dt in deg_order:
            if dt in fd_plan:
                steps.append({
                    "step": step_num,
                    "test": f"Forced degradation: {fd_plan[dt]['name']}",
                    "type": "Forced Degradation",
                    "details": fd_plan[dt]["recommended_conditions"][0]["condition"],
                })
                step_num += 1

        # Step 5: Spiking studies
        for spk_name, spk_data in spiking.items():
            steps.append({"step": step_num, "test": f"Spiking: {spk_name}", "type": "Spiking Study"})
            step_num += 1

        return steps

    def _generate_report(self, analyte, method, fd_plan, criteria, impurities, matrix) -> str:
        lines = [
            f"═══ SPECIFICITY TEST PLAN ═══",
            f"",
            f"Analyte: {analyte}",
            f"Method: {method}",
            f"Impurities: {len(impurities)} known | Matrix components: {len(matrix)}",
            f"",
            f"─── Acceptance Criteria ({method}) ───",
        ]
        for k, v in criteria.items():
            lines.append(f"  • {k.replace('_', ' ').title()}: {v}")

        lines.extend([
            f"",
            f"─── Forced Degradation Tests ({len(fd_plan)}) ───",
        ])
        for dt, data in fd_plan.items():
            lines.append(f"  • {data['name']}: {data['recommended_conditions'][0]['condition']}")

        lines.append(f"\nTotal tests planned: See full test_schedule in output.")
        return "\n".join(lines)

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            analyte = parts[0] if parts else "Unknown Analyte"
            method = parts[1] if len(parts) > 1 else "HPLC"
            impurities = [p.rstrip(';') for p in parts[2].split(';')] if len(parts) > 2 else []

            return self._run_base(analyte, impurities, method)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
