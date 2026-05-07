"""
稳定性考察方案设计工具 (Stability Study Planner)

基于 ICH Q1A-Q1E 指南，为原料药(DS)或制剂(DP)设计稳定性试验方案。
"""

import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ICH 稳定性条件定义（ICH Q1A, Q1B, Q1C, Q1E）
ICH_STABILITY_CONDITIONS = {
    "long_term": {
        "name": "Long-Term Stability",
        "ich_ref": "ICH Q1A(R2)",
        "conditions_by_zone": {
            "Zone_I_II": {"temp": "25°C ± 2°C", "RH": "60% RH ± 5%"},
            "Zone_III": {"temp": "30°C ± 2°C", "RH": "65% RH ± 5%"},
            "Zone_IVa": {"temp": "30°C ± 2°C", "RH": "75% RH ± 5%"},
            "Zone_IVb": {"temp": "25°C ± 2°C", "RH": "60% RH ± 5%"},  # Climate controlled
        },
        "min_duration_months": 12,
        "duration_renewal": 12,  # Annual data submission for registration
    },
    "accelerated": {
        "name": "Accelerated Stability",
        "ich_ref": "ICH Q1A(R2)",
        "conditions": {"temp": "40°C ± 2°C", "RH": "75% RH ± 5%"},
        "duration_months": 6,
        "time_points": [0, 1, 2, 3, 6],
    },
    "intermediate": {
        "name": "Intermediate Stability",
        "ich_ref": "ICH Q1A(R2)",
        "conditions": {"temp": "30°C ± 2°C", "RH": "65% RH ± 5%"},
        "duration_months": 6,
        "time_points": [0, 1, 2, 3, 6],
        "trigger": "Used when significant change occurs at accelerated conditions (40°C/75%RH).",
    },
}

# ICH Q1C 各气候区的长期试验时间点
LONG_TERM_TIME_POINTS = {
    "proposed_shelf_life_12m": [0, 3, 6, 9, 12],
    "proposed_shelf_life_18m": [0, 3, 6, 9, 12, 18],
    "proposed_shelf_life_24m": [0, 3, 6, 9, 12, 18, 24],
    "proposed_shelf_life_36m": [0, 3, 6, 9, 12, 18, 24, 36],
    "proposed_shelf_life_48m_plus": [0, 3, 6, 9, 12, 18, 24, 36, 48],
}

# 稳定性测试参数（基于剂型和产品类型）
STABILITY_TEST_PARAMETERS = {
    "drug_substance": {
        "general": [
            ("Appearance", "Visual inspection: color, physical form"),
            ("Assay/Potency", "HPLC/GC against reference standard"),
            ("Impurities", "Related substances by HPLC; degradation products"),
            ("Water Content", "Karl Fischer titration or loss on drying"),
            ("Dissolution" if False else "", ""),  # DS doesn't have dissolution
            ("Microbial Limits", "If applicable per pharmacopeia"),
            ("Particle Size Distribution", "If critical quality attribute"),
            ("Polymorph Form", "XRD/DSC if polymorphic risk exists"),
        ],
        "physical": [
            ("Solid form", "Microscopy, XRD"),
            ("Melting point", "DSC/TGA if relevant"),
            ("Flow properties", "For powder handling characteristics"),
        ],
    },
    "drug_product_tablet": [
        ("Appearance", "Color, shape, surface defects"),
        ("Assay", "Content of active ingredient(s)"),
        ("Related Substances/Degradants", "Impurity profile by HPLC"),
        ("Uniformity of Dosage Units", "Content uniformity / weight variation"),
        ("Dissolution", "USP apparatus I/II; multiple time points"),
        ("Water Content", "Karl Fischer"),
        ("Hardness/Friability", "Physical attributes"),
        ("Disintegration", "If not using dissolution test"),
        ("Microbial Limits", "For non-sterile products"),
    ],
    "drug_product_capsule": [
        ("Appearance", "Shell integrity, fill appearance"),
        ("Assay", "Content of active ingredient(s)"),
        ("Related Substances", "Degradant profile"),
        ("Uniformity of Dosage Units", "Content uniformity / weight variation"),
        ("Dissolution", "Apparatus I or II"),
        ("Water Content", "Karl Fischer"),
        ("Microbial Limits", "If applicable"),
    ],
    "drug_product_injection": [
        ("Appearance", "Color, clarity, particulate matter"),
        ("pH", "pH meter"),
        ("Assay", "Potency of active"),
        ("Related Substances", "Impurities and degradants"),
        ("Particulate Matter", "Light obscuration/microscopic per USP <788>"),
        ("Sterility", "Membrane filtration or direct inoculation"),
        ("Endotoxins/Bacterial Endotoxins", "LAL test"),
        ("Extractable Volume/Container Closure Integrity", "As applicable"),
        ("Preservative Efficacy", "If antimicrobial agent present"),
    ],
    "drug_product_general": [
        ("Appearance", "Visual inspection"),
        ("Identity", "Confirmatory test"),
        ("Assay", "Quantitative assay"),
        ("Impurities/Degradants", "Related substances"),
        ("Water Content/KF", "Moisture content"),
        ("Microbial Limits", "If applicable to dosage form"),
    ],
}


# ICH Q1E 显著变化定义
SIGNIFICANT_CHANGE_DEFINITION = {
    "drug_substance": [
        "5% change in assay from initial value",
        "Any degradant exceeding its acceptance criterion",
        "Failure to meet acceptance criteria for physical attributes (e.g., color, polymorph)",
        "Failure to meet acceptance criteria for pH (for solutions)",
        "Failure to meet acceptance criteria for dissolution (for capsules/tablets)",
    ],
    "drug_product": [
        "5% change in assay from initial value",
        "Any degradant exceeding its acceptance criterion",
        "Failure to meet the acceptance criteria for dissolution",
        "Failure to meet acceptance criteria for physical attributes (e.g., appearance, hardness)",
        "Failure to meet acceptance criteria for pH (for oral solutions)",
        "Failure to meet acceptance criteria for unit dose uniformity",
    ],
}


@ChemMCPManager.register_tool
class StabilityStudyPlanner(BaseTool):
    """
    稳定性考察方案设计工具。基于 ICH Q1A-Q1E 指南，
    为原料药或制剂设计完整的稳定性试验方案。
    """
    __version__ = "0.1.0"
    name = "StabilityStudyPlanner"
    func_name = "plan_stability_study"
    description = "Design a comprehensive stability study protocol for drug substance (DS) or drug product (DP) per ICH Q1A-Q1E guidelines."
    implementation_description = (
        "Generates ICH-compliant stability protocols with storage conditions, time points, "
        "test parameters, batch requirements, and acceptance criteria. Supports all ICH climate zones "
        "(I/II, III, IVa, IVb) and common dosage forms (tablet, capsule, injection, etc.)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Stability", "ICH Q1A", "ICH Q1B", "Regulatory", "Pharmaceutical Development", "CMC"]
    required_envs = []

    code_input_sig = [
        ("product_type", "str", "N/A", "Product type: 'DS' (Drug Substance) or 'DP' (Drug Product)."),
        ("dosage_form", "str", "General",
         "Dosage form (for DP): Tablet/Capsule/Injection/Oral_Solution/Cream/Ointment/etc.; "
         "or 'General' for DS."),
        ("ich_climate_zone", "str", "Zone_I_II",
         "ICH climate zone: Zone_I_II (temperate), Zone_III (hot/dry), Zone_IVa (hot/humid), Zone_IVb (hot/very humid)."),
        ("intended_shelf_life_months", "int", "24",
         "Proposed shelf life in months (determines long-term study duration)."),
        ("storage_condition_special", "str", "",
         "Special storage conditions if different from standard (e.g., 'Refrigerated 2-8°C' or 'Frozen -20°C')."),
        ("include_photostability", "bool", "True",
         "Whether to include photostability testing per ICH Q1B."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "String format: 'product_type [dosage_form] [climate_zone] [shelf_life_months]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing: study_protocol(dict with conditions/time_points/tests/batches), "
         "storage_conditions, test_parameters, batch_requirements, acceptance_criteria, report"),
    ]

    examples = [
        {
            "code_input": {
                "product_type": "DP",
                "dosage_form": "Tablet",
                "ich_climate_zone": "Zone_I_II",
                "intended_shelf_life_months": 24,
                "storage_condition_special": "",
                "include_photostability": True,
            },
            "text_input": {
                "input_params": "DP Tablet Zone_I_II 24"
            },
            "output": {
                "result": {
                    "product_type": "DP",
                    "study_protocol": {...},
                    "total_tests_estimated": "...",
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
        product_type: str,
        dosage_form: str = "General",
        ich_climate_zone: str = "Zone_I_II",
        intended_shelf_life_months: int = 24,
        storage_condition_special: str = "",
        include_photostability: bool = True,
    ) -> dict:
        """
        核心逻辑：设计稳定性试验方案

        Parameters:
            product_type: 产品类型 ('DS' 或 'DP')
            dosage_form: 剂型
            ich_climate_zone: ICH气候区
            intended_shelf_life_months: 预期有效期(月)
            storage_condition_special: 特殊储存条件
            include_photostability: 是否包含光稳定性

        Returns:
            dict: 完整的稳定性试验方案
        """
        # 输入验证和标准化
        pt = product_type.upper().strip()
        if pt not in ("DS", "DP"):
            raise ChemMCPError(f"product_type must be 'DS' or 'DP', got '{product_type}'.")

        zone = ich_climate_zone.strip()
        valid_zones = list(ICH_STABILITY_CONDITIONS["long_term"]["conditions_by_zone"].keys())
        if zone not in valid_zones:
            logger.warning(f"Unknown zone '{zone}', defaulting to Zone_I_II")
            zone = "Zone_I_II"

        # 获取储存条件
        lt_cond = ICH_STABILITY_CONDITIONS["long_term"]["conditions_by_zone"][zone]
        acc_cond = ICH_STABILITY_CONDITIONS["accelerated"]["conditions"]

        # 如果有特殊储存条件
        special_storage = None
        if storage_condition_special:
            special_storage = {
                "condition": storage_condition_special,
                "note": "Non-standard storage condition; requires justification.",
                "time_points": self._get_special_time_points(intended_shelf_life_months),
            }

        # 确定时间点
        shelf_key = self._shelf_life_to_key(intended_shelf_life_months)
        lt_time_points = LONG_TERM_TIME_POINTS.get(shelf_key, LONG_TERM_TIME_POINTS["proposed_shelf_life_24m"])
        acc_time_points = ICH_STABILITY_CONDITIONS["accelerated"]["time_points"]
        int_time_points = ICH_STABILITY_CONDITIONS["intermediate"]["time_points"]

        # 测试参数
        if pt == "DS":
            test_params_raw = STABILITY_TEST_PARAMETERS["drug_substance"]["general"]
            sig_change = SIGNIFICANT_CHANGE_DEFINITION["drug_substance"]
        else:
            df_key = f"drug_product_{dosage_form.lower()}"
            test_params_raw = STABILITY_TEST_PARAMETERS.get(df_key, STABILITY_TEST_PARAMETERS["drug_product_general"])
            sig_change = SIGNIFICANT_CHANGE_DEFINITION["drug_product"]

        # 清理空参数
        test_params = [(n, d) for n, d in test_params_raw if n]

        # 批次数要求
        batch_req = self._get_batch_requirements(pt, intended_shelf_life_months)

        # 构建完整方案
        protocol = {
            "long_term": {
                "condition": lt_cond,
                "time_points_months": lt_time_points,
                "num_time_points": len(lt_time_points),
                "duration_months": max(lt_time_points),
                "ich_reference": "ICH Q1A(R2)",
            },
            "accelerated": {
                "condition": acc_cond,
                "time_points_months": acc_time_points,
                "num_time_points": len(acc_time_points),
                "duration_months": max(acc_time_points),
                "ich_reference": "ICH Q1A(R2)",
            },
            "intermediate": {
                "condition": ICH_STABILITY_CONDITIONS["intermediate"]["conditions"],
                "time_points_months": int_time_points,
                "trigger_condition": ICH_STABILITY_CONDITIONS["intermediate"]["trigger"],
                "ich_reference": "ICH Q1A(R2)",
            },
        }
        if special_storage:
            protocol["special"] = special_storage

        # 光稳定性
        photo_plan = None
        if include_photostability:
            photo_plan = {
                "guideline": "ICH Q1B Photostability Testing of New Drug Substances and Products",
                "exposure": "Option 1: D65/ID65 lamp, ≥1.2×10^6 lux hours + ≥200 W·h/m² UV",
                "samples_to_test": ["Exposed sample", "Dark-wrapped control (aluminum foil)"],
                "evaluation": "Compare exposed vs control for appearance, assay, and degradants",
                "acceptance": "No significant change between exposed and dark-controlled samples",
            }

        # 估算总测试量
        num_batches = batch_req["minimum_batches"]
        lt_tests = len(test_params) * len(lt_time_points) * num_batches
        acc_tests = len(test_params) * len(acc_time_points) * num_batches
        total_estimated = lt_tests + acc_tests

        result = {
            "study_overview": {
                "product_type": pt,
                "dosage_form": dosage_form if pt == "DP" else "N/A",
                "climate_zone": zone,
                "proposed_shelf_life_months": intended_shelf_life_months,
                "regulatory_framework": "ICH Q1A(R2), Q1B, Q1C, Q1E",
            },
            "storage_conditions": {
                "long_term": lt_cond,
                "accelerated": acc_cond,
                "intermediate_available": True,
                "special": storage_condition_special if storage_condition_special else None,
            },
            "protocol": protocol,
            "test_parameters": [{"test": n, "method": d} for n, d in test_params],
            "batch_requirements": batch_req,
            "significant_change_definition": sig_change,
            "photostability_plan": photo_plan,
            "estimated_testing_load": {
                "long_term_tests": lt_tests,
                "accelerated_tests": acc_tests,
                "total_test_analyses": total_estimated,
                "batches": num_batches,
            },
            "report_text": self._generate_report(
                pt, dosage_form, zone, intended_shelf_life_months,
                lt_cond, acc_cond, lt_time_points, acc_time_points,
                test_params, batch_req, total_estimated
            ),
        }

        logger.info(f"Stability study plan: {pt}/{dosage_form}, {zone}, {intended_shelf_life_months}mo shelf life")
        return result

    def _shelf_life_to_key(self, months: int) -> str:
        """将有效期映射到时间点key"""
        if months <= 12:
            return "proposed_shelf_life_12m"
        elif months <= 18:
            return "proposed_shelf_life_18m"
        elif months <= 24:
            return "proposed_shelf_life_24m"
        elif months <= 36:
            return "proposed_shelf_life_36m"
        else:
            return "proposed_shelf_life_48m_plus"

    def _get_special_time_points(self, months: int) -> list:
        """特殊储存条件的时间点"""
        if "refrigerate" in storage_condition_special.lower() or "2-8" in storage_condition_special:
            return [0, 3, 6, 9, 12, 18, 24]
        elif "frozen" in storage_condition_special.lower() or "-20" in storage_condition_special:
            return [0, 3, 6, 9, 12, 18, 24]
        return [0, 3, 6, 9, 12]

    def _get_batch_requirements(self, pt: str, shelf_life: int) -> dict:
        """获取批次数要求"""
        min_batches = 3  # ICH minimum
        if shelf_life > 24:
            batches_note = (
                f"Minimum {min_batches} batches (pilot/commercial scale). "
                f"For shelf life >{shelf_life} months, consider additional data packages upon annual renewal."
            )
        else:
            batches_note = (
                f"Minimum {min_batches} batches, independent production batches, "
                f"same formulation and process as commercial scale."
            )

        return {
            "minimum_batches": min_batches,
            "batch_scale": "Pilot scale or commercial scale (representative of production)",
            "batch_selection": "Independent primary batches from same manufacturing process",
            "additional_notes": batches_note,
        }

    def _generate_report(self, pt, df, zone, shelf_lt, lt_cond, acc_cond,
                         lt_tp, acc_tp, tests, batch_req, total_tests) -> str:
        lines = [
            f"═══ STABILITY STUDY PROTOCOL ═══",
            f"",
            f"Product Type: {'Drug Substance (DS)' if pt == 'DS' else 'Drug Product (DP)'}",
            f"Dosage Form: {df}",
            f"Climate Zone: {zone}",
            f"Proposed Shelf Life: {shelf_lt} months",
            f"",
            f"─── Storage Conditions ───",
            f"  Long-term: {lt_cond['temp']}/{lt_cond['RH']}",
            f"  Accelerated: {acc_cond['temp']}/{acc_cond['RH']}",
            f"",
            f"─── Time Points ───",
            f"  Long-term (months): {lt_tp}",
            f"  Accelerated (months): {acc_tp}",
            f"",
            f"─── Test Parameters ({len(tests)} tests) ───",
        ]
        for t_name, t_method in tests[:10]:
            lines.append(f"  • {t_name}: {t_method}")
        if len(tests) > 10:
            lines.append(f"  ... and {len(tests)-10} more tests")

        lines.extend([
            f"",
            f"─── Batch Requirements ───",
            f"  Minimum batches: {batch_req['minimum_batches']}",
            f"",
            f"─── Estimated Workload ───",
            f"  Total analyses (approx): {total_tests}",
        ])
        return "\n".join(lines)

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入"""
        try:
            parts = input_params.strip().split()
            pt = parts[0].upper() if parts else "DP"
            df = parts[1] if len(parts) > 1 else "Tablet"
            zone = parts[2] if len(parts) > 2 else "Zone_I_II"
            shelf = int(parts[3]) if len(parts) > 3 else 24

            return self._run_base(pt, df, zone, shelf)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
