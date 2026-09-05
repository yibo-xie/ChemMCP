"""
Test suite for ChemMCP Tools #391-400
Analytical Chemistry QA/QC & Titration Tools
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from chemmcp.tools import (
    BackTitrationSolver,
    TitrantStandardization,
    BufferCapacityCalculator,
    MethodValidationChecklist,
    LinearityRangeValidator,
    RobustnessDoeDesigner,
    SpecificityTestDesigner,
    StabilityStudyPlanner,
    CMCDocumentationHelper,
    AuditTrailReviewer,
)


def separator(title):
    print(f"\n{'='*60}")
    print(f"  Testing: {title}")
    print(f"{'='*60}")


def print_result(label, result):
    print(f"\n  ✅ {label}")
    if isinstance(result, dict):
        for k, v in result.items():
            if k not in ("detailed_calculation", "report_text"):
                print(f"     {k}: {v}")
    else:
        print(f"     Result: {result}")


passed = 0
failed = 0


def test_tool(name, test_fn):
    global passed, failed
    try:
        test_fn()
        passed += 1
        print(f"\n  🎉 {name}: ALL TESTS PASSED")
    except Exception as e:
        failed += 1
        print(f"\n  💥 {name}: FAILED - {e}")
        import traceback
        traceback.print_exc()


# ============================================================
# Test 391: BackTitrationSolver
# ============================================================
def test_back_titration_solver():
    separator("391 BackTitrationSolver (返滴定计算)")
    tool = BackTitrationSolver()

    result = tool.run_code(
        analyte_mass=0.5000, analyte_mw=100.09, titrant_conc=0.1000,
        excess_titrant_vol=50.00, back_titrant_conc=0.1050,
        back_titrant_vol=12.50, stoich_ratio=1.0,
    )
    assert isinstance(result, dict)
    assert "analyte_purity" in result
    assert 0 < result["analyte_purity"] <= 100
    print_result("Basic back titration", {
        "purity": f"{result['analyte_purity']}%",
        "mass": f"{result['mass_analyte']}g",
        "moles": f"{result['moles_analyte']} mol",
    })

    result_text = tool.run_text("0.5000 100.09 0.1000 50.00 0.1050 12.50 1.0")
    assert abs(result_text["analyte_purity"] - result["analyte_purity"]) < 0.01
    print_result("Text interface", {"purity": f"{result_text['analyte_purity']}%"})

    try:
        tool.run_code(analyte_mass=-1.0, analyte_mw=100.0, titrant_conc=0.1,
                      excess_titrant_vol=50.0, back_titrant_conc=0.1, back_titrant_vol=10.0)
        assert False
    except Exception:
        print_result("Error handling (negative mass)", "Correctly raised error")

    try:
        tool.run_code(analyte_mass=0.5, analyte_mw=100.09, titrant_conc=0.1,
                      excess_titrant_vol=5.0, back_titrant_conc=0.5, back_titrant_vol=20.0)
        assert False
    except Exception as e:
        assert "negative" in str(e).lower() or "exceed" in str(e).lower()
        print_result("Error handling (impossible volumes)", f"Caught: {str(e)[:60]}")

    result_pure = tool.run_code(
        analyte_mass=0.5000, analyte_mw=100.09, titrant_conc=0.1000,
        excess_titrant_vol=50.00, back_titrant_conc=0.1050,
        back_titrant_vol=47.62, stoich_ratio=1.0,
    )
    print_result("Near-pure sample", {"purity": f"{result_pure['analyte_purity']}%"})


# ============================================================
# Test 392: TitrantStandardization
# ============================================================
def test_titrant_standardization():
    separator("392 TitrantStandardization (滴定剂标定计算)")
    tool = TitrantStandardization()

    result = tool.run_code(
        primary_std_mass=0.2042, primary_std_mw=105.99,
        titrant_vol=38.20, stoich_ratio=2.0,
    )
    conc = result["titrant_concentration"]
    assert 0.09 < conc < 0.12
    print_result("HCl standardization with Na2CO3", {
        "concentration": f"{conc:.4f} M",
        "moles_std": f"{result['moles_standard']:.6f} mol",
    })

    result_stats = tool.run_code(
        primary_std_mass=0.2042, primary_std_mw=105.99,
        titrant_vol=38.20, stoich_ratio=2.0,
        replicate_vols=[38.20, 38.18, 38.22, 38.15],
        confidence_level=0.95,
    )
    stats = result_stats["statistics"]
    assert stats["n_replicates"] == 4
    assert "rsd_percent" in stats
    print_result("With replicates", {
        "n": stats["n_replicates"],
        "mean_C": f"{stats['mean_concentration']:.6f} M",
        "RSD": f"{stats['rsd_percent']:.4f}%",
    })

    result_text = tool.run_text("0.2042 105.99 38.20 2.0")
    assert abs(result_text["titrant_concentration"] - conc) < 0.001
    print_result("Text interface", {"C": f"{result_text['titrant_concentration']:.4f} M"})

    result_khp = tool.run_code(
        primary_std_mass=0.5100, primary_std_mw=204.22,
        titrant_vol=24.85, stoich_ratio=1.0,
    )
    print_result("KHP -> NaOH (1:1)", {"C_NaOH": f"{result_khp['titrant_concentration']:.4f} M"})


# ============================================================
# Test 393: BufferCapacityCalculator
# ============================================================
def test_buffer_capacity_calculator():
    separator("393 BufferCapacityCalculator (缓冲容量计算)")
    tool = BufferCapacityCalculator()

    result = tool.run_code(
        total_buffer_conc=0.10, ph_initial=4.76, pka=4.76,
        delta_ph=0.1, calculation_mode="exact",
    )
    beta = result["beta"]
    expected_beta_max = 2.303 * 0.10 * 0.25
    assert abs(beta - expected_beta_max) < 0.001
    print_result("Acetate buffer at pKa", {
        "β": f"{beta:.4f} eq/(L·pH)",
        "β_max": f"{result['max_beta']:.4f}",
        "range": result["buffer_range"],
        "[HA]": f"{result['ha_fraction']*100:.1f}%",
        "[A-]": f"{result['a_minus_fraction']*100:.1f}%",
    })

    result_phos = tool.run_code(
        total_buffer_conc=0.20, ph_initial=7.40, pka=7.21,
        delta_ph=0.05, calculation_mode="exact",
    )
    print_result("Phosphate buffer pH 7.40", {
        "β": f"{result_phos['beta']:.4f}",
        "[H2PO4-]": f"{result_phos['ha_fraction']*100:.1f}%",
        "[HPO4^2-]": f"{result_phos['a_minus_fraction']*100:.1f}%",
    })

    result_off = tool.run_code(total_buffer_conc=0.10, ph_initial=3.76, pka=4.76)
    assert result_off["beta"] < result["beta"]
    print_result("Away from pKa", {"β": f"{result_off['beta']:.4f} (vs {result['beta']:.4f} at pKa)"})

    result_approx = tool.run_code(
        total_buffer_conc=0.10, ph_initial=4.76, pka=4.76,
        delta_ph=0.1, calculation_mode="approximate",
    )
    assert "approximate_beta" in result_approx
    print_result("Approximate mode", {"β_exact": f"{result['beta']:.4f}", "β_approx": f"{result_approx['approximate_beta']:.4f}"})

    result_text = tool.run_text("0.10 4.76 4.76 0.1 exact")
    assert abs(result_text["beta"] - beta) < 0.0001
    print_result("Text interface", {"β": f"{result_text['beta']:.4f}"})


# ============================================================
# Test 394: MethodValidationChecklist
# ============================================================
def test_method_validation_checklist():
    separator("394 MethodValidationChecklist (ICH Q2(R1))")
    tool = MethodValidationChecklist()

    result = tool.run_code(method_type="HPLC", analyte_name="Ibuprofen", matrix_type="Tablet")
    params = result["validation_parameters"]
    expected_params = ["Specificity", "Linearity", "Range", "Accuracy", "Precision",
                       "Detection Limit (LOD)", "Quantitation Limit (LOQ)", "Robustness"]
    for ep in expected_params:
        assert ep in params, f"Missing: {ep}"
    for pname, pdata in params.items():
        assert "description" in pdata and "acceptance_criteria" in pdata

    print_result("HPLC validation parameters", {
        "method": result["method_info"]["method_type"],
        "num_parameters": len(params),
        "total_tests": result["total_recommended_tests"],
        "params_list": list(params.keys()),
    })

    precision = params.get("Precision", {})
    assert "sub_categories" in precision
    print_result("Precision sub-categories", list(precision["sub_categories"].keys()))

    schedule = result["test_schedule"]
    assert schedule[0] == "Specificity"
    print_result("Test schedule order", schedule)

    refs = result["regulatory_references"]
    assert len(refs) > 0
    print_result("Regulatory references", refs)

    result_uv = tool.run_code(method_type="UV_Vis", analyte_name="Paracetamol solution")
    print_result("UV-Vis method", {"num_params": len(result_uv["validation_parameters"])})

    result_text = tool.run_text("GC Ethanol Assay Plasma FDA")
    assert "validation_parameters" in result_text
    print_result("Text interface (GC)", {"num_params": len(result_text["validation_parameters"])})


# ============================================================
# Test 395: LinearityRangeValidator
# ============================================================
def test_linearity_range_validator():
    separator("395 LinearityRangeValidator (线性范围验证)")
    tool = LinearityRangeValidator()

    concentrations = [50.0, 75.0, 100.0, 125.0, 150.0]
    responses = [2000*c + 100 + (i*5 - 10) for i, c in enumerate(concentrations)]

    result = tool.run_code(
        concentrations=concentrations, responses=responses,
        target_correlation=0.999, confidence_level=0.95, unit="μg/mL",
    )
    r = result["correlation"]["r"]
    assert r > 0.999
    print_result("Near-perfect linearity", {
        "r": f"{r:.6f}", "R²": f"{result['correlation']['r_squared']:.6f}",
        "slope": f"{result['regression']['slope']:.2f}",
        "intercept": f"{result['regression']['intercept']:.2f}",
        "assessment": result["linearity_assessment"],
        "Sy/x": f"{result['precision']['standard_error']:.4f}",
    })

    back_calc = result["back_calculations"]
    assert len(back_calc) == 5
    print_result("Back-calculations", [
        f"nom={bc['nominal_conc']}, calc={bc['back_calculated_conc']}, bias={bc['bias_percent']}%"
        for bc in back_calc[:3]
    ] + ["..."])

    residuals = result["residual_analysis"]
    print_result("Residuals", {
        "max_abs": residuals["max_absolute_residual"],
        "max_rel_%": residuals["max_relative_residual_percent"],
    })

    outliers = result["outliers_detected"]
    print_result("Outliers", {"found": len(outliers), "details": outliers if outliers else "None ✅"})

    assert "report_text" in result and len(result["report_text"]) > 50
    print_result("Report", f"{len(result['report_text'])} chars")

    result_text = tool.run_text("50.0,10100;75.0,15100;100.0,20100;125.0,25100;150.0,30100 0.999")
    assert result_text["correlation"]["r"] > 0.999
    print_result("Text interface", {"r": f"{result_text['correlation']['r']:.6f}"})


# ============================================================
# Test 396: RobustnessDoeDesigner
# ============================================================
def test_robustness_doe_designer():
    separator("396 RobustnessDoeDesigner (Plackett-Burman)")
    tool = RobustnessDoeDesigner()

    factors = [
        {"name": "Organic_Ratio", "low": 38.0, "high": 42.0, "unit": "%"},
        {"name": "Buffer_pH", "low": 2.8, "high": 3.2, "unit": ""},
        {"name": "Flow_Rate", "low": 0.9, "high": 1.1, "unit": "mL/min"},
        {"name": "Column_Temp", "low": 28, "high": 32, "unit": "°C"},
        {"name": "Wavelength", "low": 252, "high": 258, "unit": "nm"},
        {"name": "Injection_Vol", "low": 8, "high": 12, "unit": "μL"},
        {"name": "Gradient_Start", "low": -2, "high": 2, "unit": "%"},
    ]

    result = tool.run_code(
        factors=factors, num_runs=0, response_name="Resolution_RS1",
        include_center_point=True, center_replicates=3,
    )

    di = result["design_info"]
    assert di["num_factors"] == 7
    assert di["pb_runs"] >= 7
    assert di["total_runs"] == di["pb_runs"] + 3

    matrix = result["design_matrix"]
    assert len(matrix) == di["pb_runs"]
    print_result("PB Design summary", {
        "type": di["design_type"], "factors": di["num_factors"],
        "PB_runs": di["pb_runs"], "center_points": di["center_point_runs"],
        "total_runs": di["total_runs"],
    })
    print_result("Matrix (first 3)", matrix[:3])

    fs = result["factor_summary"]
    assert len(fs) == 7
    print_result("Factor summary (first 3)", fs[:3])

    ee = result["effect_estimation"]
    print_result("Effect formula", ee["general_formula"])

    aliasing = result["aliasing_structure"]
    print_result("Aliasing", aliasing["note"])

    ai = result["analysis_instructions"]
    assert len(ai) > 0
    print_result("Analysis steps", f"{len(ai)} steps")

    cps = result["center_points"]
    assert len(cps) == 3
    print_result("Center points", f"{len(cps)} CPs")

    result_min = tool.run_code(
        factors=[{"name": "pH", "low": 3.0, "high": 5.0},
                 {"name": "Temp", "low": 25, "high": 35, "unit": "°C"}],
        num_runs=4,
    )
    assert result_min["design_info"]["num_factors"] == 2
    print_result("Min 2 factors N=4", {"runs": result_min["design_info"]["total_runs"]})

    result_text = tool.run_text("pH,3.0,5.0;Temp,25,35,C 4")
    assert "design_matrix" in result_text
    print_result("Text interface", {"factors": result_text["design_info"]["num_factors"]})


# ============================================================
# Test 397: SpecificityTestDesigner
# ============================================================
def test_specificity_test_designer():
    separator("397 SpecificityTestDesigner (专属性测试方案设计)")
    tool = SpecificityTestDesigner()

    result = tool.run_code(
        analyte_name="Aspirin",
        known_impurities=["Salicylic acid", "Acetylsalicyloyl salicylic acid"],
        method_type="HPLC",
        matrix_components=["Starch", "Microcrystalline cellulose", "Magnesium stearate"],
        degradation_types=[], include_mass_balance=True,
    )

    fd_plan = result["forced_degradation_plan"]
    expected_deg = ["acidic", "basic", "oxidative", "thermal", "photolytic", "humidity"]
    for deg in expected_deg:
        assert deg in fd_plan, f"Missing: {deg}"
        assert len(fd_plan[deg]["recommended_conditions"]) > 0

    print_result("Forced degradation tests", {
        "analyte": result["analyte_name"], "num_tests": len(fd_plan),
        "types": list(fd_plan.keys()),
    })

    acid_detail = fd_plan["acidic"]
    print_result("Acidic degradation", {
        "conditions": [c["condition"] for c in acid_detail["recommended_conditions"]],
        "target": acid_detail["target_degradation"],
    })

    criteria = result["acceptance_criteria"]
    assert "resolution" in criteria
    print_result("Acceptance criteria (HPLC)", criteria)

    spiking = result["spiking_protocol"]
    assert len(spiking) > 0
    print_result("Spiking protocols", list(spiking.keys()))

    mb = result["mass_balance_assessment"]
    assert mb is not None and "procedure" in mb
    print_result("Mass balance", {"range": mb["acceptance_range"], "steps": len(mb["procedure"])})

    schedule = result["test_schedule"]
    assert len(schedule) > 0
    print_result("Schedule", f"{len(schedule)} steps")

    result_simple = tool.run_code(
        analyte_name="Glucose", method_type="UV_Vis",
        degradation_types=["oxidative"], include_mass_balance=False,
    )
    print_result("Simple UV-Vis", {"tests": len(result_simple["forced_degradation_plan"])})

    result_text = tool.run_text("Ibuprofen HPLC ImpA;ImpB")
    assert "forced_degradation_plan" in result_text
    print_result("Text interface", {"tests": len(result_text["forced_degradation_plan"])})


# ============================================================
# Test 398: StabilityStudyPlanner
# ============================================================
def test_stability_study_planner():
    separator("398 StabilityStudyPlanner (ICH Q1A-Q1E)")
    tool = StabilityStudyPlanner()

    result = tool.run_code(
        product_type="DP", dosage_form="Tablet",
        ich_climate_zone="Zone_I_II", intended_shelf_life_months=24,
        include_photostability=True,
    )

    overview = result["study_overview"]
    assert overview["product_type"] == "DP"

    sc = result["storage_conditions"]
    lt_cond = sc["long_term"]
    acc_cond = sc["accelerated"]
    print_result("Storage conditions", {"LT": lt_cond, "ACC": acc_cond})

    protocol = result["protocol"]
    print_result("Time points", {
        "LT (mo)": protocol["long_term"]["time_points_months"],
        "ACC (mo)": protocol["accelerated"]["time_points_months"],
    })

    tp = result["test_parameters"]
    assert len(tp) > 0
    print_result("Test params ({n})".format(n=len(tp)), [t["test"] for t in tp])

    br = result["batch_requirements"]
    assert br["minimum_batches"] >= 3
    print_result("Batch requirements", br)

    sc_def = result["significant_change_definition"]
    assert len(sc_def) > 0
    print_result("Significant change (DP)", sc_def[:3])

    photo = result["photostability_plan"]
    assert photo is not None
    print_result("Photostability", {"guideline": photo["guideline"]})

    wl = result["estimated_testing_load"]
    print_result("Workload estimate", wl)

    result_ds = tool.run_code(product_type="DS", dosage_form="General",
                              ich_climate_zone="Zone_I_II", intended_shelf_life_months=36)
    ds_tp = result_ds["protocol"]["long_term"]["time_points_months"]
    assert max(ds_tp) >= 36
    print_result("DS 36-month", {"timepoints": ds_tp, "tests": len(result_ds["test_parameters"])})

    result_z3 = tool.run_code(product_type="DP", dosage_form="Capsule",
                              ich_climate_zone="Zone_III", intended_shelf_life_months=24)
    print_result("Zone III", result_z3["storage_conditions"]["long_term"])

    result_cold = tool.run_code(product_type="DP", dosage_form="Injection",
                                ich_climate_zone="Zone_I_II", intended_shelf_life_months=18,
                                storage_condition_special="Refrigerated 2-8°C")
    assert result_cold["storage_conditions"]["special"] is not None
    print_result("Refrigerated special", result_cold["storage_conditions"]["special"])

    result_text = tool.run_text("DP Tablet Zone_I_II 24")
    assert "protocol" in result_text
    print_result("Text interface", {"shelf_life": result_text["study_overview"]["proposed_shelf_life_months"]})


# ============================================================
# Test 399: CMCDocumentationHelper
# ============================================================
def test_cmc_documentation_helper():
    separator("399 CMCDocumentationHelper (CMC ICH M4Q)")
    tool = CMCDocumentationHelper()

    result = tool.run_code(
        document_section="all", product_type="DS",
        drug_name="NovelDrugAPI", include_regulatory_refs=True, detail_level="standard",
    )

    structure = result["document_structure"]
    expected_sections = ["S.1_General_Information", "S.2_Manufacture", "S.3_Characterization",
                         "S.4_Control_of_Drug_Substance", "S.5_Reference_Standards_or_Materials",
                         "S.6_Container_Closure_System", "S.7_Stability"]
    for es in expected_sections:
        assert es in structure, f"Missing: {es}"

    print_result("CMC DS full outline", {
        "drug": result["document_info"]["drug_name"],
        "sections": len(structure),
        "names": list(structure.keys()),
        "items": result["total_content_items"],
    })

    s4 = structure.get("S.4_Control_of_Drug_Substance", {})
    if s4:
        print_result("S.4 subsections", list(s4.get("subsections", {}).keys()))

    result_s4 = tool.run_code(document_section="S.4", product_type="DS",
                               drug_name="ExampleAPI", detail_level="detailed")
    s4_detail = result_s4["document_structure"].get("S.4_Control_of_Drug_Substance", {})
    assert "title" in s4_detail
    print_result("S.4 Detailed", {"title": s4_detail["title"],
                                    "subs": list(s4_detail.get("subsections", {}).keys())})

    result_dp = tool.run_code(document_section="all", product_type="DP",
                               drug_name="WonderTablet", detail_level="outline")
    dp_sections = list(result_dp["document_structure"].keys())
    assert any("P." in s for s in dp_sections)
    print_result("CMC DP outline", {"sections": dp_sections, "items": result_dp["total_content_items"]})

    refs = result["regulatory_references"]
    assert len(refs) > 0
    print_result("Regulatory refs", {"categories": list(refs.keys()),
                                     "total": sum(len(v) for v in refs.values())})

    checklist = result["writing_checklist"]
    assert len(checklist) > 0
    print_result("Checklist", f"{len(checklist)} items")
    for item in checklist[:5]:
        print(f"     {item}")

    result_text = tool.run_text("S.2 DS MyDrug detailed")
    assert "document_structure" in result_text
    print_result("Text interface", {"sections": list(result_text["document_structure"].keys())})


# ============================================================
# Test 400: AuditTrailReviewer
# ============================================================
def test_audit_trail_reviewer():
    separator("400 AuditTrailReviewer (ALCOA+ 审计追踪检查)")
    tool = AuditTrailReviewer()

    audit_records = [
        {"timestamp": "2026-01-15T09:30:00Z", "user_id": "analyst_01", "action": "CREATE",
         "field": "Assay_Result_Batch_001", "old_value": None, "new_value": "99.5%",
         "reason": "Initial data entry after HPLC analysis"},
        {"timestamp": "2026-01-15T09:31:00Z", "user_id": "analyst_01", "action": "CREATE",
         "field": "Related_Substances_Batch_001", "old_value": None, "new_value": "Total impurities: 0.08%",
         "reason": "Initial entry of impurity profile"},
        {"timestamp": "2026-01-15T10:15:00Z", "user_id": "analyst_01", "action": "UPDATE",
         "field": "Assay_Result_Batch_001", "old_value": "99.5%", "new_value": "99.8%",
         "reason": "Correction: integration parameter adjusted"},
        {"timestamp": "2026-01-15T11:00:00Z", "user_id": "analyst_01", "action": "UPDATE",
         "field": "System_Suitability", "old_value": "Pending", "new_value": "PASS",
         "reason": "SST results entered"},
        {"timestamp": "2026-01-15T14:30:00Z", "user_id": "reviewer_02", "action": "REVIEW",
         "field": "Batch_001_Release", "old_value": None, "new_value": "APPROVED",
         "reason": "Second-person review complete"},
        {"timestamp": "2026-01-16T08:00:00Z", "user_id": "admin_sys", "action": "SYSTEM",
         "field": "Backup_Completed", "old_value": None, "new_value": "Backup_20260116_0800",
         "reason": "Scheduled automated nightly backup"},
        {"timestamp": "2026-01-17T13:45:00Z", "user_id": "qa_auditor", "action": "REVIEW",
         "field": "Audit_Trail_Review_Q1", "old_value": None, "new_value": "COMPLIANT",
         "reason": "Quarterly ALCOA+ review completed"},
    ]

    system_context = {"system_name": "LIMS v3.2", "gxp_type": "GMP",
                      "retention_years": 10, "has_backup": True}

    result = tool.run_code(audit_records=audit_records, check_level="comprehensive",
                            system_context=system_context)

    summary = result["review_summary"]
    assert summary["total_records_reviewed"] == len(audit_records)
    score = summary["overall_compliance_percent"]
    assert 0 <= score <= 100
    print_result("Overall compliance", {
        "records": summary["total_records_reviewed"], "score": f"{score:.1f}%",
        "status": summary["compliance_status"],
        "principles": summary["principles_evaluated"],
    })

    ps = result["principle_scores"]
    assert len(ps) == 9
    print_result("Principle scores (9 ALCOA+)", {
        name: f"{data['score_percent']}% ({data['items_passed']}/{data['items_total']})"
        for name, data in ps.items()
    })

    issues = result["issues_found"]
    print_result("Issues", {"count": len(issues),
                           "details": [(i["severity"], i.get("full_principle","")) for i in issues[:5]] if issues else "None ✅"})

    risk = result["risk_assessment"]
    assert "risk_level" in risk
    print_result("Risk assessment", {"level": risk["risk_level"],
                                      "summary": risk["issue_summary"],
                                      "response_time": risk["recommended_response_time"]})

    recs = result["recommendations"]
    assert len(recs) > 0
    print_result("Recommendations", [f"[{r['priority']}] {r['recommendation'][:60]}..." for r in recs[:3]])

    assert len(result["report_text"]) > 100
    print_result("Report", f"{len(result['report_text'])} chars")

    result_basic = tool.run_code(audit_records=audit_records, check_level="basic",
                                  system_context=system_context)
    basic_ps = result_basic["principle_scores"]
    assert len(basic_ps) == 5
    print_result("Basic level (5 principles)", list(basic_ps.keys()))

    result_no_backup = tool.run_code(
        audit_records=audit_records[:3], check_level="comprehensive",
        system_context={"gxp_type": "GLP", "retention_years": 5, "has_backup": False},
    )
    no_backup_issues = [i for i in result_no_backup["issues_found"] if i.get("principle") == "Enduring"]
    print_result("No backup flagged", {"enduring_issues": len(no_backup_issues),
                                         "score": f"{result_no_backup['review_summary']['overall_compliance_percent']}%"})

    try:
        tool.run_code(audit_records=[], check_level="comprehensive")
        assert False
    except Exception as e:
        assert "empty" in str(e).lower()
        print_result("Error (empty records)", f"Caught: {e}")

    result_text = tool.run_text("comprehensive")
    assert "principle_scores" in result_text
    print_result("Text interface (demo)", {
        "score": f"{result_text['review_summary']['overall_compliance_percent']}%",
        "status": result_text["review_summary"]["compliance_status"],
    })


# ============================================================
# Run all tests
# ============================================================
if __name__ == "__main__":
    print("\n" + "="*70)
    print("  ChemMCP Tools #391-400 Test Suite")
    print("  Analytical Chemistry QA/QC & Titration Tools")
    print("="*70)

    test_tool("#391 BackTitrationSolver", test_back_titration_solver)
    test_tool("#392 TitrantStandardization", test_titrant_standardization)
    test_tool("#393 BufferCapacityCalculator", test_buffer_capacity_calculator)
    test_tool("#394 MethodValidationChecklist", test_method_validation_checklist)
    test_tool("#395 LinearityRangeValidator", test_linearity_range_validator)
    test_tool("#396 RobustnessDoeDesigner", test_robustness_doe_designer)
    test_tool("#397 SpecificityTestDesigner", test_specificity_test_designer)
    test_tool("#398 StabilityStudyPlanner", test_stability_study_planner)
    test_tool("#399 CMCDocumentationHelper", test_cmc_documentation_helper)
    test_tool("#400 AuditTrailReviewer", test_audit_trail_reviewer)

    print("\n" + "="*70)
    print(f"  RESULTS: {passed} PASSED | {failed} FAILED | TOTAL: {passed+failed}")
    print("="*70)

    if failed > 0:
        sys.exit(1)
    else:
        print("\n  🎉 ALL 10 TOOLS TESTS PASSED! 🎉\n")
        sys.exit(0)
