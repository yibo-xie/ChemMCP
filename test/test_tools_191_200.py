"""
Test suite for ChemMCP tools #191-200 (Advanced Organic Chemistry & Physical Organic Tools).
Run: python -m pytest test/test_tools_191_200.py -v
Or: python test/test_tools_191_200.py
"""

import sys
sys.path.insert(0, "src")


def _run_and_check(tool_cls, kwargs, name, min_length=100, must_contain=None, must_not_contain=None):
    """Helper: run a tool and check output quality."""
    tool = tool_cls()
    result = tool.run_code(**kwargs)
    assert isinstance(result, str), f"{name}: expected str, got {type(result)}"
    assert len(result) >= min_length, f"{name}: result too short ({len(result)} < {min_length})"
    if must_contain:
        for phrase in must_contain:
            assert phrase.lower() in result.lower(), f"{name}: missing '{phrase}'"
    if must_not_contain:
        for phrase in must_not_contain:
            assert phrase.lower() not in result.lower(), f"{name}: should not contain '{phrase}'"
    print(f"  [PASS] {name} (len={len(result)})")
    return result


def test_191_ring_formation_strategy():
    """#191 RingFormationStrategy - 成环反应策略"""
    from chemmcp.tools import RingFormationStrategy
    return _run_and_check(
        RingFormationStrategy,
        dict(target_ring_size=6, starting_material_hint="diacid", constraints="stereocontrol"),
        "191 RingFormationStrategy",
        min_length=500,
        must_contain=["6", "ring", "strategy"],
    )


def test_192_asymmetric_synthesis_guide():
    """#192 AsymmetricSynthesisGuide - 不对称合成方法选择"""
    from chemmcp.tools import AsymmetricSynthesisGuide
    return _run_and_check(
        AsymmetricSynthesisGuide,
        dict(target_type="chiral alcohol", substrate_hint="ketone", constraints="high ee (>99%)"),
        "192 AsymmetricSynthesisGuide",
        min_length=500,
        must_contain=["chiral", "alcohol", "asymmet"],
    )


def test_193_total_synthesis_planner():
    """#193 TotalSynthesisPlanner - 多步合成路线规划"""
    from chemmcp.tools import TotalSynthesisPlanner
    return _run_and_check(
        TotalSynthesisPlanner,
        dict(target_molecule="aspirin", complexity_level="simple"),
        "193 TotalSynthesisPlanner",
        min_length=300,
        must_contain=["aspirin", "synth"],
    )


def test_194_pka_predictor():
    """#194 PkaPredictor - 预测有机酸的pKa值 (acetic acid ~4.76)"""
    from chemmcp.tools import PkaPredictor
    result = _run_and_check(
        PkaPredictor,
        dict(molecule="acetic acid", solvent="water"),
        "194 PkaPredictor",
        min_length=100,
        must_contain=["4.7", "pka", "acetic"],
    )
    # Verify pKa value is in reasonable range
    assert "4.76" in result or "4.8" in result or "4.7" in result, f"Expected pKa ~4.76 for acetic acid"
    return result


def test_195_hammett_sigma_lookup():
    """#195 HammettSigmaLookup - 查询Hammett σ常数 (p-NO2 should be strongly EWG)"""
    from chemmcp.tools import HammettSigmaLookup
    return _run_and_check(
        HammettSigmaLookup,
        dict(substituent="p-NO2", include_swain_lupton=True),
        "195 HammettSigmaLookup",
        min_length=300,
        must_contain=["no2", "hammett", "0.7"],
    )


def test_196_resonance_structure_generator():
    """#196 ResonanceStructureGenerator - 生成共振结构式 (carboxylate has 2+ forms)"""
    from chemmcp.tools import ResonanceStructureGenerator
    return _run_and_check(
        ResonanceStructureGenerator,
        dict(molecule="carboxylate", show_curved_arrows=True),
        "196 ResonanceStructureGenerator",
        min_length=300,
        must_contain=["resonance", "carboxylat"],
    )


def test_197_inductive_effect_analyzer():
    """#197 InductiveEffectAnalyzer - 诱导效应 (Cl is electron-withdrawing)"""
    from chemmcp.tools import InductiveEffectAnalyzer
    return _run_and_check(
        InductiveEffectAnalyzer,
        dict(molecule="chloroacetic acid", focus_property="acidity"),
        "197 InductiveEffectAnalyzer",
        min_length=200,
        must_contain=["induct", "chloro", "acid"],
    )


def test_198_hyperconjugation_explainer():
    """#198 HyperconjugationExplainer - 超共轭效应 (tert-butyl cation stable)"""
    from chemmcp.tools import HyperconjugationExplainer
    return _run_and_check(
        HyperconjugationExplainer,
        dict(molecule="tert-butyl cation", question="why is tertiary more stable?"),
        "198 HyperconjugationExplainer",
        min_length=200,
        must_contain=["hyperconjug", "tert-butyl", "stabil"],
    )


def test_199_steric_effect_analyzer():
    """#199 StericEffectAnalyzer - 位阻效应"""
    from chemmcp.tools import StericEffectAnalyzer
    return _run_and_check(
        StericEffectAnalyzer,
        dict(molecule="tert-butylcyclohexane", reaction_context="conformation"),
        "199 StericEffectAnalyzer",
        min_length=200,
        must_contain=["steric", "tert-butyl"],
    )


def test_200_conformational_analyzer():
    """#200 ConformationalAnalyzer - 构象稳定性分析 (Newman投影)"""
    from chemmcp.tools import ConformationalAnalyzer
    return _run_and_check(
        ConformationalAnalyzer,
        dict(molecule="butane", analysis_type="newman"),
        "200 ConformationalAnalyzer",
        min_length=300,
        must_contain=["butane", "newman", "conform"],
    )


if __name__ == "__main__":
    tests = [
        ("191 RingFormationStrategy", test_191_ring_formation_strategy),
        ("192 AsymmetricSynthesisGuide", test_192_asymmetric_synthesis_guide),
        ("193 TotalSynthesisPlanner", test_193_total_synthesis_planner),
        ("194 PkaPredictor", test_194_pka_predictor),
        ("195 HammettSigmaLookup", test_195_hammett_sigma_lookup),
        ("196 ResonanceStructureGenerator", test_196_resonance_structure_generator),
        ("197 InductiveEffectAnalyzer", test_197_inductive_effect_analyzer),
        ("198 HyperconjugationExplainer", test_198_hyperconjugation_explainer),
        ("199 StericEffectAnalyzer", test_199_steric_effect_analyzer),
        ("200 ConformationalAnalyzer", test_200_conformational_analyzer),
    ]

    passed = 0
    failed = 0
    results_log = []
    for name, fn in tests:
        try:
            fn()
            passed += 1
            results_log.append(f"✅ {name}")
        except Exception as e:
            failed += 1
            results_log.append(f"❌ {name}: {e}")
            print(f"  [FAIL] {name}: {e}")

    print(f"\n{'='*60}")
    print(f"Test Results: {passed}/{passed+failed} PASSED")
    for r in results_log:
        print(f"  {r}")
    if failed == 0:
        print("\n🎉 All 10 tests PASSED! Tools #191-200 are working correctly.")
    else:
        print(f"\n⚠️  {failed} test(s) need attention")
