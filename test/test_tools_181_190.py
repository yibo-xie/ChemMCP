"""
Test suite for ChemMCP tools #181-190.
Run: python -m pytest test/test_tools_181_190.py -v
Or: python test/test_tools_181_190.py
"""

import sys
sys.path.insert(0, "src")

def test_isotope_pattern_generator():
    from chemmcp.tools import IsotopePatternGenerator
    tool = IsotopePatternGenerator()
    result = tool.run_code(molecular_formula="CHCl3", charge=0, min_abundance=0.01)
    print(f"[181] IsotopePatternGenerator CHCl3: peaks={result['peak_count']}, nominal_mass={result['nominal_mass']}")
    assert result["peak_count"] >= 2, f"Expected at least 2 peaks for CHCl3, got {result['peak_count']}"
    assert result["nominal_mass"] == 118, f"Expected nominal mass 118, got {result['nominal_mass']}"
    assert len(result["peaks"]) > 0
    print(f"  Peaks: {[(p['mz'], p['abundance_pct']) for p in result['peaks']]}")
    print("  [PASS] IsotopePatternGenerator")


def test_uv_vis_predictor():
    from chemmcp.tools import UvVisPredictor
    tool = UvVisPredictor()
    result = tool.run_code(smiles="C=CC=CC(C)=C", solvent="methanol", analysis_mode="auto")
    print(f"[182] UvVisPredictor diene: chromophores={result['detected_chromophores']}, bands={len(result['predicted_absorptions'])}")
    assert "predicted_absorptions" in result
    assert len(result["predicted_absorptions"]) > 0
    for band in result["predicted_absorptions"]:
        assert "lambda_max_nm" in band
        assert 170 <= band["lambda_max_nm"] <= 600, f"λmax out of range: {band['lambda_max_nm']}"
    print(f"  Absorptions: {[(a['chromophore'], a['lambda_max_nm']) for a in result['predicted_absorptions']]}")
    print("  [PASS] UvVisPredictor")


def test_spectrum_to_structure():
    from chemmcp.tools import SpectrumToStructure
    tool = SpectrumToStructure()
    result = tool.run_code(
        molecular_formula="C6H12O",
        ir_peaks=[[1715, "strong"], [2950, "strong"]],
        nmr_peaks=[[2.1, "s", 3], [2.4, "t", 2], [1.6, "sextet", 2], [0.95, "t", 3]],
        ms_mz=100.0,
    )
    print(f"[183] SpectrumToStructure C6H12O: DOU={result['degree_of_unsaturation']}, FGs={result['suggested_functional_groups']}")
    assert result["degree_of_unsaturation"] == 1
    assert len(result["suggested_functional_groups"]) > 0
    assert "structural_hints" in result
    print("  [PASS] SpectrumToStructure")


def test_dept_interpreter():
    from chemmcp.tools import DeptInterpreter
    tool = DeptInterpreter()
    result = tool.run_code(
        dept_90_peaks=[28.0, 68.0],
        dept_135_peaks=[14.0, 22.5, 28.0, 31.5, 42.0, 68.0],
        regular_13c_peaks=[14.0, 22.5, 28.0, 31.5, 42.0, 60.0, 68.0, 210.0],
    )
    counts = result["summary_counts"]
    print(f"[184] DeptInterpreter: CH3={counts['CH3']}, CH2={counts['CH2']}, CH={counts['CH']}, Cq={counts['Cq']}")
    assert counts["CH3"] >= 1
    assert counts["Cq"] >= 1  # Should detect quaternary carbons (60 and 210 ppm)
    assert len(result["carbon_assignments"]) == 8
    print("  [PASS] DeptInterpreter")


def test_cosy_noesy_guide():
    from chemmcp.tools import CosyNoesyGuide
    tool = CosyNoesyGuide()
    result = tool.run_code(
        peaks_1d=[[7.25, "m", 5], [4.15, "q", 2], [3.65, "s", 3], [2.85, "t", 2], [1.35, "t", 3]],
        cross_peaks=[[4.15, 2.85], [2.85, 1.35]],
        experiment_type="cosy",
    )
    print(f"[185] CosyNoesyGuide: spin_systems={len(result['spin_systems'])}, inferences={len(result['structural_inferences'])}")
    assert len(result["spin_systems"]) > 0
    assert "connectivity_map" in result
    assert "recommended_next_experiments" in result
    print("  [PASS] CosyNoesyGuide")


def test_retrosynthesis_analyzer():
    from chemmcp.tools import RetrosynthesisAnalyzer
    tool = RetrosynthesisAnalyzer()
    result = tool.run_code(target_smiles="CC(=O)c1ccccc1", max_depth=2, focus_strategy="auto")
    print(f"[186] RetrosynthesisAnalyzer: steps={result['total_steps_estimate']}, strategy={result['strategy_used']}")
    assert result["total_steps_estimate"] >= 1
    assert len(result["suggested_disconnections"]) > 0
    for disc in result["suggested_disconnections"]:
        assert "bond" in disc or "disconnection" in disc
        assert "rationale" in disc
    print("  [PASS] RetrosynthesisAnalyzer")


def test_synthon_identifier():
    from chemmcp.tools import SynthonIdentifier
    tool = SynthonIdentifier()
    result = tool.run_code(target_smiles="CC(=O)c1ccccc1", disconnection_bond=-1, include_reagents=True)
    print(f"[187] SynthonIdentifier: synthons_found={result['total_synthon_pairs_found']}")
    assert result["total_synthon_pairs_found"] > 0
    assert "recommended_disconnection" in result
    for syn in result["identified_synthons"]:
        assert "synthon_plus" in syn or "synthon_minus" in syn
    print("  [PASS] SynthonIdentifier")


def test_disconnection_suggester():
    from chemmcp.tools import DisconnectionSuggester
    tool = DisconnectionSuggester()
    result = tool.run_code(target_smiles="CC(=O)Oc1ccccc1", strategy="auto", max_suggestions=3)
    print(f"[188] DisconnectionSuggester: sites={result['total_sites_found']}, top_score={result['ranked_disconnections'][0]['score'] if result['ranked_disconnections'] else 'N/A'}")
    assert result["total_sites_found"] > 0
    assert len(result["ranked_disconnections"]) > 0
    top = result["ranked_disconnections"][0]
    assert "position" in top
    assert "difficulty" in top
    print(f"  Top suggestion: {top['position']} (score={top['score']}, {top['difficulty']})")
    print("  [PASS] DisconnectionSuggester")


def test_functional_group_interconversion():
    from chemmcp.tools import FunctionalGroupInterconversion
    tool = FunctionalGroupInterconversion()
    result = tool.run_code(source_fg="alcohol", target_fg="aldehyde", context_smiles="")
    print(f"[189] FG Interconversion alcohol→aldehyde: possible={result['conversion_possible']}, pathways={result.get('pathway_count', 0)}")
    assert result["conversion_possible"] is True
    assert result["pathway_count"] >= 2  # Should have PCC, Swern, DMP at minimum
    assert "best_method" in result
    for pathway in result["pathways"]:
        assert "method" in pathway
        assert "reagents" in pathway
        assert "yield" in pathway
    print(f"  Best method: {result['best_method']} ({result['best_yield']})")
    print("  [PASS] FunctionalGroupInterconversion")


def test_carbon_chain_builder():
    from chemmcp.tools import CarbonChainBuilder
    tool = CarbonChainBuilder()
    result = tool.run_code(current_smiles="CCO", operation="elongate", target_carbon_count=0, method_preference="auto")
    print(f"[190] CarbonChainBuilder CCO: c_count={result['current_carbon_count']}, methods={len(result['suggested_methods'])}")
    assert result["current_carbon_count"] == 2
    assert result["operation"] == "elongate"
    assert len(result["suggested_methods"]) > 0
    for m in result["suggested_methods"]:
        assert "name" in m
        assert "reagents" in m
        assert "expected_yield" in m
    print(f"  Summary: {result['summary']}")
    print("  [PASS] CarbonChainBuilder")


# === Text interface tests ===
def test_text_interface_isotope():
    from chemmcp.tools import IsotopePatternGenerator
    tool = IsotopePatternGenerator()
    result = tool.run_text("CHCl3 0 0.01")
    assert result["formula"] == "CHCl3"
    print("[TEXT PASS] IsotopePatternGenerator text input")


def test_text_interface_uvvis():
    from chemmcp.tools import UvVisPredictor
    tool = UvVisPredictor()
    result = tool.run_text("C=CC=CC(C)=C methanol auto")
    assert len(result["predicted_absorptions"]) > 0
    print("[TEXT PASS] UvVisPredictor text input")


def test_text_interface_fg_conversion():
    from chemmcp.tools import FunctionalGroupInterconversion
    tool = FunctionalGroupInterconversion()
    result = tool.run_text("alcohol aldehyde")
    assert result["conversion_possible"] is True
    print("[TEXT PASS] FunctionalGroupInterconversion text input")


if __name__ == "__main__":
    tests = [
        ("IsotopePatternGenerator (#181)", test_isotope_pattern_generator),
        ("UvVisPredictor (#182)", test_uv_vis_predictor),
        ("SpectrumToStructure (#183)", test_spectrum_to_structure),
        ("DeptInterpreter (#184)", test_dept_interpreter),
        ("CosyNoesyGuide (#185)", test_cosy_noesy_guide),
        ("RetrosynthesisAnalyzer (#186)", test_retrosynthesis_analyzer),
        ("SynthonIdentifier (#187)", test_synthon_identifier),
        ("DisconnectionSuggester (#188)", test_disconnection_suggester),
        ("FunctionalGroupInterconversion (#189)", test_functional_group_interconversion),
        ("CarbonChainBuilder (#190)", test_carbon_chain_builder),
        ("Text Interface - Isotope", test_text_interface_isotope),
        ("Text Interface - UV-Vis", test_text_interface_uvvis),
        ("Text Interface - FG Conv", test_text_interface_fg_conversion),
    ]

    passed = 0
    failed = 0
    errors = []

    for name, test_fn in tests:
        try:
            test_fn()
            passed += 1
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"  [FAIL] {name}: {e}")

    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed out of {len(tests)} tests")
    if errors:
        print("\nFailed tests:")
        for name, err in errors:
            print(f"  ❌ {name}: {err[:200]}")
    print(f"{'='*60}")

    sys.exit(0 if failed == 0 else 1)
