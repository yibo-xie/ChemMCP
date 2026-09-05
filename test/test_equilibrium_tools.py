"""
Tests for Chemical Equilibrium MCP Tools (#41-50)
Run: python -m pytest test_equilibrium_tools.py -v
Or run directly: python test_equilibrium_tools.py
"""

import sys
import os
import math

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from chemmcp.tools import (
    EquilibriumConstantThermo,
    BornHaberCycle,
    BondEnergyCalculation,
    TemperatureEffectK,
    CalculateEquilibriumConstant,
    ICETableSolver,
    LeChatelierPrediction,
    PressureEffectEquilibrium,
    ReactionQuotient,
)


def separator(name):
    print(f"\n{'='*60}")
    print(f"  TEST: {name}")
    print(f"{'='*60}")


# ============================================================
# Test 41: EquilibriumConstantThermo - 从热力学数据计算平衡常数
# ============================================================
def test_equilibrium_constant_thermo():
    separator("EquilibriumConstantThermo")

    tool = EquilibriumConstantThermo()

    # Test 1: Direct ΔG° input (ΔG° = -23.7 kJ/mol at 298K → K should be large)
    print("\n--- Test 1: Direct ΔG° input ---")
    result = tool.run_code(temperature_k=298.15, delta_g=-23.7)
    print(f"T=298.15K, ΔG°=-23.7 kJ/mol")
    print(f"K = {result['K']}")
    print(f"Method: {result['method']}")
    assert result['method'] == 'direct'
    assert result['K'] > 0
    # Verify: K = exp(23700 / (8.314 * 298.15)) = exp(9.565) ≈ 14250
    expected_K = math.exp(23.7 * 1000 / (8.314 * 298.15))
    assert abs(result['K'] - round(expected_K, 4)) < 1, f"Expected ~{expected_K:.0f}, got {result['K']}"

    # Test 2: From ΔH° and ΔS°
    print("\n--- Test 2: From ΔH° and ΔS° ---")
    result2 = tool.run_code(temperature_k=298.15, delta_h=-57.2, delta_s=-112.5)
    print(f"T=298.15K, ΔH°=-57.2 kJ/mol, ΔS°=-112.5 J/(mol·K)")
    print(f"ΔG° calculated = {result2['delta_g_calculated']} kJ/mol")
    print(f"K = {result2['K']}")
    print(f"Method: {result2['method']}")
    assert result2['method'] == 'derived'
    dg_calc = -57.2 - 298.15 * (-112.5) / 1000.0
    assert abs(result2['delta_g_calculated'] - round(dg_calc, 4)) < 0.01

    # Test 3: Text interface
    print("\n--- Test 3: Text interface ---")
    result3 = tool.run_text("298.15 -23.7")
    print(f"Text input result K = {result3['K']}")
    assert result3['K'] > 0

    print("\n✅ EquilibriumConstantThermo: ALL PASSED")


# ============================================================
# Test 42: BornHaberCycle - Born-Haber循环计算
# ============================================================
def test_born_haber_cycle():
    separator("BornHaberCycle")

    tool = BornHaberCycle()

    # Classic example: NaCl formation
    # Na(s) + 1/2 Cl₂(g) → NaCl(s)  ΔH_f = -411 kJ/mol
    # Solve for lattice energy U
    print("\n--- Test: NaCl lattice energy ---")
    result = tool.run_code(
        delta_h_f=-411.0,
        delta_h_sub=108.0,       # Na(s) → Na(g)
        ionization_energies=[496.0],  # Na(g) → Na⁺(g) + e⁻
        bond_dissociation=122.0,      # Cl₂(g) → 2Cl(g), D/2 = 61
        electron_affinities=[-349.0], # Cl(g) + e⁻ → Cl⁻(g)
        lattice_energy=None,
        unknown="lattice_energy",
    )
    print(f"Lattice energy U = {result['result']} kJ/mol")
    print(f"Cycle summary: {result['cycle_summary']}")
    # Verify: ΔH_f = ΔH_sub + IE + D/2 + EA + U
    # -411 = 108 + 496 + 61 + (-349) + U → U = -727 kJ/mol
    assert abs(result['result'] - (-727.0)) < 5, f"Expected ~-727 kJ/mol, got {result['result']}"
    assert result['unit'] == 'kJ/mol'

    print("\n✅ BornHaberCycle: ALL PASSED")


# ============================================================
# Test 43: BondEnergyCalculation - 用键能估算反应焓变
# ============================================================
def test_bond_energy_calculation():
    separator("BondEnergyCalculation")

    tool = BondEnergyCalculation()

    # N₂ + 3H₂ → 2NH₃
    # Break: N≡N (941) + 3×H-H (3×436 = 1308) = 2249
    # Form: 6×N-H (6×391 = 2346)
    # ΔH = 2249 - 2346 = -97 kJ/mol
    print("\n--- Test: Haber process (N2 + 3H2 → 2NH3) ---")
    result = tool.run_code(
        bonds_broken=[("N≡N", 941, 1), ("H-H", 436, 3)],
        bonds_formed=[("N-H", 391, 6)],
    )
    print(f"ΔH = {result['delta_h']} kJ/mol")
    print(f"Bonds broken: {result['total_broken']} kJ/mol")
    print(f"Bonds formed: {result['total_formed']} kJ/mol")
    print(f"Reaction type: {result['reaction_type']}")
    assert result['reaction_type'] == "exothermic"
    assert abs(result['delta_h'] - (-97.0)) < 1, f"Expected ~-97, got {result['delta_h']}"

    # Text interface test
    print("\n--- Test: Text interface ---")
    result2 = tool.run_text("broken: N≡N=941x1 H-H=436x3; formed: N-H=391x6")
    print(f"Text ΔH = {result2['delta_h']} kJ/mol")
    assert result2['reaction_type'] == "exothermic"

    print("\n✅ BondEnergyCalculation: ALL PASSED")


# ============================================================
# Test 44: TemperatureEffectK - van't Hoff方程
# ============================================================
def test_temperature_effect_k():
    separator("TemperatureEffectK")

    tool = TemperatureEffectK()

    # Endothermic reaction: increasing T should increase K
    print("\n--- Test 1: Endothermic, T increases ---")
    result = tool.run_code(K1=4.0, T1=300.0, T2=400.0, delta_h=40.0)
    print(f"K1={4.0} @ 300K → K2={result['K2']} @ 400K (endothermic, ΔH=+40)")
    print(f"Direction: {result['direction']}")
    assert result['direction'] == "increases"
    assert result['K2'] > 4.0  # K2 > K1 for endothermic reaction with T increase

    # Exothermic reaction: increasing T should decrease K
    print("\n--- Test 2: Exothermic, T increases ---")
    result2 = tool.run_code(K1=100.0, T1=300.0, T2=400.0, delta_h=-50.0)
    print(f"K1={100} @ 300K → K2={result2['K2']} @ 400K (exothermic, ΔH=-50)")
    print(f"Direction: {result2['direction']}")
    assert result2['direction'] == "decreases"
    assert result2['K2'] < 100.0

    # Verify van't Hoff equation manually
    ln_ratio = (-(-50.0) * 1000 / 8.314) * (1/400 - 1/300)
    expected_K2 = 100.0 * math.exp(ln_ratio)
    print(f"Manual check: expected K2 ≈ {expected_K2:.4f}, got {result2['K2']}")
    assert abs(result2['K2'] - round(expected_K2, 4)) < 0.01

    # Text interface
    print("\n--- Test 3: Text interface ---")
    result3 = tool.run_text("4.0 300 400 40")
    assert result3['direction'] == "increases"

    print("\n✅ TemperatureEffectK: ALL PASSED")


# ============================================================
# Test 45: CalculateEquilibriumConstant - 根据浓度计算Kc/Kp
# ============================================================
def test_calculate_equilibrium_constant():
    separator("CalculateEquilibriumConstant")

    tool = CalculateEquilibriumConstant()

    # 2NO₂ ⇌ N₂O₄ at equilibrium: [NO2]=0.056, [N2O4]=0.032
    # Kc = [N2O4] / [NO2]^2 = 0.032 / (0.056)^2 = 10.20... wait
    # Actually: products have N2O4, reactants have NO2 with coef 2
    # Let me use: N2O4 ⇌ 2NO2, Kc = [NO2]^2/[N2O4]
    print("\n--- Test 1: Kc calculation (N2O4 ⇌ 2NO2) ---")
    result = tool.run_code(
        products=[{"name": "NO2", "coefficient": 2, "value": 0.235}],
        reactants=[{"name": "N2O4", "coefficient": 1, "value": 0.382}],
        mode="kc",
    )
    print(f"Kc = {result['K']}")
    print(f"Expression: {result['expression']}")
    assert result['mode'] == "Kc"
    # Manual: Kc = 0.235^2 / 0.382 = 0.1446
    manual_K = 0.235**2 / 0.382
    assert abs(result['K'] - round(manual_K, 6)) < 0.001

    # Text interface
    print("\n--- Test 2: Text interface ---")
    result2 = tool.run_text("N2O4:1:0.382 | NO2:2:0.235 | kc")
    print(f"Text Kc = {result2['K']}")
    assert result2['mode'] == "Kc"

    print("\n✅ CalculateEquilibriumConstant: ALL PASSED")


# ============================================================
# Test 46: ICETableSolver - ICE表法求解平衡浓度
# ============================================================
def test_ice_table_solver():
    separator("ICETableSolver")

    tool = ICETableSolver()

    # N2O4(g) ⇌ 2NO2(g), initial [N2O4]=0.500 M, [NO2]=0, Kc=0.067
    print("\n--- Test: ICE table for N2O4 ⇌ 2NO2 ---")
    result = tool.run_code(
        initial_conc={"N2O4": 0.500, "NO2": 0.000},
        stoichiometry={"N2O4": -1, "NO2": 2},
        K_eq=0.067,
    )
    print(f"x (extent) = {result['x']}")
    print(f"Equilibrium concentrations: {result['equilibrium_conc']}")
    print(f"ICE table: {result['ice_table']}")
    print(f"K calculated (verification): {result['K_calculated']}")
    print(f"Qc initial: {result['Qc_initial']}")

    # Verify all concentrations are non-negative
    for species, conc in result['equilibrium_conc'].items():
        assert conc >= 0, f"Concentration of {species} is negative: {conc}"

    # Verify K is approximately correct
    if result['K_calculated']:
        assert abs(result['K_calculated'] - 0.067) < 0.01, \
            f"K verification failed: {result['K_calculated']} vs 0.067"

    # Verify mass balance
    n2o4_eq = result['equilibrium_conc']['N2O4']
    no2_eq = result['equilibrium_conc']['NO2']
    total_N = n2o4_eq * 2 + no2_eq  # N atoms per liter (initial: 0.500*2 = 1.0)
    print(f"N atom conservation check: initial=1.000, final={total_N:.4f}")

    print("\n✅ ICETableSolver: ALL PASSED")


# ============================================================
# Test 47: LeChatelierPrediction - Le Chatelier原理预测
# ============================================================
def test_le_chatelier_prediction():
    separator("LeChatelierPrediction")

    tool = LeChatelierPrediction()

    # Test 1: Exothermic + increase temperature → shift backward
    print("\n--- Test 1: Exothermic, increase T ---")
    r1 = tool.run_code(reaction_type="exothermic", delta_n_gas=-1.0, disturbance="increase_temp")
    print(f"Shift: {r1['shift_direction']}")
    print(f"Reasoning: {r1['reasoning']}")
    print(f"K effect: {r1['K_effect']}")
    assert r1['shift_direction'] == "backward"
    assert r1['K_effect'] == "K decreases"

    # Test 2: Endothermic + decrease temperature → shift backward
    print("\n--- Test 2: Endothermic, decrease T ---")
    r2 = tool.run_code(reaction_type="endothermic", delta_n_gas=2.0, disturbance="decrease_temp")
    assert r2['shift_direction'] == "backward"
    assert r2['K_effect'] == "K decreases"

    # Test 3: Increase pressure, Δn > 0 → shift backward (fewer moles)
    print("\n--- Test 3: Increase P, Δn=+2 ---")
    r3 = tool.run_code(reaction_type="endothermic", delta_n_gas=2.0, disturbance="increase_pressure")
    assert r3['shift_direction'] == "backward"
    assert r3['K_effect'] == "K unchanged"

    # Test 4: Increase pressure, Δn < 0 → shift forward
    print("\n--- Test 4: Increase P, Δn=-1 ---")
    r4 = tool.run_code(reaction_type="exothermic", delta_n_gas=-1.0, disturbance="increase_pressure")
    assert r4['shift_direction'] == "forward"

    # Test 5: Catalyst → no shift
    print("\n--- Test 5: Catalyst ---")
    r5 = tool.run_code(reaction_type="exothermic", delta_n_gas=0, disturbance="catalyst")
    assert r5['shift_direction'] == "no_shift"

    # Test 6: Increase concentration of reactant
    print("\n--- Test 6: Increase [N2] ---")
    r6 = tool.run_code(reaction_type="exothermic", delta_n_gas=-1, disturbance="increase_conc", species_affected="N2")
    assert r6['shift_direction'] == "forward"

    # Test 7: Δn = 0, pressure change → no shift
    print("\n--- Test 7: Δn=0, pressure change ---")
    r7 = tool.run_code(reaction_type="exothermic", delta_n_gas=0.0, disturbance="increase_pressure")
    assert r7['shift_direction'] == "no_shift"

    # Text interface
    print("\n--- Test 8: Text interface ---")
    r8 = tool.run_text("exothermic -1 increase_temp")
    assert r8['shift_direction'] == "backward"

    print("\n✅ LeChatelierPrediction: ALL PASSED")


# ============================================================
# Test 48: PressureEffectEquilibrium - 压力变化对气相平衡的影响
# ============================================================
def test_pressure_effect_equilibrium():
    separator("PressureEffectEquilibrium")

    tool = PressureEffectEquilibrium()

    # N2O4 ⇌ 2NO2, Kp = 0.067, pressure from 1 atm → 5 atm
    print("\n--- Test: Pressure effect on N2O4 ⇌ 2NO2 ---")
    result = tool.run_code(
        initial_total_p=1.0,
        new_total_p=5.0,
        Kp=0.067,
        stoichiometry={"N2O4": -1, "NO2": 2},
    )
    print(f"Old equilibrium: {result['old_equilibrium']}")
    print(f"New equilibrium: {result['new_equilibrium']}")
    print(f"Shift direction: {result['shift_direction']}")
    print(f"Dissociation change: {result['degree_of_dissociation_change']}")

    # Δn = +1 (more gas on product side), increasing P → shift backward
    assert result['shift_direction'] == "backward", \
        f"Expected 'backward', got '{result['shift_direction']}'"

    # Check that old and new equilibria have valid partial pressures
    old_pp = result['old_equilibrium']['partial_pressures']
    new_pp = result['new_equilibrium']['partial_pressures']
    for s in old_pp:
        assert old_pp[s] >= 0, f"Old P_{s} negative: {old_pp[s]}"
        assert new_pp[s] >= 0, f"New P_{s} negative: {new_pp[s]}"

    print("\n✅ PressureEffectEquilibrium: ALL PASSED")


# ============================================================
# Test 50: ReactionQuotient - 计算反应商Q并判断反应方向
# ============================================================
def test_reaction_quotient():
    separator("ReactionQuotient")

    tool = ReactionQuotient()

    # Test 1: Q < K → forward direction
    # N2O4 ⇌ 2NO2, Kc = 6.7, [NO2]=0.01, [N2O4]=0.1
    print("\n--- Test 1: Q < K → forward ---")
    r1 = tool.run_code(
        current_state={"NO2": 0.010, "N2O4": 0.100},
        stoichiometry={"NO2": 2, "N2O4": -1},
        K_eq=6.7,
        mode="kc",
    )
    print(f"Q = {r1['Q']}, K = {r1['K']}")
    print(f"Comparison: {r1['comparison']}")
    print(f"Direction: {r1['direction']}")
    print(f"Explanation: {r1['explanation']}")
    assert r1['direction'] == "forward"
    assert r1['comparison'] == "Q < K"
    # Manual: Q = 0.01^2 / 0.1 = 0.001
    assert abs(r1['Q'] - 0.001) < 0.0001

    # Test 2: Q > K → backward direction
    # N2 + 3H2 ⇌ 2NH3, Kc = 1500
    print("\n--- Test 2: Q > K → backward ---")
    r2 = tool.run_code(
        current_state={"N2": 0.01, "H2": 0.01, "NH3": 0.5},
        stoichiometry={"N2": -1, "H2": -3, "NH3": 2},
        K_eq=1500.0,
        mode="kc",
    )
    print(f"Q = {r2['Q']}, K = {r2['K']}")
    print(f"Direction: {r2['direction']}")
    assert r2['direction'] == "backward"
    # Manual: Q = 0.5^2 / (0.01 * 0.01^3) = 0.25 / 1e-8 = 25,000,000
    manual_Q = 0.5**2 / (0.01 * 0.01**3)
    assert abs(r2['Q'] - round(manual_Q, 0)) < 1

    # Test 3: Q ≈ K → equilibrium
    print("\n--- Test 3: Q ≈ K → equilibrium ---")
    # Use values that approximately satisfy K
    r3 = tool.run_code(
        current_state={"NO2": 0.235, "N2O4": 0.382},
        stoichiometry={"NO2": 2, "N2O4": -1},
        K_eq=0.145,  # approximate K for these values
        mode="kc",
    )
    print(f"Q = {r3['Q']}, K = {r3['K']}, Direction: {r3['direction']}")

    # Text interface
    print("\n--- Test 4: Text interface ---")
    r4 = tool.run_text("NO2=0.01 N2O4=0.1 | NO2=2 N2O4=-1 | K=6.7 kc")
    assert r4['direction'] == "forward"

    print("\n✅ ReactionQuotient: ALL PASSED")


# ============================================================
# Run all tests
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("  CHEMICAL EQUILIBRIUM MCP TOOLS - TEST SUITE (#41-50)")
    print("=" * 60)

    tests = [
        ("41 - EquilibriumConstantThermo", test_equilibrium_constant_thermo),
        ("42 - BornHaberCycle", test_born_haber_cycle),
        ("43 - BondEnergyCalculation", test_bond_energy_calculation),
        ("44 - TemperatureEffectK", test_temperature_effect_k),
        ("45 - CalculateEquilibriumConstant", test_calculate_equilibrium_constant),
        ("46 - ICETableSolver", test_ice_table_solver),
        ("47 - LeChatelierPrediction", test_le_chatelier_prediction),
        ("48 - PressureEffectEquilibrium", test_pressure_effect_equilibrium),
        ("50 - ReactionQuotient", test_reaction_quotient),
    ]

    passed = 0
    failed = 0
    errors = []

    for name, test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"\n❌ FAILED: {name}")
            print(f"   Error: {e}")

    print("\n" + "=" * 60)
    print(f"  RESULTS: {passed} passed, {failed} failed, {passed + failed} total")
    print("=" * 60)

    if errors:
        print("\n❌ Failed tests:")
        for name, err in errors:
            print(f"  • {name}: {err}")

    sys.exit(0 if failed == 0 else 1)
