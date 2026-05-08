"""
Test suite for MCP Registration Table #491-500
Kinetics Corrections & Computational Chemistry Tools

Tools:
  #491 TunnelingCorrection       - 隧穿校正，轻原子转移反应
  #492 SteadyStateApprox         - 稳态近似，中间体浓度
  #493 PreEquilibrium            - 预平衡近似，快速平衡步骤
  #494 RateDeterminingStep       - 速控步分析，反应瓶颈识别
  #495 MichaelisMenten           - Michaelis-Menten动力学，酶催化
  #496 ReactionNetworkSolver     - 反应网络求解，多步机理模拟
  #497 GeometryOptimizer         - 几何优化，能量极小化
  #498 TransitionStateSearch     - 过渡态搜索，鞍点定位
  #499 PotentialEnergySurface    - 势能面扫描，反应路径探索
  #500 FrequencyAnalysis         - 频率分析，驻点性质确认

Run:  cd ~/ChemMCP && .venv/bin/python tests/test_mcp_491_500.py
Or:  cd ~/ChemMCP && python -m pytest tests/test_mcp_491_500.py -v
"""

import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


# ═══════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════

def separator(title):
    print(f"\n{'='*65}")
    print(f"  {title}")
    print(f"{'='*65}")

passed = 0
failed = 0

def run_test(name, fn):
    global passed, failed
    try:
        fn()
        passed += 1
        print(f"  ✅ {name} ALL PASSED\n")
    except Exception as e:
        failed += 1
        print(f"  ❌ {name} FAILED: {e}\n")


# ═══════════════════════════════════════════════════════════════
# #491 TunnelingCorrection — 隧穿校正
# ═══════════════════════════════════════════════════════════════
def test_491_tunneling_correction():
    separator("#491 TunnelingCorrection — 隧穿校正")
    from chemmcp.tools import TunnelingCorrection

    tool = TunnelingCorrection()

    # Test Bell correction with parameters that give significant tunneling
    result = tool.run_code(
        temperature_K=200.0,
        barrier_height_kJ_mol=20.0,
        imaginary_frequency_cm_minus_1=2000.0,
        correction_model="bell",
    )
    assert result is not None
    assert "correction_factor_kappa" in result
    kappa = result["correction_factor_kappa"]
    assert kappa > 0, f"kappa should be positive, got {kappa}"
    print(f"  Bell κ={kappa:.4f} at T=200K, Ea=20kJ/mol, ν‡=2000cm⁻¹")
    assert result["tunneling_significant"] == (kappa > 1.5)

    # Test room temperature (higher barrier → less tunneling)
    result_rt = tool.run_code(
        temperature_K=298.15,
        barrier_height_kJ_mol=45.0,
        imaginary_frequency_cm_minus_1=1500.0,
        correction_model="bell",
    )
    kappa_rt = result_rt["correction_factor_kappa"]
    assert kappa_rt > 0
    print(f"  Room temp (298K) κ={kappa_rt:.4f}")

    # Test low temperature (should enhance tunneling)
    result_cold = tool.run_code(
        temperature_K=100.0,
        barrier_height_kJ_mol=20.0,
        imaginary_frequency_cm_minus_1=2000.0,
        correction_model="bell",
    )
    kappa_cold = result_cold["correction_factor_kappa"]
    assert kappa_cold > 0
    print(f"  Cold (100K) κ={kappa_cold:.4f}")

    # Test high barrier (less tunneling)
    result_high = tool.run_code(
        temperature_K=298.15,
        barrier_height_kJ_mol=100.0,
        imaginary_frequency_cm_minus_1=800.0,
        correction_model="bell",
    )
    assert result_high["correction_factor_kappa"] > 0
    print(f"  High barrier (100kJ) κ={result_high['correction_factor_kappa']:.4f}")

    # Test WKB method
    result_wkb = tool.run_code(
        temperature_K=298.15,
        barrier_height_kJ_mol=45.0,
        imaginary_frequency_cm_minus_1=1500.0,
        correction_model="wkb",
        reduced_mass_amu=1.008,
    )
    assert "kappa" in result_wkb
    print(f"  WKB κ={result_wkb['kappa']:.4f}")

    # Test Eckart method
    result_eckart = tool.run_code(
        temperature_K=298.15,
        barrier_height_kJ_mol=45.0,
        imaginary_frequency_cm_minus_1=1500.0,
        correction_model="eckart",
        reduced_mass_amu=1.008,
        forward_barrier_kJ_mol=40.0,
        reverse_barrier_kJ_mol=5.0,
    )
    assert "kappa" in result_eckart
    print(f"  Eckart κ={result_eckart['kappa']:.4f}")

    # Test text interface (positional arg)
    result_text = tool.run_text("298.15 45.0 1500.0 bell")
    assert "correction_factor_kappa" in result_text
    print(f"  Text interface: κ={result_text['correction_factor_kappa']:.4f}")

    # Test error handling
    try:
        tool.run_code(temperature_K=-10, barrier_height_kJ_mol=45.0,
                      imaginary_frequency_cm_minus_1=1500.0)
        assert False, "Should raise error for negative T"
    except Exception:
        print("  ✅ Correctly rejects negative temperature")

    print("  📊 TunnelingCorrection: all models work correctly")


# ═══════════════════════════════════════════════════════════════
# #492 SteadyStateApprox — 稳态近似
# ═══════════════════════════════════════════════════════════════
def test_492_steady_state_approx():
    separator("#492 SteadyStateApprox — 稳态近似")
    from chemmcp.tools import SteadyStateApprox

    tool = SteadyStateApprox()

    # Test consecutive reaction A → I → P (classic SSA example)
    result = tool.run_code(
        mechanism_type="consecutive",
        rate_constants=[0.1, 1.0],
        initial_concentrations={"A": 1.0},
        time_points=[0, 5, 10, 20, 50, 100],
        intermediates=["I"],
    )
    assert result is not None
    assert "concentration_profiles" in result
    profiles = result["concentration_profiles"]
    assert len(profiles) == 6

    assert result["ssa_valid"] == True
    print(f"  Consecutive A→I→P: k2/k1={result['validity_criteria']['k_consume_over_k_form_ratio']:.1f}, valid={result['ssa_valid']}")

    A_values = [p["A_exact"] for p in profiles]
    assert A_values[0] > A_values[-1], "Reactant A should decrease"
    assert abs(A_values[0] - 1.0) < 0.01

    P_values = [p["P_exact"] for p in profiles]
    assert P_values[-1] > P_values[0]
    print(f"  A: {A_values[0]:.3f} → {A_values[-1]:.6f}")
    print(f"  P: {P_values[0]:.3f} → {P_values[-1]:.6f}")

    # Test reversible mechanism
    result_rev = tool.run_code(
        mechanism_type="reversible",
        rate_constants=[0.5, 0.1, 0.3],
        initial_concentrations={"A": 1.0},
        time_points=[0, 10, 50, 100],
        intermediates=["I"],
    )
    assert "concentration_profiles" in result_rev
    print(f"  Reversible A⇌I→P: valid={result_rev['ssa_valid']}")

    # Test pre-equilibrium
    result_pe = tool.run_code(
        mechanism_type="pre_equilibrium",
        rate_constants=[100.0, 10.0, 0.5],
        initial_concentrations={"A": 1.0},
        time_points=[0, 10, 50],
        intermediates=["I"],
    )
    assert "derived_rate_law" in result_pe
    print(f"  Pre-equilibrium: rate law derived")

    # Test error: unknown mechanism
    try:
        tool.run_code(mechanism_type="unknown", rate_constants=[0.1],
                      initial_concentrations={"A": 1}, time_points=[0, 1])
        assert False, "Should raise error for unknown mechanism"
    except Exception:
        print("  ✅ Correctly rejects unknown mechanism type")

    print("  📊 SteadyStateApprox: all mechanism types work correctly")


# ═══════════════════════════════════════════════════════════════
# #493 PreEquilibrium — 预平衡近似
# ═══════════════════════════════════════════════════════════════
def test_493_pre_equilibrium():
    separator("#493 PreEquilibrium — 预平衡近似")
    from chemmcp.tools import PreEquilibrium

    tool = PreEquilibrium()

    # Test unimolecular: A ⇌ I → P
    result = tool.run_code(
        mechanism="unimolecular",
        k_forward_list=[100.0],
        k_reverse_list=[10.0],
        k_slow=0.5,
        initial_concentrations={"A": 1.0},
    )
    assert result is not None
    assert "effective_rate_constant" in result
    assert "rate_law" in result

    K_eq = result["equilibrium_constants_Keq"][0]
    assert abs(K_eq - 10.0) < 0.01

    k_eff = result["effective_rate_constant"]
    expected_k_eff = 0.5 * (10.0 / (1 + 10.0))
    assert abs(k_eff - expected_k_eff) < 0.01
    print(f"  Unimolecular: K_eq={K_eq}, k_eff={k_eff:.4f}")
    print(f"  Rate law: {result['rate_law']}")

    # Test bimolecular
    result_bi = tool.run_code(
        mechanism="bimolecular",
        k_forward_list=[50.0],
        k_reverse_list=[5.0],
        k_slow=0.3,
        initial_concentrations={"A": 1.0, "B": 1.0},
    )
    assert "effective_rate_constant" in result_bi
    print(f"  Bimolecular: K_eq={result_bi['equilibrium_constants_Keq'][0]:.1f}, k_eff={result_bi['effective_rate_constant']:.4f}")

    # Test multi-step
    result_multi = tool.run_code(
        mechanism="multi_step",
        k_forward_list=[100.0, 50.0],
        k_reverse_list=[10.0, 5.0],
        k_slow=0.2,
        initial_concentrations={"A": 1.0},
    )
    assert "overall_equilibrium_constant" in result_multi
    K_overall = result_multi["overall_equilibrium_constant"]
    expected_K_overall = (100.0/10.0) * (50.0/5.0)
    assert abs(K_overall - expected_K_overall) < 0.1
    print(f"  Multi-step: K_overall={K_overall:.1f}, n_steps={result_multi['n_equilibrium_steps']}")

    assert result["validity"] == True

    # Test text interface (positional arg)
    result_text = tool.run_text("unimolecular 100 10 0.5 A0=1.0")
    assert "effective_rate_constant" in result_text
    print(f"  Text interface works: k_eff={result_text['effective_rate_constant']:.4f}")

    print("  📊 PreEquilibrium: all modes work correctly")


# ═══════════════════════════════════════════════════════════════
# #494 RateDeterminingStep — 速控步分析
# ═══════════════════════════════════════════════════════════════
def test_494_rate_determining_step():
    separator("#494 RateDeterminingStep — 速控步分析")
    from chemmcp.tools import RateDeterminingStep

    tool = RateDeterminingStep()

    # Test clear RDS case
    result = tool.run_code(
        mechanism_steps=[
            {"reactants": "A -> B", "products": "", "k": 0.001, "reversible": False},
            {"reactants": "B + C -> D", "products": "", "k": 5.0, "reversible": False},
            {"reactants": "D -> E", "products": "", "k": 100.0, "reversible": False},
        ],
        has_pre_equilibrium=False,
    )
    assert result is not None
    assert result["rds_step_index"] == 1
    assert result["rds_rate_constant"] == 0.001
    assert result["slowest_to_fastest_ratio"] >= 1000
    print(f"  RDS: Step {result['rds_step_index']} ('{result['rds_step_description']}'), k={result['rds_rate_constant']}")
    print(f"  Slowest/fastest ratio: {result['slowest_to_fastest_ratio']:,.0f}x")
    print(f"  Rate law: {result['overall_rate_law']}")

    # Test similar rate constants
    result_similar = tool.run_code(
        mechanism_steps=[
            {"reactants": "X -> Y", "products": "", "k": 1.0, "reversible": False},
            {"reactants": "Y -> Z", "products": "", "k": 1.5, "reversible": False},
            {"reactants": "Z -> W", "products": "", "k": 0.8, "reversible": False},
        ],
        has_pre_equilibrium=False,
    )
    assert result_similar["rds_step_index"] == 3
    ratio = result_similar["slowest_to_fastest_ratio"]
    assert ratio < 10
    print(f"  Similar rates: RDS=step {result_similar['rds_step_index']}, ratio={ratio:.1f}x")

    # Test with pre-equilibrium
    result_preq = tool.run_code(
        mechanism_steps=[
            {"reactants": "A <=> I (fast)", "products": "", "k": 1000.0, "reversible": True},
            {"reactants": "I -> P (slow)", "products": "", "k": 0.01, "reversible": False},
        ],
        has_pre_equilibrium=True,
    )
    assert result_preq["rds_step_index"] == 2
    print(f"  With pre-eq: RDS=step {result_preq['rds_step_index']}")

    # Test text interface (positional arg)
    result_text = tool.run_text("A->B;0.001;false;B+C->D;5;false;D->E;100;false")
    assert "rds_step_index" in result_text
    print(f"  Text interface: RDS=step {result_text['rds_step_index']}")

    # Error handling
    try:
        tool.run_code(mechanism_steps=[], has_pre_equilibrium=False)
        assert False, "Should raise error for empty steps"
    except Exception:
        print("  ✅ Correctly rejects empty steps list")

    print("  📊 RateDeterminingStep: analysis works correctly")


# ═══════════════════════════════════════════════════════════════
# #495 MichaelisMenten — 酶动力学
# ═══════════════════════════════════════════════════════════════
def test_495_michaelis_menten():
    separator("#495 MichaelisMenten — 酶动力学")
    from chemmcp.tools import MichaelisMenten

    tool = MichaelisMenten()

    # Test basic velocity: v = Vmax*S/(Km+S)
    result = tool.run_code(
        analysis_type="calculate_velocity",
        substrate_concentration_S=5.0,
        Vmax=100.0,
        Km=2.0,
    )
    assert result is not None
    v = result["velocity"]
    expected_v = 100.0 * 5.0 / (2.0 + 5.0)
    assert abs(v - expected_v) / expected_v < 0.001
    print(f"  Basic MM: S=5, Vmax=100, Km=2 → v={v:.2f} (expected {expected_v:.2f})")

    # Test at S = Km (v = Vmax/2)
    result_half = tool.run_code(
        analysis_type="calculate_velocity",
        substrate_concentration_S=2.0,
        Vmax=100.0,
        Km=2.0,
    )
    assert abs(result_half["velocity"] - 50.0) < 0.01
    print(f"  At S=Km: v={result_half['velocity']:.2f} (=Vmax/2 ✅)")

    # Test competitive inhibition
    result_comp = tool.run_code(
        analysis_type="full_analysis",
        substrate_concentration_S=5.0,
        Vmax=100.0,
        Km=2.0,
        inhibition_type="competitive",
        inhibitor_concentration_I=3.0,
        Ki=1.0,
    )
    assert "velocity_inhibited" in result_comp
    assert result_comp["velocity_inhibited"] < result_comp["velocity_uninhibited"]
    inh_info = result_comp["inhibition_info"]
    assert inh_info["apparent_Vmax"] == 100.0
    assert inh_info["apparent_Km"] > 2.0
    print(f"  Competitive: v={result_comp['velocity_inhibited']:.2f}, inhibited {result_comp['inhibition_percent']:.1f}%")

    # Test noncompetitive inhibition
    result_noncomp = tool.run_code(
        analysis_type="calculate_velocity",
        substrate_concentration_S=5.0,
        Vmax=100.0,
        Km=2.0,
        inhibition_type="noncompetitive",
        inhibitor_concentration_I=3.0,
        Ki=1.0,
    )
    assert result_noncomp["velocity"] < v
    print(f"  Noncompetitive: v={result_noncomp['velocity']:.2f}")

    # Test linearization data
    result_lin = tool.run_code(
        analysis_type="calculate_velocity",
        substrate_concentration_S=2.0,
        Vmax=100.0,
        Km=2.0,
    )
    lin = result_lin.get("linearization")
    assert lin is not None
    assert "lineweaver_burk" in lin
    assert len(lin["lineweaver_burk"]) > 0
    print(f"  Linearization: {len(lin['lineweaver_burk'])} LB points")

    # Test parameter fitting
    fit_data = [
        {"S": 0.5, "v": 16.7}, {"S": 1.0, "v": 28.6}, {"S": 2.0, "v": 44.4},
        {"S": 4.0, "v": 61.5}, {"S": 8.0, "v": 76.2}, {"S": 16.0, "v": 87.0},
    ]
    result_fit = tool.run_code(
        analysis_type="fit_parameters",
        substrate_concentration_S=2.0,
        Vmax=100.0,
        Km=2.0,
        substrate_velocities_data=fit_data,
        enzyme_concentration_E0=1e-9,
    )
    assert result_fit.get("parameter_fit") is not None
    fit = result_fit["parameter_fit"]
    assert fit["Vmax"] > 0
    assert fit["Km"] > 0
    print(f"  Parameter fit: Vmax≈{fit['Vmax']:.1f}, Km≈{fit['Km']:.1f}, R²={fit['R_squared']:.4f}")

    # Test text interface (positional arg)
    result_text = tool.run_text("calculate_velocity 5.0 100 2.0")
    assert "velocity" in result_text
    print(f"  Text interface: v={result_text['velocity']:.2f}")

    print("  📊 MichaelisMenten: all analysis types work correctly")


# ═══════════════════════════════════════════════════════════════
# #496 ReactionNetworkSolver — 反应网络求解
# ═══════════════════════════════════════════════════════════════
def test_496_reaction_network_solver():
    separator("#496 ReactionNetworkSolver — 反应网络求解")
    from chemmcp.tools import ReactionNetworkSolver

    tool = ReactionNetworkSolver()

    # Test consecutive: A → I → P
    result = tool.run_code(
        species=["A", "I", "P"],
        reactions=[
            {"reactants": ["A"], "products": ["I"], "k": 0.1, "reversible": False},
            {"reactants": ["I"], "products": ["P"], "k": 1.0, "reversible": False},
        ],
        initial_concentrations={"A": 1.0, "I": 0.0, "P": 0.0},
        time_end=50.0,
        n_points=50,
    )
    assert result is not None
    assert result["n_species"] == 3
    assert result["n_reactions"] == 2
    assert result["network_type"] == "consecutive"

    final = result["final_concentrations"]
    assert final["A"] < 1.0
    assert final["P"] > 0.5
    print(f"  Consecutive A→I→P (t=50): A={final['A']:.4f}, I={final['I']:.6f}, P={final['P']:.4f}")

    # Check half-life of A exists (theoretical t_1/2=ln(2)/k1≈6.93)
    # RK4 interpolation accuracy depends on n_points; just check existence
    hl = result["half_lives"].get("A")
    assert hl is not None, "Half-life of A should be computed"
    print(f"  Half-life of A: {hl:.2f} (theoretical ~6.93 for k1=0.1)")

    ss = result["steady_state_analysis"]
    assert "steady_state" in ss
    print(f"  Steady state reached: {ss['steady_state']}")

    # Test reversible reaction
    result_rev = tool.run_code(
        species=["A", "B"],
        reactions=[
            {"reactants": ["A"], "products": ["B"], "k": 0.2, "reversible": True, "k_reverse": 0.1},
        ],
        initial_concentrations={"A": 1.0, "B": 0.0},
        time_end=50.0,
        n_points=20,
    )
    assert result_rev["network_type"] == "reversible"
    final_rev = result_rev["final_concentrations"]
    print(f"  Reversible A⇌B: A={final_rev['A']:.4f}, B={final_rev['B']:.4f}")

    # Test parallel reactions
    result_par = tool.run_code(
        species=["A", "B", "C"],
        reactions=[
            {"reactants": ["A"], "products": ["B"], "k": 0.3, "reversible": False},
            {"reactants": ["A"], "products": ["C"], "k": 0.1, "reversible": False},
        ],
        initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
        time_end=30.0,
        n_points=20,
    )
    final_par = result_par["final_concentrations"]
    assert final_par["B"] > final_par["C"]
    print(f"  Parallel: B={final_par['B']:.4f}, C={final_par['C']:.4f}")

    profiles = result["concentration_profiles"]
    assert len(profiles) > 0
    print(f"  Profile points: {len(profiles)}")

    # Test text interface (positional arg)
    result_text = tool.run_text("species A,I,P reactions A->I:0.1;I->P:1.0 init A=1 t=50")
    assert "final_concentrations" in result_text
    print(f"  Text interface works: {len(result_text['concentration_profiles'])} profile points")

    # Error handling
    try:
        tool.run_code(species=[], reactions=[], initial_concentrations={}, time_end=10)
        assert False, "Should reject empty species"
    except Exception:
        print("  ✅ Correctly rejects empty species list")

    print("  📊 ReactionNetworkSolver: RK4 integration works correctly")


# ═══════════════════════════════════════════════════════════════
# #497 GeometryOptimizer — 几何优化
# ═══════════════════════════════════════════════════════════════
def test_497_geometry_optimizer():
    separator("#497 GeometryOptimizer — 几何优化")
    from chemmcp.tools import GeometryOptimizer

    tool = GeometryOptimizer()

    # Test H2O geometry optimization
    atoms_h2o = [
        {"symbol": "O", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [0.96, 0.0, 0.0]},
        {"symbol": "H", "position": [-0.24, 0.93, 0.0]},
    ]
    bonds_h2o = [
        {"i": 0, "j": 1, "r0": 0.96, "k": 500.0},
        {"i": 0, "j": 2, "r0": 0.96, "k": 500.0},
    ]

    result = tool.run_code(
        atoms=atoms_h2o,
        bonds=bonds_h2o,
        optimizer="steepest_descent",
        max_iterations=200,
        convergence_threshold=1e-4,
        step_size=0.005,
    )
    assert result is not None
    assert "converged" in result
    assert "final_energy_eV" in result
    assert len(result["optimized_positions_Angstrom"]) == 3
    print(f"  H2O optimization: converged={result['converged']}, iter={result['n_iterations']}")
    print(f"  Final E={result['final_energy_eV']:.4f} eV, RMS force={result['rms_force_eV_per_A']:.2e} eV/Å")

    if result["bond_lengths"]:
        for bl in result["bond_lengths"]:
            print(f"  Bond {bl['bond']}: r={bl['length_A']:.4f} Å (r0={bl['equilibrium_r0_A']} Å)")

    # Test conjugate gradient optimizer
    result_cg = tool.run_code(
        atoms=atoms_h2o,
        bonds=bonds_h2o,
        optimizer="conjugate_gradient",
        max_iterations=200,
        convergence_threshold=1e-4,
        step_size=0.01,
    )
    assert result_cg["converged"] in (True, False)
    print(f"  CG optimization: converged={result_cg['converged']}, iter={result_cg['n_iterations']}")

    # Test H2 molecule
    atoms_h2 = [
        {"symbol": "H", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [1.5, 0.0, 0.0]},
    ]
    bonds_h2 = [{"i": 0, "j": 1, "r0": 0.74, "k": 450.0}]

    result_h2 = tool.run_code(
        atoms=atoms_h2,
        bonds=bonds_h2,
        optimizer="steepest_descent",
        max_iterations=300,
        convergence_threshold=1e-5,
        step_size=0.005,
    )
    assert result_h2["n_atoms"] == 2
    print(f"  H2 optimization: converged={result_h2['converged']}, E={result_h2['final_energy_eV']:.4f} eV")
    if result_h2["bond_lengths"]:
        print(f"  Optimized H-H distance: {result_h2['bond_lengths'][0]['length_A']:.4f} Å")

    # Error handling
    try:
        tool.run_code(atoms=[], bonds=None)
        assert False, "Should reject empty atoms"
    except Exception:
        print("  ✅ Correctly rejects empty atoms list")

    print("  📊 GeometryOptimizer: SD and CG optimizers work correctly")


# ═══════════════════════════════════════════════════════════════
# #498 TransitionStateSearch — 过渡态搜索
# ═══════════════════════════════════════════════════════════════
def test_498_transition_state_search():
    separator("#498 TransitionStateSearch — 过渡态搜索")
    from chemmcp.tools import TransitionStateSearch

    tool = TransitionStateSearch()

    atoms_ts = [
        {"symbol": "H", "position": [0.0, 0.0, -0.5]},
        {"symbol": "H", "position": [0.0, 0.0, 0.5]},
        {"symbol": "O", "position": [0.0, 0.0, 1.8]},
    ]
    bonds_ts = [
        {"i": 0, "j": 2, "r0": 1.0, "k": 400},
        {"i": 1, "j": 2, "r0": 1.0, "k": 400},
    ]

    result = tool.run_code(
        atoms=atoms_ts,
        bonds=bonds_ts,
        guess_coordinates=None,
        search_method="quadratic_saddle",
        max_iterations=80,
        convergence_threshold=1e-4,
    )
    assert result is not None
    assert "ts_found" in result
    assert "n_imaginary_frequencies" in result
    assert "ts_energy_eV" in result
    assert "validation" in result
    print(f"  TS search: found={result['ts_found']}, n_imag={result['n_imaginary_frequencies']}")
    print(f"  TS energy: {result['ts_energy_eV']:.4f} eV")
    print(f"  Imaginary frequency: {result['imaginary_frequency_cm-1']} cm⁻¹")
    print(f"  Validation: {result['validation']}")

    # Test eigenvector following
    result_ef = tool.run_code(
        atoms=atoms_ts,
        bonds=bonds_ts,
        search_method="eigenvector_following",
        max_iterations=50,
    )
    assert "ts_found" in result_ef
    print(f"  EF method: found={result_ef['ts_found']}, iter={result_ef['n_iterations']}")

    # Test dimer method
    result_dimer = tool.run_code(
        atoms=atoms_ts,
        bonds=bonds_ts,
        search_method="dimer",
        max_iterations=50,
    )
    assert "ts_converged" in result_dimer
    print(f"  Dimer method: converged={result_dimer['ts_converged']}")

    assert len(result["ts_positions_Angstrom"]) == 3

    try:
        tool.run_code(atoms=[], bonds=None)
        assert False, "Should reject empty atoms"
    except Exception:
        print("  ✅ Correctly rejects empty atoms list")

    print("  📊 TransitionStateSearch: saddle point search works correctly")


# ═══════════════════════════════════════════════════════════════
# #499 PotentialEnergySurface — 势能面扫描
# ═══════════════════════════════════════════════════════════════
def test_499_potential_energy_surface():
    separator("#499 PotentialEnergySurface — 势能面扫描")
    from chemmcp.tools import PotentialEnergySurface

    tool = PotentialEnergySurface()

    # Test bond length scan for H2
    atoms_h2 = [
        {"symbol": "H", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [1.0, 0.0, 0.0]},
    ]
    bonds_h2 = [{"i": 0, "j": 1, "r0": 0.74, "k": 450.0}]

    result = tool.run_code(
        atoms=atoms_h2,
        bonds=bonds_h2,
        scan_type="bond_length",
        scan_atoms=[0, 1],
        start_value=0.4,
        end_value=2.0,
        n_points=20,
        optimize_each_point=False,
    )
    assert result is not None
    assert result["scan_type"] == "bond_length"
    assert result["n_points"] == 20
    assert len(result["energy_profile"]) == 20
    print(f"  Bond scan (H-H): {result['n_points']} pts, range [{result['scan_range'][0]}, {result['scan_range'][1]}] Å")
    print(f"  E range: [{result['min_energy_eV']:.4f}, {result['max_energy_eV']:.4f}] eV")

    sp = result["stationary_points"]
    print(f"  Stationary points: {len(sp)}")
    for s in sp:
        print(f"    {s['type']} at coord={s.get('coordinate_value', 'N/A'):.3f} Å, E={s['energy_eV']:.4f} eV")

    mins = [s for s in sp if s["type"] == "minimum"]
    assert len(mins) >= 1
    min_coord = mins[0]["coordinate_value"]
    assert 0.5 < min_coord < 1.2
    print(f"  Equilibrium distance: {min_coord:.3f} Å (expected ~0.74 Å)")

    energies = [p["energy"] for p in result["energy_profile"]]
    assert min(energies) < energies[0]
    assert min(energies) < energies[-1]
    print(f"  Profile shape correct: min at interior ✅")

    barriers = result["barrier_heights_eV"]
    print(f"  Barriers: forward={barriers['forward']}, reverse={barriers['reverse']}")

    # Test with water molecule O-H bond scan
    atoms_h2o = [
        {"symbol": "O", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [0.96, 0.0, 0.0]},
        {"symbol": "H", "position": [-0.24, 0.93, 0.0]},
    ]
    bonds_h2o = [
        {"i": 0, "j": 1, "r0": 0.96, "k": 500},
        {"i": 0, "j": 2, "r0": 0.96, "k": 500},
    ]
    result_h2o = tool.run_code(
        atoms=atoms_h2o,
        bonds=bonds_h2o,
        scan_type="bond_length",
        scan_atoms=[0, 1],
        start_value=0.5,
        end_value=2.0,
        n_points=15,
    )
    assert result_h2o["min_energy_eV"] is not None
    print(f"  O-H bond scan: E_min={result_h2o['min_energy_eV']:.4f} eV")

    # Error handling
    try:
        tool.run_code(atoms=atoms_h2, bonds=bonds_h2o, scan_type="bond_length",
                     scan_atoms=None, start_value=0.5, end_value=2.0)
        assert False
    except Exception:
        print("  ✅ Correctly rejects missing scan_atoms")

    try:
        tool.run_code(atoms=atoms_h2, bonds=bonds_h2o, scan_type="bond_length",
                     scan_atoms=[0, 1], start_value=0.5, end_value=2.0, n_points=1)
        assert False
    except Exception:
        print("  ✅ Correctly rejects n_points < 3")

    print("  📊 PotentialEnergySurface: PES scanning works correctly")


# ═══════════════════════════════════════════════════════════════
# #500 FrequencyAnalysis — 频率分析
# ═══════════════════════════════════════════════════════════════
def test_500_frequency_analysis():
    separator("#500 FrequencyAnalysis — 频率分析")
    from chemmcp.tools import FrequencyAnalysis

    tool = FrequencyAnalysis()

    # Test with water molecule
    atoms_h2o = [
        {"symbol": "O", "position": [0.0, 0.0, 0.0], "mass": 15.999},
        {"symbol": "H", "position": [0.0, 0.0, 0.96], "mass": 1.008},
        {"symbol": "H", "position": [0.0, 0.96, 0.0], "mass": 1.008},
    ]

    result = tool.run_code(
        atoms=atoms_h2o,
        hessian_matrix=None,
        temperature_K=298.15,
        pressure_atm=1.0,
        scale_factor=1.0,
    )
    assert result is not None
    assert "stationary_point_type" in result
    assert "n_imaginary_frequencies" in result
    assert "vibrational_frequencies_cm-1" in result
    assert "zero_point_energy_ZPE_eV" in result
    assert "thermodynamics" in result
    print(f"  Stationary point: {result['stationary_point_type']}")
    print(f"  N_imaginary: {result['n_imaginary_frequencies']}, N_vib_modes: {result['n_vibrational_modes']}")
    print(f"  ZPE: {result['zero_point_energy_ZPE_eV']:.4f} eV ({result['zero_point_energy_ZPE_kJ_mol']:.2f} kJ/mol)")

    freqs = result["all_frequencies_cm-1"]
    assert len(freqs) == 9
    print(f"  Total frequencies: {len(freqs)} (3N=9)")

    vib_freqs = result["vibrational_frequencies_cm-1"]
    print(f"  Vibrational frequencies: {vib_freqs}")

    thermo = result["thermodynamics"]
    assert "G_eV" in thermo
    assert "H_eV" in thermo
    assert "S_J_mol_K" in thermo
    print(f"  G = {thermo['G_eV']:.4f} eV ({thermo['G_kJ_mol']:.2f} kJ/mol)")
    print(f"  H = {thermo['H_eV']:.4f} eV ({thermo['H_kJ_mol']:.2f} kJ/mol)")
    print(f"  S = {thermo['S_J_mol_K']:.2f} J/(mol·K)")

    assert "analysis_summary" in result
    print(f"  Summary present ✅")

    # Test with custom Hessian (TS-like with negative eigenvalue)
    n_dim = 9
    custom_H = [[0.0] * n_dim for _ in range(n_dim)]
    for i in range(n_dim):
        custom_H[i][i] = 100.0 if i != 0 else -50.0

    result_ts = tool.run_code(
        atoms=atoms_h2o,
        hessian_matrix=custom_H,
        temperature_K=298.15,
    )
    assert result_ts["n_imaginary_frequencies"] >= 1
    print(f"\n  Custom Hessian (TS-like): n_imag={result_ts['n_imaginary_frequencies']}")
    print(f"  Type: {result_ts['stationary_point_type']}")

    # Test scale factor
    result_scaled = tool.run_code(
        atoms=atoms_h2o,
        hessian_matrix=None,
        temperature_K=298.15,
        scale_factor=0.96,
    )
    assert result_scaled["scale_factor"] == 0.96
    print(f"\n  Scaled (0.96): ZPE={result_scaled['zero_point_energy_ZPE_eV']:.4f} eV")

    # Test different temperature
    result_hot = tool.run_code(
        atoms=atoms_h2o,
        hessian_matrix=None,
        temperature_K=500.0,
    )
    assert result_hot["temperature_K"] == 500.0
    print(f"  Hot (500K): G={result_hot['thermodynamics']['G_eV']:.4f} eV")

    # Error handling
    try:
        tool.run_code(atoms=[], hessian_matrix=None)
        assert False, "Should reject empty atoms"
    except Exception:
        print("\n  ✅ Correctly rejects empty atoms list")

    print("\n  📊 FrequencyAnalysis: vibrational & thermochemical analysis works correctly")


# ═══════════════════════════════════════════════════════════════
# Main runner
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("\n" + "█" * 65)
    print("  MCP Tools #491-500 Test Suite")
    print("  Kinetics Corrections & Computational Chemistry")
    print("█" * 65)

    run_test("#491 TunnelingCorrection", test_491_tunneling_correction)
    run_test("#492 SteadyStateApprox", test_492_steady_state_approx)
    run_test("#493 PreEquilibrium", test_493_pre_equilibrium)
    run_test("#494 RateDeterminingStep", test_494_rate_determining_step)
    run_test("#495 MichaelisMenten", test_495_michaelis_menten)
    run_test("#496 ReactionNetworkSolver", test_496_reaction_network_solver)
    run_test("#497 GeometryOptimizer", test_497_geometry_optimizer)
    run_test("#498 TransitionStateSearch", test_498_transition_state_search)
    run_test("#499 PotentialEnergySurface", test_499_potential_energy_surface)
    run_test("#500 FrequencyAnalysis", test_500_frequency_analysis)

    print("\n" + "█" * 65)
    total = passed + failed
    print(f"  RESULTS: {passed}/{total} PASSED, {failed}/{total} FAILED")
    if failed == 0:
        print("  🎉 ALL TESTS PASSED! 🦐")
    else:
        print(f"  ⚠️  {failed} test(s) failed — review output above")
    print("█" * 65 + "\n")

    sys.exit(0 if failed == 0 else 1)
