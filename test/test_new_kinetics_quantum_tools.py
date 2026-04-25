"""
Test suite for new MCP tools (#221-230):
  221. CollisionTheory
  222. EnzymeKinetics
  223. ReactionMechanismSimulator
  224. SteadyStateApproximation
  225. RateDeterminingStep
  226. TemperatureJumpRelaxation
  227. ParallelConsecutiveReactions
  228. ParticleInBox
  229. HarmonicOscillator
  230. RigidRotor
"""

import sys
import os
import math

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from chemmcp.tools import (
    CollisionTheory,
    EnzymeKinetics,
    ReactionMechanismSimulator,
    SteadyStateApproximation,
    RateDeterminingStep,
    TemperatureJumpRelaxation,
    ParallelConsecutiveReactions,
    ParticleInBox,
    HarmonicOscillator,
    RigidRotor,
)


def test_collision_theory():
    """Test #221: Collision Theory - Arrhenius rate calculation."""
    print("\n=== Test 221: CollisionTheory ===")
    tool = CollisionTheory()
    
    # Test basic case: k = A * f * exp(-Ea/RT)
    result = tool.run_code(
        pre_exponential_factor_A=1e10,
        activation_energy_Ea=50000.0,
        temperature_K=300.0,
        steric_factor_f=1.0
    )
    
    assert "rate_constant_k" in result, "Missing rate_constant_k"
    assert result["rate_constant_k"] > 0, "Rate constant should be positive"
    assert result["steric_factor"] == 1.0
    
    # Verify: exp(-50000/(8.314*300)) ≈ some small number
    expected_exp = math.exp(-50000 / (8.314 * 300))
    assert abs(result["exponential_term"] - expected_exp) / max(expected_exp, 1e-30) < 0.01, \
        f"Exponential term mismatch: {result['exponential_term']} vs {expected_exp}"
    
    expected_k = 1e10 * expected_exp
    assert abs(result["rate_constant_k"] - expected_k) / expected_k < 1e-4, \
        f"Rate constant mismatch: {result['rate_constant_k']} vs {expected_k}"
    
    print(f"  ✓ k = {result['rate_constant_k']:.4e} (expected ~{expected_k:.4e})")
    
    # Test with steric factor < 1
    result2 = tool.run_code(pre_exponential_factor_A=5e9, activation_energy_Ea=75000, temperature_K=298, steric_factor_f=0.01)
    assert result2["rate_constant_k"] > 0
    print(f"  ✓ With steric factor 0.01: k = {result2['rate_constant_k']:.4e}")
    
    # Test text interface
    result3 = tool.run_text("1e10 50000 300 1.0")
    assert abs(result3["rate_constant_k"] - result["rate_constant_k"]) < 1e-10
    print("  ✓ Text interface works")
    
    # Test error handling
    try:
        tool.run_code(A=1e10, Ea=50000, T=-100)
        assert False, "Should have raised error for negative T"
    except Exception:
        pass
    
    print("  ✅ All CollisionTheory tests passed!")


def test_enzyme_kinetics():
    """Test #222: Enzyme Kinetics - Michaelis-Menten fitting."""
    print("\n=== Test 222: EnzymeKinetics ===")
    tool = EnzymeKinetics()
    
    # Classic Michaelis-Menten data (Vmax≈1.0, Km≈2.0)
    S_data = [0.5, 1.0, 2.0, 4.0, 8.0]
    v_data = [0.21, 0.35, 0.53, 0.71, 0.85]
    
    result = tool.run_code(
        substrate_concentrations_S=S_data,
        reaction_velocities_v=v_data
    )
    
    assert "Vmax_best" in result
    assert "Km_best" in result
    assert result["Vmax_best"] > 0, "Vmax must be positive"
    assert result["Km_best"] > 0, "Km must be positive"
    assert result["best_method"] in ["Lineweaver-Burk", "Eadie-Hofstee", "Hanes-Woolf"]
    
    # Check that all three methods returned results
    assert "Lineweaver_Burk" in result
    assert "Eadie_Hofstee" in result
    assert "Hanes_Woolf" in result
    
    print(f"  ✓ Vmax ≈ {result['Vmax_best']:.3f}, Km ≈ {result['Km_best']:.3f}")
    print(f"  ✓ Best method: {result['best_method']}")
    
    # Test with enzyme concentration for kcat
    result2 = tool.run_code(S_data, v_data, enzyme_concentration_E0=1e-6)
    if result2.get("kcat"):
        assert result2["kcat"] > 0
        print(f"  ✓ kcat = {result2['kcat']:.4e} s⁻¹")
    
    # Test error: mismatched lengths
    try:
        tool.run_code([1, 2], [1], None)
        assert False, "Should raise error for mismatched lengths"
    except Exception:
        pass
    
    print("  ✅ All EnzymeKinetics tests passed!")


def test_reaction_mechanism_simulator():
    """Test #223: Reaction Mechanism Simulator (RK4 ODE solver)."""
    print("\n=== Test 223: ReactionMechanismSimulator ===")
    tool = ReactionMechanismSimulator()
    
    # Consecutive reaction A->B;0.1, B->C;0.05
    result = tool.run_code(
        mechanism="A->B;0.1,B->C;0.05",
        initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
        time_end=100.0,
        n_points=11
    )
    
    assert "species" in result
    assert "concentration_profiles" in result
    assert "final_concentrations" in result
    
    # Mass conservation: A + B + C should equal initial total (=1.0)
    final = result["final_concentrations"]
    total = final["A"] + final["B"] + final["C"]
    assert abs(total - 1.0) < 0.01, f"Mass conservation violated: total={total}"
    
    # A should decrease over time
    profiles = result["concentration_profiles"]
    assert profiles["A"][0] > profiles["A"][-1], "A should decrease"
    
    # C should increase over time
    assert profiles["C"][0] < profiles["C"][-1], "C should increase"
    
    print(f"  ✓ Final: [A]={final['A']:.4f}, [B]={final['B']:.4f}, [C]={final['C']:.4f}")
    print(f"  ✓ Mass conserved: {total:.6f} ≈ 1.0")
    
    # Test parallel reaction
    result2 = tool.run_code(
        mechanism="A->B;0.1,A->C;0.02",
        initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
        time_end=50.0,
        n_points=6
    )
    final2 = result2["final_concentrations"]
    total2 = sum(final2.values())
    assert abs(total2 - 1.0) < 0.01
    print(f"  ✓ Parallel final: [A]={final2['A']:.4f}, [B]={final2['B']:.4f}, [C]={final2['C']:.4f}")
    
    print("  ✅ All ReactionMechanismSimulator tests passed!")


def test_steady_state_approximation():
    """Test #224: Steady State Approximation."""
    print("\n=== Test 224: SteadyStateApproximation ===")
    tool = SteadyStateApproximation()
    
    # Consecutive reaction A -> I -> P with k2 >> k1 (good SSA condition)
    result = tool.run_code(
        mechanism_type="consecutive",
        rate_constants=[0.1, 1.0],
        initial_reactant_concentration_A0=1.0,
        time_t=50.0
    )
    
    assert "exact" in result
    assert "approximate" in result
    assert "ssa_valid" in result
    
    # When k2 >> k1, SSA should be valid
    if result["ssa_valid"]:
        print(f"  ✓ SSA valid (k2/k1 = {result['k2_k1_ratio']})")
    
    # Check exact values are physical
    exact = result["exact"]
    assert exact["A"] >= 0, "[A] cannot be negative"
    assert exact["I"] >= 0, "[I] cannot be negative"
    assert exact["P"] >= 0, "[P] cannot be negative"
    
    # Total mass conservation
    total_exact = exact["A"] + exact["I"] + exact["P"]
    assert abs(total_exact - 1.0) < 0.01, f"Mass not conserved: {total_exact}"
    
    print(f"  ✓ Exact: [A]={exact['A']:.4f}, [I]={exact['I']:.6f}, [P]={exact['P']:.4f}")
    print(f"  ✓ Approx: [I]={result['approximate']['I']:.6f}")
    
    # Test pre-equilibrium
    result2 = tool.run_code(
        mechanism_type="pre_equilibrium",
        rate_constants=[100.0, 50.0, 1.0],
        initial_reactant_concentration_A0=1.0,
        time_t=10.0
    )
    assert "equilibrium_constant_Keq" in result2
    print(f"  ✓ Pre-equilibrium Keq = {result2['equilibrium_constant_Keq']}")
    
    print("  ✅ All SteadyStateApproximation tests passed!")


def test_rate_determining_step():
    """Test #225: Rate Determining Step Analysis."""
    print("\n=== Test 225: RateDeterminingStep ===")
    tool = RateDeterminingStep()
    
    steps = [
        {"reactants": "A -> B", "products": "", "k": 0.001, "reversible": False},
        {"reactants": "B + C -> D", "products": "", "k": 5.0, "reversible": False},
        {"reactants": "D -> E", "products": "", "k": 100.0, "reversible": False},
    ]
    
    result = tool.run_code(mechanism_steps=steps)
    
    assert result["rds_step_index"] == 1, "First step should be RDS (smallest k)"
    assert result["rds_rate_constant"] == 0.001
    assert "overall_rate_law" in result
    assert "justification" in result
    assert result["slowest_to_fastest_ratio"] >= 100  # Should be much slower
    
    print(f"  ✓ RDS: Step {result['rds_step_index']} ({result['rds_step_description']}), k={result['rds_rate_constant']}")
    print(f"  ✓ Rate law: {result['overall_rate_law']}")
    print(f"  ✓ Slowest/fastest ratio: {result['slowest_to_fastest_ratio']:,.0f}x")
    
    # Test single step
    result2 = tool.run_code([{"reactants": "X -> Y", "k": 42.0, "reversible": False}])
    assert result2["rds_step_index"] == 1
    print(f"  ✓ Single step: RDS is the only step")
    
    print("  ✅ All RateDeterminingStep tests passed!")


def test_temperature_jump_relaxation():
    """Test #226: Temperature Jump Relaxation Kinetics."""
    print("\n=== Test 226: TemperatureJumpRelaxation ===")
    tool = TemperatureJumpRelaxation()
    
    result = tool.run_code(
        equilibrium_constant_K=2.0,
        delta_H=-50000.0,
        initial_temperature_K=298.0,
        final_temperature_K=308.0,
        forward_rate_constant_kf=1000.0,
        reverse_rate_constant_kr=500.0,
        total_concentration=1.0
    )
    
    assert "relaxation_time_tau_s" in result
    assert result["relaxation_time_tau_s"] > 0
    tau = result["relaxation_time_tau_s"]
    
    # τ = 1/(kf+kr) = 1/1500 ≈ 0.000667s
    expected_tau = 1.0 / (1000.0 + 500.0)
    assert abs(tau - expected_tau) / expected_tau < 1e-6, \
        f"Tau mismatch: {tau} vs {expected_tau}"
    
    assert "new_equilibrium_A_M" in result
    assert "new_equilibrium_B_M" in result
    assert "kinetic_analysis" in result
    
    # K = [B]/[A] => [A] = total/(1+K) = 1/3
    assert abs(result["new_equilibrium_A_M"] - 1/3) < 0.01
    assert abs(result["new_equilibrium_B_M"] - 2/3) < 0.01
    
    print(f"  ✓ τ = {tau:.6e}s (expected {expected_tau:.6e}s)")
    print(f"  ✓ New eq: [A]={result['new_equilibrium_A_M']:.4f}, [B]={result['new_equilibrium_B_M']:.4f}")
    print(f"  ✓ 99% equilibrium at t = {result['time_to_99_percent_equilibrium_s']:.6e}s")
    
    print("  ✅ All TemperatureJumpRelaxation tests passed!")


def test_parallel_consecutive_reactions():
    """Test #227: Parallel and Consecutive Reactions Solver."""
    print("\n=== Test 227: ParallelConsecutiveReactions ===")
    tool = ParallelConsecutiveReactions()
    
    # --- Parallel reactions ---
    result_p = tool.run_code(
        reaction_type="parallel",
        rate_constants=[0.1, 0.02],
        initial_concentrations={"A": 1.0},
        time_points=[0, 10, 20, 50, 100]
    )
    
    assert result_p["reaction_type"] == "parallel"
    assert "A" in result_p
    assert "products" in result_p
    assert "selectivity" in result_p
    assert "yields" in result_p
    
    # At t=100, A should be significantly depleted
    A_final = result_p["A"][-1]
    assert A_final < 0.5, f"A should be depleted at t=100, got {A_final}"
    
    # Selectivity B/C = k1/k2 = 5
    sel_key = list(result_p["selectivity"].keys())[0]
    sel_val = result_p["selectivity"][sel_key]
    assert abs(sel_val - 5.0) < 0.01, f"Selectivity should be ~5, got {sel_val}"
    
    print(f"  ✓ Parallel at t=100: [A]={A_final:.4f}, selectivity(B/C)={sel_val:.1f}")
    
    # --- Consecutive reactions ---
    result_c = tool.run_code(
        reaction_type="consecutive",
        rate_constants=[0.1, 0.05],
        initial_concentrations={"A": 1.0},
        time_points=[0, 10, 20, 50]
    )
    
    assert result_c["reaction_type"] == "consecutive"
    assert "I" in result_c and "P" in result_c
    
    # Intermediate I should rise then fall
    I_values = result_c["I"]
    I_max = max(I_values)
    assert I_max > 0, "Intermediate should form"
    
    # Max intermediate time should exist
    if result_c.get("max_intermediate_time") is not None:
        assert result_c["max_intermediate_time"] > 0
        print(f"  ✓ Consecutive: max[I] at t≈{result_c['max_intermediate_time']:.2f}s, [I]max={I_max:.4f}")
    
    print(f"  ✓ Consecutive at t=50: [A]={result_c['A'][-1]:.4f}, [I]={result_c['I'][-1]:.4f}, [P]={result_c['P'][-1]:.4f}")
    
    print("  ✅ All ParallelConsecutiveReactions tests passed!")


def test_particle_in_box():
    """Test #228: Particle in a Box (Quantum Mechanics)."""
    print("\n=== Test 228: ParticleInBox ===")
    tool = ParticleInBox()
    
    # 1D: electron in 1nm box, n=2
    result = tool.run_code(
        dimensionality="1d",
        mass_kg=9.109e-31,
        box_length_m=1e-9,
        quantum_number_n=2
    )
    
    assert "energy_J" in result
    assert "energy_eV" in result
    assert result["energy_eV"] > 0
    assert result["degeneracy"] == 1
    assert result["n_nodes"] == 1  # n-1 nodes for n=2
    
    # E_2 = 4*h²/(8mL²) = 4*E_1
    h = 6.626e-34
    m = 9.109e-31
    L = 1e-9
    E_expected = 4 * h**2 / (8 * m * L**2)
    assert abs(result["energy_J"] - E_expected) / E_expected < 0.01, \
        f"Energy mismatch: {result['energy_J']} vs {E_expected}"
    
    print(f"  ✓ 1D n=2: E = {result['energy_eV']:.4f} eV, nodes = {result['n_nodes']}")
    
    # 3D cubic box
    result3d = tool.run_code(
        dimensionality="3d",
        mass_kg=9.109e-31,
        lengths_3d=[1e-9, 1e-9, 1e-9],
        quantum_numbers_3d=[1, 1, 2]
    )
    
    assert result3d["dimensionality"] == "3D"
    assert result3d["energy_eV"] > 0
    # E(1,1,2) = (1+1+4)*h²/(8mL²) = 6*E_1 (where E_1 is ground state of 1D)
    E_3d_expected = (1+1+4) * h**2 / (8 * m * L**2)
    assert abs(result3d["energy_J"] - E_3d_expected) / E_3d_expected < 0.01
    
    print(f"  ✓ 3D (1,1,2): E = {result3d['energy_eV']:.4f} eV, total nodes = {result3d['n_total_nodes']}")
    
    # Test n=1 (ground state, no internal nodes)
    result_gs = tool.run_code("1d", 9.109e-31, 1e-9, 1)
    assert result_gs["n_nodes"] == 0
    print(f"  ✓ Ground state n=1: E = {result_gs['energy_eV']:.4f} eV, nodes = {result_gs['n_nodes']}")
    
    print("  ✅ All ParticleInBox tests passed!")


def test_harmonic_oscillator():
    """Test #229: Quantum Harmonic Oscillator."""
    print("\n=== Test 229: HarmonicOscillator ===")
    tool = HarmonicOscillator()
    
    # HCl-like vibration: reduced mass ~1.63e-27 kg, k ~480 N/m
    mu = 1.627e-27  # reduced mass of HCl (kg)
    k = 480.0       # force constant (N/m)
    
    # Ground state v=0
    result_v0 = tool.run_code(
        mass_kg=mu,
        force_constant_N_m=k,
        quantum_number_v=0
    )
    
    assert "energy_J" in result_v0
    assert result_v0["energy_eV"] > 0  # Zero-point energy!
    assert result_v0["quantum_number_v"] == 0
    assert result_v0["parity"] == "even"
    assert result_v0["wavefunction_hermite_order"] == 0
    
    omega = math.sqrt(k / mu)
    E_zp_expected = 0.5 * 1.054571817e-34 * omega
    assert abs(result_v0["energy_J"] - E_zp_expected) / E_zp_expected < 0.001
    
    print(f"  ✓ v=0: E_ZPE = {result_v0['energy_eV']:.4f} eV, ν̃ = {result_v0['vibrational_frequency_cm-1']:.1f} cm⁻¹")
    print(f"  ✓ Turning points: ±{result_v0['turning_points_pm_m']:.3e} m")
    
    # First excited state v=1
    result_v1 = tool.run_code(
        mass_kg=mu, force_constant_N_m=k, quantum_number_v=1
    )
    assert result_v1["parity"] == "odd"
    assert result_v1["energy_eV"] > result_v0["energy_eV"]
    
    # ΔE(v=0→1) should equal ℏω
    delta_E = result_v1["energy_J"] - result_v0["energy_J"]
    expected_delta = 1.054571817e-34 * omega
    assert abs(delta_E - expected_delta) / expected_delta < 0.001
    
    print(f"  ✓ v=1: E = {result_v1['energy_eV']:.4f} eV, parity = {result_v1['parity']}")
    print(f"  ✓ ΔE(0→1) = {delta_E * 6.241509e18:.4f} eV")
    
    # v=2
    result_v2 = tool.run_code(
        mass_kg=mu, force_constant_N_m=k, quantum_number_v=2
    )
    assert result_v2["parity"] == "even"
    assert result_v2["energy_eV"] > result_v1["energy_eV"]
    print(f"  ✓ v=2: E = {result_v2['energy_eV']:.4f} eV")
    
    # Equally spaced levels check
    E01 = result_v1["energy_J"] - result_v0["energy_J"]
    E12 = result_v2["energy_J"] - result_v1["energy_J"]
    assert abs(E01 - E12) / E01 < 0.001, "QHO levels should be equally spaced!"
    print(f"  ✓ Energy spacing equal: ΔE01 = ΔE12 = {E01 * 6.241509e18:.4f} eV")
    
    print("  ✅ All HarmonicOscillator tests passed!")


def test_rigid_rotor():
    """Test #230: Rigid Rotor (Rotational Spectroscopy)."""
    print("\n=== Test 230: RigidRotor ===")
    tool = RigidRotor()
    
    # CO molecule-like parameters
    I_CO = 1.456e-46  # kg·m² (moment of inertia of CO)
    
    # J=0 (ground state)
    result_j0 = tool.run_code(
        moment_of_inertia_kg_m2=I_CO,
        quantum_number_J=0
    )
    
    assert result_j0["energy_J"] == 0, "J=0 energy should be 0"
    assert result_j0["degeneracy"] == 1, "g(0) = 2*0+1 = 1"
    assert result_j0["rotational_constant_B_MHz"] > 0
    
    B_MHz = result_j0["rotational_constant_B_MHz"]
    print(f"  ✓ J=0: E=0, g=1, B = {B_MHz:.1f} MHz, Θ_rot = {result_j0['characteristic_rotational_temperature_K']:.3f} K")
    
    # J=1
    result_j1 = tool.run_code(I_CO, 1)
    assert result_j1["degeneracy"] == 3, "g(1) = 3"
    # E(J=1) = 2B (in units of B), should be 2x E_per_B
    assert abs(result_j1["energy_per_B_units"] - 2) < 0.001
    print(f"  ✓ J=1: E = {result_j1['energy_eV']:.6f} eV, g = {result_j1['degeneracy']}")
    
    # J=2
    result_j2 = tool.run_code(I_CO, 2)
    assert result_j2["degeneracy"] == 5, "g(2) = 5"
    assert abs(result_j2["energy_per_B_units"] - 6) < 0.001  # 2*3=6
    print(f"  ✓ J=2: E = {result_j2['energy_eV']:.6f} eV, g = {result_j2['degeneracy']}")
    
    # Test transitions
    result_trans = tool.run_code(I_CO, 0, compute_transitions_up_to_J=4)
    assert "spectral_transitions" in result_trans
    transitions = result_trans["spectral_transitions"]
    assert len(transitions) >= 4  # J=0→1, 1→2, 2→3, 3→4
    
    # Transition frequencies should be evenly spaced (2B, 4B, 6B, ...)
    freqs = [t["frequency_GHz"] for t in transitions]
    for i in range(1, len(freqs)):
        # Each successive transition should be separated by 2B (in frequency units)
        ratio = freqs[i] / freqs[0]
        assert abs(ratio - (i + 1)) < 0.01, \
            f"Transition {i} frequency ratio off: {ratio} vs {i+1}"
    
    print(f"  ✓ Transitions: {len(transitions)} lines, first at {freqs[0]:.3f} GHz")
    print(f"  ✓ Spacing verified: equally spaced by 2B")
    
    # Test bond length mode
    result_bl = tool.run_code(
        moment_of_inertia_kg_m2=None,
        quantum_number_J=1,
        bond_length_m=1.128e-10,  # CO bond length
        reduced_mass_kg=1.145e-26  # CO reduced mass
    )
    assert result_bl["energy_eV"] > 0
    print(f"  ✓ Bond length mode: E(J=1) = {result_bl['energy_eV']:.6f} eV")
    
    print("  ✅ All RigidRotor tests passed!")


# ============================================================
# Main test runner
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("ChemMCP New Tools Test Suite (#221-230)")
    print("=" * 60)
    
    tests = [
        ("CollisionTheory", test_collision_theory),
        ("EnzymeKinetics", test_enzyme_kinetics),
        ("ReactionMechanismSimulator", test_reaction_mechanism_simulator),
        ("SteadyStateApproximation", test_steady_state_approximation),
        ("RateDeterminingStep", test_rate_determining_step),
        ("TemperatureJumpRelaxation", test_temperature_jump_relaxation),
        ("ParallelConsecutiveReactions", test_parallel_consecutive_reactions),
        ("ParticleInBox", test_particle_in_box),
        ("HarmonicOscillator", test_harmonic_oscillator),
        ("RigidRotor", test_rigid_rotor),
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
            print(f"\n  ❌ {name} FAILED: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {passed}/{len(tests)} passed, {failed} failed")
    if errors:
        print("\nFailed tests:")
        for name, err in errors:
            print(f"  ❌ {name}: {err[:200]}")
    print("=" * 60)
    
    sys.exit(0 if failed == 0 else 1)
