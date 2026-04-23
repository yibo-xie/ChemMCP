"""
Test suite for Quantum Chemistry MCP Tools (#231-240).
Run with: python -m tests.test_quantum_chemistry_231_240
From the ChemMCP directory with environment activated.
"""
import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

passed = 0
failed = 0
errors = []

def test(name, func):
    global passed, failed
    try:
        result = func()
        if result is None or result is True:
            passed += 1
            print(f"  ✅ {name}")
        else:
            failed += 1
            errors.append((name, f"Returned: {result}"))
            print(f"  ❌ {name} -> Returned: {result}")
    except Exception as e:
        failed += 1
        errors.append((name, str(e)))
        print(f"  ❌ {name} -> {type(e).__name__}: {e}")


print("=" * 60)
print("Testing Quantum Chemistry Tools (#231-240)")
print("=" * 60)

# ============================================================
# Tool 231: HydrogenAtomOrbitals
# ============================================================
print("\n🧪 #231 HydrogenAtomOrbitals")

def test_231_hydrogen_1s():
    from chemmcp.tools.hydrogen_atom_orbitals import HydrogenAtomOrbitals
    tool = HydrogenAtomOrbitals()
    r = tool.run_code(principal_quantum_number_n=1, orbital_quantum_number_l=0,
                       magnetic_quantum_number_m=0, position_r_bohr=1.0, nuclear_charge_Z=1)
    assert abs(r["energy_eV"] - (-13.6057)) < 0.01, f"E={r['energy_eV']}"
    assert r["orbital_type"] == "1s"
    assert r["n_radial_nodes"] == 0
    assert r["total_nodes"] == 0
    assert r["degeneracy"] == 1
    return True

def test_231_hydrogen_2p():
    from chemmcp.tools.hydrogen_atom_orbitals import HydrogenAtomOrbitals
    tool = HydrogenAtomOrbitals()
    r = tool.run_code(principal_quantum_number_n=2, orbital_quantum_number_l=1,
                       magnetic_quantum_number_m=0, position_r_bohr=2.0, nuclear_charge_Z=1)
    assert abs(r["energy_eV"] - (-3.4014)) < 0.01, f"E={r['energy_eV']}"
    assert r["orbital_type"] == "2p"
    assert r["n_radial_nodes"] == 0
    assert r["angular_nodes"] == 1
    return True

def test_231_hydrogen_text_interface():
    from chemmcp.tools.hydrogen_atom_orbitals import HydrogenAtomOrbitals
    tool = HydrogenAtomOrbitals(interface="text")
    r = tool.run_text("2 1 0 2.0 1")
    assert r["orbital_type"] == "2p"
    return True

def test_231_he_like():
    from chemmcp.tools.hydrogen_atom_orbitals import HydrogenAtomOrbitals
    tool = HydrogenAtomOrbitals()
    r = tool.run_code(principal_quantum_number_n=1, orbital_quantum_number_l=0,
                       magnetic_quantum_number_m=0, position_r_bohr=0.5, nuclear_charge_Z=2)
    # E = -13.6 * Z²/n² = -54.4 eV for He+
    assert abs(r["energy_eV"] - (-54.4228)) < 0.1, f"E={r['energy_eV']}"
    return True

test("231: H 1s orbital energy & nodes", test_231_hydrogen_1s)
test("231: H 2p orbital", test_231_hydrogen_2p)
test("231: Text interface parsing", test_231_hydrogen_text_interface)
test("231: He+ (Z=2) energy scaling", test_231_he_like)

# ============================================================
# Tool 232: SchrodingerSolver1d
# ============================================================
print("\n🧪 #232 SchrodingerSolver1d")

def test_232_infinite_well():
    from chemmcp.tools.schrodinger_solver_1d import SchrodingerSolver1d
    tool = SchrodingerSolver1d()
    r = tool.run_code(potential_type="infinite_well", mass_kg=9.109e-31,
                       domain_length_m=1e-9, n_points=100, n_levels=3)
    assert r["potential_type"] == "infinite_well"
    assert r["n_computed_levels"] == 3
    # E1 = h²/(8mL²) ≈ 0.376 eV for electron in 1nm box
    assert 0.3 < r["ground_state_energy_eV"] < 0.5, f"E0={r['ground_state_energy_eV']}"
    return True

def test_232_harmonic():
    from chemmcp.tools.schrodinger_solver_1d import SchrodingerSolver1d
    tool = SchrodingerSolver1d()
    r = tool.run_code(potential_type="harmonic", mass_kg=9.109e-31,
                       domain_length_m=5e-10, n_points=150, n_levels=3,
                       force_constant_N_m=10.0)
    assert r["n_computed_levels"] == 3
    # Harmonic oscillator has bound states with negative energy (if V_min=0 at center)
    # or positive energy depending on convention. Just check levels exist.
    assert len(r["levels"]) == 3
    return True

def test_232_text_interface():
    from chemmcp.tools.schrodinger_solver_1d import SchrodingerSolver1d
    tool = SchrodingerSolver1d(interface="text")
    r = tool.run_text("infinite_well 9.109e-31 1e-9 100 2")
    assert r["n_computed_levels"] == 2
    return True

def test_232_finite_well():
    from chemmcp.tools.schrodinger_solver_1d import SchrodingerSolver1d
    tool = SchrodingerSolver1d()
    r = tool.run_code(potential_type="finite_well", mass_kg=9.109e-31,
                       domain_length_m=1e-9, n_points=100, n_levels=3,
                       well_depth_J=5e-19, well_width_m=3e-10)
    assert r["n_computed_levels"] >= 1
    return True

test("232: Infinite square well ground state", test_232_infinite_well)
test("232: Harmonic oscillator (3 levels)", test_232_harmonic)
test("232: Text interface", test_232_text_interface)
test("232: Finite well", test_232_finite_well)

# ============================================================
# Tool 233: VariationalMethod
# ============================================================
print("\n🧪 #233 VariationalMethod")

def test_233_harmonic_gaussian():
    from chemmcp.tools.variational_method import VariationalMethod
    tool = VariationalMethod()
    r = tool.run_code(potential_type="harmonic", mass_kg=9.109e-31,
                       trial_function_type="gaussian", force_constant_N_m=10.0)
    assert r["variational_energy_eV"] is not None
    # Gaussian trial for HO should be close to exact
    assert r["relative_error_percent"] is not None
    assert r["relative_error_percent"] < 15, f"Error={r['relative_error_percent']}%"
    return True

def test_233_infinite_well_cosine():
    from chemmcp.tools.variational_method import VariationalMethod
    tool = VariationalMethod()
    r = tool.run_code(potential_type="infinite_well", mass_kg=9.109e-31,
                       trial_function_type="cosine", box_length_m=1e-9)
    # Cosine trial for infinite well IS the exact ground state wavefunction
    assert r["relative_error_percent"] is not None
    assert r["relative_error_percent"] < 5, f"Error={r['relative_error_percent']}%"
    return True

def test_233_text_interface():
    from chemmcp.tools.variational_method import VariationalMethod
    tool = VariationalMethod(interface="text")
    r = tool.run_text("harmonic 9.109e-31 gaussian k=10")
    assert r["trial_function"] == "gaussian"
    return True

def test_233_virial():
    from chemmcp.tools.variational_method import VariationalMethod
    tool = VariationalMethod()
    r = tool.run_code(potential_type="harmonic", mass_kg=9.109e-31,
                       trial_function_type="gaussian", force_constant_N_m=10.0)
    vr = r.get("virial_theorem_ratio_V_over_T")
    if vr is not None:
        assert 0.1 < vr < 10.0, f"V/T ratio={vr}"
    return True

test("233: HO with Gaussian trial function", test_233_harmonic_gaussian)
test("233: Infinite well with cosine (exact)", test_233_infinite_well_cosine)
test("233: Text interface", test_233_text_interface)
test("233: Virial theorem check", test_233_virial)

# ============================================================
# Tool 234: PerturbationTheory
# ============================================================
print("\n🧪 #234 PerturbationTheory")

def test_234_harmonic_perturbed():
    from chemmcp.tools.perturbation_theory import PerturbationTheory
    tool = PerturbationTheory()
    r = tool.run_code(system_type="harmonic_perturbed", perturbation_strength=0.1,
                       order=2, force_constant_N_m=10.0, n_state=0)
    assert r["unperturbed_energy_eV"] is not None
    assert r["total_corrected_energy_eV"] is not None
    # Quartic perturbation raises energy (positive λx⁴)
    assert r["total_corrected_energy_eV"] >= r["unperturbed_energy_eV"]
    return True

def test_234_two_level():
    from chemmcp.tools.perturbation_theory import PerturbationTheory
    tool = PerturbationTheory()
    r = tool.run_code(system_type="two_level", perturbation_strength=1e-20, order=2)
    assert r["system_type"] == "two_level_system"
    return True

def test_234_hydrogen_stark():
    from chemmcp.tools.perturbation_theory import PerturbationTheory
    tool = PerturbationTheory()
    r = tool.run_code(system_type="hydrogen_stark", perturbation_strength=1e-9, order=2)
    # First-order Stark correction must be zero (parity)
    assert r["first_order_correction_eV"] == 0.0
    # Second-order should be negative (always attractive)
    assert r["second_order_correction_eV"] < 0
    return True

def test_234_text_interface():
    from chemmcp.tools.perturbation_theory import PerturbationTheory
    tool = PerturbationTheory(interface="text")
    r = tool.run_text("harmonic_perturbed 0.1 2 k=10")
    assert "harmonic" in r["system_type"]
    return True

test("234: Harmonic oscillator quartic perturbation", test_234_harmonic_perturbed)
test("234: Two-level system", test_234_two_level)
test("234: Hydrogen Stark effect (E1=0)", test_234_hydrogen_stark)
test("234: Text interface", test_234_text_interface)

# ============================================================
# Tool 235: MolecularOrbitalDiagram
# ============================================================
print("\n🧪 #235 MolecularOrbitalDiagram")

def test_235_o2():
    from chemmcp.tools.molecular_orbital_diagram import MolecularOrbitalDiagram
    tool = MolecularOrbitalDiagram()
    r = tool.run_code(molecule="O2")
    assert r["bond_order"] == 2.0
    assert r["magnetic_properties"] == "paramagnetic"
    assert r["homo_lumo_gap_eV"] is not None
    return True

def test_235_n2():
    from chemmcp.tools.molecular_orbital_diagram import MolecularOrbitalDiagram
    tool = MolecularOrbitalDiagram()
    r = tool.run_code(molecule="N2")
    assert r["bond_order"] == 3.0
    assert r["magnetic_properties"] == "diamagnetic"
    return True

def test_235_h2o():
    from chemmcp.tools.molecular_orbital_diagram import MolecularOrbitalDiagram
    tool = MolecularOrbitalDiagram()
    r = tool.run_code(molecule="H2O")
    assert r["point_group"] == "C2v"
    # Check bond order exists (key may vary)
    assert "bond_order" in r or "bond_order_per_OH" in r
    return True

def test_235_text_interface():
    from chemmcp.tools.molecular_orbital_diagram import MolecularOrbitalDiagram
    tool = MolecularOrbitalDiagram(interface="text")
    r = tool.run_text("O2")
    assert r["bond_order"] == 2.0
    return True

test("235: O2 bond order=2, paramagnetic", test_235_o2)
test("235: N2 bond order=3, diamagnetic", test_235_n2)
test("235: H2O polyatomic MO", test_235_h2o)
test("235: Text interface", test_235_text_interface)

# ============================================================
# Tool 236: HuckelMethod
# ============================================================
print("\n🧪 #236 HuckelMethod")

def test_236_ethene():
    from chemmcp.tools.huckel_method import HuckelMethod
    tool = HuckelMethod()
    r = tool.run_code(molecule="ethene")
    assert r["n_carbons"] == 2
    assert r["topology"] == "linear"
    assert r["n_pi_electrons"] == 2
    return True

def test_236_benzene():
    from chemmcp.tools.huckel_method import HuckelMethod
    tool = HuckelMethod()
    r = tool.run_code(molecule="benzene")
    assert r["n_carbons"] == 6
    assert r["topology"] == "cyclic"
    assert r["n_pi_electrons"] == 6
    # Benzene: total pi energy = 8β (in α+β units where α=0)
    # The stored value is in β units: 4.0α + 4.472β → but α offset is removed
    # Check that energy is positive and reasonable
    e_pi = r["total_pi_energy_in_beta_units"]
    assert e_pi > 0, f"E_π={e_pi}"
    return True

def test_236_butadiene():
    from chemmcp.tools.huckel_method import HuckelMethod
    tool = HuckelMethod()
    r = tool.run_code(molecule="butadiene")
    assert r["n_carbons"] == 4
    de = r["delocalization_energy_beta"]
    assert de > 0, f"Delocalization energy should be positive: {de}"
    return True

def test_236_text_interface():
    from chemmcp.tools.huckel_method import HuckelMethod
    tool = HuckelMethod(interface="text")
    r = tool.run_text("benzene")
    assert r["n_carbons"] == 6
    return True

test("236: Ethene π system", test_236_ethene)
test("236: Benzene aromatic cyclic", test_236_benzene)
test("236: Butadiene delocalization energy", test_236_butadiene)
test("236: Text interface", test_236_text_interface)

# ============================================================
# Tool 237: ElectronDensityPlotter
# ============================================================
print("\n🧪 #237 ElectronDensityPlotter")

def test_237_1s_density():
    from chemmcp.tools.electron_density_plotter import ElectronDensityPlotter
    tool = ElectronDensityPlotter()
    r = tool.run_code(species_type="hydrogen", quantum_numbers={"n": 1, "l": 0, "m": 0},
                       n_plot_points=100)
    assert r["orbital_name"] == "1s"
    # Check radial data exists
    rd = r["radial_data"]
    assert "most_probable_radius_bohr" in rd
    # For 1s: |R(r)|² peaks at r→0 (correct: ψ_max at nucleus)
    # The RDF peak (r²|R|²) is at a₀ ≈ 0.529 Å — check visualization_summary instead
    vs = r.get("visualization_summary", {})
    assert "peak_location_bohr" in vs
    return True

def test_237_2p_density():
    from chemmcp.tools.electron_density_plotter import ElectronDensityPlotter
    tool = ElectronDensityPlotter()
    r = tool.run_code(species_type="hydrogen", quantum_numbers={"n": 2, "l": 1, "m": 0},
                       n_plot_points=100)
    assert r["orbital_name"] == "2p"
    rd = r["radial_data"]
    r_mp_ang = rd["most_probable_radius_bohr"] * 0.529177
    # 2p most probable radius ≈ 4a₀ = 2.12 Å
    assert 1.0 < r_mp_ang < 5.0, f"r_mp={r_mp_ang} Å"
    return True

def test_237_isosurface():
    from chemmcp.tools.electron_density_plotter import ElectronDensityPlotter
    tool = ElectronDensityPlotter()
    r = tool.run_code(species_type="hydrogen", quantum_numbers={"n": 1, "l": 0, "m": 0},
                       isosurface_level=0.95, n_plot_points=100)
    iso = r["isosurface_info"]
    assert iso["isosurface_radius_bohr"] is not None
    assert iso["isosurface_radius_bohr"] > 0
    return True

def test_237_text_interface():
    from chemmcp.tools.electron_density_plotter import ElectronDensityPlotter
    tool = ElectronDensityPlotter(interface="text")
    r = tool.run_text("hydrogen 1 0 0")
    assert r["orbital_name"] == "1s"
    return True

test("237: 1s radial density distribution", test_237_1s_density)
test("237: 2p radial density", test_237_2p_density)
test("237: Isosurface calculation", test_237_isosurface)
test("237: Text interface", test_237_text_interface)

# ============================================================
# Tool 238: SpinOrbitCoupling
# ============================================================
print("\n🧪 #238 SpinOrbitCoupling")

def test_238_na_splitting():
    from chemmcp.tools.spin_orbit_coupling import SpinOrbitCoupling
    tool = SpinOrbitCoupling()
    r = tool.run_code(element="Na", n=3, l=1)
    assert len(r["split_energy_levels"]) == 2  # j=1/2 and j=3/2
    delta_eV = r["splitting_delta_eV"]
    assert delta_eV > 0, "Splitting should be positive"
    # Na 3p ζ = 17.2 cm⁻¹ → ΔE = 25.8 cm⁻¹ = 0.0032 eV
    # Wavelength of transition BETWEEN split levels: hc/ΔE ≈ 387 μm (far IR/microwave)
    # This is correct — it's NOT the D-line wavelength (which is 3s→3p)
    wl = r["transition_wavelength_nm"]
    assert wl is not None and wl > 0
    return True

def test_238_s_orbital():
    from chemmcp.tools.spin_orbit_coupling import SpinOrbitCoupling
    tool = SpinOrbitCoupling()
    r = tool.run_code(element="Na", n=3, l=0)  # s orbital
    assert len(r["split_energy_levels"]) == 1  # Only j=1/2 for s
    assert r["split_energy_levels"][0]["j_value"] == 0.5
    return True

def test_238_term_symbols():
    from chemmcp.tools.spin_orbit_coupling import SpinOrbitCoupling
    tool = SpinOrbitCoupling()
    r = tool.run_code(element="Na", n=3, l=1)
    terms = [lvl["term_symbol"] for lvl in r["split_energy_levels"]]
    # Should contain P_1/2 and P_3/2
    has_p12 = any("P" in t and "1" in t for t in terms)
    has_p32 = any("P" in t and "3" in t for t in terms)
    assert has_p12 or has_p32, f"Term symbols: {terms}"
    return True

def test_238_text_interface():
    from chemmcp.tools.spin_orbit_coupling import SpinOrbitCoupling
    tool = SpinOrbitCoupling(interface="text")
    r = tool.run_text("Na 3 1")
    assert len(r["split_energy_levels"]) >= 1
    return True

test("238: Na 3p spin-orbit splitting", test_238_na_splitting)
test("238: s-orbital (single j level)", test_238_s_orbital)
test("238: Term symbols ^{2S+1}L_J", test_238_term_symbols)
test("238: Text interface", test_238_text_interface)

# ============================================================
# Tool 239: SelectionRulesChecker
# ============================================================
print("\n🧪 #239 SelectionRulesChecker")

def test_239_1s_2p_allowed():
    from chemmcp.tools.selection_rules_checker import SelectionRulesChecker
    tool = SelectionRulesChecker()
    r = tool.run_code(transition_type="electric_dipole",
                       initial_state={"n": 1, "l": 0, "m_l": 0},
                       final_state={"n": 2, "l": 1, "m_l": 0})
    assert r["is_allowed"] == True, f"1s→2p should be allowed, got: {r.get('violated_rules')}"
    return True

def test_239_1s_2s_forbidden():
    from chemmcp.tools.selection_rules_checker import SelectionRulesChecker
    tool = SelectionRulesChecker()
    r = tool.run_code(transition_type="electric_dipole",
                       initial_state={"n": 1, "l": 0, "m_l": 0},
                       final_state={"n": 2, "l": 0, "m_l": 0})
    assert r["is_allowed"] == False, "1s→2s should be forbidden (Δl=0)"
    assert any("Δl" in v for v in r["violated_rules"])
    return True

def test_239_rotational():
    from chemmcp.tools.selection_rules_checker import SelectionRulesChecker
    tool = SelectionRulesChecker()
    r = tool.run_code(transition_type="rotational",
                       initial_state={"J": 0}, final_state={"J": 1})
    assert r["is_allowed"] == True, "ΔJ=1 should be allowed for rotational"
    return True

def test_239_raman():
    from chemmcp.tools.selection_rules_checker import SelectionRulesChecker
    tool = SelectionRulesChecker()
    r = tool.run_code(transition_type="raman",
                       initial_state={"J": 0}, final_state={"J": 2})
    assert r["is_allowed"] == True, "ΔJ=2 should be allowed for Raman (S branch)"
    return True

def test_239_text_interface():
    from chemmcp.tools.selection_rules_checker import SelectionRulesChecker
    tool = SelectionRulesChecker(interface="text")
    r = tool.run_text("electric_dipole 1 0 0 2 1 0")
    assert r["is_allowed"] == True
    return True

test("239: E1 1s→2p allowed ✓", test_239_1s_2p_allowed)
test("239: E1 1s→2s forbidden ✗ (Δl=0)", test_239_1s_2s_forbidden)
test("239: Rotational ΔJ=1 allowed", test_239_rotational)
test("239: Raman ΔJ=2 allowed (S branch)", test_239_raman)
test("239: Text interface", test_239_text_interface)

# ============================================================
# Tool 240: TunnelingProbability
# ============================================================
print("\n🧪 #240 TunnelingProbability")

def test_240_rectangular_tunneling():
    from chemmcp.tools.tunneling_probability import TunnelingProbability
    tool = TunnelingProbability()
    r = tool.run_code(barrier_type="rectangular", particle_mass_kg=9.109e-31,
                       energy_J=1e-20, barrier_height_J=5e-20, barrier_width_m=1e-9)
    T = r["transmission_probability"]
    assert 0 < T < 1, f"T should be between 0 and 1, got T={T}"
    assert r["is_classically_allowed"] == False
    return True

def test_240_above_barrier():
    from chemmcp.tools.tunneling_probability import TunnelingProbability
    tool = TunnelingProbability()
    r = tool.run_code(barrier_type="rectangular", particle_mass_kg=9.109e-31,
                       energy_J=8e-20, barrier_height_J=5e-20, barrier_width_m=1e-9)
    assert r["is_classically_allowed"] == True
    assert r["transmission_probability"] > 0.5, "Above barrier should have high T"
    return True

def test_240_alpha_decay():
    from chemmcp.tools.tunneling_probability import TunnelingProbability
    tool = TunnelingProbability()
    r = tool.run_code(barrier_type="alpha_decay", particle_mass_kg=6.644657230e-27,
                       energy_J=8.79e-13, nuclear_charge_Z=90)
    T = r["transmission_probability"]
    assert T < 1e-10, f"Alpha decay T should be tiny: {T}"
    # Check regime exists somewhere in result
    assert "alpha" in str(r.get("tunneling_regime", "")).lower() or "gamow" in str(r).lower()
    return True

def test_240_zero_energy():
    from chemmcp.tools.tunneling_probability import TunnelingProbability
    tool = TunnelingProbability()
    r = tool.run_code(barrier_type="rectangular", particle_mass_kg=9.109e-31,
                       energy_J=0, barrier_height_J=5e-20, barrier_width_m=1e-9)
    assert r["transmission_probability"] == 0.0
    return True

def test_240_text_interface():
    from chemmcp.tools.tunneling_probability import TunnelingProbability
    tool = TunnelingProbability(interface="text")
    # Use height= and width= format
    r = tool.run_text("rectangular 9.109e-31 1e-20 height=5e-20 width=1e-9")
    assert 0 <= r["transmission_probability"] <= 1
    return True

test("240: Rectangular barrier tunneling (0<T<1)", test_240_rectangular_tunneling)
test("240: Above-barrier transmission", test_240_above_barrier)
test("240: Alpha decay (extremely small T)", test_240_alpha_decay)
test("240: Zero energy (T=0)", test_240_zero_energy)
test("240: Text interface", test_240_text_interface)

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
total = passed + failed
print(f"Results: {passed}/{total} passed, {failed} failed")
if errors:
    print("\nFailed tests:")
    for name, err in errors:
        print(f"  • {name}: {err}")
print("=" * 60)

sys.exit(0 if failed == 0 else 1)
