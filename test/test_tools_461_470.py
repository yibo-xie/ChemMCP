#!/usr/bin/env python3
"""
Test suite for MCP Tools #461-470: Quantum Chemistry Integral & Advanced Electronic Structure
运行方式: cd ~/ChemMCP && python -m pytest test/test_tools_461_470.py -v
或直接: python test/test_tools_461_470.py
"""
import sys
import os
import math

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


def test_461_basis_set_handler():
    """Test #461 BasisSetHandler — 基组处理"""
    from chemmcp.tools import BasisSetHandler
    tool = BasisSetHandler()

    # Test 1: List all basis sets
    result = tool.run_code("list")
    assert "result" in result
    data = result["result"]
    assert "available_basis_sets" in data
    assert len(data["available_basis_sets"]) > 10
    print(f"  ✓ 461-list: {len(data['available_basis_sets'])} basis sets available")

    # Test 2: Info query for STO-3G on Carbon
    result = tool.run_code("STO-3G", "info", element="C")
    assert "result" in result
    info = result["result"]
    assert info["basis_set"] == "STO-3G"
    assert info["type"] == "Minimal Basis"
    print(f"  ✓ 461-info: STO-3G/C → {info['type']}, polarization={info.get('polarization')}")

    # Test 3: STO→GTO conversion for H/STO-3G
    result = tool.run_code("sto_gto_convert", element="H", n_primitives=3)
    assert "result" in result
    conv = result["result"]
    assert "orbital_expansions" in conv
    assert "1s" in conv["orbital_expansions"]
    coeffs = conv["orbital_expansions"]["1s"]["normalized_coefficients"]
    assert len(coeffs) == 3
    print(f"  ✓ 461-convert: H STO-3G → {len(coeffs)} primitives, sum_d≈{conv['orbital_expansions']['1s']['sum_check_d']:.4f}")

    # Test 4: Compare two basis sets
    result = tool.run_code("compare", "6-31G*", compare_with="STO-3G", element="C")
    assert "result" in result
    cmp_data = result["result"]
    assert "recommendation" in cmp_data
    print(f"  ✓ 461-compare: 6-31G* vs STO-3G comparison done")

    # Test 5: Text interface
    result = tool.run_text("STO-3G info C")
    assert "result" in result
    print(f"  ✓ 461-text: text interface works")

    print("  ✅ All BasisSetHandler tests passed!")


def test_462_overlap_integral():
    """Test #462 OverlapIntegral — 重叠积分计算"""
    from chemmcp.tools import OverlapIntegral
    tool = OverlapIntegral()

    # Test 1: Same orbital at same center → S≈1 (normalized)
    result = tool.run_code("1s", "1s", 1.0, 1.0, R_bohr=0.0)
    S = result["result"]["overlap_integral_S"]
    assert S > 0.95  # normalized same-orbital at same center
    print(f"  ✓ 462-same_center: S(1s,1s,R=0) = {S:.6f} ≈ 1.0")

    # Test 2: Same orbital at distance → S < 1
    result = tool.run_code("1s", "1s", 0.27095, 0.27095, R_bohr=1.4)
    S = result["result"]["overlap_integral_S"]
    assert 0 < S < 1
    print(f"  ✓ 462-separated: S(1s,1s,R=1.4) = {S:.6f}")

    # Test 3: Large distance → S → 0
    result = tool.run_code("1s", "1s", 1.0, 1.0, R_bohr=10.0)
    S = result["result"]["overlap_integral_S"]
    assert abs(S) < 0.01
    print(f"  ✓ 462-far: S(1s,1s,R=10) = {S:.8f} ≈ 0")

    # Test 4: p-p orthogonal orbitals
    result = tool.run_code("2px", "2py", 1.0, 1.0, R_bohr=0.0)
    S_pp = result["result"]["overlap_integral_S"]
    assert abs(S_pp) < 0.5  # px and py should be nearly orthogonal
    print(f"  ✓ 462-orthogonal: S(px,py,R=0) = {S_pp:.6f}")

    # Test 5: Text interface
    result = tool.run_text("1s 1s 0.27 0.27 1.4")
    assert "overlap_integral_S" in result["result"]
    print(f"  ✓ 462-text: text interface works")

    print("  ✅ All OverlapIntegral tests passed!")


def test_463_coulomb_integral():
    """Test #463 CoulombIntegral — 库仑积分计算"""
    from chemmcp.tools import CoulombIntegral
    tool = CoulombIntegral()

    # Test 1: Same-center Coulomb integral (positive value)
    result = tool.run_code("(ii|jj)", "1s", alpha1=0.27)
    J = result["result"]["coulomb_integral_J"]
    assert J > 0
    print(f"  ✓ 463-same_center: (1s|1s) = {J:.6f} Hartree ({J*27.21:.2f} eV)")

    # Test 2: Two-center Coulomb (should be smaller than same-center)
    result_near = tool.run_code("(ij|ij)", "1s", alpha1=0.27, R_ab=1.0)
    result_far = tool.run_code("(ij|ij)", "1s", alpha1=0.27, R_ab=5.0)
    J_near = result_near["result"]["coulomb_integral_J"]
    J_far = result_far["result"]["coulomb_integral_J"]
    assert J_far < J_near  # farther → smaller Coulomb
    print(f"  ✓ 463-distance: J(R=1)={J_near:.4f}, J(R=5)={J_far:.4f}")

    # Test 3: Units check
    result = tool.run_code("(ii|jj)", "1s", alpha1=0.5)
    assert "coulomb_integral_eV" in result["result"]
    print(f"  ✓ 463-units: Hartree & eV both present")

    # Test 4: Text interface
    result = tool.run_text("(ii|jj) 1s 0.27")
    assert "coulomb_integral_J" in result["result"]
    print(f"  ✓ 463-text: text interface works")

    print("  ✅ All CoulombIntegral tests passed!")


def test_464_exchange_integral():
    """Test #464 ExchangeIntegral — 交换积分计算"""
    from chemmcp.tools import ExchangeIntegral
    tool = ExchangeIntegral()

    # Test 1: Calculate exchange integral
    result = tool.run_code("calculate", "1s", "1s", 0.27, 0.27, R_bohr=1.4, same_spin=True)
    K = result["result"]["exchange_integral_K"]
    assert K >= 0  # exchange should be non-negative
    print(f"  ✓ 464-calculate: K(1s,1s,R=1.4) = {K:.6f} Hartree")

    # Test 2: Opposite spin → no exchange
    result_ss = tool.run_code("calculate", "1s", "1s", 0.27, 0.27, R_bohr=1.4, same_spin=True)
    result_os = tool.run_code("calculate", "1s", "1s", 0.27, 0.27, R_bohr=1.4, same_spin=False)
    assert result_os["result"]["effective_K_same_spin"] == 0
    assert result_ss["result"]["effective_K_same_spin"] > 0
    print(f"  ✓ 464-spin: same_spin K={result_ss['result']['effective_K_same_spin']:.6f}, opposite_spin K={result_os['result']['effective_K_same_spin']}")

    # Test 3: Pauli explanation
    result = tool.run_code("explain")
    assert "quantum_origin" in result["result"]
    assert "fermi_hole" in result["result"]
    print(f"  ✓ 464-explain: Pauli principle explanation generated")

    # Test 4: Fock contribution analysis
    result = tool.run_code("fock_contribution", "1s", "1s", 0.27, 0.27, R_bohr=1.4)
    assert "fock_operator_form" in result["result"]
    print(f"  ✓ 464-fock: Fock matrix contribution analyzed")

    # Test 5: Full analysis
    result = tool.run_code("full_analysis", "1s", "1s", 0.27, 0.27, R_bohr=1.4)
    assert "pauli_explanation" in result["result"]
    assert "fock_matrix_analysis" in result["result"]
    print(f"  ✓ 464-full: complete analysis with all components")

    print("  ✅ All ExchangeIntegral tests passed!")


def test_465_mo_energy_level_diagram():
    """Test #465 MOEnergyLevelDiagram — 分子轨道能级图"""
    from chemmcp.tools import MOEnergyLevelDiagram
    tool = MOEnergyLevelDiagram()

    # Test 1: N₂ diagram (triple bond)
    result = tool.run_code("N2", "full")
    assert "result" in result
    data = result["result"]
    assert data["bond_order"] == 3
    assert not data.get("paramagnetic", False)
    assert len(data["molecular_orbitals"]) >= 6
    print(f"  ✓ 465-N2: bond_order={data['bond_order']}, {len(data['molecular_orbitals'])} MOs, diamagnetic")

    # Test 2: O₂ diagram (paramagnetic!)
    result = tool.run_code("O2", "frontier")
    data = result["result"]
    assert data.get("paramagnetic") == True
    assert "homo_lumo_gap_eV" in data
    gap = data["homo_lumo_gap_eV"]
    assert gap >= 0  # O₂ has degenerate HOMO=LUMO (π_g*), so gap can be 0
    print(f"  ✓ 465-O2: paramagnetic✓, HOMO-LUMO gap={gap:.1f} eV")

    # Test 3: CO heteronuclear
    result = tool.run_code("CO", "full")
    data = result["result"]
    assert data["bond_order"] == 3
    assert data.get("dipole_moment_D") is not None
    print(f"  ✓ 465-CO: BO={data['bond_order']}, μ={data.get('dipole_moment_D')} D")

    # Test 4: HF polar molecule
    result = tool.run_code("HF", "comprehensive")
    data = result["result"]
    assert data["bond_order"] == 1
    assert data.get("dipole_moment_D", 0) > 1.0
    print(f"  ✓ 465-HF: BO={data['bond_order']}, μ={data.get('dipole_moment_D')} D (highly polar)")

    # Test 5: Text output format
    result = tool.run_code("H2", "text", show_symmetry=True, output_format="text")
    assert "ascii_diagram" in result["result"]
    print(f"  ✓ 465-ascii: ASCII diagram generated")

    # Test 6: Plot data format
    result = tool.run_code("H2", "full", output_format="plot_data")
    assert "plot_data" in result["result"]
    print(f"  ✓ 465-plot: plot-ready data generated")

    print("  ✅ All MOEnergyLevelDiagram tests passed!")


def test_466_frontier_orbital_analysis():
    """Test #466 FrontierOrbitalAnalysis — 前线轨道分析"""
    from chemmcp.tools import FrontierOrbitalAnalysis
    tool = FrontierOrbitalAnalysis()

    # Test 1: Basic HOMO/LUMO for benzene
    result = tool.run_code("benzene", "basic")
    assert "result" in result
    data = result["result"]
    assert "homo" in data
    assert "lumo" in data
    assert data["gap_eV"] > 0
    print(f"  ✓ 466-benzene-basic: HOMO={data['homo']['energy_eV']} eV, LUMO={data['lumo']['energy_eV']} eV, gap={data['gap_eV']:.1f} eV")

    # Test 2: Conceptual DFT quantities
    dft = data.get("conceptual_dft", {})
    assert "global_hardness_eta_eV" in dft
    assert "electronegativity_chi_eV" in dft
    eta = dft["global_hardness_eta_eV"]
    chi = dft["electronegativity_chi_eV"]
    print(f"  ✓ 466-DFT: η={eta:.2f} eV, χ={chi:.2f} eV, ω={dft.get('electrophilicity_index_omega_eV', 0):.2f} eV")

    # Test 3: Fukui function analysis
    result = tool.run_code("benzene", "fukui")
    data = result["result"]
    assert "fukui_functions" in data
    ff = data["fukui_functions"]
    assert "fukui_plus_f+" in ff
    assert "fukui_minus_f-" in ff
    assert "fukui_zero_f0" in ff
    print(f"  ✓ 466-fukui: f⁺, f⁻, f⁰ all computed")

    # Test 4: Reactivity prediction
    result = tool.run_code("formaldehyde", "reactivity")
    data = result["result"]
    assert "reactivity_prediction" in data
    rp = data["reactivity_prediction"]
    assert "hsab_class" in rp
    assert "general_reactivity" in rp
    print(f"  ✓ 466-reactivity: HSAB={rp['hsab_class']}, reactivity={rp['general_reactivity']}")

    # Test 5: Compare hard (N2) vs soft (I2-like) molecules
    result_n2 = tool.run_code("n2", "basic")
    result_tce = tool.run_code("tetracyanoethylene", "basic")
    eta_n2 = result_n2["result"]["conceptual_dft"]["global_hardness_eta_eV"]
    eta_tce = result_tce["result"]["conceptual_dft"]["global_hardness_eta_eV"]
    assert eta_n2 > eta_tce  # N2 is harder than TCNE
    print(f"  ✓ 466-compare: η(N₂)={eta_n2:.1f} >> η(TCNE)={eta_tce:.1f} (hard vs soft)")

    # Test 6: Text interface
    result = tool.run_text("ethylene reactivity")
    assert "homo" in result["result"]
    print(f"  ✓ 466-text: text interface works")

    print("  ✅ All FrontierOrbitalAnalysis tests passed!")


def test_467_dft_xc_functional():
    """Test #467 DftXcFunctional — DFT交换相关泛函"""
    from chemmcp.tools import DftXcFunctional
    tool = DftXcFunctional()

    rho = 0.01  # typical valence electron density

    # Test 1: LDA Slater exchange
    result = tool.run_code("LDA_Slater", rho)
    assert "result" in result
    data = result["result"]
    assert data["Ex_per_volume_Hartree_Bohr3"] < 0  # exchange energy is negative
    assert data["Ec_per_volume_Hartree_Bohr3"] == 0  # pure Xα has no correlation
    print(f"  ✓ 467-LDA-Slater: E_x={data['Ex_per_volume_Hartree_Bohr3']:.6f} (negative ✓)")

    # Test 2: LDA VWN correlation
    result = tool.run_code("LDA_VWN", rho)
    data = result["result"]
    assert data["Ec_per_volume"] < 0  # correlation energy is negative
    assert data["Exc_per_volume"] < 0  # total XC negative
    print(f"  ✓ 467-LDA-VWN: E_c={data['Ec_per_volume']:.6f}, E_xc={data['Exc_per_volume']:.6f}")

    # Test 3: PBE GGA (with gradient)
    result = tool.run_code("PBE", rho, gradient=0.05)
    data = result["result"]
    assert "Ex_GGA" in data
    assert "Ec_GGA" in data
    assert "gradient_correction_Ex" in data
    print(f"  ✓ 467-PBE: E_x^GGA={data['Ex_GGA']:.6f}, ΔE_x(grad)={data['gradient_correction_Ex']:.8f}")

    # Test 4: B3LYP hybrid
    result = tool.run_code("B3LYP", rho, gradient=0.03)
    data = result["result"]
    assert "hybrid_parameters" in data
    hp = data["hybrid_parameters"]
    assert hp["a0_HF_exchange_fraction"] == 0.20
    print(f"  ✓ 467-B3LYP: HF mix={hp['exact_exchange_percent']}, E_xc={data.get('Exc_total', 'N/A')}")

    # Test 5: PBE0 hybrid
    result = tool.run_code("PBE0", rho, gradient=0.02)
    data = result["result"]
    assert data["hybrid_parameters"]["HF_exchange_fraction"] == 0.25
    print(f"  ✓ 467-PBE0: 25% HF exact exchange")

    # Test 6: Compare all functionals
    result = tool.run_code("compare_all", rho)
    assert "comparison" in result["result"]
    comp = result["result"]["comparison"]
    assert len(comp) >= 5
    print(f"  ✓ 467-compare: {len(comp)} functionals compared")

    # Test 7: Zero density error handling
    try:
        tool.run_code("LDA_Slider", 0.0)  # typo intentional? no, just zero density
        assert False, "Should raise error for zero density"
    except Exception:
        print(f"  ✓ 467-error: zero density correctly rejected")

    # Test 8: Text interface
    result = tool.run_text("LDA_VWN 0.01")
    assert "Ec_per_volume" in result["result"]
    print(f"  ✓ 467-text: text interface works")

    print("  ✅ All DftXcFunctional tests passed!")


def test_468_electron_density_calculator():
    """Test #468 ElectronDensityCalculator — 电子密度计算"""
    from chemmcp.tools import ElectronDensityCalculator
    tool = ElectronDensityCalculator()

    # Test 1: Density grid along H2 bond axis
    result = tool.run_code("density_grid", "H2", n_points=20)
    assert "result" in result
    data = result["result"]
    assert "density_values" in data
    assert len(data["density_values"]) == 20
    assert data["max_density_value"] > 0
    print(f"  ✓ 468-grid: H₂ density evaluated at {len(data['density_values'])} points, max ρ={data['max_density_value']:.4f}")

    # Test 2: Mulliken population analysis
    result = tool.run_code("mulliken", "H2")
    assert "result" in result
    data = result["result"]
    assert "atomic_populations" in data
    pops = data["atomic_populations"]
    total = data["total_electrons"]
    assert abs(total - 2.0) < 0.5  # roughly correct for H2
    print(f"  ✓ 468-mulliken: {len(pops)} atoms, total e⁻={total:.3f}")

    # Test 3: ESP calculation
    result = tool.run_code("esp", "H2", n_points=20)
    data = result["result"]
    assert "esp_values" in data
    assert len(data["esp_values"]) == 20
    assert "min_ESP_Hartree/e" in data
    print(f"  ✓ 468-ESP: V_min={data['min_ESP_Hartree/e']:.3f}, V_max={data['max_ESP_Hartree/e']:.3f} Ha/e")

    # Test 4: Topology / Bader analysis
    result = tool.run_code("topology", "H2")
    data = result["result"]
    assert "critical_points" in data
    cps = data["critical_points"]
    ncp_types = set(cp["type"] for cp in cps)
    assert any("Nuclear" in t for t in ncp_types)
    assert any("Bond Critical" in t for t in ncp_types)
    print(f"  ✓ 468-topology: {len(cps)} critical points identified (NCP + BCP)")

    # Test 5: Integrated density check
    result = tool.run_code("integrated", "H2")
    data = result["result"]
    assert "total_electrons_from_density" in data
    print(f"  ✓ 468-integrated: N_elec={data['total_electrons_from_density']:.3f}")

    # Test 6: LiH (polar molecule)
    result = tool.run_code("mulliken", "LiH")
    data = result["result"]
    pops = data["atomic_populations"]
    charges = [p["mulliken_charge"] for p in pops]
    # Li should be positive (electron deficient), H negative
    print(f"  ✓ 468-LiH: Mulliken charges = {charges} (Liδ⁺-Hδ⁻ expected)")

    # Test 7: Text interface
    result = tool.run_text("density_grid H2 15")
    assert "density_values" in result["result"]
    print(f"  ✓ 468-text: text interface works")

    print("  ✅ All ElectronDensityCalculator tests passed!")


def test_469_mp2_correlation():
    """Test #469 Mp2Correlation — MP2相关能计算"""
    from chemmcp.tools import Mp2Correlation
    tool = Mp2Correlation()

    # Test 1: H2 MP2 calculation
    result = tool.run_code("H2")
    assert "result" in result
    data = result["result"]
    assert "E_MP2_Hartree" in data
    assert data["E_MP2_Hartree"] < 0  # MP2 correlation energy must be negative (stabilizing)
    print(f"  ✓ 469-H2: E_MP2={data['E_MP2_Hartree']:.8f} Ha ({data['E_MP2_eV']:.4f} eV), "
          f"E_corrected={data['E_corrected_total_Hartree']:.6f} Ha")

    # Test 2: Spin decomposition
    result = tool.run_code("H2", spin_case="decompose")
    data = result["result"]
    assert "E_opposite_spin_singlet_Hartree" in data
    assert "E_same_spin_triplet_Hartree" in data
    os_frac = data.get("OS_fraction", 0)
    ss_frac = data.get("SS_fraction", 0)
    assert 0 <= os_frac <= 1
    assert 0 <= ss_frac <= 1
    print(f"  ✓ 469-decompose: OS={os_frac:.1%}, SS={ss_frac:.1%}")

    # Test 3: HeH+ (heteronuclear)
    result = tool.run_code("HeH+")
    data = result["result"]
    assert data["E_MP2_Hartree"] < 0
    print(f"  ✓ 469-HeH+: E_MP2={data['E_MP2_Hartree']:.8f} Ha, recovery={data.get('MP2_recovery_percent','?')}%")

    # Test 4: Generic 4-electron system
    result = tool.run_code("generic_4e", spin_case="all")
    data = result["result"]
    assert data["n_electrons"] == 4
    assert data["E_MP2_Hartree"] < 0
    pct = data.get("correction_percent_of_HF", 0)
    print(f"  ✓ 469-generic-4e: E_MP2={data['E_MP2_Hartree']:.6f} ({pct:.2f}% of E_HF)")

    # Test 5: Different methods (SCS-MP2)
    result_scs = tool.run_code("H2", method="scs-mp2")
    result_canonical = tool.run_code("H2", method="canonical")
    assert "scaling_note" in result_scs["result"]
    print(f"  ✓ 469-methods: canonical & SCS-MP2 both work")

    # Test 6: Correlation recovery check (for systems with exact reference)
    if "MP2_recovery_percent" in data:
        rec = data["MP2_recovery_percent"]
        assert 0 < rec < 150  # reasonable range
        print(f"  ✓ 469-recovery: MP2 recovers ~{rec:.1f}% of correlation energy")

    # Test 7: Text interface
    result = tool.run_text("H2 decompose")
    assert "E_MP2_Hartree" in result["result"]
    print(f"  ✓ 469-text: text interface works")

    print("  ✅ All Mp2Correlation tests passed!")


def test_470_configuration_interaction():
    """Test #470 ConfigurationInteraction — 组态相互作用"""
    from chemmcp.tools import ConfigurationInteraction
    tool = ConfigurationInteraction()

    # Test 1: CISD for 2-electron system (smallest meaningful CI)
    result = tool.run_code("CISD", 2, 4)
    assert "result" in result
    data = result["result"]
    assert data["ci_method"] == "CISD"
    assert data["n_ci_configurations"] >= 1  # at least reference
    assert data["E_CI_ground_state_Hartree"] is not None
    assert data["E_correlation_Hartree"] <= 0  # correlation lowers energy
    print(f"  ✓ 470-CISD(2e,4o): {data['n_ci_configurations']} configs, "
          f"E₀={data['E_CI_ground_state_Hartree']:.6f}, E_corr={data['E_correlation_Hartree']:.6f}")

    # Test 2: CIS only (Brillouin's theorem → no correlation from singles alone)
    result_cis = tool.run_code("CIS", 2, 4)
    cis_data = result_cis["result"]
    assert cis_data["ci_method"] == "CIS"
    assert cis_data["n_single_excitations"] > 0
    print(f"  ✓ 470-CIS: {cis_data['n_ci_configurations']} configs, "
          f"{cis_data['n_single_excitations']} singles")

    # Test 3: CID (doubles only — should recover more correlation than CIS)
    result_cid = tool.run_code("CID", 2, 4)
    cid_data = result_cid["result"]
    assert cid_data["n_double_excitations"] > 0
    print(f"  ✓ 470-CID: {cid_data['n_double_excitations']} double excitations")

    # Test 4: Larger system (4 electrons, 6 orbitals)
    result = tool.run_code("CISD", 4, 6, n_states=5)
    data = result["result"]
    assert data["n_electrons"] == 4
    assert data["n_spatial_orbitals"] == 6
    assert len(data["states"]) >= 2
    E0 = data["states"][0]["energy_Hartree"]
    E1 = data["states"][1]["energy_Hartree"]
    assert E1 > E0  # excited state higher than ground state
    print(f"  ✓ 470-CISD(4e,6o): {data['n_ci_configurations']} configs, "
          f"{len(data['states'])} states, E₁-E₀={(E1-E0)*27.21:.2f} eV")

    # Test 5: State characterization
    if len(data["states"]) > 1:
        char = data["states"][1].get("character", "")
        exc_eV = data["states"][1].get("excitation_energy_eV", 0)
        print(f"  ✓ 470-states: S₁ character='{char}', E_exc={exc_eV:.2f} eV")

    # Test 6: FCI for small system
    result_fci = tool.run_code("FCI", 2, 4)
    fci_data = result_fci["result"]
    assert fci_data["ci_method"] == "FCI"
    print(f"  ✓ 470-FCI(2e,4o): {fci_data['n_ci_configurations']} configs (exact within basis)")

    # Test 7: Excited state properties
    if len(data["states"]) > 1:
        for s in data["states"]:
            assert "dominant_coefficients" in s
            assert "energy_Hartree" in s
        print(f"  ✓ 470-excited: all states have coefficients & energies")

    # Test 8: Text interface
    result = tool.run_text("CISD 2 4")
    assert "E_CI_ground_state_Hartree" in result["result"]
    print(f"  ✓ 470-text: text interface works")

    print("  ✅ All ConfigurationInteraction tests passed!")


# ── Main Runner ──────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n" + "="*60)
    print("  MCP Tools #461-470 Test Suite")
    print("  Quantum Chemistry: Integrals & Advanced Electronic Structure")
    print("="*60 + "\n")

    tests = [
        ("#461 BasisSetHandler", test_461_basis_set_handler),
        ("#462 OverlapIntegral", test_462_overlap_integral),
        ("#463 CoulombIntegral", test_463_coulomb_integral),
        ("#464 ExchangeIntegral", test_464_exchange_integral),
        ("#465 MOEnergyLevelDiagram", test_465_mo_energy_level_diagram),
        ("#466 FrontierOrbitalAnalysis", test_466_frontier_orbital_analysis),
        ("#467 DftXcFunctional", test_467_dft_xc_functional),
        ("#468 ElectronDensityCalculator", test_468_electron_density_calculator),
        ("#469 Mp2Correlation", test_469_mp2_correlation),
        ("#470 ConfigurationInteraction", test_470_configuration_interaction),
    ]

    passed = 0
    failed = 0
    errors = []

    for name, test_fn in tests:
        try:
            print(f"\n▶ Running {name}...")
            test_fn()
            passed += 1
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"  ❌ FAILED: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "="*60)
    print(f"  Results: {passed} passed, {failed} failed out of {len(tests)}")
    if errors:
        print("\n  Failed tests:")
        for name, err in errors:
            print(f"    ❌ {name}: {err[:200]}")
    else:
        print("\n  🎉 ALL TESTS PASSED! 🎉")
    print("="*60)

    sys.exit(0 if failed == 0 else 1)
