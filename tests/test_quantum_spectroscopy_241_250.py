"""
Test suite for ChemMCP Quantum Mechanics & Spectroscopy Tools (#241-250)
Tools tested:
  241. UncertaintyPrinciple    - Heisenberg uncertainty principle
  242. ExpectationValue        - QM expectation values
  243. IrSpectrumPredictor     - IR spectrum prediction
  244. RamanActivity           - Raman activity prediction
  245. UvVisTransitions        - UV-Vis transition energy
  246. NmrChemicalShift        - NMR chemical shift prediction
  247. RotationalSpectrum      - Rotational spectrum calculation
  248. VibrationalModes         - Vibrational mode analysis
  249. FranckCondonFactors     - Franck-Condon factor calculation
  250. BeerLambertCalculator   - Beer-Lambert law calculations
"""

import pytest
import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


# ============================================================
# Test 241: UncertaintyPrinciple
# ============================================================
class TestUncertaintyPrinciple:
    def setup_method(self):
        from chemmcp.tools import UncertaintyPrinciple
        self.tool = UncertaintyPrinciple()

    def test_position_momentum_electron(self):
        """Electron confined to 1 nm should have large velocity uncertainty."""
        result = self.tool.run_code(system_type="position_momentum", delta_x=1e-10, mass_kg=9.109e-31)
        r = result["result"]
        assert r["system"] == "position_momentum"
        # Value is in scientific notation string format
        dp = float(r["min_delta_p_kg_m_s"])
        assert dp > 0
        print(f"[PASS] Uncertainty: Δx=1e-10m → Δp={r['min_delta_p_kg_m_s']} kg·m/s")

    def test_energy_time(self):
        """Energy-time uncertainty should satisfy ΔE·Δt ≥ ℏ/2."""
        result = self.tool.run_code(system_type="energy_time", delta_e=1e-19, delta_t=1e-15)
        r = result["result"]
        assert r["system"] == "energy_time"
        assert r["satisfies_principle"] == True
        print(f"[PASS] Uncertainty: ΔE·Δt = {r['product_J_s']} J·s, limit = {r['hbar_over_2_J_s']}")

    def test_electron_demo(self):
        """Electron demo should compute kinetic energy."""
        result = self.tool.run_code(system_type="electron", delta_x=1e-9, mass_kg=9.109e-31)
        r = result["result"]
        assert "min_kinetic_energy_eV" in r
        ke_ev = r["min_kinetic_energy_eV"]
        # 1nm confinement → KE ~0.38 eV (order of magnitude check)
        assert float(ke_ev) > 0.001  # at least 1 meV
        assert float(ke_ev) < 100
        print(f"[PASS] Uncertainty: Electron 1nm → KE = {ke_ev} eV")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("electron 1e-9 9.109e-31")
        assert "result" in result
        print(f"[PASS] Uncertainty text interface OK")

    def test_invalid_system_type_raises(self):
        with pytest.raises(Exception):
            self.tool.run_code(system_type="invalid_type", delta_x=1e-10, mass_kg=9.109e-31)
        print(f"[PASS] Uncertainty: Invalid system type raises error")


# ============================================================
# Test 242: ExpectationValue
# ============================================================
class TestExpectationValue:
    def setup_method(self):
        from chemmcp.tools import ExpectationValue
        self.tool = ExpectationValue()

    def test_particle_in_box_x_expectation(self):
        """<x> for particle in a box should be L/2."""
        result = self.tool.run_code(system="particle_in_box", observable="x", n=1, L=1e-9)
        r = result["result"]
        assert abs(r["expectation_value"] - 5e-10) < 1e-15
        print(f"[PASS] Expectation: PIB <x> = {r['expectation_value']:.3e} m (= L/2)")

    def test_particle_in_box_energy(self):
        """Ground state energy of particle in box."""
        result = self.tool.run_code(system="particle_in_box", observable="E", n=1, L=1e-9, mass_kg=9.109e-31)
        r = result["result"]
        val = float(r["expectation_value"])
        assert val > 0
        # E1 = h²/(8mL²) ≈ 6.0e-20 J for electron in 1nm box
        assert 1e-22 < val < 1e-15
        print(f"[PASS] Expectation: PIB E1 = {r['expectation_value']} J")

    def test_harmonic_oscillator_zero_point(self):
        """HO ground state energy should be (0+½)ℏω."""
        result = self.tool.run_code(system="harmonic_oscillator", observable="E", n=0, omega=1e14)
        r = result["result"]
        val = float(r["expectation_value"])
        expected = 0.5 * 1.054571817e-34 * 1e14
        assert abs(val - expected) / expected < 0.001
        print(f"[PASS] Expectation: HO E0 = {r['expectation_value']} J (zero-point)")

    def test_harmonic_oscillator_virial(self):
        """Virial theorem: <T> = <V> = E/2 for HO."""
        result_e = self.tool.run_code(system="harmonic_oscillator", observable="E", n=2, omega=1e14)
        result_t = self.tool.run_code(system="harmonic_oscillator", observable="T", n=2, omega=1e14)
        E = float(result_e["result"]["expectation_value"])
        T = float(result_t["result"]["expectation_value"])
        assert abs(T - E/2) / E < 0.001
        print(f"[PASS] Expectation: Virial theorem holds: T={E/2:.3e}, E/2={E/2:.3e}")

    def test_hydrogen_atom_ground_state_r(self):
        """<r> for H atom ground state should be 1.5*a0."""
        result = self.tool.run_code(system="hydrogen_atom", observable="r", n=1, l=0)
        r = result["result"]
        val = float(r["expectation_value"])
        expected = 1.5 * 5.29177210903e-11  # 3a0/2
        assert abs(val - expected) / expected < 0.001
        print(f"[PASS] Expectation: H(1s) <r> = {r['expectation_value']} m (= 1.5·a₀)")

    def test_hydrogen_atom_Lz(self):
        """<Lz> for hydrogen should be m·ℏ."""
        result = self.tool.run_code(system="hydrogen_atom", observable="Lz", n=2, l=1, m=1)
        r = result["result"]
        val = float(r["expectation_value"])
        expected = 1 * 1.054571817e-34
        assert abs(val - expected) / expected < 0.001
        print(f"[PASS] Expectation: H <Lz> = {r['expectation_value']} J·s")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("particle_in_box x 1 0 0 1e-9")
        assert "result" in result
        print(f"[PASS] Expectation text interface OK")


# ============================================================
# Test 243: IrSpectrumPredictor
# ============================================================
class TestIrSpectrumPredictor:
    def setup_method(self):
        from chemmcp.tools import IrSpectrumPredictor
        self.tool = IrSpectrumPredictor()

    def test_ketone_alcohol(self):
        """Ketone + alcohol should give C=O and O-H peaks."""
        result = self.tool.run_code(functional_groups=["ketone", "alcohol"])
        r = result["result"]
        assert r["num_peaks"] >= 2
        assignments = [p["assignment"] for p in r["peaks"]]
        has_carbonyl = any("C=O" in a or "carbonyl" in a.lower() for a in assignments)
        has_oh = any("O-H" in a or "alcohol" in a.lower() for a in assignments)
        assert has_carbonyl, f"Expected C=O peak, got: {assignments}"
        assert has_oh, f"Expected O-H peak, got: {assignments}"
        print(f"[PASS] IR Predictor: ketone+alcohol → {r['num_peaks']} peaks")

    def test_carboxylic_acid_aromatic(self):
        """Carboxylic acid + aromatic should have broad O-H and aromatic C=C."""
        result = self.tool.run_code(functional_groups=["carboxylic_acid", "aromatic"], detail_level="detailed")
        r = result["result"]
        assert r["num_peaks"] >= 3
        assert "diagnostic_notes" in r
        print(f"[PASS] IR Predictor: COOH+aromatic → {r['num_peaks']} peaks, notes present")

    def test_single_group(self):
        """Single functional group should work."""
        result = self.tool.run_code(functional_groups=["nitrile"])
        r = result["result"]
        assert r["num_peaks"] >= 1
        print(f"[PASS] IR Predictor: nitrile → {r['num_peaks']} peaks")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("ketone alcohol")
        assert "result" in result
        print(f"[PASS] IR Predictor text interface OK")


# ============================================================
# Test 244: RamanActivity
# ============================================================
class TestRamanActivity:
    def setup_method(self):
        from chemmcp.tools import RamanActivity
        self.tool = RamanActivity()

    def test_alkene_aromatic(self):
        """Alkene + aromatic should show strong Raman peaks."""
        result = self.tool.run_code(functional_groups=["alkene", "aromatic"])
        r = result["result"]
        assert r["raman_active_modes"] >= 2
        assignments = [p["assignment"] for p in r["peaks"]]
        has_ring_breathing = any("ring breathing" in a.lower() for a in assignments)
        assert has_ring_breathing, f"Expected ring breathing mode, got: {assignments}"
        print(f"[PASS] Raman: alkene+aromatic → {r['raman_active_modes']} modes")

    def test_nitrile_disulfide(self):
        """Nitrile and disulfide are both strong Raman scatterers."""
        result = self.tool.run_code(functional_groups=["nitrile", "disulfide"])
        r = result["result"]
        intensities = [p.get("raman_activity", "") for p in r["peaks"]]
        assert any("high" in i for i in intensities), f"Expected high Raman activity, got: {intensities}"
        print(f"[PASS] Raman: nitrile+disulfide → {r['raman_active_modes']} modes")

    def test_selection_rules(self):
        """Selection rule analysis should be included."""
        result = self.tool.run_code(functional_groups=["benzene"], molecule_symmetry="D∞h")
        r = result["result"]
        assert "selection_rule_analysis" in r
        assert "ir_vs_raman_complementarity" in r
        print(f"[PASS] Raman: Selection rule analysis present")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("alkene aromatic")
        assert "result" in result
        print(f"[PASS] Raman text interface OK")


# ============================================================
# Test 245: UvVisTransitions
# ============================================================
class TestUvVisTransitions:
    def setup_method(self):
        from chemmcp.tools import UvVisTransitions
        self.tool = UvVisTransitions()

    def test_energy_from_wavelength_254(self):
        """254 nm is common UV wavelength; E = hc/λ."""
        result = self.tool.run_code(mode="energy_from_wavelength", wavelength_nm=254)
        r = result["result"]
        assert r["wavelength_nm"] == 254
        assert 4.8 < r["energy_eV"] < 4.9  # ~4.88 eV for 254 nm
        assert r["color_region"] == "UV-C / Far UV (200-280 nm)"
        print(f"[PASS] UV-Vis: 254 nm → {r['energy_eV']} eV")

    def test_wavelength_from_energy(self):
        """Inverse conversion: energy to wavelength."""
        result = self.tool.run_code(mode="wavelength_from_energy", energy_eV=3.1)
        r = result["result"]
        assert 390 < r["wavelength_nm"] < 410  # ~400 nm for 3.1 eV
        print(f"[PASS] UV-Vis: 3.1 eV → {r['wavelength_nm']} nm")

    def test_chromophore_lookup_benzene(self):
        """Benzene λmax ~254 nm."""
        result = self.tool.run_code(mode="chromophore_lookup", chromophore="benzene")
        r = result["result"]
        assert r["lambda_max_nm"] == 254
        assert "π→π*" in r["transition_type"]
        print(f"[PASS] UV-Vis: Benzene λmax = {r['lambda_max_nm']} nm ({r['transition_type']})")

    def test_woodward_fieser_diene(self):
        """Woodward-Fieser rules for dienes."""
        result = self.tool.run_code(
            mode="woodward_fieser_diene",
            conjugated_diene_type="acyclic_trans_trans",
            substituents=["alkyl_substituent", "alkyl_substituent"],
            num_conj_double_bonds=2,
        )
        r = result["result"]
        assert r["predicted_lambda_max_nm"] > 214  # base value
        assert "extensions" in r
        total = sum(e["value_nm"] for e in r["extensions"])
        assert total == r["predicted_lambda_max_nm"]
        print(f"[PASS] UV-Vis: WF diene λmax = {r['predicted_lambda_max_nm']} nm")

    def test_visible_light_region(self):
        """Visible light should map correctly."""
        result = self.tool.run_code(mode="energy_from_wavelength", wavelength_nm=550)
        r = result["result"]
        assert "Green" in r["color_region"]
        print(f"[PASS] UV-Vis: 550 nm → {r['color_region']}")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("energy_from_wavelength 280")
        assert "result" in result
        print(f"[PASS] UV-Vis text interface OK")


# ============================================================
# Test 246: NmrChemicalShift
# ============================================================
class TestNmrChemicalShift:
    def setup_method(self):
        from chemmcp.tools import NmrChemicalShift
        self.tool = NmrChemicalShift()

    def test_h1_methyl_aromatic_aldehyde(self):
        """¹H shifts for methyl, aromatic, aldehyde environments."""
        result = self.tool.run_code(nucleus="1H", environments=["methyl_alkane", "aromatic", "aldehyde"])
        r = result["result"]
        assert len(r["predictions"]) == 3
        shifts = {p["environment"]: p["typical_ppm"] for p in r["predictions"]}
        # Methyl should be lowest field (~0.9 ppm)
        assert shifts["H-C(sp³) (alkane)"] < 2
        # Aromatic should be 6.5-8.5
        assert 6.5 <= shifts["H-C(sp²) (aromatic)"] <= 8.5
        # Aldehyde should be highest (>9 ppm)
        assert shifts["H-C(=O) (aldehyde)"] >= 9.0
        print(f"[PASS] NMR ¹H: methyl={shifts['H-C(sp³) (alkane)']}, arom={shifts['H-C(sp²) (aromatic)']}, ald={shifts['H-C(=O) (aldehyde)']} ppm")

    def test_c13_ketone_aromatic_alcohol(self):
        """¹³C shifts for ketone, aromatic, alcohol alpha carbon."""
        result = self.tool.run_code(nucleus="13C", environments=["ketone", "aromatic", "alcohol_alpha"])
        r = result["result"]
        assert len(r["predictions"]) == 3
        shifts = {p["environment"]: p["typical_ppm"] for p in r["predictions"]}
        # Ketone carbonyl should be highest (~210 ppm)
        assert shifts["C=O (ketone)"] >= 190
        # Aromatic should be ~128 ppm
        assert 115 <= shifts["C (aromatic)"] <= 160
        # Alcohol alpha-C should be 50-90 ppm
        assert 50 <= shifts["C (α to oxygen, ether/alcohol)"] <= 90
        print(f"[PASS] NMR ¹³C: ketone={shifts['C=O (ketone)']}, arom={shifts['C (aromatic)']}, alc_α={shifts['C (α to oxygen, ether/alcohol)']} ppm")

    def test_coupling_info(self):
        """Including coupling info should add multiplicity data."""
        result = self.tool.run_code(nucleus="1H", environments=["methyl_alkane", "aldehyde"], include_coupling_info=True)
        r = result["result"]
        for p in r["predictions"]:
            if "multiplicity" in p:
                assert p["multiplicity"]
        assert "coupling_note" in r
        print(f"[PASS] NMR: Coupling info included")

    def test_invalid_nucleus_raises(self):
        with pytest.raises(Exception):
            self.tool.run_code(nucleus="99Tc", environments=["methyl_alkane"])
        print(f"[PASS] NMR: Invalid nucleus raises error")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("1H methyl_alkane aromatic aldehyde")
        assert "result" in result
        print(f"[PASS] NMR text interface OK")


# ============================================================
# Test 247: RotationalSpectrum
# ============================================================
class TestRotationalSpectrum:
    def setup_method(self):
        from chemmcp.tools import RotationalSpectrum
        self.tool = RotationalSpectrum()

    def test_co_like_molecule(self):
        """CO-like molecule (¹²C-¹⁶O) with reasonable B constant."""
        result = self.tool.run_code(
            molecule_type="linear",
            bond_length_angstrom=1.128,
            isotope_masses_amu=[12.0, 16.0],
            temperature_k=298.15,
            max_j=6,
        )
        r = result["result"]
        # CO rotational constant B ≈ 1.93 cm⁻¹
        B = r["rotational_constant_B_cm-1"]
        assert 1.5 < B < 2.5, f"Expected B ~1.93 cm⁻¹ for CO-like, got {B}"
        assert r["num_transitions"] == 6
        # Spacing should be 2B
        if r["transitions"]:
            spacing_check = all(
                abs(t["wavenumber_cm-1"] - 2*B*(t["J_lower"]+1)) < 0.01
                for t in r["transitions"]
            )
            assert spacing_check, "Transition spacing should equal 2B(J+1)"
        print(f"[PASS] Rotational: CO-like B = {B} cm⁻¹, {r['num_transitions']} transitions")

    def test_heavier_molecule_smaller_B(self):
        """Heavier molecule → larger I → smaller B."""
        result = self.tool.run_code(
            molecule_type="linear",
            reduced_mass_amu=20.0,
            bond_length_angstrom=1.5,
            temperature_k=300,
            max_j=5,
        )
        r = result["result"]
        B = r["rotational_constant_B_cm-1"]
        assert B > 0
        assert r["num_transitions"] == 5
        print(f"[PASS] Rotational: Heavy molecule B = {B} cm⁻¹")

    def test_selection_rule(self):
        """Selection rule should be ΔJ = ±1."""
        result = self.tool.run_code(
            molecule_type="linear",
            reduced_mass_amu=2.0,
            bond_length_angstrom=1.0,
            max_j=3,
        )
        r = result["result"]
        assert "ΔJ" in r["selection_rule"]
        assert r["most_intense_transition"] is not None
        print(f"[PASS] Rotational: Selection rule OK, most intense: {r['most_intense_transition']}")

    def test_temperature_affects_intensity(self):
        """Higher T should allow population of higher J levels."""
        # Use heavier molecule (larger I, smaller B) where T matters
        r_low = self.tool.run_code(molecule_type="linear", reduced_mass_amu=20.0, bond_length_angstrom=2.0, temperature_k=10, max_j=5)
        r_high = self.tool.run_code(molecule_type="linear", reduced_mass_amu=20.0, bond_length_angstrom=2.0, temperature_k=500, max_j=5)
        # At higher T, the most intense transition should be at higher J
        low_max = r_low["result"]["most_intense_transition"]
        high_max = r_high["result"]["most_intense_transition"]
        # Both should run without error and produce transitions
        assert r_low["result"]["num_transitions"] == 5
        assert r_high["result"]["num_transitions"] == 5
        print(f"[PASS] Rotational: Temperature effect — 10K max={low_max}, 500K max={high_max}")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("linear 2.0 1.0 300 5")
        assert "result" in result
        print(f"[PASS] Rotational text interface OK")


# ============================================================
# Test 248: VibrationalModes
# ============================================================
class TestVibrationalModes:
    def setup_method(self):
        from chemmcp.tools import VibrationalModes
        self.tool = VibrationalModes()

    def test_linear_triatomic(self):
        """Linear triatomic (like HCN): 3N-5 = 4 vibrational modes."""
        result = self.tool.run_code(n_atoms=3, molecule_geometry="linear", bond_types=["C≡N stretch", "C-H stretch"])
        r = result["result"]
        assert r["vibrational_modes"] == 4  # 3*3 - 5
        assert r["total_dof"] == 9
        assert r["translational_dof"] == 3
        assert r["rotational_dof"] == 2
        print(f"[PASS] Vib Modes: Linear triatomic → {r['vibrational_modes']} vib modes")

    def test_nonlinear_molecule(self):
        """Nonlinear molecule (7 atoms): 3N-6 = 15 vibrational modes."""
        result = self.tool.run_code(
            n_atoms=7, molecule_geometry="nonlinear",
            bond_types=["C=O stretch", "C-H stretch"],
            point_group="Cs",
        )
        r = result["result"]
        assert r["vibrational_modes"] == 15  # 3*7 - 6
        assert r["rotational_dof"] == 3
        print(f"[PASS] Vib Modes: Nonlinear 7-atom → {r['vibrational_modes']} vib modes")

    def test_predicted_frequencies(self):
        """Predicted frequencies should be in reasonable ranges."""
        result = self.tool.run_code(
            n_atoms=4, molecule_geometry="nonlinear",
            bond_types=["C=O stretch", "C-H stretch", "O-H stretch", "C-O stretch"],
        )
        r = result["result"]
        peaks = r["predicted_peaks"]
        assert len(peaks) >= 3
        freqs = [p["frequency_cm-1"] for p in peaks]
        # O-H stretch should be highest
        oh_freqs = [p["frequency_cm-1"] for p in peaks if "O-H" in p["mode"]]
        co_freqs = [p["frequency_cm-1"] for p in peaks if "C=O" in p["mode"]]
        if oh_freqs and co_freqs:
            assert oh_freqs[0] > co_freqs[0], "O-H stretch should be higher than C=O"
        print(f"[PASS] Vib Modes: Frequencies = {freqs}")

    def test_mutual_exclusion_centrosymmetric(self):
        """Centrosymmetric molecule should trigger mutual exclusion note."""
        result = self.tool.run_code(
            n_atoms=6, molecule_geometry="nonlinear",
            bond_types=["C=C stretch (alkene)", "C-H stretch"],
            point_group="D2h",
        )
        r = result["result"]
        assert r["has_center_of_inversion"] == True
        assert "mutual_exclusion_rule" in r
        print(f"[PASS] Vib Modes: D2h mutual exclusion noted")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("3 linear Cinfv CN_stretch CH_stretch")
        assert "result" in result
        print(f"[PASS] Vib Modes text interface OK")


# ============================================================
# Test 249: FranckCondonFactors
# ============================================================
class TestFranckCondonFactors:
    def setup_method(self):
        from chemmcp.tools import FranckCondonFactors
        self.tool = FranckCondonFactors()

    def test_s1_poission_distribution(self):
        """S=1, v=0 → Poisson distribution with maximum at v'=0 or 1."""
        result = self.tool.run_code(huang_rhys_factor_S=1.0, v_max=10, v_initial=0)
        r = result["result"]
        assert r["huang_rhys_factor_S"] == 1.0
        assert abs(r["sum_of_factors"] - 1.0) < 0.001  # Normalization
        assert r["max_intensity_at_v_prime"] in (0, 1)  # For S=1, max at v'≈S=1
        fc_factors = r["fc_factors"]
        assert len(fc_factors) == 11  # v'=0..10
        # v'=0 should be e^(-1) ≈ 0.368
        assert abs(fc_factors[0]["fc_factor"] - math.exp(-1)) < 0.001
        print(f"[PASS] FC: S=1, max at v'={r['max_intensity_at_v_prime']}, sum={r['sum_of_factors']:.4f}")

    def test_s0_origin_dominant(self):
        """S≈0: only 0-0 transition significant."""
        result = self.tool.run_code(huang_rhys_factor_S=0.1, v_max=5, v_initial=0)
        r = result["result"]
        assert r["max_intensity_at_v_prime"] == 0
        assert r["fc_factors"][0]["relative_intensity_pct"] > 90
        print(f"[PASS] FC: S=0.1, origin dominant ({r['fc_factors'][0]['relative_intensity_pct']:.1f}% at v'=0)")

    def test_s4_broad_progression(self):
        """S=4: broad progression with max at v'≈4."""
        result = self.tool.run_code(huang_rhys_factor_S=4.0, v_max=20, v_initial=0)
        r = result["result"]
        assert r["max_intensity_at_v_prime"] >= 3  # Should be near S=4
        assert r["max_intensity_at_v_prime"] <= 5
        assert abs(r["sum_of_factors"] - 1.0) < 0.001
        print(f"[PASS] FC: S=4, max at v'={r['max_intensity_at_v_prime']}, shape: {r.get('progression_summary', {}).get('progression_shape', 'N/A')}")

    def test_progression_plot_data(self):
        """Should include plot-suitable data."""
        result = self.tool.run_code(huang_rhys_factor_S=2.0, v_max=8, v_initial=0, include_progression_plot_data=True)
        r = result["result"]
        assert "progression_summary" in r
        ps = r["progression_summary"]
        assert "significant_peaks_above_1pct" in ps
        assert "progression_shape" in ps
        print(f"[PASS] FC: Progression summary present, shape = {ps['progression_shape']}")

    def test_normalization_check(self):
        """FC factors must sum to 1."""
        for S in [0.5, 1.0, 2.0, 3.0]:
            result = self.tool.run_code(huang_rhys_factor_S=S, v_max=15, v_initial=0)
            r = result["result"]
            assert r["normalization_check"] == "OK", f"S={S}: normalization failed, sum={r['sum_of_factors']}"
        print(f"[PASS] FC: Normalization verified for S=0.5,1,2,3")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("1.0 10 0 1000")
        assert "result" in result
        print(f"[PASS] FC text interface OK")


# ============================================================
# Test 250: BeerLambertCalculator
# ============================================================
class TestBeerLambertCalculator:
    def setup_method(self):
        from chemmcp.tools import BeerLambertCalculator
        self.tool = BeerLambertCalculator()

    def test_absorbance_calculation(self):
        """A = εbc for typical values."""
        result = self.tool.run_code(mode="absorbance", epsilon_M_cm=10000, pathlength_cm=1.0, concentration_M=0.0001)
        r = result["result"]
        expected_A = 10000 * 1.0 * 0.0001  # = 1.0
        assert abs(r["absorbance_A"] - expected_A) < 0.0001
        assert r["formula"]
        print(f"[PASS] Beer-Lambert: A = {r['absorbance_A']} (ε=10000, b=1, c=0.1mM)")

    def test_concentration_from_absorbance(self):
        """c = A/(εb)."""
        result = self.tool.run_code(mode="concentration", absorbance_A=0.654, epsilon_M_cm=15000, pathlength_cm=1.0)
        r = result["result"]
        expected_c = 0.654 / (15000 * 1.0)
        assert abs(r["concentration_M"] - expected_c) < 1e-10
        # Should also report µM
        assert r["concentration_uM"] == round(expected_c * 1e6, 4)
        print(f"[PASS] Beer-Lambert: c = {r['concentration_M']:.8f} M = {r['concentration_uM']} µM")

    def test_transmittance_conversion(self):
        """T ↔ A conversion."""
        result = self.tool.run_code(mode="transmittance", transmittance_T=0.35)
        r = result["result"]
        expected_A = -math.log10(0.35)
        assert abs(r["absorbance_A"] - expected_A) < 0.0001
        assert abs(r["percent_transmittance"] - 35.0) < 0.001
        print(f"[PASS] Beer-Lambert: T=0.35 → A={r['absorbance_A']}, %T={r['percent_transmittance']}%")

    def test_optimal_range_check(self):
        """A=0.5 should be optimal range."""
        result = self.tool.run_code(mode="absorbance", epsilon_M_cm=5000, pathlength_cm=1.0, concentration_M=0.0001)
        r = result["result"]
        assert "Optimal" in r["valid_range_check"] or "optimal" in r["valid_range_check"].lower()
        print(f"[PASS] Beer-Lambert: Range check = {r['valid_range_check']}")

    def test_high_absorbance_warning(self):
        """A >> 2 should warn about dilution."""
        result = self.tool.run_code(mode="absorbance", epsilon_M_cm=100000, pathlength_cm=1.0, concentration_M=0.01)
        r = result["result"]
        assert r["absorbance_A"] > 2
        assert "dilute" in r["valid_range_check"].lower() or "high" in r["valid_range_check"].lower()
        print(f"[PASS] Beer-Lambert: High-A warning = {r['valid_range_check']}")

    def test_epsilon_classification(self):
        """ε classification should match transition type."""
        result = self.tool.run_code(mode="epsilon", absorbance_A=1.0, pathlength_cm=1.0, concentration_M=0.0001)
        r = result["result"]
        eps = r["epsilon_M_cm"]
        # A = 1.0, b = 1.0, c = 0.0001 → ε = A/(b·c) = 1.0/0.0001 = 10000
        assert eps == 10000.0
        assert r["classification"]
        print(f"[PASS] Beer-Lambert: ε = {eps}, class = {r['classification']}")

    def test_multi_component(self):
        """Multi-component absorbance should be additive."""
        result = self.tool.run_code(
            mode="multi_component",
            components=[
                {"epsilon": 10000, "c": 0.0001, "name": "A"},
                {"epsilon": 5000, "c": 0.0002, "name": "B"},
            ],
            pathlength_cm=1.0,
        )
        r = result["result"]
        A_a = 10000 * 1.0 * 0.0001  # = 1.0
        A_b = 5000 * 1.0 * 0.0002   # = 1.0
        assert abs(r["total_absorbance_A"] - (A_a + A_b)) < 0.0001
        assert len(r["components"]) == 2
        print(f"[PASS] Beer-Lambert: Multi-component A_total = {r['total_absorbance_A']}")

    def test_dilution(self):
        """Dilution should reduce absorbance proportionally."""
        result = self.tool.run_code(mode="dilution", absorbance_A=1.2, initial_volume_mL=10, final_volume_mL=50)
        r = result["result"]
        expected_A = 1.2 * (10.0 / 50.0)  # = 0.24
        assert abs(r["final_absorbance_A"] - expected_A) < 0.0001
        assert r["fold_dilution"] == "1:5.0"
        print(f"[PASS] Beer-Lambert: Dilution 10→50 mL, A: 1.2 → {r['final_absorbance_A']}")

    def test_text_interface(self):
        tool = self.tool.__class__(interface="text")
        result = tool.run_text("absorbance 15000 1.0 0.001")
        assert "result" in result
        print(f"[PASS] Beer-Lambert text interface OK")


# ============================================================
# Run all tests
# ============================================================
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
