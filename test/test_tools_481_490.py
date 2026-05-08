"""
Test suite for MCP tools #481-490: Spectroscopy, NMR, Kinetics & Reaction Coordinate.
Run with: uv run python -m pytest test/test_tools_481_490.py -v
"""
import pytest
import math


class TestRamanSpectrumPredictor_482:
    """Tests for #482 RamanSpectrumPredictor."""

    def setup_method(self):
        from chemmcp.tools import RamanSpectrumPredictor
        self.tool = RamanSpectrumPredictor()

    def test_basic_prediction(self):
        result = self.tool.run_code(
            functional_groups=["benzene", "alkyne"],
            smiles=None,
            detail_level="standard"
        )
        r = result["result"]
        assert "molecule_type" in r
        assert "peaks" in r
        assert len(r["peaks"]) >= 3
        # Ring breathing should be present for benzene
        assignments = [p["assignment"] for p in r["peaks"]]
        assert any("breathing" in a.lower() or "benzene" in a.lower() or "aromatic" in a.lower() for a in assignments)

    def test_thiol_detection(self):
        """Thiols have diagnostic S-H stretch ~2575 cm⁻¹ unique to Raman."""
        result = self.tool.run_code(
            functional_groups=["thiol", "alkane"],
            detail_level="standard"
        )
        r = result["result"]
        shifts = [p["shift_cm-1"] for p in r["peaks"]]
        # S-H stretch should be around 2550-2600
        assert any(2550 <= s <= 2610 for s in shifts), f"S-H stretch not found in {shifts}"

    def test_text_input(self):
        result = self.tool.run_text("ketone alkene")
        assert "result" in result

    def test_depolarization_info(self):
        result = self.tool.run_code(
            functional_groups=["alkyne"],
            detail_level="detailed"
        )
        r = result["result"]
        for peak in r["peaks"]:
            assert "depolarization_ratio" in peak
            assert "activity" in peak


class TestUvVisSpectrum_483:
    """Tests for #483 UvVisSpectrum."""

    def setup_method(self):
        from chemmcp.tools import UvVisSpectrum
        self.tool = UvVisSpectrum()

    def test_enone_woodward_fieser(self):
        result = self.tool.run_code(
            chromophore_type="enone",
            substituents=["alkyl_beta"],
            conjugation_length=2,
            solvent="ethanol",
            calculation_mode="woodward_fieser"
        )
        r = result["result"]
        assert 210 <= r["lambda_max_nm"] <= 250  # enone base 215 + beta-alkyl 12 ≈ 227
        assert r["molar_absorptivity"] > 0
        assert "transition_type" in r

    def test_aromatic_lookup(self):
        result = self.tool.run_code(
            chromophore_type="aromatic",
            substituents=["OH"],
            solvent="water",
            calculation_mode="auto"
        )
        r = result["result"]
        assert 250 <= r["lambda_max_nm"] <= 300  # Phenol ~270 nm

    def test_polyene_estimate(self):
        result = self.tool.run_code(
            chromophore_type="polyene",
            conjugation_length=5,
            calculation_mode="estimate"
        )
        r = result["result"]
        # Linear polyene formula: λ ≈ 114 + 47*5 = 349 nm (visible region)
        assert 300 <= r["lambda_max_nm"] <= 400

    def test_text_input(self):
        result = self.tool.run_text("diene OR_group 3 hexane woodward_fieser")
        assert "result" in result


class TestNmrShielding_484:
    """Tests for #484 NmrShielding."""

    def setup_method(self):
        from chemmcp.tools import NmrShielding
        self.tool = NmrShielding()

    def test_proton_with_chlorine(self):
        """CH₃-Cl type proton should be deshielded to ~3.0-3.5 ppm."""
        result = self.tool.run_code(
            nucleus="1H",
            atom_environment=["-Cl", "-CH3"],
            hybridization="sp3",
            reference_compound="TMS"
        )
        r = result["result"]
        assert 2.5 <= r["chemical_shift_ppm"] <= 4.5
        assert r["shielding_constant"] < 0  # Deshielded → negative σ

    def test_carbonyl_carbon_13c(self):
        """Carbonyl ¹³C should be highly deshielded (>160 ppm)."""
        result = self.tool.run_code(
            nucleus="13C",
            atom_environment=["carbonyl"],
            hybridization="sp2"
        )
        r = result["result"]
        assert r["chemical_shift_ppm"] >= 160

    def test_aromatic_proton(self):
        """Aromatic protons should be in 6.5-8.5 ppm range."""
        result = self.tool.run_code(
            nucleus="1H",
            atom_environment=["aromatic", "ortho_EWG"],
            hybridization="sp2"
        )
        r = result["result"]
        assert 6.0 <= r["chemical_shift_ppm"] <= 10.0

    def test_shielding_constant_relation(self):
        """δ = σ_ref − σ_sample; with TMS reference (σ_ref=0): δ = −σ"""
        result = self.tool.run_code(
            nucleus="1H",
            atom_environment=["-Cl"],
            hybridization="sp3"
        )
        r = result["result"]
        assert abs(r["chemical_shift_ppm"] + r["shielding_constant"]) < 0.01

    def test_text_input(self):
        result = self.tool.run_text("1H -NO2 -CH3 sp3 TMS standard")
        assert "result" in result


class TestSpinSpinCoupling_485:
    """Tests for #485 SpinSpinCoupling."""

    def setup_method(self):
        from chemmcp.tools import SpinSpinCoupling
        self.tool = SpinSpinCoupling()

    def test_vicinal_sp3_coupling(self):
        """³JHH(sp³) free rotation should be ~5-14 Hz, typical ~7 Hz."""
        result = self.tool.run_code(
            coupled_nuclei="HH",
            num_bonds=3,
            hybridization="sp3",
            geometry="default"
        )
        r = result["result"]
        assert 3 <= r["j_coupling_hz"] <= 16

    def test_trans_alkene_karplus(self):
        """Trans alkene ³JHH should be large (~10-19 Hz)."""
        result = self.tool.run_code(
            coupled_nuclei="HH",
            num_bonds=3,
            hybridization="sp2",
            geometry="trans",
            dihedral_angle_deg=180.0
        )
        r = result["result"]
        assert r["j_coupling_hz"] >= 10  # Trans J is large (Karplus at 180°)
        assert "karplus_analysis" in r

    def test_cis_alkene_smaller_j(self):
        """Cis alkene ³JHH should be smaller than trans."""
        result_cis = self.tool.run_code(
            coupled_nuclei="HH", num_bonds=3, hybridization="sp2",
            geometry="cis", dihedral_angle_deg=0.0
        )
        result_trans = self.tool.run_code(
            coupled_nuclei="HH", num_bonds=3, hybridization="sp2",
            geometry="trans", dihedral_angle_deg=180.0
        )
        assert result_cis["result"]["j_coupling_hz"] < result_trans["result"]["j_coupling_hz"]

    def test_one_bond_ch_coupling(self):
        """¹JCH sp² should be ~150-170 Hz."""
        result = self.tool.run_code(
            coupled_nuclei="CH",
            num_bonds=1,
            hybridization="sp2"
        )
        r = result["result"]
        assert 140 <= r["j_coupling_hz"] <= 180

    def test_karplus_curve_properties(self):
        """Karplus curve: J(0°) and J(180°) should be maxima, J(90°) minimum."""
        j0 = self.tool._karplus_calc(0, "sp3", "default")
        j90 = self.tool._karplus_calc(90, "sp3", "default")
        j180 = self.tool._karplus_calc(180, "sp3", "default")
        # Both 0° and 180° should be larger than 90°
        assert j0 > j90
        assert j180 > j90

    def test_text_input(self):
        result = self.tool.run_text("CC 2 sp2 geminal detailed")
        assert "result" in result


class TestRateLawIntegrator_486:
    """Tests for #486 RateLawIntegrator."""

    def setup_method(self):
        from chemmcp.tools import RateLawIntegrator
        self.tool = RateLawIntegrator()

    def test_first_order_decay(self):
        """First-order: [A] = [A]₀·exp(-kt). After t=0, [A]=[A]₀."""
        result = self.tool.run_code(
            order=1, k=0.001, initial_concentration=1.0, time_s=0
        )
        assert abs(result["result"]["concentration_M"] - 1.0) < 1e-6

    def test_first_order_half_life(self):
        """First-order t₁/₂ = ln(2)/k."""
        result = self.tool.run_code(order=1, k=0.693, initial_concentration=1.0)
        t_half = result["result"]["half_life_s"]
        assert abs(t_half - 1.0) < 0.02  # ln(2)/0.693 ≈ 1.001

    def test_zero_order_depletion(self):
        """Zero-order: [A] = [A]₀ - kt. Should reach zero at t=[A]₀/k."""
        result = self.tool.run_code(
            order=0, k=0.01, initial_concentration=1.0, time_s=100
        )
        assert result["result"]["concentration_M"] == 0.0  # Exactly depleted

    def test_second_order_half_life(self):
        """Second-order t₁/₂ = 1/(k·[A]₀)."""
        k = 2.0
        A0 = 0.5
        result = self.tool.run_code(order=2, k=k, initial_concentration=A0)
        expected_t_half = 1.0 / (k * A0)  # = 1.0
        assert abs(result["result"]["half_life_s"] - expected_t_half) < 0.02

    def test_target_fraction_time(self):
        """Time to reach 50% remaining should equal half-life."""
        result = self.tool.run_code(
            order=1, k=1.0, initial_concentration=1.0, target_fraction=0.5
        )
        t_frac = result["result"].get("time_to_fraction_0.5_s")
        t_half = self.tool._calc_half_life(1, 1.0, 1.0)
        assert abs(t_frac - t_half) < 0.01

    def test_fraction_remaining_consistency(self):
        """fraction_remaining = [A]/[A]₀."""
        A0 = 2.0
        result = self.tool.run_code(order=1, k=0.5, initial_concentration=A0, time_s=1.0)
        At = result["result"]["concentration_M"]
        frac = result["result"]["fraction_remaining"]
        assert abs(frac - At / A0) < 1e-6

    def test_text_input(self):
        result = self.tool.run_text("1 0.001 1.0 3600")
        assert "result" in result


class TestArrheniusCalculator_487:
    """Tests for #487 ArrheniusCalculator."""

    def setup_method(self):
        from chemmcp.tools import ArrheniusCalculator
        self.tool = ArrheniusCalculator()

    def test_calculate_k_basic(self):
        """k = A·exp(-Ea/RT) must give positive k."""
        result = self.tool.run_code(
            mode="calculate_k",
            ea_kj_mol=50.0,
            pre_exponential_A=1e11,
            temperature_k=298.15
        )
        r = result["result"]
        assert r["rate_constant_k"] > 0
        assert r["ea_kj_mol"] == 50.0

    def test_high_ea_gives_small_k(self):
        """Higher Ea → smaller k at same T."""
        r_low = self.tool.run_code(mode="calculate_k", ea_kj_mol=50, pre_exponential_A=1e11, temperature_k=298.15)["result"]
        r_high = self.tool.run_code(mode="calculate_k", ea_kj_mol=100, pre_exponential_A=1e11, temperature_k=298.15)["result"]
        assert r_low["rate_constant_k"] > r_high["rate_constant_k"]

    def test_higher_T_gives_larger_k(self):
        """Higher T → larger k at same Ea."""
        r_cold = self.tool.run_code(mode="calculate_k", ea_kj_mol=60, pre_exponential_A=1e11, temperature_k=280)["result"]
        r_hot = self.tool.run_code(mode="calculate_k", ea_kj_mol=60, pre_exponential_A=1e11, temperature_k=350)["result"]
        assert r_hot["rate_constant_k"] > r_cold["rate_constant_k"]

    def test_two_point_ea(self):
        """Two-point method should return positive Ea for normal kinetics."""
        result = self.tool.run_code(
            mode="two_point_ea",
            T1_k=300.0, k1=1e-4,
            T2_k=320.0, k2=1e-3
        )
        r = result["result"]
        assert r["ea_kj_mol"] > 0
        assert r["pre_exponential_A_estimated"] > 0

    def test_arrhenius_plot_regression(self):
        """Arrhenius plot with synthetic data should recover known parameters."""
        # Generate synthetic data with known Ea=75 kJ/mol, A=1e13
        R_local = 8.314462618
        Ea_syn = 75.0  # kJ/mol
        A_syn = 1e13
        temps = [298.15, 308.15, 318.15, 328.15, 338.15]
        ks = [A_syn * math.exp(-Ea_syn * 1000 / (R_local * T)) for T in temps]

        result = self.tool.run_code(
            mode="arrhenius_plot",
            temperatures_k=",".join(str(t) for t in temps),
            rate_constants=",".join(str(k) for k in ks)
        )
        r = result["result"]
        # Recovered Ea should be close to 75 kJ/mol
        assert abs(r["ea_kj_mol"] - Ea_syn) / Ea_syn < 0.05  # Within 5%
        assert r["r_squared"] > 0.99

    def test_text_input(self):
        result = self.tool.run_text("calculate_k 75 1e13 298.15")
        assert "result" in result


class TestEyringEquation_488:
    """Tests for #488 EyringEquation."""

    def setup_method(self):
        from chemmcp.tools import EyringEquation
        self.tool = EyringEquation()

    def test_calculate_k_from_dg(self):
        """k from ΔG‡ should be positive and match TST prefactor scaling."""
        result = self.tool.run_code(
            mode="calculate_k",
            delta_g_kj_mol=75.0,
            temperature_k=298.15
        )
        r = result["result"]
        assert r["rate_constant_s-1"] > 0
        assert r["delta_g_dd_kj_mol"] == 75.0
        assert "eyring_equation" in r

    def test_calculate_k_from_dh_ds(self):
        """From ΔH‡ and ΔS‡: ΔG‡ = ΔH‡ − TΔS‡."""
        result = self.tool.run_code(
            mode="calculate_k",
            delta_g_kj_mol=0,  # Not used when dH/dS provided
            delta_h_kj_mol=72.0,
            delta_s_j_mol_K=-10.0,
            temperature_k=298.15
        )
        r = result["result"]
        # ΔG‡ = 72 - 298.15*(-10/1000) = 72 + 2.98 = 74.98 kJ/mol
        expected_dg = 72.0 - 298.15 * (-10.0) / 1000.0
        assert abs(r["delta_g_dd_kj_mol"] - expected_dg) < 0.1
        assert r["rate_constant_s-1"] > 0

    def test_negative_ds_ordered_ts(self):
        """Negative ΔS‡ gives slower k than same ΔG‡ with ΔS‡=0."""
        r_neg_s = self.tool.run_code(
            mode="calculate_k", delta_h_kj_mol=70.0, delta_s_j_mol_K=-20.0, temperature_k=298.15
        )["result"]
        r_zero_s = self.tool.run_code(
            mode="calculate_k", delta_h_kj_mol=70.0, delta_s_j_mol_K=0.0, temperature_k=298.15
        )["result"]
        # More negative ΔS‡ → higher ΔG‡ → slower
        assert r_neg_s["rate_constant_s-1"] < r_zero_s["rate_constant_s-1"]

    def test_eyring_plot_fit(self):
        """Eyring plot should extract reasonable activation parameters."""
        # Generate data with known ΔH‡=72 kJ/mol, ΔS‡=-10 J/(mol·K)
        kB_l = 1.380649e-23
        h_l = 6.62607015e-34
        R_l = 8.314462618
        dH = 72000  # J/mol
        dS = -10.0  # J/(mol·K)

        temps = [298.15, 308.15, 318.15, 328.15, 338.15]
        ks = []
        for T in temps:
            dg = dH - T * dS
            pref = kB_l * T / h_l
            ks.append(pref * math.exp(-dg / (R_l * T)))

        result = self.tool.run_code(
            mode="eyring_plot",
            temperatures_k=",".join(str(t) for t in temps),
            rate_constants=",".join(str(k) for k in ks)
        )
        r = result["result"]
        assert abs(r["delta_h_dd_kj_mol"] - 72.0) / 72.0 < 0.05  # Within 5%
        assert r["r_squared"] > 0.99

    def test_compare_arrhenius(self):
        """Compare mode should return both Arrhenius and TST equivalent params."""
        result = self.tool.run_code(
            mode="compare_arrhenius",
            arrhenius_ea_kj_mol=75.0,
            arrhenius_A=1e13,
            temperature_k=298.15
        )
        r = result["result"]
        assert "arrhenius" in r
        assert "tst_equivalent" in r
        assert r["arrhenius"]["ea_kj_mol"] == 75.0

    def test_text_input(self):
        result = self.tool.run_text("calculate_k 80.0")
        assert "result" in result


class TestReactionCoordinate_490:
    """Tests for #490 ReactionCoordinate."""

    def setup_method(self):
        from chemmcp.tools import ReactionCoordinate
        self.tool = ReactionCoordinate()

    def test_analytical_exothermic(self):
        """Exothermic reaction: products lower than reactants."""
        result = self.tool.run_code(
            input_mode="analytical",
            reactant_energy_kj_mol=0.0,
            ts_energy_kj_mol=75.0,
            product_energy_kj_mol=-20.0,
            n_points=50
        )
        r = result["result"]
        assert r["reaction_type"] == "exothermic"
        assert r["forward_barrier_kj_mol"] == 75.0
        assert r["reverse_barrier_kj_mol"] == 95.0  # 75 - (-20) = 95
        assert r["reaction_energy_kj_mol"] == -20.0
        assert 0 < r["ts_position_normalized"] < 1
        assert len(r["irc_points"]) >= 40

    def test_analytical_endergonic(self):
        """Endergonic reaction: products higher than reactants."""
        result = self.tool.run_code(
            input_mode="analytical",
            reactant_energy_kj_mol=0.0,
            ts_energy_kj_mol=80.0,
            product_energy_kj_mol=30.0,
            n_points=30
        )
        r = result["result"]
        assert "endergonic" in r["reaction_type"]
        assert r["forward_barrier_kj_mol"] == 80.0
        assert r["reverse_barrier_kj_mol"] == 50.0  # 80 - 30 = 50
        assert r["reaction_energy_kj_mol"] == 30.0

    def test_manual_profile(self):
        """Manual energy profile should find correct TS position."""
        profile = [
            (0.0, 0.0),
            (0.2, 25.0),
            (0.4, 55.0),
            (0.6, 78.0),   # TS
            (0.8, 45.0),
            (1.0, 10.0),
        ]
        result = self.tool.run_code(
            input_mode="manual",
            energy_profile=profile,
            detail_level="detailed"
        )
        r = result["result"]
        assert r["ts_energy_kj_mol"] == 78.0
        assert r["max_energy_point"]["energy_kj_mol"] == 78.0
        assert r["max_energy_point"]["coord"] == 0.6

    def test_thermoneutral_reaction(self):
        """Thermoneutral: symmetric barriers."""
        result = self.tool.run_code(
            input_mode="analytical",
            reactant_energy_kj_mol=0.0,
            ts_energy_kj_mol=65.0,
            product_energy_kj_mol=0.0,
            n_points=20
        )
        r = result["result"]
        assert "thermoneutral" in r["reaction_type"]
        assert abs(r["forward_barrier_kj_mol"] - r["reverse_barrier_kj_mol"]) < 0.01

    def test_hammond_postulate_exothermic(self):
        """Exothermic → late TS (position > 0.5)."""
        result = self.tool.run_code(
            input_mode="analytical",
            reactant_energy_kj_mol=0.0,
            ts_energy_kj_mol=60.0,
            product_energy_kj_mol=-50.0
        )
        r = result["result"]
        assert r["ts_position_normalized"] > 0.5  # Late TS for exothermic

    def test_hammond_postulate_endergonic(self):
        """Endergonic → early TS (position < 0.5)."""
        result = self.tool.run_code(
            input_mode="analytical",
            reactant_energy_kj_mol=0.0,
            ts_energy_kj_mol=70.0,
            product_energy_kj_mol=50.0
        )
        r = result["result"]
        assert r["ts_position_normalized"] < 0.5  # Early TS for endergonic

    def test_keq_calculation(self):
        """K_eq should favor products for exothermic reactions."""
        result = self.tool.run_code(
            input_mode="analytical",
            reactant_energy_kj_mol=0.0,
            ts_energy_kj_mol=50.0,
            product_energy_kj_mol=-30.0,
            temperature_k=298.15
        )
        r = result["result"]
        keq = r.get("equilibrium_constant", {})
        assert keq.get("K_eq", 0) > 1  # Products favored

    def test_text_input_analytical(self):
        result = self.tool.run_text("analytical 0 75 -20 298.15 50 standard")
        assert "result" in result

    def test_text_input_manual(self):
        result = self.tool.run_text("manual 0,0 0.5,80 1.0,10 basic")
        assert "result" in result


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
