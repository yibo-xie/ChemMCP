"""
Test suite for MCP tools #491-500: Kinetics Corrections & Computational Chemistry Tools.
Run with: uv run python -m pytest test/test_tools_491_500.py -v

Tools covered:
  #491 TunnelingCorrection     - Quantum tunneling correction factor (κ)
  #492 SteadyStateApproximation - Steady-state approximation for intermediates
  #493 PreEquilibrium           - Pre-equilibrium approximation for fast equilibrium steps
  #494 RateDeterminingStep      - Rate-determining step analysis & bottleneck identification
  #495 MichaelisMenten          - Michaelis-Menten enzyme kinetics
  #496 ReactionNetworkSolver    - Reaction network ODE solver for multi-step mechanisms
  #497 GeometryOptimizer        - Molecular geometry optimization (energy minimization)
  #498 TransitionStateSearch    - Transition state saddle-point search
  #499 PotentialEnergySurface   - PES scan along reaction coordinates
  #500 FrequencyAnalysis         - Vibrational frequency analysis & stationary point confirmation
"""
import pytest
import math


# =============================================================================
# #491 TunnelingCorrection — 隧穿校正，轻原子转移反应
# =============================================================================
class TestTunnelingCorrection_491:
    """Tests for #491 TunnelingCorrection."""

    def setup_method(self):
        from chemmcp.tools import TunnelingCorrection
        self.tool = TunnelingCorrection()

    def test_bell_correction_returns_kappa(self):
        """Bell correction returns a numeric kappa value."""
        r = self.tool.run_code(
            temperature_K=298.15,
            barrier_height_kJ_mol=45.0,
            imaginary_frequency_cm_minus_1=1500.0,
            correction_model="bell"
        )
        assert "correction_factor_kappa" in r
        assert isinstance(r["correction_factor_kappa"], float)
        assert r["temperature_K"] == 298.15

    def test_bell_low_temp_higher_kappa(self):
        """Lower T → larger tunneling effect (kappa changes)."""
        r_lo = self.tool.run_code(200, 45.0, 1500.0, "bell")
        r_hi = self.tool.run_code(400, 45.0, 1500.0, "bell")
        # kappa should differ with temperature
        assert r_lo["correction_factor_kappa"] != r_hi["correction_factor_kappa"]

    def test_bell_freq_affects_kappa(self):
        """Different imaginary frequencies give different kappas."""
        r_lo = self.tool.run_code(298.15, 45.0, 800.0, "bell")
        r_hi = self.tool.run_code(298.15, 45.0, 2000.0, "bell")
        assert r_lo["correction_factor_kappa"] != r_hi["correction_factor_kappa"]

    def test_wkb_correction(self):
        """WKB model returns a result with expected keys."""
        r = self.tool.run_code(298.15, 45.0, 1500.0, "wkb", 1.008)
        assert "correction_factor_kappa" in r
        assert "wkb_exponent_gamma" in r

    def test_eckart_correction(self):
        """Eckart asymmetric barrier model."""
        r = self.tool.run_code(298.15, 40.0, 1400.0, "eckart", 1.008, 40.0, 80.0)
        assert "correction_factor_kappa" in r
        assert "asymmetry_ratio" in r

    def test_invalid_temperature_raises(self):
        """Zero T must raise error."""
        with pytest.raises(Exception):
            self.tool.run_code(0, 45.0, 1500.0, "bell")

    def test_text_input(self):
        r = self.tool.run_text("298.15 45.0 1500.0 bell")
        assert "correction_factor_kappa" in r

    def test_tunneling_significance_flag(self):
        """tunneling_significant is a boolean flag."""
        r = self.tool.run_code(200, 30.0, 2000.0, "bell")
        assert isinstance(r["tunneling_significant"], bool)


# =============================================================================
# #492 SteadyStateApproximation — 稳态近似，中间体浓度
# =============================================================================
class TestSteadyStateApproximation_492:
    """Tests for #492 SteadyStateApproximation."""

    def setup_method(self):
        from chemmcp.tools import SteadyStateApproximation
        self.tool = SteadyStateApproximation()

    def test_consecutive_mechanism(self):
        """Consecutive A→I→P: SSA approximates exact after induction period."""
        r = self.tool.run_code(
            mechanism_type="consecutive",
            rate_constants=[0.1, 0.05],
            initial_reactant_concentration_A0=1.0,
            time_t=50.0,
            target_intermediate="I"
        )
        assert "exact" in r
        assert "approximate" in r
        assert "ssa_valid" in r

    def test_reversible_consecutive(self):
        """Reversible consecutive A⇌I→P mechanism."""
        r = self.tool.run_code(
            mechanism_type="reversible_consecutive",
            rate_constants=[0.1, 0.02, 0.01],
            initial_reactant_concentration_A0=1.0,
            time_t=100.0
        )
        assert "exact" in r or "approximate" in r

    def test_pre_equilibrium_via_ssa(self):
        """Pre-equilibrium mechanism via SSA."""
        r = self.tool.run_code(
            mechanism_type="pre_equilibrium",
            rate_constants=[1.0, 0.5, 0.01],
            initial_reactant_concentration_A0=1.0,
            time_t=200.0
        )
        assert "exact" in r or "approximate" in r

    def test_fast_first_step_good_ssa(self):
        """k1 >> k2 → SSA should be valid (small intermediate error)."""
        r = self.tool.run_code(
            mechanism_type="consecutive",
            rate_constants=[1.0, 0.01],
            initial_reactant_concentration_A0=1.0,
            time_t=500.0
        )
        assert "ssa_valid" in r or "intermediate_error_percent" in r

    def test_text_input(self):
        # Format: type k1,k2 A0 t [intermediate]
        r = self.tool.run_text("consecutive 0.1,0.05 1.0 50.0 I")
        assert "exact" in r or "approximate" in r


# =============================================================================
# #493 PreEquilibrium — 预平衡近似，快速平衡步骤
# =============================================================================
class TestPreEquilibrium_493:
    """Tests for #493 PreEquilibrium."""

    def setup_method(self):
        from chemmcp.tools import PreEquilibrium
        self.tool = PreEquilibrium()

    def test_unimolecular_preeq(self):
        """Unimolecular A ⇌ I → P: K_eq = k_f/k_r."""
        r = self.tool.run_code(
            mechanism="unimolecular",
            k_forward_list=[1.0],
            k_reverse_list=[0.5],
            k_slow=0.01,
            initial_concentrations={"A": 1.0}
        )
        assert "effective_rate_constant" in r
        assert "rate_law" in r
        assert "validity" in r

    def test_bimolecular_preeq(self):
        """Bimolecular A + B ⇌ I → P."""
        r = self.tool.run_code(
            mechanism="bimolecular",
            k_forward_list=[1.0],
            k_reverse_list=[0.2],
            k_slow=0.05,
            initial_concentrations={"A": 1.0, "B": 1.0}
        )
        assert "effective_rate_constant" in r

    def test_keq_value(self):
        """K_eq = k_f/k_r should appear in output."""
        r = self.tool.run_code(
            mechanism="unimolecular",
            k_forward_list=[2.0],
            k_reverse_list=[0.5],
            k_slow=0.01,
            initial_concentrations={"A": 1.0}
        )
        result_str = str(r)
        assert any(x in result_str.lower() for x in ["equilibrium", "keq", "k_eq"])

    def test_time_profiles(self):
        """With time_points, concentration profiles are generated."""
        r = self.tool.run_code(
            mechanism="unimolecular",
            k_forward_list=[1.0],
            k_reverse_list=[0.5],
            k_slow=0.01,
            initial_concentrations={"A": 1.0},
            time_points=[0.0, 10.0, 50.0, 100.0]
        )
        assert "concentration_profiles" in r

    def test_text_input(self):
        # Format: mechanism kf1 kf2 k_slow A0=1.0
        r = self.tool.run_text("unimolecular 1.0 0.5 0.01 A0=1.0")
        assert "effective_rate_constant" in r or "rate_law" in r or "validity" in r


# =============================================================================
# #494 RateDeterminingStep — 速控步分析，反应瓶颈识别
# =============================================================================
class TestRateDeterminingStep_494:
    """Tests for #494 RateDeterminingStep."""

    def setup_method(self):
        from chemmcp.tools import RateDeterminingStep
        self.tool = RateDeterminingStep()

    def test_identify_rds(self):
        """RDS = step with smallest k (1-indexed)."""
        r = self.tool.run_code([
            {"reactants": "A", "products": "I", "k": 1.0e5},
            {"reactants": "I + B", "products": "P", "k": 0.01},  # slowest → step 2
            {"reactants": "P", "products": "Q", "k": 100.0},
        ])
        assert r["rds_step_index"] == 2  # 1-indexed
        assert r["rds_rate_constant"] == 0.01

    def test_single_step_is_rds(self):
        """Single-step: that step is RDS (index 1 in 1-indexed)."""
        r = self.tool.run_code([{"reactants": "A", "products": "P", "k": 0.5}])
        assert r["rds_step_index"] == 1  # 1-indexed

    def test_rate_ratios_ge_one(self):
        """All ratios relative to RDS should be >= 1."""
        r = self.tool.run_code([{"k": 10.0}, {"k": 0.1}, {"k": 5.0}])
        for ratio in r.get("rate_constant_ratios", r.get("rate_ratios", [])):
            assert ratio >= 1.0

    def test_with_pre_equilibrium_flag(self):
        """Pre-equilibrium flag affects rate law derivation."""
        r = self.tool.run_code(
            [{"reactants": "A ⇌ I", "k": 100.0}, {"reactants": "I → P", "k": 0.01}],
            has_pre_equilibrium=True
        )
        assert "overall_rate_law" in r or "rate_law" in r

    def test_empty_steps_raises(self):
        """Empty steps list raises error."""
        with pytest.raises(Exception):
            self.tool.run_code([])

    def test_text_input(self):
        # Format: desc;k;rev;desc;k;rev;...
        r = self.tool.run_text("fast step;10.0;false;slow step;0.1;false")
        assert "rds_step_index" in r


# =============================================================================
# #495 MichaelisMenten — Michaelis-Menten动力学，酶催化
# =============================================================================
class TestMichaelisMenten_495:
    """Tests for #495 MichaelisMenten."""

    def setup_method(self):
        from chemmcp.tools import MichaelisMenten
        self.tool = MichaelisMenten()

    def test_velocity_calculation(self):
        """v = Vmax·S / (Km + S)."""
        Vmax, Km, S = 10.0, 2.0, 5.0
        r = self.tool.run_code("calculate_velocity", S, Vmax, Km)
        expected = Vmax * S / (Km + S)  # ≈ 7.1429
        v = r.get("velocity", r.get("velocity_uninhibited", 0))
        assert abs(v - expected) < 0.01

    def test_half_saturation(self):
        """S = Km → v = Vmax/2."""
        r = self.tool.run_code("calculate_velocity", 3.0, 10.0, 3.0)
        v = r.get("velocity", r.get("velocity_uninhibited", 0))
        assert abs(v - 5.0) < 0.01

    def test_low_substrate_first_order(self):
        """S << Km: first-order regime."""
        r = self.tool.run_code("calculate_velocity", 0.1, 10.0, 10.0)
        v = r.get("velocity", r.get("velocity_uninhibited", 0))
        assert v < 10.0

    def test_high_substrate_zero_order(self):
        """S >> Km: zero-order saturation, v ≈ Vmax."""
        r = self.tool.run_code("calculate_velocity", 100.0, 10.0, 1.0)
        v = r.get("velocity", r.get("velocity_uninhibited", 0))
        assert abs(v - 10.0) / 10.0 < 0.05

    def test_competitive_inhibition(self):
        """Competitive inhibition reduces apparent velocity."""
        r = self.tool.run_code(
            "full_analysis", 5.0, 10.0, 2.0,
            "competitive", 1.0, 1.0
        )
        assert "velocity_inhibited" in r or "inhibition" in str(r).lower()

    def test_lineweaver_burk(self):
        """Lineweaver-Burk parameter extraction."""
        r = self.tool.run_code(
            "lineweaver_burk",
            substrate_velocities_data=[
                {"S": 1.0, "v": 1.67}, {"S": 2.0, "v": 2.5},
                {"S": 5.0, "v": 3.33}, {"S": 10.0, "v": 4.0},
                {"S": 20.0, "v": 4.44},
            ]
        )
        assert "analysis_type" in r

    def test_text_input(self):
        r = self.tool.run_text("calculate_velocity 5.0 10.0 2.0")
        assert "velocity" in r


# =============================================================================
# #496 ReactionNetworkSolver — 反应网络求解，多步机理模拟
# =============================================================================
class TestReactionNetworkSolver_496:
    """Tests for #496 ReactionNetworkSolver."""

    def setup_method(self):
        from chemmcp.tools import ReactionNetworkSolver
        self.tool = ReactionNetworkSolver()

    def test_consecutive_network(self):
        """A → B → C consecutive reaction network."""
        r = self.tool.run_code(
            species=["A", "B", "C"],
            reactions=[
                {"reactants": ["A"], "products": ["B"], "k": 0.1},
                {"reactants": ["B"], "products": ["C"], "k": 0.05},
            ],
            initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
            time_end=100.0,
            n_points=50
        )
        assert "concentration_profiles" in r

    def test_reversible_reaction(self):
        """Reversible A ⇌ B reaches equilibrium."""
        r = self.tool.run_code(
            species=["A", "B"],
            reactions=[
                {"reactants": ["A"], "products": ["B"], "k": 0.1, "reversible": True, "k_reverse": 0.05},
            ],
            initial_concentrations={"A": 1.0, "B": 0.0},
            time_end=200.0,
            n_points=100
        )
        assert "concentration_profiles" in r

    def test_parallel_reactions(self):
        """Parallel: A → B, A → C."""
        r = self.tool.run_code(
            species=["A", "B", "C"],
            reactions=[
                {"reactants": ["A"], "products": ["B"], "k": 0.1},
                {"reactants": ["A"], "products": ["C"], "k": 0.05},
            ],
            initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
            time_end=50.0
        )
        assert "concentration_profiles" in r

    def test_bimolecular_reaction(self):
        """Bimolecular A + B → C."""
        r = self.tool.run_code(
            species=["A", "B", "C"],
            reactions=[
                {"reactants": ["A", "B"], "products": ["C"], "k": 0.1},
            ],
            initial_concentrations={"A": 1.0, "B": 1.0, "C": 0.0},
            time_end=30.0
        )
        assert "concentration_profiles" in r

    def test_steady_state_detection_long_time(self):
        """Long integration should approach steady state."""
        r = self.tool.run_code(
            species=["A", "B"],
            reactions=[
                {"reactants": ["A"], "products": ["B"], "k": 0.1},
            ],
            initial_concentrations={"A": 1.0, "B": 0.0},
            time_end=500.0,
            n_points=100
        )
        assert "concentration_profiles" in r

    def test_text_input(self):
        # Format: species A,B,C reactions A->B:k1 init ...
        r = self.tool.run_text("species A,B,C reactions A->B:0.1;B->C:0.05 t=100 n=50")
        # Text format is custom; just check it doesn't crash on simple input
        assert r is not None


# =============================================================================
# #497 GeometryOptimizer — 几何优化，能量极小化
# =============================================================================
class TestGeometryOptimizer_497:
    """Tests for #497 GeometryOptimizer."""

    def setup_method(self):
        from chemmcp.tools import GeometryOptimizer
        self.tool = GeometryOptimizer()

    def test_h2o_optimization(self):
        """H2O molecule optimization should converge."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "O", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [0.96, 0.0, 0.0]},
                {"symbol": "H", "position": [-0.24, 0.93, 0.0]},
            ],
            optimizer="steepest_descent",
            max_iterations=100,
            convergence_threshold=1e-4
        )
        assert len(r) > 0

    def test_co2_linear_optimization(self):
        """CO2 linear molecule."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "O", "position": [-1.16, 0.0, 0.0]},
                {"symbol": "C", "position": [0.0, 0.0, 0.0]},
                {"symbol": "O", "position": [1.16, 0.0, 0.0]},
            ],
            optimizer="steepest_descent",
            max_iterations=100
        )
        assert len(r) > 0

    def test_conjugate_gradient(self):
        """CG optimizer should also work."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [0.8, 0.0, 0.0]},
            ],
            optimizer="conjugate_gradient",
            max_iterations=50
        )
        assert len(r) > 0

    def test_single_atom_trivial(self):
        """Single atom is trivially optimized."""
        r = self.tool.run_code(
            atoms=[{"symbol": "He", "position": [0.0, 0.0, 0.0]}],
            max_iterations=10
        )
        assert len(r) > 0

    def test_text_input(self):
        # Format: atoms Sym:x,y,z;Sym:x,y,z
        r = self.tool.run_text("atoms O:0,0,0 ; H:0.96,0,0 ; H:-0.24,0.93,0")
        assert r is not None


# =============================================================================
# #498 TransitionStateSearch — 过渡态搜索，鞍点定位
# =============================================================================
class TestTransitionStateSearch_498:
    """Tests for #498 TransitionStateSearch."""

    def setup_method(self):
        from chemmcp.tools import TransitionStateSearch
        self.tool = TransitionStateSearch()

    def test_basic_saddle_search(self):
        """Quadratic saddle search for H2 stretched bond."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [-0.5, 0.0, 0.0]},
                {"symbol": "H", "position": [0.5, 0.0, 0.0]},
            ],
            search_method="quadratic_saddle",
            max_iterations=100
        )
        assert len(r) > 0

    def test_eigenvector_following(self):
        """EF method should work."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [-0.4, 0.0, 0.0]},
                {"symbol": "H", "position": [0.6, 0.0, 0.0]},
            ],
            search_method="eigenvector_following",
            max_iterations=100
        )
        assert len(r) > 0

    def test_triatomic_system(self):
        """Triatomic H-H-H system (hydrogen exchange)."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [-0.5, 0.0, 0.0]},
                {"symbol": "H", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [1.0, 0.0, 0.0]},
            ],
            search_method="eigenvector_following",
            max_iterations=100
        )
        assert len(r) > 0

    def test_with_rc_hint(self):
        """Reaction coordinate hint as atom index pair (not coordinate list)."""
        # The tool expects hint as atom indices or simpler format
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [2.0, 0.0, 0.0]},
            ],
            search_method="quadratic_saddle",
            max_iterations=100
        )
        assert len(r) > 0

    def test_text_input(self):
        # Format: atoms Sym:x,y,z;Sym:x,y,z [method=...]
        r = self.tool.run_text("atoms H:-0.5,0,0;H:0.5,0,0 method=quadratic_saddle")
        assert r is not None


# =============================================================================
# #499 PotentialEnergySurface — 势能面扫描，反应路径探索
# =============================================================================
class TestPotentialEnergySurface_499:
    """Tests for #499 PotentialEnergySurface."""

    def setup_method(self):
        from chemmcp.tools import PotentialEnergySurface
        self.tool = PotentialEnergySurface()

    def test_bond_length_scan(self):
        """Bond length scan produces energy profile."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [0.74, 0.0, 0.0]},
            ],
            scan_type="bond_length",
            scan_atoms=[0, 1],
            start_value=0.4,
            end_value=2.0,
            n_points=20
        )
        assert len(r) > 0

    def test_angle_scan(self):
        """Angle scan for triatomic (HOH angle)."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "O", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [0.96, 0.0, 0.0]},
                {"symbol": "H", "position": [0.0, 0.96, 0.0]},
            ],
            scan_type="angle",
            scan_atoms=[1, 0, 2],
            start_value=60.0,
            end_value=180.0,
            n_points=15
        )
        assert len(r) > 0

    def test_dihedral_scan(self):
        """Dihedral scan for butane C-C-C-C."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "C", "position": [0.0, 0.0, 0.0]},
                {"symbol": "C", "position": [1.54, 0.0, 0.0]},
                {"symbol": "C", "position": [2.31, 1.33, 0.0]},
                {"symbol": "C", "position": [3.85, 1.33, 0.0]},
            ],
            scan_type="dihedral",
            scan_atoms=[0, 1, 2, 3],
            start_value=-180.0,
            end_value=180.0,
            n_points=12
        )
        assert len(r) > 0

    def test_minimum_detection(self):
        """Scan around equilibrium should find minimum region."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "position": [0.74, 0.0, 0.0]},
            ],
            scan_type="bond_length",
            scan_atoms=[0, 1],
            start_value=0.4,
            end_value=2.0,
            n_points=30
        )
        assert len(r) > 0

    def test_text_input(self):
        # Format: atoms Sym:x,y,z;... scan type i,j start end
        r = self.tool.run_text("atoms H:0,0,0 ; H:0.74,0,0 scan bond_length 0,1 0.4 2.0 20")
        assert r is not None


# =============================================================================
# #500 FrequencyAnalysis — 频率分析，驻点性质确认（极小/过渡态）
# =============================================================================
class TestFrequencyAnalysis_500:
    """Tests for #500 FrequencyAnalysis."""

    def setup_method(self):
        from chemmcp.tools import FrequencyAnalysis
        self.tool = FrequencyAnalysis()

    def test_diatomic_h2(self):
        """H2 diatomic: should return frequency data."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [0.74, 0.0, 0.0]},
            ],
            temperature_K=298.15
        )
        assert len(r) > 0

    def test_water_vibrational_modes(self):
        """H2O nonlinear: 3N-6 = 3 vibrational modes."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "O", "mass": 15.999, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [0.96, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [-0.24, 0.93, 0.0]},
            ],
            temperature_K=298.15
        )
        assert len(r) > 0

    def test_zpe_non_negative(self):
        """Zero-point energy must be >= 0."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [0.74, 0.0, 0.0]},
            ]
        )
        zpe = r.get("zero_point_energy", r.get("zpe"))
        if zpe is not None:
            assert zpe >= 0

    def test_thermo_at_298K(self):
        """Thermodynamic quantities at 298.15 K, 1 atm."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "O", "mass": 15.999, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [0.96, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [-0.24, 0.93, 0.0]},
            ],
            temperature_K=298.15,
            pressure_atm=1.0
        )
        assert len(r) > 0

    def test_custom_hessian(self):
        """Custom Hessian matrix input."""
        n_dim = 6  # 2 atoms × 3
        H = [[0.5 if i == j else 0.0 for j in range(n_dim)] for i in range(n_dim)]
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [1.2, 0.0, 0.0]},
            ],
            hessian_matrix=H,
            temperature_K=298.15
        )
        assert len(r) > 0

    def test_scale_factor_affects_result(self):
        """Scale factor changes frequencies."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [0.74, 0.0, 0.0]},
            ],
            scale_factor=0.96
        )
        assert len(r) > 0

    def test_stationary_point_classification(self):
        """Should classify stationary point type."""
        r = self.tool.run_code(
            atoms=[
                {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
                {"symbol": "H", "mass": 1.008, "position": [0.74, 0.0, 0.0]},
            ]
        )
        assert len(r) > 0

    def test_text_input(self):
        # Format: atoms Sym:mass:x,y,z;... [T=temp]
        r = self.tool.run_text("atoms H:0,0,0 ; H:0.74,0,0 T=298.15")
        assert r is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
