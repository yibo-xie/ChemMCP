import importlib

_tool_module_map = {
    "WebSearch":    "web_search",
    "MoleculeCaptioner": "molecule_captioner",
    "MoleculeGenerator": "molecule_generator",
    "PubchemSearchQA": "pubchem_search_qa",
    "PubchemSearch": "pubchem_search",  # Not registerred as an MCP tool
    "ForwardSynthesis": "forward_synthesis",
    "Retrosynthesis": "retrosynthesis",
    "MoleculeSimilarity": "molecule_similarity",
    "MoleculeWeight": "molecule_weight",
    "FunctionalGroups": "functional_groups",
    "SmilesCanonicalization": "smiles_canonicalization",
    "MoleculeAtomCount": "molecule_atom_count",
    "MoleculePrice": "molecule_price",
    "PatentCheck": "patent_check",
    "Iupac2Smiles": "iupac2smiles",
    "Smiles2Iupac": "smiles2iupac",
    "Smiles2Formula": "smiles2formula",
    "Name2Smiles": "name2smiles",
    "Selfies2Smiles": "selfies2smiles",
    "Smiles2Selfies": "smiles2selfies",
    "MoleculeSmilesCheck": "molecule_smiles_check",
    "ReactionSmilesCheck": "reaction_smiles_check",
    "Smiles2Cas": "smiles2cas",
    "SafetyCheck": "safety_check",
    "MoleculeModifier": "molecule_modifier",
    "MoleculeVisualizer": "molecule_visualizer",
    "KcKpConverter": "kc_kp_converter",

    # --- New chemistry MCP tools (Element Properties & Bonding) ---
    "GetElementInfo": "get_element_info",
    "GetElectronConfiguration": "get_electron_configuration",
    "GetOxidationStates": "get_oxidation_states",
    "GetIonizationEnergy": "get_ionization_energy",
    "GetElectronAffinity": "get_electron_affinity",
    "CompareElements": "compare_elements",
    "GetIsotopes": "get_isotopes",
    "PeriodicTrend": "periodic_trend",
    "GetElementDiscovery": "get_element_discovery",
    "ElementAbundance": "element_abundance",
    "DrawLewisStructure": "draw_lewis_structure",
    "PredictVseprGeometry": "predict_vsepr_geometry",
    "CalculateFormalCharge": "calculate_formal_charge",
    "GetBondLength": "get_bond_length",
    "GetBondEnergy": "get_bond_energy",
    "PredictHybridization": "predict_hybridization",
    "AnalyzeMolecularOrbital": "analyze_molecular_orbital",
    "GetCrystalStructure": "get_crystal_structure",
    "CalculateLatticeEnergy": "calculate_lattice_energy",
    "PredictPolarity": "predict_polarity",

    # --- Electrochemistry & Reaction Tools (#21-30) ---
    "SymmetryPointGroup": "symmetry_point_group",
    "CoordinationGeometry": "coordination_geometry",
    "BalanceEquation": "balance_equation",
    "BalanceRedox": "balance_redox",
    "IdentifyReactionType": "identify_reaction_type",
    "PredictProducts": "predict_products",
    "GetStandardPotential": "get_standard_potential",
    "CalculateCellPotential": "calculate_cell_potential",
    "NernstEquation": "nernst_equation",
    "IdentifyOxidizingAgent": "identify_oxidizing_agent",

    # --- Thermochemistry & Reaction Analysis Tools (#31-40) ---
    "AssignOxidationNumber": "assign_oxidation_number",
    "DisproportionationCheck": "disproportionation_check",
    "PrecipitationPrediction": "precipitation_prediction",
    "GasEvolutionPrediction": "gas_evolution_prediction",
    "GetStandardEnthalpy": "get_standard_enthalpy",
    "GetStandardEntropy": "get_standard_entropy",
    "GetGibbsEnergy": "get_gibbs_energy",
    "CalculateReactionEnthalpy": "calculate_reaction_enthalpy",
    "CalculateGibbsChange": "calculate_gibbs_change",
    "SpontaneityCheck": "spontaneity_check",

    # --- Chemical Equilibrium Tools (#41-50) ---
    "EquilibriumConstantThermo": "equilibrium_constant_thermo",
    "BornHaberCycle": "born_haber_cycle",
    "BondEnergyCalculation": "bond_energy_calculation",
    "TemperatureEffectK": "temperature_effect_k",
    "CalculateEquilibriumConstant": "calculate_equilibrium_constant",
    "ICETableSolver": "ice_table_solver",
    "LeChatelierPrediction": "le_chatelier_prediction",
    "PressureEffectEquilibrium": "pressure_effect_equilibrium",
    # Note: KcKpConverter already exists as #49 (convert_Kc_Kp)
    "ReactionQuotient": "reaction_quotient",

    # --- Acid-Base Equilibrium & Titration Tools (#51-60) ---
    "DegreeOfDissociation": "degree_of_dissociation",
    "CommonIonEffect": "common_ion_effect",
    "GetPka": "get_pka",
    "GetPkb": "get_pkb",
    "CalculatePH": "calculate_ph",
    "BufferPreparation": "buffer_preparation",
    "BufferCapacity": "buffer_capacity",
    "TitrationCurve": "titration_curve",
    "EquivalencePoint": "equivalence_point",
    "PolyproticAcid": "polyprotic_acid",

    # --- Solubility & Precipitation Analysis Tools (#61-70) ---
    "AmphotericSpecies": "amphoteric_species",
    "AcidBaseStrengthCompare": "acid_base_strength_compare",
    "GetKsp": "get_ksp",
    "CalculateSolubility": "calculate_solubility",
    "WillPrecipitate": "will_precipitate",
    "SelectivePrecipitation": "selective_precipitation",
    "DissolvePrecipitate": "dissolve_precipitate",
    "CommonIonSolubility": "common_ion_solubility",
    "ComplexIonSolubility": "complex_ion_solubility",
    "SolubilityRules": "solubility_rules",

    # --- Coordination Chemistry Tools (#71-80) ---
    "GetFormationConstant": "get_formation_constant",
    "NameComplex": "name_complex",
    "ComplexFromName": "complex_from_name",
    "CrystalFieldSplitting": "crystal_field_splitting",
    "PredictMagnetism": "predict_magnetism",
    "DElectronConfiguration": "d_electron_configuration",
    "SpectrochemicalSeries": "spectrochemical_series",
    "PredictColor": "predict_color",
    "IsomerTypes": "isomer_types",
    "ChelateEffect": "chelate_effect",

    # --- Inorganic Chemistry: Advanced Topics (#81-90) ---
    "JahnTellerDistortion": "jahn_teller_distortion",
    "LigandFieldDiagram": "ligand_field_diagram",
    "AlkaliMetalProperties": "alkali_metal_properties",
    "AlkalineEarthProperties": "alkaline_earth_properties",
    "HalogenProperties": "halogen_properties",
    "NobleGasCompounds": "noble_gas_compounds",
    "TransitionMetalChemistry": "transition_metal_chemistry",
    "LanthanideProperties": "lanthanide_properties",
    "ActinideProperties": "actinide_properties",
    "MainGroupTrends": "main_group_trends",

    # --- Nuclear Chemistry & General Tools (#91-100) ---
    "RadioactiveDecay": "radioactive_decay",
    "HalfLifeCalculation": "half_life_calculation",
    "NuclearEquationBalance": "nuclear_equation_balance",
    "BindingEnergy": "binding_energy",
    "MassDefect": "mass_defect",
    "DecaySeries": "decay_series",
    "GetPhysicalConstant": "get_physical_constant",
    "UnitConversion": "unit_conversion",
    "SignificantFigures": "significant_figures",
    "DimensionalAnalysis": "dimensional_analysis",

    # --- Stereochemistry & Advanced Molecular Analysis Tools (#101-110) ---
    "IupacNamer": "iupac_namer",
    "CommonNameLookup": "common_name_lookup",
    "SmilesToStructure": "smiles_to_structure",
    "StructureToSmiles": "structure_to_smiles",
    "FunctionalGroupIdentifier": "functional_group_identifier",
    "StereocenterFinder": "stereocenter_finder",
    "RsConfigurator": "rs_configurator",
    "EzConfigurator": "ez_configurator",
    "CisTransAnalyzer": "cis_trans_analyzer",
    "ConstitutionalIsomerGenerator": "constitutional_isomer_generator",

    # --- Advanced Stereochemistry & Aromaticity Tools (#111-115) ---
    "StereoisomerCounter": "stereoisomer_counter",
    "MesoCompoundChecker": "meso_compound_checker",
    "RingSystemAnalyzer": "ring_system_analyzer",
    "AromaticSystemDetector": "aromatic_system_detector",
    "TautomerGenerator": "tautomer_generator",

    # --- Reaction Mechanism Tools (#116-120) ---
    "Sn1Mechanism": "sn1_mechanism",
    "Sn2Mechanism": "sn2_mechanism",
    "E1Mechanism": "e1_mechanism",
    "E2Mechanism": "e2_mechanism",
    "ElectrophilicAddition": "electrophilic_addition",

    # --- Reaction Mechanism Tools (#121-130) ---
    "NucleophilicAddition": "nucleophilic_addition",
    "ElectrophilicAromaticSubstitution": "electrophilic_aromatic_substitution",
    "NucleophilicAromaticSubstitution": "nucleophilic_aromatic_substitution",
    "RadicalChainMechanism": "radical_chain_mechanism",
    "CarbocationRearrangement": "carbocation_rearrangement",
    "AldolMechanism": "aldol_mechanism",
    "ClaisenMechanism": "claisen_mechanism",
    "MichaelMechanism": "michael_mechanism",
    "DielsAlderMechanism": "diels_alder_mechanism",
    "GrignardMechanism": "grignard_mechanism",

    # --- Advanced Reaction Mechanism & Prediction Tools (#131-140) ---
    "WittigMechanism": "wittig_mechanism",
    "OxidationMechanism": "oxidation_mechanism",
    "ReductionMechanism": "reduction_mechanism",
    "PericyclicAnalyzer": "pericyclic_analyzer",
    "ArrowPushingValidator": "arrow_pushing_validator",
    "ReactionPredictor": "reaction_predictor",
    "RegioselectivityPredictor": "regioselectivity_predictor",
    "StereoselectivityPredictor": "stereoselectivity_predictor",
    "LeavingGroupRanker": "leaving_group_ranker",
    "NucleophilicityRanker": "nucleophilicity_ranker",

    # --- Advanced Analysis & Prediction Tools (#141-150) ---
    "BasicityVsNucleophilicity": "basicity_vs_nucleophilicity",
    "SolventSelector": "solvent_selector",
    "TemperatureAdvisor": "temperature_advisor",
    "CatalystRecommender": "catalyst_recommender",
    "ProtectingGroupSelector": "protecting_group_selector",
    "DeprotectionConditions": "deprotection_conditions",
    "CompetitionAnalyzer": "competition_analyzer",
    "KineticVsThermodynamic": "kinetic_vs_thermodynamic",
    "AcidBasePredictor": "acid_base_predictor",
    "OxidationStateCalculator": "oxidation_state_calculator",

    # --- Electron Transfer & Reaction Energy Tools (#151-153) ---
    "ElectronSinkIdentifier": "electron_sink_identifier",
    "ElectronSourceIdentifier": "electron_source_identifier",
    "ReactionEnergyEstimator": "reaction_energy_estimator",

    # --- Named Organic Reaction Tools (#154-160) ---
    "AldolReaction": "aldol_reaction",
    "ClaisenCondensation": "claisen_condensation",
    "DielsAlderReaction": "diels_alder_reaction",
    "GrignardReaction": "grignard_reaction",
    "WittigReaction": "wittig_reaction",
    "FriedelCraftsReaction": "friedel_crafts_reaction",
    "SuzukiCoupling": "suzuki_coupling",

    # --- Named Organic Reaction Tools (#161-170) ---
    "HeckReaction": "heck_reaction",
    "SonogashiraCoupling": "sonogashira_coupling",
    "BuchwaldHartwigAmination": "buchwald_hartwig_amination",
    "SwernOxidation": "swern_oxidation",
    "JonesOxidation": "jones_oxidation",
    "SharplessEpoxidation": "sharpless_epoxidation",
    "SharplessDihydroxylation": "sharpless_dihydroxylation",
    "BirchReduction": "birch_reduction",
    "WolffKishnerReduction": "wolff_kishner_reduction",
    "ClemmensenReduction": "clemmensen_reduction",

    # --- Spectroscopy & Spectral Analysis Tools (#171-180) ---
    "BeckmannRearrangement": "beckmann_rearrangement",
    "BaeyerVilligerOxidation": "baeyer_villiger_oxidation",
    "NamedReactionLookup": "named_reaction_lookup",
    "IrPeakInterpreter": "ir_peak_interpreter",
    "NmrHPredictor": "nmr_h_predictor",
    "NmrCPredictor": "nmr_c_predictor",
    "CouplingConstantAnalyzer": "coupling_constant_analyzer",
    "SplittingPatternExplainer": "splitting_pattern_explainer",
    "MassSpecFragmenter": "mass_spec_fragmenter",
    "MolecularIonCalculator": "molecular_ion_calculator",

    # --- Advanced Spectroscopy & Structure Elucidation Tools (#181-185) ---
    "IsotopePatternGenerator": "isotope_pattern_generator",
    "UvVisPredictor": "uv_vis_predictor",
    "SpectrumToStructure": "spectrum_to_structure",
    "DeptInterpreter": "dept_interpreter",
    "CosyNoesyGuide": "cosy_noesy_guide",

    # --- Retrosynthesis & Synthetic Analysis Tools (#186-190) ---
    "RetrosynthesisAnalyzer": "retrosynthesis_analyzer",
    "SynthonIdentifier": "synthon_identifier",
    "DisconnectionSuggester": "disconnection_suggester",
    "FunctionalGroupInterconversion": "functional_group_interconversion",
    "CarbonChainBuilder": "carbon_chain_builder",

    # --- Advanced Organic Chemistry & Physical Organic Tools (#191-200) ---
    "RingFormationStrategy": "ring_formation_strategy",
    "AsymmetricSynthesisGuide": "asymmetric_synthesis_guide",
    "TotalSynthesisPlanner": "total_synthesis_planner",
    "PkaPredictor": "pka_predictor",
    "HammettSigmaLookup": "hammett_sigma_lookup",
    "ResonanceStructureGenerator": "resonance_structure_generator",
    "InductiveEffectAnalyzer": "inductive_effect_analyzer",
    "HyperconjugationExplainer": "hyperconjugation_explainer",
    "StericEffectAnalyzer": "steric_effect_analyzer",
    "ConformationalAnalyzer": "conformational_analyzer",

    # --- Thermodynamics Tools (#201-210) ---
    "GibbsFreeEnergy": "gibbs_free_energy",
    "EnthalpyCalculator": "enthalpy_calculator",
    "EntropyCalculator": "entropy_calculator",
    "HeatCapacityLookup": "heat_capacity_lookup",
    "ClausiusClapeyron": "clausius_clapeyron",
    "VantHoffEquilibrium": "vant_hoff_equilibrium",
    "FugacityCalculator": "fugacity_calculator",
    "ActivityCoefficient": "activity_coefficient",
    "MaxwellRelations": "maxwell_relations",
    "CarnotEfficiency": "carnot_efficiency",

    # --- Advanced Thermodynamics Tools (#211-220) ---
    "JouleThomson": "joule_thomson",
    "ChemicalPotential": "chemical_potential",
    "PartialMolarQuantity": "partial_molar_quantity",
    "PhaseRuleAnalyzer": "phase_rule_analyzer",
    "StandardStateConverter": "standard_state_converter",
    "RateLawFitter": "rate_law_fitter",
    "ArrheniusAnalyzer": "arrhenius_analyzer",
    "HalfLifeCalculator": "half_life_calculator",
    "IntegratedRateLaw": "integrated_rate_law",
    "TransitionStateTheory": "transition_state_theory",

    # --- Chemical Kinetics & Reaction Mechanism Tools (#221-227) ---
    "CollisionTheory": "collision_theory",
    "EnzymeKinetics": "enzyme_kinetics",
    "ReactionMechanismSimulator": "reaction_mechanism_simulator",
    "SteadyStateApproximation": "steady_state_approximation",
    "RateDeterminingStep": "rate_determining_step",
    "TemperatureJumpRelaxation": "temperature_jump_relaxation",
    "ParallelConsecutiveReactions": "parallel_consecutive_reactions",

    # --- Quantum Mechanics Tools (#228-230) ---
    "ParticleInBox": "particle_in_box",
    "HarmonicOscillator": "harmonic_oscillator",
    "RigidRotor": "rigid_rotor",

    # --- Quantum Chemistry Tools (#231-240) ---
    "HydrogenAtomOrbitals": "hydrogen_atom_orbitals",
    "SchrodingerSolver1d": "schrodinger_solver_1d",
    "VariationalMethod": "variational_method",
    "PerturbationTheory": "perturbation_theory",
    "MolecularOrbitalDiagram": "molecular_orbital_diagram",
    "HuckelMethod": "huckel_method",
    "ElectronDensityPlotter": "electron_density_plotter",
    "SpinOrbitCoupling": "spin_orbit_coupling",
    "SelectionRulesChecker": "selection_rules_checker",
    "TunnelingProbability": "tunneling_probability",

    # --- Quantum Mechanics & Spectroscopy Tools (#241-250) ---
    "UncertaintyPrinciple": "uncertainty_principle",
    "ExpectationValue": "expectation_value",
    "IrSpectrumPredictor": "ir_spectrum_predictor",
    "RamanActivity": "raman_activity",
    "UvVisTransitions": "uv_vis_transitions",
    "NmrChemicalShift": "nmr_chemical_shift",
    "RotationalSpectrum": "rotational_spectrum",
    "VibrationalModes": "vibrational_modes",
    "FranckCondonFactors": "franck_condon_factors",
    "BeerLambertCalculator": "beer_lambert_calculator",

    # --- Statistical Thermodynamics & Spectroscopy Tools (#251-260) ---
    "FluorescenceLifetime": "fluorescence_lifetime",
    "StarkEffect": "stark_effect",
    "ZeemanSplitting": "zeeman_splitting",
    "SpectralLinewidth": "spectral_linewidth",
    "PartitionFunction": "partition_function",
    "BoltzmannDistribution": "boltzmann_distribution",
    "MaxwellBoltzmannSpeed": "maxwell_boltzmann_speed",
    "FermiDiracDistribution": "fermi_dirac_distribution",
    "BoseEinsteinDistribution": "bose_einstein_distribution",
    "StatisticalEntropy": "statistical_entropy",

    # --- Physical Chemistry: Statistical Mechanics & Electrochemistry Tools (#261-270) ---
    "EquipartitionTheorem": "equipartition_theorem",
    "DebyeModel": "debye_model",
    "EinsteinModel": "einstein_model",
    "EnsembleAverages": "ensemble_averages",
    # NernstEquation already exists (#265)
    "CellEmfCalculator": "cell_emf_calculator",
    "FaradayElectrolysis": "faraday_electrolysis",
    "DebyeHuckelActivity": "debye_huckel_activity",
    "ConductivityCalculator": "conductivity_calculator",
    "TafelEquation": "tafel_equation",

    # --- Electrochemistry & Surface Chemistry Tools (#271-280) ---
    "OverpotentialAnalyzer": "overpotential_analyzer",
    "PourbaixDiagramLookup": "pourbaix_diagram_lookup",
    "IonTransportNumber": "ion_transport_number",
    "ButlerVolmerKinetics": "butler_volmer_kinetics",
    "LangmuirIsotherm": "langmuir_isotherm",
    "BETSurfaceArea": "bet_surface_area",
    "FreundlichIsotherm": "freundlich_isotherm",
    "SurfaceTensionCalculator": "surface_tension_calculator",
    "ContactAngleAnalyzer": "contact_angle_analyzer",
    "GibbsAdsorption": "gibbs_adsorption",

    # --- Colloidal & Phase Equilibrium Tools (#281-290) ---
    "ColloidalStability": "colloidal_stability",
    "ZetaPotential": "zeta_potential",
    "BinaryPhaseDiagram": "binary_phase_diagram",
    "RaoultLaw": "raoult_law",
    "HenryLaw": "henrys_law",
    "BoilingPointElevation": "boiling_point_elevation",
    "FreezingPointDepression": "freezing_point_depression",
    "OsmoticPressure": "osmotic_pressure",
    "EutecticPointFinder": "eutectic_point_finder",
    "LeverRuleCalculator": "lever_rule_calculator",

    # --- New MCP Tools (#291-300) ---
    "PointGroupIdentifier": "point_group_identifier",
    "SymmetryOperations": "symmetry_operations",
    "BondOrderCalculator": "bond_order_calculator",
    "DipoleMomentEstimator": "dipole_moment_estimator",
    "HybridizationAnalyzer": "hybridization_analyzer",
    "VseprGeometry": "vsepr_geometry",
    "IdealGasCalculator": "ideal_gas_calculator",
    "VanDerWaalsGas": "van_der_waals_gas",
    "CompressibilityFactor": "compressibility_factor",
    "VirialEquation": "virial_equation",

    # --- Sample Preparation & Analytical Chemistry Tools (#301-310) ---
    "SampleDilutionCalculator": "sample_dilution_calculator",
    "StandardSolutionPrep": "standard_solution_prep",
    "ExtractionOptimizer": "extraction_optimizer",
    "DigestionProtocolSelector": "digestion_protocol_selector",
    "FiltrationGuide": "filtration_guide",
    "SPEMethodDesigner": "spe_method_designer",
    "DerivatizationReagentSelector": "derivatization_reagent_selector",
    "MatrixMatchingAdvisor": "matrix_matching_advisor",
    "SamplePreservationGuide": "sample_preservation_guide",
    "HomogenizationProtocol": "homogenization_protocol",

    # --- Analytical Chemistry Tools (#311-320) ---
    "CentrifugationCalculator": "centrifugation_calculator",
    "EvaporationEstimator": "evaporation_estimator",
    "PhAdjustmentBuffer": "ph_adjustment_buffer",
    "InternalStandardSelector": "internal_standard_selector",
    "RecoveryCalculator": "recovery_calculator",
    "UvVisWavelengthSelector": "uv_vis_wavelength_selector",
    "MolarAbsorptivityLookup": "molar_absorptivity_lookup",
    "IrSpectrumInterpreter": "ir_spectrum_interpreter",
    "FtirBaselineCorrector": "ftir_baseline_corrector",

    # --- Spectroscopy & Analytical Chemistry Tools (#321-330) ---
    "RamanShiftCalculator": "raman_shift_calculator",
    "FluorescenceQuantumYield": "fluorescence_quantum_yield",
    "ExcitationEmissionOptimizer": "excitation_emission_optimizer",
    "AasFlameSelector": "aas_flame_selector",
    "AasInterferenceChecker": "aas_interference_checker",
    "IcpOesLineSelector": "icp_oes_line_selector",
    "IcpMsIsotopeSelector": "icp_ms_isotope_selector",
    "PlasmaConditionOptimizer": "plasma_condition_optimizer",
    "XrfMatrixCorrection": "xrf_matrix_correction",
    "XrdPhaseIdentifier": "xrd_phase_identifier",

    # --- Spectroscopy & Chromatography Tools (#331-340) ---
    "NmrChemicalShiftPredictor": "nmr_chemical_shift_predictor",
    # Note: CouplingConstantAnalyzer already registered as #177
    "SpectralDeconvolution": "spectral_deconvolution",
    "HplcColumnSelector": "hplc_column_selector",
    "MobilePhaseOptimizer": "mobile_phase_optimizer",
    "RetentionTimePredictor": "retention_time_predictor",
    "PlateNumberCalculator": "plate_number_calculator",
    "ResolutionCalculator": "resolution_calculator",
    "GcOvenProgramDesigner": "gc_oven_program_designer",
    "GcCarrierGasSelector": "gc_carrier_gas_selector",

    # --- Chromatography & Analytical Chemistry Tools (#341-350) ---
    "GcColumnBleedPredictor": "gc_column_bleed_predictor",
    "IonChromatographyEluent": "ion_chromatography_eluent",
    "SecCalibrationCurve": "sec_calibration_curve",
    "PeakPurityAnalyzer": "peak_purity_analyzer",
    "SystemSuitabilityChecker": "system_suitability_checker",
    "DeadVolumeCalculator": "dead_volume_calculator",
    "VanDeemterAnalyzer": "van_deemter_analyzer",
    "CapacityFactorCalculator": "capacity_factor_calculator",
    "SelectivityFactorCalculator": "selectivity_factor_calculator",
    "MolecularIonCalculator": "molecular_ion_calculator",

    # --- Mass Spectrometry Tools (#351-360) ---
    "IsotopePatternSimulator": "isotope_pattern_simulator",
    "FragmentationPredictor": "fragmentation_predictor",
    "AdductIonIdentifier": "adduct_ion_identifier",
    "MRmTransitionOptimizer": "mrm_transition_optimizer",
    "MassAccuracyCalculator": "mass_accuracy_calculator",
    "CollisionEnergyOptimizer": "collision_energy_optimizer",
    "IonSuppressionChecker": "ion_suppression_checker",
    "ElementalCompositionCalculator": "elemental_composition_calculator",
    "MsMsSpectrumAnnotator": "msms_spectrum_annotator",
    "MatrixClusterIdentifier": "matrix_cluster_identifier",
}

__all__ = list(_tool_module_map.keys())

def __getattr__(name: str):
    if name in _tool_module_map:
        module_name = _tool_module_map.get(name)
        if module_name is None:
            raise AttributeError(f"No mapping for tool {name!r} in chemmcp.tools")
        module = importlib.import_module(f"{__name__}.{module_name}")
        try:
            return getattr(module, name)
        except AttributeError:
            raise ImportError(f"Module {module_name!r} has no attribute {name!r}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return list(__all__)
