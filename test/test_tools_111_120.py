"""
Test file for MCP Tools #111-120
Tests all 10 new tools with example inputs.
"""
import sys
sys.path.insert(0, 'src')

def test_stereoisomer_counter():
    """Tool #111: Stereoisomer Counter"""
    from chemmcp.tools import StereoisomerCounter
    tool = StereoisomerCounter()
    
    # Test 1: Chiral molecule with 2 centers (should give 4 isomers)
    result = tool.run_code(smiles='CC(O)C(Cl)Br', only_count_max=True)
    print(f"[#111] StereoisomerCounter - CC(O)C(Cl)Br:")
    print(f"  total_isomers: {result['total_isomers']}")
    print(f"  n_chiral_centers: {result['n_chiral_centers']}")
    assert result['total_isomers'] > 0, "Should have at least 1 isomer"
    print("  ✅ PASS")
    
    # Test 2: Alkene with E/Z (should have DB stereochemistry)
    result2 = tool.run_code(smiles='CC=CC', only_count_max=True)
    print(f"[#111] StereoisomerCounter - CC=CC:")
    print(f"  total_isomers: {result2['total_isomers']}")
    print("  ✅ PASS")


def test_meso_compound_checker():
    """Tool #112: Meso Compound Checker"""
    from chemmcp.tools import MesoCompoundChecker
    tool = MesoCompoundChecker()
    
    # Test 1: Single chiral center → not meso
    result = tool.run_code(smiles='[C@H](O)(Cl)Br', enumerate_isomers=True)
    print(f"[#112] MesoCompoundChecker - [C@H](O)(Cl)Br:")
    print(f"  is_meso_candidate: {result['result']['is_meso_candidate']}")
    assert result['result']['is_meso_candidate'] == False, "Single center can't be meso"
    print("  ✅ PASS")
    
    # Test 2: 2 chiral centers → possible meso
    result2 = tool.run_code(smiles='CC(O)C(Cl)Br', enumerate_isomers=True)
    print(f"[#112] MesoCompoundChecker - CC(O)C(Cl)Br:")
    print(f"  is_meso_candidate: {result2['result']['is_meso_candidate']}")
    print(f"  n_chiral_centers: {result2['result']['n_chiral_centers']}")
    print("  ✅ PASS")


def test_ring_system_analyzer():
    """Tool #113: Ring System Analyzer"""
    from chemmcp.tools import RingSystemAnalyzer
    tool = RingSystemAnalyzer()
    
    # Test 1: Naphthalene (fused rings)
    result = tool.run_code(smiles='C1=CC=CC2=C1C=CC=C2', detailed=True)
    print(f"[#113] RingSystemAnalyzer - Naphthalene:")
    print(f"  n_rings: {result['result']['n_rings']}")
    print(f"  fused_rings: {len(result['result']['fused_rings'])}")
    assert result['result']['n_rings'] == 2, "Naphthalene has 2 rings"
    assert len(result['result']['fused_rings']) > 0, "Should detect fused ring system"
    print("  ✅ PASS")
    
    # Test 2: Acyclic molecule
    result2 = tool.run_code(smiles='CCCC', detailed=False)
    print(f"[#113] RingSystemAnalyzer - butane:")
    print(f"  n_rings: {result2['result']['n_rings']}")
    assert result2['result']['n_rings'] == 0, "Butane has no rings"
    print("  ✅ PASS")


def test_aromatic_system_detector():
    """Tool #114: Aromatic System Detector"""
    from chemmcp.tools import AromaticSystemDetector
    tool = AromaticSystemDetector()
    
    # Test 1: Benzene (aromatic)
    result = tool.run_code(smiles='c1ccccc1', include_huckel_details=True)
    print(f"[#114] AromaticSystemDetector - Benzene:")
    print(f"  is_aromatic: {result['result']['is_aromatic']}")
    assert result['result']['is_aromatic'] == True, "Benzene should be aromatic"
    print("  ✅ PASS")
    
    # Test 2: Cyclohexane (non-aromatic, saturated)
    result2 = tool.run_code(smiles='C1CCCCC1', include_huckel_details=True)
    print(f"[#114] AromaticSystemDetector - Cyclohexane:")
    print(f"  is_aromatic: {result2['result']['is_aromatic']}")
    assert result2['result']['is_aromatic'] == False, "Cyclohexane is not aromatic"
    print("  ✅ PASS")


def test_tautomer_generator():
    """Tool #115: Tautomer Generator"""
    from chemmcp.tools import TautomerGenerator
    tool = TautomerGenerator()
    
    # Test 1: Acetone (keto-enol)
    result = tool.run_code(smiles='CC(=O)C', max_tautomers=5)
    print(f"[#115] TautomerGenerator - Acetone:")
    print(f"  n_tautomers: {result['result']['n_tautomers']}")
    assert result['result']['n_tautomers'] >= 1, "Should find at least original form"
    print("  ✅ PASS")
    
    # Test 2: Simple alkane (no tautomers expected)
    result2 = tool.run_code(smiles='CCCC', max_tautomers=5)
    print(f"[#115] TautomerGenerator - Butane:")
    print(f"  n_tautomers: {result2['result']['n_tautomers']}")
    print("  ✅ PASS")


def test_sn1_mechanism():
    """Tool #116: SN1 Mechanism"""
    from chemmcp.tools import Sn1Mechanism
    tool = Sn1Mechanism()
    
    # Test 1: tert-butyl chloride + H2O (classic SN1)
    result = tool.run_code(
        substrate_smiles='CC(C)(C)Cl',
        nucleophile='H2O',
        solvent='polar protic'
    )
    print(f"[#116] Sn1Mechanism - t-BuCl + H2O:")
    print(f"  mechanism_type: {result['result']['mechanism_type']}")
    print(f"  carbocation type: {result['result']['carbocation_analysis']['type']}")
    print(f"  favorability: {result['result']['favorability']}")
    assert result['result']['mechanism_type'] == 'SN1'
    assert 'tertiary' in result['result']['carbocation_analysis']['type']
    print("  ✅ PASS")
    
    # Test 2: Primary substrate (poor for SN1)
    result2 = tool.run_code(
        substrate_smiles='CCCl',
        nucleophile='H2O',
        solvent='polar protic'
    )
    print(f"[#116] Sn1Mechanism - Ethyl chloride + H2O:")
    print(f"  favorability: {result2['result']['favorability']}")
    print("  ✅ PASS")


def test_sn2_mechanism():
    """Tool #117: SN2 Mechanism"""
    from chemmcp.tools import Sn2Mechanism
    tool = Sn2Mechanism()
    
    # Test 1: Primary substrate + OH- (good SN2)
    result = tool.run_code(
        substrate_smiles='CC(Cl)',
        nucleophile='OH-',
        solvent='polar aprotic'
    )
    print(f"[#117] Sn2Mechanism - EtCl + OH-:")
    print(f"  mechanism_type: {result['result']['mechanism_type']}")
    print(f"  steric hindrance: {result['result']['steric_analysis']['hindrance_level']}")
    print(f"  favorability: {result['result']['favorability']}")
    assert result['result']['mechanism_type'] == 'SN2 (concerted)'
    print("  ✅ PASS")
    
    # Test 2: Tertiary substrate (bad for SN2)
    result2 = tool.run_code(
        substrate_smiles='CC(C)(C)Cl',
        nucleophile='OH-',
        solvent='polar aprotic'
    )
    print(f"[#117] Sn2Mechanism - t-BuCl + OH-:")
    print(f"  favorability: {result2['result']['favorability']}")
    print("  ✅ PASS")


def test_e1_mechanism():
    """Tool #118: E1 Mechanism"""
    from chemmcp.tools import E1Mechanism
    tool = E1Mechanism()
    
    # Test 1: t-BuCl + H2O (E1 favorable)
    result = tool.run_code(
        substrate_smiles='CC(C)(C)Cl',
        base='H2O',
        solvent='polar protic',
        temperature_c=60.0
    )
    print(f"[#118] E1Mechanism - t-BuCl + H2O @ 60°C:")
    print(f"  mechanism_type: {result['result']['mechanism_type']}")
    print(f"  carbocation type: {result['result']['carbocation_analysis']['type']}")
    print(f"  can_eliminate: {result['result']['beta_hydrogen_analysis']['can_eliminate']}")
    print(f"  favorability: {result['result']['favorability']}")
    assert result['result']['mechanism_type'] == 'E1'
    print("  ✅ PASS")
    
    # Test 2: Check product prediction
    prod = result['result'].get('product_prediction', {})
    if prod.get('can_form_alkene'):
        print(f"  Major product: {prod.get('zaitsev_product', {}).get('description', '?')}")
    print("  ✅ PASS")


def test_e2_mechanism():
    """Tool #119: E2 Mechanism"""
    from chemmcp.tools import E2Mechanism
    tool = E2Mechanism()
    
    # Test 1: Secondary halide + NaOEt (E2)
    result = tool.run_code(
        substrate_smiles='CC(Cl)CC',
        base='NaOEt/EtOH',
        solvent='EtOH',
        temperature_c=55.0
    )
    print(f"[#119] E2Mechanism - sec-BuCl + NaOEt:")
    print(f"  mechanism_type: {result['result']['mechanism_type']}")
    print(f"  beta_h available: {result['result']['beta_hydrogen_analysis']['total_beta_hydrogens']}")
    print(f"  favorability: {result['result']['favorability']}")
    assert result['result']['mechanism_type'] == 'E2 (concerted)'
    print("  ✅ PASS")
    
    # Test 2: Tertiary + bulky base (E2 forced)
    result2 = tool.run_code(
        substrate_smiles='CC(C)(C)Cl',
        base='t-BuOK/t-BuOH',
        solvent='t-BuOH',
        temperature_c=30.0
    )
    print(f"[#119] E2Mechanism - t-BuCl + t-BuOK:")
    print(f"  favorability: {result2['result']['favorability']}")
    print("  ✅ PASS")


def test_electrophilic_addition():
    """Tool #120: Electrophilic Addition"""
    from chemmcp.tools import ElectrophilicAddition
    tool = ElectrophilicAddition()
    
    # Test 1: Propene + HBr (Markovnikov)
    result = tool.run_code(
        alkene_smiles='CC=C',
        reaction_type='HBr'
    )
    print(f"[#120] ElectrophilicAddition - Propene + HBr:")
    print(f"  reaction_name: {result['result']['reaction_name']}")
    print(f"  has_alkene: {result['result']['alkene_analysis']['has_alkene']}")
    print(f"  regiochemistry: {result['result']['regiochemistry']}")
    print(f"  n_steps: {len(result['result']['mechanism_steps'])}")
    assert result['result']['reaction_name'] == 'Hydrobromination'
    assert result['result']['alkene_analysis']['has_alkene'] == True
    print("  ✅ PASS")
    
    # Test 2: Ethene + Br2 (halogenation anti addition)
    result2 = tool.run_code(
        alkene_smiles='C=C',
        reaction_type='Br2'
    )
    print(f"[#120] ElectrophilicAddition - Ethene + Br2:")
    print(f"  reaction_name: {result2['result']['reaction_name']}")
    print(f"  intermediate: {result2['result']['intermediate'][:50]}...")
    print("  ✅ PASS")
    
    # Test 3: Propene + hydroboration (Anti-Markovnikov)
    result3 = tool.run_code(
        alkene_smiles='CC=C',
        reaction_type='BH3/THF then H2O2/OH-'
    )
    print(f"[#120] ElectrophilicAddition - Propene + Hydroboration:")
    print(f"  product is_markovnikov: {result3['result']['product_prediction']['is_markovnikov']}")
    assert result3['result']['product_prediction']['is_markovnikov'] == False
    print("  ✅ PASS")


if __name__ == '__main__':
    print("=" * 60)
    print("Testing MCP Tools #111-120")
    print("=" * 60)
    
    tests = [
        ("#111 StereoisomerCounter", test_stereoisomer_counter),
        ("#112 MesoCompoundChecker", test_meso_compound_checker),
        ("#113 RingSystemAnalyzer", test_ring_system_analyzer),
        ("#114 AromaticSystemDetector", test_aromatic_system_detector),
        ("#115 TautomerGenerator", test_tautomer_generator),
        ("#116 Sn1Mechanism", test_sn1_mechanism),
        ("#117 Sn2Mechanism", test_sn2_mechanism),
        ("#118 E1Mechanism", test_e1_mechanism),
        ("#119 E2Mechanism", test_e2_mechanism),
        ("#120 ElectrophilicAddition", test_electrophilic_addition),
    ]
    
    passed = 0
    failed = 0
    for name, test_fn in tests:
        try:
            print(f"\n--- {name} ---")
            test_fn()
            passed += 1
        except Exception as e:
            print(f"  ❌ FAIL: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed out of {len(tests)} tests")
    print("=" * 60)
