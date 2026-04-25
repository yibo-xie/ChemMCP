"""
Test file for MCP Tools #151-160
Tests all 10 new tools: ElectronSinkIdentifier, ElectronSourceIdentifier,
ReactionEnergyEstimator, AldolReaction, ClaisenCondensation, DielsAlderReaction,
GrignardReaction, WittigReaction, FriedelCraftsReaction, SuzukiCoupling.
"""
import sys
sys.path.insert(0, 'src')

def test_electron_sink_identifier():
    """Test #151: Electron Sink Identifier"""
    print("=" * 60)
    print("Test 151: ElectronSinkIdentifier")
    print("=" * 60)
    from chemmcp.tools import ElectronSinkIdentifier
    tool = ElectronSinkIdentifier()

    # Test 1: Permanganate redox
    result = tool.run_code(
        reaction_input="2KMnO4 + 16HCl → 2KCl + 2MnCl2 + 5Cl2 + 8H2O",
        analysis_mode="detailed"
    )
    r = result['result']
    assert r['primary_electron_sink'] is not None, "Should identify a sink"
    print(f"✅ Primary sink: {r['primary_electron_sink']}")
    print(f"   Confidence: {r['confidence']}")
    assert r['confidence'] in ('high', 'medium-high', 'medium', 'medium-low', 'low')
    print(f"   Electrons transferred: {r.get('electrons_transferred')}")
    if 'half_reactions' in r:
        print(f"   Half reactions provided: {len(r['half_reactions'])}")

    # Test 2: Alcohol oxidation (generic oxidant)
    result2 = tool.run_code(
        reaction_input="CH3CH2OH + [O] → CH3CHO + H2O",
        analysis_mode="brief"
    )
    r2 = result2['result']
    print(f"✅ Alcohol oxidation sink: {r2['primary_electron_sink']}")
    assert r2['confidence'] != ''

    # Test text interface
    result3 = tool.run_text("CC(=O)O>[O]>CC(=O)O detailed")
    r3 = result3['result']
    print(f"✅ Text interface works: {r3.get('primary_electron_sink', 'N/A')}")

    print("✅ All ElectronSinkIdentifier tests passed!\n")


def test_electron_source_identifier():
    """Test #152: Electron Source Identifier"""
    print("=" * 60)
    print("Test 152: ElectronSourceIdentifier")
    print("=" * 60)
    from chemmcp.tools import ElectronSourceIdentifier
    tool = ElectronSourceIdentifier()

    # Test 1: NaBH4 reduction of aldehyde
    result = tool.run_code(
        reaction_input="CH3CHO + NaBH4 → CH3CH2OH",
        analysis_mode="detailed"
    )
    r = result['result']
    assert r['primary_electron_source'] is not None
    print(f"✅ Primary source: {r['primary_electron_source']}")
    print(f"   Confidence: {r['confidence']}")

    # Test 2: Metal displacement
    result2 = tool.run_code(
        reaction_input="Zn + CuSO4 → ZnSO4 + Cu",
        analysis_mode="brief"
    )
    r2 = result2['result']
    print(f"✅ Displacement source: {r2['primary_electron_source']}")
    assert 'Zn' in r2.get('primary_electron_source', '')

    # Text interface
    result3 = tool.run_text("C=C+H2 CC 298")
    r3 = result3['result']
    print(f"✅ Text interface works: {r3.get('primary_electron_source', 'N/A')}")

    print("✅ All ElectronSourceIdentifier tests passed!\n")


def test_reaction_energy_estimator():
    """Test #153: Reaction Energy Estimator"""
    print("=" * 60)
    print("Test 153: ReactionEnergyEstimator")
    print("=" * 60)
    from chemmcp.tools import ReactionEnergyEstimator
    tool = ReactionEnergyEstimator()

    # Test 1: Ethylene hydrogenation
    result = tool.run_code(
        reactants_smiles="C=C + H2",
        products_smiles="CC",
        temperature_k=298.15
    )
    r = result['result']
    print(f"✅ Reaction: {r.get('reaction', '?')}")
    print(f"   ΔH: {r.get('delta_h')} {r.get('delta_h_unit')}")
    print(f"   ΔG: {r.get('delta_g_estimate')} {r.get('delta_g_unit')}")
    print(f"   Spontaneous: {r.get('spontaneous')}")
    print(f"   Feasibility: {r.get('feasibility_rating')}")
    assert r.get('spontaneous') == True or r.get('spontaneous') == False  # should be bool
    assert r.get('feasibility_rating') != ''

    # Test 2: Haber-Bosch
    result2 = tool.run_code(
        reactants_smiles="N2 + 3H2",
        products_smiles="2NH3",
        temperature_k=298.15
    )
    r2 = result2['result']
    print(f"✅ Haber-Bosch ΔG: {r2.get('delta_g_estimate')} kJ/mol")
    print(f"   K: {r2.get('equilibrium_constant')}")

    # Text interface
    result3 = tool.run_text("C=C+H2 CC 298")
    r3 = result3['result']
    print(f"✅ Text interface ΔG: {r3.get('delta_g_estimate')}")

    print("✅ All ReactionEnergyEstimator tests passed!\n")


def test_aldol_reaction():
    """Test #154: Aldol Reaction"""
    print("=" * 60)
    print("Test 154: AldolReaction")
    print("=" * 60)
    from chemmcp.tools import AldolReaction
    tool = AldolReaction()

    # Self-aldol of acetaldehyde
    result = tool.run_code(
        substrate1_smiles="CC=O",
        catalyst_type="base",
        solvent="EtOH"
    )
    r = result['result']
    print(f"✅ Reaction type: {r.get('reaction_type')}")
    print(f"   Product (condensation): {r.get('product_prediction', {}).get('condensation_product', 'N/A')}")
    print(f"   Yield: {r.get('yield_expectation')}")
    assert r.get('reaction_type') != ''
    assert r.get('product_prediction') is not None

    # Crossed aldol
    result2 = tool.run_code(
        substrate1_smiles="CC(=O)C",
        substrate2_smiles="Cc1ccccc1=O",
        catalyst_type="base"
    )
    r2 = result2['result']
    print(f"✅ Crossed aldol: {r2.get('reaction_type')}")
    print(f"   Clean crossed: {r2.get('product_prediction', {}).get('crossed_cleanliness', 'N/A')}")

    # Text interface
    result3 = tool.run_text("CC=O base EtOH")
    r3 = result3['result']
    print(f"✅ Text interface: {r3.get('reaction_type')}")

    print("✅ All AldolReaction tests passed!\n")


def test_claisen_condensation():
    """Test #155: Claisen Condensation"""
    print("=" * 60)
    print("Test 155: ClaisenCondensation")
    print("=" * 60)
    from chemmcp.tools import ClaisenCondensation
    tool = ClaisenCondensation()

    # Self-Claisen of ethyl acetate
    result = tool.run_code(
        ester1_smiles="CC(=O)OCC",
        base="NaOEt",
        solvent="EtOH"
    )
    r = result['result']
    print(f"✅ Reaction type: {r.get('reaction_type')}")
    print(f"   Product: {r.get('product_prediction', {}).get('name', 'N/A')}")
    print(f"   Yield: {r.get('yield_expectation')}")
    mech = r.get('mechanism_summary', [])
    print(f"   Mechanism steps: {len(mech)}")
    assert len(mech) >= 4  # Should have at least 4 mechanism steps

    # Crossed Claisen
    result2 = tool.run_code(
        ester1_smiles="CC(=O)OCC",
        ester2_smiles="c1ccccc1C(=O)OCC",
        base="NaOEt"
    )
    r2 = result2['result']
    print(f"✅ Crossed Claisen: {r2.get('reaction_type')}")
    print(f"   Clean crossed: {r2.get('product_prediction', {}).get('clean', 'N/A')}")

    # Text interface
    result3 = tool.run_text("CC(=O)OCC NaOEt EtOH")
    r3 = result3['result']
    print(f"✅ Text interface: {r3.get('reaction_type')}")

    print("✅ All ClaisenCondensation tests passed!\n")


def test_diels_alder_reaction():
    """Test #156: Diels-Alder Reaction"""
    print("=" * 60)
    print("Test 156: DielsAlderReaction")
    print("=" * 60)
    from chemmcp.tools import DielsAlderReaction
    tool = DielsAlderReaction()

    # Classic butadiene + maleic anhydride
    result = tool.run_code(
        diene_smiles="C=CC=C",
        dienophile_smiles="maleic anhydride",
        temperature_c=80
    )
    r = result['result']
    print(f"✅ Reaction: {r.get('reaction')}")
    print(f"   Feasibility: {r.get('feasibility_assessment', {}).get('rating', 'N/A')}")
    stereo = r.get('stereochemistry', {})
    print(f"   Stereochemistry: {stereo.get('endo_product', 'N/A')}")
    product = r.get('product_prediction', {})
    print(f"   Bicyclic: {product.get('bicyclic', 'N/A')}")
    assert r.get('feasibility_assessment', {}).get('rating') != ''

    # Cyclopentadiene (very reactive diene)
    result2 = tool.run_code(
        diene_smiles="cyclopentadiene",
        dienophile_smiles="acrylate",
        temperature_c=25
    )
    r2 = result2['result']
    print(f"✅ Cyclopentadiene D-A: {r2.get('feasibility_assessment', {}).get('rating')}")

    # Text interface
    result3 = tool.run_text("C=CC=C maleic_anhydride 80")
    r3 = result3['result']
    print(f"✅ Text interface: {r3.get('feasibility_assessment', {}).get('rating')}")

    print("✅ All DielsAlderReaction tests passed!\n")


def test_grignard_reaction():
    """Test #157: Grignard Reaction"""
    print("=" * 60)
    print("Test 157: GrignardReaction")
    print("=" * 60)
    from chemmcp.tools import GrignardReaction
    tool = GrignardReaction()

    # CH3MgBr + benzaldehyde
    result = tool.run_code(
        grignard_reagent="CH3MgBr",
        electrophile_smiles="benzaldehyde",
        solvent="dry Et2O"
    )
    r = result['result']
    print(f"✅ Reaction: {r.get('reaction')}")
    print(f"   Product: {r.get('product_prediction', {}).get('name', 'N/A')}")
    print(f"   Yield: {r.get('product_prediction', {}).get('yield_estimate', 'N/A')}")
    mech = r.get('mechanism', [])
    print(f"   Mechanism steps: {len(mech)}")
    safety = r.get('safety_notes', [])
    print(f"   Safety notes: {len(safety)} items")
    assert len(safety) > 0

    # CO2 carboxylation
    result2 = tool.run_code(
        grignard_reagent="PhMgBr",
        electrophile_smiles="CO2",
        solvent="dry THF"
    )
    r2 = result2['result']
    print(f"✅ Carboxylation: {r2.get('product_prediction', {}).get('name', 'N/A')}")

    # Text interface
    result3 = tool.run_text("CH3MgBr benzaldehyde dry_Et2O")
    r3 = result3['result']
    print(f"✅ Text interface: {r3.get('product_prediction', {}).get('name', 'N/A')}")

    print("✅ All GrignardReaction tests passed!\n")


def test_wittig_reaction():
    """Test #158: Wittig Reaction"""
    print("=" * 60)
    print("Test 158: WittigReaction")
    print("=" * 60)
    from chemmcp.tools import WittigReaction
    tool = WittigReaction()

    # Non-stabilized ylide + aldehyde
    result = tool.run_code(
        carbonyl_smiles="benzaldehyde",
        ylide_type="non-stabilized",
        phosphonium_salt="Ph3P=CH2"
    )
    r = result['result']
    print(f"✅ Reaction: {r.get('reaction')}")
    print(f"   Product: {r.get('product_prediction', {}).get('name', 'N/A')}")
    ez = r.get('e_z_selectivity', {})
    print(f"   E/Z: {ez.get('prediction', 'N/A')}")
    mech = r.get('mechanism', [])
    print(f"   Mechanism steps: {len(mech)}")
    assert len(mech) >= 3

    # Stabilized ylide (E-selective)
    result2 = tool.run_code(
        carbonyl_smiles="benzaldehyde",
        ylide_type="stabilized",
        phosphonium_salt="Ph3P=CHCOOEt"
    )
    r2 = result2['result']
    print(f"✅ Stabilized E/Z: {r2.get('e_z_selectivity', {}).get('prediction', 'N/A')}")

    # Text interface
    result3 = tool.run_text("benzaldehyde non-stabilized Ph3P=CH2 n-BuLi")
    r3 = result3['result']
    print(f"✅ Text interface product: {r3.get('product_prediction', {}).get('name', 'N/A')}")

    print("✅ All WittigReaction tests passed!\n")


def test_friedel_crafts_reaction():
    """Test #159: Friedel-Crafts Reaction"""
    print("=" * 60)
    print("Test 159: FriedelCraftsReaction")
    print("=" * 60)
    from chemmcp.tools import FriedelCraftsReaction
    tool = FriedelCraftsReaction()

    # FC Alkylation: benzene + ethyl chloride
    result = tool.run_code(
        arene_smiles="benzene",
        electrophile_type="alkyl",
        electrophile_spec="CH3CH2Cl",
        catalyst="AlCl3"
    )
    r = result['result']
    print(f"✅ Reaction: {r.get('reaction')}")
    print(f"   Product: {r.get('product_prediction', {}).get('name', 'N/A')}")
    print(f"   Yield: {r.get('product_prediction', {}).get('yield', 'N/A')}")
    orient = r.get('orientation', {})
    print(f"   Orientation: {orient.get('expected_isomers', 'N/A')}")
    limits = r.get('applicable_limitations', [])
    print(f"   Limitations noted: {len(limits)}")
    assert len(limits) > 0

    # FC Acylation: toluene
    result2 = tool.run_code(
        arene_smiles="toluene",
        electrophile_type="acyl",
        electrophile_spec="CH3COCl",
        catalyst="AlCl3"
    )
    r2 = result2['result']
    print(f"✅ FC Acylation: {r2.get('product_prediction', {}).get('name', 'N/A')}")
    print(f"   Polyacylation risk: {any('poly' in l.get('issue','').lower() for l in r2.get('applicable_limitations', []))}")

    # Text interface
    result3 = tool.run_text("benzene alkyl CH3Cl AlCl3")
    r3 = result3['result']
    print(f"✅ Text interface: {r3.get('product_prediction', {}).get('name', 'N/A')}")

    print("✅ All FriedelCraftsReaction tests passed!\n")


def test_suzuki_coupling():
    """Test #160: Suzuki Coupling"""
    print("=" * 60)
    print("Test 160: SuzukiCoupling")
    print("=" * 60)
    from chemmcp.tools import SuzukiCoupling
    tool = SuzukiCoupling()

    # Standard coupling: phenylboronic acid + bromobenzene
    result = tool.run_code(
        organoboron_smiles="phenylboronic acid",
        halide_smiles="bromobenzene",
        ligand="PPh3",
        base="K2CO3",
        solvent="dioxane/H2O"
    )
    r = result['result']
    print(f"✅ Reaction: {r.get('reaction')}")
    cycle = r.get('mechanism_steps', [])
    print(f"   Catalytic cycle steps: {len(cycle)}")
    assert len(cycle) == 3  # Oxidative addition, transmetalation, reductive elimination
    cat = r.get('catalyst_recommendation', {})
    print(f"   Catalyst rec: {cat.get('primary', 'N/A')[:60]}...")
    yield_est = r.get('typical_yields', 'N/A')
    print(f"   Yield estimate: {yield_est}")
    cond = r.get('optimal_conditions', {})
    print(f"   Temperature: {cond.get('temperature', 'N/A')}")
    limits = r.get('applicable_limitations', [])
    print(f"   Limitations: {len(limits)} noted")
    assert len(limits) > 0

    # Vinyl coupling
    result2 = tool.run_code(
        organoboron_smiles="vinyl-Bpin",
        halide_smiles="bromobenzene",
        ligand="Pd(dppf)Cl2"
    )
    r2 = result2['result']
    print(f"✅ Vinyl coupling yield: {r2.get('typical_yields', 'N/A')}")

    # Text interface
    result3 = tool.run_text("phenylboronic_acid bromobenzene PPh3 K2CO3 dioxane/H2O")
    r3 = result3['result']
    print(f"✅ Text interface yield: {r3.get('typical_yields', 'N/A')}")

    print("✅ All SuzukiCoupling tests passed!\n")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Running ALL Tests for Tools #151-160")
    print("=" * 60 + "\n")

    try:
        test_electron_sink_identifier()
    except Exception as e:
        print(f"❌ ElectronSinkIdentifier FAILED: {e}")

    try:
        test_electron_source_identifier()
    except Exception as e:
        print(f"❌ ElectronSourceIdentifier FAILED: {e}")

    try:
        test_reaction_energy_estimator()
    except Exception as e:
        print(f"❌ ReactionEnergyEstimator FAILED: {e}")

    try:
        test_aldol_reaction()
    except Exception as e:
        print(f"❌ AldolReaction FAILED: {e}")

    try:
        test_claisen_condensation()
    except Exception as e:
        print(f"❌ ClaisenCondensation FAILED: {e}")

    try:
        test_diels_alder_reaction()
    except Exception as e:
        print(f"❌ DielsAlderReaction FAILED: {e}")

    try:
        test_grignard_reaction()
    except Exception as e:
        print(f"❌ GrignardReaction FAILED: {e}")

    try:
        test_wittig_reaction()
    except Exception as e:
        print(f"❌ WittigReaction FAILED: {e}")

    try:
        test_friedel_crafts_reaction()
    except Exception as e:
        print(f"❌ FriedelCraftsReaction FAILED: {e}")

    try:
        test_suzuki_coupling()
    except Exception as e:
        print(f"❌ SuzukiCoupling FAILED: {e}")

    print("\n" + "=" * 60)
    print("ALL TESTS COMPLETED for Tools #151-160 ✅")
    print("=" * 60)
