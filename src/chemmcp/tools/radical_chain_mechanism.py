"""
Radical Chain Mechanism (Tool #124)
自由基链式反应机理：引发、传递（链增长）、终止步骤。
Provides radical chain reaction mechanism analysis: initiation, propagation,
and termination steps with BDE analysis and selectivity prediction.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


# Bond dissociation energies (kJ/mol) — approximate
BDE_DATA = {
    'H3C-H': 439,  # primary C-H
    'H2CH2C-H': 423,  # secondary C-H
    'H3CH2C-H': 413,  # tertiary C-H
    'PhCH2-H': 372,   # benzylic C-H
    'H2C=CH-H': 465,  # vinylic C-H (strong!)
    'HC≡C-H': 556,    # acetylenic C-H (very strong)
    'Br-CH2CH3': 285, # C-Br (1°)
    'I-CH2CH3': 222,  # C-I (1°)
    'Cl-CH2CH3': 352, # C-Cl (1°)
    'F-CH3': 485,     # C-F (very strong)
    'HO-OH': 213,     # O-O in H2O2
    'CH3O-OCH3': 159, # O-O in peroxides
    'Br-Br': 193,     # Br-Br
    'Cl-Cl': 243,     # Cl-Cl
    'I-I': 151,       # I-I
    'H-Br': 366,      # H-Br
    'H-I': 298,       # H-I
    'H-Cl': 432,      # H-Cl
    'H-F': 567,      # H-F
}

RADICAL_INITIATORS = {
    'NBS/hv or AIBN': {'name': 'NBS/light or AIBN', 'type': 'radical initiator', 'initiation_temp': '~80°C or hv'},
    'ROOR': {'name': 'peroxide (ROOR)', 'type': 'thermal initiator', 'initiation_temp': '>80°C'},
    'hv': {'name': 'light (UV)', 'type': 'photoinitiation', 'initiation_temp': 'ambient with UV'},
    'AIBN': {'name': "AIBN (azobisisobutyronitrile)", 'type': 'thermal initiator', 'initiation_temp': '~70-80°C'},
    'Bz2O2': {'name': 'benzoyl peroxide', 'type': 'thermal initiator', 'initiation_temp': '~80°C'},
}


@ChemMCPManager.register_tool
class RadicalChainMechanism(BaseTool):
    __version__ = "0.1.0"
    name = "RadicalChainMechanism"
    func_name = 'explain_radical_chain_mechanism'
    description = "Explain radical chain reaction mechanisms: initiation (radical generation), propagation (chain-carrying steps), and termination (combination/disproportionation). Covers halogenation (allylic/benzylic selective), polymerization, and radical addition reactions."
    implementation_description = "Analyzes substrate for weak bonds (low BDE sites), identifies the most favorable hydrogen abstraction positions based on bond strength and radical stability, provides complete chain mechanism with all three stages, predicts regioselectivity of radical reactions, and evaluates the role of initiators and inhibitors."
    categories = ["Reaction"]
    tags = ["Radical", "Chain Reaction", "Initiation", "Propagation", "Termination", "BDE", "Selectivity"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES string of the substrate.'),
        ('reagent', 'str', 'NBS/hv or AIBN', 'Reagent/initiator. Options: NBS/hv or AIBN, ROOR, hv, AIBN, Bz2O2, Cl2/hv, Br2/hv, HBr/ROOR.'),
        ('reaction_type', 'str', 'auto', 'Reaction type: auto, bromination, chlorination, polymerization, addition.'),
        ('solvent', 'str', '', 'Optional solvent (e.g., CCl4, benzene).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles reagent [reaction_type] [solvent]. E.g., "CC(C)=CC NBS/hv_or_AIBN bromination".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing bde_analysis, initiation_steps, propagation_steps, termination_steps, selectivity_prediction, product_prediction, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'C=CC',
                'reagent': 'NBS/hv or AIBN',
                'reaction_type': 'bromination',
                'solvent': '',
            },
            'text_input': {'query': 'C=CC NBS/hv_or_AIBN bromination'},
            'output': {
                'result': {
                    'substrate': 'propene',
                    'reaction_type': 'Allylic Bromination (Wohl-Ziegler)',
                    'mechanism': 'radical chain',
                    'selectivity': 'allylic position favored (resonance-stabilized allylic radical, BDE ≈ 370 kJ/mol)',
                    'product': '3-bromopropene',
                    'favorability': 'excellent — NBS provides low, steady [Br•]',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'c1ccccc1C',
                'reagent': 'NBS/hv or AIBN',
                'reaction_type': 'bromination',
                'solvent': '',
            },
            'text_input': {'query': 'c1ccccc1C NBS/hv_or_AIBN bromination'},
            'output': {
                'result': {
                    'substrate': 'toluene',
                    'selectivity': 'benzylic position (BDE = 372 kJ/mol) >> other C-H bonds',
                    'product': 'benzyl bromide',
                    'radical_intermediate': 'benzyl radical (resonance-stabilized by phenyl ring)',
                    'favorability': 'excellent',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'CC=C',
                'reagent': 'HBr/ROOR',
                'reaction_type': 'addition',
                'solvent': '',
            },
            'text_input': {'query': 'CC=C HBr/ROOR addition'},
            'output': {
                'result': {
                    'substrate': 'propene',
                    'reaction_type': 'Anti-Markovnikov Hydrobromination (Kharasch effect)',
                    'mechanism': 'radical chain (peroxide effect)',
                    'product': '1-bromopropane (anti-Markovnikov)',
                    'key_step': 'Br• adds to less substituted carbon → more stable 2° radical → abstracts Br from HBr',
                    'favorability': 'good — requires peroxide or light initiation',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, reagent: str = 'NBS/hv or AIBN', reaction_type: str = 'auto', solvent: str = '') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Analyze substrate for weak bonds / abstraction sites
        substrate_analysis = self._analyze_substrate(mol)

        # 2. Identify reagent/initiator
        reagent_info = self._classify_reagent(reagent)

        # 3. Determine reaction type
        rxn_type = self._determine_reaction_type(reaction_type, reagent_info)

        # 4. Build mechanism steps
        mechanism = self._build_full_mechanism(substrate_analysis, reagent_info, rxn_type)

        # 5. Selectivity analysis
        selectivity = self._analyze_selectivity(substrate_analysis, rxn_type)

        # 6. Product prediction
        products = self._predict_products(substrate_analysis, rxn_type, selectivity)

        # 7. Favorability
        favorability = self._evaluate_favorability(substrate_analysis, reagent_info, rxn_type)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'substrate_analysis': substrate_analysis,
                'reagent_info': reagent_info,
                'reaction_type': rxn_type,
                'bde_analysis': substrate_analysis.get('bde_sites', []),
                **mechanism,
                'selectivity': selectivity,
                'product_prediction': products,
                'favorability': favorability,
                'summary': self._build_summary(rxn_type, products, selectivity, favorability),
            }
        }

        logger.info(f"Radical chain: {substrate_smiles} + {reagent} ({rxn_type}) → {favorability}")
        return result

    def _analyze_substrate(self, mol):
        """Analyze substrate for radical vulnerability."""
        sites = []

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue

            # Check for C-H bonds at this carbon
            n_h = atom.GetTotalNumHs()
            if n_h == 0:
                continue

            neighbors = atom.GetNeighbors()
            n_c_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6)

            # Determine C-H bond strength based on substitution
            has_double_bond_neighbor = any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()) and
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
                for n in neighbors if n.GetAtomicNum() == 6
            )
            is_allylic = has_double_bond_neighbor
            is_benzylic = any(n.GetIsAromatic() for n in neighbors)

            if n_c_neighbors >= 3:
                c_class = 'tertiary'; bde = 405  # ~tertiary C-H
            elif n_c_neighbors == 2:
                c_class = 'secondary'; bde = 420
            elif n_c_neighbors == 1:
                c_class = 'primary'; bde = 435
            else:
                c_class = 'methyl'; bde = 439

            if is_allylic:
                c_class += '/allylic'; bde = 365
            if is_benzylic:
                c_class += '/benzylic'; bde = 370

            # Check for vinylic/acetylenic C-H (very strong, unreactive toward radicals)
            is_vinylic = any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()) and
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() >= 2.0
                for n in neighbors if n.GetAtomicNum() == 6
            )
            if atom.GetTotalNumHs() > 0:
                for bond in mol.GetBonds():
                    if bond.GetBeginAtomIdx() == atom.GetIdx():
                        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                        if end_atom.GetAtomicNum() == 1:
                            if any(
                                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()) and
                                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() >= 2.0
                                for n in atom.GetNeighbors() if n.GetAtomicNum() == 6
                            ):
                                bde = 465; c_class = 'vinylic'

            radical_stability = self._classify_radical_stability(c_class)

            sites.append({
                'atom_idx': atom.GetIdx(),
                'carbon_class': c_class,
                'n_hydrogens': n_h,
                'estimated_bde_kjmol': bde,
                'radical_stability': radical_stability,
                'is_allylic': is_allylic,
                'is_benzylic': is_benzylic,
            })

        # Sort by BDE (lowest first = most vulnerable to H-abstraction)
        sites.sort(key=lambda s: s['estimated_bde_kjmol'])

        return {
            'abstraction_sites': sites,
            'most_vulnerable_site': sites[0] if sites else None,
            'bde_sites': sites[:5],  # top 5 weakest
        }

    def _classify_radical_stability(self, c_class):
        order = ['benzylic', 'allylic/tertiary', 'tertiary', 'allylic/secondary', 'secondary',
                 'allylic/primary', 'primary', 'methyl', 'vinylic']
        for i, pattern in enumerate(order):
            if pattern in c_class or c_class in pattern:
                stability_levels = ['exceptionally stable', 'very stable', 'stable', 'moderately stable',
                                   'somewhat stable', 'unstable', 'very unstable', 'extremely unstable', 'unreactive']
                return stability_levels[min(i, len(stability_levels)-1)]
        return 'unknown'

    def _classify_reagent(self, reagent):
        """Classify reagent/initiator."""
        r_clean = reagent.upper().replace(' ', '').replace('/', '')

        for key, info in RADICAL_INITIATORS.items():
            key_clean = key.upper().replace(' ', '').replace('/', '')
            if key_clean in r_clean or r_clean in key_clean:
                return {**info, 'key': key}

        # Special cases
        if 'HBR' in r_clean and 'ROOR' in r_clean:
            return {'name': 'HBr/peroxide (Kharasch)', 'type': 'anti-Markovnikov addition', 'initiation_temp': 'RT + peroxide'}
        if 'CL2' in r_clean:
            return {'name': 'Cl₂/light', 'type': 'chlorination', 'initiation_temp': 'hv or Δ'}
        if 'BR2' in r_clean:
            return {'name': 'Br₂/light', 'type': 'bromination', 'initiation_temp': 'hv or Δ'}

        return {'name': reagent, 'type': 'unknown', 'initiation_temp': '?'}

    def _determine_reaction_type(self, requested, reagent_info):
        if requested != 'auto':
            return requested
        rtype = reagent_info.get('type', '')
        if 'bromin' in rtype: return 'bromination'
        if 'chlorin' in rtype: return 'chlorination'
        if 'anti-markovnikov' in rtype or 'kharasch' in rtype.lower(): return 'addition'
        if 'polymer' in rtype: return 'polymerization'
        return 'bromination'  # default

    def _build_full_mechanism(self, substrate, reagent_info, rxn_type):
        """Build complete radical chain mechanism."""
        init_name = reagent_info.get('name', 'Initiator')

        # Initiation
        if rxn_type == 'bromination':
            if 'NBS' in init_name:
                initiation = [
                    {
                        'step': 'I-1',
                        'name': 'Initiator Decomposition / Br• Generation',
                        'equation': 'NBS → Br• (trace) OR ROOR → 2 RO•',
                        'details': (
                            f"Trace Br• from NBS (or RO• from peroxide/AIBN) serves as the "
                            f"chain-carrying radical initiator. Light (hv) or heat (~80°C) triggers this step."
                        ),
                    }
                ]
            else:
                initiation = [
                    {
                        'step': 'I-1',
                        'name': 'Radical Initiation',
                        'equation': 'Initiator → 2 R• (e.g., (PhCOO)₂ → 2 PhCOO• → Ph• + CO₂)',
                        'details': f"{init_name} homolytically cleaves to generate radicals.",
                    }
                ]

            propagation = [
                {
                    'step': 'P-1',
                    'name': 'Hydrogen Abstraction (Chain Transfer)',
                    'equation': f"R• + Substrate-H → R-H + Substrate•",
                    'details': (
                        f"The radical abstracts the most weakly bonded hydrogen (lowest BDE site), "
                        f"generating a substrate radical. This step determines regioselectivity."
                    ),
                },
                {
                    'step': 'P-2',
                    'name': 'Halogen Atom Transfer',
                    'equation': "Substrate• + Br₂ (or NBS) → Product-Br + Br•",
                    'details': (
                        "The substrate radical reacts with Br₂ (or NBS as Br source) to give "
                        "the brominated product and regenerate Br•, sustaining the chain."
                    ),
                },
            ]

        elif rxn_type == 'addition':
            initiation = [{
                'step': 'I-1', 'name': 'Peroxide Initiation',
                'equation': 'ROOR → 2 RO•',
                'details': 'Peroxide homolytically cleaves upon heating or UV to generate alkoxy radicals.',
            }]
            propagation = [
                {
                    'step': 'P-1', 'name': 'Radical Addition to Alkene',
                    'equation': 'RO• + HBr → ROH + Br•  (or Br• adds directly)',
                    'details': 'Br• radical adds to the less substituted end of the alkene.',
                },
                {
                    'step': 'P-2', 'name': 'Chain Transfer',
                    'equation': 'Carbon radical + HBr → Alkyl bromide + Br•',
                    'details': 'The carbon radical abstracts Br from HBr, giving anti-Markovnikov product and regenerating Br•.',
                },
            ]

        else:  # generic
            initiation = [{'step': 'I-1', 'name': 'Initiation', 'equation': 'Init → 2 R•', 'details': 'Homolytic cleavage.'}]
            propagation = [
                {'step': 'P-1', 'name': 'Propagation Step 1', 'equation': 'R• + S → Intermediate', 'details': 'Chain transfer.'},
                {'step': 'P-2', 'name': 'Propagation Step 2', 'equation': 'Intermediate + Reagent → P + R•', 'details': 'Chain regeneration.'},
            ]

        termination = [
            {
                'step': 'T-1',
                'name': 'Combination',
                'equation': 'R• + R• → R-R',
                'details': 'Two radicals combine to form a covalent bond.',
            },
            {
                'step': 'T-2',
                'name': 'Disproportionation',
                'equation': 'R• + R\'• → R-H + R\'(unsaturated)',
                'details': 'One radical abstracts H from another, giving alkane + alkene.',
            },
            {
                'step': 'T-3',
                'name': 'Other Termination Pathways',
                'equation': 'Various radical-radical coupling modes',
                'details': 'Any two radicals present can terminate the chain.',
            },
        ]

        return {
            'initiation_steps': initiation,
            'propagation_steps': propagation,
            'termination_steps': termination,
        }

    def _analyze_selectivity(self, substrate, rxn_type):
        """Analyze radical selectivity."""
        sites = substrate.get('abstraction_sites', [])
        if not sites:
            return {'no_abstraction_sites': True}

        best = sites[0]
        worst = sites[-1] if len(sites) > 1 else best

        bde_diff = worst['estimated_bde_kjmol'] - best['estimated_bde_kjmol']

        # At 25°C, each 10 kJ/mol BDE difference gives ~50× selectivity ratio
        # (from Arrhenius equation, assuming similar pre-exponential factors)
        relative_rate_best = 10 ** (bde_diff / (2.303 * 0.008314 * 298)) * 100 if bde_diff > 0 else 100
        relative_rate_worst = 100

        return {
            'most_favored_position': f"C({best['atom_idx']}) — {best['carbon_class']} C-H (BDE ≈ {best['estimated_bde_kjmol']} kJ/mol)",
            'least_favored_position': f"C({worst['atom_idx']}) — {worst['carbon_class']} C-H (BDE ≈ {worst['estimated_bde_kjmol']} kJ/mol)",
            'bde_range': f'{best["estimated_bde_kjmol"]} – {worst["estimated_bde_kjmol"]} kJ/mol',
            'selectivity_ratio_approximate': f'>{int(relative_rate_best)}:1 (best vs worst site)',
            'radical_stability_order': [f"{s['carbon_class']} (BDE={s['estimated_bde_kjmol']})" for s in sites[:5]],
        }

    def _predict_products(self, substrate, rxn_type, selectivity):
        """Predict major product(s)."""
        best = substrate.get('most_vulnerable_site')
        if not best:
            return {'error': 'No suitable abstraction sites found.'}

        if rxn_type == 'bromination':
            return {
                'major_product': f'brominated at {best["carbon_class"]} position (C{best["atom_idx"]})',
                'regiochemistry': f'Selective for lowest-BDE C-H bond ({best["radical_stability"]} radical)',
                'byproducts': 'dibromination (if excess NBS/Br₂), allylic rearrangement products possible',
            }
        elif rxn_type == 'addition':
            return {
                'major_product': 'anti-Markovnikov addition product',
                'regiochemistry': 'Radical adds to less substituted carbon → radical at more substituted position',
                'note': 'Only for HX where X is a good radical trap (Br, I; NOT Cl or F)',
            }
        return {'major_product': f'radical-substituted product at most reactive site'}

    def _evaluate_favorability(self, substrate, reagent_info, rxn_type):
        score = 3  # baseline
        best = substrate.get('most_vulnerable_site')
        if best:
            bde = best.get('estimated_bde_kjmol', 500)
            if bde < 380: score += 3  # benzylic/allylic
            elif bde < 420: score += 2  # tertiary
            elif bde < 440: score += 1  # secondary/primary

        rtype = reagent_info.get('type', '')
        if 'NBS' in rtype or 'AIBN' in rtype: score += 2
        if 'peroxide' in rtype.lower(): score += 1

        if score >= 7: return 'excellent'
        elif score >= 5: return 'good'
        elif score >= 3: return 'moderate'
        return 'possible but may need optimization'

    def _build_summary(self, rxn_type, products, selectivity, favorability):
        prod = products.get('major_product', '?')
        sel = selectivity.get('most_favored_position', '?')
        return f"Radical {rxn_type}. Major product: {prod}. Selectivity: {sel}. Favorability: {favorability}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        reag = parts[1] if len(parts) > 1 else 'NBS/hv or AIBN'
        rtype = parts[2] if len(parts) > 2 else 'auto'
        solv = parts[3] if len(parts) > 3 else ''
        return self._run_base(sub, reag, rtype, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
