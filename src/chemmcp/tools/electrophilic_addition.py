"""
Electrophilic Addition Reaction Mechanism (Tool #120)
烯烃亲电加成反应机理（如卤化、水合、氢卤化等）。
Provides electrophilic addition mechanism analysis for alkenes: halogenation, hydration, hydrohalogenation, etc.
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


# Common electrophilic addition reaction types
ADDITION_REACTIONS = {
    'HBr': {
        'name': 'Hydrobromination',
        'reagent': 'HBr',
        'electrophile': 'H⁺',
        'nucleophile': 'Br⁻',
        'regiochemistry': "Markovnikov's rule (H adds to less substituted C)",
        'stereochemistry': 'anti addition (via bromonium ion or carbocation)',
        'category': 'hydrohalogenation',
    },
    'HCl': {
        'name': 'Hydrochlorination',
        'reagent': 'HCl',
        'electrophile': 'H⁺',
        'nucleophile': 'Cl⁻',
        'regiochemistry': "Markovnikov's rule",
        'stereochemistry': 'anti or mixed depending on intermediate',
        'category': 'hydrohalogenation',
    },
    'HI': {
        'name': 'Hydroiodination',
        'reagent': 'HI',
        'electrophile': 'H⁺',
        'nucleophile': 'I⁻',
        'regiochemistry': "Markovnikov's rule",
        'stereochemistry': 'anti addition',
        'category': 'hydrohalogenation',
    },
    'H2O/H2SO4': {
        'name': 'Acid-Catalyzed Hydration',
        'reagent': 'H₂O / H₂SO₄',
        'electrophile': 'H⁺ (from acid)',
        'nucleophile': 'H₂O',
        'regiochemistry': "Markovnikov's rule → alcohol",
        'stereochemistry': 'mixed (carbocation intermediate allows rotation)',
        'category': 'hydration',
    },
    'Br2': {
        'name': 'Bromination',
        'reagent': 'Br₂',
        'electrophile': 'Br⁺ (polarizable Br₂)',
        'nucleophile': 'Br⁻',
        'regiochemistry': 'symmetric (no regioselectivity issue)',
        'stereochemistry': 'ANTI addition via cyclic bromonium ion',
        'category': 'halogenation',
    },
    'Cl2': {
        'name': 'Chlorination',
        'reagent': 'Cl₂',
        'electrophile': 'Cl⁺',
        'nucleophile': 'Cl⁻',
        'regiochemistry': 'symmetric',
        'stereochemistry': 'ANTI addition via chloronium ion',
        'category': 'halogenation',
    },
    'BH3/THF then H2O2/OH-': {
        'name': 'Hydroboration-Oxidation',
        'reagent': 'BH₃ then H₂O₂/OH⁻',
        'electrophile': 'BH₃ (electrophilic B)',
        'nucleophile': 'H⁻ (from boron)',
        'regiochemistry': "Anti-Markovnikov (OH adds to less substituted C)",
        'stereochemistry': 'SYN addition',
        'category': 'hydroboration',
    },
    'Hg(OAc)2/H2O then NaBH4': {
        'name': 'Oxymercuration-Demercuration',
        'reagent': 'Hg(OAc)₂/H₂O then NaBH₄',
        'electrophile': 'Hg(OAc)⁺',
        'nucleophile': 'H₂O',
        'regiochemistry': "Markovnikov (without rearrangement)",
        'stereochemistry': 'mixed (no rearrangement)',
        'category': 'oxymercuration',
    },
    'OsO4/NMO': {
        'name': 'Dihydroxylation (Syn)',
        'reagent': 'OsO₄ / NMO',
        'electrophile': 'OsO₄ (cyclic osmate ester)',
        'nucleophile': '-O-Os bond',
        'regiochemistry': 'symmetric',
        'stereochemistry': 'SYN addition (both OH on same face)',
        'category': 'dihydroxylation',
    },
    'KMnO4/cold_dilute': {
        'name': 'Dihydroxylation (KMnO₄)',
        'reagent': 'cold dilute KMnO₄',
        'electrophile': 'MnO₄⁻',
        'nucleophile': '-O-Mn bond',
        'regiochemistry': 'symmetric',
        'stereochemistry': 'SYN addition',
        'category': 'dihydroxylation',
    },
    'H2/Pt': {
        'name': 'Catalytic Hydrogenation',
        'reagent': 'H₂ / Pt, Pd, or Ni',
        'electrophile': 'H₂ (syn adsorbed on surface)',
        'nucleophile': 'H• (surface)',
        'regiochemistry': 'symmetric',
        'stereochemistry': 'SYN addition (cis alkane product)',
        'category': 'hydrogenation',
    },
}


@ChemMCPManager.register_tool
class ElectrophilicAddition(BaseTool):
    __version__ = "0.1.0"
    name = "ElectrophilicAddition"
    func_name = 'explain_electrophilic_addition'
    description = "Explain electrophilic addition reactions to alkenes: halogenation, hydrohalogenation, hydration, hydroboration, dihydroxylation, hydrogenation, and more. Includes mechanism steps, regiochemistry (Markovnikov/anti-Markovnikov), and stereochemistry."
    implementation_description = "Analyzes the alkene substrate structure, identifies substitution pattern and stereochemistry of the double bond, matches with the specified reagent type, and provides a complete stepwise mechanism including intermediate structures, regiochemical outcome (Markovnikov vs anti-Markovnikov), and stereochemical result (syn vs anti)."
    categories = ["Reaction"]
    tags = ["Electrophilic Addition", "Alkene", "Markovnikov", "Halogenation", "Hydration", "Mechanism"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('alkene_smiles', 'str', 'N/A', 'SMILES string of the alkene substrate.'),
        ('reaction_type', 'str', 'HBr', 'Type of electrophilic addition. Options: HBr, HCl, HI, H2O/H2SO4, Br2, Cl2, BH3/THF then H2O2/OH-, Hg(OAc)2/H2O then NaBH4, OsO4/NMO, KMnO4/cold_dilute, H2/Pt.'),
        ('solvent', 'str', '', 'Optional solvent (e.g., CCl4 for halogenation, THF for hydroboration).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: alkene_smiles reaction_type [solvent]. E.g., "CC=C HBr".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing mechanism_steps, regiochemistry, stereochemistry, predicted_product, intermediate_analysis, and Markovnikov explanation.'),
    ]
    examples = [
        {
            'code_input': {
                'alkene_smiles': 'CC=C',
                'reaction_type': 'HBr',
                'solvent': '',
            },
            'text_input': {'query': 'CC=C HBr'},
            'output': {
                'result': {
                    'alkene': 'propene (terminal, monosubstituted)',
                    'reaction': 'Hydrobromination',
                    'mechanism_type': 'electrophilic addition via carbocation',
                    'steps': [
                        {'step': 1, 'name': 'Electrophilic attack (protonation)', 'description': 'π electrons attack H⁺ → more stable carbocation forms'},
                        {'step': 2, 'name': 'Nucleophilic attack', 'description': 'Br⁻ attacks carbocation → alkyl bromide'},
                    ],
                    'regiochemistry': 'Markovnikov: H adds to CH₂ (less sub.), Br adds to CH (more sub.) → 2-bromopropane',
                    'stereochemistry': 'No new stereocenter in this case',
                    'predicted_product': '2-bromopropane (major)',
                }
            },
        },
        {
            'code_input': {
                'alkene_smiles': 'C=C',
                'reaction_type': 'Br2',
                'solvent': 'CCl4',
            },
            'text_input': {'query': 'C=C Br2'},
            'output': {
                'result': {
                    'reaction': 'Bromination',
                    'intermediate': 'cyclic bromonium ion (3-membered ring)',
                    'stereochemistry': 'ANTI addition → meso or racemic mixture depending on substituents',
                    'predicted_product': '1,2-dibromoethane',
                }
            },
        },
        {
            'code_input': {
                'alkene_smiles': 'CC=C',
                'reaction_type': 'BH3/THF then H2O2/OH-',
                'solvent': 'THF',
            },
            'text_input': {'query': 'CC=C BH3/THF then H2O2/OH-'},
            'output': {
                'result': {
                    'reaction': 'Hydroboration-Oxidation',
                    'regiochemistry': 'Anti-Markovnikov: OH adds to terminal carbon → 1-propanol',
                    'stereochemistry': 'SYN addition',
                }
            },
        },
    ]

    def _run_base(self, alkene_smiles: str, reaction_type: str = 'HBr', solvent: str = '') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(alkene_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(alkene_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Analyze alkene
        alkene_analysis = self._analyze_alkene(mol)

        # Get reaction info
        rxn_info = ADDITION_REACTIONS.get(reaction_type)
        if rxn_info is None:
            # Try fuzzy match
            rxn_info = self._fuzzy_match_reaction(reaction_type)
            if rxn_info is None:
                available = ', '.join(ADDITION_REACTIONS.keys())
                raise ChemMCPInputError(
                    f"Unknown reaction type '{reaction_type}'. Available types: {available}"
                )

        # Build mechanism steps
        steps = self._build_mechanism_steps(alkene_analysis, rxn_info)

        # Predict product
        product = self._predict_product(alkene_analysis, rxn_info)

        # Stereochemistry analysis
        stereo = self._analyze_stereochemistry(alkene_analysis, rxn_info)

        result = {
            'result': {
                'alkene_smiles': alkene_smiles,
                'alkene_analysis': alkene_analysis,
                'reaction_type': reaction_type,
                'reaction_name': rxn_info['name'],
                'reaction_category': rxn_info['category'],
                'reagent': rxn_info['reagent'],
                'mechanism_steps': steps,
                'intermediate': self._describe_intermediate(rxn_info),
                'regiochemistry': rxn_info['regiochemistry'],
                'regiochemistry_explanation': self._explain_regiochemistry(alkene_analysis, rxn_info),
                'stereochemistry': stereo,
                'product_prediction': product,
                'markovnikov_analysis': self._markovnikov_detail(alkene_analysis, rxn_info),
                'summary': self._build_summary(alkene_analysis, rxn_info, product),
            }
        }

        logger.info(f"Electrophilic addition: {alkene_smiles} + {reaction_type} → {rxn_info['name']}")
        return result

    def _analyze_alkene(self, mol):
        """Analyze alkene structure."""
        # Find double bond(s)
        double_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                atom1 = mol.GetAtomWithIdx(a1)
                atom2 = mol.GetAtomWithIdx(a2)

                # Count substituents on each carbon
                subs1 = sum(1 for n in atom1.GetNeighbors() if n.GetIdx() != a2 and n.GetAtomicNum() > 1)
                subs2 = sum(1 for n in atom2.GetNeighbors() if n.GetIdx() != a1 and n.GetAtomicNum() > 1)

                # Check existing stereochemistry
                stereo = bond.GetStereo()
                stereo_str = {Chem.BondStereo.STEREONONE: 'none', Chem.BondStereo.STEREOZ: 'Z', Chem.BondStereo.STEREOE: 'E'}.get(stereo, str(stereo))

                double_bonds.append({
                    'atom1_idx': a1, 'atom1_symbol': atom1.GetSymbol(), 'atom1_subs': subs1,
                    'atom2_idx': a2, 'atom2_symbol': atom2.GetSymbol(), 'atom2_subs': subs2,
                    'substitution': f"{subs1 + 1},{subs2 + 1}",  # total substituents per carbon
                    'configured_stereochemistry': stereo_str,
                })

        if not double_bonds:
            return {'has_alkene': False, 'n_double_bonds': 0}

        # Classify overall alkene
        db = double_bonds[0]
        total_sub = db['atom1_subs'] + db['atom2_subs']
        if total_sub == 0: alkene_class = 'unsubstituted (ethene-type)'
        elif total_sub <= 1: alkene_class = 'monosubstituted (terminal)'
        elif total_sub == 2: alkene_class = 'disubstituted'
        elif total_sub == 3: alkene_class = 'trisubstituted'
        else: alkene_class = 'tetrasubstituted'

        is_terminal = db['atom1_subs'] == 0 or db['atom2_subs'] == 0
        is_symmetric = db['atom1_subs'] == db['atom2_subs']

        return {
            'has_alkene': True,
            'n_double_bonds': len(double_bonds),
            'double_bond_info': db,
            'alkene_class': alkene_class,
            'is_terminal_alkene': is_terminal,
            'is_symmetric': is_symmetric,
            'more_substituted_side': 'C2' if db['atom2_subs'] > db['atom1_subs'] else 'C1' if db['atom1_subs'] > db['atom2_subs'] else 'equal',
        }

    def _fuzzy_match_reaction(self, rt):
        """Fuzzy match reaction type."""
        rt_lower = rt.lower().replace(' ', '')
        for key, val in ADDITION_REACTIONS.items():
            key_clean = key.lower().replace(' ', '').replace('/', '')
            if rt_lower in key_clean or key_clean in rt_lower:
                return val
        return None

    def _build_mechanism_steps(self, alkene, rxn):
        category = rxn.get('category', '')
        steps = []

        if category == 'hydrohalogenation':
            steps = [
                {
                    'step': 1, 'name': 'Electrophilic Attack (Protonation)',
                    'equation': f"Alkene + {rxn['electrophile']} → Carbocation",
                    'details': (
                        f"π electrons of the alkene attack the electrophilic {rxn['electrophile']}. "
                        f"The proton adds to the less-substituted carbon (or following Markovnikov rule), "
                        f"generating the most stable carbocation on the more-substituted carbon."
                    ),
                    'rate_relevant': True,
                },
                {
                    'step': 2, 'name': 'Nucleophilic Attack',
                    'equation': f"Carbocation + {rxn['nucleophile']} → Product",
                    'details': (
                        f"The nucleophile ({rxn['nucleophile']}) attacks the planar carbocation "
                        f"to form the addition product."
                    ),
                    'rate_relevant': False,
                },
            ]
        elif category == 'halogenation':
            steps = [
                {
                    'step': 1, 'name': 'Electrophilic Attack & Cyclic Halonium Formation',
                    'equation': f"Alkene + X₂ → Cyclic halonium ion + X⁻",
                    'details': (
                        f"The polarizable {rxn['reagent']} molecule approaches the alkene π cloud. "
                        f"One X atom acts as electrophile, forming a 3-membered cyclic halonium ion. "
                        f"This prevents free rotation and ensures anti stereochemistry."
                    ),
                    'rate_relevant': True,
                },
                {
                    'step': 2, 'name': 'Halide Ion Attack (Backside)',
                    'equation': f"Cyclic halonium ion + X⁻ → Vicinal dihalide",
                    'details': (
                        f"The halide ion attacks from the backside of the halonium ring (anti addition), "
                        f"opening the 3-membered ring to form the vicinal dihalide."
                    ),
                    'rate_relevant': False,
                },
            ]
        elif category == 'hydration':
            steps = [
                {
                    'step': 1, 'name': 'Protonation',
                    'equation': 'Alkene + H⁺ → Carbocation',
                    'details': 'Proton adds to less-substituted carbon (Markovnikov).',
                    'rate_relevant': True,
                },
                {
                    'step': 2, 'name': 'Nucleophilic Attack by Water',
                    'equation': 'Carbocation + H₂O → Oxonium ion',
                    'details': 'Water attacks carbocation.',
                    'rate_relevant': False,
                },
                {
                    'step': 3, 'name': 'Deprotonation',
                    'equation': 'Oxonium ion → Alcohol',
                    'details': 'Base removes a proton to give the neutral alcohol.',
                    'rate_relevant': False,
                },
            ]
        elif category == 'hydroboration':
            steps = [
                {
                    'step': 1, 'name': 'Concerted Syn Hydroboration',
                    'equation': 'Alkene + BH₃ → Alkylborane',
                    'details': (
                        'BH₃ adds in a concerted, syn fashion across the double bond. '
                        'Boron (electrophilic) adds to the LESS substituted carbon (steric control), '
                        'H adds to the MORE substituted carbon.'
                    ),
                    'rate_relevant': True,
                },
                {
                    'step': 2, 'name': 'Oxidation',
                    'equation': 'Alkylborane + H₂O₂/OH⁻ → Alcohol',
                    'details': 'Oxidation replaces B with OH with retention of configuration.',
                    'rate_relevant': False,
                },
            ]
        elif category == 'dihydroxylation':
            steps = [
                {
                    'step': 1, 'name': 'Cyclic Osmate/Ester Formation',
                    'equation': f'Alkene + {rxn["reagent"].split("/")[0]} → Cyclic ester',
                    'details': 'Concerted syn addition forms a cyclic osmate or manganate ester.',
                    'rate_relevant': True,
                },
                {
                    'step': 2, 'name': 'Hydrolysis',
                    'equation': 'Cyclic ester → cis-diol',
                    'details': 'Hydrolysis releases the vicinal diol with syn stereochemistry.',
                    'rate_relevant': False,
                },
            ]
        elif category == 'hydrogenation':
            steps = [
                {
                    'step': 1, 'name': 'Syn Addition on Metal Surface',
                    'equation': 'Alkene + H₂ (on catalyst) → Alkane',
                    'details': (
                        'Both H atoms add simultaneously from the same face of the alkene '
                        '(syn addition) via adsorption on the metal catalyst surface (Pt, Pd, Ni).'
                    ),
                    'rate_relevant': True,
                },
            ]
        else:
            steps = [{
                'step': 1, 'name': 'Electrophilic Addition',
                'equation': f'Alkene + {rxn["reagent"]} → Product',
                'details': f'Standard electrophilic addition with {rxn["reagent"]}.',
                'rate_relevant': True,
            }]

        return steps

    def _describe_intermediate(self, rxn):
        category = rxn.get('category', '')
        intermediates = {
            'hydrohalogenation': 'Carbocation (planar sp² hybridized). More stable carbocation determines regiochemistry.',
            'halogenation': 'Cyclic halonium ion (3-membered ring). Prevents rotation → anti addition.',
            'hydration': 'Carbocation intermediate (may undergo rearrangement).',
            'hydroboration': 'No discrete intermediate — concerted 4-center transition state.',
            'oxymercuration': 'Mercurinium ion (3-membered ring, similar to halonium). No rearrangement!',
            'dihydroxylation': 'Cyclic osmate/manganate ester (5-membered ring).',
            'hydrogenation': 'Surface-adsorbed species on metal catalyst (no isolated intermediate).',
        }
        return intermediates.get(category, 'Reaction-specific intermediate.')

    def _predict_product(self, alkene, rxn):
        if not alkene.get('has_alkene'):
            return {'error': 'No alkene found in substrate.'}

        db = alkene.get('double_bond_info', {})
        cat = rxn.get('category', '')
        more_sub = alkene.get('more_substituted_side', '')

        products = {
            'hydrohalogenation': f"{rxn['reagent']} adds across double bond → {rxn['nucleophile']} on more-substituted side ({more_sub}), H on less-substituted. Product: alkyl {rxn['reagent'].split('/')[0] if '/' in rxn['reagent'] else rxn['reagent'][:-1]}ide.",
            'halogenation': f"Vicinal di{rxn['reagent'][0].lower()}ide (anti addition).",
            'hydration': 'Alcohol (Markovnikov: OH on more-substituted carbon).',
            'hydroboration': 'Alcohol (Anti-Markovnikov: OH on less-substituted carbon).',
            'oxymercuration': 'Alcohol (Markovnikov, no rearrangement).',
            'dihydroxylation': 'Vicinal cis-diol (syn addition).',
            'hydrogenation': 'Alkane (syn addition, saturated).',
        }

        return {
            'product_type': cat,
            'description': products.get(cat, 'Addition product formed.'),
            'is_markovnikov': cat not in ('hydroboration',),
        }

    def _analyze_stereochemistry(self, alkene, rxn):
        cat = rxn.get('category', '')
        has_configured = alkene.get('double_bond_info', {}).get('configured_stereochemistry', 'none')
        is_terminal = alkene.get('is_terminal_alkene', False)
        is_symmetric = alkene.get('is_symmetric', False)

        stereo_map = {
            'hydrohalogenation': (
                'Non-stereospecific (carbocation intermediate allows rotation). '
                'If a new stereocenter is formed, racemic mixture results. '
                'No defined syn/anti relationship.'
            ),
            'halogenization': (
                'ANTI addition via cyclic halonium ion. '
                f'Two new stereocenters are formed with opposite absolute configuration. '
                f'{"Racemic pair (enantiomers)" if not is_symmetric else "Meso product"} expected.'
            ),
            'hydration': 'Not stereospecific (carbocation intermediate).',
            'hydroboration': 'SYN addition (concerted process). Both new groups add to same face.',
            'oxymercuration': 'Not stereospecific (mercurinium ion can open from either side).',
            'dihydroxylation': 'SYN addition → cis-diol (both OH on same face).',
            'hydrogenation': 'SYN addition → alkane (both H add from same face of catalyst).',
        }

        base = stereo_map.get(cat, 'Standard addition stereochemistry.')

        if is_terminal and cat in ('hydrohalogenation', 'hydration', 'hydroboration'):
            base += ' Terminal alkene → no new stereocenter at terminal carbon.'

        if has_configured != 'none':
            base += f' Starting alkene geometry: {has_configured}.'

        return base

    def _explain_regiochemistry(self, alkene, rxn):
        cat = rxn.get('category', '')
        more_sub = alkene.get('more_substituted_side', '')

        if cat in ('halogenation', 'dihydroxylation', 'hydrogenation'):
            return f"Symmetric reagent or symmetric addition → no regioselectivity issue."

        if cat == 'hydroboration':
            return (
                "Anti-Markovnikov: Steric factors direct boron to the LESS hindered (less substituted) "
                "carbon. Electronic factors also contribute (partial positive on B seeks electron density)."
            )

        # Default: Markovnikov
        return (
            f"Markovnikov's Rule: The electrophile (H⁺ or equivalent) adds to the less-substituted carbon "
            f"(which has more H atoms), generating the more stable carbocation on the more-substituted side "
            f"({more_sub}). The nucleophile then adds to this carbocation."
        )

    def _markovnikov_detail(self, alkene, rxn):
        cat = rxn.get('category', '')
        db = alkene.get('double_bond_info', {})

        detail = {
            'rule_applied': 'Markovnikov' if cat != 'hydroboration' else 'Anti-Markovnikov',
            'rationale': (
                "The electrophile adds so as to produce the most stable intermediate carbocation. "
                f"Tertiary > secondary > primary > methyl stability order applies."
            ) if cat != 'hydroboration' else (
                "Steric control directs the bulky BH₃ group to the less-hindered terminal carbon. "
                "This gives anti-Markovnikov selectivity after oxidation."
            ),
            'possible_rearrangement': 'Yes — carbocation may rearrange (hydride/alkyl shift)' if cat in ('hydrohalogenation', 'hydration') else 'No',
            'asymmetry_note': (
                f"Asymmetric alkene (C1: {db.get('atom1_subs', 0)} subs, C2: {db.get('atom2_subs', 0)} subs) "
                f"→ regioselective outcome expected."
            ) if not alkene.get('is_symmetric') else (
                "Symmetric alkene → no regioselectivity issue (both sides equivalent)."
            ),
        }
        return detail

    def _build_summary(self, alkene, rxn, product):
        return (
            f"{rxn['name']} of {alkene.get('alkene_class', '?')} alkene. "
            f"{rxn.get('regiochemistry', '')}. "
            f"Product: {product.get('description', '?')}."
        )

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 2)
        alkene = parts[0] if len(parts) > 0 else ''
        rxn_type = parts[1] if len(parts) > 1 else 'HBr'
        solv = parts[2] if len(parts) > 2 else ''
        return self._run_base(alkene, rxn_type, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
