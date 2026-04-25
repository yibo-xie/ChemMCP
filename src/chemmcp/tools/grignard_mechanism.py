"""
Grignard Mechanism (Tool #130)
格氏试剂反应机理：Grignard 试剂形成、对羰基的亲核加成、
四面体中间体、酸性后处理、反应活性顺序（甲醛->1°醇、醛->2°醇、酮->3°醇）、
限制条件（无活泼 H）。
Provides Grignard reagent reaction mechanism analysis: formation, nucleophilic addition
to carbonyls (aldehydes, ketones, esters, CO2, epoxides), tetrahedral intermediate,
acidic workup, reactivity order, and functional group compatibility.
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


# Carbonyl reactivity toward Grignard reagents (relative rates)
CARBONYL_REACTIVITY = {
    'formaldehyde':     {'rank': 1, 'product_type': '1° alcohol (primary)',   'description': 'H2C=O -> RCH2OH', 'reactivity': 'very high'},
    'aldehyde':          {'rank': 2, 'product_type': '2° alcohol (secondary)', 'description': 'R-CHO -> R-RCHOH (secondary)', 'reactivity': 'high'},
    'ketone':            {'rank': 3, 'product_type': '3° alcohol (tertiary)', 'description': 'R-CO-R -> R-R-COH-R (tertiary)', 'reactivity': 'moderate-high'},
    'ester':             {'rank': 4, 'product_type': '3° alcohol + 1° alcohol', 'description': 'R-COO-R -> R(R)-COH + R-OH (2 equiv needed)', 'reactivity': 'moderate (2 equiv RMgX needed)'},
    'acid_chloride':     {'rank': 5, 'product_type': '3° alcohol after workup', 'description': 'RCOCl -> tertiary alcohol', 'reactivity': 'high but often over-addition'},
    'CO2':               {'rank': 6, 'product_type': 'carboxylic acid',        'description': 'CO2 -> RCOOH',           'reactivity': 'good'},
    'epoxide':           {'rank': 7, 'product_type': 'alcohol (+2 carbons)',    'description': 'epoxide -> RCH2CH(OH)R', 'reactivity': 'good (regioselective at less hindered C)'},
    'nitrile':           {'rank': 8, 'product_type': 'ketone (after hydrolysis)','description': 'RC#N -> RC(=O)R',       'reactivity': 'slow, needs heat'},
    'amide':             {'rank': 9, 'product_type': 'complex (usually N/A)',   'description': 'low reactivity',      'reactivity': 'poor'},
}

INCOMPATIBLE_GROUPS = [
    ('O-H', 'hydroxyl', 'Acidic proton quenches Grignard'),
    ('N-H', 'amine/amide', 'Acidic proton (N-H) reacts with RMgX)'),
    ('-COOH', 'carboxylic acid', 'Acidic O-H proton; also C=O reacts'),
    ('SO3H', 'sulfonic acid', 'Strongly acidic — destroys Grignard'),
    ('-NO2', 'nitro group', 'Can oxidize/destroy Grignard reagent'),
    ('C=O in same mol', 'ketone/aldehyde in same molecule', 'Intramolecular reaction!'),
    ('-OXH (peroxide)', 'peroxide', 'Can cause dangerous side reactions'),
    ('active halide', '-CCl3, -CHCl2 etc.', 'May undergo Mg insertion / Wurtz coupling'),
]


@ChemMCPManager.register_tool
class GrignardMechanism(BaseTool):
    __version__ = "0.1.0"
    name = "GrignardMechanism"
    func_name = 'explain_grignard_mechanism'
    description = "Explain Grignard reagent reaction mechanism: oxidative addition of Mg into C-X bond to form RMgX, nucleophilic addition to carbonyl compounds (formaldehyde -> 1° alcohol, aldehyde -> 2°, ketone -> 3°, ester -> 3°+1°), tetrahedral intermediate, acidic workup, and functional group compatibility constraints."
    implementation_description = "Analyzes the carbonyl/substrate for Grignard suitability, classifies substrate type (aldehyde/ketone/ester/CO2/epoxide), provides complete mechanism including Grignard formation overview, nucleophilic addition step with stereochemistry, predicts product type based on carbonyl class, identifies incompatible functional groups, evaluates side reactions, and determines stoichiometry requirements."
    categories = ["Reaction"]
    tags = ["Grignard", "Organometallic", "Carbonyl Addition", "Alcohol Synthesis", "Nucleophilic Addition", "Mg"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('carbonyl_smiles', 'str', 'N/A', 'SMILES of the electrophile (aldehyde, ketone, ester, CO2 equivalent, epoxide).'),
        ('grignard_reagent', 'str', 'CH3MgBr', 'Grignard reagent formula. E.g., CH3MgBr, PhMgBr, RMgBr, vinyl-MgBr.'),
        ('solvent', 'str', 'dry Et2O or THF', 'Anhydrous solvent (must be aprotic and dry!).'),
        ('workup', 'str', 'H3O+', 'Workup conditions: H3O+, NH4Cl(aq), etc.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: carbonyl_smiles grignard_reagent [solvent] [workup]. E.g., "CC=O CH3MgBr".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing grignard_formation, substrate_classification, mechanism_steps, tetrahedral_intermediate, product_prediction, stoichiometry, compatibility_check, side_reactions, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'carbonyl_smiles': 'CC=O',
                'grignard_reagent': 'CH3MgBr',
                'solvent': 'dry Et2O',
                'workup': 'H3O+',
            },
            'text_input': {'query': 'CC=O CH3MgBr'},
            'output': {
                'result': {
                    'substrate': 'acetaldehyde (aldehyde)',
                    'grignard': 'methylmagnesium bromide (CH3MgBr)',
                    'mechanism_type': 'nucleophilic addition to aldehyde',
                    'product': '2-propanol (sec-butyl alcohol? no — isopropanol) (secondary alcohol)',
                    'stoichiometry': '1 equiv RMgX per carbonyl',
                    'key_step': 'CH3- (from CH3MgBr) attacks carbonyl C -> tetrahedral intermediate -> protonation -> alcohol',
                    'favorability': 'excellent — classic Grignard reaction',
                }
            },
        },
        {
            'code_input': {
                'carbonyl_smiles': 'CC(=O)C',
                'grignard_reagent': 'CH3MgBr',
                'solvent': 'dry THF',
                'workup': 'H3O+',
            },
            'text_input': {'query': 'CC(=O)C CH3MgBr THF'},
            'output': {
                'result': {
                    'substrate': 'acetone (ketone)',
                    'grignard': 'methylmagnesium bromide',
                    'product': 'tert-butyl alcohol (tertiary alcohol)',
                    'mechanism': 'CH3- attacks ketone C=O -> tertiary alkoxide -> tert-butanol after workup',
                    'sterics': 'Ketone is more sterically hindered than aldehyde — slower than aldehyde but still works well',
                    'favorability': 'good — ketones react readily with Grignards',
                }
            },
        },
        {
            'code_input': {
                'carbonyl_smiles': 'CC(=O)OCC',  # ethyl acetate
                'grignard_reagent': 'CH3MgBr',
                'solvent': 'dry Et2O',
                'workup': 'H3O+',
            },
            'text_input': {'query': 'ethyl_acetate CH3MgBr'},
            'output': {
                'result': {
                    'substrate': 'ethyl acetate (ester)',
                    'special_note': 'Esters react with 2 equivalents of Grignard!',
                    'step1': 'First RMgX adds -> ketone intermediate (cannot be isolated)',
                    'step2': 'Second RMgX adds to ketone -> tertiary alkoxide',
                    'product': '2-methylbutan-2-ol (tertiary alcohol) + ethanol (from -OEt elimination)',
                    'stoichiometry': '2 equiv RMgX required!',
                    'favorability': 'good — but must use excess Grignard',
                }
            },
        },
    ]

    def _run_base(self, carbonyl_smiles: str, grignard_reagent: str = 'CH3MgBr', solvent: str = 'dry Et2O or THF', workup: str = 'H3O+') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(carbonyl_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(carbonyl_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Classify substrate
        substrate_class = self._classify_substrate(mol)

        # 2. Grignard reagent analysis
        rmgx_info = self._analyze_grignard(grignard_reagent)

        # 3. Compatibility check
        compatibility = self._check_compatibility(mol)

        # 4. Build mechanism steps
        steps = self._build_mechanism(substrate_class, rmgx_info, workup)

        # 5. Tetrahedral intermediate analysis
        tetrahedral = self._analyze_tetrahedral(substrate_class)

        # 6. Product prediction
        products = self._predict_products(substrate_class, rmgx_info)

        # 7. Stoichiometry
        stoich = self._determine_stoichiometry(substrate_class)

        # 8. Side reactions
        side_rxns = self._identify_side_reactions(mol, substrate_class)

        # 9. Favorability
        favorability = self._evaluate_favorability(substrate_class, compatibility)

        result = {
            'result': {
                'substrate_smiles': carbonyl_smiles,
                'substrate_classification': substrate_class,
                'grignard_reagent': rmgx_info,
                'compatibility_check': compatibility,
                'mechanism_steps': steps,
                'tetrahedral_intermediate': tetrahedral,
                'product_prediction': products,
                'stoichiometry': stoich,
                'side_reactions': side_rxns,
                'favorability': favorability,
                'summary': self._build_summary(substrate_class, rmgx_info, products, favorability),
            }
        }

        logger.info(f"Grignard: {carbonyl_smiles} + {grignard_reagent} -> {favorability}")
        return result

    def _classify_substrate(self, mol):
        """Classify the electrophilic substrate."""
        n_carbonyl = 0
        has_epoxide = False
        has_c_n = False
        substrate_type = 'unknown'

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() >= 2.0:
                        if neighbor.GetAtomicNum() == 8:
                            n_carbonyl += 1
                            n_non_h = sum(1 for n in atom.GetNeighbors()
                                          if n.GetAtomicNum() != 1 and n.GetIdx() != neighbor.GetIdx())
                            # Check if it's an ester (C=O-O-C): carbonyl C has O neighbor that has C neighbor
                            for o_nbr in atom.GetNeighbors():
                                if o_nbr.GetAtomicNum() == 8 and o_nbr.GetIdx() != neighbor.GetIdx():
                                    for o_nbr_c in o_nbr.GetNeighbors():
                                        if o_nbr_c.GetAtomicNum() == 6 and o_nbr_c.GetIdx() != atom.GetIdx():
                                            substrate_type = 'ester'
                                            break
                                    if substrate_type == 'ester':
                                        break

                            if substrate_type != 'ester':
                                if n_non_h <= 1:
                                    substrate_type = 'aldehyde' if substrate_type == 'unknown' else substrate_type
                                else:
                                    substrate_type = 'ketone'

                        elif neighbor.GetAtomicNum() == 7 and bond.GetBondTypeAsDouble() >= 3.0:
                            has_c_n = True
                            substrate_type = 'nitrile'

            # Check for epoxide ring
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if len(ring) == 3:
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    syms = [(a.GetAtomicNum(), a.GetIdx()) for a in ring_atoms]
                    if any(s[0] == 8 for s in syms):
                        has_epoxide = True
                        substrate_type = 'epoxide'

        if n_carbonyl == 0 and not has_epoxide and not has_c_n:
            # Check for simple molecule like CO2 representation
            smiles = Chem.MolToSmiles(mol)
            if 'C(=O)=O' in smiles or 'O=C=O' in smiles:
                substrate_type = 'CO2'

        if substrate_type == 'unknown' and n_carbonyl > 0:
            substrate_type = 'carbonyl-containing'

        reactivity_info = CARBONYL_REACTIVITY.get(substrate_type, {})
        return {
            'type': substrate_type,
            'n_carbonyl_groups': n_carbonyl,
            'reactivity_rank': reactivity_info.get('rank', 99),
            'product_type': reactivity_info.get('product_type', 'unknown'),
            'reactivity': reactivity_info.get('reactivity', '?'),
        }

    def _analyze_grignard(self, reagent):
        """Analyze Grignard reagent."""
        r_clean = reagent.upper().replace(' ', '')

        organomagnesium_types = {
            'CH3MG': {'name': 'methylmagnesium bromide', 'character': 'primary, small, unhindered', 'pKa': '~ 48'},
            'PHMG':  {'name': 'phenylmagnesium bromide', 'character': 'aryl, planar, resonance-stabilized carbanion', 'pKa_CH': '~ 43'},
            'VINYL': {'name': 'vinyl Grignard', 'character': 'sp²-hybridized, retains alkene geometry', 'pKa_CH': '~ 44'},
            'ALLYL': {'name': 'allyl Grignard', 'character': 'resonance-stabilized, can rearrange', 'pKa_CH': '~ 43'},
        }

        for key, info in organomagnesium_types.items():
            if key in r_clean:
                return {**info, 'reagent_formula': reagent}

        return {
            'name': reagent,
            'character': 'organomagnesium halide (RMgX)',
            'reagent_formula': reagent,
            'pKa_CH': '~35-50 (depending on R group)',
        }

    def _check_compatibility(self, mol):
        """Check for incompatible functional groups."""
        issues = []
        compatible = True

        # Check for O-H groups
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 8:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        issues.append(('O-H', 'hydroxyl', 'Acidic proton will quench Grignard'))
                        compatible = False

            # Check for N-H
            if atom.GetAtomicNum() == 7:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        issues.append(('N-H', 'amine N-H', 'Acidic N-H proton reacts with RMgX'))
                        compatible = False

        return {
            'compatible': compatible,
            'issues': issues if issues else None,
            'note': None if compatible else (
                '⚠️ PROTECT or REMOVE acidic protons before using Grignard! '
                'Common protecting groups: TMS/TBS for OH, Boc for NH.'
            ),
        }

    def _build_mechanism(self, substrate, rmgx, workup='H3O+'):
        """Build Grignard mechanism steps."""
        r_name = rmgx.get('name', 'RMgX')
        sub_type = substrate.get('type', '')

        steps = [
            {
                'step': 0,
                'name': 'Grignard Reagent Formation (Brief Overview)',
                'equation': f'R-X + Mg (dry ether) -> R-Mg-X',
                'details': (
                    f"Metallic Mg inserts into the carbon-halogen bond via single-electron transfer (SET) mechanism "
                    f"in anhydrous diethyl ether or THF. The solvent coordinates to Mg, stabilizing the reagent. "
                    f"⚠️ Must be completely anhydrous and free of acidic protons!"
                ),
            },
            {
                'step': 1,
                'name': 'Nucleophilic Attack on Carbonyl Carbon',
                'equation': f'{r_name} + Substrate(C=O) -> Tetrahedral Alkoxide Intermediate',
                'details': (
                    f"The Grignard reagent acts as a source of carbanion (R:-), which is strongly nucleophilic. "
                    f"It attacks the electrophilic carbonyl carbon, breaking the C=O π bond and forming a new C-C σ bond. "
                    f"The carbonyl oxygen becomes an alkoxide (negatively charged)."
                ),
                'rate_determining': True,
            },
        ]

        # Ester-specific: two additions
        if sub_type == 'ester':
            steps.append({
                'step': 2,
                'name': 'Elimination (Ketone Formation)',
                'equation': 'Tetrahedral intermediate -> Ketone + -OMgX(OR\')',
                'details': (
                    "The tetrahedral intermediate collapses, expelling the alkoxide (-OR') leaving group. "
                    "This forms a ketone intermediate which is MORE reactive than the original ester toward Grignard. "
                    "The ketone cannot be isolated under normal conditions — it immediately reacts further."
                ),
            })
            steps.append({
                'step': 3,
                'name': 'Second Grignard Addition',
                'equation': f'Ketone + {r_name} -> New Tetrahedral Alkoxide (tertiary)',
                'details': "A second equivalent of RMgX adds to the ketone -> tertiary alkoxide.",
            })
            final_workup_step = 4
        else:
            final_workup_step = 2

        steps.append({
            'step': final_workup_step,
            'name': 'Acidic Workup',
            'equation': 'Alkoxide + H3O+ -> Alcohol Product',
            'details': (
                f"Dilute acid ({workup}) protonates the alkoxide to give the neutral alcohol product. "
                f"Mg salts are removed during aqueous workup."
            ),
        })

        return steps

    def _analyze_tetrahedral(self, substrate):
        """Analyze tetrahedral intermediate."""
        sub_type = substrate.get('type', '')
        prod_type = substrate.get('product_type', '?')

        return {
            'geometry': 'sp³-hybridized carbon (tetrahedral, ~109.5°)',
            'groups_attached': ['O- (alkoxide)', 'R (from Grignard)', 'original substituents from carbonyl'],
            'charge': '-1 on oxygen',
            'fate': 'Protonation upon workup -> alcohol' if sub_type != 'ester' else 'Collapse -> ketone -> second addition',
            'product_alcohol_class': prod_type,
        }

    def _predict_products(self, substrate, rmgx):
        """Predict final product(s)."""
        sub_type = substrate.get('type', '?')
        r_name = rmgx.get('name', 'R')
        prod_type = substrate.get('product_type', '?')

        products = {
            'aldehyde': f"{prod_type} — R'CH(OH)R (from {sub_type} + {r_name})",
            'ketone': f"{prod_type} — R'R''(OH)R (from ketone + {r_name})",
            'ester': f"{prod_type} + primary alcohol (from ester + 2 equiv {r_name}; -OR eliminated)",
            'formaldehyde': f"{prod_type} — RCH2OH (primary alcohol, from formaldehyde + {r_name})",
            'CO2': f"carboxylic acid RCOOH (after acidic workup; from CO2 + {r_name})",
            'epoxide': f"alcohol with chain extended by 2 carbons (anti addition across epoxide)",
            'nitrile': f"ketone R(C=O)R' (after acidic hydrolysis of imine intermediate)",
        }

        return {
            'main_product': products.get(sub_type, f'addition product from {sub_type}'),
            'alcohol_class': prod_type,
            'new_bonds_formed': '1 new C-C bond (esters: 2 C-C bonds)',
        }

    def _determine_stoichiometry(self, substrate):
        """Determine required Grignard equivalents."""
        sub_type = substrate.get('type', '')

        if sub_type == 'ester':
            return {'equiv_rm gx_needed': 2, 'reason': 'Ester -> ketone intermediate -> tertiary alcohol'}
        elif sub_type == 'acid_chloride':
            return {'equiv_rm gx_needed': 2, 'reason': 'Acid chloride -> aldehyde/ketone -> tertiary alcohol'}
        elif sub_type == 'CO2':
            return {'equiv_rm gx_needed': 1, 'reason': 'Single addition to CO2'}
        else:
            return {'equiv_rm gx_needed': 1, 'reason': 'Single addition to carbonyl'}

    def _identify_side_reactions(self, mol, substrate):
        """Identify potential side reactions."""
        reactions = []
        sub_type = substrate.get('type', '')

        if sub_type in ('aldehyde', 'ketone'):
            reactions.append({
                'reaction': 'Enolization or reduction',
                'conditions': 'Excess RMgX can act as base -> enolate formation (minor pathway)',
                'severity': 'minor',
            })

        if sub_type == 'ester':
            reactions.append({
                'reaction': 'Incomplete reaction (only 1 equiv added)',
                'conditions': 'Insufficient RMgX -> mixture of ketone + tertiary alcohol',
                'severity': 'manageable — use excess RMgX',
            })

        reactions.append({
            'reaction': 'Wurtz-type coupling',
            'conditions': 'RMgX + R\'X (halide impurity) -> R-R\' (homocoupling)',
            'severity': 'minor if pure substrates used',
        })

        reactions.append({
            'reaction': 'Quenching by moisture/air',
            'conditions': 'H2O or O2 in system -> ROH + Mg(OH)X (destroyed reagent)',
            'severity': 'catastrophic — must exclude moisture and air!',
        })

        return reactions

    def _evaluate_favorability(self, substrate, compatibility):
        score = 3
        rank = substrate.get('reactivity_rank', 99)
        if rank <= 2: score += 3
        elif rank <= 4: score += 2
        elif rank <= 6: score += 1

        if compatibility.get('compatible'): score += 2
        else: score -= 3

        if score >= 7: return 'excellent — should proceed cleanly'
        elif score >= 5: return 'good — standard conditions sufficient'
        elif score >= 3: return 'moderate — may need optimization'
        return 'problematic — check compatibility and substrate type'

    def _build_summary(self, substrate, rmgx, products, fav):
        sub_type = substrate.get('type', '?')
        r_name = rmgx.get('name', 'RMgX')
        prod = products.get('main_product', '?')
        return f"Grignard reaction: {sub_type} + {r_name} -> {prod}. Favorability: {fav}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        carb = parts[0] if len(parts) > 0 else ''
        rmgx = parts[1] if len(parts) > 1 else 'CH3MgBr'
        solv = parts[2] if len(parts) > 2 else 'dry Et2O or THF'
        wk = parts[3] if len(parts) > 3 else 'H3O+'
        return self._run_base(carb, rmgx, solv, wk)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
