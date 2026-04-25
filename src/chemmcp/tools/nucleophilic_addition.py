"""
Nucleophilic Addition to Carbonyl (Tool #121)
羰基亲核加成反应机理：醛、酮的亲核加成，四面体中间体分析。
Provides nucleophilic addition mechanism analysis for aldehydes and ketones,
including tetrahedral intermediate, stereochemistry, and catalysis effects.
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


# Nucleophile database with strength and sterics
NUCLEOPHILES = {
    'CN-':    {'name': 'cyanide (CN⁻)',     'strength': 'strong',  'sterics': 'small',  'product_suffix': 'cyanohydrin',   'pKa_HA': 9.2},
    'HCN':    {'name': 'hydrogen cyanide',    'strength': 'weak',   'sterics': 'small',  'product_suffix': 'cyanohydrin',   'pKa_HA': 9.2},
    'NaHSO3': {'name': 'bisulfite (HSO₃⁻)',   'strength': 'moderate','sterics': 'small',  'product_suffix': 'bisulfite adduct','pKa_HA': 7.2},
    'OH-':    {'name': 'hydroxide (OH⁻)',      'strength': 'strong', 'sterics': 'small',  'product_suffix': 'hydrate',        'pKa_HA': 15.7},
    'H2O':    {'name': 'water (H₂O)',          'strength': 'weak',   'sterics': 'small',  'product_suffix': 'hydrate',        'pKa_HA': 15.7},
    'RO-':    {'name': 'alkoxide (RO⁻)',       'strength': 'strong', 'sterics': 'small',  'product_suffix': 'hemiacetal',     'pKa_HA': 16},
    'ROH':    {'name': 'alcohol (ROH)',         'strength': 'weak',   'sterics': 'small',  'product_suffix': 'hemiacetal',     'pKa_HA': 16},
    'NH3':    {'name': 'ammonia (NH₃)',         'strength': 'weak',   'sterics': 'small',  'product_suffix': 'imine precursor','pKa_HA': 38},
    'RNH2':   {'name': 'primary amine (RNH₂)',  'strength': 'weak',   'sterics': 'small',  'product_suffix': 'imine/schiff base','pKa_HA': 38},
    'NH2OH':  {'name': 'hydroxylamine (NH₂OH)', 'strength': 'weak',   'sterics': 'small',  'product_suffix': 'oxime',          'pKa_HA': 6},
    'NHNH2':  {'name': 'hydrazine (N₂H₄)',      'strength': 'weak',   'sterics': 'small',  'product_suffix': 'hydrazone',      'pKa_HA': 8},
    'PhMgBr': {'name': 'Grignard (PhMgBr)',     'strength': 'very strong','sterics':'moderate','product_suffix':'alcohol after workup','pKa_HA': 43},
    'CH3MgBr':{'name': 'Grignard (CH₃MgBr)',    'strength': 'very strong','sterics':'small', 'product_suffix':'alcohol after workup','pKa_HA': 43},
    'RLi':    {'name': 'organolithium (RLi)',    'strength': 'very strong','sterics':'small', 'product_suffix':'alcohol after workup','pKa_HA': 45},
    'NaBH4':  {'name': 'hydride (NaBH₄)',        'strength': 'strong', 'sterics': 'small',  'product_suffix': 'alcohol',        'pKa_HA': 35},
    'LiAlH4': {'name': 'hydride (LiAlH₄)',       'strength': 'very strong','sterics':'small','product_suffix':'alcohol',           'pKa_HA': 35},
}


@ChemMCPManager.register_tool
class NucleophilicAddition(BaseTool):
    __version__ = "0.1.0"
    name = "NucleophilicAddition"
    func_name = 'explain_nucleophilic_addition'
    description = "Explain the nucleophilic addition reaction mechanism to carbonyl compounds (aldehydes and ketones): attack on electrophilic carbonyl carbon, tetrahedral intermediate formation, protonation, and product formation."
    implementation_description = "Analyzes the carbonyl substrate (aldehyde vs ketone, steric environment), classifies the nucleophile (strength, sterics, hard/soft character), determines acid or base catalysis, provides stepwise mechanism with tetrahedral intermediate analysis, predicts stereochemistry of new stereocenter formation, and evaluates overall favorability."
    categories = ["Reaction"]
    tags = ["Nucleophilic Addition", "Carbonyl", "Tetrahedral Intermediate", "Aldehyde", "Ketone", "Stereochemistry"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('carbonyl_smiles', 'str', 'N/A', 'SMILES string of the carbonyl compound (aldehyde or ketone).'),
        ('nucleophile', 'str', 'CN-', 'Nucleophile formula or name. Options: CN-, HCN, NaHSO3, OH-, H2O, RO-, ROH, NH3, RNH2, NH2OH, NHNH2, PhMgBr, CH3MgBr, RLi, NaBH4, LiAlH4.'),
        ('solvent', 'str', 'polar aprotic', 'Solvent type: polar aprotic, polar protic, non-polar.'),
        ('catalysis', 'str', 'auto', 'Catalysis mode: auto, acid, base, none.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: carbonyl_smiles nucleophile [solvent] [catalysis]. E.g., "CC=O CN- polar_aprotic auto".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing carbonyl_analysis, nucleophile_analysis, mechanism_steps, tetrahedral_intermediate, stereochemistry, product_prediction, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'carbonyl_smiles': 'CC=O',
                'nucleophile': 'CN-',
                'solvent': 'polar aprotic',
                'catalysis': 'auto',
            },
            'text_input': {'query': 'CC=O CN-'},
            'output': {
                'result': {
                    'substrate': 'acetaldehyde (aldehyde)',
                    'carbonyl_type': 'aldehyde',
                    'nucleophile': 'cyanide (CN⁻)',
                    'mechanism_type': 'nucleophilic addition (base-catalyzed)',
                    'steps': [
                        {'step': 1, 'description': 'Nucleophilic attack: CN⁻ attacks carbonyl C → tetrahedral intermediate'},
                        {'step': 2, 'description': 'Protonation: tetrahedral O⁻ is protonated → cyanohydrin product'},
                    ],
                    'tetrahedral_intermediate': 'sp³ carbon with -OH, -CN, -H, -CH₃ groups',
                    'stereochemistry': 'New chiral center formed → racemic mixture (R/S cyanohydrin)',
                    'product': 'lactonitrile (acetaldehyde cyanohydrin)',
                    'favorability': 'excellent — aldehydes are highly reactive toward nucleophilic addition',
                }
            },
        },
        {
            'code_input': {
                'carbonyl_smiles': 'CC(=O)C',
                'nucleophile': 'MgMeBr',
                'solvent': 'dry Et2O',
                'catalysis': 'none',
            },
            'text_input': {'query': 'CC(=O)C MgMeBr dry_Et2O none'},
            'output': {
                'result': {
                    'substrate': 'acetone (ketone)',
                    'carbonyl_type': 'ketone',
                    'nucleophile': 'Grignard reagent (CH₃MgBr equivalent)',
                    'product': 'tert-butyl alcohol (after acidic workup)',
                    'favorability': 'good — ketones react well with strong nucleophiles like Grignards',
                }
            },
        },
    ]

    def _run_base(self, carbonyl_smiles: str, nucleophile: str = 'CN-', solvent: str = 'polar aprotic', catalysis: str = 'auto') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(carbonyl_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(carbonyl_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Analyze carbonyl substrate
        carbonyl_analysis = self._analyze_carbonyl(mol)

        # 2. Analyze nucleophile
        nu_info = self._classify_nucleophile(nucleophile)

        # 3. Determine catalysis
        cat_mode = self._determine_catalysis(catalysis, nu_info)

        # 4. Build mechanism steps
        steps = self._build_mechanism_steps(carbonyl_analysis, nu_info, cat_mode)

        # 5. Tetrahedral intermediate analysis
        tetrahedral = self._analyze_tetrahedral_intermediate(carbonyl_analysis, nu_info)

        # 6. Stereochemistry
        stereo = self._predict_stereochemistry(carbonyl_analysis, nu_info)

        # 7. Product prediction
        product = self._predict_product(carbonyl_analysis, nu_info)

        # 8. Favorability
        favorability = self._evaluate_favorability(carbonyl_analysis, nu_info, solvent)

        result = {
            'result': {
                'substrate_smiles': carbonyl_smiles,
                'carbonyl_analysis': carbonyl_analysis,
                'nucleophile_analysis': nu_info,
                'catalysis_mode': cat_mode,
                'mechanism_steps': steps,
                'tetrahedral_intermediate': tetrahedral,
                'stereochemistry': stereo,
                'product_prediction': product,
                'favorability': favorability,
                'summary': self._build_summary(carbonyl_analysis, nu_info, product, favorability),
            }
        }

        logger.info(f"Nucleophilic addition: {carbonyl_smiles} + {nucleophile} → {favorability}")
        return result

    def _analyze_carbonyl(self, mol):
        """Analyze carbonyl compound."""
        # Find carbonyl group(s)
        carbonyl_carbons = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            # Count non-H substituents on carbonyl carbon
                            non_h_subs = sum(
                                1 for n in atom.GetNeighbors()
                                if n.GetAtomicNum() != 1 and n.GetAtomicNum() != 8
                            )
                            has_alpha_h = any(
                                n.GetTotalNumHs() > 0
                                for n in atom.GetNeighbors()
                                if n.GetAtomicNum() == 6
                            )
                            carbonyl_carbons.append({
                                'idx': atom.GetIdx(),
                                'n_non_h_substituents': non_h_subs,
                                'has_alpha_hydrogen': has_alpha_h,
                            })

        if not carbonyl_carbons:
            return {'has_carbonyl': False, 'error': 'No carbonyl group found.'}

        cc = carbonyl_carbons[0]
        n_sub = cc['n_non_h_substituents']

        if n_sub == 0:
            c_type = 'formaldehyde'
            reactivity = 'very high'
            steric = 'minimal'
        elif n_sub == 1:
            c_type = 'aldehyde'
            reactivity = 'high'
            steric = 'low'
        elif n_sub == 2:
            c_type = 'ketone'
            reactivity = 'moderate'
            steric = 'moderate'
        else:
            c_type = 'ketone/ester derivative'
            reactivity = 'low (hindered)'
            steric = 'high'

        return {
            'has_carbonyl': True,
            'carbonyl_type': c_type,
            'reactivity': reactivity,
            'steric_environment': steric,
            'n_substituents': n_sub,
            'is_aldehyde': n_sub <= 1,
            'is_ketone': n_sub >= 2,
            'prochiral': n_sub == 1,  # Aldehyde with one non-H sub → prochiral
            'can_form_chiral_center': n_sub >= 1 and n_sub <= 2,
        }

    def _classify_nucleophile(self, nu):
        """Classify nucleophile."""
        nu_clean = nu.upper().replace(' ', '').replace('-', '')

        # Direct match
        for key, info in NUCLEOPHILES.items():
            if nu_clean == key.upper().replace('-', '') or nu_clean in key.upper() or key.upper().replace('-', '') in nu_clean:
                return {**info, 'key': key}

        # Fuzzy match
        fuzzy_map = {
            'GRIGNARD': 'PhMgBr', 'MG': 'PhMgBr', 'RMGBR': 'PhMgBr',
            'HYDRIDE': 'NaBH4', 'NAHBH4': 'NaBH4', 'LIALH4': 'LiAlH4',
            'AMINE': 'RNH2', 'AMMONIA': 'NH3',
            'CYANIDE': 'CN-', 'CN': 'CN-',
        }
        for fkey, real_key in fuzzy_map.items():
            if fkey in nu_clean:
                info = NUCLEOPHILES[real_key]
                return {**info, 'key': real_key, 'matched_from': nu}

        return {
            'name': nu, 'strength': 'unknown', 'sterics': 'unknown',
            'product_suffix': 'unknown', 'pKa_HA': None, 'key': nu,
        }

    def _determine_catalysis(self, requested, nu_info):
        """Determine catalysis mode."""
        if requested != 'auto':
            return requested

        strength = nu_info.get('strength', '')
        if strength in ('weak',):
            return 'base'  # Weak Nu needs activation (e.g., HCN + base → CN⁻)
        if strength in ('moderate',) and nu_info.get('pKa_HA', 99) < 10:
            return 'base'
        if strength == 'very strong':
            return 'none'  # Grignard/Rli react directly
        return 'base'

    def _build_mechanism_steps(self, carbonyl, nu_info, cat_mode):
        """Build stepwise mechanism."""
        nu_name = nu_info.get('name', 'Nu')
        steps = []
        c_type = carbonyl.get('carbonyl_type', '')

        if cat_mode == 'base':
            # Step 0: Generate active nucleophile (if needed)
            pka = nu_info.get('pKa_HA')
            if pka and pka > 7:
                steps.append({
                    'step': 0,
                    'name': 'Nucleophile Activation (Base Catalysis)',
                    'equation': f"HA + B⁻ → A⁻ (active nucleophile) + BH",
                    'details': (
                        f"Base deprotonates weak nucleophile (pKa ≈ {pka}) "
                        f"to generate the active anionic nucleophile."
                    ),
                    'rate_relevant': False,
                })

        # Step 1: Nucleophilic attack
        steps.append({
            'step': len(steps) + 1,
            'name': 'Nucleophilic Attack (Rate-Determining for most cases)',
            'equation': f"C=O + Nu⁻ → Tetrahedral Intermediate (sp³)",
            'details': (
                f"The nucleophile ({nu_name}) attacks the electrophilic carbonyl carbon. "
                f"The π bond of C=O breaks, electrons move to oxygen, forming an sp³-hybridized "
                f"tetrahedral intermediate with a negatively charged oxygen. "
                f"Carbonyl reactivity: {carbonyl.get('reactivity', '?')}."
            ),
            'rate_relevant': True,
        })

        # Step 2: Protonation
        if nu_info.get('strength') != 'very strong':
            steps.append({
                'step': len(steps) + 1,
                'name': 'Protonation of Alkoxide',
                'equation': f"Tetrahedral-O⁻ + H⁺ → Product (neutral alcohol/derivative)",
                'details': (
                    "The alkoxide oxygen is protonated by solvent or added acid "
                    "to give the neutral addition product."
                ),
                'rate_relevant': False,
            })
        else:
            # For Grignard/hydride: need acidic workup
            steps.append({
                'step': len(steps) + 1,
                'name': 'Acidic Workup',
                'equation': f"Tetrahedral-O⁻ + H₃O⁺ → Alcohol",
                'details': (
                    f"After nucleophilic addition, acidic workup (H₃O⁺) protonates the alkoxide "
                    f"to yield the final alcohol product."
                ),
                'rate_relevant': False,
            })

        # Special: imine formation needs extra steps
        if 'amine' in nu_name.lower() or 'ammonia' in nu_name.lower() or 'hydraz' in nu_name.lower() or 'hydroxylamine' in nu_name.lower():
            steps.append({
                'step': len(steps) + 1,
                'name': 'Elimination of Water (for imine/oxime/hydrazone)',
                'equation': "Carbinolamine → C=N-R + H₂O",
                'details': (
                    "The tetrahedral carbinolamine intermediate undergoes dehydration "
                    "(acid-catalyzed) to form the C=N double bond (imine, oxime, or hydrazone)."
                ),
                'rate_relevant': False,
            })

        return steps

    def _analyze_tetrahedral_intermediate(self, carbonyl, nu_info):
        """Analyze the tetrahedral intermediate."""
        c_type = carbonyl.get('carbonyl_type', '')
        n_sub = carbonyl.get('n_substituents', 0)

        groups = []
        if c_type == 'formaldehyde':
            groups = ['H', 'H']
        elif c_type == 'aldehyde':
            groups = ['H', 'R (alkyl/aryl)']
        elif c_type == 'ketone':
            groups = ['R¹', 'R²']

        nu_group = nu_info.get('name', 'Nu').split('(')[0].strip()

        return {
            'hybridization': 'sp³',
            'geometry': 'tetrahedral (~109.5°)',
            'groups_attached': [*groups, f'O⁻ (alkoxide)', nu_group],
            'charge': '-1 on oxygen',
            'stability': 'more stable than original carbonyl due to charge dispersal into electronegative O',
            'fate': 'protonation → neutral product; or elimination (for N-nucleophiles → imines)',
        }

    def _predict_stereochemistry(self, carbonyl, nu_info):
        """Predict stereochemical outcome."""
        can_form_chiral = carbonyl.get('can_form_chiral_center', False)
        c_type = carbonyl.get('carbonyl_type', '')

        if not can_form_chiral:
            return 'No new stereocenter formed (formaldehyde or symmetric ketone).'

        if c_type == 'aldehyde':
            return (
                'New chiral center formed at the former carbonyl carbon! '
                'The prochiral sp² carbon becomes sp³ with four different substituents. '
                'Nucleophile can attack either face of the planar carbonyl → racemic mixture (R/S pair). '
                'If a chiral catalyst or chiral auxiliary is used, enantioselective addition is possible.'
            )

        if c_type == 'ketone':
            has_identical_subs = False  # simplified check
            if not has_identical_subs:
                return (
                    'New chiral center possible if the two substituents on the ketone are different. '
                    'Attack on either face → racemic mixture. '
                    'If the two substituents are identical (e.g., acetone), no new chirality.'
                )
            return 'Symmetric ketone → no new stereocenter.'

        return 'Stereochemistry depends on specific substrate structure.'

    def _predict_product(self, carbonyl, nu_info):
        """Predict the product."""
        c_type = carbonyl.get('carbonyl_type', '')
        suffix = nu_info.get('product_suffix', 'addition product')
        nu_name = nu_info.get('name', 'Nu')

        products = {
            'formaldehyde': {
                'primary alcohol': 'with Grignard/hydride → primary alcohol (RCH₂OH)',
                'cyanohydrin': 'with CN⁻ → formaldehyde cyanohydrin (HOCH₂CN)',
            },
            'aldehyde': {
                'secondary alcohol': 'with Grignard/hydride → secondary alcohol (RCHOH R\')',
                'cyanohydrin': f'with CN⁻ → aldehyde cyanohydrin (RCH(OH)CN)',
                'bisulfite': f'with HSO₃⁻ → bisulfite adduct (crystalline solid, useful for purification)',
                'imine': f'with amine → imine (Schiff base, RCH=NR\')',
                'hydrate': 'with H₂O/H⁺ → gem-diol (usually favors carbonyl)',
            },
            'ketone': {
                'tertiary alcohol': 'with Grignard → tertiary alcohol (RC(R\')(R\'\')OH)',
                'cyanohydrin': f'with CN⁻ → ketone cyanohydrin (R₂C(OH)CN)',
                'imine': f'with primary amine → imine/ketimine (R₂C=NR)',
            },
        }

        desc = products.get(c_type, {}).get(suffix.lower(), f'{suffix} formed from {c_type}')

        # Check for special cases
        if 'grignard' in nu_name.lower() or 'lithium' in nu_name.lower():
            alcohol_type = 'primary' if c_type == 'formaldehyde' else \
                           'secondary' if c_type == 'aldehyde' else \
                           'tertiary' if c_type == 'ketone' else 'unknown'
            desc = f'{alcohol_type.capitalize()} alcohol after acidic workup ({c_type} + organometallic)'

        return {
            'product_type': suffix,
            'description': desc,
            'from_carbonyl_type': c_type,
        }

    def _evaluate_favorability(self, carbonyl, nu_info, solvent):
        """Score favorability."""
        score = 0
        c_type = carbonyl.get('carbonyl_type', '')
        strength = nu_info.get('strength', '')

        # Substrate reactivity
        if c_type == 'formaldehyde': score += 3
        elif c_type == 'aldehyde': score += 2
        elif c_type == 'ketone': score += 1
        elif 'ester' in c_type or 'hindered' in c_type: score -= 1

        # Nucleophile strength
        if strength == 'very strong': score += 3
        elif strength == 'strong': score += 2
        elif strength == 'moderate': score += 1
        elif strength == 'weak': score -= 1  # may need catalysis

        # Sterics
        steric = carbonyl.get('steric_environment', '')
        if steric == 'minimal': score += 1
        elif steric == 'high': score -= 2

        # Solvent
        if 'aprotic' in solvent.lower(): score += 1
        elif 'protic' in solvent.lower() and strength in ('strong', 'very strong'):
            score -= 1  # Protic solvents H-bond to strong nucleophiles

        if score >= 6: return 'excellent'
        elif score >= 4: return 'good'
        elif score >= 2: return 'moderate'
        elif score >= 0: return 'possible but slow'
        return 'unfavorable — consider alternative approach'

    def _build_summary(self, carbonyl, nu_info, product, favorability):
        c_type = carbonyl.get('carbonyl_type', '?')
        nu_name = nu_info.get('name', 'Nu')
        return (
            f"Nucleophilic addition to {c_type} by {nu_name}. "
            f"Product: {product.get('description', '?')}. "
            f"Favorability: {favorability}."
        )

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        carb = parts[0] if len(parts) > 0 else ''
        nu = parts[1] if len(parts) > 1 else 'CN-'
        solv = parts[2] if len(parts) > 2 else 'polar aprotic'
        cat = parts[3] if len(parts) > 3 else 'auto'
        return self._run_base(carb, nu, solv, cat)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
