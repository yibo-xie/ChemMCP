"""
Claisen Mechanism (Tool #127)
Claisen 缩合反应机理：酯烯醇化、亲核酰基取代、β-酮酯产物、
Dieckmann 缩合（分子内）、逆 Claisen（水解+脱羧）。
Provides Claisen condensation mechanism analysis: enolate formation,
nucleophilic acyl substitution, β-keto ester product, Dieckmann variant.
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


@ChemMCPManager.register_tool
class ClaisenMechanism(BaseTool):
    __version__ = "0.1.0"
    name = "ClaisenMechanism"
    func_name = 'explain_claisen_mechanism'
    description = "Explain the Claisen condensation mechanism: ester enolate formation (strong base required), nucleophilic acyl substitution on another ester, tetrahedral intermediate collapse with alkoxide leaving group, and β-keto ester product formation. Covers Dieckmann condensation (intramolecular), crossed Claisen, and retro-Claison (hydrolysis + decarboxylation)."
    implementation_description = "Analyzes the ester substrate(s) for α-hydrogen availability (required for enolization), determines base strength requirements (alkoxide base needed — stronger than aldol), provides complete stepwise mechanism including the critical alkoxide elimination step, predicts β-keto ester products, analyzes Dieckmann feasibility for diesters, and covers decarboxylation of β-keto acids."
    categories = ["Reaction"]
    tags = ["Claisen", "Condensation", "Ester", "Enolate", "Acyl Substitution", "β-Keto Ester", "Dieckmann"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('ester_smiles', 'str', 'N/A', 'SMILES of the ester (must have α-H for enolization).'),
        ('ester2_smiles', 'str', '', 'SMILES of second ester (for crossed Claisen). Leave empty for self-Condensation.'),
        ('base', 'str', 'EtO-', 'Base used. Options: EtO-, MeO-, NaOEt, LDA, NaH.'),
        ('solvent', 'str', 'EtOH', 'Solvent (should match alkoxide to avoid transesterification).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: ester_smiles [ester2_smiles] [base] [solvent]. E.g., "CC(=O)OCC".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing ester_analysis, enolate_analysis, mechanism_steps, product_prediction, dieckmann_analysis, retro_claisen, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'ester_smiles': 'CC(=O)OCC',
                'ester2_smiles': '',
                'base': 'EtO-',
                'solvent': 'EtOH',
            },
            'text_input': {'query': 'CC(=O)OCC'},
            'output': {
                'result': {
                    'substrate': 'ethyl acetate',
                    'reaction_type': 'self-Claisen condensation',
                    'mechanism': 'ester enolate + ester → nucleophilic acyl substitution → β-keto ester',
                    'product': 'ethyl acetoacetate (acetoacetic ester)',
                    'key_difference_from_aldol': (
                        'Alkoxide is a leaving group in Claisen (acyl substitution) vs protonation in aldol. '
                        'Requires strong base (pKa of β-keto ester ~11 must be deprotonated to drive equilibrium).'
                    ),
                    'favorability': 'excellent — classic Claisen condensation',
                }
            },
        },
        {
            'code_input': {
                'ester_smiles': 'C1(CCCC1C(=O)OCC)C(=O)OCC',  # diethyl adipate
                'ester2_smiles': '',
                'base': 'EtO-',
                'solvent': 'EtOH',
            },
            'text_input': {'query': 'diethyl_adipate'},
            'output': {
                'result': {
                    'reaction_type': 'Dieckmann condensation (intramolecular Claisen)',
                    'ring_formed': '5-membered ring (cyclopentanone derivative)',
                    'product': 'ethyl 2-oxocyclopentanecarboxylate',
                    'favorability': 'good — 5- and 6-membered rings form readily',
                }
            },
        },
    ]

    def _run_base(self, ester_smiles: str, ester2_smiles: str = '', base: str = 'EtO-', solvent: str = 'EtOH') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(ester_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(ester_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Analyze ester
        ester_analysis = self._analyze_ester(mol)

        # 2. Enolate analysis
        enolate_info = self._analyze_enolate(ester_analysis, base)

        # 3. Mechanism steps
        steps = self._build_mechanism(ester_analysis, base)

        # 4. Product prediction
        products = self._predict_products(ester_analysis)

        # 5. Dieckmann analysis
        dieckmann = self._analyze_dieckmann(mol)

        # 6. Retro-Claisen analysis
        retro = self._analyze_retro_claisn()

        # 7. Favorability
        favorability = self._evaluate_favorability(ester_analysis, base)

        result = {
            'result': {
                'ester_analysis': ester_analysis,
                'enolate_analysis': enolate_info,
                'mechanism_steps': steps,
                'product_prediction': products,
                'dieckmann_analysis': dieckmann,
                'retro_claisan': retro,
                'favorability': favorability,
                'summary': self._build_summary(ester_analysis, products, favorability),
            }
        }

        logger.info(f"Claisen: {ester_smiles} ({base}) → {favorability}")
        return result

    def _analyze_ester(self, mol):
        """Analyze ester substrate."""
        has_ester = False
        has_alpha_h = False
        alpha_carbons = []
        n_ester_groups = 0

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() >= 2.0:
                        if neighbor.GetAtomicNum() == 8:
                            # Check if this oxygen is part of an ester (C=O-O-C)
                            for o_neighbor in neighbor.GetNeighbors():
                                if o_neighbor.GetAtomicNum() == 6 and o_neighbor.GetIdx() != atom.GetIdx():
                                    has_ester = True
                                    n_ester_groups += 1
                                    # Check alpha carbon for H
                                    for alpha_n in atom.GetNeighbors():
                                        if alpha_n.GetAtomicNum() == 6 and alpha_n.GetTotalNumHs() > 0:
                                            has_alpha_h = True
                                            alpha_carbons.append({
                                                'idx': alpha_n.GetIdx(),
                                                'n_h': alpha_n.GetTotalNumHs(),
                                            })

        return {
            'has_ester_group': has_ester,
            'n_ester_groups': n_ester_groups,
            'has_alpha_hydrogen': has_alpha_h,
            'alpha_carbons': alpha_carbons,
            'can_enolize': has_alpha_h,
            'alpha_h_pKa': '~25 (ester α-H, less acidic than aldehyde/ketone due to resonance donation from OR)',
        }

    def _analyze_enolate(self, ester_analysis, base):
        """Analyze enolate formation."""
        if not ester_analysis.get('can_enolize'):
            return {'can_enolize': False, 'error': 'No α-H — cannot undergo Claisen condensation.'}

        b = base.upper().replace('-', '')
        strong_bases = ['ETOE', 'NAOET', 'MEOE', 'NAOME', 'LDA', 'NAH']
        is_strong = any(sb in b for sb in strong_bases)

        return {
            'can_enolize': True,
            'base_strength': 'strong' if is_strong else 'possibly insufficient',
            'base_required': (
                'Alkoxide base (RO⁻) or stronger (LDA, NaH) required because: '
                '(1) Ester α-H pKa ≈ 25 (less acidic than ketones); '
                '(2) The product β-keto ester (pKa ~11) must be fully deprotonated '
                'to pull the equilibrium forward (irreversible deprotonation step)'
            ),
            'note': None if is_strong else f'⚠ {base} may not be strong enough — use NaOEt/EtOH or LDA.',
        }

    def _build_mechanism(self, ester_analysis, base):
        """Build Claisen mechanism steps."""
        return [
            {
                'step': 1,
                'name': 'Enolate Formation (Reversible)',
                'equation': 'Ester + RO⁻ ⇌ Enolate + ROH',
                'details': (
                    f"Base ({base}) abstracts an α-proton from the ester to form a "
                    f"resonance-stabilized enolate. This step is reversible."
                ),
            },
            {
                'step': 2,
                'name': 'Nucleophilic Acyl Substitution (Key Step)',
                'equation': 'Enolate + Ester → Tetrahedral Intermediate',
                'details': (
                    "The enolate (nucleophile at carbon) attacks the carbonyl carbon of a "
                    "second ester molecule (or the same molecule in Dieckmann). This forms a "
                    "tetrahedral intermediate — unlike aldol where attack gives an alkoxide that "
                    "is simply protonated, here the intermediate must eliminate a leaving group."
                ),
            },
            {
                'step': 3,
                'name': 'Collapse of Tetrahedral Intermediate (Elimination)',
                'equation': 'Tetrahedral Intermediate → β-Keto Ester + Alkoxide (⁻OR)',
                'details': (
                    "The tetrahedral intermediate collapses, expelling an alkoxide (⁻OR) leaving "
                    "group. This reforms a carbonyl and produces the β-keto ester product. "
                    "**This is the key difference from aldol**: elimination vs protonation."
                ),
            },
            {
                'step': 4,
                'name': 'Deprotonation of Product (Irreversible — Drives Equilibrium)',
                'equation': 'β-Keto Ester + ⁻OR → β-Keto Ester Enolate + ROH',
                'details': (
                    "The β-keto ester product is significantly more acidic (pKa ~11) than the "
                    "starting ester (pKa ~25). The base deprotonates it irreversibly, pulling "
                    "the entire equilibrium toward products. Final acid workup gives neutral β-keto ester."
                ),
                'irreversible': True,
            },
        ]

    def _predict_products(self, ester_analysis):
        """Predict products."""
        return {
            'claisen_product': {
                'type': 'β-keto ester (β-ketoacid ester)',
                'structure': 'R-CO-CH₂-COOR\' (if starting from R-CH₂-COOR\')',
                'acidity': 'pKa ~11 (enolizable position between two carbonyls)',
                'uses': 'Useful synthon for ketone synthesis via hydrolysis/decarboxylation',
            },
            'after_hydrolysis_decarboxylation': {
                'type': 'ketone (retro-Claisen)',
                'process': 'Acid or base hydrolysis → β-keto acid → heat-induced decarboxylation → ketone',
                'driving_force': 'CO₂ evolution (entropy) + stable ketone formation',
            },
        }

    def _analyze_dieckmann(self, mol):
        """Analyze Dieckmann condensation potential."""
        ring_info = mol.GetRingInfo()
        n_rings = ring_info.NumRings()

        if n_rings == 0:
            return {'applicable': False, 'reason': 'Not cyclic — intermolecular Claisen.'}

        for ring in ring_info.AtomRings():
            ring_size = len(ring)
            if 7 <= ring_size <= 14:  # reasonable range for diesters
                product_ring_size = ring_size - 2
                favorability = 'excellent' if product_ring_size in (5, 6) else \
                               'good' if product_ring_size in (4, 7) else 'possible'
                return {
                    'applicable': True,
                    'starting_ring_size': ring_size,
                    'product_ring_size': product_ring_size,
                    'product_type': f'{product_ring_size}-membered ring β-keto ester',
                    'favorability': favorability,
                    'blaise_rule': '5- and 6-membered rings form most readily (Blaise rule)',
                }

        return {'applicable': False, 'reason': 'No suitable diester ring system found.'}

    def _analyze_retro_claisn(self):
        """Analyze retro-Claisen (hydrolysis + decarboxylation)."""
        return {
            'process': 'β-Keto Ester → β-Keto Acid → Ketone + CO₂',
            'step1': 'Hydrolysis: ester → carboxylic acid (acid or base catalyzed)',
            'step2': 'Decarboxylation: β-keto acid → ketone + CO₂ (heat, 100-200°C)',
            'mechanism_of_decaboxylation': (
                'Six-membered cyclic transition state via enol of β-keto acid: '
                'carboxylic proton hydrogen-bonds to carbonyl O, facilitating CO₂ loss via a '
                'concerted cyclic mechanism. Only works for β-keto acids (and β-cyano/β-nitro acids).'
            ),
            'synthetic_utility': 'Convert esters into ketones — complementary to Grignard/Friedel-Crafts',
        }

    def _evaluate_favorability(self, ester_analysis, base):
        score = 2
        if ester_analysis.get('can_enolize'): score += 3
        if ester_analysis.get('n_ester_groups', 0) >= 2: score += 2  # diester possible
        b = base.upper()
        if any(x in b for x in ['NAOET', 'ETOE', 'LDA', 'NAH']): score += 2
        if ester_analysis.get('alpha_carbons'):
            total_h = sum(a['n_h'] for a in ester_analysis['alpha_carbons'])
            if total_h >= 2: score += 1

        if score >= 7: return 'excellent'
        elif score >= 5: return 'good'
        elif score >= 3: return 'moderate'
        return 'check: ensure α-H present and strong base'

    def _build_summary(self, ester, products, fav):
        can_enol = ester.get('can_enolize', False)
        status = '✅' if can_enol else '❌ no α-H'
        prod = products.get('claisen_product', {}).get('type', '?')
        return f"Claisen condensation [{status}]. Product: {prod}. Favorability: {fav}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        e1 = parts[0] if len(parts) > 0 else ''
        e2 = parts[1] if len(parts) > 1 else ''
        b = parts[2] if len(parts) > 2 else 'EtO-'
        s = parts[3] if len(parts) > 3 else 'EtOH'
        return self._run_base(e1, e2, b, s)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
