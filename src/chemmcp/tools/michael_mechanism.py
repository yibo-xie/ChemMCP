"""
Michael Mechanism (Tool #128)
Michael 加成反应机理：共轭（1,4-）加成 vs 直接（1,2-）加成、
软硬酸碱分析、Robinson 增环化。
Provides Michael addition mechanism analysis: conjugate (1,4-) addition to α,β-unsaturated
carbonyl compounds vs direct (1,2-) addition, HSAB analysis, and Robinson annulation context.
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


# Michael acceptor reactivity order (more reactive = better acceptor)
ACCEPTOR_REACTIVITY = {
    'nitroalkene':       ('exceptional', '-NO₂ strongly electron-withdrawing + stabilizes enolate'),
    'enone':             ('excellent', 'conjugated ketone — classic Michael acceptor'),
    'enal':              ('very good', 'conjugated aldehyde'),
    'unsaturated ester':  ('good', 'α,β-unsaturated ester (acrylate)'),
    'unsaturated nitrile':('good', 'acrylonitrile type'),
    'unsaturated amide':  ('moderate', 'α,β-unsaturated amide'),
    'unsaturated ketone': ('good', 'general enone'),
    'vinyl sulfone':     ('good', '-SO₂Ph strong EWG'),
    'vinyl phosphonate': ('moderate', '-P(O)(OR)₂ moderate EWG'),
}

DONOR_TYPES = {
    'β-dicarbonyl':   ('excellent', 'pKa ~9-13, very stable enolate'),
    'malonate':        ('excellent', 'pKa ~13, doubly stabilized'),
    'acetoacetate':    ('very good', 'pKa ~11, β-keto ester'),
    'nitroalkane':     ('good', 'pKa ~10, nitro-stabilized'),
    'nitrile':         ('moderate', 'pKa ~25, needs stronger base'),
    'aldehyde/enolate':('moderate', 'standard enolate'),
    'ketone/enolate':  ('fair', 'pKa ~20, may compete with 1,2-addition'),
}


@ChemMCPManager.register_tool
class MichaelMechanism(BaseTool):
    __version__ = "0.1.0"
    name = "MichaelMechanism"
    func_name = 'explain_michael_mechanism'
    description = "Explain the Michael (conjugate) addition mechanism: nucleophilic 1,4-addition to α,β-unsaturated carbonyl compounds via enolate donors, competition with 1,2-direct addition, HSAB (hard/soft) analysis, and Robinson annulation as an extension. Covers donor reactivity, acceptor strength, and regioselectivity control."
    implementation_description = "Analyzes the Michael donor for enolizability and acidity, classifies the Michael acceptor by electrophilicity, provides complete stepwise mechanism for conjugate addition, explains the thermodynamic vs kinetic control of 1,4- vs 1,2-addition using HSAB theory, covers Robinson annulation (Michael + aldol), and evaluates reaction conditions."
    categories = ["Reaction"]
    tags = ["Michael Addition", "Conjugate Addition", "1,4-Addition", "Enolate", "HSAB", "Robinson Annulation", "α,β-Unsaturated"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('donor_smiles', 'str', 'N/A', 'SMILES of the Michael donor (active methylene compound with acidic α-H).'),
        ('acceptor_smiles', 'str', 'N/A', 'SMILES of the Michael acceptor (α,β-unsaturated carbonyl compound).'),
        ('base', 'str', 'base', 'Base for enolization. Options: OH-, RO-, piperidine, DBU, LDA, etc.'),
        ('solvent', 'str', 'polar aprotic', 'Solvent type.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: donor_smiles acceptor_smiles [base] [solvent]. E.g., "CC(=O)CC(=O)OEt CC=CO".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing donor_analysis, acceptor_analysis, mechanism_steps, hsab_analysis, product_prediction, competition_with_12_addition, robinson_context, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'donor_smiles': 'CC(=O)CC(=O)Oc1ccccc1',  # ethyl acetoacetate
                'acceptor_smiles': 'C=CC=O',  # propenal (acrolein)
                'base': 'base',
                'solvent': 'polar aprotic',
            },
            'text_input': {'query': 'ethyl_acetoacetate acrolein'},
            'output': {
                'result': {
                    'donor': 'ethyl acetoacetate (β-keto ester, excellent donor)',
                    'acceptor': 'acrolein (α,β-unsaturated aldehyde, very good acceptor)',
                    'mechanism_type': 'Michael addition (1,4-conjugate addition)',
                    'product': '1,5-dicarbonyl compound after protonation',
                    'key_feature': 'New C-C bond formed at β-position of acceptor',
                    'favorability': 'excellent — classic Michael system',
                }
            },
        },
        {
            'code_input': {
                'donor_smiles': 'CC(=O)C(C)=O',
                'acceptor_smiles': 'CC(=O)C=Cc1ccccc1',  # benzylideneacetone
                'base': 'base',
                'solvent': 'polar aprotic',
            },
            'text_input': {'query': 'acetone benzylideneacetone'},
            'output': {
                'result': {
                    'reaction_type': 'Robinson annulation precursor (Michael + intramolecular aldol)',
                    'donor_enolate': 'acetone → enamine or enolate',
                    'acceptor': 'benzylideneacetone (extended conjugation, good acceptor)',
                    'michael_product': '2,6-heptanedione derivative',
                    'subsequent_aldol': 'Intramolecular aldol → cyclohexenone (Robinson annulation)',
                    'favorability': 'good — Robinson annulation is a powerful ring-forming strategy',
                }
            },
        },
    ]

    def _run_base(self, donor_smiles: str, acceptor_smiles: str, base: str = 'base', solvent: str = 'polar aprotic') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(donor_smiles) or not is_smiles(acceptor_smiles):
            raise ChemMCPInputError("Invalid SMILES string(s).")

        mol_donor = Chem.MolFromSmiles(donor_smiles)
        mol_acceptor = Chem.MolFromSmiles(acceptor_smiles)
        if mol_donor is None or mol_acceptor is None:
            raise ChemMCPInputError("Cannot parse SMILES.")

        # 1. Donor analysis
        donor_analysis = self._analyze_donor(mol_donor)

        # 2. Acceptor analysis
        acceptor_analysis = self._analyze_acceptor(mol_acceptor)

        # 3. HSAB analysis
        hsab = self._hsab_analysis(donor_analysis, acceptor_analysis)

        # 4. Mechanism steps
        steps = self._build_mechanism(donor_analysis, acceptor_analysis)

        # 5. Product prediction
        products = self._predict_products(donor_analysis, acceptor_analysis)

        # 6. Competition with 1,2-addition
        competition = self._analyze_competition(donor_analysis, acceptor_analysis, base)

        # 7. Robinson context
        robinson = self._robinson_analysis(donor_analysis, acceptor_analysis)

        # 8. Favorability
        favorability = self._evaluate_favorability(donor_analysis, acceptor_analysis)

        result = {
            'result': {
                'donor_analysis': donor_analysis,
                'acceptor_analysis': acceptor_analysis,
                'hsab_analysis': hsab,
                'mechanism_steps': steps,
                'product_prediction': products,
                'competition_with_12_addition': competition,
                'robinson_annulation_context': robinson,
                'favorability': favorability,
                'summary': self._build_summary(donor_analysis, acceptor_analysis, products, favorability),
            }
        }

        logger.info(f"Michael: {donor_smiles} + {acceptor_smiles} → {favorability}")
        return result

    def _analyze_donor(self, mol):
        """Analyze Michael donor."""
        has_active_methylene = False
        active_positions = []
        n_carbonyls = 0

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() >= 2.0 and neighbor.GetAtomicNum() == 8:
                        n_carbonyls += 1

                        # Check alpha carbon
                        for alpha in atom.GetNeighbors():
                            if alpha.GetAtomicNum() == 6 and alpha.GetTotalNumHs() > 0:
                                has_active_methylene = True
                                active_positions.append({
                                    'idx': alpha.GetIdx(),
                                    'n_h': alpha.GetTotalNumHs(),
                                    'between_carbonyls': sum(
                                        1 for nn in alpha.GetNeighbors()
                                        if any(mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx()) and
                                               mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx()).GetBondTypeAsDouble() >= 2.0
                                               for nnn in nn.GetNeighbors() if nnn.GetAtomicNum() == 8)
                                        if nn.GetAtomicNum() == 6
                                    ),
                                })

        # Classify donor
        if n_carbonyls >= 2 and has_active_methylene:
            dtype = 'β-dicarbonyl'; drating = 'excellent'
        elif n_carbonyls >= 1 and has_active_methylene:
            dtype = 'monocarbonyl enolate'; drating = 'moderate'
        elif has_active_methylene:
            dtype = 'simple enolizable'; drating = 'fair'
        else:
            dtype = 'no active methylene'; drating = 'poor'

        return {
            'has_active_methylene': has_active_methylene,
            'active_positions': active_positions,
            'n_carbonyl_groups': n_carbonyls,
            'donor_type': dtype,
            'donor_rating': drating,
            'can_form_enolate': has_active_methylene,
        }

    def _analyze_acceptor(self, mol):
        """Analyze Michael acceptor."""
        has_conjugated_system = False
        acceptor_info = {'has_conjugation': False}

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() == 2.0 and neighbor.GetAtomicNum() == 6:
                        # Found C=C, check if conjugated to C=O
                        for c_neighbor in atom.GetNeighbors():
                            cb = mol.GetBondBetweenAtoms(atom.GetIdx(), c_neighbor.GetIdx())
                            if cb and cb.GetBondTypeAsDouble() >= 2.0 and c_neighbor.GetAtomicNum() == 8:
                                has_conjugated_system = True
                                # Check what's attached to the β-carbon
                                beta_subs = sum(1 for n in neighbor.GetNeighbors()
                                                if n.GetAtomicNum() != 1 and n.GetIdx() != atom.GetIdx())
                                acceptor_info = {
                                    'has_conjugation': True,
                                    'type': 'α,β-unsaturated carbonyl',
                                    'beta_substitution': beta_subs,
                                    'is_extended': self._check_extended_conjugation(mol, neighbor),
                                }
                                break
                        if has_conjugated_system: break
                if has_conjugated_system: break

        return acceptor_info

    def _check_extended_conjugation(self, mol, beta_atom):
        """Check for extended conjugation beyond the enone."""
        for neighbor in beta_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(beta_atom.GetIdx(), neighbor.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 2.0 and neighbor.GetAtomicNum() == 6:
                return True  # Extended (diene/aromatic)
        return False

    def _hsab_analysis(self, donor, acceptor):
        """Hard-Soft Acid-Base analysis."""
        d_type = donor.get('donor_type', '')
        a_has_cj = acceptor.get('has_conjugation', False)

        donor_char = 'Soft nucleophile' if 'dicarbonyl' in d_type else (
                      'Borderline nucleophile' if 'enolate' in d_type else 'Hard nucleophile')
        return {
            'donor_character': (
                f"{donor_char} — enolates are softer than alkoxides due to charge delocalization"
            ),
            'acceptor_character': (
                'Soft electrophile — the β-carbon of an α,β-unsaturated carbonyl is '
                'softer than the carbonyl carbon because it is more polarizable '
                '(positive charge delocalized over two atoms via resonance)'
            ),
            'prediction': (
                'Soft-soft match favors 1,4-(Michael) addition. '
                'Hard bases (RO⁻, R₂NH) favor 1,2-addition; soft bases (R₂CuLi, enolates) favor 1,4-addition.'
            ),
            'rule_of_thumb': (
                '1,2-addition is kinetically favored (irreversible, low T); '
                '1,4-addition is thermodynamically favored (reversible, often under equilibrating conditions)'
            ),
        }

    def _build_mechanism(self, donor, acceptor):
        """Build Michael addition mechanism."""
        return [
            {
                'step': 1,
                'name': 'Donor Deprotonation / Enolate Formation',
                'equation': 'Donor-H + Base → Enolate + BH',
                'details': (
                    f"Base deprotonates the active methylene position of the donor ({donor.get('donor_type','?')}) "
                    f"to form a resonance-stabilized enolate."
                ),
            },
            {
                'step': 2,
                'name': 'Conjugate (1,4-) Nucleophilic Attack',
                'equation': 'Enolate + Acceptor → New Enolate (at acceptor α-position)',
                'details': (
                    "The enolate attacks the **β-carbon** of the α,β-unsaturated acceptor (not the carbonyl carbon!). "
                    "This is a conjugate addition: electrons flow through the π system, forming a new C-C bond at the β-position "
                    "and placing negative charge on the carbonyl oxygen (enolate of the adduct)."
                ),
                'key_distinction': 'Attack at β-carbon (soft site) distinguishes Michael from direct (1,2-) addition.',
            },
            {
                'step': 3,
                'name': 'Protonation',
                'equation': 'Adduct Enolate + H⁺ → 1,5-Dicarbonyl Product',
                'details': (
                    "The enolate intermediate is protonated to give the neutral 1,5-dicarbonyl compound. "
                    "The product has carbonyl groups separated by 3 carbons (1,5-relationship) — a versatile synthetic handle."
                ),
            },
        ]

    def _predict_products(self, donor, acceptor):
        """Predict Michael product."""
        return {
            'product_type': '1,5-dicarbonyl compound',
            'structure_description': (
                'Donor fragment bonded to β-carbon of acceptor. '
                'Result: two carbonyl groups with a 3-carbon spacer between them.'
            ),
            'further_transformations': [
                'Aldol cyclization → cyclic enone (Robinson annulation)',
                'Intramolecular aldol → 5-7 membered rings',
                'Further Michael additions → polyannulation',
                'Reduction → various saturated products',
            ],
        }

    def _analyze_competition(self, donor, acceptor, base):
        """Analyze 1,4- vs 1,2- competition."""
        b_lower = base.lower()
        hard_bases = ['oh-', 'ro-', 'nh3', 'rn h2', 'li alh4']
        soft_bases = ['r2culi', 'enolate', 'phosphine', 'thiol', 'organocopper']

        is_hard = any(hb in b_lower for hb in hard_bases)
        is_soft = any(sb in b_lower for sb in soft_bases)

        if is_hard:
            prediction = '1,2-direct addition may compete or dominate'
            reason = 'Hard bases prefer the harder carbonyl carbon (direct addition).'
        elif is_soft or 'enolate' in donor.get('donor_type', ''):
            prediction = '1,4-Michael addition favored'
            reason = 'Soft enolate/copper reagents prefer the softer β-carbon (conjugate addition).'
        else:
            prediction = 'depends on conditions'
            reason = 'Under kinetic control (low T): 1,2 may dominate. Under thermodynamic control (reversible): 1,4 dominates.'

        return {
            'prediction': prediction,
            'reason': reason,
            'controlling_factor': 'HSAB character of nucleophile + reversibility of addition',
        }

    def _robinson_analysis(self, donor, acceptor):
        """Analyze Robinson annulation potential."""
        d_carbonyls = donor.get('n_carbonyl_groups', 0)
        a_cj = acceptor.get('has_conjugation', False)

        can_do_robinson = d_carbonyls >= 1 and a_cj

        return {
            'applicable': can_do_robinson,
            'description': (
                'Robinson annulation = Michael addition + intramolecular aldol condensation. '
                'Forms a 6-membered ring cyclohexenone — one of the most powerful ring-forming reactions in organic synthesis.'
            ) if can_do_robinson else 'Not a Robinson annulation system.',
            'steps_if_applicable': [
                'Step A: Michael addition gives 1,5-dicarbonyl',
                'Step B: Base-catalyzed enolization of ketone donor fragment',
                'Step C: Intramolecular aldol attack on remaining ketone/aldehyde',
                'Step D: Dehydration → α,β-unsaturated cyclic ketone (cyclohexenone derivative)',
            ] if can_do_robinson else [],
        }

    def _evaluate_favorability(self, donor, acceptor):
        score = 2
        if donor.get('can_form_enolate'): score += 3
        if donor.get('donor_rating') == 'excellent': score += 2
        elif donor.get('donor_rating') in ('very good', 'good'): score += 1
        if acceptor.get('has_conjugation'): score += 2
        if acceptor.get('is_extended'): score += 1  # extended conjugation = good acceptor

        if score >= 7: return 'excellent'
        elif score >= 5: return 'good'
        elif score >= 3: return 'moderate'
        return 'possible but check donor acidity and acceptor activation'

    def _build_summary(self, donor, acceptor, products, fav):
        d = donor.get('donor_type', '?')
        a = acceptor.get('type', '?')
        prod = products.get('product_type', '?')
        return f"Michael addition: {d} donor + {a} acceptor → {prod}. Favorability: {fav}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        donor = parts[0] if len(parts) > 0 else ''
        acc = parts[1] if len(parts) > 1 else ''
        b = parts[2] if len(parts) > 2 else 'base'
        s = parts[3] if len(parts) > 3 else 'polar aprotic'
        return self._run_base(donor, acc, b, s)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
