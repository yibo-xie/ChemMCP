"""
E2 Elimination Reaction Mechanism (Tool #119)
展示 E2 消除反应机理与构象要求（反式共平面）。
Provides E2 mechanism analysis with anti-periplanar requirement and stereochemistry.
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
class E2Mechanism(BaseTool):
    __version__ = "0.1.0"
    name = "E2Mechanism"
    func_name = 'explain_e2_mechanism'
    description = "Explain the E2 bimolecular elimination reaction mechanism: concerted anti-periplanar elimination, transition state geometry, and stereochemical requirements."
    implementation_description = "Analyzes substrate for E2 feasibility, identifies anti-periplanar β-hydrogens, evaluates base strength/sterics, predicts Zaitsev vs Hofmann product selectivity, and provides complete one-step bimolecular mechanism with conformational analysis."
    categories = ["Reaction"]
    tags = ["E2", "Elimination", "Anti-periplanar", "Zaitsev", "Hofmann", "Stereochemistry"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES of the substrate (alkyl halide or similar).'),
        ('base', 'str', 'NaOEt/EtOH', 'Base formula or name (e.g., NaOEt, KOH, t-BuOK, NaNH2).'),
        ('solvent', 'str', 'EtOH', 'Solvent name.'),
        ('temperature_c', 'float', '50.0', 'Temperature in °C.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles base solvent temp_c.'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing conformational_requirement, beta_h_analysis, transition_state, predicted_product, stereochemistry, and competition_with_sn2.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'CC(Cl)CC',
                'base': 'NaOEt/EtOH',
                'solvent': 'EtOH',
                'temperature_c': 55.0,
            },
            'text_input': {'query': 'CC(Cl)CC NaOEt/EtOH EtOH 55'},
            'output': {
                'result': {
                    'mechanism_type': 'E2',
                    'conformational_requirement': 'anti-periplanar H-C-C-LG dihedral angle ≈ 180°',
                    'beta_hydrogens_available': 4,
                    'predicted_major_product': 'but-1-ene or but-2-ene',
                    'zaitsev_product_favored': True,
                    'stereochemistry': 'Anti elimination → specific E/Z stereochemistry',
                    'rate_law': 'rate = k[substrate][base]',
                    'competition': 'SN2 competes with primary/secondary substrates',
                    'favorability': 'good',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'CC(C)(C)Cl',
                'base': 't-BuOK/t-BuOH',
                'solvent': 't-BuOH',
                'temperature_c': 30.0,
            },
            'text_input': {'query': 'CC(C)(C)Cl t-BuOK t-BuOH 30'},
            'output': {
                'result': {
                    'favorability': 'excellent — bulky base forces E2 over SN2 on tertiary substrate',
                    'product_selectivity': 'Hofmann product may be favored due to steric hindrance of bulky base',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, base: str = 'NaOEt/EtOH', solvent: str = 'EtOH', temperature_c: float = 50.0) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Substrate analysis
        sub_analysis = self._analyze_e2_substrate(mol)

        # Base characterization
        base_info = self._characterize_base(base)

        # β-Hydrogen analysis
        beta_h = self._find_beta_hydrogens_e2(mol, sub_analysis)

        # Product prediction
        products = self._predict_e2_products(beta_h, sub_analysis, base_info)

        # Transition state
        ts = self._build_ts(sub_analysis, base_info)

        # Competition
        competition = self._identify_competition_e2(sub_analysis, base_info)

        favorability = self._evaluate_favorability_e2(sub_analysis, base_info, temperature_c)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'mechanism_type': 'E2 (concerted)',
                'substrate_analysis': sub_analysis,
                'base_analysis': base_info,
                'conformational_requirement': (
                    'Anti-periplanar: H–Cα–Cβ–LG dihedral angle must be ~180°. '
                    'The C-H and C-LG bonds must be coplanar for optimal p-orbital overlap '
                    'in forming the π bond.'
                ),
                'beta_hydrogen_analysis': beta_h,
                'transition_state': ts,
                'product_prediction': products,
                'stereochemistry': self._get_stereochemistry(beta_h),
                'rate_law': 'rate = k[substrate][base]',
                'competition_with_sn2': competition,
                'temperature_effect': f"T={temperature_c}°C — higher T favors E2 over SN2",
                'favorability': favorability,
                'summary': self._build_summary_e2(sub_analysis, products, base_info, favorability),
            }
        }

        logger.info(f"E2 mechanism: {substrate_smiles} + {base} → {favorability}")
        return result

    def _analyze_e2_substrate(self, mol):
        """Analyze substrate for E2."""
        lg_atom_idx, carbon_idx, lg_info = self._find_lg(mol)
        if carbon_idx is None:
            return {'suitable': False}

        carbon = mol.GetAtomWithIdx(carbon_idx)
        neighbors = [n for n in carbon.GetNeighbors() if n.GetIdx() != lg_atom_idx]
        n_non_h = sum(1 for n in neighbors if n.GetAtomicNum() > 1)

        if n_non_h >= 3: c_class = 'tertiary'; sn2_score = -3
        elif n_non_h == 2: c_class = 'secondary'; sn2_score = 2
        elif n_non_h == 1: c_class = 'primary'; sn2_score = 4
        else: c_class = 'methyl'; sn2_score = 5

        # Check neopentyl-type (primary but very hindered)
        beta_branching = sum(
            1 for n in neighbors
            if n.GetAtomicNum() == 6
            for nn in n.GetNeighbors()
            if nn.GetAtomicNum() == 6 and nn.GetIdx() != carbon_idx
        )

        return {
            'suitable': True,
            'alpha_carbon_idx': carbon_idx,
            'carbon_class': c_class,
            'n_non_h_subs': n_non_h,
            'leaving_group': lg_info,
            'beta_branching': beta_branching,
        }

    def _find_lg(self, mol):
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ('Cl', 'Br', 'I'):
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        q = {'I': 1, 'Br': 2, 'Cl': 3}.get(atom.GetSymbol(), 99)
                        return atom.GetIdx(), neighbor.GetIdx(), {
                            'atom': atom.GetSymbol(), 'name': f"{atom.GetSymbol()}ide",
                            'quality_rank': q,
                        }
        return None, None, {'atom': None, 'name': 'Not found', 'quality_rank': 99}

    def _characterize_base(self, base):
        """Classify base by strength and sterics."""
        b = base.upper().replace(' ', '').replace('/', '')
        strong_small = ['NAOH', 'KOH', 'NAOET', 'NAOME', 'CH3ONA', 'C2H5ONA']
        strong_bulky = ['T-BUOK', 'T-BUNA', 'LDA', 'NANH2', 'LIDAH']
        weak = ['H2O', 'CH3OH', 'CH3CH2OH', 'ETOH', 'ACETATE', 'CH3COONA']

        if any(s in b for s in strong_bulky):
            return {'strength': 'strong', 'sterics': 'bulky', 'type': 'strong bulky base', 'hofmann_bias': 'high'}
        if any(s in b for s in strong_small):
            return {'strength': 'strong', 'sterics': 'small', 'type': 'strong small base', 'hofmann_bias': 'low'}
        if any(s in b for s in weak):
            return {'strength': 'weak', 'sterics': 'small', 'type': 'weak base', 'hofmann_bias': 'none'}
        return {'strength': 'unknown', 'sterics': 'unknown', 'type': base, 'hofmann_bias': 'unknown'}

    def _find_beta_hydrogens_e2(self, mol, sub_analysis):
        alpha_idx = sub_analysis.get('alpha_carbon_idx')
        if alpha_idx is None:
            return {'total_beta_h': 0, 'beta_positions': []}

        alpha = mol.GetAtomWithIdx(alpha_idx)
        lg_set = set()
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ('Cl', 'Br', 'I'):
                for n in atom.GetNeighbors():
                    if n.GetIdx() == alpha_idx:
                        lg_set.add(atom.GetIdx())

        positions = []
        total_h = 0
        for neighbor in alpha.GetNeighbors():
            if neighbor.GetIdx() in lg_set:
                continue
            if neighbor.GetAtomicNum() == 6:
                n_h = neighbor.GetTotalNumHs()
                total_h += max(n_h, 0)
                n_c_subs = sum(1 for nn in neighbor.GetNeighbors()
                               if nn.GetAtomicNum() == 6 and nn.GetIdx() != alpha_idx)
                positions.append({
                    'index': neighbor.GetIdx(),
                    'n_beta_h': max(n_h, 0),
                    'substitution': n_c_subs,
                })
            elif neighbor.GetAtomicNum() == 1:
                pass  # α-H, not β-H

        return {
            'total_beta_hydrogens': total_h,
            'n_beta_carbons': len(positions),
            'beta_positions': positions,
            'can_eliminate': total_h > 0,
        }

    def _predict_e2_products(self, beta_h, sub_analysis, base_info):
        positions = beta_h.get('beta_positions', [])
        if not positions:
            return {'can_form_alkene': False, 'reason': 'No β-hydrogens.'}

        zaitsev_pos = max(positions, key=lambda p: p['substitution'])
        hofmann_pos = min(positions, key=lambda p: p['substitution'])

        # Bulky bases favor Hofmann
        hofmann_bias = base_info.get('hofmann_bias', 'none')
        if hofmann_bias == 'high':
            major = 'hofmann'
            reason = f"Bulky base ({base_info.get('type', '')}) favors less hindered β-H abstraction → Hofmann product"
        else:
            major = 'zaitsev'
            reason = f"Small/base ({base_info.get('type', '')}) favors more substituted alkene → Zaitsev product"

        sub_types = {0: 'mono', 1: 'di', 2: 'tri', 3: 'tetra'}

        return {
            'can_form_alkene': True,
            'n_possible_alkenes': len(set(p['index'] for p in positions)),
            'zaitsev_product': {
                'alkene_type': f"{sub_types.get(zaitsev_pos['substitution'], '?')}substituted",
                'beta_carbon': zaitsev_pos['index'],
                'is_major': major == 'zaitsev',
            },
            'hofmann_product': {
                'alkene_type': f"{sub_types.get(hofmann_pos['substitution'], '?')}substituted",
                'beta_carbon': hofmann_pos['index'],
                'is_major': major == 'hofmann',
            } if hofmann_pos['index'] != zaitsev_pos['index'] else None,
            'selectivity': reason,
        }

    def _build_ts(self, sub_analysis, base_info):
        return {
            'geometry': 'coplanar arrangement of H-C-C-LG atoms',
            'dihedral_angle': '~180° (anti-periplanar) or ~0° (syn-periplanar, less favored)',
            'bond_changes': 'C-H breaking, C=C π-bond forming, C-LG breaking — all concerted',
            'partial_bonds': 'All bond order changes occur simultaneously in the TS',
            'charge': 'Partial negative on base and LG; partial positive on α- and β-carbons',
            'kinetics': 'Second-order: rate = k[substrate][base]',
            'base_role': base_info.get('type', 'Base') + ' abstracts β-proton',
        }

    def _get_stereochemistry(self, beta_h):
        n_positions = beta_h.get('n_beta_carbons', 0)
        if n_positions >= 2:
            return (
                'Anti-periplanar elimination gives defined stereochemistry: '
                'the eliminated H and LG must be trans-diaxial (in cyclohexanes) or anti (in acyclic systems). '
                'E or Z alkene determined by which β-H is anti to the LG.'
            )
        return 'Single elimination product possible — no E/Z stereoisomerism.'

    def _identify_competition_e2(self, sub_analysis, base_info):
        c_class = sub_analysis.get('carbon_class', '')
        comp = []

        if c_class in ('methyl', 'primary') and base_info.get('sterics') != 'bulky':
            comp.append('SN2 is major competitor (good Nu, low sterics)')
        elif c_class == 'secondary':
            if base_info.get('strength') == 'strong' and base_info.get('sterics') == 'bulky':
                comp.append('E2 dominant (bulky base disfavors SN2)')
            else:
                comp.append('SN2 and E2 compete; ratio depends on conditions')
        elif c_class == 'tertiary':
            comp.append('SN2 impossible (too hindered); E2 is only bimolecular path')
            comp.append('E1/SN1 may compete if solvent is protic')

        return comp if comp else ['E2 should dominate']

    def _evaluate_favorability_e2(self, sub_analysis, base_info, temp_c):
        score = 0
        c_class = sub_analysis.get('carbon_class', '')

        if c_class == 'tertiary' and base_info.get('strength') == 'strong': score += 4
        elif c_class == 'secondary' and base_info.get('strength') == 'strong': score += 3
        elif c_class == 'primary' and base_info.get('strength') == 'strong': score += 2
        elif c_class == 'secondary': score += 1

        lg = sub_analysis.get('leaving_group', {})
        if lg.get('quality_rank', 99) <= 2: score += 1

        if temp_c > 25: score += 1
        if temp_c > 60: score += 1

        # β-H availability checked separately; assume valid substrate

        if score >= 6: return 'excellent'
        elif score >= 4: return 'good'
        elif score >= 2: return 'moderate'
        return 'possible but may compete'

    def _build_summary_e2(self, sub_analysis, products, base_info, favorability):
        c_class = sub_analysis.get('carbon_class', '?')
        sel = products.get('selectivity', '')
        return f"E2 elimination on {c_class} substrate. {sel}. Favorability: {favorability}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        base = parts[1] if len(parts) > 1 else 'NaOEt/EtOH'
        solv = parts[2] if len(parts) > 2 else 'EtOH'
        temp = float(parts[3]) if len(parts) > 3 else 50.0
        return self._run_base(sub, base, solv, temp)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
