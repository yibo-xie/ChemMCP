"""
E1 Elimination Reaction Mechanism (Tool #118)
展示 E1 消除反应机理：离子化→消除。
Provides E1 mechanism analysis with carbocation intermediate and elimination step.
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
class E1Mechanism(BaseTool):
    __version__ = "0.1.0"
    name = "E1Mechanism"
    func_name = 'explain_e1_mechanism'
    description = "Explain the E1 unimolecular elimination reaction mechanism: ionization to form carbocation intermediate, then base-induced β-hydrogen elimination to form alkene."
    implementation_description = "Analyzes substrate for E1 suitability (carbocation stability), identifies β-hydrogens available for elimination, predicts major product via Zaitsev's rule, evaluates competition with SN1, and provides complete two-step mechanism with regiochemistry and stereochemistry."
    categories = ["Reaction"]
    tags = ["E1", "Elimination", "Carbocation", "Zaitsev Rule", "Alkene Formation"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES of the substrate (alkyl halide or similar).'),
        ('base', 'str', 'H2O', 'Base formula or name (e.g., H2O, CH3OH, EtOH).'),
        ('solvent', 'str', 'polar protic', 'Solvent type.'),
        ('temperature_c', 'float', '50.0', 'Temperature in °C (higher T favors E1 over SN1).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles base solvent temp_c.'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing steps, carbocation_analysis, beta_hydrogen_analysis, predicted_product, zaitsev_analysis, and competition_with_sn1.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'CC(C)(C)Cl',
                'base': 'H2O',
                'solvent': 'polar protic',
                'temperature_c': 60.0,
            },
            'text_input': {'query': 'CC(C)(C)Cl H2O polar_protic 60'},
            'output': {
                'result': {
                    'mechanism_type': 'E1',
                    'n_steps': 2,
                    'carbocation': {'type': 'tertiary', 'stability': 'very high'},
                    'beta_hydrogens': [{'position': ..., 'count': 9}],
                    'predicted_major_product': 'isobutylene (2-methylpropene)',
                    'zaitsev_product': 'more substituted alkene (trisubstituted)',
                    'rate_law': 'rate = k[substrate]',
                    'competition': 'SN1 competes significantly; E1 favored at higher T',
                    'favorability': 'excellent',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, base: str = 'H2O', solvent: str = 'polar protic', temperature_c: float = 50.0) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Analyze substrate
        sub_analysis = self._analyze_e1_substrate(mol)

        # Find β-hydrogens
        beta_h_analysis = self._find_beta_hydrogens(mol, sub_analysis)

        # Predict products via Zaitsev/Hofmann analysis
        product_analysis = self._predict_products(beta_h_analysis, sub_analysis)

        # Build mechanism steps
        steps = self._build_e1_steps(sub_analysis, base, temperature_c)

        # Competition with SN1
        competition = self._analyze_sn1_competition(sub_analysis, base, temperature_c)

        favorability = self._evaluate_favorability(sub_analysis, temperature_c)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'mechanism_type': 'E1',
                'n_steps': len(steps),
                'carbocation_analysis': sub_analysis.get('carbocation', {}),
                'beta_hydrogen_analysis': beta_h_analysis,
                'product_prediction': product_analysis,
                'steps': steps,
                'rate_law': 'rate = k[substrate]',
                'stereochemistry': (
                    'E/Z mixture possible for unsymmetrical alkenes. '
                    'Anti-periplanar elimination preferred but not required (carbocation intermediate allows rotation).'
                ),
                'competition_with_sn1': competition,
                'temperature_effect': f"T={temperature_c}°C — higher temperature favors E1 (elimination) over SN1 (substitution) due to positive ΔS‡.",
                'favorability': favorability,
                'summary': self._build_summary(sub_analysis, product_analysis, favorability),
            }
        }

        logger.info(f"E1 mechanism: {substrate_smiles} → {favorability}")
        return result

    def _analyze_e1_substrate(self, mol):
        """Analyze substrate for E1 suitability."""
        analysis = {}
        lg_atom_idx, carbon_idx, lg_info = self._find_leaving_group(mol)
        analysis['leaving_group'] = lg_info

        if carbon_idx is None:
            analysis['suitable'] = False
            return analysis

        carbon = mol.GetAtomWithIdx(carbon_idx)
        neighbors = [n for n in carbon.GetNeighbors() if n.GetIdx() != lg_atom_idx]
        n_c_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6)

        if n_c_neighbors >= 3:
            cation_type = 'tertiary'; stability = 'very high'
        elif n_c_neighbors == 2:
            cation_type = 'secondary'; stability = 'moderate'
        elif n_c_neighbors == 1:
            cation_type = 'primary'; stability = 'poor (unlikely)'
        else:
            cation_type = 'methyl'; stability = 'none'

        # Check allylic/benzylic
        is_allylic = any(b.GetBondTypeAsDouble() == 2.0 for n in neighbors for b in n.GetBonds())
        is_benzylic = any(n.GetIsAromatic() for n in neighbors)
        if is_allylic: cation_type += '/allylic'; stability += ' + resonance'
        if is_benzylic: cation_type += '/benzylic'; stability += ' + resonance'

        analysis['carbon'] = {'index': carbon_idx, 'cation_type': cation_type}
        analysis['carbocation'] = {'type': cation_type, 'stability': stability}
        analysis['suitable'] = True
        return analysis

    def _find_leaving_group(self, mol):
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

    def _find_beta_hydrogens(self, mol, sub_analysis):
        """Find all β-carbons and their hydrogens."""
        carbon_idx = sub_analysis.get('carbon', {}).get('index')
        if carbon_idx is None:
            return {'total_beta_h': 0, 'beta_positions': []}

        alpha_carbon = mol.GetAtomWithIdx(carbon_idx)
        lg_atoms = set()
        for a in ['Cl', 'Br', 'I']:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == a:
                    for n in atom.GetNeighbors():
                        if n.GetIdx() == carbon_idx:
                            lg_atoms.add(atom.GetIdx())

        beta_carbons = []
        total_h = 0
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetIdx() in lg_atoms:
                continue
            if neighbor.GetAtomicNum() == 6:  # β-carbon
                n_h = neighbor.GetTotalNumHs()
                # Also count implicit H
                n_h = max(n_h, neighbor.GetTotalNumHs())
                total_h += n_h
                beta_carbons.append({
                    'index': neighbor.GetIdx(),
                    'symbol': 'C',
                    'n_beta_hydrogens': n_h,
                    'substitution_level': sum(
                        1 for nn in neighbor.GetNeighbors()
                        if nn.GetAtomicNum() == 6 and nn.GetIdx() != carbon_idx
                    ),
                })
            elif neighbor.GetAtomicNum() == 1:  # H directly on α-C (no β-H here)
                pass

        return {
            'total_beta_hydrogens': total_h,
            'n_beta_carbons': len(beta_carbons),
            'beta_positions': beta_carbons,
            'can_eliminate': total_h > 0,
        }

    def _predict_products(self, beta_h_analysis, sub_analysis):
        """Predict alkene products using Zaitsev's rule."""
        positions = beta_h_analysis.get('beta_positions', [])
        if not positions:
            return {'can_form_alkene': False, 'reason': 'No β-hydrogens found.'}

        # Find most substituted alkene (Zaitsev product)
        best_pos = max(positions, key=lambda p: p['substitution_level'])
        zaitsev_sub = best_pos['substitution_level']

        # Hofmann product (least substituted)
        worst_pos = min(positions, key=lambda p: p['substitution_level'])

        alkene_types = {
            0: 'monosubstituted',
            1: 'disubstituted (terminal)',
            2: 'trisubstituted (internal)',
            3: 'tetrasubstituted',
        }

        return {
            'can_form_alkene': True,
            'n_possible_alkenes': len(positions),
            'zaitsev_product': {
                'description': f'Most substituted alkene ({alkene_types.get(zaitsev_sub, "?")})',
                'beta_carbon_index': best_pos['index'],
                'substitution_level': zaitsev_sub,
                'is_major_product': True,
            },
            'hofmann_product': {
                'description': f'Least substituted alkene ({alkene_types.get(worst_pos["substitution_level"], "?")})',
                'beta_carbon_index': worst_pos['index'],
                'substitution_level': worst_pos['substitution_level'],
                'is_minor_product': True,
            } if worst_pos['index'] != best_pos['index'] else None,
            'regioselectivity': 'Zaitsev product favored' if zaitsev_sub >= 2 else 'May give mixture',
        }

    def _build_e1_steps(self, sub_analysis, base, temp_c):
        lg = sub_analysis.get('leaving_group', {})
        cat = sub_analysis.get('carbocation', {})
        return [
            {
                'step': 1,
                'name': 'Ionization (Rate-Determining Step)',
                'equation': 'R-LG → R⁺ + LG⁻',
                'details': (
                    f"Leaving group ({lg.get('name', '?')}) departs forming "
                    f"{cat.get('type', '?')} carbocation. Slow, unimolecular."
                ),
                'rate_determining': True,
            },
            {
                'step': 2,
                'name': 'Base-Mediated β-H Elimination',
                'equation': 'R⁺ + :B → Alkene + BH⁺',
                'details': (
                    f"Base ({base}) abstracts a β-proton, electron pair forms π bond. "
                    f"Fast step after carbocation formation."
                ),
                'rate_determining': False,
            },
        ]

    def _analyze_sn1_competition(self, sub_analysis, base, temp_c):
        cat = sub_analysis.get('carbocation', {})
        cat_type = cat.get('type', '')

        strong_bases = ['OH-', 'EtO-', 't-BuO-', 'NH2-', 'CH3O-']
        weak_bases = ['H2O', 'CH3OH', 'EtOH', 'AcO-']

        base_is_strong = any(b in base.upper().replace(' ', '') for b in strong_bases)
        base_is_weak = any(b in base.upper().replace(' ', '') for b in weak_bases)

        comp = []
        if 'tertiary' in cat_type:
            if base_is_weak or base.lower() in ('h2o', 'ch3oh'):
                comp.append('SN1:E1 ≈ 3:1 at moderate T (weak Nu favors substitution)')
            elif base_is_strong:
                comp.append('E2 may compete strongly (strong base)')
        elif 'secondary' in cat_type:
            comp.append('SN1, E1, SN2, E2 all possible — conditions determine outcome')

        if temp_c > 80:
            comp.append('High T favors E1 (elimination has larger +ΔS‡)')
        elif temp_c < 25:
            comp.append('Low T favors SN1 (substitution has lower ΔH‡)')

        return comp if comp else ['Both SN1 and E1 proceed through same carbocation intermediate']

    def _evaluate_favorability(self, sub_analysis, temp_c):
        cat = sub_analysis.get('carbocation', {})
        score = 0
        ct = cat.get('type', '')
        if 'tertiary' in ct: score += 3
        elif 'secondary' in ct: score += 2
        if 'benzylic' in ct: score += 2
        if 'allylic' in ct: score += 1
        if temp_c > 50: score += 1
        if temp_c > 80: score += 1
        lg = sub_analysis.get('leaving_group', {})
        if lg.get('quality_rank', 99) <= 2: score += 1

        if score >= 5: return 'excellent'
        elif score >= 3: return 'good'
        elif score >= 2: return 'moderate'
        return 'poor'

    def _build_summary(self, sub_analysis, product_analysis, favorability):
        cat = sub_analysis.get('carbocation', {})
        prod = product_analysis.get('zaitsev_product', {})
        return (
            f"E1 elimination: {cat.get('type', '?')} carbocation → "
            f"{prod.get('description', 'alkene')}. Favorability: {favorability}."
        )

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        base = parts[1] if len(parts) > 1 else 'H2O'
        solv = parts[2] if len(parts) > 2 else 'polar protic'
        temp = float(parts[3]) if len(parts) > 3 else 50.0
        return self._run_base(sub, base, solv, temp)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
