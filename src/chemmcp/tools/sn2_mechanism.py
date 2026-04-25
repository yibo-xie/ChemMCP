"""
SN2 Reaction Mechanism (Tool #117)
展示 SN2 反应的协同机理与过渡态。
Provides SN2 mechanism analysis with concerted backside attack and Walden inversion.
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
class Sn2Mechanism(BaseTool):
    __version__ = "0.1.0"
    name = "Sn2Mechanism"
    func_name = 'explain_sn2_mechanism'
    description = "Explain the SN2 nucleophilic substitution reaction mechanism: concerted backside attack, pentacoordinate transition state, and Walden inversion of configuration."
    implementation_description = "Analyzes substrate sterics to determine SN2 feasibility, identifies steric hindrance factors, evaluates nucleophile and solvent effects, and provides a complete one-step bimolecular mechanism description with stereochemical outcome."
    categories = ["Reaction"]
    tags = ["SN2", "Nucleophilic Substitution", "Walden Inversion", "Transition State", "Mechanism"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES of the substrate (alkyl halide or similar).'),
        ('nucleophile', 'str', 'OH-', 'Nucleophile formula or name (e.g., OH-, CN-, CH3O-, I-).'),
        ('solvent', 'str', 'polar aprotic', 'Solvent type: polar aprotic, polar protic, non-polar.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles nucleophile solvent. E.g., "CC(Cl) OH- polar_aprotic".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing mechanism_type, transition_state_analysis, steric_hindrance, rate_law, stereochemistry, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'CC(Cl)',
                'nucleophile': 'OH-',
                'solvent': 'polar aprotic',
            },
            'text_input': {'query': 'CC(Cl) OH- polar_aprotic'},
            'output': {
                'result': {
                    'substrate': 'ethyl chloride',
                    'mechanism_type': 'SN2',
                    'steric_hindrance': 'low (primary carbon)',
                    'transition_state': {'geometry': 'trigonal bipyramidal', 'nu_entry': 'axial (backside)', 'lg_exit': 'axial (front)'},
                    'rate_law': 'rate = k[substrate][nucleophile]',
                    'stereochemistry': 'Walden inversion (complete configuration reversal)',
                    'favorability': 'excellent',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'CC(C)(C)Cl',
                'nucleophile': 'OH-',
                'solvent': 'polar aprotic',
            },
            'text_input': {'query': 'CC(C)(C)Cl OH- polar_aprotic'},
            'output': {
                'result': {
                    'substrate': 'tert-butyl chloride',
                    'steric_hindrance': 'very high (tertiary carbon)',
                    'favorability': 'very poor — E2 elimination dominates',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, nucleophile: str = 'OH-', solvent: str = 'polar aprotic') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Analyze substrate for SN2 suitability
        sub_analysis = self._analyze_sn2_substrate(mol)

        # Build transition state analysis
        ts_analysis = self._build_transition_state(sub_analysis, nucleophile)

        # Evaluate favorability
        favorability = self._evaluate_sn2_favorability(sub_analysis, nucleophile, solvent)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'substrate_name': sub_analysis.get('common_name', 'unknown'),
                'mechanism_type': 'SN2 (concerted)',
                'steric_analysis': sub_analysis.get('steric', {}),
                'leaving_group': sub_analysis.get('leaving_group', {}),
                'transition_state': ts_analysis,
                'steps': [self._build_single_step(nucleophile, sub_analysis)],
                'rate_law': 'rate = k[substrate][nucleophile]',
                'stereochemistry': (
                    'Walden inversion: complete inversion of configuration at stereocenter. '
                    'Backside attack → umbrella flip. R → S or S → R.'
                ),
                'solvent_effect': self._get_solvent_effect(solvent),
                'nucleophile_strength': self._classify_nucleophile(nu=nucleophile),
                'competeting_reactions': self._identify_competing(sub_analysis),
                'favorability': favorability,
                'summary': self._build_summary(sub_analysis, favorability, nucleophile, solvent),
            }
        }

        logger.info(f"SN2 mechanism: {substrate_smiles} + {nucleophile} → {favorability}")
        return result

    def _analyze_sn2_substrate(self, mol):
        """Analyze substrate for SN2 suitability."""
        analysis = {}

        # Find leaving group and alpha-carbon
        lg_atom_idx, carbon_idx, lg_info = self._find_leaving_group(mol)
        analysis['leaving_group'] = lg_info

        if carbon_idx is None:
            analysis['suitable'] = False
            return analysis

        carbon = mol.GetAtomWithIdx(carbon_idx)
        neighbors = [n for n in carbon.GetNeighbors() if n.GetIdx() != lg_atom_idx]

        # Count substituents on alpha-carbon (excluding LG and H)
        n_non_h_subs = sum(1 for n in neighbors if n.GetAtomicNum() > 1)

        # Steric classification
        if n_non_h_subs == 0:
            steric_class = 'methyl'
            hindrance = 'none'
            sn2_score = 5
        elif n_non_h_subs == 1:
            steric_class = 'primary'
            hindrance = 'low'
            sn2_score = 4
        elif n_non_h_subs == 2:
            steric_class = 'secondary'
            hindrance = 'moderate'
            sn2_score = 2
        else:
            steric_class = 'tertiary'
            hindrance = 'severe'
            sn2_score = -3

        # Check beta-branching (additional steric hindrance)
        beta_branching = 0
        beta_groups = []
        for n in neighbors:
            if n.GetAtomicNum() == 6:
                for nn in n.GetNeighbors():
                    if nn.GetIdx() != carbon_idx and nn.GetAtomicNum() == 6:
                        beta_branching += 1
                        beta_groups.append(f"{nn.GetSymbol()}({nn.GetIdx()})")

        if beta_branching >= 2:
            hindrance += ' + significant β-branching'
            sn2_score -= 1

        # Check for resonance stabilization of LG (vinyl/aryl — bad for SN2)
        is_vinyl_aryl = self._check_vinyl_aryl(mol, carbon_idx)

        analysis['steric'] = {
            'alpha_carbon': carbon_idx,
            'carbon_class': steric_class,
            'n_non_h_substituents': n_non_h_subs,
            'hindrance_level': hindrance,
            'sn2_score': sn2_score,
            'beta_branching': beta_branching,
            'beta_groups': beta_groups,
            'is_vinyl_aryl': is_vinyl_aryl,
        }
        analysis['suitable'] = True
        return analysis

    def _find_leaving_group(self, mol):
        """Find leaving group."""
        lg_elements = {'Cl', 'Br', 'I', 'F'}
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in lg_elements:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        quality_map = {'I': 1, 'Br': 2, 'Cl': 3, 'F': 5}
                        return atom.GetIdx(), neighbor.GetIdx(), {
                            'atom': atom.GetSymbol(),
                            'name': f"{atom.GetSymbol()}ide",
                            'quality_rank': quality_map.get(atom.GetSymbol(), 99),
                        }
        return None, None, {'atom': None, 'name': 'Not found', 'quality_rank': 99}

    def _check_vinyl_aryl(self, mol, carbon_idx):
        """Check if alpha-carbon is part of vinyl or aryl system."""
        carbon = mol.GetAtomWithIdx(carbon_idx)
        for bond in carbon.GetBonds():
            if bond.GetBondTypeAsDouble() >= 2.0:
                return True  # double-bonded carbon (vinyl)
        if carbon.GetIsAromatic():
            return True
        return False

    def _build_transition_state(self, sub_analysis, nucleophile):
        """Describe the SN2 transition state."""
        steric = sub_analysis.get('steric', {})
        hindrance = steric.get('hindrance_level', '')

        nu_name = nucleophile

        return {
            'geometry': 'trigonal bipyramidal (TBP)',
            'nucleophile_position': 'axial — backside attack (180° from leaving group)',
            'leaving_group_position': 'axial — frontside departure',
            'equatorial_positions': 'occupied by the three substituents on α-carbon',
            'bond_order_partial': 'Nu---C partially formed (~50%), C-LG partially broken (~50%)',
            'charge_distribution': f'Partial negative on Nu and LG; partial positive on α-C',
            'energy_profile': 'Single high-energy barrier (no intermediate)',
            'steric_impact': f'Hindrance: {hindrance}. Higher hindrance raises TS energy dramatically.',
        }

    def _build_single_step(self, nucleophile, sub_analysis):
        """Build the single concerted step."""
        lg = sub_analysis.get('leaving_group', {})
        return {
            'step': 1,
            'name': 'Concerted Nucleophilic Attack & Leaving Group Departure',
            'details': (
                f"{nucleophile} attacks the α-carbon from the backside (180° from {lg.get('name', 'LG')}) "
                f"in a single concerted step. As the new Nu-C bond forms, the C-LG bond breaks simultaneously. "
                f"The transition state has trigonal bipyramidal geometry."
            ),
            'rate_determining': True,
            'reversible': False,
        }

    def _get_solvent_effect(self, solvent):
        s = solvent.lower().replace('_', ' ')
        effects = {
            'polar aprotic': 'Excellent — does NOT solvate nucleophile strongly. Best for SN2.',
            'polar protic': 'Poor — strongly solvates nucleophile via H-bonding, reducing reactivity.',
            'non-polar': 'Poor — cannot dissolve ionic nucleophiles.',
        }
        return effects.get(s, 'Unknown.')

    def _classify_nucleophile(self, nu):
        strong = ['CN-', 'I-', 'RS-', 'SH-', 'N3-', 'CH3O-', 'C2H5O-', 't-BuO-', 'NH2-', 'OH-']
        moderate = ['CH3COO-', 'F-', 'pyridine']
        weak = ['H2O', 'CH3OH', 'CH3CH2OH']
        if any(n in nu for n in strong): return 'strong nucleophile'
        if any(n in nu for n in moderate): return 'moderate nucleophile'
        if any(n in nu for n in weak): return 'weak nucleophile'
        return 'unknown'

    def _identify_competing(self, sub_analysis):
        """Identify competing reactions."""
        steric = sub_analysis.get('steric', {})
        c_class = steric.get('carbon_class', '')
        competing = []

        if c_class == 'tertiary':
            competing.append('E2 elimination (major pathway)')
            competing.append('SN1 (if good LG + protic solvent)')
        elif c_class == 'secondary':
            competing.append('E2 elimination (with strong base)')
            competing.append('SN1 (with protic solvent + weak Nu)')
        elif c_class in ('primary', 'methyl'):
            competing.append('E2 elimination (only with very strong, bulky base like t-BuOK)')
            if steric.get('beta_branching', 0) >= 2:
                competing.append('E2 favored over SN2 due to β-branching')

        return competing if competing else ['None dominant — SN2 should be the major pathway']

    def _evaluate_sn2_favorability(self, sub_analysis, nucleophile, solvent):
        steric = sub_analysis.get('steric', {})
        score = steric.get('sn2_score', 0)
        lg = sub_analysis.get('leaving_group', {})

        if lg.get('quality_rank', 99) <= 2: score += 1
        elif lg.get('quality_rank', 99) <= 4: score += 0

        nu_str = self._classify_nucleophile(nucleophile)
        if 'strong' in nu_str: score += 2
        elif 'moderate' in nu_str: score += 1

        solv = solvent.lower().replace('_', ' ')
        if 'aprotic' in solv: score += 2
        elif 'protic' in solv: score -= 1

        if steric.get('is_vinyl_aryl'): score -= 4

        if score >= 7: return 'excellent'
        elif score >= 5: return 'good'
        elif score >= 3: return 'moderate'
        elif score >= 1: return 'poor but possible'
        return 'very unlikely / use different mechanism'

    def _build_summary(self, sub_analysis, favorability, nucleophile, solvent):
        steric = sub_analysis.get('steric', {})
        parts = [
            f"SN2 reaction: {steric.get('carbon_class', '?')} substrate.",
            f"Steric hindrance: {steric.get('hindrance_level', '?')}.",
            f"Nucleophile: {nucleophile} ({self._classify_nucleophile(nucleophile)}).",
            f"Solvent: {solvent}.",
            f"Favorability: {favorability}.",
        ]
        comp = self._identify_competing(sub_analysis)
        if comp and comp[0] != 'None dominant':
            parts.append(f"Competing: {comp[0]}.")
        return " ".join(parts)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        nu = parts[1] if len(parts) > 1 else 'OH-'
        solv = parts[2] if len(parts) > 2 else 'polar aprotic'
        return self._run_base(sub, nu, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
