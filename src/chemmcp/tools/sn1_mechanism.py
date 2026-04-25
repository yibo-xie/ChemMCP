"""
SN1 Reaction Mechanism (Tool #116)
展示 SN1 反应的分步机理与中间体（碳正离子）。
Provides step-by-step SN1 mechanism analysis with carbocation intermediate.
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


# Common leaving groups and their quality (lower = better)
LEAVING_GROUPS = {
    'I': {'name': 'Iodide', 'quality': 1, 'pKa': -10},
    'Br': {'name': 'Bromide', 'quality': 2, 'pKa': -9},
    'Cl': {'name': 'Chloride', 'quality': 3, 'pKa': -7},
    'tosylate': {'name': 'Tosylate (OTs)', 'quality': 0, 'pKa': -2},
    'triflate': {'name': 'Triflate (OTf)', 'quality': -1, 'pKa': -14},
    '[F-]': {'name': 'Fluoride', 'quality': 6, 'pKa': 3.2},
    '[OH-]': {'name': 'Hydroxide', 'quality': 8, 'pKa': 15.7},
}

CARBOCATION_STABILITY = {
    'methyl': 1,
    'primary': 2,
    'secondary': 4,
    'tertiary': 8,
    'allylic': 7,
    'benzylic': 9,
    'resonance_stabilized': 8,
}


@ChemMCPManager.register_tool
class Sn1Mechanism(BaseTool):
    __version__ = "0.1.0"
    name = "Sn1Mechanism"
    func_name = 'explain_sn1_mechanism'
    description = "Explain the SN1 nucleophilic substitution reaction mechanism step-by-step: ionization to form carbocation intermediate, nucleophilic attack, and product formation."
    implementation_description = "Analyzes substrate structure to determine carbocation stability, identifies leaving groups, evaluates solvent and temperature effects, and provides a complete two-step SN1 mechanism with rate law, stereochemistry (racemization), and rearrangement possibilities."
    categories = ["Reaction"]
    tags = ["SN1", "Nucleophilic Substitution", "Carbocation", "Mechanism", "Organic Chemistry"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES of the substrate (alkyl halide or similar).'),
        ('nucleophile', 'str', 'H2O', 'Nucleophile formula or name (e.g., H2O, CH3OH, OH-, CN-).'),
        ('solvent', 'str', 'polar protic', 'Solvent type: polar protic, polar aprotic, non-polar.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles nucleophile solvent. E.g., "CC(C)(C)Cl H2O polar_protic".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing steps, carbocation_analysis, rate_law, stereochemistry, rearrangement_possibility, and full mechanism description.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'CC(C)(C)Cl',
                'nucleophile': 'H2O',
                'solvent': 'polar protic',
            },
            'text_input': {'query': 'CC(C)(C)Cl H2O polar_protic'},
            'output': {
                'result': {
                    'substrate': 'tert-butyl chloride',
                    'mechanism_type': 'SN1',
                    'n_steps': 2,
                    'carbocation': {'type': 'tertiary', 'stability': 'high', 'structure': '(CH3)3C+'},
                    'steps': [
                        {'step': 1, 'description': 'Ionization: C-Cl bond breaks heterolytically → carbocation + Cl-', 'rate_determining': True},
                        {'step': 2, 'description': 'Nucleophilic attack: H2O attacks carbocation → protonated alcohol', 'rate_determining': False},
                        {'step': 3, 'description': 'Deprotonation: Loss of H+ gives final alcohol product', 'rate_determining': False},
                    ],
                    'rate_law': 'rate = k[substrate]',
                    'stereochemistry': 'racemization (planar carbocation)',
                    'rearrangement_possible': False,
                    'favorability': 'excellent',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'CC(Cl)C',
                'nucleophile': 'CH3OH',
                'solvent': 'polar protic',
            },
            'text_input': {'query': 'CC(Cl)C CH3OH polar_protic'},
            'output': {
                'result': {
                    'substrate': 'isopropyl chloride (secondary)',
                    'carbocation': {'type': 'secondary', 'stability': 'moderate'},
                    'rearrangement_possible': True,
                    'favorability': 'good',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, nucleophile: str = 'H2O', solvent: str = 'polar protic') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Analyze substrate
        substrate_analysis = self._analyze_sn1_substrate(mol, substrate_smiles)

        # Build mechanism steps
        steps = self._build_sn1_steps(substrate_analysis, nucleophile, solvent)

        # Evaluate overall favorability
        favorability = self._evaluate_favorability(substrate_analysis, nucleophile, solvent)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'substrate_name': substrate_analysis.get('common_name', 'unknown'),
                'mechanism_type': 'SN1',
                'n_steps': len(steps),
                'carbocation_analysis': substrate_analysis.get('carbocation', {}),
                'leaving_group_analysis': substrate_analysis.get('leaving_group', {}),
                'steps': steps,
                'rate_law': 'rate = k[substrate]',
                'stereochemistry': 'Racemization via planar sp² carbocation intermediate; partial inversion if ion pair effects.',
                'rearrangement_possibility': substrate_analysis.get('can_rearrange', False),
                'rearrangement_details': substrate_analysis.get('rearrangement_info', ''),
                'solvent_effect': self._get_solvent_effect(solvent),
                'nucleophile_strength': self._classify_nucleophile(nucleophile),
                'favorability': favorability,
                'overall_summary': self._build_overall_summary(substrate_analysis, favorability, nucleophile, solvent),
            }
        }

        logger.info(f"SN1 mechanism: {substrate_smiles} + {nucleophile} → {favorability}")
        return result

    def _analyze_sn1_substrate(self, mol, smiles):
        """Analyze substrate for SN1 suitability."""
        analysis = {}

        # Identify carbon with leaving group
        lg_atom, carbon_idx, lg_info = self._find_leaving_group(mol)
        analysis['leaving_group'] = lg_info

        if carbon_idx is None:
            analysis['suitable'] = False
            analysis['reason'] = 'No suitable leaving group found.'
            return analysis

        carbon = mol.GetAtomWithIdx(carbon_idx)
        neighbors = carbon.GetNeighbors()

        # Determine carbocation type (primary/secondary/tertiary)
        n_carbon_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6 and n.GetIdx() != (lg_atom or -1))

        # Count all non-H, non-LG neighbors
        n_non_h = sum(1 for n in neighbors if n.GetAtomicNum() != 1 and n.GetIdx() != (lg_atom or -1))
        n_c_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6 and n.GetIdx() != (lg_atom or -1))

        if n_c_neighbors >= 2:
            cation_type = 'tertiary'
            stability = 'very high'
        elif n_c_neighbors == 1:
            cation_type = 'secondary'
            stability = 'moderate'
        else:
            cation_type = 'primary'
            stability = 'poor (unlikely SN1)'

        # Check for allylic/benzylic stabilization
        is_allylic = self._check_allylic(mol, carbon_idx)
        is_benzylic = self._check_benzylic(mol, carbon_idx)

        if is_allylic:
            cation_type += ' / allylic'
            stability = 'enhanced by resonance'
        if is_benzylic:
            cation_type += ' / benzylic'
            stability = 'greatly enhanced by resonance'

        analysis['carbon'] = {
            'index': carbon_idx,
            'symbol': carbon.GetSymbol(),
            'n_neighbors': n_non_h,
            'n_c_neighbors': n_c_neighbors,
            'cation_type': cation_type,
        }

        # Check for possible rearrangements (hydride shift, alkyl shift)
        can_rearrange, rearr_info = self._check_rearrangement(mol, carbon_idx, cation_type)
        analysis['can_rearrange'] = can_rearrange
        analysis['rearrangement_info'] = rearr_info

        analysis['carbocation'] = {
            'type': cation_type,
            'stability': stability,
            'is_allylic': is_allylic,
            'is_benzylic': is_benzylic,
        }

        analysis['suitable'] = True
        return analysis

    def _find_leaving_group(self, mol):
        """Find the best leaving group on the molecule."""
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in LEAVING_GROUPS:
                # Find attached carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        return atom.GetIdx(), neighbor.GetIdx(), {
                            'atom': sym,
                            'name': LEAVING_GROUPS[sym]['name'],
                            'quality_rank': LEAVING_GROUPS[sym]['quality'],
                            'pKa': LEAVING_GROUPS[sym]['pKa'],
                        }

            # Check for OTs/OTf groups (represented as O-S(=O)...)
            if sym == 'S':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        for o_neighbor in neighbor.GetNeighbors():
                            if o_neighbor.GetAtomicNum() == 6:
                                return atom.GetIdx(), o_neighbor.GetIdx(), {
                                    'atom': 'OSO2R (sulfonate ester)',
                                    'name': 'Tosylate/Triflate',
                                    'quality_rank': 0,
                                    'pKa': '< -2',
                                }

        return None, None, {'atom': None, 'name': 'Not found', 'quality_rank': 99, 'pKa': None}

    def _check_allylic(self, mol, carbon_idx):
        """Check if carbon is adjacent to a double bond."""
        carbon = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            for bond in neighbor.GetBonds():
                if bond.GetBondTypeAsDouble() == 2.0:
                    return True
        return False

    def _check_benzylic(self, mol, carbon_idx):
        """Check if carbon is attached to an aromatic ring."""
        carbon = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True
        return False

    def _check_rearrangement(self, mol, carbon_idx, cation_type):
        """Check if carbocation can rearrange to more stable form."""
        if 'tertiary' in cation_type or 'benzylic' in cation_type:
            return False, 'Already stable tertiary/benzylic carbocation.'

        carbon = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Check if this neighbor has more carbon substituents (hydride/alkyl shift source)
                n_nbr = neighbor.GetNeighbors()
                n_c_of_nbr = sum(1 for n in n_nbr if n.GetAtomicNum() == 6)
                if n_c_of_nbr >= 2:
                    return True, (
                        f"Possible hydride or alkyl shift from C({neighbor.GetIdx()}) "
                        f"(which has {n_c_of_nbr} carbon neighbors) could form a more stable carbocation."
                    )
        return False, ''

    def _build_sn1_steps(self, analysis, nucleophile, solvent):
        """Build step-by-step mechanism."""
        lg = analysis.get('leaving_group', {})
        cat = analysis.get('carbocation', {})
        steps = []

        # Step 1: Ionization (rate-determining)
        steps.append({
            'step': 1,
            'name': 'Ionization (Rate-Determining Step)',
            'equation': f'R-LG → R⁺ + LG⁻',
            'details': (
                f"The {lg.get('name', 'leaving group')} departs with the electron pair, "
                f"forming a {cat.get('type', 'carbocation')}. "
                f"This is slow and reversible. Solvent ({solvent}) stabilizes both ions."
            ),
            'rate_determining': True,
            'reversible': True,
        })

        # Step 2: Nucleophilic attack
        nu_name = self._get_nucleophile_name(nucleophile)
        steps.append({
            'step': 2,
            'name': 'Nucleophilic Attack',
            'equation': f'R⁺ + {nu_name} → R-{nu_name[:2] if len(nu_name) > 2 else nu_name}⁺ (or neutral)',
            'details': (
                f"{nu_name} attacks the planar carbocation from either face. "
                f"Fast step — not rate-determining."
            ),
            'rate_determining': False,
            'reversible': False,
        })

        # Step 3: Deprotonation (if needed)
        if nucleophile.upper() in ('H2O', 'CH3OH', 'ROH', 'HO'):
            steps.append({
                'step': 3,
                'name': 'Deprotonation',
                'equation': 'R-OH₂⁺ + B → R-OH + BH⁺',
                'details': 'A base removes the extra proton from the oxonium/carbocation intermediate.',
                'rate_determining': False,
                'reversible': False,
            })

        return steps

    def _get_solvent_effect(self, solvent):
        s = solvent.lower()
        effects = {
            'polar protic': 'Excellent — stabilizes ions, promotes ionization. Best for SN1.',
            'polar aprotic': 'Moderate — stabilizes cations but not as well as protic solvents.',
            'non-polar': 'Poor — cannot stabilize ions. SN1 very unlikely.',
        }
        return effects.get(s, 'Unknown solvent effect.')

    def _classify_nucleophile(self, nu):
        strong_nu = ['CN-', 'I-', 'RS-', 'NH2-', 'OH-', 'CH3O-', 'C2H5O-', 'N3-', 'SH-']
        weak_nu = ['H2O', 'CH3OH', 'CH3CH2OH', 'ROH', 'carboxylic acids']
        nu_upper = nu.upper().replace(' ', '')
        if any(n in nu_upper for n in strong_nu):
            return 'strong nucleophile'
        elif any(n in nu_upper for n in weak_nu):
            return 'weak nucleophile'
        return 'moderate nucleophile'

    def _get_nucleophile_name(self, nu):
        names = {
            'H2O': 'water (H₂O)', 'CH3OH': 'methanol (CH₃OH)', 'OH-': 'hydroxide (OH⁻)',
            'CN-': 'cyanide (CN⁻)', 'CH3CH2OH': 'ethanol (CH₃CH₂OH)',
            'CH3COO-': 'acetate (CH₃COO⁻)', 'I-': 'iodide (I⁻)',
        }
        return names.get(nu, nu)

    def _evaluate_favorability(self, analysis, nucleophile, solvent):
        cat = analysis.get('carbocation', {})
        lg = analysis.get('leaving_group', {})
        cat_type = cat.get('type', '')

        score = 0
        if 'tertiary' in cat_type: score += 3
        elif 'secondary' in cat_type: score += 2
        elif 'primary' in cat_type: score -= 2
        if 'benzylic' in cat_type: score += 2
        if 'allylic' in cat_type: score += 1
        if lg.get('quality_rank', 99) <= 2: score += 2
        elif lg.get('quality_rank', 99) <= 4: score += 1
        if 'polar protic' in solvent.lower(): score += 2
        elif 'polar aprotic' in solvent.lower(): score += 0
        else: score -= 1

        if score >= 6: return 'excellent'
        elif score >= 4: return 'good'
        elif score >= 2: return 'moderate (possible but may compete with E1/SN2)'
        elif score >= 0: return 'poor (slow, likely needs forcing conditions)'
        return 'very unlikely'

    def _build_overall_summary(self, analysis, favorability, nucleophile, solvent):
        cat = analysis.get('carbocation', {})
        lg = analysis.get('leaving_group', {})
        parts = [
            f"SN1 reaction of {analysis.get('common_name', 'substrate')} with {nucleophile}.",
            f"Carbocation: {cat.get('type', '?')} ({cat.get('stability', '?')}).",
            f"Leaving group: {lg.get('name', '?')}.",
            f"Favorability: {favorability}.",
        ]
        if analysis.get('can_rearrange'):
            parts.append(f"⚠ Rearrangement expected: {analysis.get('rearrangement_info', '')}")
        return " ".join(parts)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        nu = parts[1] if len(parts) > 1 else 'H2O'
        solv = parts[2] if len(parts) > 2 else 'polar protic'
        return self._run_base(sub, nu, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
