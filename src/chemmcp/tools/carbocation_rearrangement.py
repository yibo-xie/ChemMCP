"""
Carbocation Rearrangement (Tool #125)
碳正离子重排：氢迁移（hydride shift）、烷基迁移（alkyl shift）、
Wagner-Meerwein 重排、片呐醇重排、扩环/缩环。
Provides carbocation rearrangement analysis: hydride shifts, alkyl shifts,
ring expansions, Wagner-Meerwein, and pinacol rearrangement pathways.
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


CARBOCATION_STABILITY_ORDER = [
    ('resonance-stabilized (allylic/benzylic)', 10),
    ('tertiary (3°)', 8),
    ('secondary (2°)', 4),
    ('primary (1°)', 1),
    ('methyl', 0),
]


@ChemMCPManager.register_tool
class CarbocationRearrangement(BaseTool):
    __version__ = "0.1.0"
    name = "CarbocationRearrangement"
    func_name = 'explain_carbocation_rearrangement'
    description = "Explain carbocation rearrangement mechanisms: hydride (H⁻) shift, alkyl (R⁻) shift, ring expansion/contraction, Wagner-Meerwein rearrangement, and pinacol rearrangement. Predicts the most favorable rearrangement pathway to a more stable carbocation."
    implementation_description = "Analyzes the substrate structure around the initial carbocation position, identifies possible 1,2-shift donors (adjacent C-H for hydride, adjacent C-C for alkyl), evaluates the stability gain from each possible rearrangement, predicts the thermodynamically favored product, and provides detailed stepwise mechanism with driving force analysis."
    categories = ["Reaction"]
    tags = ["Carbocation", "Rearrangement", "Hydride Shift", "Alkyl Shift", "Wagner-Meerwein", "Pinacol", "Ring Expansion"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES of the molecule that forms or contains a carbocation.'),
        ('cation_position', 'str', 'auto', 'Position of initial carbocation: auto-detect, or specify atom index/description.'),
        ('rearrangement_type', 'str', 'all', 'Type to analyze: all, hydride_only, alkyl_only, ring_expansion.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles [cation_position] [rearrangement_type]. E.g., "CC(C)(C)CCl auto all".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing initial_cation_analysis, possible_shifts, best_pathway, final_product, and driving_force.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'CC(C)CCl',
                'cation_position': 'auto',
                'rearrangement_type': 'all',
            },
            'text_input': {'query': 'CC(C)CCl'},
            'output': {
                'result': {
                    'substrate': 'isobutyl chloride → 2° cation after ionization',
                    'initial_cation': {'type': 'secondary', 'stability': 'moderate'},
                    'possible_rearrangements': [
                        {'type': 'hydride shift', 'from': 'C3 (tertiary)', 'to': 'C2 (cation)', 'product_cation': 'tertiary', 'driving_force': '+1 stability level'},
                    ],
                    'most_favorable': 'hydride shift from tertiary C-H → tertiary carbocation',
                    'final_product': 'tert-butyl cation (or tert-butyl alcohol after Nu attack)',
                    'favorability': 'very likely — significant stability gain',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'C1CCCC1COH',  # cyclopentylmethanol derivative
                'cation_position': 'auto',
                'rearrangement_type': 'all',
            },
            'text_input': {'query': 'C1CCCC1COH'},
            'output': {
                'result': {
                    'initial_cation': 'primary benzylic-type cation (on side chain)',
                    'possible_rearrangements': ['ring expansion to cyclohexyl cation'],
                    'most_favorable': 'ring expansion (Wagner-Meerwein): 5-membered → 6-membered ring',
                    'final_product': 'cyclohexanone (after pinacol-like rearrangement)',
                    'favorability': 'likely — ring strain relief + stability increase',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, cation_position: str = 'auto', rearrangement_type: str = 'all') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Identify/analyze initial carbocation position
        cation_analysis = self._identify_carbocation(mol, cation_position)

        # 2. Find all possible 1,2-shifts
        possible_shifts = self._find_possible_shifts(mol, cation_analysis, rearrangement_type)

        # 3. Evaluate each shift's driving force
        evaluated_shifts = self._evaluate_shifts(possible_shifts, cation_analysis)

        # 4. Determine best pathway
        best_pathway = self._select_best_pathway(evaluated_shifts, cation_analysis)

        # 5. Final product prediction
        final_product = self._predict_final_product(cation_analysis, best_pathway)

        # 6. Favorability
        favorability = self._evaluate_favorability(cation_analysis, best_pathway)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'initial_cation_analysis': cation_analysis,
                'possible_rearrangements': evaluated_shifts,
                'best_pathway': best_pathway,
                'final_product_prediction': final_product,
                'favorability': favorability,
                'summary': self._build_summary(cation_analysis, best_pathway, final_product, favorability),
            }
        }

        logger.info(f"Carbocation rearrangement: {substrate_smiles} → {best_pathway.get('shift_type', 'none')} → {favorability}")
        return result

    def _identify_carbocation(self, mol, position):
        """Identify the initial carbocation."""
        if position != 'auto':
            try:
                idx = int(position)
                atom = mol.GetAtomWithIdx(idx)
                return self._classify_cation(atom, mol)
            except ValueError:
                pass

        # Auto-detect: find the most likely cation site
        # Look for carbon with leaving group attached
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ('Cl', 'Br', 'I'):
                        return self._classify_cation(atom, mol)

        # If no LG found, look for -OH group (could be protonated in acid)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        for nn in neighbor.GetNeighbors():
                            if nn.GetAtomicNum() == 1:  # O-H
                                return self._classify_cation(atom, mol)

        # Default: find least substituted carbon with H
        candidates = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
                n_c = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
                candidates.append((atom, n_c))

        if candidates:
            best = min(candidates, key=lambda x: x[1])
            return self._classify_cation(best[0], mol)

        return {'type': 'unknown', 'stability': 'unknown', 'position': None}

    def _classify_cation(self, atom, mol):
        """Classify carbocation type at this atom."""
        neighbors = atom.GetNeighbors()
        n_c_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6)

        # Exclude heteroatoms from substitution count for basic classification
        n_non_h_non_x = sum(1 for n in neighbors if n.GetAtomicNum() not in (1,) and n.GetSymbol() not in ('Cl','Br','I'))

        if n_non_h_non_x >= 3: ctype = 'tertiary (3°)'
        elif n_non_h_non_x == 2: ctype = 'secondary (2°)'
        elif n_non_h_non_x == 1: ctype = 'primary (1°)'
        else: ctype = 'methyl'

        # Check for resonance stabilization
        is_allylic = any(
            mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()) and
            mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() >= 2.0
            for n in neighbors if n.GetAtomicNum() == 6
        )
        is_benzylic = any(n.GetIsAromatic() for n in neighbors)

        if is_benzylic: ctype += ', benzylic'
        if is_allylic: ctype += ', allylic'

        stability_score = 8 if 'tertiary' in ctype else (4 if 'secondary' in ctype else 1)
        if is_benzylic: stability_score += 2
        if is_allylic: stability_score += 1

        return {
            'position': atom.GetIdx(),
            'type': ctype,
            'stability_score': stability_score,
            'n_carbon_neighbors': n_c_neighbors,
            'is_allylic': is_allylic,
            'is_benzylic': is_benzylic,
        }

    def _find_possible_shifts(self, mol, cation_analysis, rtype):
        """Find all possible 1,2-shifts."""
        pos = cation_analysis.get('position')
        if pos is None:
            return []

        cation_atom = mol.GetAtomWithIdx(pos)
        shifts = []
        analyzed_indices = set()

        for neighbor in cation_atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            if nidx in analyzed_indices:
                continue
            analyzed_indices.add(nidx)

            if neighbor.GetAtomicNum() != 6:
                continue

            # Hydride shift: neighbor has H atoms
            n_h = neighbor.GetTotalNumHs()
            if n_h > 0 and rtype in ('all', 'hydride_only'):
                donor_class = self._classify_donor_carbon(neighbor, mol)
                shifts.append({
                    'shift_type': 'hydride shift (H⁻)',
                    'donor_idx': nidx,
                    'donor_class': donor_class,
                    'description': f"H⁻ migrates from C({nidx}) ({donor_class}) to cationic center C({pos})",
                })

            # Alkyl shift: neighbor has carbon substituent(s) that could migrate
            if rtype in ('all', 'alkyl_only'):
                for sub_neighbor in neighbor.GetNeighbors():
                    if sub_neighbor.GetIdx() == pos:
                        continue
                    if sub_neighbor.GetAtomicNum() == 6:
                        # This carbon could migrate with its bond
                        n_sub_of_sub = sum(1 for sn in sub_neighbor.GetNeighbors()
                                           if sn.GetAtomicNum() == 6 and sn.GetIdx() != nidx)
                        migrating_group = f'C({sub_neighbor.GetIdx()})' if n_sub_of_sub == 0 else \
                                          f'CH₂-' if n_sub_of_sub <= 1 else \
                                          f'alkyl (≥C₂)'
                        shifts.append({
                            'shift_type': f'alkyl shift ({migrating_group}⁻)',
                            'donor_idx': nidx,
                            'acceptor_idx': pos,
                            'migrating_atom_idx': sub_neighbor.GetIdx(),
                            'donor_class': self._classify_donor_carbon(neighbor, mol),
                            'description': f"{migrating_group} migrates from C({nidx}) to C({pos})",
                        })

            # Ring expansion check
            ring_info = mol.GetRingInfo()
            if ring_info.NumRings() > 0:
                for ring in ring_info.AtomRings():
                    if pos in ring and nidx in ring:
                        ring_size = len(ring)
                        if ring_size <= 5:  # Small rings expand readily
                            shifts.append({
                                'shift_type': f'ring expansion ({ring_size}→{ring_size+1})',
                                'donor_idx': nidx,
                                'ring_atoms': list(ring),
                                'current_ring_size': ring_size,
                                'new_ring_size': ring_size + 1,
                                'description': (
                                    f"Ring bond migration: {ring_size}-membered ring expands "
                                    f"to {ring_size+1}-membered ring (relieves angle strain)"
                                ),
                            })

        return shifts

    def _classify_donor_carbon(self, atom, mol):
        """Classify the carbon donating H or alkyl group."""
        n_c = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
        if n_c >= 3: return 'tertiary'
        elif n_c == 2: return 'secondary'
        elif n_c == 1: return 'primary'
        return 'methyl'

    def _evaluate_shifts(self, shifts, cation_analysis):
        """Evaluate each shift's driving force."""
        current_score = cation_analysis.get('stability_score', 0)
        current_type = cation_analysis.get('type', '')

        for shift in shifts:
            stype = shift.get('shift_type', '')
            donor_class = shift.get('donor_class', '')

            # Estimate new cation stability after shift
            if 'tertiary' in donor_class or 'ring expansion' in stype:
                new_score = 8  # likely becomes tertiary
            elif 'secondary' in donor_class:
                new_score = 4
            else:
                new_score = current_score + 1  # minimal improvement

            # Ring expansion bonus
            if 'expansion' in stype:
                old_size = shift.get('current_ring_size', 6)
                if old_size <= 4: new_score += 3  # big strain relief
                elif old_size == 5: new_score += 2

            shift['stability_gain'] = new_score - current_score
            shift['new_cation_stability'] = ('very high' if new_score >= 8 else
                                             'high' if new_score >= 5 else
                                             'moderate' if new_score >= 3 else 'low')
            shift['driving_force'] = (
                f"{'Strong' if shift['stability_gain'] >= 4 else 'Moderate' if shift['stability_gain'] >= 2 else 'Weak'} driving force: "
                f"{current_type} → more stable cation (+{shift['stability_gain']} stability units)"
            )

        # Sort by stability gain (best first)
        shifts.sort(key=lambda s: s.get('stability_gain', 0), reverse=True)
        return shifts

    def _select_best_pathway(self, shifts, cation_analysis):
        """Select the most favorable rearrangement pathway."""
        if not shifts:
            return {'shift_type': 'no rearrangement', 'reason': 'No favorable 1,2-shifts available.', 'stability_gain': 0}

        best = shifts[0]
        if best.get('stability_gain', 0) <= 0:
            return {'shift_type': 'no rearrangement favored', 'reason': 'All possible shifts lead to equally or less stable cations.', 'stability_gain': 0}

        return {
            **best,
            'reason': f"Most favorable: {best.get('driving_force', '?')}",
            'will_occur': best.get('stability_gain', 0) >= 2,
        }

    def _predict_final_product(self, cation_analysis, best_pathway):
        """Predict final structure after rearrangement."""
        ptype = best_pathway.get('shift_type', '')
        if 'no rearrangement' in ptype:
            return {'product': 'No rearrangement — original cation reacts directly'}

        if 'hydride' in ptype:
            return {
                'product': 'Rearranged carbocation (more substituted)',
                'structural_change': 'H migrated to original cation center; cation moved to adjacent carbon',
                'then': 'Nucleophile attacks new (more stable) cation → rearranged product',
            }
        elif 'alkyl' in ptype:
            return {
                'product': 'Rearranged carbocation (skeletal change)',
                'structural_change': 'Alkyl group migrated; skeleton rearranged',
                'then': 'Nucleophile attacks new cation → skeletal-rearranged product',
            }
        elif 'ring expansion' in ptype:
            old_size = best_pathway.get('current_ring_size', '?')
            return {
                'product': f'Ring-expanded product ({old_size}→{old_size+1} membered)',
                'structural_change': f'{old_size}-membered ring expanded by one atom',
                'then': 'Nucleophile captures ring-expanded cation',
            }
        return {'product': 'Rearranged product'}

    def _evaluate_favorability(self, cation_analysis, best_pathway):
        gain = best_pathway.get('stability_gain', 0)
        will = best_pathway.get('will_occur', False)

        if not will:
            ctype = cation_analysis.get('type', '')
            if 'tertiary' in ctype or 'benzylic' in ctype:
                return 'already stable — rearrangement unlikely'
            return 'unstable cation but no good rearrangement path'

        if gain >= 4: return 'very likely — strong driving force'
        if gain >= 2: return 'likely — moderate driving force'
        return 'possible but may compete with direct capture'

    def _build_summary(self, cation, best, product, fav):
        init = cation.get('type', '?')
        shift = best.get('shift_type', '?')
        prod = product.get('product', '?')
        return f"Carbocation rearrangement: {init} cation → [{shift}] → {prod}. {fav}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        pos = parts[1] if len(parts) > 1 else 'auto'
        rtype = parts[2] if len(parts) > 2 else 'all'
        return self._run_base(sub, pos, rtype)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
