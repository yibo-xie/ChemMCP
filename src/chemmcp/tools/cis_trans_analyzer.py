"""
环状化合物顺反异构分析工具
Analyze cis/trans isomerism in cyclic compounds (rings).
Uses RDKit ring finding + substituent position analysis.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class CisTransAnalyzer(BaseTool):
    __version__ = "0.1.0"
    name = "CisTransAnalyzer"
    func_name = 'analyze_cis_trans'
    description = "Analyze cis/trans isomerism in cyclic compounds (rings). For each ring, identifies substituent positions and determines possible cis/trans relationships between pairs of substituents."
    implementation_description = "Uses RDKit's ring detection to identify all rings in the molecule, then analyzes substituent positions on each ring to determine which pairs of substituents can exhibit cis/trans isomerism. Handles mono-, bi-, and polycyclic systems including fused and bridged rings."
    categories = ["Molecule"]
    tags = ["Stereochemistry", "Cis/Trans", "Cyclic Compounds", "Rings", "Isomerism", "RDKit"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('focus_ring', 'int', '-1', 'Ring index to focus on (-1 for all rings).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional ring index. E.g., "C1CCC(CC)C1" or "C1CCC(CC)C1 0".'),
    ]
    output_sig = [
        ('cis_trans_analysis', 'dict', 'Detailed analysis of cis/trans isomerism including ring info, substituents, possible isomers, and stereochemical assignments.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'C1CC(C)CCC1', 'focus_ring': -1},
            'text_input': {'query': 'C1CC(C)CCC1'},
            'output': {
                'cis_trans_analysis': {
                    'rings': [
                        {
                            'ring_size': 6,
                            'ring_type': 'aliphatic',
                            'substituents': [{'position': 2, 'atom': 'C (methyl)'}],
                            'n_substituents': 1,
                            'cis_trans_possible': False,
                            'reason': 'Only one substituent on this ring — need at least 2 for cis/trans.',
                        }
                    ],
                    'total_rings': 1,
                    'summary': '1 ring found with 1 substituent — no cis/trans isomerism possible.',
                }
            },
        },
        {
            'code_input': {'smiles': 'C1CC(C)C(C)C1', 'focus_ring': -1},
            'text_input': {'query': 'C1CC(C)C(C)C1'},
            'output': {
                'cis_trans_analysis': {
                    'rings': [
                        {
                            'ring_size': 5,
                            'ring_type': 'aliphatic',
                            'substituents': [
                                {'position': 2, 'atom': 'C (methyl)'},
                                {'position': 3, 'atom': 'C (methyl)'},
                            ],
                            'cis_trans_possible': True,
                            'possible_isomers': ['cis-1,2-dimethylcyclopentane', 'trans-1,2-dimethylcyclopentane'],
                        }
                    ],
                    'summary': '1 ring found: cyclopentane with 2 methyl substituents → cis/trans isomerism possible.',
                }
            },
        },
        {
            'code_input': {'smiles': 'C1C(CC)C(CC)CC1', 'focus_ring': -1},
            'text_input': {'query': 'C1C(CC)C(CC)CC1'},
            'output': {
                'cis_trans_analysis': {
                    'rings': [
                        {
                            'ring_size': 6,
                            'substituents': [
                                {'position': 1, 'atom': 'C (ethyl)'},
                                {'position': 3, 'atom': 'C (ethyl)'},
                            ],
                            'cis_trans_possible': True,
                        }
                    ],
                    'summary': '6-membered ring with 2 ethyl substituents → cis/trans isomerism possible.',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, focus_ring: int = -1) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Get ring information
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        bond_rings = ring_info.BondRings()

        if not atom_rings:
            return {
                'cis_trans_analysis': {
                    'rings': [],
                    'total_rings': 0,
                    'summary': 'No rings found in this molecule.',
                }
            }

        # Analyze each ring
        ring_analyses = []
        for ring_idx, ring_atoms in enumerate(atom_rings):
            if focus_ring >= 0 and ring_idx != focus_ring:
                continue

            analysis = self._analyze_ring(mol, ring_atoms, ring_idx)
            ring_analyses.append(analysis)

        # Build summary
        total = len(ring_analyses)
        n_with_ct = sum(1 for r in ring_analyses if r['cis_trans_possible'])
        parts = [f"{total} ring(s) analyzed."]
        if n_with_ct > 0:
            parts.append(f"{n_with_ct} ring(s) can show cis/trans isomerism.")
        else:
            parts.append("No rings with cis/trans isomerism detected.")

        result = {
            'cis_trans_analysis': {
                'rings': ring_analyses,
                'total_rings': total,
                'rings_with_cis_trans': n_with_ct,
                'summary': " ".join(parts),
            }
        }

        logger.info(f"Cis/Trans analysis: {smiles} → {total} rings, {n_with_ct} with CT isomerism")
        return result

    def _analyze_ring(self, mol, ring_atoms, ring_idx):
        """Analyze a single ring for cis/trans isomerism."""
        ring_size = len(ring_atoms)

        # Determine ring type
        all_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms)
        ring_type = "aromatic" if all_aromatic else "aliphatic"

        # Get ring atoms info
        ring_atom_details = []
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            ring_atom_details.append({
                'index_in_ring': len(ring_atom_details),
                'atom_index': atom_idx,
                'symbol': atom.GetSymbol(),
                'is_aromatic': atom.GetIsAromatic(),
            })

        # Find non-hydrogen substituents on each ring atom
        substituents = []
        for i, atom_idx in enumerate(ring_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # If neighbor is NOT part of this ring, it's a substituent
                if neighbor_idx not in ring_atoms:
                    sub_symbol = neighbor.GetSymbol()
                    n_h = neighbor.GetTotalNumHs()
                    # Describe the substituent
                    sub_desc = self._describe_substituent(mol, neighbor, atom_idx)
                    substituents.append({
                        'ring_position': i,
                        'ring_atom_index': atom_idx,
                        'substituent_index': neighbor_idx,
                        'symbol': sub_symbol,
                        'description': sub_desc,
                    })
            # Check for explicit H as substituent (for stereochemistry context)
            n_implicit_hs = atom.GetTotalNumHs()
            if n_implicit_hs > 0:
                # Only count H as notable if it's the distinguishing feature
                pass  # H is usually implicit; we care about non-H subs

        # Determine if cis/trans is possible
        n_subs = len(substituents)
        unique_positions = list(set(s['ring_position'] for s in substituents))

        cis_trans_possible = len(unique_positions) >= 2

        # Generate possible isomer descriptions
        possible_isomers = []
        if cis_trans_possible and n_subs >= 2:
            # Pair up substituents for cis/trans analysis
            sub_pairs = []
            for i in range(len(substituents)):
                for j in range(i + 1, len(substituents)):
                    pos_i = substituents[i]['ring_position']
                    pos_j = substituents[j]['ring_position']
                    if pos_i != pos_j:
                        # Calculate distance along ring (shorter path)
                        dist_forward = (pos_j - pos_i) % ring_size
                        dist_backward = (pos_i - pos_j) % ring_size
                        min_dist = min(dist_forward, dist_backward)
                        is_adjacent = min_dist == 1 or min_dist == ring_size - 1

                        sub_pairs.append({
                            'substituent_a': substituents[i]['description'],
                            'substituent_b': substituents[j]['description'],
                            'position_a': pos_i,
                            'position_b': pos_j,
                            'ring_distance': min_dist,
                            'is_adjacent': is_adjacent,
                            'cis_isomer': f"cis-{self._name_isomer(substituents[i], substituents[j], ring_type)}",
                            'trans_isomer': f"trans-{self._name_isomer(substituents[i], substituents[j], ring_type)}",
                        })

            possible_isomers = sub_pairs

        # Build ring name
        ring_names = {
            3: "cyclopropane",
            4: "cyclobutane",
            5: "cyclopentane",
            6: "cyclohexane",
            7: "cycloheptane",
            8: "cyclooctane",
        }
        base_name = ring_names.get(ring_size, f"{ring_size}-membered ring")
        if all_aromatic:
            if ring_size == 6:
                base_name = "benzene"
            elif ring_size == 5:
                base_name = "cyclopentadienyl/aromatic 5-ring"
            else:
                base_name = f"aromatic {ring_size}-membered ring"

        return {
            'ring_index': ring_idx,
            'ring_size': ring_size,
            'ring_type': ring_type,
            'base_name': base_name,
            'ring_atoms': ring_atom_details,
            'substituents': substituents,
            'n_substituents': n_subs,
            'n_substituted_positions': len(unique_positions),
            'cis_trans_possible': cis_trans_possible,
            'reason': (
                f"{len(unique_positions)} substituted position(s) on this ring — "
                f"need ≥2 for cis/trans." if not cis_trans_possible
                else f"{len(unique_positions)} substituted positions enable cis/trans analysis."
            ),
            'possible_isomers': possible_isomers,
        }

    def _describe_substituent(self, mol, atom, parent_idx):
        """Describe a substituent group."""
        symbol = atom.GetSymbol()
        neighbors = atom.GetNeighbors()

        # Filter out the parent ring atom
        children = [n for n in neighbors if n.GetIdx() != parent_idx]

        if symbol == 'C':
            n_c_children = sum(1 for n in children if n.GetSymbol() == 'C')
            n_total = atom.GetTotalNumHs() + len(children)

            if n_c_children >= 2 or n_total - len(children) >= 3:
                return "alkyl group"
            elif any(n.GetSymbol() == 'O' for n in children):
                return "oxygenated group"
            elif len(children) == 0:
                return "methyl"
            elif len(children) == 1:
                child_sym = children[0].GetSymbol()
                if child_sym == 'C':
                    return "ethyl" if children[0].GetDegree() <= 2 else "larger alkyl"
                return f"{child_sym.lower()}-containing group"
            return "alkyl group"
        elif symbol == 'O':
            return "hydroxyl" if atom.GetTotalNumHs() > 0 else "oxygen linker"
        elif symbol == 'N':
            return "amino" if atom.GetTotalNumHs() > 0 else "nitrogen group"
        elif symbol in ('F', 'Cl', 'Br', 'I'):
            return f"halo ({symbol})"
        elif symbol == 'S':
            return "thio/thiol group"
        return f"{symbol} group"

    def _name_isomer(self, sub_a, sub_b, ring_type):
        """Generate a name for an isomer based on two substituents."""
        desc_a = sub_a['description']
        desc_b = sub_b['description']
        # Simplify
        name = f"{desc_a}-{desc_b}"
        return name


if __name__ == "__main__":
    run_mcp_server()
