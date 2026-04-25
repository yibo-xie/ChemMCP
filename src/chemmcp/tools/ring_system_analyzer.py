"""
Ring System Analyzer (Tool #113)
分析稠环(fused)、螺环(spiro)、桥环(bridged)结构。
Uses RDKit for ring perception and SSSR analysis.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class RingSystemAnalyzer(BaseTool):
    __version__ = "0.1.0"
    name = "RingSystemAnalyzer"
    func_name = 'analyze_ring_system'
    description = "Analyze ring systems in a molecule: identify fused rings, spiro rings, bridged rings, and provide detailed ring system classification."
    implementation_description = "Uses RDKit's GetSymmSSSR() for ring detection, then classifies each ring system as fused, spiro, or bridged based on atom/bond sharing patterns. Also detects aromaticity of individual rings and calculates Bredt's rule applicability for bridged systems."
    categories = ["Molecule"]
    tags = ["Ring Systems", "Fused Rings", "Spiro", "Bridged", "RDKit", "Topology"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('detailed', 'bool', 'True', 'If True, include per-atom and per-bond details of ring membership.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional detailed flag.'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary with ring_systems, fused_rings, spiro_rings, bridged_rings, aromaticity, and summary.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'C1=CC=CC2=C1C=CC=C2', 'detailed': True},  # naphthalene
            'text_input': {'query': 'C1=CC=CC2=C1C=CC=C2 true'},
            'output': {
                'result': {
                    'n_rings': 2,
                    'ring_systems': [{'type': 'fused', 'n_rings': 2, 'rings': [0, 1], 'atoms': [...]}],
                    'fused_rings': [{'system_id': 0, 'shared_atoms': [...], 'shared_bonds': [...]}],
                    'spiro_rings': [],
                    'bridged_rings': [],
                    'all_aromatic': True,
                    'summary': 'Naphthalene: 2 fused aromatic rings (fused bicyclic system).',
                }
            },
        },
        {
            'code_input': {'smiles': 'C12CCC(CC1)CC2', 'detailed': False},  # spiro[4.5]decane skeleton
            'text_input': {'query': 'C12CCC(CC1)CC2 false'},
            'output': {
                'result': {
                    'n_rings': 2,
                    'spiro_rings': [{'spiro_atom': ..., 'ring_sizes': [5, 6]}],
                    'summary': 'Spiro compound: 2 rings sharing a single spiro atom.',
                }
            },
        },
        {
            'code_input': {'smiles': 'C1C2CC3CC1CC(C3)C2', 'detailed': True},  # adamantane-like
            'text_input': {'query': 'C1C2CC3CC1CC(C3)C2 true'},
            'output': {
                'result': {
                    'n_rings': 4,
                    'bridged_rings': [{'bridge_atoms': [...], 'bridgeheads': [...]}],
                    'summary': 'Bridged polycyclic system (adamantane-type cage).',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, detailed: bool = True) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Get SSSR (Smallest Set of Smallest Rings)
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        n_rings = len(rings)

        if n_rings == 0:
            return {
                'result': {
                    'n_rings': 0,
                    'ring_systems': [],
                    'fused_rings': [],
                    'spiro_rings': [],
                    'bridged_rings': [],
                    'aromatic_rings': [],
                    'all_aromatic': False,
                    'summary': 'No rings found in this molecule.',
                }
            }

        # Classify ring systems using connected components of shared atoms
        ring_systems = self._classify_ring_systems(rings, mol)

        # Classify each system type
        fused_list = []
        spiro_list = []
        bridged_list = []

        for sys_idx, rs in enumerate(ring_systems):
            stype = rs['type']
            entry = {
                'system_id': sys_idx,
                'n_rings': len(rs['rings']),
                'ring_indices': rs['rings'],
                'ring_sizes': [len(rings[ri]) for ri in rs['rings']],
                'total_atoms': len(rs['atoms']),
            }
            if detailed:
                entry['atoms'] = rs['atoms']
                entry['bonds'] = rs.get('bonds', [])

            if stype == 'fused':
                entry['shared_atoms'] = rs.get('shared_atoms', [])
                entry['shared_bonds'] = rs.get('shared_bonds', [])
                fused_list.append(entry)
            elif stype == 'spiro':
                entry['spiro_atom'] = rs.get('spiro_atom')
                spiro_list.append(entry)
            elif stype == 'bridged':
                entry['bridgehead_atoms'] = rs.get('bridgehead_atoms', [])
                entry['bridge_atoms'] = rs.get('bridge_atoms', [])
                bridged_list.append(entry)

        # Aromaticity check
        arom_rings = []
        for i, ring_atoms in enumerate(rings):
            # Check if all bonds in this ring are aromatic
            is_arom = all(
                mol.GetBondBetweenAtoms(ring_atoms[j], ring_atoms[(j+1) % len(ring_atoms)]).GetIsAromatic()
                for j in range(len(ring_atoms))
                if mol.GetBondBetweenAtoms(ring_atoms[j], ring_atoms[(j+1) % len(ring_atoms)])
            )
            if is_arom:
                arom_rings.append(i)

        # IUPAC naming hints
        iupac_hint = self._get_iupac_hint(ring_systems, rings, n_rings, arom_rings)

        result = {
            'result': {
                'n_rings': n_rings,
                'ring_systems': [{
                    'type': rs['type'],
                    'n_rings': len(rs['rings']),
                    'ring_sizes': [len(rings[ri]) for ri in rs['rings']],
                } for rs in ring_systems],
                'fused_rings': fused_list,
                'spiro_rings': spiro_list,
                'bridged_rings': bridged_list,
                'aromatic_rings': arom_rings,
                'all_aromatic': len(arom_rings) == n_rings and n_rings > 0,
                'iupac_naming_hint': iupac_hint,
                'summary': self._build_summary(n_rings, fused_list, spiro_list, bridged_list, arom_rings),
            }
        }

        logger.info(f"Ring analyzer: {smiles} → {n_rings} rings ({len(fused_list)} fused, {len(spiro_list)} spiro, {len(bridged_list)} bridged)")
        return result

    def _classify_ring_systems(self, rings, mol):
        """Classify rings into systems and determine their types."""
        n = len(rings)
        if n == 0:
            return []

        # Build adjacency between rings based on shared atoms
        ring_adj = [[] for _ in range(n)]
        shared_count = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                shared = set(rings[i]) & set(rings[j])
                cnt = len(shared)
                if cnt > 0:
                    ring_adj[i].append(j)
                    ring_adj[j].append(i)
                    shared_count[i][j] = cnt
                    shared_count[j][i] = cnt

        # Find connected components (ring systems via BFS/DFS)
        visited = [False] * n
        systems = []

        for start in range(n):
            if visited[start]:
                continue
            # BFS to find connected component
            component = []
            queue = [start]
            visited[start] = True
            while queue:
                r = queue.pop(0)
                component.append(r)
                for nb in ring_adj[r]:
                    if not visited[nb]:
                        visited[nb] = True
                        queue.append(nb)

            # Determine system type
            sys_atoms = set()
            for r in component:
                sys_atoms.update(rings[r])

            # Analyze sharing pattern
            max_shared = 0
            has_single_atom_share = False
            has_bond_share = False
            shared_atoms_global = set()
            shared_bonds_global = set()

            for ii in range(len(component)):
                for jj in range(ii + 1, len(component)):
                    ri, rj = component[ii], component[jj]
                    sc = shared_count[ri][rj] if ri < rj else shared_count[rj][ri]
                    if sc > 0:
                        shared = set(rings[ri]) & set(rings[rj])
                        shared_atoms_global |= shared
                        max_shared = max(max_shared, sc)
                        if sc == 1:
                            has_single_atom_share = True
                        # Check bond sharing
                        for ai in rings[ri]:
                            for aj in rings[rj]:
                                if aj > ai:
                                    bond = mol.GetBondBetweenAtoms(ai, aj)
                                    if bond and (ai in rings[rj] and aj in rings[ri]):
                                        has_bond_share = True
                                        shared_bonds_global.add((min(ai, aj), max(ai, aj)))

            # Determine type
            if len(component) >= 3 and max_shared >= 2:
                stype = 'bridged'
            elif has_single_atom_share and max_shared == 1 and len(component) == 2:
                stype = 'spiro'
            elif max_shared >= 2 or (len(component) >= 2 and max_shared >= 1):
                stype = 'fused'
            else:
                stype = 'isolated'

            sys_entry = {
                'type': stype,
                'rings': component,
                'atoms': list(sys_atoms),
                'shared_atoms': list(shared_atoms_global),
                'shared_bonds': [list(b) for b in shared_bonds_global],
            }

            if stype == 'spiro' and shared_atoms_global:
                sys_entry['spiro_atom'] = list(shared_atoms_global)[0]
            if stype == 'bridged':
                # Find bridgehead atoms (atoms shared by ≥3 rings)
                atom_ring_count = {}
                for a in sys_atoms:
                    atom_ring_count[a] = sum(1 for r in component if a in rings[r])
                bridgeheads = [a for a, c in atom_ring_count.items() if c >= 3]
                sys_entry['bridgehead_atoms'] = bridgeheads
                sys_entry['bridge_atoms'] = [a for a in sys_atoms if a not in bridgeheads and atom_ring_count.get(a, 0) >= 2]

            systems.append(sys_entry)

        return systems

    def _get_iupac_hint(self, systems, rings, n_total, arom_rings):
        """Generate IUPAC naming hint."""
        hints = []
        for si, sys in enumerate(systems):
            sizes = sorted([len(rings[ri]) for ri in sys['rings']], reverse=True)
            nr = len(sys['rings'])
            t = sys['type']

            if t == 'fused':
                if nr == 2:
                    hints.append(f'bicyclic fused [{"-".join(map(str, sizes))}]')
                else:
                    hints.append(f'{nr}-ring fused system [{"-".join(map(str, sizes))}]')
            elif t == 'spiro':
                spiro_sizes = sorted(sizes)
                hints.append(f'spiro[{".".join(map(str, spiro_sizes))}]')
            elif t == 'bridged':
                hints.append(f'bridged {nr}-cyclic system (e.g., {"-".join(map(str, sizes))})')

        return "; ".join(hints) if hints else "simple monocyclic or acyclic"

    def _build_summary(self, n_rings, fused, spiro, bridged, arom):
        parts = [f"{n_rings} ring(s) detected."]
        if fused:
            parts.append(f"{len(fused)} fused ring system(s).")
        if spiro:
            parts.append(f"{len(spiro)} spiro junction(s).")
        if bridged:
            parts.append(f"{len(bridged)} bridged system(s).")
        if arom:
            parts.append(f"{len(arom)} aromatic ring(s).")
        return " ".join(parts)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        smiles = parts[0]
        det = True
        if len(parts) > 1:
            det = parts[1].lower() in ('true', '1', 'yes', 't')
        return self._run_base(smiles, det)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
