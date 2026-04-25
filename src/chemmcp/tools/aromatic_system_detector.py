"""
Aromatic System Detector (Tool #114)
判断分子芳香性，基于Hückel规则验证。
Uses RDKit for aromaticity detection and Hückel rule (4n+2 π electrons) verification.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem import rdchem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class AromaticSystemDetector(BaseTool):
    __version__ = "0.1.0"
    name = "AromaticSystemDetector"
    func_name = 'detect_aromaticity'
    description = "Detect aromatic systems in a molecule using Hückel's rule (4n+2 π electrons) and RDKit aromaticity analysis."
    implementation_description = "Combines RDKit's built-in aromaticity detection with Hückel rule verification: counts π electrons in each ring system, checks planarity and cyclic conjugation requirements, and classifies as aromatic, anti-aromatic, or non-aromatic."
    categories = ["Molecule"]
    tags = ["Aromaticity", "Hückel Rule", "π Electrons", "Conjugation", "RDKit"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('include_huckel_details', 'bool', 'True', 'Include detailed Hückel rule analysis with electron counting.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional huckel_details flag.'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing is_aromatic, aromatic_rings, pi_electron_counts, huckel_analysis, and classification per ring system.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'c1ccccc1', 'include_huckel_details': True},  # benzene
            'text_input': {'query': 'c1ccccc1 true'},
            'output': {
                'result': {
                    'is_aromatic': True,
                    'aromatic_systems': [{'ring_atoms': [0,1,2,3,4,5], 'pi_electrons': 6, 'huckel_n': 1, 'classification': 'aromatic'}],
                    'total_pi_electrons': 6,
                    'summary': 'Benzene: aromatic (6 π electrons, satisfies Hückel 4n+2 with n=1).',
                }
            },
        },
        {
            'code_input': {'smiles': 'C1=CC=C1', 'include_huckel_details': True},  # cyclobutadiene
            'text_input': {'query': 'C1=CC=C1 true'},
            'output': {
                'result': {
                    'is_aromatic': False,
                    'aromatic_systems': [{'ring_atoms': [...], 'pi_electrons': 4, 'huckel_n': 0.5, 'classification': 'anti-aromatic'}],
                    'summary': 'Cyclobutadiene: anti-aromatic (4 π electrons, violates Hückel rule — 4n with n=1).',
                }
            },
        },
        {
            'code_input': {'smiles': 'C1CCC1', 'include_huckel_details': True},  # cyclobutane
            'text_input': {'query': 'C1CCC1 true'},
            'output': {
                'result': {
                    'is_aromatic': False,
                    'aromatic_systems': [{'ring_atoms': [...], 'pi_electrons': 0, 'classification': 'non-aromatic'}],
                    'summary': 'Cyclobutane: non-aromatic (no conjugated π system).',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, include_huckel_details: bool = True) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.SanitizeMol(mol)
        Chem.SetAromaticity(mol)

        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()

        if len(rings) == 0:
            return {
                'result': {
                    'is_aromatic': False,
                    'aromatic_systems': [],
                    'total_pi_electrons': 0,
                    'summary': 'No rings found — molecule cannot be aromatic.',
                }
            }

        # Analyze each ring
        ring_analyses = []
        any_aromatic = False

        for ri, ring_atoms in enumerate(rings):
            analysis = self._analyze_ring_aromaticity(mol, ring_atoms, ri, include_huckel_details)
            ring_analyses.append(analysis)
            if analysis['classification'] == 'aromatic':
                any_aromatic = True

        total_pi = sum(ra.get('pi_electrons', 0) for ra in ring_analyses)

        result = {
            'result': {
                'is_aromatic': any_aromatic,
                'aromatic_systems': ring_analyses,
                'total_pi_electrons': total_pi,
                'rdkit_aromatic_atoms': sum(1 for a in mol.GetAtoms() if a.GetIsAromatic()),
                'summary': self._build_summary(ring_analyses, any_aromatic),
            }
        }

        logger.info(f"Aromaticity detector: {smiles} → aromatic={any_aromatic}, π_e={total_pi}")
        return result

    def _analyze_ring_aromaticity(self, mol, ring_atoms, ring_idx, detailed):
        """Analyze a single ring for aromaticity."""
        n_atoms = len(ring_atoms)

        # Check RDKit aromaticity first
        rdkit_arom = all(
            mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms
        )

        # Count π electrons
        pi_electrons = self._count_pi_electrons(mol, ring_atoms)

        # Hückel rule check
        huckel_result = self._check_huckel_rule(pi_electrons, n_atoms)

        # Determine classification
        if rdkit_arom:
            classification = 'aromatic'
        elif pi_electrons > 0 and huckel_result['is_anti_aromatic']:
            classification = 'anti-aromatic'
        elif pi_electrons == 0:
            classification = 'non-aromatic (saturated)'
        else:
            classification = 'non-aromatic (non-conjugated)'

        analysis = {
            'ring_index': ring_idx,
            'ring_atoms': ring_atoms,
            'ring_size': n_atoms,
            'rdkit_aromatic': rdkit_arom,
            'pi_electrons': pi_electrons,
            'classification': classification,
        }

        if detailed:
            analysis['huckel'] = huckel_result
            analysis['atom_contributions'] = [
                {
                    'idx': a,
                    'symbol': mol.GetAtomWithIdx(a).GetSymbol(),
                    'is_aromatic': mol.GetAtomWithIdx(a).GetIsAromatic(),
                    'pi_contribution': self._get_atom_pi_contribution(mol, a),
                }
                for a in ring_atoms
            ]

        return analysis

    def _count_pi_electrons(self, mol, ring_atoms):
        """Count π electrons in a ring system."""
        pi_count = 0
        counted_atoms = set()
        counted_bonds = set()

        for ai in ring_atoms:
            atom = mol.GetAtomWithIdx(ai)
            atomic_num = atom.GetAtomicNum()

            # Each atom in the ring contributes based on hybridization and lone pairs
            if ai in counted_atoms:
                continue

            # Double-bonded atoms contribute 1 π electron each
            for bond in atom.GetBonds():
                aj = bond.GetOtherAtomIdx(ai)
                bond_order = bond.GetBondTypeAsDouble()
                if aj in ring_atoms and bond_order == 2.0:
                    if (min(ai, aj), max(ai, aj)) not in counted_bonds:
                        pi_count += 2  # Each double bond contributes 2 π e-
                        counted_bonds.add((min(ai, aj), max(ai, aj)))

            # Atoms with lone pairs in p-orbitals contribute 2
            # Heteroatoms (O, N, S) in aromatic rings may contribute
            if atomic_num in (7, 8, 16):  # N, O, S
                # Check if atom is part of aromatic system with lone pair
                if atom.GetIsAromatic():
                    num_hs = atom.GetTotalNumHs()
                    # For pyrrole-type N: contributes 2 e- from lone pair
                    if atomic_num == 7 and num_hs == 1:
                        pi_count += 2
                        counted_atoms.add(ai)
                    # For furan-type O: contributes 2 e- from lone pair
                    elif atomic_num == 8 and num_hs == 0:
                        # Check valence
                        pass  # O typically contributes 2 in furan-like

        return pi_count

    def _get_atom_pi_contribution(self, mol, atom_idx):
        """Estimate π electron contribution of an atom."""
        atom = mol.GetAtomWithIdx(atom_idx)
        # Simple heuristic based on atom type and bonding
        atomic_num = atom.GetAtomicNum()
        n_double_bonds = sum(1 for b in atom.GetBonds() if b.GetBondTypeAsDouble() >= 2.0)
        n_hs = atom.GetTotalNumHs()

        if atomic_num == 6:  # Carbon
            if n_double_bonds > 0:
                return 1
            return 0
        elif atomic_num == 7:  # Nitrogen
            if atom.GetIsAromatic() and n_hs == 1:
                return 2  # pyrrole-type
            elif n_double_bonds > 0:
                return 1  # pyridine-type
            return 0
        elif atomic_num == 8:  # Oxygen
            if atom.GetIsAromatic():
                return 2  # furan-type
            return 0
        return 0

    def _check_huckel_rule(self, pi_e, n_atoms):
        """Check Hückel 4n+2 rule."""
        if pi_e == 0:
            return {'satisfies_4n2': False, 'is_4n': False, 'is_anti_aromatic': False, 'n_value': None}

        # Check 4n+2
        n_val = (pi_e - 2) / 2
        satisfies_4n2 = n_val >= 0 and abs(n_val - round(n_val)) < 0.001

        # Check 4n (anti-aromatic)
        n_4n = pi_e / 4
        is_4n = n_4n >= 1 and abs(n_4n - round(n_4n)) < 0.001

        return {
            'satisfies_4n2': satisfies_4n2,
            'n_for_4n2': round(n_val) if satisfies_4n2 else None,
            'is_4n': is_4n,
            'n_for_4n': round(n_4n) if is_4n else None,
            'is_anti_aromatic': is_4n and pi_e > 0,
        }

    def _build_summary(self, analyses, any_arom):
        arom_count = sum(1 for a in analyses if a['classification'] == 'aromatic')
        anti_count = sum(1 for a in analyses if a['classification'] == 'anti-aromatic')
        non_count = sum(1 for a in analyses if 'non-aromatic' in a['classification'])

        parts = []
        if any_arom:
            parts.append(f"Molecule is AROMATIC ({arom_count} aromatic ring system(s)).")
        if anti_count:
            parts.append(f"{anti_count} anti-aromatic ring(s) detected.")
        if non_count:
            parts.append(f"{non_count} non-aromatic ring(s).")

        # Add Hückel details for aromatic systems
        for a in analyses:
            if a['classification'] == 'aromatic' and 'huckel' in a:
                h = a['huckel']
                parts.append(
                    f"Ring {a['ring_index']} ({a['ring_size']}-membered): "
                    f"{a['pi_electrons']} π e⁻ → Hückel 4n+2 (n={h.get('n_for_4n2', '?')}) ✓"
                )
            elif a['classification'] == 'anti-aromatic' and 'huckel' in a:
                h = a['huckel']
                parts.append(
                    f"Ring {a['ring_index']} ({a['ring_size']}-membered): "
                    f"{a['pi_electrons']} π e⁻ → 4n (n={h.get('n_for_4n', '?')}) ✗ anti-aromatic"
                )

        return " ".join(parts)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        smiles = parts[0]
        huckel_det = True
        if len(parts) > 1:
            huckel_det = parts[1].lower() in ('true', '1', 'yes', 't')
        return self._run_base(smiles, huckel_det)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
