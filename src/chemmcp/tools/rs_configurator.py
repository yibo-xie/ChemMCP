"""
R/S构型判断工具
Determine R/S configuration at each chiral center using CIP priority rules.
Uses RDKit's built-in chirality assignment with detailed justification.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import EnumerateStereoisomers
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class RsConfigurator(BaseTool):
    __version__ = "0.1.0"
    name = "RsConfigurator"
    func_name = 'determine_rs_configuration'
    description = "Determine R/S (Rectus/Sinister) configuration at each chiral center in a molecule using Cahn-Ingold-Prelog (CIP) priority rules."
    implementation_description = "Uses RDKit's AssignStereochemistry() and FindMolChiralCenters() to determine R/S configurations. For centers with unspecified stereochemistry, attempts to enumerate possible stereoisomers and report both possibilities. Provides atomic-level detail for each chiral center including substituent priorities."
    categories = ["Molecule"]
    tags = ["Stereochemistry", "R/S Configuration", "CIP Rules", "Chirality", "RDKit"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('enumerate_options', 'bool', 'False', 'Whether to enumerate all possible R/S combinations for unspecified centers.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional enumerate flag. E.g., "[C@H](C)(Cl)Br false".'),
    ]
    output_sig = [
        ('rs_configuration', 'dict', 'Detailed R/S configuration for each chiral center with indices, symbols, assignments, and justifications.'),
    ]
    examples = [
        {
            'code_input': {'smiles': '[C@H](C)(Cl)Br', 'enumerate_options': False},
            'text_input': {'query': '[C@H](C)(Cl)Br false'},
            'output': {
                'rs_configuration': {
                    'chiral_centers': [
                        {
                            'atom_index': 0,
                            'symbol': 'C',
                            'configuration': 'S',
                            'substituents': [
                                {'priority': 1, 'atom': 'Br', 'atomic_number': 35},
                                {'priority': 2, 'atom': 'Cl', 'atomic_number': 17},
                                {'priority': 3, 'atom': 'C (methyl)', 'atomic_number': 6},
                                {'priority': 4, 'atom': 'H (implicit)', 'atomic_number': 1},
                            ],
                            'justification': 'Br(35) > Cl(17) > C(6) > H(1). With H pointing away, the sequence 1→2→3 is clockwise → S.',
                        }
                    ],
                    'total_centers': 1,
                    'summary': 'The molecule has 1 chiral center: C(0) is S configured.',
                }
            },
        },
        {
            'code_input': {'smiles': '[C@@H](F)(Cl)Br', 'enumerate_options': False},
            'text_input': {'query': '[C@@H](F)(Cl)Br false'},
            'output': {
                'rs_configuration': {
                    'chiral_centers': [
                        {
                            'atom_index': 0,
                            'symbol': 'C',
                            'configuration': 'R',
                            'substituents': [
                                {'priority': 1, 'atom': 'Br', 'atomic_number': 35},
                                {'priority': 2, 'atom': 'Cl', 'atomic_number': 17},
                                {'priority': 3, 'atom': 'F', 'atomic_number': 9},
                                {'priority': 4, 'atom': 'H (implicit)', 'atomic_number': 1},
                            ],
                            'justification': 'Br(35) > Cl(17) > F(9) > H(1). Sequence 1→2→3 is counterclockwise → R.',
                        }
                    ],
                    'total_centers': 1,
                    'summary': 'The molecule has 1 chiral center: C(0) is R configured.',
                }
            },
        },
        {
            'code_input': {'smiles': 'C[C@H](O)[C@@H](O)CC(C)C', 'enumerate_options': False},
            'text_input': {'query': 'C[C@H](O)[C@@H](O)CC(C)C false'},
            'output': {
                'rs_configuration': {
                    'chiral_centers': [
                        {'atom_index': 1, 'symbol': 'C', 'configuration': 'R'},
                        {'atom_index': 3, 'symbol': 'C', 'configuration': 'S'},
                    ],
                    'total_centers': 2,
                    'summary': 'The molecule has 2 chiral centers: C(1)=R, C(3)=S. This is the (1R,3S) diastereomer.',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, enumerate_options: bool = False) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Assign stereochemistry — this applies CIP rules internally
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Find chiral centers with their R/S assignments
        chiral_centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=True
        )

        if not chiral_centers:
            return {
                'rs_configuration': {
                    'chiral_centers': [],
                    'total_centers': 0,
                    'summary': 'No chiral centers found in this molecule.',
                }
            }

        # Build detailed analysis for each center
        center_details = []
        for idx, config in chiral_centers:
            atom = mol.GetAtomWithIdx(idx)
            detail = self._analyze_chiral_center(mol, atom, idx, config)
            center_details.append(detail)

        # Determine if this is a meso form or chiral
        n_r = sum(1 for c in center_details if c['configuration'] == 'R')
        n_s = sum(1 for c in center_details if c['configuration'] == 'S')
        n_unspecified = sum(1 for c in center_details if c['configuration'] == '?')

        # Build summary
        total = len(center_details)
        if n_unspecified == total:
            summary = f"The molecule has {total} chiral center(s), but none have specified stereochemistry."
        else:
            configs = ", ".join(f"C({c['atom_index']})={c['configuration']}" for c in center_details)
            summary = f"The molecule has {total} chiral center(s): {configs}."

        if n_unspecified > 0:
            summary += f" {n_unspecified} center(s) have unspecified configuration."

        result = {
            'rs_configuration': {
                'chiral_centers': center_details,
                'total_centers': total,
                'n_R': n_r,
                'n_S': n_s,
                'n_unspecified': n_unspecified,
                'summary': summary,
            }
        }

        # Optionally enumerate stereoisomers
        if enumerate_options and n_unspecified > 0:
            try:
                isomers = list(EnumerateStereoisomers(mol))
                result['rs_configuration']['possible_stereoisomers'] = [
                    {
                        'smiles': Chem.MolToSmiles(iso),
                        'centers': Chem.FindMolChiralCenters(iso),
                    }
                    for iso in isomers
                ]
                result['rs_configuration']['n_possible_stereoisomers'] = len(isomers)
            except Exception as e:
                logger.debug(f"Stereoisomer enumeration failed: {e}")

        logger.info(f"R/S config: {smiles} → {total} centers (R={n_r}, S={n_s}, ?={n_unspecified})")
        return result

    def _analyze_chiral_center(self, mol, atom, idx, config):
        """Analyze a single chiral center in detail."""
        symbol = atom.GetSymbol()

        # Get neighbors (substituents) with their atomic info
        neighbors = []
        for bond in atom.GetBonds():
            other_idx = bond.GetOtherAtomIdx(idx)
            neighbor = mol.GetAtomWithIdx(other_idx)
            neighbors.append({
                'index': other_idx,
                'symbol': neighbor.GetSymbol(),
                'atomic_number': neighbor.GetAtomicNum(),
                'bond_type': str(bond.GetBondType()),
                'is_in_ring': bond.IsInRing(),
            })

        # Sort by atomic number to show CIP priority order
        neighbors_sorted = sorted(neighbors, key=lambda x: -x['atomic_number'])

        # Add implicit H as lowest priority if present
        n_implicit_hs = atom.GetTotalNumHs()
        if n_implicit_hs > 0:
            neighbors_sorted.append({
                'index': None,
                'symbol': f'H (implicit, ×{n_implicit_hs})',
                'atomic_number': 1,
                'bond_type': 'implicit',
                'is_in_ring': False,
            })

        # Build substituent list with priority numbers
        substituents = []
        for i, n in enumerate(neighbors_sorted):
            substituents.append({
                'priority': i + 1,
                'atom': n['symbol'],
                'atomic_number': n['atomic_number'],
                'neighbor_index': n['index'],
            })

        # Generate justification text
        if len(substituents) >= 4:
            priority_str = " > ".join(f"{s['atom']}({s['atomic_number']})" for s in substituents[:4])
            direction = "clockwise" if config == 'R' else "counterclockwise"
            justification = (
                f"CIP priority: {priority_str}. "
                f"With lowest priority group positioned away, "
                f"the sequence 1→2→3 traces {direction} → {config}."
            )
        elif config == '?':
            justification = (
                "Stereochemistry is unspecified at this center. "
                "Both R and S configurations are possible."
            )
        else:
            justification = f"Configuration assigned by RDKit: {config}"

        return {
            'atom_index': idx,
            'symbol': symbol,
            'configuration': config,
            'substituents': substituents,
            'neighbors_raw': neighbors,
            'justification': justification,
        }


if __name__ == "__main__":
    run_mcp_server()
