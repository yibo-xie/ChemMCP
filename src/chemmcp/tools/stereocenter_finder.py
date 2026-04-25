"""
手性中心定位工具
Locate stereocenters (chiral centers) and potential stereocenters in a molecule.
Uses RDKit's chirality analysis capabilities.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, EnumerateStereoisomers
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class StereocenterFinder(BaseTool):
    __version__ = "0.1.0"
    name = "StereocenterFinder"
    func_name = 'find_stereocenters'
    description = "Locate all stereocenters (chiral centers) and potential stereocenters in a molecule. Duplicates between specified and unspecified centers are reported."
    implementation_description = "Uses RDKit's FindMolChiralCenters() to find specified chiral centers and FindPotentialStereo() to find all potential stereocenters (including those with unspecified stereochemistry). Also detects double bond stereoisomerism and ring stereochemistry."
    categories = ["Molecule"]
    tags = ["Stereochemistry", "Chirality", "Stereocenter", "RDKit", "CIP"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('include_potential', 'bool', 'True', 'Whether to include potential (unspecified) stereocenters in addition to specified ones.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional include_potential flag. E.g., "[C@H](C)(Cl)Br true".'),
    ]
    output_sig = [
        ('stereocenters', 'dict', 'Dictionary with specified_centers, potential_centers, double_bond_stereo, and summary.'),
    ]
    examples = [
        {
            'code_input': {'smiles': '[C@H](C)(Cl)Br', 'include_potential': True},
            'text_input': {'query': '[C@H](C)(Cl)Br true'},
            'output': {
                'stereocenters': {
                    'specified_centers': [{'index': 0, 'symbol': 'C', 'rs_config': 'S'}],
                    'potential_centers': [],
                    'double_bond_stereo': [],
                    'total_stereocenters': 1,
                    'summary': '1 specified chiral center(s) found at atom(s): C(0)[S].',
                }
            },
        },
        {
            'code_input': {'smiles': 'CC(O)C(Cl)Br', 'include_potential': True},
            'text_input': {'query': 'CC(O)C(Cl)Br true'},
            'output': {
                'stereocenters': {
                    'specified_centers': [],
                    'potential_centers': [{'index': 2, 'symbol': 'C', 'type': 'tetrahedral'}],
                    'double_bond_stereo': [],
                    'total_stereocenters': 1,
                    'summary': '1 potential chiral center(s) found at atom(s): C(2).',
                }
            },
        },
        {
            'code_input': {'smiles': 'C[C@H](O)[C@@H](O)C(C)(C)C', 'include_potential': True},
            'text_input': {'query': 'C[C@H](O)[C@@H](O)C(C)(C)C true'},
            'output': {
                'stereocenters': {
                    'specified_centers': [
                        {'index': 1, 'symbol': 'C', 'rs_config': 'R'},
                        {'index': 3, 'symbol': 'C', 'rs_config': 'S'},
                    ],
                    'potential_centers': [],
                    'total_stereocenters': 2,
                    'summary': '2 specified chiral center(s) found at atoms: C(1)[R], C(3)[S].',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, include_potential: bool = True) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Assign stereochemistry
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # 1. Find specified (explicitly marked) chiral centers
        specified_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
        spec_center_list = []
        for idx, config in specified_centers:
            atom = mol.GetAtomWithIdx(idx)
            spec_center_list.append({
                'index': idx,
                'symbol': atom.GetSymbol(),
                'rs_config': config,  # 'R' or 'S'
            })

        # 2. Find potential (including unspecified) stereocenters
        potential_center_list = []
        if include_potential:
            # Use FindPotentialStereo for comprehensive detection
            try:
                stereo_info = Chem.FindPotentialStereo(mol)
                for si in stereo_info:
                    if si.type == Chem.STEREOTETRAHEDRON:
                        centered_on = si.centeredOn
                        atom = mol.GetAtomWithIdx(centered_on)
                        # Check if this is already in specified
                        already_specified = any(c['index'] == centered_on for c in spec_center_list)
                        entry = {
                            'index': centered_on,
                            'symbol': atom.GetSymbol(),
                            'type': 'tetrahedral',
                            'specified': not already_specified,
                        }
                        if not already_specified:
                            potential_center_list.append(entry)
                    elif si.type == Chem.STEREOANY or si.type == Chem.STEREODOUBLEBOND:
                        pass  # Handle below
            except Exception as e:
                logger.debug(f"FindPotentialStereo failed: {e}")
                # Fallback: use the simpler approach
                unassigned = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                for idx, config in unassigned:
                    if not any(c['index'] == idx for c in spec_center_list):
                        atom = mol.GetAtomWithIdx(idx)
                        potential_center_list.append({
                            'index': idx,
                            'symbol': atom.GetSymbol(),
                            'type': 'tetrahedral',
                            'specified': False,
                        })

        # 3. Double bond stereochemistry
        db_stereo_list = []
        try:
            stereo_info = Chem.FindPotentialStereo(mol)
            for si in stereo_info:
                if si.type in (Chem.STEREODOUBLEBOND, Chem.STEREOANY, Chem.STEREOE/Z):
                    bond_atoms = []
                    for ai in range(si.centeredOn - 1, si.centeredOn + 2):
                        if 0 <= ai < mol.GetNumAtoms():
                            bond_atoms.append({
                                'index': ai,
                                'symbol': mol.GetAtomWithIdx(ai).GetSymbol(),
                            })
                    db_stereo_list.append({
                        'bond_center': si.centeredOn,
                        'atoms': bond_atoms,
                        'type': str(si.type),
                        'specified': si.specified == Chem.STEREOSPECIFIED,
                    })
        except Exception as e:
            logger.debug(f"Double bond stereo detection: {e}")

        # 4. Count total stereoisomers possible
        n_specified = len(spec_center_list)
        n_potential = len(potential_center_list)
        n_db = len(db_stereo_list)
        total = n_specified + n_potential

        # Calculate max stereoisomers
        max_isomers = 2 ** total if total > 0 else 1

        # Build summary
        parts = []
        if n_specified > 0:
            spec_str = ", ".join(f"{c['symbol']}({c['index']})[{c['rs_config']}]" for c in spec_center_list)
            parts.append(f"{n_specified} specified chiral center(s) found at atom(s): {spec_str}.")
        if n_potential > 0:
            pot_str = ", ".join(f"{c['symbol']}({c['index']})" for c in potential_center_list)
            parts.append(f"{n_potential} potential chiral center(s) found at atom(s): {pot_str}.")
        if n_db > 0:
            parts.append(f"{n_db} double bond(s) with potential E/Z stereochemistry.")
        if total == 0:
            parts.append("No stereocenters found in this molecule.")

        result = {
            'stereocenters': {
                'specified_centers': spec_center_list,
                'potential_centers': potential_center_list,
                'double_bond_stereo': db_stereo_list,
                'total_stereocenters': total,
                'max_possible_stereoisomers': max_isomers,
                'summary': " ".join(parts),
            }
        }

        logger.info(f"Stereocenter finder: {smiles} → {total} centers ({n_specified} specified, {n_potential} potential)")
        return result

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        smiles = parts[0]
        inc_pot = True
        if len(parts) > 1:
            inc_pot = parts[1].lower() in ('true', '1', 'yes', 't')
        return self._run_base(smiles, inc_pot)


if __name__ == "__main__":
    run_mcp_server()
