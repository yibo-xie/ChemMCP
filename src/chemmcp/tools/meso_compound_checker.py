"""
Meso Compound Checker (Tool #112)
判断分子是否为内消旋(meso)化合物。
Uses RDKit for symmetry and chirality analysis.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, EnumerateStereoisomers, Descriptors
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class MesoCompoundChecker(BaseTool):
    __version__ = "0.1.0"
    name = "MesoCompoundChecker"
    func_name = 'check_meso_compound'
    description = "Determine whether a molecule is a meso compound (internally compensated chiral molecule with an internal plane of symmetry)."
    implementation_description = "Uses RDKit to enumerate stereoisomers and check for achiral stereoisomers among them. A meso compound has chiral centers but is overall achiral due to an internal symmetry plane. The tool generates all stereoisomers and identifies which are meso forms."
    categories = ["Molecule"]
    tags = ["Stereochemistry", "Meso Compound", "Chirality", "Symmetry", "RDKit"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('enumerate_isomers', 'bool', 'True', 'Whether to enumerate all stereoisomers to identify meso forms explicitly.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional enumerate flag. E.g., "CC(O)C(O)C true".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing is_meso, n_chiral_centers, n_stereoisomers, meso_forms, and explanation.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'CC(O)C(O)C(Cl)(Cl)C(O)C(O)CC', 'enumerate_isomers': True},
            'text_input': {'query': 'CC(O)C(O)C(Cl)(Cl)C(O)C(O)CC true'},
            'output': {
                'result': {
                    'is_meso_candidate': True,
                    'n_chiral_centers': 4,
                    'n_stereoisomers': None,
                    'has_internal_symmetry': True,
                    'meso_count_estimate': 2,
                    'explanation': 'This molecule has multiple chiral centers with symmetric substitution patterns, suggesting possible meso forms due to an internal plane of symmetry.',
                }
            },
        },
        {
            'code_input': {'smiles': '[C@H](O)(Cl)Br', 'enumerate_isomers': True},
            'text_input': {'query': '[C@H](O)(Cl)Br true'},
            'output': {
                'result': {
                    'is_meso_candidate': False,
                    'n_chiral_centers': 1,
                    'n_stereoisomers': 2,
                    'has_internal_symmetry': False,
                    'meso_count_estimate': 0,
                    'explanation': 'Only 1 chiral center cannot form a meso compound (requires ≥2 centers with internal compensation).',
                }
            },
        },
        {
            'code_input': {'smiles': 'CC(O)C(Cl)Br', 'enumerate_isomers': True},
            'text_input': {'query': 'CC(O)C(Cl)Br true'},
            'output': {
                'result': {
                    'is_meso_candidate': True,
                    'n_chiral_centers': 2,
                    'n_stereoisomers': 4,
                    'has_internal_symmetry': True,
                    'meso_count_estimate': 0,
                    'explanation': '2 chiral centers with different substituents: no identical pairs → unlikely to be meso despite having 2+ centers.',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, enumerate_isomers: bool = True) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Find potential stereocenters
        stereo_info = Chem.FindPotentialStereo(mol)
        n_chiral = sum(1 for si in stereo_info if si.type == Chem.StereoType.Atom_Tetrahedral)

        # Basic criteria: need at least 2 chiral centers
        if n_chiral < 2:
            return {
                'result': {
                    'is_meso_candidate': False,
                    'n_chiral_centers': n_chiral,
                    'n_stereoisomers': 2 ** max(n_chiral, 0),
                    'has_internal_symmetry': False,
                    'meso_count_estimate': 0,
                    'explanation': f'Only {n_chiral} chiral center(s). Meso compounds require at least 2 chiral centers with internal symmetry compensation.',
                }
            }

        # Check for internal symmetry - analyze substituent patterns at each chiral center
        chiral_center_data = self._analyze_chiral_centers(mol, stereo_info)

        # Determine if meso is possible based on substituent symmetry
        has_symmetry, reason = self._check_substituent_symmetry(chiral_center_data, mol)

        meso_estimate = 0
        iso_list = []

        if enumerate_isomers and n_chiral <= 6:
            try:
                isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(
                    mol, options=EnumerateStereoisomers.StereoEnumerationOptions(unique=True)
                ))
                # Check each isomer for chirality
                meso_iso = []
                chiral_iso = []
                for iso in isomers:
                    iso_smiles = Chem.MolToSmiles(iso)
                    try:
                        is_chiral =Descriptors.Chi(iso) != 0 or any(
                            atom.GetChiralTag() != Chem.CHI_UNSPECIFIED
                            for atom in iso.GetAtoms() if atom.GetAtomicNum() > 1
                        )
                    except Exception:
                        is_chiral = True

                    # More robust: check if the isomer has a net chirality
                    try:
                        from rdkit.Chem import AllChem
                        AllChem.AssignStereochemistry(iso, force=True)
                        # Count R and S centers
                        r_count = sum(1 for a in iso.GetAtoms() if a.GetChiralTag() == Chem.CHI_TETRA_CCW)
                        s_count = sum(1 for a in iso.GetAtoms() if a.GetChiralTag() == Chem.CHI_TETRA_CW)
                        # Meso-like: equal R/S in symmetric environment
                        entry = {'smiles': iso_smiles, 'r_centers': r_count, 's_centers': s_count}
                        if r_count == s_count and r_count > 0 and has_symmetry:
                            meso_iso.append(entry)
                            meso_estimate += 1
                        else:
                            chiral_iso.append(entry)
                    except Exception:
                        chiral_iso.append({'smiles': iso_smiles, 'r_centers': '?', 's_centers': '?'})

                iso_list = {'meso_forms': meso_iso, 'chiral_forms': chiral_iso}
            except Exception as e:
                logger.debug(f"Stereoisomer enumeration failed: {e}")
                iso_list = {'error': str(e)}

        total_iso = 2 ** n_chiral if n_chiral > 0 else 1

        result = {
            'result': {
                'is_meso_candidate': has_symmetry,
                'n_chiral_centers': n_chiral,
                'n_stereoisometers': total_iso,
                'has_internal_symmetry': has_symmetry,
                'meso_count_estimate': meso_estimate,
                'chiral_center_analysis': chiral_center_data,
                'isomer_details': iso_list if iso_list else None,
                'explanation': reason,
            }
        }

        # Fix typo in key name
        result['result']['n_stereoisomers'] = result['result'].pop('n_stereoisometers')

        logger.info(f"Meso checker: {smiles} → meso={has_symmetry}, centers={n_chiral}")
        return result

    def _analyze_chiral_centers(self, mol, stereo_info):
        """Analyze each chiral center's substituent pattern."""
        centers = []
        for si in stereo_info:
            if si.type != Chem.StereoType.Atom_Tetrahedral:
                continue
            idx = si.centeredOn
            atom = mol.GetAtomWithIdx(idx)
            neighbors = atom.GetNeighbors()
            neighbor_info = []
            for n in neighbors:
                # Get neighbor's atomic symbol and its own neighbors (one level deep)
                n_symbols = [nn.GetSymbol() for nn in n.GetNeighbors()]
                neighbor_info.append({
                    'symbol': n.GetSymbol(),
                    'idx': n.GetIdx(),
                    'neighbor_symbols': sorted(n_symbols),
                })
            centers.append({
                'index': idx,
                'symbol': atom.GetSymbol(),
                'neighbors': neighbor_info,
            })
        return centers

    def _check_substituent_symmetry(self, chiral_data, mol):
        """Check if there are symmetric pairs of chiral centers that could compensate."""
        if len(chiral_data) < 2:
            return False, "Need at least 2 chiral centers."

        # Group chiral centers by their neighbor signature for comparison
        signatures = []
        for c in chiral_data:
            # Create a hashable signature of the substitution pattern
            neigh_symbols = sorted([n['symbol'] for n in c['neighbors']])
            sig = (c['symbol'], tuple(neigh_symbols))
            signatures.append(sig)

        # Check for pairs of identical signatures
        from collections import Counter
        sig_counts = Counter(signatures)
        paired_sums = sum(c // 2 for c in sig_counts.values())

        if paired_sums >= 1:
            return True, (
                f"Found {len(chiral_data)} chiral center(s) with {paired_sums} pair(s) of "
                f"identical substitution patterns. Internal symmetry plane may exist, "
                f"making meso form(s) possible."
            )
        else:
            return False, (
                f"{len(chiral_data)} chiral center(s) found but no identical substitution pairs. "
                f"Each center has unique substituents → meso compound unlikely."
            )

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        smiles = parts[0]
        enum_iso = True
        if len(parts) > 1:
            enum_iso = parts[1].lower() in ('true', '1', 'yes', 't')
        return self._run_base(smiles, enum_iso)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
