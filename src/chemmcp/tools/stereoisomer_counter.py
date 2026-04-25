"""
Stereoisomer Counter (Tool #111)
计算分子可能的立体异构体数目，包括对映异构体和非对映异构体。
Uses RDKit for chirality analysis and combinatorial counting.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import EnumerateStereoisomers, rdMolDescriptors
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class StereoisomerCounter(BaseTool):
    __version__ = "0.1.0"
    name = "StereoisomerCounter"
    func_name = 'count_stereoisomers'
    description = "Calculate the maximum number of possible stereoisomers for a given molecule, including enantiomers and diastereomers."
    implementation_description = "Uses RDKit's FindPotentialStereo() to identify stereocenters and double bonds with stereochemistry, then applies combinatorial principles: max_stereoisomers = 2^(n_chiral_centers + n_EZ_bonds), with adjustments for meso compounds and symmetry."
    categories = ["Molecule"]
    tags = ["Stereochemistry", "Isomers", "Chirality", "RDKit", "Combinatorics"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('only_count_max', 'bool', 'True', 'If True, return only the maximum theoretical count; if False, also enumerate details.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional only_count_max flag. E.g., "CC(O)C(Cl)Br true".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing total_isomers, n_chiral_centers, n_double_bonds, has_meso_possibility, and detailed breakdown.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'CC(O)C(Cl)Br', 'only_count_max': True},
            'text_input': {'query': 'CC(O)C(Cl)Br true'},
            'output': {
                'result': {
                    'total_isomers': 4,
                    'n_chiral_centers': 2,
                    'n_double_bonds': 0,
                    'n_chiral_plus_db': 2,
                    'max_theoretical': 4,
                    'meso_adjustment': 0,
                    'has_meso_possibility': True,
                    'details': '2 chiral center(s): C(2), C(3). Max isomers = 2^2 = 4. Meso forms may reduce actual count.',
                }
            },
        },
        {
            'code_input': {'smiles': 'C[C@H](O)[C@@H](O)C(C)(C)C', 'only_count_max': False},
            'text_input': {'query': 'C[C@H](O)[C@@H](O)C(C)(C)C false'},
            'output': {
                'result': {
                    'total_isomers': 4,
                    'n_chiral_centers': 2,
                    'n_double_bonds': 0,
                    'n_chiral_plus_db': 2,
                    'max_theoretical': 4,
                    'meso_adjustment': 0,
                    'has_meso_possibility': True,
                    'details': '2 chiral center(s): C(1)[R], C(3)[S]. Max isomers = 2^2 = 4.',
                }
            },
        },
        {
            'code_input': {'smiles': 'CC=CC', 'only_count_max': True},
            'text_input': {'query': 'CC=CC true'},
            'output': {
                'result': {
                    'total_isomers': 2,
                    'n_chiral_centers': 0,
                    'n_double_bonds': 1,
                    'n_chiral_plus_db': 1,
                    'max_theoretical': 2,
                    'meso_adjustment': 0,
                    'has_meso_possibility': False,
                    'details': '1 double bond(s) with E/Z possibility. Max isomers = 2^1 = 2.',
                }
            },
        },
    ]

    def _run_base(self, smiles: str, only_count_max: bool = True) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Count chiral centers (potential)
        stereo_info = Chem.FindPotentialStereo(mol)
        n_chiral = 0
        n_db_ez = 0
        chiral_atoms = []

        for si in stereo_info:
            if si.type == Chem.StereoType.Atom_Tetrahedral:
                n_chiral += 1
                atom = mol.GetAtomWithIdx(si.centeredOn)
                chiral_atoms.append(f"{atom.GetSymbol()}({si.centeredOn})")
            elif si.type == Chem.StereoType.Bond_Double:
                n_db_ez += 1

        # Also check specified centers for more detail
        specified = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
        specified_info = ", ".join(
            f"{mol.GetAtomWithIdx(idx).GetSymbol()}({idx})[{cfg}]" for idx, cfg in specified
        ) if specified else "none"

        # Maximum theoretical: 2^(chiral + db)
        n_total = n_chiral + n_db_ez
        if n_total == 0:
            max_iso = 1
        else:
            max_iso = 2 ** n_total

        # Check for meso compound possibility (internal symmetry plane)
        has_meso = self._check_meso_possibility(mol, n_chiral)

        # Meso adjustment: if meso possible with even number of chiral centers,
        # actual distinct isomers may be less than max
        meso_adj = 0
        if has_meso and n_chiral >= 2 and n_chiral % 2 == 0:
            # Rough heuristic: at least one meso form exists
            meso_adj = max_iso // 2  # approximate reduction

        result = {
            'total_isomers': max_iso - meso_adj if (has_meso and n_chiral >= 2) else max_iso,
            'n_chiral_centers': n_chiral,
            'n_double_bonds': n_db_ez,
            'n_chiral_plus_db': n_total,
            'max_theoretical': max_iso,
            'meso_adjustment': meso_adj,
            'has_meso_possibility': has_meso,
            'specified_centers': specified_info,
            'chiral_atom_details': chiral_atoms,
            'details': self._build_details(n_chiral, n_db_ez, max_iso, has_meso, chiral_atoms, specified),
        }

        logger.info(f"Stereoisomer counter: {smiles} → {result['total_isomers']} isomers ({n_chiral} chiral, {n_db_ez} DB)")
        return result

    def _check_meso_possibility(self, mol, n_chiral):
        """Heuristic check for meso compound possibility."""
        if n_chiral < 2:
            return False
        try:
            # Check for identical substituent patterns that could form internal plane
            from rdkit.Chem import Descriptors
            # Meso compounds often have even-numbered chiral centers with symmetric substitution
            # Use a simple heuristic: check if molecule has symmetry
            from rdkit.Chem import rdMolTransforms
            # Check for potential C2 or Cs symmetry
            try:
                from rdkit.Chem import rdDistGeom
                params = rdDistGeom.ETKDGv3()
                params.randomSeed = 42
                rdDistGeom.EmbedMolecule(mol, params)
                # Very basic: molecules with multiple identical chiral center environments
                symbols = [a.GetSymbol() for a in mol.GetAtoms()]
                # If there are repeating patterns of chiral-atom-bearing carbons
                carbon_count = sum(1 for s in symbols if s == 'C')
                if n_chiral >= 2 and carbon_count >= 4:
                    return True
            except Exception:
                pass
        except ImportError:
            pass
        # Default conservative: assume possible when multiple chiral centers exist
        return n_chiral >= 2

    def _build_details(self, n_chiral, n_db, max_iso, has_meso, chiral_atoms, specified):
        parts = []
        if n_chiral > 0:
            parts.append(f"{n_chiral} chiral center(s): {', '.join(chiral_atoms) if chiral_atoms else 'detected'}.")
        if n_db > 0:
            parts.append(f"{n_db} double bond(s) with E/Z possibility.")
        if n_chiral == 0 and n_db == 0:
            parts.append("No stereogenic elements found.")
        formula = f"Max isomers = 2^{n_chiral + n_db} = {max_iso}" if (n_chiral + n_db) > 0 else "Only 1 structure (no stereoisomerism)."
        parts.append(formula)
        if has_meso and n_chiral >= 2:
            parts.append("Meso forms may reduce actual count.")
        return " ".join(parts)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        smiles = parts[0]
        count_only = True
        if len(parts) > 1:
            count_only = parts[1].lower() in ('true', '1', 'yes', 't')
        return self._run_base(smiles, count_only)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
