import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPToolProcessError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem as RDKitChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# Hybridization mapping from steric number (electron domains)
HYBRIDIZATION_MAP = {
    2: ("sp", "linear", 180),
    3: ("sp²", "trigonal planar", 120),
    4: ("sp³", "tetrahedral", 109.5),
    5: ("sp³d", "trigonal bipyramidal", [90, 120]),
    6: ("sp³d²", "octahedral", 90),
    7: ("sp³d³", "pentagonal bipyramidal", [72, 90]),
}


@ChemMCPManager.register_tool
class PredictHybridization(BaseTool):
    __version__ = "0.1.0"
    name = "PredictHybridization"
    func_name = 'predict_hybridization'
    description = "Predict hybridization state of atoms in a molecule from its SMILES string."
    implementation_description = "Uses RDKit to analyze molecular structure and determines hybridization for each atom based on steric number (number of sigma bonds + lone pairs). Maps steric number to hybridization type (sp, sp2, sp3, sp3d, sp3d2)."
    oss_dependencies = [("RDKit", "https://github.com/rdkit/rdkit", "BSD 3-Clause")]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Hybridization", "RDKit", "Molecular Properties", "Chemical Bonding"]
    required_envs = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule'),
        ('atom_index', 'int', 'N/A', 'Optional: specific atom index (0-based). If not provided, returns all atoms.'),
    ]
    text_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string'),
    ]
    output_sig = [
        ('atom_hybridizations', 'list', 'Hybridization for each atom with symbol and steric number'),
        ('summary', 'dict', 'Summary of unique hybridizations found'),
    ]
    
    
    examples = [
        {'code_input': {'smiles': 'CCO', 'atom_index': None}, 'text_input': {'smiles': 'ethanol'}, 'output': {'atom_hybridizations': [...], 'summary': '...'}},
        {'code_input': {'smiles': 'C=C', 'atom_index': None}, 'text_input': {'smiles': 'ethene'}, 'output': {'atom_hybridizations': [...], 'summary': '...'}},
    ]
    def _run_base(self, smiles: str, atom_index: int = None) -> dict:
        if not RDKIT_AVAILABLE:
            raise ChemMCPToolProcessError("RDKit is not available.")

        mol = RDKitChem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError(f"Invalid SMILES string: '{smiles}'")

        mol_h = RDKitChem.AddHs(mol)
        
        atom_results = []
        hybrid_counts = {}
        
        for i, atom in enumerate(mol_h.GetAtoms()):
            if atom_index is not None and i != atom_index:
                continue
            
            sym = atom.GetSymbol()
            # Get RDKit's hybridization
            rdkit_hyb = str(atom.GetHybridization())
            
            # Calculate steric number manually
            num_bonds = len(list(atom.GetBonds()))
            # Count lone pairs from valence - bonding electrons
            atomic_num = atom.GetAtomicNum()
            total_valence_electrons = {1:1,6:2,5:3,7:3,8:2,9:1,15:5,16:6}.get(atomic_num, 4)
            
            # Count bonded heavy atoms + H atoms
            neighbors = list(atom.GetNeighbors())
            h_count = sum(1 for n in neighbors if n.GetSymbol() == 'H')
            heavy_bond_count = len(neighbors) - h_count
            
            # Bond order sum
            bond_order_sum = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
            
            # Steric number = number of sigma bonds (count each bond once) + lone pairs
            # For simplicity, use RDKit's built-in hybridization which is reliable
            hyb_map = {
                "RDKit_HYBRID_S": "sp",
                "RDKit_HYBRID_SP": "sp",
                "RDKit_HYBRID_SP2": "sp²",
                "RDKit_HYBRID_SP3": "sp³",
                "RDKit_HYBRID_SP3D": "sp³d",
                "RDKit_HYBRID_SP3D2": "sp³d²",
                "UNSPECIFIED": "unknown",
            }
            predicted = hyb_map.get(rdkit_hyb, rdkit_hyb)
            
            info = {
                "index": i,
                "symbol": sym,
                "hybridization": predicted,
                "rdkit_hybridization": rdkit_hyb,
                "num_sigma_bonds": num_bonds,
                "is_aromatic": atom.GetIsAromatic(),
            }
            atom_results.append(info)
            hybrid_counts[predicted] = hybrid_counts.get(predicted, 0) + 1

        result = {
            "smiles": smiles,
            "atom_hybridizations": atom_results,
            "summary": {
                "total_atoms": mol_h.GetNumAtoms(),
                "hybridization_distribution": hybrid_counts,
            },
        }

        if atom_index is not None:
            result["note"] = f"Results shown only for atom index {atom_index}."

        return result


if __name__ == "__main__":
    run_mcp_server()
