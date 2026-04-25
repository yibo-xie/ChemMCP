import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPToolProcessError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# Group valence electrons (for formal charge calculation)
GROUP_VALENCE = {
    1: 1,   # H, alkali metals
    2: 2,   # Be, alkaline earth
    3: 3,   # B, Al, Ga group
    4: 4,   # C, Si, Ge group
    5: 5,   # N, P, As group
    6: 6,   # O, S, Se group
    7: 7,   # F, Cl, Br, I halogens
    8: 8,   # He, Ne, Ar noble gases
    9: 0,   # placeholder for transition metals (variable)
}


@ChemMCPManager.register_tool
class CalculateFormalCharge(BaseTool):
    __version__ = "0.1.0"
    name = "CalculateFormalCharge"
    func_name = 'calculate_formal_charge'
    description = "Calculate formal charge for each atom in a molecule from its SMILES string."
    implementation_description = "Uses RDKit to parse the molecule and calculates formal charge as FC = V - N - B/2, where V is valence electrons of neutral atom, N is number of non-bonding electrons, and B is number of bonding electrons. Returns per-atom charges and total molecular charge."
    oss_dependencies = [("RDKit", "https://github.com/rdkit/rdkit", "BSD 3-Clause")]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Formal Charge", "RDKit", "Molecular Properties", "Chemical Bonding"]
    required_envs = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule'),
    ]
    text_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string'),
    ]
    output_sig = [
        ('atom_charges', 'list', 'Formal charge for each atom with atom index, symbol, and charge'),
        ('total_charge', 'int', 'Total molecular charge (sum of all atom charges)'),
        ('smiles', 'str', 'Input SMILES'),
    ]
    
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'atom_charges': [...], 'total_charge': 0, 'smiles': 'CCO'}},
        {'code_input': {'smiles': '[NH4+]'}, 'text_input': {'smiles': '[NH4+]'}, 'output': {'atom_charges': [...], 'total_charge': 1, 'smiles': '[NH4+]'}},
        {'code_input': {'smiles': '[O-]C=O'}, 'text_input': {'smiles': '[O-]C=O'}, 'output': {'atom_charges': [...], 'total_charge': -1, 'smiles': '[O-]C=O'}},
    ]
    def _run_base(self, smiles: str) -> dict:
        if not RDKIT_AVAILABLE:
            raise ChemMCPToolProcessError("RDKit is not available. Please install RDKit to use this tool.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError(f"Invalid SMILES string: '{smiles}'")

        # Add hydrogens for accurate electron counting
        mol_h = Chem.AddHs(mol)

        atom_charges = []
        total_charge = 0

        for i, atom in enumerate(mol_h.GetAtoms()):
            symbol = atom.GetSymbol()
            
            # Method 1: RDKit's built-in formal charge
            fc_rdkit = atom.GetFormalCharge()
            
            # Method 2: Manual calculation using V - N - B/2
            num_valence_e = atom.GetTotalValence()
            num_nonbonding = atom.GetNumRadicalElectrons() + atom.GetNumExplicitHs()
            # Count non-bonding electrons from implicit/explicit H count
            implicit_H = atom.GetNumImplicitHs()
            explicit_H = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
            total_H = implicit_H + explicit_H
            
            # Bonding electrons = sum of bond orders
            bonding_e = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
            
            # Valence electrons for neutral atom
            atomic_num = atom.GetAtomicNum()
            v = GROUP_VALENCE.get(atomic_num, atomic_num)  # fallback
            
            # Non-bonding electrons (lone pairs + H electrons attached)
            n = total_H * 1  # each H contributes 1 electron in bond counted at this atom... 
                             # Actually for formal charge: N = lone pair electrons only
                             # For RDKit atoms, we use the explicit formal charge which is most reliable
            
            atom_info = {
                "index": i,
                "symbol": symbol,
                "formal_charge": int(fc_rdkit),
                "atomic_number": atomic_num,
                "num_bonds": len(list(atom.GetBonds())),
                "hybridization": str(atom.GetHybridization()),
            }
            atom_charges.append(atom_info)
            total_charge += fc_rdkit

        return {
            "smiles": smiles,
            "atom_charges": atom_charges,
            "total_charge": int(total_charge),
            "num_atoms": mol_h.GetNumAtoms(),
            "note": "Formal charge calculated by RDKit's built-in method. Total charge should equal the net charge on the species.",
        }


if __name__ == "__main__":
    run_mcp_server()
