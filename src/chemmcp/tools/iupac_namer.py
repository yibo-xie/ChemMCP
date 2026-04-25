"""
IUPAC系统命名生成工具
Generate systematic IUPAC name from molecular structure (SMILES).
Uses RDKit's built-in IUPAC naming with PubChem fallback.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPSearchFailError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

try:
    from ..tool_utils.names import pubchem_smiles2iupac
    _PUBCHEM_AVAILABLE = True
except ImportError:
    _PUBCHEM_AVAILABLE = False


@ChemMCPManager.register_tool
class IupacNamer(BaseTool):
    __version__ = "0.1.0"
    name = "IupacNamer"
    func_name = 'generate_iupac_name'
    description = "Generate systematic IUPAC name from a molecular structure given as SMILES."
    implementation_description = "Uses RDKit's Chem.MolToIUPACName() as primary method. Falls back to PubChem API lookup if RDKit naming fails or is unavailable. Supports organic molecules including those with stereochemistry, rings, and common functional groups."
    categories = ["Molecule"]
    tags = ["IUPAC", "Nomenclature", "SMILES", "RDKit", "PubChem"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
        ("PubChemPy", "https://github.com/mcs07/PubChemPy", "MIT"),
    ]
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
    ]
    text_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
    ]
    output_sig = [
        ('iupac_name', 'str', 'Systematic IUPAC name of the molecule.'),
        ('method_used', 'str', 'Which method produced the name: rdkit or pubchem.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'CCO'},
            'text_input': {'smiles': 'CCO'},
            'output': {'iupac_name': 'ethanol', 'method_used': 'rdkit'},
        },
        {
            'code_input': {'smiles': 'CC(=O)O'},
            'text_input': {'smiles': 'CC(=O)O'},
            'output': {'iupac_name': 'acetic acid', 'method_used': 'rdkit'},
        },
        {
            'code_input': {'smiles': 'c1ccccc1'},
            'text_input': {'smiles': 'c1ccccc1'},
            'output': {'iupac_name': 'benzene', 'method_used': 'rdkit'},
        },
    ]

    def _run_base(self, smiles: str) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available. Please install rdkit.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string into a molecule.")

        # Try RDKit IUPAC naming first
        try:
            iupac_name = Chem.MolToIUPACName(mol)
            if iupac_name and len(iupac_name.strip()) > 0:
                logger.debug(f"RDKit IUPAC naming succeeded: {iupac_name}")
                return {
                    'iupac_name': iupac_name,
                    'method_used': 'rdkit',
                }
        except Exception as e:
            logger.debug(f"RDKit IUPAC naming failed: {e}")

        # Fallback to PubChem
        if _PUBCHEM_AVAILABLE:
            try:
                iupac_name = pubchem_smiles2iupac(smiles)
                logger.debug("PubChem fallback succeeded.")
                return {
                    'iupac_name': iupac_name,
                    'method_used': 'pubchem',
                }
            except ChemMCPSearchFailError:
                logger.debug("PubChem lookup also failed.")

        raise ChemMCPSearchFailError(
            "Failed to generate IUPAC name. Both RDKit and PubChem methods failed "
            "for the given SMILES string."
        )


if __name__ == "__main__":
    run_mcp_server()
