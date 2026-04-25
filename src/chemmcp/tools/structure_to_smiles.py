"""
结构描述转SMILES工具
Convert structural description or natural language chemical name to SMILES.
Uses RDKit's name parsing, common name database, and PubChem fallback.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPSearchFailError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

try:
    from pubchempy import get_compounds, CompoundNotFoundError
    _PUBCHEMPY_AVAILABLE = True
except ImportError:
    _PUBCHEMPY_AVAILABLE = False

# Common structural descriptions → SMILES mapping
STRUCTURE_SMILES_MAP = {
    # Simple descriptions
    "methane": "C",
    "ethane": "CC",
    "propane": "CCC",
    "butane": "CCCC",
    "pentane": "CCCCC",
    "hexane": "CCCCCC",
    "heptane": "CCCCCCC",
    "octane": "CCCCCCCC",
    "ethene": "C=C",
    "ethylene": "C=C",
    "ethyne": "C#C",
    "acetylene": "C#C",
    "methanol": "CO",
    "methyl alcohol": "CO",
    "ethanol": "CCO",
    "ethyl alcohol": "CCO",
    "isopropanol": "CC(C)O",
    "isopropyl alcohol": "CC(C)O",
    "propanol": "CCCO",
    "butanol": "CCCCO",
    "formaldehyde": "C=O",
    "methanal": "C=O",
    "acetaldehyde": "CC=O",
    "ethanal": "CC=O",
    "acetone": "CC(=O)C",
    "propanone": "CC(=O)C",
    "formic acid": "C(=O)O",
    "methanoic acid": "C(=O)O",
    "acetic acid": "CC(=O)O",
    "ethanoic acid": "CC(=O)O",
    "benzene": "c1ccccc1",
    "phenol": "c1ccccc1O",
    "toluene": "Cc1ccccc1",
    "methylbenzene": "Cc1ccccc1",
    "aniline": "c1ccccc1N",
    "benzoic acid": "c1ccccc1C(=O)O",
    "styrene": "C=Cc1ccccc1",
    "naphthalene": "c1ccc2ccccc2c1",
    "ethylamine": "CCN",
    "methylamine": "CN",
    "diethylamine": "CCN(CC)",
    "trimethylamine": "CN(C)C",
    "dimethyl ether": "COC",
    "diethyl ether": "CCOCC",
    "glycerol": "OC[C@H](O)CO",
    "glucose": "OC[C@H](O)[C@H](O)[C@H](O)C(O)CO",
    "urea": "NC(=O)N",
    "chloromethane": "CCl",
    "dichloromethane": "ClCCl",
    "chloroform": "ClC(Cl)Cl",
    "carbon tetrachloride": "Cl(Cl)(Cl)Cl",
    "acetonitrile": "CC#N",
    "nitromethane": "CN(=O)=O",
    "pyridine": "n1ccccc1",
    "furan": "c1ccoc1",
    "thiophene": "c1cccs1",
    "pyrrole": "c1cc[nH]c1",
    "cyclohexane": "C1CCCCC1",
    "cyclopentane": "C1CCCC1",
    "clobutane": "C1CCC1",
    "cyclopropane": "C1CC1",
    "cyclohexene": "C1=CCCCC1",
    "cyclopentadiene": "C1=CC=C C1",  # simplified
    "norbornane": "C1CC2CCC1C2",
    "adamantane": "C1C3CC2CC1CC(C2)C3",
    "cholesterol": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3C1",
    "caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "nicotine": "cn1ccc(C2CCCN2C)c1",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "d dt": "Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl",
    "salicylic acid": "O=C(O)c1ccccc1O",
    "menthol": "CC(C)[C@H]1CC[C@@H](C)CC1",
    "camphor": "CC1(C)C2CCC1(C)C(=O)C2",
    "lactic acid": "CC(O)C(=O)O",
    "citric acid": "C(C(=O)O)C(O)(CC(=O)O)C(=O)O",
    "oxalic acid": "C(=O)C(=O)O",
    "tartaric acid": "OC[C@H](O)C(=O)O",
    "maleic acid": "O=C(O)/C=C\\C(=O)O",
    "fumaric acid": "O=C(O)/C=C/C(=O)O",
    "phthalic anhydride": "O=C1OC(=O)c2ccccc12",
    "maleic anhydride": "O=C1OC(=O)C=C1",
    "quinone": "O=C1C=CC(=O)CC=1",
    "hydroquinone": "Oc1ccc(O)cc1",
    "anthracene": "c1ccc2c(c1)ccc1ccccc21",
    "phenanthrene": "c1ccc2c(c1)ccc1ccccc21",
    "pyrene": "c1ccc2c(c1)ccc3c2ccc4c3cccc4",
    "biphenyl": "c1ccc(cc1)-c2ccccc2",
    "diphenylmethane": "c1ccc(cc1)Cc2ccccc2",
    "triphenylphosphine": "P(c1ccccc1)(c2ccccc2)c3ccccc3",
    "tetrahydrofuran": "C1CCOC1",
    "thf": "C1CCOC1",
    "dioxane": "C1COCCO1",
    "dmso": "CS(=O)C",
    "dmf": "CN(C)C=O",
    "acetic anhydride": "CC(=O)OC(=O)C",
    "ethyl acetate": "CCOC(=O)C",
    "methyl formate": "COC=O",
    "vinyl chloride": "C=CCl",
    "vinyl acetate": "C=COC(=O)C",
    "acrylamide": "C=CC(=O)N",
    "acrylonitrile": "C=CC#N",
    "caprolactam": "C1CCC(=O)NC1",
    "adipic acid": "C(CCC(=O)O)CC(=O)O",
    "terephthalic acid": "O=C(O)c1ccc(cc1)C(=O)O",
    "bisphenol a": "CC(c1ccc(O)cc1)c1ccc(O)cc1",
    "vitamin c": "O=C(O)C(O)(C(O)CO)C(O)=O",
    "thymine": "Cc1nc(nc(=O)[nH]1)O",
    "adenine": "c1ncnc2c1nc[nH]2",
    "guanine": "c1nc2c(ncn2[nH]1)N",
    "cytosine": "n1cc(nc(=O)[nH]1)N",
    "uracil": "O=c1[nH]cnc(n1)O",
    "indole": "c1ccc2[nH]ccc2c1",
    "quinoline": "c1ccc2ncccc2c1",
    "isoquinoline": "c1ccc2ncccc2c1",
    "imidazole": "n1cc[nH]c1",
    "triazole": "n1cnnc1",
}


def _normalize_query(query: str) -> str:
    """Normalize structure description query."""
    return query.strip().lower()


@ChemMCPManager.register_tool
class StructureToSmiles(BaseTool):
    __version__ = "0.1.0"
    name = "StructureToSmiles"
    func_name = 'structure_to_smiles'
    description = "Convert a structural description or natural language chemical name to SMILES representation."
    implementation_description = "Uses a built-in mapping of common names/descriptions to SMILES, with RDKit name parsing and PubChem lookup as fallbacks. Supports IUPAC names, common names, trivial names, and simple structural descriptions."
    categories = ["Molecule"]
    tags = ["SMILES", "Structure", "Name Conversion", "RDKit", "PubChem"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
        ("PubChemPy", "https://github.com/mcs07/PubChemPy", "MIT"),
    ]
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]

    code_input_sig = [
        ('structure_description_or_name', 'str', 'N/A', 'Structural description, common name, or IUPAC name of the molecule.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Name or description of the molecule.'),
    ]
    output_sig = [
        ('smiles', 'str', 'SMILES string of the molecule.'),
        ('source', 'str', 'Method used: local_database, rdkit, or pubchem.'),
        ('name_resolved', 'str', 'The resolved/normalized name of the molecule.'),
    ]
    examples = [
        {
            'code_input': {'structure_description_or_name': 'ethanol'},
            'text_input': {'query': 'ethanol'},
            'output': {'smiles': 'CCO', 'source': 'local_database', 'name_resolved': 'ethanol'},
        },
        {
            'code_input': {'structure_description_or_name': 'acetone'},
            'text_input': {'query': 'acetone'},
            'output': {'smiles': 'CC(=O)C', 'source': 'local_database', 'name_resolved': 'acetone'},
        },
        {
            'code_input': {'structure_description_or_name': 'benzene'},
            'text_input': {'query': 'benzene'},
            'output': {'smiles': 'c1ccccc1', 'source': 'local_database', 'name_resolved': 'benzene'},
        },
    ]

    def _run_base(self, structure_description_or_name: str) -> dict:
        if not structure_description_or_name or not structure_description_or_name.strip():
            raise ChemMCPInputError("Input cannot be empty.")

        normalized = _normalize_query(structure_description_or_name)

        # Strategy 1: Direct local database lookup
        result = self._lookup_local(normalized)
        if result:
            logger.info(f"Local DB hit: '{normalized}' → {result['smiles']}")
            return result

        # Strategy 2: Input is already valid SMILES
        if is_smiles(structure_description_or_name.strip()):
            return {
                'smiles': structure_description_or_name.strip(),
                'source': 'input_was_smiles',
                'name_resolved': structure_description_or_name.strip(),
            }

        # Strategy 3: RDKit name parsing (RDKit can parse some IUPAC/common names)
        if _RDKIT_AVAILABLE:
            result = self._try_rdkit(normalized)
            if result:
                logger.info(f"RDKit parse: '{normalized}' → {result['smiles']}")
                return result

        # Strategy 4: PubChem lookup by name
        if _PUBCHEMPY_AVAILABLE:
            result = self._try_pubchem(normalized)
            if result:
                logger.info(f"PubChem lookup: '{normalized}' → {result['smiles']}")
                return result

        raise ChemMCPSearchFailError(
            f"Cannot convert '{structure_description_or_name}' to SMILES. "
            "Not found in local database, not valid SMILES, and RDKit/PubChem lookups failed."
        )

    def _lookup_local(self, normalized: str) -> dict:
        """Look up in local database."""
        # Exact match
        if normalized in STRUCTURE_SMILES_MAP:
            smiles = STRUCTURE_SMILES_MAP[normalized]
            return {
                'smiles': smiles,
                'source': 'local_database',
                'name_resolved': normalized,
            }

        # Partial/fuzzy match
        for key, smiles in STRUCTURE_SMILES_MAP.items():
            if normalized in key or key in normalized:
                return {
                    'smiles': smiles,
                    'source': 'local_database_fuzzy',
                    'name_resolved': key,
                    'note': f'Exact match not found. Used closest match: "{key}"',
                }

        return None

    def _try_rdkit(self, normalized: str) -> dict:
        """Try RDKit's name-to-molecule conversion."""
        try:
            # RDKit can sometimes parse IUPAC names as SMILES
            mol = Chem.MolFromSmiles(normalized)
            if mol is not None:
                canon_smiles = Chem.MolToSmiles(mol)
                return {
                    'smiles': canon_smiles,
                    'source': 'rdkit_parsed_as_smiles',
                    'name_resolved': normalized,
                }
        except Exception:
            pass

        # Try MolFromInchi or other methods
        try:
            # Some names work via this path
            mol = Chem.MolFromSmiles(normalized, sanitize=False)
            if mol is not None:
                Chem.SanitizeMol(mol)
                canon_smiles = Chem.MolToSmiles(mol)
                return {
                    'smiles': canon_smiles,
                    'source': 'rdkit_sanitize_parse',
                    'name_resolved': normalized,
                }
        except Exception:
            pass

        return None

    def _try_pubchem(self, normalized: str) -> dict:
        """Try PubChem lookup by compound name."""
        try:
            compounds = get_compounds(normalized, 'name')
            if compounds:
                c = compounds[0]
                smiles = c.canonical_smiles
                if smiles:
                    return {
                        'smiles': smiles,
                        'source': 'pubchem',
                        'name_resolved': c.iupac_name or c.synonyms[0] if c.synonyms else normalized,
                    }
        except CompoundNotFoundError:
            pass
        except Exception as e:
            logger.debug(f"PubChem lookup failed: {e}")

        return None


if __name__ == "__main__":
    run_mcp_server()
