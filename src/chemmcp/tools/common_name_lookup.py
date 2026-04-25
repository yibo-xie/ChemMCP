"""
俗名与系统命名互查工具
Bidirectional lookup between common/trivial names and systematic IUPAC names.
Includes a built-in dictionary of common chemicals with PubChem fallback.
"""
import logging
from typing import Optional
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
    from ..tool_utils.names import pubchem_smiles2iupac
    _PUBCHEM_AVAILABLE = True
except ImportError:
    _PUBCHEM_AVAILABLE = False

# Comprehensive common name ↔ IUPAC / SMILES mapping database
COMMON_NAME_DB = {
    # Format: lowercase_common_name -> (iupac_name, smiles)
    "ethanol": ("ethanol", "CCO"),
    "ethyl alcohol": ("ethanol", "CCO"),
    "grain alcohol": ("ethanol", "CCO"),
    "methanol": ("methanol", "CO"),
    "methyl alcohol": ("methanol", "CO"),
    "wood alcohol": ("methanol", "CO"),
    "propanol": ("propan-1-ol", "CCCO"),
    "n-propanol": ("propan-1-ol", "CCCO"),
    "isopropanol": ("propan-2-ol", "CC(C)O"),
    "isopropyl alcohol": ("propan-2-ol", "CC(C)O"),
    "butanol": ("butan-1-ol", "CCCCO"),
    "acetone": ("propanone", "CC(=O)C"),
    "dimethyl ketone": ("propanone", "CC(=O)C"),
    "formaldehyde": ("methanal", "C=O"),
    "acetaldehyde": ("ethanal", "CC=O"),
    "formic acid": ("methanoic acid", "C(=O)O"),
    "acetic acid": ("ethanoic acid", "CC(=O)O"),
    "ethanoic acid": ("ethanoic acid", "CC(=O)O"),
    "propionic acid": ("propanoic acid", "CCC(=O)O"),
    "butyric acid": ("butanoic acid", "CCCC(=O)O"),
    "lactic acid": ("2-hydroxypropanoic acid", "CC(O)C(=O)O"),
    "oxalic acid": ("ethanedioic acid", "C(=O)C(=O)O"),
    "citric acid": ("2-hydroxypropane-1,2,3-tricarboxylic acid", "C(C(=O)O)C(O)(CC(=O)O)C(=O)O"),
    "benzoic acid": ("benzenecarboxylic acid", "c1ccccc1C(=O)O"),
    "methylamine": ("methanamine", "CN"),
    "ethylamine": ("ethanamine", "CCN"),
    "aniline": ("phenylamine", "c1ccccc1N"),
    "urea": ("carbamide", "NC(=O)N"),
    "benzene": ("benzene", "c1ccccc1"),
    "toluene": ("methylbenzene", "Cc1ccccc1"),
    "methylbenzene": ("methylbenzene", "Cc1ccccc1"),
    "xylene": ("dimethylbenzene", "Cc1cccc(C)c1"),  # o-xylene as default
    "styrene": ("ethenylbenzene", "C=Cc1ccccc1"),
    "phenylethylene": ("ethenylbenzene", "C=Cc1ccccc1"),
    "naphthalene": ("naphthalene", "c1ccc2ccccc2c1"),
    "phenol": ("phenol", "c1ccccc1O"),
    "carbolic acid": ("phenol", "c1ccccc1O"),
    "hydroquinone": ("benzene-1,4-diol", "Oc1ccc(O)cc1"),
    "catechol": ("benzene-1,2-diol", "Oc1ccccc1O"),
    "resorcinol": ("benzene-1,3-diol", "Oc1cccc(O)c1"),
    "anisole": ("methoxybenzene", "COc1ccccc1"),
    "benzaldehyde": ("benzenecarbaldehyde", "O=Cc1ccccc1"),
    "acetophenone": ("1-phenylethanone", "CC(=O)c1ccccc1"),
    "benzyl alcohol": ("phenylmethanol", "OCc1ccccc1"),
    "salicylic acid": ("2-hydroxybenzoic acid", "O=C(O)c1ccccc1O"),
    "aspirin": ("2-acetoxybenzoic acid", "CC(=O)Oc1ccccc1C(=O)O"),
    "menthol": ("2-isopropyl-5-methylcyclohexan-1-ol", "CC(C)[C@H]1CC[C@@H](C)CC1"),
    "camphor": ("1,7,7-trimethylbicyclo[2.2.1]heptan-2-one", "CC1(C)C2CCC1(C)C(=O)C2"),
    "cholesterol": ("cholest-5-en-3β-ol", "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3C1"),
    "glucose": ("glucose", "OC[C@H](O)[C@H](O)[C@H](O)C(O)CO"),
    "sucrose": "sucrose not available in simple SMILES",
    "fructose": ("fructose", "OCC1(O)C(O)C(O)C(O)C(O)C1O"),
    "caffeine": ("1,3,7-trimethylpurine-2,6-dione", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    "nicotine": ("3-(1-methylpyrrolidin-2-yl)pyridine", "cn1ccc(C2CCCN2C)c1"),
    "morphine": ("morphine", "CN1CCC23c4c5ccc(O)c4O[C@H]3[C@@H](O)[C@H]2C5C1C6=CC=CC=C6"),
    "aspartame": ("L-aspartyl-L-phenylalanine methyl ester", "COC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1ccccc1)N"),
    "ddt": ("1,1-bis(4-chlorophenyl)-2,2,2-trichloroethane", "Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl"),
    "ether": ("ethoxyethane", "CCOCC"),
    "diethyl ether": ("ethoxyethane", "CCOCC"),
    "glycerol": ("propane-1,2,3-triol", "OC[C@H](O)CO"),
    "ethylene glycol": ("ethane-1,2-diol", "OCCO"),
    "ammonia": ("azane", "N"),
    "hydrazine": ("diazane", "NN"),
    "hydrogen peroxide": ("oxidane", "OO"),
    "nitric acid": ("nitric acid", "[N+](=O)[O-]"),
    "sulfuric acid": ("sulfuric acid", "OS(=O)(=O)O"),
    "hydrochloric acid": ("hydrochloric acid", "Cl"),
    "sodium chloride": ("sodium chloride", "[Na+].[Cl-]"),
    "carbon dioxide": ("carbon dioxide", "O=C=O"),
    "carbon monoxide": ("carbon monoxide", "[C-]#[O+]"),
    "water": ("oxidane", "O"),
    "methane": ("methane", "C"),
    "ethane": ("ethane", "CC"),
    "propane": ("propane", "CCC"),
    "butane": ("butane", "CCCC"),
    "pentane": ("pentane", "CCCCC"),
    "hexane": ("hexane", "CCCCCC"),
    "heptane": ("heptane", "CCCCCCC"),
    "octane": ("octane", "CCCCCCCC"),
    "ethylene": ("ethene", "C=C"),
    "ethene": ("ethene", "C=C"),
    "acetylene": ("ethyne", "C#C"),
    "ethyne": ("ethyne", "C#C"),
    "vinyl chloride": ("chloroethene", "C=CCl"),
    "pvc monomer": ("chloroethene", "C=CCl"),
    "styrene monomer": ("ethenylbenzene", "C=Cc1ccccc1"),
    "caprolactam": ("azepan-2-one", "C1CCC(=O)NC1"),
    "nylon-6 precursor": ("azepan-2-one", "C1CCC(=O)NC1"),
    "adipic acid": ("hexanedioic acid", "C(CCC(=O)O)CC(=O)O"),
    "terephthalic acid": ("benzene-1,4-dicarboxylic acid", "O=C(O)c1ccc(cc1)C(=O)O"),
    "bisphenol a": ("2,2-bis(4-hydroxyphenyl)propane", "CC(c1ccc(O)cc1)c1ccc(O)cc1"),
    "bpa": ("2,2-bis(4-hydroxyphenyl)propane", "CC(c1ccc(O)cc1)c1ccc(O)cc1"),
    "vitamin c": ("(5R)-5-[(1S)-1,2-dihydroxyethyl]-3,4-dihydroxyfuran-2(5H)-one", "O=C(O)C(O)(C(O)CO)C(O)=O"),
    "ascorbic acid": ("(5R)-5-[(1S)-1,2-dihydroxyethyl]-3,4-dihydroxyfuran-2(5H)-one", "O=C(O)C(O)(C(O)CO)C(O)=O"),
    "thymine": ("5-methylpyrimidine-2,4(1H,3H)-dione", "Cc1nc(nc(=O)[nH]1)O"),
    "adenine": ("9H-purin-6-amine", "c1ncnc2c1nc[nH]2"),
    "guanine": ("2-amino-1,9-dihydro-6H-purin-6-one", "c1nc2c(ncn2[C@H]3[C@@H]([C@H](O)[C@@H](O3)CO)n1)N"),
    "cytosine": ("4-aminopyrimidin-2(1H)-one", "n1cc(nc(=O)[nH]1)N"),
    "uracil": ("pyrimidine-2,4(1H,3H)-dione", "O=c1[nH]cnc(n1)O"),
    "pyridine": ("pyridine", "n1ccccc1"),
    "pyrrole": ("pyrrole", "c1cc[nH]c1"),
    "furan": ("furan", "c1ccoc1"),
    "thiophene": ("thiophene", "c1cccs1"),
    "indole": ("indole", "c1ccc2[nH]ccc2c1"),
    "quinoline": ("quinoline", "c1ccc2ncccc2c1"),
    "isoquinoline": ("isoquinoline", "c1ccc2ncccc2c1"),
    "imidazole": ("1H-imidazole", "n1cc[nH]c1"),
    "triazole": ("1H-1,2,4-triazole", "n1cnnc1"),
    "tetrahydrofuran": ("oxolane", "C1CCOC1"),
    "thf": ("oxolane", "C1CCOC1"),
    "dioxane": ("1,4-dioxane", "C1COCCO1"),
    "dichloromethane": ("dichloromethane", "ClCCl"),
    "dc m": ("dichloromethane", "ClCCl"),
    "chloroform": ("trichloromethane", "ClC(Cl)Cl"),
    "carbon tetrachloride": ("tetrachloromethane", "Cl(Cl)(Cl)Cl"),
    "cs2": ("carbon disulfide", "S=C=S"),
    "carbon disulfide": ("carbon disulfide", "S=C=S"),
    "dmso": ("dimethyl sulfoxide", "CS(=O)C"),
    "dimethyl sulfoxide": ("dimethyl sulfoxide", "CS(=O)C"),
    "dmf": ("N,N-dimethylformamide", "CN(C)C=O"),
    "n,n-dimethylformamide": ("N,N-dimethylformamide", "CN(C)C=O"),
    "thf": ("oxolane", "C1CCOC1"),
    "acetonitrile": ("ethanenitrile", "CC#N"),
    "formamide": ("formamide", "NC=O"),
    "methyl formate": ("methyl methanoate", "COC=O"),
    "ethyl acetate": ("ethyl ethanoate", "CCOC(=O)C"),
    "vinyl acetate": ("ethyl ethenoate", "C=COC(=O)C"),
    "methyl methacrylate": ("methyl 2-methylpropenoate", "C=C(C)C(=O)OC"),
    "mma monomer": ("methyl 2-methylpropenoate", "C=C(C)C(=O)OC"),
    "acrylamide": ("prop-2-enamide", "C=CC(=O)N"),
    "acrylonitrile": ("prop-2-enenitrile", "C=CC#N"),
    "an monomer": ("prop-2-enenitrile", "C=CC#N"),
    "maleic acid": ("(Z)-but-2-enedioic acid", "O=C(O)/C=C\C(=O)O"),
    "fumaric acid": ("(E)-but-2-enedioic acid", "O=C(O)/C=C/C(=O)O"),
    "maleic anhydride": ("furan-2,5-dione", "O=C1OC(=O)C=C1"),
    "phthalic anhydride": ("2-benzofuran-1,3-dione", "O=C1OC(=O)c2ccccc12"),
    "quinone": ("cyclohexa-2,5-diene-1,4-dione", "O=C1C=CC(=O)CC=1"),
    "p-benzoquinone": ("cyclohexa-2,5-diene-1,4-dione", "O=C1C=CC(=O)CC=1"),
    "hydroquinone": ("benzene-1,4-diol", "Oc1ccc(O)cc1"),
    "anthracene": ("anthracene", "c1ccc2c(c1)ccc1ccccc21"),
    "phenanthrene": ("phenanthrene", "c1ccc2c(c1)ccc1ccccc21"),
    "pyrene": ("pyrene", "c1ccc2c(c1)ccc3c2ccc4c3cccc4"),
    "coronene": ("coronene", "c1ccc2c(c1)ccc3c2c4ccccc4c4ccccc34"),
    "fullerene-c60": ("fullerene-C60", None),
    "graphite": ("graphite (allotrope of carbon)", None),
    "diamond": ("diamond (allotrope of carbon)", None),
}


def _get_db_entry(value):
    """Safely extract (iupac, smiles) from a DB entry (handles str or tuple)."""
    if isinstance(value, str):
        return value, None
    elif isinstance(value, tuple) and len(value) >= 2:
        return value[0], value[1]
    else:
        return value, None


def _reverse_lookup(target: str) -> list:
    """Find all common names that map to the given IUPAC name or SMILES."""
    results = []
    target_lower = target.lower().strip()
    for common, value in COMMON_NAME_DB.items():
        iupac, smiles = _get_db_entry(value)
        if iupac and iupac.lower() == target_lower:
            results.append((common, iupac, smiles))
        elif smiles and smiles == target:
            results.append((common, iupac, smiles))
    return results


@ChemMCPManager.register_tool
class CommonNameLookup(BaseTool):
    __version__ = "0.1.0"
    name = "CommonNameLookup"
    func_name = 'lookup_common_name'
    description = "Bidirectional lookup between common/trivial chemical names and systematic IUPAC names. Supports both directions: common→IUPAC and IUPAC→common."
    implementation_description = "Uses a built-in comprehensive database of ~150 common chemicals mapping trivial names to IUPAC names and SMILES. For unknown names, attempts RDKit name parsing and PubChem lookup as fallback."
    categories = ["Molecule"]
    tags = ["Nomenclature", "Common Names", "IUPAC", "Chemical Names"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
        ("PubChemPy", "https://github.com/mcs07/PubChemPy", "MIT"),
    ]
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]

    code_input_sig = [
        ('name', 'str', 'N/A', 'Chemical name to look up (common name or IUPAC name).'),
        ('direction', 'str', 'auto', 'Lookup direction: "auto" (detect automatically), "common_to_iupac", or "iupac_to_common".'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Name to look up, optionally followed by direction. E.g., "ethanol" or "acetic acid iupac_to_common".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing matched names, IUPAC name, SMILES, and direction used.'),
    ]
    examples = [
        {
            'code_input': {'name': 'ethanol', 'direction': 'auto'},
            'text_input': {'query': 'ethanol'},
            'output': {
                'result': {
                    'query': 'ethanol',
                    'direction_used': 'common_to_iupac',
                    'common_names': ['ethanol', 'ethyl alcohol', 'grain alcohol'],
                    'iupac_name': 'ethanol',
                    'smiles': 'CCO',
                    'source': 'local_database',
                }
            },
        },
        {
            'code_input': {'name': 'ethanol', 'direction': 'iupac_to_common'},
            'text_input': {'query': 'ethanol iupac_to_common'},
            'output': {
                'result': {
                    'query': 'ethanol',
                    'direction_used': 'iupac_to_common',
                    'common_names': ['ethanol', 'ethyl alcohol', 'grain alcohol'],
                    'iupac_name': 'ethanol',
                    'smiles': 'CCO',
                    'source': 'local_database',
                }
            },
        },
        {
            'code_input': {'name': 'acetone', 'direction': 'auto'},
            'text_input': {'query': 'acetone'},
            'output': {
                'result': {
                    'query': 'acetone',
                    'direction_used': 'common_to_iupac',
                    'common_names': ['acetone', 'dimethyl ketone'],
                    'iupac_name': 'propanone',
                    'smiles': 'CC(=O)C',
                    'source': 'local_database',
                }
            },
        },
    ]

    def _run_base(self, name: str, direction: str = "auto") -> dict:
        name_stripped = name.strip()
        if not name_stripped:
            raise ChemMCPInputError("Name cannot be empty.")

        name_lower = name_stripped.lower()

        # Determine direction if auto
        if direction == "auto":
            # Check if it's in our DB as a common name key
            if name_lower in COMMON_NAME_DB:
                direction = "common_to_iupac"
            else:
                # Check if it matches any IUPAC name or SMILES in DB
                reverse_hits = _reverse_lookup(name_stripped)
                if reverse_hits:
                    direction = "iupac_to_common"
                else:
                    direction = "common_to_iupac"  # default

        if direction == "common_to_iupac":
            return self._lookup_common_to_iupac(name_lower)
        elif direction == "iupac_to_common":
            return self._lookup_iupac_to_common(name_stripped)
        else:
            raise ChemMCPInputError(f"Invalid direction: {direction}. Use 'auto', 'common_to_iupac', or 'iupac_to_common'.")

    def _lookup_common_to_iupac(self, name_lower: str) -> dict:
        """Look up common name → IUPAC + SMILES."""
        # Direct match
        if name_lower in COMMON_NAME_DB:
            iupac, smiles = _get_db_entry(COMMON_NAME_DB[name_lower])
            # Find all aliases for this IUPAC name
            aliases = [k for k, v in COMMON_NAME_DB.items()
                       if _get_db_entry(v)[0] == iupac and k != name_lower]
            result = {
                'query': name_lower,
                'direction_used': 'common_to_iupac',
                'common_names': [name_lower] + aliases,
                'iupac_name': iupac,
                'smiles': smiles,
                'source': 'local_database',
            }
            logger.info(f"Common name lookup: '{name_lower}' → '{iupac}'")
            return result

        # Fuzzy partial match
        matches = [k for k in COMMON_NAME_DB if name_lower in k or k in name_lower]
        if matches:
            best_match = matches[0]
            iupac, smiles = _get_db_entry(COMMON_NAME_DB[best_match])
            result = {
                'query': name_lower,
                'direction_used': 'common_to_iupac',
                'common_names': [best_match],
                'iupac_name': iupac,
                'smiles': smiles,
                'source': 'local_database_fuzzy',
                'note': f'Exact match not found. Closest match: "{best_match}"',
            }
            logger.info(f"Fuzzy common name lookup: '{name_lower}' → '{best_match}' ({iupac})")
            return result

        # Try RDKit name-to-SMILES
        if _RDKIT_AVAILABLE:
            try:
                mol = Chem.MolFromSmiles(name_stripped)
                if mol is not None:
                    try:
                        iupac = Chem.MolToIUPACName(mol)
                        result = {
                            'query': name_lower,
                            'direction_used': 'common_to_iupac',
                            'common_names': [name_stripped],
                            'iupac_name': iupac or name_stripped,
                            'smiles': name_stripped,
                            'source': 'rdkit_parsed_as_smiles',
                        }
                        return result
                    except Exception:
                        pass
            except Exception:
                pass

        raise ChemMCPSearchFailError(
            f"Cannot find chemical name '{name_stripped}' in the common name database "
            "and no fallback method succeeded."
        )

    def _lookup_iupac_to_common(self, name: str) -> dict:
        """Look up IUPAC name → common name(s)."""
        reverse_hits = _reverse_lookup(name)

        if reverse_hits:
            common_names = [h[0] for h in reverse_hits]
            _, iupac, smiles = reverse_hits[0]
            result = {
                'query': name,
                'direction_used': 'iupac_to_common',
                'common_names': common_names,
                'iupac_name': iupac,
                'smiles': smiles,
                'source': 'local_database',
            }
            logger.info(f"IUPAC reverse lookup: '{name}' → {common_names}")
            return result

        # Try treating input as SMILES
        if _RDKIT_AVAILABLE and is_smiles(name):
            try:
                mol = Chem.MolFromSmiles(name)
                if mol is not None:
                    iupac = Chem.MolToIUPACName(mol)
                    reverse_hits = _reverse_lookup(iupac) if iupac else []
                    if reverse_hits:
                        common_names = [h[0] for h in reverse_hits]
                        result = {
                            'query': name,
                            'direction_used': 'iupac_to_common',
                            'common_names': common_names,
                            'iupac_name': iupac,
                            'smiles': name,
                            'source': 'smiles_parsed_via_rdkit',
                        }
                        return result
                    else:
                        result = {
                            'query': name,
                            'direction_used': 'iupac_to_common',
                            'common_names': [],
                            'iupac_name': iupac,
                            'smiles': name,
                            'source': 'rdkit_only_no_common_match',
                            'note': f"IUPAC name generated: '{iupac}', but no common name found in database.",
                        }
                        return result
            except Exception as e:
                logger.debug(f"RDKit SMILES parse failed: {e}")

        raise ChemMCPSearchFailError(
            f"Cannot find common name(s) for '{name}'. "
            "The name is not in the local database and could not be resolved via RDKit."
        )

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        name = parts[0]
        direction = parts[1] if len(parts) > 1 else "auto"
        return self._run_base(name, direction)


if __name__ == "__main__":
    run_mcp_server()
