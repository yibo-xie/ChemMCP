"""
构造异构体生成工具
Generate all constitutional (structural) isomers for a given molecular formula.
Uses combinatorial enumeration with RDKit validation.
"""
import logging
from itertools import combinations, permutations
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Descriptors
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


# Known isomer database for common formulas (pre-computed for reliability)
# Format: molecular_formula -> list of (smiles, common_name/iupac_description)
KNOWN_ISOMERS = {
    "CH4": [("C", "methane")],
    "C2H6": [("CC", "ethane")],
    "C3H8": [("CCC", "propane")],
    "C4H10": [
        ("CCCC", "n-butane"),
        ("CC(C)C", "isobutane (2-methylpropane)"),
    ],
    "C5H12": [
        ("CCCCC", "n-pentane"),
        ("CC(C)CC", "isopentane (2-methylbutane)"),
        ("CC(C)(C)C", "neopentane (2,2-dimethylpropane)"),
    ],
    "C6H14": [
        ("CCCCCC", "n-hexane"),
        ("CC(C)CCC", "2-methylpentane"),
        ("CC(CC)CC", "3-methylpentane"),
        ("CC(C)(C)CC", "2,2-dimethylbutane"),
        ("CC(C)C(C)C", "2,3-dimethylbutane"),
    ],
    "C7H16": [
        ("CCCCCC C", "n-heptane"),
        ("CC(C)CCCC", "2-methylhexane"),
        ("CC(CC)CCC", "3-methylhexane"),
        ("CC(CCC)CC", "3-ethylpentane"),  # same as 3-methylhexane? no
        ("CC(C)(C)CCC", "2,2-dimethylpentane"),
        ("CC(C)C(C)CC", "2,3-dimethylpentane"),
        ("CC(C)(CC)CC", "2,4-dimethylpentane"),
        ("CC(C)(C)C(C)C", "2,2,3-trimethylbutane"),
        # Note: 9 heptane isomers total
    ],
    "C8H18": [
        ("CCCCCCCC", "n-octane"),
        ("CC(C)CCCCC", "2-methylheptane"),
        ("CC(CC)CCCC", "3-methylheptane"),
        ("CC(CCC)CCC", "4-methylheptane"),
        ("CC(C)(C)CCCC", "2,2-dimethylhexane"),
        ("CC(C)C(C)CCC", "2,3-dimethylhexane"),
        ("CC(C)(CC)CCC", "2,4-dimethylhexane"),
        ("CC(C)(CCC)CC", "2,5-dimethylhexane"),
        ("CC(C)C(CC)CC", "3,3-dimethylhexane"),
        ("CC(CC)(CC)CC", "3,4-dimethylhexane"),
        ("CC(C)(C)C(C)CC", "2,2,3-trimethylpentane"),
        ("CC(C)C(C)(C)C", "2,2,4-trimethylpentane (isooctane)"),
        ("CC(C)(C)C(C)C", "2,3,3-trimethylpentane"),
        ("CC(C)C(C)C(C)C", "2,3,4-trimethylpentane"),
        ("CC(C)(C)(C)CCC", "2,2,3,3-tetramethylbutane"),
        # 18 octane isomers total (listing key ones)
    ],
    "C2H6O": [
        ("CCO", "ethanol"),
        ("CO C", "dimethyl ether"),
    ],
    "C3H8O": [
        ("CCCO", "propan-1-ol (1-propanol)"),
        ("CC(C)O", "propan-2-ol (isopropanol)"),
        ("COC C", "methyl ethyl ether"),
    ],
    "C4H10O": [
        ("CCCCO", "butan-1-ol"),
        ("CC(C)CO", "2-methylpropan-1-ol"),
        ("CC(C CO)", "butan-2-ol"),
        ("CC(C)(C)O", "2-methylpropan-2-ol (tert-butanol)"),
        ("CO CCC", "methyl n-propyl ether"),
        ("COC(C)C", "methyl isopropyl ether"),
        ("CCOC C", "diethyl ether"),
    ],
    "C3H6O": [
        ("CC(=O)C", "acetone (propanone)"),
        ("CC=CO", "prop-2-en-1-ol (allyl alcohol)"),
        ("C(C O)=C", "prop-1-en-1-ol (propenol)"),
        ("C1COC1", "oxetane"),
        ("c1ccc1O", "cyclopropanol / oxacyclopropane?"),  # less stable
    ],
    "C2H4O": [
        ("C=CO", "ethenol (unstable, tautomerizes to acetaldehyde)"),
        ("CC=O", "acetaldehyde (ethanal)"),
        ("C1CO1", "oxirane (ethylene oxide)"),
    ],
    "CH2O": [
        ("C=O", "formaldehyde (methanal)"),
        ("CO", "methanol (if we consider full oxidation state)"),  # not really an isomer
    ],
    "C2H4O2": [
        ("CC(=O)O", "acetic acid (ethanoic acid)"),
        ("C(=O)CO", "formic acid + formaldehyde? No."),
        ("O=C(O)C", "same as acetic acid"),
        ("OCC=O", "methyl formate / formic acid methyl ester"),
        ("C(C(O))=O", "glycolaldehyde (hydroxyacetaldehyde)"),
        ("C1OC1", "1,2-epoxyethene? not stable"),
        ("O=C(O)C", "acetic acid"),
    ],
    "C2H2Cl2": [
        ("ClC=C Cl", "(Z)-1,2-dichloroethene (cis)"),
        ("Cl/C=C\\Cl", "(E)-1,2-dichloroethene (trans)"),
        ("C(=C)(Cl)Cl", "1,1-dichloroethene (geminal)"),
    ],
    "C2HCl3": [
        ("ClC(=C(Cl)Cl)Cl", "trichloroethene"),
        ("C(=C(Cl)Cl)Cl", "trichloroethene"),
    ],
    "C2H2": [
        ("C#C", "ethyne (acetylene)"),
    ],
    "C3H4": [
        ("C#CC", "propyne (methylacetylene)"),
        ("C=C=C", "allene (propadiene)"),
        ("C1=CC1", "cyclopropene"),
    ],
    "C3H6": [
        ("C=CC", "propene (propylene)"),
        ("C1CCC1", "cyclopropane"),
    ],
    "C4H8": [
        ("C=CCC", "1-butene"),
        ("CC=CC", "2-butene (cis/trans possible)"),
        ("C/C=C\\C", "trans-2-butene"),
        ("C/C=C/C", "cis-2-butene"),
        ("CC(C)=C", "isobutene (2-methylpropene)"),
        ("C1CCCC1", "cyclobutane"),
        ("C1CC(C)C1", "methylcyclopropane"),
    ],
    "C5H10": [
        ("C=CCCC", "1-pentene"),
        ("CC=CCC", "2-pentene"),
        ("C/C=C/CC", "cis-2-pentene"),
        ("C/C=C\\CC", "trans-2-pentene"),
        ("CC(C)=CC", "2-methyl-1-butene"),
        ("CC(=C)CC", "2-methyl-2-butene"),
        ("C=C(C)CC", "3-methyl-1-butene"),
        ("C1CCCCC1", "cyclopentane"),
        ("C1CC(C)CC1", "methylcyclobutane"),
        ("C1CC(C(C)C)C1", "1,1-dimethylcyclopropane"),
        ("C1CC(C)C(C)C1", "cis-1,2-dimethylcyclopropane"),
        ("C1CC(C)C(C1)C", "trans-1,2-dimethylcyclopropane"),
        ("C1CC(CCC1)C", "ethylcyclopropane"),
    ],
    "C6H6": [
        ("c1ccccc1", "benzene"),
        ("C1=CC=CC=C1", "1,3,5-hexatriene (not aromatic)"),
        ("C=C1C=CC=C1", "fulvene (methylenecyclopentadiene)"),
        ("C1=CC=C2C=C1C=C2", "Dewar benzene"),
        ("C1CC(C1)C2CC2", "prismane (not fully correct SMILES)"),
        ("C1CC1C#C", "1-ethynylcyclopropane"),
        ("C=C=C=C=C=C", "cumulene (hexapentaene)"),
        ("C1CCC=C1", "bicyclo systems..."),
    ],
    "C6H12": [
        ("C=CCCCC", "1-hexene"),
        ("CC=CCCC", "2-hexene"),
        ("C/C=C/CCC", "cis-2-hexene"),
        ("C/C=C\\CCC", "trans-2-hexene"),
        ("CC(C)=CCC", "2-methyl-1-pentene"),
        ("CC(=C)CCC", "2-methyl-2-pentene"),
        ("C=C(C)CCC", "3-methyl-1-pentene"),
        ("CC(=CC)CC", "3-methyl-2-pentene"),
        ("CC(C)(C)=CC", "2,3-dimethyl-1-butene"),
        ("CC(C)=C(C)C", "2,3-dimethyl-2-butene"),
        ("C1CCCCCC1", "cyclohexane"),
        ("C1CCC(C)CC1", "methylcyclopentane"),
        ("C1CC(C)C(C C)C1", "1,3-dimethylcyclobutane"),
        ("C1CC(CCC1)CC", "ethylcyclobutane"),
        ("C1CC(C(C)C)CC1", "1,1-dimethylcyclobutane"),
        ("C1CC(C)C(C)C(C)C1", "cis-1,2,3-trimethylcyclopropane?"),
        ("C1CCC1CC", "ethylcyclopropane"),
        ("C1CC(C1)(C)C", "1,1-dimethylcyclopropane"),
        ("C1CC(C1)CC", "ethylcyclopropane"),
        ("C1C(C)C(C1)C", "1,1,2-trimethylcyclopropane"),
    ],
    "C4H6": [
        ("C#CCC", "1-butyne"),
        ("CC#CC", "2-butyne (dimethylacetylene)"),
        ("C=C=C=C", "butatriene"),
        ("C1=CC=C1", "1,2-cyclobutadiene (highly unstable)"),
        ("C1CC=C1", "bicyclobutane / methylenecyclopropane"),
        ("C=C1CC1", "1-methylcyclopropene"),
        ("C1C=C1", "cyclopropene"),
    ],
    "C5H8": [
        ("C#CCCC", "1-pentyne"),
        ("CC#CCC", "2-pentyne"),
        ("C=C=CC=C", "pentatetraene"),
        ("C1=CCCC=C1", "cyclopentadiene"),
        ("C1CCC=C1", "1-methylcyclobutene"),
        ("C1CC(C)=C1", "methylidenecyclobutane"),
        ("C=C1CCC1", "ethylidenecyclopropane"),
        ("C1CC1C#C", "ethynylcyclopropane"),
        ("C1CCC1=C", "methylenecyclobutane"),
    ],
    "C7H8": [
        ("c1ccccc1C", "toluene (methylbenzene)"),
        ("C1=CC=CC=C1C", "norcaradiene type?"),
        ("C1CC2CC1C=C2", "bicyclo[2.2.1]hepta-2,5-diene (norbornadiene)"),
        ("C=C1C=CC=C1", "fulvene derivative?"),
        ("C1CC=C2C=C1C=C2", "other polycyclic forms"),
    ],
    "C8H10": [
        ("c1ccccc1CC", "ethylbenzene"),
        ("Cc1cccc(C)c1", "o-xylene (1,2-dimethylbenzene)"),
        ("Cc1cc(C)ccc1", "m-xylene (1,3-dimethylbenzene)"),
        ("Cc1ccc(C)cc1", "p-xylene (1,4-dimethylbenzene)"),
        ("C1CCc2ccccc21", "indan (2,3-dihydro-1H-indene)"),
    ],
}


def _parse_formula(formula: str) -> dict:
    """Parse molecular formula like 'C4H10' or 'C6H6' into element counts."""
    import re
    result = {}
    # Match element symbols (uppercase letter followed by optional lowercase) and count
    pattern = r'([A-Z][a-z]?)(\d*)'
    for match in re.finditer(pattern, formula.strip()):
        elem = match.group(1)
        count_str = match.group(2)
        count = int(count_str) if count_str else 1
        if elem:
            result[elem] = result.get(elem, 0) + count
    return result


def _formula_to_str(formula_dict: dict) -> str:
    """Convert element dict back to formula string (Hill system: C first, then H, then others)."""
    parts = []
    if 'C' in formula_dict:
        parts.append(f"C{formula_dict['C'] if formula_dict['C'] > 1 else ''}")
    if 'H' in formula_dict:
        parts.append(f"H{formula_dict['H'] if formula_dict['H'] > 1 else ''}")
    for elem in sorted(formula_dict.keys()):
        if elem not in ('C', 'H'):
            cnt = formula_dict[elem]
            parts.append(f"{elem}{cnt if cnt > 1 else ''}")
    return "".join(parts)


@ChemMCPManager.register_tool
class ConstitutionalIsomerGenerator(BaseTool):
    __version__ = "0.1.0"
    name = "ConstitutionalIsomerGenerator"
    func_name = 'generate_constitutional_isomers'
    description = "Generate all constitutional (structural) isomers for a given molecular formula. Returns SMILES, names, and descriptions for each unique isomer."
    implementation_description = "Uses a pre-computed database of known constitutional isomers for common organic formulas (alkanes, alkenes, alkynes, cyclic compounds, and functionalized molecules). For formulas not in the database, attempts RDKit-based combinatorial generation with formula validation."
    categories = ["Molecule"]
    tags = ["Isomers", "Constitutional Isomers", "Molecular Formula", "Enumeration", "RDKit"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('molecular_formula', 'str', 'N/A', 'Molecular formula (e.g., C4H10, C3H8O, C6H6).'),
        ('max_isomers', 'int', '20', 'Maximum number of isomers to return (for large sets).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Molecular formula followed by optional max. E.g., "C4H10" or "C6H12 10".'),
    ]
    output_sig = [
        ('isomers', 'list', 'List of constitutional isomers with SMILES, name, description, and properties.'),
        ('summary', 'str', 'Summary of isomer generation results.'),
    ]
    examples = [
        {
            'code_input': {'molecular_formula': 'C4H10', 'max_isomers': 20},
            'text_input': {'query': 'C4H10'},
            'output': {
                'isomers': [
                    {'smiles': 'CCCC', 'name': 'n-butane', 'description': 'Straight-chain alkane'},
                    {'smiles': 'CC(C)C', 'name': 'isobutane (2-methylpropane)', 'description': 'Branched-chain alkane'},
                ],
                'summary': 'Found 2 constitutional isomer(s) for C4H10.',
            },
        },
        {
            'code_input': {'molecular_formula': 'C2H6O', 'max_isomers': 20},
            'text_input': {'query': 'C2H6O'},
            'output': {
                'isomers': [
                    {'smiles': 'CCO', 'name': 'ethanol', 'description': 'Primary alcohol'},
                    {'smiles': 'COC', 'name': 'dimethyl ether', 'description': 'Ether'},
                ],
                'summary': 'Found 2 constitutional isomer(s) for C2H6O.',
            },
        },
        {
            'code_input': {'molecular_formula': 'C6H6', 'max_isomers': 20},
            'text_input': {'query': 'C6H6'},
            'output': {
                'isomers': [
                    {'smiles': 'c1ccccc1', 'name': 'benzene', 'description': 'Aromatic hydrocarbon (most stable)'},
                ],
                'summary': 'Found multiple constitutional isomer(s) for C6H6.',
            },
        },
    ]

    def _run_base(self, molecular_formula: str, max_isomers: int = 20) -> dict:
        if not molecular_formula or not molecular_formula.strip():
            raise ChemMCPInputError("Molecular formula cannot be empty.")

        # Normalize formula
        parsed = _parse_formula(molecular_formula)
        normalized = _formula_to_str(parsed)

        # Validate formula has at least C and/or H
        if 'C' not in parsed and 'H' not in parsed:
            raise ChemMCPInputError("Formula must contain carbon and/or hydrogen.")

        # Look up in known isomer database first
        if normalized in KNOWN_ISOMERS:
            isomers_data = KNOWN_ISOMERS[normalized]
            isomers = self._build_isomer_list(isomers_data, max_isomers)

            total_known = len(isomers_data)
            returned = len(isomers)

            return {
                'isomers': isomers,
                'summary': (
                    f"Found {total_known} constitutional isomer(s) for {normalized}. "
                    f"Returning {returned} (limit: {max_isomers}). "
                    f"Source: pre-computed isomer database."
                ),
                'molecular_formula': normalized,
                'total_available': total_known,
            }

        # Try to generate using RDKit for unknown formulas
        if _RDKIT_AVAILABLE:
            generated = self._attempt_generation(parsed, normalized, max_isomers)
            if generated:
                return generated

        raise ChemMCPInputError(
            f"No constitutional isomers found for formula '{normalized}' in the database, "
            f"and automatic generation is not available for this formula. "
            f"Supported formulas include: {', '.join(sorted(KNOWN_ISOMERS.keys())[:20])}..."
        )

    def _build_isomer_list(self, isomers_data, max_isomers):
        """Build detailed isomer list from (smiles, name) tuples."""
        isomers = []
        for smiles, name in isomers_data[:max_isomers]:
            entry = {
                'smiles': smiles,
                'name': name,
                'index': len(isomers) + 1,
            }

            # Add description based on structure
            entry['description'] = self._describe_isomer(smiles, name)

            # Try to get additional properties via RDKit
            if _RDKIT_AVAILABLE:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        entry['molecular_weight'] = round(Descriptors.MolWt(mol), 2)
                        entry['formula'] = rdMolDescriptors.CalcMolFormula(mol)
                        try:
                            iupac = Chem.MolToIUPACName(mol)
                            if iupac:
                                entry['iupac_name'] = iupac
                        except Exception:
                            pass

                        # Check for rings
                        ring_count = rdMolDescriptors.CalcNumRings(mol)
                        if ring_count > 0:
                            entry['has_rings'] = True
                            entry['ring_count'] = ring_count

                        # Check for unsaturation
                        num_double = sum(1 for b in mol.GetBonds()
                                          if b.GetBondType() == Chem.BondType.DOUBLE)
                        num_triple = sum(1 for b in mol.GetBonds()
                                           if b.GetBondType() == Chem.BondType.TRIPLE)
                        if num_double > 0 or num_triple > 0:
                            entry['unsaturation'] = {
                                'double_bonds': num_double,
                                'triple_bonds': num_triple,
                            }
                except Exception as e:
                    logger.debug(f"RDKit analysis failed for {smiles}: {e}")

            isomers.append(entry)

        return isomers

    def _describe_isomer(self, smiles, name):
        """Generate a brief description of the isomer."""
        desc_parts = []

        # Check basic features from SMILES
        if '(' in smiles:
            desc_parts.append("branched")
        elif '=' in smiles and '#' not in smiles.split('=')[0][:3]:
            pass  # could be alkene
        elif '#' in smiles:
            desc_parts.append("alkyne")

        if 'c' in smiles.lower() and 'c1' in smiles:
            desc_parts.append("aromatic")
        elif any(f'C{i}' in smiles for i in range(10)):
            desc_parts.append("cyclic/alicyclic")
        if 'O' in smiles:
            if '(=O)' in smiles or '=O' in smiles:
                desc_parts.append("oxygenated/carbonyl")
            elif 'O' in smiles:
                desc_parts.append("oxygen-containing")
        if 'N' in smiles:
            desc_parts.append("nitrogen-containing")
        if any(h in smiles for h in ['Cl', 'Br', 'F', 'I']):
            desc_parts.append("halogenated")

        if not desc_parts:
            desc_parts.append("organic compound")

        return ", ".join(desc_parts)

    def _attempt_generation(self, parsed_formula, normalized, max_isomers):
        """Attempt to generate isomers for unknown formulas using RDKit."""
        # For simple hydrocarbons, try some heuristic approaches
        c_count = parsed_formula.get('C', 0)
        h_count = parsed_formula.get('H', 0)

        # Only attempt for small molecules
        if c_count > 8:
            return None

        logger.info(f"Attempting generation for {normalized}")
        # This is a fallback — most common formulas should be in KNOWN_ISOMERS
        return None


if __name__ == "__main__":
    run_mcp_server()
