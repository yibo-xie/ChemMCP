"""
E/Z构型判断工具
Determine E/Z configuration of double bonds using CIP priority rules.
Uses RDKit's double bond stereochemistry analysis.
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
class EzConfigurator(BaseTool):
    __version__ = "0.1.0"
    name = "EzConfigurator"
    func_name = 'determine_ez_configuration'
    description = "Determine E/Z (Entgegen/Zusammen) configuration of double bonds in a molecule using Cahn-Ingold-Prelog (CIP) priority rules."
    implementation_description = "Uses RDKit's FindPotentialStereo() and bond stereochemistry analysis to identify double bonds capable of E/Z isomerism, then determines E or Z assignment for each based on CIP priorities of substituents on each side of the double bond."
    categories = ["Molecule"]
    tags = ["Stereochemistry", "E/Z Configuration", "Double Bond", "CIP Rules", "RDKit"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('enumerate_options', 'bool', 'False', 'Whether to enumerate all possible E/Z combinations.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional enumerate flag.'),
    ]
    output_sig = [
        ('ez_configuration', 'dict', 'Detailed E/Z configuration for each double bond with atom indices, assignments, and CIP justifications.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'C/C=C/C', 'enumerate_options': False},
            'text_input': {'query': 'C/C=C/C false'},
            'output': {
                'ez_configuration': {
                    'double_bonds': [
                        {
                            'bond_index_info': 'between atoms 1 and 2',
                            'atoms': ['C', 'C'],
                            'configuration': 'Z',
                            'left_substituents': [{'atom': 'C (methyl)', 'priority': 1}, {'atom': 'H', 'priority': 2}],
                            'right_substituents': [{'atom': 'C (methyl)', 'priority': 1}, {'atom': 'H', 'priority': 2}],
                            'justification': 'Higher priority groups on same side → Z (zusammen).',
                        }
                    ],
                    'total_double_bonds': 1,
                    'summary': '1 double bond found: C(1)=C(2) has Z configuration.',
                }
            },
        },
        {
            'code_input': {'smiles': 'C/C=C\\C', 'enumerate_options': False},
            'text_input': {'query': 'C/C=C\\C false'},
            'output': {
                'ez_configuration': {
                    'double_bonds': [
                        {
                            'configuration': 'E',
                            'justification': 'Higher priority groups on opposite sides → E (entgegen).',
                        }
                    ],
                    'total_double_bonds': 1,
                    'summary': '1 double bond found with E configuration.',
                }
            },
        },
        {
            'code_input': {'smiles': 'CC=CC', 'enumerate_options': False},
            'text_input': {'query': 'CC=CC false'},
            'output': {
                'ez_configuration': {
                    'double_bonds': [
                        {
                            'configuration': 'unspecified',
                            'can_show_ez_isomerism': True,
                            'justification': 'Double bond exists but E/Z not specified in input SMILES. Both E and Z isomers are possible.',
                        }
                    ],
                    'summary': '1 double bond capable of E/Z isomerism found (unspecified).',
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

        # Assign stereochemistry
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Find all double bonds
        double_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                db_info = self._analyze_double_bond(mol, bond)
                if db_info:
                    double_bonds.append(db_info)

        if not double_bonds:
            return {
                'ez_configuration': {
                    'double_bonds': [],
                    'total_double_bonds': 0,
                    'summary': 'No double bonds found in this molecule.',
                }
            }

        total = len(double_bonds)
        configs_str = ", ".join(
            f"{db.get('bond_atoms_str', f'db_{i}')}: {db['configuration']}"
            for i, db in enumerate(double_bonds)
        )
        summary = f"{total} double bond(s) analyzed: {configs_str}."

        result = {
            'ez_configuration': {
                'double_bonds': double_bonds,
                'total_double_bonds': total,
                'summary': summary,
            }
        }

        # Optionally enumerate stereoisomers
        if enumerate_options:
            try:
                isomers = list(EnumerateStereoisomers(mol))
                result['ez_configuration']['possible_stereoisomers'] = [
                    {
                        'smiles': Chem.MolToSmiles(iso),
                    }
                    for iso in isomers
                ]
                result['ez_configuration']['n_possible_stereoisomers'] = len(isomers)
            except Exception as e:
                logger.debug(f"Stereoisomer enumeration failed: {e}")

        logger.info(f"E/Z config: {smiles} → {total} double bonds")
        return result

    def _analyze_double_bond(self, mol, bond):
        """Analyze a single double bond for E/Z configuration."""
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        begin_atom = mol.GetAtomWithIdx(begin_idx)
        end_atom = mol.GetAtomWithIdx(end_idx)

        # Get stereo info from the bond
        stereo = bond.GetStereo()
        bond_atoms_str = f"{begin_atom.GetSymbol()}({begin_idx})={end_atom.GetSymbol()}({end_idx})"

        # Get substituents on each carbon of the double bond
        left_subs = self._get_substituents(mol, begin_atom, exclude_idx=end_idx)
        right_subs = self._get_substituents(mol, end_atom, exclude_idx=begin_idx)

        # Check if E/Z isomerism is possible (need 2 different subs on each side)
        left_symbols = set(s['symbol'] for s in left_subs)
        right_symbols = set(s['symbol'] for s in right_subs)
        can_ez = len(left_subs) >= 2 and len(right_subs) >= 2

        # Determine configuration
        config = None
        justification = ""

        if stereo == Chem.BondStereo.STEREOZ or stereo == Chem.BondStereo.STEREOCIS:
            config = "Z"
            justification = (
                f"Higher priority groups on both sides of the {bond_atoms_str} double bond "
                f"are on the same side → Z (zusammen, 'together')."
            )
        elif stereo == Chem.BondStereo.STEREOE or stereo == Chem.BondStereo.STEREOTRANS:
            config = "E"
            justification = (
                f"Higher priority groups on both sides of the {bond_atoms_str} double bond "
                f"are on opposite sides → E (entgegen, 'opposite')."
            )
        else:
            config = "unspecified"
            if can_ez:
                justification = (
                    f"The {bond_atoms_str} double bond can show E/Z isomerism, "
                    f"but stereochemistry was not specified in the input. "
                    f"Both E and Z isomers are possible."
                )
            else:
                justification = (
                    f"The {bond_atoms_str} double bond does not have distinct E/Z isomers "
                    f"(one or both sides have identical substituents)."
                )

        # Sort substituents by atomic number for display
        left_sorted = sorted(left_subs, key=lambda x: -x['atomic_number'])
        right_sorted = sorted(right_subs, key=lambda x: -x['atomic_number'])

        return {
            'bond_index_info': f'between atoms {begin_idx} and {end_idx}',
            'bond_atoms_str': bond_atoms_str,
            'atoms': [begin_atom.GetSymbol(), end_atom.GetSymbol()],
            'begin_atom_idx': begin_idx,
            'end_atom_idx': end_idx,
            'configuration': config,
            'stereo_tag': str(stereo),
            'can_show_ez_isomerism': can_ez,
            'left_substituents': left_sorted,
            'right_substituents': right_sorted,
            'justification': justification,
        }

    def _get_substituents(self, mol, atom, exclude_idx):
        """Get substituent information for an atom, excluding the double-bond partner."""
        subs = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() == exclude_idx:
                continue
            sub_symbol = neighbor.GetSymbol()
            # Check if it's part of a larger group
            n_h = neighbor.GetTotalNumHs()
            if n_h > 0:
                display = f"{sub_symbol}H{n_h}" if n_h > 1 else f"{sub_symbol}H"
            else:
                # Check for larger group
                group_desc = self._describe_group(mol, neighbor, atom.GetIdx())
                display = group_desc if group_desc else sub_symbol
            subs.append({
                'index': neighbor.GetIdx(),
                'symbol': display,
                'base_symbol': sub_symbol,
                'atomic_number': neighbor.GetAtomicNum(),
            })

        # Add implicit H
        n_implicit = atom.GetTotalNumHs() - sum(1 for s in subs if 'H' in s['base_symbol'])
        if n_implicit > 0:
            subs.append({
                'index': None,
                'symbol': 'H' if n_implicit == 1 else f'H×{n_implicit}',
                'base_symbol': 'H',
                'atomic_number': 1,
            })

        return subs

    def _describe_group(self, mol, start_atom, exclude_idx):
        """Briefly describe a substituent group."""
        symbol = start_atom.GetSymbol()
        if symbol == 'C':
            neighbors = [n for n in start_atom.GetNeighbors()
                         if n.GetIdx() != exclude_idx]
            n_c = sum(1 for n in neighbors if n.GetSymbol() == 'C')
            if n_c >= 2:
                return f"C (group)"
            elif any(n.GetSymbol() == 'O' for n in neighbors):
                return "C (oxygenated)"
            return "C (alkyl)"
        return None


if __name__ == "__main__":
    run_mcp_server()
