"""
官能团识别工具（增强版）
Identify ALL functional groups in a molecule comprehensively.
Enhanced version with positions, counts, and SMARTS pattern details.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

# Comprehensive functional group definitions with SMARTS patterns
# Format: (group_name, smarts_pattern, category, description)
FUNCTIONAL_GROUPS_DB = [
    # === Oxygen-containing groups ===
    ("Alcohol (hydroxyl)", "[OH]", "oxygen", "Hydroxyl group (-OH) attached to sp³ carbon"),
    ("Phenol (phenolic -OH)", "[Ox2H][Cx1]", "oxygen", "Hydroxyl directly bonded to aromatic ring"),
    ("Aldehyde", "[CX3H1](=O)[#6]", "oxygen", "Carbonyl with at least one H (terminal C=O)"),
    ("Ketone", "[#6][CX3](=O)[#6]", "oxygen", "Carbonyl between two carbons (internal C=O)"),
    ("Carboxylic acid", "[CX3](=O)[OX2H1]", "oxygen", "-COOH group"),
    ("Ester", "[#6][CX3](=O)[OX2H0][#6]", "oxygen", "-COO- linkage between carbon atoms"),
    ("Acid chloride", "[CX3](=O)[Cl]", "oxygen", "-COCl group (acyl halide)"),
    ("Acid anhydride", "[#6][CX3](=O)[OX2][CX3](=O)[#6]", "oxygen", "-CO-O-CO- linkage"),
    ("Ether (acyclic)", "[OD2]([#6])[#6]", "oxygen", "C-O-C where O is bonded to two carbons"),
    ("Ether (cyclic)", "[OD2]1[CX4][CX4]1", "oxygen", "Oxygen in a ring (e.g., THF, dioxane)"),
    ("Peroxide", "[OX2][OX2]", "oxygen", "O-O single bond (ROOR)"),
    ("Hydroperoxide", "[OX2][OX2H]", "oxygen", "-OOH group"),
    ("Acetal/Ketal", "[OX2]([#6])[#6][OX2]([#6])", "oxygen", "Two oxygens on same carbon (C(OR)2)"),
    ("Hemiacetal/Hemiketal", "[OX2H]([#6])[OX2]([#6])", "oxygen", "One OH and one OR on same carbon"),

    # === Nitrogen-containing groups ===
    ("Primary amine", "[NX3;H2;!$(NC=O);!$(N[CD3])]", "nitrogen", "-NH2 group (not amide)"),
    ("Secondary amine", "[NX3;H1;!$(NC=O)][#6]", "nitrogen", ">NH group (not amide)"),
    ("Tertiary amine", "[NX3;H0;!$(NC=O)]([#6])([#6])[#6]", "nitrogen", "N bonded to 3 carbons"),
    ("Amide", "[NX3][CX3](=[OX1])", "nitrogen", "-CONH- or -CONR2 group"),
    ("Imide", "[NX3][CX3](=[OX1])[CX3](=[OX1])", "nitrogen", "Nitrogen between two carbonyls"),
    ("Nitrile (cyano)", "[C-]#[N+]", "nitrogen", "-C≡N triple bond"),
    ("Isocyanide", "[CX2+]#[N-]", "nitrogen", "-N⁺≡C⁻ isocyanide"),
    ("Nitro", "[N+](=O)[O-]", "nitrogen", "-NO2 group"),
    ("Nitroso", "[N](=O)", "nitrogen", "-N=O group"),
    ("Azo", "[N]=[N]", "nitrogen", "-N=N- double bond"),
    ("Diazo", "[N]=[N+]=[N-]", "nitrogen", "-N=N⁺=N⁻ diazo group"),
    ("Hydrazine", "[NX3][NX3]", "nitrogen", "-NH-NH- group"),
    ("Imine (Schiff base)", "[NX2]=[#6]", "nitrogen", "C=N double bond (not aromatic)"),
    ("Oxime", "[NX3]=[OX1]", "nitrogen", "=N-OH group"),
    ("Hydrazone", "[NX3]=[NX3]", "nitrogen", "=N-N< group"),
    ("Enamine", "[NX3][#6]=[#6]", "nitrogen", "N-C=C group"),
    ("Pyridine-like N", "n", "nitrogen_heterocycle", "Aromatic nitrogen in ring (pyridine-type)"),
    ("Pyrrole-like N", "[nH]", "nitrogen_heterocycle", "Aromatic NH in ring (pyrrole-type)"),
    ("Amidine", "[NX3]=[CX3][NX3]", "nitrogen", "-C(=NH)-NH2 group"),
    ("Urea/Carbamate", "[NX3][CX3](=[OX1])[NX3]", "nitrogen", "-NH-CO-NH- group"),

    # === Sulfur-containing groups ===
    ("Thiol (mercaptan)", "[SH]", "sulfur", "-SH group"),
    ("Thioether (sulfide)", "[SD2]([#6])[#6]", "sulfur", "C-S-C linkage"),
    ("Disulfide", "[SD2][SD2]", "sulfur", "-S-S- disulfide bridge"),
    ("Sulfoxide", "[S;D3](=O)([#6])[#6]", "sulfur", "-S(=O)- sulfoxide"),
    ("Sulfone", "[S;D4](=O)(=O)([#6])[#6]", "sulfur", "-SO2- sulfone"),
    ("Sulfonic acid", "[S](=O)(=O)[O]", "sulfur", "-SO3H group"),
    ("Sulfonate ester", "[S](=O)(=O)[O][#6]", "sulfur", "-SO3R sulfonate ester"),
    ("Sulfonyl chloride", "[S](=O)(=O)Cl", "sulfur", "-SO2Cl group"),
    ("Sulfonamide", "[S](=O)(=O)[NX3]", "sulfur", "-SO2NH2 group"),
    ("Thiocarbonyl", "[#6]=[S]", "sulfur", "C=S double bond"),
    ("Thiocarboxylic acid", "[CX3](=[SX1])[OX2H1]", "sulfur", "-C(=S)OH group"),
    ("Isothiocyanate", "[N]=C=S", "sulfur", "-N=C=S group"),

    # === Halogen groups ===
    ("Fluorine substituent", "[F]", "halogen", "-F fluorine atom"),
    ("Chlorine substituent", "[Cl]", "halogen", "-Cl chlorine atom"),
    ("Bromine substituent", "[Br]", "halogen", "-Br bromine atom"),
    ("Iodine substituent", "[I]", "halogen", "-I iodine atom"),
    ("Polyfluoro", "[#6](F)(F)F", "halogen", "-CF3 trifluoromethyl group"),
    ("Perhaloalkyl", "[#6]([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]", "halogen", "Carbon with multiple halogens"),

    # === Carbon-carbon multiple bonds ===
    ("Alkene (C=C)", "[#6]=[#6;!$([#6]=[#6][#6]=[#6])]", "unsaturation", "Carbon-carbon double bond (non-aromatic)"),
    ("Alkyne (C≡C)", "[#6]#[#6]", "unsaturation", "Carbon-carbon triple bond"),
    ("Cumulene", "[#6]=[#6]=[#6]", "unsaturation", "Consecutive double bonds (C=C=C)"),
    ("Allene", "[#6]=[#6]=[#6]", "unsaturation", "Allene (C=C=C) substructure"),

    # === Aromatic systems ===
    ("Benzene ring", "c1ccccc1", "aromatic", "6-membered aromatic carbon ring"),
    ("Phenyl group", "c1ccccc1", "aromatic", "Benzene ring as substituent"),
    ("Fused benzene (naphthalene)", "c1ccc2ccccc2c1", "aromatic_polycyclic", "Two fused benzene rings"),
    ("Heteroaromatic (5-membered)", "[a][a][a][a][a]", "aromatic_heterocycle", "5-membered aromatic heterocycle"),
    ("Heteroaromatic (6-membered)", "[a][a][a][a][a][a]", "aromatic_heterocycle", "6-membered aromatic heterocycle"),

    # === Carbonyl derivatives ===
    ("Ketene", "[#6]=[C]=[OX1]", "carbonyl", "C=C=O ketene group"),
    ("β-lactam (4-ring amide)", "[C]1[C][C]([N]1)[CX3](=[OX1])", "heterocycle", "4-membered cyclic amide"),
    ("γ-lactam (5-ring amide)", "[C]1CC[C]([N]1)[CX3](=[OX1])", "heterocycle", "5-membered cyclic amide"),
    ("δ-lactam (6-ring amide)", "[C]1CCC[C]([N]1)[CX3](=[OX1])", "heterocycle", "6-membered cyclic amide"),
    ("β-lactone (4-ring ester)", "[C]1[C][C]([O]1)[CX3](=[OX1])", "heterocycle", "4-membered cyclic ester"),
    ("γ-lactone (5-ring ester)", "[C]1CC[C]([O]1)[CX3](=[OX1])", "heterocycle", "5-membered cyclic ester"),
    ("δ-lactone (6-ring ester)", "[C]1CCC[C]([O]1)[CX3](=[OX1])", "heterocycle", "6-membered cyclic ester"),

    # === Phosphorus groups ===
    ("Phosphine", "[P]([#6])([#6])[#6]", "phosphorus", "Tertiary phosphine PR3"),
    ("Phosphate ester", "[P](=O)([O])([O])[O]", "phosphorus", "Phosphate -OP(O)(O)O- group"),
    ("Phosphonate", "[P]([#6])(=O)([O])[O]", "phosphorus", "Phosphonate R-PO(OR)2 group"),

    # === Other notable groups ===
    ("Epoxide (oxirane)", "[C]1[C][O]1", "heterocycle", "3-membered oxygen-containing ring"),
    ("Oxirane", "[C]1[C][O]1", "heterocycle", "3-membered epoxide ring"),
    ("Aziridine", "[C]1[C][N]1", "heterocycle", "3-membered nitrogen-containing ring"),
    ("Cyclopropane", "[C]1[C][C]1", "ring_strain", "3-membered carbon ring (high strain)"),
    ("Cyclobutane", "[C]1[C][C][C]1", "ring_strain", "4-membered carbon ring"),
    ("Quaternary carbon", "[#6]([#6])([#6])([#6])[#6]", "structural", "Carbon bonded to 4 other carbons"),
    ("Gem-dihalide", "[#6]([#9,#17,#35,#53])([F,Cl,Br,I])", "structural", "Carbon with two halogens"),
    ("α,β-unsaturated carbonyl", "[#6]=[#6][CX3](=[OX1])", "conjugated", "Conjugated enone system"),
    ("1,3-Diene", "[#6]=[#6][#6]=[#6]", "conjugated", "Conjugated diene system"),
    ("Benzylic position", "[#6][c]1ccccc1", "structural", "CH2/C directly attached to aromatic ring"),
    ("Allylic position", "[#6][#6]=[#6]", "structural", "Carbon adjacent to C=C"),
    ("Propargylic position", "[#6][#6]#[#6]", "structural", "Carbon adjacent to C≡C"),
]


@ChemMCPManager.register_tool
class FunctionalGroupIdentifier(BaseTool):
    __version__ = "0.1.0"
    name = "FunctionalGroupIdentifier"
    func_name = 'identify_functional_groups'
    description = "Comprehensively identify all functional groups in a molecule, including group names, counts, atom positions, SMARTS patterns used, and categorization."
    implementation_description = "Uses a comprehensive database of ~80+ functional group SMARTS patterns covering oxygen-, nitrogen-, sulfur-, halogen-containing groups, unsaturated systems, aromatic systems, carbonyl derivatives, phosphorus groups, and structural motifs. Each pattern is matched against the molecule using RDKit's substructure search."
    categories = ["Molecule"]
    tags = ["Functional Groups", "SMILES", "RDKit", "SMARTS", "Molecular Analysis"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
    ]
    text_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
    ]
    output_sig = [
        ('functional_groups', 'dict', 'Detailed dictionary of all identified functional groups with names, counts, positions, and categories.'),
        ('summary', 'str', 'Human-readable summary of functional groups found.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'CCO'},
            'text_input': {'smiles': 'CCO'},
            'output': {
                'functional_groups': {'groups': [{'name': 'Alcohol (hydroxyl)', 'count': 1, 'category': 'oxygen'}]},
                'summary': 'This molecule contains 1 functional group(s): Alcohol (hydroxyl).',
            },
        },
        {
            'code_input': {'smiles': 'CC(=O)O'},
            'text_input': {'smiles': 'CC(=O)O'},
            'output': {
                'functional_groups': {'groups': [{'name': 'Carboxylic acid', 'count': 1, 'category': 'oxygen'}]},
                'summary': 'This molecule contains 1 functional group(s): Carboxylic acid.',
            },
        },
        {
            'code_input': {'smiles': 'c1ccccc1'},
            'text_input': {'smiles': 'c1ccccc1'},
            'output': {
                'functional_groups': {'groups': [{'name': 'Benzene ring', 'count': 1, 'category': 'aromatic'}]},
                'summary': 'This molecule contains 1 functional group(s): Benzene ring.',
            },
        },
    ]

    def _run_base(self, smiles: str) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Scan all functional group patterns
        detected_groups = []
        for fg_name, smarts, category, description in FUNCTIONAL_GROUPS_DB:
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is None:
                    continue

                matches = mol.GetSubstructMatches(pattern)
                if matches:
                    # Get matched atom indices for each occurrence
                    match_details = []
                    for match in matches:
                        match_details.append({
                            'atom_indices': list(match),
                            'atom_symbols': [mol.GetAtomWithIdx(i).GetSymbol() for i in match],
                        })

                    detected_groups.append({
                        'name': fg_name,
                        'smarts_pattern': smarts,
                        'category': category,
                        'description': description,
                        'count': len(matches),
                        'matches': match_details,
                    })
            except Exception as e:
                logger.debug(f"Pattern matching failed for {fg_name}: {e}")
                continue

        # Group by category for summary
        by_category = {}
        for g in detected_groups:
            cat = g['category']
            if cat not in by_category:
                by_category[cat] = []
            by_category[cat].append(g)

        # Build human-readable summary
        if detected_groups:
            group_names = [f"{g['name']} ({g['count']})" for g in detected_groups]
            summary = f"This molecule contains {len(detected_groups)} functional group type(s): {', '.join(group_names)}."
        else:
            summary = "No known functional groups were identified in this molecule."

        result = {
            'functional_groups': {
                'total_types': len(detected_groups),
                'total_occurrences': sum(g['count'] for g in detected_groups),
                'groups': detected_groups,
                'by_category': {k: [{'name': g['name'], 'count': g['count'], 'description': g['description']}
                                   for g in v] for k, v in by_category.items()},
            },
            'summary': summary,
            'smiles_input': smiles,
        }

        logger.info(f"FG identification: {smiles} → {len(detected_groups)} types found")
        return result


if __name__ == "__main__":
    run_mcp_server()
