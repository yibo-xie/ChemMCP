"""
SMILES转结构式描述工具
Convert SMILES string to human-readable structural formula description.
Analyzes molecule using RDKit and generates detailed text description.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, Draw, rdmolfiles, rdmolops
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@ChemMCPManager.register_tool
class SmilesToStructure(BaseTool):
    __version__ = "0.1.0"
    name = "SmilesToStructure"
    func_name = 'smiles_to_structure'
    description = "Convert a SMILES string into a human-readable structural formula description including atom connectivity, bond types, rings, and functional groups."
    implementation_description = "Parses the SMILES string using RDKit, then analyzes the molecular graph to produce a structured text description: atom list with hybridization, bond connectivity, ring systems, functional groups, and overall molecular characteristics."
    categories = ["Molecule"]
    tags = ["SMILES", "Structure", "Molecular Description", "RDKit", "Visualization"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule.'),
        ('detail_level', 'str', 'medium', 'Level of detail: "brief", "medium", or "detailed".'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'SMILES string followed by optional detail level. E.g., "CCO medium".'),
    ]
    output_sig = [
        ('structure_description', 'str', 'Human-readable structural formula description of the molecule.'),
        ('molecular_formula', 'str', 'Molecular formula (e.g., C2H6O).'),
        ('molecular_weight', 'float', 'Molecular weight in g/mol.'),
    ]
    examples = [
        {
            'code_input': {'smiles': 'CCO', 'detail_level': 'medium'},
            'text_input': {'query': 'CCO medium'},
            'output': {
                'structure_description': 'Ethanol: A 2-carbon chain with a hydroxyl group on the terminal carbon. Structure: CH3-CH2-OH.',
                'molecular_formula': 'C2H6O',
                'molecular_weight': 46.07,
            },
        },
        {
            'code_input': {'smiles': 'c1ccccc1', 'detail_level': 'brief'},
            'text_input': {'query': 'c1ccccc1 brief'},
            'output': {
                'structure_description': 'Benzene: An aromatic 6-membered carbon ring (C6H6) with alternating double bonds.',
                'molecular_formula': 'C6H6',
                'molecular_weight': 78.11,
            },
        },
        {
            'code_input': {'smiles': 'CC(=O)O', 'detail_level': 'detailed'},
            'text_input': {'query': 'CC(=O)O detailed'},
            'output': {
                'structure_description': 'Acetic acid: A 2-carbon carboxylic acid. The methyl group (CH3-) is bonded to a carboxyl group (-COOH). Structure: CH3-C(=O)-OH.',
                'molecular_formula': 'C2H4O2',
                'molecular_weight': 60.05,
            },
        },
    ]

    def _run_base(self, smiles: str, detail_level: str = "medium") -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # Basic properties
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)

        # Generate description based on detail level
        if detail_level == "brief":
            desc = self._describe_brief(mol)
        elif detail_level == "detailed":
            desc = self._describe_detailed(mol)
        else:
            desc = self._describe_medium(mol)

        logger.info(f"Structure description for {smiles} ({detail_level})")

        return {
            'structure_description': desc,
            'molecular_formula': formula,
            'molecular_weight': round(mw, 2),
        }

    def _get_atom_info(self, mol):
        """Get basic atom information."""
        atoms = []
        for atom in mol.GetAtoms():
            atoms.append({
                'idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'hybridization': str(atom.GetHybridization()),
                'formal_charge': atom.GetFormalCharge(),
                'num_hs': atom.GetTotalNumHs(),
                'is_aromatic': atom.GetIsAromatic(),
                'degree': atom.GetDegree(),
            })
        return atoms

    def _get_bond_info(self, mol):
        """Get bond information."""
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                'begin': bond.GetBeginAtomIdx(),
                'end': bond.GetEndAtomIdx(),
                'type': str(bond.GetBondType()),
                'is_in_ring': bond.IsInRing(),
            })
        return bonds

    def _detect_ring_systems(self, mol):
        """Detect ring information."""
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        ring_descs = []
        for ring in rings:
            size = len(ring)
            atoms_in_ring = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
            # Check aromaticity
            is_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
            ring_type = "aromatic" if is_aromatic else "aliphatic"
            ring_descs.append({
                'size': size,
                'atoms': atoms_in_ring,
                'type': ring_type,
            })
        return ring_descs

    def _describe_brief(self, mol):
        """Brief one-line description."""
        try:
            iupac = Chem.MolToIUPACName(mol)
        except Exception:
            iupac = None

        formula = rdMolDescriptors.CalcMolFormula(mol)
        ring_info = mol.GetRingInfo()
        n_rings = ring_info.NumRings()

        parts = []
        if iupac:
            parts.append(f"{iupac}")
        parts.append(f"Molecular formula: {formula}")
        if n_rings > 0:
            parts.append(f"Contains {n_rings} ring(s)")

        return ". ".join(parts) + "."

    def _describe_medium(self, mol):
        """Medium detail description."""
        atoms = self._get_atom_info(mol)
        bonds = self._get_bond_info(mol)
        rings = self._detect_ring_systems(mol)

        # Build atom-bond connectivity summary
        atom_symbols = [a['symbol'] for a in atoms]
        n_atoms = len(atoms)
        n_bonds = len(bonds)
        n_carbons = sum(1 for a in atoms if a['symbol'] == 'C')
        n_heteroatoms = sum(1 for a in atoms if a['symbol'] not in ('C', 'H'))

        # Build skeleton description
        lines = []

        # Try IUPAC name first
        try:
            iupac = Chem.MolToIUPACName(mol)
            if iupac:
                lines.append(f"{iupac}:")
        except Exception:
            pass

        # Main chain / core structure
        if n_carbons > 0 and n_carbons == n_atoms - n_heteroatoms or n_heteroatoms <= 2:
            # Likely organic molecule - describe as carbon skeleton
            carbon_chain = self._describe_carbon_skeleton(mol, atoms, bonds)
            if carbon_chain:
                lines.append(carbon_chain)

        # Ring info
        if rings:
            ring_parts = []
            for r in rings:
                ring_parts.append(f"{r['size']}-membered {r['type']} ring ({''.join(r['atoms'])})")
            lines.append(f"Contains {len(rings)} ring(s): {'; '.join(ring_parts)}")

        # Heteroatom features
        hetero_features = self._describe_heteroatoms(atoms, bonds)
        if hetero_features:
            lines.append(hetero_features)

        # Condensed structural formula
        condensed = self._condensed_formula(mol, atoms, bonds)
        if condensed:
            lines.append(f"Condensed structure: {condensed}")

        return "\n".join(lines) if lines else f"Molecule with {n_atoms} atoms and {n_bonds} bonds."

    def _describe_detailed(self, mol):
        """Detailed multi-section description."""
        atoms = self._get_atom_info(mol)
        bonds = self._get_bond_info(mol)
        rings = self._detect_ring_systems(mol)

        sections = []

        # Section 1: Overview
        try:
            iupac = Chem.MolToIUPACName(mol)
        except Exception:
            iupac = None

        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = round(Descriptors.MolWt(mol), 2)
        overview = f"Molecule: {iupac or 'Unknown'}\nFormula: {formula} | MW: {mw} g/mol | Atoms: {len(atoms)} | Bonds: {len(bonds)}"
        sections.append(overview)

        # Section 2: Atom-by-atom table
        atom_lines = ["\nAtoms:"]
        for a in atoms:
            h_str = f"H{a['num_hs']}" if a['num_hs'] > 0 else ""
            arom_str = " (aromatic)" if a['is_aromatic'] else ""
            charge_str = f" (charge={a['formal_charge']})" if a['formal_charge'] != 0 else ""
            atom_lines.append(
                f"  [{a['idx']}] {a['symbol']}{h_str} — "
                f"hyb:{a['hybridization'].split('.')[-1]}, "
                f"degree:{a['degree']}{arom_str}{charge_str}"
            )
        sections.append("\n".join(atom_lines))

        # Section 3: Bond table
        bond_lines = ["\nBonds:"]
        for b in bonds:
            a1_sym = mol.GetAtomWithIdx(b['begin']).GetSymbol()
            a2_sym = mol.GetAtomWithIdx(b['end']).GetSymbol()
            ring_str = " [ring]" if b['is_in_ring'] else ""
            bond_lines.append(
                f"  {a1_sym}({b['begin']})={b['type'].split('.')[-1]}={a2_sym}({b['end']}){ring_str}"
            )
        sections.append("\n".join(bond_lines))

        # Section 4: Ring systems
        if rings:
            ring_lines = ["\nRing Systems:"]
            for i, r in enumerate(rings):
                bond_atoms = ''.join(r['atoms'])
                ring_lines.append(
                    f"  Ring {i+1}: {r['size']}-membered {r['type']} — atoms: {bond_atoms}"
                )
            sections.append("\n".join(ring_lines))

        # Section 5: Functional group detection (basic)
        fg_detected = self._detect_functional_groups_simple(mol, atoms, bonds)
        if fg_detected:
            sections.append("\nFunctional Groups:\n  " + fg_detected)

        # Section 6: Condensed formula
        condensed = self._condensed_formula(mol, atoms, bonds)
        if condensed:
            sections.append(f"\nCondensed Structural Formula: {condensed}")

        return "\n".join(sections)

    def _describe_carbon_skeleton(self, mol, atoms, bonds):
        """Describe the carbon skeleton of an organic molecule."""
        carbons = [(i, a) for i, a in enumerate(atoms) if a['symbol'] == 'C']
        if not carbons:
            return None

        # Simple approach: build adjacency from bonds
        adj = {i: [] for i, _ in carbons}
        for b in bonds:
            if b['begin'] in adj and b['end'] in adj:
                adj[b['begin']].append(b['end'])
                adj[b['end']].append(b['begin'])

        # Describe based on carbon count
        n_c = len(carbons)
        c_desc = f"{n_c}-carbon"

        # Check for chains vs rings
        ring_info = mol.GetRingInfo()
        carbons_in_rings = set()
        for ring in ring_info.AtomRings():
            ring_carbons = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
            carbons_in_rings.update(ring_carbons)

        n_ring_c = len(carbons_in_rings & set(i for i, _ in carbons))
        n_chain_c = n_c - n_ring_c

        parts = []
        if n_ring_c > 0:
            parts.append(f"{n_ring_c} cyclic carbon(s)")
        if n_chain_c > 0:
            parts.append(f"{n_chain_c}-carbon chain")

        # Add substituent info
        hetero_atoms = [(i, a) for i, a in enumerate(atoms) if a['symbol'] not in ('C', 'H')]
        if hetero_atoms:
            substs = []
            for i, a in hetero_atoms:
                # Find what this heteroatom is bonded to
                neighbors = [b['end'] if b['begin'] == i else b['begin'] for b in bonds
                             if b['begin'] == i or b['end'] == i]
                c_neighbors = [n for n in neighbors if mol.GetAtomWithIdx(n).GetSymbol() == 'C']
                if c_neighbors:
                    substs.append(f"{a['symbol']}")
            if substs:
                parts.append(f"containing {', '.join(substs)}")

        return "A " + " with ".join(parts) + ("." if not parts[-1].endswith(".") else "")

    def _describe_heteroatoms(self, atoms, bonds):
        """Describe heteroatom features."""
        hetero = [(i, a) for i, a in enumerate(atoms) if a['symbol'] not in ('C', 'H')]
        if not hetero:
            return ""

        features = []
        for i, a in hetero:
            sym = a['symbol']
            hs = a['num_hs']
            if sym == 'O':
                if hs == 1:
                    features.append("hydroxyl group (-OH)")
                elif hs == 0:
                    # Could be carbonyl, ether, ester...
                    features.append("oxygen-containing group (possibly C=O, -O-, COOH)")
            elif sym == 'N':
                if hs >= 1:
                    features.append(f"amino group (-NH{'x' if hs > 1 else ''}{hs})")
                else:
                    features.append("nitrogen group (possibly cyano, nitro, amide)")
            elif sym == 'S':
                features.append("sulfur-containing group")
            elif sym in ('F', 'Cl', 'Br', 'I'):
                features.append(f"{sym} halogen substituent")
            elif sym == 'P':
                features.append("phosphorus group")

        return "Features: " + "; ".join(features) + "." if features else ""

    def _condensed_formula(self, mol, atoms, bonds):
        """Generate a simplified condensed structural formula."""
        # This is a simplified version — RDKit's MolToSmiles gives us something close
        try:
            return Chem.MolToSmiles(mol)
        except Exception:
            return None

    def _detect_functional_groups_simple(self, mol, atoms, bonds):
        """Simple functional group detection for detailed output."""
        groups = []
        smarts_patterns = {
            'alcohol': '[OH]',
            'aldehyde': '[CX3H1](=O)[#6]',
            'ketone': '[#6][CX3](=O)[#6]',
            'carboxylic acid': '[#6](=O)[OX2H1]',
            'ester': '[#6][CX3](=O)[OX2H0][#6]',
            'amine': '[NX3;H2,H1;!$(NC=O)]',
            'amide': '[NX3][CX3](=[OX1])',
            'nitrile': '[C-]#[N+]',
            'nitro': '[N+](=O)[O-]',
            'alkene': '[#6]=[#6]',
            'alkyne': '[#6]#[#6]',
            'ether': '[#6][OX2][#6]',
            'halide': '[F,Cl,Br,I]',
            'thiol': '[SH]',
            'phenyl': 'c1ccccc1',
        }

        for name, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matches = mol.GetSubstructMatches(pattern)
                groups.append(f"{name} ({len(matches)} found)")

        return ", ".join(groups) if groups else None

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split(None, 1)
        smiles = parts[0]
        detail = parts[1] if len(parts) > 1 else "medium"
        if detail not in ("brief", "medium", "detailed"):
            detail = "medium"
        return self._run_base(smiles, detail)


if __name__ == "__main__":
    run_mcp_server()
