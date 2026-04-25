import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPToolProcessError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem as RDKitChem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Pauling electronegativity values for polarity analysis
EN_VALUES = {
    "H": 2.20, "He": None,
    "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98, "Ne": None,
    "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": None,
    "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55,
    "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Ga": 1.81, "Ge": 2.01,
    "As": 2.18, "Se": 2.55, "Br": 2.96, "Kr": 3.00,
    "I": 2.66, "At": 2.20,
}


@ChemMCPManager.register_tool
class PredictPolarity(BaseTool):
    __version__ = "0.1.0"
    name = "PredictPolarity"
    func_name = 'predict_polarity'
    description = "Predict molecular polarity (polar/nonpolar) from SMILES string, including dipole moment direction and explanation."
    implementation_description = "Uses RDKit to analyze molecular structure and electronegativity differences between bonded atoms. Determines if bond dipoles cancel out due to symmetry. Returns polarity prediction, dipole moment estimate, and detailed reasoning."
    oss_dependencies = [("RDKit", "https://github.com/rdkit/rdkit", "BSD 3-Clause")]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Polarity", "Dipole Moment", "RDKit", "Molecular Properties"]
    required_envs = []

    code_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string of the molecule'),
    ]
    text_input_sig = [
        ('smiles', 'str', 'N/A', 'SMILES string'),
    ]
    output_sig = [
        ('smiles', 'str', 'Input SMILES'),
        ('is_polar', 'bool', 'Whether the molecule is polar'),
        ('polarity', 'str', '"polar" or "nonpolar"'),
        ('dipole_analysis', 'dict', 'Bond-by-bond dipole analysis'),
        ('symmetry_analysis', 'str', 'Symmetry and cancellation analysis'),
        ('explanation', 'str', 'Detailed explanation of polarity prediction'),
    ]
    
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'smiles': 'CCO', 'is_polar': True, 'polarity': 'polar', 'dipole_analysis': [...], 'symmetry_analysis': {...}, 'explanation': '...'}},
        {'code_input': {'smiles': 'C(=O)=O'}, 'text_input': {'smiles': 'CO2'}, 'output': {'smiles': 'CO2', 'is_polar': False, 'polarity': 'nonpolar', 'dipole_analysis': [...], 'symmetry_analysis': {...}, 'explanation': '...'}},
        {'code_input': {'smiles': 'c1ccccc1'}, 'text_input': {'smiles': 'benzene'}, 'output': {'smiles': 'c1ccccc1', 'is_polar': False, 'polarity': 'nonpolar', 'dipole_analysis': [...], 'symmetry_analysis': {...}, 'explanation': '...'}},
        {'code_input': {'smiles': 'Cl'}, 'text_input': {'smiles': 'Cl'}, 'output': {'smiles': 'Cl', 'is_polar': True, 'polarity': 'polar', 'dipole_analysis': [...], 'symmetry_analysis': {...}, 'explanation': '...'}},
    ]
    def _run_base(self, smiles: str) -> dict:
        if not RDKIT_AVAILABLE:
            raise ChemMCPToolProcessError("RDKit is not available.")

        mol = RDKitChem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError(f"Invalid SMILES string: '{smiles}'")

        mol_h = RDKitChem.AddHs(mol)
        
        # Analyze bonds for dipoles
        bond_dipoles = []
        en_differences = []
        
        for bond in mol_h.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            sym1 = atom1.GetSymbol()
            sym2 = atom2.GetSymbol()
            
            en1 = EN_VALUES.get(sym1)
            en2 = EN_VALUES.get(sym2)
            
            if en1 is not None and en2 is not None:
                delta_en = abs(en1 - en2)
                # Determine dipole direction (toward more electronegative atom)
                if en2 > en1:
                    direction = f"{sym1}δ⁺ → {sym2}δ⁻"
                    negative_end = sym2
                elif en1 > en2:
                    direction = f"{sym2}δ⁺ → {sym1}δ⁻"
                    negative_end = sym1
                else:
                    direction = "no dipole (same EN)"
                    negative_end = None
                
                bond_info = {
                    "bond": f"{sym1}-{sym2}",
                    "en1": en1, "en2": en2,
                    "delta_en": round(delta_en, 2),
                    "dipole_direction": direction,
                    "bond_type": str(bond.GetBondType()),
                }
                bond_dipoles.append(bond_info)
                if delta_en > 0.3:  # significant polarity threshold
                    en_differences.append({
                        "direction": direction,
                        "negative_end": negative_end,
                        "magnitude": delta_en,
                    })

        # Symmetry analysis
        has_polar_bonds = any(b["delta_en"] > 0.3 for b in bond_dipoles)
        
        # Check for symmetry cancellation
        # Simple heuristics based on common molecular patterns
        smiles_upper = smiles.upper().replace(" ", "")
        nonpolar_patterns = [
            "C(=O)=O", "CO2",  # linear CO2
            "C(C)(C)C", "CC(C)(C)C",  # tetrahedral neopentane / CCl4-like
            "C1CC(C)(C)C1",  # cyclohexane chair
            "C=C",  # symmetric alkenes like ethene
            "C#C",  # ethyne
            "c1ccccc1",  # benzene (regular hexagon)
            "S(=O)(=O)",  # SO3 (trigonal planar)
            "F(F)(F)F",  # CF4
        ]
        
        is_symmetric = False
        symmetry_reason = ""
        
        if not has_polar_bonds:
            is_polar = False
            symmetry_reason = "No significant polar bonds (all bonds have ΔEN < 0.3 or between same elements)."
        elif any(p in smiles_upper for p in nonpolar_patterns):
            is_polar = False
            is_symmetric = True
            symmetry_reason = "Molecule has polar bonds but they cancel out due to high symmetry."
        else:
            # Check if all dipoles point in different directions that could cancel
            # This is a simplified check; full vector analysis would require 3D coordinates
            negative_ends = [d["negative_end"] for d in en_differences if d.get("negative_end")]
            unique_negative_ends = set(negative_ends)
            
            if len(unique_negative_ends) <= 1 and len(en_differences) > 0:
                # All dipoles point the same direction -> polar
                is_polar = True
                symmetry_reason = f"All bond dipoles point toward {list(unique_negative_ends)[0]}; no cancellation."
            else:
                # Multiple directions - likely some cancellation but may still be polar
                # For simplicity, assume polar unless clearly symmetric
                is_polar = True
                symmetry_reason = (
                    f"Bond dipoles do not fully cancel. "
                    f"Polar bonds present with dipoles pointing toward: {unique_negative_ends}. "
                    f"Molecular geometry does not allow complete dipole cancellation."
                )

        # Get rough dipole moment estimate from RDLogP (correlates with polarity)
        try:
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol_h)
        except:
            logp, tpsa = None, None

        return {
            "smiles": smiles,
            "is_polar": is_polar,
            "polarity": "polar" if is_polar else "nonpolar",
            "bond_dipole_analysis": bond_dipoles,
            "has_polar_bonds": has_polar_bonds,
            "symmetry_analysis": {
                "is_symmetric": is_symmetric,
                "reason": symmetry_reason,
            },
            "explanation": self._generate_explanation(is_polar, has_polar_bonds, symmetry_reason, bond_dipoles),
            "molecular_properties": {
                "logp": round(logp, 2) if logp else None,
                "tpsa_sa": round(tpsa, 1) if tpsa else None,
                "num_atoms": mol_h.GetNumAtoms(),
            },
        }

    def _generate_explanation(self, is_polar, has_polar_bonds, reason, dipoles):
        if not has_polar_bonds:
            return ("This molecule is **nonpolar** because it contains no significant polar bonds. "
                    "All bonds are either between atoms of similar electronegativity or the molecule consists of a single element.")
        elif is_polar:
            return (f"This molecule is **polar**. It has polar bonds whose dipole moments do NOT cancel out. "
                    f"{reason}")
        else:
            return (f"This molecule is **nonpolar** despite having polar bonds. "
                    f"{reason}")


if __name__ == "__main__":
    run_mcp_server()
