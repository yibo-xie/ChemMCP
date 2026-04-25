"""
Bonding data module: bond lengths (pm), bond energies (kJ/mol), covalent radii.
Data from CRC Handbook / standard chemistry references.
"""

from typing import Dict, Optional, Tuple, List

# Average bond lengths in picometers (pm)
# Key: (element1, element2, bond_type) -> length_pm
BOND_LENGTHS: Dict[Tuple[str, str, str], float] = {
    # Single bonds
    ("H", "H", "single"): 74,
    ("H", "C", "single"): 109,
    ("H", "N", "single"): 101,
    ("H", "O", "single"): 96,
    ("H", "F", "single"): 92,
    ("H", "Si", "single"): 148,
    ("H", "S", "single"): 134,
    ("H", "P", "single"): 144,
    ("H", "Cl", "single"): 127,
    ("H", "Br", "single"): 141,
    ("H", "I", "single"): 161,
    ("C", "C", "single"): 154,
    ("C", "N", "single"): 147,
    ("C", "O", "single"): 143,
    ("C", "F", "single"): 135,
    ("C", "Si", "single"): 186,
    ("C", "P", "single"): 187,
    ("C", "S", "single"): 182,
    ("C", "Cl", "single"): 177,
    ("C", "Br", "single"): 194,
    ("C", "I", "single"): 214,
    ("N", "N", "single"): 145,
    ("N", "O", "single"): 136,
    ("N", "F", "single"): 136,
    ("N", "Cl", "single"): 175,
    ("N", "Br", "single"): 190,
    ("O", "O", "single"): 148,
    ("O", "F", "single"): 142,
    ("O", "Si", "single"): 166,
    ("O", "P", "single"): 163,
    ("O", "S", "single"): 151,
    ("O", "Cl", "single"): 170,
    ("O", "Br", "single"): 184,
    ("O", "I", "single"): 204,
    ("F", "F", "single"): 142,
    ("F", "Si", "single"): 156,
    ("F", "P", "single"): 156,
    ("F", "S", "single"): 158,
    ("F", "Cl", "single"): 166,
    ("F", "Br", "single"): 178,
    ("Si", "Si", "single"): 234,
    ("Si", "Cl", "single"): 202,
    ("Si", "Br", "single"): 215,
    ("P", "P", "single"): 221,
    ("P", "Cl", "single"): 204,
    ("P", "Br", "single"): 222,
    ("S", "S", "single"): 205,
    ("S", "Cl", "single"): 207,
    ("S", "Br", "single"): 225,
    ("Cl", "Cl", "single"): 199,
    ("Cl", "Br", "single"): 214,
    ("Br", "Br", "single"): 228,
    ("I", "I", "single"): 267,
    ("C", "B", "single"): 159,
    ("B", "H", "single"): 119,
    ("B", "F", "single"): 131,
    ("B", "Cl", "single"): 175,
    ("B", "O", "single"): 136,
    ("C", "Metal", "single"): 200,  # generic C-metal
    # Double bonds
    ("C", "C", "double"): 134,
    ("C", "N", "double"): 129,
    ("C", "O", "double"): 123,
    ("C", "S", "double"): 160,
    ("N", "N", "double"): 125,
    ("N", "O", "double"): 120,
    ("O", "O", "double"): 121,
    ("S", "O", "double"): 143,
    ("P", "O", "double"): 150,
    # Triple bonds
    ("C", "C", "triple"): 120,
    ("C", "N", "triple"): 116,
    ("C", "P", "triple"): 154,
    ("N", "N", "triple"): 110,
    # Aromatic (average)
    ("C", "C", "aromatic"): 139,
    ("C", "N", "aromatic"): 134,
    ("C", "O", "aromatic"): 136
}

# Bond dissociation energies in kJ/mol (average values at 298K)
# Key: bond_spec -> energy_kJ_mol
BOND_ENERGIES: Dict[str, float] = {
    # Single bonds
    "H-H": 436,
    "H-C": 413,
    "H-N": 391,
    "H-O": 463,
    "H-F": 567,
    "H-Si": 318,
    "H-S": 363,
    "H-P": 322,
    "H-Cl": 432,
    "H-Br": 366,
    "H-I": 299,
    "C-C": 347,
    "C-N": 305,
    "C-O": 358,
    "C-F": 485,
    "C-Si": 301,
    "C-P": 264,
    "C-S": 272,
    "C-Cl": 339,
    "C-Br": 276,
    "C-I": 238,
    "C-B": 356,
    "N-N": 163,
    "N-O": 201,
    "N-F": 270,
    "N-Cl": 200,
    "O-O": 146,
    "O-F": 190,
    "O-Si": 452,
    "O-P": 335,
    "O-S": 265,
    "O-Cl": 203,
    "O-Br": 243,
    "F-F": 158,
    "F-Si": 565,
    "F-P": 490,
    "F-S": 284,
    "F-Cl": 255,
    "Si-Si": 226,
    "Si-Cl": 381,
    "Si-O": 452,
    "P-P": 200,
    "P-Cl": 331,
    "S-S": 226,
    "S-Cl": 255,
    "S-Br": 243,
    "S-F": 284,
    "Cl-Cl": 242,
    "Cl-Br": 218,
    "Br-Br": 193,
    "I-I": 151,
    "B-H": 389,
    "B-F": 613,
    "B-Cl": 456,
    "B-O": 536,
    "B-C": 356,
    # Double bonds
    "C=C": 614,
    "C=N": 615,
    "C=O": 799,  # CO2 / ketone average; formaldehyde C=O is 745
    "C=S": 577,
    "N=N": 418,
    "N=O": 607,
    "O=O": 498,
    "P=O": 535,
    "S=O": 522,
    # Triple bonds
    "C≡C": 839,
    "C≡N": 891,
    "N≡N": 945,
    "C≡P": 552,
}

# Covalent radii in pm (Cordero et al., 2008)
COVALENT_RADII: Dict[str, float] = {
    "H": 31, "He": 28, "Li": 128, "Be": 96, "B": 84, "C": 76, "N": 71,
    "O": 66, "F": 57, "Ne": 58, "Na": 166, "Mg": 141, "Al": 121, "Si": 111,
    "P": 107, "S": 105, "Cl": 102, "Ar": 106, "K": 178, "Ca": 176, "Sc": 170,
    "Ti": 160, "V": 153, "Cr": 139, "Mn": 139, "Fe": 132, "Co": 126, "Ni": 124,
    "Cu": 132, "Zn": 122, "Ga": 122, "Ge": 120, "As": 119, "Se": 120, "Br": 120,
    "Kr": 116, "Rb": 195, "Sr": 190, "Y": 162, "Zr": 148, "Nb": 146, "Mo": 138,
    "Tc": 136, "Ru": 136, "Rh": 134, "Pd": 131, "Ag": 136, "Cd": 136, "In": 142,
    "Sn": 139, "Sb": 139, "Te": 138, "I": 139, "Xe": 140, "Cs": 224, "Ba": 198,
    "La": 172, "Ce": 166, "Pr": 164, "Nd": 164, "Pm": 164, "Sm": 163, "Eu": 166,
    "Gd": 165, "Tb": 164, "Dy": 163, "Ho": 162, "Er": 162, "Tm": 162, "Yb": 170,
    "Lu": 162, "Hf": 150, "Ta": 146, "W": 139, "Re": 137, "Os": 137, "Ir": 136,
    "Pt": 136, "Au": 136, "Hg": 132, "Tl": 145, "Pb": 146, "Bi": 148, "Po": 140,
    "At": 150, "Rn": 150, "Fr": 226, "Ra": 202, "Ac": 176, "Th": 170, "Pa": 164,
    "U": 168, "Np": 164, "Pu": 164, "Am": 166, "Cm": 166, "Bk": 166, "Cf": 166,
    "Es": 166, "Fm": 164, "Md": 164, "No": 164, "Lr": 162,
}

# Born exponents for lattice energy calculation
BORN_EXPONENTS: Dict[str, int] = {
    "He": 5, "Ne": 7, "Ar": 9, "Kr": 10, "Xe": 12,
    "Li+": 7, "Na+": 9, "K+": 10, "Rb+": 11, "Cs+": 12,
    "Be2+": 6, "Mg2+": 8, "Ca2+": 9, "Sr2+": 10, "Ba2+": 11,
    "Al3+": 9, "Sc3+": 9, "Y3+": 10, "La3+": 10,
    "Ti4+": 9, "Zr4+": 10, "Hf4+": 10,
    "O2-": 7, "S2-": 9, "F-": 7, "Cl-": 9, "Br-": 10, "I-": 11,
}

# Madelung constants for common crystal structure types
MADELUNG_CONSTANTS: Dict[str, float] = {
    "NaCl (rock salt)": 1.74756,
    "CsCl": 1.76267,
    "ZnS (zinc blende)": 1.63806,
    "ZnS (wurtzite)": 1.641,
    "CaF2 (fluorite)": 2.51939,
    "TiO2 (rutile)": 2.408,
    "Al2O3 (corundum)": 4.1719,
}


def get_bond_length(elem1: str, elem2: str, bond_type: str = "single") -> Optional[float]:
    """Get bond length in pm. Returns None if not found."""
    e1, e2 = elem1.capitalize(), elem2.capitalize()
    key = (e1, e2, bond_type)
    rkey = (e2, e1, bond_type)
    return BOND_LENGTHS.get(key) or BOND_LENGTHS.get(rkey)


def get_bond_energy(bond_spec: str) -> Optional[float]:
    """Get bond energy in kJ/mol from a bond spec like 'C-C', 'C=C', 'C≡C'."""
    return BOND_ENERGIES.get(bond_spec)


def get_covalent_radius(element: str) -> Optional[float]:
    """Get covalent radius in pm."""
    el = element.capitalize()
    if len(el) == 1:
        el = el.upper()
    elif el.lower() in ("fe", "co", "ni", "cu", "zn", "se"):
        el = {"fe":"Fe","co":"Co","ni":"Ni","cu":"Cu","zn":"Zn","se":"Se"}.get(el.lower(), el.title())
    return COVALENT_RADII.get(el)


def estimate_bond_length(elem1: str, elem2: str) -> float:
    """Estimate bond length as sum of covalent radii."""
    r1 = get_covalent_radius(elem1)
    r2 = get_covalent_radius(elem2)
    if r1 and r2:
        return round(r1 + r2, 1)
    raise ValueError(f"Cannot estimate bond length for {elem1}-{elem2}: missing covalent radius data")
