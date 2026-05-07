import logging
import math
from typing import Optional, List, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_AMU_TO_KG = 1.66054e-27  # kg per amu
_H = 6.62607015e-34  # J·s (Planck constant)
_C = 2.99792458e8  # m/s (speed of light)
_NA = 6.02214076e23  # mol^-1 (Avogadro's number)

# Conversion: rotational constant B = h / (8 * pi^2 * I * c) in cm^-1
def _I_to_rotational_constant(I_kg_m2: float) -> float:
    """Convert moment of inertia (kg·m^2) to rotational constant (cm^-1)."""
    if I_kg_m2 <= 0:
        return float('inf')
    return _H / (8.0 * math.pi ** 2 * I_kg_m2 * _C) * 100  # m^-1 to cm^-1


@ChemMCPManager.register_tool
class MomentOfInertia(BaseTool):
    """
    转动惯量计算工具 — 计算分子的主转动惯量、转动常数和分子分类。
    
    支持：双原子分子、线性三原子分子、对称/不对称陀螺分子。
    """
    __version__ = "0.1.0"
    name = "MomentOfInertia"
    func_name = "calculate_moment_of_inertia"
    description = "Calculate moment of inertia for molecules, determine rotational constants and molecular rotor classification."
    implementation_description = "Computes principal moments of inertia from atomic masses and geometry, classifies the molecule as spherical top, linear, symmetric top, or asymmetric top, and returns rotational constants A, B, C in cm^-1."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Moment of Inertia", "Rotational Spectroscopy", "Molecular Geometry", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("molecule_type", "str", "N/A", "Molecule type: 'diatomic', 'linear_polyatomic', 'spherical_top', 'symmetric_top', 'asymmetric_top'"),
        ("masses", "list", "N/A", "List of atomic masses in atomic mass units (amu), e.g., [12.0, 16.0]"),
        ("geometry", "dict", "N/A", "Geometry dict with keys depending on type. For diatomic: {'bond_length_A': 1.13}. For linear: {'bond_lengths_A': [1.13], 'angles': []}. For tops: {'coordinates': [(x,y,z), ...] in Angstroms}"),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: molecule_type|masses|geometry_json. Example: 'diatomic|[12,16]|{\"bond_length_A\":1.13}'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing: Ia, Ib, Ic (kg·m^2), A, B, C (cm^-1), classification, reduced_mass (amu) if diatomic."),
    ]

    examples = [
        {
            "code_input": {
                "molecule_type": "diatomic",
                "masses": [12.0, 16.0],
                "geometry": {"bond_length_A": 1.13},
            },
            "text_input": {
                "input_str": "diatomic|[12,16]|{\"bond_length_A\":1.13}"
            },
            "output": {
                "result": {
                    "classification": "linear",
                    "I_perp": "<value>",
                    "I_parallel": 0,
                    "B": "<value>",
                    "reduced_mass": 6.857,
                }
            },
        },
        {
            "code_input": {
                "molecule_type": "spherical_top",
                "masses": [12.0, 1.008, 1.008, 1.008, 1.008],
                "geometry": {"coordinates": [(0, 0, 0), (0, 0, 1.09), (0, 1.88, -0.36), (0, -1.88, -0.36), (1.21, 0, -0.36)]},
            },
            "text_input": {
                "input_str": "spherical_top|[12,1.008,1.008,1.008,1.008]|{\"coordinates\":[[0,0,0],[0,0,1.09],[0,1.88,-0.36],[0,-1.88,-0.36],[1.21,0,-0.36]]}"
            },
            "output": {
                "result": {
                    "classification": "spherical top",
                    "Ia": "<value>",
                    "Ib": "<value>",
                    "Ic": "<value>",
                    "A": "<value>",
                    "B": "<value>",
                    "C": "<value>",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_diatomic(self, masses: List[float], geometry: dict) -> dict:
        """Diatomic molecule: I = mu * r^2"""
        m1, m2 = masses[0], masses[1]
        r_angstrom = geometry.get("bond_length_A", 1.0)
        r_m = r_angstrom * 1e-10
        mu_amu = (m1 * m2) / (m1 + m2)
        mu_kg = mu_amu * _AMU_TO_KG
        I = mu_kg * r_m ** 2
        B_cm = _I_to_rotational_constant(I)
        return {
            "classification": "linear",
            "I_parallel": 0.0,
            "I_perp": I,
            "Ia": 0.0,
            "Ib": I,
            "Ic": I,
            "A": float('inf'),
            "B": round(B_cm, 4),
            "C": round(B_cm, 4),
            "reduced_mass_amu": round(mu_amu, 4),
            "bond_length_angstrom": r_angstrom,
        }

    def _calc_linear_polyatomic(self, masses: List[float], geometry: dict) -> dict:
        """Linear polyatomic molecule."""
        coords = geometry.get("coordinates", [])
        bond_lengths = geometry.get("bond_lengths_A", [])
        if not coords and bond_lengths:
            # Build coordinates along z-axis
            coords = [[0.0, 0.0, 0.0]]
            z = 0.0
            for bl in bond_lengths:
                z += bl
                coords.append([0.0, 0.0, z])
        return self._calc_from_coords(masses, coords, is_linear=True)

    def _calc_from_coords(self, masses: List[float], coords: list, is_linear: bool = False) -> dict:
        """General calculation from Cartesian coordinates in Angstroms."""
        n = len(masses)
        if n != len(coords):
            raise ChemMCPError(f"Number of masses ({n}) doesn't match number of coordinates ({len(coords)}).")

        # Center of mass
        total_mass = sum(masses)
        com = [0.0, 0.0, 0.0]
        for i in range(n):
            for j in range(3):
                com[j] += masses[i] * coords[i][j]
        com = [c / total_mass for c in com]

        # Inertia tensor
        tensor = [[0.0] * 3 for _ in range(3)]
        for i in range(n):
            m = masses[i] * _AMU_TO_KG
            x = (coords[i][0] - com[0]) * 1e-10
            y = (coords[i][1] - com[1]) * 1e-10
            z = (coords[i][2] - com[2]) * 1e-10
            r2 = x*x + y*y + z*z
            tensor[0][0] += m * (r2 - x*x)
            tensor[1][1] += m * (r2 - y*y)
            tensor[2][2] += m * (r2 - z*z)
            tensor[0][1] -= m * x * y
            tensor[0][2] -= m * x * z
            tensor[1][2] -= m * y * z

        tensor[1][0] = tensor[0][1]
        tensor[2][0] = tensor[0][2]
        tensor[2][1] = tensor[1][2]

        # Eigenvalues of inertia tensor → principal moments
        eigenvalues = self._symmetric_eigen3x3(tensor)
        Ia, Ib, Ic = sorted(eigenvalues)

        # Classification
        eps = 1e-50
        if is_linear:
            classification = "linear"
        elif abs(Ia - Ib) < eps * max(Ia, 1) and abs(Ib - Ic) < eps * max(Ib, 1):
            classification = "spherical top"
        elif abs(Ia - Ib) < eps * max(Ia, 1):
            classification = "prolate symmetric top"  # Ia ≈ Ib < Ic
        elif abs(Ib - Ic) < eps * max(Ib, 1):
            classification = "oblate symmetric top"   # Ia < Ib ≈ Ic
        else:
            classification = "asymmetric top"

        A_cm = _I_to_rotational_constant(Ia) if Ia > eps else float('inf')
        B_cm = _I_to_rotational_constant(Ib) if Ib > eps else float('inf')
        C_cm = _I_to_rotational_constant(Ic) if Ic > eps else float('inf')

        return {
            "classification": classification,
            "Ia": round(Ia, 40),
            "Ib": round(Ib, 40),
            "Ic": round(Ic, 40),
            "A": round(A_cm, 4) if A_cm != float('inf') else A_cm,
            "B": round(B_cm, 4) if B_cm != float('inf') else B_cm,
            "C": round(C_cm, 4) if C_cm != float('inf') else C_cm,
        }

    def _symmetric_eigen3x3(self, M: list) -> tuple:
        """Compute eigenvalues of a symmetric 3x3 matrix using analytical method."""
        a, b, c = M[0][0], M[1][1], M[2][2]
        d, e, f = M[0][1], M[0][2], M[1][2]

        p = -(a + b + c)
        q = a*b + b*c + c*a - d*d - e*e - f*f
        r = -(a*b*c + 2*d*e*f - a*f*f - b*e*e - c*d*d)

        # Solve cubic: x^3 + p*x^2 + q*x + r = 0
        Q = (3*q - p*p) / 9.0
        R = (9*p*q - 27*r - 2*p**3) / 54.0
        D = Q**3 + R**3

        sqrt_D = math.sqrt(abs(D)) if D >= 0 else math.sqrt(-D)
        theta = math.acos(R / math.sqrt(-Q**3)) if Q < 0 else 0.0

        if D <= 0:
            sq = math.sqrt(-Q)
            root1 = 2*sq*math.cos(theta/3.0) - p/3.0
            root2 = 2*sq*math.cos((theta+2*math.pi)/3.0) - p/3.0
            root3 = 2*sq*math.cos((theta+4*math.pi)/3.0) - p/3.0
        else:
            S = math.cbrt(R + sqrt_D)
            T = math.cbrt(R - sqrt_D)
            root1 = (S + T) - p/3.0
            real_part = -(S+T)/2.0 - p/3.0
            imag_part = (math.sqrt(3)/2.0)*(S-T)
            root2 = real_part
            root3 = real_part

        return (root1, root2, root3)

    def _run_base(self, molecule_type: str, masses: List[float], geometry: dict) -> dict:
        """Core logic: calculate moment of inertia based on molecule type."""
        mol_type = molecule_type.lower().strip()

        if mol_type == "diatomic":
            return self._calc_diatomic(masses, geometry)
        elif mol_type == "linear_polyatomic":
            return self._calc_linear_polyatomic(masses, geometry)
        elif mol_type in ("spherical_top", "symmetric_top", "asymmetric_top"):
            coords = geometry.get("coordinates", [])
            return self._calc_from_coords(masses, coords)
        else:
            raise ChemMCPError(
                f"Unsupported molecule type: '{molecule_type}'. "
                f"Supported: 'diatomic', 'linear_polyatomic', 'spherical_top', 'symmetric_top', 'asymmetric_top'."
            )

    def _run_text(self, input_str: str) -> dict:
        """Parse text input like: diatomic|[12,16]|{'bond_length_A':1.13}"""
        try:
            parts = input_str.split("|", 2)
            if len(parts) < 3:
                raise ValueError("Expected format: molecule_type|masses|geometry_json")
            
            mol_type = parts[0].strip()
            import json
            masses = json.loads(parts[1].replace("'", '"'))
            geometry = json.loads(parts[2].replace("'", '"'))
            
            return self._run_base(mol_type, masses, geometry)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'molecule_type|[m1,m2,...]|<json_geometry>'")
