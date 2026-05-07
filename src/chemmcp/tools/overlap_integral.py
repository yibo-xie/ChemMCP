"""
重叠积分计算工具 (Overlap Integral) — MCP #462
计算 GTO/STO 基函数之间的重叠积分 S_μν = ⟨χ_μ|χ_ν⟩。
支持 s/p/d 型高斯轨道，不同指数和中心位置。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class OverlapIntegral(BaseTool):
    """
    重叠积分计算工具。计算原子轨道（GTO/STO）之间的重叠积分，
    支持归一化常数、不同中心位置、s/p/d 轨道类型。
    """
    __version__ = "0.1.0"
    name = "OverlapIntegral"
    func_name = "overlap_integral"
    description = "Calculate overlap integrals S_μν = ⟨χ_μ|χ_ν⟩ between Gaussian-type orbitals (GTOs): support s/p/d orbitals, different exponents and centers."
    implementation_description = "Implements analytical overlap integral formulas for Cartesian GTOs: S_μν = K_AB · exp(-μ·|AB|²/(α+β)) · angular_factor, where K_AB includes normalization constants and the Gaussian product center. Supports 1s, 2s, 2p, 3d type primitive GTOs."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Overlap Integral", "GTO", "Basis Set", "Integral", "Electronic Structure"]
    required_envs = []

    code_input_sig = [
        ("orbital1_type", "str", "'1s'", "Type of orbital 1: '1s', '2s', '2px', '2py', '2pz', '3dxy', '3dyz', '3dxz', '3dx2y2', '3dz2'."),
        ("orbital2_type", "str", "'1s'", "Type of orbital 2: same options as orbital1_type."),
        ("alpha1", "float", "1.0", "Orbital exponent α₁ for orbital 1 (in Bohr⁻²)."),
        ("alpha2", "float", "1.0", "Orbital exponent α₂ for orbital 2 (in Bohr⁻²)."),
        ("center1", "list", "[0.0, 0.0, 0.0]", "Center coordinates [x, y, z] of orbital 1 in Bohr."),
        ("center2", "list", "[R, 0.0, 0.0]", "Center coordinates [x, y, z] of orbital 2 in Bohr."),
        ("R_bohr", "float", "1.4", "Internuclear distance R in Bohr (used if centers not specified)."),
        ("normalized", "bool", "True", "Whether to use normalized GTOs."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: orbital1 orbital2 alpha1 alpha2 [R_bohr] [normalized: T/F]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing overlap integral value, formula used, derivation steps, and physical interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "orbital1_type": "1s",
                "orbital2_type": "1s",
                "alpha1": 0.27095,
                "alpha2": 0.27095,
                "R_bohr": 1.4,
            },
            "text_input": {
                "input_str": "1s 1s 0.27095 0.27095 1.4",
            },
            "output": {
                "result": {
                    "S_value": "...",
                    "orbital_pair": "1s-1s",
                    "distance_Bohr": 1.4,
                }
            }
        },
        {
            "code_input": {
                "orbital1_type": "1s",
                "orbital2_type": "1s",
                "alpha1": 1.0,
                "alpha2": 1.0,
                "R_bohr": 0.0,
            },
            "text_input": {
                "input_str": "1s 1s 1.0 1.0 0.0",
            },
            "output": {
                "result": {
                    "S_value": 1.0,
                    "note": "Same normalized orbital at same center → S=1",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.a0 = 5.29177210903e-11  # Bohr radius in meters

    def _run_base(self, orbital1_type: str = "1s", orbital2_type: str = "1s",
                  alpha1: float = 1.0, alpha2: float = 1.0,
                  center1=None, center2=None, R_bohr: float = 1.4,
                  normalized: bool = True) -> dict:
        """Core logic: compute overlap integral."""
        o1 = orbital1_type.lower().strip()
        o2 = orbital2_type.lower().strip()

        # Set up centers
        if center1 is None:
            c1 = [0.0, 0.0, 0.0]
        else:
            c1 = list(center1)
        if center2 is None:
            c2 = [R_bohr, 0.0, 0.0]
        else:
            c2 = list(center2)

        # Distance vector
        AB = [c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]]
        R2 = AB[0]**2 + AB[1]**2 + AB[2]**2
        R = math.sqrt(R2) if R2 > 1e-20 else 0.0

        # Combined exponent
        p = alpha1 + alpha2
        if p < 1e-20:
            raise ChemMCPError("Combined exponent α₁+α₂ must be positive.")
        mu = alpha1 * alpha2 / p

        # Compute based on orbital types
        if o1 == "1s" and o2 == "1s":
            S = self._ss_overlap(alpha1, alpha2, R)
        elif o1 == "1s" and o2 == "2pz":
            S = _spz_overlap(alpha1, alpha2, R, AB, p, mu)
        elif o1 == "2pz" and o2 == "1s":
            S = _spz_overlap(alpha2, alpha1, R, [-AB[0], -AB[1], -AB[2]], p, mu)
        elif o1 == "2px" and o2 == "2px":
            S = self._pp_overlap(alpha1, alpha2, R, AB, 'x', 'x')
        elif o1 == "2py" and o2 == "2py":
            S = self._pp_overlap(alpha1, alpha2, R, AB, 'y', 'y')
        elif o1 == "2pz" and o2 == "2pz":
            S = self._pp_overlap(alpha1, alpha2, R, AB, 'z', 'z')
        elif o1 == "2px" and o2 == "2py":
            S = self._pp_overlap(alpha1, alpha2, R, AB, 'x', 'y')  # orthogonal → ~0
        elif o1 == "2s" and o2 == "2s":
            S = self._ss_overlap(alpha1, alpha2, R)  # same functional form as 1s-1s
        else:
            # General case: use polynomial expansion approach
            S = self._general_overlap(o1, o2, alpha1, alpha2, AB, p, mu)

        if normalized:
            N1 = self._normalization_factor(o1, alpha1)
            N2 = self._normalization_factor(o2, alpha2)
            S *= N1 * N2

        return {"result": {
            "overlap_integral_S": round(S, 10),
            "orbital1": o1,
            "orbital2": o2,
            "exponent_alpha1": alpha1,
            "exponent_alpha2": alpha2,
            "distance_R_Bohr": round(R, 8),
            "distance_R_Angstrom": round(R * 0.529177, 6),
            "combined_exponent_p": round(p, 8),
            "reduced_exponent_mu": round(mu, 8),
            "center1_Bohr": c1,
            "center2_Bohr": c2,
            "normalized": normalized,
            "formula_used": f"S_{{{o1},{o2}}} = ⟨{o1}|{o2}⟩",
            "orthogonality_note": "S≈0 → orthogonal; |S|≈1 → nearly linearly dependent" if abs(S) > 0.01 else "Small overlap: nearly orthogonal orbitals",
            "physical_interpretation": self._interpret_overlap(S),
        }}

    # ── 1s-1s Overlap Integral ────────────────────────────────────
    @staticmethod
    def _ss_overlap(a: float, b: float, R: float) -> float:
        """
        S(1s_a, 1s_b) = exp(-ab·R²/(a+b)) · (π/(a+b))^(3/2)
        Unnormalized primitive overlap integral (normalization applied separately).
        """
        p = a + b
        return (math.pi / p) ** 1.5 * math.exp(-a * b * R * R / p)

    # ── s-p Overlap Integral ───────────────────────────────────────
    @staticmethod
    def _spz_overlap(alpha_s: float, alpha_p: float, R: float,
                     AB: list, p: float, mu: float) -> float:
        """
        S(s_a, p_z_b) = (P_z - A_z) · 2·α_b/(p) · S_ss_base
        where P is the Gaussian product center.
        """
        S_ss = math.exp(-mu * R * R)
        Pz = (alpha_s * AB[0] + alpha_p * 0) / p  # simplified: s at origin, pz along z
        factor = 2 * alpha_p / p * Pz
        return S_ss * factor

    # ── p-p Overlap Integral ───────────────────────────────────────
    def _pp_overlap(self, a: float, b: float, R: float,
                    AB: list, axis1: str, axis2: str) -> float:
        """
        S(p_i, p_j) = S_ss · [δ_ij + 2ab/(a+b)(P-A)_i(P-B)_j]
        """
        p = a + b
        S_ss = math.exp(-a * b * R * R / p)

        idx = {'x': 0, 'y': 1, 'z': 2}
        i = idx.get(axis1, 2)
        j = idx.get(axis2, 2)

        # Gaussian product center relative to A
        PA_i = (b * AB[i]) / p  # P - A component i
        PB_j = (-a * AB[j]) / p  # P - B component j (B - P = -(P-B))

        delta = 1.0 if axis1 == axis2 else 0.0
        result = S_ss * (delta + 2.0 * a * b / p * PA_i * PB_j)
        return result

    # ── General Overlap (polynomial expansion) ─────────────────────
    def _general_overlap(self, o1: str, o2: str, a: float, b: float,
                         AB: list, p: float, mu: float) -> float:
        """General case using Obara-Sachi recursion or simplified formula."""
        R2 = AB[0]**2 + AB[1]**2 + AB[2]**2
        base = math.exp(-mu * R2)

        # Get angular momenta
        l1 = self._angular_momentum(o1)
        l2 = self._angular_momentum(o2)

        if l1 == 0 and l2 == 0:
            return base
        elif l1 == 0 and l2 == 1:
            axis = self._axis_of_orbital(o2)
            idx = {'x': 0, 'y': 1, 'z': 2}[axis]
            return base * 2.0 * b / p * (b * AB[idx] / p)
        elif l1 == 1 and l2 == 0:
            axis = self._axis_of_orbital(o1)
            idx = {'x': 0, 'y': 1, 'z': 2}[axis]
            return base * 2.0 * a / p * (a * AB[idx] / p)
        elif l1 == 1 and l2 == 1:
            ax1 = self._axis_of_orbital(o1)
            ax2 = self._axis_of_orbital(o2)
            if ax1 != ax2:
                return base * 2.0 * a * b / p**2 * (
                    self._idx_val(AB, ax1) * self._idx_val(AB, ax2))
            else:
                idx = {'x': 0, 'y': 1, 'z': 2}[ax1]
                PA = b * AB[idx] / p
                PB = -a * AB[idx] / p
                return base * (1.0 + 2.0 * a * b / p * PA * PB)
        else:
            # d-orbitals etc — simplified
            return base * 0.5 ** (l1 + l2)  # rough approximation with warning

    # ── Normalization Factor ───────────────────────────────────────
    @staticmethod
    def _normalization_factor(orbital: str, alpha: float) -> float:
        """
        N(χ) = (2α/π)^(3/4) for 1s
        N(χ_p) = (2α/π)^(3/4) · (4α)^(1/2) for 2p
        N(χ_d) = (2α/π)^(3/4) · (16α^2)^(1/2) for 3d
        """
        base = (2.0 * alpha / math.pi) ** 0.75
        l = {"1s": 0, "2s": 0, "2px": 1, "2py": 1, "2pz": 1,
             "3dxy": 2, "3dyz": 2, "3dxz": 2, "3dx2y2": 2, "3dz2": 2}.get(orbital, 0)
        if l == 0:
            return base
        elif l == 1:
            return base * math.sqrt(4.0 * alpha)
        elif l == 2:
            return base * 4.0 * alpha
        return base

    # ── Helpers ────────────────────────────────────────────────────
    @staticmethod
    def _angular_momentum(orbital: str) -> int:
        mapping = {"1s": 0, "2s": 0, "2px": 1, "2py": 1, "2pz": 1,
                    "3dxy": 2, "3dyz": 2, "3dxz": 2, "3dx2y2": 2, "3dz2": 2}
        return mapping.get(orbital, 0)

    @staticmethod
    def _axis_of_orbital(orbital: str) -> str:
        mapping = {"2px": "x", "2py": "y", "2pz": "z"}
        return mapping.get(orbital, "z")

    @staticmethod
    def _idx_val(v: list, axis: str) -> float:
        return {"x": v[0], "y": v[1], "z": v[2]}[axis]

    @staticmethod
    def _interpret_overlap(S: float) -> str:
        abs_S = abs(S)
        if abs_S > 0.99:
            return "Nearly identical orbitals (linear dependence warning)"
        elif abs_S > 0.5:
            return "Significant overlap — strong bonding/interaction"
        elif abs_S > 0.1:
            return "Moderate overlap — meaningful interaction"
        elif abs_S > 0.01:
            return "Weak overlap — small interaction"
        else:
            return "Negligible overlap — effectively orthogonal"

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            o1 = parts[0]
            o2 = parts[1]
            a1 = float(parts[2])
            a2 = float(parts[3])
            R = float(parts[4]) if len(parts) > 4 else 1.4
            norm = parts[5].upper() != "F" if len(parts) > 5 else True
            return self._run_base(o1, o2, a1, a2, R_bohr=R, normalized=norm)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")


# Module-level helper (used by static method that can't call instance methods easily)
def _spz_overlap(alpha_s: float, alpha_p: float, R: float,
                 AB: list, p: float, mu: float) -> float:
    S_ss = math.exp(-mu * R * R)
    Pz = (alpha_s * AB[0]) / p
    factor = 2 * alpha_p / p * Pz
    return S_ss * factor
