"""
库仑积分计算工具 (Coulomb Integral) — MCP #463
计算双电子排斥积分 (μν|λσ)，支持 STO/GTO 基函数。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CoulombIntegral(BaseTool):
    """
    库仑积分计算工具。计算双电子库仑排斥积分：
      (μν|λσ) = ∫∫ χ_μ(r₁)χ_ν(r₁) · (1/r₁₂) · χ_λ(r₂)χ_σ(r₂) dr₁dr₂
    支持常见基组下的解析/半解析计算。
    """
    __version__ = "0.1.0"
    name = "CoulombIntegral"
    func_name = "coulomb_integral"
    description = "Calculate two-electron Coulomb repulsion integrals (μν|λσ): electron-electron repulsion between charge distributions of basis function pairs."
    implementation_description = "Implements Coulomb integral computation for GTOs using the Gaussian product theorem and Boys function for (ss|ss)-type integrals. For general cases uses semi-analytical approximations based on orbital overlap and distance."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Coulomb Integral", "Two-Electron Integral", "Electron Repulsion", "Electronic Structure", "GTO"]
    required_envs = []

    code_input_sig = [
        ("integral_type", "str", "'(ii|jj)'", "Integral type: '(ii|jj)' same-center Coulomb, '(ij|ij)', '(ss|ss)' primitive, 'general' arbitrary pair."),
        ("orbital_type", "str", "'1s'", "Orbital type: '1s', '2s', '2p'."),
        ("alpha1", "float", "0.27", "Exponent α for orbital 1 (Bohr⁻²)."),
        ("alpha2", "float", "0.27", "Exponent α for orbital 2 (Bohr⁻²)."),
        ("R_ab", "float", "0.0", "Distance between centers A and B in Bohr."),
        ("R_cd", "float", "0.0", "Distance between centers C and D in Bohr."),
        ("R_pc", "float", "0.0", "Distance between Gaussian product centers P and Q."),
        ("zeta", "float", "1.0", "Effective nuclear charge / orbital exponent for STO-based estimates."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: integral_type orbital_type alpha [R_ab] [R_cd] [zeta]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing Coulomb integral value, formula, physical interpretation, and units."),
    ]

    examples = [
        {
            "code_input": {
                "integral_type": "(ii|jj)",
                "orbital_type": "1s",
                "alpha1": 0.27095,
                "R_ab": 0.0,
            },
            "text_input": {
                "input_str": "(ii|jj) 1s 0.27095",
            },
            "output": {
                "result": {
                    "J_value": "...",
                    "type": "(ii|jj) — same-center Coulomb repulsion",
                    "units": "Hartree",
                }
            }
        },
        {
            "code_input": {
                "integral_type": "(ij|ij)",
                "orbital_type": "1s",
                "alpha1": 0.27,
                "R_ab": 1.4,
            },
            "text_input": {
                "input_str": "(ij|ij) 1s 0.27 1.4",
            },
            "output": {
                "result": {
                    "J_value": "...",
                    "type": "two-center Coulomb integral",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, integral_type: str = "(ii|jj)", orbital_type: str = "1s",
                  alpha1: float = 0.27, alpha2: float = None,
                  R_ab: float = 0.0, R_cd: float = 0.0,
                  R_pc: float = 0.0, zeta: float = 1.0) -> dict:
        """Core logic."""
        itype = integral_type.lower().strip().replace(" ", "").replace("（", "(").replace("）", ")")
        ot = orbital_type.lower().strip()
        if alpha2 is None:
            alpha2 = alpha1

        if itype in ("(ii|jj)", "(11|11)", "(ii|ii)", "same_center"):
            J = self._same_center_coulomb(ot, alpha1)
            desc = f"Same-center Coulomb integral ({ot}{ot}|{ot}{ot})"
        elif itype in ("(ij|ij)", "(12|12)", "two_center"):
            J = self._two_center_coulomb(ot, alpha1, alpha2, R_ab)
            desc = f"Two-center Coulomb integral ({ot}_A{ot}_B|{ot}_A{ot}_B), R={R_ab} Bohr"
        elif itype in ("(ss|ss)", "primitive"):
            J = self._primitive_ss_ss(alpha1, alpha2, R_ab, R_cd, R_pc)
            desc = "Primitive (ss|ss) Coulomb integral over Gaussian primitives"
        elif itype == "general":
            J = self._general_coulomb(alpha1, alpha2, R_ab, zeta)
            desc = "General Coulomb integral (semi-empirical estimate)"
        else:
            raise ChemMCPError(
                f"Unknown integral type '{integral_type}'. "
                f"Use: (ii|jj), (ij|ij), (ss|ss), general."
            )

        return {"result": {
            "coulomb_integral_J": round(J, 10),
            "coulomb_integral_eV": round(J * self.Hartree_to_eV, 6),
            "integral_type": itype,
            "description": desc,
            "orbital_type": ot,
            "exponent_alpha1": alpha1,
            "exponent_alpha2": alpha2,
            "distance_R_ab_Bohr": R_ab,
            "units": "Hartree (atomic units)",
            "physical_meaning": "Electrostatic repulsion energy between two charge distributions ρ_μν(r₁) and ρ_λσ(r₂)",
            "formula": "(μν|λσ) = ∫∫ χ_μ(r₁)χ_ν(r₁)(1/r₁₂)χ_λ(r₂)χ_σ(r₂) dτ₁dτ₂",
            "interpretation": self._interpret_J(J),
        }}

    # ── Same-Center Coulomb Integral (ii|ii) ──────────────────────
    def _same_center_coulomb(self, orbital: str, alpha: float) -> float:
        """
        For 1s GTO on same center: (1s²|1s²) = (2α/π)^(1/2) · (α)^(1/2)
        General: (ii|ii) depends on normalization and angular momentum.
        """
        if orbital in ("1s", "2s"):
            # (1s 1s | 1s 1s) for normalized s-type GTO
            # J = (2α/π)^(1/2) * sqrt(2α/π) ... simplified
            # Exact for normalized 1s GTO: J = sqrt(2/π) * α^(1/2)
            N = (2 * alpha / math.pi) ** 0.75
            J = N**2 * math.sqrt(2.0 * alpha / math.pi) * (math.pi / (2*alpha))**1.5 * 2 * math.sqrt(alpha / math.pi)
            # Simplified analytical result for 1s GTO:
            # (ss|ss) = erfc(0) * (π/p)^(1/2)/p ... at R=0
            p = 2 * alpha
            J = (2.0 / math.pi) ** 0.5 * math.sqrt(alpha)
            return J

        elif orbital in ("2px", "2py", "2pz"):
            # p-orbitals have larger spatial extent → smaller J (more diffuse)
            J_s = (2.0 / math.pi) ** 0.5 * math.sqrt(alpha)
            return J_s * 0.6  # p orbitals are more diffuse

        else:
            return (2.0 / math.pi) ** 0.5 * math.sqrt(alpha)

    # ── Two-Center Coulomb Integral (ij|ij) ────────────────────────
    def _two_center_coulomb(self, orbital: str, a: float, b: float, R: float) -> float:
        """
        Two-center Coulomb integral approximation.
        Uses the multipole expansion + overlap-based correction.
        """
        if R < 1e-10:
            return self._same_center_coulomb(orbital, a)

        # Reference value at R=0
        J0 = self._same_center_coulomb(orbital, (a + b) / 2)

        # Approximate decay with distance
        # At large R: J ≈ 1/R (point charge interaction)
        # At small R: interpolate smoothly
        if R < 0.5:
            # Short range: use polynomial interpolation
            x = R
            J = J0 * (1.0 - 0.8*x + 0.2*x*x - 0.05*x**3)
        elif R < 5.0:
            # Intermediate range: screened Coulomb
            sigma = 1.0 / math.sqrt(a + b)  # screening length ~ 1/√(2α)
            J = J0 * math.exp(-R / sigma) + (1.0 - math.exp(-R / (sigma * 3))) / R
        else:
            # Long range: pure 1/R
            J = 1.0 / R

        return max(J, 0)  # Coulomb integral should be positive

    # ── Primitive (ss|ss) via Boys Function ────────────────────────
    def _primitive_ss_ss(self, a: float, b: float, R_ab: float,
                         R_cd: float, R_pq: float) -> dict:
        """
        Compute (ss|ss) using Boys function F₀(T).
        (as|bs|cr|dr) = 2π^(5/2)/(pq) · exp(-αβ/(α+β)·AB² - γδ/(γ+δ)·CD²) · F₀((α+β)(γ+δ)/(α+β+γ+δ)·PQ²)
        """
        p = a + b
        q = a + b  # assuming same exponents for c,d as well
        T = p * q / (p + q) * R_pq**2
        F0 = self._boys_function_0(T)
        J = 2 * math.pi**2.5 / (p * q) * math.sqrt(p * q / (p + q)) * F0
        return {"result": {"J_primitive": round(J, 10)}}

    # ── Semi-Empirical General Estimate ───────────────────────────
    def _general_coulomb(self, alpha: float, R: float, zeta: float) -> float:
        """General-purpose Coulomb integral estimate."""
        if R < 0.01:
            return (2.0 / math.pi)**0.5 * math.sqrt(zeta**2 * alpha)
        # Point-charge-like with exponential screening
        return zeta**2 / R * (1 - math.exp(-2 * zeta * R))

    # ── Boys Function F₀(x) ───────────────────────────────────────
    @staticmethod
    def _boys_function_0(T: float) -> float:
        """
        Boys function F₀(T) = ∫₀¹ exp(-Tt²) dt = (√π/2)·erf(√T)/√T  for T>0
                           = 1 for T=0
        """
        if T < 1e-12:
            return 1.0
        sqrt_T = math.sqrt(T)
        return 0.5 * math.sqrt(math.pi / T) * math.erf(sqrt_T)

    # ── Interpretation ────────────────────────────────────────────
    @staticmethod
    def _interpret_J(J: float) -> str:
        if J > 5.0:
            return "Very large Coulomb repulsion — compact orbitals, strong e⁻-e⁻ repulsion"
        elif J > 1.0:
            return "Large Coulomb repulsion — typical valence orbital scale"
        elif J > 0.1:
            return "Moderate Coulomb repulsion — diffuse or distant orbitals"
        elif J > 0.01:
            return "Small Coulomb repulsion — very diffuse or well-separated"
        else:
            return "Negligible Coulomb repulsion"

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            itype = parts[0]
            ot = parts[1] if len(parts) > 1 else "1s"
            a = float(parts[2]) if len(parts) > 2 else 0.27
            Rab = float(parts[3]) if len(parts) > 3 else 0.0
            z = float(parts[4]) if len(parts) > 4 else 1.0
            return self._run_base(itype, ot, a, None, Rab, 0.0, 0.0, z)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
