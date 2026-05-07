"""
交换积分计算工具 (Exchange Integral) — MCP #464
计算交换积分 K_ij = (ij|ji)，解释 Pauli 原理效应。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ExchangeIntegral(BaseTool):
    """
    交换积分计算工具。计算非经典交换积分 K_ij = (ij|ji)，
    解释 Pauli 不相容原理的量子力学来源、同自旋交换能贡献，
    以及在 Hartree-Fock 方程中的作用。
    """
    __version__ = "0.1.0"
    name = "ExchangeIntegral"
    func_name = "exchange_integral"
    description = "Calculate exchange integrals K_ij = (ij|ji): Pauli principle effect, same-spin electron exchange energy, and Hartree-Fock Fock matrix contribution."
    implementation_description = "Computes exchange integrals for GTO basis functions. Explains the quantum mechanical origin of exchange as a consequence of antisymmetry (Pauli). Shows how K_ij enters the Fock matrix: F_μν = H_μν + Σ_λσ P_λσ[(μν|λσ) - 0.5(μλ|νσ)]."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Exchange Integral", "Pauli Principle", "Hartree-Fock", "Electronic Structure", "Fock Matrix"]
    required_envs = []

    code_input_sig = [
        ("calculation", "str", "'calculate'", "Type: 'calculate' (compute K), 'explain' (Pauli effect explanation), 'fock_contribution' (Fock matrix term), 'full_analysis' (complete analysis)."),
        ("orbital_i", "str", "'1s'", "Orbital i type: '1s', '2s', '2px', '2py', '2pz'."),
        ("orbital_j", "str", "'1s'", "Orbital j type: same options."),
        ("alpha_i", "float", "0.27", "Exponent α_i for orbital i (Bohr⁻²)."),
        ("alpha_j", "float", "0.27", "Exponent α_j for orbital j (Bohr⁻²)."),
        ("overlap_S", "float", "None", "Pre-computed overlap integral S_ij (if None, computed automatically)."),
        ("coulomb_J", "float", "None", "Pre-computed Coulomb integral J_ij (if None, estimated)."),
        ("R_bohr", "float", "1.4", "Distance between orbital centers in Bohr."),
        ("same_spin", "bool", "True", "Same spin? Exchange only contributes for same-spin electrons."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: calculation orbital_i orbital_j [alpha_i alpha_j R_bohr] [same_spin T/F]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing exchange integral value, Pauli effect analysis, Fock matrix contribution, and physical interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "calculation": "calculate",
                "orbital_i": "1s",
                "orbital_j": "1s",
                "alpha_i": 0.27,
                "alpha_j": 0.27,
                "R_bohr": 1.4,
                "same_spin": True,
            },
            "text_input": {
                "input_str": "calculate 1s 1s 0.27 0.27 1.4",
            },
            "output": {
                "result": {
                    "K_value": "...",
                    "K_eV": "...",
                    "same_spin": True,
                    "note": "Exchange integral for same-spin electrons",
                }
            }
        },
        {
            "code_input": {
                "calculation": "explain",
            },
            "text_input": {
                "input_str": "explain",
            },
            "output": {
                "result": {
                    "explanation": "...",
                    "pauli_principle_origin": "...",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, calculation: str = "calculate",
                  orbital_i: str = "1s", orbital_j: str = "1s",
                  alpha_i: float = 0.27, alpha_j: float = 0.27,
                  overlap_S: float = None, coulomb_J: float = None,
                  R_bohr: float = 1.4, same_spin: bool = True) -> dict:
        """Core logic."""
        calc = calculation.lower().strip()

        if calc == "explain":
            return self._explain_pauli()
        elif calc in ("calculate", "compute"):
            return self._calculate_exchange(
                orbital_i, orbital_j, alpha_i, alpha_j,
                overlap_S, coulomb_J, R_bohr, same_spin)
        elif calc == "fock_contribution":
            return self._fock_contribution(
                orbital_i, orbital_j, alpha_i, alpha_j, R_bohr, same_spin)
        elif calc == "full_analysis":
            result = self._calculate_exchange(
                orbital_i, orbital_j, alpha_i, alpha_j,
                overlap_S, coulomb_J, R_bohr, same_spin)["result"]
            explain = self._explain_pauli()["result"]
            fock = self._fock_contribution(
                orbital_i, orbital_j, alpha_i, alpha_j, R_bohr, same_spin)["result"]
            return {"result": {**result, **{"pauli_explanation": explain}, **{"fock_matrix_analysis": fock}}}
        else:
            raise ChemMCPError(
                f"Unknown calculation '{calculation}'. "
                f"Use: calculate, explain, fock_contribution, full_analysis."
            )

    # ── Calculate Exchange Integral ────────────────────────────────
    def _calculate_exchange(self, oi: str, oj: str, ai: float, aj: float,
                            S: float, J: float, R: float, ss: bool) -> dict:
        """Compute K_ij = (ij|ji)."""

        # Compute overlap if not provided
        if S is None:
            p = ai + aj
            mu = ai * aj / p
            S = math.exp(-mu * R * R)

        # Estimate Coulomb integral if not provided
        if J is None:
            J = self._estimate_coulomb(ai, aj, R)

        # Exchange integral approximation:
        # For GTOs on different centers: K ≈ S² · J_eff
        # This captures the key physics: exchange requires overlap
        # The exact formula depends on the specific orbitals involved
        K = S * S * J * self._angular_factor(oi, oj)

        # If not same spin → no exchange contribution (in restricted HF)
        effective_K = K if ss else 0.0

        return {"result": {
            "exchange_integral_K": round(K, 10),
            "effective_K_same_spin": round(effective_K, 10),
            "K_eV": round(K * self.Hartree_to_eV, 6),
            "overlap_S_ij": round(S, 10),
            "coulomb_J_estimate": round(J, 10),
            "ratio_K_over_J": round(K / J if abs(J) > 1e-15 else 0, 6),
            "same_spin": ss,
            "contributes_to_Fock": ss,
            "orbital_i": oi,
            "orbital_j": oj,
            "distance_R_Bohr": round(R, 6),
            "formula": "K_ij = (ij|ji) = ∫∫ χ_i(r₁)χ_j(r₁)(1/r₁₂)χ_j(r₂)χ_i(r₂)dτ₁dτ₂",
            "approximation_used": "K ≈ S² · J · angular_factor (valid for well-separated GTO centers)",
            "physical_meaning": (
                "Exchange integral arises from the antisymmetry requirement of the total wave function "
                "(Pauli principle). It has NO classical analog — it is a purely quantum mechanical effect."
            ),
            "interpretation": self._interpret_K(K, S, ss),
        }}

    # ── Pauli Principle Explanation ────────────────────────────────
    def _explain_pauli(self) -> dict:
        return {"result": {
            "title": "Exchange Integral & Pauli Principle",
            "quantum_origin": (
                "The exchange integral K_ij originates from the antisymmetry of the fermionic wave function. "
                "For a two-electron system: Ψ(x₁,x₂) = (1/√2)[χ_i(1)χ_j(2) - χ_j(1)χ_i(2)]\n\n"
                "The minus sign (from particle exchange) creates an additional term in the expectation value "
                "of the Hamiltonian that has NO classical counterpart."
            ),
            "mathematical_derivation": (
                "⟨Ψ|Ĥ|Ψ⟩ = h_ii + h_jj + J_ij - K_ij\n\n"
                "where:\n"
                "• J_ij = (ij|ij) — Coulomb (classical electrostatic repulsion)\n"
                "• K_ij = (ij|ji) — Exchange (purely quantum, from antisymmetry)\n\n"
                "Key insight: K_ij > 0 means the exchange LOWERS the energy for parallel spins "
                "(because it's subtracted: E = ... + J - K), favoring aligned spins."
            ),
            "fermi_hole": (
                "The exchange effect creates a 'Fermi hole' around each electron — a region of reduced "
                "probability of finding another same-spin electron. This is NOT a physical repulsion but "
                "a correlation built into the wave function's structure.\n\n"
                "The Fermi hole depth: ρ_x(r₂) ≈ -|χ_i(r₂)|² |χ_j(r₁)|² / ρ(r₁) for r₂ near r₁"
            ),
            "consequences": [
                "Exchange stabilizes parallel-spin configurations (origin of ferromagnetism)",
                "K_ij decays rapidly with orbital overlap (short-range effect)",
                "Opposite-spin electrons have no (direct) exchange interaction in HF theory",
                "Exchange is responsible for the 'exchange-correlation' hole in DFT",
            ],
            "in_hartree_fock": (
                "In the Fock operator: F_μν(α) = H_μν + Σ_ν P_λσ(α)[(μν|λσ) - (μλ|νσ)_same_spin]\n\n"
                "The exchange term -(μλ|νσ) makes the Fock matrix non-linear in the density P, "
                "requiring iterative SCF solution."
            ),
        }}

    # ── Fock Matrix Contribution Analysis ──────────────────────────
    def _fock_contribution(self, oi: str, oj: str, ai: float, aj: float,
                           R: float, ss: bool) -> dict:
        p = ai + aj
        mu = ai * aj / p
        S = math.exp(-mu * R * R)
        J = self._estimate_coulomb(ai, aj, R)
        K = S * S * J * self._angular_factor(oi, oj)

        return {"result": {
            "title": "Fock Matrix Exchange Contribution",
            "fock_operator_form": "F_μν = H_μν(core) + Σ_λσ P_λσ[(μν|λσ) - δ_{σ_λ,σ_σ}(μλ|νσ)]",
            "coulomb_part_G_μν": round(J, 8),
            "exchange_part_minus_K_μν": round(-K if ss else 0, 8),
            "net_two_electron_contribution": round(J - (K if ss else 0), 8),
            "is_nonlinear": True,
            "scf_required": True,
            "reason_for_iteration": (
                "The Fock matrix depends on the density matrix P through both J and K terms. "
                "Since P is constructed from the eigenvectors of F itself, an iterative "
                "(self-consistent) procedure is required."
            ),
            "convergence_criterion": "‖P_new - P‖ < threshold (typically 10⁻⁶ to 10⁻⁸)",
        }}

    # ── Helper: Estimate Coulomb ───────────────────────────────────
    @staticmethod
    def _estimate_coulomb(a: float, b: float, R: float) -> float:
        if R < 0.01:
            return (2.0 / math.pi)**0.5 * math.sqrt((a+b)/2)
        return math.sqrt(a*b) * (1.0/R + 0.5) * math.exp(-(a+b)*R*R/4)

    # ── Helper: Angular Factor ─────────────────────────────────────
    @staticmethod
    def _angular_factor(oi: str, oj: str) -> float:
        """Angular momentum factor for exchange integral."""
        s_orbitals = {"1s", "2s"}
        p_orbitals = {"2px", "2py", "2pz"}

        if oi in s_orbitals and oj in s_orbitals:
            return 1.0
        elif (oi in s_orbitals and oj in p_orbitals) or (oi in p_orbitals and oj in s_orbitals):
            return 0.5
        elif oi in p_orbitals and oj in p_orbitals:
            if oi == oj:
                return 0.3
            else:
                return 0.05  # orthogonal p-orbitals → tiny exchange
        return 0.5

    # ── Interpretation ────────────────────────────────────────────
    @staticmethod
    def _interpret_K(K: float, S: float, ss: bool) -> str:
        if not ss:
            return "No direct exchange contribution (opposite spins in restricted HF)"
        if K > 0.5:
            return f"Large exchange (K={K:.4f}) — strong Fermi hole, significant stabilization of parallel-spin configuration"
        elif K > 0.1:
            return f"Moderate exchange (K={K:.4f}) — noticeable spin-pairing effect"
        elif K > 0.01:
            return f"Small exchange (K={K:.4f}) — weak overlap-dependent exchange"
        else:
            return f"Negligible exchange (K={K:.4f}) — orbitals have minimal overlap"

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            calc = parts[0]
            oi = parts[1] if len(parts) > 1 else "1s"
            oj = parts[2] if len(parts) > 2 else "1s"
            ai = float(parts[3]) if len(parts) > 3 else 0.27
            aj = float(parts[4]) if len(parts) > 4 else 0.27
            R = float(parts[5]) if len(parts) > 5 else 1.4
            ss = parts[6].upper() != "F" if len(parts) > 6 else True
            return self._run_base(calc, oi, oj, ai, aj, None, None, R, ss)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
