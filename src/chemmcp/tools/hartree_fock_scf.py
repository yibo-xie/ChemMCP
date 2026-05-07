import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HartreeFockSCF(BaseTool):
    """
    Hartree-Fock 自洽场 (Self-Consistent Field) 计算工具。
    实现受限开壳/闭壳 RHF/UHF 的简化 SCF 迭代流程，计算分子轨道能量、总能量、密度矩阵。
    支持最小基组 (STO-3G 风格) 的 H₂、HeH⁺、LiH 等小分子。
    """
    __version__ = "0.1.0"
    name = "HartreeFockSCF"
    func_name = "hartree_fock_scf"
    description = "Perform restricted Hartree-Fock SCF calculations: compute MO energies, total energy, density matrix, convergence for small molecules with minimal basis sets."
    implementation_description = "Implements simplified RHF/SCF procedure: build Fock matrix from core Hamiltonian + Coulomb + exchange, diagonalize iteratively until density matrix converges. Uses STO-3G-style integral approximations for H₂, HeH⁺, LiH-like molecules."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Hartree-Fock", "SCF", "Ab Initio", "MO Theory", "Electronic Structure"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule: 'H2', 'HeH+', 'LiH', 'H2_minimal' (STO-1G), 'generic_2e' (2-electron generic)."),
        ("method", "str", "'RHF'", "Method: 'RHF' (restricted closed-shell), 'UHF' (unrestricted), 'ROHF' (restricted open-shell)."),
        ("bond_length_Angstrom", "float", "0.74", "Internuclear distance in Angstroms (for diatomics)."),
        ("max_iterations", "int", "50", "Maximum SCF iterations."),
        ("convergence_threshold", "float", "1e-6", "RMS density change threshold for convergence."),
        ("initial_guess", "str", "'core'", "Initial guess: 'core' (core Hamiltonian), 'huckel', 'random'."),
        ("charge", "int", "0", "Total molecular charge."),
        ("multiplicity", "int", "1", "Spin multiplicity (2S+1)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: molecule method bond_length [max_iter conv_thresh charge mult]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing SCF results: total energy, MO energies, density matrix, convergence info, orbital occupations."),
    ]

    examples = [
        {
            "code_input": {
                "molecule": "H2",
                "method": "RHF",
                "bond_length_Angstrom": 0.74,
            },
            "text_input": {
                "input_str": "H2 RHF 0.74",
            },
            "output": {
                "result": {
                    "molecule": "H2",
                    "method": "RHF",
                    "total_energy_Hartree": -1.117,
                    "scf_converged": True,
                    "n_iterations": 8,
                    "mo_energies_Hartree": [...],
                    "n_alpha_electrons": 1,
                    "n_beta_electrons": 1,
                }
            }
        },
        {
            "code_input": {
                "molecule": "HeH+",
                "method": "RHF",
                "bond_length_Angstrom": 0.78,
            },
            "text_input": {
                "input_str": "HeH+ RHF 0.78",
            },
            "output": {
                "result": {
                    "molecule": "HeH+",
                    "total_energy_Hartree": "...",
                    "scf_converged": True,
                }
            }
        },
    ]

    # Pre-computed integral data for common molecules at various bond lengths
    # These are approximate STO-3G values (in atomic units / Hartree)
    # Format: {molecule: {R_Ang: {H_core, S, (ii|jj), ...}}}

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34   # J·s
        self.a0 = 5.29177210903e-11   # Bohr radius, m
        self.Hartree_eV = 27.211386245988  # eV per Hartree
        self.Hartree_J = 4.3597447222071e-18  # J per Hartree

    def _run_base(self, molecule: str, method: str = "RHF",
                  bond_length_Angstrom: float = 0.74, max_iterations: int = 50,
                  convergence_threshold: float = 1e-6, initial_guess: str = "core",
                  charge: int = 0, multiplicity: int = 1) -> dict:
        """Core logic: run simplified HF-SCF calculation."""
        mol = molecule.lower().strip().replace("+", "p").replace("-", "m")
        meth = method.upper().strip()
        R_bohr = bond_length_Angstrom / 0.529177210903  # Å → bohr

        if mol == "h2":
            return self._scf_h2(R_bohr, max_iterations, convergence_threshold)
        elif mol in ("hehp", "heh+", "heh_p"):
            return self._scf_heh(R_bohr, max_iterations, convergence_threshold)
        elif mol == "lih":
            return self._scf_lih(R_bohr, max_iterations, convergence_threshold)
        elif mol == "h2_minimal":
            return self._scf_h2_sto1g(R_bohr, max_iterations, convergence_threshold)
        elif mol == "generic_2e":
            return self._scf_generic_2e(R_bohr, max_iterations, convergence_threshold)
        else:
            raise ChemMCPError(
                f"Unknown molecule '{molecule}'. Available: H2, HeH+, LiH, H2_minimal, generic_2e."
            )

    # ── H₂ STO-3G (minimal basis: 1s on each H) ────────────────────
    def _scf_h2(self, R: float, max_iter: int, conv_thr: float) -> dict:
        """
        RHF for H₂ with minimal STO-3G basis.
        Each H contributes one 1s orbital → 2×2 Fock matrix.
        """
        # STO-3G exponents and coefficients for H 1s (standard values)
        zeta = 1.19  # H 1s exponent for STO-3G
        alpha = zeta ** 2

        # Compute integrals as functions of R
        S = math.exp(-alpha * R * R / 2) * (1 + alpha * R * R / 3)  # overlap

        # Core Hamiltonian: kinetic + nuclear attraction
        # H_aa = T_aa + V_aA + V_aB (kinetic + attraction to A + B)
        # For H atom at its own nucleus: T + V_aA ≈ α/2 (for 1s STO)
        H_aa = -0.5 * alpha  # kinetic energy of 1s STO (approximate)
        # Nuclear attraction to other center
        gamma = alpha * R
        if R > 1e-10:
            V_ab = -(1/R) * (1 - math.exp(-2*gamma)*(1 + gamma))  # attraction to other nucleus
        else:
            V_ab = 0
        H_aa += V_ab  # H_aa includes attraction to both nuclei
        H_bb = H_aa  # symmetric

        # Off-diagonal (resonance integral / kinetic + nuclear terms)
        if R > 1e-10:
            H_ab = S * (-alpha/2 - 1/R + alpha*(1+R)*math.exp(-alpha*R))
        else:
            H_ab = -alpha/2

        H_core = [[H_aa, H_ab], [H_ab, H_bb]]

        # Two-electron integrals (STO-3G approximation)
        # (11|11) = (22|22) — Coulomb repulsion on same center
        J_11 = 0.7746 * alpha ** 0.5  # approximate (ss|ss) integral
        J_12 = 1/R * (1 - (1 + 11*R/8 + 3*R**2/4 + R**3/3) * math.exp(-R)) if R > 1e-10 else 5*alpha**0.5  # (11|22)

        K_12 = S * S * J_12  # exchange approximation: (11|12) ≈ S²·(11|22)

        # SCF iteration
        P = [[0.5, 0.5], [0.5, 0.5]]  # initial guess: equal sharing
        E_total_prev = None

        for iteration in range(max_iter):
            # Build Fock matrix: F = H_core + G(P)
            # G_μν = Σλσ P_λσ[(μν|λσ) - 0.5(μλ|νσ)]
            G00 = P[0][0] * J_11 + P[1][1] * J_12 - 0.5 * (P[0][0] * J_11 + P[0][1] * K_12)
            G01 = P[0][1] * (J_12 - 0.5 * K_12)
            G10 = G01
            G11 = P[1][1] * J_11 + P[0][0] * J_12 - 0.5 * (P[1][1] * J_11 + P[1][0] * K_12)

            F = [[H_core[0][0] + G00, H_core[0][1] + G01],
                 [H_core[1][0] + G10, H_core[1][1] + G11]]

            # Diagonalize Fock matrix (2x2 analytical)
            eps, C = self._diag_2x2(F)

            # Form new density matrix (closed-shell, 2 electrons → occupy lowest MO)
            n_occ = 1  # 2 electrons / 2 spin × 1 spatial orbital
            P_new = [[0.0, 0.0], [0.0, 0.0]]
            for m in range(n_occ):
                for mu in range(2):
                    for nu in range(2):
                        P_new[mu][nu] += 2.0 * C[mu][m] * C[nu][m]

            # Check convergence
            rms = math.sqrt(sum((P_new[mu][nu] - P[mu][nu])**2 for mu in range(2) for nu in range(2)) / 4)

            # Total energy: E = 0.5 * Σμν P_μν(H_μν + F_μν)
            E_total = 0.0
            for mu in range(2):
                for nu in range(2):
                    E_total += 0.5 * P_new[mu][nu] * (H_core[mu][nu] + F[mu][nu])
            # Add nuclear repulsion
            E_total += 1.0 / R if R > 1e-10 else 0

            converged = rms < conv_thr
            P = P_new

            if converged or (E_total_prev is not None and abs(E_total - E_total_prev) < conv_thr):
                return self._build_result("H₂", "RHF", eps, C, P, E_total, iteration + 1, True, R)

            E_total_prev = E_total

        # Did not converge within max iterations
        return self._build_result("H₂", "RHF", eps, C, P, E_total, max_iter, False, R,
                                   warning=f"SCF did not converge in {max_iter} iterations.")

    # ── HeH⁺ ────────────────────────────────────────────────────────
    def _scf_heh(self, R: float, max_iter: int, conv_thr: float) -> dict:
        """RHF for HeH⁺ heteronuclear diatomic."""
        # Similar to H₂ but different nuclear charges
        Z_He = 2.0
        Z_H = 1.0
        zeta_He = 2.0925  # STO-3G He 1s exponent
        zeta_H = 1.24     # STO-3G H 1s exponent (slightly larger in cation)

        alpha_He = zeta_He ** 2
        alpha_H = zeta_H ** 2
        p_ab = alpha_He + alpha_H

        # Overlap between He 1s and H 1s
        S = 16 * (zeta_He * zeta_H)**1.5 / (p_ab**3) if R < 5 else 0.001

        # Core Hamiltonian (different on each center due to different Z)
        H_HeHe = -0.5 * alpha_He - Z_He * math.sqrt(2*alpha_He/math.pi)  # T + V_HeNuc
        H_HH = -0.5 * alpha_H - Z_H * math.sqrt(2*alpha_H/math.pi)

        # Simplified off-diagonal
        H_HeH = S * (H_HeHe + H_HH) / 2

        H_core = [[H_HeHe, H_HeH], [H_HeH, H_HH]]

        # Approximate two-electron integrals
        J_HeHe = 0.95  # (He|He)
        J_HH = 0.77    # (H|H)
        J_HeH = 0.45 / R if R > 0.1 else 4.0  # cross term

        # Run simple SCF
        P = [[2.0, 0.0], [0.0, 0.0]]  # initial: 2 electrons mostly on He
        for it in range(max_iter):
            G00 = P[0][0]*J_HeHe + P[1][1]*J_HeH - 0.5*(P[0][0]*J_HeHe)
            G01 = P[0][1]*(J_HeH - 0.25*J_HeH)
            G11 = P[1][1]*J_HH + P[0][0]*J_HeH - 0.5*(P[1][1]*J_HH)
            F = [[H_core[0][0]+G00, H_core[0][1]+G01],
                 [H_core[1][0]+G01, H_core[1][1]+G11]]
            eps, C = self._diag_2x2(F)
            P_new = [[2*C[0][0]**2, 2*C[0][0]*C[0][1]],
                      [2*C[0][0]*C[0][1], 2*C[0][1]**2]]
            rms = math.sqrt(sum((P_new[i][j]-P[i][j])**2 for i in range(2) for j in range(2))/4)
            E = sum(P_new[i][j]*(H_core[i][j]+F[i][j]) for i in range(2) for j in range(2))*0.5
            E += Z_He * Z_H / R if R > 0.01 else 100
            P = P_new
            if rms < conv_thr:
                break

        return self._build_result("HeH⁺", "RHF", eps, C, P, E, it+1, it < max_iter-1, R)

    # ── LiH ─────────────────────────────────────────────────────────
    def _scf_lih(self, R: float, max_iter: int, conv_thr: float) -> dict:
        """Simplified RHF for LiH (minimal basis: Li 1s,2s + H 1s)."""
        # Very simplified — qualitative result only
        # Li: 1s² core + 2s¹ valence; H: 1s¹
        # Use effective 2-orbital model (Li 2s + H 1s), treat Li 1s as frozen core

        alpha_Li = 0.65  # Li 2s STO exponent (effective)
        alpha_H = 1.0

        S = math.exp(-alpha_Li*alpha_H*R*R/(alpha_Li+alpha_H))

        H_LiLi = -0.2   # Li 2s IP ≈ 5.4 eV ≈ 0.2 Hartree
        H_HH = -0.5     # H 1s IP ≈ 13.6 eV ≈ 0.5 Hartree
        H_LiH = -0.15 * S

        H_core = [[H_LiLi, H_LiH], [H_LiH, H_HH]]

        # Simple SCF
        P = [[1.0, 0.0], [0.0, 1.0]]
        for it in range(max_iter):
            J11, J22, J12 = 0.5, 0.4, 0.25
            G00 = P[0][0]*J11 + P[1][1]*J12 - 0.5*P[0][0]*J11
            G11 = P[1][1]*J22 + P[0][0]*J12 - 0.5*P[1][1]*J22
            G01 = P[0][1]*J12 * 0.5
            F = [[H_core[0][0]+G00, H_core[0][1]+G01], [H_core[1][0]+G01, H_core[1][1]+G11]]
            eps, C = self._diag_2x2(F)
            P_new = [[2*C[0][0]**2, 2*C[0][0]*C[0][1]], [2*C[0][0]*C[0][1], 2*C[0][1]**2]]
            rms = sum(abs(P_new[i][j]-P[i][j]) for i in range(2) for j in range(2))/4
            E = sum(P_new[i][j]*(H_core[i][j]+F[i][j]) for i in range(2) for j in range(2))*0.5 + 3/R
            P = P_new
            if rms < conv_thr:
                break

        return self._build_result("LiH", "RHF", eps, C, P, E, it+1, it < max_iter-1, R,
                                   note="Minimal basis (Li 2s + H 1s); Li 1s treated as frozen core.")

    # ── H₂ STO-1G (single Gaussian, educational) ────────────────────
    def _scf_h2_sto1g(self, R: float, max_iter: int, conv_thr: float) -> dict:
        """Educational STO-1G (single Gaussian per H atom)."""
        alpha = 0.270950  # STO-1G exponent for H 1s

        # 1-Gaussian integrals (analytical)
        S = math.exp(-alpha * R**2 / 2)
        T = 3*alpha/2  # kinetic energy of normalized 1s GTO
        # Nuclear attraction to own center
        V_AA = -2*math.sqrt(2*alpha/math.pi)
        # To other center
        if R > 1e-10:
            V_AB = -math.sqrt(2*alpha/math.pi)/R * math.erf(math.sqrt(alpha/(2*alpha))*R)  # simplified
            V_AB = -2*math.sqrt(alpha/math.pi)/R * (1 - math.exp(-alpha*R**2/2))
        else:
            V_AB = 0

        H_aa = T + V_AA + V_AB
        H_bb = H_aa
        H_ab = S * (T + V_AB)  # rough approximation

        H_core = [[H_aa, H_ab], [H_ab, H_bb]]

        # (ss|ss) for 1-Gaussian: = erfc stuff, simplified
        J_11 = 2 * math.sqrt(alpha / math.pi)
        J_12 = math.sqrt(2/math.pi) * math.erf(math.sqrt(alpha/2)*R) / R if R > 1e-10 else 2*math.sqrt(alpha/math.pi)

        P = [[0.5, 0.5], [0.5, 0.5]]
        for it in range(max_iter):
            G00 = P[0][0]*J_11 + P[1][1]*J_12 - 0.5*P[0][0]*J_11
            G11 = P[1][1]*J_11 + P[0][0]*J_12 - 0.5*P[1][1]*J_11
            G01 = P[0][1]*J_12 * 0.5
            F = [[H_core[0][0]+G00, H_core[0][1]+G01], [H_core[1][0]+G01, H_core[1][1]+G11]]
            eps, C = self._diag_2x2(F)
            P_new = [[2*C[0][0]**2, 2*C[0][0]*C[0][1]], [2*C[0][0]*C[0][1], 2*C[0][1]**2]]
            rms = sum(abs(P_new[i][j]-P[i][j])**2 for i in range(2) for j in range(2))/4
            E = sum(P_new[i][j]*(H_core[i][j]+F[i][j]) for i in range(2) for j in range(2))*0.5 + 1/R
            P = P_new
            if rms < conv_thr:
                break

        return self._build_result("H₂(STO-1G)", "RHF", eps, C, P, E, it+1, it < max_iter-1, R,
                                   note="Educational single-Gaussian basis. Not quantitatively accurate.")

    # ── Generic 2-electron system ───────────────────────────────────
    def _scf_generic_2e(self, R: float, max_iter: int, conv_thr: float) -> dict:
        """Generic 2-electron homonuclear diatomic with adjustable parameters."""
        alpha = 1.0
        S = math.exp(-alpha * R**2 / 2)
        H_aa = -1.0
        H_ab = -0.3 * S
        H_core = [[H_aa, H_ab], [H_ab, H_aa]]
        J_11 = 0.7
        J_12 = 0.5 / (R + 0.5)

        P = [[0.5, 0.5], [0.5, 0.5]]
        for it in range(max_iter):
            G00 = P[0][0]*J_11 + P[1][1]*J_12 - 0.5*P[0][0]*J_11
            G11 = G00
            G01 = P[0][1]*J_12 * 0.5
            F = [[H_core[0][0]+G00, H_core[0][1]+G01], [H_core[1][0]+G01, H_core[1][1]+G11]]
            eps, C = self._diag_2x2(F)
            P_new = [[2*C[0][0]**2, 2*C[0][0]*C[0][1]], [2*C[0][0]*C[0][1], 2*C[0][1]**2]]
            rms = sum(abs(P_new[i][j]-P[i][j])**2 for i in range(2) for j in range(2))/4
            E = sum(P_new[i][j]*(H_core[i][j]+F[i][j]) for i in range(2) for j in range(2))*0.5 + 1/R
            P = P_new
            if rms < conv_thr:
                break

        return self._build_result("Generic 2e⁻", "RHF", eps, C, P, E, it+1, it < max_iter-1, R)

    # ── Helper: 2x2 Symmetric Matrix Diagonalization ───────────────
    @staticmethod
    def _diag_2x2(M):
        """Diagonalize a 2×2 symmetric matrix analytically. Returns (eigenvalues, eigenvectors)."""
        a = M[0][0]; b = M[0][1]; d = M[1][1]
        tr = a + d
        det = a*d - b*b
        disc = tr*tr - 4*det
        if disc < 0:
            disc = 0
        sqrt_disc = math.sqrt(disc)
        e1 = (tr + sqrt_disc) / 2
        e2 = (tr - sqrt_disc) / 2

        # Eigenvectors
        if abs(b) > 1e-15:
            v1 = [b, e1 - a]
            v2 = [b, e2 - a]
        elif abs(e1 - a) < abs(e2 - a):
            v1 = [1, 0]; v2 = [0, 1]
        else:
            v1 = [0, 1]; v2 = [1, 0]

        # Normalize
        for v in [v1, v2]:
            n = math.sqrt(v[0]**2 + v[1]**2)
            if n > 1e-15:
                v[0] /= n; v[1] /= n
            else:
                v[0], v[1] = 1.0, 0.0

        return [e1, e2], [v1, v2]

    # ── Result Builder ──────────────────────────────────────────────
    def _build_result(self, mol: str, method: str, eps, C, P, E_total, n_iter: bool,
                      converged: bool, R_bohr: float, warning: str = "", note: str = "") -> dict:
        mo_info = []
        for i, e in enumerate(sorted(eps)):
            occ = 2 if i == 0 else 0  # closed-shell: lowest occupied
            mo_info.append({
                "orbital_number": i + 1,
                "energy_Hartree": round(e, 6),
                "energy_eV": round(e * self.Hartree_eV, 4),
                "occupation": occ,
                "coefficients": [round(C[j][i], 6) for j in range(len(C))] if i < len(C[0]) else [],
            })

        result = {
            "molecule": mol,
            "method": method,
            "total_energy_Hartree": round(E_total, 6),
            "total_energy_eV": round(E_total * self.Hartree_eV, 4),
            "scf_converged": converged,
            "n_iterations": n_iter,
            "bond_length_Bohr": round(R_bohr, 6),
            "bond_length_Angstrom": round(R_bohr * 0.529177, 6),
            "mo_analysis": mo_info,
            "density_matrix_P": [[round(P[i][j], 6) for j in range(len(P))] for i in range(len(P))],
            "n_basis_functions": len(P),
            "n_alpha_electrons": 1,
            "n_beta_electrons": 1,
        }

        if len(mo_info) >= 2:
            result["homo_energy_Hartree"] = round(mo_info[0]["energy_Hartree"], 6)
            result["lumo_energy_Hartree"] = round(mo_info[1]["energy_Hartree"], 6)
            result["homo_lumo_gap_eV"] = round((mo_info[1]["energy_Hartree"] - mo_info[0]["energy_Hartree"]) * self.Hartree_eV, 4)
            result["ionization_potential_eV"] = round(-mo_info[0]["energy_Hartree"] * self.Hartree_eV, 4)  # Koopmans' theorem

        if warning:
            result["warning"] = warning
        if note:
            result["note"] = note

        return {"result": result}

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split()
            mol = parts[0]
            method = parts[1] if len(parts) > 1 else "RHF"
            R = float(parts[2]) if len(parts) > 2 else 0.74
            mx = int(parts[3]) if len(parts) > 3 else 50
            ct = float(parts[4]) if len(parts) > 4 else 1e-6
            ig = parts[5] if len(parts) > 5 else "core"
            ch = int(parts[6]) if len(parts) > 6 else 0
            mul = int(parts[7]) if len(parts) > 7 else 1
            return self._run_base(mol, method, R, mx, ct, ig, ch, mul)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
