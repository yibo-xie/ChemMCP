import logging
import math
from typing import Optional, List, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SlaterDeterminant(BaseTool):
    """
    Slater 行列式构建工具。
    构建反对称多电子波函数（Slater 行列式），计算行列式值、归一化、
    自旋轨道占据、激发组态（单/双激发），以及简单矩阵元计算。
    """
    __version__ = "0.1.0"
    name = "SlaterDeterminant"
    func_name = "slater_determinant_build"
    description = "Construct antisymmetric many-electron wave functions via Slater determinants: build, normalize, evaluate, and compute matrix elements for electronic configurations."
    implementation_description = "Implements Slater determinant construction for N-electron systems: evaluates determinant value, normalization, spin-orbital occupancy, excited configurations (singles/doubles), and one- and two-electron matrix elements using the Slater-Condon rules."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Slater Determinant", "Antisymmetry", "Wave Function", "Electronic Structure", "Pauli Principle"]
    required_envs = []

    code_input_sig = [
        ("calculation", "str", "N/A", "Type: 'build' (construct determinant), 'evaluate' (numerical value), 'normalize', 'excited_config' (generate excitations), 'matrix_element' (Slater-Condon rules)."),
        ("system", "str", "'generic'", "System: 'helium', 'lithium', 'beryllium', 'h2_minimal', 'generic_Ne'."),
        ("n_electrons", "int", "2", "Number of electrons."),
        ("spin_orbitals", "list", "None", "List of occupied spin-orbital indices (0-based) for custom configuration. If None, uses ground state filling."),
        ("excitation_level", "int", "1", "Excitation level: 1 (single), 2 (double), 3 (triple)."),
        ("coordinates", "list", "None", "Evaluation points [[x1,y1,z1], [x2,y2,z2], ...] in Bohr for numerical evaluation."),
        ("basis_type", "str", "'sto-1g'", "Basis type for orbital evaluation: 'sto-1g', 'hydrogenic', 'gaussian', 'planewave'."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: calculation system [n_electrons] [excitation_level]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing Slater determinant info: orbital list, determinant value, normalization, excitation configurations, matrix element values."),
    ]

    examples = [
        {
            "code_input": {
                "calculation": "build",
                "system": "helium",
                "n_electrons": 2,
            },
            "text_input": {
                "input_str": "build helium 2",
            },
            "output": {
                "result": {
                    "system": "Helium",
                    "n_electrons": 2,
                    "occupied_spin_orbitals": ["1sα", "1sβ"],
                    "determinant_form": "|1sα(1) 1sβ(2)| / √2",
                    "antisymmetric": True,
                    "singlet_state": True,
                }
            }
        },
        {
            "code_input": {
                "calculation": "excited_config",
                "system": "beryllium",
                "n_electrons": 4,
                "excitation_level": 1,
            },
            "text_input": {
                "input_str": "excited_config beryllium 4 1",
            },
            "output": {
                "result": {
                    "n_singles_excitations": "...",
                    "excited_configurations": [...],
                }
            }
        },
        {
            "code_input": {
                "calculation": "matrix_element",
                "system": "generic_Ne",
                "n_electrons": 2,
            },
            "text_input": {
                "input_str": "matrix_element generic_Ne 2",
            },
            "output": {
                "result": {
                    "one_electron_elements": [...],
                    "two_electron_coulomb": "...",
                    "two_electron_exchange": "...",
                    "slater_condon_rules_applied": True,
                }
            }
        },
    ]

    # Standard orbital ordering: spatial × spin
    # Spin-orbital index: χ_i(r,σ) = φ_k(r)·α(σ) or β(σ)
    ORBITAL_NAMES = {
        0: "1sα",   1: "1sβ",
        2: "2sα",   3: "2sβ",
        4: "2pxα",  5: "2pxβ",
        6: "2pyα",  7: "2pyβ",
        8: "2pzα",  9: "2pzβ",
        10: "3sα", 11: "3sβ",
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34
        self.a0 = 5.29177210903e-11  # Bohr radius, m

    def _run_base(self, calculation: str, system: str = "generic_Ne",
                  n_electrons: int = 2, spin_orbitals=None,
                  excitation_level: int = 1, coordinates=None,
                  basis_type: str = "sto-1g") -> dict:
        """Core logic."""
        calc = calculation.lower().strip()
        sys = system.lower().strip()

        if calc == "build":
            return self._build_determinant(sys, n_electrons, spin_orbitals)
        elif calc == "evaluate":
            return self._evaluate_determinant(sys, n_electrons, spin_orbitals, coordinates, basis_type)
        elif calc == "normalize":
            return self._normalization(sys, n_electrons)
        elif calc in ("excited_config", "excited_configuration"):
            return self._excited_configs(sys, n_electrons, excitation_level)
        elif calc in ("matrix_element", "slater_condon"):
            return self._matrix_element(sys, n_electrons)
        else:
            raise ChemMCPError(
                f"Unknown calculation '{calculation}'. "
                f"Use: build, evaluate, normalize, excited_config, matrix_element."
            )

    # ── Build Slater Determinant ────────────────────────────────────
    def _build_determinant(self, sys: str, n_elec: int, custom_occ=None) -> dict:
        """Build the Slater determinant configuration."""
        if custom_occ is not None:
            occ = sorted(custom_occ)
        else:
            occ = list(range(n_elec))  # ground state: fill lowest

        orbital_list = [self.ORBITAL_NAMES.get(i, f"χ_{i}") for i in occ]

        # Determine term symbol from electron count
        n_alpha = sum(1 for i in occ if i % 2 == 0)  # even index → α
        n_beta = n_elec - n_alpha
        S_spin = abs(n_alpha - n_beta) / 2.0
        mult = int(2 * S_spin + 1)

        # Construct symbolic representation
        if n_elec <= 6:
            det_str = " |" + "  ".join(orbital_list) + "| "
        else:
            det_str = f" |χ₁ χ₂ ... χ_{n_elec}| "

        result = {
            "system": sys.capitalize() if sys != "generic_ne" else "Generic",
            "n_electrons": n_elec,
            "n_alpha_electrons": n_alpha,
            "n_beta_electrons": n_beta,
            "total_spin_S": round(S_spin, 4),
            "spin_multiplicity": mult,
            "term_symbol": self._guess_term_symbol(n_alpha, n_beta),
            "occupied_spin_orbitals": orbital_list,
            "occupied_indices": occ,
            "determinant_representation": f"Ψ = (1/√{n_elec}!) · det[χ_i(x_j)]{det_str}",
            "antisymmetric_property": "Ψ(...,x_i,...,x_j,...) = -Ψ(...,x_j,...,x_i,...)",
            "pauli_principle_satisfied": len(occ) == len(set(occ)),
        }

        # Add system-specific info
        if sys == "helium":
            result["electronic_configuration"] = "1s²"
            result["ground_state_term"] = "¹S₀"
        elif sys == "lithium":
            result["electronic_configuration"] = "1s² 2s¹"
            result["ground_state_term"] = "²S_{1/2}"
        elif sys == "beryllium":
            result["electronic_configuration"] = "1s² 2s²"
            result["ground_state_term"] = "¹S₀"
        elif sys == "h2_minimal":
            result["electronic_configuration"] = "σg(1s)²"
            result["ground_state_term"] = "¹Σ_g⁺"

        return {"result": result}

    # ── Evaluate Determinant Numerically ───────────────────────────
    def _evaluate_determinant(self, sys: str, n_elec: int, custom_occ, coords, basis: str) -> dict:
        """Evaluate Slater determinant at given coordinates."""
        if coords is None or len(coords) < n_elec:
            raise ChemMCPError(f"Need at least {n_elec} coordinate points for {n_elec} electrons.")

        if custom_occ is not None:
            occ = sorted(custom_occ)
        else:
            occ = list(range(n_elec))

        # Build matrix D[i][j] = χ_{occ[i]}(r_j, σ_j)
        N = n_elec
        D = []
        for i in range(N):  # rows = spin-orbitals
            row = []
            orb_idx = occ[i]
            spatial_idx = orb_idx // 2  # which spatial orbital
            spin = "alpha" if orb_idx % 2 == 0 else "beta"

            for j in range(N):  # columns = electrons
                r = coords[j]
                x, y, z = r[0], r[1], r[2]

                # Evaluate spatial part
                phi_val = self._eval_spatial(spatial_idx, x, y, z, basis)

                # Spin part (simplified: assume all alpha for eval, or use provided spin)
                sigma_j = j % 2  # alternate spin assignment
                spin_val = 1.0 if (sigma_j == 0 and spin == "alpha") or (sigma_j == 1 and spin == "beta") else 0.0

                row.append(phi_val * spin_val)
            D.append(row)

        # Compute determinant
        det_val = self._det(D)
        norm_factor = 1.0 / math.sqrt(math.factorial(N))

        return {"result": {
            "system": sys,
            "n_electrons": N,
            "determinant_value": round(det_val * norm_factor, 8),
            "raw_determinant": round(det_val, 8),
            "normalization_factor": round(norm_factor, 8),
            "matrix_D": [[round(v, 6) for v in row] for row in D],
            "evaluation_points_Bohr": coords,
            "basis": basis,
        }}

    # ── Normalization Check ────────────────────────────────────────
    def _normalization(self, sys: str, n_elec: int) -> dict:
        """Check and explain normalization of Slater determinant."""
        N = n_elec
        norm = 1.0 / math.sqrt(math.factorial(N))

        # For orthonormal spin-orbitals ⟨χ_i|χ_j⟩ = δ_ij:
        # ∫ |Ψ|² dτ = 1/N! Σ_p δ_p ... = 1
        ortho_norm_proof = (
            "For orthonormal spin-orbitals:\n"
            f"  ∫|Ψ|² dτ₁...dτ_{N} = (1/{N}!) Σ_P (-1)^P Π_i ⟨χ_i|χ_P(i)⟩\n"
            f"                       = (1/{math.factorial(N)}) · {math.factorial(N)} · 1 = 1  ✓"
        )

        return {"result": {
            "system": sys,
            "n_electrons": N,
            "normalized": True,
            "normalization_constant": round(norm, 10),
            "factorial_denominator": math.factorial(N),
            "proof": ortho_norm_proof,
            "condition": "Spin-orbitals must be orthonormal for exact normalization.",
        }}

    # ── Excited Configurations ─────────────────────────────────────
    def _excited_configs(self, sys: str, n_elec: int, exc_level: int) -> dict:
        """Generate excited configurations (singles/doubles/triples)."""
        occ = list(range(n_elec))  # ground state occupied
        virt = list(range(n_elec, min(n_elec + 10, len(self.ORBITAL_NAMES))))  # virtual orbitals

        excitations = []

        if exc_level == 1:
            # Singles: one electron moved from occ → virt
            for i in occ:
                for a in virt:
                    excitations.append({
                        "type": "single_excitation",
                        "notation": f"|iⁱᵃ⟩",
                        "hole_index": i,
                        "hole_orbital": self.ORBITAL_NAMES.get(i, f"χ_{i}"),
                        "particle_index": a,
                        "particle_orbital": self.ORBITAL_NAMES.get(a, f"χ_{a}"),
                        "excitation_energy_approx": f"ε_{a} - ε_{i}",
                    })
        elif exc_level == 2:
            # Doubles: two electrons moved
            for ii in range(len(occ)):
                for jj in range(ii+1, len(occ)):
                    for aa in range(len(virt)):
                        for bb in range(aa+1, len(virt)):
                            excitations.append({
                                "type": "double_excitation",
                                "notation": f"|ijⁱʲᵃᵇ⟩",
                                "holes": [occ[ii], occ[jj]],
                                "particles": [virt[aa], virt[bb]],
                                "hole_orbitals": [self.ORBITAL_NAMES.get(occ[ii]), self.ORBITAL_NAMES.get(occ[jj])],
                                "particle_orbitals": [self.ORBITAL_NAMES.get(virt[aa]), self.ORBITAL_NAMES.get(virt[bb])],
                            })
                            # Also same-spin doubles
                            if len(virt) >= 2:
                                pass  # already covered by loop
        elif exc_level == 3:
            # Triples
            for i in occ[:min(3, len(occ))]:
                for a in virt[:min(3, len(virt))]:
                    excitations.append({
                        "type": "triple_excitation",
                        "notation": f"|...^{a}...⟩",
                        "description": f"Triple excitation involving hole {i} → particle {a}",
                    })

        total_singles = len(occ) * len(virt) if exc_level == 1 else 0
        total_doubles = (len(occ)*(len(occ)-1)//2) * (len(virt)*(len(virt)-1)//2) if exc_level == 2 else 0

        return {"result": {
            "system": sys,
            "n_electrons": n_elec,
            "n_occupied": len(occ),
            "n_virtual": len(virt),
            "excitation_level": exc_level,
            "n_excitations_generated": len(excitations),
            "total_possible_singles": total_singles,
            "total_possible_doubles": total_doubles,
            "excitations": excitations[:20],  # limit output size
            "truncated": len(excitations) > 20,
        }}

    # ── Matrix Elements (Slater-Condon Rules) ───────────────────────
    def _matrix_element(self, sys: str, n_elec: int) -> dict:
        """Compute one- and two-electron matrix elements using Slater-Condon rules."""
        occ = list(range(n_elec))

        # One-electron operator: ⟨Φ|ĥ|Φ⟩ = Σ_i h_ii
        one_el_elements = []
        for i in occ:
            orbital_name = self.ORBITAL_NAMES.get(i, f"χ_{i}")
            one_el_elements.append({
                "orbital": orbital_name,
                "index": i,
                "integral": f"h_{i}{i}",
                "contribution": "⟨i|h|i⟩",
            })

        # Two-electron Coulomb: J_ij = (ii|jj)
        coulomb_integrals = []
        exchange_integrals = []
        for i_idx in range(len(occ)):
            for j_idx in range(i_idx, len(occ)):
                i = occ[i_idx]
                j = occ[j_idx]
                orb_i = self.ORBITAL_NAMES.get(i, f"χ_{i}")
                orb_j = self.ORBITAL_NAMES.get(j, f"χ_{j}")

                if i == j:
                    # Self-interaction cancels in HF but present in general
                    coulomb_integrals.append({
                        "pair": f"{orb_i}-{orb_j}",
                        "J_integral": f"({i}{i}|{j}{j})",
                        "type": "coulomb_self",
                    })
                else:
                    coulomb_integrals.append({
                        "pair": f"{orb_i}-{orb_j}",
                        "J_integral": f"({i}{i}|{j}{j})",
                        "type": "coulomb",
                    })

                    # Exchange integral K_ij = (ij|ji) — only for same-spin pairs
                    same_spin = (i % 2) == (j % 2)
                    if same_spin:
                        exchange_integrals.append({
                            "pair": f"{orb_i}-{orb_j}",
                            "K_integral": f"({i}{j}|{j}{i})",
                            "type": "exchange",
                            "same_spin": True,
                        })

        # Total energy expression
        E_expr = "E = Σ_i h_ii + ½ Σ_ij [(ii|jj) - (ij|ji)]_same_spin"

        return {"result": {
            "system": sys,
            "n_electrons": n_elec,
            "method": "Slater-Condon Rules",
            "one_electron_sum": one_el_elements,
            "two_electron_coulomb_J": coulomb_integrals,
            "two_electron_exchange_K": exchange_integrals,
            "total_energy_expression": E_expr,
            "slater_condon_rule_1": "⟨Φ|ĥ|Φ⟩ = Σ_i h_ii  (sum over occupied spin-orbitals)",
            "slater_condon_rule_2a": "⟨Φ|ĝ|Φ⟩ = Σ_{i≤j} [(ii|jj) - (ij|ji)δ_{σ_i,σ_j}]",
            "slater_condon_rule_2b": "Singles: ⟨Φ^i_a|ĥ|Φ⟩ = h_ai  (if spins match)",
            "slater_condon_rule_2c": "Doubles: ⟨Φ^ij_ab|ĝ|Φ⟩ = (ab|ij) - (aj|ib)  (CIS-like)",
            "n_coulomb_terms": len(coulomb_integrals),
            "n_exchange_terms": len(exchange_integrals),
        }}

    # ── Spatial Orbital Evaluation Helpers ──────────────────────────
    @staticmethod
    def _eval_spatial(spatial_idx: int, x: float, y: float, z: float, basis: str) -> float:
        """Evaluate spatial orbital at position (x,y,z) in Bohr."""
        r = math.sqrt(x*x + y*y + z*z)

        if basis == "sto-1g":
            # Single Gaussian: φ(r) = (2α/π)^{3/4} exp(-αr²)
            exponents = {0: 0.270950, 1: 0.270950, 2: 0.13009, 3: 0.13009,
                        4: 0.02730, 5: 0.02730, 6: 0.02730, 7: 0.02730}
            alpha = exponents.get(spatial_idx, 0.28)
            norm = (2*alpha/math.pi)**0.75
            val = norm * math.exp(-alpha * r*r)

            # Add angular dependence for p orbitals
            if spatial_idx in (4, 5):  # px
                val *= x if r > 1e-10 else 1.0
            elif spatial_idx in (6, 7):  # py
                val *= y if r > 1e-10 else 1.0
            elif spatial_idx in (8, 9):  # pz
                val *= z if r > 1e-10 else 1.0

            return val

        elif basis == "hydrogenic":
            # Hydrogen-like orbitals
            Z_eff = {0: 2.0, 1: 2.0, 2: 1.0, 3: 1.0}.get(spatial_idx, 1.0)
            n_dict = {0: 1, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
            n = n_dict.get(spatial_idx, 1)
            l_dict = {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1}
            l = l_dict.get(spatial_idx, 0)
            rho = Z_eff * r / n
            val = math.exp(-rho)  # very simplified radial part
            if l == 1:
                if spatial_idx in (4, 5): val *= x/r if r > 1e-10 else 0
                elif spatial_idx in (6, 7): val *= y/r if r > 1e-10 else 0
                elif spatial_idx in (8, 9): val *= z/r if r > 1e-10 else 0
            return val

        else:  # gaussian default
            alpha = 0.28
            return (2*alpha/math.pi)**0.75 * math.exp(-alpha * r*r)

    @staticmethod
    def _det(M):
        """Compute determinant of square matrix (recursive cofactor expansion)."""
        n = len(M)
        if n == 1:
            return M[0][0]
        if n == 2:
            return M[0][0]*M[1][1] - M[0][1]*M[1][0]
        if n == 3:
            return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                  - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                  + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
        # General case — LU decomposition would be better but keep it simple
        det = 0.0
        for j in range(n):
            minor = [[M[row][col] for col in range(n) if col != j] for row in range(1, n)]
            det += ((-1)**j) * M[0][j] * SlaterDeterminant._det(minor)
        return det

    @staticmethod
    def _guess_term_symbol(n_alpha: int, n_beta: int) -> str:
        """Guess atomic term symbol from occupation."""
        S = abs(n_alpha - n_beta) / 2.0
        mult = int(2*S + 1)
        L_map = {0: "S", 1: "P", 2: "D", 3: "F"}
        # Very rough guess based on subshell filling
        if n_alpha + n_beta <= 2:
            L_sym = "S"
        elif n_alpha + n_beta <= 4:
            L_sym = "P" if (n_alpha + n_beta) % 2 == 1 else "S"
        else:
            L_sym = "S"
        return f"^{mult}{L_sym}"

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        try:
            parts = input_str.strip().split()
            calc = parts[0]
            sys = parts[1] if len(parts) > 1 else "generic_Ne"
            ne = int(parts[2]) if len(parts) > 2 else 2
            exlvl = int(parts[3]) if len(parts) > 3 else 1
            return self._run_base(calc, sys, ne, None, exlvl)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
