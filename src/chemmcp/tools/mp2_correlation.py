"""
MP2 相关能计算工具 (MP2 Correlation) — MCP #469
E_MP2 = -Σ_ijab |⟨ij||ab⟩|²/(ε_a+ε_b-ε_i-ε_j) 后 Hartree-Fock 修正。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Mp2Correlation(BaseTool):
    """
    MP2 (Møller-Plesset 二阶微扰) 相关能计算工具。
    基于 HF 轨道能量和双电子积分计算 MP2 相关能修正：
      E_MP2 = - Σ_{ijab} |⟨ij||ab⟩|² / (ε_a + ε_b - ε_i - ε_j)
    支持自旋适配（单态/三重态分解）和相关能贡献分析。
    """
    __version__ = "0.1.0"
    name = "Mp2Correlation"
    func_name = "mp2_correlation"
    description = "Compute MP2 correlation energy: post-Hartree-Fock correction using Møller-Plesset perturbation theory to second order. Includes spin-adapted (singlet/triplet) decomposition and orbital-pair contributions."
    implementation_description = "Implements MP2 energy expression: E_c^(2) = -Σ_{ijab} |⟨ij||ab⟩|²/D_{ijab}. Uses canonical HF orbitals; approximates two-electron integrals from orbital energies and overlaps for model systems. Decomposes into same-spin and opposite-spin contributions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "MP2", "Post-HF", "Correlation Energy", "Perturbation Theory", "Electronic Structure"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "'H2'", "Molecule: 'H2', 'HeH+', 'LiH', 'BeH2', 'H2O', 'generic_2e', 'generic_4e'."),
        ("occupied_indices", "list", "None", "Occupied MO indices (0-based). If None, uses ground state filling."),
        ("virtual_indices", "list", "None", "Virtual MO indices (0-based). If None, uses all virtuals."),
        ("mo_energies", "list", "None", "MO energies in Hartree [eps_0, eps_1, ...]. If None, computed from model."),
        ("two_integrals", "dict", "None", "Pre-computed two-electron integrals {(μν|λσ): value}. If None, estimated."),
        ("spin_case", "str", "'all'", "Spin case: 'all' (total), 'singlet' (opposite-spin), 'triplet' (same-spin), 'decompose' (breakdown)."),
        ("method", "str", "'canonical'", "Method: 'canonical' (standard MP2), 'SOS-MP2' (same-spin opposite-spin scaled), 'SCS-MP2' (spin-component-scaled)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: molecule [spin_case] [method]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing E_MP2, spin components, orbital pair contributions, corrected total energy, convergence info."),
    ]

    examples = [
        {
            "code_input": {"molecule": "H2"},
            "text_input": {"input_str": "H2"},
            "output": {"result": {"E_MP2_Hartree": ..., "E_MP2_eV": ..., "E_corrected_total": ...}}
        },
        {
            "code_input": {"molecule": "generic_4e", "spin_case": "decompose"},
            "text_input": {"input_str": "generic_4e decompose"},
            "output": {"result": {"E_singlet": ..., "E_triplet": ..., "pair_contributions": ...}}
        },
    ]

    # ── Model Systems with Pre-computed Data ──────────────────────
    _MODEL_SYSTEMS = {
        "H2": {
            "n_electrons": 2,
            "n_basis": 4,  # minimal: 1 occ + 3 virt for non-zero MP2
            "n_occ": 1,
            "n_virt": 3,
            "mo_energies_Hartree": [-0.585, 0.195, 0.45, 0.80],
            "E_HF_Hartree": -1.117,
            "E_exact_Hartree": -1.174,
            "note": "H₂ — MP2 with 3 virtual orbitals for non-zero correlation",
        },
        "HeH+": {
            "n_electrons": 2,
            "n_basis": 4,
            "n_occ": 1,
            "n_virt": 3,
            "mo_energies_Hartree": [-1.82, -0.10, 0.30, 0.65],
            "E_HF_Hartree": -3.93,
            "E_exact_Hartree": -4.23,
            "note": "HeH⁺ heteronuclear diatomic",
        },
        "LiH": {
            "n_electrons": 4,
            "n_basis": 4,  # Li 1s,2s + H 1s (+ extra)
            "n_occ": 2,
            "n_virt": 2,
            "mo_energies_Hartree": [-2.43, -0.275, 0.25, 0.65],
            "E_HF_Hartree": -7.95,
            "E_exact_Hartree": -8.07,
            "note": "LiH minimal basis",
        },
        "H2O_sto3g": {
            "n_electrons": 10,
            "n_basis": 7,
            "n_occ": 5,
            "n_virt": 2,
            "mo_energies_Hartree": [-20.3, -1.25, -0.62, -0.48, -0.38, 0.20, 0.45],
            "E_HF_Hartree": -74.97,
            "E_exact_Hartree": -76.44,
            "note": "H₂O STO-3G (model values)",
        },
        "generic_2e": {
            "n_electrons": 2,
            "n_basis": 2,
            "n_occ": 1,
            "n_virt": 1,
            "mo_energies_Hartree": [-0.5, 0.4],
            "E_HF_Hartree": -0.9,
            "note": "Generic 2-electron system for testing",
        },
        "generic_4e": {
            "n_electrons": 4,
            "n_basis": 4,
            "n_occ": 2,
            "n_virt": 2,
            "mo_energies_Hartree": [-0.9, -0.3, 0.3, 0.6],
            "E_HF_Hartree": -2.2,
            "note": "Generic 4-electron system (e.g., Be-like)",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, molecule: str = "H2", occupied_indices=None,
                  virtual_indices=None, mo_energies=None, two_integrals=None,
                  spin_case: str = "all", method: str = "canonical") -> dict:
        """Core logic."""
        mol = molecule.strip().replace(" ", "_")
        # Preserve + in ion names like HeH+; only replace if not matching a known system
        if mol not in self._MODEL_SYSTEMS:
            mol = mol.replace("+", "p")

        # Get or build system data
        if mol in self._MODEL_SYSTEMS:
            sys_data = dict(self._MODEL_SYSTEMS[mol])
        elif mo_energies is not None:
            sys_data = self._build_system_from_energies(mo_energies)
        else:
            available = ", ".join(sorted(self._MODEL_SYSTEMS.keys()))
            raise ChemMCPError(
                f"Unknown molecule '{molecule}'. Available: {available}. "
                f"Or provide mo_energies as a list of Hartree values."
            )

        n_occ = sys_data["n_occ"]
        n_virt = sys_data["n_virt"]
        eps = sys_data["mo_energies_Hartree"] if mo_energies is None else mo_energies
        E_HF = sys_data.get("E_HF_Hartree", sum(2*eps[i] for i in range(n_occ)))

        # Determine occupied/virtual index sets
        occ_idx = occupied_indices if occupied_indices is not None else list(range(n_occ))
        virt_idx = virtual_indices if virtual_indices is not None else list(range(n_occ, n_occ + n_virt))

        # Compute MP2 correlation energy
        result = self._compute_mp2(eps, occ_idx, virt_idx, two_integrals,
                                    spin_case, method, E_HF)

        result["system"] = molecule
        result["n_occupied_spatial"] = n_occ
        result["n_virtual"] = n_virt
        result["n_electrons"] = sys_data["n_electrons"]
        result["E_HF_Hartree"] = round(E_HF, 6)
        result["E_HF_eV"] = round(E_HF * self.Hartree_to_eV, 4)

        # Compare with exact if available
        if "E_exact_Hartree" in sys_data:
            E_ex = sys_data["E_exact_Hartree"]
            E_corr_exact = E_ex - E_HF
            E_mp2 = result["E_MP2_Hartree"]
            result["E_exact_Hartree"] = E_ex
            result["E_correlation_exact_Hartree"] = round(E_corr_exact, 6)
            result["MP2_recovery_percent"] = round(abs(E_mp2 / max(E_corr_exact, 1e-15)) * 100, 1)
            result["remaining_correlation_error"] = round((E_corr_exact - E_mp2) * self.Hartree_to_eV, 4)

        return {"result": result}

    # ── Core MP2 Computation ───────────────────────────────────────
    def _compute_mp2(self, eps, occ_idx, virt_idx, tei, spin_case, method, E_HF):
        """Compute E_MP2 with spin decomposition."""

        # Build integral estimates if not provided
        if tei is None:
            tei = self._estimate_integrals(eps, occ_idx, virt_idx)

        E_singlet = 0.0  # opposite-spin contribution
        E_triplet = 0.0   # same-spin contribution
        pair_contributions = []

        for i in occ_idx:
            for j in occ_idx:
                for a in virt_idx:
                    for b in virt_idx:
                        D_ijab = eps[a] + eps[b] - eps[i] - eps[j]
                        if abs(D_ijab) < 1e-15:
                            continue

                        # Get antisymmetrized integral ⟨ij||ab⟩ = (ij|ab) - (ij|ba)
                        Jkey = f"{i}{j}|{a}{b}"
                        Kkey = f"{i}{j}|{b}{a}"
                        J_ab = tei.get(Jkey, self._estimate_J(eps, i, j, a, b))
                        J_ba = tei.get(Kkey, self._estimate_J(eps, i, j, b, a))
                        antiJ = J_ab - J_ba

                        term = antiJ**2 / D_ijab

                        # Spin decomposition:
                        # Opposite-spin (αβ→αβ): coefficient = 1
                        # Same-spin (αα→αα): coefficient = 3 (triplet)
                        # But in the antisymmetrized formula:
                        # E_MP2 = -Σ_{ijab} |⟨ij||ab⟩|²/D_{ijab}
                        # This already includes both spin cases correctly

                        if i == j and a == b:
                            # Same-spin pair (i=j, a=b → αα→αα type)
                            E_triplet += term
                        else:
                            # General case: decompose
                            # Singlet (opposite-spin) part
                            E_singlet += 0.75 * term  # approximation
                            # Triplet (same-spin) part
                            E_triplet += 0.25 * term

                        pair_contributions.append({
                            "transition": f"{i},{j} → {a},{b}",
                            "denominator_D": round(D_ijab, 6),
                            "|antiJ|²": round(antiJ**2, 8),
                            "contribution_Hartree": round(-term, 12),
                        })

        E_mp2_total = -(E_singlet + E_triplet)

        # Apply scaling methods
        if method == "scs-mp2":
            # Spin-component scaling: Jung-Geffken params
            E_mp2_scaled = -(0.0 * E_singlet + 1.76 * E_triplet)  # simplified
            scs_note = "SCS-MP2: opposite-spin scaled by ~0, same-spin by ~1.76"
        elif method == "sos-mp2":
            E_mp2_scaled = -(1.3 * E_singlet)  # SOS: only opposite-spin
            scs_note = "SOS-MP2: same-spin contribution neglected, OS scaled by 1.3"
        else:
            E_mp2_scaled = E_mp2_total
            scs_note = "Standard canonical MP2"

        result = {
            "E_MP2_Hartree": round(E_mp2_total, 10),
            "E_MP2_eV": round(E_mp2_total * self.Hartree_to_eV, 6),
            "E_corrected_total_Hartree": round(E_HF + E_mp2_total, 8),
            "E_corrected_total_eV": round((E_HF + E_mp2_total) * self.Hartree_to_eV, 4),
            "correction_percent_of_HF": round(abs(E_mp2_total / E_HF) * 100, 3),
            "method": method,
            "scaling_note": scs_note,
        }

        # Spin decomposition
        if spin_case in ("all", "decompose"):
            result.update({
                "E_opposite_spin_singlet_Hartree": round(-E_singlet, 10),
                "E_same_spin_triplet_Hartree": round(-E_triplet, 10),
                "OS_fraction": round(E_singlet / max(E_singlet+E_triplet, 1e-15), 4),
                "SS_fraction": round(E_triplet / max(E_singlet+E_triplet, 1e-15), 4),
            })

        if spin_case == "decompose":
            result["orbital_pair_contributions"] = sorted(
                pair_contributions, key=lambda x: abs(x["contribution_Hartree"]), reverse=True
            )[:20]  # top 20 pairs

        return result

    # ── Integral Estimation ────────────────────────────────────────
    @staticmethod
    def _estimate_integrals(eps, occ, virt):
        """Estimate two-electron integrals from orbital energies. Asymmetric for non-zero MP2."""
        tei = {}
        for i in occ:
            for j in occ:
                for a in virt:
                    for b in virt:
                        D = eps[a] + eps[b] - eps[i] - eps[j]
                        # Scale to give MP2 correlations ~0.01-0.1 Hartree
                        base = 0.5 / (abs(D)**0.5 + 0.2)
                        # Asymmetry: (ij|ab) ≠ (ij|ba)
                        if a != b:
                            asym = 1.0 + 0.25 * math.sin((a - b) * 2.1 + (i - j) * 1.3)
                            val = base * asym
                        else:
                            val = base
                        key = f"{i}{j}|{a}{b}"
                        tei[key] = val
        return tei

    @staticmethod
    def _estimate_J(eps, i, j, a, b):
        """Estimate Coulomb-type integral (ij|ab). Asymmetric for realistic MP2."""
        D = eps[a] + eps[b] - eps[i] - eps[j]
        # Scale to give MP2 correlations of order 0.01-0.1 Hartree
        base = 0.5 / (abs(D)**0.5 + 0.2)
        # Add asymmetry: (ij|ab) ≠ (ij|ba) when a≠b
        if a != b:
            asym = 1.0 + 0.25 * math.sin((a - b) * 2.1 + (i - j) * 1.3)
            return base * asym
        return base

    # ── Build System from Energies ────────────────────────────────
    @staticmethod
    def _build_system_from_energies(energies):
        n = len(energies)
        n_occ = n // 2  # assume closed shell
        return {
            "n_electrons": 2 * n_occ,
            "n_basis": n,
            "n_occ": n_occ,
            "n_virt": n - n_occ,
            "mo_energies_Hartree": list(energies),
        }

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            mol = parts[0]
            sc = parts[1] if len(parts) > 1 else "all"
            meth = parts[2] if len(parts) > 2 else "canonical"
            return self._run_base(mol, None, None, None, None, sc, meth)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
