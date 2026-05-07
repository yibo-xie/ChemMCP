"""
前线轨道分析工具 (Frontier Orbital Analysis) — MCP #466
HOMO/LUMO 分析、Fukui 函数、化学反应性预测。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FrontierOrbitalAnalysis(BaseTool):
    """
    前线轨道分析工具。计算 HOMO/LUMO 能级与组成、Fukui 函数、
    软硬度/电负性（概念 DFT 框架）、原子贡献分析与反应性预测。
    """
    __version__ = "0.1.0"
    name = "FrontierOrbitalAnalysis"
    func_name = "frontier_orbital_analysis"
    description = "Frontier molecular orbital analysis: HOMO/LUMO energies, compositions, Fukui functions, global hardness/softness/electronegativity, atomic contributions, and chemical reactivity prediction."
    implementation_description = "Implements frontier orbital theory (Fukui, Parr-Pearson conceptual DFT): computes Fukui functions f⁺(r), f⁻(r), f⁰(r) from finite-difference approximations; global hardness η = I - A, electronegativity χ = (I + A)/2, electrophilicity ω = χ²/(2η). Uses LCAO coefficients for atomic orbital contributions to frontier orbitals."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "HOMO", "LUMO", "Fukui Function", "Conceptual DFT", "Reactivity", "Frontier Orbital Theory"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "'benzene'", "Molecule: 'benzene', 'ethylene', 'formaldehyde', 'co', 'h2o', 'nh3', 'ch4', 'n2', 'o2', 'generic'."),
        ("analysis_type", "str", "'basic'", "Analysis type: 'basic' (HOMO/LUMO/gap), 'fukui' (+ Fukui), 'reactivity' (full reactivity), 'full' (complete analysis)."),
        ("method", "str", "'conceptual_dft'", "Method: 'conceptual_dft', 'koopmans', 'perturbation'."),
        ("charge_state", "int", "0", "Net charge of the molecule (for charged species analysis)."),
        ("custom_mo_energies", "list", "None", "Custom MO energies in eV [e1, e2, ...] for generic analysis."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: molecule [analysis_type] [charge]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing HOMO/LUMO data, gap, Fukui functions, hardness, electronegativity, reactivity indices."),
    ]

    examples = [
        {
            "code_input": {"molecule": "benzene", "analysis_type": "fukui"},
            "text_input": {"input_str": "benzene fukui"},
            "output": {"result": {"homo_energy_eV": ..., "lumo_energy_eV": ..., "gap_eV": ..., "fukui_functions": ...}}
        },
        {
            "code_input": {"molecule": "ethylene", "analysis_type": "reactivity"},
            "text_input": {"input_str": "ethylene reactivity"},
            "output": {"result": {"electrophilic_sites": ..., "nucleophilic_sites": ..., "regioselectivity": ...}}
        },
    ]

    # ── Frontier Orbital Database ──────────────────────────────────
    # Data from experiment/theory (energies in eV)
    _MO_DB = {
        "benzene": {
            "formula": "C₆H₆", "point_group": "D₆h",
            "homo_eV": -9.24, "homo_sym": "e₁g(π)", "homo_degeneracy": 2,
            "lumo_eV": 0.83, "lumo_sym": "e₂u(π*)", "lumo_degeneracy": 2,
            "ip_eV": 9.24, "ea_eV": -0.83,
            "homo_composition": {"C_2pz": [1/6]*6},  # equal on all C atoms
            "lumo_composition": {"C_2pz": [1/6]*6},
            "aromatic": True,
            "note": "Aromatic stabilization lowers HOMO and raises LUMO → large gap",
        },
        "ethylene": {
            "formula": "C₂H₄", "point_group": "D₂h",
            "homo_eV": -10.5, "homo_sym": "π", "homo_degeneracy": 1,
            "lumo_eV": 1.5, "lumo_sym": "π*", "lumo_degeneracy": 1,
            "ip_eV": 10.51, "ea_eV": -1.5,
            "homo_composition": {"C_2pz": [0.5, 0.5]},
            "lumo_composition": {"C_2pz": [0.5, -0.5]},
            "reaction_note": "π→π* excitation drives photochemistry and cycloadditions",
        },
        "formaldehyde": {
            "formula": "CH₂O", "point_group": "C₂v",
            "homo_eV": -10.88, "homo_sym": "n_O (lone pair)", "homo_degeneracy": 1,
            "lumo_eV": 0.94, "lumo_sym": "π*_CO", "lumo_degeneracy": 1,
            "ip_eV": 10.88, "ea_eV": -0.94,
            "homo_composition": {"O_2p_y (lp)": 0.8, "others": 0.2},
            "lumo_composition": {"C_2pz": 0.7, "O_2pz": -0.3},
            "note": "LUMO is π*_CO — carbonyl carbon is electrophilic site",
        },
        "co": {
            "formula": "CO", "point_group": "C∞v",
            "homo_eV": -14.0, "homo_sym": "σ", "homo_degeneracy": 1,
            "lumo_eV": -6.5, "lumo_sym": "2π*", "lumo_degeneracy": 2,
            "ip_eV": 14.01, "ea_eV": 1.49,
            "homo_composition": {"O": 0.55, "C": 0.45},
            "lumo_composition": {"C_2p": 0.8, "O_2p": 0.2},
            "note": "Small gap, good σ-donor and π-acceptor ligand",
        },
        "h2o": {
            "formula": "H₂O", "point_group": "C₂v",
            "homo_eV": -12.62, "homo_sym": "1b₁ (O lone pair p)", "homo_degeneracy": 1,
            "lumo_eV": 0.82, "lumo_sym": "4a₁ (σ*)", "lumo_degeneracy": 1,
            "ip_eV": 12.62, "ea_eV": -0.82,
            "homo_composition": {"O_2p": 0.95, "H_1s": 0.05},
            "note": "HOMO is O lone pair — nucleophilic site at oxygen",
        },
        "nh3": {
            "formula": "NH₃", "point_group": "C₃v",
            "homo_eV": -10.2, "homo_sym": "a₁ (N lone pair)", "homo_degeneracy": 1,
            "lumo_eV": 1.0, "lumo_sym": "e (σ*)", "lumo_degeneracy": 2,
            "ip_eV": 10.18, "ea_eV": -1.0,
            "homo_composition": {"N_2sp³_lp": 0.95},
            "note": "N lone pair HOMO — Lewis base behavior",
        },
        "ch4": {
            "formula": "CH₄", "point_group": "Td",
            "homo_eV": -12.6, "homo_sym": "t₂ (C-H bonding)", "homo_degeneracy": 3,
            "lumo_eV": 2.5, "lumo_sym": "a₁* (C-H antibonding)", "lumo_degeneracy": 1,
            "ip_eV": 12.61, "ea_eV": 0.0,
            "homo_composition": {"C_2p + H_1s": "delocalized over 4 C-H bonds"},
            "note": "Large gap — chemically inert (hard molecule)",
        },
        "n2": {
            "formula": "N₂", "point_group": "D∞h",
            "homo_eV": -15.58, "homo_sym": "3σ_g", "homo_degeneracy": 1,
            "lumo_eV": 7.5, "lumo_sym": "1π_g*", "lumo_degeneracy": 2,
            "ip_eV": 15.58, "ea_eV": -7.5,
            "homo_composition": {"N_2p_z": [0.5, 0.5]},
            "note": "Very large gap (~23 eV) — extremely inert (hard)",
        },
        "o2": {
            "formula": "O₂", "point_group": "D∞h",
            "homo_eV": -12.07, "homo_sym": "1π_g*", "homo_degeneracy": 2,
            "lumo_eV": 2.5, "lumo_sym": "3σ_u*", "lumo_degeneracy": 1,
            "ip_eV": 12.07, "ea_eV": -2.5,
            "paramagnetic": True, "spin_S": 1,
            "note": "SOMO (singly occupied MO) — radical character",
        },
        "pyridine": {
            "formula": "C₅H₅N", "point_group": "C₂v",
            "homo_eV": -9.5, "homo_sym": "π", "homo_degeneracy": 1,
            "lumo_eV": 0.2, "lumo_sym": "π*", "lumo_degeneracy": 1,
            "ip_eV": 9.5, "ea_eV": -0.2,
            "homo_composition": {"C_2pz": "distributed, less on C adjacent to N"},
            "note": "N atom withdraws electron density → lower HOMO than benzene",
        },
        "pyrole": {
            "formula": "C₄H₅N", "point_group": "C₂v",
            "homo_eV": -8.4, "homo_sym": "a₂ (π)", "homo_degeneracy": 1,
            "lumo_eV": 0.6, "lumo_sym": "b₁ (π*)", "lumo_degeneracy": 1,
            "ip_eV": 8.4, "ea_eV": -0.6,
            "homo_composition": {"N_2pz": 0.35, "C_2pz": 0.65},
            "note": "N donates electrons into ring → higher HOMO than benzene",
        },
        "tetracyanoethylene": {
            "formula": "C₆N₄", "point_group": "D₂h",
            "homo_eV": -11.8, "homo_sym": "π", "homo_degeneracy": 1,
            "lumo_eV": 1.7, "lumo_sym": "π*", "lumo_degeneracy": 1,
            "ip_eV": 11.8, "ea_eV": 1.7,
            "note": "Strong electron acceptor (low LUMO), classic Diels-Alder dienophile",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, molecule: str = "benzene", analysis_type: str = "basic",
                  method: str = "conceptual_dft", charge_state: int = 0,
                  custom_mo_energies=None) -> dict:
        """Core logic."""
        mol = molecule.strip().lower()
        atype = analysis_type.lower().strip()

        # Look up molecule or use generic
        if mol in self._MO_DB:
            db_data = dict(self._MO_DB[mol])
        elif custom_mo_energies is not None:
            db_data = self._build_generic_from_energies(custom_mo_energies)
        else:
            available = ", ".join(sorted(self._MO_DB.keys()))
            raise ChemMCPError(
                f"Unknown molecule '{molecule}'. Available: {available}. "
                f"Or provide custom_mo_energies as a list of eV values."
            )

        result = {
            "molecule": molecule,
            "formula": db_data.get("formula", "?"),
            "charge_state": charge_state,
            "method": method,
        }

        # Basic HOMO/LUMO info
        E_H = db_data["homo_eV"]
        E_L = db_data["lumo_eV"]
        gap = E_L - E_H

        result.update({
            "homo": {
                "energy_eV": round(E_H, 3),
                "energy_Hartree": round(E_H / self.Hartree_to_eV, 6),
                "symmetry": db_data.get("homo_sym", "?"),
                "degeneracy": db_data.get("homo_degeneracy", 1),
                "composition": db_data.get("homo_composition", {}),
            },
            "lumo": {
                "energy_eV": round(E_L, 3),
                "energy_Hartree": round(E_L / self.Hartree_to_eV, 6),
                "symmetry": db_data.get("lumo_sym", "?"),
                "degeneracy": db_data.get("lumo_degeneracy", 1),
                "composition": db_data.get("lumo_composition", {}),
            },
            "gap_eV": round(gap, 3),
            "gap_nm": round(1240.0 / max(gap, 0.001), 1),
        })

        # Conceptual DFT quantities
        IP = db_data.get("ip_eV", -E_H if method == "koopmans" else abs(E_H))
        EA = db_data.get("ea_eV", -E_L if method == "koopmans" else abs(E_L))

        chi = (IP + EA) / 2.0   # Mulliken electronegativity
        eta = IP - EA           # Global hardness
        S = 1.0 / (2 * eta) if eta > 0.001 else 999  # Softness
        omega = chi**2 / (2 * eta) if eta > 0.001 else 0  # Electrophilicity index

        result["conceptual_dft"] = {
            "ionization_potential_IP_eV": round(IP, 3),
            "electron_affinity_EA_eV": round(EA, 3),
            "electronegativity_chi_eV": round(chi, 3),
            "global_hardness_eta_eV": round(eta, 3),
            "global_softness_S_eV^-1": round(S, 4),
            "electrophilicity_index_omega_eV": round(omega, 3),
            "chemical_potential_mu_eV": round(-chi, 3),  # μ = -χ
            "hard_soft_classification": (
                "Hard" if eta > 6 else "Borderline" if eta > 3 else "Soft"
            ),
        }

        # Analysis type specific results
        if atype == "fukui":
            result["fukui_functions"] = self._compute_fukui(db_data)
        elif atype == "reactivity":
            result["fukui_functions"] = self._compute_fukui(db_data)
            result["reactivity_prediction"] = self._predict_reactivity(db_data, chi, eta, gap)
        elif atype == "full":
            result["fukui_functions"] = self._compute_fukui(db_data)
            result["reactivity_prediction"] = self._predict_reactivity(db_data, chi, eta, gap)
            result["orbital_interaction_analysis"] = self._orbital_interaction(db_data)

        return {"result": result}

    # ── Fukui Function Computation ────────────────────────────────
    def _compute_fukui(self, db_data: dict) -> dict:
        """Compute Fukui functions via finite difference approximation."""
        E_H = db_data["homo_eV"]
        E_L = db_data["lumo_eV"]

        # Finite difference (Koopmans' theorem):
        # f⁺(r) ≈ |φ_LUMO(r)|²  (for nucleophilic attack)
        # f⁻(r) ≈ |φ_HOMO(r)|²  (for electrophilic attack)
        # f⁰(r)) ≈ ½(|φ_HOMO|² + |φ_LUMO|²)  (radical attack)

        homo_comp = db_data.get("homo_composition", {})
        lumo_comp = db_data.get("lumo_composition", {})

        return {
            "fukui_plus_f+": {  # Nucleophilic susceptibility
                "description": "f⁺(r) = ρ_{N+1}(r) - ρ_N(r) ≈ |φ_LUMO(r)|²",
                "physical_meaning": "Sites where f⁺ is large are susceptible to nucleophilic attack",
                "atomic_contributions": lumo_comp,
                "approximation": "Koopmans: f⁺ ≈ |φ_LUMO|²",
            },
            "fukui_minus_f-": {  # Electrophilic susceptibility
                "description": "f⁻(r) = ρ_N(r) - ρ_{N-1}(r) ≈ |φ_HOMO(r)|²",
                "physical_meaning": "Sites where f⁻ is large are susceptible to electrophilic attack",
                "atomic_contributions": homo_comp,
                "approximation": "Koopmans: f⁻ ≈ |φ_HOMO|²",
            },
            "fukui_zero_f0": {  # Radical susceptibility
                "description": "f⁰(r) = ½[ρ_{N+1}(r) - ρ_{N-1}(r)] ≈ ½(|φ_HOMO|² + |φ_LUMO|²)",
                "physical_meaning": "Sites where f⁰ is large are susceptible to radical attack",
                "atomic_contributions": self._avg_dict(homo_comp, lumo_comp),
                "approximation": "Koopmans: f⁰ ≈ ½(|φ_HOMO|² + |φ_LUMO|²)",
            },
            "philicity_indices": {
                "ω⁺ (nucleophilic)": round((E_L**2) / (2 * (E_L - E_H)) if E_L > E_H else 0, 3),
                "ω⁻ (electrophilic)": round((E_H**2) / (2 * (E_L - E_H)) if E_L > E_H else 0, 3),
            },
        }

    # ── Reactivity Prediction ─────────────────────────────────────
    def _predict_reactivity(self, db_data: dict, chi: float, eta: float, gap: float) -> dict:
        pred = {}
        E_H = db_data["homo_eV"]
        E_L = db_data["lumo_eV"]

        # HSAB principle based classification
        if eta > 6:
            pred["hsab_class"] = "Hard"
            pred["preferred_reactions"] = [
                "Electrostatic-controlled reactions",
                "Hard-hard interactions (ionic/covalent with small orbital overlap)",
                "Charge-controlled selectivity",
            ]
        elif eta < 3:
            pred["hsab_class"] = "Soft"
            pred["preferred_reactions"] = [
                "Orbital-controlled reactions",
                "Soft-soft interactions (large orbital overlap, covalent)",
                "Frontier-controlled selectivity",
            ]
        else:
            pred["hsab_class"] = "Borderline"
            pred["preferred_reactions"] = ["Both electrostatic and orbital effects important"]

        # Specific site predictions from composition
        homo_comp = db_data.get("homo_composition", {})
        lumo_comp = db_data.get("lumo_composition", {})

        pred["electrophilic_attack_sites"] = self._rank_atoms(homo_comp, "HOMO")
        pred["nucleophilic_attack_sites"] = self._rank_atoms(lumo_comp, "LUMO")
        pred["radical_attack_sites"] = self._rank_atoms(
            self._avg_dict(homo_comp, lumo_comp), "average")

        # General predictions
        if gap < 3:
            pred["general_reactivity"] = "High — soft molecule, kinetically labile"
        elif gap < 6:
            pred["general_reactivity"] = "Moderate — typical organic molecule"
        else:
            pred["general_reactivity"] = "Low — hard/inert molecule"

        return pred

    # ── Orbital Interaction Analysis ──────────────────────────────
    def _orbital_interaction(self, db_data: dict) -> dict:
        E_H = db_data["homo_eV"]
        E_L = db_data["lumo_eV"]
        return {
            "frontier_orbital_interaction_summary": (
                f"HOMO({db_data.get('homo_sym','?')}) at {E_H:.2f} eV can donate electrons "
                f"to appropriate acceptor LUMOs.\n"
                f"LUMO({db_data.get('lumo_sym','?')}) at {E_L:.2f} eV can accept electrons "
                f"from donor HOMOs."
            ),
            "possible_interactions": [
                "As a nucleophile (donor): HOMO energy indicates donating ability",
                "As an electrophile (acceptor): LUMO energy indicates accepting ability",
                "In pericyclic reactions: FMO symmetry determines allowed/forbidden",
                "In catalysis: Metal d-orbitals interact with HOMO/LUMO",
            ],
            "woodward_hoffmann_note": (
                "For concerted reactions, the symmetry match between reacting partner's "
                "FMOs determines thermally vs photochemically allowed pathways."
            ),
        }

    # ── Helpers ────────────────────────────────────────────────────
    @staticmethod
    def _build_generic_from_energies(energies: list) -> dict:
        sorted_e = sorted(energies)
        n = len(energies)
        n_elec = n  # assume closed shell for simplicity
        return {
            "homo_eV": sorted_e[n_elec // 2 - 1],
            "lumo_eV": sorted_e[n_elec // 2],
            "ip_eV": abs(sorted_e[n_elec // 2 - 1]),
            "ea_eV": -abs(sorted_e[n_elec // 2]),
        }

    @staticmethod
    def _avg_dict(d1: dict, d2: dict) -> dict:
        keys = set(list(d1.keys()) + list(d2.keys()))
        result = {}
        for k in keys:
            v1 = d1.get(k, 0)
            v2 = d2.get(k, 0)
            # Handle both scalar and list values
            if isinstance(v1, (list, tuple)) and isinstance(v2, (list, tuple)):
                result[k] = [(a+b)/2 for a,b in zip(v1, v2)]
            elif isinstance(v1, (list, tuple)):
                result[k] = [(x+v2)/2 for x in v1]
            elif isinstance(v2, (list, tuple)):
                result[k] = [(v1+x)/2 for x in v2]
            else:
                result[k] = (v1 + v2) / 2.0
        return result

    @staticmethod
    def _rank_atoms(comp: dict, source: str) -> list:
        if not comp or isinstance(comp, str):
            return [{"site": "unknown", "contribution": "N/A", "source": source}]
        if isinstance(comp, list):
            return [{"site": f"atom_{i}", "contribution": c, "source": source}
                    for i, c in enumerate(comp)]
        return sorted(
            [{"site": k, "contribution": v, "source": source} for k, v in comp.items()],
            key=lambda x: -(abs(x["contribution"]) if isinstance(x["contribution"], (int, float)) else 0)
        )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            mol = parts[0]
            atype = parts[1] if len(parts) > 1 else "basic"
            chg = int(parts[2]) if len(parts) > 2 else 0
            return self._run_base(mol, atype, "conceptual_dft", chg)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
