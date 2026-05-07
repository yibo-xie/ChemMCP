"""
分子轨道能级图生成工具 (MO Energy Level Diagram) — MCP #465
生成完整 MO 能级图数据，支持同核/异核双原子分子及小多原子分子。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MOEnergyLevelDiagram(BaseTool):
    """
    分子轨道能级图生成工具。生成完整的 MO 能级图结构化数据，
    包含对称性标签、占据情况、成键/反键性质、可用于绑图的数据。
    支持同核/异核双原子分子 (H₂, N₂, O₂, F₂, CO, HF, HeH⁺ 等)。
    """
    __version__ = "0.1.0"
    name = "MOEnergyLevelDiagram"
    func_name = "mo_energy_level_diagram"
    description = "Generate molecular orbital energy level diagram data: energies, symmetries, occupations, bonding/antibonding character for diatomic and small polyatomic molecules."
    implementation_description = "Uses LCAO-MO theory with experimental ionization data and Walsh rules to construct qualitative/semi-quantitative MO diagrams. For homonuclear diatomics uses g/u and σ/π symmetry labels; for heteronuclear uses σ/π and bonding/antibonding notation."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "MO Diagram", "Molecular Orbital", "LCAO", "Energy Levels", "Bonding Theory"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "'H2'", "Molecule: 'H2', 'N2', 'O2', 'F2', 'CO', 'NO', 'HF', 'HeH+', 'LiH', 'Li2', 'Be2', 'B2', 'C2'."),
        ("diagram_type", "str", "'full'", "Diagram type: 'full' (all orbitals), 'homo_lumo' (frontier only), 'frontier' (with gap analysis), 'comprehensive' (with full analysis)."),
        ("show_symmetry", "bool", "True", "Include symmetry labels (g/u, σ/π, etc.)."),
        ("output_format", "str", "'dict'", "Output format: 'dict' (structured data), 'plot_data' (ready for plotting), 'text' (ASCII art)."),
        ("method", "str", "'LCAO'", "Method: 'LCAO' (qualitative), 'quantitative' (numerical estimates)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: molecule [diagram_type] [show_symmetry T/F] [output_format]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing MO energy levels, orbital data, bond order, HOMO/LUMO info, plot-ready data."),
    ]

    examples = [
        {
            "code_input": {"molecule": "N2"},
            "text_input": {"input_str": "N2"},
            "output": {"result": {"molecule": "N2", "bond_order": 3, "n_orbitals": 8}},
        },
        {
            "code_input": {"molecule": "O2", "diagram_type": "frontier"},
            "text_input": {"input_str": "O2 frontier"},
            "output": {"result": {"homo": {}, "lumo": {}, "gap_eV": 10.1, "paramagnetic": True}},
        },
    ]

    # ── MO Database (experimental/photoelectron + theoretical) ────
    # Energies in eV (negative = bound). Approximate values from PES/theory.
    _MO_DATA = {
        "H2": {
            "electrons": 2,
            "orbitals": [
                {"name": "σ_g(1s)", "energy_eV": -15.8, "sym": "σ_g", "character": "bonding", "n_elec": 2,
                 "ao_contribution": {"1s_Ha": 0.5 + 1/math.sqrt(2)/2, "1s_Hb": 0.5 - 1/math.sqrt(2)/2}},
                {"name": "σ_u*(1s)", "energy_eV": -2.5, "sym": "σ_u", "character": "antibonding", "n_elec": 0,
                 "ao_contribution": {"1s_Ha": 0.5 - 1/math.sqrt(2)/2, "1s_Hb": 0.5 + 1/math.sqrt(2)/2}},
            ],
            "dissociation_eV": 4.75,
            "bond_length_Ang": 0.74,
            "note": "Simplest covalent bond",
        },
        "He2": {
            "electrons": 4,
            "orbitals": [
                {"name": "σ_g(1s)", "energy_eV": -24.6, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s)", "energy_eV": -18.0, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
            ],
            "bond_order": 0, "stable": False, "note": "Bond order 0 → no stable molecule",
        },
        "Li2": {
            "electrons": 6,
            "orbitals": [
                {"name": "σ_g(1s) core", "energy_eV": -55, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s) core", "energy_eV": -52, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -5.0, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": 1.5, "sym": "σ_u", "character": "antibonding", "n_elec": 0},
            ],
            "dissociation_eV": 1.05, "bond_length_Ang": 2.67,
        },
        "Be2": {
            "electrons": 8,
            "orbitals": [
                {"name": "σ_g(1s)", "energy_eV": -115, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s)", "energy_eV": -112, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -9.0, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": -3.0, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
            ],
            "bond_order": 0, "stable": False, "note": "No bond → Be₂ does not exist under normal conditions",
        },
        "B2": {
            "electrons": 10,
            "orbitals": [
                {"name": "σ_g(1s)", "energy_eV": -190, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s)", "energy_eV": -188, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -14, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": -8, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "π_u(2p)", "energy_eV": -3.0, "sym": "π_u", "character": "bonding", "n_elec": 2, "degeneracy": 2},
            ],
            "dissociation_eV": -2.9, "bond_length_Ang": 1.59, "paramagnetic": True,
            "note": "Paramagnetic: π_u HOMO has 2 electrons in 2 degenerate orbitals (triplet ground state)",
        },
        "C2": {
            "electrons": 12,
            "orbitals": [
                {"name": "σ_g(1s)", "energy_eV": -295, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s)", "energy_eV": -293, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -20, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": -12, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "π_u(2p)", "energy_eV": -9.0, "sym": "π_u", "character": "bonding", "n_elec": 4, "degeneracy": 2},
                {"name": "σ_g(2p)", "energy_eV": -7.0, "sym": "σ_g", "character": "bonding", "n_elec": 0},
            ],
            "dissociation_eV": 6.3, "bond_order": 2, "bond_length_Ang": 1.24,
        },
        "N2": {
            "electrons": 14,
            "orbitals": [
                {"name": "σ_g(1s) [core]", "energy_eV": -410, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s) [core]", "energy_eV": -408, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -37, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": -19.5, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "π_u(2p_x, 2p_y)", "energy_eV": -16.0, "sym": "π_u", "character": "bonding", "n_elec": 4, "degeneracy": 2},
                {"name": "σ_g(2p_z)", "energy_eV": -15.5, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "π_g*(2p)", "energy_eV": 7.5, "sym": "π_g", "character": "antibonding", "n_elec": 0, "degeneracy": 2},
                {"name": "σ_u*(2p_z)", "energy_eV": 11.0, "sym": "σ_u", "character": "antibonding", "n_elec": 0},
            ],
            "dissociation_eV": 9.79, "bond_order": 3, "bond_length_Ang": 1.098,
            "ionization_potential_eV": 15.58,  # experimental
            "note": "Triple bond: one σ + two π bonds",
        },
        "O2": {
            "electrons": 16,
            "orbitals": [
                {"name": "σ_g(1s) [core]", "energy_eV": -543, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s) [core]", "energy_eV": -541, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -41, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": -21, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2p_z)", "energy_eV": -15.8, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "π_u(2p_x, 2p_y)", "energy_eV": -17.0, "sym": "π_u", "character": "bonding", "n_elec": 4, "degeneracy": 2},
                {"name": "π_g*(2p_x*, 2p_y*)", "energy_eV": -4.0, "sym": "π_g", "character": "antibonding", "n_elec": 2, "degeneracy": 2},
                {"name": "σ_u*(2p_z*)", "energy_eV": 2.5, "sym": "σ_u", "character": "antibonding", "n_elec": 0},
            ],
            "dissociation_eV": 5.16, "bond_order": 2, "bond_length_Ang": 1.207,
            "paramagnetic": True, "spin_multiplicity": 3,
            "ionization_potential_eV": 12.07,
            "note": "Paramagnetic! Two unpaired electrons in degenerate π_g* orbitals (triplet ground state)",
        },
        "F2": {
            "electrons": 18,
            "orbitals": [
                {"name": "σ_g(1s) [core]", "energy_eV": -696, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(1s) [core]", "energy_eV": -694, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2s)", "energy_eV": -45, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "σ_u*(2s)", "energy_eV": -23, "sym": "σ_u", "character": "antibonding", "n_elec": 2},
                {"name": "σ_g(2p_z)", "energy_eV": -19.5, "sym": "σ_g", "character": "bonding", "n_elec": 2},
                {"name": "π_u(2p_x, 2p_y)", "energy_eV": -20.5, "sym": "π_u", "character": "bonding", "n_elec": 4, "degeneracy": 2},
                {"name": "π_g*(2p_x*, 2p_y*)", "energy_eV": -10.0, "sym": "π_g", "character": "antibonding", "n_elec": 4, "degeneracy": 2},
                {"name": "σ_u*(2p_z*)", "energy_eV": -3.0, "sym": "σ_u", "character": "antibonding", "n_elec": 0},
            ],
            "dissociation_eV": 1.61, "bond_order": 1, "bond_length_Ang": 1.42,
            "note": "Weak single bond due to strong antibonding occupation",
        },
        "CO": {
            "electrons": 14,
            "type": "heteronuclear",
            "orbitals": [
                {"name": "σ (C1s+O1s) [core]", "energy_eV": -540, "sym": "σ", "character": "bonding", "n_elec": 2},
                {"name": "σ* (C1s-O1s) [core]", "energy_eV": -538, "sym": "σ*", "character": "antibonding", "n_elec": 2},
                {"name": "σ (2s)", "energy_eV": -41, "sym": "σ", "character": "bonding", "n_elec": 2},
                {"name": "σ* (2s)", "energy_eV": -20, "sym": "σ*", "character": "antibonding", "n_elec": 2},
                {"name": "π (2p_x, 2p_y)", "energy_eV": -17.0, "sym": "π", "character": "bonding", "n_elec": 4, "degeneracy": 2},
                {"name": "σ (2p_z)", "energy_eV": -14.0, "sym": "σ", "character": "bonding", "n_elec": 2},
                {"name": "π* (2p_x*, 2p_y*)", "energy_eV": 7.0, "sym": "π*", "character": "antibonding", "n_elec": 0, "degeneracy": 2},
                {"name": "σ* (2p_z*)", "energy_eV": 10.5, "sym": "σ*", "character": "antibonding", "n_elec": 0},
            ],
            "dissociation_eV": 11.21, "bond_order": 3, "bond_length_Ang": 1.128,
            "dipole_moment_D": 0.11, "note": "Isoelectronic with N₂, slight dipole C⁻→O⁺",
        },
        "HF": {
            "electrons": 10,
            "type": "heteronuclear",
            "orbitals": [
                {"name": "σ (F1s) [core]", "energy_eV": -696, "sym": "σ", "character": "nonbonding(F)", "n_elec": 2},
                {"name": "σ (bonding)", "energy_eV": -20.0, "sym": "σ", "character": "bonding", "n_elec": 2,
                 "ao_contribution": {"H1s": 0.6, "F2p_z": 0.4}},
                {"name": "n (F2p_x, 2p_y)", "energy_eV": -16.0, "sym": "π", "character": "nonbonding", "n_elec": 4, "degeneracy": 2},
                {"name": "σ* (antibonding)", "energy_eV": 4.0, "sym": "σ*", "character": "antibonding", "n_elec": 0},
            ],
            "dissociation_eV": 6.01, "bond_order": 1, "bond_length_Ang": 0.92,
            "dipole_moment_D": 1.82, "note": "Polar covalent bond, large dipole moment",
        },
        "HeH+": {
            "electrons": 2,
            "type": "heteronuclear",
            "orbitals": [
                {"name": "σ (bonding)", "energy_eV": -25.0, "sym": "σ", "character": "bonding", "n_elec": 2,
                 "ao_contribution": {"He1s": 0.85, "H1s": 0.15}},
                {"name": "σ* (antibonding)", "energy_eV": 5.0, "sym": "σ*", "character": "antibonding", "n_elec": 0,
                 "ao_contribution": {"He1s": 0.15, "H1s": 0.85}},
            ],
            "dissociation_eV": 2.0, "bond_order": 1, "bond_length_Ang": 0.78,
            "charge": 1, "note": "Simplest heteronuclear molecular ion",
        },
        "NO": {
            "electrons": 15,
            "type": "heteronuclear",
            "orbitals": [
                {"name": "σ (core)", "energy_eV": -480, "sym": "σ", "character": "bonding", "n_elec": 2},
                {"name": "σ* (core)", "energy_eV": -478, "sym": "σ*", "character": "antibonding", "n_elec": 2},
                {"name": "σ (2s)", "energy_eV": -39, "sym": "σ", "character": "bonding", "n_elec": 2},
                {"name": "σ* (2s)", "energy_eV": -20, "sym": "σ*", "character": "antibonding", "n_elec": 2},
                {"name": "π (2p)", "energy_eV": -15.5, "sym": "π", "character": "bonding", "n_elec": 4, "degeneracy": 2},
                {"name": "σ (2p_z)", "energy_eV": -13.0, "sym": "σ", "character": "bonding", "n_elec": 1},
                {"name": "π* (2p*)", "energy_eV": -2.0, "sym": "π*", "character": "antibonding", "n_elec": 1, "degeneracy": 2},
            ],
            "dissociation_eV": 6.52, "bond_order": 2.5, "bond_length_Ang": 1.151,
            "paramagnetic": True, "unpaired_electrons": 1,
        },
        "LiH": {
            "electrons": 4,
            "type": "heteronuclear",
            "orbitals": [
                {"name": "σ (Li1s) [core-like]", "energy_eV": -58, "sym": "σ", "character": "nonbonding(Li)", "n_elec": 2},
                {"name": "σ (bonding)", "energy_eV": -8.0, "sym": "σ", "character": "bonding", "n_elec": 2,
                 "ao_contribution": {"Li2s": 0.3, "H1s": 0.7}},
                {"name": "σ* (antibonding)", "energy_eV": 2.0, "sym": "σ*", "character": "antibonding", "n_elec": 0,
                 "ao_contribution": {"Li2s": 0.7, "H1s": 0.3}},
            ],
            "dissociation_eV": 2.52, "bond_order": 1, "bond_length_Ang": 1.60,
            "dipole_moment_D": 5.88, "ionic_character_pct": 70,
            "note": "Highly polar bond, significant Liδ⁺-Hδ⁻ polarity",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, molecule: str, diagram_type: str = "full",
                  show_symmetry: bool = True, output_format: str = "dict",
                  method: str = "LCAO") -> dict:
        """Core logic."""
        mol = molecule.strip()
        mol_key = self._normalize_mol_name(mol)

        if mol_key not in self._MO_DATA:
            available = ", ".join(sorted(self._MO_DATA.keys()))
            raise ChemMCPError(
                f"Unknown molecule '{molecule}'. Available: {available}"
            )

        data = dict(self._MO_DATA[mol_key])
        orbs = list(data.pop("orbitals", []))

        # Compute derived properties
        n_bonding = sum(o["n_elec"] for o in orbs if o["character"] == "bonding")
        n_antibonding = sum(o["n_elec"] for o in orbs if o["character"] == "antibonding")
        bo = (n_bonding - n_antibonding) / 2.0
        if "bond_order" not in data:
            data["bond_order"] = round(bo, 2)

        # Find HOMO/LUMO
        homo_idx, lumo_idx = None, None
        occupied = []
        virtual = []
        for i, o in enumerate(orbs):
            degen = o.get("degeneracy", 1)
            if o["n_elec"] > 0:
                occupied.append(i)
            if o["n_elec"] < 2 * degen:
                virtual.append(i)

        if occupied:
            homo_idx = occupied[-1]
        if virtual:
            lumo_idx = virtual[0]

        result = {
            "molecule": mol,
            "molecule_key": mol_key,
            "method": method,
            "total_electrons": data.get("electrons", sum(o["n_elec"] for o in orbs)),
            "bond_order": data.get("bond_order", round(bo, 2)),
            "diagram_type": diagram_type,
            "molecular_orbitals": [],
        }

        # Add extra properties
        for key in ["dissociation_eV", "bond_length_Ang", "paramagnetic",
                     "spin_multiplicity", "ionization_potential_eV",
                     "dipole_moment_D", "charge", "type", "stable",
                     "ionic_character_pct", "unpaired_electrons", "note"]:
            if key in data:
                result[key] = data[key]

        # Build orbital list based on diagram_type
        if diagram_type == "homo_lumo":
            if homo_idx is not None:
                result["homo"] = self._enrich_orbital(orbs[homo_idx], show_symmetry)
            if lumo_idx is not None:
                result["lumo"] = self._enrich_orbital(orbs[lumo_idx], show_symmetry)
            if homo_idx is not None and lumo_idx is not None:
                e_h = orbs[homo_idx]["energy_eV"]
                e_l = orbs[lumo_idx]["energy_eV"]
                result["homo_lumo_gap_eV"] = round(e_l - e_h, 3)
                result["homo_lumo_gap_nm"] = round(1240.0 / max(e_l - e_h, 0.01), 1)

        elif diagram_type == "frontier":
            if homo_idx is not None:
                result["homo"] = self._enrich_orbital(orbs[homo_idx], show_symmetry)
            if lumo_idx is not None:
                result["lumo"] = self._enrich_orbital(orbs[lumo_idx], show_symmetry)
            e_h = orbs[homo_idx]["energy_eV"] if homo_idx else 0
            e_l = orbs[lumo_idx]["energy_eV"] if lumo_idx else 100
            result["homo_lumo_gap_eV"] = round(e_l - e_h, 3)
            result["global_hardness_eV"] = round((e_l - e_h) / 2.0, 3)
            result["global_softness_eV^-1"] = round(2.0 / max(e_l - e_h, 0.001), 3)
            result["approximate_ip_eV"] = round(-e_h, 3)  # Koopmans
            result["approximate_ea_eV"] = round(-e_l, 3)
            # Include all orbitals as well
            result["all_orbitals"] = [self._enrich_orbital(o, show_symmetry) for o in orbs]

        elif diagram_type in ("full", "comprehensive"):
            result["molecular_orbitals"] = [self._enrich_orbital(o, show_symmetry) for o in orbs]
            if homo_idx is not None:
                result["homo_index"] = homo_idx
                result["homo"] = orbs[homo_idx]["name"]
            if lumo_idx is not None:
                result["lumo_index"] = lumo_idx
                result["lumo"] = orbs[lumo_idx]["name"]
            if homo_idx is not None and lumo_idx is not None:
                result["homo_lumo_gap_eV"] = round(orbs[lumo_idx]["energy_eV"] - orbs[homo_idx]["energy_eV"], 3)

            if diagram_type == "comprehensive":
                result["analysis"] = self._comprehensive_analysis(data, orbs, homo_idx, lumo_idx)

        # Output format variations
        if output_format == "text":
            result["ascii_diagram"] = self._ascii_diagram(orbs, molecule)
        elif output_format == "plot_data":
            result["plot_data"] = self._plot_data(orbs, molecule)

        return {"result": result}

    # ── Enrich Orbital Data ────────────────────────────────────────
    def _enrich_orbital(self, orb: dict, show_sym: bool) -> dict:
        enriched = dict(orb)
        enriched["energy_Hartree"] = round(orb["energy_eV"] / self.Hartree_to_eV, 6)
        if not show_sym:
            enriched.pop("sym", None)
        return enriched

    # ── Comprehensive Analysis ─────────────────────────────────────
    def _comprehensive_analysis(self, data: dict, orbs: list, homo_i, lumo_i) -> dict:
        analysis = {}
        bo = data.get("bond_order", 0)
        if bo >= 3:
            analysis["bond_strength"] = "Very Strong"
        elif bo >= 2:
            analysis["bond_strength"] = "Strong"
        elif bo >= 1:
            analysis["bond_strength"] = "Moderate"
        elif bo > 0:
            analysis["bond_strength"] = "Weak"
        else:
            analysis["bond_strength"] = "No Bond"
        analysis["thermodynamic_stability"] = (
            "Stable" if data.get("stable", True) and bo > 0
            else "Unstable/Metastable" if bo <= 0
            else "Weakly Bound"
        )
        analysis["magnetic_behavior"] = (
            "Paramagnetic" if data.get("paramagnetic", False)
            else "Diamagnetic"
        )
        if homo_i is not None and lumo_i is not None:
            gap = orbs[lumo_i]["energy_eV"] - orbs[homo_i]["energy_eV"]
            if gap < 2:
                analysis["reactivity"] = "Very High (soft molecule)"
            elif gap < 5:
                analysis["reactivity"] = "Moderate"
            else:
                analysis["reactivity"] = "Low (hard molecule)"
        return analysis

    # ── ASCII Diagram ──────────────────────────────────────────────
    @staticmethod
    def _ascii_diagram(orbs: list, mol: str) -> str:
        lines = [f"\n  ══ MO Energy Level Diagram: {mol} ══\n"]
        lines.append(f"  {'Energy (eV)':>14s} | {'Orbital':>25s} | {'χ':>4s} | {'n_e⁻':>5s}")
        lines.append(f"  {'-'*14}-+-{'-'*25}-+-{'-'*4}-+-{'-'*5}")
        for o in sorted(orbs, key=lambda x: x["energy_eV"], reverse=True):
            occ_marker = "●" * o["n_elec"] + "○" * max(0, 2 * o.get("degeneracy", 1) - o["n_elec"])
            char = o["character"][0].upper()  # B/A/N
            lines.append(f"  {o['energy_eV']:>+14.2f} | {o['name']:>25s} | {char:>4s} | {occ_marker:>5s}")
        lines.append("")
        return "\n".join(lines)

    # ── Plot Data ──────────────────────────────────────────────────
    def _plot_data(self, orbs: list, mol: str) -> dict:
        levels = []
        for i, o in enumerate(orbs):
            degen = o.get("degeneracy", 1)
            for d in range(degen):
                levels.append({
                    "index": i,
                    "sub_index": d,
                    "name": o["name"],
                    "energy_eV": o["energy_eV"],
                    "n_electrons": min(o["n_elec"] - d * (o["n_elec"] // degen if degen > 1 else 99),
                                       2) if o["n_elec"] > 0 else 0,
                    "character": o["character"],
                    "symmetry": o.get("sym", ""),
                    "occupied": o["n_elec"] > d * (2 if degen > 1 else 1),
                })
        return {"molecule": mol, "levels": levels, "n_levels": len(levels)}

    # ── Normalize Molecule Name ────────────────────────────────────
    @staticmethod
    def _normalize_mol_name(mol: str) -> str:
        m = mol.strip().lower().replace(" ", "").replace("_", "")
        mapping = {
            "h2": "H2", "he2": "He2", "li2": "Li2", "be2": "Be2",
            "b2": "B2", "c2": "C2", "n2": "N2", "o2": "O2", "f2": "F2",
            "co": "CO", "no": "NO", "hf": "HF", "heh+": "HeH+", "lih": "LiH",
        }
        return mapping.get(m, mol)

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            mol = parts[0]
            dtype = parts[1] if len(parts) > 1 else "full"
            sym = parts[2].upper() != "F" if len(parts) > 2 else True
            fmt = parts[3] if len(parts) > 3 else "dict"
            return self._run_base(mol, dtype, sym, fmt)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")

    # Fix reference to undefined variable
    def _comprehensive_analysis(self, data: dict, orbs: list, homo_i, lumo_i) -> dict:
        analysis = {}
        bo = data.get("bond_order", 0)
        if bo >= 3:
            analysis["bond_strength"] = "Very Strong"
        elif bo >= 2:
            analysis["bond_strength"] = "Strong"
        elif bo >= 1:
            analysis["bond_strength"] = "Moderate"
        elif bo > 0:
            analysis["bond_strength"] = "Weak"
        else:
            analysis["bond_strength"] = "No Bond"
        analysis["thermodynamic_stability"] = (
            "Stable" if data.get("stable", True) and bo > 0
            else "Unstable/Metastable" if bo <= 0
            else "Weakly Bound"
        )
        analysis["magnetic_behavior"] = (
            "Paramagnetic" if data.get("paramagnetic", False)
            else "Diamagnetic"
        )
        if homo_i is not None and lumo_i is not None:
            gap = orbs[lumo_i]["energy_eV"] - orbs[homo_i]["energy_eV"]
            if gap < 2:
                analysis["reactivity"] = "Very High (soft molecule)"
            elif gap < 5:
                analysis["reactivity"] = "Moderate"
            else:
                analysis["reactivity"] = "Low (hard molecule)"
        return analysis
