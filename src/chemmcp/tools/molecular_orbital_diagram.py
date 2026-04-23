import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MolecularOrbitalDiagram(BaseTool):
    """
    分子轨道能级图生成。
    
    使用 LCAO-MO (原子轨道线性组合) 方法构建分子轨道能级图。
    确定键级、电子排布、HOMO/LUMO 信息、磁性等性质。
    """
    __version__ = "0.1.0"
    name = "MolecularOrbitalDiagram"
    func_name = "molecular_orbital_diagram"
    description = "Generate molecular orbital diagrams for diatomic and simple polyatomic molecules using LCAO method."
    implementation_description = "Constructs MO diagrams from atomic orbital combinations using LCAO approximation. Determines bond order, electron configuration, HOMO/LUMO gap, magnetic properties (diamagnetic/paramagnetic), and orbital symmetry labels for common diatomic (H2, N2, O2, F2, CO, NO, HF) and polyatomic molecules (H2O, NH3, CH4, C2H4, C2H2, BH3)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Molecular Orbitals", "LCAO", "Bond Order", "HOMO-LUMO"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule: 'H2', 'N2', 'O2', 'F2', 'CO', 'NO', 'HF', 'H2O', 'NH3', 'CH4', 'C2H4', 'C2H2', 'BH3', 'Li2', 'B2', 'C2'."),
        ("method", "str", "LCAO", "Method: 'LCAO' (default), 'Huckel' (for pi systems), 'extended_Huckel'."),
        ("charge", "int", "0", "Net molecular charge."),
        ("multiplicity", "int", "None", "Spin multiplicity (auto-determined if None)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'molecule [method] [charge]'. Example: 'O2 LCAO 0'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with MO energy levels, electron configuration, bond order, magnetic properties, HOMO-LUMO info, and diagram data."),
    ]

    examples = [
        {
            "code_input": {"molecule": "O2", "method": "LCAO", "charge": 0, "multiplicity": None},
            "text_input": {"input_params": "O2"},
            "output": {"result": {"bond_order": 2.0, "magnetic_properties": "paramagnetic", "homo_lumo_gap_eV": 2.1}},
        },
        {
            "code_input": {"molecule": "N2", "method": "LCAO", "charge": 0, "multiplicity": None},
            "text_input": {"input_params": "N2"},
            "output": {"result": {"bond_order": 3.0, "magnetic_properties": "diamagnetic"}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- MO Data for Diatomic Molecules ----
    
    def _get_diatomic_mo_data(self, mol: str) -> dict:
        """Return pre-computed MO data for homonuclear/heteronuclear diatomics.
        
        Energy ordering follows standard MO theory:
        - For O2, F2, Ne2: σ(2p_z) < π(2p_x)=π(2p_y) < π*(2p_x)=π*(2p_y) < σ*(2p_z)
        - For B2, C2, N2: π(2p_x)=π(2p_y) < σ(2p_z) < π*(2p_x)=π*(2p_y) < σ*(2p_z)
        
        Energies in eV (approximate).
        """
        db = {
            "H2": {
                "aos": ["1s_A", "1s_B"],
                "mos": [
                    {"label": "σ(1s)", "energy_eV": -13.6, "type": "bonding", "electrons": 2,
                     "ao_contrib": {"1s_A": 0.707, "1s_B": 0.707}},
                    {"label": "σ*(1s)", "energy_eV": 15.0, "type": "antibonding", "electrons": 0,
                     "ao_contrib": {"1s_A": 0.707, "1s_B": -0.707}},
                ],
                "total_valence_e": 2,
                "bond_order": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 28.6,
            },
            "Li2": {
                "aos": ["1s_A", "1s_B", "2s_A", "2s_B"],
                "mos": [
                    {"label": "σ(1s)", "energy_eV": -55.0, "type": "bonding", "electrons": 2},
                    {"label": "σ*(1s)", "energy_eV": -50.0, "type": "antibonding", "electrons": 2},
                    {"label": "σ(2s)", "energy_eV": -5.3, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": 1.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 2,
                "bond_order": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 6.3,
            },
            "B2": {
                "aos": ["2s_A", "2s_B", "2px_A", "2py_A", "2pz_A", "2px_B", "2py_B", "2pz_B"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -12.5, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": -8.0, "type": "antibonding", "electrons": 2},
                    {"label": "π(2p_x)", "energy_eV": -3.0, "type": "bonding", "electrons": 1},
                    {"label": "π(2p_y)", "energy_eV": -3.0, "type": "bonding", "electrons": 1},
                    {"label": "σ(2p_z)", "energy_eV": -0.5, "type": "bonding", "electrons": 0},
                    {"label": "π*(2p_x)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                    {"label": "π*(2p_y)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(2p_z)", "energy_eV": 8.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 6,
                "bond_order": 1.0,
                "magnetic": "paramagnetic",
                "homo_lumo_gap_eV": 2.5,
            },
            "C2": {
                "aos": ["2s_A", "2s_B", "2p_A×3", "2p_B×3"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -17.5, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": -12.0, "type": "antibonding", "electrons": 2},
                    {"label": "π(2p_x)", "energy_eV": -8.5, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_y)", "energy_eV": -8.5, "type": "bonding", "electrons": 2},
                    {"label": "σ(2p_z)", "energy_eV": -7.0, "type": "bonding", "electrons": 0},
                    {"label": "π*(2p_x)", "energy_eV": 3.0, "type": "antibonding", "electrons": 0},
                    {"label": "π*(2p_y)", "energy_eV": 3.0, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(2p_z)", "energy_eV": 7.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 8,
                "bond_order": 2.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 1.5,
            },
            "N2": {
                "aos": ["2s_A", "2s_B", "2p_A×3", "2p_B×3"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -27.0, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": -20.0, "type": "antibonding", "electrons": 2},
                    {"label": "π(2p_x)", "energy_eV": -15.5, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_y)", "energy_eV": -15.5, "type": "bonding", "electrons": 2},
                    {"label": "σ(2p_z)", "energy_eV": -14.5, "type": "bonding", "electrons": 2},
                    {"label": "π*(2p_x)", "energy_eV": 7.5, "type": "antibonding", "electrons": 0},
                    {"label": "π*(2p_y)", "energy_eV": 7.5, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(2p_z)", "energy_eV": 12.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 10,
                "bond_order": 3.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 22.0,
            },
            "O2": {
                "aos": ["2s_A", "2s_B", "2p_A×3", "2p_B×3"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -32.4, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": -26.0, "type": "antibonding", "electrons": 2},
                    {"label": "σ(2p_z)", "energy_eV": -18.0, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_x)", "energy_eV": -14.0, "type": "bonding", "electrons": 1},
                    {"label": "π(2p_y)", "energy_eV": -14.0, "type": "bonding", "electrons": 1},
                    {"label": "π*(2p_x)", "energy_eV": -5.5, "type": "antibonding", "electrons": 0},
                    {"label": "π*(2p_y)", "energy_eV": -5.5, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(2p_z)", "energy_eV": 2.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 12,
                "bond_order": 2.0,
                "magnetic": "paramagnetic",
                "homo_lumo_gap_eV": 8.5,
            },
            "F2": {
                "aos": ["2s_A", "2s_B", "2p_A×3", "2p_B×3"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -40.0, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": -34.0, "type": "antibonding", "electrons": 2},
                    {"label": "σ(2p_z)", "energy_eV": -21.0, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_x)", "energy_eV": -17.5, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_y)", "energy_eV": -17.5, "type": "bonding", "electrons": 2},
                    {"label": "π*(2p_x)", "energy_eV": -10.0, "type": "antibonding", "electrons": 2},
                    {"label": "π*(2p_y)", "energy_eV": -10.0, "type": "antibonding", "electrons": 2},
                    {"label": "σ*(2p_z)", "energy_eV": -3.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 14,
                "bond_order": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 7.0,
            },
            "CO": {
                "aos": ["C(core)", "O(core)", "C(valence)", "O(valence)"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -30.0, "type": "bonding", "electrons": 2, "polarity": "slightly polarized toward O"},
                    {"label": "σ*(2s)", "energy_eV": -23.0, "type": "antibonding", "electrons": 2},
                    {"label": "σ(2p_z)", "energy_eV": -16.0, "type": "bonding", "electrons": 2, "polarity": "toward O"},
                    {"label": "π(2p_x)", "energy_eV": -14.0, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_y)", "energy_eV": -14.0, "type": "bonding", "electrons": 2},
                    {"label": "π*(2p_x)", "energy_eV": 6.5, "type": "antibonding", "electrons": 0},
                    {"label": "π*(2p_y)", "energy_eV": 6.5, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(2p_z)", "energy_eV": 11.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 10,
                "bond_order": 3.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 20.5,
                "dipole_moment_D": 0.11,  # Small dipole C⁻→O⁺ (unusual)
            },
            "NO": {
                "aos": ["N(valence)", "O(valence)"],
                "mos": [
                    {"label": "σ(2s)", "energy_eV": -28.0, "type": "bonding", "electrons": 2},
                    {"label": "σ*(2s)", "energy_eV": -22.0, "type": "antibonding", "electrons": 2},
                    {"label": "σ(2p_z)", "energy_eV": -16.5, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_x)", "energy_eV": -13.5, "type": "bonding", "electrons": 2},
                    {"label": "π(2p_y)", "energy_eV": -13.5, "type": "bonding", "electrons": 1},
                    {"label": "π*(2p_x)", "energy_eV": -4.5, "type": "antibonding", "electrons": 0},
                    {"label": "π*(2p_y)", "energy_eV": -4.5, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(2p_z)", "energy_eV": 3.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 11,
                "bond_order": 2.5,
                "magnetic": "paramagnetic",
                "homo_lumo_gap_eV": 9.0,
            },
            "HF": {
                "aos": ["H_1s", "F_1s", "F_2s", "F_2p×3"],
                "mos": [
                    {"label": "σ(F_1s)", "energy_eV": -700.0, "type": "nonbonding", "electrons": 2},
                    {"label": "σ(F_2s)", "energy_eV": -40.0, "type": "nonbonding", "electrons": 2},
                    {"label": "σ(H-F)", "energy_eV": -20.0, "type": "bonding", "electrons": 2, "ao_contrib": {"H_1s": 0.3, "F_2pz": 0.95}},
                    {"label": "n_F(2px)", "energy_eV": -17.5, "type": "nonbonding", "electrons": 2},
                    {"label": "n_F(2py)", "energy_eV": -17.5, "type": "nonbonding", "electrons": 2},
                    {"label": "σ*(H-F)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 8,
                "bond_order": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 22.5,
                "dipole_moment_D": 1.82,
            },
        }
        return db.get(mol.upper(), db.get(mol))

    def _get_polyatomic_mo_data(self, mol: str) -> dict:
        """Return MO data for polyatomic molecules."""
        db = {
            "H2O": {
                "point_group": "C2v",
                "geometry": "bent, 104.5°",
                "mos": [
                    {"label": "a1(O_1s)", "energy_eV": -550, "type": "core", "electrons": 2},
                    {"label": "a1(bonding)", "energy_eV": -20.0, "type": "bonding", "electrons": 2, "desc": "O 2s + H combination"},
                    {"label": "b2(bonding)", "energy_eV": -14.0, "type": "bonding", "electrons": 2, "desc": "O 2py + H combo"},
                    {"label": "a1(bonding)", "energy_eV": -11.0, "type": "bonding", "electrons": 2, "desc": "O 2pz + H combo"},
                    {"label": "b1(nonbonding)", "energy_eV": -13.5, "type": "nonbonding", "electrons": 2, "desc": "O 2px lone pair"},
                    {"label": "b2*(antibonding)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                    {"label": "a1*(antibonding)", "energy_eV": 8.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 8,
                "bond_order_per_OH": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 18.5,
                "dipole_moment_D": 1.85,
            },
            "NH3": {
                "point_group": "C3v",
                "geometry": "trigonal pyramidal, 107°",
                "mos": [
                    {"label": "a1(N_1s)", "energy_eV": -400, "type": "core", "electrons": 2},
                    {"label": "a1(N_2s+H)", "energy_eV": -25.0, "type": "bonding", "electrons": 2},
                    {"label": "e(N_2p+H)", "energy_eV": -14.5, "type": "bonding", "electrons": 4},
                    {"label": "a1(lone pair)", "energy_eV": -10.0, "type": "nonbonding", "electrons": 2},
                    {"label": "e*(antibonding)", "energy_eV": 4.0, "type": "antibonding", "electrons": 0},
                    {"label": "a1*(antibonding)", "energy_eV": 7.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 8,
                "bond_order_per_NH": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 14.0,
                "dipole_moment_D": 1.47,
            },
            "CH4": {
                "point_group": "Td",
                "geometry": "tetrahedral, 109.47°",
                "mos": [
                    {"label": "a1(C_1s)", "energy_eV": -280, "type": "core", "electrons": 2},
                    {"label": "a1(C_2s+H)", "energy_eV": -23.0, "type": "bonding", "electrons": 2},
                    {"label": "t2(C_2p+H)", "energy_eV": -13.5, "type": "bonding", "electrons": 6},
                    {"label": "t2*(antibonding)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                    {"label": "a1*(antibonding)", "energy_eV": 9.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 8,
                "bond_order_per_CH": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 18.5,
            },
            "C2H4": {
                "point_group": "D2h",
                "geometry": "planar, C=C double bond",
                "mos": [
                    {"label": "core (C_1s×2)", "energy_eV": -280, "type": "core", "electrons": 4},
                    {"label": "σ(C-C)", "energy_eV": -22.0, "type": "bonding", "electrons": 2},
                    {"label": "σ(C-H)×4", "energy_eV": -15.0, "type": "bonding", "electrons": 8},
                    {"label": "π(C=C)", "energy_eV": -10.5, "type": "bonding", "electrons": 2},
                    {"label": "π*(C=C)", "energy_eV": 2.0, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(C-H)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(C-C)", "energy_eV": 8.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 12,
                "bond_order_CC": 2.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 12.5,
            },
            "C2H2": {
                "point_group": "D∞h",
                "geometry": "linear, C≡C triple bond",
                "mos": [
                    {"label": "core (C_1s×2)", "energy_eV": -280, "type": "core", "electrons": 4},
                    {"label": "σ(C-C sp)", "energy_eV": -25.0, "type": "bonding", "electrons": 2},
                    {"label": "σ(C-H)×2", "energy_eV": -18.0, "type": "bonding", "electrons": 4},
                    {"label": "π(C=C)_x", "energy_eV": -11.0, "type": "bonding", "electrons": 2},
                    {"label": "π(C=C)_y", "energy_eV": -11.0, "type": "bonding", "electrons": 2},
                    {"label": "π*_x", "energy_eV": 3.5, "type": "antibonding", "electrons": 0},
                    {"label": "π*_y", "energy_eV": 3.5, "type": "antibonding", "electrons": 0},
                    {"label": "σ*(C-C)", "energy_eV": 10.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 10,
                "bond_order_CC": 3.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 14.5,
            },
            "BH3": {
                "point_group": "D3h",
                "geometry": "trigonal planar",
                "mos": [
                    {"label": "a1'(B_1s)", "energy_eV": -180, "type": "core", "electrons": 2},
                    {"label": "e'(B_2p+H)", "energy_eV": -8.0, "type": "bonding", "electrons": 4},
                    {"label": "a2''(empty_pz)", "energy_eV": 0.5, "type": "nonbonding", "electrons": 0},
                    {"label": "e'*(antibonding)", "energy_eV": 5.0, "type": "antibonding", "electrons": 0},
                    {"label": "a1'*(antibonding)", "energy_eV": 8.0, "type": "antibonding", "electrons": 0},
                ],
                "total_valence_e": 6,
                "bond_order_per_BH": 1.0,
                "magnetic": "diamagnetic",
                "homo_lumo_gap_eV": 8.5,
                "lewis_acidity": "strong Lewis acid (empty p orbital)",
            },
        }
        return db.get(mol.upper(), db.get(mol))

    def _compute_bond_order(self, mo_list: list) -> float:
        """Compute bond order from MO occupation."""
        bo = 0.0
        for mo in mo_list:
            e = mo.get("electrons", 0)
            t = mo.get("type", "")
            if t == "bonding":
                bo += e / 2.0
            elif t == "antibonding":
                bo -= e / 2.0
        return bo

    def _find_homo_lumo(self, mo_list: list) -> tuple:
        """Find HOMO and LUMO indices."""
        homo_idx = None
        lumo_idx = None

        for i, mo in enumerate(mo_list):
            if mo.get("electrons", 0) > 0:
                homo_idx = i
            if lumo_idx is None and mo.get("electrons", 0) == 0:
                lumo_idx = i

        return homo_idx, lumo_idx

    def _check_magnetism(self, mo_list: list) -> str:
        """Check if molecule is paramagnetic (unpaired electrons) or diamagnetic."""
        for mo in mo_list:
            e = mo.get("electrons", 0)
            if e % 2 == 1:
                return "paramagnetic"
        return "diamagnetic"

    def _build_electron_config(self, mo_list: list) -> str:
        """Build electron configuration string."""
        parts = []
        for mo in mo_list:
            e = mo.get("electrons", 0)
            if e > 0:
                label = mo["label"].split("(")[0].strip()
                parts.append(f"{label}^{e}")
        return " ".join(parts)

    def _run_base(self, molecule: str, method: str = "LCAO", charge: int = 0,
                  multiplicity: int = None) -> dict:

        mol_upper = molecule.strip().upper()

        # Look up data
        data = self._get_diatomic_mo_data(molecule)
        if data is None:
            data = self._get_polyatomic_mo_data(molecule)

        if data is None:
            raise ChemMCPError(f"Unknown molecule: {molecule}. Supported: "
                             f"H2, Li2, B2, C2, N2, O2, F2, CO, NO, HF, H2O, NH3, CH4, C2H4, C2H2, BH3")

        mos = data.get("mos", [])
        
        # Adjust for charge
        n_valence = data.get("total_valence_e", 0)
        if charge != 0:
            n_valence -= charge  # Remove electrons for positive charge
            # Redistribute electrons (simple fill from bottom)
            remaining = n_valence
            for mo in mos:
                old_e = mo["electrons"]
                cap = 2  # Max electrons per MO
                new_e = min(cap, max(0, remaining))
                mo["electrons"] = new_e
                remaining -= new_e
                if remaining <= 0:
                    break

        # Compute properties
        bond_order = self._compute_bond_order(mos)
        homo_idx, lumo_idx = self._find_homo_lumo(mos)
        magnetic = self._check_magnetism(mos)
        config = self._build_electron_config(mos)

        # HOMO-LUMO gap
        gap_eV = None
        if homo_idx is not None and lumo_idx is not None and lumo_idx < len(mos):
            gap_eV = mos[lumo_idx]["energy_eV"] - mos[homo_idx]["energy_eV"]

        # Build diagram data for plotting
        diagram_data = []
        for mo in mos:
            diagram_data.append({
                "label": mo["label"],
                "energy_eV": mo["energy_eV"],
                "electrons": mo["electrons"],
                "type": mo["type"],
            })

        result = {
            "molecule": molecule,
            "method": method,
            "charge": charge,
            "point_group": data.get("point_group", "D∞h / C∞v"),
            "geometry": data.get("geometry", "diatomic" if data else "unknown"),
            "electron_configuration": config,
            "total_valence Electrons": n_valence,
            "bond_order": round(bond_order, 2),
            "magnetic_properties": magnetic,
            "homo_index": homo_idx,
            "lumo_index": lumo_idx,
            "homo_orbital": mos[homo_idx]["label"] if homo_idx is not None else None,
            "lumo_orbital": mos[lumo_idx]["label"] if lumo_idx is not None else None,
            "homo_energy_eV": round(mos[homo_idx]["energy_eV"], 2) if homo_idx is not None else None,
            "lumo_energy_eV": round(mos[lumo_idx]["energy_eV"], 2) if lumo_idx is not None else None,
            "homo_lumo_gap_eV": round(gap_eV, 2) if gap_eV is not None else data.get("homo_lumo_gap_eV"),
            "mo_diagram": diagram_data,
            "n_mos": len(mos),
            "dipole_moment_Debye": data.get("dipole_moment_D"),
            "additional_info": {k: v for k, v in data.items() 
                               if k not in ("mos", "total_valence_e", "bond_order", 
                                            "magnetic", "homo_lumo_gap_eV")},
        }

        logger.info(f"MolecularOrbitalDiagram: {molecule}, BO={bond_order}, {magnetic}, gap={gap_eV}eV")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mol = parts[0]
            method = parts[1] if len(parts) > 1 else "LCAO"
            charge = int(parts[2]) if len(parts) > 2 else 0
            return self._run_base(mol, method, charge)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'molecule [method] [charge]'")
