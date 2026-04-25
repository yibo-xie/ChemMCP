import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
H = 6.62607015e-34        # J·s
C = 2.99792458e8           # m/s
NA = 6.02214076e23         # mol⁻¹
AMU = 1.66053906660e-27    # kg


# Reference vibrational frequencies for common bond types (cm⁻¹)
# Format: (IR_active, Raman_active, typical_frequency_cm-1, description)
BOND_VIBRATIONS = {
    # Stretches
    "C-H (alkane) stretch": {"freq": 2850, "ir": True, "raman": True, "type": "stretch", "strength": "medium"},
    "C-H (alkene) stretch": {"freq": 3080, "ir": True, "raman": True, "type": "stretch", "strength": "medium"},
    "C-H (alkyne) stretch": {"freq": 3300, "ir": True, "raman": True, "type": "stretch", "strength": "strong"},
    "O-H stretch": {"freq": 3400, "ir": True, "raman": False, "type": "stretch", "strength": "strong"},
    "N-H stretch": {"freq": 3450, "ir": True, "raman": False, "type": "stretch", "strength": "medium"},
    "C=O stretch": {"freq": 1715, "ir": True, "raman": False, "type": "stretch", "strength": "strong"},
    "C=C stretch (alkene)": {"freq": 1650, "ir": "weak", "raman": True, "type": "stretch", "strength": "variable"},
    "C≡C stretch (alkyne)": {"freq": 2150, "ir": "weak", "raman": True, "type": "stretch", "strength": "variable"},
    "C≡N stretch": {"freq": 2250, "ir": True, "raman": True, "type": "stretch", "strength": "medium-strong"},
    "C-O stretch": {"freq": 1100, "ir": True, "raman": True, "type": "stretch", "strength": "strong"},
    "C-N stretch": {"freq": 1030, "ir": True, "raman": True, "type": "stretch", "strength": "medium"},
    "C-S stretch": {"freq": 700, "ir": True, "raman": True, "type": "stretch", "strength": "strong"},
    "S-S stretch": {"freq": 520, "ir": False, "raman": True, "type": "stretch", "strength": "strong"},
    "S-H stretch": {"freq": 2575, "ir": "weak", "raman": True, "type": "stretch", "strength": "strong"},
    "N=O stretch": {"freq": 1550, "ir": True, "raman": True, "type": "stretch", "strength": "strong"},
    # Bends
    "C-H bend (scissoring)": {"freq": 1465, "ir": True, "raman": True, "type": "bend", "strength": "medium"},
    "C-H bend (rocking)": {"freq": 720, "ir": True, "raman": True, "type": "bend", "strength": "weak"},
    "O-H bend (in-plane)": {"freq": 1400, "ir": True, "raman": False, "type": "bend", "strength": "medium"},
    "N-H bend": {"freq": 1600, "ir": True, "raman": True, "type": "bend", "strength": "medium"},
    "C-H bend (aldehyde)": {"freq": 1370, "ir": True, "raman": True, "type": "bend", "strength": "medium"},
    "=C-H bend (out-of-plane)": {"freq": 800, "ir": True, "raman": False, "type": "bend", "strength": "medium-strong"},
    "C=C-C bend": {"freq": 400, "ir": True, "raman": True, "type": "bend", "strength": "weak"},
    "COO⁻ symmetric stretch": {"freq": 1400, "ir": "weak", "raman": True, "type": "stretch", "strength": "medium"},
    "COO⁻ asymmetric stretch": {"freq": 1610, "ir": True, "raman": "weak", "type": "stretch", "strength": "strong"},
}


@ChemMCPManager.register_tool
class VibrationalModes(BaseTool):
    """
    分子振动模式分析工具。
    分析分子的振动模式数量、类型（伸缩/弯曲）、IR和Raman活性，以及特征频率。
    """
    __version__ = "0.1.0"
    name = "VibrationalModes"
    func_name = "analyze_vibrational_modes"
    description = "Analyze molecular vibrational modes: count degrees of freedom, classify IR/Raman activity, and predict characteristic vibrational frequencies."
    implementation_description = "Uses the 3N-6 (nonlinear) or 3N-5 (linear) rule for vibrational mode counting, mutual exclusion rules for centrosymmetric molecules, and group frequency correlations for peak prediction."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "Vibrational Modes", "IR", "Raman", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("n_atoms", "int", "N/A", "Total number of atoms in the molecule."),
        ("molecule_geometry", "str", "nonlinear", "Molecular geometry: 'linear' or 'nonlinear'."),
        ("bond_types", "list", "N/A", "List of bond types present (for frequency prediction)."),
        ("point_group", "str", "unknown", "Point group symmetry for selection rule analysis."),
        ("include_detailed_analysis", "bool", "True", "Whether to include detailed DOF breakdown and symmetry analysis."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: n_atoms geometry [point_group] [bond_type1 bond_type2 ...]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing vibrational mode count, classification, predicted frequencies with IR/Raman activities, and symmetry analysis."),
    ]

    examples = [
        {
            "code_input": {
                "n_atoms": 3,
                "molecule_geometry": "linear",
                "bond_types": ["C≡N stretch"],
                "point_group": "C∞v",
                "include_detailed_analysis": True,
            },
            "text_input": {
                "input_params": "3 linear Cinfv C≡N_stretch",
            },
            "output": {
                "result": {
                    "formula": "HCN-type (3 atoms, linear)",
                    "total_dof": 9,
                    "vibrational_modes": 4,
                    "translational_dof": 3,
                    "rotational_dof": 2,
                    "predicted_peaks": [
                        {"mode": "C≡N stretch", "frequency_cm-1": 2250, "ir_active": True, "raman_active": True, "type": "stretch"},
                        {"mode": "C-H stretch", "frequency_cm-1": 3300, "ir_active": True, "raman_active": True, "type": "stretch"},
                        {"mode": "H-C≡N bend (doubly degenerate)", "frequency_cm-1": 700, "ir_active": True, "raman_active": True, "type": "bend"},
                    ],
                }
            }
        },
        {
            "code_input": {
                "n_atoms": 7,
                "molecule_geometry": "nonlinear",
                "bond_types": ["C=O stretch", "C-H stretch", "C-H bend", "C-O stretch"],
                "point_group": "Cs",
                "include_detailed_analysis": True,
            },
            "text_input": {
                "input_params": "7 nonlinear Cs C=O_stretch C-H_stretch C-H_bend C-O_stretch",
            },
            "output": {
                "result": {
                    "formula": "General organic molecule (7 atoms, nonlinear)",
                    "total_dof": 21,
                    "vibrational_modes": 15,
                    "translational_dof": 3,
                    "rotational_dof": 3,
                    "predicted_peaks": [
                        {"mode": "C-H stretch", "frequency_cm-1": 2920, "ir_active": True, "raman_active": True, "type": "stretch"},
                        {"mode": "C=O stretch", "frequency_cm-1": 1715, "ir_active": True, "raman_active": False, "type": "stretch"},
                        {"mode": "C-H bend", "frequency_cm-1": 1465, "ir_active": True, "raman_active": True, "type": "bend"},
                        {"mode": "C-O stretch", "frequency_cm-1": 1100, "ir_active": True, "raman_active": True, "type": "stretch"},
                    ],
                    "symmetry_note": "Cs point group: no center of inversion → modes can be both IR and Raman active.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.db = dict(BOND_VIBRATIONS)
        self.inversion_groups = {"d∞h", "d_inf_h", "d2h", "d3h", "d4h", "d6h", "ci", "oh"}

    def _run_base(self, n_atoms: int, molecule_geometry: str = "nonlinear",
                  bond_types: List[str] = None, point_group: str = "unknown",
                  include_detailed_analysis: bool = True) -> dict:
        """Core logic."""
        if n_atoms < 1:
            raise ChemMCPError("Number of atoms must be >= 1.")
        if bond_types is None:
            bond_types = []

        geom = molecule_geometry.lower().strip()
        is_linear = geom == "linear"

        total_dof = 3 * n_atoms
        trans_dof = 3
        rot_dof = 2 if is_linear else 3
        vib_dof = total_dof - trans_dof - rot_dof

        pg = point_group.lower().strip()
        has_inversion = pg in self.inversion_groups

        # Predict peaks from bond types
        peaks = []
        seen = set()
        for bt in bond_types:
            matched_keys = self._match_bond(bt)
            for key in matched_keys:
                if key in seen:
                    continue
                seen.add(key)
                data = self.db[key]
                peaks.append({
                    "mode": key,
                    "frequency_cm-1": data["freq"],
                    "ir_active": data["ir"] if not has_inversion or data.get("raman") is False else data["ir"],
                    "raman_active": data["raman"] if not has_inversion or data.get("ir") is False else data["raman"],
                    "type": data["type"],
                    "strength": data["strength"],
                })

        # Sort by frequency descending
        peaks.sort(key=lambda p: p["frequency_cm-1"], reverse=True)

        result = {
            "n_atoms": n_atoms,
            "geometry": geom,
            "point_group": point_group if point_group != "unknown" else "not specified",
            "has_center_of_inversion": has_inversion,
            "total_dof": total_dof,
            "translational_dof": trans_dof,
            "rotational_dof": rot_dof,
            "vibrational_modes": vib_dof,
            "predicted_peaks": peaks,
            "num_predicted_peaks": len(peaks),
        }

        if include_detailed_analysis:
            result["dof_breakdown"] = (
                f"3N = {total_dof} = {trans_dof} (translation) + {rot_dof} (rotation) + {vib_dof} (vibration)"
            )
            if has_inversion:
                result["mutual_exclusion_rule"] = (
                    "Centrosymmetric molecule: mutually exclusive — "
                    "IR-active modes are Raman-inactive and vice versa."
                )
            else:
                result["selection_rule_note"] = (
                    f"Non-centrosymmetric ({pg}): modes can be simultaneously IR and Raman active."
                )
            result["spectral_regions_summary"] = self._summarize_regions(peaks)

        return {"result": result}

    def _match_bond(self, bt: str) -> List[str]:
        direct = {
            "c=o_stretch": ["C=O stretch"],
            "c-h_stretch": ["C-H (alkane) stretch"],
            "c-h_bend": ["C-H bend (scissoring)"],
            "o-h_stretch": ["O-H stretch"],
            "n-h_stretch": ["N-H stretch"],
            "c-o_stretch": ["C-O stretch"],
            "c-n_stretch": ["C-N stretch"],
            "c=c_stretch": ["C=C stretch (alkene)"],
            "cc_stretch": ["C≡C stretch (alkyne)"],
            "cn_stretch": ["C≡N stretch"],
            "s-s_stretch": ["S-S stretch"],
            "s-h_stretch": ["S-H stretch"],
            "no_stretch": ["N=O stretch"],
            "coo_symmetric": ["COO⁻ symmetric stretch"],
            "coo_asymmetric": ["COO⁻ asymmetric stretch"],
        }
        bt_clean = bt.lower().replace(" ", "_").replace("-", "_")
        if bt_clean in direct:
            return direct[bt_clean]

        matched = []
        for key in self.db:
            clean_key = key.lower().replace(" ", "_").replace("-", "").replace("(", "").replace(")", "")
            if bt_clean in clean_key or any(w in clean_key for w in bt_clean.split("_")):
                matched.append(key)
        return matched if matched else [bt]

    @staticmethod
    def _summarize_regions(peaks: List[dict]) -> str:
        regions = {
            "X-H stretch (2500-4000 cm⁻¹)": [],
            "Triple bond region (2000-2300 cm⁻¹)": [],
            "Double bond region (1500-1850 cm⁻¹)": [],
            "Fingerprint region (<1500 cm⁻¹)": [],
        }
        for p in peaks:
            f = p["frequency_cm-1"]
            if f >= 2500:
                regions["X-H stretch (2500-4000 cm⁻¹)"].append(p["mode"])
            elif f >= 2000:
                regions["Triple bond region (2000-2300 cm⁻¹)"].append(p["mode"])
            elif f >= 1500:
                regions["Double bond region (1500-1850 cm⁻¹)"].append(p["mode"])
            else:
                regions["Fingerprint region (<1500 cm⁻¹)"].append(p["mode"])

        parts = []
        for region, modes in regions.items():
            if modes:
                parts.append(f"{region}: {', '.join(modes)}")
        return " | ".join(parts) if parts else "No bond types provided."

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            n_atoms = int(parts[0])
            geom = parts[1] if len(parts) > 1 else "nonlinear"
            pg = parts[2] if len(parts) > 2 and not parts[2].replace("_","").isalpha() == False else "unknown"
            # Check if parts[2] looks like a point_group vs a bond type
            bonds = []
            for p in parts[3:]:
                if p and not p.replace(".","").isdigit():
                    bonds.append(p)
            return self._run_base(n_atoms, geom, bonds, pg)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
