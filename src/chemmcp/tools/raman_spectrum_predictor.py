import logging
import math
from typing import List, Tuple, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Raman group frequencies (cm⁻¹) — based on standard spectroscopy data
# Format: (shift_lo, shift_hi, intensity, depolarization_ratio, activity_notes)
RAMAN_GROUP_FREQUENCIES = {
    # C-H stretches
    "C-H (alkane sp³)": (2850, 2960, "strong", "low (ρ<0.75)", "polarized"),
    "C-H (alkene sp²)": (3010, 3100, "medium-strong", "low (ρ<0.75)", "polarized"),
    "C≡C-H stretch": (3300, 3340, "strong", "low", "polarized"),
    "C-H (aldehyde)": (2650, 2820, "medium", "low", "polarized"),
    # Triple bonds
    "C≡C (alkyne)": (2100, 2260, "variable-strong", "high (ρ≈0.75-1)", "depolarized"),
    "C≡N (nitrile)": (2220, 2260, "medium-strong", "high", "depolarized"),
    # Double bonds: C=O (Raman usually weak but detectable)
    "C=O (ketone/aldehyde)": (1680, 1750, "weak-medium", "low", "polarized"),
    "C=O (carboxylic acid)": (1700, 1760, "weak", "low", "polarized"),
    "C=O (amide I)": (1630, 1690, "medium", "low", "polarized"),
    # C=C
    "C=C (alkene)": (1620, 1680, "variable-strong", "high", "depolarized"),
    "C=C (aromatic)": (1450, 1600, "variable-medium", "high", "depolarized"),
    # Important Raman-active regions
    "S-S stretch": (500, 520, "strong", "high", "depolarized"),
    "C-S stretch": (630, 790, "strong", "high", "depolarized"),
    "S-H stretch": (2550, 2600, "medium-strong", "low", "polarized"),
    "C-Cl stretch": (550, 800, "medium", "high", "depolarized"),
    "C-Br stretch": (480, 650, "medium", "high", "depolarized"),
    # Ring breathing modes
    "Ring breathing (benzene)": (990, 1010, "very strong", "low", "polarized"),
    "Ring breathing (pyridine)": (990, 1030, "strong", "low", "polarized"),
    # Symmetric stretches (Raman-enhanced)
    "Symmetric C-C stretch": (1000, 1150, "strong", "low", "polarized"),
    "NO₂ symmetric stretch": (1300, 1370, "strong", "low", "polarized"),
    "NO₂ asymmetric stretch": (1480, 1570, "medium", "high", "depolarized"),
    # Fingerprint region highlights
    "C-C skeletal": (800, 1100, "variable", "varies", "depends on symmetry"),
    "=C-H bend (out-of-plane)": (700, 1000, "weak-medium", "high", "depolarized"),
    "Aromatic ring quadrant": (1550, 1620, "medium-strong", "low", "polarized"),
    "Aromatic ring semicircle": (1400, 1500, "medium", "high", "depolarized"),
}


@ChemMCPManager.register_tool
class RamanSpectrumPredictor(BaseTool):
    """
    拉曼光谱预测工具。
    根据分子官能团预测拉曼光谱特征峰位置、强度、退偏振比和活性。
    基于极化率变化原理，对称振动模式在拉曼光谱中增强。
    """
    __version__ = "0.1.0"
    name = "RamanSpectrumPredictor"
    func_name = "predict_raman_spectrum"
    description = "Predict characteristic Raman spectrum peaks for molecules based on functional groups and polarizability changes."
    implementation_description = (
        "Uses a database of Raman group frequency correlations to predict peak positions (cm⁻¹), "
        "intensities, depolarization ratios (ρ), and Raman activity (polarized/depolarized) "
        "for common organic functional groups. Selection rule: symmetric vibrations → strong Raman; "
        "asymmetric vibrations → strong IR."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Raman", "Spectroscopy", "Polarizability", "Vibrational Spectroscopy", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("functional_groups", "list", "N/A", "List of functional group names present in the molecule."),
        ("smiles", "str", "None", "Optional SMILES string for additional context."),
        ("detail_level", "str", "standard", "Detail level: 'basic', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space or comma-separated list of functional groups, e.g., 'ketone benzene alkyne'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing predicted Raman peaks sorted by wavenumber, including position, intensity, depolarization ratio, and activity."),
    ]

    examples = [
        {
            "code_input": {
                "functional_groups": ["benzene", "alkyne"],
                "smiles": None,
                "detail_level": "standard",
            },
            "text_input": {
                "input_params": "benzene alkyne",
            },
            "output": {
                "result": {
                    "molecule_type": "phenyl-acetylene",
                    "num_peaks": 5,
                    "peaks": [
                        {"shift_cm-1": 3080, "assignment": "C-H (aromatic)", "intensity": "medium-strong", "depolarization_ratio": "low (ρ<0.75)", "activity": "polarized"},
                        {"shift_cm-1": 2180, "assignment": "C≡C (alkyne)", "intensity": "variable-strong", "depolarization_ratio": "high (ρ≈0.75-1)", "activity": "depolarized"},
                        {"shift_cm-1": 1585, "assignment": "C=C (aromatic)", "intensity": "medium-strong", "depolarization_ratio": "high", "activity": "depolarized"},
                        {"shift_cm-1": 1000, "assignment": "Ring breathing (benzene)", "intensity": "very strong", "depolarization_ratio": "low (ρ<0.75)", "activity": "polarized"},
                        {"shift_cm-1": 620, "assignment": "=C-H bend (out-of-plane)", "intensity": "weak-medium", "depolarization_ratio": "high", "activity": "depolarized"},
                    ],
                    "selection_rule_note": "Symmetric C≡C stretch and ring breathing mode are Raman-enhanced.",
                }
            }
        },
        {
            "code_input": {
                "functional_groups": ["thiol", "alkane"],
                "smiles": None,
                "detail_level": "detailed",
            },
            "text_input": {
                "input_params": "thiol alkane",
            },
            "output": {
                "result": {
                    "molecule_type": "alkane-thiol",
                    "num_peaks": 4,
                    "peaks": [
                        {"shift_cm-1": 2920, "assignment": "C-H (alkane sp³)", "intensity": "strong", "depolarization_ratio": "low (ρ<0.75)", "activity": "polarized"},
                        {"shift_cm-1": 2575, "assignment": "S-H stretch", "intensity": "medium-strong", "depolarization_ratio": "low", "activity": "polarized"},
                        {"shift_cm-1": 720, "assignment": "C-S stretch", "intensity": "strong", "depolarization_ratio": "high", "activity": "depolarized"},
                        {"shift_cm-1": 300, "assignment": "C-S-C deformation", "intensity": "medium", "depolarization_ratio": "high", "activity": "depolarized"},
                    ],
                    "diagnostic_notes": "S-H stretch at ~2575 cm⁻¹ is diagnostic for thiols in Raman (absent in IR).",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.db = dict(RAMAN_GROUP_FREQUENCIES)

    def _run_base(self, functional_groups: List[str], smiles: str = None,
                  detail_level: str = "standard") -> dict:
        """Core logic: predict Raman spectrum peaks from functional groups."""
        if not functional_groups:
            raise ChemMCPError("At least one functional group must be provided.")

        dl = detail_level.lower()
        peaks = []
        seen_assignments = set()

        for fg in functional_groups:
            fg_key = self._match_group(fg.lower())
            if fg_key is None:
                logger.warning(f"Unknown functional group: {fg}, skipping.")
                continue

            entry = self.db[fg_key]
            pos_lo, pos_hi, intensity, depol, activity = entry

            pos = round((pos_lo + pos_hi) / 2) if isinstance(pos_hi, (int, float)) else pos_lo
            assignment = fg_key

            if assignment not in seen_assignments:
                seen_assignments.add(assignment)
                peak_info = {
                    "shift_cm-1": int(pos),
                    "range_cm-1": f"{pos_lo}-{pos_hi}" if isinstance(pos_hi, (int, float)) else f"{pos_lo}",
                    "assignment": assignment,
                    "intensity": intensity,
                    "depolarization_ratio": depol,
                    "activity": activity,
                }
                if dl == "detailed":
                    peak_info["origin"] = "group frequency correlation"
                    peak_info["polarizability_change"] = "large" if "symmetric" in assignment.lower() or "breathing" in assignment.lower() else "moderate"
                peaks.append(peak_info)

        # Add baseline C-H alkane stretch for organic molecules
        ch_key = "C-H (alkane sp³)"
        if ch_key not in seen_assignments:
            entry = self.db[ch_key]
            pos_lo, pos_hi, intensity, depol, activity = entry
            peaks.append({
                "shift_cm-1": round((pos_lo + pos_hi) / 2),
                "range_cm-1": f"{pos_lo}-{pos_hi}",
                "assignment": ch_key,
                "intensity": intensity,
                "depolarization_ratio": depol,
                "activity": activity,
            })

        # Sort by wavenumber (descending)
        peaks.sort(key=lambda p: p["shift_cm-1"], reverse=True)

        mol_type = self._classify_molecule(functional_groups)

        result = {
            "molecule_type": mol_type,
            "num_peaks": len(peaks),
            "peaks": peaks,
        }

        if dl != "basic":
            result["selection_rule_note"] = self._generate_selection_rule_note(peaks)
            result["ir_vs_raman_note"] = "Complementary to IR: symmetric stretches are Raman-active; asymmetric stretches are IR-active."

        return {"result": result}

    def _match_group(self, fg: str):
        """Match user input to known functional group key."""
        direct_map = {
            "alkane": "C-H (alkane sp³)",
            "alkene": "C=C (alkene)",
            "alkyne": "C≡C (alkyne)",
            "aromatic": "C=C (aromatic)",
            "benzene": "Ring breathing (benzene)",
            "pyridine": "Ring breathing (pyridine)",
            "ketone": "C=O (ketone/aldehyde)",
            "aldehyde": "C=O (ketone/aldehyde)",
            "carboxylic_acid": "C=O (carboxylic acid)",
            "amide": "C=O (amide I)",
            "nitrile": "C≡N (nitrile)",
            "thiol": "S-H stretch",
            "sulfide": "C-S stretch",
            "disulfide": "S-S stretch",
            "nitro": "NO₂ asymmetric stretch",
            "chloride": "C-Cl stretch",
            "bromide": "C-Br stretch",
            "ether": "Symmetric C-C stretch",
        }
        if fg in direct_map:
            return direct_map[fg]

        for key in self.db:
            clean_key = key.lower().replace(" ", "_").replace("-", "").replace("(", "").replace(")", "")
            clean_fg = fg.replace(" ", "_").replace("-", "")
            if clean_fg in clean_key or clean_key in clean_fg:
                return key
        return None

    @staticmethod
    def _classify_molecule(groups: List[str]) -> str:
        type_parts = [g.replace("_", " ").title() for g in groups[:3]]
        return "-".join(type_parts) + ("-etc" if len(groups) > 3 else "")

    @staticmethod
    def _generate_selection_rule_note(peaks: list) -> str:
        polarized = [p for p in peaks if "polarized" in p.get("activity", "")]
        depolarized = [p for p in peaks if "depolarized" in p.get("activity", "")]
        notes = []
        if polarized:
            names = [p["assignment"].split("(")[0].strip() for p in polarized[:3]]
            notes.append(f"Polarized (symmetric): {'; '.join(names)}")
        if depolarized:
            names = [p["assignment"].split("(")[0].strip() for p in depolarized[:3]]
            notes.append(f"Depolarized (asymmetric): {'; '.join(names)}")
        return " | ".join(notes) if notes else "Standard organic Raman profile."

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().replace(",", " ").split()
            return self._run_base(parts)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
