import logging
import math
from typing import List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Raman activity reference data for common functional groups
# Based on standard Raman spectroscopy references
RAMAN_GROUP_DATA = {
    # (frequency range cm⁻¹, typical intensity, Raman activity level, notes)
    "O-H stretch": {"range": (3200, 3750), "intensity": "weak", "activity": "low", "notes": "Generally weak in Raman"},
    "N-H stretch": {"range": (3300, 3500), "intensity": "medium", "activity": "moderate", "notes": "Amide I strong in Raman"},
    "C-H stretch (alkane)": {"range": (2850, 2960), "intensity": "medium-strong", "activity": "moderate", "notes": ""},
    "C≡C stretch": {"range": (2100, 2260), "intensity": "strong", "activity": "high", "notes": "Symmetric triple bonds are Raman-active"},
    "C≡N stretch": {"range": (2220, 2260), "intensity": "strong", "activity": "high", "notes": "Nitrile group gives sharp Raman peak"},
    "C=O stretch": {"range": (1650, 1780), "intensity": "weak-medium", "activity": "low-moderate", "notes": "Carbonyls generally weak in Raman unless conjugated"},
    "C=C stretch (alkene)": {"range": (1620, 1680), "intensity": "medium-strong", "activity": "moderate-high", "notes": "Symmetric alkenes are Raman-active"},
    "C=C aromatic": {"range": (1550, 1620), "intensity": "strong", "activity": "high", "notes": "Aromatic ring breathing modes are very Raman-active"},
    "S-S stretch": {"range": (500, 550), "intensity": "strong", "activity": "high", "notes": "Disulfide bonds are Raman markers for proteins"},
    "S-H stretch": {"range": (2550, 2600), "intensity": "strong", "activity": "high", "notes": "Unique region with no interference"},
    "C-S stretch": {"range": (630, 790), "intensity": "strong", "activity": "high", "notes": ""},
    "C-Cl stretch": {"range": (550, 800), "intensity": "medium", "activity": "moderate", "notes": ""},
    "Ring breathing (benzene)": {"range": (990, 1010), "intensity": "very strong", "activity": "very high", "notes": "Monosubstituted benzene ring breathing ~992 cm⁻¹"},
    "Ring breathing (pyridine)": {"range": (990, 1030), "intensity": "very strong", "activity": "very high", "notes": "~991 cm⁻¹"},
    "C=C-C bend": {"range": (200, 600), "intensity": "variable", "activity": "variable", "notes": "Fingerprint region"},
}


@ChemMCPManager.register_tool
class RamanActivity(BaseTool):
    """
    拉曼活性判断和频率预测工具。
    判断分子或官能团的拉曼活性，预测特征拉曼位移峰位置和强度。
    """
    __version__ = "0.1.0"
    name = "RamanActivity"
    func_name = "predict_raman_activity"
    description = "Predict Raman spectroscopy activity and characteristic peak positions/frequencies for molecules and functional groups."
    implementation_description = "Uses molecular symmetry analysis rules (mutual exclusion rule) and a database of group frequency correlations to predict Raman-active vibrational modes."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "Raman", "Vibrational Spectroscopy", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("functional_groups", "list", "N/A", "List of functional groups or bond types to analyze."),
        ("molecule_symmetry", "str", "unknown", "Point group symmetry of molecule (e.g., 'C2v', 'D∞h', 'Td'). Helps determine IR/Raman mutual exclusion."),
        ("include_selection_rules", "bool", "True", "Whether to include symmetry-based selection rule analysis."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated list of functional groups, e.g., 'alkene nitrile aromatic'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing predicted Raman peaks, activity assessment, selection rule analysis, and comparison with IR."),
    ]

    examples = [
        {
            "code_input": {
                "functional_groups": ["alkene", "aromatic"],
                "molecule_symmetry": "unknown",
                "include_selection_rules": True,
            },
            "text_input": {
                "input_params": "alkene aromatic",
            },
            "output": {
                "result": {
                    "molecule_type": "alkene-aromatic",
                    "raman_active_modes": 4,
                    "peaks": [
                        {"position_cm-1": 3065, "assignment": "C-H stretch (aromatic)", "intensity": "medium-strong", "raman_activity": "high"},
                        {"position_cm-1": 1600, "assignment": "C=C aromatic quadrant stretch", "intensity": "strong", "raman_activity": "high"},
                        {"position_cm-1": 1000, "assignment": "Ring breathing mode", "intensity": "very strong", "raman_activity": "very_high"},
                        {"position_cm-1": 1650, "assignment": "C=C stretch (alkene)", "intensity": "medium-strong", "raman_activity": "moderate-high"},
                    ],
                    "selection_rule_notes": "For molecules with center of inversion, mutually exclusive: IR-active modes are Raman-inactive and vice versa.",
                    "ir_vs_raman_complementarity": "Raman excels at symmetric vibrations (C=C, ring breathing); IR excels at dipole-changing modes (C=O, O-H).",
                }
            }
        },
        {
            "code_input": {
                "functional_groups": ["nitrile", "disulfide"],
                "molecule_symmetry": "unknown",
                "include_selection_rules": True,
            },
            "text_input": {
                "input_params": "nitrile disulfide",
            },
            "output": {
                "result": {
                    "molecule_type": "nitrile-disulfide",
                    "raman_active_modes": 2,
                    "peaks": [
                        {"position_cm-1": 2240, "assignment": "C≡N stretch (nitrile)", "intensity": "strong", "raman_activity": "high"},
                        {"position_cm-1": 525, "assignment": "S-S stretch (disulfide)", "intensity": "strong", "raman_activity": "high"},
                    ],
                    "diagnostic_note": "Both C≡N and S-S stretches are excellent Raman markers due to large polarizability changes.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.db = dict(RAMAN_GROUP_DATA)
        # Symmetry groups with center of inversion → mutual exclusion
        self.inversion_groups = {"d∞h", "d_inf_h", "d2h", "d3h", "d4h", "d6h", "ci", "oh"}

    def _run_base(self, functional_groups: List[str], molecule_symmetry: str = "unknown",
                  include_selection_rules: bool = True) -> dict:
        """Core logic: predict Raman activity."""
        if not functional_groups:
            raise ChemMCPError("At least one functional group must be provided.")

        sym = molecule_symmetry.lower().strip()
        has_inversion = sym in self.inversion_groups

        peaks = []
        seen = set()

        for fg in functional_groups:
            matched_keys = self._match_group(fg.lower())
            for key in matched_keys:
                if key in seen:
                    continue
                seen.add(key)
                data = self.db[key]
                rlo, rhi = data["range"]
                pos = round((rlo + rhi) / 2)

                peaks.append({
                    "position_cm-1": pos,
                    "range_cm-1": f"{rlo}-{rhi}",
                    "assignment": key,
                    "intensity": data["intensity"],
                    "raman_activity": data["activity"],
                    "notes": data["notes"],
                })

        # Sort by position descending
        peaks.sort(key=lambda p: p["position_cm-1"], reverse=True)

        result = {
            "molecule_type": "-".join(functional_groups[:3]),
            "raman_active_modes": len(peaks),
            "peaks": peaks,
        }

        if include_selection_rules:
            result["selection_rule_analysis"] = self._selection_rules(sym, has_inversion)
            result["ir_vs_raman_complementarity"] = (
                "Raman: best for symmetric, nonpolar bonds (C≡C, C=C, S-S, ring breathing). "
                "IR: best for polar bonds with dipole changes (C=O, O-H, N-H). "
                "Combined use provides complete vibrational information."
            )

        return {"result": result}

    def _match_group(self, fg: str) -> List[str]:
        """Match user input to known Raman groups."""
        direct = {
            "alkene": ["C=C stretch (alkene)", "C-H stretch (alkane)"],
            "aromatic": ["C=C aromatic", "Ring breathing (benzene)"],
            "nitrile": ["C≡N stretch"],
            "alkyne": ["C≡C stretch"],
            "disulfide": ["S-S stretch"],
            "thiol": ["S-H stretch"],
            "thioether": ["C-S stretch"],
            "alcohol": ["O-H stretch"],
            "amine": ["N-H stretch"],
            "ketone": ["C=O stretch"],
            "halide": ["C-Cl stretch"],
            "benzene": ["C=C aromatic", "Ring breathing (benzene)"],
            "pyridine": ["Ring breathing (pyridine)"],
        }
        if fg in direct:
            return direct[fg]

        matched = []
        for key in self.db:
            clean_key = key.lower().replace(" ", "_").replace("-", "").replace("(", "").replace(")", "")
            if fg.replace(" ", "_").replace("-", "") in clean_key or clean_key in fg.replace(" ", "_"):
                matched.append(key)
        return matched if matched else [fg]

    def _selection_rules(self, sym: str, has_inversion: bool) -> dict:
        """Generate symmetry-based selection rule analysis."""
        analysis = {
            "molecule_symmetry": sym if sym != "unknown" else "not specified",
            "has_center_of_inversion": has_inversion,
            "mutual_exclusion_applies": has_inversion,
        }

        if has_inversion:
            analysis["rule"] = (
                "Mutual Exclusion Rule: For centrosymmetric molecules, "
                "vibrations that are IR-active are Raman-inactive, and vice versa."
            )
            analysis["implication"] = (
                "A mode cannot be both IR and Raman active simultaneously. "
                "Complete vibrational characterization requires both techniques."
            )
        elif sym != "unknown":
            analysis["rule"] = (
                f"For point group {sym.upper()}, some modes may be both IR and Raman active. "
                "Check character table for specific mode symmetries."
            )
            analysis["implication"] = "IR and Raman provide complementary but overlapping information."
        else:
            analysis["rule"] = "Specify molecular symmetry for detailed selection rule analysis."
            analysis["implication"] = "General prediction based on group frequencies only."

        return analysis

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().replace(",", " ").split()
            return self._run_base(parts)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
