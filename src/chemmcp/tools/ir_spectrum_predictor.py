import logging
import math
from typing import List, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Reference IR group frequencies (cm⁻¹) — common functional groups
# Based on standard organic spectroscopy data (Silverstein, Pavia, etc.)
IR_GROUP_FREQUENCIES = {
    # O-H stretches
    "O-H (alcohol, free)": (3600, 3750, "sharp", "broad" if False else "sharp"),
    "O-H (alcohol, H-bonded)": (3200, 3600, "broad", "broad"),
    "O-H (carboxylic acid)": (2500, 3300, "very broad", "very broad"),
    # N-H stretches
    "N-H (primary amine)": (3400, 3500, "doublet", "medium"),
    "N-H (secondary amine)": (3300, 3400, "singlet", "weak"),
    "N-H (amide)": (3100, 3500, "doublet/singlet", "medium"),
    # C-H stretches
    "C-H (alkane sp³)": (2850, 2960, "medium", "medium"),
    "C-H (alkene sp²)": (3010, 3100, "medium", "medium"),
    "C-H (alkyne sp)": (3300, "sharp", "strong"),
    "C-H (aldehyde)": (2650, 2820, "doublet (Fermi)", "medium"),
    # Triple bonds
    "C≡C (alkyne)": (2100, 2260, "weak/absent if symmetrical", "variable"),
    "C≡N (nitrile)": (2220, 2260, "sharp", "medium-strong"),
    # Double bonds: C=O
    "C=O (saturated aliphatic)": (1710, 1745, "sharp", "strong"),
    "C=O (α,β-unsaturated)": (1680, 1710, "sharp", "strong"),
    "C=O (conjugated ketone)": (1665, 1685, "sharp", "strong"),
    "C=O (carboxylic acid)": (1710, 1760, "sharp", "strong"),
    "C=O (ester)": (1735, 1750, "sharp", "strong"),
    "C=O (amide I band)": (1630, 1690, "sharp", "strong"),
    "C=O (anhydride)": (1740, 1820, "doublet", "strong"),
    "C=O (acid chloride)": (1770, 1815, "sharp", "strong"),
    "C=O (aldehyde)": (1720, 1740, "sharp", "strong"),
    # C=C
    "C=C (alkene)": (1620, 1680, "variable", "variable"),
    "C=C (aromatic)": (1450, 1600, "2-4 peaks", "variable"),
    # Other important regions
    "N-H bend (amine/amide)": (1550, 1640, "medium", "medium"),
    "C-H bend (aldehyde)": (1340, 1390, "medium-strong", "medium"),
    "C-H bend (alkane)": (1465, 1370, "two peaks", "medium"),
    "C-O (alcohol/ether)": (1000, 1260, "strong", "strong"),
    "C-O (ester)": (1300, 1000, "two peaks", "strong"),
    "C-O (carboxylic acid)": (1210, 1320, "strong", "strong"),
    "NO₂ (symmetric)": (1300, 1370, "strong", "strong"),
    "NO₂ (asymmetric)": (1500, 1570, "strong", "strong"),
    "C-F": (1000, 1400, "strong", "strong"),
    "C-Cl": (600, 800, "strong", "strong"),
    "C-Br": (500, 600, "strong", "medium"),
    "≡C-H bend (alkyne)": (620, 700, "strong", "medium"),
    "=C-H bend (alkene)": (650, 1000, "medium", "medium"),
    "aromatic out-of-plane C-H bend": (690, 900, "characteristic pattern", "strong"),
}


@ChemMCPManager.register_tool
class IrSpectrumPredictor(BaseTool):
    """
    红外光谱特征峰预测工具。
    根据分子中的官能团预测红外光谱的特征吸收峰位置、强度和形状。
    """
    __version__ = "0.1.0"
    name = "IrSpectrumPredictor"
    func_name = "predict_ir_spectrum"
    description = "Predict characteristic infrared (IR) absorption peaks for molecules based on their functional groups."
    implementation_description = "Uses a comprehensive database of group frequency correlations to predict IR peak positions (wavenumber in cm⁻¹), intensities, and peak shapes for common organic functional groups."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "IR", "Infrared", "Functional Groups", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("functional_groups", "list", "N/A", "List of functional group names present in the molecule."),
        ("smiles", "str", "None", "Optional SMILES string for additional context."),
        ("detail_level", "str", "standard", "Detail level: 'basic', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space or comma-separated list of functional groups, e.g., 'ketone alcohol alkene'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing predicted IR peaks sorted by wavenumber, including position, intensity, shape, and assignment."),
    ]

    examples = [
        {
            "code_input": {
                "functional_groups": ["ketone", "alcohol"],
                "smiles": None,
                "detail_level": "standard",
            },
            "text_input": {
                "input_params": "ketone alcohol",
            },
            "output": {
                "result": {
                    "molecule_type": "hydroxy-ketone",
                    "num_peaks": 4,
                    "peaks": [
                        {"position_cm-1": 3400, "assignment": "O-H stretch (H-bonded alcohol)", "intensity": "broad strong", "shape": "broad"},
                        {"position_cm-1": 2920, "assignment": "C-H stretch (alkane)", "intensity": "medium", "shape": "sharp"},
                        {"position_cm-1": 1715, "assignment": "C=O stretch (saturated ketone)", "intensity": "strong", "shape": "sharp"},
                        {"position_cm-1": 1150, "assignment": "C-O stretch (alcohol)", "intensity": "strong", "shape": "broad"},
                    ],
                    "fingerprint_region_note": "Check 1500-400 cm⁻¹ fingerprint region for unique identification.",
                }
            }
        },
        {
            "code_input": {
                "functional_groups": ["carboxylic_acid", "aromatic"],
                "smiles": None,
                "detail_level": "detailed",
            },
            "text_input": {
                "input_params": "carboxylic_acid aromatic",
            },
            "output": {
                "result": {
                    "molecule_type": "aromatic carboxylic acid",
                    "num_peaks": 6,
                    "peaks": [
                        {"position_cm-1": 2900, "assignment": "O-H stretch (carboxylic acid dimer)", "intensity": "very broad", "shape": "very broad"},
                        {"position_cm-1": 3060, "assignment": "C-H stretch (aromatic)", "intensity": "medium", "shape": "sharp"},
                        {"position_cm-1": 1700, "assignment": "C=O stretch (carboxylic acid)", "intensity": "strong", "shape": "sharp"},
                        {"position_cm-1": 1600, "assignment": "C=C stretch (aromatic ring)", "intensity": "variable", "shape": "sharp"},
                        {"position_cm-1": 1460, "assignment": "C=C stretch (aromatic ring)", "intensity": "variable", "shape": "sharp"},
                        {"position_cm-1": 1250, "assignment": "C-O stretch (carboxylic acid)", "intensity": "strong", "shape": "broad"},
                    ],
                    "diagnostic_notes": "Very broad O-H (2500-3300) + strong C=O ~1700 is diagnostic for carboxylic acids.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.db = dict(IR_GROUP_FREQUENCIES)

    def _run_base(self, functional_groups: List[str], smiles: str = None,
                  detail_level: str = "standard") -> dict:
        """Core logic: predict IR spectrum peaks from functional groups."""
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
            pos_lo, pos_hi, shape, intensity = self._parse_entry(entry)

            pos = round((pos_lo + pos_hi) / 2) if isinstance(pos_hi, (int, float)) else pos_lo
            assignment = fg_key

            if assignment not in seen_assignments:
                seen_assignments.add(assignment)
                peak_info = {
                    "position_cm-1": int(pos),
                    "range_cm-1": f"{pos_lo}-{pos_hi}" if isinstance(pos_hi, (int, float)) else f"{pos_lo}",
                    "assignment": assignment,
                    "intensity": intensity,
                    "shape": shape,
                }
                if dl == "detailed":
                    peak_info["origin"] = "group frequency correlation"
                    peak_info["reliability"] = "high" if "C=O" in assignment or "O-H" in assignment or "C≡" in assignment else "moderate"
                peaks.append(peak_info)

        # Always add C-H alkane stretch as baseline for organic molecules
        ch_key = "C-H (alkane sp³)"
        if ch_key not in seen_assignments and any(
            "alkane" not in g.lower() and "alkene" not in g.lower() and "alkyne" not in g.lower()
            for g in functional_groups
        ):
            entry = self.db[ch_key]
            pos_lo, pos_hi, shape, intensity = self._parse_entry(entry)
            peaks.append({
                "position_cm-1": round((pos_lo + pos_hi) / 2),
                "range_cm-1": f"{pos_lo}-{pos_hi}",
                "assignment": ch_key,
                "intensity": intensity,
                "shape": shape,
            })

        # Sort by wavenumber (descending)
        peaks.sort(key=lambda p: p["position_cm-1"], reverse=True)

        mol_type = self._classify_molecule(functional_groups)

        result = {
            "molecule_type": mol_type,
            "num_peaks": len(peaks),
            "peaks": peaks,
        }

        if dl != "basic":
            result["diagnostic_notes"] = self._generate_diagnostic_notes(peaks, functional_groups)
            result["fingerprint_region_note"] = (
                "Check 1500-400 cm⁻¹ fingerprint region for unique molecular identification."
            )

        return {"result": result}

    def _match_group(self, fg: str):
        """Match user input to known functional group key."""
        direct_map = {
            "ketone": "C=O (saturated aliphatic)",
            "aldehyde": "C=O (aldehyde)",
            "alcohol": "O-H (alcohol, H-bonded)",
            "carboxylic_acid": "O-H (carboxylic acid)",
            "carboxylic acid": "O-H (carboxylic acid)",
            "amine": "N-H (primary amine)",
            "primary_amine": "N-H (primary amine)",
            "secondary_amine": "N-H (secondary amine)",
            "amide": "N-H (amide)",
            "alkene": "C=C (alkene)",
            "alkyne": "C≡C (alkyne)",
            "aromatic": "C=C (aromatic)",
            "benzene": "C=C (aromatic)",
            "ester": "C=O (ester)",
            "nitrile": "C≡N (nitrile)",
            "anhydride": "C=O (anhydride)",
            "acid_chloride": "C=O (acid chloride)",
            "nitro": "NO₂ (asymmetric)",
            "halide": "C-Cl",
            "ether": "C-O (alcohol/ether)",
        }
        if fg in direct_map:
            return direct_map[fg]

        # Fuzzy match
        for key in self.db:
            clean_key = key.lower().replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
            clean_fg = fg.replace(" ", "_").replace("-", "_")
            if clean_fg in clean_key or clean_key in clean_fg:
                return key
        return None

    @staticmethod
    def _parse_entry(entry):
        """Parse a DB entry into (lo, hi, shape, intensity)."""
        if len(entry) == 4:
            lo, hi, shape, intensity = entry
        elif len(entry) == 3:
            lo, hi, intensity = entry
            shape = "unknown"
        else:
            lo, intensity = entry[0], entry[-1]
            hi = lo
            shape = "unknown"
        return lo, hi, shape, intensity

    @staticmethod
    def _classify_molecule(groups: List[str]) -> str:
        gs = [g.lower().replace(" ", "_") for g in groups]
        type_parts = []
        for g in groups:
            type_parts.append(g.replace("_", " ").title())
        return "-".join(type_parts[:3]) + ("-etc" if len(groups) > 3 else "")

    @staticmethod
    def _generate_diagnostic_notes(peaks: List[dict], groups: List[str]) -> str:
        notes = []
        strong_peaks = [p for p in peaks if "strong" in p.get("intensity", "")]
        if strong_peks := [p["assignment"] for p in strong_peaks]:
            notes.append(f"Strongest diagnostic peaks: {'; '.join(strong_peks[:3])}")
        if any("O-H" in g for g in groups):
            notes.append("Broad O-H stretch indicates hydrogen bonding — check for intermolecular association.")
        if any("C=O" in p["assignment"] for p in peaks):
            carbonyl = [p for p in peaks if "C=O" in p["assignment"]]
            if carbonyl:
                notes.append(f"Carbonyl region: {carbonyl[0]['range_cm-1']} cm⁻¹ — check conjugation effects.")
        return " | ".join(notes) if notes else "Standard organic IR profile."

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().replace(",", " ").split()
            return self._run_base(parts)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
