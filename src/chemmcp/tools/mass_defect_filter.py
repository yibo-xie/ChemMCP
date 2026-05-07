import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MassDefectFilter(BaseTool):
    """
    Mass defect filtering analysis for mass spectrometry data.
    Filters peaks by mass defect (fractional mass) window to reduce chemical noise.
    """
    __version__ = "0.1.0"
    name = "MassDefectFilter"
    func_name = "mass_defect_filter"
    description = "Perform mass defect filtering analysis on MS peak lists. Filters peaks by their fractional mass (mass defect) within a specified window to isolate compounds of interest and reduce background noise."
    implementation_description = "Calculates mass defect (MD = nominal mass - exact mass, or fractional part of exact mass) for each m/z value, then filters peaks whose mass defect falls outside the user-specified window. Supports Kendrick mass defect (KMD) analysis with customizable base."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Mass Spectrometry", "Mass Defect", "Kendrick Mass Defect", "Peak Filtering", "Data Analysis"]
    required_envs = []

    code_input_sig = [
        ("mz_list", "list", "N/A", "List of m/z values (floats) to filter."),
        ("intensity_list", "list", "N/A", "List of corresponding intensity values (floats). Must match mz_list length."),
        ("md_low", "float", "-0.2", "Lower bound of mass defect window (e.g., -0.2)."),
        ("md_high", "float", "0.5", "Upper bound of mass defect window (e.g., 0.5)."),
        ("mode", "str", "fractional", "Mass defect mode: 'fractional' (exact - floor) or 'kendrick' (KMD with base)."),
        ("base", "float", "14.01565", "Base for Kendrick mass defect calculation (default = CH2 = 14.01565). Default: 14.01565."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Space-separated: m1,i1 m2,i2 ... md_low md_high [mode] [base]. Example: '200.123,1000 201.234,500 -0.1 0.3 fractional' or pipe/comma-separated m/z,intensity pairs."),
    ]

    output_sig = [
        ("filtered_peaks", "list", "List of [mz, intensity] pairs that passed the mass defect filter."),
        ("rejected_peaks", "list", "List of [mz, intensity, mass_defect] pairs that were rejected."),
        ("mass_defects", "list", "Calculated mass defect for each input peak."),
        ("filter_stats", "dict", "Statistics: total_peaks, passed, rejected, pass_rate (%)."),
        ("md_window", "tuple", "The applied mass defect window (low, high)."),
        ("mode_used", "str", "The mass defect mode used for calculation."),
    ]

    examples = [
        {
            "code_input": {
                "mz_list": [200.1567, 300.2345, 400.1123, 500.3456, 150.9876],
                "intensity_list": [10000, 50000, 8000, 25000, 120000],
                "md_low": 0.05,
                "md_high": 0.25,
            },
            "text_input": {"input_string": "200.1567,10000 300.2345,50000 400.1123,8000 500.3456,25000 150.9876,120000 0.05 0.25"},
            "output": {
                "filtered_peaks": [[300.2345, 50000], [400.1123, 8000]],
                "rejected_peaks": [[200.1567, 10000, 0.1567], [500.3456, 25000, 0.3456], [150.9876, 120000, -0.0124]],
                "filter_stats": {"total_peaks": 5, "passed": 2, "rejected": 3, "pass_rate": 40.0},
                "md_window": (0.05, 0.25),
            }
        },
        {
            "code_input": {
                "mz_list": [100.0750, 200.1250, 300.1850],
                "intensity_list": [5000, 10000, 7500],
                "md_low": 0.06,
                "md_high": 0.20,
                "mode": "kendrick",
                "base": 14.01565,
            },
            "text_input": {"input_string": "100.0750,5000 200.1250,10000 300.1850,7500 0.06 0.20 kendrick"},
            "output": {
                "filtered_peaks": [[200.1250, 10000]],
                "rejected_peaks": [[100.0750, 5000, 0.0482], [300.1850, 7500, 0.2134]],
                "filter_stats": {"total_peaks": 3, "passed": 1, "rejected": 2, "pass_rate": 33.33},
                "md_window": (0.06, 0.20),
                "mode_used": "kendrick",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_mass_defect(self, mz: float, mode: str = "fractional", base: float = 14.01565) -> float:
        """Calculate mass defect for a given m/z value."""
        if mode == "fractional":
            return mz - math.floor(mz)
        elif mode == "kendrick":
            # KMD = round(mz * (base / nominal_base)) - mz * (base / nominal_base)
            # where nominal_base is the nearest integer of base
            nominal_base = round(base)
            km = mz * (base / nominal_base)
            return round(km) - km
        else:
            raise ChemMCPError(f"Unknown mass defect mode: {mode}. Use 'fractional' or 'kendrick'.")

    def _run_base(
        self,
        mz_list: List[float],
        intensity_list: List[float],
        md_low: float = -0.2,
        md_high: float = 0.5,
        mode: str = "fractional",
        base: float = 14.01565,
    ) -> dict:
        """Perform mass defect filtering on MS peak data."""
        if len(mz_list) != len(intensity_list):
            raise ChemMCPError("mz_list and intensity_list must have the same length.")
        if len(mz_list) == 0:
            raise ChemMCPError("Input peak lists cannot be empty.")

        if md_low > md_high:
            raise ChemMCPError(f"md_low ({md_low}) must be <= md_high ({md_high}).")

        # Calculate mass defects
        mass_defects = []
        filtered_peaks = []
        rejected_peaks = []

        for mz, intensity in zip(mz_list, intensity_list):
            md = self._calc_mass_defect(mz, mode, base)
            mass_defects.append(round(md, 6))

            if md_low <= md <= md_high:
                filtered_peaks.append([round(mz, 6), intensity])
            else:
                rejected_peaks.append([round(mz, 6), intensity, round(md, 6)])

        total = len(mz_list)
        passed = len(filtered_peaks)
        rejected = len(rejected_peaks)
        pass_rate = round(100.0 * passed / total, 2) if total > 0 else 0.0

        return {
            "filtered_peaks": filtered_peaks,
            "rejected_peaks": rejected_peaks,
            "mass_defects": mass_defects,
            "filter_stats": {
                "total_peaks": total,
                "passed": passed,
                "rejected": rejected,
                "pass_rate": pass_rate,
            },
            "md_window": (md_low, md_high),
            "mode_used": mode,
        }

    def _run_text(self, input_string: str) -> dict:
        """Parse text input string."""
        parts = input_string.strip().split()

        # Parse m/z,intensity pairs
        mz_list = []
        intensity_list = []
        idx = 0
        while idx < len(parts):
            p = parts[idx]
            if "," in p and "." in p.replace(",", "").replace("-", ""):
                sub = p.split(",")
                if len(sub) == 2:
                    try:
                        mz_list.append(float(sub[0]))
                        intensity_list.append(float(sub[1]))
                        idx += 1
                        continue
                    except ValueError:
                        pass
            break

        if not mz_list:
            raise ChemMCPError(
                f"Could not parse m/z,intensity pairs from '{input_string}'. "
                f"Format: 'mz1,int1 mz2,int2 ... [md_low] [md_high] [mode] [base]'"
            )

        md_low = float(parts[idx]) if idx < len(parts) else -0.2
        idx += 1
        md_high = float(parts[idx]) if idx < len(parts) else 0.5
        idx += 1
        mode = parts[idx] if idx < len(parts) else "fractional"
        idx += 1
        base = float(parts[idx]) if idx < len(parts) else 14.01565

        return self._run_base(mz_list, intensity_list, md_low, md_high, mode, base)
