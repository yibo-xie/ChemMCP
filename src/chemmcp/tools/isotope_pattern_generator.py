"""
Isotope Pattern Generator - generates theoretical isotope distribution patterns
for a given molecular formula using polynomial expansion algorithm.
"""

import logging
import re
from collections import defaultdict
from itertools import product as iter_product
from typing import Dict, List, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common isotope data: element -> [(mass, abundance), ...]
ISOTOPE_DATA: Dict[str, List[Tuple[float, float]]] = {
    "H":  [(1.007825, 0.999885), (2.014102, 0.000115)],
    "C":  [(12.000000, 0.9893), (13.003355, 0.0107)],
    "N":  [(14.003074, 0.99632), (15.000109, 0.00368)],
    "O":  [(15.994915, 0.99757), (16.999132, 0.00038), (17.999160, 0.00205)],
    "F":  [(18.998403, 1.0)],
    "Na": [(22.989770, 1.0)],
    "Si": [(27.976927, 0.92223), (28.976495, 0.04685), (29.973770, 0.03092)],
    "P":  [(30.973762, 1.0)],
    "S":  [(31.972071, 0.9499), (32.971458, 0.0075), (33.967867, 0.0425), (34.969034, 0.0001)],
    "Cl": [(34.968853, 0.7576), (36.965903, 0.2424)],
    "Br": [(78.918337, 0.5069), (80.916291, 0.4931)],
    "I":  [(126.904473, 1.0)],
    "B":  [(10.012937, 0.199), (11.009305, 0.801)],
    "Mg": [(23.985042, 0.7899), (24.985837, 0.1000), (25.982593, 0.1101)],
    "Ca": [(39.962591, 0.96941), (41.958618, 0.00647), (42.958767, 0.00135), (43.955481, 0.02086), (44.956186, 0.00187), (47.952534, 0.00004)],
    "Fe": [(53.939610, 0.05845), (55.934938, 0.91754), (56.935394, 0.02119), (57.933276, 0.00282)],
    "Cu": [(62.929598, 0.69017), (64.927791, 0.30983)],
    "Zn": [(63.929142, 0.4917), (65.926033, 0.2767), (66.927127, 0.0410), (67.924844, 0.1850), (69.925319, 0.0061)],
    "K":  [(38.963707, 0.932581), (39.963999, 0.000117), (40.961826, 0.067302)],
}


def _parse_formula(formula: str) -> Dict[str, int]:
    """Parse molecular formula like 'C6H12O6' into {'C': 6, 'H': 12, 'O': 6}."""
    pattern = r"([A-Z][a-z]?)(\d*)"
    matches = re.findall(pattern, formula)
    result: Dict[str, int] = {}
    for elem, count in matches:
        if not elem:
            continue
        result[elem] = result.get(elem, 0) + (int(count) if count else 1)
    return result


@ChemMCPManager.register_tool
class IsotopePatternGenerator(BaseTool):
    __version__      = "0.1.0"
    name             = "IsotopePatternGenerator"
    func_name        = "generate_isotope_pattern"
    description      = "Generate theoretical isotope distribution pattern for a given molecular formula."
    implementation_description = "Uses polynomial expansion algorithm with natural isotope abundances to compute mass-to-charge ratios and relative intensities."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Isotopes", "Mass Spectrometry", "Molecular Formula"]
    required_envs    = []

    code_input_sig   = [
        ("molecular_formula", "str", "N/A", "Molecular formula string (e.g., 'C6H12O6', 'CHCl3')."),
        ("charge", "int", "0", "Charge of the ion (e.g., +1, -1)."),
        ("min_abundance", "float", "0.01", "Minimum relative abundance threshold (0-1). Peaks below this are filtered out."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'molecular_formula charge min_abundance'. Example: 'C6H12O6 0 0.01'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: formula, monoisotopic_mass, charge, peaks (list of {mz, abundance_pct, intensity}), peak_count, nominal_mass."),
    ]

    examples         = [
        {
            "code_input": {"molecular_formula": "CHCl3", "charge": 0, "min_abundance": 0.01},
            "text_input": {"input_params": "CHCl3 0 0.01"},
            "output": {
                "result": {
                    "formula": "CHCl3",
                    "monoisotopic_mass": 117.9,
                    "charge": 0,
                    "peak_count": 3,
                    "nominal_mass": 118,
                    "peaks": [
                        {"mz": 117.9, "abundance_pct": 100.0, "intensity": 1.0},
                        {"mz": 119.9, "abundance_pct": 96.3, "intensity": 0.963},
                        {"mz": 121.9, "abundance_pct": 31.1, "intensity": 0.311},
                    ]
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, molecular_formula: str, charge: int = 0, min_abundance: float = 0.01) -> dict:
        """Core logic: compute isotope pattern via convolution."""
        elements = _parse_formula(molecular_formula)
        if not elements:
            raise ChemMCPError(f"Cannot parse molecular formula: '{molecular_formula}'")

        # Build isotope distributions per element (with multiplicities)
        # Each element's contribution: convolve its own isotopes (count times)
        total_dist: Dict[float, float] = {0.0: 1.0}  # mass -> probability

        for elem, count in elements.items():
            if elem not in ISOTOPE_DATA:
                raise ChemMCPError(f"No isotope data for element '{elem}'. Available: {sorted(ISOTOPE_DATA.keys())}")

            isotopes = ISOTOPE_DATA[elem]
            # Single element distribution
            elem_dist: Dict[float, float] = {}
            for mass, ab in isotopes:
                elem_dist[mass] = ab

            # Raise to power `count` by repeated convolution
            powered = self._power_dist(elem_dist, count)

            # Convolve into total
            total_dist = self._convolve(total_dist, powered)

        if not total_dist:
            raise ChemMCPError("Isotope pattern computation produced no results.")

        # Normalize and convert to peaks
        max_prob = max(total_dist.values())
        peaks = []
        for mass in sorted(total_dist.keys()):
            prob = total_dist[mass] / max_prob
            if prob * 100 >= min_abundance * 100:
                mz = mass - charge * 0.00054858  # electron mass correction
                peaks.append({
                    "mz": round(mz, 4),
                    "abundance_pct": round(prob * 100, 2),
                    "intensity": round(prob, 6),
                })

        if not peaks:
            raise ChemMCPError("No peaks above minimum abundance threshold.")

        return {
            "formula": molecular_formula,
            "monoisotopic_mass": round(peaks[0]["mz"], 4),
            "charge": charge,
            "peak_count": len(peaks),
            "nominal_mass": int(round(peaks[0]["mz"])),
            "peaks": peaks,
        }

    @staticmethod
    def _convolve(d1: Dict[float, float], d2: Dict[float, float]) -> Dict[float, float]:
        """Convolve two mass->probability distributions."""
        result: Dict[float, float] = defaultdict(float)
        for m1, p1 in d1.items():
            for m2, p2 in d2.items():
                result[round(m1 + m2, 6)] += p1 * p2
        return dict(result)

    @staticmethod
    def _power_dist(dist: Dict[float, float], power: int) -> Dict[float, float]:
        """Raise a distribution to a positive integer power via repeated squaring."""
        if power <= 0:
            return {0.0: 1.0}
        result: Dict[float, float] = {0.0: 1.0}
        base = dist.copy()
        while power > 0:
            if power % 2 == 1:
                result = IsotopePatternGenerator._convolve(result, base)
            base = IsotopePatternGenerator._convolve(base, base)
            power //= 2
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            formula = parts[0]
            charge = int(parts[1]) if len(parts) > 1 else 0
            min_ab = float(parts[2]) if len(parts) > 2 else 0.01
            return self._run_base(formula, charge, min_ab)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'formula charge min_abundance'")
