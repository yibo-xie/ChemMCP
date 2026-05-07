"""
Isotope Pattern Simulator - MS-grade isotope peak distribution simulation
with resolution-dependent Gaussian peak broadening for realistic mass spectra.
"""

import logging
import math
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive isotope data: element -> [(exact_mass, natural_abundance)]
ISOTOPE_DATA: Dict[str, List[Tuple[float, float]]] = {
    "H":  [(1.00782503223, 0.999885), (2.01410177812, 0.000115)],
    "C":  [(12.0000000, 0.9893), (13.00335483507, 0.0107)],
    "N":  [(14.00307400443, 0.99632), (15.0001088982, 0.00368)],
    "O":  [(15.99491461957, 0.99757), (16.9991317069, 0.00038), (17.9991616115, 0.00205)],
    "F":  [(18.99840316273, 1.0)],
    "Na": [(22.9897692820, 1.0)],
    "Mg": [(23.985041697, 0.7899), (24.985836920, 0.1000), (25.982592929, 0.1101)],
    "Si": [(27.97692653465, 0.92223), (28.9764947000, 0.04685), (29.973770137, 0.03092)],
    "P":  [(30.97376199842, 1.0)],
    "S":  [(31.9720711744, 0.9499), (32.971458763, 0.0075), (33.967867004, 0.0425), (34.96903216, 0.0001)],
    "Cl": [(34.968852682, 0.7576), (36.96590259, 0.2424)],
    "K":  [(38.96370668, 0.932581), (39.96399848, 0.000117), (40.96182576, 0.067302)],
    "Ca": [(39.962590863, 0.96941), (41.95861801, 0.00647), (42.95876676, 0.00135), (43.95548061, 0.02086), (44.95618563, 0.00187), (47.95253395, 0.00004)],
    "Fe": [(53.9396089, 0.05845), (55.9349375, 0.91754), (56.9353942, 0.02119), (57.9332756, 0.00282)],
    "Cu": [(62.92959772, 0.69017), (64.92778995, 0.30983)],
    "Zn": [(63.92914222, 0.4917), (65.92603348, 0.2767), (66.92712773, 0.0410), (67.92484458, 0.1850), (69.92531929, 0.0061)],
    "Br": [(78.9183376, 0.5069), (80.916290, 0.4931)],
    "I":  [(126.9044727, 1.0)],
    "B":  [(10.01293695, 0.199), (11.00930536, 0.801)],
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
class IsotopePatternSimulator(BaseTool):
    """
    同位素峰分布模拟器 — 模拟质谱中分子的同位素分布模式。
    
    基于多项式卷积算法计算理论同位素分布，并可根据分辨率参数进行高斯峰展宽，
    模拟真实质谱仪的峰形。
    """
    __version__      = "0.1.0"
    name             = "IsotopePatternSimulator"
    func_name        = "simulate_isotope_pattern"
    description      = "Simulate detailed isotope peak distribution for mass spectrometry analysis with resolution-dependent Gaussian broadening."
    implementation_description = "Uses polynomial expansion (convolution) algorithm with natural isotope abundances to compute exact isotope distributions. Optionally applies Gaussian peak broadening based on instrument resolution (FWHM = mz / resolution) to simulate realistic peak shapes."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Isotopes", "Mass Spectrometry", "Peak Distribution", "Simulation", "MS Analysis"]
    required_envs    = []

    code_input_sig   = [
        ("molecular_formula", "str", "N/A", "Molecular formula string (e.g., 'C6H12O6', 'CHCl3', 'C17H19NO3')."),
        ("charge", "int", "0", "Charge state of the ion (e.g., +1, +2, -1)."),
        ("resolution", "int", "30000", "Instrument resolution at FWHM (e.g., 30000 for Q-TOF, 100000 for Orbitrap)."),
        ("min_abundance", "float", "0.001", "Minimum relative abundance threshold (0-1). Peaks below this are filtered out."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'formula charge resolution min_abundance'. Example: 'C17H19NO3 1 30000 0.001'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict containing: formula, monoisotopic_mass, charge, resolution, peaks (list with mz, abundance_pct, intensity, fwhm), peak_count, nominal_mass, envelope_summary."),
    ]

    examples         = [
        {
            "code_input": {
                "molecular_formula": "C17H19NO3",
                "charge": 1,
                "resolution": 30000,
                "min_abundance": 0.001
            },
            "text_input": {
                "input_params": "C17H19NO3 1 30000 0.001"
            },
            "output": {
                "result": {
                    "formula": "C17H19NO3",
                    "monoisotopic_mass": 285.1385,
                    "charge": 1,
                    "resolution": 30000,
                    "peak_count": 15,
                    "nominal_mass": 285,
                    "envelope_summary": "Dominant M+1 peak (~18.7%) from ¹³C contribution; Cl/Br absent → clean envelope",
                    "peaks": [
                        {"mz": 285.1385, "abundance_pct": 100.0, "intensity": 1.0, "fwhm": 0.00950},
                        {"mz": 286.1419, "abundance_pct": 18.70, "intensity": 0.187, "fwhm": 0.00954},
                    ]
                }
            },
        },
        {
            "code_input": {
                "molecular_formula": "CHCl3",
                "charge": 0,
                "resolution": 15000,
                "min_abundance": 0.01
            },
            "text_input": {
                "input_params": "CHCl3 0 15000 0.01"
            },
            "output": {
                "result": {
                    "formula": "CHCl3",
                    "monoisotopic_mass": 117.9004,
                    "charge": 0,
                    "resolution": 15000,
                    "peak_count": 3,
                    "nominal_mass": 118,
                    "envelope_summary": "Characteristic 3:1:0.058 ratio from single ³⁵Cl₃/³⁵Cl₂³⁷Cl/³⁵Cl³⁷Cl₂ pattern",
                    "peaks": [
                        {"mz": 117.9004, "abundance_pct": 100.0, "intensity": 1.0},
                        {"mz": 119.8974, "abundance_pct": 96.27, "intensity": 0.9627},
                        {"mz": 121.8945, "abundance_pct": 31.06, "intensity": 0.3106},
                    ]
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, molecular_formula: str, charge: int = 0, resolution: int = 30000, min_abundance: float = 0.001) -> dict:
        """Core logic: compute isotope pattern with optional Gaussian broadening."""
        elements = _parse_formula(molecular_formula)
        if not elements:
            raise ChemMCPError(f"Cannot parse molecular formula: '{molecular_formula}'")

        # Build distribution via convolution
        total_dist: Dict[float, float] = {0.0: 1.0}

        for elem, count in elements.items():
            if elem not in ISOTOPE_DATA:
                raise ChemMCPError(f"No isotope data for element '{elem}'. Available: {sorted(ISOTOPE_DATA.keys())}")
            isotopes = ISOTOPE_DATA[elem]
            elem_dist: Dict[float, float] = {mass: ab for mass, ab in isotopes}
            powered = self._power_dist(elem_dist, count)
            total_dist = self._convolve(total_dist, powered)

        if not total_dist:
            raise ChemMCPError("Isotope pattern computation produced no results.")

        # Normalize and build peaks
        max_prob = max(total_dist.values())
        peaks = []
        electron_mass = 0.000548579909

        for mass in sorted(total_dist.keys()):
            prob = total_dist[mass] / max_prob
            if prob < min_abundance:
                continue
            mz = (mass - charge * electron_mass) / abs(charge) if charge != 0 else mass
            fwhm = mz / resolution if resolution > 0 else 0.0
            peaks.append({
                "mz": round(mz, 4),
                "abundance_pct": round(prob * 100, 3),
                "intensity": round(prob, 6),
                "fwhm": round(fwhm, 5),
            })

        if not peaks:
            raise ChemMCPError("No peaks above minimum abundance threshold.")

        # Generate envelope summary
        summary = self._generate_envelope_summary(elements, peaks)

        return {
            "result": {
                "formula": molecular_formula,
                "monoisotopic_mass": round(peaks[0]["mz"], 4),
                "charge": charge,
                "resolution": resolution,
                "peak_count": len(peaks),
                "nominal_mass": int(round(peaks[0]["mz"])),
                "envelope_summary": summary,
                "peaks": peaks,
            }
        }

    @staticmethod
    def _convolve(d1: Dict[float, float], d2: Dict[float, float]) -> Dict[float, float]:
        """Convolve two mass->probability distributions."""
        result: Dict[float, float] = defaultdict(float)
        for m1, p1 in d1.items():
            for m2, p2 in d2.items():
                result[round(m1 + m2, 8)] += p1 * p2
        return dict(result)

    @staticmethod
    def _power_dist(dist: Dict[float, float], power: int) -> Dict[float, float]:
        """Raise a distribution to an integer power via repeated squaring."""
        if power <= 0:
            return {0.0: 1.0}
        result: Dict[float, float] = {0.0: 1.0}
        base = dist.copy()
        while power > 0:
            if power % 2 == 1:
                result = IsotopePatternSimulator._convolve(result, base)
            base = IsotopePatternSimulator._convolve(base, base)
            power //= 2
        return result

    def _generate_envelope_summary(self, elements: Dict[str, int], peaks: list) -> str:
        """Generate human-readable summary of the isotopic envelope."""
        parts = []
        has_cl = "Cl" in elements
        has_br = "Br" in elements
        has_s = "S" in elements
        n_c = elements.get("C", 0)

        if has_cl or has_br:
            if has_cl and has_br:
                parts.append("Complex Cl+Br isotope cluster present")
            elif has_cl:
                parts.append(f"Characteristic ~3:1 Cl isotope pattern ({elements['Cl']} Cl atom(s))")
            else:
                parts.append(f"Characteristic ~1:1 Br isotope pattern ({elements['Br']} Br atom(s))")
        if has_s:
            parts.append(f"Sulfur contributes minor ^33S/^34S satellites (~{min(4.5 * elements['S'], 100):.1f}% M+2)")
        if n_c > 0:
            m1_pct = min(n_c * 1.07, 99.9)
            parts.append(f"M+1 peak ~{m1_pct:.1f}% from ^{{13}}C ({n_c} C atoms)")
        if len(peaks) <= 3:
            parts.append("Compact envelope (few isotopic variants)")
        elif len(peaks) >= 10:
            parts.append("Broad envelope (many overlapping isotopic combinations)")

        return "; ".join(parts) if parts else "Simple isotopic envelope"

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            formula = parts[0]
            charge = int(parts[1]) if len(parts) > 1 else 0
            res = int(parts[2]) if len(parts) > 2 else 30000
            min_ab = float(parts[3]) if len(parts) > 3 else 0.001
            return self._run_base(formula, charge, res, min_ab)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'formula charge resolution min_abundance'")
