import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SplittingPatternExplainer(BaseTool):
    """
    NMR 裂分峰型解释工具。
    基于 n+1 规则和 Pascal 三角形解释峰裂分模式。
    """
    __version__ = "0.1.0"
    name = "SplittingPatternExplainer"
    func_name = "explain_splitting_pattern"
    description = "Explain NMR peak splitting patterns using the n+1 rule and Pascal's triangle. Returns pattern name, number of peaks, intensity ratios, and visual diagram."
    implementation_description = "Implements the n+1 rule for NMR first-order splitting analysis. Uses Pascal's triangle for intensity ratios. Covers singlet through multiplet patterns up to n=15 neighbors, including common abbreviations (s, d, t, q, quint, sext, sept) and complex cases."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["NMR", "Splitting Pattern", "n+1 Rule", "Pascal Triangle", "Spectroscopy", "Multiplicity"]
    required_envs = []

    code_input_sig = [
        ("n_neighbors", "int", "N/A", "Number of equivalent neighboring protons (n)."),
        ("J_coupling_hz", "float", "7.0", "Coupling constant J in Hz (for peak spacing; default 7.0 Hz typical)."),
        ("show_diagram", "bool", "True", "Whether to include ASCII art diagram of the splitting pattern."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'n_neighbors [J_hz]'. Example: '3 7.2'"),
    ]

    output_sig = [
        ("pattern_info", "dict", "Complete splitting pattern information including name, peaks, intensities, diagram, and notes."),
    ]

    examples = [
        {
            "code_input": {"n_neighbors": 3, "J_coupling_hz": 7.0, "show_diagram": True},
            "text_input": {"input_params": "3"},
            "output": {
                "pattern_info": {
                    "pattern_name": "quartet",
                    "n_peaks": 4,
                    "intensity_ratio": [1, 3, 3, 1],
                }
            },
        },
        {
            "code_input": {"n_neighbors": 6, "J_coupling_hz": 7.0, "show_diagram": True},
            "text_input": {"input_params": "6"},
            "output": {
                "pattern_info": {
                    "pattern_name": "septet",
                    "n_peaks": 7,
                    "intensity_ratio": [1, 6, 15, 20, 15, 6, 1],
                }
            },
        },
    ]

    # Pascal's triangle rows (precomputed up to n=20)
    # Row n contains coefficients for (a+b)^n
    _PASCAL_TRIANGLE = [
        [1],                          # n=0: singlet
        [1, 1],                        # n=1: doublet
        [1, 2, 1],                     # n=2: triplet
        [1, 3, 3, 1],                  # n=3: quartet
        [1, 4, 6, 4, 1],               # n=4: quintet
        [1, 5, 10, 10, 5, 1],          # n=5: sextet
        [1, 6, 15, 20, 15, 6, 1],      # n=6: septet
        [1, 7, 21, 35, 35, 21, 7, 1],   # n=7: octet
        [1, 8, 28, 56, 70, 56, 28, 8, 1],  # n=8: nonet
        [1, 9, 36, 84, 126, 126, 84, 36, 9, 1],  # n=9: decet
        [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1],  # n=10
        [1, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1],  # n=11
        [1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1],  # n=12
        [1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1],  # n=13
        [1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364, 91, 14, 1],  # n=14
        [1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1],  # n=15
        [1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560, 120, 16, 1],  # n=16
        [1, 17, 136, 680, 2380, 6188, 12376, 19448, 24310, 24310, 19448, 12376, 6188, 2380, 680, 136, 17, 1],  # n=17
        [1, 18, 153, 816, 3060, 8568, 18564, 31824, 43758, 48620, 43758, 31824, 18564, 8568, 3060, 816, 153, 18, 1],  # n=18
        [1, 19, 171, 969, 3876, 11628, 27132, 50388, 75582, 92378, 92378, 75582, 50388, 27132, 11628, 3876, 969, 171, 19, 1],  # n=19
        [1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756, 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1],  # n=20
    ]

    # Abbreviation mapping
    _PATTERN_NAMES = {
        0: ("singlet", "s"),
        1: ("doublet", "d"),
        2: ("triplet", "t"),
        3: ("quartet", "q"),
        4: ("quintet", "quint" if False else "p"),  # p for pentet
        5: ("sextet", "sex" if False else "sext"),
        6: ("septet", "sep"),
        7: ("octet", "o"),
        8: ("nonet", "non"),
        9: ("decet", "dec"),
        10: ("undecet", "und"),
        11: ("dodecet", "ddo"),
        12: ("tridecet", "tri"),
        13: ("tetradecet", "tt"),
        14: ("pentadecet", "pd"),
        15: ("hexadecet", "hd"),
        16: ("heptadecet", "hp"),
        17: ("octadecet", "od"),
        18: ("nonadecet", "nd"),
        19: ("viginticet", "vi"),
        20: ("unviginticet", "uv"),
    }

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _generate_ascii_diagram(self, ratios: list, j_hz: float, max_width: int = 50) -> str:
        """Generate ASCII art of splitting pattern."""
        n_peaks = len(ratios)
        total = sum(ratios)

        if total == 0:
            return ""

        lines = []
        lines.append(f"\n  {'Peak':^8} | {'Intensity':^10} | {'Relative':^10} | {'Diagram'}")
        lines.append(f"  {'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*max_width}")

        max_ratio = max(ratios)
        for i, r in enumerate(ratios):
            rel_pct = round(r / total * 100, 1)
            bar_len = int(r / max_ratio * max_width)
            bar = "█" * bar_len
            offset = (j_hz * i) if n_peaks > 1 else 0
            offset_str = f"+{offset:.1f}Hz" if offset > 0 else "center"

            lines.append(
                f"  {i+1:^8} | {r:^10} | {rel_pct:>8}%  | {bar}"
            )

        # Draw simplified stick spectrum below
        lines.append(f"\n  Stick Spectrum (spacing = J = {j_hz} Hz):")
        lines.append("  " + "-" * (max_width + 10))

        # Create a simple line-based spectrum
        spectrum_line = []
        scale = max_width / (n_peaks - 1) if n_peaks > 1 else max_width
        for i in range(n_peaks):
            pos = int(i * scale) if n_peaks > 1 else int(max_width / 2)
            height = int((ratios[i] / max_ratio) * 6)
            spectrum_line.append((pos, height))

        for row in range(6, 0, -1):
            line = "  "
            for pos, h in spectrum_line:
                bar_start = 0
                for p2, h2 in spectrum_line:
                    if p2 == pos and h >= row:
                        line += " " * (pos - len(line) + 1) + "│"
                        break
                    elif p2 == pos and h < row:
                        line += " " * (pos - len(line) + 1) + "·"
                        break
                else:
                    continue
            lines.append(line.rstrip())

        return "\n".join(lines)

    def _get_complex_pattern_notes(self, n: int) -> str:
        """Get additional notes for complex splitting patterns."""
        notes = []

        if n == 0:
            notes.append("No neighboring protons → no coupling → single sharp peak.")
            notes.append("Common for: isolated CH₃ (no adjacent H), OH/NH₂ (exchange-broadened), quaternary C-H.")
        elif n == 1:
            notes.append("One neighboring proton splits signal into two equal peaks.")
            notes.append("Common in: -CH-CH₃ groups, trans/cis olefinic protons.")
        elif n == 2:
            notes.append("Two equivalent neighbors give 1:2:3 triplet pattern.")
            notes.append("Classic example: CH₃-CH₂-X (CH₃ appears as triplet).")
        elif n == 3:
            notes.append("Three equivalent neighbors give 1:3:3:1 quartet.")
            notes.append("Classic example: CH₃-CH₂-X (CH₂ appears as quartet from CH₃ coupling).")
        elif n >= 10:
            notes.append("Large n values approximate a binomial distribution envelope.")
            notes.append("At very large n, the pattern approaches a Gaussian/Poisson shape.")
            notes.append("In practice, such patterns often appear as broad multiplets due to overlapping couplings.")

        # General rules
        notes.append("")
        notes.append("KEY RULES:")
        notes.append("• n+1 rule applies to FIRST-ORDER systems (Δν >> J)")
        notes.append("• Only EQUIVALENT protons are counted together")
        notes.append("• Non-equivalent neighbors produce more complex patterns (doublet of doublets, etc.)")
        notes.append("• Coupling is mutual: if Hᵃ splits Hᵇ into n+1 peaks, Hᵇ also splits Hᵃ")

        return "\n".join(notes)

    def _run_base(self, n_neighbors: int, J_coupling_hz: float = 7.0,
                  show_diagram: bool = True) -> dict:
        """
        Explain splitting pattern.

        Args:
            n_neighbors: Number of equivalent neighboring protons
            J_coupling_hz: Coupling constant in Hz
            show_diagram: Whether to generate ASCII diagram

        Returns:
            Dict with full pattern information
        """
        if not isinstance(n_neighbors, int) or n_neighbors < 0:
            raise ChemMCPError("n_neighbors must be a non-negative integer.")

        if n_neighbors > 20:
            raise ChemMCPError("n_neighbors must be ≤ 20. For larger values, use approximation methods.")

        if J_coupling_hz <= 0:
            raise ChemMCPError("J_coupling_hz must be positive.")

        n = n_neighbors
        n_peaks = n + 1
        ratios = self._PASCAL_TRIANGLE[n]
        full_name, abbrev = self._PATTERN_NAMES[n]

        # Calculate peak positions (relative to center)
        positions = [(i - n / 2) * J_coupling_hz for i in range(n_peaks)]

        # Total relative intensity
        total_intensity = sum(ratios)
        normalized = [round(r / total_intensity, 4) for r in ratios]

        result = {
            "pattern_info": {
                "n_neighbors": n,
                "n_peaks": n_peaks,
                "pattern_name": full_name,
                "abbreviation": abbrev,
                "rule_applied": f"n+1 rule: {n} equivalent neighbor(s) → {n_peaks} peak(s)",
                "intensity_ratio": ratios,
                "normalized_intensities": normalized,
                "total_coefficient_sum": total_intensity,
                "peak_positions_hz": [round(p, 2) for p in positions],
                "coupling_constant_J_hz": J_coupling_hz,
                "peak_spacing_hz": J_coupling_hz,
                "total_pattern_width_hz": round(n * J_coupling_hz, 2),
                "pascal_row": f"Pascal triangle row {n}: ({', '.join(map(str, ratios))})",
                "notes": self._get_complex_pattern_notes(n),
            }
        }

        if show_diagram:
            result["pattern_info"]["ascii_diagram"] = self._generate_ascii_diagram(ratios, J_coupling_hz)

        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'n_neighbors [J_hz]'")

        try:
            n = int(parts[0])
        except ValueError:
            raise ChemMCPError(f"Invalid n_neighbors: '{parts[0]}' must be an integer")

        j_hz = float(parts[1]) if len(parts) > 1 else 7.0

        return self._run_base(n, j_hz)
