"""
d电子在晶体场轨道中的排布工具
d-electron orbital distribution in crystal field (ASCII diagrams).
"""
import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DElectronConfiguration(BaseTool):
    """
    展示d电子在晶体场分裂轨道中的排布（ASCII图示）。
    支持八面体、四面体，以及高/低自旋选项。
    """
    __version__ = "0.1.0"
    name = "DElectronConfiguration"
    func_name = "d_electron_configuration"
    description = "Display d-electron orbital distribution in crystal field splitting with ASCII art diagrams. Shows t2g/eg (octahedral) or e/t2 (tetrahedral) filling for high-spin and low-spin configurations."
    implementation_description = "Generates detailed ASCII orbital filling diagrams based on CFT rules. For octahedral d4-d7, shows both spin states. Includes electron count per orbital set, total electrons, and spin multiplicity."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Crystal Field", "Electron Configuration", "Orbital Diagram", "CFT"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal ion, e.g., 'Fe2+', 'Co3+', 'Cr3+'."),
        ("geometry", "str", "octahedral", "'octahedral' or 'tetrahedral'."),
        ("field_strength", "str", "both", "'strong', 'weak', or 'both' (shows both spin states where applicable)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion geometry [field_strength]', e.g., 'Fe2+ octahedral both' or 'Mn2+ tetrahedral'."),
    ]

    output_sig = [
        ("metal_ion", "str", "Metal ion analyzed."),
        ("d_count", "int", "Total number of d electrons."),
        ("geometry", "str", "Coordination geometry."),
        ("orbital_diagram", "str", "ASCII/text orbital diagram with electron occupancy."),
        ("configuration_summary", "dict", "Summary of electron counts per orbital set, unpaired electrons, multiplicity."),
        ("spin_state_info", "str", "Description of which spin state(s) are possible and which is favored."),
        ("filling_sequence", "str", "Step-by-step electron filling order explanation."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Fe2+",
                "geometry": "octahedral",
                "field_strength": "both",
            },
            "text_input": {
                "query": "Fe2+ octahedral both"
            },
            "output": {
                "metal_ion": "Fe2+",
                "d_count": 6,
                "geometry": "Octahedral",
                "orbital_diagram": "=== HIGH-SPIN ===\n...\n=== LOW-SPIN ===\n...",
                "configuration_summary": {"high_spin": {"t2g": 4, "eg": 2, "unpaired": 4}, "low_spin": {"t2g": 6, "eg": 0, "unpaired": 0}},
                "spin_state_info": "Weak field → high-spin; Strong field → low-spin",
                "filling_sequence": "e⁻¹→t2g, e⁻²→t2g, ..., e⁻⁴→t2g, e⁻⁵→t2g/eg(decision point), ...",
            }
        },
        {
            "code_input": {
                "metal_ion": "Cr3+",
                "geometry": "octahedral",
                "field_strength": "weak",
            },
            "text_input": {
                "query": "Cr3+ octahedral"
            },
            "output": {
                "metal_ion": "Cr3+",
                "d_count": 3,
                "geometry": "Octahedral",
                "orbital_diagram": "ASCII diagram showing t2g³ eg⁰",
                "configuration_summary": {"t2g": 3, "eg": 0, "unpaired": 3, "multiplicity": "4 (quartet)"},
                "spin_state_info": "Only one configuration exists for d3 (no spin crossover)",
                "filling_sequence": "e⁻¹→t2g(↑), e⁻²→t2g(↑), e⁻³→t2g(↑)",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize d-count database."""
        self._d_counts = {
            "ti3+": 1, "v3+": 2, "cr3+": 3, "cr2+": 4,
            "mn3+": 4, "mn2+": 5, "fe3+": 5, "fe2+": 6,
            "co3+": 6, "co2+": 7, "ni2+": 8, "ni3+": 7,
            "cu2+": 9, "cu+": 10, "zn2+": 10, "ag+": 10,
            "pt2+": 8, "pd2+": 8, "au3+": 8,
            "v2+": 3, "ti2+": 2,
        }

    def _get_d_count(self, metal_ion: str) -> int:
        key = metal_ion.lower().replace(" ", "")
        if key not in self._d_counts:
            raise ChemMCPError(
                f"Unknown metal ion '{metal_ion}'. Known: "
                f"{', '.join(sorted(set(k.upper() for k in self._d_counts.keys())))}"
            )
        return self._d_counts[key]

    def _make_oct_diagram(self, n: int, strong: bool = None) -> str:
        """Generate ASCII orbital diagram for octahedral field."""
        lines = []

        # Determine configuration
        if strong is True or (strong is None and n not in (4, 5, 6, 7)):
            # Low-spin or unique config
            configs = {
                0: [(0, "_")], 1: [(1, "↑")], 2: [(2, "↑↑")], 3: [(3, "↑↑↑")],
                4: [(4, "↓↓↓↓")], 5: [(5, "↓↓↓↓↓")], 6: [(6, "↓↓↓↓↓↓")],
                7: [(6, "↓↓↓↓↓↓"), (1, "↑")], 8: [(6, "↓↓↓↓↓↓"), (2, "↑↑")],
                9: [(6, "↓↓↓↓↓↓"), (3, "↑↑↑")], 10: [(6, "↓↓↓↓↓↓"), (4, "↓↓↓↓")],
            }
            if strong is False:
                # High-spin forced
                configs = {
                    0: [(0, "_")], 1: [(1, "↑")], 2: [(2, "↑↑")], 3: [(3, "↑↑↑")],
                    4: [(3, "↑↑↑"), (1, "↑")], 5: [(3, "↑↑↑"), (2, "↑↑")],
                    6: [(4, "↑↑↑↑"), (2, "↑↑")], 7: [(5, "↑↑↑↑↑"), (2, "↑↑")],
                    8: [(6, "↑↑↑↑↑↑"), (2, "↑↑")], 9: [(6, "↑↑↑↑↑↑"), (3, "↑↑↑")],
                    10: [(6, "↓↓↓↓↓↓"), (4, "↓↓↓↓")],
                }
            cfg = configs[n]
        elif strong is False:
            cfg = {
                0: [(0, "_")], 1: [(1, "↑")], 2: [(2, "↑↑")], 3: [(3, "↑↑↑")],
                4: [(3, "↑↑↑"), (1, "↑")], 5: [(3, "↑↑↑"), (2, "↑↑")],
                6: [(4, "↑↑↑↑"), (2, "↑↑")], 7: [(5, "↑↑↑↑↑"), (2, "↑↑")],
                8: [(6, "↑↑↑↑↑↑"), (2, "↑↑")], 9: [(6, "↑↑↑↑↑↑"), (3, "↑↑↑")],
                10: [(6, "↓↓↓↓↓↓"), (4, "↓↓↓↓")],
            }[n]
        else:
            # Both: return dual diagram
            return self._make_oct_dual(n)

        # Build single diagram
        t2g_count = cfg[0][0] if cfg else 0
        eg_count = cfg[1][0] if len(cfg) > 1 else 0

        lines.append("     Energy")
        lines.append("       ↑")
        lines.append("   ┌───────────┐")
        lines.append("   │  eg (dx²-y², dz²)  │")
        eg_occ = self._format_orbitals(eg_count, 2, strong=True if n in (4,5,6,7) and strong else False)
        lines.append(f"   │   {eg_occ}   │")
        lines.append("   ├───────────┤  Δo")
        lines.append("   │  t2g (dxy, dxz, dyz)│")
        t2g_occ = self._format_orbitals(t2g_count, 3, strong=(strong == True))
        lines.append(f"   │   {t2g_occ}   │")
        lines.append("   └───────────┘")

        label = "LOW-SPIN" if strong else ("HIGH-SPIN" if strong is False else "")
        if label:
            lines.insert(0, f"  === {label} CONFIGURATION (d^{n}) ===")
        else:
            lines.insert(0, f"  === OCTAHEDRAL d^{n} ===")

        return "\n".join(lines)

    def _make_oct_dual(self, n: int) -> str:
        """Generate side-by-side comparison of HS and LS for d4-d7."""
        hs_configs = {
            4: (3, 1), 5: (3, 2), 6: (4, 2), 7: (5, 2),
        }
        ls_configs = {
            4: (4, 0), 5: (5, 0), 6: (6, 0), 7: (6, 1),
        }

        hs_t2g, hs_eg = hs_configs[n]
        ls_t2g, ls_eg = ls_configs[n]

        lines = []
        lines.append(f"  ╔══════════════════════════════════════════════════════╗")
        lines.append(f"  ║         OCTAHEDRAL d^{n} — SPIN STATE COMPARISON        ║")
        lines.append(f"  ╠════════════════════╦══════════════════════════════════╣")
        lines.append(f"  ║   HIGH-SPIN (weak)  ║      LOW-SPIN (strong)           ║")
        lines.append(f"  ╠════════════════════╬══════════════════════════════════╣")

        # Header row
        lines.append(f"  ║       Energy ↑      ║       Energy ↑                  ║")
        lines.append(f"  ║  ┌────────────┐    ║  ┌────────────┐                 ║")
        lines.append(f"  ║  │ eg          │    ║  │ eg          │                 ║")
        hs_eg_s = self._format_orbitals(hs_eg, 2, False)
        ls_eg_s = self._format_orbitals(ls_eg, 2, True)
        lines.append(f"  ║  │ {hs_eg_s:^12s} │    ║  │ {ls_eg_s:^12s} │                 ║")
        lines.append(f"  ║  ├────────────┤ Δo  ║  ├────────────┤ Δo              ║")
        lines.append(f"  ║  │ t2g         │    ║  │ t2g         │                 ║")
        hs_t2g_s = self._format_orbitals(hs_t2g, 3, False)
        ls_t2g_s = self._format_orbitals(ls_t2g, 3, True)
        lines.append(f"  ║  │ {hs_t2g_s:^12s} │    ║  │ {ls_t2g_s:^12s} │                 ║")
        lines.append(f"  ║  └────────────┘    ║  └────────────┘                 ║")

        hs_unp = hs_t2g % 3 + min(hs_eg, 2) if hs_t2g <= 3 else (hs_t2g - 3) + hs_eg
        # Recalculate properly
        hs_unp = self._count_unpaired_hs(n)
        ls_unp = self._count_unpaired_ls(n)

        lines.append(f"  ╠════════════════════╬══════════════════════════════════╣")
        lines.append(f"  ║ t2g^{hs_t2g} eg^{hs_eg}  {hs_unp}e⁻ unpaired ║ t2g^{ls_t2g} eg^{ls_eg}  {ls_unp}e⁻ unpaired       ║")
        lines.append(f"  ║ μso={self._mu(hs_unp):.2f} BM             ║ μso={self._mu(ls_unp):.2f} BM                    ║")
        lines.append(f"  ╚════════════════════╩══════════════════════════════════╝")
        return "\n".join(lines)

    def _format_orbitals(self, count: int, max_orbs: int, paired_mode: bool = False, strong: bool = False) -> str:
        """Format orbital occupancy string."""
        # Use strong as alias for paired_mode
        if strong:
            paired_mode = strong
        if count == 0:
            return "  ".join(["__"] * max_orbs)

        if paired_mode:
            # All paired first
            pairs = count // 2
            singles = count % 2
            result = []
            for i in range(max_orbs):
                if i < pairs:
                    result.append("↓↑")
                elif i < pairs + singles:
                    result.append("↑ ")
                else:
                    result.append("__")
            return " ".join(result)
        else:
            # Fill one per orbital first (Hund's rule)
            result = []
            remaining = count
            # First pass: one electron each
            for i in range(max_orbs):
                if remaining > 0:
                    result.append("↑ ")
                    remaining -= 1
                else:
                    result.append("__")
            # Second pass: pair up
            for i in range(max_orbs):
                if remaining > 0:
                    result[i] = "↓↑"
                    remaining -= 1
            return " ".join(result)

    def _count_unpaired_hs(self, n: int) -> int:
        """Count unpaired electrons in high-spin octahedral."""
        hs_data = {0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:4, 7:3, 8:2, 9:1, 10:0}
        return hs_data.get(n, 0)

    def _count_unpaired_ls(self, n: int) -> int:
        """Count unpaired in low-spin octahedral."""
        ls_data = {0:0, 1:1, 2:2, 3:3, 4:2, 5:1, 6:0, 7:1, 8:2, 9:1, 10:0}
        return ls_data.get(n, 0)

    def _mu(self, n: int) -> float:
        import math
        if n == 0:
            return 0.0
        return round(math.sqrt(n * (n + 2)), 2)

    def _make_tet_diagram(self, n: int) -> str:
        """Generate ASCII orbital diagram for tetrahedral field."""
        # Tetrahedral: e (lower), t2 (higher), always high-spin
        tet_fill = {
            0: (0, 0), 1: (1, 0), 2: (2, 0), 3: (2, 1),
            4: (2, 2), 5: (3, 2), 6: (4, 2), 7: (4, 3),
            8: (4, 4), 9: (4, 4), 10: (4, 4),
        }
        e_c, t2_c = tet_fill[n]

        lines = []
        lines.append(f"  === TETRAHEDRAL d^{n} (always high-spin) ===")
        lines.append("       Energy ↑")
        lines.append("   ┌────────────────┐")
        lines.append("   │  t2 (dxy,dxz,dyz)│")
        t2_s = self._format_orbitals(t2_c, 3, False)
        lines.append(f"   │  {t2_s}")
        lines.append("   ├────────────────┤  Δt")
        lines.append("   │  e  (dz²,dx²-y²) │")
        e_s = self._format_orbitals(e_c, 2, False)
        lines.append(f"   │  {e_s}")
        lines.append("   └────────────────┘")

        unpaired = self._count_unpaired_tet(n)
        lines.append(f"  e^{e_c} t2^{t2_c}, {unpaired} unpaired e⁻, μso = {self._mu(unpaired):.2f} BM")
        return "\n".join(lines)

    def _count_unpaired_tet(self, n: int) -> int:
        tet_u = {0:0, 1:1, 2:2, 3:3, 4:4, 5:3, 6:2, 7:3, 8:2, 9:1, 10:0}
        return tet_u.get(n, 0)

    def _run_base(self, metal_ion: str, geometry: str = "octahedral",
                  field_strength: str = "both") -> dict:
        """Generate d-electron orbital diagram."""
        geo = geometry.lower()
        if geo not in ("octahedral", "tetrahedral"):
            raise ChemMCPError("Geometry must be 'octahedral' or 'tetrahedral'.")

        field = field_strength.lower()
        n = self._get_d_count(metal_ion)

        if geo == "octahedral":
            if field == "both":
                diagram = self._make_oct_diagram(n, None)
                show_both = n in (4, 5, 6, 7)
                spin_info = (
                    f"For d{n}: weak-field ligands → high-spin; strong-field ligands → low-spin"
                    if show_both else
                    f"d{n} has only one possible arrangement (no spin crossover)"
                )
                summary = {
                    "high_spin": {
                        "t2g": {4:3, 5:3, 6:4, 7:5}.get(n, min(n, 3)),
                        "eg": {4:1, 5:2, 6:2, 7:2}.get(n, max(n-3, 0)),
                        "unpaired": self._count_unpaired_hs(n),
                        "multiplicity": self._count_unpaired_hs(n) + 1,
                    },
                    "low_spin": {
                        "t2g": {4:4, 5:5, 6:6, 7:6}.get(n, min(n, 6)),
                        "eg": {4:0, 5:0, 6:0, 7:1}.get(n, max(n-6, 0)),
                        "unpaired": self._count_unpaired_ls(n),
                        "multiplicity": self._count_unpaired_ls(n) + 1,
                    },
                } if show_both else {
                    "only": {
                        "t2g": min(n, 3) if n <= 3 else ({8:6, 9:6, 10:6}.get(n, n)),
                        "eg": 0 if n <= 3 else ({8:2, 9:3, 10:4}.get(n, n-3)),
                        "unpaired": {0:0, 1:1, 2:2, 3:3, 8:2, 9:1, 10:0}.get(n, n),
                        "multiplicity": {0:1, 1:2, 2:3, 3:4, 8:3, 9:2, 10:1}.get(n, n+1),
                    }
                }
            elif field == "strong":
                diagram = self._make_oct_diagram(n, True)
                t2g_val = min(n, 6) if n in (4,5,6,7) else (min(n, 3) if n <= 3 else {8:6, 9:6, 10:6}[n])
                eg_val = 0 if n in (4,5,6) else (1 if n == 7 else ({8:2, 9:3, 10:4}[n]))
                spin_info = f"Low-spin configuration enforced (strong field)"
                summary = {"t2g": t2g_val, "eg": eg_val, "unpaired": self._count_unpaired_ls(n), "multiplicity": self._count_unpaired_ls(n)+1}
            else:
                diagram = self._make_oct_diagram(n, False)
                t2g_v = min(n, 3) if n <= 3 else {4:3, 5:3, 6:4, 7:5, 8:6, 9:6, 10:6}[n]
                eg_v = 0 if n <= 3 else {4:1, 5:2, 6:2, 7:2, 8:2, 9:3, 10:4}[n]
                spin_info = f"High-spin configuration enforced (weak field)"
                summary = {"t2g": t2g_v, "eg": eg_v, "unpaired": self._count_unpaired_hs(n), "multiplicity": self._count_unpaired_hs(n)+1}

            filling = self._filling_sequence_oct(n)
        else:
            diagram = self._make_tet_diagram(n)
            spin_info = "Tetrahedral complexes are always high-spin (Δt ≈ 4/9 Δo is small)"
            tet_f = {0:(0,0), 1:(1,0), 2:(2,0), 3:(2,1), 4:(2,2), 5:(3,2), 6:(4,2), 7:(4,3), 8:(4,4), 9:(4,4), 10:(4,4)}
            e_c, t2_c = tet_f[n]
            summary = {"e": e_c, "t2": t2_c, "unpaired": self._count_unpaired_tet(n), "multiplicity": self._count_unpaired_tet(n)+1}
            filling = self._filling_sequence_tet(n)

        logger.info(f"d-electron config: {metal_ion} d^{n} {geo} {field}")

        return {
            "metal_ion": metal_ion,
            "d_count": n,
            "geometry": geo.capitalize(),
            "orbital_diagram": diagram,
            "configuration_summary": summary,
            "spin_state_info": spin_info,
            "filling_sequence": filling,
        }

    def _filling_sequence_oct(self, n: int) -> str:
        steps = []
        for i in range(1, n + 1):
            if i <= 3:
                steps.append(f"e⁻^{i} → t2g (Hund's rule, parallel spins)")
            elif i == 4:
                steps.append(f"e⁻⁴ → decision point: t2g (LS, pair up) vs eg (HS, Hund's)")
            elif i == 5:
                steps.append(f"e⁻⁵ → decision point: t2g (LS) vs eg (HS)")
            elif i == 6:
                steps.append(f"e⁻⁶ → decision point: t2g (LS, full) vs eg (HS)")
            elif i == 7:
                steps.append(f"e⁻⁷ → t2g (LS, full) + eg¹ OR t2g⁵ + eg² (HS)")
            else:
                orb = "t2g" if i <= (6 + (n-6)//2 * 0 + 999) else "eg"
                actual = "t2g" if i <= 6 else "eg"
                steps.append(f"e⁻^{i} → {actual}")
        return "; ".join(steps)

    def _filling_sequence_tet(self, n: int) -> str:
        steps = []
        for i in range(1, n + 1):
            if i <= 2:
                steps.append(f"e⁻^{i} → e (lower in tetrahedral)")
            elif i <= 5:
                steps.append(f"e⁻^{i} → t2 (Hund's rule)")
            else:
                actual = "e" if i <= 4 else "t2"
                steps.append(f"e⁻^{i} → {actual} (always HS)")
        return "; ".join(steps)

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Format: 'metal_ion [geometry] [field_strength]'. Example: 'Fe2+ octahedral both'")
        metal = parts[0]
        geo = parts[1] if len(parts) > 1 else "octahedral"
        field = parts[2] if len(parts) > 2 else "both"
        return self._run_base(metal, geo, field)
