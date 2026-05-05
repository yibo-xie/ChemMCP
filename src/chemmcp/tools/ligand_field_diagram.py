import logging
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LigandFieldDiagram(BaseTool):
    """
    生成配体场能级图（Crystal Field / Ligand Field Theory）。
    支持八面体(Oh)、四面体(Td)、平面正方形(D4h)几何构型，
    以及不同配体强度（强场/弱场）下的 d 轨道分裂。
    """
    __version__ = "0.1.0"
    name = "LigandFieldDiagram"
    func_name = "generate_ligand_field_diagram"
    description = "Generate ligand/crystal field splitting diagrams for octahedral, tetrahedral, and square planar geometries with d-orbital energy levels."
    implementation_description = "Uses crystal field theory (CFT) data to compute relative d-orbital energies for Oh/Td/D4h geometries, including Δo, Δt values, electron filling orders for weak/strong field cases, and spectrochemical series-based Δ estimates."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Crystal Field", "Ligand Field", "d-orbitals", "Coordination Chemistry", "Energy Diagram"]
    required_envs = []

    code_input_sig = [
        ("geometry", "str", "N/A", "Geometry: 'octahedral', 'tetrahedral', or 'square_planar'."),
        ("d_count", "int", "N/A", "Number of d electrons (0-10)."),
        ("field_strength", "str", "weak", "Field strength: 'weak' (high-spin) or 'strong' (low-spin)."),
        ("ligand", "str", "H2O", "Ligand name (used to estimate Δ from spectrochemical series)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'geometry d_count field_strength [ligand]'. Example: 'octahedral 6 strong NH3'."),
    ]

    output_sig = [
        ("diagram_data", "dict", "Complete diagram data: orbital energies, electron configuration, CFSE, spin state, magnetic moment."),
        ("ascii_diagram", "str", "ASCII art representation of the energy level diagram."),
        ("explanation", "str", "Detailed explanation of the splitting pattern and electron filling."),
    ]

    examples = [
        {
            "code_input": {
                "geometry": "octahedral",
                "d_count": 6,
                "field_strength": "strong",
                "ligand": "NH3",
            },
            "text_input": {
                "input_params": "octahedral 6 strong NH3"
            },
            "output": {
                "diagram_data": {
                    "geometry": "octahedral",
                    "delta_o_dq": 1.0,
                    "orbitals": {
                        "eg": {"energy_dq": 0.6, "orbitals": ["dx2-y2", "dz2"], "electrons": 0},
                        "t2g": {"energy_dq": -0.4, "orbitals": ["dxy", "dxz", "dyz"], "electrons": 6},
                    },
                    "electron_config": "(t2g)^6(eg)^0",
                    "spin_state": "low_spin",
                    "unpaired_electrons": 0,
                    "cfse_dq": -2.4,
                    "magnetic_moment_bm": 0.0,
                },
                "ascii_diagram": "      eg (dx²-y², dz²)  ━━━━━ +0.6Δo ━━━━\n                              ↑↓\n                               \n      t2g (dxy,dxz,dyz) ━━━ -0.4Δo ━━━\n                          ↑↓   ↑↓   ↑↓",
                "explanation": "In an octahedral strong-field complex with d^6 (e.g., [Co(NH3)6]3+), the large Δo forces all 6 electrons to pair in the lower t2g set, giving a low-spin diamagnetic configuration with CFSE = -2.4Δo.",
            }
        },
    ]

    # Spectrochemical series (approximate Δo values relative to H2O ≈ 1.0)
    SPECTROCHEMICAL_SERIES = {
        "I-": 0.7, "Br-": 0.76, "SCN-": 0.77, "Cl-": 0.78, "N3-": 0.83,
        "F-": 0.9, "urea": 0.91, "OH-": 0.83, "ox2-": 0.97, "H2O": 1.0,
        "NCS-": 1.02, "CH3CN": 1.0, "py": 1.07, "NH3": 1.23, "en": 1.28,
        "bipy": 1.28, "phen": 1.34, "NO2-": 1.4, "PPh3": 1.3,
        "CN-": 1.7, "CO": 1.9,
    }

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize orbital data."""
        pass

    def _get_delta_factor(self, ligand: str) -> float:
        """Get Δ factor from spectrochemical series."""
        ligand_key = None
        for k in self.SPECTROCHEMICAL_SERIES:
            if k.lower() == ligand.lower().strip():
                ligand_key = k
                break
        if ligand_key:
            return self.SPECTROCHEMICAL_SERIES[ligand_key]
        # Default for unknown ligands
        return 1.0

    def _fill_orbitals_oh(self, d_n: int, field: str) -> dict:
        """Fill d orbitals for octahedral geometry."""
        is_strong = field.lower() == "strong"
        # In Oh: t2g (3 orbitals, -0.4Δ each), eg (2 orbitals, +0.6Δ each)
        # Pairing energy P is typically ~0.8Δo for the boundary

        t2g_electrons = 0
        eg_electrons = 0
        unpaired = 0
        spin_type = ""

        if is_strong:
            # Strong field: fill t2g completely first
            t2g_electrons = min(d_n, 6)
            eg_electrons = max(0, d_n - 6)
            unpaired = max(0, d_n % 5) if d_n <= 6 else (d_n - 6)
            # Actually for strong field:
            # d0-d3: all unpaired in t2g; d4-d6: paired in t2g; d7: t2g^6 eg^1; d8: t2g^6 eg^2
            if d_n <= 3:
                unpaired = d_n
            elif d_n <= 6:
                unpaired = 0
            elif d_n == 7:
                unpaired = 1
            elif d_n == 8:
                unpaired = 2
            elif d_n == 9:
                unpaired = 1
            else:
                unpaired = 0
            spin_type = "low_spin" if d_n in (4, 5, 6, 7) else ("high_spin" if d_n <= 3 else "")
        else:
            # Weak field (high spin): follow Hund's rule
            # Fill t2g with one e- each first (max 3), then eg with one e- each (max 2), then pair
            remaining = d_n
            # First round: one electron per orbital (t2g: 3, eg: 2)
            if remaining > 0:
                t2g_first = min(remaining, 3)
                t2g_electrons += t2g_first
                unpaired += t2g_first
                remaining -= t2g_first
            if remaining > 0:
                eg_first = min(remaining, 2)
                eg_electrons += eg_first
                unpaired += eg_first
                remaining -= eg_first
            # Second round: pair up (decrement unpaired as we fill already-occupied orbitals)
            if remaining > 0:
                t2g_second = min(remaining, 6 - t2g_electrons)
                unpaired -= t2g_second
                t2g_electrons += t2g_second
                remaining -= t2g_second
            if remaining > 0:
                eg_second = min(remaining, 4 - eg_electrons)
                unpaired -= eg_second
                eg_electrons += eg_second
            spin_type = "high_spin"

        cfse = t2g_electrons * (-0.4) + eg_electrons * 0.6  # in units of Δo

        return {
            "t2g_e": t2g_electrons,
            "eg_e": eg_electrons,
            "unpaired": unpaired,
            "spin": spin_type,
            "cfse_dq": round(cfse, 4),
            "config": f"(t2g)^{t2g_electrons}(eg)^{eg_electrons}",
        }

    def _fill_orbitals_td(self, d_n: int) -> dict:
        """Fill d orbitals for tetrahedral geometry (always high-spin)."""
        # Td: e (2 orbitals, -0.6Δt), t2 (3 orbitals, +0.4Δt); Δt ≈ 4/9 Δo
        e_electrons = 0
        t2_electrons = 0
        unpaired = 0
        remaining = d_n

        # Always high-spin in Td (Δt is small)
        # Fill e first (lower in Td!), then t2
        if remaining > 0:
            e_first = min(remaining, 2)
            e_electrons += e_first
            unpaired += e_first
            remaining -= e_first
        if remaining > 0:
            t2_first = min(remaining, 3)
            t2_electrons += t2_first
            unpaired += t2_first
            remaining -= t2_first
        # Pairing
        if remaining > 0:
            t2_second = min(remaining, 6 - t2_electrons)
            t2_electrons += t2_second
            remaining -= t2_second
        if remaining > 0:
            e_second = min(remaining, 4 - e_electrons)
            e_electrons += e_second

        cfse = e_electrons * (-0.6) + t2_electrons * 0.4  # in units of Δt

        return {
            "e_e": e_electrons,
            "t2_e": t2_electrons,
            "unpaired": unpaired,
            "spin": "high_spin",
            "cfse_dqt": round(cfse, 4),
            "config": f"(e)^{e_electrons}(t2)^{t2_electrons}",
        }

    def _fill_orbitals_sq(self, d_n: int) -> dict:
        """Fill d orbitals for square planar (D4h) geometry."""
        # D4h order (lowest to highest): dxz/dyz < dxy < dz2 < dx2-y2
        # Energies relative to barycenter (in units of Δo):
        # dxz,dyz: -4.28Dq; dxy: -2.28Dq; dz2: 0.86Dq; dx2-y2: 12.28Dq
        # Simplified: most d^8 are square planar (low-spin type)
        # For simplicity, we use a simplified model
        orbitals_order = [("dxz", -4.28), ("dyz", -4.28), ("dxy", -2.28),
                          ("dz2", 0.86), ("dx2_y2", 12.28)]
        electrons_per_orb = [0] * 5
        unpaired = 0
        remaining = d_n

        # Strong-field-like filling (square planar implies strong field)
        for i, (name, energy) in enumerate(orbitals_order):
            if remaining <= 0:
                break
            electrons_per_orb[i] = min(remaining, 2)
            if electrons_per_orb[i] == 1:
                unpaired += 1
            remaining -= electrons_per_orb[i]

        cfse = sum(electrons_per_orb[i] * orbitals_order[i][1] for i in range(5))

        config_parts = []
        orbital_names = ["dxz", "dyz", "dxy", "dz2", "dx2_y2"]
        for i, n in enumerate(electrons_per_orb):
            if n > 0:
                config_parts.append(f"{orbital_names[i]}^{n}")

        return {
            "electrons": {orbital_names[i]: electrons_per_orb[i] for i in range(5)},
            "unpaired": unpaired,
            "spin": "low_spin" if d_n >= 7 else "unknown",
            "cfse_dq": round(cfse, 4),
            "config": "(" + ")(".join(config_parts) + ")" if config_parts else "()",
        }

    def _make_ascii_oh(self, data: dict, d_n: int, field: str) -> str:
        """Generate ASCII diagram for octahedral geometry."""
        t2g_e = data["t2g_e"]
        eg_e = data["eg_e"]

        # Build electron symbols for t2g
        t2g_symbols = self._orbital_symbols(t2g_e, 3)
        eg_symbols = self._orbital_symbols(eg_e, 2)

        lines = []
        lines.append("    Octahedral (Oh) Crystal Field Splitting")
        lines.append(f"    d^{d_n} | {'Strong Field (Low Spin)' if field == 'strong' else 'Weak Field (High Spin)'}")
        lines.append("")
        lines.append("           Energy")
        lines.append("              │")
        lines.append(f"    ┌─────────┼──────────────┐")
        lines.append(f"    │         │  eg (dx²-y², dz²)  ← +0.6Δo  (+3/5Δo)")
        lines.append(f"    │    {' '.join(eg_symbols)}             │")
        lines.append(f"    ├─────────┼──────────────┤")
        lines.append(f"    │                  Δo     │")
        lines.append(f"    │         │                     │")
        lines.append(f"    │  t2g (dxy, dxz, dyz)  ← -0.4Δo  (-2/5Δo)")
        lines.append(f"    │   {'  '.join(t2g_symbols)}          │")
        lines.append(f"    └─────────┴──────────────┘")
        lines.append("")
        lines.append(f"    Config: {data['config']}")
        lines.append(f"    Unpaired e⁻: {data['unpaired']}  |  CFSE: {data['cfse_dq']:.2f}Δo")

        return "\n".join(lines)

    def _make_ascii_td(self, data: dict, d_n: int) -> str:
        """Generate ASCII diagram for tetrahedral geometry."""
        e_e = data["e_e"]
        t2_e = data["t2_e"]

        e_symbols = self._orbital_symbols(e_e, 2)
        t2_symbols = self._orbital_symbols(t2_e, 3)

        lines = []
        lines.append("    Tetrahedral (Td) Crystal Field Splitting")
        lines.append(f"    d^{d_n} | High Spin (always)")
        lines.append("")
        lines.append("           Energy")
        lines.append("              │")
        lines.append(f"    ┌─────────┼──────────────────┐")
        lines.append(f"    │         │  t2 (dxy,dxz,dyz)  ← +0.4Δt  (+2/5Δt)")
        lines.append(f"    │   {'  '.join(t2_symbols)}            │")
        lines.append(f"    ├─────────┼──────────────────┤")
        lines.append(f"    │                   Δt       │  (Δt ≈ 4/9 Δo)")
        lines.append(f"    │         │                      │")
        lines.append(f"    │      e (dz2,dx2-y2)  ← -0.6Δt  (-3/5Δt)")
        lines.append(f"    │       {' '.join(e_symbols)}       │")
        lines.append(f"    └─────────┴──────────────────┘")
        lines.append("")
        lines.append(f"    Config: {data['config']}")
        lines.append(f"    Unpaired e⁻: {data['unpaired']}  |  CFSE: {data['cfse_dqt']:.2f}Δt")

        return "\n".join(lines)

    def _make_ascii_sq(self, data: dict, d_n: int) -> str:
        """Generate ASCII diagram for square planar geometry."""
        el = data["electrons"]
        orb_names = ["dxz", "dyz", "dxy", "dz2", "dx²-y²"]
        energies = ["-4.28Dq", "-4.28Dq", "-2.28Dq", "+0.86Dq", "+12.28Dq"]

        lines = []
        lines.append("    Square Planar (D4h) Crystal Field Splitting")
        lines.append(f"    d^{d_n}")
        lines.append("")
        lines.append("           Energy")
        lines.append("              │")
        for i in range(4, -1, -1):
            sym = self._orbital_symbols(el.get(orb_names[i].replace("²", "2").replace("-", "_"), 0) if isinstance(el, dict) else list(el.values())[i], 1)[0]
            bar = "━" * (18 - len(energies[i]))
            lines.append(f"    ┌─────────┼{'─' * 25}┐")
            lines.append(f"    │    {sym}  │  {orb_names[i]:8s}  ← {energies[i]:>8s} │")
        lines.append(f"    └─────────┴{'─' * 25}┘")
        lines.append("")
        lines.append(f"    Config: {data['config']}")
        lines.append(f"    Unpaired e⁻: {data['unpaired']}  |  CFSE: {data['cfse_dq']:.2f}Dq")

        return "\n".join(lines)

    @staticmethod
    def _orbital_symbols(n_electrons: int, n_orbitals: int) -> List[str]:
        """Generate electron occupancy symbols for orbitals."""
        symbols = []
        remaining = n_electrons
        for i in range(n_orbitals):
            if remaining >= 2:
                symbols.append("↑↓")
                remaining -= 2
            elif remaining == 1:
                symbols.append("↑ ")
                remaining -= 1
            else:
                symbols.append("  ")
        return symbols

    def _run_base(self, geometry: str, d_count: int, field_strength: str = "weak", ligand: str = "H2O") -> dict:
        """Core logic: generate ligand field diagram."""
        geometry = geometry.lower().strip()
        if geometry not in ("octahedral", "tetrahedral", "square_planar"):
            raise ChemMCPError("Geometry must be 'octahedral', 'tetrahedral', or 'square_planar'.")

        if not isinstance(d_count, int) or d_count < 0 or d_count > 10:
            raise ChemMCPError("d_count must be an integer between 0 and 10.")

        field = field_strength.lower().strip()
        delta_factor = self._get_delta_factor(ligand)

        if geometry == "octahedral":
            data = self._fill_orbitals_oh(d_count, field)
            ascii_diag = self._make_ascii_oh(data, d_count, field)
            explanation = self._explain_oh(d_count, field, ligand, delta_factor)
            diagram_data = {
                "geometry": "octahedral",
                "delta_o_relative": delta_factor,
                "orbitals": {
                    "eg": {"energy_dq": 0.6, "orbitals": ["dx2-y2", "dz2"], "electrons": data["eg_e"]},
                    "t2g": {"energy_dq": -0.4, "orbitals": ["dxy", "dxz", "dyz"], "electrons": data["t2g_e"]},
                },
                "electron_config": data["config"],
                "spin_state": data["spin"],
                "unpaired_electrons": data["unpaired"],
                "cfse_dq": data["cfse_dq"],
                "magnetic_moment_bm": round((data["unpaired"] * (data["unpaired"] + 2)) ** 0.5, 2) if data["unpaired"] > 0 else 0.0,
            }
        elif geometry == "tetrahedral":
            data = self._fill_orbitals_td(d_count)
            ascii_diag = self._make_ascii_td(data, d_count)
            explanation = self._explain_td(d_count, ligand)
            diagram_data = {
                "geometry": "tetrahedral",
                "delta_t_relative": delta_factor * 4.0 / 9.0,
                "orbitals": {
                    "e": {"energy_dqt": -0.6, "orbitals": ["dz2", "dx2-y2"], "electrons": data["e_e"]},
                    "t2": {"energy_dqt": 0.4, "orbitals": ["dxy", "dxz", "dyz"], "electrons": data["t2_e"]},
                },
                "electron_config": data["config"],
                "spin_state": data["spin"],
                "unpaired_electrons": data["unpaired"],
                "cfse_dqt": data["cfse_dqt"],
                "magnetic_moment_bm": round((data["unpaired"] * (data["unpaired"] + 2)) ** 0.5, 2) if data["unpaired"] > 0 else 0.0,
            }
        else:  # square_planar
            data = self._fill_orbitals_sq(d_count)
            ascii_diag = self._make_ascii_sq(data, d_count)
            explanation = self._explain_sq(d_count, ligand)
            diagram_data = {
                "geometry": "square_planar",
                "orbitals": data.get("electrons", {}),
                "electron_config": data["config"],
                "spin_state": data["spin"],
                "unpaired_electrons": data["unpaired"],
                "cfse_dq": data["cfse_dq"],
                "magnetic_moment_bm": round((data["unpaired"] * (data["unpaired"] + 2)) ** 0.5, 2) if data["unpaired"] > 0 else 0.0,
            }

        return {
            "diagram_data": diagram_data,
            "ascii_diagram": ascii_diag,
            "explanation": explanation,
        }

    def _explain_oh(self, d_n: int, field: str, ligand: str, delta_f: float) -> str:
        """Generate explanation for octahedral case."""
        spin_txt = "low-spin" if field == "strong" else "high-spin"
        examples_map = {
            (6, "strong"): "[Co(NH3)6]3+, [Fe(CN)6]4-, [Co(en)3]3+",
            (6, "weak"): "[Fe(H2O)6]2+, [CoF6]3-",
            (3, "any"): "[Cr(NH3)6]3+, [Cr(H2O)6]3+",
            (9, "any"): "[Cu(H2O)6]2+, [Cu(NH3)6]2+",
            (5, "weak"): "[Mn(H2O)6]2+, [Fe(H2O)6]3+",
            (5, "strong"): "[Fe(CN)6]3-, [Mn(CN)6]3-",
        }
        ex = examples_map.get((d_n, field)) or examples_map.get((d_n, "any")) or "various complexes"
        return (
            f"In octahedral crystal field, the five degenerate d-orbitals split into lower-energy t2g "
            f"(dxy, dxz, dyz, -0.4Δo) and higher-energy eg (dx²-y², dz², +0.6Δo). "
            f"For d^{d_n} with {ligand} (Δ factor ≈ {delta_f:.2f}), the {spin_txt} configuration results. "
            f"The spectrochemical series ranks {ligand} as a {'strong' if delta_f > 1.2 else 'intermediate' if delta_f > 1.0 else 'weak'}-field ligand. "
            f"Examples: {ex}."
        )

    def _explain_td(self, d_n: int, ligand: str) -> str:
        return (
            f"In tetrahedral crystal field, d-orbitals split into lower e (dz², dx²-y², -0.6Δt) and higher t2 "
            f"(dxy, dxz, dyz, +0.4Δt). Note that Δt ≈ 4/9 Δo, so tetrahedral complexes are almost always high-spin. "
            f"For d^{d_n} with {ligand}, all electrons occupy orbitals following Hund's rule before pairing."
        )

    def _explain_sq(self, d_n: int, ligand: str) -> str:
        return (
            f"Square planar geometry can be derived from extreme octahedral elongation (Jahn-Teller) or from "
            f"strong-field d^8 configurations. The d-orbital energies (in Dq): dxz/dyz (-4.28) < dxy (-2.28) < "
            f"dz² (+0.86) < dx²-y² (+12.28). Most common for d^8 ions like Ni(II), Pd(II), Pt(II), Au(III) "
            f"with strong-field ligands like CN-, PR3, or NH3."
        )

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need at least geometry and d_count.")
            geo = parts[0]
            d_n = int(parts[1])
            field = parts[2] if len(parts) > 2 else "weak"
            ligand = parts[3] if len(parts) > 3 else "H2O"
            return self._run_base(geo, d_n, field, ligand)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'geometry d_count [field_strength] [ligand]'")
