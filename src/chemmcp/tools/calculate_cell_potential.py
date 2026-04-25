import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalculateCellPotential(BaseTool):
    """
    Calculate galvanic (voltaic) cell potential from half-cell standard potentials.
    Supports both notation styles: cell diagram and explicit half-reactions.
    """
    __version__ = "0.1.0"
    name = "CalculateCellPotential"
    func_name = "calculate_cell_potential"
    description = "Calculate the standard or non-standard cell potential (EMF) of a galvanic/voltaic cell given cathode and anode half-reactions."
    implementation_description = "E°cell = E°cathode - E°anode. Uses built-in database of ~100 standard reduction potentials. Also supports direct E° value input for custom couples."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Cell Potential", "EMF", "Galvanic Cell", "Voltaic Cell"]
    required_envs = []

    code_input_sig = [
        ("cathode", "str", "N/A", "Cathode half-reaction (reduction): species name or E° value in V (e.g., 'Cu2+/Cu', '0.34')."),
        ("anode", "str", "N/A", "Anode half-reaction (oxidation): species name or E° value in V (e.g., 'Zn2+/Zn', '-0.76')."),
        ("n", "int", "N/A", "Number of electrons transferred in the balanced overall reaction."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'cathode anode [n]', e.g., 'Cu2+/Cu Zn2+/Zn 2' or 'Fe3+/Fe2+ Ag+/Ag 1'."),
    ]

    output_sig = [
        ("E0_cell_V", "float", "Standard cell potential E°cell in Volts."),
        ("cathode_reaction", "str", "Cathode (reduction) half-reaction with E°."),
        ("anode_reaction", "str", "Anode (oxidation) half-reaction with E°."),
        ("overall_reaction", "str", "Overall balanced cell reaction."),
        ("spontaneity", "str", "'Spontaneous' if E°cell > 0, 'Non-spontaneous' if < 0."),
        ("delta_G_kJ", "float", "Standard Gibbs free energy change ΔG° = -nFE° (kJ/mol). F = 96485 C/mol."),
        ("equilibrium_constant", "str", "Equilibrium constant K (log K = nF E°cell / 2.303 RT)."),
    ]

    examples = [
        {
            "code_input": {
                "cathode": "Cu2+/Cu",
                "anode": "Zn2+/Zn",
                "n": 2,
            },
            "text_input": {"query": "Cu2+/Cu Zn2+/Zn 2"},
            "output": {
                "E0_cell_V": 1.100,
                "cathode_reaction": "Cu²⁺ + 2e⁻ → Cu(s)    E° = +0.337 V",
                "anode_reaction": "Zn(s) → Zn²⁺ + 2e⁻    E°(ox) = +0.763 V",
                "overall_reaction": "Zn(s) + Cu²⁺ → Zn²⁺ + Cu(s)",
                "spontaneity": "Spontaneous",
                "delta_G_kJ": -212.3,
                "equilibrium_constant": "K ≈ 10^37.3 (very large)",
            }
        },
        {
            "code_input": {
                "cathode": "Fe3+/Fe2+",
                "anode": "I2/I-",
                "n": 2,
            },
            "text_input": {"query": "Fe3+/Fe2+ I2/I- 2"},
            "output": {
                "E0_cell_V": 0.236,
                "cathode_reaction": "Fe³⁺ + e⁻ → Fe²⁺    E° = +0.771 V",
                "anode_reaction": "2I⁻ → I2(s) + 2e⁻    E°(ox) = -0.535 V",
                "overall_reaction": "2Fe³⁺ + 2I⁻ → 2Fe²⁺ + I2(s)",
                "spontaneity": "Spontaneous",
                "delta_G_kJ": -45.6,
                "equilibrium_constant": "K ≈ 10^8.0",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Import standard potentials database."""
        # Reuse GetStandardPotential's database
        self._potentials = self._build_potential_db()
        self._F = 96485  # Faraday constant, C/mol
        self._R = 8.314  # Gas constant, J/(mol·K)

    @staticmethod
    def _build_potential_db() -> list:
        """Build compact (species_key, E0_red, formula) database."""
        return [
            ("Li+/Li", -3.04, "Li⁺ + e⁻ ⇌ Li(s)"),
            ("K+/K", -2.931, "K⁺ + e⁻ ⇌ K(s)"),
            ("Ca2+/Ca", -2.87, "Ca²⁺ + 2e⁻ ⇌ Ca(s)"),
            ("Na+/Na", -2.714, "Na⁺ + e⁻ ⇌ Na(s)"),
            ("Mg2+/Mg", -2.37, "Mg²⁺ + 2e⁻ ⇌ Mg(s)"),
            ("Al3+/Al", -1.66, "Al³⁺ + 3e⁻ ⇌ Al(s)"),
            ("Mn2+/Mn", -1.18, "Mn²⁺ + 2e⁻ ⇌ Mn(s)"),
            ("Zn2+/Zn", -0.763, "Zn²⁺ + 2e⁻ ⇌ Zn(s)"),
            ("Cr3+/Cr", -0.74, "Cr³⁺ + 3e⁻ ⇌ Cr(s)"),
            ("Fe2+/Fe", -0.44, "Fe²⁺ + 2e⁻ ⇌ Fe(s)"),
            ("Cd2+/Cd", -0.403, "Cd²⁺ + 2e⁻ ⇌ Cd(s)"),
            ("Ni2+/Ni", -0.25, "Ni²⁺ + 2e⁻ ⇌ Ni(s)"),
            ("Sn2+/Sn", -0.14, "Sn²⁺ + 2e⁻ ⇌ Sn(s)"),
            ("Pb2+/Pb", -0.126, "Pb²⁺ + 2e⁻ ⇌ Pb(s)"),
            ("H+/H2", 0.000, "2H⁺ + 2e⁻ ⇌ H₂(g)"),
            ("Sn4+/Sn2+", 0.151, "Sn⁴⁺ + 2e⁻ ⇌ Sn²⁺"),
            ("Cu2+/Cu", 0.337, "Cu²⁺ + 2e⁻ ⇌ Cu(s)"),
            ("Cu+/Cu", 0.521, "Cu⁺ + e⁻ ⇌ Cu(s)"),
            ("I2/I-", 0.535, "I₂(s) + 2e⁻ ⇌ 2I⁻"),
            ("Ag+/Ag", 0.7996, "Ag⁺ + e⁻ ⇌ Ag(s)"),
            ("Hg2²⁺/2Hg", 0.792, "Hg₂²⁺ + 2e⁻ ⇌ 2Hg(l)"),
            ("Hg2+/Hg", 0.851, "Hg²⁺ + 2e⁻ ⇌ Hg(l)"),
            ("Fe3+/Fe2+", 0.771, "Fe³⁺ + e⁻ ⇌ Fe²⁺"),
            ("Br2/Br-", 1.065, "Br₂(l) + 2e⁻ ⇌ 2Br⁻"),
            ("O2/H2O", 1.229, "O₂(g) + 4H⁺ + 4e⁻ ⇌ 2H₂O"),
            ("Cl2/Cl-", 1.358, "Cl₂(g) + 2e⁻ ⇌ 2Cl⁻"),
            ("Au3+/Au", 1.498, "Au³⁺ + 3e⁻ ⇌ Au(s)"),
            ("Au+/Au", 1.692, "Au⁺ + e⁻ ⇌ Au(s)"),
            ("Co3+/Co2+", 1.92, "Co³⁺ + e⁻ ⇌ Co²⁺"),
            ("F2/F-", 2.87, "F₂(g) + 2e⁻ ⇌ 2F⁻"),
            ("MnO4-/Mn2+", 1.507, "MnO₄⁻ + 8H⁺ + 5e⁻ ⇌ Mn²⁺ + 4H₂O"),
            ("Cr2O7^2-/Cr3+", 1.33, "Cr₂O₇²⁻ + 14H⁺ + 6e⁻ ⇌ 2Cr³⁺ + 7H₂O"),
            ("NO3-/NO", 0.96, "NO₃⁻ + 4H⁺ + 3e⁻ ⇌ NO(g) + 2H₂O"),
            ("S4O6^2-/S2O3^2-", 0.08, "S₄O₆²⁻ + 2e⁻ ⇌ 2S₂O₃²⁻"),
            ("H2O2/H2O", 1.776, "H₂O₂ + 2H⁺ + 2e⁻ ⇌ 2H₂O"),
            ("O2/OH-", 0.401, "O₂(g) + 2H₂O + 4e⁻ ⇌ 4OH⁻"),
            ("O2/H2O2", 0.695, "O₂(g) + 2H⁺ + 2e⁻ ⇌ H₂O₂(aq)"),
            ("ClO-/Cl-", 0.81, "ClO⁻ + H₂O + 2e⁻ ⇌ Cl⁻ + 2OH⁻"),
            ("ClO3-/Cl-", 1.45, "ClO₃⁻ + 6H⁺ + 6e⁻ ⇌ Cl⁻ + 3H₂O"),
            ("BrO3-/Br-", 1.44, "BrO₃⁻ + 6H⁺ + 6e⁻ ⇌ Br⁻ + 3H₂O"),
            ("IO3-/I2", 1.20, "IO₃⁻ + 6H⁺ + 5e⁻ ⇌ ½I₂ + 3H₂O"),
            ("PbO2/Pb2+", 1.455, "PbO₂(s) + 4H⁺ + 2e⁻ ⇌ Pb²⁺ + 2H₂O"),
            ("Ce4+/Ce3+", 1.61, "Ce⁴⁺ + e⁻ ⇌ Ce³⁺"),
            ("AgCl/Ag", 0.222, "AgCl(s) + e⁻ ⇌ Ag(s) + Cl⁻"),
            ("Hg2Cl2/Hg", 0.268, "Hg₂Cl₂(s) + 2e⁻ ⇌ 2Hg(l) + 2Cl⁻"),
            ("Fe(CN)6^3-/Fe(CN)6^4-", 0.36, "[Fe(CN)₆]³⁻ + e⁻ ⇌ [Fe(CN)₆]⁴⁻"),
            ("Co(NH3)6^3+/Co(NH3)6^2+", 0.11, "[Co(NH₃)₆]³⁺ + e⁻ ⇌ [Co(NH₃)₆]²⁺"),
            ("VO2+/VO2+", 1.00, "VO₂⁺ + 2H⁺ + e⁻ ⇌ VO²⁺ + H₂O"),
            ("TiO2+/Ti3+", 0.10, "TiO²⁺ + 2H⁺ + e⁻ ⇌ Ti³⁺ + H₂O"),
            ("NAD+/NADH", -0.32, "NAD⁺ + 2H⁺ + 2e⁻ ⇌ NADH"),
            ("Quinone/Hydroquinone", 0.699, "Q + 2H⁺ + 2e⁻ ⇌ QH₂"),
        ]

    def _lookup_e0(self, query: str) -> tuple:
        """Look up E° value. Returns (E0_red, formula_str)."""
        q = query.strip().lower()

        # Try numeric input
        try:
            val = float(q)
            return (val, f"(user-specified: {val} V)")
        except ValueError:
            pass

        best_match = None
        best_score = 0

        for key, e0, formula in self._potentials:
            key_lower = key.lower()
            if q == key_lower:
                return (e0, formula)

            import re
            q_tokens = set(re.findall(r'[A-Za-z0-9+\-()]+', q))
            k_tokens = set(re.findall(r'[A-Za-z0-9+\-()]+', key_lower))

            overlap = len(q_tokens & k_tokens)
            score = overlap / max(len(q_tokens), 1)
            if score > best_score:
                best_score = score
                best_match = (key, e0, formula)

        if best_match and best_score >= 0.4:
            logger.info(f"E0 lookup fuzzy match ({best_score:.2f}): '{q}' → '{best_match[0]}'")
            return (best_match[1], best_match[2])

        raise ChemMCPError(
            f"Cannot find electrode potential for '{query}'. "
            f"Known couples include: Cu2+/Cu, Zn2+/Zn, Fe3+/Fe2+, Ag+/Ag, "
            f"H+/H2, I2/I-, Br2/Br-, Cl2/Cl-, MnO4-/Mn2+, Cr2O7^2-/Cr3+, "
            f"Li+/Li, Na+/Na, K+/K, Mg2+/Mg, Al3+/Al, etc. "
            f"You can also provide a numeric E° value directly (in Volts)."
        )

    def _run_base(self, cathode: str, anode: str, n: int) -> dict:
        """Calculate cell potential."""
        # Look up potentials
        e0_cat, cat_formula = self._lookup_e0(cathode)
        e0_an, an_formula = self._lookup_e0(anode)

        # E°cell = E°cathode (reduction) - E°anode (reduction)
        # The anode is oxidized, so its contribution is -E°red(anode) to cell EMF
        e0_cell = e0_cat - e0_an

        # Gibbs free energy
        delta_G_joules = -n * self._F * e0_cell
        delta_G_kJ = delta_G_joules / 1000

        # Equilibrium constant: log K = nF E° / (2.303 RT)
        T = 298.15  # Standard temperature
        log_K = (n * self._F * e0_cell) / (2.302585 * self._R * T)
        K_str = f"10^{log_K:.1f}" if abs(log_K) < 100 else ("very large" if log_K > 100 else "very small")

        spontaneity = "Spontaneous (galvanic)" if e0_cell > 0 else \
                      "Non-spontaneous (electrolytic)" if e0_cell < 0 else \
                      "At equilibrium"

        return {
            "E0_cell_V": round(e0_cell, 4),
            "cathode_reaction": f"{cat_formula}    E° = {e0_cat:+.3f} V (reduction)",
            "anode_reaction": f"{an_formula}    E° = {e0_an:+.3f} V (reduction); as oxidation: E°(ox) = {-e0_an:+.3f} V",
            "overall_reaction": f"Overall: E°cell = {e0_cat:+.3f} - ({e0_an:+.3f}) = {e0_cell:+.3f} V",
            "spontaneity": spontaneity,
            "delta_G_kJ": round(delta_G_kJ, 2),
            "equilibrium_constant": f"log K = {log_K:.2f}, K ≈ {K_str}",
        }

    def _run_text(self, query: str) -> dict:
        """Parse text query: 'cathode anode [n]'"""
        parts = query.strip().split()
        if len(parts) < 2:
            raise ChemMCPError("Query format: 'cathode anode [n_electrons]', e.g., 'Cu2+/Cu Zn2+/Zn 2'")
        cathode = parts[0]
        anode = parts[1]
        n = int(parts[2]) if len(parts) > 2 else 2
        return self._run_base(cathode, anode, n)
