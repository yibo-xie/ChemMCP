"""
螯合效应分析工具
Chelate effect analysis: compare stability of chelated vs non-chelated complexes.
"""
import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ChelateEffect(BaseTool):
    """
    分析螯合效应：比较螯合配合物与单齿配体配合物的稳定性差异。
    计算熵贡献、有效摩尔浓度、环大小影响等。
    """
    __version__ = "0.1.0"
    name = "ChelateEffect"
    func_name = "chelate_effect"
    description = "Analyze the chelate effect: compare thermodynamic stability of chelating ligands vs equivalent monodentate ligands, including entropy contribution and ring size effects."
    implementation_description = "Uses literature data for common chelate systems (en vs NH3, EDTA vs acetate, oxalate vs formate). Computes ΔΔG° from ΔΔH° and TΔΔS°, explains effective molarity concept, analyzes 5- vs 6-membered chelate rings."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Chelate Effect", "Thermodynamics", "Stability", "Entropy"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal ion, e.g., 'Cu2+', 'Ni2+', 'Co2+', 'Ca2+'."),
        ("chelating_ligand", "str", "N/A", "Chelating ligand name, e.g., 'en', 'EDTA', 'oxalate', 'acac', 'phen'."),
        ("monodentate_analog", "str", "NH3", "Equivalent monodentate ligand for comparison (default: NH3 for en)."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (default 25°C = 298.15 K)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion chelating_ligand [monodentate_analog] [temperature_K]', e.g., 'Cu2+ en' or 'Ni2+ EDTA acetate'."),
    ]

    output_sig = [
        ("metal_ion", "str", "Metal ion analyzed."),
        ("chelating_system", "dict", "Details of the chelating ligand system (denticity, ring size, log β)."),
        ("monodentate_system", "dict", "Details of the equivalent monodentate system (log β)."),
        ("stability_comparison", "dict", "log β difference, K ratio, ΔΔG° (kJ/mol)."),
        ("entropy_contribution", "dict", "TΔS° contribution to the chelate effect (kJ/mol)."),
        ("enthalpy_contribution", "dict", "ΔH° contribution if available (kJ/mol)."),
        ("effective_molarity", "str", "Effective molarity explanation."),
        ("ring_analysis", "dict", "Chelate ring size analysis (5- vs 6-membered preferred)."),
        ("explanation", "str", "Comprehensive explanation of the chelate effect for this system."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Cu2+",
                "chelating_ligand": "en",
                "monodentate_analog": "NH3",
                "temperature_k": 298.15,
            },
            "text_input": {
                "query": "Cu2+ en"
            },
            "output": {
                "metal_ion": "Cu2+",
                "chelating_system": {"ligand": "ethylenediamine (en)", "denticity": 2, "ring_size": "5-membered", "log_beta_2": 20.8},
                "monodentate_system": {"ligand": "ammonia (NH3)", "denticity": 1, "log_beta_2": 12.67, "n_ligands": 2},
                "stability_comparison": {"delta_log_beta": 8.13, "K_ratio": "1.35 × 10⁸", "delta_delta_G_kjmol": -46.4},
                "entropy_contribution": {"T_delta_S": "-20 to -30", "dominant_factor": True},
                "ring_analysis": {"ring_size": "5-membered", "stability": "optimal (minimal angle strain)"},
                "enthalpy_contribution": {"delta_H": "small (chelate effect is mostly entropic)"},
                "effective_molarity": "When one end of en binds, [NH2] local concentration at Cu2+ is very high",
                "explanation": "[Cu(en)2]2+ is ~10⁸ times more stable than [Cu(NH3)4]2+. The chelate effect is primarily entropic: one molecule of en replaces two NH3 → increase in number of free particles → +ΔS.",
            }
        },
        {
            "code_input": {
                "metal_ion": "Ca2+",
                "chelating_ligand": "EDTA",
                "monodentate_analog": "acetate",
                "temperature_k": 298.15,
            },
            "text_input": {
                "query": "Ca2+ EDTA acetate"
            },
            "output": {
                "metal_ion": "Ca2+",
                "chelating_system": {"ligand": "EDTA", "denticity": 6, "ring_size": "five 5-membered rings", "log_beta_1": 10.69},
                "monodentate_system": {"ligand": "acetate (OAc-)", "denticity": 1, "log_beta_2": "~2.0"},
                "stability_comparison": {"delta_log_beta": "~8.7", "K_ratio": "~5 × 10⁸", "delta_delta_G_kjmol": -49.6},
                "entropy_contribution": {"T_delta_S": "dominant", "explanation": "1 EDTA replaces 6 monodentates → +5 particles released"},
                "ring_analysis": {"ring_size": "five 5-membered rings", "stability_assessment": "optimal for EDTA"},
                "enthalpy_contribution": {"delta_H": "varies, often small or favorable"},
                "effective_molarity": "EDTA's 6 donor groups create huge effective molarity advantage",
                "explanation": "EDTA forms extremely stable complexes even with Ca²⁺ (a hard, class A metal). This is why EDTA is used in water softening, titrations, and medicine.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize chelate effect database with literature data."""
        # Chelate formation constants (log β) at 25°C
        # Format: (metal, chelator) → {log_beta, denticity, notes}
        self._chelate_data = {
            # Ethylenediamine (en) systems
            ("cu2+", "en"): {"log_beta": [10.71, 20.85], "denticity": 2, "ring_size": 5, "notes": "[Cu(en)]2+ blue, [Cu(en)2]2+ deep purple"},
            ("ni2+", "en"): {"log_beta": [7.96, 14.28, 18.77], "denticity": 2, "ring_size": 5, "notes": "[Ni(en)3]2+ violet"},
            ("zn2+", "en"): {"log_beta": [5.92, 11.07, 13.93], "denticity": 2, "ring_size": 5, "notes": ""},
            ("co2+", "en"): {"log_beta": [5.95, 10.80, 13.82], "denticity": 2, "ring_size": 5, "notes": ""},
            ("fe2+", "en"): {"log_beta": [4.34, 7.65, 9.70], "denticity": 2, "ring_size": 5, "notes": ""},
            ("cr3+", "en"): {"log_beta": [None, None, None], "denticity": 2, "ring_size": 5, "notes": "kinetically inert, slow substitution"},
            ("ag+", "en"): {"log_beta": [None, None], "denticity": 2, "ring_size": 6.5, "notes": "prefers linear geometry"},

            # Oxalate (ox) systems
            ("cu2+", "oxalate"): {"log_beta": [4.5, 8.3], "denticity": 2, "ring_size": 5, "notes": "[Cu(ox)2]2-"},
            ("fe3+", "oxalate"): {"log_beta": [9.4, 16.2, 20.2], "denticity": 2, "ring_size": 5, "notes": "[Fe(ox)3]3- green"},
            ("fe2+", "oxalate"): {"log_beta": [4.2, 7.3, 9.5], "denticity": 2, "ring_size": 5, "notes": ""},
            ("co3+", "oxalate"): {"log_beta": [None, None, ~20], "denticity": 2, "ring_size": 5, "notes": ""},

            # Acetylacetonate (acac) systems
            ("cu2+", "acac"): {"log_beta": [None, 16.0], "denticity": 2, "ring_size": 6, "notes": "[Cu(acac)2 volatile"},
            ("co3+", "acac"): {"log_beta": [None, None, 15], "denticity": 2, "ring_size": 6, "notes": "[Co(acac)3]"},
            ("cr3+", "acac"): {"log_beta": [None, None, ~14], "denticity": 2, "ring_size": 6, "notes": "[Cr(acac)3]"},

            # Phenanthroline (phen) systems
            ("fe2+", "phen"): {"log_beta": [5.85, 11.35, 21.35], "denticity": 2, "ring_size": 5, "notes": "ferroin indicator (red)"},
            ("cu2+", "phen"): {"log_beta": [9.1, 16.0, 21.0], "denticity": 2, "ring_size": 5, "notes": ""},
            ("ni2+", "phen"): {"log_beta": [8.8, 17.1, 24.5], "denticity": 2, "ring_size": 5, "notes": ""},

            # EDTA systems
            ("ca2+", "edta"): {"log_beta": [10.69], "denticity": 6, "ring_size": "5×5-membered", "notes": "water hardness"},
            ("mg2+", "edta"): {"log_beta": [8.64], "denticity": 6, "ring_size": "5×5-membered", "notes": ""},
            ("fe3+", "edta"): {"log_beta": [25.1], "denticity": 6, "ring_size": "5×5-membered", "notes": "very stable"},
            ("cu2+", "edta"): {"log_beta": [18.80], "denticity": 6, "ring_size": "5×5-membered", "notes": "deep blue"},
            ("zn2+", "edta"): {"log_beta": [16.50], "denticity": 6, "ring_size": "5×5-membered", "notes": ""},
            ("ni2+", "edta"): {"log_beta": [18.56], "denticity": 6, "ring_size": "5×5-membered", "notes": ""},
            ("co2+", "edta"): {"log_beta": [16.26], "denticity": 6, "ring_size": "5×5-membered", "notes": ""},
            ("pb2+", "edta"): {"log_beta": [18.3], "denticity": 6, "ring_size": "5×5-membered", "notes": ""},
            ("cd2+", "edta"): {"log_beta": [16.46], "denticity": 6, "ring_size": "5×5-membered", "notes": ""},
            ("al3+", "edta"): {"log_beta": [16.5], "denticity": 6, "ring_size": "5×5-membered", "notes": "slow formation"},
        }

        # Monodentate reference data (log β for n equivalents)
        self._mono_data = {
            "nh3": {
                "cu2+": [4.15, 7.65, 10.54, 12.67],
                "ni2+": [2.80, 5.04, 6.87, 7.91, 8.22, 8.30],
                "zn2+": [2.37, 4.81, 7.31, 9.46],
                "co2+": [2.11, 3.85, 4.94, 5.50, 5.73, 5.22],
                "ag+": [3.36, 7.23],
            },
            "acetate": {
                "cu2+": [2.22, 3.63],  # approximate
                "ca2+": [1.23, 1.82],
                "mg2+": [1.27, 2.03],
            },
            "h2o": {  # not really a comparison but for context
                "cu2+": [0],  # aqua is reference
            },
        }

        # Chelate ring size stability order
        self._ring_stability = {
            3: ("very unstable", "high angle strain, rarely observed"),
            4: ("unstable", "significant strain, rare except Pd/Pt"),
            5: ("most stable", "near-ideal bond angles, optimal for most metals"),
            6: ("stable", "slightly less than 5 for saturated rings; good for conjugated/acac"),
            7: ("less stable", "increasing strain, entropy penalty decreases"),
            8: ("unstable", "large rings often more flexible but lower effective concentration"),
        }

        # R = 8.314 J/(mol·K)
        self.R = 8.314

    def _get_chelate_info(self, metal: str, chelator: str) -> dict:
        """Look up chelate data."""
        key = (metal.lower(), chelator.lower())
        if key not in self._chelate_data:
            raise ChemMCPError(
                f"No chelate data for ({metal}, {chelator}). Available: "
                f"en, oxalate, acac, phen, EDTA with metals: Cu2+, Ni2+, Zn2+, Co2+, Fe2+/3+, "
                f"Ca2+, Mg2+, Pb2+, Cd2+, Al3+."
            )
        return self._chelate_data[key]

    def _get_mono_log_beta(self, metal: str, mono: str, n: int) -> float:
        """Get cumulative log β for n monodentate ligands."""
        key = mono.lower()
        if key not in self._mono_data:
            # Estimate: typical log K1 for weak monodentate ≈ 1-2
            return float(n * 1.5)  # rough estimate
        mdata = self._mono_data[key]
        mkey = metal.lower()
        if mkey in mdata:
            vals = mdata[mkey]
            if n <= len(vals):
                return vals[n - 1]
            return vals[-1] + (n - len(vals)) * 0.5  # extrapolate
        return float(n * 1.5)

    def _run_base(self, metal_ion: str, chelating_ligand: str,
                  monodentate_analog: str = "NH3",
                  temperature_k: float = 298.15) -> dict:
        """Analyze chelate effect for a given system."""
        T = temperature_k
        metal = metal_ion.lower().replace(" ", "")
        chel = chelating_ligand.lower().strip()
        mono = monodentate_analog.upper().strip()

        # Get chelate data
        chel_info = self._get_chelate_info(metal, chel)
        dent = chel_info["denticity"]
        beta_list = chel_info["log_beta"]
        ring_size_raw = chel_info["ring_size"]

        # Use highest coordination number available
        max_n_chel = len(beta_list)
        log_beta_chel = beta_list[-1]
        if log_beta_chel is None or (isinstance(log_beta_chel, str) and "none" in str(log_beta_chel).lower()):
            # Find first non-None value
            for v in beta_list:
                if v is not None and not (isinstance(v, str)):
                    log_beta_chel = v
                    break

        # Number of monodentate ligands for equivalent denticity
        n_mono = dent * max_n_chel  # each chelator provides 'dent' donor atoms
        # But we should compare same total donor atoms
        # If en (bidentate) × 2 = 4 N donors → compare with 4 NH3
        total_donor_atoms = dent * max_n_chel
        log_beta_mono = self._get_mono_log_beta(metal, mono, total_donor_atoms)

        # Stability comparison
        delta_log_beta = log_beta_chel - log_beta_mono
        K_ratio = 10 ** delta_log_beta
        delta_G = -delta_log_beta * 2.303 * self.R * T / 1000  # kJ/mol

        # Entropy estimation
        # When 1 chelator (dentate) replaces 'dent' monodentates:
        # Δn_particles = (dent - 1) more free particles → +ΔS
        # TΔS ≈ (dent - 1) × (10-20) kJ/mol at 298K (rough)
        delta_n = dent - 1  # net particle change per chelator
        Ts_delta_S_est = delta_n * max_n_chel * 15.0  # rough estimate kJ/mol

        # Ring analysis
        if isinstance(ring_size_raw, int):
            ring_stability, ring_note = self._ring_stability.get(ring_size_raw, ("unknown", ""))
        else:
            ring_stability = "multiple rings (see notes)"
            ring_note = ring_size_raw

        # Effective molarity explanation
        em_explanation = (
            f"Effective molarity: when one end of a multidentate ligand binds to the metal, "
            f"the local concentration of the other donor group(s) near the metal is very high "
            f"(often 10⁻²–10⁴ M), making the subsequent binding step intramolecularly favored. "
            f"For a bidentate ligand like {chelating_ligand}, this effectively gives a huge rate "
            f"and equilibrium advantage over two separate monodentate ligands."
        )

        # Build explanation
        chel_name = chelating_ligand.upper() if len(chelating_ligand) <= 4 else chelating_ligand.capitalize()
        explanation = (
            f"**Chelate Effect Analysis: [{metal_ion}({chel_name}){max_n_chel}] vs [{metal_ion}({mono}){total_donor_atoms}]**\n\n"
            f"**Key Data:**\n"
            f"• Chelating ligand ({chel_name}): log β{max_n_chel} = {log_beta_chel}\n"
            f"• Monodentate ({mono}), {total_donor_atoms} equiv: log β ≈ {log_beta_mono:.2f}\n"
            f"• Δlog β = {delta_log_beta:.2f}\n"
            f"• Stability ratio: Kchel/Kmono = 10^{delta_log_beta:.2f} = {self._format_large(K_ratio)}\n"
            f"• ΔΔG° = {delta_G:.1f} kJ/mol (more negative = more favorable)\n\n"
            f"**Why?**\n"
            f"The chelate effect is primarily **entropic**: replacing {total_donor_atoms} molecules of {mono} "
            f"with {max_n_chel} molecule(s) of {chel_name} increases the number of free particles in solution "
            f"(→ +ΔS, → more negative ΔG at constant T).\n\n"
            f"**Ring Size:** {ring_size_raw}-membered chelate ring(s) — {ring_stability}. {ring_note}\n\n"
            f"{em_explanation}"
        )

        logger.info(f"Chelate effect: {metal_ion} + {chelating_ligand} vs {mono}: Δlogβ={delta_log_beta:.2f}")

        return {
            "metal_ion": metal_ion,
            "chelating_system": {
                "ligand": f"{chelating_ligand}",
                "denticity": dent,
                "ring_size": str(ring_size_raw),
                "coordination_number": max_n_chel,
                "log_beta": log_beta_chel,
                "notes": chel_info.get("notes", ""),
            },
            "monodentate_system": {
                "ligand": mono,
                "denticity": 1,
                "n_ligands": total_donor_atoms,
                "log_beta": round(log_beta_mono, 2),
            },
            "stability_comparison": {
                "delta_log_beta": round(delta_log_beta, 2),
                "K_ratio": self._format_large(K_ratio),
                "delta_delta_G_kjmol": round(delta_G, 1),
            },
            "entropy_contribution": {
                "estimated_T_delta_S_kjmol": round(Ts_delta_S_est, 1),
                "particle_change_per_chelator": f"+{delta_n} free particles",
                "is_dominant": True,
            },
            "ring_analysis": {
                "ring_size": str(ring_size_raw),
                "stability_assessment": ring_stability,
                "note": ring_note,
                "preference": "5-membered > 6 > 7 > 4 (for most transition metals)",
            },
            "effective_molarity": em_explanation,
            "explanation": explanation,
        }

    @staticmethod
    def _format_large(x: float) -> str:
        """Format large numbers in scientific notation."""
        if abs(x) >= 1e6:
            return f"{x:.2e}"
        elif abs(x) >= 1000:
            return f"{x:.1f}"
        else:
            return f"{x:.2f}"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        if len(parts) < 2:
            raise ChemMCPError("Format: 'metal_ion chelating_ligand [monodentate_analog] [temperature_K]'. Example: 'Cu2+ en' or 'Ca2+ EDTA acetate'")
        metal = parts[0]
        chel = parts[1]
        mono = parts[2] if len(parts) > 2 else "NH3"
        temp = float(parts[3]) if len(parts) > 3 else 298.15
        return self._run_base(metal, chel, mono, temp)
