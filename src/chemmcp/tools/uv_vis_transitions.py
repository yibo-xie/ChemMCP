import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
H = 6.62607015e-34        # J·s, Planck constant
C = 2.99792458e8           # m/s, speed of light
NA = 6.02214076e23         # mol⁻¹, Avogadro's number
EV_TO_J = 1.602176634e-19  # J per eV


# Common chromophores and their UV-Vis absorption data
# (lambda_max in nm, epsilon M⁻¹cm⁻¹, transition type, solvent)
CHROMOPHORE_DATA = {
    # π → π* transitions
    "ethylene (C=C)": {"lambda_nm": 165, "epsilon": 15000, "transition": "π→π*", "notes": "Isolated alkene"},
    "1,3-butadiene": {"lambda_nm": 217, "epsilon": 21000, "transition": "π→π*", "notes": "Conjugated diene"},
    "isoprene unit": {"lambda_nm": 222, "epsilon": 23700, "transition": "π→π*", "notes": "Each additional C=C adds ~30-40 nm"},
    "β-carotene (11 conj.)": {"lambda_nm": 452, "epsilon": 152000, "transition": "π→π*", "notes": "Highly conjugated polyene"},
    "benzene": {"lambda_nm": 254, "epsilon": 250, "transition": "π→π* (forbidden)", "notes": "Weak due to symmetry"},
    "naphthalene": {"lambda_nm": 275, "epsilon": 6500, "transition": "π→π*", "notes": ""},
    "anthracene": {"lambda_nm": 356, "epsilon": 9000, "transition": "π→π*", "notes": ""},
    "phenol": {"lambda_nm": 270, "epsilon": 1450, "transition": "π→π*", "notes": "-OH auxochrome: red shift + bathochromic"},
    "aniline": {"lambda_nm": 280, "epsilon": 1430, "transition": "π→π*", "notes": "-NH₂ auxochrome"},
    "nitrobenzene": {"lambda_nm": 320, "epsilon": 12500, "transition": "π→π*", "notes": "-NO₂ strong acceptor"},
    "styrene": {"lambda_nm": 244, "transition": "π→π*", "notes": "Ph-CH=CH₂ conjugation"},
    "enone (α,β-unsat. ketone)": {"lambda_nm": 215, "epsilon": 5000, "transition": "π→π*", "notes": "Woodward-Fieser rules apply"},
    "enone (extended)": {"lambda_nm": 254, "transition": "π→π*", "notes": "With ring extension"},
    # n → π* transitions
    "ketone (saturated)": {"lambda_nm": 280, "epsilon": 15, "transition": "n→π*", "notes": "Weak, symmetry-forbidden"},
    "aldehyde": {"lambda_nm": 290, "epsilon": 12, "transition": "n→π*", "notes": ""},
    "carboxylic acid": {"lambda_nm": 204, "epsilon": 60, "transition": "n→π*", "notes": ""},
    "ester": {"lambda_nm": 207, "epsilon": 70, "transition": "n→π*", "notes": ""},
    "amide": {"lambda_nm": 205, "epsilon": 60, "transition": "n→π*", "notes": ""},
    "nitro group (-NO₂)": {"lambda_nm": 270, "epsilon": 20, "transition": "n→π*", "notes": "Also strong π→π* ~320 nm"},
    "azobenzene": {"lambda_nm": 440, "transition": "n→π*", "notes": "N=N chromophore"},
    # Charge transfer
    "iodine (I₂) in CCl₄": {"lambda_nm": 520, "transition": "σ→σ* / CT", "notes": "Visible color (violet)"},
    "KMnO₄ (aq)": {"lambda_nm": 525, "transition": "LMCT", "notes": "Intense purple color"},
    "[Cu(NH₃)₄]²⁺": {"lambda_nm": 600, "transition": "d-d", "notes": "Deep blue complex"},
    "[Cu(H₂O)₆]²⁺": {"lambda_nm": 800, "transition": "d-d", "notes": "Light blue"},
}


@ChemMCPManager.register_tool
class UvVisTransitions(BaseTool):
    """
    紫外-可见光谱跃迁能量计算工具。
    计算电子跃迁能量、预测吸收波长、应用Woodward-Fieser规则估算共轭体系λmax。
    """
    __version__ = "0.1.0"
    name = "UvVisTransitions"
    func_name = "calculate_uv_vis_transitions"
    description = "Calculate UV-Visible electronic transition energies, predict absorption wavelengths, and estimate λmax for conjugated systems using Woodward-Fieser rules."
    implementation_description = "Uses Planck's relation E=hc/λ to convert between wavelength and energy, plus Woodward-Fieser empirical rules for dienes and enones, and a chromophore database for common organic/inorganic absorbers."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Spectroscopy", "UV-Vis", "Electronic Transitions", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Calculation mode: 'energy_from_wavelength', 'wavelength_from_energy', 'woodward_fieser_diene', 'woodward_fieser_enone', or 'chromophore_lookup'."),
        ("wavelength_nm", "float", "0", "Wavelength in nm (for energy calculation)."),
        ("energy_eV", "float", "0", "Energy in eV (for wavelength calculation)."),
        ("chromophore", "str", "None", "Chromophore name for database lookup."),
        ("conjugated_diene_type", "str", "None", "Diene type for Woodward-Fieser: 'acyclic_trans_trans', 'acyclic_cis_trans', 'homoannular', 'heteroannular'."),
        ("substituents", "list", "[]", "List of substituent descriptors for Woodward-Fieser corrections."),
        ("num_conj_double_bonds", "int", "0", "Number of conjugated double bonds (for extension correction)."),
        ("base_lambda_nm", "float", "0", "Base λmax value (optional override)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Mode-specific parameters. E.g., 'energy_from_wavelength 254' or 'chromophore_lookup benzene'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing calculated energies/wavelengths, transition assignments, and spectroscopic interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "energy_from_wavelength",
                "wavelength_nm": 254,
                "energy_eV": 0,
                "chromophore": None,
                "conjugated_diene_type": None,
                "substituents": [],
                "num_conj_double_bonds": 0,
                "base_lambda_nm": 0,
            },
            "text_input": {
                "input_params": "energy_from_wavelength 254",
            },
            "output": {
                "result": {
                    "mode": "energy_from_wavelength",
                    "wavelength_nm": 254,
                    "energy_eV": round(H * C / (254e-9) / EV_TO_J, 4),
                    "energy_J": round(H * C / (254e-9), 6),
                    "wavenumber_cm-1": round(1e7 / 254, 2),
                    "frequency_hz": round(C / (254e-9), 4),
                    "color_region": "UV-C / near UV-B",
                    "interpretation": "254 nm is in the ultraviolet region; commonly used for HPLC detection and DNA absorbance.",
                }
            }
        },
        {
            "code_input": {
                "mode": "chromophore_lookup",
                "wavelength_nm": 0,
                "energy_eV": 0,
                "chromophore": "benzene",
                "conjugated_diene_type": None,
                "substituents": [],
                "num_conj_double_bonds": 0,
                "base_lambda_nm": 0,
            },
            "text_input": {
                "input_params": "chromophore_lookup benzene",
            },
            "output": {
                "result": {
                    "chromophore": "benzene",
                    "lambda_max_nm": 254,
                    "epsilon_M_cm": 250,
                    "transition_type": "π→π* (forbidden)",
                    "energy_eV": round(H * C / (254e-9) / EV_TO_J, 4),
                    "notes": "Weak absorption due to symmetry-forbidden transition.",
                }
            }
        },
        {
            "code_input": {
                "mode": "woodward_fieser_diene",
                "wavelength_nm": 0,
                "energy_eV": 0,
                "chromophore": None,
                "conjugated_diene_type": "acyclic_trans_trans",
                "substituents": ["alkyl_substituent", "alkyl_substituent"],
                "num_conj_double_bonds": 2,
                "base_lambda_nm": 0,
            },
            "text_input": {
                "input_params": "woodward_fieser_diene acyclic_trans_trans alkyl alkyl 2",
            },
            "output": {
                "result": {
                    "method": "Woodward-Fieser Rules for Dienes",
                    "base_value_nm": 214,
                    "extensions": [
                        {"item": "Acyclic trans-trans diene base", "value_nm": 214},
                        {"item": "Alkyl substituent × 2", "value_nm": 10},
                        {"item": "Double bond extending conjugation × 1", "value_nm": 30},
                    ],
                    "predicted_lambda_max_nm": 254,
                    "estimated_epsilon": "~21000",
                    "region": "UV region",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = H
        self.c = C
        self.na = NA
        self.eV_to_j = EV_TO_J

    def _run_base(self, mode: str, wavelength_nm: float = 0.0, energy_eV: float = 0.0,
                  chromophore: str = None, conjugated_diene_type: str = None,
                  substituents: List[str] = None, num_conj_double_bonds: int = 0,
                  base_lambda_nm: float = 0.0) -> dict:
        """Core logic."""
        m = mode.lower().replace("-", "_")
        if substituents is None:
            substituents = []

        if m == "energy_from_wavelength":
            return self._calc_energy(wavelength_nm)
        elif m == "wavelength_from_energy":
            return self._calc_wavelength(energy_eV)
        elif m == "chromophore_lookup":
            return self._lookup_chromophore(chromophore)
        elif m in ("woodward_fieser_diene", "woodward_fieser_diene"):
            return self._wf_diene(conjugated_diene_type, substituents, num_conj_double_bonds)
        elif m in ("woodward_fieser_enone", "woodward_fieser_enone"):
            return self._wf_enone(conjugated_diene_type, substituents, num_conj_double_bonds, base_lambda_nm)
        else:
            raise ChemMCPError(
                f"Unknown mode '{mode}'. Use: 'energy_from_wavelength', 'wavelength_from_energy', "
                "'chromophore_lookup', 'woodward_fieser_diene', or 'woodward_fieser_enone'."
            )

    def _calc_energy(self, wl_nm: float) -> dict:
        """Convert wavelength to energy."""
        if wl_nm <= 0:
            raise ChemMCPError("Wavelength must be positive.")
        wl_m = wl_nm * 1e-9
        energy_J = self.h * self.c / wl_m
        energy_eV = energy_J / self.eV_to_j
        wavenumber = 1e7 / wl_nm
        freq_hz = self.c / wl_m

        region = self._spectral_region(wl_nm)

        return {"result": {
            "mode": "energy_from_wavelength",
            "wavelength_nm": wl_nm,
            "energy_eV": round(energy_eV, 4),
            "energy_J": round(energy_J, 6),
            "energy_kJ_per_mol": round(energy_J * self.na / 1000, 2),
            "wavenumber_cm-1": round(wavenumber, 2),
            "frequency_THz": round(freq_hz / 1e12, 2),
            "color_region": region,
            "interpretation": f"{wl_nm} nm → {round(energy_eV, 2)} eV ({round(energy_J*self.na/1000, 1)} kJ/mol). {self._region_desc(region)}.",
        }}

    def _calc_wavelength(self, energy_eV: float) -> dict:
        """Convert energy to wavelength."""
        if energy_eV <= 0:
            raise ChemMCPError("Energy must be positive.")
        energy_J = energy_eV * self.eV_to_j
        wl_m = self.h * self.c / energy_J
        wl_nm = wl_m * 1e9
        wavenumber = 1e7 / wl_nm
        region = self._spectral_region(wl_nm)

        return {"result": {
            "mode": "wavelength_from_energy",
            "energy_eV": energy_eV,
            "energy_J": round(energy_J, 6),
            "wavelength_nm": round(wl_nm, 2),
            "wavenumber_cm-1": round(wavenumber, 2),
            "color_region": region,
            "interpretation": f"{energy_eV} eV → {round(wl_nm, 1)} nm ({region}). {self._region_desc(region)}.",
        }}

    def _lookup_chromophore(self, name: str) -> dict:
        """Look up chromophore data."""
        if not name:
            raise ChemMCPError("Chromophore name must be provided.")

        # First try exact match
        name_lower = name.lower()
        if name_lower in CHROMOPHORE_DATA:
            best_match = name_lower
        else:
            # Try prefix match (exact name at start of key)
            best_match = None
            best_score = 0
            for key in CHROMOPHORE_DATA:
                kl = key.lower()
                # Prefer exact or prefix match
                if kl == name_lower:
                    best_match = key
                    break
                elif kl.startswith(name_lower + " ") or name_lower in kl:
                    if len(key) > best_score:
                        best_match = key
                        best_score = len(key)

        if best_match is None:
            raise ChemMCPError(f"Chromophore '{name}' not found in database.")

        data = CHROMOPHORE_DATA[best_match]
        lam = data["lambda_nm"]
        energy_eV = self.h * self.c / (lam * 1e-9) / self.eV_to_j

        return {"result": {
            "chromophore": best_match,
            "lambda_max_nm": lam,
            "epsilon_M_cm": data.get("epsilon", "unknown"),
            "transition_type": data.get("transition", "unknown"),
            "energy_eV": round(energy_eV, 4),
            "energy_kJ_per_mol": round(energy_eV * 96.485, 2),
            "wavenumber_cm-1": round(1e7 / lam, 1),
            "spectral_region": self._spectral_region(lam),
            "notes": data.get("notes", ""),
        }}

    def _wf_diene(self, diene_type: str, substituents: List[str], n_ext: int) -> dict:
        """Woodward-Fieser rules for conjugated dienes."""
        # Base values (nm)
        bases = {
            "acyclic_trans_trans": 214,
            "acyclic_cis_trans": 253,
            "homoannular": 253,
            "heteroannular": 214,
        }

        if not diene_type or diene_type.lower() not in bases:
            raise ChemMCPError(
                f"Invalid diene type '{diene_type}'. Choose from: {list(bases.keys())}"
            )

        dt = diene_type.lower()
        total = bases[dt]
        extensions = [{"item": f"Base value: {diene_type}", "value_nm": bases[dt]}]

        # Substituent increments
        inc = {
            "alkyl_substituent": 5,
            "ring_residue": 5,
            "exocyclic_double_bond": 5,
            "halogen": 5,
            "or_substituent": 6,
            "oac_substituent": 6,
            "sr_substituent": 30,
            "nr2_substituent": 60,
        }

        for sub in substituents:
            sl = sub.lower().strip()
            if sl in inc:
                total += inc[sl]
                extensions.append({"item": f"{sub} (+{inc[sl]} nm)", "value_nm": inc[sl]})
            else:
                extensions.append({"item": f"{sub} (unrecognized, +0)", "value_nm": 0})

        # Double bond extension
        ext_val = 30 * max(0, n_ext - 2)  # Each additional double bond beyond diene
        if ext_val > 0:
            total += ext_val
            extensions.append({"item": f"Extending conjugation ×{max(0, n_ext - 2)} (+{ext_val} nm)", "value_nm": ext_val})

        return {"result": {
            "method": "Woodward-Fieser Rules for Conjugated Dienes",
            "diene_type": diene_type,
            "base_value_nm": bases[dt],
            "extensions": extensions,
            "predicted_lambda_max_nm": total,
            "estimated_epsilon": "~" + str(10000 + len(substituents) * 3000),
            "region": self._spectral_region(total),
            "disclaimer": "Woodward-Fieser rules are empirical; actual values may vary ±5-10 nm depending on solvent.",
        }}

    def _wf_enone(self, enone_type: str, substituents: List[str], n_ext: int, base: float) -> dict:
        """Woodward-Fieser rules for α,β-unsaturated ketones (enones)."""
        base_val = base if base > 0 else 215  # Acyclic enone default
        total = base_val
        extensions = [{"item": f"Base value (enone)", "value_nm": base_val}]

        # Enone substituent increments (nm)
        inc = {
            "alpha_alkyl": 10,
            "beta_alkyl": 12,
            "gamma_alkyl_delta_or_higher": 5,
            "beta_exocyclic_double_bond": 5,
            "alpha_beta_gamma_delta_alkyl_in_ring_a": 7,
            "alpha_beta_alkyl_in_ring_b": 10,
            "homodienone_alpha_beta_gamma_delta_alkyl": 34,
            "hydroxyl_alpha": 35,
            "hydroxyl_beta": 30,
            "hydroxyl_gamma_delta": 50,
            "alkoxy_alpha": 35,
            "alkoxy_beta": 30,
            "acyloxy_alpha": 10,
            "halogen_alpha": 15,
            "halogen_beta": 12,
            "ring_a_six_membered": 39,
            "ring_b_five_membered": 12,
            "ring_b_six_membered": 5,
            "ring_b_seven_plus": 15,
        }

        for sub in substituents:
            sl = sub.lower().strip()
            if sl in inc:
                total += inc[sl]
                extensions.append({"item": f"{sub} (+{inc[sl]} nm)", "value_nm": inc[sl]})

        ext_val = 30 * max(0, n_ext - 2)
        if ext_val > 0:
            total += ext_val
            extensions.append({"item": f"Extended conjugation (+{ext_val} nm)", "value_nm": ext_val})

        return {"result": {
            "method": "Woodward-Fieser Rules for α,β-Unsaturated Ketones (Enones)",
            "base_value_nm": base_val,
            "extensions": extensions,
            "predicted_lambda_max_nm": total,
            "region": self._spectral_region(total),
            "disclaimer": "Empirical rules; accuracy ±5-15 nm.",
        }}

    @staticmethod
    def _spectral_region(wl_nm: float) -> str:
        if wl_nm < 200:
            return "Vacuum UV (<200 nm)"
        elif wl_nm < 280:
            return "UV-C / Far UV (200-280 nm)"
        elif wl_nm < 315:
            return "UV-B (280-315 nm)"
        elif wl_nm < 400:
            return "UV-A / Near UV (315-400 nm)"
        elif wl_nm < 420:
            return "Violet (~400-420 nm)"
        elif wl_nm < 495:
            return "Blue (420-495 nm)"
        elif wl_nm < 570:
            return "Green (495-570 nm)"
        elif wlnm := wl_nm < 590:
            return "Yellow (570-590 nm)"
        elif wl_nm < 620:
            return "Orange (590-620 nm)"
        elif wl_nm < 750:
            return "Red (620-750 nm)"
        else:
            return "Infrared (>750 nm)"

    @staticmethod
    def _region_desc(region: str) -> str:
        descs = {
            "Vacuum UV (<200 nm)": "High-energy vacuum ultraviolet; absorbed by air and quartz.",
            "UV-C / Far UV (200-280 nm)": "Germicidal UV; protein/DNA absorption region.",
            "UV-B (280-315 nm)": "Causes sunburn; partially absorbed by ozone layer.",
            "UV-A / Near UV (315-400 nm)": "Near UV; causes tanning; many organic chromophores absorb here.",
            "Violet (~400-420 nm)": "Visible violet; S→T transitions often appear here.",
            "Blue (420-495 nm)": "Visible blue; charge-transfer bands often appear here.",
            "Green (495-570 nm)": "Visible green; [Cu(NH₃)₄]²⁺ absorbs here.",
            "Yellow (570-590 nm)": "Visible yellow; complementary to blue.",
            "Orange (590-620 nm)": "Visible orange; complementary to blue-green.",
            "Red (620-750 nm)": "Visible red; d-d transitions of some metal complexes.",
            "Infrared (>750 nm)": "Infrared region; beyond electronic transitions into vibrational.",
        }
        return descs.get(region, "")

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            mode = parts[0]

            if mode == "energy_from_wavelength":
                return self._run_base(mode, wavelength_nm=float(parts[1]))
            elif mode == "wavelength_from_energy":
                return self._run_base(mode, energy_eV=float(parts[1]))
            elif mode == "chromophore_lookup":
                return self._run_base(mode, chromophore=" ".join(parts[1:]))
            elif mode in ("woodward_fieser_diene",):
                dt = parts[1] if len(parts) > 1 else "acyclic_trans_trans"
                subs = parts[2:-1] if len(parts) > 3 else []
                n_ext = int(parts[-1]) if len(parts) > 2 else 2
                return self._run_base(mode, conjugated_diene_type=dt, substituents=subs, num_conj_double_bonds=n_ext)
            else:
                raise ValueError(f"Unknown mode: {mode}")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
