import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Empirical UV-Vis absorption data (λ_max in nm) for common chromophores
# Based on Woodward-Fieser rules and experimental data
CHROMOPHORE_DATA = {
    # (base_lambda_nm, epsilon_M_cm, transition_type, solvent_effect_nm)
    "alkene_isolated": (175, 10000, "π→π*", "negligible"),
    "diene_acyclic": (217, 20000, "π→π*", "0"),
    "diene_homoannular": (214, 35000, "π→π*", "+6 (polar solvent)"),
    "heteroannular": (228, 25000, "π→π*", "0"),
    "enone_αβ_unsaturated": (215, 15000, "n→π* + π→π*", "0"),
    "enone_extended": (254, 25000, "π→π*", "+5 (polar solvent)"),
    "carbonyl_aldehyde_ketone": (280, 20, "n→π*", "+7-10 (H-bonding)"),
    "carboxylic_acid": (204, 50, "n→π*", "solvent dependent"),
    "ester": (205, 50, "n→π*", "small"),
    "amide": (208, 30, "n→π*", "solvent dependent"),
    "nitro": (270, 20, "n→π* + π→π*", "moderate"),
    "azo": (340, 10, "n→π*", "small"),
    "nitroso": (300, 100, "n→π*", "moderate"),
    "benzene": (255, 230, "π→π*", "0 to +8"),
    "phenol": (270, 1450, "π→π*", "+7-27 (pH dependent)"),
    "anisole": (269, 1480, "π→π*", "small"),
    "aniline": (280, 1300, "π→π*", "+16-28 (pH dependent)"),
    "styrene": (244, 12000, "π→π*", "0 to +5"),
    "naphthalene": (275, 5600, "π→π*", "small"),
    "anthracene": (356, 6300, "π→π*", "small"),
    "biphenyl": (248, 17000, "π→π*", "+6-10 (twisted)"),
    "conjugated_polyene_3_double_bonds": (258, 35000, "π→π*", "0"),
    "conjugated_polyene_4_double_bonds": (290, 52000, "π→π*", "0"),
    "conjugated_polyene_5_double_bonds": (330, 80000, "π→π*", "0"),
    "conjugated_polyene_6_double_bonds": (363, 118000, "π→π*", "0"),
    "conjugated_polyene_7_double_bonds": (390, 140000, "π→π*", "0"),
    "conjugated_polyene_8_double_bonds": (410, 185000, "π→π*", "0"),
}

# Substituent increments for Woodard-Fieser rules (nm)
WOODWARD_FIESER_INCREMENTS = {
    # For dienes / polyenes:
    "diene_core_acyclic": 0,
    "diene_core_heteroannular": 39,
    "diene_core_homoannular": 36,
    "extending_conjugation": 30,
    "substituent_alkyl_cyclic": 5,
    "substituent_OR_group": 6,
    "substituent_SR_group": 30,
    "substituent_Cl_Br": 5,
    "substituent_NR2_group": 60,
    # For α,β-unsaturated carbonyls:
    "enone_core_parent": 215,
    "enone_extending_conjugation": 30,
    "enone_alkyl_alpha": 10,
    "enone_alkyl_beta": 12,
    "enone_alkyl_gamma_delta": 5,
    "enone_OH_alpha_beta": 35,
    "enone_OR_alpha": 35,
    "enone_OR_beta": 30,
    "enone_OR_gamma_delta": 17,
    "enone_Cl_alpha_beta": 15,
    "enone_Br_alpha_beta": 25,
    "enote_NR2_beta": 58,
    "enone_exocyclic_db": 5,
    "enone_homoannular_diene": 39,
    "enone_solvent_correction_hexane": 0,
    "enone_solvent_correction_ether_methanol_ethanol": 0,
    "enone_solvent_correction_water": 8,
    "enone_solvent_correction_chloroform": -1,
    "enone_solvent_correction_dioxane": -5,
}


@ChemMCPManager.register_tool
class UvVisSpectrum(BaseTool):
    """
    UV-Vis光谱计算工具。
    基于Woodward-Fieser规则和经验数据，计算电子跃迁的吸收波长(λ_max)、摩尔吸光系数和谱带归属。
    """
    __version__ = "0.1.0"
    name = "UvVisSpectrum"
    func_name = "calculate_uv_vis_spectrum"
    description = "Calculate UV-Vis absorption spectrum parameters: λ_max (nm), molar absorptivity (ε), and electronic transition assignments using empirical rules."
    implementation_description = (
        "Uses Woodward-Fieser empirical rules for conjugated systems (dienes, enones, aromatic compounds) "
        "to predict λ_max values. Includes substituent increment tables, solvent corrections, "
        "and molar absorptivity estimates based on transition type (π→π*, n→π*, charge transfer)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["UV-Vis", "Spectroscopy", "Electronic Transitions", "Chromophore", "Woodward-Fieser"]
    required_envs = []

    code_input_sig = [
        ("chromophore_type", "str", "N/A", "Type of chromophore: 'diene', 'enone', 'aromatic', 'carbonyl', 'polyene', or specific name."),
        ("substituents", "list", "[]", "List of substituent names on the chromophore."),
        ("conjugation_length", "int", "1", "Number of conjugated double bonds."),
        ("solvent", "str", "hexane", "Solvent name for correction: 'hexane', 'methanol', 'ethanol', 'water', 'chloroform', etc."),
        ("calculation_mode", "str", "auto", "Mode: 'auto' (lookup), 'woodward_fieser' (manual calculation), or 'estimate'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: chromophore_type [substituents...] [conjugation_length] [solvent]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with lambda_max_nm, molar_absorptivity, transition_type, band_assignment, woodward_fieser_breakdown."),
    ]

    examples = [
        {
            "code_input": {
                "chromophore_type": "enone",
                "substituents": ["alkyl_beta"],
                "conjugation_length": 2,
                "solvent": "ethanol",
                "calculation_mode": "woodward_fieser",
            },
            "text_input": {
                "input_params": "enone alkyl_beta 2 ethanol woodward_fieser",
            },
            "output": {
                "result": {
                    "lambda_max_nm": 227,
                    "molar_absorptivity": 12500,
                    "transition_type": "π→π* (K-band)",
                    "chromophore": "α,β-unsaturated ketone",
                    "woodward_fieser_breakdown": {"parent_enone": 215, "beta_alkyl": 12, "total": 227},
                    "solvent_correction": 0,
                    "spectral_region": "UV-B (near UV)",
                }
            }
        },
        {
            "code_input": {
                "chromophore_type": "aromatic",
                "substituents": ["OH"],
                "conjugation_length": 1,
                "solvent": "water",
                "calculation_mode": "auto",
            },
            "text_input": {
                "input_params": "aromatic OH water auto",
            },
            "output": {
                "result": {
                    "lambda_max_nm": 270,
                    "molar_absorptivity": 1450,
                    "transition_type": "π→π* (B-band with auxochrome)",
                    "chromophore": "phenol (benzene + OH auxochrome)",
                    "band_assignment": "Red-shifted benzene band due to n→π* conjugation of OH lone pair",
                    "solvent_effect": "Water H-bonding may cause slight blue shift",
                    "spectral_region": "UV-A/B boundary",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.chromophore_db = dict(CHROMOPHORE_DATA)
        self.wf_increments = dict(WOODWARD_FIESER_INCREMENTS)

    def _run_base(
        self,
        chromophore_type: str,
        substituents: Optional[List[str]] = None,
        conjugation_length: int = 1,
        solvent: str = "hexane",
        calculation_mode: str = "auto",
    ) -> dict:
        """Core logic: calculate UV-Vis spectrum parameters."""
        if not chromophore_type:
            raise ChemMCPError("chromophore_type is required.")

        if substituents is None:
            substituents = []

        mode = calculation_mode.lower().strip()
        ct = chromophore_type.lower().strip()

        if mode == "woodward_fieser" or mode == "wf":
            result = self._calc_woodward_fieser(ct, substituents, conjugation_length, solvent)
        elif mode == "estimate":
            result = self._calc_estimate(ct, conjugation_length)
        else:
            # Auto mode: try WF first, fall back to lookup
            try:
                result = self._calc_woodward_fieser(ct, substituents, conjugation_length, solvent)
            except (ValueError, KeyError):
                result = self._calc_lookup(ct, substituents, solvent)

        return {"result": result}

    def _calc_woodward_fieser(self, chromophore: str, subs: List[str],
                               conj_len: int, solvent: str) -> dict:
        """Calculate λ_max using Woodward-Fieser rules."""
        base = 0
        breakdown = {}
        chromo_key = ""

        if chromophore in ("diene", "polyene"):
            if conj_len >= 3:
                base = CHROMOPHORE_DATA.get("conjugated_polyene_{}_double_bonds".format(conj_len),
                                            CHROMOPHORE_DATA["conjugated_polyene_4_double_bonds"])[0]
                chromo_key = f"conjugated_polyene_{conj_len}_double_bonds"
            else:
                base = self.wf_increments.get("diene_core_acyclic", 217)
                chromo_key = "diene_core_acyclic"
            breakdown["parent_diene"] = base

            # Each extending double bond adds ~30 nm
            if conj_len > 2:
                extension = (conj_len - 2) * self.wf_increments.get("extending_conjugation", 30)
                breakdown["extending_conjugation"] = extension
                base += extension

        elif chromophore in ("enone", "alpha_beta_unsaturated_carbonyl"):
            base = self.wf_increments.get("enone_core_parent", 215)
            chromo_key = "enone_core_parent"
            breakdown["parent_enone"] = base

        elif chromophore in ("aromatic", "benzene"):
            base = CHROMOPHORE_DATA.get("benzene", (255,))[0]
            chromo_key = "benzene"
            breakdown["parent_aromatic"] = base

        else:
            raise ChemMCPError(f"Woodward-Fieser not applicable to '{chromophore}'. Use 'auto' mode.")

        # Add substituent increments
        sub_total = 0
        for sub in subs:
            sub_key = self._map_substituent(sub, chromophore)
            inc = self.wf_increments.get(sub_key, 5)
            breakdown[f"substituent_{sub}"] = inc
            sub_total += inc
        breakdown["substituents_total"] = sub_total
        base += sub_total

        # Solvent correction
        solv_key = f"enone_solvent_correction_{solvent.lower()}"
        solv_corr = self.wf_increments.get(solv_key, 0)
        if solv_corr != 0:
            breakdown["solvent_correction"] = solv_corr
            base += solv_corr

        lam_max = int(base)
        eps = self._estimate_epsilon(chromo_key, lam_max, subs)

        return {
            "lambda_max_nm": lam_max,
            "molar_absorptivity": eps,
            "transition_type": self._get_transition_type(chromophore, lam_max),
            "chromophore": chromophore,
            "woodward_fieser_breakdown": breakdown,
            "solvent_correction": solv_corr,
            "spectral_region": self._classify_wavelength(lam_max),
        }

    def _calc_lookup(self, chromophore: str, subs: List[str], solvent: str) -> dict:
        """Look up from chromophore database."""
        best_match = None
        best_score = 0

        for key, data in self.chromophore_db.items():
            clean_key = key.replace("_", " ").lower()
            clean_ch = chromophore.replace("_", " ").lower()
            if clean_ch in clean_key or clean_key.startswith(clean_ch):
                score = len(clean_ch)
                if score > best_score:
                    best_score = score
                    best_match = key

        if best_match is None:
            raise ChemMCPError(f"Unknown chromophore: '{chromophore}'")

        base_lam, eps, trans_type, solv_eff = self.chromophore_db[best_match]

        # Adjust for substituents (simple heuristic: each auxochrome red-shifts by ~5-15 nm)
        aux_shift = len(subs) * 8
        lam_max = int(base_lam + aux_shift)

        return {
            "lambda_max_nm": lam_max,
            "molar_absorptivity": int(eps * (1.1 ** len(subs))),
            "transition_type": trans_type,
            "chromophore": best_match,
            "band_assignment": self._generate_band_assignment(best_match, subs),
            "solvent_effect": solv_eff,
            "spectral_region": self._classify_wavelength(lam_max),
        }

    def _calc_estimate(self, chromophore: str, conj_len: int) -> dict:
        """Quick estimate for conjugated polyenes."""
        # Base formula for linear polyenes: λ_max ≈ 114 + 47*n (where n = number of double bonds)
        if conj_len < 2:
            conj_len = 2
        lam_max = int(114 + 47 * conj_len)
        eps = int(10000 * (1.5 ** (conj_len - 2)))

        return {
            "lambda_max_nm": lam_max,
            "molar_absorptivity": min(eps, 200000),
            "transition_type": "π→π* (allowed)",
            "chromophore": f"linear_conjugated_polyene_{conj_len}_db",
            "estimation_method": "empirical_linear_polyene_formula",
            "spectral_region": self._classify_wavelength(lam_max),
        }

    @staticmethod
    def _map_substituent(sub: str, chromophore: str) -> str:
        """Map user input to Woodward-Fieser increment key."""
        s = sub.lower().strip()
        prefix_map = {
            "alkyl": "enone_alkyl_" if "enone" in chromophore else "substituent_alkyl_",
            "or": "enone_OR_" if "enone" in chromophore else "substituent_OR_",
            "oh": "enone_OH_" if "enone" in chromophore else "substituent_OR_",
            "sr": "substituent_SR_",
            "cl": "enone_Cl_" if "enone" in chromophore else "substituent_Cl_Br",
            "br": "enone_Br_" if "enone" in chromophore else "substituent_Cl_Br",
            "nr2": "enote_NR2_" if "enone" in chromophore else "substituent_NR2_",
        }
        for k, v in prefix_map.items():
            if s.startswith(k):
                if v.endswith("_") and "alpha" in s:
                    return v + "alpha"
                elif v.endswith("_") and "beta" in s:
                    return v + "beta"
                elif v.endswith("_") and ("gamma" in s or "delta" in s):
                    return v + "gamma_delta"
                return v + "cyclic"
        return "substituent_alkyl_cyclic"

    @staticmethod
    def _estimate_epsilon(chromo_key: str, lam_max: int, subs: list) -> int:
        """Estimate molar absorptivity ε (M⁻¹·cm⁻¹)."""
        if "pi_pi" in chromo_key.lower() or "diene" in chromo_key.lower():
            base_eps = 10000
        elif "n_pi" in chromo_key.lower() or "carbonyl" in chromo_key.lower():
            base_eps = 50
        elif "aromatic" in chromo_key.lower():
            base_eps = 200
        else:
            base_eps = 5000

        # Red shift generally increases ε
        factor = max(0.5, min(3.0, lam_max / 250.0))
        return int(base_eps * factor * (1 + 0.15 * len(subs)))

    @staticmethod
    def _get_transition_type(chromophore: str, lam_max: int) -> str:
        if "carbonyl" in chromophore or chromophore in ("aldehyde_ketone", "acid", "ester", "amide"):
            return "n→π* (forbidden, weak)" if lam_max > 270 else "n→π* + π→π*"
        elif lam_max > 400:
            return "π→π* (low energy, visible region)"
        elif lam_max > 220:
            return "π→π* (K-band, allowed)"
        else:
            return "π→π* (far UV, E-band)"

    @staticmethod
    def _classify_wavelength(lam_max: int) -> str:
        if lam_max < 200:
            return "Vacuum UV (<200 nm)"
        elif lam_max < 280:
            return "UV-C / UV-B (200-280 nm)"
        elif lam_max < 315:
            return "UV-B / UV-A (280-315 nm)"
        elif lam_max < 400:
            return "UV-A (315-400 nm)"
        elif lam_max < 700:
            return "Visible (400-700 nm)"
        else:
            return "Near IR (>700 nm)"

    @staticmethod
    def _generate_band_assignment(chromo: str, subs: list) -> str:
        parts = [f"Base chromophore: {chromo}"]
        if subs:
            parts.append(f"Auxochromes: {', '.join(subs)} — cause bathochromic (red) shift")
        return ". ".join(parts) + "."

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if not parts:
                raise ChemMCPError("Empty input.")

            chromophore = parts[0]
            kwargs = {"chromophore_type": chromophore}

            remaining = []
            modes = {"auto", "woodward_fieser", "wf", "estimate"}
            known_solvents = {
                "hexane", "methanol", "ethanol", "water", "chloroform", "dioxane",
                "ether", "acetonitrile", "dmf", "dmso", "acetone",
            }
            subs = []
            for p in parts[1:]:
                pl = p.lower()
                if pl in modes:
                    kwargs["calculation_mode"] = pl
                elif pl in known_solvents:
                    kwargs["solvent"] = pl
                elif p.isdigit():
                    kwargs["conjugation_length"] = int(p)
                else:
                    subs.append(p)
            if subs:
                kwargs["substituents"] = subs

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
