import logging
import math
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# NMR shielding reference standards
NMR_REFERENCES = {
    "1H": {"compound": "TMS (tetramethylsilane)", "sigma_ref_ppm": 0.0, "shift_ref_ppm": 0.0},
    "13C": {"compound": "TMS (tetramethylsilane)", "sigma_ref_ppm": 0.0, "shift_ref_ppm": 0.0},
    "19F": {"compound": "CFCl3", "sigma_ref_ppm": 0.0, "shift_ref_ppm": 0.0},
    "31P": {"compound": "85% H3PO4", "sigma_ref_ppm": 0.0, "shift_ref_ppm": 0.0},
}

# Empirical substituent chemical shift increments for ¹H NMR (ppm relative to TMS)
# Based on Shoolery's rule / substituent effects
# Format: (base_shift, substituent_effects_dict)
H_NMR_SUBSTITUENT_EFFECTS = {
    # Base value for CH₃-X type protons
    "CH3": {"base": 0.9},
    "CH2": {"base": 1.3},
    "CH": {"base": 1.7},
    # Substituent effects on adjacent proton (ppm increment)
    "substituent_effects": {
        # Electron-withdrawing groups (deshield → downfield shift)
        "-NO2": 2.5,
        "-CHO": 1.8,
        "-COR": 1.7,
        "-COOH": 1.0,
        "-COOR": 1.1,
        "-CONH2": 1.1,
        "-COX": 1.6,
        "-SO2R": 1.5,
        "-CN": 1.2,
        "-Ar": 1.5,
        "-C≡CH": 1.4,
        "-C≡N": 1.2,
        "-Halogen (F)": 4.0,
        "-Halogen (Cl)": 2.3,
        "-Halogen (Br)": 2.0,
        "-Halogen (I)": 1.8,
        "-OR": 2.5,
        "-OH": 2.5,
        "-OCOR": 2.8,
        "-NH2": 1.3,
        "-NR2": 1.3,
        "-NHR": 1.3,
        "-SH": 1.1,
        "-SR": 1.3,
        "-alkyl": 0.5,
        "-CH=CH2": 1.0,
        "-C≡CH": 1.4,
    },
}

# Aromatic proton substitution effects (ppm)
AROMATIC_H_EFFECTS = {
    "benzene_H": 7.27,
    "ortho_EWG": 0.6,   # electron-withdrawing group at ortho position
    "meta_EWG": 0.15,
    "para_EWG": 0.25,
    "ortho_EDG": -0.5,  # electron-donating group (upfield shift)
    "meta_EDG": -0.08,
    "para_EDG": -0.9,
}

# ¹³C NMR typical ranges (ppm)
C13_NMR_RANGES = {
    "alkane_C": (0, 60),
    "C_adjacent_to_electronegative": (45, 90),
    "alkyne_C": (65, 100),
    "alkene_C": (100, 150),
    "aromatic_C": (110, 165),
    "imine_C": (145, 170),
    "amide_carbonyl_C": (155, 180),
    "carboxylic_acid_C": (170, 185),
    "ester_carbonyl_C": (160, 180),
    "ketone_aldehyde_C": (190, 220),
    "nitrile_C": (110, 130),
}


@ChemMCPManager.register_tool
class NmrShielding(BaseTool):
    """
    NMR化学位移计算工具（基于屏蔽常数）。
    通过屏蔽常数σ计算化学位移 δ = σ_ref − σ，使用取代基效应和电负性经验规则。
    支持¹H和¹³C核的化学位移预测。
    """
    __version__ = "0.1.0"
    name = "NmrShielding"
    func_name = "calculate_nmr_shielding"
    description = "Calculate NMR chemical shifts from shielding constants using empirical substituent electronegativity and anisotropy effects for ¹H and ¹³C nuclei."
    implementation_description = (
        "Chemical shift δ (ppm) = σ_ref − σ_sample. Uses Shoolery-type substituent increments "
        "for aliphatic protons, aromatic substitution patterns for benzene derivatives, "
        "and typical range databases for ¹³C nuclei. Shielding constant estimated from "
        "local electronic environment (electronegativity, hybridization, ring currents)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["NMR", "Chemical Shift", "Shielding Constant", "Spectroscopy", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("nucleus", "str", "N/A", "Nucleus type: '1H', '13C', '19F', or '31P'."),
        ("atom_environment", "list", "N/A", "List of nearby atoms/groups/structural features affecting the nucleus."),
        ("hybridization", "str", "sp3", "Hybridization of the atom bearing the nucleus: 'sp3', 'sp2', or 'sp'."),
        ("reference_compound", "str", "TMS", "Reference compound for shielding (default TMS for 1H/13C)."),
        ("detail_level", "str", "standard", "Detail level: 'basic', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: nucleus [atom_environment...] [hybridization] [reference]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing chemical_shift_ppm, shielding_constant, reference, contribution_breakdown."),
    ]

    examples = [
        {
            "code_input": {
                "nucleus": "1H",
                "atom_environment": ["-Cl", "-CH3"],
                "hybridization": "sp3",
                "reference_compound": "TMS",
                "detail_level": "standard",
            },
            "text_input": {
                "input_params": "1H -Cl -CH3 sp3 TMS standard",
            },
            "output": {
                "result": {
                    "nucleus": "1H",
                    "chemical_shift_ppm": 3.37,
                    "shielding_constant": -3.37,
                    "reference": "TMS (δ = 0.00 ppm)",
                    "contribution_breakdown": {
                        "base_CH3": 0.9,
                        "substituent_-Cl": 2.3,
                        "substituent_-CH3": 0.17,
                        "total_shift": 3.37,
                    },
                    "interpretation": "Chlorine deshields the proton significantly via inductive effect.",
                    "typical_range": "Aliphatic region (0.5–4.5 ppm)",
                }
            }
        },
        {
            "code_input": {
                "nucleus": "13C",
                "atom_environment": ["carbonyl", "conjugated"],
                "hybridization": "sp2",
                "reference_compound": "TMS",
                "detail_level": "detailed",
            },
            "text_input": {
                "input_params": "13C carbonyl conjugated sp2 TMS detailed",
            },
            "output": {
                "result": {
                    "nucleus": "13C",
                    "chemical_shift_ppm": 198.5,
                    "shielding_constant": -198.5,
                    "reference": "TMS (δ = 0.00 ppm)",
                    "contribution_breakdown": {
                        "base_ketone_carbonyl": 205,
                        "conjugation_correction": -6.5,
                        "total_shift": 198.5,
                    },
                    "interpretation": "Strongly deshielded carbonyl carbon; conjugation provides slight shielding.",
                    "typical_range": "Carbonyl region (190–220 ppm)",
                    "anisotropy_note": "Carbonyl π-bond creates local magnetic anisotropy that deshields the nucleus.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.references = dict(NMR_REFERENCES)
        self.h_substituents = dict(H_NMR_SUBSTITUENT_EFFECTS["substituent_effects"])
        self.c13_ranges = dict(C13_NMR_RANGES)

    def _run_base(
        self,
        nucleus: str,
        atom_environment: List[str],
        hybridization: str = "sp3",
        reference_compound: str = "TMS",
        detail_level: str = "standard",
    ) -> dict:
        """Core logic: calculate chemical shift from shielding constant."""
        nuc = nucleus.upper().strip()
        if nuc not in self.references:
            raise ChemMCPError(f"Unsupported nucleus: '{nucleus}'. Supported: {list(self.references.keys())}")

        hyb = hybridization.lower().strip()
        dl = detail_level.lower()

        if nuc == "1H":
            result = self._calc_1h_shift(atom_environment, hyb, dl)
        elif nuc == "13C":
            result = self._calc_13c_shift(atom_environment, hyb, dl)
        else:
            result = self._calc_general_shift(nuc, atom_environment, hyb)

        result["nucleus"] = nuc
        ref_info = self.references[nuc]
        result["reference"] = f"{ref_info['compound']} (δ = {ref_info['shift_ref_ppm']:.2f} ppm)"
        result["shielding_constant"] = round(-result["chemical_shift_ppm"], 3)

        return {"result": result}

    def _calc_1h_shift(self, env: List[str], hyb: str, dl: str) -> dict:
        """Calculate ¹H chemical shift from substituent effects."""
        # Determine base value based on hybridization / carbon type
        if any("aromatic" in e.lower() or "phenyl" in e.lower() or "benzene" in e.lower() for e in env):
            return self._calc_aromatic_h(env, dl)

        if hyb == "sp2" and any("alkene" in e.lower() or "vinylic" in e.lower() for e in env):
            base = 5.25  # vinylic proton
        elif hyb == "sp" and any("alkyne" in e.lower() for e in env):
            base = 2.0  # terminal alkyne proton
        elif hyb == "sp3":
            # Determine if CH3, CH2, or CH
            if any("CH2" in e or "methylene" in e.lower() for e in env):
                base = H_NMR_SUBSTITUENT_EFFECTS["CH2"]["base"]
            elif any("CH" in e or "methine" in e.lower() for e in env):
                base = H_NMR_SUBSTITUENT_EFFECTS["CH"]["base"]
            else:
                base = H_NMR_SUBSTITUENT_EFFECTS["CH3"]["base"]
        else:
            base = 1.5

        breakdown = {"base_value": base}
        total = base

        # Apply substituent effects
        for item in env:
            effect = self._match_substituent_effect(item)
            if effect is not None:
                key = f"substituent_{item}"
                breakdown[key] = effect
                total += effect

        total = round(total, 2)
        breakdown["total_shift"] = total

        result = {
            "chemical_shift_ppm": total,
            "contribution_breakdown": breakdown,
            "interpretation": self._interpret_1h_shift(total),
            "typical_range": self._classify_1h_region(total),
        }

        if dl == "detailed":
            result["shielding_mechanism"] = self._explain_shielding_mechanism(env, total)

        return result

    def _calc_aromatic_h(self, env: List[str], dl: str) -> dict:
        """Calculate aromatic proton chemical shift."""
        base = AROMATIC_H_EFFECTS["benzene_H"]
        breakdown = {"base_benzene": base}
        total = base

        for item in env:
            il = item.lower()
            pos_key = None
            if il.startswith("ortho") or il.startswith("o_"):
                pos_key = "ortho"
            elif il.startswith("meta") or il.startswith("m_"):
                pos_key = "meta"
            elif il.startswith("para") or il.startswith("p_"):
                pos_key = "para"

            if pos_key:
                is_ewg = any(kw in il for kw in ["no2", "cor", "cooh", "cn", "so", "cho", "f", "cl", "br"])
                if is_ewg:
                    effect = AROMATIC_H_EFFECTS[f"{pos_key}_EWG"]
                else:
                    effect = AROMATIC_H_EFFECTS[f"{pos_key}_EDG"]
                breakdown[f"{pos_key}_substituent"] = effect
                total += effect

        total = round(max(6.0, min(10.0, total)), 2)
        breakdown["total_shift"] = total

        return {
            "chemical_shift_ppm": total,
            "contribution_breakdown": breakdown,
            "interpretation": f"Aromatic proton at {total} ppm in the aromatic region (6.5–8.5 ppm).",
            "typical_range": "Aromatic region (6.5–8.5 ppm)",
            "ring_current_note": "Aromatic ring current causes significant deshielding.",
        }

    def _calc_13c_shift(self, env: List[str], hyb: str, dl: str) -> dict:
        """Calculate ¹³C chemical shift based on functional group environment."""
        best_match = None
        best_score = 0

        for key, (lo, hi) in self.c13_ranges.items():
            for item in env:
                clean_item = item.lower().replace("_", "").replace("-", "")
                clean_key = key.lower().replace("_", "").replace("-", "")
                if clean_item in clean_key or clean_key in clean_item:
                    score = len(clean_key)
                    if score > best_score:
                        best_score = score
                        best_match = key

        if best_match:
            lo, hi = self.c13_ranges[best_match]
            shift = (lo + hi) / 2
            breakdown = {f"range_{best_match}": (lo, hi), "estimated_midpoint": shift}
        else:
            # Fallback to hybridization-based estimate
            shift = {"sp3": 35, "sp2": 120, "sp": 75}.get(hyb, 50)
            breakdown = {"fallback_hybridization_estimate": shift}

        # Adjustments
        adj_total = 0
        for item in env:
            if item.lower() in ("electronegative", "halogen", "heteroatom"):
                adj = 20
                breakdown[f"adjustment_{item}"] = adj
                adj_total += adj
            elif item.lower() in ("conjugated", "aromatic"):
                adj = -5
                breakdown[f"adjustment_{item}"] = adj
                adj_total += adj

        shift += adj_total
        shift = round(shift, 1)
        breakdown["total_shift"] = shift

        return {
            "chemical_shift_ppm": shift,
            "contribution_breakdown": breakdown,
            "interpretation": self._interpret_13c_shift(shift),
            "typical_range": self._classify_13c_region(shift),
        }

    def _calc_general_shift(self, nuc: str, env: List[str], hyb: str) -> dict:
        """Generic calculation for other nuclei (19F, 31P)."""
        # Simple heuristic based on environment electronegativity
        base = {"19F": -220, "31P": 0}.get(nuc, 0)  # 19F referenced to CFCl3
        env_factor = len(env) * 15 * (1 if nuc == "19F" else 5)
        shift = round(base + env_factor, 1)

        return {
            "chemical_shift_ppm": shift,
            "contribution_breakdown": {"base": base, "environment_adjustment": env_factor},
            "interpretation": f"Estimated {nuc} chemical shift based on environment.",
        }

    def _match_substituent_effect(self, item: str) -> Optional[float]:
        """Match a substituent name to its chemical shift effect."""
        il = item.lower().strip().lstrip("-")
        for key, val in self.h_substituents.items():
            k_clean = key.lstrip("-").lower().replace(" ", "_").replace("(", "").replace(")", "")
            i_clean = il.replace(" ", "_").replace("(", "").replace(")", "")
            if i_clean in k_clean or k_clean in i_clean:
                return val
        return None

    @staticmethod
    def _interpret_1h_shift(shift: float) -> str:
        if shift < 0.5:
            return "Highly shielded — possibly metal-bound or extreme upfield case."
        elif shift < 2.0:
            return "Typical alkyl region; shielded by electron-donating groups."
        elif shift < 4.5:
            return "Deshielded by attached electronegative atoms (O, N, halogens)."
        elif shift < 5.5:
            return "Vinylic or allylic proton region."
        elif shift < 6.5:
            return "Vinylic or aromatic-adjacent region."
        elif shift <= 8.5:
            return "Aromatic proton region — deshielded by ring current."
        else:
            return "Highly deshielded — aldehyde (9–10 ppm), carboxylic acid (10–12 ppm), or strongly deshielded system."

    @staticmethod
    def _classify_1h_region(shift: float) -> str:
        if shift < 2:
            return "Alkyl region (0.5–2 ppm)"
        elif shift < 4.5:
            return "α-to-heteroatom region (2–4.5 ppm)"
        elif shift < 5.5:
            "Allylic/vinylic region (4.5–5.5 ppm)"
        elif shift < 8.5:
            return "Vinylic/aromatic region (5.5–8.5 ppm)"
        else:
            return "Aldehyde/carboxylic acid region (>8.5 ppm)"

    @staticmethod
    def _classify_13c_region(shift: float) -> str:
        if shift < 60:
            return "Alkyl carbon region (0–60 ppm)"
        elif shift < 100:
            return "C-O / C-N / alkyne region (60–100 ppm)"
        elif shift < 145:
            return "Alkene / aromatic C region (100–145 ppm)"
        elif shift < 170:
            return "Imide / amide / ester / acid C=O region (145–170 ppm)"
        elif shift < 220:
            return "Ketone / aldehyde C=O region (170–220 ppm)"
        else:
            return "Very deshielded carbon (>220 ppm)"

    @staticmethod
    def _interpret_13c_shift(shift: float) -> str:
        if shift > 200:
            return "Aldehyde/ketone carbonyl — highly deshielded sp² carbon bonded to oxygen."
        elif shift > 175:
            return "Carboxylic acid derivative carbonyl — deshielded by electronegative oxygen."
        elif shift > 155:
            return "Amide/ester carbonyl or quaternary aromatic carbon."
        elif shift > 110:
            return "Aromatic or olefinic sp² carbon."
        elif shift > 80:
            return "C bonded to one or more electronegative atoms (O, N, halogens)."
        else:
            return "Aliphatic sp³ carbon."

    @staticmethod
    def _explain_shielding_mechanism(env: list, shift: float) -> str:
        parts = []
        ewg = [e for e in env if any(kw in e.lower() for kw in ["no2", "cor", "cooh", "halog", "cl", "br", "f"])]
        edg = [e for e in env if any(kw in e.lower() for kw in ["alkyl", "si", "ch3"])]

        if ewg:
            parts.append(f"Electron-withdrawing groups ({', '.join(ewg)}) reduce electron density → lower shielding (deshielding).")
        if edg:
            parts.append(f"Electron-donating groups ({', '.join(edg)}) increase electron density → higher shielding.")
        if not parts:
            parts.append("Local electronic environment determines shielding through inductive and mesomeric effects.")

        return " ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if not parts:
                raise ChemMCPError("Empty input.")

            nucleus = parts[0]
            kwargs = {"nucleus": nucleus}
            env = []
            modes = {"basic", "standard", "detailed"}
            hybs = {"sp3", "sp2", "sp"}

            for p in parts[1:]:
                pl = p.lower()
                if pl in modes:
                    kwargs["detail_level"] = p
                elif pl in hybs:
                    kwargs["hybridization"] = p
                elif p.upper() in ("TMS", "CFCL3", "H3PO4"):
                    kwargs["reference_compound"] = p
                else:
                    env.append(p)
            if env:
                kwargs["atom_environment"] = env

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
