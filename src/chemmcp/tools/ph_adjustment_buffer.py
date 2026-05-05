"""
pH Adjustment & Buffer Preparation — pH调节缓冲液配制计算
Henderson-Hasselbalch 方程、酸碱滴定计算、缓冲液配制方案
"""
import logging
import math
from typing import Optional, List, Dict, Any, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见缓冲体系数据 ────────────────────────────────────────────────
BUFFER_SYSTEMS: Dict[str, dict] = {
    "acetate": {
        "weak_acid": "Acetic acid (CH₃COOH)",
        "conjugate_base": "Sodium acetate (CH₃COONa)",
        "pKa": 4.76,
        "mw_acid": 60.05, "mw_salt": 82.03,
        "effective_range": (3.76, 5.76),
        "notes": "Common for pH 4-5.5; compatible with most enzymes.",
    },
    "phosphate": {
        "weak_acid": "NaH₂PO₄ (sodium dihydrogen phosphate)",
        "conjugate_base": "Na₂HPO₄ (disodium hydrogen phosphate)",
        "pKa": 7.21,
        "mw_acid": 120.0, "mw_salt": 141.96,
        "effective_range": (6.21, 8.21),
        "notes": "Most common biological buffer; may inhibit some enzymes at high conc.",
    },
    "tris": {
        "weak_acid": "Tris base (tris(hydroxymethyl)aminomethane)",
        "conjugate_base": "Tris-HCl",
        "pKa": 8.06,
        "mw_acid": 121.14, "mw_salt": 157.6,
        "effective_range": (7.06, 9.06),
        "notes": "Very common in biochemistry; temperature-sensitive pKa (-0.031/°C).",
    },
    "carbonate": {
        "weak_acid": "NaHCO₃ (sodium bicarbonate)",
        "conjugate_base": "Na₂CO₃ (sodium carbonate)",
        "pKa": 10.33,
        "mw_acid": 84.01, "mw_salt": 105.99,
        "effective_range": (9.33, 11.33),
        "notes": "For alkaline range; releases CO₂; use fresh.",
    },
    "citrate": {
        "weak_acid": "Citric acid / NaH₂Citrate",
        "conjugate_base": "Na₂HCitrate / Na₃Citrate",
        "pKa": 3.13,  # using pKa2 for main buffer region
        "pKa_list": [3.13, 4.76, 6.40],
        "mw_acid": 210.14, "mw_salt": 294.10,
        "effective_range": (2.13, 6.40),
        "notes": "Wide range (pKa1-3); good chelating properties.",
    },
    "borate": {
        "weak_acid": "Boric acid (H₃BO₃)",
        "conjugate_base": "Borax (Na₂B₄O₇)",
        "pKa": 9.24,
        "mw_acid": 61.83, "mw_salt": 381.37,
        "effective_range": (8.24, 10.24),
        "notes": "Good for alkaline electrophoresis buffers.",
    },
    "ammonium": {
        "weak_acid": "NH₄Cl (ammonium chloride)",
        "conjugate_base": "NH₃·H₂O (aqueous ammonia)",
        "pKa": 9.25,
        "mw_acid": 53.49, "mw_salt": 17.03,
        "effective_range": (8.25, 10.25),
        "notes": "Volatile buffer; useful when buffer must be removed later.",
    },
    "formate": {
        "weak_acid": "Formic acid (HCOOH)",
        "conjugate_base": "Sodium formate (HCOONa)",
        "pKa": 3.75,
        "mw_acid": 46.03, "mw_salt": 68.01,
        "effective_range": (2.75, 4.75),
        "notes": "Good for acidic HPLC mobile phases; volatile.",
    },
    "succinate": {
        "weak_acid": "Succinic acid / NaHSuccinate",
        "conjugate_base": "Na₂Succinate",
        "pKa": 5.64,  # pKa2
        "pKa_list": [4.21, 5.64],
        "mw_acid": 118.09, "mw_salt": 162.12,
        "effective_range": (4.64, 6.64),
        "notes": "Good for pH ~5-6 mesophilic enzyme work.",
    },
    "hepes": {
        "weak_acid": "HEPES free acid",
        "conjugate_base": "HEPES sodium salt",
        "pKa": 7.55,
        "mw_acid": 238.30, "mw_salt": 260.31,
        "effective_range": (6.55, 8.55),
        "notes": "Good's buffer; minimal metal binding; cell culture compatible.",
    },
    "mes": {
        "weak_acid": "MES free acid",
        "conjugate_base": "MES sodium salt",
        "pKa": 6.15,
        "mw_ac": 195.24, "mw_salt": 213.14,
        "effective_range": (5.15, 7.15),
        "notes": "Good's buffer; excellent for pH 5.5-6.7; low membrane permeability.",
    },
}


@ChemMCPManager.register_tool
class PhAdjustmentBuffer(BaseTool):
    """
    pH 调节与缓冲液配制工具。
    基于 Henderson-Hasselbalch 方程：pH = pKa + log([A⁻]/[HA])
    支持缓冲液配制、pH调节所需酸/碱量计算、稀释计算。
    """
    __version__ = "0.1.0"
    name = "PhAdjustmentBuffer"
    func_name = "calculate_ph_buffer"
    description = "Calculate pH adjustment requirements and buffer preparation recipes using the Henderson-Hasselbalch equation."
    implementation_description = "Implements Henderson-Hasselbalch equation for buffer ratio calculation, acid/base titration for pH adjustment, and complete buffer preparation protocols with built-in data for common buffer systems."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["pH", "Buffer", "Henderson-Hasselbalch", "Titration", "Solution Preparation", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "buffer_prep", "Mode: 'buffer_prep', 'ph_adjust', 'dilution', 'system_select'."),
        ("target_ph", "float", "7.0", "Target pH value."),
        ("concentration_M", "float", "0.1", "Total buffer concentration in mol/L (M)."),
        ("volume_mL", "float", "1000", "Final buffer volume in mL."),
        ("buffer_system", "str", "phosphate", "Buffer system name (e.g., 'phosphate', 'acetate', 'tris')."),
        ("current_ph", "float", "7.0", "Current pH (for ph_adjust mode)."),
        ("solution_volume_mL", "float", "1000", "Solution volume to adjust (for ph_adjust mode)."),
        ("acid_type", "str", "hcl", "Acid type: 'hcl', 'h2so4', 'ch3cooh', etc."),
        ("base_type", "str", "naoh", "Base type: 'naoh', 'nahco3', 'nh3', etc."),
        ("acid_concentration_M", "float", "1.0", "Stock acid concentration (M)."),
        ("base_concentration_M", "float", "1.0", "Stock base concentration (M)."),
        ("temperature_C", "float", "25.0", "Temperature in °C (affects Tris pKa)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., 'buffer_prep 7.4 0.1 1000 phosphate' or 'ph_adjust 7.0 7.4 500 hcl 1.0'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with calculated amounts, recipe steps, and notes."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "buffer_prep",
                "target_ph": 7.4,
                "concentration_M": 0.1,
                "volume_mL": 1000,
                "buffer_system": "phosphate",
                "current_ph": 7.0,
                "solution_volume_mL": 1000,
                "acid_type": "hcl",
                "base_type": "naoh",
                "acid_concentration_M": 1.0,
                "base_concentration_M": 1.0,
                "temperature_C": 25.0,
            },
            "text_input": {
                "input_params": "buffer_prep 7.4 0.1 1000 phosphate",
            },
            "output": {
                "result": {
                    "mode": "buffer_prep",
                    "note": "Henderson-Hasselbalch calculation result.",
                }
            }
        },
        {
            "code_input": {
                "mode": "ph_adjust",
                "target_ph": 7.4,
                "current_ph": 7.0,
                "solution_volume_mL": 500,
                "acid_type": "",
                "base_type": "naoh",
                "base_concentration_M": 1.0,
                "volume_mL": 1000,
                "concentration_M": 0.1,
                "buffer_system": "phosphate",
                "temperature_C": 25.0,
            },
            "text_input": {
                "input_params": "ph_adjust 7.0 7.4 500 naoh 1.0",
            },
            "output": {
                "result": {
                    "mode": "ph_adjust",
                    "note": "Calculated amount of base needed to raise pH.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def henderson_hasselbalch(pKa: float, target_ph: float) -> tuple:
        """Calculate [A⁻]/[HA] ratio from HH equation."""
        ratio = 10 ** (target_ph - pKa)
        frac_A = ratio / (1 + ratio)
        frac_HA = 1 / (1 + ratio)
        return frac_A, frac_HA

    def _get_system(self, name: str) -> dict:
        key = name.lower().strip()
        if key not in BUFFER_SYSTEMS:
            available = ", ".join(sorted(BUFFER_SYSTEMS.keys()))
            raise ChemMCPError(
                f"Unknown buffer system '{name}'. Available: {available}"
            )
        return BUFFER_SYSTEMS[key]

    def _run_base(self, mode: str = "buffer_prep", target_ph: float = 7.0,
                  concentration_M: float = 0.1, volume_mL: float = 1000.0,
                  buffer_system: str = "phosphate", current_ph: float = 7.0,
                  solution_volume_mL: float = 1000.0, acid_type: str = "hcl",
                  base_type: str = "naoh", acid_concentration_M: float = 1.0,
                  base_concentration_M: float = 1.0,
                  temperature_C: float = 25.0) -> dict:

        m = mode.lower().strip()

        if m == "system_select":
            return self._select_system(target_ph)

        elif m == "buffer_prep":
            return self._calc_buffer_prep(
                target_ph, concentration_M, volume_mL, buffer_system, temperature_C)

        elif m in ("ph_adjust", "adjust"):
            return self._calc_ph_adjust(
                current_ph, target_ph, solution_volume_mL,
                acid_type, base_type, acid_concentration_M, base_concentration_M)

        elif m in ("dilution", "dilute"):
            return self._calc_dilution(concentration_M, volume_mL)

        else:
            raise ChemMCPError(
                f"Unknown mode '{mode}'. Use: 'buffer_prep', 'ph_adjust', "
                "'dilution', or 'system_select'."
            )

    def _select_system(self, target_ph: float) -> dict:
        """Recommend best buffer system(s) for a given pH."""
        suitable = []
        for name, info in BUFFER_SYSTEMS.items():
            lo, hi = info["effective_range"]
            if lo <= target_ph <= hi:
                distance = abs(target_ph - info["pKa"])
                suitable.append({
                    "system": name,
                    "pKa": info["pKa"],
                    "distance_from_pKa": round(distance, 2),
                    "effective_range": f"{lo:.2f} - {hi:.2f}",
                    "notes": info.get("notes", ""),
                })
        suitable.sort(key=lambda x: x["distance_from_pKa"])

        if not suitable:
            # Find closest
            all_sys = []
            for name, info in BUFFER_SYSTEMS.items():
                lo, hi = info["effective_range"]
                mid = (lo + hi) / 2
                all_sys.append((name, abs(mid - target_ph), info))
            all_sys.sort(key=lambda x: x[1])
            closest_name, _, closest_info = all_sys[0]
            return {"result": {
                "mode": "system_select",
                "target_ph": target_ph,
                "warning": f"No standard buffer system covers pH {target_ph}. Closest recommendation:",
                "recommended": closest_name,
                **closest_info,
            }}

        return {"result": {
            "mode": "system_select",
            "target_ph": target_ph,
            "suitable_systems": suitable[:5],
            "best_choice": suitable[0]["system"] if suitable else None,
            "note": "Choose the system whose pKa is closest to target pH (±1 unit optimal).",
        }}

    def _calc_buffer_prep(self, target_ph: float, conc: float, vol: float,
                          system_name: str, temp_C: float) -> dict:
        sys_info = self._get_system(system_name)
        pKa = sys_info["pKa"]

        # Temperature correction for Tris
        effective_pKa = pKa
        if system_name.lower() == "tris":
            effective_pKa = pKa - 0.031 * (temp_C - 25.0)

        frac_A, frac_HA = self.henderson_hasselbalch(effective_pKa, target_ph)

        total_mol = conc * vol / 1000  # L
        mol_A = total_mol * frac_A
        mol_HA = total_mol * frac_HA

        mw_salt = sys_info.get("mw_salt", 0)
        mw_acid = sys_info.get("mw_acid", 0)

        mass_A_g = mol_A * mw_salt if mw_salt else 0
        mass_HA_g = mol_HA * mw_acid if mw_acid else 0

        lo, hi = sys_info["effective_range"]

        return {"result": {
            "mode": "buffer_prep",
            "buffer_system": system_name,
            "target_ph": target_ph,
            "total_concentration_M": conc,
            "final_volume_mL": vol,
            "pKa_used": round(effective_pKa, 3),
            "temperature_correction_applied": system_name.lower() == "tris",
            "temperature_C": temp_C,
            "ratio_base_to_acid": round(frac_A / frac_HA, 3) if frac_HA > 0 else float('inf'),
            "fraction_conjugate_base_pct": round(frac_A * 100, 2),
            "fraction_weak_acid_pct": round(frac_HA * 100, 2),
            "moles_conjugate_base": round(mol_A, 6),
            "moles_weak_acid": round(mol_HA, 6),
            "mass_conjugate_base_g": round(mass_A_g, 4) if mass_A_g else "N/A",
            "mass_weak_acid_g": round(mass_HA_g, 4) if mass_HA_g else "N/A",
            "effective_range": f"{lo:.2f} - {hi:.2f}",
            "in_range": lo <= target_ph <= hi,
            "warning": "" if lo <= target_ph <= hi
            else f"⚠ Target pH {target_ph} outside effective range ({lo:.2f}-{hi:.2f}). Buffer capacity will be reduced.",
            "preparation_protocol": self._prep_protocol(system_name, mass_HA_g, mass_A_g, vol, target_ph),
            "notes": sys_info.get("notes", ""),
        }}

    def _calc_ph_adjust(self, current_ph: float, target_ph: float,
                        vol_mL: float, acid_type: str, base_type: str,
                        acid_conc: float, base_conc: float) -> dict:
        """
        Estimate volume of acid/base needed for pH adjustment.
        Simplified model: assumes dilute aqueous solution without significant buffering.
        For buffered solutions, this gives a rough estimate.
        """
        delta_ph = target_ph - current_ph
        abs_delta = abs(delta_ph)

        if abs_delta < 0.01:
            return {"result": {
                "mode": "ph_adjust",
                "current_ph": current_ph,
                "target_ph": target_ph,
                "delta_ph": delta_ph,
                "message": "pH is already at target. No adjustment needed.",
            }}

        # Rough estimate: for unbuffered water, ~0.05 mL of 1M NaOH per 100mL per pH unit
        # This is highly approximate — actual value depends on buffering capacity
        base_ml_per_ph_unit_per_100ml = 0.4  # empirical, order of magnitude
        acid_ml_per_ph_unit_per_100ml = 0.4

        if delta_ph > 0:
            # Need to RAISE pH → add base
            est_vol_mL = (abs_delta * base_ml_per_ph_unit_per_100ml *
                         (vol_mL / 100) / base_conc)
            reagent = base_type.upper()
            conc = base_conc
        else:
            # Need to LOWER pH → add acid
            est_vol_mL = (abs_delta * acid_ml_per_ph_unit_per_100ml *
                         (vol_mL / 100) / acid_conc)
            reagent = acid_type.upper() if acid_type else "HCl"
            conc = acid_conc

        return {"result": {
            "mode": "ph_adjust",
            "current_ph": current_ph,
            "target_ph": target_ph,
            "delta_ph": round(delta_ph, 3),
            "direction": "raise (add base)" if delta_ph > 0 else "lower (add acid)",
            "reagent": reagent,
            "reagent_concentration_M": conc,
            "estimated_volume_mL": round(est_vol_mL, 3),
            "estimated_volume_uL": round(est_vol_mL * 1000, 1),
            "solution_volume_mL": vol_mL,
            "caution": (
                "This is an ESTIMATE for weakly buffered solutions. "
                "Add gradually while monitoring with pH meter. "
                "Strongly buffered solutions will require more reagent."
            ),
            "protocol": self._adjust_protocol(current_ph, target_ph, est_vol_mL, reagent, conc),
        }}

    def _calc_dilution(self, conc: float, vol: float) -> dict:
        """Calculate stock solution preparation for making working concentration."""
        common_stocks = [0.5, 1.0, 2.0, 5.0, 10.0]
        results = []
        for stock_c in common_stocks:
            if stock_c > conc:
                dilution_factor = stock_c / conc
                stock_vol_for_1L = 1000 / dilution_factor
                results.append({
                    "stock_concentration_M": stock_c,
                    "dilution_factor": round(dilution_factor, 1),
                    "stock_volume_for_1L_mL": round(stock_vol_for_1L, 2),
                    "water_volume_for_1L_mL": round(1000 - stock_vol_for_1L, 2),
                })

        return {"result": {
            "mode": "dilution",
            "working_concentration_M": conc,
            "final_volume_mL": vol,
            "stock_options": results,
            "recommendation": f"For {conc}M working solution, prepare a {results[0]['stock_concentration_M']}M stock first." if results else "Concentration too high for typical stocks.",
        }}

    @staticmethod
    def _prep_protocol(sys_name: str, mass_acid: float, mass_base: float,
                       vol: float, target_ph: float) -> List[str]:
        return [
            f"1. Weigh {mass_acid:.3f} g of weak acid component.",
            f"2. Weigh {mass_base:.3f} g of conjugate base component.",
            f"3. Dissolve both in ~{vol*0.8:.0f} mL of deionized water.",
            f"4. Adjust final volume to {vol:.0f} mL with deionized water.",
            f"5. Verify pH with calibrated meter (target: {target_ph}). Fine-tune if needed.",
            f"6. Filter sterilize (0.22 μm) if required.",
            f"7. Store at 4°C; label with date, concentration, and pH.",
        ]

    @staticmethod
    def _adjust_protocol(curr_ph: float, tgt_ph: float, est_vol: float,
                         reagent: str, conc: float) -> List[str]:
        direction = "raise" if tgt_ph > curr_ph else "lower"
        return [
            f"1. Calibrate pH meter with fresh standard buffers (pH 4, 7, 10).",
            f"2. Place solution on magnetic stirrer; immerse pH electrode.",
            f"3. Prepare {reagent} ({conc} M) in a burette or graduated pipette.",
            f"4. Add reagent dropwise while stirring (~{est_vol:.1f} mL estimated total).",
            f"5. Monitor pH after each addition; slow down as approaching {tgt_ph}.",
            f"6. Stop at pH {tgt_ph}; record actual volume used.",
            f"7. If overshoot, back-titrate with opposite reagent.",
        ]

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            mode = parts[0]
            if mode == "buffer_prep":
                ph = float(parts[1]) if len(parts) > 1 else 7.0
                c = float(parts[2]) if len(parts) > 2 else 0.1
                v = float(parts[3]) if len(parts) > 3 else 1000
                sys = parts[4] if len(parts) > 4 else "phosphate"
                return self._run_base(mode, target_ph=ph, concentration_M=c,
                                       volume_mL=v, buffer_system=sys)
            elif mode in ("ph_adjust", "adjust"):
                curr = float(parts[1]) if len(parts) > 1 else 7.0
                tgt = float(parts[2]) if len(parts) > 2 else 7.4
                v = float(parts[3]) if len(parts) > 3 else 500
                reagent = parts[4] if len(parts) > 4 else "naoh"
                rc = float(parts[5]) if len(parts) > 5 else 1.0
                return self._run_base(mode, current_ph=curr, target_ph=tgt,
                                       solution_volume_mL=v, base_type=reagent,
                                       base_concentration_M=rc)
            elif mode == "system_select":
                ph = float(parts[1]) if len(parts) > 1 else 7.0
                return self._run_base(mode, target_ph=ph)
            elif mode in ("dilution", "dilute"):
                c = float(parts[1]) if len(parts) > 1 else 0.1
                v = float(parts[2]) if len(parts) > 2 else 1000
                return self._run_base(mode, concentration_M=c, volume_mL=v)
            else:
                raise ValueError(f"Unknown text mode: {mode}")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
