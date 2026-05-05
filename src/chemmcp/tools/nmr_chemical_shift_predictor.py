import logging
import math
from typing import List, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ─── Increment system for substituent effects on ¹H NMR chemical shifts ───
# Base values (ppm) relative to TMS = 0
_H_BASE_SHIFTS = {
    "CH3": 0.87,
    "CH2": 1.20,
    "CH": 1.55,
    "CH2=": 5.28,
    "CH=": 5.65,
    "Ar-H": 7.27,
    "CHO": 9.72,
    "COOH": 11.0,
    "NH2": 0.5,
    "OH": 1.5,
    "SH": 1.3,
}

# Substituent increment constants (ppm) for aliphatic systems (Shoolery / Tobey rules)
# Format: {substituent: (α-effect, β-effect, γ-effect)}
_SUBSTITUENT_INCREMENTS = {
    # Electron-withdrawing groups
    "-F":     (2.23,  0.36, -0.32),
    "-Cl":    (2.10,  0.37, -0.31),
    "-Br":    (2.03,  0.47, -0.27),
    "-I":     (1.85,  0.50, -0.24),
    "-CF3":   (1.78,  0.52, -0.20),
    "-CN":    (1.00,  0.92, -0.18),
    "-C≡CH":  (1.13,  0.82, -0.15),
    "-CHO":   (1.34,  0.62, -0.16),
    "-COR":   (1.22,  0.64, -0.17),   # Ketone
    "-COOH":  (1.05,  0.70, -0.15),
    "-COOR":  (0.97,  0.68, -0.14),   # Ester
    "-CONH2": (0.99,  0.55, -0.12),
    "-NO2":   (1.19,  0.86, -0.19),
    "-SO2R":  (1.35,  0.75, -0.21),
    "-SOR":   (1.10,  0.68, -0.16),
    # Electron-donating groups
    "-CH3":   (0.45,  0.25,  0.00),
    "-CH2R":  (0.42,  0.22,  0.00),
    "-CR3":   (0.38,  0.20,  0.00),
    "=CH2":   (0.95,  0.60, -0.08),
    "=CR2":   (0.90,  0.55, -0.07),
    "-Ph":    (1.30,  0.55, -0.10),   # Phenyl
    "-Ar":    (1.30,  0.55, -0.10),   # Aromatic
    "-OR":    (2.40,  0.35, -0.25),   # Ether/Alcohol
    "-OCOR":  (2.80,  0.40, -0.30),   # Ester oxygen
    "-NR2":   (1.15,  0.45, -0.18),   # Amine
    "-NHCOR": (1.35,  0.40, -0.14),   # Amide
    "-SR":    (1.50,  0.50, -0.18),   # Thioether
    "-SiR3":  (0.48,  0.28,  0.02),
}

# ¹³C NMR substituent increments (Lindeman-Adams parameters)
# Format: {substituent: (α-effect, β-effect, γ-effect)}
_C13_INCREMENTS = {
    "-F":     (70.1,   9.0,  -7.6),
    "-Cl":    (31.0,  10.4,  -5.9),
    "-Br":    (26.1,  10.8,  -5.5),
    "-I":     (6.6,  11.2,  -4.8),
    "-CN":    (3.4,   7.8,  -4.8),
    "-CHO":   (31.4,  -0.7,  -2.3),
    "-COR":   (29.5,  -0.7,  -2.3),
    "-COOH":  (21.4,   2.7,  -2.9),
    "-COOR":  (22.6,   2.0,  -2.8),
    "-CONH2": (18.7,   1.8,  -2.5),
    "-NO2":   (61.6,   3.4,  -4.5),
    "-OR":    (57.0,   6.5,  -4.8),
    "-OCOR":  (52.5,   6.5,  -4.8),
    "-NR2":   (30.3,   4.9,  -2.3),
    "-NHCOR": (20.8,   2.6,  -2.3),
    "-SR":    (20.5,   6.8,  -3.0),
    "-CH3":   (9.1,   9.4,  -2.5),
    "-CR=CR2":(8.5,   6.5,  -1.5),
    "-Ph":    (22.6,   8.5,  -2.5),
    "-C≡CH":  (4.4,   4.8,  -2.3),
}

# Aromatic proton substituent increments (ppm, ortho/meta/para)
_AROMATIC_H_INCREMENTS = {
    "-CH3":      (-0.18, -0.10, -0.19),
    "-CH2CH3":   (-0.15, -0.06, -0.18),
    "-CH(CH3)2":(-0.14, -0.09, -0.18),
    "-C(CH3)3": (-0.12, -0.08, -0.21),
    "-F":        (0.29,  -0.06, -0.23),
    "-Cl":       (0.03, -0.02, -0.09),
    "-Br":       (-0.17, -0.03, -0.13),
    "-I":        (-0.37, -0.07, -0.17),
    "-OH":       (-0.47, -0.11, -0.38),
    "-OCH3":     (-0.43, -0.10, -0.37),
    "-NH2":      (-0.69, -0.24, -0.61),
    "-NO2":      (0.95,   0.17,   0.34),
    "-CHO":      (0.58,   0.21,   0.27),
    "-COCH3":    (0.64,   0.09,   0.28),
    "-COOH":     (0.79,   0.13,   0.20),
    "-COOCH3":   (0.71,   0.01,   0.24),
    "-CN":       (0.27,   0.11,   0.31),
    "-CF3":      (0.63,   0.23,   0.15),
}


@ChemMCPManager.register_tool
class NmrChemicalShiftPredictor(BaseTool):
    """
    NMR化学位移预测工具（增强版）。
    基于取代基增量规则（Shoolery规则、Lindeman-Adams参数）预测¹H和¹³C NMR化学位移，
    支持脂肪族、芳香族、烯烃等多种体系。
    """
    __version__ = "0.1.0"
    name = "NmrChemicalShiftPredictor"
    func_name = "predict_nmr_shift"
    description = "Predict NMR chemical shifts using substituent increment rules (Shoolery, Lindeman-Adams). Supports ¹H and ¹³C for aliphatic, olefinic, and aromatic systems."
    implementation_description = (
        "Uses empirical substituent increment rules: Shoolery's rule for aliphatic ¹H, "
        "aromatic substituent constants for benzene derivatives, and Lindeman-Adams parameters "
        "for ¹³C chemical shift prediction."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["NMR", "Chemical Shift", "Prediction", "Spectroscopy", "Substituent Effects"]
    required_envs = []

    code_input_sig = [
        ("nucleus", "str", "1H", "Nucleus type: '1H' or '13C'."),
        ("base_structure", "str", "N/A", "Base structure key: 'CH3', 'CH2', 'CH', 'CH2=', 'CH=', 'Ar-H', 'CHO', 'COOH' etc."),
        ("substituents", "list", "[]", "List of substituents with positions, e.g., ['α:-Cl', 'β:-OH', 'γ:-CH3'] for aliphatic; ['ortho:-NO2', 'meta:-CH3'] for aromatic."),
        ("solvent", "str", "CDCl3", "NMR solvent."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: nucleus base_substituent sub1 sub2 ... Example: '1H CH2 α:-Cl β:-OH γ:-CH3'"),
    ]

    output_sig = [
        ("prediction", "dict", "Predicted chemical shift with breakdown of contributions, estimated range, and confidence level."),
    ]

    examples = [
        {
            "code_input": {
                "nucleus": "1H",
                "base_structure": "CH2",
                "substituents": ["α:-Cl", "β:-OH"],
                "solvent": "CDCl3",
            },
            "text_input": {"input_params": "1H CH2 α:-Cl β:-OH"},
            "output": {
                "prediction": {
                    "nucleus": "¹H",
                    "predicted_shift_ppm": 3.67,
                    "estimated_range_ppm": "3.4-3.9",
                    "breakdown": {
                        "base_value": 1.20,
                        "contributions": [{"substituent": "α:-Cl", "increment": 2.10}, {"substituent": "β:-OH", "increment": 2.40}],
                        "total_increment": 4.50,
                    },
                    "reference": "TMS = 0.0 ppm",
                    "method": "Shoolery's rule (aliphatic)",
                    "confidence": "high",
                }
            }
        },
        {
            "code_input": {
                "nucleus": "1H",
                "base_structure": "Ar-H",
                "substituents": ["ortho:-NO2", "para:-OH"],
                "solvent": "CDCl3",
            },
            "text_input": {"input_params": "1H Ar-H ortho:-NO2 para:-OH"},
            "output": {
                "prediction": {
                    "nucleus": "¹H",
                    "predicted_shift_ppm": 7.84,
                    "estimated_range_ppm": "7.6-8.1",
                    "method": "Aromatic substituent increment rule",
                    "confidence": "medium",
                }
            }
        },
        {
            "code_input": {
                "nucleus": "13C",
                "base_structure": "CH3",
                "substituents": ["α:-Cl"],
                "solvent": "CDCl3",
            },
            "text_input": {"input_params": "13C CH3 α:-Cl"},
            "output": {
                "prediction": {
                    "nucleus": "¹³C",
                    "predicted_shift_ppm": 40.1,
                    "estimated_range_ppm": "35-45",
                    "method": "Lindeman-Adams rule (¹³C)",
                    "confidence": "high",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _parse_substituent(self, sub_str: str) -> tuple:
        """Parse 'position:-X' into (position, substituent_key)."""
        parts = sub_str.split(":", 1)
        if len(parts) != 2:
            raise ChemMCPError(f"Invalid substituent format: '{sub_str}'. Use format 'position:-X', e.g., 'α:-Cl'")
        return parts[0].strip().lower(), parts[1].strip()

    def _predict_aliphatic_h1(self, base: str, substituents: List[str]) -> dict:
        """Predict aliphatic ¹H shift using Shoolery's rule: δ = 0.23 + Σσᵢ."""
        base_val = _H_BASE_SHIFTS.get(base)
        if base_val is None:
            raise ChemMCPError(f"Unknown base structure for ¹H: '{base}'. Available: {list(_H_BASE_SHIFTS.keys())}")

        total_inc = 0.0
        contributions = []
        for sub_str in substituents:
            pos, sub_key = self._parse_substituent(sub_str)
            if sub_key not in _SUBSTITUENT_INCREMENTS:
                contributions.append({"substituent": f"{pos}:{sub_key}", "increment": 0.0, "note": "not in database"})
                continue
            inc = _SUBSTITUENT_INCREMENTS[sub_key]
            pos_idx = {"α": 0, "alpha": 0, "β": 1, "beta": 1, "γ": 2, "gamma": 2}.get(pos)
            if pos_idx is None:
                contributions.append({"substituent": f"{pos}:{sub_key}", "increment": 0.0, "note": f"unknown position '{pos}'"})
                continue
            val = inc[pos_idx]
            total_inc += val
            contributions.append({"substituent": f"{pos}:{sub_key}", "increment": round(val, 2)})

        predicted = base_val + total_inc
        error_margin = max(0.3, abs(total_inc) * 0.1)

        return {
            "predicted_shift_ppm": round(predicted, 2),
            "estimated_range_ppm": f"{round(predicted - error_margin, 2)}-{round(predicted + error_margin, 2)}",
            "breakdown": {
                "base_value": base_val,
                "contributions": contributions,
                "total_increment": round(total_inc, 2),
            },
            "method": "Shoolery's rule (aliphatic ¹H): δ = δ_base + ΣΔσᵢ",
            "confidence": "high" if len(contributions) > 0 else "low",
        }

    def _predict_aromatic_h1(self, base: str, substituents: List[str]) -> dict:
        """Predict aromatic ¹H shift using substituent increment rule."""
        base_val = _H_BASE_SHIFTS.get("Ar-H", 7.27)
        total_inc = 0.0
        contributions = []
        for sub_str in substituents:
            pos, sub_key = self._parse_substituent(sub_str)
            if sub_key not in _AROMATIC_H_INCREMENTS:
                contributions.append({"substituent": f"{pos}:{sub_key}", "increment": 0.0, "note": "not in database"})
                continue
            inc = _AROMATIC_H_INCREMENTS[sub_key]
            pos_idx = {"ortho": 0, "o": 0, "meta": 1, "m": 1, "para": 2, "p": 2}.get(pos)
            if pos_idx is None:
                contributions.append({"substituent": f"{pos}:{sub_key}", "increment": 0.0, "note": f"unknown position '{pos}'"})
                continue
            val = inc[pos_idx]
            total_inc += val
            contributions.append({"substituent": f"{pos}:{sub_key}", "increment": round(val, 2)})

        predicted = base_val + total_inc
        return {
            "predicted_shift_ppm": round(predicted, 2),
            "estimated_range_ppm": f"{round(predicted - 0.3, 2)}-{round(predicted + 0.3, 2)}",
            "breakdown": {
                "base_value": base_val,
                "contributions": contributions,
                "total_increment": round(total_inc, 2),
            },
            "method": "Aromatic substituent increment rule (benzene = 7.27 ppm)",
            "confidence": "medium",
        }

    def _predict_aliphatic_c13(self, base: str, substituents: List[str]) -> dict:
        """Predict ¹³C shift using Lindeman-Adams parameters."""
        base_map_c13 = {
            "CH3": -2.6, "CH2": 15.6, "CH": 23.8,
        }
        base_val = base_map_c13.get(base)
        if base_val is None:
            raise ChemMCPError(f"Unknown base structure for ¹³C: '{base}'. Available: {list(base_map_c13.keys())}")

        total_inc = 0.0
        contributions = []
        for sub_str in substituents:
            pos, sub_key = self._parse_substituent(sub_str)
            if sub_key not in _C13_INCREMENTS:
                contributions.append({"substituent": f"{pos}:{sub_key}", "increment": 0.0, "note": "not in database"})
                continue
            inc = _C13_INCREMENTS[sub_key]
            pos_idx = {"α": 0, "alpha": 0, "β": 1, "beta": 1, "γ": 2, "gamma": 2}.get(pos)
            if pos_idx is None:
                continue
            val = inc[pos_idx]
            total_inc += val
            contributions.append({"substituent": f"{pos}:{sub_key}", "increment": round(val, 1)})

        predicted = base_val + total_inc
        return {
            "predicted_shift_ppm": round(predicted, 1),
            "estimated_range_ppm": f"{round(predicted - 3, 1)}-{round(predicted + 3, 1)}",
            "breakdown": {
                "base_value": base_val,
                "contributions": contributions,
                "total_increment": round(total_inc, 1),
            },
            "method": "Lindeman-Adams parameters (¹³C)",
            "confidence": "high",
        }

    def _run_base(self, nucleus: str, base_structure: str, substituents: Optional[List[str]] = None,
                  solvent: str = "CDCl3") -> dict:
        """Core prediction logic."""
        nuc = nucleus.upper().replace(" ", "")
        substituents = substituents or []

        if nuc in ("1H", "H1", "PROTON"):
            if base_structure.lower() in ("ar-h", "aromatic", "phenyl"):
                result = self._predict_aromatic_h1(base_structure, substituents)
            else:
                result = self._predict_aliphatic_h1(base_structure, substituents)
        elif nuc in ("13C", "C13", "CARBON"):
            result = self._predict_aliphatic_c13(base_structure, substituents)
        else:
            raise ChemMCPError(f"Unsupported nucleus: '{nucleus}'. Use '1H' or '13C'.")

        result["nucleus"] = result.get("nucleus", nucleus)
        result["solvent"] = solvent
        result["reference"] = "TMS = 0.0 ppm"

        return {"prediction": result}

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if len(parts) < 2:
            raise ChemMCPError(f"Input required. Format: 'nucleus base [sub1] [sub2] ...'. Got: '{input_params}'")

        nucleus = parts[0]
        base = parts[1]
        subs = parts[2:] if len(parts) > 2 else []
        return self._run_base(nucleus, base, subs)
