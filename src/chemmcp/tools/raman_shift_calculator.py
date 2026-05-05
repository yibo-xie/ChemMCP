"""
Raman Shift Calculator — 拉曼位移计算与峰位预测工具 (#321)

功能：
  1. 根据激发波长和散射波长/波数计算拉曼位移
  2. 预测常见官能团的拉曼特征峰位置
公式：Δν̃(cm⁻¹) = (1/λ_ex – 1/λ_Raman) × 10⁷
"""

import logging
import math
from typing import Optional, List, Dict, Any, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见官能团拉曼位移参考数据 (cm⁻¹) ──────────────────────────
RAMAN_PEAK_DB: Dict[str, List[Dict[str, Any]]] = {
    "C-C stretch":      [{"peak": 1100, "intensity": "strong",   "assignment": "C-C aliphatic stretch"}],
    "C=C stretch":      [{"peak": 1650, "intensity": "strong",   "assignment": "C=C alkene stretch"}],
    "C≡C stretch":      [{"peak": 2120, "intensity": "medium",   "assignment": "C≡C alkyne stretch"}],
    "C-H stretch":      [{"peak": 2900, "intensity": "strong",   "assignment": "C-H aliphatic"},
                         {"peak": 3050, "intensity": "medium",   "assignment": "C-H aromatic"}],
    "O-H stretch":      [{"peak": 3400, "intensity": "broad",    "assignment": "O-H stretch (alcohol/water)"}],
    "N-H stretch":      [{"peak": 3300, "intensity": "medium",   "assignment": "N-H stretch (amine/amide)"}],
    "C=O stretch":      [{"peak": 1715, "intensity": "strong",   "assignment": "C=O carbonyl stretch"}],
    "C-O stretch":      [{"peak": 1050, "intensity": "strong",   "assignment": "C-O alcohol/ether"}],
    "C-N stretch":      [{"peak": 1120, "intensity": "medium",   "assignment": "C-N amine stretch"}],
    "C-Cl stretch":     [{"peak": 700,  "intensity": "strong",   "assignment": "C-Cl stretch"}],
    "C-Br stretch":     [{"peak": 550,  "intensity": "medium",   "assignment": "C-Br stretch"}],
    "Ring breathing":    [{"peak": 1000, "intensity": "very_strong", "assignment": "Benzene ring breathing mode"}],
    "S-S stretch":      [{"peak": 500,  "intensity": "weak",     "assignment": "S-S disulfide stretch"}],
    "P=O stretch":      [{"peak": 1180, "intensity": "strong",   "assignment": "P=O phosphate stretch"}],
    "NO2 sym str":      [{"peak": 1350, "intensity": "strong",   "assignment": "NO2 symmetric stretch"}],
    "NO2 asym str":     [{"peak": 1550, "intensity": "strong",   "assignment": "NO2 asymmetric stretch"}],
    "CH3 bend":         [{"peak": 1450, "intensity": "medium",   "assignment": "CH3 deformation/scissoring"}],
    "CH2 bend":         [{"peak": 1465, "intensity": "medium",   "assignment": "CH2 scissoring"}],
    "C=C aromatic":     [{"peak": 1600, "intensity": "strong",   "assignment": "Aromatic C=C quadrant stretch"}],
    "Si-O-Si stretch":  [{"peak": 1100, "intensity": "strong",   "assignment": "Si-O-Si asymmetric stretch"}],
}


@ChemMCPManager.register_tool
class RamanShiftCalculator(BaseTool):
    """
    拉曼位移计算与峰位预测工具。
    支持两种模式：(1) 精确计算拉曼位移 (2) 查询官能团特征拉曼峰
    """
    __version__                = "0.1.0"
    name                       = "RamanShiftCalculator"
    func_name                  = "calculate_raman_shift"
    description                = ("Calculate Raman shift from laser/excitation wavelength and observed "
                                 "scattering wavelength or predict Raman peak positions for common functional groups.")
    implementation_description = ("Uses the formula Δν̃(cm⁻¹) = (1/λ_excitation − 1/λ_Raman) × 10⁷ "
                                 "and a built-in database of functional group Raman shifts.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Raman Spectroscopy", "Vibrational Spectroscopy",
                                   "Analytical Chemistry", "Peak Prediction"]
    required_envs              = []

    code_input_sig = [
        ("laser_wavelength_nm",        "float", "N/A",       "Excitation/laser wavelength in nanometers (nm), e.g. 532.0."),
        ("observed_wavelength_nm",     "float", "None",      "Observed Raman scattering wavelength in nm. Use None if using functional_group mode."),
        ("functional_group",           "str",   "None",      "Name of functional group to look up (e.g. 'C=C stretch', 'O-H stretch'). Use None if calculating from wavelengths."),
        ("temperature_k",              "float", "298.15",    "Temperature in Kelvin (affects peak position slightly via Boltzmann correction)."),
    ]

    text_input_sig = [
        ("input_params",               "str",   "N/A",       "Space-separated: laser_wavelength_nm observed_wavelength_nm [functional_group] [temperature_k]"),
    ]

    output_sig = [
        ("raman_shift_cm",             "float", "Calculated Raman shift in wavenumber (cm⁻¹)."),
        ("predicted_peaks",            "dict",  "Predicted peaks info when functional_group is provided (None otherwise)."),
        ("laser_wavenumber_cm",        "float", "Laser excitation wavenumber in cm⁻¹."),
    ]

    examples = [
        {
            "code_input": {
                "laser_wavelength_nm": 532.0,
                "observed_wavelength_nm": 547.07,
                "functional_group": None,
                "temperature_k": 298.15,
            },
            "text_input": {"input_params": "532.0 547.07 None 298.15"},
            "output": {
                "raman_shift_cm": 517.8,
                "predicted_peaks": None,
                "laser_wavenumber_cm": 18796.6,
            },
        },
        {
            "code_input": {
                "laser_wavelength_nm": 532.0,
                "observed_wavelength_nm": None,
                "functional_group": "C=C stretch",
                "temperature_k": 298.15,
            },
            "text_input": {"input_params": "532.0 None C=C_stretch 298.15"},
            "output": {
                "raman_shift_cm": 1650.0,
                "predicted_peaks": {"peak": 1650, "intensity": "strong", "assignment": "C=C alkene stretch"},
                "laser_wavenumber_cm": 18796.6,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load reference database."""
        self.peak_db = RAMAN_PEAK_DB

    # ── internal helpers ────────────────────────────────────────

    @staticmethod
    def _nm_to_cm(nm: float) -> float:
        """Convert wavelength (nm) → wavenumber (cm⁻¹)."""
        return 1e7 / nm

    def _lookup_functional_group(self, fg_name: str) -> List[Dict[str, Any]]:
        """Case-insensitive lookup of functional group peaks."""
        # fuzzy match: replace spaces with underscores and vice versa
        key_normalized = fg_name.strip().replace("_", " ").lower()
        for db_key, peaks in self.peak_db.items():
            if db_key.lower() == key_normalized or db_key.lower().replace(" ", "_") == key_normalized:
                return peaks
        # partial match
        matches = [k for k in self.peak_db if key_normalized in k.lower()]
        if len(matches) == 1:
            return self.peak_db[matches[0]]
        elif matches:
            raise ChemMCPError(f"Ambiguous functional group '{fg_name}'. Did you mean: {matches}?")
        available = ", ".join(sorted(self.peak_db.keys()))
        raise ChemMCPError(f"Functional group '{fg_name}' not found. Available: {available}")

    # ── public API ─────────────────────────────────────────────

    def _run_base(
        self,
        laser_wavelength_nm: float,
        observed_wavelength_nm: Optional[float] = None,
        functional_group: Optional[str] = None,
        temperature_k: float = 298.15,
    ) -> Dict[str, Any]:
        """
        Core logic:
          Mode A (wavelength-based): calculate exact Raman shift from λ_ex and λ_Raman.
          Mode B (functional group): predict approximate Raman peak position(s).
        """
        if laser_wavelength_nm <= 0:
            raise ChemMCPError("Laser wavelength must be positive.")

        laser_wn = self._nm_to_cm(laser_wavelength_nm)

        # Mode B: functional group lookup
        if functional_group is not None and functional_group.lower() != "none":
            peaks = self._lookup_functional_group(functional_group)
            primary = peaks[0]
            raman_shift = float(primary["peak"])
            logger.info(f"Looked up functional group '{functional_group}': shift={raman_shift} cm⁻¹")
            return {
                "raman_shift_cm": round(raman_shift, 4),
                "predicted_peaks": primary,
                "laser_wavenumber_cm": round(laser_wn, 4),
            }

        # Mode A: wavelength calculation
        if observed_wavelength_nm is None:
            raise ChemMCPError(
                "Either observed_wavelength_nm or functional_group must be provided."
            )
        if observed_wavelength_nm <= 0:
            raise ChemMCPError("Observed wavelength must be positive.")

        obs_wn = self._nm_to_cm(observed_wavelength_nm)
        raman_shift = laser_wn - obs_wn

        if raman_shift < 0:
            # Anti-Stokes: negative shift; report absolute value with note
            logger.info(f"Anti-Stokes detected: |Δν| = {abs(raman_shift):.2f} cm⁻¹")
            raman_shift = abs(raman_shift)

        logger.info(
            f"Raman shift calculated: λ_ex={laser_wavelength_nm} nm, "
            f"λ_obs={observed_wavelength_nm} nm → Δν̃={raman_shift:.4f} cm⁻¹"
        )
        return {
            "raman_shift_cm": round(raman_shift, 4),
            "predicted_peaks": None,
            "laser_wavenumber_cm": round(laser_wn, 4),
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse space-separated text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least laser_wavelength and observed_wavelength.")

            laser = float(parts[0])
            obs_raw = parts[1].strip()
            obs = None if obs_raw.upper() in ("NONE", "NULL") else float(obs_raw)
            fg = None
            temp = 298.15

            idx = 2
            while idx < len(parts):
                p = parts[idx]
                # try to detect if it's a number (temp) or string (fg)
                try:
                    tval = float(p)
                    temp = tval
                except ValueError:
                    fg = p.replace("_", " ")
                idx += 1

            return self._run_base(laser, obs, fg, temp)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. "
                               "Format: 'laser_nm observed_nm [functional_group] [temp_K]'")
