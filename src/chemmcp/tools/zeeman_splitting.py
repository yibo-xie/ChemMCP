import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ZeemanSplitting(BaseTool):
    """
    Zeeman 效应 — 外磁场中原子能级/谱线分裂分析。
    
    支持两种模式：
      - normal（正常 Zeeman 效应）：自旋 S=0，分裂为 3 条线 (π 和 σ±)
      - anomalous（反常 Zeeman 效应）：S≠0，使用 Landé g 因子
    """
    __version__                 = "0.1.0"
    name                        = "ZeemanSplitting"
    func_name                   = "analyze_zeeman_splitting"
    description                 = "Analyze Zeeman effect spectral line splitting in an external magnetic field (normal and anomalous)."
    implementation_description  = "Normal Zeeman: ΔE = μB·m·B, splits into 3 lines. Anomalous: ΔE = μB·gJ·mJ·B using Landé g-factor. Returns transition energies and wavelengths."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Zeeman Effect", "Atomic Physics", "Magnetic Field", "Spectroscopy"]
    required_envs               = []

    code_input_sig   = [
        ("magnetic_field_T",          "float",  "N/A",     "Magnetic field strength in Tesla."),
        ("orbital_angular_momentum_l","int",    "N/A",     "Orbital angular momentum quantum number L."),
        ("total_spin_s",             "float",  "N/A",     "Total spin quantum number S."),
        ("transition_wavelength_nm", "float",  "None",    "Transition wavelength in nm (optional, for spectral line analysis)."),
        ("lande_g_factor",           "float",  "None",    "Landé g-factor (if None, auto-calculated from L, S, J)."),
        ("total_angular_momentum_j", "float",  "None",    "Total angular momentum quantum number J. If None, computed as L+S."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'B(T) L S [wavelength_nm] [g_J] [J]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with zeeman_type, lande_g, splitting_components, transition_lines."),
    ]

    examples         = [
        {
            "code_input": {
                "magnetic_field_T":            1.0,
                "orbital_angular_momentum_l":  1,
                "total_spin_s":               1,
                "total_angular_momentum_j":    2.0,
                "transition_wavelength_nm":   None,
                "lande_g_factor":             None,
            },
            "text_input": {
                "input_params":               "1.0 1 1",
            },
            "output": {
                "result": {
                    "magnetic_field_T": 1.0,
                    "L": 1,
                    "S": 1,
                    "J": 2.0,
                    "zeeman_type": "anomalous",
                    "lande_g_factor": 1.25,
                    "splitting_components": [
                        {"mJ": -2, "delta_E_eV": -5.798e-4},
                        {"mJ": -1, "delta_E_eV": -2.899e-4},
                        {"mJ": 0,  "delta_E_eV": 0.0},
                        {"mJ": 1,  "delta_E_eV": 2.899e-4},
                        {"mJ": 2,  "delta_E_eV": 5.798e-4},
                    ],
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """物理常数"""
        self.mu_B = 9.2740100783e-24   # J/T (Bohr magneton)
        self.eV   = 1.602176634e-19     # J/eV
        self.h    = 6.62607015e-34       # J·s
        self.c    = 2.99792458e8         # m/s

    def _compute_lande_g(self, L: int, S: float, J: float) -> float:
        """Compute Landé g-factor: g_J = 1 + [J(J+1) + S(S+1) - L(L+1)] / [2J(J+1)]"""
        if J == 0:
            return 1.0  # limit case; actually undefined but conventionally set
        
        numerator = J * (J + 1) + S * (S + 1) - L * (L + 1)
        denominator = 2.0 * J * (J + 1)
        
        if abs(denominator) < 1e-15:
            return 1.0
        
        return 1.0 + numerator / denominator

    def _run_base(
        self,
        magnetic_field_T: float,
        orbital_angular_momentum_l: int,
        total_spin_s: float,
        transition_wavelength_nm: Optional[float] = None,
        lande_g_factor: Optional[float] = None,
        total_angular_momentum_j: Optional[float] = None,
    ) -> dict:
        """Core logic for Zeeman splitting."""
        B  = float(magnetic_field_T)
        L  = int(orbital_angular_momentum_l)
        S  = float(total_spin_s)
        
        if B < 0:
            raise ChemMCPError("Magnetic field must be non-negative.")
        if L < 0:
            raise ChemMCPError("L must be >= 0.")
        if S < 0:
            raise ChemMCPError("S must be >= 0.")

        # Determine J
        if total_angular_momentum_j is not None:
            J = float(total_angular_momentum_j)
        else:
            J = L + S  # default to maximum J (could also be |L-S| ... L+S)

        # Determine g-factor
        if lande_g_factor is not None:
            gJ = float(lande_g_factor)
        else:
            gJ = self._compute_lande_g(L, S, J)

        # Classify Zeeman type
        is_normal = (abs(S) < 1e-10)  # S ≈ 0 → normal Zeeman
        zeeman_type = "normal" if is_normal else "anomalous"

        # Generate splitting components
        components = []
        mJ_min = int(-abs(J))
        mJ_max = int(abs(J))
        
        for mJ in range(mJ_min, mJ_max + 1):
            delta_E_J = self.mu_B * gJ * mJ * B
            delta_E_eV = delta_E_J / self.eV
            
            comp = {
                "mJ": mJ,
                "delta_E_J": round(delta_E_J, 25),
                "delta_E_eV": round(delta_E_eV, 12),
            }
            
            # Add wavelength shift if transition wavelength given
            if transition_wavelength_nm is not None:
                wl = float(transition_wavelength_nm)
                E_photon = (self.h * self.c) / (wl * 1e-9)  # J
                E_photon_eV = E_photon / self.eV
                
                if E_photon > 0:
                    frac_shift = delta_E_eV / E_photon_eV
                    wl_shifted_nm = wl / (1.0 + frac_shift)
                    comp["wavelength_shifted_nm"] = round(wl_shifted_nm, 6)
                    comp["wavelength_delta_nm"] = round(wl_shifted_nm - wl, 6)
            
            components.append(comp)

        result = {
            "magnetic_field_T": B,
            "L": L,
            "S": S,
            "J": J,
            "zeeman_type": zeeman_type,
            "lande_g_factor": round(gJ, 8),
            "mu_B_J_per_T": self.mu_B,
            "splitting_components": components,
            "num_split_levels": len(components),
            "total_splitting_range_eV": round(
                abs(components[-1]["delta_E_eV"] - components[0]["delta_E_eV"]), 12
            ) if len(components) >= 2 else 0.0,
        }

        # Transition line analysis
        if transition_wavelength_nm is not None:
            wl = float(transition_wavelength_nm)
            lines = []
            for comp in components:
                if "wavelength_shifted_nm" in comp:
                    lines.append({
                        "mJ": comp["mJ"],
                        "wavelength_nm": comp["wavelength_shifted_nm"],
                        "polarization": self._get_polarization(mJ, is_normal),
                    })
            result["transition_wavelength_original_nm"] = wl
            result["zeeman_lines"] = lines

        logger.info(f"Zeeman {zeeman_type}: B={B}T, L={L}, S={S}, J={J}, gJ={gJ}, {len(components)} levels")
        return result

    def _get_polarization(self, mJ: int, is_normal: bool) -> str:
        """Determine polarization of Zeeman component."""
        if is_normal:
            if mJ == 0:
                return "π (Δm=0)"
            elif mJ > 0:
                return "σ⁺ (Δm=+1)"
            else:
                return "σ⁻ (Δm=-1)"
        else:
            if mJ == 0:
                return "π (Δm=0)"
            elif mJ > 0:
                return "σ⁺ (Δm>0)"
            else:
                return "σ⁻ (Δm<0)"

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need at least 'B L S' params.")
            
            B = float(parts[0])
            L = int(parts[1])
            S = float(parts[2])
            wl = float(parts[3]) if len(parts) > 3 and parts[3].lower() != "none" else None
            gJ = float(parts[4]) if len(parts) > 4 and parts[4].lower() != "none" else None
            J = float(parts[5]) if len(parts) > 5 and parts[5].lower() != "none" else None
            
            return self._run_base(B, L, S, wl, gJ, J)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'B(T) L S [wavelength_nm] [gJ] [J]'")
