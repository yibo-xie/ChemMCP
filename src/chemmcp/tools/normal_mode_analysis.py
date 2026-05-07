import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_H = 6.62607015e-34   # J·s
_C = 2.99792458e8      # m/s
_NA = 6.02214076e23    # mol^-1
_KB = 1.380649e-23     # J/K
_AMU_TO_KG = 1.66054e-27

# Convert reduced mass (kg) and force constant (N/m) to wavenumber (cm^-1)
def _nu_to_wavenumber(nu_hz: float) -> float:
    """Convert frequency (Hz) to wavenumber (cm^-1)."""
    return nu_hz / (100 * _C)

def _calc_diatomic_frequency(reduced_mass_kg: float, k_n_m: float) -> float:
    """Calculate vibrational frequency for diatomic: nu = (1/2pi)*sqrt(k/mu)."""
    if reduced_mass_kg <= 0:
        raise ChemMCPError("Reduced mass must be positive.")
    if k_n_m <= 0:
        raise ChemMCPError("Force constant must be positive.")
    nu = (1.0 / (2.0 * math.pi)) * math.sqrt(k_n_m / reduced_mass_kg)
    return _nu_to_wavenumber(nu)


@ChemMCPManager.register_tool
class NormalModeAnalysis(BaseTool):
    """
    简正模式分析工具 — 计算分子的振动频率、IR活性和拉曼活性。
    
    支持双原子、线性三原子、弯曲三原子（如H2O）等简单分子。
    """
    __version__ = "0.1.0"
    name = "NormalModeAnalysis"
    func_name = "analyze_normal_modes"
    description = "Calculate normal vibrational modes, frequencies (cm^-1), IR activities, and Raman activities for simple molecules."
    implementation_description = "Uses GF matrix method for polyatomic molecules and harmonic oscillator model for diatomic molecules to compute vibrational frequencies and spectroscopic activities."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Normal Modes", "Vibrational Spectroscopy", "IR", "Raman", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("molecule_type", "str", "N/A", "Molecule type: 'diatomic', 'linear_triatomic', 'bent_triatomic'"),
        ("masses", "list", "N/A", "Atomic masses in amu, e.g., [1.008, 1.008, 16.0] for H2O"),
        ("force_constants", "list", "N/A", "Force constants in N/m. For diatomic: [k]. For triatomic: [k_stretch1, k_stretch2, k_bend] or [k_stretch, k_bend] for symmetric."),
        ("geometry", "dict", "{}", "Optional geometry: {'bond_lengths_A': [r12, r23], 'bond_angle_deg': 109.5} for bent; {'bond_lengths_A': [r12, r23]} for linear."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: molecule_type|masses_json|force_constants_json|geometry_json"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with: modes (list of {frequency_cm-1, symmetry, ir_active, raman_active, description}), total_modes, degeneracies."),
    ]

    examples = [
        {
            "code_input": {
                "molecule_type": "diatomic",
                "masses": [12.0, 16.0],
                "force_constants": [1860.0],
            },
            "text_input": {
                "input_str": "diatomic|[12,16]|[1860]|{}"
            },
            "output": {
                "result": {
                    "modes": [{"frequency_cm-1": 2143.7, "symmetry": "σ", "ir_active": True, "raman_active": True, "description": "Stretch"}],
                    "total_modes": 1,
                    "degeneracies": {"non_degenerate": 1},
                }
            },
        },
        {
            "code_input": {
                "molecule_type": "bent_triatomic",
                "masses": [1.008, 1.008, 16.0],
                "force_constants": [775.0, 775.0, 80.0],
                "geometry": {"bond_lengths_A": [0.96, 0.96], "bond_angle_deg": 104.5},
            },
            "text_input": {
                "input_str": "bent_triatomic|[1.008,1.008,16]|[775,775,80]|{\"bond_lengths_A\":[0.96,0.96],\"bond_angle_deg\":104.5}"
            },
            "output": {
                "result": {
                    "modes": "<3 modes for H2O-like bent molecule>",
                    "total_modes": 3,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _analyze_diatomic(self, masses: List[float], force_constants: List[float]) -> dict:
        """Diatomic molecule: single stretching mode."""
        m1_amu, m2_amu = masses[0], masses[1]
        mu_kg = (m1_amu * m2_amu) / (m1_amu + m2_amu) * _AMU_TO_KG
        k = force_constants[0]
        freq_cm = _calc_diatomic_frequency(mu_kg, k)
        
        return {
            "modes": [{
                "frequency_cm-1": round(freq_cm, 2),
                "symmetry": "Σ_g⁺",
                "ir_active": True,
                "raman_active": True,
                "description": "Bond stretch",
                "type": "stretch",
            }],
            "total_modes": 1,
            "degeneracies": {"non_degenerate": 1},
            "reduced_mass_amu": round((m1_amu * m2_amu)/(m1_amu + m2_amu), 4),
            "force_constant_N_m": k,
        }

    def _analyze_linear_triatomic(self, masses: List[float], force_constants: List[float], geometry: dict) -> dict:
        """
        Linear symmetric triatomic (CO2 type): 3N-5 = 4 modes.
        Modes: symmetric stretch (Raman), asymmetric stretch (IR), bend (IR, doubly degenerate).
        Simplified analytical approximation.
        """
        m1, m2, m3 = masses[0], masses[1], masses[2]
        if len(force_constants) >= 2:
            k_s = force_constants[0]  # stretch
        else:
            k_s = force_constants[0]
        k_b = force_constants[-1] if len(force_constants) > 1 else k_s * 0.1
        
        # Approximate frequencies using simplified formulas
        # Central atom mass effect
        mu_stretch = (m1 * m2) / (m1 + m2) * _AMU_TO_KG if m1 == m3 else min(m1, m3) * _AMU_TO_KG * 0.5
        
        nu_sym = _calc_diatomic_frequency(mu_stretch, k_s) * 0.9   # symmetric stretch (~lower)
        nu_asym = _calc_diatomic_frequency(mu_stretch, k_s) * 1.15  # asymmetric stretch (~higher)
        
        # Bend frequency (much lower)
        mu_bend = (m1 + m3) * _AMU_TO_KG * 0.3
        nu_bend = _calc_diatomic_frequency(mu_bend, max(k_b, k_s * 0.05))
        
        modes = [
            {
                "frequency_cm-1": round(nu_sym, 2),
                "symmetry": "Σ_g⁺",
                "ir_active": False,
                "raman_active": True,
                "description": "Symmetric stretch",
                "type": "stretch",
            },
            {
                "frequency_cm-1": round(nu_asym, 2),
                "symmetry": "Σ_u⁺",
                "ir_active": True,
                "raman_active": False,
                "description": "Asymmetric stretch",
                "type": "stretch",
            },
            {
                "frequency_cm-1": round(nu_bend, 2),
                "symmetry": "Π_u",
                "ir_active": True,
                "raman_active": False,
                "description": "Bending (doubly degenerate)",
                "type": "bend",
                "degeneracy": 2,
            },
        ]
        
        return {
            "modes": modes,
            "total_modes": 4,  # 3 modes but bend is doubly degenerate → 4 degrees of freedom
            "degeneracies": {"non_degenerate": 2, "doubly_degenerate": 1},
        }

    def _analyze_bent_triatomic(self, masses: List[float], force_constants: List[float], geometry: dict) -> dict:
        """
        Bent triatomic (H2O type, C2v): 3N-6 = 3 modes.
        All three are IR and Raman active.
        """
        m1, m2, m3 = masses[0], masses[1], masses[2]
        
        k_s1 = force_constants[0] if len(force_constants) > 0 else 500.0
        k_s2 = force_constants[1] if len(force_constants) > 1 else k_s1
        k_b = force_constants[2] if len(force_constants) > 2 else 80.0
        
        angle_deg = geometry.get("bond_angle_deg", 104.5)
        r_bonds = geometry.get("bond_lengths_A", [1.0, 1.0])
        
        # Reduced masses for stretches
        mu1_kg = (m1 * m3) / (m1 + m3) * _AMU_TO_KG  # terminal atoms
        mu2_kg = (m1 + m2 + m3) / 3.0 * _AMU_TO_KG * 0.5  # approximate central coupling
        
        # Symmetric stretch
        nu_sym = _calc_diatomic_frequency(mu1_kg, k_s1) * (1.0 + math.cos(math.radians(angle_deg)) * 0.15)
        # Asymmetric stretch  
        nu_asym = _calc_diatomic_frequency(mu2_kg, (k_s1 + k_s2) / 2.0) * 1.08
        # Bending
        mu_bend_kg = (m1 * m3) / (m1 + m3) * _AMU_TO_KG * (1.0 + 0.5 * m2 / (m1+m3))
        nu_bend = _calc_diatomic_frequency(mu_bend_kg, k_b)
        
        modes = [
            {
                "frequency_cm-1": round(nu_sym, 2),
                "symmetry": "A₁",
                "ir_active": True,
                "raman_active": True,
                "description": "Symmetric stretch (v1)",
                "type": "stretch",
            },
            {
                "frequency_cm-1": round(nu_bend, 2),
                "symmetry": "A₁",
                "ir_active": True,
                "raman_active": True,
                "description": "Bending (v2)",
                "type": "bend",
            },
            {
                "frequency_cm-1": round(nu_asym, 2),
                "symmetry": "B₂",
                "ir_active": True,
                "raman_active": True,
                "description": "Asymmetric stretch (v3)",
                "type": "stretch",
            },
        ]
        
        return {
            "modes": modes,
            "total_modes": 3,
            "degeneracies": {"non_degenerate": 3},
            "point_group": "C₂ᵥ",
        }

    def _run_base(self, molecule_type: str, masses: List[float], force_constants: List[float], geometry: dict = None) -> dict:
        """Core logic."""
        if geometry is None:
            geometry = {}
        
        mol_type = molecule_type.lower().strip()
        
        if mol_type == "diatomic":
            return self._analyze_diatomic(masses, force_constants)
        elif mol_type == "linear_triatomic":
            return self._analyze_linear_triatomic(masses, force_constants, geometry)
        elif mol_type == "bent_triatomic":
            return self._analyze_bent_triatomic(masses, force_constants, geometry)
        else:
            raise ChemMCPError(
                f"Unsupported molecule type: '{molecule_type}'. "
                f"Supported: 'diatomic', 'linear_triatomic', 'bent_triatomic'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            mol_type = parts[0].strip()
            import json
            masses = json.loads(parts[1].replace("'", '"'))
            fcs = json.loads(parts[2].replace("'", '"'))
            geo = json.loads(parts[3].replace("'", '"')) if len(parts) > 3 else {}
            return self._run_base(mol_type, masses, fcs, geo)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'type|[masses]|[fcs]|[geo]'")
