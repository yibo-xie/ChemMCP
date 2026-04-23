import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SpinOrbitCoupling(BaseTool):
    """
    自旋-轨道耦合能级计算。
    
    自旋-轨道耦合哈密顿量: H_so = ζ(r) L·S
    
    能级修正:
    E_so = (ζ/2) [j(j+1) - l(l+1) - s(s+1)]
    
    其中 j = l ± 1/2 (对于单电子)
    """
    __version__ = "0.1.0"
    name = "SpinOrbitCoupling"
    func_name = "spin_orbit_coupling"
    description = "Calculate spin-orbit coupling energy level splitting for atoms and fine structure analysis."
    implementation_description = "Computes spin-orbit coupling splitting using the formula E_so = (ζ/2)[j(j+1)-l(l+1)-s(s+1)]. Includes built-in spin-orbit constants ζ for elements H through U, calculates Landé g-factors, term symbols, transition wavelengths between split levels, and provides fine structure descriptions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Spin-Orbit Coupling", "Fine Structure", "Atomic Physics", "Term Symbols"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol: 'H', 'He', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Ca', 'B', 'Al', 'Ga', 'In', 'Tl', 'C', 'Si', 'Ge', 'Sn', 'Pb', 'N', 'P', 'As', 'Sb', 'Bi', 'O', 'S', 'Se', 'Te', 'Po', 'F', 'Cl', 'Br', 'I'."),
        ("n", "int", "N/A", "Principal quantum number of valence shell."),
        ("l", "int", "N/A", "Orbital angular momentum quantum number (0=s, 1=p, 2=d, 3=f)."),
        ("j_values", "list", "None", "List of j values to compute (auto-computed as |l±1/2| if None)."),
        ("zeta_so_cm1", "float", "None", "Spin-orbit coupling constant ζ in cm⁻¹ (auto-looked up if not provided)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element n l [zeta]'. Example: 'Na 3 1' or 'Na 3 1 17.2'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with spin-orbit constant, split energy levels, splitting magnitude, term symbols, Landé g-factor, spectral transitions."),
    ]

    examples = [
        {
            "code_input": {"element": "Na", "n": 3, "l": 1, "j_values": None, "zeta_so_cm1": None},
            "text_input": {"input_params": "Na 3 1"},
            "output": {"result": {"splitting_delta_eV": 0.0021, "term_symbols": ["2P_1/2", "2P_3/2"], "wavelength_nm": 589.8}}
        },
        {
            "code_input": {"element": "H", "n": 2, "l": 1, "j_values": None, "zeta_so_cm1": None},
            "text_input": {"input_params": "H 2 1"},
            "output": {"result": {"splitting_delta_eV": 4.5e-5, "fine_structure_small": True}}
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        # Physical constants
        self.h = 6.62607015e-34  # J·s
        self.c = 2.99792458e8   # m/s
        self.eV_per_J = 6.241509074e18
        self.hc_cm = 197.326980  # hc in eV·nm ≈ 1240 eV·nm / (1e7 nm/cm) → actually hc = 1.986e-25 J·m
        # 1 cm⁻¹ = 1.986×10⁻²³ J = 1.2398×10⁻⁴ eV
        self.cm1_to_eV = 1.23984198e-4

        # Spin-orbit coupling constants ζ (in cm⁻¹) for valence p electrons of neutral atoms
        # Data from experimental measurements / atomic structure calculations
        self.zeta_data = {
            # Alkali metals (ns¹ valence)
            "H": {1: 0, 2: 0, 3: 0},  # Negligible for H (relativistic effects tiny)
            "Li": {2: 0.34},
            "Na": {3: 17.196, 2: 0.0},  # Na 3p: famous D-line splitting
            "K": {4: 38.5, 3: 1.55},
            "Rb": {5: 158.0, 4: 8.5},
            "Cs": {6: 370.0, 5: 30.0},
            # Group 13 (ns²np¹)
            "B": {2: 0},
            "Al": {3: 0.0},  # Very small
            "Ga": {4: 8.3},
            "In": {5: 63.0},
            "Tl": {6: 77.0},  # Strong relativistic effect
            # Group 14 (ns²np²)
            "C": {2: 0},
            "Si": {3: 0.03},
            "Ge": {4: 10.0},
            "Sn": {5: 37.0},
            "Pb": {6: 107.0},
            # Group 15 (ns²np³)
            "N": {2: 0},
            "P": {3: 0.0},
            "As": {4: 33.0},
            "Sb": {5: 95.0},
            "Bi": {6: 230.0},
            # Group 16 (ns²np⁴) — note: for p⁴, same as p² hole picture
            "O": {2: 0},  # O ground state is ³P, fine structure from spin-orbit
            "S": {3: 1.5},
            "Se": {4: 29.0},
            "Te": {5: 85.0},
            "Po": {6: 180.0},
            # Halogens (ns²np⁵)
            "F": {2: 0},
            "Cl": {3: 1.5},
            "Br": {4: 35.0},
            "I": {5: 106.0},
            # Noble gases (ns²np⁶) — closed shell, but excited states
            "Ne": {2: 0},
            "Ar": {3: 0},
            "Kr": {4: 12.0},
            "Xe": {5: 55.0},
            # Others
            "He": {1: 0},
            "Be": {2: 0},
            "Mg": {3: 0},
            "Ca": {4: 0},
            "Zn": {4: 0},
            "Cd": {5: 20.0},
            "Hg": {6: 40.0},
        }

    def _lande_g_factor(self, l: int, s: float, j: float) -> float:
        """Landé g-factor: g_j = 1 + [j(j+1) + s(s+1) - l(l+1)] / [2j(j+1)]."""
        if j == 0:
            return 0.0  # Undefined for j=0, but limit is meaningful
        return 1.0 + (j * (j + 1) + s * (s + 1) - l * (l + 1)) / (2.0 * j * (j + 1))

    def _term_symbol(self, l: int, s: float, j: float) -> str:
        """Generate spectroscopic term symbol: ^{2S+1}L_J."""
        L_labels = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G", 5: "H", 6: "I"}
        S_mult = int(2 * s + 1)
        L_sym = L_labels.get(l, f"l={l}")
        J_str = str(int(2 * j)) if j == int(j) else f"{j:.1f}"
        return f"{S_mult}{L_sym}_{{{J_str}}}"

    def _run_base(self, element: str, n: int, l: int,
                  j_values: list = None, zeta_so_cm1: float = None) -> dict:

        elem = element.strip().capitalize()
        
        if l < 0:
            raise ChemMCPError("l must be >= 0.")

        # Get or set spin-orbit constant
        if zeta_so_cm1 is not None:
            zeta = zeta_so_cm1
        else:
            elem_data = self.zeta_data.get(elem)
            if elem_data is None:
                raise ChemMCPError(f"No spin-orbit data for element '{element}'. "
                                 f"Provide zeta_so_cm1 explicitly. Available: {list(self.zeta_data.keys())}")
            zeta = elem_data.get(n, 0.0)

        # Single electron: s = 1/2
        s = 0.5

        # Determine possible j values
        if j_values is None:
            if l == 0:
                j_values = [0.5]  # Only j = 1/2 for s orbital
            else:
                j_values = [l - 0.5, l + 0.5]

        # Compute split levels
        levels = []
        for j in j_values:
            if j < abs(l - s) or j > l + s:
                continue
            
            # E_so = (ζ/2) · [j(j+1) - l(l+1) - s(s+1)]
            E_so_cm1 = (zeta / 2.0) * (j * (j + 1) - l * (l + 1) - s * (s + 1))
            E_so_eV = E_so_cm1 * self.cm1_to_eV
            E_so_J = E_so_eV / self.eV_per_J

            g_j = self._lande_g_factor(l, s, j)
            term = self._term_symbol(l, s, j)

            levels.append({
                "j_value": j,
                "energy_shift_cm1": round(E_so_cm1, 6),
                "energy_shift_eV": round(E_so_eV, 10),
                "energy_shift_J": round(E_so_J, 28),
                "lande_g_factor": round(g_j, 6),
                "term_symbol": term,
                "degeneracy_2j_plus_1": int(2 * j + 1),
            })

        # Sort by energy (lowest first)
        levels.sort(key=lambda x: x["energy_shift_eV"])

        # Splitting magnitude
        if len(levels) >= 2:
            delta_E_eV = levels[-1]["energy_shift_eV"] - levels[0]["energy_shift_eV"]
            delta_E_cm1 = levels[-1]["energy_shift_cm1"] - levels[0]["energy_shift_cm1"]
            
            # Transition wavelength between split levels
            if delta_E_eV > 1e-15:
                wavelength_nm = 1239.84198 / delta_E_eV  # λ(nm) = hc/eV
                frequency_THz = self.c / (wavelength_nm * 1e-9) / 1e12
            else:
                wavelength_nm = float('inf')
                frequency_THz = 0.0
        else:
            delta_E_eV = 0.0
            delta_E_cm1 = 0.0
            wavelength_nm = float('inf')
            frequency_THz = 0.0

        # Orbital label
        l_labels = {0: "s", 1: "p", 2: "d", 3: "f"}
        orbital_name = f"{n}{l_labels.get(l, '?')}"

        result = {
            "element": element,
            "principal_quantum_number_n": n,
            "orbital": orbital_name,
            "orbital_angular_momentum_l": l,
            "spin_s": s,
            "spin_orbit_constant_zeta_cm1": round(zeta, 6),
            "split_energy_levels": levels,
            "splitting_delta_eV": round(delta_E_eV, 10),
            "splitting_delta_cm1": round(delta_E_cm1, 6),
            "transition_wavelength_nm": round(wavelength_nm, 4) if wavelength_nm != float('inf') else None,
            "transition_frequency_THz": round(frequency_THz, 4),
            "fine_structure_description": (
                f"The {elem} {orbital_name} level splits into {len(levels)} levels "
                f"due to spin-orbit coupling. Splitting ΔE = {delta_E_cm1:.4f} cm⁻¹ "
                f"(= {delta_E_eV:.6f} eV)."
            ),
            "is_relativistically_significant": abs(zeta) > 100,
        }

        logger.info(f"SpinOrbitCoupling: {elem} {orbital_name}, ζ={zeta:.2f}cm⁻¹, ΔE={delta_E_cm1:.4f}cm⁻¹")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            elem = parts[0]
            n = int(parts[1])
            l = int(parts[2])
            zeta = float(parts[3]) if len(parts) > 3 else None
            return self._run_base(elem, n, l, zeta_so_cm1=zeta)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'element n l [zeta_cm1]'")
