import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants (SI)
_R_J = 8.314462618        # J/(mol·K)
_H = 6.62607015e-34       # J·s
_C = 2.99792458e8          # m/s
_KB = 1.380649e-23         # J/K
_NA = 6.02214076e23        # mol^-1
_AMU_TO_KG = 1.66054e-27   # kg

# Planck constant / Boltzmann constant ratio for rotational partition function
_H_OVER_KB = _H / _KB      # K·s


@ChemMCPManager.register_tool
class AdvancedPartitionFunction(BaseTool):
    """
    高级配分函数计算工具 — 计算分子的完整配分函数（平动、转动、振动、电子）。
    
    用于统计热力学计算，支持双原子和多原子分子。
    """
    __version__ = "0.1.0"
    name = "AdvancedPartitionFunction"
    func_name = "calculate_partition_function"
    description = "Calculate full molecular partition function: translational, rotational, vibrational, and electronic contributions with thermodynamic functions."
    implementation_description = "Computes q_trans via particle-in-box model, q_rot via rigid rotor (linear/non-linear), q_vib via harmonic oscillator, q_elec from ground state degeneracy. Derives U, H, S, G, A, Cv from partition functions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Partition Function", "Statistical Mechanics", "Thermodynamics", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("molecular_mass_amu", "float", "N/A", "Molecular mass in atomic mass units."),
        ("molecule_type", "str", "diatomic", "'diatomic' or 'nonlinear_polyatomic'"),
        ("volume_m3", "float", "0.02445", "System volume in m^3 (default: 1 atm, 298K → ~0.02445 m^3)."),
        ("rotational_constants_cm", "list", "[]", "Rotational constants [B] for diatomic or [A, B, C] for nonlinear in cm^-1. Empty=auto from geometry."),
        ("vibrational_frequencies_cm", "list", "[]", "Vibrational frequencies in cm^-1 (wavenumbers)."),
        ("electronic_degeneracy", "int", "1", "Ground state electronic degeneracy g0."),
        ("sigma", "int", "1", "Symmetry number (1 for heteronuclear, 2 for homonuclear diatomic)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: T|mass|type|volume|[rot_consts]|[vib_freqs]|g0|sigma"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with q_total, q_trans, q_rot, q_vib, q_elec, and derived thermodynamic quantities (U, S, Cv, H, G, A) per mole."),
    ]

    examples = [
        {
            "code_input": {
                "temperature_k": 298.15,
                "molecular_mass_amu": 28.0,
                "molecule_type": "diatomic",
                "rotational_constants_cm": [1.93],
                "vibrational_frequencies_cm": [2359.6],
                "electronic_degeneracy": 1,
                "sigma": 2,
            },
            "text_input": {
                "input_str": "298.15|28|diatomic|0.02445|[1.93]|[2359.6]|1|2"
            },
            "output": {
                "result": {
                    "q_trans": "<value>",
                    "q_rot": "<value>",
                    "q_vib": "<value>",
                    "q_elec": 1,
                    "q_total": "<value>",
                    "S_total_J_mol_K": "<value>",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _q_trans(self, T: float, M_amu: float, V: float) -> float:
        """Translational partition function: q_trans = (2*pi*M*kB*T/h^2)^(3/2) * V."""
        M_kg = M_amu * _AMU_TO_KG
        Lambda = _H / math.sqrt(2 * math.pi * M_kg * _KB * T)  # thermal de Broglie wavelength
        q_t = V / (Lambda ** 3)
        return q_t

    def _q_rot_diatomic(self, T: float, B_cm: float, sigma: int) -> float:
        """Rigid rotor partition function for linear molecule: q_rot = T / (sigma * theta_rot)."""
        if B_cm <= 0:
            return 1.0
        # theta_rot = hc*B / k_B  (in Kelvin), with B in cm^-1
        theta_rot = (_H * _C * 100 * B_cm) / _KB  # 100 converts cm^-1 to m^-1
        if T < 0.1 * theta_rot:
            # Low-T approximation: sum of first few J levels
            return 1.0 + 3.0 * math.exp(-2.0 * theta_rot / T) + 5.0 * math.exp(-6.0 * theta_rot / T)
        return T / (sigma * theta_rot)

    def _q_rot_nonlinear(self, T: float, rot_consts_cm: List[float], sigma: int) -> float:
        """Rotational partition function for nonlinear top."""
        A_cm, B_cm, C_cm = rot_consts_cm[0], rot_consts_cm[1], rot_consts_cm[2]
        theta_A = (_H * _C * 100 * A_cm) / _KB
        theta_B = (_H * _C * 100 * B_cm) / _KB
        theta_C = (_H * _C * 100 * C_cm) / _KB
        q_r = (math.sqrt(math.pi) / sigma) * (T ** 1.5) / math.sqrt(theta_A * theta_B * theta_C)
        return q_r

    def _q_vib(self, T: float, freqs_cm: List[float]) -> float:
        """Vibrational partition function (harmonic oscillator): q_vib = prod_i 1/(1-exp(-h*c*nu_i/kT))."""
        if not freqs_cm:
            return 1.0
        q_v = 1.0
        for nu_cm in freqs_cm:
            if nu_cm <= 0:
                continue
            x = (_H * _C * 100 * nu_cm) / (_KB * T)  # h*c*nu_bar / (k_B*T), dimensionless
            if x > 500:  # essentially zero-point only
                continue
            q_v *= 1.0 / (1.0 - math.exp(-x))
        return q_v

    def _q_elec(self, g0: int) -> float:
        """Electronic partition function (ground state only)."""
        return float(g0)

    def _thermodynamics_from_q(self, T: float, q_t: float, q_r: float, q_v: float, q_e: float,
                                 M_amu: float, molecule_type: str, rot_consts_cm: list,
                                 vib_freqs_cm: list, sigma: int) -> dict:
        """Derive thermodynamic quantities from partition functions."""
        q_tot = q_t * q_r * q_v * q_e
        
        # Translational contributions to U, Cv, S (per mole)
        # U_trans = 3/2 * R * T
        U_trans = 1.5 * _R_J * T
        Cv_trans = 1.5 * _R_J
        
        # S_trans = R * ln(q_e/N_A) + 5/2*R  (Sackur-Tetrode, simplified)
        # Using: S_trans/R = ln(q_trans/N_A) + 5/2
        N_A = _NA
        S_trans = _R_J * (math.log(max(q_t / N_A, 1e-300)) + 2.5)
        
        # Rotational contributions
        if molecule_type == "diatomic" and rot_consts_cm:
            U_rot = _R_J * T  # 2 degrees of freedom → RT
            Cv_rot = _R_J
            # S_rot = R * (ln(q_rot) + 1) for linear at high T
            S_rot = _R_J * (math.log(max(q_r, 1e-300)) + 1.0)
        elif molecule_type == "nonlinear_polyatomic" and len(rot_consts_cm) >= 3:
            U_rot = 1.5 * _R_J * T  # 3 degrees of freedom
            Cv_rot = 1.5 * _R_J
            S_rot = _R_J * (math.log(max(q_r, 1e-300)) + 1.5)
        else:
            U_rot = 0.0
            Cv_rot = 0.0
            S_rot = 0.0
        
        # Vibrational contributions (Einstein model per mode)
        U_vib = 0.0
        Cv_vib = 0.0
        S_vib = 0.0
        for nu_cm in vib_freqs_cm:
            if nu_cm <= 0:
                continue
            x = (_H * _C * 100 * nu_cm) / (_KB * T)
            if x > 500:
                continue
            exp_x = math.exp(-x)
            exp_x_1 = 1.0 - exp_x
            
            # U_vib_mode = R * theta_v * (1/2 + 1/(exp(x)-1)), but we report thermal part
            # U_vib_thermal = R * T * x / (exp(x)-1)
            U_vib += _R_J * T * x / (math.exp(x) - 1.0)
            
            # Cv_vib_mode = R * x^2 * exp(x) / (exp(x)-1)^2
            ex = math.exp(x)
            Cv_vib += _R_J * (x ** 2) * ex / ((ex - 1.0) ** 2)
            
            # S_vib_mode = R * (x/(exp(x)-1) - ln(1-exp(-x)))
            S_vib += _R_J * (x / (ex - 1.0) - math.log(max(exp_x_1, 1e-300)))
        
        # Totals
        U_total = U_trans + U_rot + U_vib  # per mole, excluding ZPE
        H_total = U_total + _R_J * T  # H = U + pV = U + RT for ideal gas
        Cv_total = Cv_trans + Cv_rot + Cv_vib
        S_total = S_trans + S_rot + S_vib
        G_total = H_total - T * S_total
        A_total = U_total - T * S_total
        
        return {
            "q_total": round(q_tot, 4),
            "q_trans": round(q_t, 4),
            "q_rot": round(q_r, 4),
            "q_vib": round(q_v, 4),
            "q_elec": q_e,
            "U_J_mol": round(U_total, 2),
            "H_J_mol": round(H_total, 2),
            "Cv_J_mol_K": round(Cv_total, 2),
            "S_J_mol_K": round(S_total, 2),
            "G_J_mol": round(G_total, 2),
            "A_J_mol": round(A_total, 2),
            "U_trans_J_mol": round(U_trans, 2),
            "U_rot_J_mol": round(U_rot, 2),
            "U_vib_J_mol": round(U_vib, 2),
            "S_trans_J_mol_K": round(S_trans, 2),
            "S_rot_J_mol_K": round(S_rot, 2),
            "S_vib_J_mol_K": round(S_vib, 2),
        }

    def _run_base(self, temperature_k: float, molecular_mass_amu: float, molecule_type: str = "diatomic",
                  volume_m3: float = 0.02445, rotational_constants_cm: List[float] = None,
                  vibrational_frequencies_cm: List[float] = None, electronic_degeneracy: int = 1,
                  sigma: int = 1) -> dict:
        """Core logic."""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive.")
        if molecular_mass_amu <= 0:
            raise ChemMCPError("Molecular mass must be positive.")
        
        if rotational_constants_cm is None:
            rotational_constants_cm = []
        if vibrational_frequencies_cm is None:
            vibrational_frequencies_cm = []
        
        # Calculate each component
        q_t = self._q_trans(temperature_k, molecular_mass_amu, volume_m3)
        
        mol_type = molecule_type.lower().strip()
        if mol_type == "diatomic":
            B_cm = rotational_constants_cm[0] if rotational_constants_cm else 1.0
            q_r = self._q_rot_diatomic(temperature_k, B_cm, sigma)
        elif mol_type == "nonlinear_polyatomic":
            if len(rotational_constants_cm) < 3:
                raise ChemMCPError("Nonlinear polyatomic requires 3 rotational constants [A, B, C] in cm^-1.")
            q_r = self._q_rot_nonlinear(temperature_k, rotational_constants_cm, sigma)
        else:
            raise ChemMCPError(f"Unsupported molecule type: '{molecule_type}'. Use 'diatomic' or 'nonlinear_polyatomic'.")
        
        q_v = self._q_vib(temperature_k, vibrational_frequencies_cm)
        q_e = self._q_elec(electronic_degeneracy)
        
        # Derive thermodynamic quantities
        thermo = self._thermodynamics_from_q(
            temperature_k, q_t, q_r, q_v, q_e,
            molecular_mass_amu, mol_type, rotational_constants_cm,
            vibrational_frequencies_cm, sigma
        )
        
        logger.info(f"Partition function calculated: q_tot={thermo['q_total']}")
        return thermo

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            T = float(parts[0].strip())
            mass = float(parts[1].strip())
            mol_type = parts[2].strip()
            V = float(parts[3].strip()) if len(parts) > 3 else 0.02445
            import json
            rot = json.loads(parts[4]) if len(parts) > 4 else []
            vib = json.loads(parts[5]) if len(parts) > 5 else []
            g0 = int(parts[6]) if len(parts) > 6 else 1
            sig = int(parts[7]) if len(parts) > 7 else 1
            return self._run_base(T, mass, mol_type, V, rot, vib, g0, sig)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T|mass|type|volume|[rot]|[vib]|g0|sigma'")
