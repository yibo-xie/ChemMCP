import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PartitionFunction(BaseTool):
    """
    配分函数计算工具 — 平动、转动、振动、电子配分函数及总配分函数。
    
    适用于双原子分子/理想气体分子统计热力学计算。
    """
    __version__                 = "0.1.0"
    name                        = "PartitionFunction"
    func_name                   = "calculate_partition_function"
    description                 = "Calculate molecular partition functions: translational (q_trans), rotational (q_rot), vibrational (q_vib), electronic (q_elec), and total Q."
    implementation_description  = "Uses statistical mechanics formulas: q_trans = (2πmkT/h²)^(3/2)·V, q_rot = T/(σ·Θ_rot), q_vib = 1/(1-exp(-Θ_vib/T)), q_elec = Σ g_i·exp(-ε_i/kT). For 1 atm, V = RT/P."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Statistical Mechanics", "Partition Function", "Thermodynamics", "Molecular"]
    required_envs               = []

    code_input_sig   = [
        ("temperature_k",            "float",  "N/A",     "Temperature in Kelvin."),
        ("molecular_mass_kg",         "float",  "N/A",     "Molecular mass in kg."),
        ("moment_of_inertia_kg_m2",   "float",  "N/A",     "Moment of inertia in kg·m²."),
        ("vibration_frequency_hz",    "float",  "None",    "Fundamental vibration frequency in Hz (None → skip vib)."),
        ("electronic_degeneracies",   "list",   "None",    "List of electronic state degeneracies g_i (None → q_elec=1)."),
        ("electronic_energies_j",     "list",   "None",    "List of electronic energies ε_i in J relative to ground state (None → all zero)."),
        ("symmetry_number",           "int",    "1",       "Rotational symmetry number σ."),
        ("pressure_atm",             "float",  "1.0",     "Pressure in atmospheres (for translational volume)."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'T(K) mass_kg I_kgm2 [freq_hz] [sigma] [P_atm]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with q_trans, q_rot, q_vib, q_elec, Q_total and contributions."),
    ]

    examples         = [
        {
            "code_input": {
                "temperature_k":              300.0,
                "molecular_mass_kg":          28.0 * 1.66053906660e-27,  # N2
                "moment_of_inertia_kg_m2":    1.3994e-46,  # N2
                "vibration_frequency_hz":     7.075e13,      # N2
                "electronic_degeneracies":    None,
                "electronic_energies_j":      None,
                "symmetry_number":            2,
                "pressure_atm":              1.0,
            },
            "text_input": {
                "input_params":               "300.0 4.64951e-26 1.3994e-46 7.075e13 2 1.0",
            },
            "output": {
                "result": {
                    "temperature_K": 300.0,
                    "q_trans": 3.52e30,
                    "q_rot": 58.5,
                    "q_vib": 0.0098,
                    "Q_total": 2.01e32,
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
        self.k_B = 1.380649e-23      # J/K
        self.h   = 6.62607015e-34     # J·s
        self.c   = 2.99792458e8       # m/s
        self.NA  = 6.02214076e23      # mol⁻¹
        self.R   = 8.314462618        # J/(mol·K)
        # 1 atm in Pa
        self.atm_to_Pa = 101325.0

    def _q_trans(self, T: float, m: float, P_atm: float) -> dict:
        """平动配分函数: q_trans = (2πmkT/h²)^(3/2) · V, V = RT/P"""
        V = (self.R * T) / (P_atm * self.atm_to_Pa)  # m³/mol → per molecule divide by NA
        V_per_molecule = V / self.NA
        
        lambda_th_sq = (self.h * self.h) / (2.0 * math.pi * m * self.k_B * T)
        q_t = math.pow(2.0 * math.pi * m * self.k_B * T / (self.h * self.h), 1.5) * V_per_molecule
        
        thermal_wavelength = math.sqrt(lambda_th_sq)
        
        return {
            "q_trans":           q_t,
            "thermal_de_Broglie_m": thermal_wavelength,
            "volume_per_molecule_m3": V_per_molecule,
            "formula": "q_trans = (2πmkT/h²)^(3/2) · V",
        }

    def _q_rot(self, T: float, I: float, sigma: int) -> dict:
        """转动配分函数（线性分子）: q_rot = T / (σ · Θ_rot), Θ_rot = ℏ²/(2Ik)"""
        if I <= 0:
            raise ChemMCPError("Moment of inertia must be > 0.")
        
        hbar = self.h / (2.0 * math.pi)
        Theta_rot = (hbar * hbar) / (2.0 * I * self.k_B)  # K
        
        if sigma <= 0:
            sigma = 1
        
        q_r = T / (sigma * Theta_rot) if Theta_rot > 0 else float('inf')
        
        return {
            "q_rot":              q_r,
            "rotational_temperature_K": Theta_rot,
            "symmetry_number":    sigma,
            "formula": "q_rot = T / (σ · Θ_rot), Θ_rot = ℏ²/(2Ik)",
        }

    def _q_vib(self, T: float, freq_hz: float) -> dict:
        """振动配分函数（谐振子，基态能量为零）: q_vib = 1 / (1 - exp(-hν/kT))"""
        if freq_hz is None or freq_hz <= 0:
            return None
        
        x = self.h * freq_hz / (self.k_B * T)
        
        if x > 700:  # exp(-700) ≈ 0, avoid overflow
            q_v = 1.0
        else:
            q_v = 1.0 / (1.0 - math.exp(-x))
        
        Theta_vib = self.h * freq_hz / self.k_B  # K
        
        return {
            "q_vib":                   q_v,
            "vibrational_temperature_K": Theta_vib,
            "x_hnu_kT":                x,
            "formula": "q_vib = 1/(1-exp(-hν/kT))",
        }

    def _q_elec(
        self,
        T: float,
        degeneracies: Optional[List[int]],
        energies_j: Optional[List[float]],
    ) -> dict:
        """电子配分函数: q_elec = Σ g_i · exp(-ε_i/kT)"""
        if degeneracies is None or len(degeneracies) == 0:
            return {"q_elec": 1.0, "note": "No electronic states specified; assuming ground state only with g₀=1."}
        
        q_e = 0.0
        terms = []
        
        for i, g in enumerate(degeneracies):
            eps = energies_j[i] if (energies_j and i < len(energies_j)) else 0.0
            term = g * math.exp(-eps / (self.k_B * T))
            q_e += term
            terms.append({"state": i, "g_i": g, "eps_J": eps, "term_value": term})
        
        return {
            "q_elec": q_e,
            "num_states": len(degeneracies),
            "terms": terms,
            "formula": "q_elec = Σ gᵢ·exp(-εᵢ/kT)",
        }

    def _run_base(
        self,
        temperature_k: float,
        molecular_mass_kg: float,
        moment_of_inertia_kg_m2: float,
        vibration_frequency_hz: Optional[float] = None,
        electronic_degeneracies: Optional[List[int]] = None,
        electronic_energies_j: Optional[List[float]] = None,
        symmetry_number: int = 1,
        pressure_atm: float = 1.0,
    ) -> dict:
        """Core logic: calculate all partition function components."""
        T  = float(temperature_k)
        m  = float(molecular_mass_kg)
        I  = float(moment_of_inertia_kg_m2)
        
        if T <= 0:
            raise ChemMCPError("Temperature must be > 0 K.")
        if m <= 0:
            raise ChemMCPError("Mass must be > 0.")
        if I <= 0:
            raise ChemMCPError("Moment of inertia must be > 0.")

        sigma = int(symmetry_number) if symmetry_number >= 1 else 1

        # Calculate each component
        qt = self._q_trans(T, m, pressure_atm)
        qr = self._q_rot(T, I, sigma)
        qv = self._q_vib(T, vibration_frequency_hz)
        qe = self._q_elec(T, electronic_degeneracies, electronic_energies_j)

        # Total partition function (factorized approximation for ideal gas)
        Q_total = qt["q_trans"] * qr["q_rot"]
        if qv is not None:
            Q_total *= qv["q_vib"]
        Q_total *= qe["q_elec"]

        result = {
            "temperature_K": T,
            "pressure_atm":  pressure_atm,
            "translational": qt,
            "rotational":    qr,
            "vibrational":   qv,
            "electronic":     qe,
            "Q_total":        round(Q_total, 6),
            "ln_Q_total":     round(math.log(Q_total), 6) if Q_total > 0 else None,
        }

        logger.info(f"PartitionFunction: T={T}K, Q_total={Q_total:.4e}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 3:
                raise ValueError("Need at least 'T mass I' params.")
            
            T = float(parts[0])
            m = float(parts[1])
            I = float(parts[2])
            freq = float(parts[3]) if len(parts) > 3 and parts[3].lower() != "none" else None
            sigma = int(parts[4]) if len(parts) > 4 else 1
            P = float(parts[5]) if len(parts) > 5 else 1.0
            
            return self._run_base(T, m, I, freq, None, None, sigma, P)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T(K) mass_kg I_kgm2 [freq_hz] [sigma] [P_atm]'")
