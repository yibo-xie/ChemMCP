import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class StatisticalEntropy(BaseTool):
    """
    统计熵计算工具 — S = k_B · ln W。
    
    支持多种计算模式：
      1. 直接给定微观状态数 W（Boltzmann 熵公式）
      2. N 个粒子分配到 g 个能级/状态（W = g^N / (n₁! n₂! ...), 使用 Stirling 近似）
      3. Sackur-Tetrode 方程（理想单原子气体平动熵）
    """
    __version__                 = "0.1.0"
    name                        = "StatisticalEntropy"
    func_name                   = "calculate_statistical_entropy"
    description                 = "Calculate statistical entropy S = k_B·ln W. Supports direct W, particle occupation numbers, and Sackur-Tetrode equation for ideal gases."
    implementation_description  = "Mode 1: S=k·ln(W). Mode 2: W=N!/(Πnᵢ!) with Stirling approx. Mode 3: Sackur-Tetrode for monatomic ideal gas translational entropy."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Entropy", "Statistical Mechanics", "Thermodynamics", "Microstates", "Sackur-Tetrode"]
    required_envs               = []

    code_input_sig   = [
        ("mode",                     "str",    "direct",   "Calculation mode: 'direct', 'occupation', or 'sackur_tetrode'."),
        ("microstate_count_W",       "int",    "None",    "[mode=direct] Total number of microstates W."),
        ("number_of_particles_N",    "int",    "None",    "[mode=direct|occupation] Number of particles N."),
        ("occupation_numbers",       "list",   "None",    "[mode=occupation] List of occupation numbers [n1, n2, ...]."),
        ("temperature_k",            "float",  "None",    "[mode=sackur_tetrode] Temperature in Kelvin."),
        ("molecular_mass_kg",        "float",  "None",    "[mode=sackur_tetrode] Molecular mass in kg."),
        ("volume_m3",                "float",  "None",    "[mode=sackur_tetrode] Volume in m³ (default from NkT/P at 1atm)."),
        ("pressure_atm",             "float",  "1.0",     "[mode=sackur_tetrode] Pressure in atm (for default volume)."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Mode-specific format:\n"
                                          "  direct: 'direct W [N]'\n"
                                          "  occupation: 'occupation N n1 n2 n3 ...'\n"
                                          "  sackur_tetrode: 'sackur T M_kg [V_m3] [P_atm]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with entropy_S_J/K, entropy_S_eV/K, ln_W, and mode-specific details."),
    ]

    examples         = [
        {
            "code_input": {
                "mode":                  "direct",
                "microstate_count_W":    1024,
                "number_of_particles_N": 10,
                "occupation_numbers":     None,
                "temperature_k":         None,
                "molecular_mass_kg":      None,
                "volume_m3":              None,
                "pressure_atm":           None,
            },
            "text_input": {
                "input_params":           "direct 1024 10",
            },
            "output": {
                "result": {
                    "mode": "direct",
                    "W": 1024,
                    "ln_W": 6.9315,
                    "entropy_J_per_K": 9.57e-23,
                    "entropy_eV_per_K": 5.97e-4,
                }
            },
        },
        {
            "code_input": {
                "mode":               "sackur_tetrode",
                "microstate_count_W": None,
                "number_of_particles_N":None,
                "occupation_numbers":   None,
                "temperature_k":      298.15,
                "molecular_mass_kg":   4.0026e-3 / 6.022e23,  # He atom
                "volume_m3":           None,
                "pressure_atm":        1.0,
            },
            "text_input": {
                "input_params":        "sackur_tetrode 298.15 6.65e-27",
            },
            "output": {
                "result": {
                    "mode": "sackur_tetrode",
                    "entropy_J_per_K": 126.2,
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
        self.k_B = 1.380649e-23     # J/K
        self.h   = 6.62607015e-34    # J·s
        self.NA  = 6.02214076e23     # mol⁻¹
        self.eV  = 1.602176634e-19   # J/eV
        self.pi  = math.pi
        self.atm_to_Pa = 101325.0

    def _mode_direct(self, W: int, N: Optional[int] = None) -> dict:
        """Mode 1: Direct S = k_B · ln(W)"""
        if W <= 0:
            raise ChemMCPError("W must be a positive integer.")
        
        ln_W = math.log(float(W))
        S_J = self.k_B * ln_W
        S_eV = S_J / self.eV
        
        result = {
            "mode":              "direct",
            "W":                 int(W),
            "ln_W":              round(ln_W, 12),
            "entropy_J_per_K":   round(S_J, 25),
            "entropy_eV_per_K":  round(S_eV, 16),
            "entropy_J_mol_K":   round(S_J * self.NA, 4),
            "N":                 N,
        }
        
        if N is not None and N > 0:
            # Per-particle entropy
            result["entropy_per_particle_J_per_K"] = round(S_J / N, 28)
            # Information entropy per particle in bits: s = (ln W/N) / ln(2)
            result["information_entropy_bits_per_particle"] = round(ln_W / (N * math.log(2)), 10)
        
        return result

    def _mode_occupation(self, N: int, occupation_numbers: List[int]) -> dict:
        """Mode 2: S = k_B · ln(W) where W = N!/(n₁!·n₂!·...·nₖ!)
        
        Using Stirling's approximation: ln(n!) ≈ n·ln(n) - n
        So: ln(W) ≈ N·ln(N) - N - Σ[nᵢ·ln(nᵢ) - nᵢ]
                   = N·ln(N) - Σ nᵢ·ln(nᵢ)
        """
        if N <= 0:
            raise ChemMCPError("Number of particles N must be > 0.")
        
        if not occupation_numbers:
            raise ChemMCPError("Occupation numbers list cannot be empty.")
        
        ns = [int(n) for n in occupation_numbers]
        
        if sum(ns) != N:
            raise ChemMCPError(f"Sum of occupation numbers ({sum(ns)}) must equal N ({N}).")
        
        # Use Stirling approximation for ln(W)
        # Handle n_i = 0: ln(0!) = 0 by convention
        sum_ni_ln_ni = 0.0
        for ni in ns:
            if ni < 0:
                raise ChemMCPError("Occupation numbers cannot be negative.")
            if ni > 0:
                sum_ni_ln_ni += ni * math.log(float(ni))
            # ni == 0 contributes 0
        
        ln_W = N * math.log(float(N)) - sum_ni_ln_ni
        S_J = self.k_B * ln_W
        S_eV = S_J / self.eV
        
        return {
            "mode":               "occupation",
            "N":                  N,
            "num_states":         len(ns),
            "occupation_numbers": ns,
            "ln_W_stirling":      round(ln_W, 12),
            "W_approx":           round(math.exp(min(ln_W, 700)), 4) if ln_W < 700 else float('inf'),
            "entropy_J_per_K":    round(S_J, 25),
            "entropy_eV_per_K":   round(S_eV, 16),
            "entropy_J_mol_K":    round(S_J * self.NA, 4),
            "formula_used":       "ln(W) ≈ N·ln(N) - Σnᵢ·ln(nᵢ)  (Stirling approximation)",
        }

    def _mode_sackur_tetrode(
        self,
        temperature_k: float,
        molecular_mass_kg: float,
        volume_m3: Optional[float] = None,
        pressure_atm: float = 1.0,
    ) -> dict:
        """Mode 3: Sackur-Tetrode equation for monatomic ideal gas.
        
        S = Nk · [ln(V/N · (4πmE/(3Nh²))^(3/2)) + 5/2]
          = Nk · [ln(V/N · (2πmkT/h²)^(3/2)) + 5/2]
        
        For 1 mole: S = R · [ln((V/N_A)·(2πmkT/h²)^(3/2)) + 5/2]
        Or per molecule: S/k = ln[(V/N)·(2πmkT/h²)^(3/2)] + 5/2
        """
        T = float(temperature_k)
        m = float(molecular_mass_kg)
        
        if T <= 0:
            raise ChemMCPError("Temperature must be > 0 K.")
        if m <= 0:
            raise ChemMCPError("Mass must be > 0.")

        # Default volume from ideal gas law at given pressure
        if volume_m3 is None:
            V = (self.NA * self.k_B * T) / (pressure_atm * self.atm_to_Pa)
        else:
            V = float(volume_m3)

        # Thermal de Broglie wavelength cubed: Λ³ = (h²/(2πmkT))^(3/2)
        Lambda_cubed = (self.h ** 2 / (2.0 * self.pi * m * self.k_B * T)) ** 1.5
        
        # S/k per molecule
        S_over_k = math.log(V / (self.NA * Lambda_cubed) if Lambda_cubed > 0 else 1.0) + 2.5
        
        # For N_A molecules (1 mole)
        S_mol_J = self.NA * self.k_B * S_over_k
        # Per molecule
        S_molecule_J = self.k_B * S_over_k
        
        return {
            "mode":                  "sackur_tetrode",
            "temperature_K":         T,
            "molecular_mass_kg":     m,
            "volume_m3":             round(V, 6),
            "pressure_atm":          pressure_atm,
            "thermal_deBroglie_lambda_m": round(Lambda_cubed ** (1.0/3.0), 20),
            "S_over_k_per_molecule": round(S_over_k, 12),
            "entropy_J_mol_K":       round(S_mol_J, 4),
            "entropy_J_per_K":       round(S_mol_J, 4),
            "entropy_eV_mol_K":      round(S_mol_J / (self.NA * self.eV), 8),
            "formula_used":          "S = Nk·[ln(V/N·(2πmkT/h²)^(3/2)) + 5/2]",
        }

    def _run_base(
        self,
        mode: str = "direct",
        microstate_count_W: Optional[int] = None,
        number_of_particles_N: Optional[int] = None,
        occupation_numbers: Optional[List[int]] = None,
        temperature_k: Optional[float] = None,
        molecular_mass_kg: Optional[float] = None,
        volume_m3: Optional[float] = None,
        pressure_atm: float = 1.0,
    ) -> dict:
        """Core logic: dispatch to appropriate mode."""
        mode = mode.lower().strip()
        
        if mode == "direct":
            if microstate_count_W is None:
                raise ChemMCPError("Mode 'direct' requires microstate_count_W.")
            return self._mode_direct(microstate_count_W, number_of_particles_N)
        
        elif mode == "occupation":
            if number_of_particles_N is None or occupation_numbers is None:
                raise ChemMCPError("Mode 'occupation' requires both N and occupation_numbers.")
            return self._mode_occupation(number_of_particles_N, occupation_numbers)
        
        elif mode in ("sackur_tetrode", "sackur"):
            if temperature_k is None or molecular_mass_kg is None:
                raise ChemMCPError("Mode 'sackur_tetrode' requires both temperature_k and molecular_mass_kg.")
            return self._mode_sackur_tetrode(temperature_k, molecular_mass_kg, volume_m3, pressure_atm)
        
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Choose from: 'direct', 'occupation', 'sackur_tetrode'.")

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 'mode ...' params.")
            
            mode = parts[0].lower()
            
            if mode == "direct":
                W = int(parts[1])
                N = int(parts[2]) if len(parts) > 2 else None
                return self._run_base("direct", W, N)
            
            elif mode == "occupation":
                N = int(parts[1])
                occ_nums = [int(p) for p in parts[2:]]
                return self._run_base("occupation", None, N, occ_nums)
            
            elif mode in ("sackur", "sackur_tetrode"):
                T = float(parts[1])
                m = float(parts[2])
                V = float(parts[3]) if len(parts) > 3 and parts[3].lower() != "none" else None
                P = float(parts[4]) if len(parts) > 4 else 1.0
                return self._run_base("sackur_tetrode", None, None, None, T, m, V, P)
            
            else:
                raise ValueError(f"Unknown mode: {mode}")
                
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}.\n"
                f"Formats:\n"
                f"  'direct W [N]'\n"
                f"  'occupation N n1 n2 ...'\n"
                f"  'sackur_tetrode T(K) M(kg) [V(m3)] [P(atm)]'"
            )
