import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TemperatureJumpRelaxation(BaseTool):
    """
    温度跳跃弛豫动力学分析（T-Jump Method）。
    突然改变温度使体系偏离平衡，然后监测恢复平衡的弛豫过程。
    
    弛豫时间 τ = 1/(kf + kr) 对于简单反应 A<=>B
    """
    __version__ = "0.1.0"
    name = "TemperatureJumpRelaxation"
    func_name = "temp_jump_relaxation"
    description = "Analyze temperature-jump (T-jump) relaxation kinetics for studying fast reaction dynamics."
    implementation_description = "Calculates relaxation time constant (τ), new equilibrium position after T-jump, and describes approach to equilibrium. Supports simple reversible reactions A<=>B and provides kinetic analysis of the relaxation process."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Relaxation", "Temperature Jump", "Fast Reactions"]
    required_envs = []

    code_input_sig = [
        ("equilibrium_constant_K", "float", "N/A", "Equilibrium constant at the NEW (final) temperature."),
        ("delta_H", "float", "N/A", "Enthalpy change of reaction in J/mol."),
        ("initial_temperature_K", "float", "N/A", "Initial temperature before jump in Kelvin."),
        ("final_temperature_K", "float", "N/A", "Final temperature after T-jump in Kelvin."),
        ("forward_rate_constant_kf", "float", "N/A", "Forward rate constant at final T (optional if τ given)."),
        ("reverse_rate_constant_kr", "float", "N/A", "Reverse rate constant at final T (optional if τ given)."),
        ("relaxation_time_tau", "float", "N/A", "Directly provide relaxation time τ in seconds (optional; calculated from kf+kr if not given)."),
        ("initial_A_concentration", "float", "0.0", "Initial concentration of A before perturbation (in M or same units as total)."),
        ("total_concentration", "float", "1.0", "Total concentration [A]+[B] (default 1.0 M)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'K delta_H T_initial T_final kf kr [tau] [A0] [total]'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with relaxation time tau, new equilibrium concentrations, displacement decay profile, and kinetic analysis."),
    ]

    examples = [
        {
            "code_input": {
                "equilibrium_constant_K": 2.0,
                "delta_H": -50000.0,
                "initial_temperature_K": 298.0,
                "final_temperature_K": 308.0,
                "forward_rate_constant_kf": 1000.0,
                "reverse_rate_constant_kr": 500.0,
                "relaxation_time_tau": None,
                "initial_A_concentration": 0.0,
                "total_concentration": 1.0,
            },
            "text_input": {
                "input_params": "2.0 -50000 298 308 1000 500",
            },
            "output": {
                "result": {
                    "relaxation_time_tau_s": 0.000667,
                    "new_eq_A": 0.3333,
                    "new_eq_B": 0.6667,
                    "reaction_type": "A<=>B",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314  # J/(mol·K)

    def _run_base(self, equilibrium_constant_K: float, delta_H: float,
                  initial_temperature_K: float, final_temperature_K: float,
                  forward_rate_constant_kf: float = None,
                  reverse_rate_constant_kr: float = None,
                  relaxation_time_tau: float = None,
                  initial_A_concentration: float = 0.0,
                  total_concentration: float = 1.0) -> dict:
        
        if final_temperature_K <= 0 or initial_temperature_K <= 0:
            raise ChemMCPError("Temperatures must be positive in Kelvin.")
        if equilibrium_constant_K <= 0:
            raise ChemMCPError("Equilibrium constant must be positive.")

        # Calculate relaxation time from rate constants if not provided
        if relaxation_time_tau is None:
            if forward_rate_constant_kf is None or reverse_rate_constant_kr is None:
                raise ChemMCPError(
                    "Must provide either relaxation_time_tau OR both forward_rate_constant_kf "
                    "and reverse_rate_constant_kr.")
            # For A<=>B: τ = 1 / (kf + kr)
            tau = 1.0 / (forward_rate_constant_kf + reverse_rate_constant_kr)
        else:
            tau = relaxation_time_tau

        # New equilibrium position at final T
        K_new = equilibrium_constant_K  # User provides K at final T
        A_eq = total_concentration / (1 + K_new)
        B_eq = K_new * total_concentration / (1 + K_new)

        # Displacement from equilibrium (Δ[A])
        # If initial_A_concentration is 0, assume system was at old equilibrium
        # For demonstration, compute displacement
        delta_T = final_temperature_K - initial_temperature_K
        
        # van't Hoff to estimate old K (approximate)
        if abs(delta_T) > 0.1 and abs(delta_H) > 0:
            # ln(K2/K1) = -ΔH/R * (1/T2 - 1/T1)
            ln_ratio = -delta_H / self.R * (1.0/final_temperature_K - 1.0/initial_temperature_K)
            K_old = K_new / math.exp(ln_ratio)
            A_old_eq = total_concentration / (1 + K_old) if K_old > 0 else total_concentration
        else:
            A_old_eq = A_eq

        # Initial displacement
        if initial_A_concentration > 0:
            A_init = initial_A_concentration
        else:
            A_init = A_old_eq

        displacement = A_init - A_eq

        # Characteristic times: at t = τ, exp(-1) ≈ 37% of displacement remains
        # At t = 5τ, < 1% remains
        result = {
            "relaxation_time_tau_s": round(tau, 10),
            "new_equilibrium_A_M": round(A_eq, 8),
            "new_equilibrium_B_M": round(B_eq, 8),
            "displacement_from_equilibrium_M": round(displacement, 8),
            "temperature_jump_K": round(delta_T, 2),
            "delta_H_J_mol": delta_H,
            "equilibrium_constant_K_final": K_new,
            "fraction_remaining_at_tau": round(math.exp(-1), 4),  # ~0.368
            "fraction_remaining_at_5tau": round(math.exp(-5), 4),  # ~0.007
            "time_to_99_percent_equilibrium_s": round(5 * tau, 10),
            "kinetic_analysis": (
                f"After T-jump of {delta_T:+.1f}K, the system relaxes with τ={tau:.4e}s. "
                f"New equilibrium: [A]={A_eq:.4f}, [B]={B_eq:.4f}. "
                f"The displacement decays exponentially: Δ(t) = Δ₀·exp(-t/τ)."
            ),
        }

        if forward_rate_constant_kf is not None:
            result["forward_rate_constant_kf"] = forward_rate_constant_kf
        if reverse_rate_constant_kr is not None:
            result["reverse_rate_constant_kr"] = reverse_rate_constant_kr

        logger.info(f"TemperatureJumpRelaxation: τ={tau:.4e}s, ΔT={delta_T}K")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            K = float(parts[0])
            dH = float(parts[1])
            T_init = float(parts[2])
            T_final = float(parts[3])
            kf = float(parts[4]) if len(parts) > 4 else None
            kr = float(parts[5]) if len(parts) > 5 else None
            tau = float(parts[6]) if len(parts) > 6 else None
            A0 = float(parts[7]) if len(parts) > 7 else 0.0
            total = float(parts[8]) if len(parts) > 8 else 1.0
            return self._run_base(K, dH, T_init, T_final, kf, kr, tau, A0, total)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'K dH Ti Tf kf kr [tau] [A0] [total]'")
