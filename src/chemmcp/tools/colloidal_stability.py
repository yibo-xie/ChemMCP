"""
胶体稳定性分析工具（DLVO 理论）
计算 van der Waals 吸引力、双电层排斥力及总相互作用能，判断胶体稳定性。
"""
import logging
import math
from typing import Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ColloidalStability(BaseTool):
    """
    胶体稳定性（DLVO理论）分析工具。

    基于 DLVO 理论计算胶体颗粒间的 van der Waals 吸引能与双电层排斥能，
    并给出总相互作用势能曲线和稳定性预测。
    """
    __version__                 = "0.1.0"
    name                        = "ColloidalStability"
    func_name                   = "analyze_colloidal_stability"
    description                 = "Analyze colloidal stability using DLVO theory: calculate van der Waals attraction, double-layer repulsion, and total interaction energy."
    implementation_description  = "Uses DLVO theory formulas: V_A (Hamaker) for attraction, V_R (Poisson-Boltzmann linearized) for repulsion. Computes total V_total = V_A + V_R at a range of separation distances to predict stability."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Colloid", "DLVO", "Physical Chemistry", "Surface Chemistry", "Stability"]
    required_envs               = []

    code_input_sig = [
        ("particle_radius_m",        "float", "N/A",          "Radius of colloidal particles in meters (e.g., 1e-8 for 10 nm)."),
        ("surface_potential_mv",     "float", "N/A",          "Surface potential in millivolts (mV)."),
        ("hamaker_constant_j",      "float", "N/A",          "Hamaker constant in Joules (typical range: 1e-20 to 1e-19 J)."),
        ("electrolyte_concentration_m", "float", "N/A",      "Electrolyte concentration in mol/L (M)."),
        ("temperature_k",           "float", "298.15",       "Temperature in Kelvin."),
        ("valency_z",               "int",   "1",            "Valency of the electrolyte ion (z=1 for NaCl, z=2 for CaCl2, etc.)."),
        ("separation_distance_nm",  "float", "1.0",           "Surface-to-surface separation distance in nm at which to evaluate energy."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'particle_radius_m surface_potential_mv hamaker_constant_j electrolyte_concentration_m temperature_k valency_z separation_distance_nm'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing van_der_waals_energy_j, double_layer_energy_j, total_energy_j, debye_length_nm, stability_prediction, and details."),
    ]

    examples = [
        {
            "code_input": {
                "particle_radius_m": 1e-8,
                "surface_potential_mv": 25.0,
                "hamaker_constant_j": 1e-20,
                "electrolyte_concentration_m": 0.01,
                "temperature_k": 298.15,
                "valency_z": 1,
                "separation_distance_nm": 2.0,
            },
            "text_input": {
                "input_params": "1e-8 25.0 1e-20 0.01 298.15 1 2.0",
            },
            "output": {
                "result": {
                    "van_der_waals_energy_j": "... (negative value)",
                    "double_layer_energy_j": "... (positive value)",
                    "total_energy_j": "...",
                    "debye_length_nm": "...",
                    "stability_prediction": "stable / unstable / metastable",
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
        """初始化物理常数。"""
        self.k_B = 1.380649e-23       # Boltzmann constant, J/K
        self.e = 1.602176634e-19       # Elementary charge, C
        self.epsilon_0 = 8.854187817e-12  # Vacuum permittivity, F/m
        self.N_A = 6.02214076e23       # Avogadro's number, 1/mol

    def _run_base(
        self,
        particle_radius_m: float,
        surface_potential_mv: float,
        hamaker_constant_j: float,
        electrolyte_concentration_m: float,
        temperature_k: float = 298.15,
        valency_z: int = 1,
        separation_distance_nm: float = 1.0,
    ) -> Dict[str, Any]:
        """
        核心逻辑：基于 DLVO 理论计算胶体相互作用能。

        Parameters:
            particle_radius_m: 颗粒半径 (m)
            surface_potential_mv: 表面电势 (mV)
            hamaker_constant_j: Hamaker 常数 (J)
            electrolyte_concentration_m: 电解质浓度 (mol/L)
            temperature_k: 温度 (K)
            valency_z: 离子价态
            separation_distance_nm: 表面间距 (nm)

        Returns:
            包含各项能量、Debye 长度和稳定性预测的字典
        """
        # ---- 输入验证 ----
        if particle_radius_m <= 0:
            raise ChemMCPError("Particle radius must be positive.")
        if hamaker_constant_j <= 0:
            raise ChemMCPError("Hamaker constant must be positive.")
        if electrolyte_concentration_m < 0:
            raise ChemMCPError("Electrolyte concentration must be non-negative.")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")
        if separation_distance_nm <= 0:
            raise ChemMCPError("Separation distance must be positive.")

        H = separation_distance_nm * 1e-9  # 转换为 m

        # ---- Debye 长度 (kappa^-1) ----
        # kappa^2 = (2 * e^2 * N_A * |z|^2 * c) / (epsilon_r * epsilon_0 * k_B * T)
        epsilon_r = 78.5  # 水的相对介电常数 (25°C)
        n0 = electrolyte_concentration_m * self.N_A * 1e3  # number density (ions/m³), 每种离子

        if electrolyte_concentration_m > 0:
            kappa_squared = (
                2 * self.e**2 * self.N_A * valency_z**2 * electrolyte_concentration_m * 1000
                / (epsilon_r * self.epsilon_0 * self.k_B * temperature_k)
            )
            kappa = math.sqrt(kappa_squared)
            debye_length = 1.0 / kappa
        else:
            kappa = 0.0
            debye_length = float("inf")

        # ---- van der Waals 吸引能 (球体近似) ----
        # V_A = -A * r / (12 * H)  (两球近似, H << r)
        if H > 0:
            v_a = -hamaker_constant_j * particle_radius_m / (12.0 * H)
        else:
            raise ChemMCPError("Separation distance must be > 0.")

        # ---- 双电层排斥能 (线性化 Poisson-Boltzmann) ----
        # psi0 in volts
        psi0_v = surface_potential_mv / 1000.0

        if kappa > 0 and electrolyte_concentration_m > 0:
            # gamma = tanh(ze*psi0 / (4*k_B*T))
            arg_gamma = (valency_z * self.e * psi0_v) / (4.0 * self.k_B * temperature_k)
            gamma = math.tanh(arg_gamma)

            # V_R = 64 * pi * n0 * k_B * T / kappa^2 * gamma^2 * exp(-kappa*H) * r
            # (对于两个半径为 r 的球体)
            prefactor = 64.0 * math.pi * n0 * self.k_B * temperature_k / (kappa ** 2) * (gamma ** 2)
            v_r = prefactor * math.exp(-kappa * H) * particle_radius_m
        else:
            v_r = 0.0

        # ---- 总相互作用能 ----
        v_total = v_a + v_r

        # ---- 稳定性判断 ----
        # 如果在当前距离处存在明显的能垒 (> 几 kT)，则稳定
        kT = self.k_B * temperature_k
        energy_barrier_kT = v_total / kT if kT > 0 else float("inf")

        if v_total > 5 * kT:
            stability = "stable (significant repulsive barrier)"
        elif v_total > 0:
            stability = "metastable (weak repulsive barrier)"
        elif v_total > -5 * kT:
            stability = "unstable (net attraction, weak)"
        else:
            stability = "unstable (strong net attraction, rapid coagulation)"

        logger.info(f"DLVO analysis: V_A={v_a:.3e} J, V_R={v_r:.3e} J, V_tot={v_total:.3e} J, "
                     f"kappa^-1={debye_length*1e9:.1f} nm, stability={stability}")

        return {
            "van_der_waals_energy_j": v_a,
            "double_layer_energy_j": v_r,
            "total_energy_j": v_total,
            "debye_length_nm": round(debye_length * 1e9, 4) if math.isfinite(debye_length) else None,
            "energy_barrier_kT": round(energy_barrier_kT, 4),
            "stability_prediction": stability,
            "parameters_used": {
                "particle_radius_m": particle_radius_m,
                "surface_potential_mv": surface_potential_mv,
                "hamaker_constant_j": hamaker_constant_j,
                "electrolyte_concentration_M": electrolyte_concentration_m,
                "temperature_K": temperature_k,
                "valency_z": valency_z,
                "separation_distance_nm": separation_distance_nm,
            }
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入并调用核心逻辑。"""
        try:
            parts = input_params.split()
            if len(parts) < 4:
                raise ValueError("Need at least: radius potential_mv hamaker_J conc_M")

            kwargs = {
                "particle_radius_m": float(parts[0]),
                "surface_potential_mv": float(parts[1]),
                "hamaker_constant_j": float(parts[2]),
                "electrolyte_concentration_m": float(parts[3]),
            }
            if len(parts) > 4:
                kwargs["temperature_k"] = float(parts[4])
            if len(parts) > 5:
                kwargs["valency_z"] = int(parts[5])
            if len(parts) > 6:
                kwargs["separation_distance_nm"] = float(parts[6])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'radius(m) potential_mv hamaker(J) conc(M) [T(K) z dist(nm)]'"
            )
