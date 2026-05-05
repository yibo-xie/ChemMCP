"""
Zeta 电位与双电层计算工具
根据电泳迁移率计算 Zeta 电位，并计算 Debye 长度（双电层厚度）。
"""
import logging
import math
from typing import Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ZetaPotential(BaseTool):
    """
    Zeta 电位与双电层计算工具。

    根据电泳迁移率（electrophoretic mobility）计算 Zeta 电位，
    同时计算 Debye 长度和双电层厚度。
    支持 Smoluchowski 极限（大颗粒，κa >> 1）和 Hückel 极限（小颗粒，κa << 1）。
    """
    __version__                 = "0.1.0"
    name                        = "ZetaPotential"
    func_name                   = "calculate_zeta_potential"
    description                 = "Calculate zeta potential from electrophoretic mobility, and compute Debye length (double-layer thickness)."
    implementation_description  = "Uses Smoluchowski (κa>>1) or Hückel (κa<<1) approximation to convert electrophoretic mobility to zeta potential. Also computes the Debye length from electrolyte concentration."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Zeta Potential", "Double Layer", "Colloid", "Electrochemistry", "Surface Chemistry"]
    required_envs               = []

    code_input_sig = [
        ("electrophoretic_mobility", "float", "N/A",     "Electrophoretic mobility in m²/(V·s). Typical range: 1e-8 to 5e-8."),
        ("temperature_k",           "float", "298.15",   "Temperature in Kelvin."),
        ("viscosity_pa_s",          "float", "0.000894", "Dynamic viscosity of solvent in Pa·s (water at 25°C ≈ 0.000894)."),
        ("dielectric_constant",     "float", "78.5",      "Relative dielectric constant of solvent (water ≈ 78.5 at 25°C)."),
        ("electrolyte_concentration_m", "float", "0.01",  "Electrolyte concentration in mol/L (M), for Debye length calculation."),
        ("valency_z",               "int",   "1",        "Valency of symmetrical electrolyte (e.g., z=1 for NaCl, z=2 for MgSO4)."),
        ("particle_radius_m",       "float", "1e-7",      "Particle radius in meters, used to determine Smoluchowski vs Hückel regime."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A",
         "Space-separated string: 'mobility(m²/Vs) [T(K)] [viscosity(Pa·s)] [epsilon_r] [conc(M)] [z] [radius(m)]'"),
    ]

    output_sig = [
        ("result", "dict",
         "Dictionary containing zeta_potential_mv, debye_length_nm, double_layer_thickness_nm, regime_used, and details."),
    ]

    examples = [
        {
            "code_input": {
                "electrophoretic_mobility": 3.5e-8,
                "temperature_k": 298.15,
                "viscosity_pa_s": 0.000894,
                "dielectric_constant": 78.5,
                "electrolyte_concentration_m": 0.01,
                "valency_z": 1,
                "particle_radius_m": 1e-7,
            },
            "text_input": {
                "input_params": "3.5e-8 298.15 0.000894 78.5 0.01 1 1e-7",
            },
            "output": {
                "result": {
                    "zeta_potential_mv": "... (mV value)",
                    "debye_length_nm": "... (nm)",
                    "regime_used": "Smoluchowski or Hückel or intermediate",
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
        electrophoretic_mobility: float,
        temperature_k: float = 298.15,
        viscosity_pa_s: float = 0.000894,
        dielectric_constant: float = 78.5,
        electrolyte_concentration_m: float = 0.01,
        valency_z: int = 1,
        particle_radius_m: float = 1e-7,
    ) -> Dict[str, Any]:
        """核心逻辑：计算 Zeta 电位和双电层参数。"""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")
        if viscosity_pa_s <= 0:
            raise ChemMCPError("Viscosity must be positive.")
        if dielectric_constant <= 0:
            raise ChemMCPError("Dielectric constant must be positive.")
        if particle_radius_m <= 0:
            raise ChemMCPError("Particle radius must be positive.")

        # ---- 计算 Debye 长度 ----
        epsilon = dielectric_constant * self.epsilon_0

        if electrolyte_concentration_m > 0:
            kappa_squared = (
                2 * self.e**2 * self.N_A * valency_z**2 * electrolyte_concentration_m * 1000
                / (epsilon * self.k_B * temperature_k)
            )
            kappa = math.sqrt(kappa_squared)
            debye_length = 1.0 / kappa
        else:
            kappa = 0.0
            debye_length = float("inf")

        # ---- 判断 Smoluchowski vs Hückel vs 中间区域 ----
        kappa_a = kappa * particle_radius_m if kappa > 0 else 0.0
        epsilon_eff = dielectric_constant * self.epsilon_0

        if kappa_a > 100:
            regime = "Smoluchowski (κa >> 1)"
            zeta_v = viscosity_pa_s * electrophoretic_mobility / epsilon_eff
        elif kappa_a < 1:
            regime = "Hückel (κa << 1)"
            zeta_v = 3.0 * viscosity_pa_s * electrophoretic_mobility / (2.0 * epsilon_eff)
        else:
            # 中间区域：基于 log10(κa) 的线性插值近似 Henry 函数
            log_ka = math.log10(kappa_a) if kappa_a > 0 else -10
            t = max(0.0, min(1.0, (log_ka + 1.0) / 3.0))  # map log_ka ∈ [-1,2] → t ∈ [0,1]
            f_henry = 1.5 - 0.5 * t  # f: 1.5 → 1.0
            f_henry = max(1.0, min(1.5, f_henry))
            regime = f"intermediate Henry (κa={kappa_a:.1f}, f≈{f_henry:.2f})"
            zeta_v = viscosity_pa_s * electrophoretic_mobility / (f_henry * epsilon_eff)

        zeta_mv = zeta_v * 1000.0

        logger.info(f"Zeta potential: {zeta_mv:.2f} mV, Debye length: {debye_length*1e9:.2f} nm, "
                     f"regime: {regime}")

        return {
            "zeta_potential_mv": round(zeta_mv, 4),
            "zeta_potential_v": round(zeta_v, 10),
            "debye_length_nm": round(debye_length * 1e9, 4) if math.isfinite(debye_length) else None,
            "double_layer_thickness_nm": round(debye_length * 1e9, 4) if math.isfinite(debye_length) else None,
            "kappa_a": round(kappa_a, 4),
            "regime_used": regime,
            "parameters_used": {
                "electrophoretic_mobility_m2_Vs": electrophoretic_mobility,
                "temperature_K": temperature_k,
                "viscosity_Pa_s": viscosity_pa_s,
                "dielectric_constant": dielectric_constant,
                "electrolyte_concentration_M": electrolyte_concentration_m,
                "valency_z": valency_z,
                "particle_radius_m": particle_radius_m,
            }
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """解析文本输入并调用核心逻辑。"""
        try:
            parts = input_params.split()
            if len(parts) < 1:
                raise ValueError("Need at least: electrophoretic_mobility")

            kwargs = {"electrophoretic_mobility": float(parts[0])}
            if len(parts) > 1: kwargs["temperature_k"] = float(parts[1])
            if len(parts) > 2: kwargs["viscosity_pa_s"] = float(parts[2])
            if len(parts) > 3: kwargs["dielectric_constant"] = float(parts[3])
            if len(parts) > 4: kwargs["electrolyte_concentration_m"] = float(parts[4])
            if len(parts) > 5: kwargs["valency_z"] = int(parts[5])
            if len(parts) > 6: kwargs["particle_radius_m"] = float(parts[6])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {str(e)}. "
                f"Format: 'mobility(m²/Vs) [T(K)] [viscosity(Pa·s)] [epsilon_r] [conc(M)] [z] [radius(m)]'"
            )
