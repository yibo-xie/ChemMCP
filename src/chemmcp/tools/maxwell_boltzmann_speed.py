import logging
import math
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MaxwellBoltzmannSpeed(BaseTool):
    """
    Maxwell-Boltzmann 速率分布计算工具。
    
    计算最概然速率、平均速率、方均根速率，以及速率分布概率密度函数 (PDF) 数据。
    
    PDF: f(v) = 4π · (m/2πkT)^(3/2) · v² · exp(-mv²/2kT)
    """
    __version__                 = "0.1.0"
    name                        = "MaxwellBoltzmannSpeed"
    func_name                   = "calculate_mb_speed_distribution"
    description                 = "Calculate Maxwell-Boltzmann speed distribution: most probable, average, RMS speeds and full PDF curve data."
    implementation_description  = "Uses f(v) = 4π(m/2πkT)^(3/2)·v²·exp(-mv²/2kT). Computes characteristic speeds and generates distribution curve data for plotting."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Maxwell-Boltzmann", "Kinetic Theory", "Speed Distribution", "Statistical Mechanics"]
    required_envs               = []

    code_input_sig   = [
        ("temperature_k",            "float",  "N/A",     "Temperature in Kelvin."),
        ("molar_mass_kg_per_mol",    "float",  "N/A",     "Molar mass in kg/mol."),
        ("speed_m_s",                "float",  "None",    "Specific speed in m/s to evaluate PDF at (optional)."),
        ("v_max_factor",             "float",  "3.0",     "Multiple of v_mp for curve range max (v_max = factor * v_mp)."),
        ("num_points",               "int",    "200",     "Number of points in distribution curve."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'T(K) M_kg/mol [speed_m/s] [v_max_factor] [num_points]'"),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with v_mp, v_avg, v_rms, pdf_at_speed, distribution_curve."),
    ]

    examples         = [
        {
            "code_input": {
                "temperature_k":           300.0,
                "molar_mass_kg_per_mol":   0.032,   # O2
                "speed_m_s":               None,
                "v_max_factor":            3.0,
                "num_points":              200,
            },
            "text_input": {
                "input_params":            "300.0 0.032",
            },
            "output": {
                "result": {
                    "temperature_K": 300.0,
                    "molar_mass_kg_mol": 0.032,
                    "v_mp_m_s": 394.89,
                    "v_avg_m_s": 444.69,
                    "v_rms_m_s": 483.56,
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
        self.k_B = 1.380649e-23   # J/K
        self.NA = 6.02214076e23   # mol⁻¹
        self.pi = math.pi

    def _pdf(self, v: float, a: float) -> float:
        """Maxwell-Boltzmann speed PDF: f(v) = 4π·(a²/√π)·v³? No:
        
        f(v) = 4π · (m/(2πkT))^(3/2) · v² · exp(-mv²/(2kT))
        
        Let a = √(m/(2kT)), then f(v) = 4π · (a³/√π) · v² · exp(-a²v²)? No.
        
        Actually: f(v) = 4π · (m/(2πkT))^(3/2) · v² · exp(-mv²/(2kT))
        Let α = m/(2kT), then f(v) = 4π · (α/π)^(3/2) · v² · exp(-αv²)
        """
        if v < 0:
            return 0.0
        # a = sqrt(m / (2 * k_B * T)) is embedded in the prefactor
        # Use direct formula with precomputed alpha_over_pi
        return 4.0 * self.pi * (a ** 1.5) / math.sqrt(self.pi) * (v ** 2) * math.exp(-(a * v * v))

    def _run_base(
        self,
        temperature_k: float,
        molar_mass_kg_per_mol: float,
        speed_m_s: Optional[float] = None,
        v_max_factor: float = 3.0,
        num_points: int = 200,
    ) -> dict:
        """Core logic."""
        global a
        T = float(temperature_k)
        M = float(molar_mass_kg_per_mol)
        
        if T <= 0:
            raise ChemMCPError("Temperature must be > 0 K.")
        if M <= 0:
            raise ChemMCPError("Molar mass must be > 0.")

        # Molecular mass (per molecule)
        m = M / self.NA
        
        # Key parameter: α = m / (2kT)
        alpha = m / (2.0 * self.k_B * T)
        a = alpha  # for use in _pdf... but we need it as global or pass through
        
        # --- Characteristic speeds ---
        # Most probable: v_mp = √(2kT/m) = √(2RT/M)
        v_mp = math.sqrt(2.0 * self.k_B * T / m)
        
        # Average: v_avg = √(8kT/πm) = √(8RT/πM)
        v_avg = math.sqrt(8.0 * self.k_B * T / (math.pi * m))
        
        # Root-mean-square: v_rms = √(3kT/m) = √(3RT/M)
        v_rms = math.sqrt(3.0 * self.k_B * T / m)

        # --- PDF at specific speed ---
        pdf_val = None
        if speed_m_s is not None:
            v = float(speed_m_s)
            prefactor = 4.0 * math.pi * (alpha / math.pi) ** 1.5
            pdf_val = prefactor * (v ** 2) * math.exp(-alpha * v * v)

        # --- Distribution curve data ---
        v_max = v_max_factor * v_mp
        dv = v_max / num_points if num_points > 0 else v_max
        prefactor = 4.0 * math.pi * (alpha / math.pi) ** 1.5
        
        curve_data = []
        for i in range(num_points + 1):
            v = i * dv
            fv = prefactor * (v ** 2) * math.exp(-alpha * v * v)
            curve_data.append({
                "speed_m_s": round(v, 4),
                "pdf_s_per_m": round(fv, 10),
            })

        result = {
            "temperature_K":          T,
            "molar_mass_kg_per_mol":  M,
            "molecular_mass_kg":      round(m, 25),
            "v_mp_m_s":              round(v_mp, 4),
            "v_avg_m_s":             round(v_avg, 4),
            "v_rms_m_s":             round(v_rms, 4),
            "ratio_vavg_vmp":         round(v_avg / v_mp, 6),
            "ratio_vrms_vmp":         round(v_rms / v_mp, 6),
            "pdf_at_speed":           {"speed_m_s": speed_m_s, "pdf_value": round(pdf_val, 12)} if speed_m_s is not None else None,
            "distribution_curve":     {
                "v_min_m_s": 0.0,
                "v_max_m_s": round(v_max, 4),
                "num_points": len(curve_data),
                "peak_pdf_value": round(prefactor * (v_mp ** 2) * math.exp(-alpha * v_mp * v_mp), 12),
                "data": curve_data[:30],  # First 30 points
                "total_points_computed": len(curve_data),
            },
        }

        logger.info(f"MB Speed: T={T}K, M={M}kg/mol => v_mp={v_mp:.2f}, v_avg={v_avg:.2f}, v_rms={v_rms:.2f} m/s")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 'T M' params.")
            
            T = float(parts[0])
            M = float(parts[1])
            v = float(parts[2]) if len(parts) > 2 and parts[2].lower() != "none" else None
            vf = float(parts[3]) if len(parts) > 3 else 3.0
            np_ = int(parts[4]) if len(parts) > 4 else 200
            
            return self._run_base(T, M, v, vf, np_)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T(K) M(kg/mol) [speed] [v_max_factor] [num_points]'")
