"""
Centrifugation Calculator — 离心参数计算工具
RCF/RPM 转换、沉降时间估算（Stokes 定律）
"""
import logging
import math
from typing import Optional, List, Dict, Any, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CentrifugationCalculator(BaseTool):
    """
    离心参数计算器：RCF 与 RPM 双向转换、沉降时间估算。
    公式：RCF = 1.118 × 10⁻⁵ × r × RPM²
    沉降时间基于 Stokes 定律：t = [9η ln(r₂/r₁)] / [2ω²(ρₚ - ρₘ)rₚ²]
    """
    __version__ = "0.1.0"
    name = "CentrifugationCalculator"
    func_name = "calculate_centrifugation"
    description = "Calculate centrifugation parameters: RCF↔RPM conversion, sedimentation time estimation using Stokes' law."
    implementation_description = "Implements RCF=1.118e-5×r×RPM² for RCF/RPM bidirectional conversion, and Stokes-law-based sedimentation time estimation accounting for particle density, medium viscosity, and rotor geometry."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Centrifugation", "RCF", "RPM", "Sample Preparation", "Separation", "Stokes Law"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "rcf_from_rpm", "Mode: 'rcf_from_rpm', 'rpm_from_rcf', 'sedimentation_time', 'comparison_table'."),
        ("rpm", "float", "0", "Rotor speed in RPM (revolutions per minute)."),
        ("radius_cm", "float", "8.0", "Effective rotor radius in cm (average of min and max radius)."),
        ("rcf", "float", "0", "Relative centrifugal force (×g)."),
        ("particle_diameter_um", "float", "1.0", "Particle diameter in μm (for sedimentation calculation)."),
        ("particle_density_g_mL", "float", "1.1", "Particle density in g/mL."),
        ("medium_density_g_mL", "float", "1.0", "Medium/solvent density in g/mL (water ≈ 1.0)."),
        ("medium_viscosity_cP", "float", "1.0", "Medium dynamic viscosity in cP (centipoise; water@20°C ≈ 1.0)."),
        ("meniscus_radius_cm", "float", "5.0", "Meniscus (liquid surface) distance from rotation axis in cm."),
        ("pellet_radius_cm", "float", "10.0", "Pellet/bottom of tube distance from rotation axis in cm."),
        ("sedimentation_path_mm", "float", "50.0", "Sedimentation path length in mm (for quick estimate mode)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated params: e.g., 'rcf_from_rpm 10000 8.0' or 'rpm_from_rcf 5000 7'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with calculated values, formulas used, and practical recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "rcf_from_rpm",
                "rpm": 10000,
                "radius_cm": 8.0,
                "rcf": 0,
                "particle_diameter_um": 1.0,
                "particle_density_g_mL": 1.1,
                "medium_density_g_mL": 1.0,
                "medium_viscosity_cP": 1.0,
                "meniscus_radius_cm": 5.0,
                "pellet_radius_cm": 10.0,
                "sedimentation_path_mm": 50.0,
            },
            "text_input": {
                "input_params": "rcf_from_rpm 10000 8.0",
            },
            "output": {
                "result": {
                    "mode": "rcf_from_rpm",
                    "rpm": 10000,
                    "radius_cm": 8.0,
                    "rcf_xg": round(1.118e-5 * 8.0 * 10000**2, 2),
                    "formula": "RCF = 1.118×10⁻⁵ × r × RPM²",
                    "note": "Typical lab centrifuge range: 100-150000×g.",
                }
            }
        },
        {
            "code_input": {
                "mode": "rpm_from_rcf",
                "rpm": 0,
                "radius_cm": 7.0,
                "rcf": 5000,
                "particle_diameter_um": 1.0,
                "particle_density_g_mL": 1.1,
                "medium_density_g_mL": 1.0,
                "medium_viscosity_cP": 1.0,
                "meniscus_radius_cm": 4.5,
                "pellet_radius_cm": 9.5,
                "sedimentation_path_mm": 50.0,
            },
            "text_input": {
                "input_params": "rpm_from_rcf 5000 7",
            },
            "output": {
                "result": {
                    "mode": "rpm_from_rcf",
                    "target_rcf_xg": 5000,
                    "radius_cm": 7.0,
                    "required_rpm": round(math.sqrt(5000 / (1.118e-5 * 7.0)), 1),
                    "formula": "RPM = √(RCF / (1.118×10⁻⁵ × r))",
                }
            }
        },
        {
            "code_input": {
                "mode": "sedimentation_time",
                "rpm": 12000,
                "radius_cm": 8.5,
                "rcf": 0,
                "particle_diameter_um": 0.5,
                "particle_density_g_mL": 1.15,
                "medium_density_g_mL": 1.0,
                "medium_viscosity_cP": 1.0,
                "meniscus_radius_cm": 5.0,
                "pellet_radius_cm": 11.0,
                "sedimentation_path_mm": 60.0,
            },
            "text_input": {
                "input_params": "sedimentation_time 12000 8.5 0.5 1.15 1.0 1.0 5.0 11.0",
            },
            "output": {
                "result": {
                    "mode": "sedimentation_time",
                    "note": "Estimated sedimentation time based on Stokes law approximation.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ── Core calculations ──────────────────────────────────────────────

    @staticmethod
    def calc_rcf(rpm: float, radius_cm: float) -> float:
        """RCF = 1.118 × 10⁻⁵ × r(cm) × RPM²"""
        return 1.118e-5 * radius_cm * rpm ** 2

    @staticmethod
    def calc_rpm(rcf: float, radius_cm: float) -> float:
        """RPM = √(RCF / (1.118 × 10⁻⁵ × r))"""
        if rcf <= 0:
            raise ChemMCPError("RCF must be positive.")
        if radius_cm <= 0:
            raise ChemMCPError("Radius must be positive.")
        return math.sqrt(rcf / (1.118e-5 * radius_cm))

    def _estimate_sedimentation_time(self, rpm: float, radius_cm: float,
                                     d_um: float, rho_p: float, rho_m: float,
                                     eta_cP: float, r_meniscus: float,
                                     r_pellet: float) -> dict:
        """
        Stokes-law-based sedimentation time estimation.
        t = [9·η·ln(r₂/r₁)] / [2·ω²·(ρₚ-ρₘ)·rₚ²]
        Returns time in minutes.
        """
        if rho_p <= rho_m:
            raise ChemMCPError(
                f"Particle density ({rho_p} g/mL) must exceed medium density "
                f"({rho_m} g/mL) for sedimentation to occur."
            )
        if d_um <= 0 or eta_cP <= 0:
            raise ChemMCPError("Particle diameter and viscosity must be positive.")

        omega = rpm * 2 * math.pi / 60  # rad/s
        r_p = d_um * 1e-6 / 2  # particle radius in meters
        eta_Pa_s = eta_cP * 1e-3  # cP → Pa·s
        rho_diff = (rho_p - rho_m) * 1e6  # g/mL → kg/m³ (×1000 for both)

        # Avoid division by zero
        if omega == 0:
            raise ChemMCPError("RPM must be > 0.")

        numerator = 9 * eta_Pa_s * math.log(r_pellet / r_meniscus)
        denominator = 2 * (omega ** 2) * rho_diff * (r_p ** 2)
        t_seconds = abs(numerator / denominator)
        t_minutes = t_seconds / 60

        return {
            "estimated_time_min": round(t_minutes, 2),
            "estimated_time_hr": round(t_minutes / 60, 2),
            "omega_rad_s": round(omega, 3),
            "particle_radius_m": r_p,
            "density_diff_kg_m3": round(rho_diff, 1),
            "stokes_parameters": {
                "omega_rad_per_s": round(omega, 3),
                "viscosity_Pa_s": eta_Pa_s,
                "log_radius_ratio": round(math.log(r_pellet / r_meniscus), 4),
            },
            "recommendation": self._time_recommendation(t_minutes),
        }

    @staticmethod
    def _time_recommendation(t_min: float) -> str:
        if t_min < 1:
            return f"Very fast sedimentation (~{t_min*60:.0f}s). Consider shorter spin or lower speed."
        elif t_min < 10:
            return f"Normal sedimentation: ~{t_min:.1f}min is appropriate."
        elif t_min < 60:
            return f"Moderate time needed: ~{t_min:.0f}min. Ensure rotor is balanced."
        elif t_min < 240:
            return f"Long spin required: ~{t_min/60:.1f}hr. Consider higher RPM or ultracentrifuge."
        else:
            return f"Very long (>4h). Use ultracentrifuge or increase particle size (flocculation)."

    def _run_base(self, mode: str = "rcf_from_rpm", rpm: float = 0.0,
                  radius_cm: float = 8.0, rcf: float = 0.0,
                  particle_diameter_um: float = 1.0,
                  particle_density_g_mL: float = 1.1,
                  medium_density_g_mL: float = 1.0,
                  medium_viscosity_cP: float = 1.0,
                  meniscus_radius_cm: float = 5.0,
                  pellet_radius_cm: float = 10.0,
                  sedimentation_path_mm: float = 50.0) -> dict:

        m = mode.lower().strip().replace("-", "_").replace(" ", "_")

        if m in ("rcf_from_rpm", "rcf"):
            if rpm <= 0:
                raise ChemMCPError("RPM must be positive.")
            if radius_cm <= 0:
                raise ChemMCPError("Radius must be positive.")
            rcf_val = self.calc_rcf(rpm, radius_cm)
            return {"result": {
                "mode": "rcf_from_rpm",
                "rpm": rpm,
                "radius_cm": radius_cm,
                "rcf_xg": round(rcf_val, 2),
                "rcf_kg": round(rcf_val / 1000, 4),
                "formula": "RCF = 1.118×10⁻⁵ × r × RPM² "
                           f"= 1.118×10⁻⁵ × {radius_cm} × {rpm}² = {round(rcf_val, 2)} ×g",
                "classification": self._classify_rcf(rcf_val),
                "safety_note": "Always balance tubes symmetrically opposite each other.",
            }}

        elif m in ("rpm_from_rcf", "rpm"):
            rpm_val = self.calc_rpm(rcf, radius_cm)
            actual_rcf = self.calc_rcf(rpm_val, radius_cm)
            return {"result": {
                "mode": "rpm_from_rcf",
                "target_rcf_xg": rcf,
                "radius_cm": radius_cm,
                "required_rpm": round(rpm_val, 1),
                "actual_rcf_at_set_rpm": round(actual_rcf, 2),
                "formula": f"RPM = √({rcf} / (1.118×10⁻⁵ × {radius_cm})) = {round(rpm_val, 1)}",
                "nearest_setting": self._nearest_rpm(rpm_val),
            }}

        elif m in ("sedimentation_time", "sedimentation", "time"):
            sed_result = self._estimate_sedimentation_time(
                rpm, radius_cm, particle_diameter_um,
                particle_density_g_mL, medium_density_g_mL,
                medium_viscosity_cP, meniscus_radius_cm, pellet_radius_cm
            )
            rcf_val = self.calc_rcf(rpm, radius_cm)
            return {"result": {
                "mode": "sedimentation_time",
                "rpm": rpm,
                "radius_cm": radius_cm,
                "current_rcf_xg": round(rcf_val, 1),
                "particle_diameter_um": particle_diameter_um,
                "particle_density_g_mL": particle_density_g_mL,
                "medium_density_g_mL": medium_density_g_mL,
                **sed_result,
            }}

        elif m in ("comparison_table", "table", "compare"):
            return self._build_comparison_table(radius_cm)

        else:
            raise ChemMCPError(
                f"Unknown mode '{mode}'. "
                "Use: 'rcf_from_rpm', 'rpm_from_rcf', 'sedimentation_time', or 'comparison_table'."
            )

    def _build_comparison_table(self, radius_cm: float) -> dict:
        """Build a reference table of common RPM ↔ RCF values."""
        rpms = [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000, 14000,
                16000, 18000, 20000, 25000, 30000, 35000, 40000, 50000, 60000,
                70000, 80000, 90000, 100000, 120000, 150000]
        rows = []
        for r in rpms:
            rows.append({"rpm": r, "rcf_xg": round(self.calc_rcf(r, radius_cm), 1)})
        return {"result": {
            "mode": "comparison_table",
            "radius_cm": radius_cm,
            "table": rows,
            "note": f"All values computed at r = {radius_cm} cm from rotation axis.",
        }}

    @staticmethod
    def _classify_rcf(rcf: float) -> str:
        if rcf < 100:
            return "Low-speed (clinical/mini centrifuge)"
        elif rcf < 10000:
            return "General-purpose benchtop centrifuge"
        elif rcf < 60000:
            return "High-speed centrifuge"
        elif rcf < 150000:
            return "Ultracentrifuge"
        else:
            return "Analytical ultracentrifuge"

    @staticmethod
    def _nearest_rpm(rpm: float) -> str:
        """Suggest nearest practical dial setting."""
        hundreds = round(rpm / 100) * 100
        thousands = round(rpm / 1000) * 1000
        return f"Set to nearest {hundreds} RPM (or {thousands} RPM if coarse dial)"

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            mode = parts[0]
            if mode in ("rcf_from_rpm", "rcf"):
                rpm = float(parts[1])
                r = float(parts[2]) if len(parts) > 2 else 8.0
                return self._run_base(mode, rpm=rpm, radius_cm=r)
            elif mode in ("rpm_from_rcf", "rpm"):
                rcf = float(parts[1])
                r = float(parts[2]) if len(parts) > 2 else 8.0
                return self._run_base(mode, rcf=rcf, radius_cm=r)
            elif mode in ("sedimentation_time", "sedimentation", "time"):
                rpm = float(parts[1])
                r = float(parts[2]) if len(parts) > 2 else 8.0
                d = float(parts[3]) if len(parts) > 3 else 1.0
                rho_p = float(parts[4]) if len(parts) > 4 else 1.1
                rho_m = float(parts[5]) if len(parts) > 5 else 1.0
                eta = float(parts[6]) if len(parts) > 6 else 1.0
                rm = float(parts[7]) if len(parts) > 7 else 5.0
                rp = float(parts[8]) if len(parts) > 8 else 10.0
                return self._run_base(mode, rpm=rpm, radius_cm=r, particle_diameter_um=d,
                                       particle_density_g_mL=rho_p, medium_density_g_mL=rho_m,
                                       medium_viscosity_cP=eta, meniscus_radius_cm=rm,
                                       pellet_radius_cm=rp)
            elif mode in ("comparison_table", "table"):
                r = float(parts[1]) if len(parts) > 1 else 8.0
                return self._run_base(mode, radius_cm=r)
            else:
                raise ValueError(f"Unknown text mode: {mode}")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
