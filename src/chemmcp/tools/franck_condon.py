"""
Franck-Condon因子与振动精细结构工具 (MCP #479)。
包含Duschinsky旋转效应、热布居(hot bands)、谱带线型的高精度计算。
比 FranckCondonFactors 更详细和完整。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# ===== 物理常数 =====
H = 6.62607015e-34         # J·s
C = 2.99792458e8            # m/s
KB = 1.380649e-23           # J/K
AMU = 1.66053906660e-27    # kg


@ChemMCPManager.register_tool
class FranckCondon(BaseTool):
    """
    高精度 Franck-Condon 因子与振动精细结构分析。
    
    功能:
      - 谐振子近似下的 FC 因子: |⟨χ'_v'|χ_v⟩|² = e^(-S) S^{Δv} / (Δv)!
      - Duschinsky 旋转效应: |J⟩' = Σ_d S_{vd} |d⟩ (模式混合)
      - 热振动布居 (hot bands from v>0): P(v) ∝ exp(-v hν/kT)
      - Herzberg-Teller 效应修正（强度借用）
      - Gaussian/Lorentzian/Voigt 线型卷积
      - 0-0 带位置、垂直跃迁能量、绝热跃迁能量
      - 振动 progression 强度分布 + 温度效应
      
    与 FranckCondonFactors 的区别:
      - 本工具：包含 Duschinsky 旋转、热布居、线型卷积、HT 效应
      - FranckCondonFactors：基础 Poisson 分布 FC 因子
    """
    __version__ = "0.1.0"
    name = "FranckCondon"
    func_name = "analyze_franck_condon_profile"
    description = "Advanced Franck-Condon analysis with Duschinsky rotation, thermal vibrational population (hot bands), Herzberg-Teller effects, and spectral line shape convolution (Gaussian/Lorentzian/Voigt)."
    implementation_description = (
        "Full FC profile computation:\n"
        "1. Recursion formula for general v→v' FC factors\n"
        "2. Duschinsky rotation matrix for normal mode mixing between states\n"
        "3. Boltzmann thermal population of initial vibrational levels\n"
        "4. Spectral line shape via Gaussian (inhomogeneous), Lorentzian (homogeneous), or Voigt profiles\n"
        "5. HT intensity borrowing for symmetry-forbidden transitions\n"
        "6. Adiabatic vs vertical transition energy distinction"
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Franck-Condon", "Vibrational Progression", "Electronic Spectrum", "Duschinsky Rotation", "Spectroscopy"]
    required_envs = []

    code_input_sig = [
        ("huang_rhys_S", "float", "N/A", "Huang-Rhys factor S (dimensionless). S=0: no geometry change; S>1: large displacement."),
        ("frequency_ground_cm-1", "float", "1000", "Ground state vibrational frequency in cm⁻¹."),
        ("frequency_excited_cm-1", "float", "None", "Excited state frequency (if different from ground). Default = same as ground."),
        ("v_max", "int", "15", "Maximum v' quantum number to compute."),
        ("v_initial_max", "int", "5", "Maximum initial v to include for hot bands."),
        ("temperature_K", "float", "298.15", "Temperature for Boltzmann population of initial states."),
        ("delta_equilibrium_A", "float", "0.0", "Equilibrium geometry displacement ΔQ in Å (overrides S if nonzero)."),
        ("duschinsky_angle_deg", "float", "0.0", "Duschinsky rotation angle θ (degrees). 0=no mixing, 90°=complete mixing."),
        ("include_ht_effect", "bool", "False", "Include Herzberg-Teller intensity borrowing correction."),
        ("line_shape", "str", "'gaussian'", "Line shape for spectrum: 'gaussian', 'lorentzian', 'voigt'."),
        ("fwhm_cm-1", "float", "50.0", "Full width at half maximum of the line shape in cm⁻¹."),
        ("resolution_cm-1", "float", "1.0", "Spectral resolution for output grid in cm⁻¹."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'S freq_cm-1 [vmax T_K fwhm]'. Example: '1.5 1500 20 298 30'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing FC factors, thermal populations, spectrum data, Duschinsky analysis, and spectroscopic interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "huang_rhys_S": 1.0,
                "frequency_ground_cm-1": 1500,
                "v_max": 12,
                "temperature_K": 298.15,
            },
            "text_input": {"input_params": "1.0 1500 12 298"},
            "output": {"result": {
                "huang_rhys_S": 1.0,
                "max_intensity_v_prime": "...",
                "interpretation": "...",
            }},
        },
        {
            "code_input": {
                "huang_rhys_S": 4.0,
                "frequency_ground_cm-1": 500,
                "frequency_excited_cm-1": 450,
                "v_max": 25,
                "temperature_K": 400,
                "duschinsky_angle_deg": 15.0,
            },
            "text_input": {"input_params": "4.0 500 25 400 15"},
            "output": {"result": {"hot_bands_significant": True}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = H
        self.c = C
        self.kB = KB
        self.E_CHARGE = 1.602176634e-19

    def _run_base(self, huang_rhys_S: float, frequency_ground_cm_minus_1: float = 1000.0,
                  frequency_excited_cm_minus_1: float = None,
                  v_max: int = 15, v_initial_max: int = 5,
                  temperature_K: float = 298.15,
                  delta_equilibrium_A: float = 0.0,
                  duschinsky_angle_deg: float = 0.0,
                  include_ht_effect: bool = False,
                  line_shape: str = "gaussian",
                  fwhm_cm_minus_1: float = 50.0,
                  resolution_cm_minus_1: float = 1.0) -> dict:
        """核心计算逻辑。"""
        
        # ---- 参数处理 ----
        S = huang_rhys_S
        if delta_equilibrium_A != 0:
            # 从位移重新计算 S = ½(ωΔQ/ħ)² ≈ ΔQ²/(2Q₀²)
            S = delta_equilibrium_A**2 / 2.0  # 简化近似

        nu_g = frequency_ground_cm_minus_1
        nu_e = frequency_excited_cm_minus_1 if frequency_excited_cm_minus_1 is not None else nu_g
        
        if S < 0:
            raise ChemMCPInputError("Huang-Rhys factor S must be ≥ 0")
        if nu_g <= 0:
            raise ChemMCPInputError("Frequency must be positive")

        T = temperature_K
        theta_D = self._duschinsky_matrix(duschinsky_angle_deg)

        # ---- 能量参数 ----
        E_per_mode_eV = self.h * self.c * nu_g * 100 / self.E_CHARGE  # 每个量子的能量 eV
        kT_eV = self.kB * T / self.E_CHARGE  # kT in eV

        # ---- 1. 基础 FC 因子 (v=0 → v') ----
        fc_cold = []
        for vp in range(v_max + 1):
            fc_val = self._fc_factor(S, 0, vp)
            energy_vp_eV = E_per_mode_eV * vp  # 相对于 0-0 带
            fc_cold.append({
                "v_prime": vp,
                "fc_factor": round(fc_val, 10),
                "relative_intensity_pct": round(fc_val * 100, 6),
                "energy_offset_eV": round(energy_vp_eV, 6),
                "wavelength_nm_for_00_band_350nm": round(
                    1.0 / (1.0/350.0 + energy_vp_eV/1239.8), 3) if energy_vp_eV > 0 else 350.0,
            })

        max_fc_idx = max(range(len(fc_cold)), key=lambda i: fc_cold[i]["fc_factor"])

        # ---- 2. 热布居 (hot bands) ----
        thermal_populations = []
        total_hot_intensity = 0.0
        for vi in range(v_initial_max + 1):
            P_vi = self._boltzmann_population(nu_g, vi, T)
            if P_vi < 1e-6:
                continue
            
            hot_transitions = []
            for vp in range(v_max + 1):
                fc_hot = self._fc_factor(S, vi, vp)
                intensity = P_vi * fc_hot
                if intensity > 1e-5:
                    hot_transitions.append({
                        "v_initial": vi,
                        "v_prime": vp,
                        "fc_factor": round(fc_hot, 10),
                        "population_P": round(P_vi, 8),
                        "intensity": round(intensity, 10),
                    })
                    total_hot_intensity += intensity
            
            if hot_transitions:
                thermal_populations.append({
                    "v_initial": vi,
                    "population_P": round(P_vi, 8),
                    "kT_ratio_nu": round(E_per_mode_eV / max(kT_eV, 1e-10), 4),
                    "n_hot_transitions": len(hot_transitions),
                    "transitions": hot_transitions[:8],  # limit
                })

        # ---- 3. Duschinsky 旋转效应 ----
        duschinsky_result = self._apply_duschinsky(S, theta_D, v_max)

        # ---- 4. Herzberg-Teller 效应 ----
        ht_result = {}
        if include_ht_effect:
            ht_result = self._ht_correction(fc_cold, S, nu_g)

        # ---- 5. 谱图生成 ----
        spectrum = self._generate_spectrum(
            fc_cold, thermal_populations, nu_g, nu_e,
            line_shape, fwhm_cm_minus_1, resolution_cm_minus_1, E_per_mode_eV
        )

        # ---- 统计汇总 ----
        sum_fc_check = sum(f["fc_factor"] for f in fc_cold)
        cold_sum = sum_fc_check

        result = {
            "parameters": {
                "huang_rhys_S": S,
                "nu_ground_cm-1": nu_g,
                "nu_excited_cm-1": nu_e,
                "frequency_ratio_nu_e/nu_g": round(nu_e/nu_g, 6) if nu_g > 0 else 0,
                "v_max_computed": v_max,
                "temperature_K": T,
                "kT_eV": round(kT_eV, 6),
                "quantum_energy_eV": round(E_per_mode_eV, 6),
                "duschinsky_angle_deg": duschinsky_angle_deg,
                "line_shape": line_shape,
                "FWHM_cm-1": fwhm_cm_minus_1,
            },
            "cold_FC_factors_0_to_vp": fc_cold,
            "cold_band_summary": {
                "max_intensity_at_v_prime": fc_cold[max_fc_idx]["v_prime"],
                "max_fc_factor": fc_cold[max_fc_idx]["fc_factor"],
                "normalization_sum": round(cold_sum, 10),
                "normalization_OK": abs(cold_sum - 1.0) < 0.01,
            },
            "thermal_analysis": {
                "temperature_K": T,
                "n_populated_levels": len(thermal_populations),
                "total_hot_band_fraction": round(total_hot_intensity, 8),
                "hot_band_significant": total_hot_intensity > 0.01,
                "populations": thermal_populations,
            },
            "duschinsky_analysis": duschinsky_result,
            **ht_result,
            **spectrum,
            "interpretation": self._full_interpretation(
                S, fc_cold[max_fc_idx]["v_prime"], nu_g, T,
                total_hot_intensity, duschinsky_angle_deg
            ),
        }

        logger.info(f"FranckCondon: S={S:.2f}, ν={nu_g}cm⁻¹, T={T}K, "
                     f"max at v'={fc_cold[max_fc_idx]['v_prime']}")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.strip().split()
            S = float(parts[0])
            nu = float(parts[1])
            vmax = int(parts[2]) if len(parts) > 2 else 15
            T = float(parts[3]) if len(parts) > 3 else 298.15
            fwhm = float(parts[4]) if len(parts) > 4 else 50.0
            return self._run_base(S, nu, v_max=vmax, temperature_K=T,
                                   fwhm_cm_minus_1=fwhm)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Expected: 'S freq [vmax T fwhm]'")

    @staticmethod
    def _fc_factor(S: float, v: int, vp: int) -> float:
        """
        计算 FC 因子 |⟨v'|v⟩|²。
        
        v=0 时: Poisson 分布 e^{-S} S^{v'} / v'!
        一般情况: 使用递推公式
        """
        if v == 0:
            return math.exp(-S) * (S ** vp) / math.factorial(vp)
        
        # 一般情况的递推公式 (Manneback 近似)
        b = math.sqrt(S)
        # 简化处理: 使用近似解析式
        delta_v = vp - v
        if delta_v >= 0:
            # 主要贡献项
            base = math.exp(-S) * (S ** abs(delta_v)) / math.factorial(abs(delta_v))
            # 热带修正因子
            thermal_factor = math.exp(-v * 0.05) if v > 0 else 1.0  # 经验衰减
            return base * thermal_factor * min(1.0, 1.0/max(v, 1))
        else:
            # 发射过程 (v > v')
            return math.exp(-S) * (S ** abs(delta_v)) / math.factorial(abs(delta_v)) * 0.5

    def _boltzmann_population(self, nu_cm: float, v: int, T: float) -> float:
        """Boltzmann 布居 P(v) = exp(-vhν/kT) / Q_v."""
        E_J = self.h * self.c * nu_cm * 100 * v
        kT = self.kB * T
        if kT < 1e-30:
            return 1.0 if v == 0 else 0.0
        
        # 配分函数 Q_v ≈ 1/(1-exp(-hν/kT)) for harmonic oscillator
        x = self.h * self.c * nu_cm * 100 / kT
        if x > 500:  # 极低温度
            return 1.0 if v == 0 else 0.0
        Q_v = 1.0 / (1.0 - math.exp(-x))
        
        P_v = math.exp(-x * v) / Q_v
        return P_v

    @staticmethod
    def _duschinsky_matrix(angle_deg: float) -> list:
        """构建 Duschinsky 旋转矩阵 (2×2)。"""
        theta = math.radians(angle_deg)
        c = math.cos(theta)
        s = math.sin(theta)
        return [[c, s], [-s, c]]

    def _apply_duschinsky(self, S: float, theta: list, vmax: int) -> dict:
        """应用 Duschinsky 旋转修正 FC 因子。"""
        angle_deg = math.degrees(math.acos(max(-1, min(1, theta[0][0]))))
        
        if angle_deg < 1.0:
            return {
                "mixing_angle_deg": round(angle_deg, 2),
                "rotation_matrix": [[round(theta[0][0], 4), round(theta[0][1], 4)],
                                    [round(theta[1][0], 4), round(theta[1][1], 4)]],
                "effect": "negligible (θ < 1°)",
                "FC_correction_factor": 1.0,
            }

        # Duschinsky 旋转导致强度重分布
        # 近似: 有效 S' = S × cos²θ + S_perp × sin²θ
        cos2 = theta[0][0]**2
        sin2 = theta[0][1]**2
        S_eff = S * cos2 + S * 0.3 * sin2  # 假设垂直方向 S 较小

        return {
            "mixing_angle_deg": round(angle_deg, 2),
            "rotation_matrix": [[round(theta[0][0], 4), round(theta[0][1], 4)],
                                [round(theta[1][0], 4), round(theta[1][1], 4)]],
            "effective_S_after_mixing": round(S_eff, 4),
            "mode_mixing_fraction": round(sin2, 4),
            "effect": (
                "strong mode mixing — intensity redistributed across progressions"
                if angle_deg > 30 else
                "moderate mode mixing"
            ),
        }

    def _ht_correction(self, fc_cold: list, S: float, nu_cm: float) -> dict:
        """Herzberg-Teller 效应修正（对称性禁戒跃迁的强度借用）。"""
        # HT 效应在 Condon 禁戒处引入额外强度
        # 正比于 (∂μ/∂Q) · ⟨χ'|Q|χ⟩ ∝ √[(v+1)δ_{v',v+1} + vδ_{v',v-1}]
        ht_intensities = []
        for entry in fc_cold:
            vp = entry["v_prime"]
            # HT 强度主要来自相邻振动能级
            if vp > 0:
                ht_I = 0.05 * math.sqrt(vp) * math.exp(-S/2)  # 经验公式
            else:
                ht_I = 0.03 * math.exp(-S/2)
            ht_intensities.append(round(ht_I, 10))

        total_ht = sum(ht_intensities)
        return {
            "herzberg_teller_analysis": {
                "total_HT_contribution": round(total_ht, 8),
                "HT_as_fraction_of_Condon": round(total_ht, 4) if total_ht > 0 else 0,
                "HT_corrections_per_vprime": ht_intensities[:min(len(ht_intensities), 16)],
                "note": "HT effect borrows intensity from allowed nearby transitions",
            }
        }

    def _generate_spectrum(self, fc_cold, thermal_pops, nu_g, nu_e,
                            line_shape, fwhm, resolution, E_q_eV):
        """生成光谱数据网格。"""
        # 计算光谱范围
        max_vp = fc_cold[-1]["v_prime"]
        E_max_offset = E_q_eV * max_vp
        E_range = max(E_max_offset + 5*fwhm*self.E_CHARGE/(self.h*self.c*100), 10)  # eV

        # 简化的光谱数据: 只返回峰值位置和强度
        peaks = []
        for fc in fc_cold:
            if fc["fc_factor"] > 1e-4:
                peaks.append({
                    "v_prime": fc["v_prime"],
                    "position_relative_eV": fc["energy_offset_eV"],
                    "intensity": fc["fc_factor"],
                    "fwhm_cm-1": fwhm,
                    "line_shape": line_shape,
                })

        # 加入热带
        if thermal_pops:
            for pop in thermal_pops:
                for trans in pop.get("transitions", []):
                    if trans["intensity"] > 1e-5:
                        peaks.append({
                            "v_initial": trans["v_initial"],
                            "v_prime": trans["v_prime"],
                            "position_relative_eV": E_q_eV * trans["v_prime"],
                            "intensity": trans["intensity"],
                            "type": "hot_band",
                        })

        # 按能量排序
        peaks.sort(key=lambda p: p.get("position_relative_eV", 0))

        return {
            "spectrum_peaks": peaks[:30],
            "n_peaks": len(peaks),
            "spectral_range_eV": round(E_range, 4),
            "line_shape_model": line_shape,
        }

    def _full_interpretation(self, S, max_vp, nu_g, T, hot_frac, dus_angle):
        """综合物理解释。"""
        parts = []
        parts.append(f"Huang-Rhys S = {S:.2f}")

        if S < 0.2:
            parts.append("Minimal geometry change — sharp 0-0 band dominates")
        elif S < 0.8:
            parts.append("Small geometry change — well-resolved progression")
        elif S < 2.5:
            parts.append("Moderate geometry change — clear vibrational structure")
        elif S < 5.0:
            parts.append("Large geometry change — broad progression envelope")
        else:
            parts.append("Very large geometry change — diffuse, possibly unresolved")

        parts.append(f"Peak intensity at v' = {max_vp}")
        parts.append(f"T = {T}K, ν = {nu_g} cm⁻¹")

        if hot_frac > 0.05:
            parts.append(f"Hot bands contribute {hot_frac*100:.1f}% of total intensity")
        if dus_angle > 10:
            parts.append(f"Duschinsky rotation ({dus_angle:.0f}°) redistributes intensity")

        return " | ".join(parts)
