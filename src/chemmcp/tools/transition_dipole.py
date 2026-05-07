"""
跃迁偶极矩计算工具 (MCP #477)。
计算电子/振动跃迁的跃迁偶极矩、振子强度、Einstein系数和辐射寿命。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# ===== 物理常数 =====
H = 6.62607015e-34         # J·s
HBAR = 1.054571817e-34     # J·s
C = 2.99792458e8            # m/s
E_CHARGE = 1.602176634e-19  # C
ME = 9.1093837015e-31      # kg (电子质量)
EPSILON_0 = 8.8541878128e-12 # F/m
DEBYE = 3.33564e-30        # C·m
NA = 6.02214076e23          # mol⁻¹
EV_TO_J = E_CHARGE          # J/eV


@ChemMCPManager.register_tool
class TransitionDipole(BaseTool):
    """
    跃迁偶极矩与光谱强度计算。
    
    功能:
      - 跃迁偶极矩 |μ_if| = |⟨ψ_f|μ̂|ψ_i⟩|
      - 振子强度 f = (8π² m_e ν / (3h e²)) |μ_if|²
      - Einstein 自发发射系数 A_if = (16π³ ν³ / (3ε₀ h c³)) |μ_if|²
      - 吸收截面 σ(ν) = (π e² / (ε₀ m_e c)) f · lineshape
      - 辐射寿命 τ = 1/A_if
      - 支持常见原子/分子跃迁数据库查询
    """
    __version__ = "0.1.0"
    name = "TransitionDipole"
    func_name = "calculate_transition_dipole"
    description = "Calculate transition dipole moment, oscillator strength, absorption cross-section, Einstein A/B coefficients, and radiative lifetime for electronic and vibrational transitions."
    implementation_description = (
        "Computes transition dipole moment μ_if from wavefunction overlap or experimental data. "
        "Derives all related spectroscopic quantities: f-number, σ, ε, A_if, B_if, τ. "
        "Includes database of common transitions for reference. "
        "Key relations: f ∝ ν|μ|², A ∝ ν³|μ|², σ ∝ πe²f/(ε₀mc)·g(ν)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Transition Dipole", "Oscillator Strength", "Spectroscopy", "Electronic Transitions", "Einstein Coefficients"]
    required_envs = []

    code_input_sig = [
        ("transition_energy_eV", "float", "N/A", "Transition energy in eV. Alternatively use transition_energy_cm-1."),
        ("transition_energy_cm-1", "float", "None", "Transition energy in wavenumbers (cm⁻¹). Overrides eV if provided."),
        ("dipole_moment_D", "float", "None", "Transition dipole moment magnitude in Debye. Provide this OR compute from quantum numbers."),
        ("dipole_moment_Cm", "float", "None", "Transition dipole moment in C·m (alternative unit)."),
        ("oscillator_strength_f", "float", "None", "If known directly: provide f-value instead of μ."),
        ("initial_state", "dict", "None", "Initial state info: {'n':int, 'l':int, 'type':'atomic'/'molecular'}. For atomic H-like calculation."),
        ("final_state", "dict", "None", "Final state info (same format as initial_state)."),
        ("temperature_K", "float", "298.15", "Temperature for Boltzmann population factors."),
        ("compute_all", "bool", "True", "Whether to compute all derived quantities (A, B, σ, ε, τ)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'energy_eV mu_D'. Example: '3.4 0.5' or 'energy_cm=25000 mu_D=1.2' or lookup 'H Lyman-alpha'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing transition dipole moment, oscillator strength, cross-section, Einstein coefficients, lifetime, and analysis."),
    ]

    examples = [
        {
            "code_input": {
                "transition_energy_eV": 10.2,
                "dipole_moment_D": 2.5,
            },
            "text_input": {"input_params": "10.2 2.5"},
            "output": {"result": {
                "transition_energy_eV": 10.2,
                "oscillator_strength_f": "...",
                "Einstein_A_s-1": "...",
                "lifetime_ns": "...",
            }},
        },
        {
            "code_input": {
                "initial_state": {"n": 1, "l": 0},
                "final_state": {"n": 2, "l": 1},
            },
            "text_input": {"input_params": "H 1s->2p"},
            "output": {"result": {"transition_type": "H atom 1s→2p"}},
        },
    ]

    # ===== 常见跃迁数据库 =====
    TRANSITION_DB = {
        # 原子跃迁 (能量eV, μ/D, f, 描述)
        "H Lyman-alpha":   (10.20, 1.50, 0.416, "H 1s → 2p, UV 121.6nm"),
        "H Lyman-beta":    (12.09, 0.80, 0.079, "H 1s → 3p, UV 102.6nm"),
        "H Balmer-alpha":  (1.89,  1.02, 0.641, "H 2s/2p → 3s/d, Vis 656.3nm (Hα)"),
        "Na D-line":       (2.10,  2.98, 0.645, "Na 3s → 3p, Yellow 589nm doublet"),
        "K D-line":        (1.61,  2.85, 0.698, "K 4s → 4p, IR 766.5nm"),
        "Rb D1":           (1.59,  3.12, 0.694, "Rb 5s → 5p₁/₂, 794.98nm"),
        "Cs D1":           (1.44,  3.25, 0.718, "Cs 6s → 6p₁/₂, 894.35nm"),
        "He 1¹S₀→1¹P₁":   (21.22, 0.60, 0.276, "He resonance line, 58.43nm"),
        "Ne 3s→3p":       (19.0,  0.95, 0.30, "Ne red line ~640nm region"),
        "Hg 253.7nm":      (4.88,  0.32, 0.027, "Hg 6s²6p² → 6s²6p6d, UV resonance"),
        # 分子跃迁
        "formaldehyde n→π*": (3.65, 0.45, 0.002, "H₂CO S₀→S₁, ~340nm, weak n→π*"),
        "benzene π→π*":     (4.72, 2.10, 0.089, "C₆H₆ ¹A₁g→¹B₂u, ~260nm, symmetry-forbidden"),
        "acetone n→π*":     (4.41, 0.38, 0.0015, "(CH₃)₂CO, ~280nm"),
        "formaldehyde π→π*": (7.8, 1.55, 0.08, "H₂CO S₀→S₂, ~160nm"),
        "HCl V=0→1":        (0.36,  0.12, 3e-5, "HCl fundamental vibration, 2886cm⁻¹"),
        "CO V=0→1":         (0.27,  0.10, 2e-5, "CO stretch, 2143cm⁻¹"),
        "I₂ B←X":           (1.97,  1.05, 0.034, "I₂ visible transition, ~500nm"),
        "O₂ b¹Σ_g⁺←X³Σ_g⁻": (1.96, 0.001, 2e-9, "O₂ atmospheric band, magnetic dipole"),
        "chlorophyll-a Qy":  (1.83,  3.80, 0.10, "Chl-a red band ~680nm, strong"),
        "retinal (rhodopsin)": (2.38, 4.52, 0.18, "11-cis retinal S₀→S₁, ~520nm"),
        "DNA base π→π*":    (4.50, 1.80, 0.065, "Nucleobase average, ~260nm"),
        "tryptophan ¹La":    (4.31, 2.85, 0.11, "Trp indole S₀→S₁, ~280nm"),
        "tyrosine":          (4.42, 1.62, 0.055, "Tyr phenol S₀→S₁, ~275nm"),
        "phenylalanine":     (4.78, 1.48, 0.046, "Phe benzene S₀→S₁, ~260nm"),
        "Rhodamine 6G S₀→S₁": (2.33, 5.20, 0.24, "Fluorescent dye, ~530nm"),
        "fluorescein":        (2.63, 4.85, 0.21, "Common fluorophore, ~470nm"),
        "GFP chromophore":   (2.67, 3.92, 0.14, "Green fluorescent protein, ~465nm excitation"),
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = H
        self.hbar = HBAR
        self.c = C
        self.e = E_CHARGE
        self.me = ME
        self.eps0 = EPSILON_0
        self.D = DEBYE

    def _run_base(self, transition_energy_eV: float = None,
                  transition_energy_cm_minus_1: float = None,
                  dipole_moment_D: float = None, dipole_moment_Cm: float = None,
                  oscillator_strength_f: float = None,
                  initial_state: dict = None, final_state: dict = None,
                  temperature_K: float = 298.15,
                  compute_all: bool = True) -> dict:
        """
        核心计算逻辑。
        """
        # ---- 确定跃迁能量 ----
        if transition_energy_cm_minus_1 is not None:
            nu_cm = transition_energy_cm_minus_1
            nu_eV = nu_cm * self.c * 100 * self.h / self.e  # cm⁻¹ → eV
            nu_Hz = nu_cm * self.c * 100                     # Hz
            wavelength_nm = 1e7 / nu_cm if nu_cm > 0 else 0   # nm
        elif transition_energy_eV is not None:
            nu_eV = transition_energy_eV
            nu_Hz = nu_eV * self.e / self.h                   # Hz
            nu_cm = nu_Hz / (self.c * 100)                    # cm⁻¹
            wavelength_nm = self.c / nu_Hz * 1e9              # nm
        else:
            raise ChemMCPInputError("Must provide either transition_energy_eV or transition_energy_cm-1")

        if nu_eV <= 0:
            raise ChemMCPInputError(f"Transition energy must be positive, got {nu_eV} eV")

        result_base = {
            "transition_energy_eV": round(nu_eV, 6),
            "transition_energy_cm-1": round(nu_cm, 4),
            "wavelength_nm": round(wavelength_nm, 4),
            "frequency_THz": round(nu_Hz / 1e12, 4),
            "temperature_K": temperature_K,
        }

        # ---- 确定跃迁偶极矩 ----
        if oscillator_strength_f is not None:
            # 从 f 反推 μ
            f = oscillator_strength_f
            if nu_Hz > 0:
                mu_sq_Cm2 = (3 * self.h * self.e**2 / (8 * math.pi**2 * self.me * nu_Hz)) * f
                mu_mag_Cm = math.sqrt(max(0, mu_sq_Cm2))
                mu_mag_D = mu_mag_Cm / self.D
            else:
                mu_mag_D = 0.0; mu_mag_Cm = 0.0
        elif dipole_moment_Cm is not None:
            mu_mag_Cm = dipole_moment_Cm
            mu_mag_D = mu_mag_Cm / self.D
            f = self._compute_f(mu_mag_Cm, nu_Hz)
        elif dipole_moment_D is not None:
            mu_mag_D = dipole_moment_D
            mu_mag_Cm = mu_mag_D * self.D
            f = self._compute_f(mu_mag_Cm, nu_Hz)
        elif initial_state is not None and final_state is not None:
            # 尝试从量子数计算（类氢原子近似）
            mu_mag_D, f = self._estimate_from_states(initial_state, final_state, nu_eV)
            mu_mag_Cm = mu_mag_D * self.D
        else:
            raise ChemMCPInputError(
                "Must provide one of: dipole_moment_D, dipole_moment_Cm, "
                "oscillator_strength_f, or initial+final state quantum numbers"
            )

        result_base["transition_dipole_moment_D"] = round(mu_mag_D, 6)
        result_base["transition_dipole_moment_Cm"] = round(mu_mag_Cm, 30)
        result_base["oscillator_strength_f"] = round(f, 8)

        # ---- 计算所有派生量 ----
        derived = {}
        if compute_all and nu_Hz > 0:
            derived = self._compute_derived(mu_mag_Cm, nu_eV, nu_Hz, nu_cm, wavelength_nm)

        # ---- 数据库查找 ----
        db_match = self._find_db_match(nu_eV, mu_mag_D, f)

        result = {
            **result_base,
            "derived_quantities": derived,
            "database_reference": db_match,
            "analysis": self._analyze_transition(f, mu_mag_D, nu_eV),
        }

        logger.info(f"TransitionDipole: E={nu_eV:.3f}eV, λ={wavelength_nm:.1f}nm, "
                     f"μ={mu_mag_D:.3f}D, f={f:.6f}")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            s = input_params.strip()
            # 检查是否为数据库查询
            for key in self.TRANSITION_DB:
                if key.lower() in s.lower():
                    return self._lookup_transition(key)

            # 数值输入: energy(eV) mu(D)
            parts = s.split()
            if len(parts) >= 2:
                E = float(parts[0])
                mu = float(parts[1])
                return self._run_base(transition_energy_eV=E, dipole_moment_D=mu)
            else:
                raise ChemMCPInputError(f"Cannot parse: {s}. Use 'E_eV mu_D' or a transition name.")
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")

    def _lookup_transition(self, name: str) -> dict:
        """从数据库查询跃迁。"""
        if name not in self.TRANSITION_DB:
            raise ChemMCPError(f"Transition '{name}' not found in database. "
                             f"Available: {list(self.TRANSITION_DB.keys())[:10]}...")
        E, mu, f, desc = self.TRANSITION_DB[name]
        return self._run_base(transition_energy_eV=E, dipole_moment_D=mu)

    @staticmethod
    def _compute_f(mu_Cm: float, nu_Hz: float) -> float:
        """从 μ 和 ν 计算振子强度 f。"""
        H = 6.62607015e-34
        ME = 9.1093837015e-31
        EC = 1.602176634e-19
        if nu_Hz <= 0:
            return 0.0
        return (8 * math.pi**2 * ME * nu_Hz / (3 * H * EC**2)) * mu_Cm**2

    def _compute_derived(self, mu_Cm: float, nu_eV: float, nu_Hz: float,
                         nu_cm: float, wl_nm: float) -> dict:
        """计算所有派生光谱量。"""
        mu_sq = mu_Cm**2

        # Einstein A (自发发射速率) [s⁻¹]
        # A = (16π³ν³)/(3ε₀hc³) × |μ|²
        A_if = (16 * math.pi**3 * nu_Hz**3) / (3 * self.eps0 * self.h * self.c**3) * mu_sq

        # Einstein B [J⁻¹m³s⁻²] 或 [m³/J·s²]
        # B = (π/(3ε₀ħ²)) × |μ|²
        B_if = math.pi / (3 * self.eps0 * self.hbar**2) * mu_sq

        # 辐射寿命 τ = 1/A (如果有多个通道需求和)
        tau_s = 1.0 / A_if if A_if > 0 else float('inf')
        tau_ns = tau_s * 1e9

        # 吸收截面 σ [cm²] (峰值，假设 δ 函数线型)
        # σ_peak = (πe²/(ε₀m_ec)) × f ≈ 0.0265 × f [cm²] (当 E[eV] 时需要修正)
        f_val = self._compute_f(mu_Cm, nu_Hz) if nu_Hz > 0 else 0.0
        sigma_cm2 = (math.pi * self.e**2 / (self.eps0 * self.me * self.c)) * f_val * 1e4  # m²→cm²

        # 摩尔消光系数 ε [M⁻¹cm⁻¹]
        # ε = N_A σ / (1000 ln(10)) ≈ 8.73×10²¹ σ[cm²]
        epsilon_Mcm = NA * sigma_cm2 / (1000 * math.log(10))

        # 积分吸收强度 S = ∫ε dν̃ [M⁻¹cm⁻²]
        S_integrated = epsilon_Mcm * (self._typical_width_cm(nu_cm))  # 近似

        return {
            "Einstein_A_s-1": round(A_if, 4),
            "Einstein_B_J-1m3s-2": round(B_if, 10),
            "radiative_lifetime_s": f"{tau_s:.4e}" if tau_s < 1e10 else "∞",
            "radiative_lifetime_ns": round(tau_ns, 4) if tau_ns < 1e12 else float('inf'),
            "absorption_cross_section_peak_cm2": round(sigma_cm2, 4),
            "molar_absorptivity_epsilon_M-1cm-1": round(epsilon_Mcm, 4),
            "integrated_absorption_M-1cm-2": round(S_integrated, 4),
            "peak_absorption_OD_for_1uM_1cm": round(epsilon_Mcm * 1e-6, 4),  # OD for 1μM, 1cm path
        }

    def _estimate_from_states(self, init: dict, final: dict, nu_eV: float):
        """从类氢原子量子数估算跃迁偶极矩。"""
        ni = init.get("n", 1); li = init.get("l", 0)
        nf = final.get("n", 2); lf = final.get("l", 1)

        # 类氢原子跃迁偶极矩近似
        # |⟨n'l'|r|nl⟩| ≈ a₀ × f_nl,n'l'
        a0_B = 0.529177  # Bohr radius Å
        if li == 0 and lf == 1:  # s→p
            if ni == 1 and nf == 2:
                mu_D = 1.29  # H 1s→2p 精确值
            else:
                r_avg = a0_B * (ni*nf * math.sqrt(ni**2 + nf**2) / 2) * 0.8
                mu_D = r_avg * self.e * self.a0 / self.D  # e·Å → D
        else:
            mu_D = 0.5  # 默认估计

        f_val = self._compute_f(mu_D * self.D, nu_eV * self.e / self.h)
        return mu_D, f_val

    def _find_db_match(self, nu_eV: float, mu_D: float, f: float) -> dict:
        """在数据库中查找最接近的跃迁。"""
        best_name = None
        best_score = float('inf')
        for name, (db_E, db_mu, db_f, desc) in self.TRANSITION_DB.items():
            score = abs(db_E - nu_eV)/max(db_E, 0.1) + abs(db_mu - mu_D)/max(db_mu, 0.01) + abs(db_f - f)/max(abs(f), 0.001)
            if score < best_score:
                best_score = score
                best_name = name

        if best_name and best_score < 2.0:
            E, mu, f_db, desc = self.TRANSITION_DB[best_name]
            return {
                "closest_match": best_name,
                "description": desc,
                "db_energy_eV": E,
                "db_mu_D": mu,
                "db_f": f_db,
                "match_quality_score": round(best_score, 4),
            }
        return {"closest_match": None, "note": "No close match in database"}

    @staticmethod
    def _analyze_transition(f: float, mu_D: float, nu_eV: float) -> dict:
        """分析跃迁特征。"""
        if f > 1.0:
            strength = "very strongly allowed (super-allowed)"
        elif f > 0.1:
            strength = "strongly allowed (fully allowed electric dipole)"
        elif f > 0.01:
            strength = "moderately allowed"
        elif f > 0.001:
            strength = "weakly allowed (partially forbidden)"
        elif f > 1e-6:
            strength = "forbidden (spin/orbit/symmetry forbidden, gains intensity via vibronic coupling)"
        else:
            strength = "strictly forbidden (magnetic dipole, electric quadrupole, or higher multipole)"

        region = (
            "γ-ray" if nu_eV > 1e6 else
            "X-ray" if nu_eV > 1e3 else
            "vacuum UV" if nu_eV > 10 else
            "UV" if nu_eV > 3.1 else
            "visible" if nu_eV > 1.6 else
            "near-IR" if nu_eV > 0.5 else
            "mid-IR" if nu_eV > 0.05 else
            "far-IR/microwave" if nu_eV > 0.001 else
            "radio"
        )

        return {
            "strength_classification": strength,
            "spectral_region": region,
            "allowedness": "allowed" if f > 0.01 else "weakly allowed" if f > 1e-4 else "forbidden",
        }

    @staticmethod
    def _typical_width_cm(nu_cm: float) -> float:
        """估算典型谱带半高全宽 (cm⁻¹)。"""
        if nu_cm > 20000:  # UV
            return 3000
        elif nu_cm > 15000:  # near UV
            return 2000
        elif nu_cm > 10000:  # visible
            return 1500
        elif nu_cm > 5000:   # near IR
            return 500
        else:  # IR
            return 20
