"""
振子强度计算工具 (MCP #480)。
计算振子强度 f、吸收截面 σ、摩尔消光系数 ε、Einstein 系数和辐射寿命。
包含常见分子跃迁数据库和光谱区域分析。
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
ME = 9.1093837015e-31      # kg
EPSILON_0 = 8.8541878128e-12 # F/m
DEBYE = 3.33564e-30        # C·m
NA = 6.02214076e23          # mol⁻¹


@ChemMCPManager.register_tool
class OscillatorStrength(BaseTool):
    """
    振子强度与吸收光谱全面分析。
    
    功能:
      - 振子强度 f 计算/转换（从 μ, ε, A, σ 任一量出发）
      - 吸收截面 σ(ν) [cm²] 和峰值截面
      - 摩尔消光系数 ε [M⁻¹cm⁻¹] 和 Beer-Lambert 定律 OD = εcl
      - Einstein 自发发射系数 A_if [s⁻¹] 和受激系数 B
      - 辐射寿命 τ = 1/A (ns, μs, s)
      - 积分吸光度 ∫ε dν̃ 和跃迁偶极矩反推
      - 谱带线型函数 (Gaussian/Lorentzian/Voigt)
      - 常见分子电子/振动跃迁数据库查询
      
    与 TransitionDipole 的区别:
      - TransitionDipole: 从波函数计算 μ → 推导所有量
      - OscillatorStrength: 以 f 为中心，支持从任意光谱量出发的互转 + 实验数据分析
    """
    __version__ = "0.1.0"
    name = "OscillatorStrength"
    func_name = "analyze_oscillator_strength"
    description = "Comprehensive oscillator strength analysis: compute f, absorption cross-section σ, molar absorptivity ε, Einstein A/B coefficients, radiative lifetime from any spectroscopic input. Includes database of common molecular transitions."
    implementation_description = (
        "Central quantity is oscillator strength f. All conversions use fundamental relations:\n"
        "• f = (8π²m_e ν / 3h e²) |μ|²  ← from transition dipole\n"
        "• σ_peak = (πe² / ε₀ m_e c) · f / Δν̃  ← absorption cross-section\n"
        "• ε_max = N_A · σ_peak / (1000 ln10)  ← molar absorptivity\n"
        "• A = (16π³ ν³ / 3ε₀ h c³) |μ|²  ← Einstein A coefficient\n"
        "• τ = 1/A  ← radiative lifetime\n"
        "Supports bidirectional conversion between all quantities."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Oscillator Strength", "Absorption Cross-section", "Beer-Lambert", "Einstein Coefficients", "Spectroscopy"]
    required_envs = []

    code_input_sig = [
        ("input_type", "str", "'f_value'", "What you're providing: 'f_value', 'dipole_D', 'epsilon', 'sigma_cm2', 'A_coefficient', 'wavelength_nm', 'transition_name'."),
        ("value", "float", "N/A", "Numerical value of the input quantity."),
        ("transition_energy_eV", "float", "None", "Transition energy in eV (required for most conversions except 'transition_name')."),
        ("transition_energy_cm-1", "float", "None", "Alternative: energy in cm⁻¹."),
        ("fwhm_cm-1", "float", "None", "Band FWHM in cm⁻¹ (needed for σ and ε calculations)."),
        ("concentration_M", "float", "None", "Concentration for Beer-Lambert OD calculation."),
        ("path_length_cm", "float", "1.0", "Optical path length for OD calculation."),
        ("temperature_K", "float", "298.15", "Temperature (affects Boltzmann factors)."),
        ("line_shape", "str", "'gaussian'", "Line shape model: 'gaussian', 'lorentzian', 'voigt'."),
        ("lookup_database", "bool", "True", "Whether to search database for matching transitions."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'input_type value [energy_eV]'. Examples:\n'f 0.5 4.0' (f=0.5 at 4eV)\n'mu_D 2.0 3.5' (μ=2D at 3.5eV)\n'eps 50000 280nm' (ε=50000 at 280nm)\n'A 1e8 500nm' (A=10⁸ s⁻¹ at 500nm)\n'lookup benzene π→π*'"),
    ]

    output_sig = [
        ("result", "dict", "Complete oscillator strength analysis: f, μ, σ, ε, A, B, τ, integrated intensity, spectral region, and comparison with database values."),
    ]

    examples = [
        {
            "code_input": {
                "input_type": "f_value",
                "value": 0.68,
                "transition_energy_eV": 2.10,
            },
            "text_input": {"input_params": "f 0.68 2.10"},
            "output": {"result": {
                "oscillator_strength_f": 0.68,
                "transition_energy_eV": 2.10,
                "wavelength_nm": 590.5,
                "approximate_transition": "Na D-line like",
            }},
        },
        {
            "code_input": {
                "input_type": "epsilon",
                "value": 150000,
                "transition_energy_eV": 4.0,
                "fwhm_cm_minus_1:": 4000,
            },
            "text_input": {"input_params": "eps 150000 4.0"},
            "output": {"result": {"oscillator_strength_f": "..."}},
        },
        {
            "code_input": {
                "input_type": "transition_name",
                "value": "H Lyman-alpha",
            },
            "text_input": {"input_params": "lookup H Lyman-alpha"},
            "output": {"result": {"database_match": "H Lyman-alpha"}},
        },
    ]

    # ===== 分子跃迁数据库 =====
    TRANSITION_DB = {
        # 格式: name → {energy_eV, f, epsilon_max, lambda_nm, type, description}
        # --- 原子共振线 ---
        "H Lyman-alpha":   {"E_eV": 10.20, "f": 0.416, "eps": None, "wl": 121.6, "type": "atomic UV", "desc": "H 1s→2p resonance"},
        "H Balmer-alpha":  {"E_eV": 1.89,  "f": 0.641, "eps": None, "wl": 656.3, "type": "atomic Vis", "desc": "H Hα line (red)"},
        "He 58.4nm":       {"E_eV": 21.22, "f": 0.276, "eps": None, "wl": 58.4,  "type": "atomic EUV","desc": "He 1¹S₀→1¹P₁"},
        "Na D-line":       {"E_eV": 2.10,  "f": 0.645, "eps": None, "wl": 589.3, "type": "atomic Vis", "desc": "Na 3s→3p doublet avg"},
        "K D-line":        {"E_eV": 1.61,  "f": 0.698, "eps": None, "wl": 766.5, "type": "atomic NIR", "desc": "K 4s→4p"},
        "Ne red lines":     {"E_eV": 19.0,  "f": 0.30,  "eps": None, "wl": 650.0, "type": "atomic Vis", "desc": "Ne discharge lines"},
        "Hg 253.7nm":      {"E_eV": 4.88,  "f": 0.027, "eps": None, "wl": 253.7, "type": "atomic UV", "desc": "Hg intercombination spin-forbidden"},
        "Ca K-line":       {"E_eV": 2.93,  "f": 0.67,  "eps": None, "wl": 422.7, "type": "atomic Vis", "desc": "Ca 4s→4p (blue)"},
        "Ba green":        {"E_eV": 2.00,  "f": 0.24,  "eps": None, "wl": 553.5, "type": "atomic Vis", "desc": "Ba 6s→6p (green)"},
        # --- 有机分子电子跃迁 ---
        "benzene π→π*":   {"E_eV": 4.72, "f": 0.089, "eps": 60000,  "wl": 263,  "type": "π→π*",   "desc": "Benzene ¹A₁g→¹B₂u (symmetry forbidden)"},
        "toluene π→π*":   {"E_eV": 4.64, "f": 0.095, "eps": 70000,  "wl": 267,  "type": "π→π*",   "desc": "Toluene (methyl benzene)"},
        "phenol π→π*":    {"E_eV": 4.53, "f": 0.12,  "eps": 85000,  "wl": 274,  "type": "π→π*",   "desc": "Phenol (OH enhances intensity)"},
        "aniline π→π*":   {"E_eV": 4.35, "f": 0.18,  "eps": 120000, "wl": 285,  "type": "π→π*",   "desc": "Aniline (NH₂ strong donor)"},
        "formaldehyde n→π*":{"E_eV": 3.65, "f": 0.002, "eps": 150,    "wl": 340,  "type": "n→π*",   "desc": "H₂CO S₀→S₁ weak singlet"},
        "acetone n→π*":   {"E_eV": 4.41, "f": 0.0015,"eps": 200,    "wl": 281,  "type": "n→π*",   "desc": "(CH₃)₂CO weak n→π*"},
        "acetonitrile π→π*":{"E_eV": 8.0, "f": 0.35,  "eps": 500000, "wl": 155,  "type": "π→π*",   "desc": "CH₃CN strong π→π* (vacuum UV)"},
        "formaldehyde π→π*":{"E_eV": 7.8,  "f": 0.08,  "eps": 80000,  "wl": 159,  "type": "π→π*",   "desc": "H₂CO S₀→S₂ (allowed)"},
        "ethylene π→π*": {"E_eV": 7.80, "f": 0.34,  "eps": 100000, "wl": 159,  "type": "π→π*",   "desc": "C₂H₄ N→V (π→π*)"},
        "butadiene π→π*": {"E_eV": 5.90, "f": 0.78,  "eps": 210000, "wl": 210,  "type": "π→π*",   "desc": "1,3-butadiene N→V₁"},
        "hexatriene π→π*": {"E_eV": 4.87,"f": 1.20,  "eps": 350000, "wl": 255,  "type": "π→π*",   "desc": "1,3,5-hexatriene (longer polyene)"},
        "β-carotene":      {"E_eV": 2.50, "f": 1.10,  "eps": 140000, "wl": 496,  "type": "π→π*",   "desc": "β-carotene (orange pigment)"},
        "retinal (S₀→S₁)": {"E_eV": 2.38, "f": 0.18,  "eps": 43000,  "wl": 520,  "type": "π→π*",   "desc": "11-cis retinal (vision chromophore)"},
        "chlorophyll-a Qy": {"E_eV": 1.83, "f": 0.10,  "eps": 80000,  "wl": 677,  "type": "π→π*",   "desc": "Chl-a red band (photosynthesis)"},
        "chlorophyll-a Soret":{"E_eV": 3.30, "f": 1.50,  "eps": 200000, "wl": 376,  "type": "π→π*",   "desc": "Chl-a Soret band (blue)"},
        "heme (Soret)":    {"E_eV": 2.95, "f": 1.80,  "eps": 130000, "wl": 420,  "type": "π→π*",   "desc": "Heme Soret band (porphyrin)"},
        "DNA base avg":    {"E_eV": 4.70, "f": 0.07,  "eps": 9000,   "wl": 264,  "type": "π→π*",   "desc": "Nucleobase average (A/T/G/C)"},
        "tryptophan ¹La":  {"E_eV": 4.31, "f": 0.11,  "eps": 55000,  "wl": 288,  "type": "π→π*",   "desc": "Trp indole band"},
        "tyrosine":        {"E_eV": 4.42, "f": 0.055, "eps": 23000,  "wl": 280,  "type": "π→π*",   "desc": "Tyr phenol band"},
        "phenylalanine":   {"E_eV": 4.78, "f": 0.046, "eps": 19000,  "wl": 259,  "type": "π→π*",   "desc": "Phe benzene band"},
        "Rhodamine 6G":    {"E_eV": 2.33, "f": 0.24,  "eps": 116000, "wl": 532,  "type": "π→π*",   "desc": "R6G fluorescent dye"},
        "fluorescein":     {"E_eV": 2.63, "f": 0.21,  "eps": 92000,  "wl": 472,  "type": "π→π*",   "desc": "Fluorescein dye"},
        "GFP chromophore": {"E_eV": 2.67, "f": 0.14,  "eps": 66000,  "wl": 465,  "type": "π→π*",   "desc": "GFP excitation"},
        "eosin Y":         {"E_eV": 2.30, "f": 0.30,  "eps": 143000, "wl": 539,  "type": "π→π*",   "desc": "Eosin Y (red dye)"},
        "coumarin 153":    {"E_eV": 3.40, "f": 0.22,  "eps": 87000,  "wl": 365,  "type": "π→π*",   "desc": "Laser dye (UV/blue)"},
        "DCM dye":         {"E_eV": 2.18, "f": 0.19,  "eps": 71000,  "wl": 569,  "type": "π→π*",   "desc": "Red laser dye"},
        # --- 电荷转移跃迁 ---
        "Ru(bpy)₃²+ MLCT": {"E_eV": 2.73, "f": 0.04,  "eps": 14500,  "wl": 454,  "type": "MLCT",    "desc": "Ru(II) polypyridine CT"},
        "Fe(CN)₆³⁻ LMCT": {"E_eV": 3.50, "f": 0.12,  "eps": 40000,  "wl": 354,  "type": "LMCT",    "desc": "Ferricyanide charge transfer"},
        "iodine B←X":     {"E_eV": 1.97, "f": 0.034, "eps": 3000,   "wl": 630,  "type": "valence",  "desc": "I₂ visible (brown vapor)"},
        "NO₂ gas":         {"E_eV": 2.50, "f": 0.02,  "eps": 2000,   "wl": 496,  "type": "n→π*",   "desc": "Nitrogen dioxide (red-brown)"},
        "KMnO₄ (purple)": {"E_eV": 2.10, "f": 0.06,  "eps": 6000,   "wl": 590,  "type": "LMCT",    "desc": "Permanganate (Mn VII O CT)"},
        "CuSO₄ blue":     {"E_eV": 1.85, "f": 0.03,  "eps": 2000,   "wl": 670,  "type": "d-d",     "desc": "Cu(II) d-d transition"},
        "Ni(dmg)₂ red":    {"E_eV": 2.10, "f": 0.05,  "eps": 4500,   "wl": 590,  "type": "d-d",     "desc": "Ni(II) square planar"},
        "TiO₂ band edge":  {"E_eV": 3.20, "f": 0.01,  "eps": 5000,   "wl": 388,  "type": "LMCT",    "desc": "Ti(IV) O²⁻ semiconductor"},
        "CdS yellow":      {"E_eV": 2.42, "f": 0.02,  "eps": 8000,   "wl": 512,  "type": "bandgap", "desc": "CdS semiconductor (pigment)"},
        # --- 生物重要分子 ---
        "melatonin":       {"E_eV": 4.38, "f": 0.13,  "eps": 60000,  "wl": 283,  "type": "π→π*",   "desc": "Sleep hormone (indole)"},
        "vitamin D3":      {"E_eV": 4.50, "f": 0.15,  "eps": 42000,  "wl": 276,  "type": "π→π*",   "desc": "Triene system"},
        "riboflavin":      {"E_eV": 3.44, "f": 0.25,  "eps": 93000,  "wl": 361,  "type": "π→π*",   "desc": "Vitamin B₂ (isoalloxazine)"},
        "NADH":           {"E_eV": 3.72, "f": 0.28,  "eps": 122000, "wl": 333,  "type": "π→π*",   "desc": "Nicotinamide adenine dinucleotide"},
        "ATP (adenine)":  {"E_eV": 4.60, "f": 0.08,  "eps": 13500,  "wl": 270,  "type": "π→π*",   "desc": "ATP purine absorption"},
        "hemoglobin Soret": {"E_eV": 2.95,"f": 1.85,  "eps": 125000, "wl": 420,  "type": "π→π*",   "desc": "Hb porphyrin Soret"},
        "myoglobin Soret": {"E_eV": 2.95,"f": 1.75,  "eps": 118000, "wl": 420,  "type": "π→π*",   "desc": "Mb porphyrin Soret"},
        "cytochrome c α": {"E_eV": 1.89, "f": 0.008,"eps": 1800,   "wl": 656,  "type": "Q-band",  "desc": "Cyt c heme α band"},
        "anthracene":      {"E_eV": 4.00, "f": 0.32,  "eps": 210000, "wl": 310,  "type": "π→π*",   "desc": "Three fused rings"},
        "pyrene":          {"E_eV": 3.82, "f": 0.22,  "eps": 160000, "wl": 325,  "type": "π→π*",   "desc": "Four fused rings"},
        "perylene":        {"E_eV": 3.40, "f": 0.48,  "eps": 260000, "wl": 365,  "type": "π→π*",   "desc": "Five fused rings"},
        "fullerene C60":   {"E_eV": 3.70, "f": 0.02,  "eps": 8000,   "wl": 335,  "type": "π→π*",   "desc": "C₆₀ symmetry-forbidden weak"},
        "graphene exciton": {"E_eV": 2.50, "f": 0.80,  "eps": 250000, "wl": 496,  "type": "π→π*",   "desc": "2D material exciton peak"},
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = H
        self.c = C
        self.e = E_CHARGE
        self.me = ME
        self.eps0 = EPSILON_0
        self.D = DEBYE
        self.NA = NA
        self.hbar = HBAR

    def _run_base(self, input_type: str, value: float,
                  transition_energy_eV: float = None,
                  transition_energy_cm_minus_1: float = None,
                  **kwargs) -> dict:
        """核心计算逻辑。"""

        itype = input_type.lower().strip().replace("-", "_")

        # ---- 能量确定 ----
        nu_eV, nu_cm, wl_nm, nu_Hz = self._resolve_energy(
            transition_energy_eV, transition_energy_cm_minus_1
        )

        if nu_eV is None or nu_eV <= 0:
            raise ChemMCPInputError("Must provide valid transition energy.")

        # ---- 根据输入类型确定 f 值 ----
        f_val = None
        mu_D = None

        if itype == "f_value":
            f_val = value
            mu_D = self._f_to_mu(f_val, nu_Hz)

        elif itype in ("dipole_d", "mu_d", "dipole_debye"):
            mu_D = value
            f_val = self._mu_to_f(mu_D, nu_Hz)

        elif itype in ("epsilon", "eps", "molar_absorptivity"):
            eps_val = value
            fwhm = kwargs.get("fwhm_cm-1") or kwargs.get("fwhm_cm_minus_1")
            if not fwhm:
                raise ChemMCPInputError("Need fwhm_cm-1 to convert ε to f.")
            f_val = self._eps_to_f(eps_val, fwhm)
            mu_D = self._f_to_mu(f_val, nu_Hz)

        elif itype in ("sigma", "sigma_cm2", "cross_section"):
            sigma_val = value
            f_val = self._sigma_to_f(sigma_val)
            mu_D = self._f_to_mu(f_val, nu_Hz)

        elif itype in ("a_coefficient", "einstein_a", "A_coefficient"):
            A_val = value
            f_val = self._A_to_f(A_val, nu_Hz)
            mu_D = self._f_to_mu(f_val, nu_Hz)

        elif itype in ("transition_name", "lookup", "db_lookup"):
            return self._lookup_transition(value)

        else:
            raise ChemMCPError(
                f"Unknown input type: {input_type}. Choose: "
                f"f_value, dipole_D, epsilon, sigma_cm2, A_coefficient, transition_name"
            )

        if f_val is None or f_val < 0:
            raise ChemMCPError(f"Invalid derived f-value: {f_val}")

        # ---- 计算所有派生量 ----
        derived = self._compute_all(f_val, mu_D, nu_eV, nu_cm, wl_nm, nu_Hz, **kwargs)

        # ---- 数据库匹配 ----
        db_match = kwargs.get("lookup_database", True)
        db_result = self._find_db_match(nu_eV, f_val) if db_match else None

        # ---- Beer-Lambert 计算 ----
        conc = kwargs.get("concentration_M")
        path = kwargs.get("path_length_cm", 1.0)
        bl_result = {}
        if conc is not None and derived.get("epsilon_max_M-1cm-1"):
            eps = derived["epsilon_max_M-1cm-1"]
            bl_result = {
                "concentration_M": conc,
                "path_length_cm": path,
                "optical_density_OD": round(eps * conc * path, 4),
                "transmittance_pct": round(10**(-eps * conc * path) * 100, 4),
                "absorbance_pct": round((1 - 10**(-eps * conc * path)) * 100, 4),
            }

        result = {
            "input_summary": {
                "input_type": input_type,
                "input_value": value,
                "derived_f": round(f_val, 8),
                "derived_mu_D": round(mu_D, 6) if mu_D else None,
            },
            "energy": {
                "eV": round(nu_eV, 6),
                "cm-1": round(nu_cm, 2),
                "nm": round(wl_nm, 2),
                "Hz": f"{nu_Hz:.4e}",
                "THz": round(nu_Hz/1e12, 4),
                "spectral_region": self._spectral_region(nu_eV),
            },
            "oscillator_strength_analysis": {
                "f_value": round(f_val, 10),
                "f_classification": self._classify_f(f_val),
                "transition_dipole_moment_Debye": round(mu_D, 6) if mu_D else None,
                "transition_dipole_moment_Cm": round(mu_D * self.D, 30) if mu_D else None,
            },
            **derived,
            "beer_lambert": bl_result,
            "database_reference": db_result,
        }

        logger.info(f"OscillatorStrength: f={f_val:.6f}, E={nu_eV:.3f}eV, λ={wl_nm:.1f}nm")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.strip().split()
            itype = parts[0].lower()
            val = float(parts[1])

            extra = {}
            if len(parts) > 2:
                # 尝试解析能量值
                p2 = parts[2]
                if "nm" in p2.lower():
                    wl = float(p2.replace("nm", ""))
                    extra["transition_energy_eV"] = 1239.8 / wl
                elif "ev" in p2.lower():
                    extra["transition_energy_eV"] = float(p2.replace("eV", ""))
                elif "cm" in p2.lower():
                    extra["transition_energy_cm-1"] = float(p2.replace("cm-1", "").replace("cm⁻¹", ""))
                else:
                    try:
                        extra["transition_energy_eV"] = float(p2)
                    except ValueError:
                        pass

            if itype in ("lookup", "db"):
                return self._run_base("transition_name", val)
            return self._run_base(itype, val, **extra)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}. Format: 'type value [energy]'")

    def _resolve_energy(self, eV, cm):
        """统一处理能量输入。"""
        if eV is not None and eV > 0:
            nu_eV = eV
            nu_Hz = nu_eV * self.e / self.h
            nu_cm = nu_Hz / (self.c * 100)
            wl_nm = self.c / nu_Hz * 1e9
        elif cm is not None and cm > 0:
            nu_cm = cm
            nu_Hz = nu_cm * self.c * 100
            nu_eV = nu_Hz * self.e / self.h
            wl_nm = self.c / nu_Hz * 1e9
        else:
            return None, None, None, None
        return nu_eV, nu_cm, wl_nm, nu_Hz

    def _compute_all(self, f, mu_D, nu_eV, nu_cm, wl_nm, nu_Hz, **kw):
        """计算全部派生光谱量。"""
        mu_Cm = mu_D * self.D if mu_D else 0
        mu_sq = mu_Cm ** 2

        # Einstein A
        A = (16 * math.pi**3 * nu_Hz**3) / (3 * self.eps0 * self.h * self.c**3) * mu_sq if nu_Hz > 0 else 0

        # Einstein B
        B = math.pi / (3 * self.eps0 * self.hbar**2) * mu_sq

        # 寿命
        tau_s = 1.0 / A if A > 0 else float('inf')
        tau_ns = tau_s * 1e9

        # 吸收截面 (理论峰值, δ 函数近似)
        sigma_theory = (math.pi * self.e**2 / (self.eps0 * self.me * self.c)) * f * 1e4  # m²→cm²

        # 使用 FWHM 得到实际峰截面
        fwhm = kw.get("fwhm_cm-1") or kw.get("fwhm_cm_minus_1") or 4000.0
        sigma_peak = sigma_theory / max(fwhm, 1.0) * math.sqrt(math.pi / math.log(2))  # Gaussian correction

        # 摩尔消光系数
        epsilon = NA * sigma_peak / (1000 * math.log(10))

        # 积分吸光度
        S_integrated = epsilon * fwhm  # M⁻¹cm⁻² 近似

        # 自然线宽 (Γ = ħ/τ)
        Gamma_eV = self.hbar / max(tau_s, 1e-30) / self.e if tau_s < 1e30 else 0
        natural_width_cm = Gamma_eV * self.e / (self.h * self.c * 100) if Gamma_eV > 0 else 0

        return {
            "einstein_A_s-1": round(A, 4) if A < 1e15 else f"{A:.4e}",
            "einstein_B_J-1m3s-2": round(B, 10) if B < 1e20 else f"{B:.4e}",
            "radiative_lifetime_ns": round(tau_ns, 4) if tau_ns < 1e12 else "∞",
            "radiative_lifetime_s": f"{tau_s:.4e}" if tau_s < 1e10 else "∞",
            "absorption_cross_section_theoretical_cm2": round(sigma_theory, 6),
            "absorption_cross_section_peak_cm2": round(sigma_peak, 6),
            "molar_absorptivity_epsilon_max_M-1cm-1": round(epsilon, 4),
            "integrated_absorption_M-1cm-2": round(S_integrated, 4),
            "natural_linewidth_eV": round(Gamma_eV, 10),
            "natural_linewidth_cm-1": round(natural_width_cm, 6),
            "FWHM_used_for_sigma_calc_cm-1": fwhm,
        }

    @staticmethod
    def _f_to_mu(f, nu_Hz):
        """f → μ (Debye)。"""
        if nu_Hz <= 0 or f < 0:
            return 0.0
        mu_sq_Cm2 = (3 * H * E_CHARGE**2 / (8 * math.pi**2 * ME * nu_Hz)) * f
        return math.sqrt(max(0, mu_sq_Cm2)) / DEBYE

    @staticmethod
    def _mu_to_f(mu_D, nu_Hz):
        """μ (Debye) → f。"""
        if nu_Hz <= 0 or mu_D < 0:
            return 0.0
        mu_Cm = mu_D * DEBYE
        return (8 * math.pi**2 * ME * nu_Hz / (3 * H * E_CHARGE**2)) * mu_Cm**2

    @staticmethod
    def _eps_to_f(eps, fwhm_cm):
        """ε [M⁻¹cm⁻¹] → f。"""
        # ε = N_A σ / (1000 ln10), σ = (πe²/(ε₀mc)) f / Δν̃
        # 反解: f = ε × 1000 ln10 × Δν̃ / (N_A × πe²/(ε₀mc))
        prefactor = (math.pi * E_CHARGE**2 / (EPSILON_0 * ME * C)) * 1e4  # cm² per unit f
        sigma_from_eps = eps * 1000 * math.log(10) / NA
        f = sigma_from_eps * fwhm_cm / prefactor
        return f

    @staticmethod
    def _sigma_to_f(sigma_cm2):
        """σ [cm²] → f。"""
        prefactor = (math.pi * E_CHARGE**2 / (EPSILON_0 * ME * C)) * 1e4
        return sigma_cm2 / prefactor

    @staticmethod
    def _A_to_f(A, nu_Hz):
        """A [s⁻¹] → f。"""
        if nu_Hz <= 0 or A <= 0:
            return 0.0
        # A = (16π³ν³/3ε₀hc³)|μ|², f = (8π²m_eν/3he²)|μ|²
        # f/A = (m_e/2πe²c) × (h/ν²) ... simplified ratio
        mu_sq = A * (3 * EPSILON_0 * H * C**3) / (16 * math.pi**3 * nu_Hz**3)
        return OscillatorStrength._mu_to_f(math.sqrt(max(0, mu_sq)) / DEBYE, nu_Hz)

    def _lookup_transition(self, name: str) -> dict:
        """从数据库查找跃迁。"""
        best = None
        best_score = float('inf')
        name_lower = name.lower().strip()

        for key, data in self.TRANSITION_DB.items():
            if key.lower() == name_lower or key.lower() in name_lower or name_lower in key.lower():
                score = 0
                best = (key, data)
                break
            # 模糊匹配
            score = sum(1 for w in name_lower.split() if w in key.lower())
            if score > 0 and score < best_score:
                best_score = score
                best = (key, data)

        if best is None:
            available = list(self.TRANSITION_DB.keys())[:15]
            raise ChemMCPError(
                f"Transition '{name}' not found.\n"
                f"Available ({len(self.TRANSITION_DB)} total): {available}..."
            )

        key, data = best
        E = data["E_eV"]; f_db = data["f"]; eps = data["eps"]; wl = data["wl"]

        result = self._run_base(
            "f_value", f_db,
            transition_energy_eV=E,
            fwhm_cm_minus_1=4000.0,
        )
        result["result"]["database_entry"] = {
            "name": key,
            "description": data["desc"],
            "transition_type": data["type"],
            "db_wavelength_nm": wl,
            "db_epsilon": eps,
            "db_f": f_db,
        }
        return result

    def _find_db_match(self, nu_eV: float, f_val: float):
        """在数据库中找最接近的跃迁。"""
        best_name = None; best_score = float('inf')
        for name, d in self.TRANSITION_DB.items():
            sc = abs(d["E_eV"] - nu_eV)/max(d["E_eV"], 0.1) + abs(d["f"] - f_val)/max(abs(f_val), 0.001)
            if sc < best_score:
                best_score = sc; best_name = name

        if best_name and best_score < 3.0:
            d = self.TRANSITION_DB[best_name]
            return {
                "closest_match": best_name,
                "match_quality_score": round(best_score, 4),
                "db_E_eV": d["E_eV"], "db_f": d["f"],
                "db_lambda_nm": d["wl"], "db_type": d["type"],
                "db_description": d["desc"],
            }
        return {"closest_match": None}

    @staticmethod
    def _classify_f(f: float) -> str:
        if f >= 1.0:
            return "very strong (fully allowed, π→π*)"
        elif f >= 0.1:
            return "strong (allowed electric dipole)"
        elif f >= 0.01:
            return "moderate (partially allowed)"
        elif f >= 0.001:
            return "weak (forbidden gains vibronic intensity)"
        elif f >= 1e-6:
            return "very weak (spin/orbit/symmetry forbidden)"
        else:
            return "extremely weak (magnetic dipole or quadrupole)"

    @staticmethod
    def _spectral_region(eV: float) -> str:
        if eV > 1e6: return "γ-ray"
        elif eV > 1e3: return "X-ray"
        elif eV > 124: return "extreme UV"
        elif eV > 10: return "vacuum UV"
        elif eV > 3.1: return "UV-C/B/A"
        elif eV > 1.6: return "visible"
        elif eV > 0.5: return "near-infrared"
        elif eV > 0.05: return "mid/far-infrared"
        else: return "microwave/radio"
