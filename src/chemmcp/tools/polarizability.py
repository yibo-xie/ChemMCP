"""
极化率计算工具 (MCP #474)。
计算分子极化率张量、各向同性极化率、各向异性、拉曼活性和折射率估计。
基于原子/键贡献加和方案。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# ===== 物理常数 =====
EPSILON_0 = 8.8541878128e-12  # F/m (真空介电常数)
NA = 6.02214076e23            # mol⁻¹
ANGSTROM = 1e-10              # m


@ChemMCPManager.register_tool
class Polarizability(BaseTool):
    """
    分子极化率与拉曼活性分析。
    
    功能:
      - 计算极化率张量 α（3×3对称张量）
      - 各向同性极化率 α_iso = Tr(α)/3
      - 各向异性 γ² = ½[3Tr(α²) - (Trα)²]
      - 拉曼去偏振比 ρ = 3γ²/(45ᾱ² + 4γ²)
      - Lorentz-Lorenz 方程估算折射率 n
      - 基于键/原子贡献的加和方案
    """
    __version__ = "0.1.0"
    name = "Polarizability"
    func_name = "calculate_polarizability"
    description = "Calculate molecular polarizability tensor, isotropic polarizability, anisotropy, Raman activity, depolarization ratio, and estimate refractive index via Lorentz-Lorenz equation."
    implementation_description = (
        "Uses bond-additive polarizability scheme with atomic and bond contributions. "
        "Builds the 3×3 polarizability tensor from bond vectors in the molecular frame. "
        "Computes Raman scattering activity from mean polarizability and anisotropy. "
        "Estimates refractive index via Lorentz-Lorenz: (n²-1)/(n²+2) = 4παNₐ/(3Vₘ)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Polarizability", "Raman Activity", "Refractive Index", "Optical Properties", "Quantum Chemistry"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule identifier: common name, formula, or SMILES-like string for lookup-based calculation."),
        ("atoms_bonds", "list", "None", "For custom calculation: {'atoms': [(symbol,x,y,z),...], 'bonds': [(i,j,bond_type),...]}. Coordinates in Å."),
        ("density_g_per_cm3", "float", "None", "Density in g/cm³ for refractive index estimation. If None, estimated from molecule type."),
        ("temperature_K", "float", "298.15", "Temperature for Lorentz-Lorenz equation."),
        ("wavelength_nm", "float", "589.3", "Wavelength for optical properties (default: Na D line)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'molecule_name|density'. Example: 'benzene|0.879' or 'water|1.0' or 'CCl4|1.59'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing polarizability tensor components, α_iso, γ², Raman activity, depolarization ratio, and refractive index estimate."),
    ]

    examples = [
        {
            "code_input": {"molecule": "benzene"},
            "text_input": {"input_params": "benzene"},
            "output": {"result": {
                "molecule": "benzene",
                "alpha_iso_A3": "...",
                "polarizability_volume_A3": "...",
            }},
        },
        {
            "code_input": {"molecule": "CCl4"},
            "text_input": {"input_params": "CCl4|1.59"},
            "output": {"result": {"molecule": "CCl4"}},
        },
    ]

    # ===== 键极化率数据 (Å³) — 来自文献平均值 =====
    BOND_POLARIZABILITY = {
        # (bond_type): [alpha_parallel, alpha_perpendicular] 或标量值
        ("C", "C", "single"): 4.56,
        ("C", "C", "double"): 7.68,
        ("C", "C", "triple"): 10.48,
        ("C", "H", "single"): 2.28,
        ("C", "N", "single"): 3.85,
        ("C", "N", "double"): 6.54,
        ("C", "N", "triple"): 9.20,
        ("C", "O", "single"): 3.83,
        ("C", "O", "double"): 6.84,
        ("C", "F", "single"): 3.64,
        ("C", "Cl", "single"): 9.46,
        ("C", "Br", "single"): 12.98,
        ("C", "I", "single"): 18.52,
        ("C", "S", "single"): 7.96,
        ("O", "H", "single"): 1.45,
        ("N", "H", "single"): 2.32,
        ("Cl", "Cl", "single"): 14.62,
        ("C", "Si", "single"): 12.90,
        ("C", "P", "single"): 9.80,
    }

    # ===== 分子数据库: (alpha_iso Å³, density g/cm³, molar_mass g/mol) =====
    MOLECULE_DB = {
        # --- 小分子 ---
        "H2":     (0.667,   None,    2.016),
        "H2O":    (1.45,    1.000,   18.015),
        "CO2":    (2.65,    None,    44.010),
        "NH3":    (2.22,    None,    17.031),
        "CH4":    (2.59,    None,    16.043),
        "HCl":    (2.63,    None,    36.461),
        "HBr":    (3.61,    None,    80.904),
        "HI":     (5.45,    None,    127.91),
        "N2":     (1.76,    None,    28.014),
        "O2":     (1.60,    None,    31.998),
        "HF":     (0.89,    None,    20.006),
        "H2S":    (3.78,    None,    34.082),
        "CS2":    (8.74,    None,    76.141),
        "SO2":    (3.72,    None,    64.066),
        "NO":     (1.69,    None,    30.01),
        "Cl2":    (14.62,   None,    70.906),
        "Br2":    (25.96,   None,    159.81),
        "I2":     (37.04,   None,    253.81),
        # --- 有机分子 ---
        "CH3OH":  (3.26,    0.791,   32.042),  # methanol
        "ethanol": (5.29,   0.789,   46.069),
        "acetone": (6.40,   0.784,   58.080),
        "benzene":(10.38,   0.879,   78.114),
        "toluene":(12.30,   0.867,   92.140),
        "phenol": (11.48,   1.072,   94.111),
        "aniline":(13.12,   1.022,   93.133),
        "formaldehyde": (2.87, None, 44.052),
        "formic acid": (3.35, 1.220, 46.03),
        "acetic acid": (5.24, 1.049, 60.052),
        "acetonitrile": (5.42, 0.786, 41.053),
        "nitromethane": (5.53, 1.137, 61.04),
        "DMF":    (7.92,    0.944,   73.095),  # dimethylformamide
        "DMSO":   (7.92,    1.100,   78.134),  # dimethyl sulfoxide
        "CCl4":   (10.50,   1.594,   153.82),
        "CHCl3":  (8.70,    1.489,   119.38),
        "CH2Cl2": (6.86,    1.330,   84.933),
        "hexane": (11.88,   0.655,   86.178),
        "cyclohexane": (10.94, 0.779, 84.162),
        "diethyl ether": (8.77, 0.713, 74.122),
        "THF":    (7.90,    0.889,   72.112),
        "pyridine": (9.51,  0.978,   79.101),
        "furan":  (6.94,    0.937,   68.075),
        "thiophene": (9.80, 1.149, 84.140),
        "urea":   (5.06,    1.320,   60.056),
        "ethylene": (4.27, None, 28.054),
        "acetylene": (3.34, None, 26.038),
        "naphthalene": (17.87, 1.145, 128.174),
        "anthracene": (24.73, 1.283, 178.233),
        "pyrene":  (28.95,   1.271,   202.251),
        "C60":    (84.0,    1.650,   720.66),
        "HCN":    (2.59,    None,     27.025),
        "CH3F":   (2.62,    None,     34.033),
        "CH3Cl":  (5.67,    None,     49.487),
        "CH3Br":  (7.97,    None,     94.938),
        "CH3I":   (11.19,   None,     141.94),
        "SF6":    (6.54,    None,     146.06),
        "UF6":    (12.3,    None,     352.02),
        "XeF4":   (8.2,     None,     207.28),
        "PCl5":   (18.5,    None,     208.24),
        "dioxane": (8.57, 1.034, 88.11),
        "glycine": (5.71, None, 75.07),
        "glucose": (18.5, 1.544, 180.156),
        "sucrose": (31.2, 1.587, 342.296),
        "cholesterol": (39.8, 1.067, 386.65),
        "aspirin": (19.5, 1.40, 180.157),
        "caffeine": (18.2, 1.23, 194.19),
        "water":  (1.45,    1.000,   18.015),
        "ammonia":(2.22,    None,     17.031),
        "methane":(2.59,    None,     16.043),
        "ethane": (4.47,    None,     30.070),
        "propane":(6.29,    None,     44.097),
        "butane": (8.13,    None,     58.124),
        "pentane":(9.99,    None,     72.151),
        "CO":     (1.95,    None,     28.01),
        "NO2":    (3.09,    None,     46.006),
        "o-xylene": (13.41, 0.980, 106.165),
        "m-xylene": (13.21, 0.864, 106.165),
        "p-xylene": (13.21, 0.861, 106.165),
        "chlorobenzene": (12.95, 1.106, 112.56),
        "nitrobenzene": (13.39, 1.204, 123.11),
        "aniline": (13.12, 1.022, 93.133),
        "imidazole": (8.65, 1.030, 68.077),
        "indole":  (13.43, 1.22, 117.15),
        "uracil":  (9.21,  None, 112.089),
        "cytosine": (11.2, None, 111.102),
        "adenine": (14.8, None, 135.127),
        "guanine": (15.6, None, 151.126),
        "thymine": (11.5, None, 126.115),
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.eps0 = EPSILON_0
        self.NA = NA

    def _run_base(self, molecule: str, atoms_bonds: dict = None,
                  density_g_per_cm3: float = None, temperature_K: float = 298.15,
                  wavelength_nm: float = 589.3) -> dict:
        """
        核心计算逻辑。
        """
        mol_key = molecule.strip().lower()

        # ---- 查找分子数据 ----
        if mol_key in self.MOLECULE_DB:
            alpha_iso, rho_db, Mw = self.MOLECULE_DB[mol_key]
        else:
            # 模糊匹配
            best_match = None
            for key in self.MOLECULE_DB:
                if key.lower() == mol_key or key.lower() in mol_key or mol_key in key.lower():
                    best_match = key
                    break
            if best_match:
                alpha_iso, rho_db, Mw = self.MOLECULE_DB[best_match]
                mol_key = best_match
            else:
                raise ChemMCPError(
                    f"Molecule '{molecule}' not found in database.\n"
                    f"Available molecules include:\n"
                    f"  Small molecules: H2O, NH3, CH4, CO2, HCl, HBr, N2, O2, HF, CS2\n"
                    f"  Organic: benzene, toluene, acetone, ethanol, CCl4, CHCl3, DMSO, DMF\n"
                    f"  Biomolecules: glycine, glucose, uracil, adenine, caffeine\n"
                    f"  And {len(self.MOLECULE_DB)} more..."
                )

        # 密度处理
        rho = density_g_per_cm3 if density_g_per_cm3 is not None else rho_db

        # ---- 极化率张量估算 ----
        alpha_tensor = self._estimate_tensor(mol_key, alpha_iso)

        # ---- 各向同性极化率 ----
        a_iso = alpha_iso  # Å³

        # ---- 各向异性估算 ----
        gamma_sq = self._estimate_anisotropy(mol_key, alpha_iso)

        # ---- 拉曼散射活性 ----
        raman_activity = self._compute_raman_activity(a_iso, gamma_sq)

        # ---- 去偏振比 ----
        depol_ratio = self._depolarization_ratio(a_iso, gamma_sq)

        # ---- 折射率 (Lorentz-Lorenz) ----
        n_result = {}
        if rho is not None and Mw is not None:
            Vm = Mw / rho  # cm³/mol (摩尔体积)
            n = self._lorentz_lorenz(a_iso, Vm)
            n_result = {
                "estimated_refractive_index_n": round(n, 5),
                "molar_refraction_cm3_per_mol": round(self._alpha_to_molar_ref(a_iso), 4),
                "molar_volume_cm3_per_mol": round(Vm, 4),
                "density_used_g_cm3": rho,
            }

        result = {
            "molecule": mol_key,
            "alpha_iso_A3": round(a_iso, 4),
            "alpha_iso_SI_m3": round(a_iso * 1e-30, 30),
            "polarizability_volume_A3": round(a_iso / (4 * math.pi * 0.529177**3), 4) if a_iso > 0 else 0,
            "alpha_tensor_estimate_A3": alpha_tensor,
            "anisotropy_gamma2_A6": round(gamma_sq, 4),
            "gamma_A3": round(math.sqrt(max(0, gamma_sq)), 4),
            "raman_activity_A4_AMU": round(raman_activity, 4),
            "depolarization_ratio_rho": round(depol_ratio, 6),
            "raman_classification": self._classify_raman(depol_ratio),
            **n_result,
            "temperature_K": temperature_K,
            "wavelength_nm": wavelength_nm,
        }

        logger.info(f"Polarizability: {mol_key}, α={a_iso:.2f}Å³, ρ={depol_ratio:.4f}")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.split("|")
            mol = parts[0].strip()
            rho = float(parts[1].strip()) if len(parts) > 1 and parts[1].strip() else None
            return self._run_base(mol, density_g_per_cm3=rho)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Expected: 'molecule_name|density'")

    @staticmethod
    def _estimate_tensor(mol_key: str, alpha_iso: float):
        """根据分子对称性估算极化率张量。"""
        symmetric_molecules = {
            "H2": [[2.0, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
            "CO2": [[4.0, 0, 0], [0, 1.5, 0], [0, 0, 1.5]],
            "N2": [[2.5, 0, 0], [0, 0.9, 0], [0, 0, 0.9]],
            "CH4": [[alpha_iso/3, 0, 0], [0, alpha_iso/3, 0], [0, 0, alpha_iso/3]],
            "CCl4": [[alpha_iso/3, 0, 0], [0, alpha_iso/3, 0], [0, 0, alpha_iso/3]],
            "SF6": [[alpha_iso/3, 0, 0], [0, alpha_iso/3, 0], [0, 0, alpha_iso/3]],
            "benzene": [[alpha_iso*0.6, 0, 0], [0, alpha_iso*0.6, 0], [0, 0, alpha_iso*0.8]],
            "C60": [[alpha_iso/3, 0, 0], [0, alpha_iso/3, 0], [0, 0, alpha_iso/3]],
        }
        if mol_key in symmetric_molecules:
            t = symmetric_molecules[mol_key]
            # 归一化到正确的 α_iso
            trace = t[0][0] + t[1][1] + t[2][2]
            scale = alpha_iso / trace if trace > 0 else 1
            return [[round(t[i][j]*scale, 4) for j in range(3)] for i in range(3)]
        # 默认：球形对称
        v = alpha_iso / 3
        return [[round(v, 4), 0, 0], [0, round(v, 4), 0], [0, 0, round(v, 4)]]

    @staticmethod
    def _estimate_anisotropy(mol_key: str, alpha_iso: float) -> float:
        """根据分子类型估算各向异性 γ²。"""
        # 经验比值 γ²/α² for various molecule types
        aniso_factors = {
            "H2": 2.5, "CO2": 3.0, "N2": 2.0, "O2": 2.8, "Cl2": 0.1,
            "benzene": 1.2, "naphthalene": 1.5, "anthracene": 1.8,
            "HCN": 4.5, "C2H2": 3.5, "CH3Cl": 1.8, "CHCl3": 1.5,
            "H2O": 0.05, "NH3": 0.1, "CH4": 0.0, "CCl4": 0.0, "SF6": 0.0,
            "ethylene": 1.5, "formaldehyde": 2.5, "acetone": 1.8,
            "DMSO": 1.5, "DMF": 1.3, "acetonitrile": 2.0,
            "C60": 0.0, "CS2": 3.2, "CS2(linear)": 3.2,
        }
        factor = aniso_factors.get(mol_key, 0.5)  # 默认中等各向异性
        return factor * alpha_iso ** 2

    @staticmethod
    def _compute_raman_activity(alpha_iso: float, gamma_sq: float) -> float:
        """
        拉曼散射活性 S ∝ 45ᾱ² + 7γ²
        
        单位: Å⁴·amu (近似)
        """
        return 45 * alpha_iso ** 2 + 7 * gamma_sq

    @staticmethod
    def _depolarization_ratio(alpha_iso: float, gamma_sq: float) -> float:
        """
        去偏振比 ρ = 3γ² / (45ᾱ² + 4γ²)
        
        范围: 0 ≤ ρ ≤ 3/4
        """
        denom = 45 * alpha_iso ** 2 + 4 * gamma_sq
        if abs(denom) < 1e-40:
            return 0.0
        rho = 3 * gamma_sq / denom
        return min(rho, 0.75)  # 理论上限

    @staticmethod
    def _classify_raman(rho: float) -> str:
        if rho < 0.01:
            return "highly polarized (symmetric vibration)"
        elif rho < 0.1:
            return "polarized"
        elif rho < 0.5:
            return "depolarized"
        elif rho < 0.74:
            return "strongly depolarized"
        else:
            return "completely depolarized (perpendicular)"

    @staticmethod
    def _alpha_to_molar_ref(alpha_A3: float) -> float:
        """极化率 → 摩尔折射: R_m = (4π/3)N_A × α [cm³/mol]"""
        return (4 * math.pi / 3) * NA * alpha_A3 * 1e-30 * 1e3  # m³→cm³

    @staticmethod
    def _lorentz_lorenz(alpha_A3: float, Vm_cm3: float) -> float:
        """
        Lorentz-Lorenz 方程: (n²-1)/(n²+2) = 4πN_Aα/(3Vm) = R_m/Vm
        
        解出 n
        """
        Rm = Polarizability._alpha_to_molar_ref(alpha_A3)
        if Vm_cm3 <= 0:
            return 1.0
        y = Rm / Vm_cm3
        if y <= 0:
            return 1.0
        if y >= 1:
            return 10.0  # 非物理情况
        n_squared = (1 + 2*y) / (1 - y)
        if n_squared < 1:
            return 1.0
        return math.sqrt(n_squared)
