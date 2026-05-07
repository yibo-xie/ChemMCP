"""
偶极矩计算工具 (MCP #473)。
从原子坐标和部分电荷计算分子的偶极矩矢量 μ = Σ qᵢ·rᵢ。
与 DipoleMomentEstimator（数据库查询）不同，本工具从坐标直接计算。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# ===== 物理常数 =====
DEBYE = 3.33564e-30   # C·m (1 Debye)
ANGSTROM = 1e-10       # m


@ChemMCPManager.register_tool
class DipoleMoment(BaseTool):
    """
    从原子坐标和部分电荷计算分子偶极矩。
    
    功能:
      - 计算偶极矩矢量 μ = Σ qᵢ(rᵢ - r_c) （相对于质心/电荷中心）
      - 输出大小（Debye）和方向
      - 计算四极矩张量作为高阶近似
      - 支持不同参考点（几何中心、电荷中心、质心）
      
    与 DipoleMomentEstimator 的区别：
      - 本工具：从坐标+电荷**计算**
      - DipoleMomentEstimator：从数据库**查询**实验值
    """
    __version__ = "0.1.0"
    name = "DipoleMoment"
    func_name = "calculate_dipole_moment"
    description = "Calculate molecular dipole moment vector μ = Σqᵢrᵢ from atomic partial charges and coordinates. Returns magnitude in Debye, direction cosines, and quadrupole tensor."
    implementation_description = (
        "Computes dipole moment as μ = Σ qᵢ(rᵢ - r_ref) where r_ref can be the center of mass, "
        "charge center, or origin. Also computes traceless electric quadrupole moment tensor. "
        "All coordinates internally converted from Å to SI (meters)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Dipole Moment", "Molecular Properties", "Electrostatics", "Polarity", "Quantum Chemistry"]
    required_envs = []

    code_input_sig = [
        ("atoms", "list", "N/A", "List of atoms: [(symbol, partial_charge_e, x_A, y_A, z_A), ...]. Coordinates in Ångströms."),
        ("reference_point", "str", "'charge_center'", "Reference point for dipole calculation: 'origin', 'mass_center', 'charge_center', or custom (x,y,z)."),
        ("custom_reference", "list", "None", "Custom reference point [x, y, z] in Å if reference_point='custom'."),
        ("compute_quadrupole", "bool", "True", "Whether to compute the electric quadrupole moment tensor."),
        ("unit", "str", "'Debye'", "Output unit for dipole magnitude: 'Debye' or 'C·m'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'atom_data|reference'. Example: 'O,-0.834,0,0,-0.96;H,0.417,0.757,0.586,0.19;H,0.417,-0.757,0.586,0.19|charge_center'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing dipole moment vector, magnitude, direction, and optional quadrupole tensor."),
    ]

    examples = [
        {
            "code_input": {
                "atoms": [
                    ("O", -0.834, 0.0, 0.0, -0.96),
                    ("H", 0.417, 0.757, 0.586, 0.19),
                    ("H", 0.417, -0.757, 0.586, 0.19),
                ],
                "reference_point": "charge_center",
            },
            "text_input": {"input_params": "O,-0.834,0,0,-0.96;H,0.417,0.757,0.586,0.19;H,0.417,-0.757,0.586,0.19|charge_center"},
            "output": {"result": {
                "n_atoms": 3,
                "dipole_magnitude_Debye": "...",
                "molecule": "H2O-like",
            }},
        },
        {
            "code_input": {
                "atoms": [
                    ("C", 0.0, 0.0, 0.0, 0.0),
                    ("Cl", 0.2, 1.76, 0.0, 0.0),
                    ("Cl", 0.2, -1.76, 0.0, 0.0),
                    ("Cl", 0.2, 0.0, 1.76, 0.0),
                    ("Cl", -0.2, 0.0, -1.76, 0.0),
                ],
                "reference_point": "origin",
            },
            "text_input": {"input_params": "C,0,0,0,0;Cl,0.2,1.76,0,0;Cl,0.2,-1.76,0,0;Cl,0.2,0,1.76,0;Cl,-0.2,0,-1.76,0|origin"},
            "output": {"result": {"dipole_magnitude_Debye": 0.0}},
        },
    ]

    # 原子质量 (amu)
    ATOMIC_MASS = {
        "H": 1.008, "He": 4.003, "Li": 6.941, "Be": 9.012, "B": 10.81,
        "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
        "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.086, "P": 30.974,
        "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.098, "Br": 79.904,
        "I": 126.90, "Fe": 55.845, "Cu": 63.546, "Zn": 65.38,
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.D = DEBYE
        self.a0 = ANGSTROM

    def _run_base(self, atoms: list, reference_point: str = "charge_center",
                  custom_reference: list = None, compute_quadrupole: bool = True,
                  unit: str = "Debye") -> dict:
        """
        核心计算：偶极矩 μ = Σ qᵢ(rᵢ - r_ref)
        """
        if not atoms:
            raise ChemMCPInputError("Atoms list cannot be empty.")

        n_atoms = len(atoms)

        # 解析数据
        symbols = []
        charges_e = []
        positions_A = []  # Å
        masses = []
        for at in atoms:
            sym = at[0]; q = float(at[1])
            x, y, z = float(at[2]), float(at[3]), float(at[4])
            symbols.append(sym)
            charges_e.append(q)
            positions_A.append((x, y, z))
            masses.append(self.ATOMIC_MASS.get(sym, 12.0))

        # ---- 确定参考点 ----
        if reference_point == "origin":
            ref = (0.0, 0.0, 0.0)
        elif reference_point == "mass_center":
            M = sum(masses)
            ref = tuple(sum(m * p[i] for m, p in zip(masses, positions_A)) / M for i in range(3))
        elif reference_point == "charge_center":
            Q = sum(charges_e)
            if abs(Q) < 1e-10:
                # 中性分子用电荷加权中心
                ref = (0.0, 0.0, 0.0)
            else:
                ref = tuple(sum(q * p[i] for q, p in zip(charges_e, positions_A)) / Q for i in range(3))
        elif reference_point == "custom" and custom_reference:
            ref = (float(custom_reference[0]), float(custom_reference[1]), float(custom_reference[2]))
        else:
            ref = (0.0, 0.0, 0.0)

        # ---- 计算偶极矩 ----
        # μ = Σ qᵢ(rᵢ - r_ref)，单位：e·Å → 转换为 C·m → Debye
        mux, muy, muz = 0.0, 0.0, 0.0
        for q, pos in zip(charges_e, positions_A):
            dx = pos[0] - ref[0]
            dy = pos[1] - ref[1]
            dz = pos[2] - ref[2]
            mux += q * dx
            muy += q * dy
            muz += q * dz

        # e·Å → C·m: × (1.602e-19 C/e) × (1e-10 m/Å) = 1.602e-29 C·m/e·Å
        mux_Cm = mux * 1.602176634e-29
        muy_Cm = muy * 1.602176634e-29
        muz_Cm = muz * 1.602176634e-29

        mu_mag_Cm = math.sqrt(mux_Cm**2 + muy_Cm**2 + muz_Cm**2)
        mu_mag_D = mu_mag_Cm / DEBYE

        # 方向余弦
        if mu_mag_Cm > 1e-40:
            cos_alpha = mux_Cm / mu_mag_Cm
            cos_beta = muy_Cm / mu_mag_Cm
            cos_gamma = muz_Cm / mu_mag_Cm
        else:
            cos_alpha = cos_beta = cos_gamma = 0.0

        result = {
            "n_atoms": n_atoms,
            "molecular_formula": self._guess_formula(symbols),
            "total_charge_e": round(sum(charges_e), 6),
            "reference_point_type": reference_point,
            "reference_point_A": (round(ref[0], 6), round(ref[1], 6), round(ref[2], 6)),
            "dipole_vector_eA": (round(mux, 6), round(muy, 6), round(muz, 6)),
            "dipole_vector_Cm": (round(mux_Cm, 30), round(muy_Cm, 30), round(muz_Cm, 30)),
            "dipole_magnitude_Debye": round(mu_mag_D, 6),
            "dipole_magnitude_Cm2": round(mu_mag_Cm, 30),
            "direction_cosines": (round(cos_alpha, 6), round(cos_beta, 6), round(cos_gamma, 6)),
            "polarity_classification": self._classify_polarity(mu_mag_D),
        }

        # ---- 四极矩计算 ----
        if compute_quadrupole:
            quad = self._compute_quadrupole(charges_e, positions_A, ref)
            result["quadrupole_tensor"] = quad

        logger.info(f"DipoleMoment: {self._guess_formula(symbols)}, μ={mu_mag_D:.3f}D")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.split("|")
            atom_part = parts[0].strip()
            ref = parts[1].strip() if len(parts) > 1 else "charge_center"

            atoms = []
            for a_str in atom_part.split(";"):
                vals = a_str.strip().split(",")
                atoms.append((vals[0], float(vals[1]), float(vals[2]), float(vals[3]), float(vals[4])))

            return self._run_base(atoms, ref)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")

    def _compute_quadrupole(self, charges_e, positions_A, ref):
        """
        计算电四极矩张量 Θₐᵦ = ½Σqᵢ[3(rᵢₐ - refₐ)(rᵢᵦ - refᵦ) - |rᵢ-ref|²δₐᵦ]
        单位：e·Å²
        """
        Theta = [[0.0]*3 for _ in range(3)]
        for q, pos in zip(charges_e, positions_A):
            dr = [pos[i] - ref[i] for i in range(3)]
            r2 = dr[0]**2 + dr[1]**2 + dr[2]**2
            for a in range(3):
                for b in range(3):
                    Theta[a][b] += q * (3*dr[a]*dr[b] - r2 * (1 if a == b else 0))
        for a in range(3):
            for b in range(3):
                Theta[a][b] *= 0.5
        return {
            "quadrupole_tensor_eA2": [[round(Theta[a][b], 6) for b in range(3)] for a in range(3)],
            "trace_eA2": round(Theta[0][0] + Theta[1][1] + Theta[2][2], 6),  # 应该≈0（无迹）
        }

    @staticmethod
    def _classify_polarity(mu_D: float) -> str:
        if mu_D < 0.1:
            return "nonpolar"
        elif mu_D < 1.0:
            return "slightly polar"
        elif mu_D < 2.5:
            return "polar"
        elif mu_D < 4.0:
            return "highly polar"
        else:
            return "very highly polar (ionic character)"

    @staticmethod
    def _guess_formula(symbols) -> str:
        """简单猜测分子式。"""
        from collections import Counter
        c = Counter(symbols)
        parts = []
        for elem in sorted(c.keys()):
            if c[elem] == 1:
                parts.append(elem)
            else:
                parts.append(f"{elem}{c[elem]}")
        return "".join(parts)
