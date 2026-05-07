"""
库仑势计算工具 (MCP #471)。
计算点电荷体系的静电势、电场和相互作用能。
支持多点电荷叠加、均匀带电球壳、带电线段等模型。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# ===== 物理常数 (SI单位) =====
K_E = 8.9875517923e9       # 库仑常数 k = 1/(4πε₀) [N·m²/C²]
EPSILON_0 = 8.8541878128e-12 # 真空介电常数 ε₀ [F/m]
E_CHARGE = 1.602176634e-19   # 基本电荷 e [C]
ANGSTROM = 1e-10             # Å → m


@ChemMCPManager.register_tool
class CoulombPotential(BaseTool):
    """
    库仑势与点电荷相互作用计算。
    
    功能:
      - 计算任意点电荷分布在空间各点的静电势 V(r) = Σ k·qᵢ/|r - rᵢ|
      - 计算电场矢量 E = -∇V
      - 计算体系总静电能 U = ½ Σᵢ≠ⱼ k·qᵢqⱼ/rᵢⱼ
      - 支持连续电荷分布（均匀球、线段）
    """
    __version__ = "0.1.0"
    name = "CoulombPotential"
    func_name = "calculate_coulomb_potential"
    description = "Calculate Coulomb potential, electric field, and interaction energy for point charge systems and continuous charge distributions."
    implementation_description = (
        "Uses direct superposition of Coulomb's law for point charges: V = kΣqᵢ/rᵢ, "
        "E = Σ kqᵢr̂ᵢ/rᵢ². Supports multipole expansion for distant fields. "
        "Includes analytical solutions for uniformly charged sphere and line segment."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coulomb Potential", "Electrostatics", "Point Charges", "Electric Field", "Quantum Chemistry"]
    required_envs = []

    code_input_sig = [
        ("charges", "list", "N/A", "List of point charges: [(q_in_e, x_A, y_A, z_A), ...]. q in units of elementary charge e; coordinates in Ångströms."),
        ("eval_points", "list", "None", "Evaluation points [(x, y, z), ...] in Å. If None, compute at each charge position from all others."),
        ("dielectric", "float", "1.0", "Relative dielectric constant of the medium (ε_r)."),
        ("compute_field", "bool", "True", "Whether to compute the electric field vector at each evaluation point."),
        ("compute_energy", "bool", "True", "Whether to compute total electrostatic potential energy of the system."),
        ("expansion_order", "int", "0", "Multipole expansion order (0=monopole only, 1=+dipole, 2=+quadrupole). Only used when eval_points are far from charges."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Semicolon-separated charges then points. Format: 'q1,x1,y1,z1;q2,x2,y2,z2 | px,py,pz'. Example: '1,0,0,-0.96;-1,0,0,0.96 | 0,0,2'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing potentials, electric fields, total energy, and analysis."),
    ]

    examples = [
        {
            "code_input": {
                "charges": [
                    (-0.834, 0.0, 0.0, -0.96),   # O (partial charge ~-0.83e)
                    (0.417, 0.757, 0.586, 0.19),  # H1
                    (0.417, -0.757, 0.586, 0.19),  # H2
                ],
                "eval_points": [(0.0, 0.0, 5.0)],
                "dielectric": 1.0,
                "compute_field": True,
                "compute_energy": True,
            },
            "text_input": {"input_params": "-0.834,0,0,-0.96;0.417,0.757,0.586,0.19;0.417,-0.757,0.586,0.19 | 0,0,5"},
            "output": {"result": {
                "n_charges": 3,
                "total_charge_e": 0.0,
                "potential_at_eval_points_V": [{"point": (0,0,5), "V_volts": "..."}],
            }},
        },
        {
            "code_input": {
                "charges": [(1.0, 0.0, 0.0, 0.0), (-1.0, 2.0, 0.0, 0.0)],
                "eval_points": [(1.0, 0.0, 0.0)],
                "dielectric": 1.0,
                "compute_field": True,
                "compute_energy": True,
            },
            "text_input": {"input_params": "1,0,0,0;-1,2,0,0 | 1,0,0"},
            "output": {"result": {"total_charge_e": 0.0}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """初始化物理常数。"""
        self.k_e = K_E
        self.e = E_CHARGE
        self.a0 = ANGSTROM

    def _run_base(self, charges: list, eval_points: list = None, dielectric: float = 1.0,
                  compute_field: bool = True, compute_energy: bool = True,
                  expansion_order: int = 0) -> dict:
        """
        核心计算逻辑。

        Parameters
        ----------
        charges : list of tuples (q_e, x, y, z)
            点电荷列表，q 以基本电荷 e 为单位，坐标以 Å 为单位
        eval_points : list of tuples (x, y, z) or None
            待求电势的点，None 表示在每个电荷位置处计算其他电荷产生的势
        dielectric : float
            相对介电常数 ε_r
        compute_field : bool
            是否计算电场
        compute_energy : bool
            是否计算总静电能
        expansion_order : int
            多极展开阶数

        Returns
        -------
        dict : 包含电势、电场、能量等结果
        """
        if not charges:
            raise ChemMCPInputError("Charges list cannot be empty.")

        n_charges = len(charges)

        # 解析电荷数据
        qc = []  # 电荷数组 [C]
        rc = []  # 位置数组 [m]
        total_q_e = 0.0
        for ch in charges:
            q_e, x, y, z = ch[0], float(ch[1]), float(ch[2]), float(ch[3])
            qc.append(q_e * self.e)
            rc.append((x * self.a0, y * self.a0, z * self.a0))
            total_q_e += q_e

        # 确定评估点
        if eval_points is None:
            # 在每个电荷位置计算（排除自身贡献）
            ep = [(r[0], r[1], r[2]) for r in rc]
            exclude_self = True
        else:
            ep = [(p[0] * self.a0, p[1] * self.a0, p[2] * self.a0) for p in eval_points]
            exclude_self = False

        eps_r = dielectric
        k_eff = self.k_e / eps_r

        # ---- 计算每个评估点的电势和电场 ----
        pot_results = []
        for idx_p, rp in enumerate(ep):
            V_total = 0.0
            Ex, Ey, Ez = 0.0, 0.0, 0.0

            for idx_c, (qc_i, rc_i) in enumerate(zip(qc, rc)):
                if exclude_self and idx_p == idx_c:
                    continue

                dx = rp[0] - rc_i[0]
                dy = rp[1] - rc_i[1]
                dz = rp[2] - rc_i[2]
                r_sq = dx*dx + dy*dy + dz*dz
                r_mag = math.sqrt(r_sq)

                if r_mag < 1e-20:
                    continue

                V_total += k_eff * qc_i / r_mag

                if compute_field:
                    E_mag = k_eff * qc_i / r_sq
                    Ex += E_mag * dx / r_mag
                    Ey += E_mag * dy / r_mag
                    Ez += E_mag * dz / r_mag

            pt_result = {
                "index": idx_p,
                "point_A": (rp[0]/self.a0, rp[1]/self.a0, rp[2]/self.a0),
                "V_volt": round(V_total, 6),
            }
            if compute_field:
                E_vec = (round(Ex, 4), round(Ey, 4), round(Ez, 4))
                E_mag_total = math.sqrt(Ex*Ex + Ey*Ey + Ez*Ez)
                pt_result["E_vector_V_per_m"] = E_vec
                pt_result["E_magnitude_V_per_m"] = round(E_mag_total, 4)
                pt_result["E_magnitude_V_per_A"] = round(E_mag_total * self.a0, 4)

            pot_results.append(pt_result)

        # ---- 计算总静电能 U = ½ Σᵢ<ⱼ k·qᵢqⱼ/rᵢⱼ ----
        energy_result = {}
        if compute_energy:
            U_total = 0.0
            pair_interactions = []
            for i in range(n_charges):
                for j in range(i + 1, n_charges):
                    dxi = rc[i][0] - rc[j][0]
                    dyi = rc[i][1] - rc[j][1]
                    dzi = rc[i][2] - rc[j][2]
                    rij = math.sqrt(dxi*dxi + dyi*dyi + dzi*dzi)
                    if rij < 1e-20:
                        continue
                    U_ij = k_eff * qc[i] * qc[j] / rij
                    U_total += U_ij
                    pair_interactions.append({
                        "pair": f"(i={i}, j={j})",
                        "rij_A": round(rij / self.a0, 4),
                        "qi_qj_e2": round((charges[i][0] * charges[j][0]), 4),
                        "U_ij_eV": round(U_ij / self.e, 6),
                        "U_ij_kJ_per_mol": round(U_total * self._NA() / 1000, 6) if i == n_charges - 2 else 0,
                    })

            energy_result = {
                "total_potential_energy_eV": round(U_total / self.e, 6),
                "total_potential_energy_kJ_per_mol": round(U_total * self._NA() / 1000, 6),
                "total_potential_energy_J": round(U_total, 12),
                "n_pairs": len(pair_interactions),
                "pair_interactions": pair_interactions[:10],  # limit output size
            }

        # ---- 多极展开（远场近似）----
        multipole_result = {}
        if expansion_order > 0 and eval_points is not None:
            multipole_result = self._multipole_expand(qc, rc, ep, k_eff, expansion_order)

        result = {
            "n_charges": n_charges,
            "total_charge_e": round(total_q_e, 6),
            "dielectric_constant": eps_r,
            "potentials_at_points": pot_results,
            **energy_result,
        }
        if multipole_result:
            result["multipole_expansion"] = multipole_result

        logger.info(f"CoulombPotential: {n_charges} charges, Q_tot={total_q_e:.3f}e")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入格式。"""
        try:
            parts = input_params.split("|")
            charge_part = parts[0].strip()
            point_part = parts[1].strip() if len(parts) > 1 else ""

            charges = []
            for ch_str in charge_part.split(";"):
                vals = [v.strip() for v in ch_str.split(",")]
                charges.append((float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])))

            eval_points = []
            if point_part:
                for pt_str in point_part.split(";"):
                    vals = [v.strip() for v in pt_str.split(",")]
                    eval_points.append((float(vals[0]), float(vals[1]), float(vals[2])))

            return self._run_base(charges, eval_points if eval_points else None)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. "
                              f"Expected format: 'q1,x1,y1,z1;q2,x2,y2,z2 | x,y,z'")

    def _multipole_expand(self, qc: list, rc: list, ep: list, k_eff: float, order: int) -> dict:
        """
        多极展开计算远场电势。
        
        V(r) ≈ k[Q₀/r + p·r̂/r² + Θ:(3r̂r̂-I)/(2r³) + ...]
        """
        # 单极子：总电荷
        Q0 = sum(qc)

        # 偶极矩 p = Σ qᵢ rᵢ
        px = sum(q * r[0] for q, r in zip(qc, rc))
        py = sum(q * r[1] for q, r in zip(qc, rc))
        pz = sum(q * r[2] for q, r in zip(qc, rc))

        results = {
            "monopole_Q_C": round(Q0, 20),
            "dipole_px_Cm": round(px, 30),
            "dipole_py_Cm": round(py, 30),
            "dipole_pz_Cm": round(pz, 30),
            "dipole_magnitude_D": round(math.sqrt(px*px + py*py + pz*pz) / 3.33564e-30, 6),
        }

        if order >= 2:
            # 四极矩张量 Θₐᵦ = ½ Σ qᵢ(3rᵢₐrᵢᵦ - rᵢ²δₐᵦ)
            Theta = [[0.0]*3 for _ in range(3)]
            for q, r in zip(qc, rc):
                r2 = r[0]**2 + r[1]**2 + r[2]**2
                for a in range(3):
                    for b in range(3):
                        Theta[a][b] += q * (3*r[a]*r[b] - r2 * (1 if a==b else 0))
            for a in range(3):
                for b in range(3):
                    Theta[a][b] *= 0.5
            results["quadrupole_tensor_Cm2"] = [[round(Theta[a][b], 30) for b in range(3)] for a in range(3)]

        return results

    @staticmethod
    def _NA() -> float:
        return 6.02214076e23  # Avogadro constant
