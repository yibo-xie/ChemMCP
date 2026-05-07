"""
多极展开工具 (MCP #475)。
计算电荷分布的多极矩（单极子、偶极、四极、八极）和远场相互作用能。
"""
import logging
import math
from typing import List, Tuple, Optional, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# ===== 物理常数 =====
K_E = 8.9875517923e9       # N·m²/C²
E_CHARGE = 1.602176634e-19   # C
ANGSTROM = 1e-10             # m
NA = 6.02214076e23          # mol⁻¹


@ChemMCPManager.register_tool
class MultipoleExpansion(BaseTool):
    """
    电荷分布的多极展开与长程相互作用。
    
    功能:
      - 计算多极矩: 单极 Q₀, 偶极 p, 四极 Θ, 八极 Ω
      - 远场电势展开: V(r) = k[Q₀/r + p·r̂/r² + Θ:(3r̂r̂-I)/(2r³) + ...]
      - 两个电荷分布之间的相互作用能
      - 多极-诱导多极相互作用 (诱导能)
    """
    __version__ = "0.1.0"
    name = "MultipoleExpansion"
    func_name = "multipole_expansion_analysis"
    description = "Compute multipole moments (monopole, dipole, quadrupole, octupole) of a charge distribution and long-range interaction potentials between two distributions."
    implementation_description = (
        "Expands the electrostatic potential of a charge distribution in powers of 1/r: "
        "V(r) = kQ₀/r + kp·r̂/r² + kΘ:(3r̂r̂-I)/(2r³) + kΩ·(tensor)/r⁴ + ... "
        "Computes interaction energy U = Q₁φ₀ + p₁·E₀ + (1/3)Θ₁:∇E₀ + ... "
        "All calculations use SI units internally."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Multipole Expansion", "Electrostatics", "Long-range Interaction", "Quantum Chemistry", "Moments"]
    required_envs = []

    code_input_sig = [
        ("charges1", "list", "N/A", "Distribution 1 point charges: [(q_e, x_A, y_A, z_A), ...]. q in e units, coords in Å."),
        ("charges2", "list", "None", "Distribution 2 for interaction calculation. Same format as charges1. If None, only compute multipoles of dist-1."),
        ("max_order", "int", "3", "Maximum multipole order: 0=monopole, 1=+dipole, 2=+quadrupole, 3=+octupole."),
        ("separation_vector_A", "list", "[0, 0, 10]", "Vector from dist-1 center to dist-2 center [x,y,z] in Å, for interaction energy."),
        ("include_interaction_energy", "bool", "True", "Whether to compute interaction energy between two distributions."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'charges1|charges2|order|separation'. Example: '1,0,0,0;-1,2,0,0|1,5,0,0;-1,7,0,0|3|0,0,5'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing all multipole moments, potential expansion coefficients, and interaction energies."),
    ]

    examples = [
        {
            "code_input": {
                "charges1": [(1.0, 0.0, 0.0, 0.0), (-1.0, 1.0, 0.0, 0.0)],
                "charges2": [(1.0, 3.0, 0.0, 0.0), (-1.0, 4.0, 0.0, 0.0)],
                "max_order": 3,
                "separation_vector_A": [2.0, 0.0, 0.0],
            },
            "text_input": {"input_params": "1,0,0,0;-1,1,0,0|1,3,0,0;-1,4,0,0|3|2,0,0"},
            "output": {"result": {"n_charges_dist1": 2, "total_charge_1_e": 0.0}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.k = K_E
        self.e = E_CHARGE
        self.a0 = ANGSTROM

    def _run_base(self, charges1: list, charges2: list = None,
                  max_order: int = 3, separation_vector_A: list = None,
                  include_interaction_energy: bool = True) -> dict:
        """
        核心计算逻辑。
        """
        if not charges1:
            raise ChemMCPInputError("charges1 cannot be empty.")

        # ---- 解析分布1 ----
        q1, r1 = self._parse_charges(charges1)
        # 质心作为原点
        cm1 = self._center_of_mass(q1, r1)
        r1_shifted = [(r[0]-cm1[0], r[1]-cm1[1], r[2]-cm1[2]) for r in r1]

        # ---- 计算分布1的多极矩 ----
        mult1 = self._compute_multipoles(q1, r1_shifted, max_order)

        result = {
            "distribution_1": {
                "n_charges": len(charges1),
                "center_of_mass_A": tuple(round(x, 6) for x in cm1),
                **mult1,
            }
        }

        # ---- 分布2（如果提供）----
        if charges2 is not None:
            q2, r2 = self._parse_charges(charges2)
            cm2 = self._center_of_mass(q2, r2)
            r2_shifted = [(r[0]-cm2[0], r[1]-cm2[1], r[2]-cm2[2]) for r in r2]
            mult2 = self._compute_multipoles(q2, r2_shifted, max_order)

            result["distribution_2"] = {
                "n_charges": len(charges2),
                "center_of_mass_A": tuple(round(x, 6) for x in cm2),
                **mult2,
            }

            # ---- 相互作用能 ----
            if include_interaction_energy and separation_vector_A is not None:
                R = [separation_vector_A[i] * self.a0 for i in range(3)]
                R_mag = math.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
                if R_mag < 1e-20:
                    raise ChemMCPInputError("Separation vector cannot be zero.")

                R_hat = [R[0]/R_mag, R[1]/R_mag, R[2]/R_mag]

                inter = self._interaction_energy(mult1, mult2, R, R_hat, R_mag, max_order)
                result["interaction"] = {
                    "separation_A": [round(x, 4) for x in separation_vector_A],
                    "separation_m": round(R_mag, 20),
                    **inter,
                }

        logger.info(f"MultipoleExpansion: dist1={len(charges1)} charges, order={max_order}")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.split("|")
            c1 = self._parse_charge_str(parts[0].strip())
            c2 = self._parse_charge_str(parts[1].strip()) if len(parts) > 1 and parts[1].strip() else None
            order = int(parts[2].strip()) if len(parts) > 2 else 3
            sep = [float(x) for x in parts[3].strip().split(",")] if len(parts) > 3 and parts[3].strip() else [0, 0, 10]
            return self._run_base(c1, c2, order, sep)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")

    @staticmethod
    def _parse_charges(charges):
        """解析电荷列表为SI单位。"""
        qc = []  # C
        rc = []  # m
        for ch in charges:
            q_e = float(ch[0])
            qc.append(q_e * E_CHARGE)
            rc.append((float(ch[1])*ANGSTROM, float(ch[2])*ANGSTROM, float(ch[3])*ANGSTROM))
        return qc, rc

    @staticmethod
    def _parse_charge_str(s):
        """从字符串解析电荷列表。"""
        charges = []
        for item in s.split(";"):
            vals = item.strip().split(",")
            charges.append((float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])))
        return charges

    @staticmethod
    def _center_of_mass(qc, rc):
        """计算质心（假设单位质量）。"""
        n = len(rc)
        return tuple(sum(r[i] for r in rc) / n for i in range(3))

    def _compute_multipoles(self, qc, rc, max_order):
        """
        计算各阶多极矩。
        
        返回包含各阶矩的字典，SI单位。
        """
        result = {}

        # ===== 0阶：单极子 Q₀ = Σ qᵢ =====
        Q0 = sum(qc)
        result["monopole_Q_C"] = round(Q0, 22)
        result["monopole_Q_e"] = round(Q0 / self.e, 6)

        if max_order < 1:
            return result

        # ===== 1阶：偶极矩 p = Σ qᵢ rᵢ =====
        px = sum(q * r[0] for q, r in zip(qc, rc))
        py = sum(q * r[1] for q, r in zip(qc, rc))
        pz = sum(q * r[2] for q, r in zip(qc, rc))
        p_mag = math.sqrt(px*px + py*py + pz*pz)
        result["dipole_px_Cm"] = round(px, 30)
        result["dipole_py_Cm"] = round(py, 30)
        result["dipole_pz_Cm"] = round(pz, 30)
        result["dipole_magnitude_D"] = round(p_mag / 3.33564e-30, 6)

        if max_order < 2:
            return result

        # ===== 2阶：四极矩张量 Θₐᵦ = ½Σ qᵢ(3rᵢₐrᵢᵦ - rᵢ²δₐᵦ) =====
        Theta = [[0.0]*3 for _ in range(3)]
        for q, r in zip(qc, rc):
            r2 = r[0]**2 + r[1]**2 + r[2]**2
            for a in range(3):
                for b in range(3):
                    Theta[a][b] += q * (3*r[a]*r[b] - r2*(1 if a==b else 0))
        for a in range(3):
            for b in range(3):
                Theta[a][b] *= 0.5
        result["quadrupole_tensor_Cm2"] = [[round(Theta[a][b], 30) for b in range(3)] for a in range(3)]
        # 四极矩大小（Frobenius范数）
        theta_frob = math.sqrt(sum(Theta[a][b]**2 for a in range(3) for b in range(3)))
        result["quadrupole_magnitude_Cm2"] = round(theta_frob, 30)

        if max_order < 3:
            return result

        # ===== 3阶：八极矩张量 Ωₐᵦ𝚌 =====
        Omega = [[[0.0]*3 for _ in range(3)] for _ in range(3)]
        for q, r in zip(qc, rc):
            rx, ry, rz = r[0], r[1], r[2]
            r2 = rx*rx + ry*ry + rz*rz
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        ra = [rx, ry, rz]
                        Omega[a][b][c] += q * (
                            15*ra[a]*ra[b]*ra[c]
                            - 3*r2*(ra[a]*(1 if b==c else 0) + ra[b]*(1 if a==c else 0) + ra[c]*(1 if a==b else 0))
                        )
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    Omega[a][b][c] /= 6  # 归一化因子
        omega_frob = math.sqrt(sum(Omega[a][b][c]**2 for a in range(3) for b in range(3) for c in range(3)))
        result["octupole_tensor_Cm3"] = [[[round(Omega[a][b][c], 35) for c in range(3)] for b in range(3)] for a in range(3)]
        result["octupole_magnitude_Cm3"] = round(omega_frob, 35)

        return result

    def _interaction_energy(self, m1, m2, R, R_hat, R_mag, order):
        """
        计算两个多极分布的相互作用能。
        
        U_int = k[Q₁Q₂/R + (Q₁p₂·R̂ - Q₂p₁·R̂)/R² + ...]
        """
        k = self.k
        U_total = 0.0
        terms = []

        # 零阶-零阶: Q₁Q₂/R
        Q1 = m1.get("monopole_Q_C", 0)
        Q2 = m2.get("monopole_Q_C", 0)
        if abs(Q1) > 1e-30 or abs(Q2) > 1e-30:
            U00 = k * Q1 * Q2 / R_mag
            U_total += U00
            terms.append({"order": "(0,0)", "description": "charge-charge", "U_eV": round(U00/self.e, 8)})

        if order < 1:
            pass
        else:
            # 偶极相关项
            p1 = [m1.get("dipole_px_Cm", 0), m1.get("dipole_py_Cm", 0), m1.get("dipole_pz_Cm", 0)]
            p2 = [m2.get("dipole_px_Cm", 0), m2.get("dipole_py_Cm", 0), m2.get("dipole_pz_Cm", 0)]

            # charge-dipole: Q₁(p₂·R̂)/R²
            p2_dot_R = sum(p2[i]*R_hat[i] for i in range(3))
            U_q1p2 = -k * Q1 * p2_dot_R / (R_mag**2)
            U_total += U_q1p2
            terms.append({"order": "(0,1)", "description": "Q1-p2·R̂/R²", "U_eV": round(U_q1p2/self.e, 10)})

            p1_dot_R = sum(p1[i]*R_hat[i] for i in range(3))
            U_q2p1 = k * Q2 * p1_dot_R / (R_mag**2)
            U_total += U_q2p1
            terms.append({"order": "(1,0)", "description": "Q2·p1·R̂/R²", "U_eV": round(U_q2p1/self.e, 10)})

            # dipole-dipole: [p₁·p₂ - 3(p₁·R̂)(p₂·R̂)] / R³
            p1p2 = sum(p1[i]*p2[i] for i in range(3))
            U_dd = k * (p1p2 - 3*p1_dot_R*p2_dot_R) / (R_mag**3)
            U_total += U_dd
            terms.append({"order": "(1,1)", "description": "dipole-dipole Keesom", "U_eV": round(U_dd/self.e, 12)})

        if order >= 2:
            # 粗略估算四极贡献（数量级）
            theta1 = m1.get("quadrupole_magnitude_Cm2", 0)
            theta2 = m2.get("quadrupole_magnitude_Cm2", 0)
            U_qq_est = k * theta1 * theta2 / (R_mag**5) * 0.1  # 几何因子 ~0.1
            U_total += U_qq_est
            terms.append({"order": "(2,2)", "description": "quadrupole-quadrupole (est.)", "U_eV": round(U_qq_est/self.e, 14)})

        return {
            "total_interaction_eV": round(U_total / self.e, 10),
            "total_interaction_kJ_per_mol": round(U_total * NA / 1000, 8),
            "terms": terms,
            "convergence_note": "Higher-order terms decay as 1/R^(l1+l2+1)",
        }
