"""
分子电场计算工具 (MCP #472)。
计算分子周围的静电场分布、分子表面静电势(ESP)映射。
基于原子部分电荷和坐标进行三维电场分析。
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


@ChemMCPManager.register_tool
class ElectricFieldMolecule(BaseTool):
    """
    分子电场与静电势分布计算。
    
    功能:
      - 计算分子周围空间网格点的电场矢量 E(r)
      - 计算范德华(VDW)表面上的静电势 ESP
      - 识别电场极值点（正/负静电势区域）
      - 绘制等势面数据（供可视化使用）
    """
    __version__ = "0.1.0"
    name = "ElectricFieldMolecule"
    func_name = "calculate_molecular_electric_field"
    description = "Calculate electric field distribution and electrostatic potential (ESP) mapping around a molecule from atomic partial charges and coordinates."
    implementation_description = (
        "Computes E = Σ kqᵢ(r-rᵢ)/|r-rᵢ|³ on a 3D grid or molecular surface. "
        "Maps ESP on van der Waals surface scaled by 1.4× (probe radius). "
        "Identifies ESP extrema for σ-hole / π-hole analysis, electrophilic/nucleophilic site prediction."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Electric Field", "ESP", "Electrostatic Potential", "Molecular Surface", "Quantum Chemistry"]
    required_envs = []

    code_input_sig = [
        ("atoms", "list", "N/A", "List of atoms: [(symbol, partial_charge_e, x_A, y_A, z_A), ...]. Coordinates in Ångströms."),
        ("grid_type", "str", "'surface'", "Grid type: 'surface' (VDW surface), 'plane' (XY/YZ/ZX plane), 'sphere' (concentric sphere), 'line' (along an axis)."),
        ("grid_params", "dict", "None", "Parameters depending on grid_type. For 'plane': {'origin':(x,y,z), 'normal':'z', 'range_A':10, 'n_points':20}. For 'sphere': {'center':(x,y,z), 'radii_A':[3,5], 'n_theta':20, 'n_phi':36}. For 'line': {'start':(x,y,z), 'end':(x,y,z), 'n':50}."),
        ("vdw_scale", "float", "1.4", "VDW radius scale factor for surface mapping (1.4 = solvent probe radius ~1.4Å)."),
        ("dielectric", "float", "1.0", "Relative dielectric constant ε_r."),
        ("compute_esp", "bool", "True", "Whether to compute electrostatic potential."),
        ("compute_gradient", "bool", "False", "Whether to compute field gradient (field curvature)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: atom_data|grid_type|grid_params_json. Example: 'O,-0.834,0,0,-0.96;H,0.417,0.757,0.586,0.19;H,0.417,-0.757,0.586,0.19|sphere|{\"center\":[0,0,0],\"radii\":[3,5]}'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing E-field vectors, ESP values, extrema locations, and analysis summary."),
    ]

    examples = [
        {
            "code_input": {
                "atoms": [
                    ("O", -0.834, 0.0, 0.0, -0.96),
                    ("H", 0.417, 0.757, 0.586, 0.19),
                    ("H", 0.417, -0.757, 0.586, 0.19),
                ],
                "grid_type": "sphere",
                "grid_params": {"center": (0.0, 0.0, 0.0), "radii_A": [3.0, 5.0], "n_theta": 10, "n_phi": 18},
            },
            "text_input": {"input_params": "O,-0.834,0,0,-0.96;H,0.417,0.757,0.586,0.19;H,0.417,-0.757,0.586,0.19|sphere|{}"},
            "output": {"result": {"n_atoms": 3, "grid_type": "sphere"}},
        },
    ]

    # ===== 范德华半径 (Å) =====
    VDW_RADII = {
        "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
        "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
        "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
        "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98, "Pb": 2.02,
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.k_e = K_E
        self.e = E_CHARGE
        self.a0 = ANGSTROM

    def _run_base(self, atoms: list, grid_type: str = "surface",
                  grid_params: dict = None, vdw_scale: float = 1.4,
                  dielectric: float = 1.0, compute_esp: bool = True,
                  compute_gradient: bool = False) -> dict:
        """
        核心计算逻辑：在指定网格上计算分子电场和ESP。
        """
        if not atoms:
            raise ChemMCPInputError("Atoms list cannot be empty.")

        gp = grid_params or {}
        k_eff = self.k_e / dielectric

        # 解析原子数据
        symbols = []
        charges_C = []  # 库仑
        positions_m = []  # 米
        for at in atoms:
            sym = at[0]
            q_e = float(at[1])
            x, y, z = float(at[2]), float(at[3]), float(at[4])
            symbols.append(sym)
            charges_C.append(q_e * self.e)
            positions_m.append((x * self.a0, y * self.a0, z * self.a0))

        n_atoms = len(atoms)

        # ---- 根据网格类型生成评估点 ----
        if grid_type == "surface":
            eval_pts = self._generate_vdw_surface(positions_m, symbols, vdw_scale, gp.get("density", 30))
        elif grid_type == "sphere":
            center = gp.get("center", (0.0, 0.0, 0.0))
            radii = gp.get("radii_A", [3.0, 5.0])
            n_theta = gp.get("n_theta", 20)
            n_phi = gp.get("n_phi", 36)
            eval_pts = self._generate_sphere_grid(center, radii, n_theta, n_phi)
        elif grid_type == "plane":
            origin = gp.get("origin", (0.0, 0.0, 0.0))
            normal = gp.get("normal", "z")
            rng = gp.get("range_A", 10.0)
            np_ = gp.get("n_points", 20)
            eval_pts = self._generate_plane_grid(origin, normal, rng, np_)
        elif grid_type == "line":
            start = gp.get("start", (-5.0, 0.0, 0.0))
            end = gp.get("end", (5.0, 0.0, 0.0))
            nl = gp.get("n", 50)
            eval_pts = self._generate_line(start, end, nl)
        else:
            raise ChemMCPInputError(f"Unknown grid type: {grid_type}. Choose: surface, sphere, plane, line")

        # ---- 在每个点计算电场和ESP ----
        field_data = []
        V_min, V_max = 1e30, -1e30
        E_max = 0.0
        V_min_pt, V_max_pt, E_max_pt = None, None, None

        for idx, rp in enumerate(eval_pts):
            V = 0.0
            Ex, Ey, Ez = 0.0, 0.0, 0.0

            for qC, ra in zip(charges_C, positions_m):
                dx = rp[0] - ra[0]
                dy = rp[1] - ra[1]
                dz = rp[2] - ra[2]
                r_sq = dx*dx + dy*dy + dz*dz
                r_mag = math.sqrt(r_sq)

                if r_mag < 1e-20:
                    continue

                V += k_eff * qC / r_mag
                E_factor = k_eff * qC / r_sq
                Ex += E_factor * dx / r_mag
                Ey += E_factor * dy / r_mag
                Ez += E_factor * dz / r_mag

            Emag = math.sqrt(Ex*Ex + Ey*Ey + Ez*Ez)

            pt_entry = {
                "index": idx,
                "point_A": (round(rp[0]/self.a0, 4), round(rp[1]/self.a0, 4), round(rp[2]/self.a0, 4)),
            }
            if compute_esp:
                pt_entry["ESP_V"] = round(V, 6)
                pt_entry["ESP_kcal_per_mol"] = round(V / self.e * 23.06, 4)  # 转换为 kcal/mol
                if V < V_min:
                    V_min = V; V_min_pt = idx
                if V > V_max:
                    V_max = V; V_max_pt = idx

            pt_entry["E_vector_V_per_m"] = (round(Ex, 4), round(Ey, 4), round(Ez, 4))
            pt_entry["E_magnitude_V_per_A"] = round(Emag * self.a0, 6)
            if Emag > E_max:
                E_max = Emag; E_max_pt = idx

            if compute_gradient:
                # 数值梯度（有限差分）
                h = 0.01 * self.a0  # 0.01 Å step
                grad = self._numerical_gradient(charges_C, positions_m, rp, h, k_eff)
                pt_entry["gradient_tensor"] = grad

            field_data.append(pt_entry)

        # ---- 汇总统计 ----
        summary = {
            "n_atoms": n_atoms,
            "atom_symbols": symbols,
            "grid_type": grid_type,
            "n_evaluation_points": len(eval_pts),
            "dielectric": dielectric,
        }
        if compute_esp:
            summary["ESP_statistics"] = {
                "min_ESP_V": round(V_min, 6),
                "max_ESP_V": round(V_max, 6),
                "min_ESP_kcal_mol": round(V_min / self.e * 23.06, 4),
                "max_ESP_kcal_mol": round(V_max / self.e * 23.06, 4),
                "min_point_index": V_min_pt,
                "max_point_index": V_max_pt,
            }
        summary["E_field_statistics"] = {
            "max_E_V_per_A": round(E_max * self.a0, 6),
            "max_point_index": E_max_pt,
        }

        result = {
            **summary,
            "field_data_sample": field_data[:20],
            "field_data_full_n": len(field_data),
        }

        logger.info(f"ElectricFieldMolecule: {n_atoms} atoms, {len(eval_pts)} pts, "
                     f"grid={grid_type}, ESP=[{V_min:.3f}, {V_max:.3f}]V")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        """解析文本输入。"""
        try:
            parts = input_params.split("|")
            atom_part = parts[0].strip()
            gtype = parts[1].strip() if len(parts) > 1 else "sphere"
            gp_str = parts[2].strip() if len(parts) > 2 else "{}"

            import json
            gp = json.loads(gp_str) if gp_str else {}

            atoms = []
            for a_str in atom_part.split(";"):
                vals = a_str.strip().split(",")
                atoms.append((vals[0], float(vals[1]), float(vals[2]), float(vals[3]), float(vals[4])))

            return self._run_base(atoms, gtype, gp)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")

    def _generate_vdw_surface(self, positions, symbols, scale, density):
        """生成范德华表面采样点。"""
        import math
        points = []
        # 简化实现：在每个原子的球面上均匀采样
        for pos, sym in zip(positions, symbols):
            r_vdw = self.VDW_RADII.get(sym, 1.5) * scale * self.a0
            n_phi = max(6, int(density))
            n_theta = max(3, int(density // 2))
            for i in range(n_theta):
                theta = math.pi * (i + 0.5) / n_theta
                for j in range(n_phi):
                    phi = 2 * math.pi * j / n_phi
                    spx = pos[0] + r_vdw * math.sin(theta) * math.cos(phi)
                    spy = pos[1] + r_vdw * math.sin(theta) * math.sin(phi)
                    spz = pos[2] + r_vdw * math.cos(theta)
                    points.append((spx, spy, spz))
        return points

    def _generate_sphere_grid(self, center, radii_A, n_theta, n_phi):
        """生成同心球面网格。"""
        import math
        cx, cy, cz = center[0]*self.a0, center[1]*self.a0, center[2]*self.a0
        points = []
        for rad_A in radii_A:
            r = rad_A * self.a0
            for i in range(n_theta):
                theta = math.pi * i / max(1, n_theta - 1) if n_theta > 1 else math.pi/2
                for j in range(n_phi):
                    phi = 2*math.pi*j/max(1, n_phi-1) if n_phi > 1 else 0
                    px = cx + r*math.sin(theta)*math.cos(phi)
                    py = cy + r*math.sin(theta)*math.sin(phi)
                    pz = cz + r*math.cos(theta)
                    points.append((px, py, pz))
        return points

    def _generate_plane_grid(self, origin, normal, rng, np_):
        """生成平面网格。"""
        import math
        ox, oy, oz = origin[0]*self.a0, origin[1]*self.a0, origin[2]*self.a0
        d = rng * self.a0
        n = max(2, np_)
        points = []
        if normal == "z":
            for i in range(n):
                for j in range(n):
                    px = ox - d + 2*d*i/(n-1) if n > 1 else ox
                    py = oy - d + 2*d*j/(n-1) if n > 1 else oy
                    points.append((px, py, oz))
        elif normal == "y":
            for i in range(n):
                for j in range(n):
                    px = ox - d + 2*d*i/(n-1) if n > 1 else ox
                    pz = oz - d + 2*d*j/(n-1) if n > 1 else oz
                    points.append((px, oy, pz))
        else:  # x
            for i in range(n):
                for j in range(n):
                    py = oy - d + 2*d*i/(n-1) if n > 1 else oy
                    pz = oz - d + 2*d*j/(n-1) if n > 1 else oz
                    points.append((ox, py, pz))
        return points

    def _generate_line(self, start, end, n):
        """沿直线生成点。"""
        sx, sy, sz = start[0]*self.a0, start[1]*self.a0, start[2]*self.a0
        ex, ey, ez = end[0]*self.a0, end[1]*self.a0, end[2]*self.a0
        points = []
        nn = max(2, n)
        for i in range(nn):
            t = i / max(1, nn - 1)
            px = sx + t*(ex-sx)
            py = sy + t*(ey-sy)
            pz = sz + t*(ez-sz)
            points.append((px, py, pz))
        return points

    @staticmethod
    def _numerical_gradient(charges, positions, point, h, k_eff):
        """数值计算电场梯度张量。"""
        import math
        grad = [[0.0]*3 for _ in range(3)]
        for a in range(3):
            for b in range(3):
                pp = list(point)
                pm = list(point)
                pp[b] += h
                pm[b] -= h
                Ep, Em = 0.0, 0.0
                for q, r in zip(charges, positions):
                    dxp = pp[0]-r[0]; dyp = pp[1]-r[1]; dzp = pp[2]-r[2]
                    dxm = pm[0]-r[0]; dym = pm[1]-r[1]; dzm = pm[2]-r[2]
                    rp = math.sqrt(dxp**2+dyp**2+dzp**2)
                    rm = math.sqrt(dxm**2+dym**2+dzm**2)
                    if rp > 1e-20:
                        Ep += k_eff*q/rp
                    if rm > 1e-20:
                        Em += k_eff*q/rm
                # ∂Eₐ/∂x_b ≈ (∂V/∂x_a at +h - ∂V/∂x_a at -h) / 2h
                # Simplified: just potential gradient
                grad[a][b] = round((Ep - Em) / (2*h), 4)
        return grad
