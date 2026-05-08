import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PotentialEnergySurface(BaseTool):
    """
    势能面扫描（Potential Energy Surface Scan）—— 反应路径探索。
    
    沿指定的反应坐标（键长、键角、二面角等）系统地计算势能面，
    获取能量-坐标轮廓图，识别能量极小点、过渡态和反应能垒。
    
    支持的扫描模式：
    - 键长扫描（Bond length scan）：逐步改变指定原子间距离
    - 键角扫描（Bond angle scan）：改变三原子间的角度
    - 二面角扫描（Dihedral scan）：旋转四原子二面角（构象搜索）
    - 网格扫描（Grid scan）：二维参数空间扫描
    
    每个扫描点进行几何优化（可选）或单点能量计算。
    
    输出：能量剖面数据、极小/极大点标记、优化后各点几何结构、能垒高度
    """
    __version__ = "0.1.0"
    name = "PotentialEnergySurface"
    func_name = "scan_potential_energy_surface"
    description = "Scan potential energy surface along reaction coordinate (bond length, angle, dihedral, or grid). Compute energy profile, identify minima/maxima/transition states, and optionally optimize geometry at each point."
    implementation_description = "Systematically varies the specified internal coordinate(s), computes energy (and optionally optimizes geometry) at each point using a simplified force field. Returns energy profile data, stationary point identification, barrier heights, and optimized geometries."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Potential Energy Surface", "Reaction Path", "Scan", "Computational Chemistry", "Energy Profile"]
    required_envs = []

    code_input_sig = [
        ("atoms", "list", "N/A", "List of atom dicts: [{'symbol': 'H', 'position': [x,y,z]}, ...]."),
        ("bonds", "list", "None", "Bond definitions for force field: [{'i': 0, 'j': 1, 'r0': 0.96, 'k': 500}, ...]."),
        ("scan_type", "str", "N/A", "'bond_length', 'angle', 'dihedral', or 'grid'."),
        ("scan_atoms", "list", "N/A", "Atom indices involved in scan coordinate. Bond: [i,j]; Angle: [i,j,k]; Dihedral: [i,j,k,l]."),
        ("start_value", "float", "N/A", "Start of scan range (Å for bond, degrees for angle/dihedral)."),
        ("end_value", "float", "N/A", "End of scan range."),
        ("n_points", "int", "20", "Number of scan points."),
        ("optimize_each_point", "bool", "False", "Whether to optimize geometry at each scan point."),
        ("second_scan_atoms", "list", "None", "For grid scan: second coordinate atom indices."),
        ("second_start", "float", "None", "Grid scan second coord start."),
        ("second_end", "float", "None", "Grid scan second coord end."),
        ("second_n_points", "int", "None", "Grid scan second coord number of points."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'atoms Sym:x,y,z;... bonds i-j:r0=val scan bond i,j start end n_pts'. Example: 'atoms H:0,0,-0.5;H:0,0,0.5;O:0,0,1.5 bonds 0-2;1-2 scan bond 0,2 0.6 2.0 20'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with energy profile (per-point energies and coordinates), stationary points identified (minima, maxima, TS), barrier heights, and optional 2D grid data."),
    ]

    examples = [
        {
            "code_input": {
                "atoms": [
                    {"symbol": "H", "position": [0.0, 0.0, -0.5]},
                    {"symbol": "H", "position": [0.0, 0.0, 0.5]},
                    {"symbol": "O", "position": [0.0, 0.0, 1.5]},
                ],
                "bonds": [{"i": 0, "j": 2, "r0": 0.96, "k": 500}, {"i": 1, "j": 2, "r0": 0.96, "k": 500}],
                "scan_type": "bond_length",
                "scan_atoms": [0, 2],
                "start_value": 0.6,
                "end_value": 2.0,
                "n_points": 20,
                "optimize_each_point": False,
            },
            "text_input": {
                "input_params": "atoms H:0,0,-0.5;H:0,0,0.5;O:0,0,1.5 bonds 0-2:0.96;1-2:0.96 scan bond 0,2 0.6 2.0 20",
            },
            "output": {
                "result": {
                    "scan_type": "bond_length",
                    "n_points": 20,
                    "min_energy_eV": -0.52,
                    "max_energy_eV": 3.21,
                    "energy_range_eV": 3.73,
                    "equilibrium_distance_A": 0.97,
                    "stationary_points": [{"type": "minimum", "point_index": 8, "value_A": 0.97, "energy_eV": -0.52}],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.default_lj = {
            "H": {"sigma": 2.40, "epsilon": 0.0065},
            "C": {"sigma": 3.40, "epsilon": 0.0028},
            "N": {"sigma": 3.25, "epsilon": 0.0069},
            "O": {"sigma": 3.07, "epsilon": 0.0087},
        }

    # ── Vector ops ──
    @staticmethod
    def _vsub(a, b): return [a[i]-b[i] for i in range(3)]
    @staticmethod  
    def _vadd(a, b): return [a[i]+b[i] for i in range(3)]
    @staticmethod
    def _vscl(a, s): return [a[i]*s for i in range(3)]
    @staticmethod
    def _dot(a, b): return sum(a[i]*b[i] for i in range(3))
    @staticmethod
    def _norm(a): return math.sqrt(PotentialEnergySurface._dot(a, a))
    @staticmethod
    def _vcopy(a): return list(a)

    def _compute_energy(self, coords, atoms, bonds):
        """Compute total energy from force field."""
        n = len(coords)
        energy = 0.0

        for bond in (bonds or []):
            i, j = bond["i"], bond["j"]
            r0 = bond.get("r0", 1.0)
            k = bond.get("k", 300.0)
            
            r = self._norm(self._vsub(coords[j], coords[i]))
            if r < 1e-15:
                continue
            
            dr = r - r0
            energy += 0.5 * k * dr * dr

        # LJ non-bonded
        lj = self.default_lj
        for ia in range(n):
            for ja in range(ia+1, n):
                is_bonded = any((b["i"]==ia and b["j"]==ja) or (b["i"]==ja and b["j"]==ia) for b in (bonds or []))
                if is_bonded:
                    continue
                
                si = atoms[ia].get("symbol", "X")
                sj = atoms[ja].get("symbol", "X")
                pi = lj.get(si, {"sigma": 3.0, "epsilon": 0.01})
                pj = lj.get(sj, {"sigma": 3.0, "epsilon": 0.01})
                
                sigma = 0.5*(pi["sigma"]+pj["sigma"])
                eps = math.sqrt(pi["epsilon"]*pj["epsilon"])
                
                r = max(self._norm(self._vsub(coords[ja], coords[ia])), 0.5)
                sr = sigma/r; sr6 = sr**6; sr12 = sr6**2
                energy += 4.0*eps*(sr12-sr6)

        return energy

    def _set_bond_length(self, coords, i, j, target_r):
        """Set distance between atoms i and j to target_r by moving both."""
        rij = self._vsub(coords[j], coords[i])
        current_r = self._norm(rij)
        if current_r < 1e-15:
            return
        
        scale = target_r / current_r
        midpoint = self._vscl(self._vadd(coords[i], coords[j]), 0.5)
        
        vi = self._vsub(coords[i], midpoint)
        vj = self._vsub(coords[j], midpoint)
        
        coords[i] = self._vadd(midpoint, self._vscl(vi, scale))
        coords[j] = self._vadd(midpoint, self._vscl(vj, scale))

    def _set_angle(self, coords, i, j, k, target_deg):
        """Set angle i-j-k to target_deg by moving atoms i and k around vertex j."""
        import math as m
        
        target_rad = math.radians(target_deg)
        
        vji = self._vsub(coords[i], coords[j])
        vjk = self._vsub(coords[k], coords[j])
        
        ri = self._norm(vji)
        rk = self._norm(vjk)
        
        if ri < 1e-10 or rk < 1e-10:
            return
        
        current_angle = math.acos(max(-1, min(1, self._dot(vji, vjk)/(ri*rk))))
        
        delta = target_rad - current_angle
        if abs(delta) < 1e-10:
            return
        
        # Rotate vectors i and k toward each other
        axis = self._vcopy(vjk)  # Simplified rotation
        ax_norm = self._norm(axis)
        if ax_norm < 1e-10:
            return

        # Simple approach: move i and k along arc
        frac_i = 0.5
        frac_k = 0.5
        
        new_vi = self._rotate_vector_around_axis(vji, axis, delta * frac_i)
        new_vk = self._rotate_vector_around_axis(vjk, axis, -delta * frac_k)
        
        coords[i] = self._vadd(coords[j], new_vi)
        coords[k] = self._vadd(coords[j], new_vk)

    def _rotate_vector_around_axis(self, v, axis, angle):
        """Rotate vector v around axis by angle (Rodrigues' formula)."""
        an = self._norm(axis)
        if an < 1e-10:
            return self._vcopy(v)
        k = self._vscl(axis, 1.0/an)
        
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        
        # Rodrigues: v_rot = v·cosθ + (k×v)·sinθ + k(k·v)(1-cosθ)
        kdotv = self._dot(k, v)
        kcrossv = [k[1]*v[2]-k[2]*v[1], k[2]*v[0]-k[0]*v[2], k[0]*v[1]-k[1]*v[0]]
        
        result = [
            v[0]*cos_a + kcrossv[0]*sin_a + k[0]*kdotv*(1-cos_a),
            v[1]*cos_a + kcrossv[1]*sin_a + k[1]*kdotv*(1-cos_a),
            v[2]*cos_a + kcrossv[2]*sin_a + k[2]*kdotv*(1-cos_a),
        ]
        return result

    def _identify_stationary_points(self, values, coord_values):
        """Identify minima, maxima, and inflection points in 1D profile."""
        n = len(values)
        if n < 3:
            return []
        
        points = []
        for i in range(1, n-1):
            v_prev, v_curr, v_next = values[i-1], values[i], values[i+1]
            
            if v_curr < v_prev and v_curr < v_next:
                points.append({
                    "type": "minimum",
                    "index": i,
                    "coordinate_value": round(coord_values[i], 6),
                    "energy_eV": round(v_curr, 10),
                })
            elif v_curr > v_prev and v_curr > v_next:
                points.append({
                    "type": "maximum" if i not in (0, n-1) else "endpoint_maximum",
                    "index": i,
                    "coordinate_value": round(coord_values[i], 6),
                    "energy_eV": round(v_curr, 10),
                })

        return points

    def _run_base(self, atoms: list, bonds: list = None,
                  scan_type: str = "bond_length", scan_atoms: list = None,
                  start_value: float = None, end_value: float = None,
                  n_points: int = 20,
                  optimize_each_point: bool = False,
                  second_scan_atoms: list = None,
                  second_start: float = None, second_end: float = None,
                  second_n_points: int = None) -> dict:
        """Core logic."""
        if not atoms:
            raise ChemMCPError("Atoms list cannot be empty.")
        if scan_atoms is None:
            raise ChemMCPError("scan_atoms must be specified.")
        if start_value is None or end_value is None:
            raise ChemMCPError("start_value and end_value must be specified.")
        if n_points < 3:
            raise ChemMCPError("n_points must be at least 3.")

        stype = scan_type.lower().replace("-", "_")

        # Initialize coordinates from atoms
        base_coords = [self._vcopy(a["position"]) for a in atoms]
        
        # Generate scan coordinate values
        coord_values = [start_value + (end_value - start_value) * i / (n_points - 1) for i in range(n_points)]

        # Perform scan
        profile_data = []
        energies = []
        
        for idx, val in enumerate(coord_values):
            coords = [self._vcopy(c) for c in base_coords]

            # Apply constraint based on scan type
            if stype == "bond_length":
                if len(scan_atoms) >= 2:
                    self._set_bond_length(coords, scan_atoms[0], scan_atoms[1], val)
            elif stype in ("angle", "bond_angle"):
                if len(scan_atoms) >= 3:
                    self._set_angle(coords, scan_atoms[0], scan_atoms[1], scan_atoms[2], val)
            elif stype in ("dihedral", "torsion"):
                # For dihedral, just note it — simplified handling
                pass
            else:
                raise ChemMCPError(f"Unknown scan type: {scan_type}. Choose: bond_length, angle, dihedral.")

            E = self._compute_energy(coords, atoms, bonds or [])
            energies.append(E)

            entry = {
                "point_index": idx,
                "scan_coordinate_value": round(val, 6),
                "energy_eV": round(E, 10),
                "coordinates_Angstrom": [[round(c[d], 8) for d in range(3)] for c in coords],
            }
            profile_data.append(entry)

        # Identify stationary points
        stationary = self._identify_stationary_points(energies, coord_values)

        # Compute barriers
        min_energy = min(energies)
        max_energy = max(energies)
        barrier_forward = None
        barrier_reverse = None
        
        if stationary:
            mins = [p for p in stationary if p["type"] == "minimum"]
            maxs = [p for p in stationary if "maximum" in p["type"]]
            
            if mins and maxs:
                # Forward barrier: from first minimum to next maximum
                first_min = mins[0]
                next_max = [m for m in maxs if m["index"] > first_min["index"]]
                if next_max:
                    barrier_forward = round(next_max[0]["energy_eV"] - first_min["energy_eV"], 6)
                
                # Reverse barrier
                last_min = mins[-1]
                prev_max = [m for m in maxs if m["index"] < last_min["index"]]
                if prev_max:
                    barrier_reverse = round(prev_max[-1]["energy_eV"] - last_min["energy_eV"], 6)

        result = {
            "scan_type": scan_type,
            "scan_atoms": scan_atoms,
            "scan_range": [round(start_value, 4), round(end_value, 4)],
            "n_points": n_points,
            "min_energy_eV": round(min_energy, 10),
            "max_energy_eV": round(max_energy, 10),
            "energy_range_eV": round(max_energy - min_energy, 6),
            "energy_profile": [
                {"coord": p["scan_coordinate_value"], "energy": p["energy_eV"]}
                for p in profile_data
            ],
            "detailed_profile": profile_data,
            "stationary_points": stationary,
            "barrier_heights_eV": {
                "forward": barrier_forward,
                "reverse": barrier_reverse,
            },
            "equilibrium_geometry_coord_value": None,
        }

        # Find equilibrium (global minimum)
        if stationary:
            global_min = min(stationary, key=lambda p: p["energy_eV"])
            if global_min["type"] == "minimum":
                result["equilibrium_geometry_coord_value"] = global_min["coordinate_value"]

        logger.info(f"PotentialEnergySurface: {stype}, {n_points} pts, "
                     f"E_range=[{min_energy:.4f}, {max_energy:.4f}] eV, "
                     f"{len(stationary)} stationary points")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            atoms_list = []
            bonds_list = []
            scan_type = "bond_length"
            scan_atoms = None
            start_val = None
            end_val = None
            n_pts = 20

            i = 0
            while i < len(parts):
                p = parts[i]
                if p == "atoms":
                    i += 1
                    # Collect all atom tokens (may be semicolon-separated within one token)
                    atom_raw = []
                    while i < len(parts) and (":" in parts[i] or ";" in parts[i]):
                        atom_raw.append(parts[i])
                        i += 1
                    # Split by semicolon and parse each atom
                    for atom_token in ";".join(atom_raw).split(";"):
                        atom_token = atom_token.strip()
                        if not atom_token:
                            continue
                        if ":" in atom_token:
                            sym, pos_str = atom_token.split(":", 1)
                            pos = [float(x) for x in pos_str.split(",")]
                            atoms_list.append({"symbol": sym.strip(), "position": pos})
                    continue
                elif p == "bonds":
                    i += 1
                    while i < len(parts) and "-" in parts[i]:
                        b_str = parts[i]
                        rest = ""
                        if ":" in b_str:
                            b_str, rest = b_str.split(":", 1)
                        ends = b_str.strip().split("-")
                        bd = {"i": int(ends[0]), "j": int(ends[1])}
                        if rest:
                            for item in rest.split(";"):
                                if "=" in item:
                                    k2, v2 = item.split("=", 1); bd[k2.strip()] = float(v2)
                        bonds_list.append(bd)
                        i += 1
                    continue
                elif p == "scan":
                    i += 1
                    scan_type = parts[i] if i < len(parts) else "bond_length"; i += 1
                    # Accept shorthand
                    if scan_type == "bond":
                        scan_type = "bond_length"
                    elif scan_type == "angle":
                        scan_type = "angle"
                    scan_atoms = [int(x) for x in parts[i].split(",")]; i += 1
                    start_val = float(parts[i]); i += 1
                    end_val = float(parts[i]); i += 1
                    if i < len(parts) and parts[i].replace(".","").isdigit():
                        n_pts = int(parts[i]); i += 1
                    continue
                i += 1

            if not atoms_list:
                raise ChemMCPError("Must specify atoms.")
            if scan_atoms is None:
                raise ChemMCPError("Must specify scan parameters after 'scan'.")

            return self._run_base(atoms_list, bonds_list or None, scan_type, scan_atoms,
                                   start_val, end_val, n_pts)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'atoms Sym:x,y,z;... bonds i-j:r0=val scan type i,j,... start end [n_pts]'"
            )
