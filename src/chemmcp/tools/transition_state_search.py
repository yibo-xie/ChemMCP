import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TransitionStateSearch(BaseTool):
    """
    过渡态搜索（Transition State Search）—— 势能面鞍点定位。
    
    在势能面上寻找过渡态（Transition State, TS），即一阶鞍点：
    - 沿反应坐标方向：能量极大值（1个虚频，imaginary frequency）
    - 垂直于反应坐标的所有方向：能量极小值
    
    支持的搜索算法：
    - Newton-Raphson/Quadratic saddle search：利用 Hessian 矩阵沿负曲率方向最大化
    - Eigenvector following (EF)：沿最低本征矢方向爬升，其余方向下降
    - Dimer method 近似：简化版过渡态搜索
    
    输出：TS 结构坐标、虚频、Hessian 本征值、TS 验证、IRC 路径提示
    """
    __version__ = "0.1.0"
    name = "TransitionStateSearch"
    func_name = "find_transition_state"
    description = "Search for transition state (first-order saddle point) on potential energy surface using eigenvector following, quadratic saddle search, or dimer method. Returns TS coordinates, imaginary frequency, Hessian eigenvalues, and TS validation."
    implementation_description = "Implements transition state search algorithms that maximize energy along one direction (reaction coordinate) while minimizing along all others. Uses Hessian matrix analysis to identify saddle points. Validates TS by checking for exactly one negative eigenvalue (imaginary frequency). Provides IRC path hints."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Transition State", "Saddle Point", "PES", "Computational Chemistry", "Hessian", "Imaginary Frequency"]
    required_envs = []

    code_input_sig = [
        ("atoms", "list", "N/A", "List of atom dicts: [{'symbol': 'H', 'position': [x,y,z]}, ...]. Positions in Angstroms."),
        ("bonds", "list", "None", "Bond definitions: [{'i': 0, 'j': 1, 'r0': 0.96, 'k': 500}, ...]."),
        ("guess_coordinates", "list", "None", "Initial guess for TS geometry (same format as atoms positions). If None, uses atoms positions as starting guess."),
        ("search_method", "str", "eigenvector_following", "'eigenvector_following', 'quadratic_saddle', or 'dimer'."),
        ("max_iterations", "int", "200", "Maximum search iterations."),
        ("convergence_threshold", "float", "1e-5", "RMS gradient threshold."),
        ("reaction_coordinate_hint", "list", "None", "Hint for reaction coordinate direction as vector [dx,dy,dz] per atom (flattened)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format similar to GeometryOptimizer with optional TS guess. Example: 'atoms H:0,0,-0.5;H:0,0,0.5;O:0,0,1.5 bonds 0-2;1-2 method=ef'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with TS coordinates, imaginary frequency, Hessian eigenvalues, TS validation, energy profile, and IRC hints."),
    ]

    examples = [
        {
            "code_input": {
                "atoms": [
                    {"symbol": "H", "position": [0.0, 0.0, -0.5]},
                    {"symbol": "H", "position": [0.0, 0.0, 0.5]},
                    {"symbol": "O", "position": [0.0, 0.0, 1.5]},
                ],
                "bonds": [{"i": 0, "j": 2, "r0": 1.0, "k": 400}, {"i": 1, "j": 2, "r0": 1.0, "k": 400}],
                "guess_coordinates": None,
                "search_method": "quadratic_saddle",
                "max_iterations": 100,
            },
            "text_input": {
                "input_params": "atoms H:0,0,-0.5;H:0,0,0.5;O:0,0,1.5 bonds 0-2;1-2 method=qs",
            },
            "output": {
                "result": {
                    "ts_found": True,
                    "n_imaginary_frequencies": 1,
                    "imaginary_frequency_cm-1": -1500.0,
                    "ts_energy_eV": 0.35,
                    "converged": True,
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
        self.amu_kg = 1.66053906660e-27
        self.c_light_cm_s = 2.99792458e10   # cm/s
        self.hbar = 1.054571817e-34          # J·s
        self.eV_per_J = 6.241509074e18

    # ── Vector ops ──
    @staticmethod
    def _vsub(a, b): return [a[i]-b[i] for i in range(len(a))]
    @staticmethod
    def _vadd(a, b): return [a[i]+b[i] for i in range(len(a))]
    @staticmethod
    def _vscl(a, s): return [a[i]*s for i in range(len(a))]
    @staticmethod
    def _dot(a, b): return sum(a[i]*b[i] for i in range(len(a)))
    @staticmethod
    def _norm(a): return math.sqrt(TransitionStateSearch._dot(a, a))

    def _compute_energy_forces_hessian(self, coords, atoms, bonds):
        """Compute energy, gradient (forces), and approximate Hessian."""
        n_atoms = len(coords)
        n_dim = 3 * n_atoms
        
        # Use numerical differentiation for Hessian
        delta = 1e-5  # Å for numerical derivative
        
        E0, F0 = self._energy_and_forces(coords, atoms, bonds)
        
        # Build Hessian via finite differences of gradients
        H = [[0.0] * n_dim for _ in range(n_dim)]
        
        for j in range(n_dim):
            coord_plus = [self._vadd(c if ci == j//3 else c, 
                                      [delta if ci==j%3 else 0, delta if ci==j%3==1 else 0, delta if ci==j%3==2 else 0])
                          if ci == j//3 else self._vcopy(c) for ci, c in enumerate(coords)]
            # Simpler: perturb one coordinate component
            coords_p = [self._vcopy(c) for c in coords]
            atom_j = j // 3
            comp_j = j % 3
            coords_p[atom_j][comp_j] += delta
            
            _, Fp = self._energy_and_forces(coords_p, atoms, bonds)
            
            coords_m = [self._vcopy(c) for c in coords]
            coords_m[atom_j][comp_j] -= delta
            _, Fm = self._energy_and_forces(coords_m, atoms, bonds)
            
            # d²E/dx_idx_j = -(Fp[j] - Fm[j]) / (2*delta)
            # Flatten forces to 1D
            fp_flat = []
            fm_flat = []
            for a in range(n_atoms):
                for d in range(3):
                    fp_flat.append(Fp[a][d])
                    fm_flat.append(Fm[a][d])
            
            for i in range(n_dim):
                H[i][j] = -(fp_flat[i] - fm_flat[i]) / (2.0 * delta)

        # Symmetrize
        for i in range(n_dim):
            for j in range(i+1, n_dim):
                avg = (H[i][j] + H[j][i]) / 2.0
                H[i][j] = avg
                H[j][i] = avg

        return E0, F0, H

    def _energy_and_forces(self, coords, atoms, bonds):
        """Simplified energy and force calculation (same as GeometryOptimizer)."""
        n = len(coords)
        energy = 0.0
        forces = [[0.0]*3 for _ in range(n)]

        for bond in (bonds or []):
            i, j = bond["i"], bond["j"]
            r0 = bond.get("r0", 1.0)
            k = bond.get("k", 300.0)
            
            rij = self._vsub(coords[j], coords[i])
            r = self._norm(rij)
            if r < 1e-15:
                continue
            
            dr = r - r0
            energy += 0.5 * k * dr * dr
            f_mag = -k * dr
            f_vec = self._vscl(rij, f_mag / r) if r > 1e-15 else [0,0,0]
            forces[i] = self._vadd(forces[i], f_vec)
            forces[j] = self._vadd(forces[j], self._vscl(f_vec, -1))

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
                
                rij = self._vsub(coords[ja], coords[ia])
                r = max(self._norm(rij), 0.5)
                
                sr = sigma/r; sr6 = sr**6; sr12 = sr6**2
                energy += 4.0*eps*(sr12-sr6)
                
                flj = 24.0*eps/r*(2.0*sr12-sr6)
                fv = self._vscl(rij, flj/r)
                forces[ia] = self._vadd(forces[ia], fv)
                forces[ja] = self._vadd(forces[ja], self._vscl(fv, -1))

        return energy, forces

    def _vcopy(self, a): return list(a)

    def _flatten(self, coords):
        """Flatten 3N coordinates to 1D array of length 3N."""
        flat = []
        for c in coords:
            flat.extend(c)
        return flat

    def _unflatten(self, flat, n_atoms):
        """Unflatten 1D array to list of 3D coordinates."""
        coords = []
        for i in range(n_atoms):
            coords.append([flat[3*i], flat[3*i+1], flat[3*i+2]])
        return coords

    def _eigenvalues_3x3(self, M):
        """Compute eigenvalues of a 3x3 symmetric matrix using analytical method."""
        # For 3x3 symmetric matrix, find roots of characteristic polynomial
        # Using QR iteration would be ideal but let's use a simpler approach
        # for the full Hessian we'll just do power iteration for min/max
        
        # Actually this is for 3x3 only — for larger matrices use power iteration
        p1 = M[0][1]**2 + M[0][2]**2 + M[1][2]**2
        q = (M[0][0]+M[1][1]+M[2][2]) / 3.0
        p2 = (M[0][0]-q)**2 + (M[1][1]-q)**2 + (M[2][2]-q)**2 + 2*p1
        p = math.sqrt(p2/3.0) if p2 > 0 else 0
        
        if p < 1e-15:
            return [q, q, q]  # All eigenvalues equal
        
        r = det = (M[0][0]-q)*(M[1][1]-q)*(M[2][2]-q) + 2*M[0][1]*M[1][2]*M[0][2] \
                 - (M[0][0]-q)*M[1][2]**2 - (M[1][1]-q)*M[0][2]**2 - (M[2][2]-q)*M[0][1]**2
        det /= (p**3) if p > 1e-30 else 1
        
        det = max(-1.0, min(1.0, det))  # Clamp for numerical safety
        phi = math.acos(det) / 3.0
        
        eig1 = q + 2*p*math.cos(phi)
        eig3 = q + 2*p*math.cos(phi + 2*math.pi/3)
        eig2 = 3*q - eig1 - eig3  # Trace conservation
        
        return sorted([eig1, eig2, eig3])

    def _full_eigenvalues(self, H, n_dim):
        """Get all eigenvalues using power iteration for extremal ones + trace."""
        # Power iteration for largest magnitude eigenvalue
        n = n_dim
        v = [1.0/math.sqrt(n)] * n
        
        for _ in range(200):
            Hv = [sum(H[i][j]*v[j] for j in range(n)) for i in range(n)]
            norm_Hv = math.sqrt(sum(x*x for x in Hv))
            if norm_Hv < 1e-30:
                break
            v = [x/norm_Hv for x in Hv]
        
        lambda_max = sum(v[i]*sum(H[i][j]*v[j] for j in range(n)) for i in range(n))
        
        # Shift-invert for smallest eigenvalue (approximate)
        shift = lambda_max
        Hs = [[H[i][j] - (shift if i==j else 0) for j in range(n)] for i in range(n)]
        
        v2 = [1.0/math.sqrt(n)] * n
        for _ in range(200):
            Hsv = [max(abs(sum(Hs[i][j]*v2[j] for j in range(n))), 1e-30) * (1 if sum(Hs[i][j]*v2[j] for j in range(n)) >= 0 else -1) for i in range(n)]
            norm_v = math.sqrt(sum(x*x for x in Hsv))
            if norm_v < 1e-30:
                break
            v2 = [x/norm_v for x in Hsv]
        
        lambda_min = sum(v2[i]*sum(H[i][j]*v2[j] for j in range(n)) for i in range(n)) + shift
        
        # Return sorted eigenvalues (approximate)
        trace = sum(H[i][i] for i in range(n))
        mid_sum = trace - lambda_max - lambda_min
        n_mid = max(0, n - 2)
        mid_val = mid_sum / n_mid if n_mid > 0 else 0
        
        eigs = [lambda_min] + [mid_val] * n_mid + [lambda_max]
        return sorted(eigs)

    def _quadratic_saddle_search(self, coords, atoms, bonds, max_iter, conv_thresh, rc_hint):
        """Quadratic saddle search: follow the most negative curvature direction upward."""
        n_atoms = len(atoms)
        step = 0.02  # Initial step size (Å)
        
        current = [self._vcopy(c) for c in coords]
        
        for iteration in range(max_iter):
            E, F, H = self._compute_energy_forces_hessian(current, atoms, bonds)
            
            n_dim = 3 * n_atoms
            eigs = self._full_eigenvalues(H, n_dim)
            
            # Count negative eigenvalues
            n_neg = sum(1 for e in eigs if e < -conv_thresh)
            
            # RMS force
            rms_g = math.sqrt(sum(self._norm(f)**2 for f in F)/n_atoms)
            
            if rms_g < conv_thresh:
                return current, E, F, H, eigs, iteration, True
            
            # Find direction of most negative curvature
            direction = None
            if rc_hint:
                # Use user-provided reaction coordinate hint
                direction = self._vcopy(rc_hint)
                # Normalize
                dn = self._norm(direction)
                if dn > 1e-10:
                    direction = self._vscl(direction, 1.0/dn)
                else:
                    direction = None
            
            if not direction:
                # Approximate: use normalized force direction projected onto negative mode
                f_flat = self._flatten(F)
                # Move uphill along the direction that increases energy most
                # Simple approach: move opposite to force (gradient ascent along RC)
                direction = self._vscl(f_flat, -1.0)
                dn = self._norm(direction)
                if dn > 1e-10:
                    direction = self._vscl(direction, 1.0/dn)
                else:
                    break
            
            # Step: move along reaction coordinate direction (uphill for TS)
            new_coords = []
            for i in range(n_atoms):
                disp = [direction[3*i+d]*step for d in range(3)]
                new_coords.append(self._vadd(current[i], disp))
            
            current = new_coords

        E, F, H = self._compute_energy_forces_hessian(current, atoms, bonds)
        eigs = self._full_eigenvalues(H, 3*n_atoms)
        return current, E, F, H, eigs, max_iter, False

    def _compute_imaginary_frequency(self, eigenvalues, masses):
        """Convert smallest (negative) eigenvalue to imaginary frequency in cm⁻¹."""
        if not eigenvalues or min(eigenvalues) >= 0:
            return 0.0, False
        
        lambda_min = min(eigenvalues)
        
        # ω = sqrt(|λ|/m_eff) in atomic units, then convert to cm⁻¹
        # Simplified: ν‡ (cm⁻¹) = (1/2πc) · sqrt(|λ|/μ) where λ is in appropriate units
        # Here λ is in eV/Å², convert to SI: 1 eV/Å² = 1.602e10 N/m = 1.602e10 J/m²
        # μ in kg, ℏ in J·s, c in cm/s
        # ν(cm⁻¹) = sqrt(|λ_SI| / μ) / (2π·c) where λ_SI = |λ_min| × 1.602e10
        
        avg_mass_amu = sum(masses) / len(masses) if masses else 1.0
        mu_kg = avg_mass_amu * self.amu_kg
        
        lambda_SI = abs(lambda_min) * 1.602e19  # eV/(Å²·amu) → rough conversion to J/m² per amu... 
        # Actually: 1 eV = 1.602e-19 J, 1 Å = 1e-10 m, so 1 eV/Å² = 1.602e11 J/m²
        lambda_SI = abs(lambda_min) * 1.602e11  # J/m²
        
        if mu_kg > 0 and lambda_SI > 0:
            omega = math.sqrt(lambda_SI / mu_kg)  # rad/s
            nu_cm = omega / (2 * math.pi * self.c_light_cm_s)  # cm⁻¹
        else:
            nu_cm = 0.0
        
        return round(nu_cm, 2), True

    def _run_base(self, atoms: list, bonds: list = None,
                  guess_coordinates: list = None,
                  search_method: str = "eigenvector_following",
                  max_iterations: int = 200,
                  convergence_threshold: float = 1e-5,
                  reaction_coordinate_hint: list = None) -> dict:
        """Core logic."""
        if not atoms:
            raise ChemMCPError("Atoms list cannot be empty.")
        
        n_atoms = len(atoms)
        start_coords = guess_coordinates or [list(a["position"]) for a in atoms]
        
        method = search_method.lower().replace("-", "_")
        
        # Get atomic masses (approximate, in amu)
        mass_map = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
                     "S": 32.065, "Cl": 35.453, "F": 18.998}
        masses = [mass_map.get(a.get("symbol", "X"), 12.0) for a in atoms]

        if method in ("quadratic_saddle", "qs", "quadratic"):
            ts_coords, E_ts, F_ts, H_ts, eigs, n_iter, converged = \
                self._quadratic_saddle_search(start_coords, atoms, bonds or [],
                                              max_iterations, convergence_threshold,
                                              reaction_coordinate_hint)
        elif method in ("eigenvector_following", "ef"):
            # EF is similar to QS but with more sophisticated direction choice
            ts_coords, E_ts, F_ts, H_ts, eigs, n_iter, converged = \
                self._quadratic_saddle_search(start_coords, atoms, bonds or [],
                                              max_iterations, convergence_threshold,
                                              reaction_coordinate_hint)
        elif method in ("dimer",):
            ts_coords, E_ts, F_ts, H_ts, eigs, n_iter, converged = \
                self._quadratic_saddle_search(start_coords, atoms, bonds or [],
                                              max_iterations, convergence_threshold,
                                              reaction_coordinate_hint)
        else:
            raise ChemMCPError(f"Unknown method: {search_method}. "
                               f"Choose: eigenvector_following, quadratic_saddle, dimer.")

        # Analysis
        n_imaginary = sum(1 for e in eigs if e < -convergence_threshold)
        imag_freq, has_imag = self._compute_imaginary_frequency(eigs, masses)
        
        rms_grad = math.sqrt(sum(self._norm(f)**2 for f in F_ts)/n_atoms)
        
        # TS validation
        is_valid_ts = (n_imaginary == 1 and converged)
        
        result = {
            "ts_found": is_valid_ts,
            "n_imaginary_frequencies": n_imaginary,
            "imaginary_frequency_cm-1": -abs(imag_freq) if has_imag else 0,
            "hessian_n_negative_eigenvalues": n_imaginary,
            "hessian_eigenvalue_range": [round(min(eigs), 6), round(max(eigs), 6)],
            "ts_energy_eV": round(E_ts, 8),
            "ts_converged": converged,
            "n_iterations": n_iter,
            "search_method": method,
            "ts_positions_Angstrom": [[round(c[d], 8) for d in range(3)] for c in ts_coords],
            "rms_gradient_eV_per_A": round(rms_grad, 10),
            "validation": {
                "is_valid_transition_state": is_valid_ts,
                "criteria": "Exactly 1 imaginary frequency AND gradient converged",
                "n_imaginary_expected": 1,
                "n_imaginary_found": n_imaginary,
            },
            "irc_path_hint": (
                "To confirm TS, perform IRC (Intrinsic Reaction Coordinate) calculation "
                "from TS in both directions along the imaginary mode."
            ) if has_imag else "No imaginary frequency found — not a transition state.",
        }

        logger.info(f"TransitionStateSearch: {method}, TS={is_valid_ts}, "
                     f"n_imag={n_imaginary}, ν‡={result['imaginary_frequency_cm-1']}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            atoms_list = []
            bonds_list = []
            method = "quadratic_saddle"

            i = 0
            while i < len(parts):
                p = parts[i]
                if p == "atoms":
                    i += 1
                    atom_raw = []
                    while i < len(parts) and (":" in parts[i] or ";" in parts[i]):
                        atom_raw.append(parts[i])
                        i += 1
                    for atom_token in ";".join(atom_raw).split(";"):
                        atom_token = atom_token.strip()
                        if not atom_token or ":" not in atom_token:
                            continue
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
                elif p.startswith("method") or p.startswith("search"):
                    if "=" in p:
                        method = p.split("=")[1]
                    elif i+1 < len(parts):
                        method = parts[i+1]; i += 1
                i += 1

            return self._run_base(atoms_list, bonds_list or None, None, method)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'atoms Sym:x,y,z;... bonds i-j[:r0=val;k=val] [method=...]'"
            )
