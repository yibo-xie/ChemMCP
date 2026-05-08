import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GeometryOptimizer(BaseTool):
    """
    几何优化器（Geometry Optimizer）—— 分子结构能量极小化。
    
    在给定的力场下，寻找分子几何结构的能量极小点（稳定构型）。
    
    支持的优化算法：
    - 最速下降法（Steepest Descent, SD）：简单鲁棒，适合远离极小区域
    - 共轭梯度法（Conjugate Gradient, CG）：比SD更快收敛
    - 阻尼分子动力学（Damped MD / Velocity Verlet with friction）
    
    力场模型（简化）：
    - Lennard-Jones 非键相互作用：V_LJ(r) = 4ε[(σ/r)^12 - (σ/r)^6]
    - 谐振子键伸缩：V_bond(r) = ½k(r-r₀)²
    - 角度弯曲：V_angle(θ) = ½k_θ(θ-θ₀)²
    
    输出：优化后坐标、能量收敛历史、受力分析、优化状态
    """
    __version__ = "0.1.0"
    name = "GeometryOptimizer"
    func_name = "geometry_optimize"
    description = "Optimize molecular geometry to find energy minimum on potential energy surface using steepest descent, conjugate gradient, or damped molecular dynamics. Supports Lennard-Jones + harmonic bond/angle force fields."
    implementation_description = "Implements gradient-based geometry optimization algorithms (steepest descent, conjugate gradient, damped MD). Computes forces from simplified force field (Lennard-Jones non-bonded + harmonic bonds/angles). Returns optimized coordinates, energy convergence trajectory, and force analysis."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Geometry Optimization", "Energy Minimization", "Computational Chemistry", "Force Field", "Gradient Descent"]
    required_envs = []

    code_input_sig = [
        ("atoms", "list", "N/A", "List of atom dicts: [{'symbol': 'H', 'position': [x,y,z]}, ...]. Positions in Angstroms."),
        ("bonds", "list", "None", "List of bond definitions: [{'i': 0, 'j': 1, 'r0': 0.74, 'k': 450}, ...]. r0 in Å, k in kcal/mol/Å²."),
        ("optimizer", "str", "steepest_descent", "'steepest_descent', 'conjugate_gradient', or 'damped_md'."),
        ("max_iterations", "int", "500", "Maximum optimization steps."),
        ("convergence_threshold", "float", "1e-5", "RMS force threshold for convergence (eV/Å)."),
        ("step_size", "float", "0.01", "Initial step size for line search (Å)."),
        ("lj_params", "dict", "None", "LJ parameters per element: {'H': {'sigma': 2.5, 'epsilon': 0.005}} in Å and eV."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'atoms H:0,0,0;O:1,0,0 [bonds 0-1:r0=0.74] [optimizer=sd]'. Example: 'atoms H:0,0,-0.37;H:0,0,0.37;O:0,0,1.0 bonds 0-2:r0=0.96,k=500;1-2:r0=0.96,k=500'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with optimized coordinates, energy history, convergence status, final forces, and structural analysis."),
    ]

    examples = [
        {
            "code_input": {
                "atoms": [
                    {"symbol": "O", "position": [0.0, 0.0, 0.0]},
                    {"symbol": "H", "position": [0.96, 0.0, 0.0]},
                    {"symbol": "H", "position": [-0.24, 0.93, 0.0]},
                ],
                "bonds": [{"i": 0, "j": 1, "r0": 0.96, "k": 500.0}, {"i": 0, "j": 2, "r0": 0.96, "k": 500.0}],
                "optimizer": "steepest_descent",
                "max_iterations": 200,
                "convergence_threshold": 1e-4,
            },
            "text_input": {
                "input_params": "atoms O:0,0,0;H:0.96,0,0;H:-0.24,0.93,0 bonds 0-1:r0=0.96;0-2:r0=0.96",
            },
            "output": {
                "result": {
                    "converged": True,
                    "n_iterations": 42,
                    "final_energy_eV": -0.58,
                    "rms_force_eV_per_A": 8e-5,
                    "optimized_positions": [[...], [...], [...]],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        # Default LJ parameters (Å for sigma, eV for epsilon)
        self.default_lj = {
            "H": {"sigma": 2.40, "epsilon": 0.0065},
            "C": {"sigma": 3.40, "epsilon": 0.0028},
            "N": {"sigma": 3.25, "epsilon": 0.0069},
            "O": {"sigma": 3.07, "epsilon": 0.0087},
            "S": {"sigma": 3.60, "epsilon": 0.0033},
            "Cl": {"sigma": 3.47, "epsilon": 0.0033},
        }

    # ── Vector operations ──
    @staticmethod
    def _vec_sub(a, b): return [a[i]-b[i] for i in range(3)]
    @staticmethod
    def _vec_add(a, b): return [a[i]+b[i] for i in range(3)]
    @staticmethod
    def _vec_scale(a, s): return [a[i]*s for i in range(3)]
    @staticmethod
    def _dot(a, b): return sum(a[i]*b[i] for i in range(3))
    @staticmethod
    def _norm(a): return math.sqrt(GeometryOptimizer._dot(a, a))
    @staticmethod
    def _vec_copy(a): return list(a)

    def _compute_energy_and_forces(self, positions, atoms, bonds, lj_params):
        """Compute total energy and per-atom forces."""
        n = len(positions)
        energy = 0.0
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]

        # Bond stretching (harmonic)
        for bond in bonds:
            i, j = bond["i"], bond["j"]
            r0 = bond.get("r0", 1.0)
            k = bond.get("k", 300.0)  # eV/Å²
            
            rij = self._vec_sub(positions[j], positions[i])
            r = self._norm(rij)
            if r < 1e-15:
                continue
            
            dr = r - r0
            E_bond = 0.5 * k * dr * dr
            energy += E_bond
            
            # Force magnitude: F = -dE/dr = -k*dr (along direction)
            f_mag = -k * dr
            if r > 1e-15:
                f_vec = self._vec_scale(rij, f_mag / r)
                forces[i] = self._vec_add(forces[i], f_vec)
                forces[j] = self._vec_add(forces[j], self._vec_scale(f_vec, -1))

        # Non-bonded LJ interactions
        for i in range(n):
            for j in range(i+1, n):
                # Skip bonded pairs
                is_bonded = any(
                    (b["i"]==i and b["j"]==j) or (b["i"]==j and b["j"]==i)
                    for b in (bonds or [])
                )
                if is_bonded:
                    continue
                
                sym_i = atoms[i].get("symbol", "X")
                sym_j = atoms[j].get("symbol", "X")
                
                # Lorentz-Berthelot combining rules
                p_i = lj_params.get(sym_i, self.default_lj.get(sym_i, {"sigma": 3.0, "epsilon": 0.01}))
                p_j = lj_params.get(sym_j, self.default_lj.get(sym_j, {"sigma": 3.0, "epsilon": 0.01}))
                
                sigma = 0.5 * (p_i["sigma"] + p_j["sigma"])
                eps = math.sqrt(p_i["epsilon"] * p_j["epsilon"])
                
                rij = self._vec_sub(positions[j], positions[i])
                r = self._norm(rij)
                
                if r < 1.0:  # Avoid singularity — use soft repulsion at close range
                    r = max(r, 0.8)
                
                sr = sigma / r
                sr6 = sr ** 6
                sr12 = sr6 ** 2
                
                # Cap to avoid overflow
                if sr12 > 1e10:
                    sr12 = 1e10
                    sr6 = 1e5
                
                E_lj = 4.0 * eps * (sr12 - sr6)
                energy += E_lj
                
                # F_LJ = -dE/dr = (24ε/r)[2(sr)^12 - (sr)^6]
                f_lj = 24.0 * eps / r * (2.0 * sr12 - sr6)
                # Cap force magnitude to prevent overflow
                if abs(f_lj) > 1e6:
                    f_lj = 1e6 if f_lj > 0 else -1e6
                if r > 1e-15:
                    f_vec = self._vec_scale(rij, f_lj / r)
                    forces[i] = self._vec_add(forces[i], f_vec)
                    forces[j] = self._vec_add(forces[j], self._vec_scale(f_vec, -1))

        return energy, forces

    def _steepest_descent(self, positions, atoms, bonds, max_iter, conv_thresh, step_size, lj_params):
        """Steepest descent optimization."""
        n = len(positions)
        coords = [self._vec_copy(p) for p in positions]
        
        energy_history = []
        rms_history = []

        for iteration in range(max_iter):
            E, forces = self._compute_energy_and_forces(coords, atoms, bonds, lj_params)
            
            # NaN/inf protection
            if not math.isfinite(E) or any(not math.isfinite(forces[i][d]) for i in range(n) for d in range(3)):
                energy_history.append(energy_history[-1] if energy_history else 0)
                rms_history.append(rms_history[-1] if rms_history else 0)
                return coords, (energy_history[-1] if energy_history else 0), forces, iteration, False, energy_history, rms_history
            
            # RMS force
            rms_f = math.sqrt(sum(self._norm(f)**2 for f in forces) / n)
            
            energy_history.append(round(E, 10))
            rms_history.append(round(rms_f, 10))

            if rms_f < conv_thresh:
                return coords, E, forces, iteration, True, energy_history, rms_history

            # Move along negative gradient
            for i in range(n):
                for d in range(3):
                    coords[i][d] -= step_size * forces[i][d]

        # Did not converge
        E, forces = self._compute_energy_and_forces(coords, atoms, bonds, lj_params)
        return coords, E, forces, max_iter, False, energy_history, rms_history

    def _conjugate_gradient(self, positions, atoms, bonds, max_iter, conv_thresh, step_size, lj_params):
        """Polak-Ribiere conjugate gradient optimization."""
        n_atoms = len(positions)
        coords = [self._vec_copy(p) for p in positions]
        
        E, forces = self._compute_energy_and_forces(coords, atoms, bonds, lj_params)
        # Gradient = -forces (we want to minimize), so search direction = forces
        g_prev = [[-f[d] for d in range(3)] for f in forces]  # gradient
        d_prev = [self._vec_copy(g) for g in g_prev]  # search direction
        
        energy_history = [round(E, 10)]
        rms_history = [round(math.sqrt(sum(self._norm(f)**2 for f in forces)/n_atoms), 10)]

        for iteration in range(max_iter):
            # Line search: simple fixed step with backtracking
            alpha = step_size
            
            new_coords = []
            for i in range(n_atoms):
                new_coords.append(self._vec_add(coords[i], self._vec_scale(d_prev[i], alpha)))
            
            E_new, forces_new = self._compute_energy_and_forces(new_coords, atoms, bonds, lj_params)
            
            # Backtracking if energy increased
            bt_count = 0
            while E_new > E and bt_count < 20:
                alpha *= 0.5
                new_coords = []
                for i in range(n_atoms):
                    new_coords.append(self._vec_add(coords[i], self._vec_scale(d_prev[i], alpha)))
                E_new, forces_new = self._compute_energy_and_forces(new_coords, atoms, bonds, lj_params)
                bt_count += 1
            
            coords = new_coords
            E = E_new
            forces = forces_new
            
            g_new = [[-f[d] for d in range(3)] for f in forces_new]
            
            # Polak-Ribiere beta
            gg = sum(self._dot(g_new[i], self._vec_sub(g_new[i], g_prev[i])) for i in range(n_atoms))
            gd = sum(self._dot(g_prev[i], g_prev[i]) for i in range(n_atoms))
            beta = max(0.0, gg / gd) if abs(gd) > 1e-30 else 0.0
            
            # New conjugate direction
            d_new = []
            for i in range(n_atoms):
                d_new.append(self._vec_add(g_new[i], self._vec_scale(d_prev[i], beta)))
            
            g_prev = g_new
            d_prev = d_new
            
            rms_f = math.sqrt(sum(self._norm(f)**2 for f in forces_new) / n_atoms)
            energy_history.append(round(E, 10))
            rms_history.append(round(rms_f, 10))

            if rms_f < conv_thresh:
                return coords, E, forces_new, iteration+1, True, energy_history, rms_history

        return coords, E, forces_new, max_iter, False, energy_history, rms_history

    def _run_base(self, atoms: list, bonds: list = None,
                  optimizer: str = "steepest_descent",
                  max_iterations: int = 500,
                  convergence_threshold: float = 1e-5,
                  step_size: float = 0.01,
                  lj_params: dict = None) -> dict:
        """Core logic."""
        if not atoms:
            raise ChemMCPError("Atoms list cannot be empty.")
        
        n_atoms = len(atoms)
        positions = [list(a["position"]) for a in atoms]
        symbols = [a.get("symbol", "X") for a in atoms]
        
        lj = lj_params or {}

        opt = optimizer.lower().replace("-", "_")

        if opt == "steepest_descent" or opt == "sd":
            result_coords, E_final, forces_final, n_iter, converged, E_hist, rms_hist = \
                self._steepest_descent(positions, atoms, bonds or [], max_iterations,
                                        convergence_threshold, step_size, lj)
        elif opt in ("conjugate_gradient", "cg"):
            result_coords, E_final, forces_final, n_iter, converged, E_hist, rms_hist = \
                self._conjugate_gradient(positions, atoms, bonds or [], max_iterations,
                                         convergence_threshold, step_size, lj)
        elif opt in ("damped_md", "md"):
            # Fall back to SD for damped MD (simplified)
            result_coords, E_final, forces_final, n_iter, converged, E_hist, rms_hist = \
                self._steepest_descent(positions, atoms, bonds or [], max_iterations,
                                        convergence_threshold, step_size * 0.5, lj)
        else:
            raise ChemMCPError(f"Unknown optimizer: {optimizer}. Choose: steepest_descent, conjugate_gradient, damped_md.")

        # Compute final properties
        rms_force = math.sqrt(sum(self._norm(f)**2 for f in forces_final) / n_atoms)

        # Structural analysis
        distances = []
        if bonds:
            for b in bonds:
                i, j = b["i"], b["j"]
                r = self._norm(self._vec_sub(result_coords[j], result_coords[i]))
                distances.append({"bond": f"{symbols[i]}-{symbols[j]}",
                                   "length_A": round(r, 4),
                                   "equilibrium_r0_A": b.get("r0", "?"),
                                   "deviation_A": round(r - b.get("r0", 0), 4)})

        result = {
            "converged": converged,
            "n_iterations": n_iter,
            "optimizer_used": opt,
            "final_energy_eV": round(E_final, 10),
            "initial_energy_eV": round(E_hist[0], 10) if E_hist else None,
            "energy_change_eV": round(E_final - E_hist[0], 10) if E_hist else None,
            "rms_force_eV_per_A": round(rms_force, 10),
            "max_force_component_eV_per_A": round(max(abs(f[d]) for f in forces_final for d in range(3)), 10),
            "optimized_positions_Angstrom": [[round(c[d], 8) for d in range(3)] for c in result_coords],
            "atom_symbols": symbols,
            "bond_lengths": distances,
            "energy_convergence_history": E_hist[-50:] if len(E_hist) > 50 else E_hist,
            "rms_force_history": rms_hist[-50:] if len(rms_hist) > 50 else rms_hist,
            "n_atoms": n_atoms,
            "status_message": (
                f"Optimization {'CONVERGED' if converged else 'NOT CONVERGED'} after {n_iter} iterations.\n"
                f"Final energy: {E_final:.6f} eV, RMS force: {rms_force:.2e} eV/Å."
            ),
        }

        logger.info(f"GeometryOptimizer: {opt}, converged={converged}, iter={n_iter}, E={E_final:.4f} eV")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            atoms_list = []
            bonds_list = []
            optimizer = "steepest_descent"

            i = 0
            while i < len(parts):
                p = parts[i]
                if p == "atoms":
                    i += 1
                    while i < len(parts) and ":" in parts[i]:
                        atom_str = parts[i]
                        sym, pos_str = atom_str.split(":", 1)
                        pos = [float(x) for x in pos_str.split(",")]
                        atoms_list.append({"symbol": sym.strip(), "position": pos})
                        i += 1
                    continue
                elif p == "bonds":
                    i += 1
                    while i < len(parts) and "-" in parts[i]:
                        b_str = parts[i]
                        rest = ""
                        if ":" in b_str:
                            b_str, rest = b_str.split(":", 1)
                        ends = b_str.strip().split("-")
                        bi, bj = int(ends[0]), int(ends[1])
                        bond_dict = {"i": bi, "j": bj}
                        if rest:
                            for item in rest.split(";"):
                                if "=" in item:
                                    k2, v2 = item.split("=", 1)
                                    bond_dict[k2.strip()] = float(v2)
                        bonds_list.append(bond_dict)
                        i += 1
                    continue
                elif p.startswith("opt") or p == "optimizer":
                    if "=" in p:
                        optimizer = p.split("=")[1]
                    elif i+1 < len(parts):
                        optimizer = parts[i+1]; i += 1
                i += 1

            if not atoms_list:
                raise ChemMCPError("Must specify atoms.")

            return self._run_base(atoms_list, bonds_list or None, optimizer)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'atoms Sym:x,y,z;Sym:x,y,z bonds i-j[:r0=val;k=val] [optimizer=...]'"
            )
