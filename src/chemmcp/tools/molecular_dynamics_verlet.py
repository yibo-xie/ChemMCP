import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Default Lennard-Jones parameters (reduced units)
_LJ_EPSILON_DEFAULT = 1.0  # ε: well depth
_LJ_SIGMA_DEFAULT = 1.0   # σ: finite distance at which potential is zero


def _distance(r1: List[float], r2: List[float], box: float) -> List[float]:
    """Compute minimum-image displacement vector with PBC."""
    dr = [r2[i] - r1[i] for i in range(3)]
    for i in range(3):
        if dr[i] > box / 2:
            dr[i] -= box
        elif dr[i] < -box / 2:
            dr[i] += box
    return dr


def _vec_norm(v: List[float]) -> float:
    return math.sqrt(sum(x ** 2 for x in v))


@ChemMCPManager.register_tool
class MolecularDynamicsVerlet(BaseTool):
    """
    Verlet 积分器工具。
    分子动力学模拟核心，使用位置 Verlet 算法进行时间演化。
    """
    __version__ = "0.1.0"
    name = "MolecularDynamicsVerlet"
    func_name = "verlet_integrate"
    description = "Position Verlet integrator for molecular dynamics simulation with Lennard-Jones, harmonic, or Coulomb force fields."
    implementation_description = "Implements the position Verlet algorithm: r(t+dt)=r+v*dt+0.5*a*dt², compute new a, then v(t+dt)=v+0.5*(a_old+a_new)*dt. Supports Lennard-Jones (12-6), harmonic bond, and Coulomb potentials with periodic boundary conditions and minimum image convention."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Molecular Dynamics", "Verlet Integration", "Simulation", "Computational Chemistry", "Physics", "Time Evolution"]
    required_envs = []

    code_input_sig = [
        ("positions", "list", "N/A", "List of [x,y,z] coordinates for each atom."),
        ("velocities", "list", "N/A", "List of [vx,vy,vz] velocities for each atom."),
        ("masses", "list", "N/A", "List of atomic masses (float)."),
        ("force_function", "str", "lj", "Force field type: 'lj' (Lennard-Jones), 'harmonic', 'coulomb'."),
        ("dt", "float", "0.001", "Time step size."),
        ("n_steps", "int", "1000", "Number of integration steps."),
        ("box_size", "float", "10.0", "Cubic box size for periodic boundary conditions."),
        ("ff_params", "dict", "{}", "Force field parameters dict (epsilon, sigma, k, r0, q, etc.)."),
        ("trajectory_sample_every", "int", "10", "Save trajectory frame every N steps."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all MD parameters."),
    ]

    output_sig = [
        ("trajectory", "list", "List of position snapshots (each is list of [x,y,z] per atom)."),
        ("final_positions", "list", "Final positions after integration."),
        ("final_velocities", "list", "Final velocities after integration."),
        ("total_energy_over_time", "list", "[step, KE, PE, Total] at each sampled step."),
        ("kinetic_energy_final", "float", "Final kinetic energy."),
        ("potential_energy_final", "float", "Final potential energy."),
        ("total_energy_final", "float", "Final total energy."),
        ("temperature", "float", "Instantaneous temperature from KE."),
        ("n_atoms", "int", "Number of atoms."),
        ("n_steps_completed", "int", "Total steps completed."),
        ("force_field_type", "str", "Force field used."),
        ("energy_drift", "float", "(E_final - E_initial) / |E_initial| relative drift."),
    ]

    examples = [
        {
            "code_input": {
                "positions": [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
                "velocities": [[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]],
                "masses": [1.0, 1.0],
                "force_function": "lj",
                "dt": 0.005,
                "n_steps": 100,
                "box_size": 10.0,
            },
            "text_input": {"params_str": '{"positions":[[0,0,0],[1.5,0,0]],"velocities":[[0,0,0],[0,0.5,0]],"masses":[1,1],"force_function":"lj","dt":0.005,"n_steps":100,"box_size":10}'},
            "output": {
                "n_atoms": 2,
                "n_steps_completed": 100,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- Force computation ----

    @staticmethod
    def _compute_forces_lj(positions: List[List[float]], box: float, eps: float, sig: float) -> tuple:
        """Compute LJ forces and potential energy."""
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe = 0.0
        rcut = min(box / 2.0, 3.0 * sig)  # cutoff

        for i in range(n):
            for j in range(i + 1, n):
                dr = _distance(positions[i], positions[j], box)
                r = _vec_norm(dr)
                if r < 1e-10 or r > rcut:
                    continue
                # LJ: V = 4ε[(σ/r)^12 - (σ/r)^6]
                # F = -dV/dr * (dr_vec/r)
                sr = sig / r
                sr6 = sr ** 6
                sr12 = sr6 ** 2
                f_mag = 24.0 * eps * (2.0 * sr12 - sr6) / r  # magnitude along r direction
                pe += 4.0 * eps * (sr12 - sr6)

                # Apply forces (Newton's 3rd law)
                fi = [f_mag * dr[k] / r for k in range(3)]
                fj = [-fi[k] for k in range(3)]
                for k in range(3):
                    forces[i][k] += fi[k]
                    forces[j][k] += fj[k]

        return forces, pe

    @staticmethod
    def _compute_forces_harmonic(positions: List[List[float]], box: float, k: float, r0: float) -> tuple:
        """Harmonic spring between all pairs: V = 0.5*k*(r-r0)²."""
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                dr = _distance(positions[i], positions[j], box)
                r = _vec_norm(dr)
                if r < 1e-10:
                    continue
                f_mag = k * (r - r0)
                pe += 0.5 * k * (r - r0) ** 2
                fi = [f_mag * dr[k] / r for k in range(3)]
                fj = [-fi[k] for k in range(3)]
                for kk in range(3):
                    forces[i][kk] += fi[kk]
                    forces[j][kk] += fj[kk]
        return forces, pe

    @staticmethod
    def _compute_forces_coulomb(positions: List[List[float]], box: float, charges: List[float], ke: float = 332.0) -> tuple:
        """Coulomb interaction: V = ke*q1*q2/r."""
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                dr = _distance(positions[i], positions[j], box)
                r = _vec_norm(dr)
                if r < 1e-10:
                    continue
                f_mag = ke * charges[i] * charges[j] / (r * r)
                pe += ke * charges[i] * charges[j] / r
                fi = [f_mag * dr[k] / r for k in range(3)]
                fj = [-fi[k] for k in range(3)]
                for kk in range(3):
                    forces[i][kk] += fi[kk]
                    forces[j][kk] += fj[kk]
        return forces, pe

    def _get_forces(self, positions, ff_type, box, params):
        """Dispatch to correct force routine."""
        if ff_type == "lj":
            eps = params.get("epsilon", _LJ_EPSILON_DEFAULT)
            sig = params.get("sigma", _LJ_SIGMA_DEFAULT)
            return self._compute_forces_lj(positions, box, eps, sig)
        elif ff_type == "harmonic":
            k = params.get("k", 100.0)
            r0 = params.get("r0", 1.0)
            return self._compute_forces_harmonic(positions, box, k, r0)
        elif ff_type == "coulomb":
            charges = params.get("charges", [1.0, -1.0])
            ke = params.get("ke", 332.0)
            return self._compute_forces_coulomb(positions, box, charges, ke)
        else:
            raise ChemMCPError(f"Unknown force_function '{ff_type}'. Use: lj, harmonic, coulomb.")

    @staticmethod
    def _kinetic_energy(velocities, masses):
        ke = 0.0
        for v, m in zip(velocities, masses):
            ke += 0.5 * m * sum(vi ** 2 for vi in v)
        return ke

    @staticmethod
    def _temperature_from_ke(ke, n_atoms, kb_local=1.0):
        """T = 2*KE / (3*N*kb) for 3D system."""
        dof = 3 * n_atoms
        return 2.0 * ke / (dof * kb_local) if dof > 0 else 0.0

    def _run_base(
        self,
        positions: List[List[float]],
        velocities: List[List[float]],
        masses: List[float],
        force_function: str = "lj",
        dt: float = 0.001,
        n_steps: int = 1000,
        box_size: float = 10.0,
        ff_params: Optional[dict] = None,
        trajectory_sample_every: int = 10,
    ) -> dict:
        """Core logic: Position Verlet integration."""
        n = len(positions)
        if n != len(velocities) or n != len(masses):
            raise ChemMCPError(f"Length mismatch: positions={n}, velocities={len(velocities)}, masses={len(masses)}")
        if n == 0:
            raise ChemMCPError("No atoms provided.")
        if dt <= 0:
            raise ChemMCPError(f"dt must be positive, got {dt}.")
        if n_steps < 1:
            raise ChemMCPError(f"n_steps must be >= 1, got {n_steps}.")

        params = ff_params or {}
        ff = force_function.lower().strip()
        box = box_size
        sample_every = max(trajectory_sample_every, 1)

        # Deep copy initial state
        pos = [[p[i] for i in range(3)] for p in positions]
        vel = [[v[i] for i in range(3)] for v in velocities]

        # Initial forces
        forces, pe = self._get_forces(pos, ff, box, params)

        # Accelerations
        acc = [[forces[i][j] / masses[i] for j in range(3)] for i in range(n)]

        trajectory = []
        energy_history = []
        ke_init = self._kinetic_energy(vel, masses)
        e_total_init = ke_init + pe

        for step in range(n_steps):
            # Position Verlet Step 1: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt²
            for i in range(n):
                for j in range(3):
                    pos[i][j] += vel[i][j] * dt + 0.5 * acc[i][j] * dt * dt
                    # PBC wrap
                    while pos[i][j] < 0:
                        pos[i][j] += box
                    while pos[i][j] >= box:
                        pos[i][j] -= box

            # Compute new forces at new positions
            forces_new, pe_new = self._get_forces(pos, ff, box, params)
            acc_new = [[forces_new[i][j] / masses[i] for j in range(3)] for i in range(n)]

            # Position Verlet Step 2: v(t+dt) = v(t) + 0.5*(a(t)+a(t+dt))*dt
            for i in range(n):
                for j in range(3):
                    vel[i][j] += 0.5 * (acc[i][j] + acc_new[i][j]) * dt

            forces, pe, acc = forces_new, pe_new, acc_new

            # Sample trajectory & energies
            if step % sample_every == 0 or step == n_steps - 1:
                traj_frame = [[pos[i][j] for j in range(3)] for i in range(n)]
                trajectory.append(traj_frame)
                ke_now = self._kinetic_energy(vel, masses)
                e_tot = ke_now + pe
                energy_history.append([step, round(ke_now, 12), round(pe, 12), round(e_tot, 12)])

        ke_final = self._kinetic_energy(vel, masses)
        e_total_final = ke_final + pe
        drift = (e_total_final - e_total_init) / abs(e_total_init) if abs(e_total_init) > 1e-30 else 0.0
        temp = self._temperature_from_ke(ke_final, n)

        logger.info(
            f"MD Verlet ({ff}): {n_steps} steps, {n} atoms, "
            f"E_init={e_total_init:.6g}, E_final={e_total_final:.6g}, drift={drift:.2e}"
        )

        return {
            "trajectory": trajectory,
            "final_positions": [[pos[i][j] for j in range(3)] for i in range(n)],
            "final_velocities": [[vel[i][j] for j in range(3)] for i in range(n)],
            "total_energy_over_time": energy_history,
            "kinetic_energy_final": round(ke_final, 12),
            "potential_energy_final": round(pe, 12),
            "total_energy_final": round(e_total_final, 12),
            "temperature": round(temp, 8),
            "n_atoms": n,
            "n_steps_completed": n_steps,
            "force_field_type": ff,
            "energy_drift": round(drift, 12),
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
