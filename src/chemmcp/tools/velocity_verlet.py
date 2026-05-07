import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _dist_vec(r1: List[float], r2: List[float], box: float) -> List[float]:
    """Minimum-image displacement with PBC."""
    dr = [r2[i] - r1[i] for i in range(3)]
    for i in range(3):
        if dr[i] > box / 2:
            dr[i] -= box
        elif dr[i] < -box / 2:
            dr[i] += box
    return dr


def _norm(v: List[float]) -> float:
    return math.sqrt(sum(x ** 2 for x in v))


@ChemMCPManager.register_tool
class VelocityVerlet(BaseTool):
    """
    速度 Verlet 算法工具。
    保守系统时间演化，使用速度 Verlet（类蛙跳）积分算法，具有更好的能量守恒性。
    """
    __version__ = "0.1.0"
    name = "VelocityVerlet"
    func_name = "velocity_verlet_integrate"
    description = "Velocity Verlet (leapfrog) integrator for conservative system time evolution with energy conservation monitoring."
    implementation_description = "Implements the velocity Verlet algorithm: v(t+dt/2)=v(t)+0.5*a*dt, r(t+dt)=r(t)+v(t+dt/2)*dt, compute a(t+dt), v(t+dt)=v(t+dt/2)+0.5*a(t+dt)*dt. Symplectic integrator with excellent long-term energy conservation. Supports LJ, harmonic, Coulomb force fields."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Velocity Verlet", "Molecular Dynamics", "Symplectic Integrator", "Simulation", "Energy Conservation", "Physics"]
    required_envs = []

    code_input_sig = [
        ("positions", "list", "N/A", "List of [x,y,z] coordinates for each atom."),
        ("velocities", "list", "N/A", "List of [vx,vy,vz] velocities for each atom."),
        ("masses", "list", "N/A", "List of atomic masses (float)."),
        ("force_function", "str", "lj", "Force field type: 'lj', 'harmonic', 'coulomb'."),
        ("dt", "float", "0.001", "Time step size."),
        ("n_steps", "int", "1000", "Number of integration steps."),
        ("box_size", "float", "10.0", "Cubic box size for PBC."),
        ("ff_params", "dict", "{}", "Force field parameters."),
        ("trajectory_sample_every", "int", "10", "Save trajectory every N steps."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("trajectory", "list", "Position snapshots over time."),
        ("final_positions", "list", "Final positions after integration."),
        ("final_velocities", "list", "Final velocities after integration."),
        ("energies_over_time", "list", "[step, KE, PE, Total] at each sample."),
        ("kinetic_energy_final", "float", "Final kinetic energy."),
        ("potential_energy_final", "float", "Final potential energy."),
        ("total_energy_final", "float", "Final total energy."),
        ("temperature", "float", "Instantaneous temperature from KE."),
        ("energy_drift", "float", "Relative energy drift (E_final-E_init)/|E_init|."),
        ("max_energy_deviation", "float", "Maximum |E(t)-E(0)|/|E(0)| during simulation."),
        ("n_atoms", "int", "Number of atoms."),
        ("n_steps_completed", "int", "Total steps completed."),
        ("force_field_type", "str", "Force field used."),
        ("time_array", "list", "Time values at each sampled step."),
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
                "energy_drift": 0.0,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- Force routines ----

    @staticmethod
    def _lj_forces(positions, box, eps=1.0, sig=1.0):
        n = len(positions)
        forces = [[0.0] * 3 for _ in range(n)]
        pe = 0.0
        rcut = min(box / 2.0, 3.0 * sig)
        for i in range(n):
            for j in range(i + 1, n):
                dr = _dist_vec(positions[i], positions[j], box)
                r = _norm(dr)
                if r < 1e-10 or r > rcut:
                    continue
                sr = sig / r
                sr6 = sr ** 6
                sr12 = sr6 ** 2
                f_mag = 24.0 * eps * (2.0 * sr12 - sr6) / r
                pe += 4.0 * eps * (sr12 - sr6)
                fk = [f_mag * dr[k] / r for k in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]
        return forces, pe

    @staticmethod
    def _harmonic_forces(positions, box, k=100.0, r0=1.0):
        n = len(positions)
        forces = [[0.0] * 3 for _ in range(n)]
        pe = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                dr = _dist_vec(positions[i], positions[j], box)
                r = _norm(dr)
                if r < 1e-10:
                    continue
                f_mag = k * (r - r0)
                pe += 0.5 * k * (r - r0) ** 2
                fk = [f_mag * dr[k] / r for k in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]
        return forces, pe

    @staticmethod
    def _coulomb_forces(positions, box, charges, ke=332.0):
        n = len(positions)
        forces = [[0.0] * 3 for _ in range(n)]
        pe = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                dr = _dist_vec(positions[i], positions[j], box)
                r = _norm(dr)
                if r < 1e-10:
                    continue
                f_mag = ke * charges[i] * charges[j] / (r * r)
                pe += ke * charges[i] * charges[j] / r
                fk = [f_mag * dr[k] / r for k in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]
        return forces, pe

    def _compute_forces(self, pos, ff_type, box, params):
        if ff_type == "lj":
            return self._lj_forces(pos, box, params.get("epsilon", 1.0), params.get("sigma", 1.0))
        elif ff_type == "harmonic":
            return self._harmonic_forces(pos, box, params.get("k", 100.0), params.get("r0", 1.0))
        elif ff_type == "coulomb":
            return self._coulomb_forces(pos, box, params.get("charges", [1.0, -1.0]), params.get("ke", 332.0))
        else:
            raise ChemMCPError(f"Unknown force_function '{ff_type}'. Use: lj, harmonic, coulomb.")

    @staticmethod
    def _ke(velocities, masses):
        return sum(0.5 * m * sum(v ** 2 for v in vel) for m, vel in zip(masses, velocities))

    @staticmethod
    def _temp(ke, n_atoms, kb=1.0):
        dof = 3 * n_atoms
        return 2.0 * ke / (dof * kb) if dof > 0 else 0.0

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
        """Core logic: Velocity Verlet integration."""
        n = len(positions)
        if n != len(velocities) or n != len(masses):
            raise ChemMCPError(f"Length mismatch: pos={n}, vel={len(velocities)}, mass={len(masses)}")
        if n == 0:
            raise ChemMCPError("No atoms provided.")
        if dt <= 0:
            raise ChemMCPError(f"dt must be positive, got {dt}.")
        if n_steps < 1:
            raise ChemMCPError(f"n_steps must be >= 1, got {n_steps}.")

        params = ff_params or {}
        ff = force_function.lower().strip()
        box = box_size
        samp = max(trajectory_sample_every, 1)

        # Deep copy initial state
        pos = [[p[i] for i in range(3)] for p in positions]
        vel = [[v[i] for i in range(3)] for v in velocities]

        # Initial acceleration
        forces, pe = self._compute_forces(pos, ff, box, params)
        acc = [[forces[i][j] / masses[i] for j in range(3)] for i in range(n)]

        trajectory = []
        energies = []
        ke_init = self._ke(vel, masses)
        e_init = ke_init + pe
        max_dev = 0.0

        for step in range(n_steps):
            # Step 1: v(t + dt/2) = v(t) + 0.5 * a(t) * dt
            for i in range(n):
                for j in range(3):
                    vel[i][j] += 0.5 * acc[i][j] * dt

            # Step 2: r(t + dt) = r(t) + v(t + dt/2) * dt
            for i in range(n):
                for j in range(3):
                    pos[i][j] += vel[i][j] * dt
                    while pos[i][j] < 0:
                        pos[i][j] += box
                    while pos[i][j] >= box:
                        pos[i][j] -= box

            # Step 3: Compute a(t + dt) from new positions
            forces, pe = self._compute_forces(pos, ff, box, params)
            acc_new = [[forces[i][j] / masses[i] for j in range(3)] for i in range(n)]

            # Step 4: v(t + dt) = v(t + dt/2) + 0.5 * a(t + dt) * dt
            for i in range(n):
                for j in range(3):
                    vel[i][j] += 0.5 * acc_new[i][j] * dt

            acc = acc_new

            # Sample
            if step % samp == 0 or step == n_steps - 1:
                traj_frame = [[pos[i][j] for j in range(3)] for i in range(n)]
                trajectory.append(traj_frame)
                ke_now = self._ke(vel, masses)
                e_tot = ke_now + pe
                dev = abs(e_tot - e_init) / abs(e_init) if abs(e_init) > 1e-30 else 0.0
                max_dev = max(max_dev, dev)
                energies.append([step, round(ke_now, 12), round(pe, 12), round(e_tot, 12)])

        ke_final = self._ke(vel, masses)
        e_final = ke_final + pe
        drift = (e_final - e_init) / abs(e_init) if abs(e_init) > 1e-30 else 0.0
        temp = self._temp(ke_final, n)
        times = [step * dt for step in range(n_steps) if step % samp == 0 or step == n_steps - 1]

        logger.info(
            f"VelocityVerlet ({ff}): {n_steps} steps, {n} atoms, "
            f"E_init={e_init:.6g}, E_final={e_final:.6g}, drift={drift:.2e}, max_dev={max_dev:.2e}"
        )

        return {
            "trajectory": trajectory,
            "final_positions": [[pos[i][j] for j in range(3)] for i in range(n)],
            "final_velocities": [[vel[i][j] for j in range(3)] for i in range(n)],
            "energies_over_time": energies,
            "kinetic_energy_final": round(ke_final, 12),
            "potential_energy_final": round(pe, 12),
            "total_energy_final": round(e_final, 12),
            "temperature": round(temp, 8),
            "energy_drift": round(drift, 14),
            "max_energy_deviation": round(max_dev, 14),
            "n_atoms": n,
            "n_steps_completed": n_steps,
            "force_field_type": ff,
            "time_array": [round(t, 8) for t in times],
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
