import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


def _dr_pbc(r1, r2, box):
    """Minimum image displacement vector with PBC."""
    dr = [r2[i] - r1[i] for i in range(3)]
    hb = box / 2.0
    for i in range(3):
        if dr[i] > hb:
            dr[i] -= box
        elif dr[i] < -hb:
            dr[i] += box
    return dr


def _abs_vec(v):
    return math.sqrt(sum(x * x for x in v))


@ChemMCPManager.register_tool
class ForceFieldCalculator(BaseTool):
    """
    力场计算工具。
    分子力学能量与力的计算，支持 Lennard-Jones、Morse、谐振子、Buckingham 势函数。
    """
    __version__ = "0.1.0"
    name = "ForceFieldCalculator"
    func_name = "calculate_force_field"
    description = "Force field calculation for molecular mechanics: potential energy and atomic forces using LJ, Morse, harmonic, or Buckingham potentials."
    implementation_description = "Computes pairwise molecular mechanics interactions: Lennard-Jones (12-6), Morse potential, harmonic bond stretch, and Buckingham (exp-6). Returns total PE, per-atom forces (negative gradient), virial tensor components, and stress tensor diagonal."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Force Field", "Molecular Mechanics", "Potential Energy", "Forces", "Computational Chemistry", "LJ", "Morse"]
    required_envs = []

    code_input_sig = [
        ("positions", "list", "N/A", "List of [x,y,z] coordinates for each atom."),
        ("atom_types", "list", "N/A", "List of atom type identifiers (int or str) for each atom."),
        ("force_field_type", "str", "lj", "Force field: 'lj', 'morse', 'harmonic', 'buckingham'."),
        ("parameters", "dict", "N/A", "Force field parameters dict. See description for required keys per type."),
        ("box_size", "float", "None", "Cubic box size for PBC (None = no PBC)."),
        ("cutoff", "float", "None", "Interaction cutoff distance (None = auto from ff_type)."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all force field parameters."),
    ]

    output_sig = [
        ("total_potential_energy", "float", "Total potential energy of the system."),
        ("forces", "list", "List of [fx,fy,fz] force vectors per atom (F = -∇V)."),
        ("per_atom_energies", "list", "Energy contribution per atom (half of pair energy assigned to each)."),
        ("virial", "float", "Virial Σ r_ij · F_ij (for pressure calculation)."),
        ("stress_tensor_diagonal", "list", "[σ_xx, σ_yy, σ_zz] diagonal stress components."),
        ("n_atoms", "int", "Number of atoms."),
        ("n_pairs", "int", "Number of interacting pairs within cutoff."),
        ("force_field_type", "str", "Force field used."),
        ("pairwise_details", "list", "Optional list of (i,j,r,V_pair,F_mag) for each pair."),
        ("parameters_used", "dict", "Parameters actually used in the computation."),
    ]

    examples = [
        {
            "code_input": {
                "positions": [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
                "atom_types": [0, 1],
                "force_field_type": "lj",
                "parameters": {"epsilon": 1.0, "sigma": 1.0},
            },
            "text_input": {"params_str": '{"positions":[[0,0,0],[1.2,0,0]],"atom_types":[0,1],"force_field_type":"lj","parameters":{"epsilon":1,"sigma":1}}'},
            "output": {
                "total_potential_energy": 0.0,
                "n_atoms": 2,
                "n_pairs": 1,
            },
        },
        {
            "code_input": {
                "positions": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
                "atom_types": [0, 1, 2],
                "force_field_type": "harmonic",
                "parameters": {"k": 200.0, "r0": 1.0},
            },
            "text_input": {"params_str": '{"positions":[[0,0,0],[1,0,0],[2.5,0,0]],"atom_types":[0,1,2],"force_field_type":"harmonic","parameters":{"k":200,"r0":1}}'},
            "output": {
                "n_pairs": 3,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- Individual potentials ----

    @staticmethod
    def _compute_lj(positions, params, box, cutoff):
        """Lennard-Jones: V=4ε[(σ/r)^12-(σ/r)^6], F=-dV/dr."""
        eps = params.get("epsilon", 1.0)
        sig = params.get("sigma", 1.0)
        rc = cutoff or min(3.0 * sig, box / 2.0 if box else 3.0 * sig)
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe_total = 0.0
        virial = 0.0
        pair_details = []
        n_pairs = 0

        for i in range(n):
            for j in range(i + 1, n):
                dr = _dr_pbc(positions[i], positions[j], box) if box else [positions[j][c] - positions[i][c] for c in range(3)]
                r = _abs_vec(dr)
                if r < 1e-10 or r > rc:
                    continue

                sr = sig / r
                sr6 = sr ** 6
                sr12 = sr6 ** 2
                V = 4.0 * eps * (sr12 - sr6)
                # dV/dr = 4ε[-12σ^12/r^13 + 6σ^7/r^7]
                # F_mag = -dV/dr (along r direction)
                f_mag = 24.0 * eps * (2.0 * sr12 - sr6) / r

                fk = [f_mag * dr[c] / r for c in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]

                virial += f_mag * r
                pe_total += V
                n_pairs += 1
                pair_details.append((i, j, round(r, 8), round(V, 12), round(f_mag, 10)))

        return pe_total, forces, virial, n_pairs, pair_details

    @staticmethod
    def _compute_morse(positions, params, box, cutoff):
        """Morse: V=D_e*[1-exp(-a(r-r0))]^2."""
        De = params.get("De", 1.0)
        a_morse = params.get("a", 2.0)
        r0 = params.get("r0", 1.0)
        rc = cutoff or (r0 + 5.0 / a_morse)  # ~5 decay lengths
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe_total = 0.0
        virial = 0.0
        pair_details = []
        n_pairs = 0

        for i in range(n):
            for j in range(i + 1, n):
                dr = _dr_pbc(positions[i], positions[j], box) if box else [positions[j][c] - positions[i][c] for c in range(3)]
                r = _abs_vec(dr)
                if r < 1e-10 or r > rc:
                    continue

                exp_ar = math.exp(-a_morse * (r - r0))
                V = De * (1.0 - exp_ar) ** 2
                # dV/dr = 2*De*(1-exp(-a(r-r0)))*a*exp(-a(r-r0))
                #       = 2*a*De*exp(-a(r-r0))*(1-exp(-a(r-r0)))
                dVdr = 2.0 * a_morse * De * exp_ar * (1.0 - exp_ar)
                f_mag = -dVdr  # F = -dV/dr

                fk = [f_mag * dr[c] / r for c in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]

                virial += (-dVdr) * r
                pe_total += V
                n_pairs += 1
                pair_details.append((i, j, round(r, 8), round(V, 12), round(f_mag, 10)))

        return pe_total, forces, virial, n_pairs, pair_details

    @staticmethod
    def _compute_harmonic(positions, params, box, cutoff):
        """Harmonic bond: V=0.5*k*(r-r0)²."""
        k = params.get("k", 100.0)
        r0 = params.get("r0", 1.0)
        rc = cutoff or (r0 + 3.0 * math.sqrt(params.get("k", 100.0) / 10.0)) if False else None
        # For harmonic, no natural cutoff; use large default
        rc = rc or 1e10
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe_total = 0.0
        virial = 0.0
        pair_details = []
        n_pairs = 0

        for i in range(n):
            for j in range(i + 1, n):
                dr = _dr_pbc(positions[i], positions[j], box) if box else [positions[j][c] - positions[i][c] for c in range(3)]
                r = _abs_vec(dr)
                if r < 1e-10 or r > rc:
                    continue

                dr_r0 = r - r0
                V = 0.5 * k * dr_r0 ** 2
                f_mag = -k * dr_r0  # F = -dV/dr = -k*(r-r0)

                fk = [f_mag * dr[c] / r for c in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]

                virial += k * dr_r0 * r
                pe_total += V
                n_pairs += 1
                pair_details.append((i, j, round(r, 8), round(V, 12), round(f_mag, 10)))

        return pe_total, forces, virial, n_pairs, pair_details

    @staticmethod
    def _compute_buckingham(positions, params, box, cutoff):
        """Buckingham: V=A*exp(-Br) - C/r^6."""
        A = params.get("A", 1000.0)
        B = params.get("B", 3.6)
        C = params.get("C", 1.0)
        rc = cutoff or 10.0
        n = len(positions)
        forces = [[0.0, 0.0, 0.0] for _ in range(n)]
        pe_total = 0.0
        virial = 0.0
        pair_details = []
        n_pairs = 0

        for i in range(n):
            for j in range(i + 1, n):
                dr = _dr_pbc(positions[i], positions[j], box) if box else [positions[j][c] - positions[i][c] for c in range(3)]
                r = _abs_vec(dr)
                if r < 0.5 or r > rc:  # Buckingham diverges at small r
                    continue

                exp_br = math.exp(-B * r)
                V = A * exp_br - C / (r ** 6)
                # dV/dr = -A*B*exp(-Br) + 6*C/r^7
                dVdr = -A * B * exp_br + 6.0 * C / (r ** 7)
                f_mag = -dVdr

                fk = [f_mag * dr[c] / r for c in range(3)]
                for c in range(3):
                    forces[i][c] += fk[c]
                    forces[j][c] -= fk[c]

                virial += (-dVdr) * r
                pe_total += V
                n_pairs += 1
                pair_details.append((i, j, round(r, 8), round(V, 12), round(f_mag, 10)))

        return pe_total, forces, virial, n_pairs, pair_details

    def _run_base(
        self,
        positions: List[List[float]],
        atom_types: List,
        force_field_type: str = "lj",
        parameters: dict = None,
        box_size: Optional[float] = None,
        cutoff: Optional[float] = None,
    ) -> dict:
        """Core logic: compute force field energy and forces."""
        n = len(positions)
        if n == 0:
            raise ChemMCPError("No atoms provided.")
        if len(atom_types) != n:
            raise ChemMCPError(f"atom_types length ({len(atom_types)}) != positions length ({n}).")
        if not parameters:
            raise ChemMCPError("parameters dict is required.")

        ff = force_field_type.lower().strip()
        box = box_size

        if ff == "lj":
            pe, forces, virial, npairs, pdetails = self._compute_lj(positions, parameters, box, cutoff)
        elif ff == "morse":
            pe, forces, virial, npairs, pdetails = self._compute_morse(positions, parameters, box, cutoff)
        elif ff == "harmonic":
            pe, forces, virial, npairs, pdetails = self._compute_harmonic(positions, parameters, box, cutoff)
        elif ff == "buckingham":
            pe, forces, virial, npairs, pdetails = self._compute_buckingham(positions, parameters, box, cutoff)
        else:
            raise ChemMCPError(f"Unknown force_field_type '{ff}'. Use: lj, morse, harmonic, buckingham.")

        # Per-atom energies (half to each atom)
        per_atom_e = [0.0] * n
        for i, j, r, V, fm in pdetails:
            per_atom_e[i] += V / 2.0
            per_atom_e[j] += V / 2.0

        # Stress tensor diagonal (diagonal approximation): σ_αα = -(1/V)(PE_αα + W_αα)
        # Simplified: use virial / volume
        vol = box ** 3 if box else 1.0
        stress_diag = [-virial / vol / 3.0] * 3 if vol > 0 else [0.0, 0.0, 0.0]

        logger.info(
            f"ForceFieldCalculator ({ff}): {n} atoms, {npairs} pairs, "
            f"PE={pe:.6g}, |F_max|={max(_abs_vec(f) for f in forces):.4g}"
        )

        return {
            "total_potential_energy": round(pe, 14),
            "forces": [[round(f[c], 12) for c in range(3)] for f in forces],
            "per_atom_energies": [round(e, 14) for e in per_atom_e],
            "virial": round(virial, 12),
            "stress_tensor_diagonal": [round(s, 12) for s in stress_diag],
            "n_atoms": n,
            "n_pairs": npairs,
            "force_field_type": ff,
            "pairwise_details": pdetails[:50],  # limit output size
            "parameters_used": dict(parameters),
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
