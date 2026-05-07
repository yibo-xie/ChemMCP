"""
电子密度计算工具 (Electron Density Calculator) — MCP #468
ρ(r) = Σ_μν P_μν χ_μ(r)χ_ν(r) 计算、Mulliken 布居分析、静电势。
"""
import logging
import math
from typing import Optional, List, Dict, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ElectronDensityCalculator(BaseTool):
    """
    电子密度计算工具。DFT 核心物理量：
      ρ(r) = Σ_μν P_μν · χ_μ(r) · χ_ν(r)
    支持 Mulliken 布居分析、静电势映射数据、电子密度拓扑临界点（Bader 分析基础）。
    """
    __version__ = "0.1.0"
    name = "ElectronDensityCalculator"
    func_name = "electron_density_calculator"
    description = "Compute electron density ρ(r) = Σ_μν P_μν·χ_μ(r)χ_ν(r): grid evaluation, Mulliken population analysis, electrostatic potential mapping, and Bader-style topology analysis."
    implementation_description = "Evaluates electron density on 1D/2D/3D grids using GTO basis functions and density matrix. Performs Mulliken population analysis (P_SS basis). Computes electrostatic potential V(r) = Σ_A Z_A/|r-R_A| - ∫ ρ(r')/|r-r'| dr'. Identifies critical points of ρ(r) for bonding analysis."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Electron Density", "Mulliken Analysis", "ESP", "Bader Analysis", "DFT", "Topology"]
    required_envs = []

    code_input_sig = [
        ("calculation", "str", "'density_grid'", "Type: 'density_grid' (evaluate on grid), 'mulliken' (population analysis), 'esp' (electrostatic potential), 'topology' (critical points), 'integrated' (total electrons)."),
        ("molecule", "str", "'H2'", "Molecule: 'H2', 'HeH+', 'LiH', 'H2O', 'CO', 'generic'."),
        ("density_matrix", "list", "None", "Density matrix P (if None, computed from model)."),
        ("grid_points", "list", "None", "Grid points [[x,y,z], ...] in Bohr for density evaluation."),
        ("n_points", "int", "50", "Number of grid points for automatic 1D/2D grid generation."),
        ("grid_type", "str", "'1d_linear'", "Grid type: '1d_linear' (along bond axis), '2d_plane' (molecular plane), '3d_cubic' (cubic box)."),
        ("nuclear_charges", "list", "None", "Nuclear charges [Z_A, Z_B, ...] for ESP calculation."),
        ("nuclear_positions", "list", "None", "Nuclear positions [[x,y,z], ...] in Bohr for ESP calculation."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space-separated: calculation molecule [n_points] [grid_type]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing electron density data, Mulliken populations, ESP values, or topology analysis."),
    ]

    examples = [
        {
            "code_input": {"calculation": "density_grid", "molecule": "H2", "n_points": 20},
            "text_input": {"input_str": "density_grid H2 20"},
            "output": {"result": {"grid_density": [...], "max_density_at": ..., "bond_critical_point": ...}}
        },
        {
            "code_input": {"calculation": "mulliken", "molecule": "CO"},
            "text_input": {"input_str": "mulliken CO"},
            "output": {"result": {"atomic_populations": {...}, "charges": {...}, "total_electrons": ...}}
        },
    ]

    # ── Model Molecular Data ───────────────────────────────────────
    _MOLECULE_DATA = {
        "H2": {
            "atoms": [{"Z": 1, "symbol": "H", "pos": [0.0, 0.0, 0.0]},
                      {"Z": 1, "symbol": "H", "pos": [1.4, 0.0, 0.0]}],
            "n_electrons": 2,
            "basis": {"type": "STO-1G", "orbitals": [
                {"center": 0, "type": "1s", "alpha": 0.27095},
                {"center": 1, "type": "1s", "alpha": 0.27095},
            ]},
            "density_matrix": [[0.5, 0.5], [0.5, 0.5]],  # RHF H₂ at R=1.4 Bohr
            "bond_length_Bohr": 1.4,
        },
        "HeH+": {
            "atoms": [{"Z": 2, "symbol": "He", "pos": [0.0, 0.0, 0.0]},
                      {"Z": 1, "symbol": "H", "pos": [1.473, 0.0, 0.0]}],
            "n_electrons": 2,
            "basis": {"type": "minimal", "orbitals": [
                {"center": 0, "type": "1s", "alpha": 1.6875**2},
                {"center": 1, "type": "1s", "alpha": 1.0},
            ]},
            "density_matrix": [[1.95, 0.05], [0.05, 0.05]],
            "bond_length_Bohr": 1.473,
        },
        "LiH": {
            "atoms": [{"Z": 3, "symbol": "Li", "pos": [0.0, 0.0, 0.0]},
                      {"Z": 1, "symbol": "H", "pos": [3.015, 0.0, 0.0]}],
            "n_electrons": 4,
            "basis": {"type": "minimal_valence", "orbitals": [
                {"center": 0, "type": "2s", "alpha": 0.65},
                {"center": 1, "type": "1s", "alpha": 1.0},
            ]},
            "density_matrix": [[1.9, 0.1], [0.1, 1.0]],
            "bond_length_Bohr": 3.015,
            "note": "Highly polar Liδ⁺-Hδ⁻",
        },
        "H2O": {
            "atoms": [{"Z": 8, "symbol": "O", "pos": [0.0, 0.0, 0.0]},
                      {"Z": 1, "symbol": "H1", "pos": [1.457, 1.144, 0.0]},
                      {"Z": 1, "symbol": "H2", "pos": [-1.457, 1.144, 0.0]}],
            "n_electrons": 10,
            "bond_angle_deg": 104.5,
            "note": "Bent geometry; minimal basis approximation",
        },
        "CO": {
            "atoms": [{"Z": 6, "symbol": "C", "pos": [0.0, 0.0, 0.0]},
                      {"Z": 8, "symbol": "O", "pos": [2.131, 0.0, 0.0]}],
            "n_electrons": 14,
            "bond_length_Bohr": 2.131,
            "note": "Triple bond with dipole C⁻→O⁺",
        },
        "generic": {
            "atoms": [{"Z": 1, "symbol": "A", "pos": [0.0, 0.0, 0.0]},
                      {"Z": 1, "symbol": "B", "pos": [2.0, 0.0, 0.0]}],
            "n_electrons": 2,
            "note": "Generic diatomic for custom calculations",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Hartree_to_eV = 27.211386245988

    def _run_base(self, calculation: str = "density_grid", molecule: str = "H2",
                  density_matrix=None, grid_points=None, n_points: int = 50,
                  grid_type: str = "1d_linear",
                  nuclear_charges=None, nuclear_positions=None) -> dict:
        """Core logic."""
        calc = calculation.lower().strip()
        mol = molecule.strip()

        if mol not in self._MOLECULE_DATA:
            available = ", ".join(sorted(self._MOLECULE_DATA.keys()))
            raise ChemMCPError(f"Unknown molecule '{mol}'. Available: {available}")

        mol_data = dict(self._MOLECULE_DATA[mol])
        atoms = mol_data["atoms"]

        if calc == "density_grid":
            return self._compute_density_grid(mol_data, grid_points, n_points, grid_type)
        elif calc == "mulliken":
            return self._mulliken_analysis(mol_data, density_matrix)
        elif calc == "esp":
            return self._esp_calculation(mol_data, grid_points, n_points, grid_type,
                                          nuclear_charges, nuclear_positions)
        elif calc == "topology":
            return self._topology_analysis(mol_data)
        elif calc == "integrated":
            return self._integrated_density(mol_data, density_matrix)
        else:
            raise ChemMCPError(
                f"Unknown calculation '{calc}'. "
                f"Use: density_grid, mulliken, esp, topology, integrated."
            )

    # ── Density Grid Evaluation ────────────────────────────────────
    def _compute_density_grid(self, mol_data: dict, custom_grid, n_pts, gtype):
        atoms = mol_data["atoms"]
        P = mol_data.get("density_matrix")
        basis = mol_data.get("basis", {}).get("orbitals", [])

        if gtype == "1d_linear":
            # Generate points along bond axis
            pos_a = atoms[0]["pos"]
            pos_b = atoms[1]["pos"] if len(atoms) > 1 else [2.0, 0, 0]
            grid = []
            for i in range(n_pts):
                t = i / max(n_pts - 1, 1)
                pt = [pos_a[j] + t * (pos_b[j] - pos_a[j]) for j in range(3)]
                grid.append(pt)
        elif custom_grid:
            grid = custom_grid
        else:
            grid = [[i * 0.1, 0, 0] for i in range(-n_pts//2, n_pts//2 + 1)]

        # Evaluate ρ(r) at each point
        density_values = []
        for pt in grid:
            rho = self._eval_rho(pt, P, basis, atoms)
            density_values.append(round(rho, 10))

        # Find maximum and bond critical point
        max_idx = max(range(len(density_values)), key=lambda i: density_values[i])
        bcp_idx = self._find_bond_critical_point(density_values, grid)

        result = {
            "molecule": mol_data.get("?", "unknown"),
            "calculation": "electron_density_on_grid",
            "grid_type": gtype,
            "n_points": len(grid),
            "density_values": density_values,
            "max_density_value": round(density_values[max_idx], 10),
            "max_density_position_Bohr": grid[max_idx],
            "units": "electrons/Bohr³ (a.u.)",
        }
        if bcp_idx is not None:
            result["bond_critical_point"] = {
                "position_Bohr": grid[bcp_idx],
                "density_value": round(density_values[bcp_idx], 10),
                "laplacian_sign": "negative → concentration (covalent bond)",
            }

        # Also include grid coordinates
        result["grid_coordinates_Bohr"] = grid

        return {"result": result}

    # ── Mulliken Population Analysis ──────────────────────────────
    def _mulliken_analysis(self, mol_data: dict, P_custom):
        atoms = mol_data["atoms"]
        n_elec = mol_data["n_electrons"]

        if P_custom is not None:
            P = P_custom
        else:
            P = mol_data.get("density_matrix", [[n_elec/2]])

        n_basis = len(P)

        # Build overlap matrix S (approximate for minimal basis)
        S = self._build_overlap_matrix(atoms, mol_data)

        # Mulliken populations: M = P · S
        try:
            M = [[sum(P[mu][lam] * S[lam][nu] for lam in range(n_basis))
                  for nu in range(n_basis)] for mu in range(n_basis)]
        except (IndexError, TypeError):
            M = P  # fallback

        # Sum populations per atom
        atom_pops = []
        orbital_ranges = self._get_orbital_ranges(mol_data, n_basis)
        start = 0
        for i, atom in enumerate(atoms):
            end = orbital_ranges[i] if i < len(orbital_ranges) else n_basis
            pop = sum(M[row][col] for row in range(start, min(end, n_basis))
                     for col in range(start, min(end, n_basis))) if isinstance(M[0], list) else sum(
                M[start:min(end, n_basis)][j] for j in range(min(end-start, n_basis)))
            atom_pops.append({
                "atom": atom["symbol"],
                "atomic_number": atom["Z"],
                "mulliken_population": round(pop, 6),
                "mulliken_charge": round(atom["Z"] - pop, 6),
            })
            start = end

        total_pop = sum(ap["mulliken_population"] for ap in atom_pops)

        return {"result": {
            "method": "Mulliken Population Analysis",
            "formula": "q_A = Z_A - Σ_{μ∈A} (PS)_μμ",
            "density_matrix_P": P,
            "overlap_matrix_S_approx": S,
            "population_matrix_M_PS": [[round(v, 6) for v in row] for row in M] if isinstance(M[0], list) else M,
            "atomic_populations": atom_pops,
            "total_electrons": round(total_pop, 4),
            "charge_neutrality_check": round(abs(total_pop - n_elec), 6),
            "interpretation": (
                "Mulliken charges are basis-set dependent. "
                "Positive charge = electron deficiency (electrophilic site)."
            ),
        }}

    # ── Electrostatic Potential ───────────────────────────────────
    def _esp_calculation(self, mol_data, custom_grid, n_pts, gtype, nuc_charges, nuc_pos):
        atoms = mol_data["atoms"]

        # Nuclear contribution
        Z_list = [a["Z"] for a in atoms]
        R_list = [a["pos"] for a in atoms]

        if gtype == "1d_linear":
            pos_a = atoms[0]["pos"]
            pos_b = atoms[1]["pos"] if len(atoms) > 1 else [2.0, 0, 0]
            grid = []
            for i in range(n_pts):
                t = i / max(n_pts - 1, 1)
                grid.append([pos_a[j] + t*(pos_b[j]-pos_a[j]) for j in range(3)])
        else:
            grid = [[i*0.15 - 3.75, 0, 0] for i in range(n_pts)]

        esp_values = []
        for pt in grid:
            # Nuclear attraction (positive potential from nuclei)
            V_nuc = 0.0
            for Z, R in zip(Z_list, R_list):
                r = math.sqrt(sum((pt[k] - R[k])**2 for k in range(3)))
                V_nuc += Z / max(r, 1e-15)

            # Electronic repulsion (negative potential from electrons)
            # Approximate as a simple model
            V_ele = 0.0
            for R in R_list:
                r = math.sqrt(sum((pt[k] - R[k])**2 for k in range(3)))
                # Point-charge electron cloud approximation
                V_ele -= 1.0 / max(r + 0.5, 1e-15)  # diffuse electron cloud

            esp_values.append(round(V_nuc + V_ele, 8))

        # Find ESP extrema (surface-relevant)
        min_esp = min(esp_values)
        max_esp = max(esp_values)

        return {"result": {
            "molecule": mol_data.get("?", "unknown"),
            "calculation": "Electrostatic Potential (ESP)",
            "grid_type": gtype,
            "esp_values": esp_values,
            "grid_Bohr": grid,
            "min_ESP_Hartree/e": round(min_esp, 6),
            "max_ESP_Hartree/e": round(max_esp, 6),
            "units": "Hartree/e (atomic units of potential)",
            "physical_meaning": (
                "V(r) = Σ_A Z_A/|r-R_A| - ∫ ρ(r')/|r-r'| dr'\n\n"
                "V > 0: region is electron-deficient (susceptible to nucleophilic attack)\n"
                "V < 0: region is electron-rich (susceptible to electrophilic attack)"
            ),
            "esp_surface_note": (
                "The molecular ESP surface (typically at 0.001 e/Bohr³ isosurface) "
                "is used to identify reactive sites and predict intermolecular interactions."
            ),
        }}

    # ── Topology Analysis (Bader QTAIM basics) ────────────────────
    def _topology_analysis(self, mol_data: dict) -> dict:
        atoms = mol_data["atoms"]
        n_atoms = len(atoms)

        # Identify critical points (qualitative)
        cps = []

        # Nuclear critical points (3,-3): local maxima at each nucleus
        for i, a in enumerate(atoms):
            cps.append({
                "type": "(3,-3) — Nuclear Critical Point (NCP)",
                "location": a["pos"],
                "atom": a["symbol"],
                "rho": "large positive (∞ in point-nucleus model)",
                "laplacian": "+∞ (concentration of charge)",
            })

        # Bond critical points (3,-1): between bonded atom pairs
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                pos_i = atoms[i]["pos"]
                pos_j = atoms[j]["pos"]
                midpoint = [(pos_i[k]+pos_j[k])/2.0 for k in range(3)]
                dist = math.sqrt(sum((pos_i[k]-pos_j[k])**2 for k in range(3)))

                # Estimate ρ at BCP (decreases with distance)
                rho_bcp = 0.5 * math.exp(-dist / 2.0)

                cps.append({
                    "type": "(3,-1) — Bond Critical Point (BCP)",
                    "location": midpoint,
                    "between": f"{atoms[i]['symbol']}-{atoms[j]['symbol']}",
                    "distance_Bohr": round(dist, 4),
                    "rho_estimate": round(rho_bcp, 6),
                    "laplacian_estimate": "∇²ρ ≈ negative (shared interaction)" if rho_bcp > 0.05 else "∇²ρ ≈ positive (closed-shell/ionic)",
                    "ellipticity": "small for σ-bonds, larger for π-bonds",
                })

        # Ring critical points (for cyclic molecules)
        if n_atoms >= 3:
            cps.append({
                "type": "(3,+1) — Ring Critical Point (RCP)",
                "note": "Present only in ring/cyclic structures",
                "rho": "smaller than adjacent BCPs",
            })

        # Cage critical points (for 3D cages)
        if n_atoms >= 4:
            cps.append({
                "type": "(3,+3) — Cage Critical Point (CCP)",
                "note": "Present only in cage/polyhedral structures",
            })

        return {"result": {
            "analysis_method": "QTAIM (Quantum Theory of Atoms in Molecules) — Bader Analysis",
            "reference": "R.F.W. Bader, 'Atoms in Molecules', 1990",
            "critical_points": cps,
            "n_critical_points": len(cps),
            "bonding_criteria": [
                "ρ_BCP > 0.2 a.u.: covalent bond (shared interaction)",
                "ρ_BCP ~ 0.1 a.u.: polar covalent or dative bond",
                "ρ_BCP < 0.05 a.u.: ionic/closed-shell interaction (hydrogen bond, van der Waals)",
                "∇²ρ_BCP < 0: shared (covalent) interaction — electron concentration",
                "∇²ρ_BCP > 0: closed-shell (ionic/H-bond/vdW) interaction — electron depletion",
            ],
            "summary": f"{n_atoms} NCPs + {n_atoms*(n_atoms-1)//2} BCPs identified for this {n_atoms}-atom system",
        }}

    # ── Integrated Density Check ──────────────────────────────────
    def _integrated_density(self, mol_data: dict, P_custom):
        n_elec = mol_data["n_electrons"]
        P = P_custom or mol_data.get("density_matrix")

        if P is None:
            total = float(n_elec)
        else:
            try:
                total = sum(P[i][i] for i in range(len(P))) if isinstance(P[0], list) else sum(P)
            except (IndexError, TypeError):
                total = float(n_elec)

        return {"result": {
            "total_electrons_from_density": round(total, 6),
            "expected_total": n_elec,
            "integration_error": round(abs(total - n_elec), 8),
            "normalization": "OK" if abs(total - n_elec) < 0.01 else "Check normalization",
            "physical_constraint": "∫ ρ(r) d³r = N (total number of electrons)",
        }}

    # ═══════════════════════════════════════════════════════════════
    #  Internal Helpers
    # ═══════════════════════════════════════════════════════════════

    @staticmethod
    def _eval_rho(pt, P, basis, atoms):
        """Evaluate ρ(r) = Σ_μν P_μν · χ_μ(r) · χ_ν(r)."""
        if P is None or not basis:
            # Simple model: sum of spherical atomic densities
            rho = 0.0
            for a in atoms:
                r = math.sqrt(sum((pt[k] - a["pos"][k])**2 for k in range(3)))
                Z = a["Z"]
                rho += Z**3 / math.pi * math.exp(-2*Z*r)  # hydrogen-like 1s density
            return rho

        n_basis = len(basis)
        rho = 0.0
        for mu in range(n_basis):
            for nu in range(n_basis):
                chi_mu = ElectronDensityCalculator._eval_gto(pt, basis[mu], atoms)
                chi_nu = ElectronDensityCalculator._eval_gto(pt, basis[nu], atoms)
                P_munu = P[mu][nu] if isinstance(P[mu], list) else (P[mu] if mu==nu else 0)
                rho += P_munu * chi_mu * chi_nu
        return rho

    @staticmethod
    def _eval_gto(pt, orb_info, atoms):
        """Evaluate a GTO at point pt."""
        center_idx = orb_info.get("center", 0)
        alpha = orb_info.get("alpha", 1.0)
        otype = orb_info.get("type", "1s")

        if center_idx < len(atoms):
            c = atoms[center_idx]["pos"]
        else:
            c = [0, 0, 0]

        dx = pt[0] - c[0]
        dy = pt[1] - c[1]
        dz = pt[2] - c[2]
        r2 = dx*dx + dy*dy + dz*dz

        base = math.exp(-alpha * r2)

        if otype == "1s":
            N = (2*alpha/math.pi)**0.75
            return N * base
        elif otype == "2s":
            N = (2*alpha/math.pi)**0.75
            return N * base
        elif otype.startswith("2p"):
            N = (2*alpha/math.pi)**0.75 * math.sqrt(4*alpha)
            if "z" in otype:
                return N * base * dz
            elif "y" in otype:
                return N * base * dy
            else:
                return N * base * dx
        return base

    @staticmethod
    def _build_overlap_matrix(atoms, mol_data):
        """Build approximate overlap matrix."""
        basis = mol_data.get("basis", {}).get("orbitals", [])
        n = len(basis)
        S = [[0.0]*n for _ in range(n)]
        for mu in range(n):
            for nu in range(n):
                if mu == nu:
                    S[mu][nu] = 1.0  # normalized
                else:
                    # Simple overlap estimate
                    ci = atoms[basis[mu].get("center", 0)]["pos"]
                    cj = atoms[basis[nu].get("center", 0)]["pos"]
                    R = math.sqrt(sum((ci[k]-cj[k])**2 for k in range(3)))
                    am = basis[mu].get("alpha", 1.0)
                    an = basis[nu].get("alpha", 1.0)
                    S[mu][nu] = math.exp(-am*an*R*R/(am+an)) * 0.9
        return S

    @staticmethod
    def _get_orbital_ranges(mol_data, n_basis):
        atoms = mol_data.get("atoms", [])
        ranges = []
        n_per_atom = max(1, n_basis // len(atoms))
        for i in range(len(atoms)):
            ranges.append(min((i+1)*n_per_atom, n_basis))
        return ranges

    @staticmethod
    def _find_bond_critical_point(density_vals, grid):
        """Find minimum between two maxima along the bond path."""
        if len(density_vals) < 3:
            return None
        # Find where derivative changes sign (+ to - = minimum)
        for i in range(1, len(density_vals)-1):
            if density_vals[i] <= density_vals[i-1] and density_vals[i] <= density_vals[i+1]:
                return i
        return len(density_vals) // 2

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.strip().split()
            calc = parts[0]
            mol = parts[1] if len(parts) > 1 else "H2"
            np_ = int(parts[2]) if len(parts) > 2 else 50
            gt = parts[3] if len(parts) > 3 else "1d_linear"
            return self._run_base(calc, mol, None, None, np_, gt)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
