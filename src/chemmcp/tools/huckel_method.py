import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HuckelMethod(BaseTool):
    """
    休克尔分子轨道法计算π电子体系。
    
    构建休克尔矩阵: 对角元 = α (库仑积分), 邻接非对角元 = β (共振积分)
    求解本征值 ε_k 和本征向量 c_ki
    
    π电子总能量: E_π = Σ n_k ε_k
    离域能 = E_π - E_localized (定域化参考能量)
    """
    __version__ = "0.1.0"
    name = "HuckelMethod"
    func_name = "huckel_method"
    description = "Hückel MO theory for conjugated π-electron systems: linear/cyclic polyenes, aromatic compounds, and radicals."
    implementation_description = "Builds Hückel Hamiltonian matrix (α on diagonal, β for bonded neighbors). Solves for orbital energies ε_k and LCAO coefficients c_ki. Computes total π energy, delocalization/resonance energy, charge densities, bond orders, frontier orbitals (HOMO/LUMO), and spectral transition estimates."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Chemistry", "Hückel Method", "Pi Electrons", "Conjugated Systems", "Aromaticity"]
    required_envs = []

    code_input_sig = [
        ("molecule", "str", "N/A", "Molecule: 'ethene', 'butadiene', 'hexatriene', 'benzene', 'naphthalene', 'cyclobutadiene', 'allyl', 'pentadienyl', 'cyclopentadienyl'."),
        ("n_carbons", "int", "None", "Number of carbon atoms for linear polyene (auto-set from molecule name if None)."),
        ("topology", "str", "None", "'linear' or 'cyclic' (auto-determined from molecule if None)."),
        ("charge", "int", "0", "Net charge (affects electron count)."),
        ("ionization_state", "str", "neutral", "'neutral', 'cation', or 'anion' (for radical species)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'molecule [n_carbons] [topology] [charge]'. Example: 'benzene cyclic 0'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with orbital energies, MO coefficients, electron configuration, total π energy, delocalization energy, charge densities, bond orders, frontier orbital info."),
    ]

    examples = [
        {
            "code_input": {"molecule": "butadiene", "n_carbons": None, "topology": None, "charge": 0, "ionization_state": "neutral"},
            "text_input": {"input_params": "butadiene"},
            "output": {"result": {"n_carbons": 4, "total_pi_energy_alpha_beta": "4α + 4.472β", "delocalization_energy_beta": 0.472}},
        },
        {
            "code_input": {"molecule": "benzene", "n_carbons": None, "topology": None, "charge": 0, "ionization_state": "neutral"},
            "text_input": {"input_params": "benzene"},
            "output": {"result": {"n_carbons": 6, "resonance_energy_beta": 2.0, "aromatic": True}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _build_huckel_matrix(self, n: int, topology: str, connectivity: list = None) -> List[List[float]]:
        """Build n×n Hückel matrix H_ij.
        
        H_ii = α = 0 (in units of α)
        H_ij = β = 1 if atoms i,j are bonded, else 0
        We use units where α=0, β=1 → eigenvalues are in units of β (offset by α).
        """
        if connectivity:
            # Custom connectivity
            adj = [[0.0] * n for _ in range(n)]
            for i, j in connectivity:
                if 0 <= i < n and 0 <= j < n:
                    adj[i][j] = adj[j][i] = 1.0
            return adj

        H = [[0.0] * n for _ in range(n)]

        if topology == "cyclic":
            for i in range(n):
                H[i][i] = 0.0  # α = 0 in α+β units
                H[i][(i + 1) % n] = 1.0  # β = 1
                H[(i + 1) % n][i] = 1.0
        elif topology == "branched":
            # Default linear for now
            for i in range(n):
                H[i][i] = 0.0
                if i > 0:
                    H[i][i-1] = H[i-1][i] = 1.0
        else:  # linear
            for i in range(n):
                H[i][i] = 0.0
                if i > 0:
                    H[i][i-1] = H[i-1][i] = 1.0

        return H

    def _solve_eigenvalues(self, H: List[List[float]]) -> tuple:
        """
        Solve eigenvalue problem for symmetric real matrix using Jacobi/QR iteration.
        
        Returns (eigenvalues_sorted, eigenvectors) where eigenvectors[v][i] = coefficient of AO i in MO v.
        Uses power iteration + deflation for simplicity.
        """
        n = len(H)

        # Special case: analytic solutions for uniform chains
        # For linear chain: ε_k = α + 2β·cos(kπ/(n+1)), k=1..n
        # For cyclic: ε_k = α + 2β·cos(2kπ/n), k=0..n-1
        
        # Check if it's a standard form (tridiagonal with ones on off-diagonals)
        is_standard_linear = all(
            abs(H[i][j] - (1.0 if abs(i-j) == 1 else 0.0)) < 1e-10
            for i in range(n) for j in range(n)
        ) and all(abs(H[i][i]) < 1e-10 for i in range(n))

        is_cyclic = (
            all(abs(H[i][i]) < 1e-10 for i in range(n)) and
            all(abs(H[i][(i+1)%n] - 1.0) < 1e-10 for i in range(n)) and
            all(abs(H[(i+1)%n][i] - 1.0) < 1e-10 for i in range(n)) and
            sum(1 for i in range(n) for j in range(n) if abs(H[i][j]) > 0.5) == 2 * n
        )

        if is_standard_linear:
            eigenvalues = []
            eigenvectors = []
            for k in range(1, n + 1):
                eps = 2.0 * math.cos(k * math.pi / (n + 1))
                eigenvalues.append(eps)
                vec = [math.sin(k * i * math.pi / (n + 1)) for i in range(1, n + 1)]
                norm = math.sqrt(sum(v * v for v in vec))
                vec = [v / norm for v in vec]
                eigenvectors.append(vec)
            return eigenvalues, eigenvectors

        if is_cyclic:
            eigenvalues = []
            eigenvectors = []
            for k in range(n):
                eps = 2.0 * math.cos(2.0 * k * math.pi / n)
                eigenvalues.append(eps)
                vec = []
                phase = 2.0 * math.pi * k / n
                for i in range(n):
                    if k == 0:
                        val = 1.0 / math.sqrt(n)
                    elif 2 * k == n:
                        val = ((-1) ** i) / math.sqrt(n)
                    else:
                        val = math.cos(i * phase) / math.sqrt(n / 2.0)
                    vec.append(val)
                # Normalize
                norm = math.sqrt(sum(v * v for v in vec))
                vec = [v / norm for v in vec]
                eigenvectors.append(vec)
            return eigenvalues, eigenvectors

        # General case: numerical diagonalization via Jacobi-like iteration
        # Use power iteration for dominant eigenvalue + deflation
        eigenvalues = []
        eigenvectors = []
        H_work = [row[:] for row in H]

        for _mode in range(n):
            # Power iteration
            v = [1.0 / math.sqrt(n)] * n
            for _iter in range(3000):
                w = [sum(H_work[i][j] * v[j] for j in range(n)) for i in range(n)]
                nw = math.sqrt(sum(x * x for x in w))
                if nw < 1e-30:
                    break
                v_new = [x / nw for x in w]
                diff = max(abs(v_new[i] - v[i]) for i in range(n))
                v = v_new
                if diff < 1e-12:
                    break
            
            # Rayleigh quotient
            ev = sum(v[i] * sum(H_work[i][j] * v[j] for j in range(n)) for i in range(n))
            eigenvalues.append(ev)
            eigenvectors.append(list(v))

            # Deflation
            if _mode < n - 1:
                for i in range(n):
                    for j in range(n):
                        H_work[i][j] -= ev * v[i] * v[j]

        # Sort by eigenvalue
        paired = sorted(zip(eigenvalues, eigenvectors), key=lambda p: p[0])
        return [p[0] for p in paired], [p[1] for p in paired]

    def _get_molecule_params(self, mol: str) -> dict:
        """Get molecule parameters: n_carbons, topology, connectivity."""
        params = {
            "ethene": {"n": 2, "topo": "linear"},
            "ethylene": {"n": 2, "topo": "linear"},
            "butadiene": {"n": 4, "topo": "linear"},
            "hexatriene": {"n": 6, "topo": "linear"},
            "octatetraene": {"n": 8, "topo": "linear"},
            "benzene": {"n": 6, "topo": "cyclic"},
            "cyclobutadiene": {"n": 4, "topo": "cyclic"},
            "cyclooctatetraene": {"n": 8, "topo": "cyclic"},
            "allyl": {"n": 3, "topo": "linear"},  # allyl radical/cation/anion
            "allyl_radical": {"n": 3, "topo": "linear"},
            "pentadienyl": {"n": 5, "topo": "linear"},
            "cyclopentadienyl": {"n": 5, "topo": "cyclic"},
            "naphthalene": {"n": 10, "topo": None, "conn": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,0),(4,9)]},
            "fulvene": {"n": 6, "topo": None, "conn": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(5,2)]},  # 5-membered ring + exocyclic
        }
        key = mol.lower().strip()
        return params.get(key, params.get(mol))

    def _count_pi_electrons(self, n: int, ion_state: str, charge: int) -> int:
        """Count π electrons based on system type."""
        base_e = n  # Each carbon contributes 1 π electron in neutral state
        
        if ion_state == "cation":
            base_e -= 1
        elif ion_state == "anion":
            base_e += 1
        
        base_e -= charge  # Positive charge removes electrons
        
        return max(0, base_e)

    def _run_base(self, molecule: str, n_carbons: int = None, topology: str = None,
                  charge: int = 0, ionization_state: str = "neutral") -> dict:

        # Get molecule parameters
        mol_params = self._get_molecule_params(molecule)
        
        if mol_params:
            n = n_carbons if n_carbons else mol_params["n"]
            topo = topology if topology else mol_params["topo"]
            conn = mol_params.get("conn")
        else:
            if n_carbons is None:
                raise ChemMCPError(f"Unknown molecule '{molecule}' and n_carbons not specified.")
            n = n_carbons
            topo = topology if topology else "linear"
            conn = None

        # Build and solve Hückel matrix
        H = self._build_huckel_matrix(n, topo, conn)
        energies, coeffs = self._solve_eigenvalues(H)

        # Count π electrons
        n_pi_e = self._count_pi_electrons(n, ionization_state, charge)

        # Fill electrons (Aufbau principle): 2 per orbital (spin up/down)
        remaining = n_pi_e
        mo_data = []
        for idx in range(len(energies)):
            e_occ = min(2, max(0, remaining))
            mo_data.append({
                "orbital_index": idx,
                "energy_alpha_plus_x_beta": round(energies[idx], 8),
                "energy_description": f"α + {energies[idx]:.4f}β",
                "coefficients": [round(c, 6) for c in coeffs[idx]],
                "electrons": e_occ,
                "occupation": "filled" if e_occ == 2 else ("singly" if e_occ == 1 else "empty"),
            })
            remaining -= e_occ

        # Total π energy
        E_pi = sum(md["electrons"] * md["energy_alpha_plus_x_beta"] for md in mo_data)
        E_pi_str = f"{E_pi:.4f}α + {E_pi:.4f}β" if E_pi >= 0 else f"{E_pi:.4f}α {E_pi:.4f}β"

        # Delocalization energy (vs localized double bonds)
        n_dbonds = n_pi_e // 2  # Number of localized π bonds
        E_loc = n_dbonds * 2.0  # Each localized bond contributes 2α + 2β → relative to α: 2β each
        E_deloc = E_pi - E_loc

        # Resonance energy for aromatic systems
        resonance = None
        aromatic = False
        if topo == "cyclic" and n_pi_e == n and n % 2 == 0:  # Hückel's rule: 4n+2
            # Check 4n+2 rule
            huckel_num = (n_pi_e - 2) / 4
            if huckel_num == int(huckel_num) and n_pi_e >= 6:
                aromatic = True
            # For benzene: RE = 2β (3 double bonds × 2β = 6β vs actual 8β)
            if molecule.lower() in ("benzene",):
                resonance = 2.0 * abs(E_pi / n if n > 0 else 1)  # Approximate

        # Charge densities: q_i = Σ_k n_k · c_ki²
        q_densities = [0.0] * n
        for md in mo_data:
            nk = md["electrons"]
            for i in range(n):
                q_densities[i] += nk * md["coefficients"][i] ** 2

        # Bond orders: P_ij = Σ_k n_k · c_ki · c_kj
        bond_orders = {}
        for i in range(n):
            for j in range(i + 1, n):
                # Check if bonded
                is_bonded = False
                if conn:
                    is_bonded = (i, j) in conn or (j, i) in conn
                else:
                    is_bonded = (abs(i - j) == 1) or (topo == "cyclic" and (j == (i+1) % n or i == (j+1) % n))
                
                if is_bonded:
                    P_ij = 0.0
                    for md in mo_data:
                        nk = md["electrons"]
                        P_ij += nk * md["coefficients"][i] * md["coefficients"][j]
                    bond_orders[(i, j)] = round(P_ij, 6)

        # Frontier orbitals
        homo_idx = None
        lumo_idx = None
        for i, md in enumerate(mo_data):
            if md["electrons"] > 0:
                homo_idx = i
            if lumo_idx is None and md["electrons"] == 0:
                lumo_idx = i

        # Spectral transition estimate (HOMO→LUMO)
        gap_beta = None
        if homo_idx is not None and lumo_idx is not None:
            gap_beta = mo_data[lumo_idx]["energy_alpha_plus_x_beta"] - mo_data[homo_idx]["energy_alpha_plus_x_beta"]
            # |β| ≈ 2.5 eV for C=C bonds
            gap_eV_est = abs(gap_beta) * 2.5

        result = {
            "molecule": molecule,
            "n_carbons": n,
            "topology": topo,
            "charge": charge,
            "ionization_state": ionization_state,
            "n_pi_electrons": n_pi_e,
            "orbital_energies": mo_data,
            "total_pi_energy_str": E_pi_str,
            "total_pi_energy_in_beta_units": round(E_pi, 4),
            "delocalization_energy_beta": round(E_deloc, 4),
            "resonance_energy_beta": round(resonance, 4) if resonance is not None else None,
            "is_aromatic": aromatic,
            "huckel_4n2_rule_satisfied": aromatic,
            "charge_densities": [round(q, 4) for q in q_densities],
            "bond_orders": bond_orders,
            "frontier_orbital_info": {
                "homo_index": homo_idx,
                "lumo_index": lumo_idx,
                "homo_energy_beta": round(mo_data[homo_idx]["energy_alpha_plus_x_beta"], 4) if homo_idx is not None else None,
                "lumo_energy_beta": round(mo_data[lumo_idx]["energy_alpha_plus_x_beta"], 4) if lumo_idx is not None else None,
                "gap_beta": round(gap_beta, 4) if gap_beta is not None else None,
                "estimated_gap_eV": round(gap_eV_est, 2) if gap_beta is not None else None,
            },
        }

        logger.info(f"HuckelMethod: {molecule} (n={n}, {topo}), E_π={E_pi_str}, deloc={E_deloc:.4f}β")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mol = parts[0]
            kwargs = {}
            if len(parts) > 1:
                try:
                    kwargs["n_carbons"] = int(parts[1])
                except ValueError:
                    kwargs["topology"] = parts[1]
            if len(parts) > 2:
                try:
                    kwargs["charge"] = int(parts[2])
                except ValueError:
                    kwargs["topology"] = parts[2]
            if len(parts) > 3:
                kwargs["ionization_state"] = parts[3]
            return self._run_base(mol, **kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'molecule [n|topo] [charge] [ion_state]'")
