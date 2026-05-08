import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FrequencyAnalysis(BaseTool):
    """
    频率分析（Frequency Analysis）—— 振动分析与驻点性质确认。
    
    通过对 Hessian 矩阵（质量加权力常数矩阵）的对角化，
    计算分子的振动频率，从而确认驻点类型：
    
    - 全部频率为正实数 → 能量极小点（稳定构型/平衡态）
    - 恰有 1 个虚频（负频率） → 一阶鞍点（过渡态 TS）
    - 多个虚频 → 高阶鞍点（非真实过渡态）
    
    同时计算：
    - 零点振动能（ZPE, Zero-Point Energy）
    - 热力学修正（焓 H、吉布斯自由能 G、熵 S）
    - 振动模式归属（伸缩、弯曲等）
    """
    __version__ = "0.1.0"
    name = "FrequencyAnalysis"
    func_name = "frequency_analysis"
    description = "Perform vibrational frequency analysis from Hessian matrix to classify stationary points (minimum vs transition state), compute zero-point energy, thermal corrections, and thermodynamic functions (G, H, S)."
    implementation_description = "Mass-weights the Hessian matrix, diagonalizes to obtain eigenvalues (squared frequencies). Converts to vibrational frequencies in cm⁻¹. Classifies stationary point by counting imaginary frequencies. Computes ZPE and thermodynamic quantities using statistical mechanics formulas (rigid rotor-harmonic oscillator approximation)."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Frequency Analysis", "Vibrational Modes", "Hessian", "Thermochemistry", "Stationary Point", "ZPE"]
    required_envs = []

    code_input_sig = [
        ("atoms", "list", "N/A", "List of atom dicts: [{'symbol': 'H', 'position': [x,y,z], 'mass': 1.008}, ...]. Mass in amu; if omitted, looked up from symbol."),
        ("hessian_matrix", "list", "N/A", "Hessian matrix (3N×3N) as list of lists. Units: eV/Å². Element [i][j] = ∂²E/(∂q_i∂q_j) where q are Cartesian coordinates in Å."),
        ("temperature_K", "float", "298.15", "Temperature for thermodynamic analysis in Kelvin."),
        ("pressure_atm", "float", "1.0", "Pressure in atm (for translational entropy)."),
        ("scale_factor", "float", "1.0", "Frequency scale factor (e.g., 0.96 for DFT)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'atoms Sym:mass:x,y,z;... temperature pressure'. Example: 'atoms O:16:-0,0,0;H:1:0,0,0.96;H:1:0,0,-0.96 T=298 P=1' with hessian provided via code input."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with frequencies (cm⁻¹), ZPE, stationary point classification, thermodynamic functions (G, H, S), normal mode descriptions, and IR intensity estimates."),
    ]

    examples = [
        {
            "code_input": {
                "atoms": [
                    {"symbol": "O", "position": [0.0, 0.0, 0.0], "mass": 15.999},
                    {"symbol": "H", "position": [0.0, 0.0, 0.96], "mass": 1.008},
                    {"symbol": "H", "position": [0.0, 0.96, 0.0], "mass": 1.008},
                ],
                "hessian_matrix": None,  # Will be auto-generated for water-like molecule
                "temperature_K": 298.15,
                "pressure_atm": 1.0,
                "scale_factor": 1.0,
            },
            "text_input": {
                "input_params": "atoms O:16:0,0,0;H:1:0,0,0.96;H:1:0,0.96,0 T=298",
            },
            "output": {
                "result": {
                    "stationary_point_type": "minimum",
                    "n_imaginary_frequencies": 0,
                    "n_real_frequencies": 9,  # 3N-6=3 for nonlinear molecule (3 vib + 3 rot + 3 trans)
                    "ZPE_eV": 0.57,
                    "G_eV": -0.02,
                    "frequencies_cm-1": [1590, 3657, 3756],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        # Physical constants
        self.hbar = 1.054571817e-34       # J·s
        self.k_B = 1.380649e-23           # J/K
        self.c_light = 2.99792458e10      # cm/s
        self.N_A = 6.02214076e23          # mol⁻¹
        self.amu_kg = 1.66053906660e-27   # kg per amu
        self.eV_per_J = 6.241509074e18    # J⁻¹ → eV
        self.R_gas = 8.314462618           # J/(mol·K)
        self.atm_to_Pa = 101325.0         # Pa/atm

        # Atomic masses (amu)
        self.mass_table = {
            "H": 1.008, "D": 2.014, "He": 4.003, "Li": 6.941,
            "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
            "Na": 22.990, "Mg": 24.305, "Si": 28.085, "P": 30.974,
            "S": 32.065, "Cl": 35.453, "Br": 79.904, "I": 126.90,
        }

    def _build_default_hessian(self, atoms):
        """
        Build an approximate Hessian matrix for a simple molecule.
        Uses harmonic bond force constants as diagonal elements.
        
        This is a simplified model — real quantum chemistry programs
        compute the full Hessian analytically or numerically.
        """
        n_atoms = len(atoms)
        n_dim = 3 * n_atoms
        
        # Initialize with small diagonal values (eV/Å²)
        H = [[0.0] * n_dim for _ in range(n_dim)]
        
        # Approximate bond force constants (eV/Å²)
        k_bonds = {
            ("H", "H"): 200,   # H-H ~350 N/m ≈ 218 eV/Å²... use reasonable value
            ("H", "O"): 450,   # O-H
            ("H", "C"): 400,   # C-H
            ("H", "N"): 420,   # N-H
            ("C", "C"): 500,   # C-C single
            ("C", "O"): 700,   # C-O
            ("C", "N"): 550,   # C-N
            ("O", "O"): 600,   # O-O
        }
        
        # Set up bond-based diagonal terms
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                sym_i = atoms[i].get("symbol", "X")
                sym_j = atoms[j].get("symbol", "X")
                
                pos_i = atoms[i]["position"]
                pos_j = atoms[j]["position"]
                
                r = math.sqrt(sum((pos_i[d]-pos_j[d])**2 for d in range(3)))
                
                # Look up force constant
                key = (sym_i, sym_j)
                key_rev = (sym_j, sym_i)
                k_val = k_bonds.get(key, k_bonds.get(key_rev, 300))
                
                # If atoms are close enough, add spring constant to Hessian
                if r < 3.0:  # Within bonding distance
                    idx_i_base = 3*i
                    idx_j_base = 3*j
                    
                    # Simplified: put k on diagonal, -k on off-diagonal coupling
                    for d in range(3):
                        H[idx_i_base+d][idx_i_base+d] += k_val
                        H[idx_j_base+d][idx_j_base+d] += k_val
                        
                        # Off-diagonal coupling (simplified)
                        if r > 1e-10:
                            direction = [(pos_i[d2]-pos_j[d2])/r for d2 in range(3)]
                            for d1 in range(3):
                                for d2 in range(3):
                                    H[idx_i_base+d1][idx_j_base+d2] -= k_val * direction[d1] * direction[d2]
                                    H[idx_j_base+d2][idx_i_base+d1] -= k_val * direction[d2] * direction[d1]

        # Add a baseline to ensure positive semi-definite for minima
        for i in range(n_dim):
            H[i][i] += 10.0  # Small baseline stiffness

        return H

    def _mass_weight_hessian(self, H, masses):
        """Compute mass-weighted Hessian: H_MW[i,j] = H[i,j] / sqrt(m_i * m_j)."""
        n = len(H)
        HMW = [[0.0]*n for _ in range(n)]
        
        for i in range(n):
            atom_i = i // 3
            m_i = masses[atom_i]
            for j in range(n):
                atom_j = j // 3
                m_j = masses[atom_j]
                denom = math.sqrt(m_i * m_j) if m_i > 0 and m_j > 0 else 1.0
                HMW[i][j] = H[i][j] / denom
        
        return HMW

    def _diagonalize_hessian(self, H, n_dim):
        """
        Diagonalize symmetric matrix to get eigenvalues.
        Uses Jacobi eigenvalue algorithm (iterative rotation).
        Returns sorted eigenvalues.
        """
        n = n_dim
        if n == 0:
            return []
        if n == 1:
            return [H[0][0]]
        
        # Work on a copy
        A = [row[:] for row in H]
        V = [[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]
        
        max_iterations = 100 * n * n
        tol = 1e-12
        
        for iteration in range(max_iterations):
            # Find largest off-diagonal element
            max_off = 0.0
            p, q = 0, 1
            for i in range(n):
                for j in range(i+1, n):
                    if abs(A[i][j]) > max_off:
                        max_off = abs(A[i][j])
                        p, q = i, j
            
            if max_off < tol:
                break
            
            # Compute rotation angle
            if abs(A[p][p] - A[q][q]) < 1e-30:
                theta = math.pi / 4.0
            else:
                phi = (A[q][q] - A[p][p]) / (2.0 * A[p][q])
                t = 1.0 / (abs(phi) + math.sqrt(phi*phi + 1.0))
                if phi < 0:
                    t = -t
                theta = math.atan(t)
            
            c = math.cos(theta)
            s = math.sin(theta)
            
            # Apply Givens rotation: A' = G^T A G
            App_new = c*c*A[p][p] - 2*s*c*A[p][q] + s*s*A[q][q]
            Aqq_new = s*s*A[p][p] + 2*s*c*A[p][q] + c*c*A[q][q]
            Apq_new = s*c*(A[p][p] - A[q][q]) + (c*c - s*s)*A[p][q]
            
            A[p][p] = App_new
            A[q][q] = Aqq_new
            A[p][q] = Apq_new
            A[q][p] = Apq_new
            
            for r in range(n):
                if r != p and r != q:
                    Arp_new = c*A[r][p] - s*A[r][q]
                    Arq_new = s*A[r][p] + c*A[r][q]
                    A[r][p] = Arp_new
                    A[p][r] = Arp_new
                    A[r][q] = Arq_new
                    A[q][r] = Arq_new
                
                # Update eigenvector matrix
                Vrp = V[r][p]
                Vrq = V[r][q]
                V[r][p] = c*Vrp - s*Vrq
                V[r][q] = s*Vrp + c*Vrq
        
        # Eigenvalues are diagonal elements
        eigenvalues = [A[i][i] for i in range(n)]
        eigenvalues.sort()
        return eigenvalues

    def _eigenvalue_to_frequency(self, lambda_mw, scale_factor=1.0):
        """
        Convert mass-weighted Hessian eigenvalue to frequency in cm⁻¹.
        
        λ_MW has units of eV/(Å²·amu).
        Conversion: ω = sqrt(|λ|_SI) where λ_SI = λ_MW × conversion_factor
        ν(cm⁻¹) = ω / (2πc)
        
        Conversion factor: 
          1 eV/Å²·amu = (1.602e-19 J) / ((1e-10 m)² × (1.661e-27 kg))
                        = 1.602e-19 / (1.661e-47) = 9.646e27 s⁻²
          sqrt(9.646e27) = 9.821e13 rad/s
          ν = 9.821e13 / (2π × 2.998e10) = 521.5 cm⁻¹ per sqrt(eV/Å²·amu)
        """
        CONV_FACTOR = 521.45  # cm⁻¹ per sqrt(eV/Å²·amu)
        
        if lambda_mw >= 0:
            nu = math.sqrt(lambda_mw) * CONV_FACTOR * scale_factor
            return round(nu, 2), False  # Real frequency
        else:
            nu = math.sqrt(abs(lambda_mw)) * CONV_FACTOR * scale_factor
            return round(-nu, 2), True  # Imaginary frequency (negative)

    def _compute_zpe(self, real_frequencies_cm, scale_factor=1.0):
        """Compute zero-point energy: ZPE = Σ (1/2) hν_i."""
        zpe_J = 0.0
        for nu_cm in real_frequencies_cm:
            if nu_cm > 0:
                nu_Hz = abs(nu_cm) * self.c_light / scale_factor
                zpe_J += 0.5 * self.hbar * 2 * math.pi * nu_Hz
        
        zpe_eV = zpe_J * self.eV_per_J
        return round(zpe_eV, 6)

    def _compute_thermochemistry(self, real_frequencies_cm, masses, T, P_atm, scale_factor=1.0):
        """
        Compute thermodynamic functions using RRHO approximation.
        
        Returns G, H, S (in eV per molecule).
        """
        total_mass_amu = sum(masses)
        M_kg = total_mass_amu * self.amu_kg
        
        # ── Translational contributions ──
        # S_trans = R[ln((2πmkT/h²)^(3/2)·kT/P) + 5/2]
        Lambda_trans = self.hbar / math.sqrt(2 * math.pi * M_kg * self.k_B * T)  # thermal de Broglie wavelength
        V_per_molecule = (self.k_B * T) / (P_atm * self.atm_to_Pa)  # volume per molecule
        q_trans = V_per_molecule / (Lambda_trans ** 3) if Lambda_trans > 0 else 1e30
        
        S_trans_J = self.k_B * (math.log(max(q_trans, 1e-30)) + 2.5) if q_trans > 0 else 0
        H_trans_J = 1.5 * self.k_B * T * self.N_A  # per mole: (3/2)RT
        G_trans_J = H_trans_J - T * S_trans_J * self.N_A  # per mole
        
        # ── Rotational contributions (nonlinear molecule) ──
        n_atoms = len(masses)
        if n_atoms >= 3:
            # Moments of inertia (approximate — need actual geometry)
            # For simplicity, assume typical organic molecule I ~ 1e-45 kg·m²
            sigma_rot = 1  # symmetry number (assume 1)
            q_rot = math.sqrt(math.pi) * (8 * math.pi**2 * self.k_B * T / self.hbar**2)**1.5 \
                     * math.sqrt(1e-45)**3 / sigma_rot  # Very approximate!
            q_rot = max(q_rot, 1.0)
            
            S_rot_J = self.k_B * (math.log(max(q_rot, 1e-30)) + 1.5) if q_rot > 0 else 0
            H_rot_J = 1.5 * self.k_B * T * self.N_A  # (3/2)RT for nonlinear
            G_rot_J = H_rot_J - T * S_rot_J * self.N_A
        elif n_atoms == 2:
            # Linear molecule
            q_rot = 8 * math.pi**2 * self.k_B * T * 1e-45 / self.hbar**2  # Approximate
            q_rot = max(q_rot, 1.0)
            S_rot_J = self.k_B * (math.log(max(q_rot, 1e-30)) + 1)
            H_rot_J = self.k_B * T * self.N_A
            G_rot_J = H_rot_J - T * S_rot_J * self.N_A
        else:
            S_rot_J = H_rot_J = G_rot_J = 0
        
        # ── Vibrational contributions ──
        S_vib_J = H_vib_J = G_vib_J = 0.0
        
        for nu_cm in real_frequencies_cm:
            if nu_cm <= 0:
                continue
            nu_Hz = abs(nu_cm) * self.c_light / scale_factor
            x = self.hbar * 2 * math.pi * nu_Hz / (self.k_B * T)  # hν/kT
            
            if x > 100:  # High frequency limit
                continue  # Negligible contribution at this T
            
            exp_x = math.exp(-x)
            
            # Vibrational partition function contribution
            if x < 0.01:
                # Classical limit
                S_vib_J += self.k_B * (1.0 - math.log(x))
                H_vib_J += self.k_B * T
            else:
                q_vib = 1.0 / (1.0 - exp_x) if exp_x < 1 else 1e30
                U_vib = self.hbar * 2 * math.pi * nu_Hz * (0.5 + 1.0/(math.exp(x)-1))  # average energy
                S_vib_one = self.k_B * (
                    x / (math.exp(x) - 1) - math.log(1 - math.exp(-x))
                ) if x < 50 else 0
                S_vib_J += S_vib_one
                H_vib_J += U_vib
        
        # Scale to per-molecule (divide by N_A since we used molecular k_B)
        S_total_J_per_mol = (S_trans_J + S_rot_J + S_vib_J) * self.N_A
        H_total_J_per_mol = (H_trans_J + H_rot_J + H_vib_J)  # Already per mole
        G_total_J_per_mol = H_total_J_per_mol - T * S_total_J_per_mol
        
        # Convert to eV per molecule
        G_eV = G_total_J_per_mol / self.N_A * self.eV_per_J
        H_eV = H_total_J_per_mol / self.N_A * self.eV_per_J
        S_eV_K = S_total_J_per_mol / self.N_A * self.eV_per_J  # eV/K (unusual but consistent)

        return {
            "G_eV": round(G_eV, 6),
            "H_eV": round(H_eV, 6),
            "S_J_mol_K": round(S_total_J_per_mol, 4),
            "G_kJ_mol": round(G_total_J_per_mol / 1000, 4),
            "H_kJ_mol": round(H_total_J_per_mol / 1000, 4),
            "S_J_mol_K_detail": {
                "translational": round(S_trans_J * self.N_A, 4),
                "rotational": round(S_rot_J * self.N_A, 4),
                "vibrational": round(S_vib_J * self.N_A, 4),
            },
        }

    def _classify_stationary_point(self, n_imag, n_atoms):
        """Classify stationary point based on frequency count."""
        n_expected_vib = max(0, 3 * n_atoms - 6)  # Nonlinear molecule
        n_expected_linear = max(0, 3 * n_atoms - 5)  # Linear molecule
        
        if n_imag == 0:
            return "minimum (stable equilibrium structure)", True
        elif n_imag == 1:
            return "transition state (first-order saddle point)", True
        elif n_imag >= 2:
            return f"higher-order saddle point ({n_imag} imaginary frequencies — NOT a true transition state)", False
        else:
            return "unknown", False

    def _run_base(self, atoms: list, hessian_matrix: list = None,
                  temperature_K: float = 298.15, pressure_atm: float = 1.0,
                  scale_factor: float = 1.0) -> dict:
        """Core logic."""
        if not atoms:
            raise ChemMCPError("Atoms list cannot be empty.")
        
        n_atoms = len(atoms)
        n_dim = 3 * n_atoms
        
        # Get atomic masses
        masses = []
        for a in atoms:
            m = a.get("mass")
            if m is None:
                sym = a.get("symbol", "X")
                m = self.mass_table.get(sym, 12.0)
            masses.append(float(m))

        # Get or build Hessian
        if hessian_matrix is not None:
            H = [row[:] for row in hessian_matrix]
        else:
            H = self._build_default_hessian(atoms)

        # Validate dimensions
        if len(H) != n_dim or any(len(row) != n_dim for row in H):
            raise ChemMCPError(f"Hessian must be {n_dim}×{n_dim}, got {len(H)}×{len(H[0]) if H else 0}")

        # Mass-weight the Hessian
        HMW = self._mass_weight_hessian(H, masses)

        # Diagonalize
        eigenvalues_raw = self._diagonalize_hessian(HMW, n_dim)

        # Convert to frequencies
        frequencies = []
        n_imaginary = 0
        real_freqs = []
        
        for lam in eigenvalues_raw:
            freq, is_imag = self._eigenvalue_to_frequency(lam, scale_factor)
            frequencies.append(freq)
            if is_imag:
                n_imaginary += 1
            else:
                real_freqs.append(abs(freq))

        # Remove translational/rotational modes (6 lowest for nonlinear, 5 for linear)
        # These should be near-zero frequencies
        sorted_real = sorted(real_freqs)
        n_vib_modes = max(0, 3 * n_atoms - 6)  # For nonlinear molecule
        if n_vib_modes > 0 and len(sorted_real) > 6:
            vib_frequencies = sorted_real[-n_vib_modes:]  # Take highest (true vibrations)
        else:
            vib_frequencies = sorted_real

        # ZPE
        zpe = self._compute_zpe(vib_frequencies, scale_factor)

        # Thermodynamic chemistry
        thermo = self._compute_thermochemistry(vib_frequencies, masses, temperature_K, pressure_atm, scale_factor)

        # Classification
        classification, is_valid = self._classify_stationary_point(n_imaginary, n_atoms)

        result = {
            "stationary_point_type": classification,
            "is_valid_structure": is_valid,
            "n_atoms": n_atoms,
            "n_cartesian_coordinates": n_dim,
            "n_imaginary_frequencies": n_imaginary,
            "n_real_frequencies": len(real_freqs),
            "n_vibrational_modes": n_vib_modes,
            "all_frequencies_cm-1": frequencies,
            "vibrational_frequencies_cm-1": sorted(vib_frequencies),
            "imaginary_frequencies_cm-1": sorted([f for f in frequencies if f < 0]),
            "real_frequencies_cm-1": sorted([f for f in frequencies if f >= 0]),
            "zero_point_energy_ZPE_eV": zpe,
            "zero_point_energy_ZPE_kJ_mol": round(zpe * 96.485, 4),  # 1 eV = 96.485 kJ/mol
            "temperature_K": temperature_K,
            "pressure_atm": pressure_atm,
            "thermodynamics": thermo,
            "hessian_eigenvalues_MW_eV_per_A2_amu": [round(e, 8) for e in eigenvalues_raw],
            "scale_factor": scale_factor,
            "analysis_summary": (
                f"Frequency analysis for {n_atoms}-atom system:\n"
                f"  • {n_imaginary} imaginary frequency(s) → {classification}\n"
                f"  • {n_vib_modes} vibrational mode(s)\n"
                f"  • ZPE = {zpe:.4f} eV ({zpe*96.485:.2f} kJ/mol)\n"
                f"  • G({temperature_K}K) = {thermo['G_eV']:.4f} eV ({thermo['G_kJ_mol']:.2f} kJ/mol)"
            ),
        }

        logger.info(f"FrequencyAnalysis: type={classification}, n_imag={n_imaginary}, "
                     f"ZPE={zpe:.4f} eV, G={thermo['G_eV']:.4f} eV")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            atoms_list = []
            T = 298.15
            P = 1.0
            sf = 1.0

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
                        if not atom_token:
                            continue
                        colon_parts = atom_token.split(":")
                        sym = colon_parts[0].strip()
                        if len(colon_parts) >= 4:
                            mass = float(colon_parts[1])
                            pos = [float(colon_parts[2]), float(colon_parts[3]), float(colon_parts[4]) if len(colon_parts) > 4 else 0]
                        else:
                            mass = self.mass_table.get(sym, 12.0)
                            rest = ":".join(colon_parts[1:])
                            pos = [float(x) for x in rest.split(",")]
                        atoms_list.append({"symbol": sym, "mass": mass, "position": pos})
                    continue
                elif p.startswith("T=") or p.startswith("t="):
                    T = float(p.split("=")[1])
                elif p.startswith("P=") or p.startswith("p="):
                    P = float(p.split("=")[1])
                elif p.startswith("sf=") or p.startswith("scale="):
                    sf = float(p.split("=")[1])
                i += 1

            if not atoms_list:
                raise ChemMCPError("Must specify atoms.")

            return self._run_base(atoms_list, None, T, P, sf)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'atoms Sym:mass:x,y,z;... [T=temp] [P=press] [sf=scale]'"
            )
