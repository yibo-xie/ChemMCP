import logging
import math
from typing import Optional, List, Dict
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ElectronDensityPlotter(BaseTool):
    """
    电子密度分布可视化。
    
    计算原子和分子的电子密度分布 |ψ(r)|²。
    生成径向分布函数 D(r) = 4πr²|R(r)|²。
    确定最概然半径、<r>、<r²> 以及等密度面信息。
    """
    __version__ = "0.1.0"
    name = "ElectronDensityPlotter"
    func_name = "electron_density_plotter"
    description = "Calculate and visualize electron density distribution for hydrogen-like atoms and quantum systems."
    implementation_description = "Computes electron probability density |ψ(r)|² for hydrogen-like atomic orbitals using analytical wavefunctions. Generates radial distribution function D(r)=4πr²R(r)², finds most probable radius, expectation values ⟨r⟩ and ⟨r²⟩, cumulative probability distribution, isosurface radii, and provides visualization data points for plotting density profiles."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Electron Density", "Visualization", "Atomic Orbitals", "Radial Distribution"]
    required_envs = []

    code_input_sig = [
        ("species_type", "str", "N/A", "'hydrogen', 'helium', 'hydrogen_molecule', 'particle_in_box', 'harmonic_oscillator', 'custom_hydrogen_like'."),
        ("quantum_numbers", "dict", "N/A", "Dictionary {n, l, m} for atomic orbitals (e.g., {'n':2,'l':1,'m':0})."),
        ("position_grid", "dict", "None", "Grid parameters: {r_min_bohr, r_max_bohr, n_points}. Auto-set if None."),
        ("Z", "int", "1", "Nuclear charge for hydrogen-like ions."),
        ("isosurface_level", "float", "None", "Fraction of max density for isosurface (default=0.95)."),
        ("n_plot_points", "int", "200", "Number of radial grid points for output data."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'species n l m [Z] [n_points]'. Example: 'hydrogen 2 1 0 1 200'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with density data arrays, most probable radius, mean radius, radial distribution data, isosurface info, visualization data."),
    ]

    examples = [
        {
            "code_input": {
                "species_type": "hydrogen",
                "quantum_numbers": {"n": 1, "l": 0, "m": 0},
                "position_grid": None,
                "Z": 1,
                "isosurface_level": None,
                "n_plot_points": 200,
            },
            "text_input": {"input_params": "hydrogen 1 0 0"},
            "output": {"result": {"most_probable_radius_angstrom": 0.529, "orbital": "1s"}},
        },
        {
            "code_input": {
                "species_type": "hydrogen",
                "quantum_numbers": {"n": 2, "l": 1, "m": 0},
                "position_grid": None,
                "Z": 1,
                "isosurface_level": None,
                "n_plot_points": 200,
            },
            "text_input": {"input_params": "hydrogen 2 1 0"},
            "output": {"result": {"most_probable_radius_angstrom": 2.116, "orbital": "2p"}},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.a0 = 5.29177210903e-11  # Bohr radius in meters
        self.a0_angstrom = 0.529177210903  # Bohr in Å

    def _factorial(self, n: int) -> int:
        if n <= 1:
            return 1
        r = 1
        for i in range(2, n + 1):
            r *= i
        return r

    def _double_factorial(self, n: int) -> int:
        if n <= 0:
            return 1
        r = 1
        while n > 0:
            r *= n
            n -= 2
        return r

    def _associated_laguerre(self, n: int, k: int, x: float) -> float:
        """Associated Laguerre polynomial L_n^k(x)."""
        if n == 0:
            return 1.0
        if n == 1:
            return float(k + 1 - x)
        L_prev2 = 1.0
        L_prev1 = float(k + 1 - x)
        for i in range(2, n + 1):
            L_curr = ((2 * i + k - 1 - x) * L_prev1 - (i + k - 1) * L_prev2) / i
            L_prev2 = L_prev1
            L_prev1 = L_curr
        return L_prev1

    def _radial_wavefunction(self, n: int, l: int, rho: float, Z: int) -> float:
        """
        R_nl(r) where ρ = 2Zr/(na₀).
        
        R_nl = sqrt((2Z/na₀)³ · (n-l-1)!/(2n[(n+l)!]³)) · ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2)
        """
        nf = self._factorial(n - l - 1)
        nlf = self._factorial(n + l)
        prefactor = math.sqrt((2.0 * Z / (n ** 3)) ** 3 * nf / (2.0 * n * (nlf ** 3)))
        return prefactor * (rho ** l) * self._associated_laguerre(n - l - 1, 2 * l + 1, rho) * math.exp(-rho / 2.0)

    def _angular_part(self, l: int, m: int, theta: float) -> float:
        """Simplified angular part magnitude |Y_lm(θ, φ)| (averaged over φ)."""
        # Use Legendre polynomial P_l^|m|(cos θ) approximation
        # For simplicity, return characteristic values
        if l == 0:
            return 1.0 / math.sqrt(4.0 * math.pi)
        elif l == 1:
            if m == 0:
                return math.sqrt(3.0 / (4.0 * math.pi)) * abs(math.cos(theta))
            else:
                return math.sqrt(3.0 / (8.0 * math.pi)) * abs(math.sin(theta))
        elif l == 2:
            if m == 0:
                return math.sqrt(5.0 / (16.0 * math.pi)) * abs(3.0 * math.cos(theta)**2 - 1.0)
            elif abs(m) == 1:
                return math.sqrt(15.0 / (8.0 * math.pi)) * abs(math.sin(theta) * math.cos(theta))
            else:
                return math.sqrt(15.0 / (32.0 * math.pi)) * math.sin(theta)**2
        return 1.0

    def _generate_radial_data(self, n: int, l: int, Z: int, n_points: int,
                                r_max_bohr: float) -> dict:
        """Generate radial density data on a grid."""
        r_min = 0.01  # Avoid r=0 singularity
        dr = (r_max_bohr - r_min) / (n_points - 1)

        r_vals = []
        density_vals = []
        rdf_vals = []  # Radial distribution function
        R_vals = []

        max_density = 0.0
        max_rdf = 0.0
        mp_radius = r_min
        mp_rdf_radius = r_min

        for i in range(n_points):
            r = r_min + i * dr
            rho = 2.0 * Z * r / n

            R = self._radial_wavefunction(n, l, rho, Z)
            prob = R * R  # |R(r)|²
            rdf = r * r * prob * 4.0 * math.pi  # D(r) = 4πr²|R|² (in units where a₀=1)

            r_vals.append(round(r, 6))
            R_vals.append(round(R, 12))
            density_vals.append(round(prob, 12))
            rdf_vals.append(round(rdf, 12))

            if prob > max_density:
                max_density = prob
                mp_radius = r
            if rdf > max_rdf:
                max_rdf = rdf
                mp_rdf_radius = r

        # Cumulative probability
        cum_prob = [0.0]
        for i in range(1, len(density_vals)):
            # ∫₀ʳ 4πr'|R(r')|² dr' ≈ trapezoidal rule
            integral = cum_prob[-1] + 0.5 * dr * (rdf_vals[i] + rdf_vals[i-1])
            cum_prob.append(round(integral, 10))

        # Normalize cumulative to 1 at large r
        if cum_prob[-1] > 0:
            cum_prob = [c / cum_prob[-1] for c in cum_prob]

        return {
            "r_bohr": r_vals,
            "density_R_squared": density_vals,
            "radial_distribution_D_r": rdf_vals,
            "radial_wavefunction_R": R_vals,
            "cumulative_probability": cum_prob,
            "most_probable_radius_bohr": round(mp_radius, 4),
            "most_probable_rdf_radius_bohr": round(mp_rdf_radius, 4),
            "max_density": round(max_density, 10),
        }

    def _compute_expectation_radii(self, n: int, l: int, Z: int, n_points: int,
                                     r_max_bohr: float) -> dict:
        """Compute <r>, <r²>, Δr by numerical integration."""
        r_min = 0.001
        dr = r_max_bohr / n_points
        
        sum_r = 0.0
        sum_r2 = 0.0
        norm = 0.0

        for i in range(n_points):
            r = r_min + i * dr
            rho = 2.0 * Z * r / n
            R = self._radial_wavefunction(n, l, rho, Z)
            prob = R * R * 4.0 * math.pi * r * r  # Volume element included
            
            norm += prob * dr
            sum_r += r * prob * dr
            sum_r2 += r * r * prob * dr

        if norm > 1e-30:
            mean_r = sum_r / norm
            mean_r2 = sum_r2 / norm
        else:
            mean_r = r_max_bohr / 2.0
            mean_r2 = mean_r ** 2

        delta_r = math.sqrt(max(0, mean_r2 - mean_r ** 2))

        return {
            "mean_radius_bohr": round(mean_r, 6),
            "mean_radius_sq_bohr2": round(mean_r2, 6),
            "uncertainty_delta_r_bohr": round(delta_r, 6),
            "mean_radius_angstrom": round(mean_r * self.a0_angstrom, 6),
        }

    def _find_isosurface(self, n: int, l: int, Z: int, level: float, n_points: int,
                          r_max_bohr: float) -> dict:
        """Find isosurface radius where density drops to given fraction of max."""
        r_min = 0.001
        dr = r_max_bohr / n_points

        max_dens = 0.0
        # First pass: find max
        for i in range(n_points):
            r = r_min + i * dr
            rho = 2.0 * Z * r / n
            R = self._radial_wavefunction(n, l, rho, Z)
            d = R * R
            if d > max_dens:
                max_dens = d

        target = level * max_dens
        iso_radius = None

        # Second pass: find where density crosses target (from outside in)
        for i in range(n_points - 1, -1, -1):
            r = r_min + i * dr
            rho = 2.0 * Z * r / n
            R = self._radial_wavefunction(n, l, rho, Z)
            if R * R >= target:
                iso_radius = r
                break

        return {
            "isosurface_fraction": level,
            "isosurface_radius_bohr": round(iso_radius, 4) if iso_radius else None,
            "isosurface_radius_angstrom": round(iso_radius * self.a0_angstrom, 4) if iso_radius else None,
            "max_density_value": round(max_dens, 10),
        }

    def _run_base(self, species_type: str, quantum_numbers: dict,
                  position_grid: dict = None, Z: int = 1,
                  isosurface_level: float = None, n_plot_points: int = 200) -> dict:

        n = quantum_numbers.get("n", 1)
        l = quantum_numbers.get("l", 0)
        m = quantum_numbers.get("m", 0)

        # Validation
        if n < 1:
            raise ChemMCPError("Principal quantum number n must be >= 1.")
        if l < 0 or l >= n:
            raise ChemMCPError(f"Must have 0 <= l < n (got l={l}, n={n}).")
        if Z < 1:
            raise ChemMCPError("Nuclear charge Z must be >= 1.")

        # Orbital label
        orbital_labels = {0: "s", 1: "p", 2: "d", 3: "f"}
        orbital_name = f"{n}{orbital_labels.get(l, '?')}"

        # Set up radial grid
        if position_grid:
            r_max = position_grid.get("r_max_bohr", n * n * 6.0 / Z)
            n_pts = position_grid.get("n_points", n_plot_points)
        else:
            r_max = n * n * 6.0 / Z  # Extends to ~6× most probable radius
            n_pts = n_plot_points

        # Generate all data
        radial_data = self._generate_radial_data(n, l, Z, n_pts, r_max)
        exp_values = self._compute_expectation_radii(n, l, Z, n_pts, r_max)

        iso_level = isosurface_level if isosurface_level else 0.95
        iso_info = self._find_isosurface(n, l, Z, iso_level, n_pts, r_max)

        # Node structure
        n_radial_nodes = n - l - 1
        n_angular_nodes = l

        result = {
            "species_type": species_type,
            "quantum_numbers": {"n": n, "l": l, "m": m},
            "orbital_name": orbital_name,
            "nuclear_charge_Z": Z,
            "node_structure": {
                "n_radial_nodes": n_radial_nodes,
                "n_angular_nodes": n_angular_nodes,
                "total_nodes": n - 1,
            },
            "radial_data": radial_data,
            "expectation_radii": exp_values,
            "isosurface_info": iso_info,
            "visualization_summary": {
                "n_plot_points": n_pts,
                "r_range_bohr": [round(radial_data["r_bohr"][0], 4), round(radial_data["r_bohr"][-1], 4)],
                "peak_location_bohr": radial_data["most_probable_radius_bohr"],
                "rdf_peak_location_bohr": radial_data["most_probable_rdf_radius_bohr"],
            },
        }

        logger.info(f"ElectronDensityPlotter: {orbital_name}, r_mp={radial_data['most_probable_radius_bohr']:.3f}a₀")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            species = parts[0]
            n = int(parts[1])
            l = int(parts[2])
            m = int(parts[3]) if len(parts) > 3 else 0
            Z = int(parts[4]) if len(parts) > 4 else 1
            npts = int(parts[5]) if len(parts) > 5 else 200
            return self._run_base(species, {"n": n, "l": l, "m": m}, Z=Z, n_plot_points=npts)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'species n l m [Z] [n_points]'")
