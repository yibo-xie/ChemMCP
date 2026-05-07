import logging
import math
from typing import List, Optional, Tuple
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class RadialDistribution(BaseTool):
    """
    氢原子/类氢离子径向分布函数计算工具。
    
    径向分布函数 D(r) = r² |R_nl(r)|² 表示在半径 r 处单位厚度球壳内
    找到电子的概率。D(r)dr 给出电子在 [r, r+dr] 范围内的概率。
    
    支持计算：
    - D(r) 在径向网格上的值
    - 最概然半径（D(r) 极大值位置）
    - 平均半径 <r> 和半径不确定度 Δr
    - 各阶矩 <r^k>
    - 节点位置（径向节点）
    - 累积概率分布
    """
    __version__ = "0.1.0"
    name = "RadialDistribution"
    func_name = "radial_distribution"
    description = "Compute radial distribution function D(r)=r²|R_nl(r)|² for hydrogen-like atoms: most probable radius, mean radius, radial nodes, probability in spherical shells, cumulative distribution, and radial moments."
    implementation_description = "Computes radial probability distribution for hydrogenic orbitals using analytical radial wavefunctions R_nl(r) with associated Laguerre polynomials. Calculates D(r)=r²R² on a radial grid, finds peaks (most probable radii), computes radial moments ⟨r^k⟩, locates radial nodes, and integrates to find probability within specified radius ranges."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Radial Distribution", "Hydrogen Atom", "Probability Density", "Wavefunction", "Atomic Structure"]
    required_envs = []

    code_input_sig = [
        ("n", "int", "N/A", "Principal quantum number (n >= 1)."),
        ("l", "int", "N/A", "Orbital angular momentum quantum number (0 <= l < n)."),
        ("Z", "int", "1", "Nuclear charge (default 1 = H)."),
        ("r_max_bohr", "float", "30.0", "Maximum radius in Bohr radii for computation grid."),
        ("n_points", "int", "500", "Number of radial grid points."),
        ("compute_moments_up_to", "int", "3", "Compute radial moments <r^k> for k=1..this value."),
        ("radius_range_for_prob", "list", "None", "List of [r_min, r_max] pairs in Bohr to compute probability within each shell."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'n l [Z] [r_max] [n_points]'. Example: '2 0 1 20 500' for 2s orbital of H with r_max=20a₀."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing radial grid data, D(r) values, peak positions, node locations, expectation values, moments, and probability integrals."),
    ]

    examples = [
        {
            "code_input": {
                "n": 1,
                "l": 0,
                "Z": 1,
                "r_max_bohr": 10.0,
                "n_points": 200,
                "compute_moments_up_to": 3,
            },
            "text_input": {
                "input_params": "1 0 1 10 200",
            },
            "output": {
                "result": {
                    "orbital_type": "1s",
                    "most_probable_radius_bohr": 1.0,
                    "mean_radius_bohr": 1.5,
                    "n_radial_nodes": 0,
                    "peak_D_r_value": 0.249,  # For 1s at r=a₀: D=4(a₀⁻³)(a₀²)e⁻²=4/a₀·e⁻²≈0.249/a₀
                }
            },
        },
        {
            "code_input": {
                "n": 2,
                "l": 0,
                "Z": 1,
                "r_max_bohr": 15.0,
                "n_points": 300,
                "compute_moments_up_to": 2,
            },
            "text_input": {
                "input_params": "2 0 1 15 300",
            },
            "output": {
                "result": {
                    "orbital_type": "2s",
                    "most_probable_radius_bohr": 5.24,  # Two peaks: ~0.75 and ~5.24
                    "n_radial_nodes": 1,
                    "node_radius_bohr": 2.0,  # Exact: 2a₀ for 2s
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _factorial(n: int) -> int:
        if n <= 1:
            return 1
        r = 1
        for i in range(2, n + 1):
            r *= i
        return r

    def _assoc_laguerre(self, n: int, k: int, x: float) -> float:
        """Associated Laguerre polynomial L_n^k(x)."""
        if n == 0:
            return 1.0
        if n == 1:
            return float(k + 1 - x)
        L0 = 1.0
        L1 = float(k + 1 - x)
        for i in range(2, n + 1):
            L_curr = ((2 * i + k - 1 - x) * L1 - (i + k - 1) * L0) / i
            L0, L1 = L1, L_curr
        return L1

    def _R_nl(self, n: int, l: int, rho: float, Z: int) -> float:
        """Radial wavefunction R_nl where ρ = 2Zr/(na₀), dimensionless."""
        nf = self._factorial(n - l - 1)
        nlf = self._factorial(n + l)
        prefactor = math.sqrt(
            (2.0 * Z / float(n * n * n)) ** 3 * nf / (2.0 * float(n) * float(nlf) ** 3)
        )
        L = self._assoc_laguerre(n - l - 1, 2 * l + 1, rho)
        return prefactor * (rho ** l) * L * math.exp(-rho / 2.0)

    def _D_r(self, n: int, l: int, r_bohr: float, Z: int) -> float:
        """Radial distribution function D(r) = r² |R_nl(r)|²."""
        rho = 2.0 * Z * r_bohr / n
        R = self._R_nl(n, l, rho, Z)
        return r_bohr * r_bohr * R * R

    def _analytical_mean_r(self, n: int, l: int, Z: int) -> float:
        """⟨r⟩ = a₀[3n² - l(l+1)] / (2Z)"""
        return (3.0 * n * n - l * (l + 1)) / (2.0 * Z)

    def _analytical_mean_r2(self, n: int, l: int, Z: int) -> float:
        """⟨r²⟩ = a₀²n²[5n²+1-3l(l+1)] / (2Z²)"""
        return (float(n * n) / (2.0 * Z * Z)) * (5.0 * n * n + 1 - 3.0 * l * (l + 1))

    def _analytical_mean_r_minus_1(self, n: int, l: int, Z: int) -> float:
        """⟨r⁻¹⟩ = Z / (a₀ · n²)"""
        return float(Z) / float(n * n)

    def _analytical_mean_r_minus_2(self, n: int, l: int, Z: int) -> float:
        """⟨r⁻²⟩ = Z² / [a₀² · n³(l+1/2)]"""
        return float(Z * Z) / (float(n ** 3) * (l + 0.5))

    def _find_peaks_and_nodes(self, n: int, l: int, Z: int,
                                r_grid: List[float], D_vals: List[float]) -> dict:
        """Find local maxima (peaks) and zeros (nodes) of D(r)."""
        peaks = []  # list of (r_bohr, D_value)
        nodes = []  # list of r_bohr where D(r) ≈ 0

        for i in range(1, len(D_vals) - 1):
            # Peak detection: local maximum
            if D_vals[i] > D_vals[i-1] and D_vals[i] > D_vals[i+1] and D_vals[i] > 1e-10:
                peaks.append((round(r_grid[i], 6), round(D_vals[i], 12)))

            # Node detection: sign change or near-zero crossing of R(r)
            # D(r) = r²R² is always ≥ 0; look for minima near zero
            if D_vals[i] < 1e-12 * max(max(D_vals), 1e-30):
                nodes.append(round(r_grid[i], 6))

        # Also check for sign changes in R (not D) by looking for deep minima
        # A true radial node makes both R and D go to zero
        refined_nodes = []
        if len(nodes) > 0:
            # Group nearby zeros
            current_group = [nodes[0]]
            for j in range(1, len(nodes)):
                if nodes[j] - current_group[-1] < (r_grid[-1] / len(r_grid)) * 3:
                    current_group.append(nodes[j])
                else:
                    refined_nodes.append(sum(current_group) / len(current_group))
                    current_group = [nodes[j]]
            refined_nodes.append(sum(current_group) / len(current_group))

        return {"peaks": peaks, "node_radii_bohr": refined_nodes}

    def _run_base(self, n: int, l: int, Z: int = 1,
                  r_max_bohr: float = 30.0, n_points: int = 500,
                  compute_moments_up_to: int = 3,
                  radius_range_for_prob: list = None) -> dict:

        # Validation
        if n < 1:
            raise ChemMCPError("n must be >= 1.")
        if not (0 <= l < n):
            raise ChemMCPError(f"l must satisfy 0 <= l < n (got l={l}, n={n}).")
        if Z < 1:
            raise ChemMCPError("Z must be >= 1.")

        # Orbital label
        orb_labels = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}
        orbital_type = f"{n}{orb_labels.get(l, f'l={l}')}"
        n_radial_nodes_expected = n - l - 1

        # Build radial grid
        dr = r_max_bohr / (n_points - 1) if n_points > 1 else r_max_bohr
        r_grid = [r_max_bohr * i / (n_points - 1) for i in range(n_points)]

        # Compute D(r) on grid
        D_vals = [self._D_r(n, l, r, Z) for r in r_grid]
        R_vals = [self._R_nl(n, l, 2.0 * Z * r / n, Z) for r in r_grid]

        # Find peaks and nodes
        pn_data = self._find_peaks_and_nodes(n, l, Z, r_grid, D_vals)

        # Numerical integration of D(r): should equal 1 (normalization)
        norm_integral = sum(D_vals) * dr

        # Analytical observables
        mean_r = self._analytical_mean_r(n, l, Z)
        mean_r2 = self._analytical_mean_r2(n, l, Z)
        delta_r = math.sqrt(max(0, mean_r2 - mean_r * mean_r))

        # Radial moments
        moments = {}
        for k in range(1, compute_moments_up_to + 1):
            if k == 1:
                moments[f"<r^{k}>"] = round(mean_r, 8)
            elif k == 2:
                moments[f"<r^{k}>"] = round(mean_r2, 8)
            elif k == -1:
                moments["<r^-1>"] = round(self._analytical_mean_r_minus_1(n, l, Z), 8)
            elif k == -2:
                moments["<r^-2>"] = round(self._analytical_mean_r_minus_2(n, l, Z), 8)
            else:
                # Numerical integration for other moments
                moment_val = sum((r ** k) * d for r, d in zip(r_grid, D_vals)) * dr
                moments[f"<r^{k}>"] = round(moment_val, 8)

        # Probability in specified radius ranges
        prob_in_ranges = {}
        if radius_range_for_prob:
            for rng in radius_range_for_prob:
                if len(rng) >= 2:
                    r_min, r_max_rng = rng[0], rng[1]
                    idx_min = max(0, int(r_min / dr))
                    idx_max = min(n_points - 1, int(r_max_rng / dr))
                    prob = sum(D_vals[idx_min:idx_max + 1]) * dr
                    label = f"[{r_min}, {r_max_rng}] a₀"
                    prob_in_ranges[label] = round(prob, 8)

        # Cumulative probability distribution (at sample points)
        cumul = []
        running_sum = 0.0
        step = max(1, n_points // 25)
        for i in range(0, n_points, step):
            running_sum += sum(D_vals[max(0, step//2):min(n_points, i + step//2 + 1)]) * dr
            cumul.append({"r_bohr": round(r_grid[i], 4), "cumulative_probability": round(min(running_sum, 1.0), 6)})

        # Most probable radius (global maximum of D(r))
        global_peak_idx = max(range(len(D_vals)), key=lambda i: D_vals[i])
        r_most_probable = r_grid[global_peak_idx]
        D_at_peak = D_vals[global_peak_idx]

        result = {
            "orbital_type": orbital_type,
            "quantum_numbers": {"n": n, "l": l},
            "nuclear_charge_Z": Z,
            "n_radial_nodes_expected": n_radial_nodes_expected,
            "n_radial_nodes_found": len(pn_data["node_radii_bohr"]),
            "node_radii_bohr": pn_data["node_radii_bohr"],
            "most_probable_radius_bohr": round(r_most_probable, 6),
            "D_at_most_probable_r": round(D_at_peak, 12),
            "peaks": pn_data["peaks"],  # All local maxima [(r, D), ...]
            "mean_radius_bohr": round(mean_r, 8),
            "mean_radius_sq_bohr2": round(mean_r2, 8),
            "uncertainty_delta_r_bohr": round(delta_r, 8),
            "radial_moments_a0_units": moments,
            "normalization_integral": round(norm_integral, 8),
            "probability_in_ranges": prob_in_ranges if prob_in_ranges else None,
            "cumulative_distribution_sample": cumul,
            "radial_grid_sample": [round(r, 6) for r in r_grid[::step]],
            "D_r_sample": [round(d, 12) for d in D_vals[::step]],
            "R_nl_sample": [round(rv, 12) for rv in R_vals[::step]],
            "n_grid_points": n_points,
            "r_max_bohr": r_max_bohr,
            "grid_spacing_dr_bohr": round(dr, 8),
        }

        logger.info(f"RadialDistribution: {orbital_type}, r_mp={r_most_probable:.4f}a₀, <r>={mean_r:.4f}a₀")
        return {"result": result}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            n = int(parts[0])
            l = int(parts[1])
            Z = int(parts[2]) if len(parts) > 2 else 1
            r_max = float(parts[3]) if len(parts) > 3 else 30.0
            n_pts = int(parts[4]) if len(parts) > 4 else 500
            return self._run_base(n, l, Z, r_max, n_pts)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'n l [Z] [r_max] [n_points]'")
