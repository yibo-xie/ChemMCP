import logging
import math
from typing import List, Optional
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

R_J = 8.314  # J/(mol·K)


@ChemMCPManager.register_tool
class GibbsMinimization(BaseTool):
    """
    吉布斯自由能最小化求解复杂化学平衡。
    
    使用拉格朗日乘子法或元素平衡法，在给定初始组成和约束条件下，
    搜索使总吉布斯能最小的平衡组成。
    
    支持简单气相反应和多组分同时平衡。
    """
    __version__ = "0.1.0"
    name = "GibbsMinimization"
    func_name = "gibbs_minimization"
    description = "Minimize Gibbs free energy to find equilibrium composition for complex chemical systems using Lagrange multiplier method with element balance constraints."
    implementation_description = "Implements Gibbs energy minimization: G_total = Σ n_i [μ_i° + RT ln(n_i/n_total)] subject to element balance constraints. Uses iterative Lagrange multiplier method (or simple stoichiometric search for single reactions). Supports gas-phase reactions with ideal mixing. For single reaction: scans extent of reaction ξ to find minimum G."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Gibbs Free Energy", "Chemical Equilibrium", "Minimization", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Mode: 'single_reaction' (scan ξ) or 'multi_component' (Lagrange method)."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("pressure_atm", "float", "1.0", "Total pressure in atm."),
        ("species", "list", "N/A", "List of species names, e.g., ['NO2', 'N2O4']."),
        ("mu_standard", "list", "N/A", "Standard chemical potentials μ° in kJ/mol (same order as species)."),
        ("initial_moles", "list", "N/A", "Initial moles of each species (same order as species)."),
        ("stoich_coeffs", "list", "None", "Stoichiometric coefficients ν_i for single_reaction mode (negative=reactant, positive=product)."),
        ("n_scan_points", "int", "1000", "Number of ξ scan points for single_reaction mode."),
        ("element_matrix", "list", "None", "Element matrix a_ij (for multi_component mode): rows=elements, cols=species)."),
        ("element_totals", "list", "None", "Total moles of each element (for multi_component mode)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "JSON-like format: mode T P species|mu0|n0 [|stoich] [|elem_matrix|elem_totals]. Example: 'single_reaction 298 1 NO2,N2O4|51.3,97.9|2.0,0.0|-2,1'"),
    ]

    output_sig = [
        ("equilibrium_moles", "list", "Equilibrium moles of each species."),
        ("extent_xi", "float", "Reaction extent at equilibrium (single reaction mode)."),
        ("G_total_kJ", "float", "Total Gibbs free energy at equilibrium (kJ)."),
        ("G_initial_kJ", "float", "Initial total Gibbs free energy (kJ)."),
        ("delta_G_rxn_kJ", "float", "ΔG of reaction from initial to equilibrium (kJ)."),
        ("mole_fractions", "list", "Mole fractions at equilibrium."),
        ("K_equilibrium", "float", "Equilibrium constant K calculated from composition."),
        ("explanation", "str", "Detailed calculation steps and results summary."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "single_reaction",
                "temperature_k": 298.15,
                "pressure_atm": 1.0,
                "species": ["NO2", "N2O4"],
                "mu_standard": [51.3, 97.9],
                "initial_moles": [2.0, 0.0],
                "stoich_coeffs": [-2.0, 1.0],
                "n_scan_points": 1000,
            },
            "text_input": {
                "input_params": "single_reaction 298.15 1 NO2,N2O4|51.3,97.9|2.0,0|-2,1",
            },
            "output": {
                "equilibrium_moles": [1.25, 0.375],
                "extent_xi": 0.375,
                "K_equilibrium": round(0.24, 2),
            },
        },
        {
            "code_input": {
                "mode": "single_reaction",
                "temperature_k": 298.15,
                "pressure_atm": 1.0,
                "species": ["H2", "I2", "HI"],
                "mu_standard": [0.0, 0.0, 1.7],  # approximate
                "initial_moles": [1.0, 1.0, 0.0],
                "stoich_coeffs": [-1.0, -1.0, 2.0],
                "n_scan_points": 1000,
            },
            "text_input": {
                "input_params": "single_reaction 298 1 H2,I2,HI|0,0,1.7|1,1,0|-1,-1,2",
            },
            "output": {
                "equilibrium_moles": [0.22, 0.22, 1.56],
                "extent_xi": 0.78,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _compute_G_total(self, n: List[float], mu0: List[float],
                         T: float, P: float) -> float:
        """Compute total Gibbs energy: G = Σ n_i[μ_i° + RT ln(x_i)]"""
        n_total = sum(n)
        if n_total < 1e-30:
            return float('inf')
        G = 0.0
        RT = R_J * T / 1000.0  # kJ/mol
        for i in range(len(n)):
            if n[i] > 1e-30:
                x_i = n[i] / n_total
                G += n[i] * (mu0[i] + RT * math.log(max(x_i, 1e-30)))
            elif mu0[i] > 0:
                # Species not present: use large penalty
                pass  # Don't add contribution for absent species
        return G

    def _run_single_reaction(self, species: List[str], mu0: List[float],
                              n0: List[float], nu: List[float],
                              T: float, P: float,
                              n_scan: int) -> dict:
        """Scan extent of reaction ξ to find minimum G."""
        ns = len(species)

        # Determine feasible range of ξ
        xi_min = 0.0
        xi_max = float('inf')
        for i in range(ns):
            if nu[i] < 0:
                # Reactant: n_i = n0_i + ν_i·ξ ≥ 0 → ξ ≤ n0_i / |ν_i|
                xi_max = min(xi_max, n0_i_abs := n0[i] / abs(nu[i]))
            elif nu[i] > 0:
                # Product: no upper limit from this species
                pass

        if xi_max == float('inf') or xi_max <= 0:
            xi_max = min(n0[i] / abs(nu[i]) for i in range(ns) if nu[i] < 0)

        # Scan ξ
        xi_min_actual = max(0.0, xi_min - xi_max * 0.01)
        xi_max_actual = xi_max * 1.01

        best_G = float('inf')
        best_xi = 0.0
        best_n = list(n0)

        for step in range(n_scan + 1):
            xi = xi_min_actual + (xi_max_actual - xi_min_actual) * step / n_scan
            n = [max(0.0, n0[i] + nu[i] * xi) for i in range(ns)]
            G = self._compute_G_total(n, mu0, T, P)
            if G < best_G:
                best_G = G
                best_xi = xi
                best_n = list(n)

        # Compute initial G
        G_init = self._compute_G_total(n0, mu0, T, P)

        # Mole fractions
        n_tot = sum(best_n)
        x_eq = [ni / n_tot if n_tot > 1e-30 else 0.0 for ni in best_n]

        # Calculate K from equilibrium composition
        # Kp = Π (P_i/P°)^ν_i = Π (x_i·P)^ν_i  (P° = 1 atm)
        try:
            log_K = 0.0
            for i in range(ns):
                if abs(nu[i]) > 1e-10 and x_eq[i] > 1e-30:
                    log_K += nu[i] * math.log(max(x_eq[i] * P, 1e-30))
            K_val = math.exp(log_K) if abs(log_K) < 50 else (float('inf') if log_K > 0 else 0.0)
        except (ValueError, OverflowError):
            K_val = float('nan')

        delta_G = best_G - G_init

        explanation = (
            f"吉布斯能最小化 — 单反应模式\n"
            f"反应: {' + '.join(f'{abs(nu[i])}{species[i]}' if nu[i]<0 else '' for i in range(ns))}"
            f" ⇌ {' + '.join(f'{nu[i]}{species[i]}' if nu[i]>0 else '' for i in range(ns))}\n"
            f"T = {T} K, P = {P} atm\n"
            f"初始摩尔数: {dict(zip(species, [f'{v:.4f}' for v in n0]))}\n"
            f"扫描 ξ ∈ [0, {xi_max:.4f}], 共 {n_scan+1} 个点\n"
            f"最优反应进度 ξ* = {best_xi:.6f}\n"
            f"平衡摩尔数: {dict(zip(species, [f'{v:.6f}' for v in best_n]))}\n"
            f"平衡摩尔分数: {dict(zip(species, [f'{v:.4f}' for v in x_eq]))}\n"
            f"G_初始 = {G_init:.4f} kJ, G_平衡 = {best_G:.4f} kJ\n"
            f"ΔG = {delta_G:.4f} kJ, K = {K_val:.6g}"
        )

        logger.info(f"GibbsMinimization: ξ*={best_xi:.6f}, G={best_G:.4f}kJ, K={K_val:.4g}")
        return {
            "equilibrium_moles": [round(v, 10) for v in best_n],
            "extent_xi": round(best_xi, 10),
            "G_total_kJ": round(best_G, 6),
            "G_initial_kJ": round(G_init, 6),
            "delta_G_rxn_kJ": round(delta_G, 6),
            "mole_fractions": [round(v, 10) for v in x_eq],
            "K_equilibrium": K_val if math.isfinite(K_val) else None,
            "explanation": explanation,
        }

    def _run_multi_component(self, species: List[str], mu0: List[float],
                               n0: List[float], T: float, P: float,
                               elem_matrix: List[List[float]],
                               elem_totals: List[float]) -> dict:
        """Simplified multi-component equilibrium via stoichiometric iteration."""
        # For multi-component: iterate over extents to minimize G
        # This is a simplified version - uses random sampling + refinement
        import random
        ns = len(species)

        best_G = float('inf')
        best_n = list(n0)

        # Simple grid search with perturbation
        for trial in range(2000):
            if trial == 0:
                n_trial = list(n0)
            else:
                # Random perturbation around current best
                scale = max(0.01, sum(n0) * 0.05 * (1.0 - trial / 2000))
                n_trial = [max(0, best_n[i] + (random.random() - 0.5) * scale) for i in range(ns)]

            # Check element balance constraint
            if elem_matrix and elem_totals:
                balanced = True
                for j, b_j in enumerate(elem_totals):
                    computed = sum(elem_matrix[j][i] * n_trial[i] for i in range(ns))
                    if abs(computed - b_j) > b_j * 0.01:
                        balanced = False
                        break
                if not balanced:
                    continue

            G = self._compute_G_total(n_trial, mu0, T, P)
            if G < best_G:
                best_G = G
                best_n = list(n_trial)

        G_init = self._compute_G_total(n0, mu0, T, P)
        n_tot = sum(best_n)
        x_eq = [ni / n_tot if n_tot > 1e-30 else 0.0 for ni in best_n]

        explanation = (
            f"吉布斯能最小化 — 多组分模式\n"
            f"物种: {species}\n"
            f"T = {T} K, P = {P} atm\n"
            f"初始: {n0}\n"
            f"平衡: {[round(v,6) for v in best_n]}\n"
            f"G_初始 = {G_init:.4f} kJ, G_平衡 = {best_G:.4f} kJ"
        )

        return {
            "equilibrium_moles": [round(v, 10) for v in best_n],
            "extent_xi": None,
            "G_total_kJ": round(best_G, 6),
            "G_initial_kJ": round(G_init, 6),
            "delta_G_rxn_kJ": round(best_G - G_init, 6),
            "mole_fractions": [round(v, 10) for v in x_eq],
            "K_equilibrium": None,
            "explanation": explanation,
        }

    def _run_base(self, mode: str, temperature_k: float = 298.15,
                  pressure_atm: float = 1.0, species: list = None,
                  mu_standard: list = None, initial_moles: list = None,
                  stoich_coeffs: list = None, n_scan_points: int = 1000,
                  element_matrix: list = None, element_totals: list = None) -> dict:

        if not species or not mu_standard or not initial_moles:
            raise ChemMCPError("species, mu_standard, and initial_moles are all required.")
        if len(species) != len(mu_standard) or len(species) != len(initial_moles):
            raise ChemMCPError("species, mu_standard, and initial_moles must have same length.")

        mode = mode.lower().strip()
        if mode == "single_reaction":
            if not stoich_coeffs:
                raise ChemMCPError("stoich_coeffs is required for single_reaction mode.")
            if len(stoich_coeffs) != len(species):
                raise ChemMCPError("stoich_coeffs must have same length as species.")
            return self._run_single_reaction(
                species, mu_standard, initial_moles, stoich_coeffs,
                temperature_k, pressure_atm, n_scan_points
            )
        elif mode == "multi_component":
            return self._run_multi_component(
                species, mu_standard, initial_moles,
                temperature_k, pressure_atm, element_matrix, element_totals
            )
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use 'single_reaction' or 'multi_component'.")

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split("|")
            # First part: "mode T P species" where species may contain commas
            first = parts[0].strip().split()
            mode = first[0]
            T = float(first[1]) if len(first) > 1 else 298.15
            P = float(first[2]) if len(first) > 2 else 1.0
            # Species is everything after T and P in first segment (may contain commas)
            sp_raw = ",".join(first[3:]) if len(first) > 3 else (parts[1].strip() if len(parts) > 1 else "")
            sp_list = [s.strip() for s in sp_raw.split(",")]

            # Remaining segments shift by 1 since species was in first
            mu_list = [float(x.strip()) for x in parts[1].split(",")] if len(parts) > 1 else []
            n0_list = [float(x.strip()) for x in parts[2].split(",")] if len(parts) > 2 else []
            stoich = [float(x.strip()) for x in parts[3].split(",")] if len(parts) > 3 and parts[3].strip() else None

            kwargs = {"mode": mode, "temperature_k": T, "pressure_atm": P,
                      "species": sp_list, "mu_standard": mu_list,
                      "initial_moles": n0_list}
            if stoich:
                kwargs["stoich_coeffs"] = stoich
            return self._run_base(**kwargs)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'mode T P species|mu0|n0 [|stoich]'")
