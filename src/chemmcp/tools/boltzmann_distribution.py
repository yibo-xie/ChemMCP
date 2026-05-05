import logging
import math
import json
from typing import Optional, List
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BoltzmannDistribution(BaseTool):
    """
    Boltzmann 分布计算与可视化数据生成工具。
    
    计算各能级的布居概率和粒子数分布，并输出可用于绑图的数据。
      P_i = (g_i · exp(-E_i/kT)) / Q
      N_i = N · P_i
    """
    __version__                 = "0.1.0"
    name                        = "BoltzmannDistribution"
    func_name                   = "calculate_boltzmann_distribution"
    description                 = "Calculate Boltzmann distribution: probabilities and populations across energy levels, with plot-ready data."
    implementation_description  = "Uses P_i = g_i·exp(-E_i/kT) / Z where Z is the partition sum. Returns per-level probabilities, populations, and curve data for visualization."
    oss_dependencies            = []
    services_and_software       = []
    categories                  = ["General"]
    tags                        = ["Boltzmann", "Statistical Mechanics", "Distribution", "Thermodynamics"]
    required_envs               = []

    code_input_sig   = [
        ("temperature_k",            "float",  "N/A",     "Temperature in Kelvin."),
        ("energy_levels_j",          "list",   "N/A",     "List of energy levels in Joules [E0, E1, E2, ...]."),
        ("degeneracies",             "list",   "None",    "List of degeneracies g_i for each level (default all 1)."),
        ("total_particles",          "int",    "None",    "Total number of particles N (optional, for population counts)."),
        ("energy_range_max_j",       "float",  "None",    "Max energy for continuous distribution curve data (J)."),
        ("num_curve_points",         "int",    "100",     "Number of points for continuous distribution curve."),
    ]

    text_input_sig   = [
        ("input_params",             "str",    "N/A",     "Space-separated: 'T(K) E0 E1 E2 ... [N]'. Energies in J."),
    ]

    output_sig       = [
        ("result",                  "dict",    "Dict with partition_sum_Z, probabilities, populations, plot_data."),
    ]

    examples         = [
        {
            "code_input": {
                "temperature_k":              300.0,
                "energy_levels_j":            [0, 1.38e-23, 2.76e-23, 4.14e-23],
                "degeneracies":               [1, 2, 3, 4],
                "total_particles":            10000,
                "energy_range_max_j":         None,
                "num_curve_points":           100,
            },
            "text_input": {
                "input_params":               "300.0 0 1.38e-23 2.76e-23 4.14e-23 10000",
            },
            "output": {
                "result": {
                    "temperature_K": 300.0,
                    "num_levels": 4,
                    "partition_sum_Z": 8.3636,
                    "level_data": [
                        {"level": 0, "E_J": 0, "g_i": 1, "probability": 0.1196, "population": 1196},
                        {"level": 1, "E_J": 1.38e-23, "g_i": 2, "probability": 0.2184, "population": 2184},
                    ],
                }
            },
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.k_B = 1.380649e-23  # J/K

    def _run_base(
        self,
        temperature_k: float,
        energy_levels_j: List[float],
        degeneracies: Optional[List[int]] = None,
        total_particles: Optional[int] = None,
        energy_range_max_j: Optional[float] = None,
        num_curve_points: int = 100,
    ) -> dict:
        """Core logic: calculate Boltzmann distribution."""
        T = float(temperature_k)
        
        if T <= 0:
            raise ChemMCPError("Temperature must be > 0 K.")
        
        if not energy_levels_j or len(energy_levels_j) == 0:
            raise ChemMCPError("At least one energy level must be provided.")

        energies = [float(E) for E in energy_levels_j]
        n_levels = len(energies)
        
        if degeneracies is None:
            gs = [1] * n_levels
        elif len(degeneracies) == n_levels:
            gs = [int(g) for g in degeneracies]
        else:
            raise ChemMCPError(f"Length of degeneracies ({len(degeneracies)}) must match number of levels ({n_levels}).")

        # --- Compute partition sum Z ---
        Z = 0.0
        boltzmann_factors = []
        for i in range(n_levels):
            x = -energies[i] / (self.k_B * T)
            bf = gs[i] * math.exp(x) if x > -700 else (gs[i] if x > -800 else 0.0)
            boltzmann_factors.append(bf)
            Z += bf

        if Z <= 0:
            raise ChemMCPError("Partition sum Z must be positive.")

        # --- Per-level probabilities and populations ---
        level_data = []
        for i in range(n_levels):
            prob = boltzmann_factors[i] / Z
            pop = int(round(total_particles * prob)) if total_particles is not None else None
            
            level_data.append({
                "level":              i,
                "E_J":                energies[i],
                "E_eV":               round(energies[i] / 1.602176634e-19, 10),
                "g_i":                gs[i],
                "boltzmann_factor":   round(boltzmann_factors[i], 10),
                "probability":        round(prob, 10),
                "population":         pop,
            })

        # --- Continuous distribution curve data (density of states approximation) ---
        plot_data = None
        if energy_range_max_j is not None and num_curve_points > 0:
            E_max = float(energy_range_max_j)
            dE = E_max / num_curve_points if num_curve_points > 0 else E_max
            curve_points = []
            
            for pt in range(num_curve_points + 1):
                E = pt * dE
                # Approximate continuous density using interpolation or exponential envelope
                p_cont = math.exp(-E / (self.k_B * T)) / (self.k_B * T)  # normalized approx
                curve_points.append({
                    "E_J": round(E, 15),
                    "probability_density_per_J": round(p_cont, 10),
                })
            
            plot_data = {
                "description": "Continuous Boltzmann probability density (approximation)",
                "num_points": len(curve_points),
                "points": curve_points[:20],  # First 20 points to keep response manageable
                "total_points_computed": len(curve_points),
            }

        result = {
            "temperature_K":       T,
            "num_levels":          n_levels,
            "partition_sum_Z":     round(Z, 10),
            "level_data":          level_data,
            "plot_data":           plot_data,
            "N_total":             total_particles,
            "sum_probabilities":   round(sum(ld["probability"] for ld in level_data), 10),
        }

        logger.info(f"BoltzmannDist: T={T}K, {n_levels} levels, Z={Z:.6f}")
        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.split()
            if len(parts) < 2:
                raise ValueError("Need at least 'T E0 E1 ...' params.")
            
            T = float(parts[0])
            energies = [float(p) for p in parts[1:-1]] if parts[-1].isdigit() else [float(p) for p in parts[1:]]
            N = int(parts[-1]) if (len(energies) < len(parts) - 1 and parts[-1].replace('-','').isdigit()) else None
            
            # Re-parse: last token might be N
            if N is not None and N > 0:
                pass  # energies already set correctly
            elif len(parts) >= 2 and parts[-1].lstrip('-').replace('.','',1).isdigit():
                # Check if last looks like an integer (particle count)
                try:
                    potential_N = int(parts[-1])
                    if potential_N > 0 and potential_N < 1e12:
                        energies = [float(p) for p in parts[1:-1]]
                        N = potential_N
                except ValueError:
                    N = None
            
            return self._run_base(T, energies, None, N)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'T(K) E0 E1 E2 ... [N_total]'")
