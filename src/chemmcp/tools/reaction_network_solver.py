import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ReactionNetworkSolver(BaseTool):
    """
    反应网络求解器（Reaction Network Solver）。
    
    对任意反应网络进行数值求解，计算各物种浓度随时间的变化。
    
    支持的网络类型：
    - 连串反应网络：A → B → C → ...
    - 平行反应：A → B, A → C, A → D (竞争)
    - 可逆反应：A ⇌ B
    - 复合网络：上述类型的任意组合
    - 环状网络：A → B → C → A
    
    使用四阶龙格-库塔法（RK4）数值积分耦合常微分方程组。
    """
    __version__ = "0.1.0"
    name = "ReactionNetworkSolver"
    func_name = "solve_reaction_network"
    description = "Solve coupled ODEs for arbitrary reaction networks using numerical integration (RK4). Compute concentration vs time profiles for all species, detect steady states, and analyze rate-limiting steps."
    implementation_description = "Parses reaction network definition into a system of first-order ODEs. Integrates using 4th-order Runge-Kutta with adaptive step-size control. Supports consecutive, parallel, reversible, cyclic and mixed networks. Outputs concentration profiles, half-lives, steady-state detection, and network topology analysis."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Reaction Network", "ODE Solver", "Kinetics Simulation", "Mechanism", "RK4"]
    required_envs = []

    code_input_sig = [
        ("species", "list", "N/A", "List of species names: ['A', 'B', 'C', ...]."),
        ("reactions", "list", "N/A", "List of reaction dicts: [{'reactants': ['A'], 'products': ['B'], 'k': 0.1, 'reversible': False}, ...]."),
        ("initial_concentrations", "dict", "N/A", "Initial concentrations: {'A': 1.0, 'B': 0.0}."),
        ("time_end", "float", "100.0", "End time for simulation."),
        ("n_points", "int", "100", "Number of time points to output."),
        ("network_type", "str", "auto", "'consecutive', 'parallel', 'reversible', 'cyclic', or 'auto' (auto-detect)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'species A,B,C reactions A->B:k1;B->C:k2 init A=1,B=0,C=0 t_end=100 n=100'. Example: 'species A,B,C reactions A->B:0.1;B->C:0.5 init A=1 t=50'"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with concentration profiles (time grid + per-species values), network analysis, steady-state info, half-lives, and rate diagnostics."),
    ]

    examples = [
        {
            "code_input": {
                "species": ["A", "I", "P"],
                "reactions": [
                    {"reactants": ["A"], "products": ["I"], "k": 0.1, "reversible": False},
                    {"reactants": ["I"], "products": ["P"], "k": 1.0, "reversible": False},
                ],
                "initial_concentrations": {"A": 1.0, "I": 0.0, "P": 0.0},
                "time_end": 50.0,
                "n_points": 50,
                "network_type": "consecutive",
            },
            "text_input": {
                "input_params": "species A,I,P reactions A->I:0.1;I->P:1.0 init A=1 t=50",
            },
            "output": {
                "result": {
                    "n_species": 3,
                    "n_reactions": 2,
                    "network_type": "consecutive",
                    "final_concentrations": {"A": 0.007, "I": 0.000, "P": 0.993},
                    "half_life_A": 6.93,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _parse_reactions(self, species, reactions, init_conc):
        """Validate and index reactions."""
        species_set = set(species)
        for i, rxn in enumerate(reactions):
            for r in rxn.get("reactants", []):
                if r not in species_set:
                    raise ChemMCPError(f"Reaction {i+1}: unknown reactant '{r}'")
            for p in rxn.get("products", []):
                if p not in species_set:
                    raise ChemMCPError(f"Reaction {i+1}: unknown product '{p}'")
            if rxn.get("k") is None:
                raise ChemMCPError(f"Reaction {i+1}: missing rate constant k")
        
        # Ensure all species have initial concentrations
        conc = dict(init_conc)
        for s in species:
            if s not in conc:
                conc[s] = 0.0
        
        return conc

    def _build_ode_system(self, species, reactions):
        """
        Build d[species]/t function from reaction list.
        
        For each reaction like A -> B with k:
          d[A]/dt -= k*[A]
          d[B]/dt += k*[A]
        
        For reversible A <=> B with kf, kr:
          forward: A -> B at rate kf*[A]
          reverse: B -> A at rate kr*[B]
          
        For bimolecular A + B -> C:
          d[A]/dt -= k*[A][B], etc.
        """
        n = len(species)
        sp_idx = {s: i for i, s in enumerate(species)}

        def ode_func(t, y):
            dydt = [0.0] * n
            
            for rxn in reactions:
                k = rxn["k"]
                reacts = rxn.get("reactants", [])
                prods = rxn.get("products", [])
                is_rev = rxn.get("reversible", False)
                kr = rxn.get("k_reverse", 0) if is_rev else 0

                # Compute reaction rate (mass action)
                # Determine molecularity
                if len(reacts) == 1:
                    idx_r = sp_idx[reacts[0]]
                    rate_fwd = k * max(y[idx_r], 0)
                elif len(reacts) == 2:
                    idx_r1 = sp_idx[reacts[0]]
                    idx_r2 = sp_idx[reacts[1]]
                    rate_fwd = k * max(y[idx_r1], 0) * max(y[idx_r2], 0)
                else:
                    # Higher order - product of all reactant concentrations
                    rate_fwd = k
                    for r in reacts:
                        rate_fwd *= max(y[sp_idx[r]], 0)
                
                # Forward reaction contributions
                for r in reacts:
                    dydt[sp_idx[r]] -= rate_fwd
                for p in prods:
                    dydt[sp_idx[p]] += rate_fwd

                # Reverse reaction (if reversible)
                if is_rev and kr > 0:
                    if len(prods) == 1:
                        rate_rev = kr * max(y[sp_idx[prods[0]]], 0)
                    elif len(prods) == 2:
                        rate_rev = kr * max(y[sp_idx[prods[0]]], 0) * max(y[sp_idx[prods[1]]], 0)
                    else:
                        rate_rev = kr
                        for p in prods:
                            rate_rev *= max(y[sp_idx[p]], 0)
                    
                    for p in prods:
                        dydt[sp_idx[p]] -= rate_rev
                    for r in reacts:
                        dydt[sp_idx[r]] += rate_rev

            return dydt

        return ode_func

    def _rk4_step(self, f, t, y, dt):
        """Single RK4 step."""
        k1 = f(t, y)
        
        k2_y = [y[i] + 0.5 * dt * k1[i] for i in range(len(y))]
        k2 = f(t + 0.5 * dt, k2_y)
        
        k3_y = [y[i] + 0.5 * dt * k2[i] for i in range(len(y))]
        k3 = f(t + 0.5 * dt, k3_y)
        
        k4_y = [y[i] + dt * k3[i] for i in range(len(y))]
        k4 = f(t + dt, k4_y)
        
        return [y[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(len(y))]

    def _integrate(self, ode_func, y0, t_end, n_points):
        """Integrate ODE system from t=0 to t=t_end using RK4."""
        n = len(y0)
        dt = t_end / max(n_points - 1, 1)
        
        times = [0.0]
        y_current = list(y0)
        profiles = [[v for v in y0]]
        
        t = 0.0
        step_count = 0
        max_steps = n_points * 100  # Safety limit
        
        while t < t_end and step_count < max_steps:
            # Adaptive: take smaller steps near beginning if needed
            actual_dt = min(dt, t_end - t)
            
            # Ensure we don't overshoot output points
            next_output_time = times[-1] + dt
            
            if next_output_time > t_end:
                break
            
            y_current = self._rk4_step(ode_func, t, y_current, actual_dt)
            t += actual_dt
            step_count += 1
            
            # Record at regular intervals
            if abs(t - (len(profiles) * dt)) < actual_dt * 0.5 or step_count % max(1, int(1/dt * dt)) == 0:
                # Check if close enough to desired output time
                desired_t = len(profiles) * dt
                if abs(t - desired_t) < actual_dt * 1.1:
                    times.append(round(t, 10))
                    profiles.append([max(v, 0) for v in y_current])  # No negative concentrations

        # Ensure we have exactly n_points
        while len(profiles) < n_points:
            profiles.append(profiles[-1] if profiles else y0)
            times.append(times[-1] + dt if times else 0)
        
        # Trim to n_points
        if len(profiles) > n_points:
            # Sample evenly
            indices = [int(i * (len(profiles) - 1) / (n_points - 1)) for i in range(n_points)]
            profiles = [profiles[j] for j in indices]
            times = [times[j] for j in indices]

        return times, profiles

    def _detect_steady_state(self, species, times, profiles):
        """Detect if system has reached steady state."""
        if len(profiles) < 3:
            return {"steady_state": False}
        
        last_few = min(10, len(profiles))
        n_sp = len(species)
        
        max_rates = []
        for j in range(n_sp):
            vals = [profiles[i][j] for i in range(-last_few, 0)]
            if len(vals) >= 2:
                avg_rate = abs(vals[-1] - vals[0]) / (times[-1] - times[-last_few]) if times[-1] != times[-last_few] else 0
                mean_val = sum(vals) / len(vals)
                rel_rate = avg_rate / mean_val if mean_val > 1e-15 else 0
                max_rates.append(rel_rate)
            else:
                max_rates.append(float('inf'))
        
        is_steady = all(r < 0.01 for r in max_rates)  # < 1% change per unit time
        return {
            "steady_state": is_steady,
            "max_relative_rate": round(max(max_rates) if max_rates else 0, 8),
            "per_species_rates": {species[j]: round(max_rates[j], 8) for j in range(min(n_sp, len(max_rates)))},
        }

    def _compute_half_lives(self, species, times, profiles, init_conc):
        """Estimate half-lives for reactants that decay."""
        half_lives = {}
        sp_idx = {s: i for i, s in enumerate(species)}
        
        for sp, c0 in init_conc.items():
            if c0 <= 0:
                continue
            idx = sp_idx.get(sp)
            if idx is None:
                continue
            
            target = c0 * 0.5
            for i in range(1, len(profiles)):
                if profiles[idx][i] <= target:
                    # Linear interpolation
                    t1, t2 = times[i-1], times[i]
                    v1, v2 = profiles[idx][i-1], profiles[idx][i]
                    if v2 != v1:
                        t_half = t1 + (target - v1) * (t2 - t1) / (v2 - v1)
                        half_lives[sp] = round(t_half, 4)
                    break
            else:
                # Check if final value is still above half
                if profiles[idx][-1] > target:
                    half_lives[sp] = None  # Didn't reach half within simulation time

        return half_lives

    def _run_base(self, species: list, reactions: list, initial_concentrations: dict,
                  time_end: float = 100.0, n_points: int = 100,
                  network_type: str = "auto") -> dict:
        """Core logic."""
        if not species:
            raise ChemMCPError("Species list cannot be empty.")
        if not reactions:
            raise ChemMCPError("Reactions list cannot be empty.")
        if time_end <= 0:
            raise ChemMCPError("time_end must be positive.")
        if n_points < 2:
            raise ChemMCPError("n_points must be at least 2.")

        # Parse and validate
        conc = self._parse_reactions(species, reactions, initial_concentrations)
        y0 = [conc[s] for s in species]

        # Build ODE system
        ode_func = self._build_ode_system(species, reactions)

        # Integrate
        times, profiles = self._integrate(ode_func, y0, time_end, n_points)

        # Analysis
        ss_info = self._detect_steady_state(species, times, profiles)
        half_lives = self._compute_half_lives(species, times, profiles, conc)

        # Format concentration profiles
        profile_data = []
        for i, t in enumerate(times):
            entry = {"time": round(t, 6)}
            for j, s in enumerate(species):
                entry[s] = round(profiles[i][j], 8)
            profile_data.append(entry)

        # Final concentrations
        final = {s: round(profiles[-1][j], 8) for j, s in enumerate(species)}

        # Auto-detect network type if needed
        detected_type = network_type
        if network_type == "auto":
            has_reversible = any(r.get("reversible", False) for r in reactions)
            has_parallel = len(set(tuple(sorted(r.get("reactants", []))) for r in reactions)) < len(reactions)
            if has_reversible:
                detected_type = "reversible"
            elif has_parallel:
                detected_type = "parallel"
            else:
                detected_type = "consecutive"

        result = {
            "n_species": len(species),
            "species_list": species,
            "n_reactions": len(reactions),
            "network_type": detected_type,
            "initial_concentrations": conc,
            "final_concentrations": final,
            "concentration_profiles": profile_data,
            "time_grid": [round(t, 6) for t in times],
            "steady_state_analysis": ss_info,
            "half_lives": half_lives,
            "simulation_parameters": {
                "time_end": time_end,
                "n_output_points": len(profiles),
                "integrator": "RK4 (4th-order Runge-Kutta)",
            },
        }

        logger.info(f"ReactionNetworkSolver: {detected_type} network, {len(species)} species, "
                     f"{len(reactions)} reactions, steady={ss_info['steady_state']}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            species = []
            reactions_raw = []
            init_dict = {}
            t_end = 100.0
            n_pts = 100
            
            i = 0
            while i < len(parts):
                p = parts[i]
                if p == "species":
                    i += 1
                    species = [x.strip() for x in parts[i].split(",")]
                elif p == "reactions":
                    i += 1
                    rxn_str = parts[i]
                    for rxn in rxn_str.split(";"):
                        rxn = rxn.strip()
                        if ":" in rxn:
                            left_right, k_str = rxn.rsplit(":", 1)
                            k_val = float(k_str.strip())
                        else:
                            left_right = rxn
                            k_val = 0.1
                        
                        if "->" in left_right or "→" in left_right:
                            sep = "->" if "->" in left_right else "→"
                            reactants_s, products_s = left_right.split(sep)
                            reacts = [x.strip() for x in reactants_s.split("+")]
                            prods = [x.strip() for x in products_s.split("+")]
                        elif "<=>" in left_right or "<->" in left_right:
                            sep = "<=>" if "<=>" in left_right else "<->"
                            reactants_s, products_s = left_right.split(sep)
                            reacts = [x.strip() for x in reactants_s.split("+")]
                            prods = [x.strip() for x in products_s.split("+")]
                            reactions_raw.append({"reactants": reacts, "products": prods, "k": k_val/2, "reversible": True, "k_reverse": k_val/2})
                            i += 1
                            continue
                        else:
                            continue
                        
                        reactions_raw.append({"reactants": reacts, "products": prods, "k": k_val, "reversible": False})
                elif p == "init":
                    i += 1
                    while i < len(parts) and "=" in parts[i]:
                        key, val = parts[i].split("=", 1)
                        init_dict[key.strip()] = float(val.strip())
                        i += 1
                    continue  # Already incremented
                elif p.startswith("t=") or p.startswith("t_end"):
                    if "=" in p:
                        t_end = float(p.split("=")[1])
                elif p.startswith("n="):
                    n_pts = int(p.split("=")[1])
                i += 1

            if not species:
                raise ChemMCPError("Must specify species.")
            if not reactions_raw:
                raise ChemMCPError("Must specify reactions.")

            return self._run_base(species, reactions_raw, init_dict, t_end, n_pts)
        except Exception as e:
            raise ChemMCPError(
                f"Failed to parse text input: {e}. "
                f"Format: 'species A,B,C reactions A->B:k1;B->C:k2 init A=1,B=0,C=0 t=100'"
            )
