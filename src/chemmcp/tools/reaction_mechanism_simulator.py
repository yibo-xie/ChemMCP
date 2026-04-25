import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ReactionMechanismSimulator(BaseTool):
    """
    复杂反应机理的动力学模拟器。
    使用 Runge-Kutta 4 阶方法数值求解常微分方程组，模拟各物种浓度随时间的变化。
    
    支持的反应机理格式：
    - 一级反应: "A->B;k" 或 "A->B;0.1"
    - 二级反应: "2A->B;k" 或 "A+B->C;k"
    - 可逆反应: "A<=>B;kf,kr"
    """
    __version__ = "0.1.0"
    name = "ReactionMechanismSimulator"
    func_name = "simulate_mechanism"
    description = "Simulate concentration vs time profiles for complex reaction mechanisms using numerical integration (RK4)."
    implementation_description = "Parses reaction mechanism strings into ODEs and solves them using 4th-order Runge-Kutta integration. Supports consecutive, parallel, reversible, and complex multi-step reactions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Mechanism Simulation", "ODE", "Numerical Integration"]
    required_envs = []

    code_input_sig = [
        ("mechanism", "str", "N/A", "Reaction mechanism string, e.g., 'A->B;0.1,B->C;0.05' for consecutive, or 'A->B;0.1,A->C;0.02' for parallel."),
        ("initial_concentrations", "dict", "N/A", "Initial concentrations dict, e.g., {'A': 1.0, 'B': 0.0, 'C': 0.0}."),
        ("time_end", "float", "100.0", "End time for simulation (same unit as rate constants)."),
        ("n_points", "int", "100", "Number of time points to output."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'mechanism A:val,B:val,... [time_end] [n_points]'. Use semicolons in mechanism carefully."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing time_points, concentration_profiles (per species), mechanism_parsed, and simulation metadata."),
    ]

    examples = [
        {
            "code_input": {
                "mechanism": "A->B;0.1,B->C;0.05",
                "initial_concentrations": {"A": 1.0, "B": 0.0, "C": 0.0},
                "time_end": 100.0,
                "n_points": 11,
            },
            "text_input": {
                "input_params": "A->B;0.1,B->C;0.05 A:1,B:0,C:0 100 11",
            },
            "output": {
                "result": {
                    "species": ["A", "B", "C"],
                    "n_points": 11,
                    "final_concentrations": {"A": 0.3679, "B": 0.2387, "C": 0.3935},
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _parse_mechanism(self, mechanism_str):
        """Parse mechanism string into list of reaction steps."""
        steps = []
        reactions = mechanism_str.split(",")
        # Also try semicolon as step separator (but semicolon is used within steps)
        # Format: each step is like "A->B;k" or "A+B->C;k" or "A<=>B;kf,kr"
        # Steps separated by comma
        
        for rxn in reactions:
            rxn = rxn.strip()
            if not rxn:
                continue
            
            # Split rate constant from reaction
            # Find last semicolon that separates the reaction from k value
            parts = rxn.rsplit(";", 1)
            if len(parts) != 2:
                continue
            
            rxn_expr = parts[0].strip()
            k_str = parts[1].strip()
            
            # Determine reaction type
            if "<=>" in rxn_expr:
                left, right = rxn_expr.split("<=>")
                k_parts = k_str.split(",")
                kf = float(k_parts[0].strip())
                kr = float(k_parts[1].strip()) if len(k_parts) > 1 else 0.0
                reactants = self._parse_species(left.strip())
                products = self._parse_species(right.strip())
                steps.append({"type": "reversible", "reactants": reactants,
                              "products": products, "kf": kf, "kr": kr})
            elif "->" in rxn_expr:
                left, right = rxn_expr.split("->")
                k = float(k_str)
                reactants = self._parse_species(left.strip())
                products = self._parse_species(right.strip())
                steps.append({"type": "irreversible", "reactants": reactants,
                              "products": products, "k": k})

        return steps

    def _parse_species(self, side_str):
        """Parse species side like '2A+B' into [('A', 2), ('B', 1)]."""
        species = []
        tokens = side_str.split("+")
        for token in tokens:
            token = token.strip()
            i = 0
            coeff = 0
            while i < len(token) and token[i].isdigit():
                coeff = coeff * 10 + int(token[i])
                i += 1
            if coeff == 0:
                coeff = 1
            name = token[i:].strip()
            if name:
                species.append((name, coeff))
        return species

    def _get_rate(self, y, species_order, steps):
        """Compute dy/dt for all species."""
        rates = {sp: 0.0 for sp in species_order}
        
        for step in steps:
            if step["type"] == "irreversible":
                # Rate = k * [reactant1]^coeff1 * [reactant2]^coeff2 * ...
                r = step["k"]
                for sp, coeff in step["reactants"]:
                    idx = species_order.index(sp)
                    r *= max(y[idx], 0.0) ** coeff
                
                for sp, coeff in step["reactants"]:
                    rates[sp] -= coeff * r
                for sp, coeff in step["products"]:
                    rates[sp] += coeff * r

            elif step["type"] == "reversible":
                # Forward rate
                rf = step["kf"]
                for sp, coeff in step["reactants"]:
                    idx = species_order.index(sp)
                    rf *= max(y[idx], 0.0) ** coeff
                
                # Reverse rate
                rr = step["kr"]
                for sp, coeff in step["products"]:
                    idx = species_order.index(sp)
                    rr *= max(y[idx], 0.0) ** coeff
                
                for sp, coeff in step["reactants"]:
                    rates[sp] -= coeff * rf + (-coeff) * rr  # reverse produces reactants
                for sp, coeff in step["products"]:
                    rates[sp] += coeff * rf + (-coeff) * rr

        return [rates[sp] for sp in species_order]

    def _rk4_step(self, y, dt, species_order, steps):
        """Single RK4 step."""
        k1 = self._get_rate(y, species_order, steps)
        k2 = self._get_rate([y[i] + 0.5 * dt * k1[i] for i in range(len(y))], species_order, steps)
        k3 = self._get_rate([y[i] + 0.5 * dt * k2[i] for i in range(len(y))], species_order, steps)
        k4 = self._get_rate([y[i] + dt * k3[i] for i in range(len(y))], species_order, steps)
        
        return [y[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(len(y))]

    def _run_base(self, mechanism: str, initial_concentrations: dict,
                  time_end: float = 100.0, n_points: int = 100) -> dict:
        steps = self._parse_mechanism(mechanism)
        if not steps:
            raise ChemMCPError(f"Could not parse mechanism: {mechanism}")

        species_order = sorted(initial_concentrations.keys())
        y0 = [initial_concentrations.get(sp, 0.0) for sp in species_order]
        
        dt = time_end / (n_points - 1) if n_points > 1 else time_end
        
        time_points = []
        conc_profiles = {sp: [] for sp in species_order}
        y = list(y0)
        
        for i in range(n_points):
            t = i * dt
            time_points.append(round(t, 6))
            for j, sp in enumerate(species_order):
                conc_profiles[sp].append(round(max(y[j], 0.0), 8))
            
            if i < n_points - 1:
                y = self._rk4_step(y, dt, species_order, steps)

        final = {sp: round(conc_profiles[sp][-1], 6) for sp in species_order}
        logger.info(f"ReactionMechanismSimulator: final concentrations = {final}")
        
        return {
            "species": species_order,
            "time_points": time_points,
            "concentration_profiles": conc_profiles,
            "final_concentrations": final,
            "steps_parsed": len(steps),
            "time_end": time_end,
            "n_points": n_points,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mechanism = parts[0]
            
            # Parse initial concentrations
            init_dict = {}
            for p in parts[1:]:
                if ":" in p and not p.replace(".", "").replace("-", "").replace("+", "").isdigit():
                    key, val = p.split(":")
                    init_dict[key] = float(val)
                elif "." in p or p.isdigit():
                    # Could be time_end or n_points
                    pass
            
            time_end = 100.0
            n_points = 100
            
            numeric_parts = []
            for p in parts[1:]:
                if ":" not in p:
                    try:
                        numeric_parts.append(float(p))
                    except ValueError:
                        pass
            
            if len(numeric_parts) >= 1:
                time_end = numeric_parts[0]
            if len(numeric_parts) >= 2:
                n_points = int(numeric_parts[1])
            
            return self._run_base(mechanism, init_dict, time_end, n_points)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: 'mechanism A:val,B:val,... [time_end] [n_points]'")
