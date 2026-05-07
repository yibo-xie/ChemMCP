import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R = 8.314462618     # J/(mol·K)
_NA = 6.02214076e23   # mol^-1

@ChemMCPManager.register_tool
class ChemicalPotentialAdvanced(BaseTool):
    """
    高级化学势计算工具 — 多组分平衡、活度系数、偏摩尔量。
    
    支持理想/非理想溶液、气体混合物、逸度-化学势关系。
    """
    __version__ = "0.1.0"
    name = "ChemicalPotentialAdvanced"
    func_name = "calculate_chemical_potential"
    description = "Calculate chemical potential (μ) for pure substances and mixtures, activity-based calculations for non-ideal systems, partial molar quantities, and multi-component equilibrium conditions."
    implementation_description = "μ_i = μ_i° + RT·ln(a_i) where a_i = γ_i·x_i for solutions or a_i = f_i/P° for gases. Supports ideal gas mixtures (μ_i = μ_i° + RT·ln(P_i/P°)), ideal solutions (μ_i = μ_i° + RT·ln(x_i)), and regular solution model with activity coefficients."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chemical Potential", "Activity", "Fugacity", "Multi-component Equilibrium", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("calculation_type", "str", "N/A", "'pure_substance', 'ideal_gas_mixture', 'ideal_solution', 'regular_solution', 'gibbs_duhem', or 'equilibrium_condition'"),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        # For pure substance
        ("mu_standard", "float", "N/A", "Standard chemical potential μ° in J/mol."),
        ("pressure_atm", "float", "1.0", "Pressure in atm (for gas)."),
        # For mixture / solution
        ("mole_fractions", "list", "N/A", "List of mole fractions [x1, x2, ...]."),
        ("mu_standard_list", "list", "N/A", "List of standard potentials [μ1°, μ2°, ...] in J/mol."),
        ("activity_coefficients", "list", "[]", "Activity coefficients [γ1, γ2, ...]. Default all 1 (ideal)."),
        # For regular solution
        ("omega_J_mol", "float", "0.0", "Regular solution interaction parameter Ω in J/mol (W parameter)."),
        # For Gibbs-Duhem check
        ("total_G_J_mol", "float", "N/A", "Total Gibbs energy per mole."),
        # For equilibrium condition
        ("reaction_mu_products", "list", "N/A", "Chemical potentials of products at reaction conditions."),
        ("reaction_mu_reactants", "list", "N/A", "Chemical potentials of reactants at reaction conditions."),
        ("stoich_coeffs_prod", "list", "N/A", "Stoichiometric coefficients of products."),
        ("stoich_coeffs_react", "list", "N/A", "Stoichiometric coefficients of reactants."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: calc_type|T|[params...]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with chemical potential(s), activities, mixing properties, and equilibrium assessment."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_type": "ideal_gas_mixture",
                "temperature_k": 300.0,
                "mole_fractions": [0.3, 0.7],
                "pressure_atm": 2.0,
            },
            "text_input": {
                "input_str": "ideal_gas_mixture|300||[0.3,0.7]|2"
            },
            "output": {
                "result": {
                    "chemical_potentials_J_mol": ["<value>", "<value>"],
                    "partial_pressures_atm": [0.6, 1.4],
                    "total_pressure_atm": 2.0,
                }
            },
        },
        {
            "code_input": {
                "calculation_type": "regular_solution",
                "temperature_k": 298.15,
                "mole_fractions": [0.5, 0.5],
                "omega_J_mol": 8000.0,
            },
            "text_input": {
                "input_str": "regular_solution|298.15|[0.5,0.5]||8000"
            },
            "output": {
                "result": {
                    "activity_coefficients": ["<value>", "<value>"],
                    "chemical_potentials_J_mol": ["<value>", "<value>"],
                    "excess_Gibbs_J_mol": "<value>",
                    "miscibility": "partially miscible" or "fully miscible",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _pure_substance(self, T: float, mu0: float, P: float, P0: float = 1.0) -> dict:
        """μ(T,P) = μ°(T) + RT·ln(P/P°), assuming ideal gas."""
        if P <= 0:
            raise ChemMCPError("Pressure must be positive.")
        if P0 <= 0:
            raise ChemMCPError("Reference pressure must be positive.")
        
        mu = mu0 + _R * T * math.log(P / P0)
        
        return {
            "calculation_type": "pure_substance",
            "temperature_K": T,
            "standard_potential_J_mol": mu0,
            "pressure_atm": P,
            "reference_pressure_atm": P0,
            "chemical_potential_J_mol": round(mu, 4),
            "RT_ln_term_J_mol": round(_R * T * math.log(P / P0), 4),
        }

    def _ideal_gas_mixture(self, T: float, x_list: List[float], P_total: float, mu0_list: List[float] = None) -> dict:
        """μ_i = μ_i° + RT·ln(P_i/P°) where P_i = x_i · P_total."""
        n = len(x_list)
        total = sum(x_list)
        xs = [xi / total for xi in x_list]
        
        mus = []
        partial_ps = []
        for i, x in enumerate(xs):
            Pi = x * P_total
            partial_ps.append(round(Pi, 6))
            
            if mu0_list is not None and i < len(mu0_list):
                mu_i = mu0_list[i] + _R * T * math.log(max(Pi, 1e-300))
            else:
                # Report relative to pure component at P_total
                mu_i = _R * T * math.log(max(Pi / 1.0, 1e-300))  # relative to 1 atm standard state
            mus.append(round(mu_i, 4))
        
        # Mixing Gibbs energy: ΔG_mix = RT·Σ x_i·ln(x_i)
        dG_mix = _R * T * sum(xi * math.log(max(xi, 1e-300)) for xi in xs)
        
        return {
            "calculation_type": "ideal_gas_mixture",
            "temperature_K": T,
            "total_pressure_atm": P_total,
            "mole_fractions": [round(x, 6) for x in xs],
            "partial_pressures_atm": partial_ps,
            "chemical_potentials_relative_J_mol": mus,
            "delta_G_mixing_J_mol": round(dG_mix, 4),
            "n_components": n,
        }

    def _ideal_solution(self, T: float, x_list: List[float], mu0_list: List[float]) -> dict:
        """μ_i = μ_i° + RT·ln(x_i). Raoult's law behavior."""
        n = len(x_list)
        total = sum(x_list)
        xs = [xi / total for xi in x_list]
        
        if len(mu0_list) != n:
            raise ChemMCPError(f"Need {n} standard potentials, got {len(mu0_list)}.")
        
        mus = []
        activities = []
        for i, x in enumerate(xs):
            a_i = x  # ideal: γ=1
            activities.append(round(a_i, 6))
            mu_i = mu0_list[i] + _R * T * math.log(max(a_i, 1e-300))
            mus.append(round(mu_i, 4))
        
        dG_mix = _R * T * sum(xi * math.log(max(xi, 1e-300)) for xi in xs)
        dS_mix = -_R * sum(xi * math.log(max(xi, 1e-300)) for xi in xs)
        dH_mix = 0.0  # ideal solution
        
        return {
            "calculation_type": "ideal_solution",
            "temperature_K": T,
            "mole_fractions": [round(x, 6) for x in xs],
            "activities": activities,
            "activity_coefficients": [1.0] * n,
            "chemical_potentials_J_mol": mus,
            "delta_G_mixing_J_mol": round(dG_mix, 4),
            "delta_S_mixing_J_mol_K": round(dS_mix, 4),
            "delta_H_mixing_J_mol": 0.0,
            "n_components": n,
        }

    def _regular_solution(self, T: float, x_list: List[float], omega: float, mu0_list: List[float] = None) -> dict:
        """Regular solution model: ln(γ_1) = (Ω/RT)(1-x_1)², ln(γ_2) = (Ω/RT)(1-x_2)²."""
        n = len(x_list)
        total = sum(x_list)
        xs = [xi / total for xi in x_list]
        
        gammas = []
        activities = []
        mus = []
        
        for i, x in enumerate(xs):
            ln_gamma = (omega / (_R * T)) * ((1.0 - x) ** 2)
            gamma = math.exp(ln_gamma)
            gammas.append(round(gamma, 6))
            a_i = gamma * x
            activities.append(round(a_i, 6))
            
            if mu0_list is not None and i < len(mu0_list):
                mu_i = mu0_list[i] + _R * T * math.log(max(a_i, 1e-300))
            else:
                mu_i = _R * T * math.log(max(a_i, 1e-300))
            mus.append(round(mu_i, 4))
        
        # Excess Gibbs energy: G^E = Ω·x1·x2 (per mole for binary)
        G_excess = omega * xs[0] * xs[1] if n == 2 else omega * sum(xs[i]*xs[j] for i in range(n) for j in range(i+1,n))
        
        # Ideal mixing part
        dG_mix_ideal = _R * T * sum(xi * math.log(max(xi, 1e-300)) for xi in xs)
        dG_mix_total = dG_mix_ideal + G_excess
        
        # Check critical miscibility temperature for binary
        T_critical = omega / (2 * _R) if n == 2 else None
        miscibility = "fully miscible" if T > (T_critical or 0) else "phase separation likely"
        
        return {
            "calculation_type": "regular_solution",
            "temperature_K": T,
            "omega_J_mol": omega,
            "mole_fractions": [round(x, 6) for x in xs],
            "activity_coefficients": gammas,
            "activities": activities,
            "chemical_potentials_J_mol": mus,
            "excess_Gibbs_J_mol": round(G_excess, 4),
            "ideal_delta_G_mixing_J_mol": round(dG_mix_ideal, 4),
            "total_delta_G_mixing_J_mol": round(dG_mix_total, 4),
            "critical_temperature_K": round(T_critical, 2) if T_critical else None,
            "miscibility_assessment": miscibility,
            "n_components": n,
        }

    def _gibbs_duhem_check(self, T: float, x_list: List[float], mu_list: List[float], total_G: float) -> dict:
        """Gibbs-Duhem check: Σ x_i·dμ_i = 0 at constant T,P. Verify consistency."""
        n = len(x_list)
        total_x = sum(x_list)
        xs = [xi / total_x for xi in x_list]
        
        # Check: G = Σ x_i·μ_i
        G_from_mu = sum(xi * mui for xi, mui in zip(xs, mu_list))
        diff = abs(G_from_mu - total_G)
        consistent = diff < max(abs(total_G), 1.0) * 0.01  # within 1%
        
        return {
            "calculation_type": "gibbs_duhem_check",
            "temperature_K": T,
            "given_total_G_J_mol": total_G,
            "computed_G_from_mu_J_mol": round(G_from_mu, 4),
            "difference_J_mol": round(diff, 4),
            "difference_percent": round(diff / max(abs(total_G), 1e-10) * 100, 4),
            "consistent": consistent,
            "note": "Gibbs-Duhem requires Σx_i·μ_i = G at constant T,P.",
        }

    def _equilibrium_condition(self, mu_prods: List[float], mu_reacts: List[float],
                                coeffs_prod: List[float], coeffs_react: List[float]) -> dict:
        """At equilibrium: Σ ν_prod·μ_prod = Σ ν_react·μ_react."""
        mu_prod_side = sum(c * m for c, m in zip(coeffs_prod, mu_prods))
        mu_react_side = sum(c * m for c, m in zip(coeffs_react, mu_reacts))
        delta_mu = mu_prod_side - mu_react_side
        
        eps = 1e-6
        if abs(delta_mu) < eps:
            status = "at_equilibrium"
        elif delta_mu < 0:
            status = "spontaneous_forward"
        else:
            status = "spontaneous_reverse"
        
        return {
            "calculation_type": "equilibrium_condition",
            "sum_nu_mu_products_J_mol": round(mu_prod_side, 4),
            "sum_nu_mu_reactants_J_mol": round(mu_react_side, 4),
            "delta_mu_reaction_J_mol": round(delta_mu, 4),
            "equilibrium_status": status,
            "driving_direction": "forward" if delta_mu < eps else ("reverse" if delta_mu > eps else "none"),
        }

    def _run_base(self, calculation_type: str, temperature_k: float = 298.15,
                  mu_standard: float = None, pressure_atm: float = 1.0,
                  mole_fractions: List[float] = None, mu_standard_list: List[float] = None,
                  activity_coefficients: List[float] = None, omega_J_mol: float = 0.0,
                  total_G_J_mol: float = None, reaction_mu_products: List[float] = None,
                  reaction_mu_reactants: List[float] = None,
                  stoich_coeffs_prod: List[float] = None, stoich_coeffs_react: List[float] = None) -> dict:
        ct = calculation_type.lower().strip()
        
        if ct == "pure_substance":
            if mu_standard is None:
                raise ChemMCPError("Need mu_standard for pure substance.")
            return self._pure_substance(temperature_k, mu_standard, pressure_atm)
        
        elif ct == "ideal_gas_mixture":
            if mole_fractions is None:
                raise ChemMCPError("Need mole_fractions.")
            return self._ideal_gas_mixture(temperature_k, mole_fractions, pressure_atm, mu_standard_list)
        
        elif ct == "ideal_solution":
            if mole_fractions is None or mu_standard_list is None:
                raise ChemMCPError("Need mole_fractions and mu_standard_list for ideal solution.")
            return self._ideal_solution(temperature_k, mole_fractions, mu_standard_list)
        
        elif ct == "regular_solution":
            if mole_fractions is None:
                raise ChemMCPError("Need mole_fractions for regular solution.")
            return self._regular_solution(temperature_k, mole_fractions, omega_J_mol, mu_standard_list)
        
        elif ct == "gibbs_duhem":
            if any(x is None for x in [mole_fractions, mu_standard_list, total_G_J_mol]):
                raise ChemMCPError("Need mole_fractions, mu_standard_list (as actual μ values), and total_G_J_mol.")
            return self._gibbs_duhem_check(temperature_k, mole_fractions, mu_standard_list, total_G_J_mol)
        
        elif ct == "equilibrium_condition":
            if any(x is None for x in [reaction_mu_products, reaction_mu_reactants, stoich_coeffs_prod, stoich_coeffs_react]):
                raise ChemMCPError("Need reaction_mu_products, reaction_mu_reactants, stoich_coeffs_prod, stoich_coeffs_react.")
            return self._equilibrium_condition(reaction_mu_products, reaction_mu_reactants, stoich_coeffs_prod, stoich_coeffs_react)
        
        else:
            raise ChemMCPError(
                f"Unknown type: '{calculation_type}'. "
                f"Options: 'pure_substance', 'ideal_gas_mixture', 'ideal_solution', "
                f"'regular_solution', 'gibbs_duhem', 'equilibrium_condition'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            ct = parts[0].strip()
            T = float(parts[1]) if len(parts) > 1 and parts[1].strip() else 298.15
            import json
            def _parse_float(s):
                s = (s or '').strip().strip("'\"")
                if not s:
                    return None
                return float(s)
            def _parse_list(s):
                s = (s or '').strip().strip("'\"")
                if not s:
                    return None
                v = json.loads(s)
                return v if isinstance(v, list) else None
            # Position-aware parsing based on _run_base signature
            mu0 = _parse_float(parts[2]) if len(parts) > 2 else None
            P = _parse_float(parts[3]) if len(parts) > 3 else None
            P = P if P is not None else 1.0
            xs = _parse_list(parts[4]) if len(parts) > 4 else None
            mu0s = _parse_list(parts[5]) if len(parts) > 5 else None
            gammas = _parse_list(parts[6]) if len(parts) > 6 else []
            omega = _parse_float(parts[7]) if len(parts) > 7 else 0.0
            Gtot = _parse_float(parts[8]) if len(parts) > 8 else None
            mu_p = _parse_list(parts[9]) if len(parts) > 9 else None
            mu_r = _parse_list(parts[10]) if len(parts) > 10 else None
            cp = _parse_list(parts[11]) if len(parts) > 11 else None
            cr = _parse_list(parts[12]) if len(parts) > 12 else None
            return self._run_base(ct, T, mu0, P, xs, mu0s, gammas, omega, Gtot, mu_p, mu_r, cp, cr)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
