import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R_J = 8.314462618     # J/(mol·K)
_R_Latm = 0.082057366   # L·atm/(mol·K)
_NA = 6.02214076e23     # mol^-1

@ChemMCPManager.register_tool
class EquationOfState(BaseTool):
    """
    状态方程求解工具 - 理想气体、范德华方程、维里方程。

    求解P-V-T关系,计算压缩因子Z和偏离理想气体的程度。
    """
    __version__ = "0.1.0"
    name = "EquationOfState"
    func_name = "solve_equation_of_state"
    description = "Solve equations of state: ideal gas law, van der Waals equation, and virial equation. Calculate compressibility factor Z and deviations from ideal gas behavior."
    implementation_description = "Ideal: PV=nRT. Van der Waals: (P + a(n/V)2)(V - nb) = nRT. Virial: Z = 1 + B/Vm + C/Vm2 + ... Uses Newton-Raphson for solving V from P and T for non-ideal equations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Equation of State", "Van der Waals", "Virial", "Compressibility", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("equation_type", "str", "N/A", "'ideal_gas', 'van_der_waals', or 'virial'"),
        ("solve_for", "str", "V", "Variable to solve for: 'P', 'V', 'T', or 'Z'"),
        ("pressure_atm", "float", "N/A", "Pressure in atm."),
        ("volume_L", "float", "N/A", "Volume in L."),
        ("temperature_k", "float", "N/A", "Temperature in Kelvin."),
        ("n_moles", "float", "1.0", "Amount of substance in moles."),
        # van der Waals parameters
        ("a_param", "float", "0.0", "van der Waals a parameter in atm·L2/mol2."),
        ("b_param", "float", "0.0", "van der Waals b parameter in L/mol."),
        # Virial coefficients
        ("virial_B", "float", "0.0", "Second virial coefficient B in L/mol."),
        ("virial_C", "float", "0.0", "Third virial coefficient C in L2/mol2."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: eq_type|solve_for|P|V|T|n|a|b|B_virial|C_virial"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with solved value(s), compressibility factor Z, deviation from ideality, and intermediate values."),
    ]

    examples = [
        {
            "code_input": {
                "equation_type": "ideal_gas",
                "solve_for": "V",
                "pressure_atm": 1.0,
                "temperature_k": 298.15,
                "n_moles": 1.0,
            },
            "text_input": {
                "input_str": "ideal_gas|V|1||298.15|1"
            },
            "output": {
                "result": {
                    "volume_L": 24.45,
                    "compressibility_factor_Z": 1.0,
                    "deviation_percent": 0.0,
                }
            },
        },
        {
            "code_input": {
                "equation_type": "van_der_waals",
                "solve_for": "P",
                "volume_L": 10.0,
                "temperature_k": 300.0,
                "n_moles": 1.0,
                "a_param": 1.39,
                "b_param": 0.0391,
            },
            "text_input": {
                "input_str": "van_der_waals|P|10|300|1|1|1.39|0.0391"
            },
            "output": {
                "result": {
                    "pressure_atm": "<value>",
                    "compressibility_factor_Z": "<value>",
                }
            },
        },
        {
            "code_input": {
                "equation_type": "virial",
                "solve_for": "Z",
                "pressure_atm": 10.0,
                "temperature_k": 300.0,
                "virial_B": -0.02,
            },
            "text_input": {
                "input_str": "virial|Z|10||300|1|||-0.02"
            },
            "output": {
                "result": {
                    "compressibility_factor_Z": "<value>",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _ideal_gas(self, solve_for: str, P: float, V: float, T: float, n: float) -> dict:
        """PV = nRT."""
        R = _R_Latm

        if solve_for.upper() == "P":
            if V <= 0:
                raise ChemMCPError("Volume must be positive.")
            val = n * R * T / V
            Z = 1.0
            return {"solved_variable": "P", "pressure_atm": round(val, 4), "volume_L": V, "temperature_K": T, "compressibility_factor_Z": 1.0, "deviation_from_ideal_pct": 0.0}

        elif solve_for.upper() == "V":
            if P <= 0:
                raise ChemMCPError("Pressure must be positive.")
            val = n * R * T / P
            return {"solved_variable": "V", "volume_L": round(val, 4), "pressure_atm": P, "temperature_K": T, "compressibility_factor_Z": 1.0, "deviation_from_ideal_pct": 0.0}

        elif solve_for.upper() == "T":
            if P <= 0 or V <= 0:
                raise ChemMCPError("P and V must be positive.")
            val = P * V / (n * R)
            return {"solved_variable": "T", "temperature_K": round(val, 4), "pressure_atm": P, "volume_L": V, "compressibility_factor_Z": 1.0, "deviation_from_ideal_pct": 0.0}

        elif solve_for.upper() == "Z":
            return {"solved_variable": "Z", "compressibility_factor_Z": 1.0, "note": "Ideal gas always has Z=1."}

        else:
            raise ChemMCPError(f"Cannot solve for '{solve_for}'. Options: P, V, T, Z.")

    def _vdw(self, solve_for: str, P: float, V: float, T: float, n: float, a: float, b: float) -> dict:
        """(P + a(n/V)2)(V - nb) = nRT."""
        R = _R_Latm

        if solve_for.upper() == "P":
            if V <= n * b:
                raise ChemMCPError(f"Volume {V} L must exceed nb = {n*b} L for vdW equation.")
            P_calc = n * R * T / (V - n * b) - a * (n ** 2) / (V ** 2)
            Vm = V / n
            Z = P_calc * Vm / (R * T)
            P_ideal = n * R * T / V
            deviation = ((P_calc - P_ideal) / P_ideal * 100) if P_ideal > 0 else 0

            return {
                "solved_variable": "P", "pressure_atm": round(P_calc, 6),
                "volume_L": V, "temperature_K": T, "n_moles": n,
                "compressibility_factor_Z": round(Z, 6),
                "ideal_pressure_atm": round(P_ideal, 4),
                "deviation_from_ideal_pct": round(deviation, 4),
                "a_atm_L2_mol2": a, "b_L_mol": b,
            }

        elif solve_for.upper() == "V":
            # Newton-Raphson to solve for V
            if P <= 0:
                raise ChemMCPError("Pressure must be positive.")

            # Initial guess: ideal gas volume
            V_guess = n * R * T / P
            V_guess = max(V_guess, n * b * 1.5)

            for iteration in range(100):
                term1 = P + a * n**2 / V_guess**2
                term2 = V_guess - n * b
                f = term1 * term2 - n * R * T

                # f'(V) = P - a*n2/V3 * (V-nb) + (P + a*n2/V2)
                f_prime = term1 + term2 * (-2 * a * n**2 / V_guess**3) + P

                if abs(f_prime) < 1e-20:
                    break
                delta = f / f_prime
                V_new = V_guess - delta

                if V_new <= n * b:
                    V_new = n * b * 1.01

                if abs(delta) < 1e-8 * abs(V_new):
                    V_guess = V_new
                    break
                V_guess = V_new

            Vm = V_guess / n
            Z = P * Vm / (R * T)
            V_ideal = n * R * T / P
            deviation = ((V_guess - V_ideal) / V_ideal * 100) if V_ideal > 0 else 0

            return {
                "solved_variable": "V", "volume_L": round(V_guess, 6),
                "pressure_atm": P, "temperature_K": T, "n_moles": n,
                "compressibility_factor_Z": round(Z, 6),
                "ideal_volume_L": round(V_ideal, 4),
                "deviation_from_ideal_pct": round(deviation, 4),
                "a_atm_L2_mol2": a, "b_L_mol": b,
            }

        elif solve_for.upper() == "Z":
            # Need V to compute Z
            if V <= 0:
                raise ChemMCPError("Need volume to compute Z for vdW gas.")
            result = self._vdw("P", P, V, T, n, a, b)
            return {**result, "solved_variable": "Z"}

        else:
            raise ChemMCPError(f"Cannot solve for '{solve_for}' with vdW. Options: P, V, Z.")

    def _virial(self, solve_for: str, P: float, V: float, T: float, n: float,
                B: float, C: float) -> dict:
        """Virial equation: Z = 1 + B/Vm + C/Vm2, where Z = PVm/RT."""
        R = _R_Latm

        if solve_for.upper() == "Z" or solve_for.upper() == "P":
            if V is not None and V > 0:
                # Volume-based virial: Z = 1 + B/Vm + C/Vm²
                Vm = V / n
                Z = 1.0 + B / Vm + C / (Vm ** 2)
                P_calc = Z * R * T / Vm
                deviation = (Z - 1.0) * 100

                return {
                    "solved_variable": "Z" if solve_for.upper() == "Z" else "P",
                    "compressibility_factor_Z": round(Z, 6),
                    "pressure_atm": round(P_calc, 6) if solve_for.upper() == "P" else P,
                    "volume_L": V, "molar_volume_L_mol": round(Vm, 4),
                    "temperature_K": T, "n_moles": n,
                    "deviation_from_ideal_pct": round(deviation, 4),
                    "B_L_mol": B, "C_L2_mol2": C,
                }
            elif P is not None and P > 0:
                # Pressure-based virial: Z ≈ 1 + B'·P, where B' = B/(RT)
                B_prime = B / (R * T)
                Z = 1.0 + B_prime * P + (C / ((R*T)**2)) * P**2
                Vm_ideal = R * T / P
                Vm = Z * Vm_ideal
                deviation = (Z - 1.0) * 100
                return {
                    "solved_variable": "Z",
                    "compressibility_factor_Z": round(Z, 6),
                    "pressure_atm": P,
                    "estimated_molar_volume_L_mol": round(Vm, 4),
                    "ideal_molar_volume_L_mol": round(Vm_ideal, 4),
                    "temperature_K": T,
                    "deviation_from_ideal_pct": round(deviation, 4),
                    "B_L_mol": B, "C_L2_mol2": C,
                    "note": "Approximate Z from pressure expansion (valid at low-moderate P).",
                }

        elif solve_for.upper() == "V":
            if P <= 0:
                raise ChemMCPError("Pressure must be positive.")
            # Solve: P*Vm/(RT) = 1 + B/Vm + C/Vm2
            # Rearranged: P*Vm3/(RT) = Vm2 + B*Vm + C
            # Cubic in Vm: (P/RT)*Vm3 - Vm2 - B*Vm - C = 0
            coef_a = P / (R * T)
            coef_b = -1.0
            coef_c = -B
            coef_d = -C

            roots = self._solve_cubic(coef_a, coef_b, coef_c, coef_d)
            # Pick the real positive root closest to ideal gas Vm
            Vm_ideal = R * T / P
            best_Vm = None
            best_diff = float('inf')
            for root in roots:
                if isinstance(root, (int, float)) and root > 0 and root > abs(B) * 0.1:
                    diff = abs(root - Vm_ideal)
                    if diff < best_diff:
                        best_diff = diff
                        best_Vm = root

            if best_Vm is None:
                raise ChemMCPError("No physically meaningful solution found for virial equation.")

            V = best_Vm * n
            Z = 1.0 + B / best_Vm + C / (best_Vm ** 2)
            deviation = (Z - 1.0) * 100

            return {
                "solved_variable": "V", "volume_L": round(V, 6),
                "molar_volume_L_mol": round(best_Vm, 6),
                "pressure_atm": P, "temperature_K": T, "n_moles": n,
                "compressibility_factor_Z": round(Z, 6),
                "ideal_molar_volume_L_mol": round(Vm_ideal, 4),
                "deviation_from_ideal_pct": round(deviation, 4),
                "B_L_mol": B, "C_L2_mol2": C,
            }

        else:
            raise ChemMCPError(f"Cannot solve for '{solve_for}' with virial. Options: P, V, Z.")

    def _solve_cubic(self, a: float, b: float, c: float, d: float) -> list:
        """Solve ax3 + bx2 + cx + d = 0, return all roots (real and complex)."""
        if abs(a) < 1e-30:
            # Quadratic
            disc = b**2 - 4*c*d
            if disc >= 0:
                sq = math.sqrt(disc)
                return [(-b + sq)/(2*c), (-b - sq)/(2*c)]
            else:
                return [complex(-b/(2*c), math.sqrt(-disc)/(2*c))]

        # Normalize: x3 + (b/a)x2 + (c/a)x + d/a = 0
        p = b / a
        q = c / a
        r = d / a

        # Depressed cubic: t3 + pt + q = 0, where x = t - p/3
        p_dep = q - p**2 / 3.0
        q_dep = r - p*q/3.0 + 2*p**3/27.0

        discriminant = (q_dep/2)**2 + (p_dep/3)**3

        if discriminant >= 0:
            s = math.cbrt(-q_dep/2 + math.sqrt(discriminant))
            t_val = math.cbrt(-q_dep/2 - math.sqrt(discriminant))
            root1 = s + t_val - p/3
            # Check for remaining real roots
            if discriminant > 1e-15:
                return [root1]
            else:
                # Triple or double root case
                if abs(s - t_val) < 1e-12:
                    root2 = -s - p/3
                    return [root1, root2, root2]
                return [root1]
        else:
            # Three real roots
            m_phi = 2 * math.sqrt(-p_dep / 3)
            phi = math.acos(-q_dep / (2 * math.sqrt(-(p_dep/3)**3)))
            root1 = m_phi * math.cos(phi/3) - p/3
            root2 = m_phi * math.cos((phi + 2*math.pi)/3) - p/3
            root3 = m_phi * math.cos((phi + 4*math.pi)/3) - p/3
            return [root1, root2, root3]

    def _run_base(self, equation_type: str, solve_for: str = "V", pressure_atm: float = None,
                  volume_L: float = None, temperature_k: float = None, n_moles: float = 1.0,
                  a_param: float = 0.0, b_param: float = 0.0,
                  virial_B: float = 0.0, virial_C: float = 0.0) -> dict:
        eq = equation_type.lower().strip()
        sf = solve_for.strip().upper()

        if eq == "ideal_gas":
            return self._ideal_gas(sf, pressure_atm, volume_L, temperature_k, n_moles)
        elif eq == "van_der_waals" or eq == "vanderwaals" or eq == "vdw":
            return self._vdw(sf, pressure_atm, volume_L, temperature_k, n_moles, a_param, b_param)
        elif eq == "virial":
            return self._virial(sf, pressure_atm, volume_L, temperature_k, n_moles, virial_B, virial_C)
        else:
            raise ChemMCPError(
                f"Unknown equation type: '{equation_type}'. "
                f"Options: 'ideal_gas', 'van_der_waals', 'virial'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            eq = parts[0].strip()
            sf = parts[1].strip() if len(parts) > 1 else "V"
            P = float(parts[2]) if len(parts) > 2 and parts[2].strip() else None
            V = float(parts[3]) if len(parts) > 3 and parts[3].strip() else None
            T = float(parts[4]) if len(parts) > 4 and parts[4].strip() else None
            n = float(parts[5]) if len(parts) > 5 else 1.0
            a = float(parts[6]) if len(parts) > 6 and parts[6].strip() else 0.0
            b = float(parts[7]) if len(parts) > 7 and parts[7].strip() else 0.0
            Bv = float(parts[8]) if len(parts) > 8 and parts[8].strip() else 0.0
            Cv = float(parts[9]) if len(parts) > 9 and parts[9].strip() else 0.0
            return self._run_base(eq, sf, P, V, T, n, a, b, Bv, Cv)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'eq_type|sf|P|V|T|n|a|b|B|C'")
