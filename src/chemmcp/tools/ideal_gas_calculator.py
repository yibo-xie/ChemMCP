import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IdealGasCalculator(BaseTool):
    """
    理想气体状态方程计算工具 (MCP #297)。
    基于 PV = nRT 进行多种计算：
    - 已知3个变量求第4个（P, V, n, T）
    - 气体密度和摩尔质量计算
    - 组合气体定律（P1V1/T1 = P2V2/T2）
    - 分压计算（Dalton定律）
    R = 0.0821 L·atm/(K·mol) 或 8.314 J/(K·mol)
    """
    __version__ = "0.1.0"
    name = "IdealGasCalculator"
    func_name = "calculate_ideal_gas"
    description = "Ideal gas law calculations: PV=nRT solver, density, molar mass, combined gas law, and Dalton's law of partial pressures."
    implementation_description = (
        "Implements the ideal gas equation PV=nRT with R=0.08206 L·atm/(K·mol). "
        "Supports solving for any variable given the other three, plus density, "
        "molar mass from density, combined gas law, and partial pressure calculations."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Ideal Gas", "Thermodynamics", "Physical Chemistry", "PV=nRT"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "solve", "Calculation mode: 'solve' (PV=nRT), 'density', 'molar_mass', 'combined', 'partial_pressure'."),
        ("P", "float", "None", "Pressure in atm. Provide 3 of {P, V, n, T} for solve mode."),
        ("V", "float", "None", "Volume in L."),
        ("n", "float", "None", "Amount of substance in mol."),
        ("T", "float", "None", "Temperature in K."),
        # Additional params for other modes
        ("mass_g", "float", "None", "Mass in grams (for density/molar_mass modes)."),
        ("molar_mass", "float", "None", "Molar mass in g/mol."),
        ("P1", "float", "None", "Initial pressure for combined gas law (atm)."),
        ("V1", "float", "None", "Initial volume for combined gas law (L)."),
        ("T1", "float", "None", "Initial temperature for combined gas law (K)."),
        ("P2", "float", "None", "Final pressure for combined gas law (atm)."),
        ("V2", "float", "None", "Final volume for combined gas law (L)."),
        ("T2", "float", "None", "Final temperature for combined gas law (K)."),
        ("mole_fractions", "list", "None", "List of mole fractions for partial pressure mode."),
        ("P_total", "float", "None", "Total pressure for partial pressure mode (atm)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Natural language query or parameter string, e.g., 'P=1 V=22.4 n=1 T=?', 'density at STP CO2', 'combined P1=1 V1=10 T1=300 P2=2 V2=5 T2=?'."),
    ]

    output_sig = [
        ("result", "float", "Calculated value(s)."),
        ("unit", "str", "Unit of the result."),
        ("step_by_step", "str", "Step-by-step solution explanation."),
        ("given", "dict", "Input parameters that were provided."),
    ]

    examples = [{'code_input': {'mode': 'solve', 'P': 1.0, 'V': 22.4, 'n': 1.0, 'T': None, 'P1': 'N/A', 'P2': 'N/A', 'P_total': 'N/A', 'T1': 'N/A', 'T2': 'N/A', 'V1': 'N/A', 'V2': 'N/A', 'mass_g': 'N/A', 'molar_mass': 'N/A', 'mole_fractions': 'N/A'}, 'text_input': {'query': 'P=1atm V=22.4L n=1mol T=?'}, 'output': {'result': 273.15, 'given': 'N/A', 'step_by_step': 'N/A', 'unit': 'N/A'}}, {'code_input': {'mode': 'density', 'T': 273.15, 'P': 1.0, 'molar_mass': 44.01, 'P1': 'N/A', 'P2': 'N/A', 'P_total': 'N/A', 'T1': 'N/A', 'T2': 'N/A', 'V': 'N/A', 'V1': 'N/A', 'V2': 'N/A', 'mass_g': 'N/A', 'mole_fractions': 'N/A', 'n': 'N/A'}, 'text_input': {'query': 'density of CO2 at STP'}, 'output': {'result': 1.96, 'given': 'N/A', 'step_by_step': 'N/A', 'unit': 'N/A'}}, {'code_input': {'mode': 'combined', 'P1': 1.0, 'V1': 10.0, 'T1': 300.0, 'P2': 2.0, 'V2': 5.0, 'T2': None, 'P': 'N/A', 'P_total': 'N/A', 'T': 'N/A', 'V': 'N/A', 'mass_g': 'N/A', 'molar_mass': 'N/A', 'mole_fractions': 'N/A', 'n': 'N/A'}, 'text_input': {'query': 'combined: P1=1 V1=10 T1=300 P2=2 V2=5 T2=?'}, 'output': {'result': 300.0, 'given': 'N/A', 'step_by_step': 'N/A', 'unit': 'N/A'}}]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 0.08206  # L·atm/(K·mol)

    def _run_base(self, mode: str = "solve", **kwargs) -> dict:
        mode = mode.lower().strip()
        if mode == "solve":
            return self._solve_pv_nrt(**kwargs)
        elif mode == "density":
            return self._calc_density(**kwargs)
        elif mode == "molar_mass":
            return self._calc_molar_mass(**kwargs)
        elif mode == "combined":
            return self._combined_law(**kwargs)
        elif mode == "partial_pressure":
            return self._partial_pressure(**kwargs)
        else:
            raise ChemMCPInputError(f"Unknown mode '{mode}'. Use: solve, density, molar_mass, combined, partial_pressure.")

    def _solve_pv_nrt(self, P=None, V=None, n=None, T=None, **kwargs) -> dict:
        """Solve PV = nRT for the missing variable."""
        provided = {k: v for k, v in [("P", P), ("V", V), ("n", n), ("T", T)] if v is not None}
        if len(provided) != 3:
            raise ChemMCPInputError(f"Provide exactly 3 of P, V, n, T. Got: {list(provided.keys())}")

        R = self.R
        steps = f"Using ideal gas law: PV = nRT, R = {R} L·atm/(K·mol)\n"

        if P is None:
            result = n * R * T / V
            unit = "atm"
            steps += f"P = nRT/V = ({n})({R})({T})/({V}) = {round(result, 4)} atm"
        elif V is None:
            result = n * R * T / P
            unit = "L"
            steps += f"V = nRT/P = ({n})({R})({T})/({P}) = {round(result, 4)} L"
        elif n is None:
            result = P * V / (R * T)
            unit = "mol"
            steps += f"n = PV/(RT) = ({P})({V})/(({R})({T})) = {round(result, 6)} mol"
        elif T is None:
            result = P * V / (n * R)
            unit = "K"
            steps += f"T = PV/(nR) = ({P})({V})/(({n})({R})) = {round(result, 2)} K"
        else:
            raise ChemMCPInputError("All values provided; nothing to solve.")

        # Verification
        check = P if P is not None else result
        check_V = V if V is not None else result
        check_n = n if n is not None else result
        check_T = T if T is not None else result

        return {
            "result": round(result, 6),
            "unit": "K",
            "step_by_step": "",
            "given": {k: v for k, v in [("P", P), ("V", V), ("n", n), ("T", T)] if v is not None},
            "verification": f"Check: PV/nT = ({check}×{check_V})/({check_n}×{check_T}) should equal R",
        }

    def _calc_density(self, P=None, V=None, n=None, T=None, mass_g=None, molar_mass=None, **kwargs) -> dict:
        """Calculate gas density ρ = PM/(RT) = m/V."""
        if P is not None and T is not None and molar_mass is not None:
            # ρ = PM/(RT)
            rho = P * molar_mass / (self.R * T)
            steps = f"Density ρ = PM/(RT) = ({P})({molar_mass})/({self.R}×{T}) = {round(rho, 4)} g/L"
            return {"result": round(rho, 4)}
        elif mass_g is not None and V is not None:
            # ρ = m/V
            rho = mass_g / V
            steps = f"Density ρ = m/V = ({mass_g})/({V}) = {round(rho, 4)} g/L"
            return {"result": round(rho, 4)}
        elif molar_mass is not None and V is not None and n is not None:
            # M = m/n, ρ = m/V = nM/V
            rho = n * molar_mass / V
            steps = f"Density ρ = nM/V = ({n})({molar_mass})/({V}) = {round(rho, 4)} g/L"
            return {"result": round(rho, 4)}
        else:
            raise ChemMCPInputError(
                "For density mode provide:\n"
                "  (1) P + T + molar_mass → ρ = PM/(RT)\n"
                "  (2) mass_g + V → ρ = m/V\n"
                "  (3) n + molar_mass + V → ρ = nM/V"
            )

    def _calc_molar_mass(self, P=None, V=None, T=None, mass_g=None, **kwargs) -> dict:
        """Calculate molar mass from ideal gas data: M = mRT/(PV)"""
        if all(v is not None for v in [P, V, T, mass_g]):
            M = mass_g * self.R * T / (P * V)
            steps = f"Molar mass M = mRT/(PV) = ({mass_g})({self.R})({T})/(({P})({V})) = {round(M, 2)} g/mol"
            return {"result": round(M, 2)}
        else:
            raise ChemMCPInputError("For molar_mass mode provide: P, V, T, mass_g")

    def _combined_law(self, P1=None, V1=None, T1=None, P2=None, V2=None, T2=None, **kwargs) -> dict:
        """Combined gas law: P1V1/T1 = P2V2/T2"""
        known_1 = sum(1 for x in [P1, V1, T1] if x is not None)
        known_2 = sum(1 for x in [P2, V2, T2] if x is not None)

        if known_1 < 2 or known_2 < 2:
            raise ChemMCPInputError("Provide at least 2 of (P1,V1,T1) and 2 of (P2,V2,T2)")

        # Calculate left side
        left = (P1 or 1) * (V1 or 1) / (T1 or 1)
        actual_left = None
        if P1 and V1 and T1:
            actual_left = P1 * V1 / T1

        steps = f"Combined gas law: P₁V₁/T₁ = P₂V₂/T₂\n"
        given_1 = {k: v for k, v in [("P1", P1), ("V1", V1), ("T1", T1)] if v is not None}
        given_2 = {k: v for k, v in [("P2", P2), ("V2", V2), ("T2", T2)] if v is not None}

        if T2 is None and P2 and V2 and actual_left:
            result = P2 * V2 / actual_left
            unit = "K"
            steps += f"{actual_left:.4f} = ({P2})({V2})/T₂\nT₂ = ({P2}×{V2})/{actual_left:.4f} = {round(result, 2)} K"
        elif P2 is None and V2 and T2 and actual_left:
            result = actual_left * T2 / V2
            unit = "atm"
            steps += f"{actual_left:.4f} = P₂({V2})/{T2}\nP₂ = {actual_left:.4f}×{T2}/{V2} = {round(result, 4)} atm"
        elif V2 is None and P2 and T2 and actual_left:
            result = actual_left * T2 / P2
            unit = "L"
            steps += f"{actual_left:.4f} = ({P2})V₂/{T2}\nV₂ = {actual_left:.4f}×{T2}/{P2} = {round(result, 4)} L"
        else:
            raise ChemMCPInputError("Need to solve for exactly one unknown among P2, V2, T2")

        return {"result": round(result, 6),
                "unit": unit,
                "step_by_step": steps,
                "given_state1": given_1, "given_state2": given_2}

    def _partial_pressure(self, P_total=None, mole_fractions=None, n_total=None, ni=None, **kwargs) -> dict:
        """Dalton's law: Pi = xi × P_total = (ni/ntotal) × P_total"""
        if P_total is not None and mole_fractions is not None:
            results = []
            steps = f"Dalton's law of partial pressures: Pi = xi × P_total\nP_total = {P_total} atm\n\n"
            for i, xi in enumerate(mole_fractions):
                Pi = xi * P_total
                results.append({"component": i+1, "mole_fraction": xi, "partial_pressure": round(Pi, 4)})
                steps += f"  Component {i+1}: P{i+1} = {xi} × {P_total} = {round(Pi, 4)} atm\n"
            return {"result": results}
        elif P_total is not None and n_total is not None and ni is not None:
            xi = ni / n_total
            Pi = xi * P_total
            steps = f"Pi = (ni/ntotal) × P_total = ({ni}/{n_total}) × {P_total} = {round(Pi, 4)} atm"
            return {"result": round(Pi, 4), "mole_fraction": round(xi, 4)}
        else:
            raise ChemMCPInputError(
                "For partial_pressure mode provide:\n"
                "  (1) P_total + mole_fractions (list)\n"
                "  (2) P_total + n_total + ni"
            )

    def _run_text(self, query: str) -> dict:
        q = query.strip()
        # Try to parse simple expressions like "P=1 V=22.4 n=1 T=?"
        try:
            kwargs = {}
            for token in q.replace(",", " ").split():
                if "=" in token:
                    key, val = token.split("=", 1)
                    key = key.strip().upper()
                    if key in ("P1", "V1", "T1", "P2", "V2", "T2"):
                        kwargs[key.lower()] = float(val)
                    elif key == "MODE":
                        kwargs["mode"] = val.strip()
                    elif key in ("P", "V", "N", "T"):
                        if key == "N":
                            kwargs["n"] = float(val)
                        else:
                            kwargs[key.lower()] = float(val)
            if len(kwargs) >= 2:
                return self._run_base(**kwargs)
        except (ValueError, IndexError):
            pass

        raise ChemMCPInputError(
            f"Could not parse query '{q}'. Use format like:\n"
            f"  'P=1 V=22.4 n=1 T=?' for PV=nRT solve\n"
            f"  'P=1 T=273 M=44' for density\n"
            f"  'P1=1 V1=10 T1=300 P2=2 V2=5 T2=?' for combined law"
        )

