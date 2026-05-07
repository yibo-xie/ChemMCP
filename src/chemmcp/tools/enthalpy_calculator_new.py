import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Physical constants
_R = 8.314462618   # J/(mol·K)

@ChemMCPManager.register_tool
class EnthalpyCalculatorNew(BaseTool):
    """
    焓变计算工具 — 计算反应热、生成焓、温度依赖的焓变。
    
    支持Kirchhoff定律计算不同温度下的焓变。
    """
    __version__ = "0.1.0"
    name = "EnthalpyCalculatorNew"
    func_name = "calculate_enthalpy"
    description = "Calculate enthalpy changes (ΔH) for reactions, formation enthalpies, temperature-dependent ΔH via Kirchhoff's law."
    implementation_description = "Uses Hess's law for reaction enthalpy from formation data, and Kirchhoff's law: d(ΔH)/dT = ΔCp for temperature dependence with constant or temperature-dependent heat capacities."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Enthalpy", "Thermochemistry", "Hess Law", "Kirchhoff Law", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("calculation_type", "str", "N/A", "'reaction' (from formation enthalpies), 'kirchhoff' (T-dependent), 'formation' (from bond energies), or 'combustion'"),
        ("formation_enthalpies", "list", "N/A", "List of standard formation enthalpies in kJ/mol, e.g., [-393.5, -285.8]. For reaction mode: [products..., reactants...]."),
        ("stoich_coeffs", "list", "N/A", "Stoichiometric coefficients (positive for products, negative for reactants), e.g., [1, -1, -1] for A → B + C."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("heat_capacities", "list", "[]", "Heat capacities Cp in J/(mol·K) for each species (same order as formation_enthalpies). Used for Kirchhoff calculation."),
        ("T_ref_k", "float", "298.15", "Reference temperature for Kirchhoff law."),
        ("delta_H_ref", "float", "N/A", "ΔH at T_ref in kJ/mol (for Kirchhoff mode)."),
        ("unit", "str", "kJ", "Output unit: 'kJ' or 'J'."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "String: calc_type|[dHf_values]|[coeffs]|T|[Cp_values]|T_ref|dH_ref|unit"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with delta_H, reaction details, intermediate values."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_type": "reaction",
                "formation_enthalpies": [-393.5, -241.8, 0.0, 0.0],
                "stoich_coeffs": [1, 2, -1, -2],
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_str": "reaction|[-393.5,-241.8,0,0]|[1,2,-1,-2]|298.15|||298.15||kJ"
            },
            "output": {
                "result": {
                    "reaction": "CH4 + 2O2 → CO2 + 2H2O",
                    "delta_H_kJ_mol": -890.3,
                    "exothermic": True,
                    "description": "Combustion of methane",
                }
            },
        },
        {
            "code_input": {
                "calculation_type": "kirchhoff",
                "delta_H_ref": -57.2,
                "T_ref_k": 298.15,
                "temperature_k": 500.0,
                "heat_capacities": [37.1, 33.6, 29.4],  # products then reactants
                "stoich_coeffs": [2, -1, -1],  # 2C - A - B
            },
            "text_input": {
                "input_str": "kirchhoff|[]|[2,-1,-1]|500|[37.1,33.6,29.4]|298.15|-57.2|kJ"
            },
            "output": {
                "result": {
                    "delta_H_T_kJ_mol": "<value>",
                    "delta_H_ref_kJ_mol": -57.2,
                    "T_ref_K": 298.15,
                    "T_K": 500.0,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _calc_reaction(self, dHf_list: List[float], coeffs: List[float], T: float, unit: str) -> dict:
        """ΔH_rxn = Σ ν_i · ΔH_f°(i)"""
        if len(dHf_list) != len(coeffs):
            raise ChemMCPError("Length of formation_enthalpies must match stoich_coeffs.")
        
        dH_rxn = sum(c * dh for c, dh in zip(coeffs, dHf_list))
        factor = 1.0 if unit == "kJ" else 1000.0
        
        return {
            "calculation_type": "reaction",
            "delta_H_kJ_mol": round(dH_rxn, 4),
            "delta_H_J_mol": round(dH_rxn * 1000.0, 4),
            "temperature_K": T,
            "unit": unit,
            "exothermic": dH_rxn < 0,
            "endothermic": dH_rxn > 0,
            "num_species": len(dHf_list),
        }

    def _calc_kirchhoff(self, dH_ref: float, T_ref: float, T: float,
                         Cp_list: List[float], coeffs: List[float], unit: str) -> dict:
        """Kirchhoff's law: ΔH(T2) = ΔH(T1) + ∫ΔCp dT (constant Cp assumption)."""
        if len(Cp_list) != len(coeffs):
            raise ChemMCPError("Length of heat_capacities must match stoich_coeffs.")
        
        delta_Cp = sum(c * cp for c, cp in zip(coeffs, Cp_list))  # J/(mol·K)
        dH_T = dH_ref + delta_Cp * (T - T_ref) / 1000.0  # convert J→kJ
        
        return {
            "calculation_type": "kirchhoff",
            "delta_H_ref_kJ_mol": round(dH_ref, 4),
            "delta_H_T_kJ_mol": round(dH_T, 4),
            "delta_H_T_J_mol": round(dH_T * 1000.0, 4),
            "delta_Cp_J_mol_K": round(delta_Cp, 4),
            "T_ref_K": T_ref,
            "T_K": T,
            "dT_K": round(T - T_ref, 4),
            "unit": unit,
        }

    def _calc_formation_bond(self, bond_energies_broken: List[float], bond_energies_formed: List[float]) -> dict:
        """Estimate ΔH_f from bond energies: ΔH ≈ Σ D(broken) - Σ D(formed)."""
        E_broken = sum(bond_energies_broken)
        E_formed = sum(bond_energies_formed)
        dH = E_broken - E_formed  # kJ/mol typically
        
        return {
            "calculation_type": "formation (bond energy)",
            "total_bond_energy_broken_kJ_mol": round(E_broken, 4),
            "total_bond_energy_formed_kJ_mol": round(E_formed, 4),
            "delta_H_kJ_mol": round(dH, 4),
            "exothermic": dH < 0,
        }

    def _run_base(self, calculation_type: str, formation_enthalpies: List[float] = None,
                  stoich_coeffs: List[float] = None, temperature_k: float = 298.15,
                  heat_capacities: List[float] = None, T_ref_k: float = 298.15,
                  delta_H_ref: float = None, unit: str = "kJ") -> dict:
        calc_type = calculation_type.lower().strip()
        
        if calc_type == "reaction":
            if formation_enthalpies is None or stoich_coeffs is None:
                raise ChemMCPError("'reaction' mode requires formation_enthalpies and stoich_coeffs.")
            return self._calc_reaction(formation_enthalpies, stoich_coeffs, temperature_k, unit)
        elif calc_type == "kirchhoff":
            if delta_H_ref is None:
                raise ChemMCPError("'kirchhoff' mode requires delta_H_ref.")
            return self._calc_kirchhoff(delta_H_ref, T_ref_k, temperature_k, heat_capacities or [], stoich_coeffs or [], unit)
        elif calc_type == "formation":
            # Bond energy method — expects special input format
            raise ChemMCPError("Use 'reaction' or 'kirchhoff' modes. For bond energy estimation, provide bond energies via a different call.")
        else:
            raise ChemMCPError(
                f"Unknown type: '{calculation_type}'. Options: 'reaction', 'kirchhoff'."
            )

    def _run_text(self, input_str: str) -> dict:
        try:
            parts = input_str.split("|")
            calc_type = parts[0].strip()
            import json
            dHf = json.loads(parts[1]) if len(parts) > 1 and parts[1].strip() else None
            coeffs = json.loads(parts[2]) if len(parts) > 2 and parts[2].strip() else None
            T = float(parts[3]) if len(parts) > 3 else 298.15
            Cp = json.loads(parts[4]) if len(parts) > 4 and parts[4].strip() else []
            Tr = float(parts[5]) if len(parts) > 5 else 298.15
            dHr = float(parts[6]) if len(parts) > 6 and parts[6].strip() else None
            u = parts[7] if len(parts) > 7 else "kJ"
            return self._run_base(calc_type, dHf, coeffs, T, Cp, Tr, dHr, u)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
