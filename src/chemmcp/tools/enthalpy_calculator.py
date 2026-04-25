import logging
from typing import Dict, Optional
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Standard enthalpies of formation ΔHf° (kJ/mol) at 298 K for common substances
_STANDARD_FORMATION_ENTHALPY = {
    # Gases (g)
    "H2(g)": 0.0, "O2(g)": 0.0, "N2(g)": 0.0, "C(graphite)": 0.0,
    "CO(g)": -110.5, "CO2(g)": -393.5, "H2O(g)": -241.8,
    "NH3(g)": -46.1, "CH4(g)": -74.8, "C2H6(g)": -84.7,
    "C2H4(g)": 52.4, "C2H2(g)": 226.7, "HCl(g)": -92.3,
    "SO2(g)": -296.8, "SO3(g)": -395.7, "NO(g)": 90.25,
    "NO2(g)": 33.2, "H2S(g)": -20.6,
    # Liquids (l)
    "H2O(l)": -285.8, "C6H6(l)": 49.0, "CH3OH(l)": -238.7,
    "C2H5OH(l)": -277.6, "C6H12O6(glucose,l)": -1273.3,
    # Solids (s)
    "CaCO3(s)": -1206.9, "CaO(s)": -635.5, "Fe2O3(s)": -824.2,
    "NaCl(s)": -411.1, "AgCl(s)": -127.0, "Al2O3(s)": -1675.7,
}

# Standard combustion enthalpies (kJ/mol) at 298 K
_COMBUSTION_ENTHALPY = {
    "H2": -285.8, "C(graphite)": -393.5, "S(s)": -296.8,
    "CH4": -890.3, "C2H6": -1559.8, "C2H4": -1411.0,
    "C2H2": -1299.5, "CH3OH": -726.1, "C2H5OH": -1366.8,
    "C6H12O6(glucose)": -2803.0,
}


@ChemMCPManager.register_tool
class EnthalpyCalculator(BaseTool):
    """
    计算反应焓变（标准生成焓法、燃烧焓法）。
    """
    __version__ = "0.1.0"
    name = "EnthalpyCalculator"
    func_name = "calculate_reaction_enthalpy"
    description = "Calculate reaction enthalpy change using standard formation enthalpies or combustion enthalpies."
    implementation_description = "Supports two modes: 'formation' uses ΣΔHf°(products) - ΣΔHf°(reactants); 'combustion' uses ΣΔHc°(reactants) - ΣΔHc°(products). Includes built-in data for common substances."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Enthalpy", "Thermochemistry", "Reaction Heat"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "formation", "Calculation mode: 'formation' or 'combustion'."),
        ("species_data", "str", "N/A", "JSON string of species with stoichiometry and values. Format: '{\"species\": \"coeff:value,...\"}'. For formation mode: values are auto-looked up; for combustion mode too."),
        ("custom_values", "str", "{}", "Optional JSON dict to override built-in enthalpy values. Format: '{\"H2O(l)\": -285.8}' in kJ/mol."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: mode species_json [custom_values_json]. Example: 'formation \"{\\\"CO2(g)\\\":1,\\\"H2O(l)\\\":2,\\\"C3H8(g)\\\":-1,\\\"O2(g)\\\":-5}\"'"),
    ]

    output_sig = [
        ("delta_h", "float", "Reaction enthalpy change ΔH in kJ/mol."),
        ("explanation", "str", "Step-by-step calculation explanation."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "formation",
                "species_data": '{"CO2(g)": 3, "H2O(l)": 4, "C3H8(g)": -1, "O2(g)": -5}',
                "custom_values": "{}",
            },
            "text_input": {
                "input_params": "formation {\"CO2(g)\":3,\"H2O(l)\":4,\"C3H8(g)\":-1,\"O2(g)\":-5}",
            },
            "output": {
                "delta_h": -2220.0,
                "explanation": "Combustion of propane: C3H8 + 5O2 → 3CO2 + 4H2O; ΔH = [3(-393.5)+4(-285.8)] - [-103.85+5(0)] = -2220.0 kJ/mol",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._formation_db = dict(_STANDARD_FORMATION_ENTHALPY)
        self._combustion_db = dict(_COMBUSTION_ENTHALPY)

    def _run_base(self, mode: str, species_data: str, custom_values: str = "{}") -> dict:
        """Core logic: calculate reaction enthalpy."""
        import json

        mode = mode.lower().strip()
        if mode not in ("formation", "combustion"):
            raise ChemMCPError(f"Invalid mode '{mode}'. Use 'formation' or 'combustion'.")

        try:
            species = json.loads(species_data)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON in species_data: {species_data}")

        try:
            overrides = json.loads(custom_values) if custom_values else {}
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON in custom_values: {custom_values}")

        # Apply overrides
        db = self._formation_db if mode == "formation" else self._combustion_db
        working_db = dict(db)
        working_db.update(overrides)

        # Calculate
        products_sum = 0.0
        reactants_sum = 0.0
        details = []

        for sp, coeff in species.items():
            sp_key = sp.strip()
            if sp_key not in working_db:
                raise ChemMCPError(f"Unknown substance '{sp_key}' not in database. Provide it via custom_values.")
            value = working_db[sp_key]
            contribution = coeff * value
            if coeff > 0:
                products_sum += contribution
                details.append(f"  +{coeff}×{sp_key}: {coeff} × ({value}) = {contribution:.2f}")
            else:
                reactants_sum += contribution
                details.append(f"  {coeff}×{sp_key}: {coeff} × ({value}) = {contribution:.2f}")

        if mode == "formation":
            delta_h = products_sum - reactants_sum
            formula_desc = "ΣΔHf°(products) - ΣΔHf°(reactants)"
        else:
            delta_h = reactants_sum - products_sum
            formula_desc = "ΣΔHc°(reactants) - ΣΔHc°(products)"

        detail_str = "\n".join(details)
        explanation = f"{formula_desc}\n{detail_str}\nΔH = {delta_h:.2f} kJ/mol"

        logger.info(f"Enthalpy calculation (mode={mode}): ΔH = {delta_h:.2f} kJ/mol")
        return {"delta_h": round(delta_h, 2), "explanation": explanation}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split(None, 2)
            mode = parts[0]
            species_data = parts[1] if len(parts) > 1 else "{}"
            custom_values = parts[2] if len(parts) > 2 else "{}"
            return self._run_base(mode, species_data, custom_values)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'mode species_json [custom_values_json]'")
