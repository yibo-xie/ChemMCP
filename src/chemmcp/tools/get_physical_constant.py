import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


PHYSICAL_CONSTANTS = {
    "NA":          {"name": "Avogadro's Number", "value": 6.02214076e23,  "unit": "mol⁻¹",       "symbol": "Nₐ"},
    "R":           {"name": "Ideal Gas Constant",   "value": 8.314462618,    "unit": "J/(mol·K)",   "symbol": "R"},
    "R_Latm":      {"name": "Gas Constant (L·atm)", "value": 0.08205736608096596, "unit": "L·atm/(mol·K)", "symbol": "R"},
    "F":           {"name": "Faraday Constant",     "value": 96485.33212,     "unit": "C/mol",       "symbol": "F"},
    "h":           {"name": "Planck Constant",      "value": 6.62607015e-34,  "unit": "J·s",         "symbol": "h"},
    "hbar":        {"name": "Reduced Planck Constant", "value": 1.054571817e-34, "unit": "J·s",     "symbol": "ℏ"},
    "c":           {"name": "Speed of Light",       "value": 299792458,       "unit": "m/s",         "symbol": "c"},
    "e_charge":    {"name": "Elementary Charge",    "value": 1.602176634e-19, "unit": "C",           "symbol": "e"},
    "me":          {"name": "Electron Mass",        "value": 9.1093837015e-31,"unit": "kg",          "symbol": "mₑ"},
    "mp":          {"name": "Proton Mass",          "value": 1.67262192369e-27,"unit":"kg",           "symbol": "mₚ"},
    "mn":          {"name": "Neutron Mass",         "value": 1.6749274986e-27,"unit":"kg",            "symbol": "mₙ"},
    "amu":         {"name": "Atomic Mass Unit",     "value": 1.66053906660e-27,"unit":"kg",           "symbol": "u"},
    "G":           {"name": "Gravitational Constant","value": 6.67430e-11,    "unit": "m³/(kg·s²)",  "symbol": "G"},
    "k_B":         {"name": "Boltzmann Constant",   "value": 1.380649e-23,    "unit": "J/K",         "symbol": "k_B"},
    "sigma":       {"name": "Stefan-Boltzmann Const","value": 5.670374419e-8, "unit": "W/(m²·K⁴)",  "symbol": "σ"},
    "epsilon_0":   {"name": "Vacuum Permittivity",  "value": 8.8541878128e-12,"unit":"F/m",            "symbol": "ε₀"},
    "mu_0":        {"name": "Vacuum Permeability",  "value": 1.25663706212e-6,"unit":"N/A²",          "symbol": "μ₀"},
    "N_A":         {"name": "Avogadro's Number",    "value": 6.02214076e23,  "unit": "mol⁻¹",       "symbol": "Nₐ"},
    "eV":          {"name": "Electron Volt",        "value": 1.602176634e-19, "unit": "J",           "symbol": "eV"},
    "amu_MeV":     {"name": "AMU in MeV/c²",       "value": 931.49410242,    "unit": "MeV/c²",      "symbol": "u→MeV"},
    "Rydberg":     {"name": "Rydberg Constant",     "value": 1.0973731568160e7,"unit":"m⁻¹",          "symbol": "R∞"},
    "Bohr_radius": {"name": "Bohr Radius",          "value": 5.29177210903e-11,"unit":"m",            "symbol": "a₀"},
    "std_pressure":{"name": "Standard Pressure (atm)","value": 1.01325e5,       "unit": "Pa",          "symbol": "P°"},
    "std_volume":  {"name": "Molar Volume (STP)",   "value": 22.41396954,      "unit": "L/mol",       "symbol": "Vm"},
    "calorie":     {"name": "Thermochemical Calorie","value": 4.184,            "unit": "J",           "symbol": "cal"},
}


@ChemMCPManager.register_tool
class GetPhysicalConstant(BaseTool):
    """
    物理常数查询工具。
    查询常用物理常数（阿伏伽德罗常数、气体常数、法拉第常数、普朗克常数等）。
    """
    __version__      = "0.1.0"
    name             = "GetPhysicalConstant"
    func_name        = "get_physical_constant"
    description      = "Query physical constants used in chemistry and physics: Avogadro number, gas constant, Faraday constant, Planck constant, etc."
    implementation_description = "Built-in database of fundamental physical constants with values, units, symbols, and descriptions."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Physical Constants", "Reference Data", "Fundamental"]
    required_envs    = []

    code_input_sig   = [
        ("constant_name", "str", "N/A", "Name or symbol of the constant to query (e.g., 'NA', 'R', 'F', 'h', 'c', 'all')."),
    ]

    text_input_sig   = [
        ("query_str", "str", "N/A", "Constant name or symbol as string."),
    ]

    output_sig       = [
        ("constant_name", "str", "The queried constant identifier."),
        ("full_name", "str", "Full descriptive name of the constant."),
        ("value", "float", "Numerical value of the constant."),
        ("unit", "str", "Unit of the constant."),
        ("symbol", "str", "Symbol notation of the constant."),
        ("description", "str", "Brief description and common usage context."),
    ]

    examples         = [
        {
            "code_input": {"constant_name": "NA"},
            "text_input": {"query_str": "NA"},
            "output": {
                "constant_name": "NA",
                "full_name": "Avogadro's Number",
                "value": 6.02214076e23,
                "unit": "mol⁻¹",
                "symbol": "Nₐ",
                "description": "Number of constituent particles (atoms, molecules, etc.) per mole.",
            }
        },
        {
            "code_input": {"constant_name": "all"},
            "text_input": {"query_str": "all"},
            "output": {
                "constant_name": "All Constants",
                "full_name": "All Available Physical Constants",
                "value": 28,
                "unit": "constants",
                "symbol": "-",
                "description": f"Available constants: {', '.join(sorted(PHYSICAL_CONSTANTS.keys()))}",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, constant_name: str) -> dict:
        """Query a physical constant by name or symbol."""
        key = constant_name.strip()

        if key.lower() == "all":
            result_list = []
            for ck, cv in sorted(PHYSICAL_CONSTANTS.items()):
                result_list.append({
                    "constant_name": ck,
                    "full_name": cv["name"],
                    "value": cv["value"],
                    "unit": cv["unit"],
                    "symbol": cv["symbol"],
                    "description": self._get_desc(ck),
                })
            return {
                "constant_name": "All Constants",
                "full_name": f"All {len(PHYSICAL_CONSTANTS)} Physical Constants",
                "value": len(PHYSICAL_CONSTANTS),
                "unit": "constants",
                "symbol": "-",
                "description": f"Available: {', '.join(sorted(PHYSICAL_CONSTANTS.keys()))}",
                "_all_constants": result_list,
            }

        # Try exact match
        if key in PHYSICAL_CONSTANTS:
            return self._format_constant(key, PHYSICAL_CONSTANTS[key])

        # Case-insensitive match
        for ck in PHYSICAL_CONSTANTS:
            if ck.lower() == key.lower():
                return self._format_constant(ck, PHYSICAL_CONSTANTS[ck])

        # Partial match
        matches = [k for k in PHYSICAL_CONSTANTS if key.lower() in k.lower()]
        if len(matches) == 1:
            return self._format_constant(matches[0], PHYSICAL_CONSTANTS[matches[0]])
        elif len(matches) > 1:
            raise ChemMCPError(f"Multiple matches for '{key}': {matches}. Be more specific.")

        available = ", ".join(sorted(PHYSICAL_CONSTANTS.keys()))
        raise ChemMCPError(f"Constant '{key}' not found. Available: {available}")

    def _format_constant(self, key: str, data: dict) -> dict:
        return {
            "constant_name": key,
            "full_name": data["name"],
            "value": data["value"],
            "unit": data["unit"],
            "symbol": data["symbol"],
            "description": self._get_desc(key),
        }

    def _get_desc(self, key: str) -> str:
        descs = {
            "NA": "Number of particles per mole. Used in stoichiometry, molar mass calculations.",
            "R": "Ideal gas law constant. PV = nRT. Used in thermodynamics, gas laws.",
            "F": "Charge per mole of electrons. Used in electrochemistry: ΔG° = -nFE°.",
            "h": "Quantum of action. E = hν. Fundamental in quantum mechanics.",
            "c": "Universal speed limit. Used in relativity, E=mc².",
            "k_B": "Relates temperature to energy. E = k_B·T. Used in statistical mechanics.",
            "e_charge": "Fundamental electric charge. Charge of one proton/electron magnitude.",
            "amu": "Standard atomic mass unit. 1/12 of C-12 atom mass.",
        }
        return descs.get(key, "A fundamental physical constant used in chemistry and physics.")

    def _run_text(self, query_str: str) -> dict:
        return self._run_base(query_str.strip())
