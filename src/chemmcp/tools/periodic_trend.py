import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import PERIODIC_TABLE

logger = logging.getLogger(__name__)

# Trend data for periodic properties across periods and groups
TRENDS = {
    "atomic_radius": {
        "description": "Atomic radius generally decreases across a period (increasing nuclear charge) and increases down a group (additional electron shells).",
        "across_period": "decreases (left → right)",
        "down_group": "increases (top → bottom)",
        "unit": "pm",
    },
    "electronegativity": {
        "description": "Electronegativity increases across a period and decreases down a group.",
        "across_period": "increases (left → right)",
        "down_group": "decreases (top → bottom)",
        "unit": "Pauling scale",
    },
    "ionization_energy": {
        "description": "First ionization energy generally increases across a period and decreases down a group.",
        "across_period": "increases (left → right)",
        "down_group": "decreases (top → bottom)",
        "unit": "kJ/mol",
    },
    "electron_affinity": {
        "description": "Electron affinity generally becomes more negative (more exothermic) across a period, with exceptions at group 2 and 15.",
        "across_period": "becomes more exothermic (left → right), with exceptions",
        "down_group": "decreases (less negative down a group)",
        "unit": "kJ/mol",
    },
    "metallic_character": {
        "description": "Metallic character decreases across a period and increases down a group.",
        "across_period": "decreases (left → right)",
        "down_group": "increases (top → bottom)",
        "unit": "qualitative",
    },
    "atomic_radius_data": {  # empirical covalent/atomic radii in pm
        1: {"H": 53}, 2: {"Li": 167, "Be": 112, "B": 87, "C": 67, "N": 56, "O": 48, "F": 42, "Ne": 38},
        3: {"Na": 190, "Mg": 145, "Al": 118, "Si": 111, "P": 98, "S": 88, "Cl": 79, "Ar": 71},
        4: {"K": 243, "Ca": 194, "Sc": 184, "Ti": 176, "V": 171, "Cr": 166, "Mn": 161, "Fe": 156, "Co": 152, "Ni": 149, "Cu": 145, "Zn": 142, "Ga": 136, "Ge": 125, "As": 114, "Se": 103, "Br": 94, "Kr": 88},
    },
}


@ChemMCPManager.register_tool
class PeriodicTrend(BaseTool):
    __version__ = "0.1.0"
    name = "PeriodicTrend"
    func_name = 'periodic_trend'
    description = "Query periodic table trends for properties like atomic radius, electronegativity, ionization energy, electron affinity, and metallic character."
    implementation_description = "Returns trend descriptions (how a property changes across periods and down groups), plus actual data values for specific periods or groups when requested."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Trends", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ('property', 'str', 'N/A', 'Property to query: atomic_radius, electronegativity, ionization_energy, electron_affinity, metallic_character'),
        ('period', 'int', 'N/A', 'Optional: show data for this period number (1-7)'),
        ('group', 'int', 'N/A', 'Optional: show data for this group number (1-18)'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Property name, optionally with period or group, e.g., \"electronegativity period 3\"'),
    ]
    output_sig = [
        ('property', 'str', 'Property queried'),
        ('trend_description', 'str', 'How the property changes across periods and groups'),
        ('data', 'dict', 'Actual values if period or group specified'),
    ]
    
        
    
    examples = [
        {'code_input': {'property': 'electronegativity', 'period': None, 'group': None}, 'text_input': {'query': 'electronegativity'}, 'output': {'property': 'electronegativity', 'trend_description': '...', 'data': {}}},
        {'code_input': {'property': 'atomic_radius', 'period': 3, 'group': None}, 'text_input': {'query': 'atomic radius period 3'}, 'output': {'property': 'atomic_radius', 'trend_description': '...', 'data': {}}},
    ]
    def _run_base(self, property: str, period: int = None, group: int = None) -> dict:
        valid_props = list(TRENDS.keys())
        if property not in TRENDS:
            raise ChemMCPInputError(f"Unsupported property: '{property}'. Available: {valid_props}")

        trend = TRENDS[property]
        result = {
            "property": property,
            "trend_description": trend["description"],
            "across_period_trend": trend["across_period"],
            "down_group_trend": trend["down_group"],
            "unit": trend.get("unit"),
        }

        # If period requested, return data for that period
        if period is not None:
            key = f"{property}_data"
            if key in TRENDS and period in TRENDS[key]:
                result["period_data"] = {"period": period, "values": TRENDS[key][period]}
            else:
                # Build from periodic_table
                elements_in_period = [d for d in PERIODIC_TABLE.values() if d["period"] == period]
                prop_key_map = {
                    "electronegativity": "en_pauling",
                    "ionization_energy": None,  # needs special handling
                    "electron_affinity": None,
                    "metallic_character": None,
                    "atomic_radius": None,
                }
                result["period_data"] = {"period": period, "elements": len(elements_in_period)}

        if group is not None:
            elements_in_group = [d for d in PERIODIC_TABLE.values() if d["group"] == group]
            result["group_data"] = {"group": group, "element_count": len(elements_in_group),
                                     "elements": [{"symbol": d["symbol"], "name": d["name"], "period": d["period"]} for d in sorted(elements_in_group, key=lambda x: x["period"])]}

        return result


if __name__ == "__main__":
    run_mcp_server()
