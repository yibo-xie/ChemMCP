import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)

# Isotope data: element -> list of (mass_number, abundance_percent, stability, half_life)
ISOTOPE_DATA: dict = {
    "H":  [(1, 99.9885, "stable", None), (2, 0.0115, "stable", None), (3, None, "radioactive", "12.32 y")],
    "He": [(3, 0.000137, "stable", None), (4, 99.999863, "stable", None), (6, None, "radioactive", "806.7 ms")],
    "Li": [(6, 7.59, "stable", None), (7, 92.41, "stable", None), (8, None, "radioactive", "838 ms")],
    "Be": [(9, 100, "stable", None), (10, None, "radioactive", "1.51e6 y")],
    "B":  [(10, 19.9, "stable", None), (11, 80.1, "stable", None)],
    "C":  [(12, 98.93, "stable", None), (13, 1.07, "stable", None), (14, None, "radioactive", "5730 y")],
    "N":  [(14, 99.634, "stable", None), (15, 0.366, "stable", None), (13, None, "radioactive", "9.97 min")],
    "O":  [(16, 99.757, "stable", None), (17, 0.038, "stable", None), (18, 0.205, "stable", None)],
    "F":  [(19, 100, "stable", None), (18, None, "radioactive", "109.77 min")],
    "Ne": [(20, 90.48, "stable", None), (21, 0.27, "stable", None), (22, 9.25, "stable", None)],
    "Na": [(23, 100, "stable", None), (22, None, "radioactive", "2.602 y")],
    "Mg": [(24, 78.99, "stable", None), (25, 10.00, "stable", None), (26, 11.01, "stable", None)],
    "Al": [(27, 100, "stable", None), (26, None, "radioactive", "7.17e5 y")],
    "Si": [(28, 92.23, "stable", None), (29, 4.68, "stable", None), (30, 3.09, "stable", None)],
    "P":  [(31, 100, "stable", None), (32, None, "radioactive", "14.26 d")],
    "S":  [(32, 94.99, "stable", None), (33, 0.75, "stable", None), (34, 4.25, "stable", None), (36, 0.01, "radioactive", "1.3e8 y")],
    "Cl": [(35, 75.76, "stable", None), (37, 24.24, "stable", None)],
    "Ar": [(36, 0.334, "stable", None), (38, 0.063, "stable", None), (40, 99.603, "stable", None)],
    "K":  [(39, 93.2581, "stable", None), (40, 0.0117, "radioactive", "1.248e9 y"), (41, 6.7302, "stable", None)],
    "Ca": [(40, 96.941, "stable", None), (42, 0.647, "stable", None), (43, 0.135, "stable", None), (44, 2.086, "stable", None), (46, 0.004, "stable", None), (48, 0.187, "radioactive", "6.4e18 y")],
    "Fe": [(54, 5.845, "stable", None), (56, 91.754, "stable", None), (57, 2.119, "stable", None), (58, 0.282, "stable", None)],
    "Cu": [(63, 69.17, "stable", None), (65, 30.83, "stable", None)],
    "Zn": [(64, 48.63, "stable", None), (66, 27.90, "stable", None), (67, 4.10, "stable", None), (68, 18.75, "stable", None), (70, 0.62, "stable", None)],
    "Br": [(79, 50.69, "stable", None), (81, 49.31, "stable", None)],
    "Sr": [(84, 0.56, "stable", None), (86, 9.86, "stable", None), (87, 7.00, "stable", None), (88, 82.58, "stable", None)],
    "Ag": [(107, 51.84, "stable", None), (109, 48.16, "stable", None)],
    "I":  [(127, 100, "stable", None), (129, None, "radioactive", "1.57e7 y")],
    "Ba": [(130, 0.106, "stable", None), (132, 0.101, "stable", None), (134, 2.417, "stable", None), (135, 6.592, "stable", None), (136, 7.854, "stable", None), (137, 11.232, "stable", None), (138, 71.698, "stable", None)],
    "Pb": [(204, 1.4, "stable", None), (206, 24.1, "stable", None), (207, 22.1, "stable", None), (208, 52.4, "stable", None)],
    "U":  [(234, 0.0055, "radioactive", "2.455e5 y"), (235, 0.72, "radioactive", "7.04e8 y"), (238, 99.2745, "radioactive", "4.468e9 y")],
}


@ChemMCPManager.register_tool
class GetIsotopes(BaseTool):
    __version__ = "0.1.0"
    name = "GetIsotopes"
    func_name = 'get_isotopes'
    description = "Get isotope information for an element including mass number, natural abundance, stability, and half-life."
    implementation_description = "Uses a built-in isotope database with data for common elements. Returns mass number, percent abundance, stability status, and half-life for radioactive isotopes."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Isotopes", "Nuclear Chemistry"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., C, U, Fe)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('element', 'str', 'Element symbol'),
        ('isotopes', 'list', 'List of isotope data dictionaries'),
        ('isotope_count', 'int', 'Number of isotopes in database'),
    ]
    
    examples = [
        {'code_input': {'element': 'C'}, 'text_input': {'element': 'C'}, 'output': {'element': 'C', 'isotopes': [{...}, {...}, {...}], 'isotope_count': 3}},
        {'code_input': {'element': 'U'}, 'text_input': {'element': 'U'}, 'output': {'element': 'U', 'isotopes': [{...}, {...}, {...}], 'isotope_count': 3}},
    ]
    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        sym = data["symbol"]
        if sym not in ISOTOPE_DATA:
            return {
                "element": sym,
                "isotopes": [],
                "isotope_count": 0,
                "note": f"Isotope data not available for {sym} in this database.",
            }
        isotopes = []
        for mass_num, abundance, stability, half_life in ISOTOPE_DATA[sym]:
            iso = {
                "mass_number": mass_num,
                "abundance_percent": round(abundance, 4) if abundance is not None else None,
                "stability": stability,
                "half_life": half_life,
            }
            if stability == "stable":
                iso["isotope_label"] = f"{sym}-{mass_num}"
            else:
                iso["isotope_label"] = f"{sym}-{mass_num} (radioactive)"
            isotopes.append(iso)

        stable_count = sum(1 for i in isotopes if i["stability"] == "stable")
        return {
            "element": sym,
            "isotopes": isotopes,
            "isotope_count": len(isotopes),
            "stable_isotopes": stable_count,
            "radioactive_isotopes": len(isotopes) - stable_count,
        }


if __name__ == "__main__":
    run_mcp_server()
