import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)

# Abundance data: crust (ppm), ocean (mg/L), universe (relative to Si=1e6)
ABUNDANCE_DATA: dict = {
    "H":  {"crust": 1400, "ocean": 107800, "universe": 2.8e10},
    "He": {"crust": 0.008, "ocean": 0.000007, "universe": 2.3e9},
    "Li": {"crust": 17, "ocean": 0.18, "universe": 60},
    "Be": {"crust": 2.0, "ocean": 0.0000056, "universe": 0.73},
    "B":  {"crust": 10, "ocean": 4.44, "universe": 10},
    "C":  {"crust": 200, "ocean": 28, "universe": 1e7},
    "N":  {"crust": 19, "ocean": 0.5, "universe": 3.1e6},
    "O":  {"crust": 461000, "ocean": 857000, "universe": 2.38e7},
    "F":  {"crust": 585, "ocean": 1.3, "universe": 840},
    "Ne": {"crust": None, "ocean": 0.00012, "universe": 1.34e6},
    "Na": {"crust": 23600, "ocean": 10770, "universe": 6e4},
    "Mg": {"crust": 23300, "ocean": 1290, "universe": 1e6},
    "Al": {"crust": 82300, "ocean": 0.002, "universe": 8.5e4},
    "Si": {"crust": 282000, "ocean": 2.98, "universe": 1e6},
    "P":  {"crust": 1050, "ocean": 0.06, "universe": 8400},
    "S":  {"crust": 350, "ocean": 905, "universe": 5e5},
    "Cl": {"crust": 145, "ocean": 19345, "universe": 5200},
    "Ar": {"crust": None, "ocean": 0.45, "universe": 1e5},
    "K":  {"crust": 20900, "ocean": 392, "universe": 3700},
    "Ca": {"crust": 41000, "ocean": 412, "universe": 6200},
    "Ti": {"crust": 5650, "ocean": 0.001, "universe": 2400},
    "V":  {"crust": 120, "ocean": 0.0025, "universe": 290},
    "Cr": {"crust": 102, "ocean": 0.0003, "universe": 1400},
    "Mn": {"crust": 950, "ocean": 0.002, "universe": 930},
    "Fe": {"crust": 56300, "ocean": 0.002, "universe": 1e6},
    "Co": {"crust": 25, "ocean": 0.00002, "universe": 2300},
    "Ni": {"crust": 84, "ocean": 0.0056, "universe": 8.4e4},
    "Cu": {"crust": 60, "ocean": 0.25, "universe": 480},
    "Zn": {"crust": 70, "ocean": 0.00494, "universe": 1300},
    "Br": {"crust": 2.4, "ocean": 67.3, "universe": 13},
    "Rb": {"crust": 90, "ocean": 0.12, "universe": 7},
    "Sr": {"crust": 370, "ocean": 7.9, "universe": 41},
    "Ag": {"crust": 0.075, "ocean": 0.00003, "universe": 0.55},
    "I":  {"crust": 0.45, "ocean": 0.06, "universe": 1},
    "Au": {"crust": 0.0024, "ocean": 0.000004, "universe": 0.14},
    "Pb": {"crust": 14, "ocean": 0.003, "universe": 1},
    "U":  {"crust": 2.7, "ocean": 0.0033, "universe": 0.01},
}


@ChemMCPManager.register_tool
class ElementAbundance(BaseTool):
    __version__ = "0.1.0"
    name = "ElementAbundance"
    func_name = 'element_abundance'
    description = "Query element abundance in Earth's crust, oceans, and universe."
    implementation_description = "Returns abundance data from geochemical surveys and astrophysical measurements. Crust abundance in ppm (parts per million by mass), ocean concentration in mg/L, cosmic abundance relative to silicon = 10^6."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Abundance", "Geochemistry", "Astrochemistry"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., O, Fe, Au)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('element', 'str', 'Element symbol'),
        ('crust_abundance_ppm', 'float', 'Abundance in Earth\'s crust (ppm)'),
        ('ocean_concentration_mg_L', 'float', 'Concentration in ocean water (mg/L)'),
        ('cosmic_abundance', 'float', 'Cosmic abundance (Si=1e6)'),
        ('rank_note', 'str', 'Ranking context'),
    ]
    
        
    examples = [
        {'code_input': {'element': 'Fe'}, 'text_input': {'element': 'Fe'}, 'output': {'element': 'Fe', 'crust_abundance_ppm': 56300, 'ocean_concentration_mg_L': None, 'cosmic_abundance': None, 'rank_note': '...'}},
        {'code_input': {'element': 'Au'}, 'text_input': {'element': 'Au'}, 'output': {'element': 'Au', 'crust_abundance_ppm': 0.004, 'ocean_concentration_mg_L': None, 'cosmic_abundance': None, 'rank_note': '...'}},
    ]
    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        sym = data["symbol"]
        if sym not in ABUNDANCE_DATA:
            return {
                "element": sym,
                "note": f"Abundance data not available for {sym} in this database.",
            }
        ab = ABUNDANCE_DATA[sym]
        rank_notes = {
            "O": "Most abundant element in Earth's crust (~46% by mass).",
            "Si": "Second most abundant element in Earth's crust.",
            "Al": "Third most abundant element in Earth's crust.",
            "Fe": "Fourth most abundant element in Earth's crust; most abundant transition metal.",
            "Au": "Rare precious metal; ~0.0024 ppm in crust.",
            "He": "Second most abundant element in the universe but rare on Earth (light, escapes gravity).",
        }
        result = {
            "element": sym,
            "crust_abundance_ppm": ab["crust"],
            "ocean_concentration_mg_L": ab["ocean"],
            "cosmic_abundance": ab["universe"],
            "rank_note": rank_notes.get(sym, ""),
        }
        # Add interpretation
        if ab["crust"] is not None and ab["crust"] > 10000:
            result["abundance_category"] = "major element (>10,000 ppm)"
        elif ab["crust"] is not None and ab["crust"] > 100:
            result["abundance_category"] = "minor/trace element"
        elif ab["crust"] is not None:
            result["abundance_category"] = "ultra-trace element (<100 ppm)"
        else:
            result["abundance_category"] = "noble gas / not applicable to crust"
        return result


if __name__ == "__main__":
    run_mcp_server()
