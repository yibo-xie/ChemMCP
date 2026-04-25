import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)

# Ionization energies in kJ/mol (NIST data)
IONIZATION_ENERGIES: dict = {
    "H":  [1312],
    "He": [2372, 5250],
    "Li": [520, 7298, 11815],
    "Be": [900, 1757, 14849],
    "B":  [801, 2427, 3660],
    "C":  [1086, 2353, 4621, 6223],
    "N":  [1402, 2856, 4579, 7475, 9445],
    "O":  [1314, 3388, 5300, 7469, 10990],
    "F":  [1681, 3374, 6050, 8408, 11222],
    "Ne": [2081, 3952, 6122, 9370, 12178],
    "Na": [496, 4563, 6913, 9544, 13354, 16613, 20117, 25496],
    "Mg": [738, 1451, 7733, 10543, 13630, 18020, 21711, 25661],
    "Al": [577, 1817, 2745, 11578, 14842, 18379, 23329, 27466],
    "Si": [787, 1577, 3232, 4356, 16091, 19800, 23780, 29294],
    "P":  [1012, 1907, 2914, 4964, 6274, 21267, 25431, 29883],
    "S":  [1000, 2252, 3357, 4557, 7004, 8496, 27107, 31720],
    "Cl": [1251, 2298, 3822, 5159, 6543, 9362, 11018, 33604],
    "Ar": [1521, 2666, 3931, 5771, 7237, 8782, 12000, 13842],
    "K":  [419, 3052, 4420, 5878, 7976, 9590, 11343, 14944],
    "Ca": [590, 1145, 4912, 6491, 8153, 10496, 12273, 14206],
    "Sc": [633, 1235, 2389, 7091, 8843, 10679, 13400, 15301],
    "Ti": [658, 1310, 2652, 4175, 9581, 11534, 14400, 15692],
    "V":  [650, 1414, 2830, 4507, 6294, 12363, 14530, 16731],
    "Cr": [653, 1592, 2987, 4743, 6702, 8745, 15456, 17820],
    "Mn": [717, 1509, 3248, 4940, 6990, 9220, 11500, 19000],
    "Fe": [762, 1561, 2957, 5290, 7240, 9560, 12060, 14580],
    "Co": [760, 1648, 3232, 4950, 7670, 9840, 12400, 14760],
    "Ni": [737, 1753, 3395, 5300, 7280, 10400, 12800, 15600],
    "Cu": [745, 1958, 3555, 5536, 7700, 9900, 13400, 16000],
    "Zn": [906, 1733, 2836, 4037, 5530, 6965, 9523, 13000],
    "Ga": [579, 1979, 2963, 6180],
    "Ge": [762, 1537, 3302, 4411, 9020],
    "As": [947, 1798, 2735, 4837, 6043, 12310],
    "Se": [941, 2045, 2974, 4144, 6590, 7880, 14990],
    "Br": [1140, 2104, 3500, 4560, 5760, 8580, 9940, 18600],
    "Kr": [1351, 2350, 3565, 5070, 6240, 7570, 11159, 12600],
    "Rb": [403, 2633, 3900, 5080, 6850, 8140, 9570, 13120],
    "Sr": [550, 1064, 4138, 5500, 6910, 8760, 10243, 12800],
    "Y":  [600, 1180, 1980, 5970, 7430, 8970, 11230, 12700],
    "Zr": [640, 1340, 2210, 3313, 7759, 9500],
    "Nb": [652, 1380, 2421, 3700, 4877, 9847, 12100],
    "Mo": [685, 1558, 2621, 4480, 5255, 6640, 12190, 13900],
    "Tc": [702, 1470, 2850],
    "Ru": [711, 1617, 2747],
    "Rh": [720, 1740, 2997],
    "Pd": [804, 1870, 3177],
    "Ag": [731, 2073, 3361],
    "Cd": [867, 1631, 3616],
    "In": [558, 1821, 2704, 5210],
    "Sn": [709, 1412, 2943, 3931, 7466],
    "Sb": [834, 1592, 2440, 4260, 5400, 6925, 9500],
    "Te": [869, 1790, 2698, 3610, 5668, 6822, 13200],
    "I":  [1008, 1846, 3180],
    "Xe": [1170, 2047, 3099],
    "Cs": [376, 2235, 3400],
    "Ba": [503, 965, 3600],
    "La": [538, 1067, 1850, 4819, 5940],
    "Ce": [534, 1050, 1949, 3547, 5325],
    "Pr": [527, 1015, 2086, 3761, 5551],
    "Nd": [533, 1040, 2130, 3900],
    "Pm": [538, 1050, 2150],
    "Sm": [544, 1070, 2260],
    "Eu": [547, 1085, 2400],
    "Gd": [593, 1160, 1997],
    "Tb": [565, 1112, 2124],
    "Dy": [573, 1130, 2200],
    "Ho": [581, 1140, 2204],
    "Er": [589, 1151, 2194],
    "Tm": [597, 1190, 1900],
    "Yb": [603, 1175, 2417],
    "Lu": [524, 1340, 2022],
    "Hf": [659, 1440, 2250],
    "Ta": [761, 1500],
    "W":  [770, 1700],
    "Re": [760, 1260],
    "Os": [840, 1600],
    "Ir": [880, 1600],
    "Pt": [870, 1791],
    "Au": [890, 1980],
    "Hg": [1007, 1810],
    "Tl": [589, 1971],
    "Pb": [716, 1450, 3082, 4160],
    "Bi": [703, 1610, 2466, 4370],
    "Po": [812,]
}


@ChemMCPManager.register_tool
class GetIonizationEnergy(BaseTool):
    __version__ = "0.1.0"
    name = "GetIonizationEnergy"
    func_name = 'get_ionization_energy'
    description = "Get ionization energy data for an element (first through available ionizations)."
    implementation_description = "Uses NIST-standard ionization energy data in kJ/mol. Returns all available successive ionization energies from IE1 upward."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Ionization Energy", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., Na, Fe)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('element', 'str', 'Element symbol'),
        ('ionization_energies', 'dict', 'IE1, IE2, ... in kJ/mol'),
        ('unit', 'str', 'Unit of measurement'),
    ]
    examples = [
        {'code_input': {'element': 'Na'}, 'text_input': {'element': 'Na'}, 'output': {'element': 'Na', 'ionization_energies': {'IE1': 496, 'IE2': 4563}, 'unit': 'kJ/mol'}},
    ]

    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        sym = data["symbol"]
        if sym not in IONIZATION_ENERGIES:
            return {"element": sym, "ionization_energies": {}, "unit": "kJ/mol", "note": f"Data not available for {sym}"}
        ie_list = IONIZATION_ENERGIES[sym]
        ie_dict = {f"IE{i+1}": val for i, val in enumerate(ie_list)}
        return {
            "element": sym,
            "ionization_energies": ie_dict,
            "unit": "kJ/mol",
            "count": len(ie_list),
        }


if __name__ == "__main__":
    run_mcp_server()
