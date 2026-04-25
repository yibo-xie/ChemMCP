import logging
import math
import re

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

ATOMIC_MASSES = {
    "n": 1.008665, "p": 1.007276, "H-1": 1.007825, "H-2": 2.014102,
    "He-3": 3.016029, "He-4": 4.002603, "Li-7": 7.016004, "Be-9": 9.012183,
    "C-12": 12.000000, "C-13": 13.003355, "N-14": 14.003074, "O-16": 15.994915,
    "Fe-56": 55.934938, "U-235": 235.043930, "U-238": 238.050788,
}
U_TO_MEV = 931.49410242


@ChemMCPManager.register_tool
class MassDefect(BaseTool):
    __version__ = "0.1.0"
    name = "MassDefect"
    func_name = "calculate_mass_defect"
    description = "Calculate nuclear mass defect and its energy equivalent for a given nuclide."
    implementation_description = "Computes mass defect as constituent mass minus actual mass, converts to energy via E=Delta_m*c^2."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Nuclear Chemistry", "Mass Defect", "Mass-Energy Equivalence"]
    required_envs = []

    code_input_sig = [
        ("nuclide", "str", "N/A", "Nuclide identifier (e.g., 'Fe-56', 'He-4')."),
        ("custom_mass_amu", "float", "0", "Custom atomic mass in u (0=use database)."),
    ]

    text_input_sig = [
        ("nuclide_str", "str", "N/A", "Nuclide string (e.g., 'Fe-56')."),
    ]

    output_sig = [
        ("nuclide", "str", "The input nuclide."),
        ("proton_count", "int", "Number of protons Z."),
        ("neutron_count", "int", "Number of neutrons A-Z."),
        ("constituent_mass_u", "float", "Total mass of separate nucleons in u."),
        ("actual_nuclear_mass_u", "float", "Actual measured nuclear mass in u."),
        ("mass_defect_u", "float", "Mass defect in u."),
        ("mass_defect_kg", "str", "Mass defect in kg."),
        ("energy_equivalent_mev", "float", "Energy equivalent in MeV."),
        ("energy_equivalent_joules", "str", "Energy equivalent in Joules."),
        ("packing_fraction", "float", "Packing fraction = Delta_m/A * 10^4."),
    ]

    examples = [
        {
            "code_input": {"nuclide": "He-4", "custom_mass_amu": 0.0},
            "text_input": {"nuclide_str": "He-4"},
            "output": {
                "nuclide": "He-4",
                "proton_count": 2,
                "neutron_count": 2,
                "constituent_mass_u": 4.032980,
                "actual_nuclear_mass_u": 4.002603,
                "mass_defect_u": 0.030377,
                "mass_defect_kg": "5.044e-29",
                "energy_equivalent_mev": 28.30,
                "energy_equivalent_joules": "4.535e-12",
                "packing_fraction": 75.94,
            }
        },
        {
            "code_input": {"nuclide": "U-235", "custom_mass_amu": 0.0},
            "text_input": {"nuclide_str": "U-235"},
            "output": {
                "nuclide": "U-235",
                "proton_count": 92,
                "neutron_count": 143,
                "constituent_mass_u": 236.972395,
                "actual_nuclear_mass_u": 235.043930,
                "mass_defect_u": 1.928465,
                "mass_defect_kg": "3.203e-27",
                "energy_equivalent_mev": 1796.48,
                "energy_equivalent_joules": "2.879e-10",
                "packing_fraction": 82.07,
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _parse_nuclide(self, s):
        m = re.match(r'^([A-Z][a-z]?)-?(\d+)$', s.strip())
        if not m:
            raise ChemMCPError(f"Cannot parse nuclide '{s}'.")
        return m.group(1), int(m.group(2))

    def _get_z(self, symbol):
        elem_z = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
                   "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,
                   "Fe":26,"Ni":28,"Zn":30,"Sr":38,"Zr":40,"Ba":56,"Pb":82,"U":92}
        z = elem_z.get(symbol)
        if z is None:
            raise ChemMCPError(f"Unknown element: '{symbol}'")
        return z

    def _run_base(self, nuclide, custom_mass_amu=0.0):
        symbol, a = self._parse_nuclide(nuclide)
        z = self._get_z(symbol)
        n_neutrons = a - z
        if custom_mass_amu > 0:
            mass = custom_mass_amu
        elif nuclide in ATOMIC_MASSES:
            mass = ATOMIC_MASSES[nuclide]
        else:
            raise ChemMCPError(f"No mass data for '{nuclide}'.")
        m_h = ATOMIC_MASSES["H-1"]
        m_n_val = ATOMIC_MASSES["n"]
        constituent_mass = z * m_h + n_neutrons * m_n_val
        dm_u = constituent_mass - mass
        dm_kg = dm_u * 1.66053906660e-27
        energy_mev = dm_u * U_TO_MEV
        energy_j = energy_mev * 1.602176634e-13
        packing_frac = (dm_u / a) * 1e4 if a > 0 else 0
        return {
            "nuclide": f"{symbol}-{a}",
            "proton_count": z,
            "neutron_count": n_neutrons,
            "constituent_mass_u": round(constituent_mass, 6),
            "actual_nuclear_mass_u": round(mass, 6),
            "mass_defect_u": round(dm_u, 6),
            "mass_defect_kg": f"{dm_kg:.3e}",
            "energy_equivalent_mev": round(energy_mev, 4),
            "energy_equivalent_joules": f"{energy_j:.3e}",
            "packing_fraction": round(packing_frac, 4),
        }

    def _run_text(self, nuclide_str):
        return self._run_base(nuclide_str.strip())
