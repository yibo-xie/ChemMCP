import logging
import math
import re

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

C_LIGHT = 2.99792458e8
AMU_TO_KG = 1.66053906660e-27
EV_TO_JOULE = 1.602176634e-19
MEV_TO_JOULE = EV_TO_JOULE * 1e6

ATOMIC_MASSES = {
    "n": 1.008665, "p": 1.007276, "H-1": 1.007825, "H-2": 2.014102,
    "He-3": 3.016029, "He-4": 4.002603, "Li-6": 6.015123, "Li-7": 7.016004,
    "Be-9": 9.012183, "C-12": 12.000000, "C-13": 13.003355, "N-14": 14.003074,
    "O-16": 15.994915, "Fe-56": 55.934938, "U-235": 235.043930, "U-238": 238.050788,
}
U_TO_MEV = 931.49410242


@ChemMCPManager.register_tool
class BindingEnergy(BaseTool):
    __version__ = "0.1.0"
    name = "BindingEnergy"
    func_name = "calculate_binding_energy"
    description = "Calculate nuclear binding energy and binding energy per nucleon for a given nuclide."
    implementation_description = "Uses E=Delta_m*c^2 with atomic mass data."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Nuclear Chemistry", "Binding Energy", "Mass-Energy Equivalence"]
    required_envs = []

    code_input_sig = [
        ("nuclide", "str", "N/A", "Nuclide identifier (e.g., 'Fe-56', 'He-4')."),
        ("custom_mass_amu", "float", "0", "Custom atomic mass in u (0=use database)."),
        ("output_unit", "str", "MeV", "Output unit: 'MeV', 'J', or 'eV'."),
    ]

    text_input_sig = [
        ("nuclide_str", "str", "N/A", "Nuclide string, optionally with unit: 'Fe-56 MeV'."),
    ]

    output_sig = [
        ("nuclide", "str", "The input nuclide."),
        ("mass_number", "int", "Mass number A."),
        ("atomic_number", "int", "Atomic number Z."),
        ("neutron_number", "int", "Neutron number N=A-Z."),
        ("nuclear_mass_u", "float", "Nuclear mass in u."),
        ("mass_defect_u", "float", "Mass defect in u."),
        ("binding_energy_mev", "float", "Total binding energy in MeV."),
        ("binding_energy_per_nucleon_mev", "float", "Binding energy per nucleon in MeV."),
        ("binding_energy_joules", "float", "Total binding energy in Joules."),
        ("unit", "str", "Unit used for primary output."),
    ]

    examples = [
        {
            "code_input": {"nuclide": "Fe-56", "custom_mass_amu": 0.0, "output_unit": "MeV"},
            "text_input": {"nuclide_str": "Fe-56"},
            "output": {
                "nuclide": "Fe-56",
                "mass_number": 56,
                "atomic_number": 26,
                "neutron_number": 30,
                "nuclear_mass_u": 55.934938,
                "mass_defect_u": 0.528462,
                "binding_energy_mev": 492.27,
                "binding_energy_per_nucleon_mev": 8.79,
                "binding_energy_joules": 7.886e-11,
                "unit": "MeV",
            }
        },
        {
            "code_input": {"nuclide": "He-4", "custom_mass_amu": 0.0, "output_unit": "MeV"},
            "text_input": {"nuclide_str": "He-4"},
            "output": {
                "nuclide": "He-4",
                "mass_number": 4,
                "atomic_number": 2,
                "neutron_number": 2,
                "nuclear_mass_u": 4.002603,
                "mass_defect_u": 0.030377,
                "binding_energy_mev": 28.30,
                "binding_energy_per_nucleon_mev": 7.07,
                "binding_energy_joules": 4.535e-12,
                "unit": "MeV",
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
        elem_z = {
            "H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
            "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,
            "Fe":26,"Ni":28,"Zn":30,"Sr":38,"Zr":40,"Ba":56,"Pb":82,"U":92,
        }
        z = elem_z.get(symbol)
        if z is None:
            raise ChemMCPError(f"Unknown element: '{symbol}'")
        return z

    def _run_base(self, nuclide, custom_mass_amu=0.0, output_unit="MeV"):
        symbol, a = self._parse_nuclide(nuclide)
        z = self._get_z(symbol)
        n = a - z
        if custom_mass_amu > 0:
            mass = custom_mass_amu
        elif nuclide in ATOMIC_MASSES:
            mass = ATOMIC_MASSES[nuclide]
        else:
            raise ChemMCPError(f"No mass data for '{nuclide}'.")
        m_h = ATOMIC_MASSES["H-1"]
        m_n = ATOMIC_MASSES["n"]
        mass_constituents = z * m_h + n * m_n
        dm = mass_constituents - mass
        eb_mev = dm * U_TO_MEV
        eb_per_nucleon = eb_mev / a if a > 0 else 0
        eb_joules = eb_mev * MEV_TO_JOULE
        return {
            "nuclide": f"{symbol}-{a}",
            "mass_number": a,
            "atomic_number": z,
            "neutron_number": n,
            "nuclear_mass_u": round(mass, 6),
            "mass_defect_u": round(dm, 6),
            "binding_energy_mev": round(eb_mev, 4),
            "binding_energy_per_nucleon_mev": round(eb_per_nucleon, 4),
            "binding_energy_joules": round(eb_joules, 20),
            "unit": output_unit,
        }

    def _run_text(self, nuclide_str):
        parts = nuclide_str.strip().split()
        nuclide = parts[0]
        unit = parts[1] if len(parts) > 1 else "MeV"
        return self._run_base(nuclide, output_unit=unit)
