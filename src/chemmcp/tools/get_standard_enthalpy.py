import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetStandardEnthalpy(BaseTool):
    """
    Query standard enthalpy of formation (dHf) at 298.15 K.
    Built-in database of ~150 common species from NIST/CRC.
    """
    __version__ = "0.1.0"
    name = "GetStandardEnthalpy"
    func_name = "get_standard_enthalpy"
    description = "Query standard enthalpy of formation (dHf) for chemical species at 298 K and 1 bar."
    implementation_description = "Built-in thermodynamic database of dHf values (kJ/mol) from NIST/CRC Handbook at standard conditions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermochemistry", "Enthalpy", "Standard Data", "Thermodynamics"]
    required_envs = []

    code_input_sig = [
        ("species", "str", "N/A", "Chemical species name or formula, e.g., 'H2O', 'CO2', 'NH3', 'Fe2O3'."),
    ]

    text_input_sig = [
        ("species_str", "str", "N/A", "Chemical species name or formula string."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with species, delta_hf (kJ/mol), unit, state, temperature."),
    ]

    examples = [
        {
            "code_input": {"species": "H2O"},
            "text_input": {"species_str": "H2O"},
            "output": {
                "result": {
                    "species": "H2O",
                    "delta_hf_kj_per_mol": -285.83,
                    "unit": "kJ/mol",
                    "state": "liquid",
                    "temperature_k": 298.15,
                }
            },
        },
        {
            "code_input": {"species": "CO2"},
            "text_input": {"species_str": "CO2"},
            "output": {
                "result": {
                    "species": "CO2",
                    "delta_hf_kj_per_mol": -393.51,
                    "unit": "kJ/mol",
                    "state": "gas",
                    "temperature_k": 298.15,
                }
            },
        },
        {
            "code_input": {"species": "NH3"},
            "text_input": {"species_str": "NH3"},
            "output": {
                "result": {
                    "species": "NH3",
                    "delta_hf_kj_per_mol": -45.94,
                    "unit": "kJ/mol",
                    "state": "gas",
                    "temperature_k": 298.15,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._db = {}
        # Elements - dHf = 0
        elems0 = ["H2(g)", "O2(g)", "N2(g)", "C(graphite)", "S(s)", "Cl2(g)",
                  "Br2(l)", "I2(s)", "F2(g)", "Na(s)", "Fe(s)", "Zn(s)",
                  "Cu(s)", "Ag(s)", "Ca(s)", "Mg(s)", "Al(s)", "Si(s)"]
        for e in elems0:
            self._db[e] = (0.0, "element")

        # Inorganic compounds: (dHf kJ/mol, state)
        compounds = [
            # Water
            ("H2O(l)", -285.83, "liquid"), ("H2O(g)", -241.82, "gas"), ("H2O", -285.83, "liquid"),
            # Hydrogen compounds
            ("HCl(g)", -92.31, "gas"), ("HBr(g)", -36.40, "gas"), ("HI(g)", 26.48, "gas"),
            ("H2S(g)", -20.63, "gas"), ("HF(g)", -273.30, "gas"),
            ("NH3(g)", -45.94, "gas"), ("NH3", -45.94, "gas"),
            ("CH4(g)", -74.81, "gas"), ("CH4", -74.81, "gas"),
            ("C2H6(g)", -84.68, "gas"), ("C2H4(g)", 52.47, "gas"), ("C2H2(g)", 227.39, "gas"),
            # Oxides
            ("CO(g)", -110.53, "gas"), ("CO", -110.53, "gas"),
            ("CO2(g)", -393.51, "gas"), ("CO2", -393.51, "gas"),
            ("SO2(g)", -296.84, "gas"), ("SO2", -296.84, "gas"),
            ("SO3(g)", -395.72, "gas"),
            ("NO(g)", 90.25, "gas"), ("NO", 90.25, "gas"),
            ("NO2(g)", 33.18, "gas"), ("N2O(g)", 82.05, "gas"), ("N2O4(g)", 9.16, "gas"),
            ("P4O10(s)", -2984.03, "solid"),
            ("Na2O(s)", -414.22, "solid"),
            # Metal oxides
            ("MgO(s)", -601.60, "solid"), ("MgO", -601.60, "solid"),
            ("CaO(s)", -635.09, "solid"), ("CaO", -635.09, "solid"),
            ("Al2O3(s)", -1675.70, "solid"), ("Al2O3", -1675.70, "solid"),
            ("Fe2O3(s)", -824.20, "solid"), ("Fe2O3", -824.20, "solid"),
            ("Fe3O4(s)", -1118.38, "solid"), ("FeO(s)", -272.00, "solid"),
            ("CuO(s)", -157.28, "solid"), ("Cu2O(s)", -168.60, "solid"),
            ("ZnO(s)", -348.28, "solid"), ("Ag2O(s)", -31.05, "solid"),
            ("PbO2(s)", -277.40, "solid"), ("PbO(s)", -217.32, "solid"),
            ("SiO2(s)", -910.86, "solid"), ("SiO2", -910.86, "solid"),
            ("Cr2O3(s)", -1139.70, "solid"), ("MnO2(s)", -520.03, "solid"),
            ("MnO(s)", -385.22, "solid"), ("TiO2(s)", -944.00, "solid"),
            ("BaO(s)", -548.10, "solid"),
            # Peroxides
            ("H2O2(l)", -187.78, "liquid"), ("H2O2(g)", -136.11, "gas"),
            ("Na2O2(s)", -510.87, "solid"),
            # Acids
            ("H2SO4(l)", -814.00, "liquid"), ("HNO3(l)", -174.10, "liquid"),
            ("H3PO4(s)", -1279.01, "solid"),
            # Hydroxides
            ("NaOH(s)", -425.61, "solid"), ("KOH(s)", -424.72, "solid"),
            ("Ca(OH)2(s)", -986.09, "solid"), ("Ba(OH)2(s)", -944.70, "solid"),
            ("Al(OH)3(s)", -1277.00, "solid"), ("Fe(OH)3(s)", -823.00, "solid"),
            ("Fe(OH)2(s)", -569.00, "solid"), ("Mg(OH)2(s)", -924.54, "solid"),
            # Halides
            ("NaCl(s)", -411.15, "solid"), ("NaCl", -411.15, "solid"),
            ("KCl(s)", -436.75, "solid"), ("AgCl(s)", -127.07, "solid"),
            ("AgBr(s)", -100.37, "solid"), ("AgI(s)", -61.84, "solid"),
            ("CaCl2(s)", -795.80, "solid"), ("MgCl2(s)", -641.32, "solid"),
            ("AlCl3(s)", -704.20, "solid"), ("FeCl2(s)", -342.67, "solid"),
            ("FeCl3(s)", -399.49, "solid"), ("CuCl2(s)", -220.10, "solid"),
            ("ZnCl2(s)", -415.05, "solid"), ("PbCl2(s)", -359.41, "solid"),
            # Carbonates
            ("CaCO3(s)", -1206.92, "solid"), ("CaCO3", -1206.92, "solid"),
            ("Na2CO3(s)", -1130.68, "solid"), ("NaHCO3(s)", -950.81, "solid"),
            ("BaCO3(s)", -1216.29, "solid"), ("MgCO3(s)", -1095.79, "solid"),
            ("KHCO3(s)", -963.19, "solid"),
            # Sulfates
            ("CaSO4(s)", -1434.52, "solid"), ("BaSO4(s)", -1473.19, "solid"),
            ("Na2SO4(s)", -1387.08, "solid"), ("K2SO4(s)", -1437.79, "solid"),
            ("CuSO4(s)", -771.36, "solid"), ("FeSO4(s)", -828.40, "solid"),
            ("Al2(SO4)3(s)", -3440.69, "solid"), ("MgSO4(s)", -1284.91, "solid"),
            ("ZnSO4(s)", -982.82, "solid"),
            # Nitrates
            ("NaNO3(s)", -467.85, "solid"), ("KNO3(s)", -494.63, "solid"),
            ("Ca(NO3)2(s)", -938.21, "solid"), ("NH4NO3(s)", -365.56, "solid"),
            ("Ba(NO3)2(s)", -992.07, "solid"),
            # Phosphates
            ("Ca3(PO4)2(s)", -4120.82, "solid"),
            # Sulfides
            ("FeS(s)", -100.02, "solid"), ("ZnS(s)", -205.98, "solid"),
            ("PbS(s)", -98.28, "solid"), ("CuS(s)", -53.06, "solid"),
            ("CdS(s)", -161.90, "solid"), ("HgS(s)", -58.16, "solid"),
            ("Na2S(s)", -364.79, "solid"), ("MnS(s)", -226.20, "solid"),
            # Ammonium salts
            ("NH4Cl(s)", -314.43, "solid"), ("(NH4)2SO4(s)", -1180.85, "solid"),
            ("NH4HCO3(s)", -849.26, "solid"),
            # Other inorganic
            ("SiCl4(g)", -657.01, "gas"), ("PCl3(g)", -287.02, "gas"),
            ("PCl5(g)", -374.89, "gas"), ("SF6(g)", -1209.00, "gas"),
            ("BF3(g)", -1135.61, "gas"), ("CS2(l)", 89.70, "liquid"),
            ("HCN(g)", 135.14, "gas"), ("PH3(g)", 5.44, "gas"),
            # Aqueous ions
            ("H+(aq)", 0.0, "aqueous"), ("OH-(aq)", -230.02, "aqueous"),
            ("Na+(aq)", -240.12, "aqueous"), ("K+(aq)", -252.14, "aqueous"),
            ("Ag+(aq)", 105.58, "aqueous"), ("Ca2+(aq)", -542.96, "aqueous"),
            ("Mg2+(aq)", -466.85, "aqueous"), ("Ba2+(aq)", -537.64, "aqueous"),
            ("Fe2+(aq)", -89.12, "aqueous"), ("Fe3+(aq)", -48.53, "aqueous"),
            ("Cu2+(aq)", 64.77, "aqueous"), ("Zn2+(aq)", -153.89, "aqueous"),
            ("Al3+(aq)", -531.00, "aqueous"),
            ("Cl-(aq)", -167.08, "aqueous"), ("Br-(aq)", -121.41, "aqueous"),
            ("I-(aq)", -55.19, "aqueous"),
            ("SO4^2-(aq)", -909.27, "aqueous"), ("CO3^2-(aq)", -677.14, "aqueous"),
            ("NO3-(aq)", -205.00, "aqueous"), ("NH4+(aq)", -132.51, "aqueous"),
            # Organic
            ("C6H12O6(s, glucose)", -1273.02, "solid"),
            ("C12H22O11(s, sucrose)", -2226.09, "solid"),
            ("C2H5OH(l)", -277.69, "liquid"), ("CH3OH(l)", -238.66, "liquid"),
            ("HCHO(g)", -108.57, "gas"), ("HCOOH(l)", -424.72, "liquid"),
            ("C6H6(l)", 49.04, "liquid"), ("C6H12O6(s)", -1273.02, "solid"),
            # Other
            ("O3(g)", 142.67, "gas"),
        ]
        for entry in compounds:
            self._db[entry[0]] = (entry[1], entry[2])

        # Aliases
        self._aliases = {
            "water": "H2O(l)", "water vapor": "H2O(g)", "steam": "H2O(g)",
            "ammonia": "NH3(g)", "methane": "CH4(g)", "ethane": "C2H6(g)",
            "ethene": "C2H4(g)", "ethyne": "C2H2(g)",
            "carbon monoxide": "CO(g)", "carbon dioxide": "CO2(g)",
            "sulfur dioxide": "SO2(g)", "sulfur trioxide": "SO3(g)",
            "nitric oxide": "NO(g)", "nitrogen dioxide": "NO2(g)",
            "hydrogen chloride": "HCl(g)", "hydrogen sulfide": "H2S(g)",
            "hydrogen fluoride": "HF(g)",
            "table salt": "NaCl(s)", "salt": "NaCl(s)",
            "lime": "CaO(s)", "quicklime": "CaO(s)",
            "iron oxide": "Fe2O3(s)", "rust": "Fe2O3(s)",
            "quartz": "SiO2(s)", "sand": "SiO2(s)",
            "corundum": "Al2O3(s)", "alumina": "Al2O3(s)",
            "calcite": "CaCO3(s)", "limestone": "CaCO3(s)",
            "gypsum": "CaSO4(s)",
            "ethanol": "C2H5OH(l)", "methanol": "CH3OH(l)",
            "benzene": "C6H6(l)",
            "formaldehyde": "HCHO(g)", "formic acid": "HCOOH(l)",
            "acetic acid": "CH3COOH(l)",
            "glucose": "C6H12O6(s, glucose)", "sucrose": "C12H22O11(s, sucrose)",
            "hydrogen peroxide": "H2O2(l)", "ozone": "O3(g)",
        }

    def _run_base(self, species: str) -> dict:
        species = species.strip()
        if not species:
            raise ChemMCPError("Species name cannot be empty.")

        result = self._lookup(species)
        if result:
            return result

        lower_key = species.lower().strip()
        if lower_key in self._aliases:
            canonical = self._aliases[lower_key]
            result = self._lookup(canonical)
            if result:
                result["query"] = species
                result["canonical_name"] = canonical
                return result

        raise ChemMCPError(
            f"Species '{species}' not found in the standard enthalpy database. "
            f"The database contains ~150 common elements and compounds."
        )

    def _run_text(self, species_str: str) -> dict:
        return self._run_base(species_str.strip())

    def _lookup(self, key: str):
        if key in self._db:
            dh, state = self._db[key]
            return {
                "species": key,
                "delta_hf_kj_per_mol": dh,
                "unit": "kJ/mol",
                "state": state,
                "temperature_k": 298.15,
                "pressure": "1 bar (standard pressure)",
                "note": f"Standard enthalpy of formation of {key} ({state}) at 298.15 K.",
            }

        import re
        stripped = re.sub(r'\([a-z]+\)$', '', key)
        if stripped in self._db:
            dh, state = self._db[stripped]
            return {
                "species": key,
                "delta_hf_kj_per_mol": dh,
                "unit": "kJ/mol",
                "state": state,
                "temperature_k": 298.15,
                "pressure": "1 bar (standard pressure)",
            }

        return None
