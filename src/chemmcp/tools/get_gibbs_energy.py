import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetGibbsEnergy(BaseTool):
    """
    Query standard Gibbs free energy of formation (dGf) at 298.15 K.
    Built-in database of ~150 common species from NIST/CRC.
    """
    __version__ = "0.1.0"
    name = "GetGibbsEnergy"
    func_name = "get_gibbs_energy"
    description = "Query standard Gibbs free energy of formation (dGf) for chemical species at 298 K and 1 bar."
    implementation_description = "Built-in thermodynamic database of dGf values (kJ/mol) from NIST/CRC Handbook at standard conditions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Gibbs Energy", "Thermodynamics", "Standard Data", "Thermochemistry"]
    required_envs = []

    code_input_sig = [
        ("species", "str", "N/A", "Chemical species name or formula, e.g., 'H2O', 'CO2', 'NH3', 'Fe2O3'."),
    ]

    text_input_sig = [
        ("species_str", "str", "N/A", "Chemical species name or formula string."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with species, delta_gf (kJ/mol), unit, state, temperature."),
    ]

    examples = [
        {
            "code_input": {"species": "H2O"},
            "text_input": {"species_str": "H2O"},
            "output": {
                "result": {
                    "species": "H2O",
                    "delta_gf_kj_per_mol": -237.14,
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
                    "delta_gf_kj_per_mol": -394.36,
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
                    "delta_gf_kj_per_mol": -16.40,
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
        # Elements - dGf = 0
        elems0 = ["H2(g)", "O2(g)", "N2(g)", "C(graphite)", "S(s)", "Cl2(g)",
                  "Br2(l)", "I2(s)", "F2(g)", "Na(s)", "Fe(s)", "Zn(s)",
                  "Cu(s)", "Ag(s)", "Ca(s)", "Mg(s)", "Al(s)", "Si(s)"]
        for e in elems0:
            self._db[e] = (0.0, "element")

        # Compounds: (dGf kJ/mol, state)
        compounds = [
            # Water
            ("H2O(l)", -237.14, "liquid"), ("H2O(g)", -228.57, "gas"), ("H2O", -237.14, "liquid"),
            # Hydrogen compounds
            ("HCl(g)", -95.30, "gas"), ("HBr(g)", -53.45, "gas"), ("HI(g)", 1.70, "gas"),
            ("H2S(g)", -33.56, "gas"), ("HF(g)", -275.46, "gas"),
            ("NH3(g)", -16.40, "gas"), ("NH3", -16.40, "gas"),
            ("CH4(g)", -50.75, "gas"), ("CH4", -50.75, "gas"),
            ("C2H6(g)", -32.82, "gas"), ("C2H4(g)", 68.15, "gas"), ("C2H2(g)", 209.20, "gas"),
            # Oxides
            ("CO(g)", -137.16, "gas"), ("CO", -137.16, "gas"),
            ("CO2(g)", -394.36, "gas"), ("CO2", -394.36, "gas"),
            ("SO2(g)", -300.19, "gas"), ("SO2", -300.19, "gas"),
            ("SO3(g)", -371.06, "gas"),
            ("NO(g)", 86.55, "gas"), ("NO", 86.55, "gas"),
            ("NO2(g)", 51.26, "gas"), ("N2O(g)", 104.17, "gas"), ("N2O4(g)", 97.82, "gas"),
            ("P4O10(s)", -2769.94, "solid"),
            ("Na2O(s)", -378.00, "solid"),
            # Metal oxides
            ("MgO(s)", -569.33, "solid"), ("MgO", -569.33, "solid"),
            ("CaO(s)", -603.54, "solid"), ("CaO", -603.54, "solid"),
            ("Al2O3(s)", -1582.27, "solid"), ("Al2O3", -1582.27, "solid"),
            ("Fe2O3(s)", -742.24, "solid"), ("Fe2O3", -742.24, "solid"),
            ("Fe3O4(s)", -1015.38, "solid"), ("FeO(s)", -255.20, "solid"),
            ("CuO(s)", -129.71, "solid"), ("Cu2O(s)", -146.00, "solid"),
            ("ZnO(s)", -320.52, "solid"), ("Ag2O(s)", -11.21, "solid"),
            ("PbO(s)", -187.89, "solid"), ("PbO2(s)", -218.99, "solid"),
            ("SiO2(s)", -856.64, "solid"), ("SiO2", -856.64, "solid"),
            ("Cr2O3(s)", -1058.09, "solid"), ("MnO2(s)", -465.14, "solid"),
            ("MnO(s)", -361.50, "solid"), ("TiO2(s)", -888.60, "solid"),
            ("BaO(s)", -520.30, "solid"),
            # Peroxides
            ("H2O2(l)", -120.35, "liquid"), ("H2O2(g)", -105.48, "gas"),
            ("Na2O2(s)", -447.69, "solid"),
            # Acids
            ("H2SO4(l)", -689.90, "liquid"), ("HNO3(l)", -80.71, "liquid"),
            ("H3PO4(s)", -1124.28, "solid"),
            # Bases
            ("NaOH(s)", -379.53, "solid"), ("KOH(s)", -379.08, "solid"),
            ("Ca(OH)2(s)", -898.49, "solid"),
            ("Al(OH)3(s)", -1147.00, "solid"), ("Fe(OH)3(s)", -696.50, "solid"),
            ("Mg(OH)2(s)", -833.58, "solid"),
            # Halides
            ("NaCl(s)", -384.14, "solid"), ("NaCl", -384.14, "solid"),
            ("KCl(s)", -408.56, "solid"), ("AgCl(s)", -109.79, "solid"),
            ("AgBr(s)", -96.90, "solid"), ("AgI(s)", -66.19, "solid"),
            ("CaCl2(s)", -748.10, "solid"), ("MgCl2(s)", -591.79, "solid"),
            ("AlCl3(s)", -628.80, "solid"), ("FeCl2(s)", -302.30, "solid"),
            ("FeCl3(s)", -333.99, "solid"), ("CuCl2(s)", -179.90, "solid"),
            ("ZnCl2(s)", -369.43, "solid"),
            # Carbonates
            ("CaCO3(s)", -1128.76, "solid"), ("CaCO3", -1128.76, "solid"),
            ("Na2CO3(s)", -1044.44, "solid"), ("NaHCO3(s)", -851.87, "solid"),
            ("BaCO3(s)", -1134.41, "solid"), ("MgCO3(s)", -1012.11, "solid"),
            # Sulfates
            ("CaSO4(s)", -1321.74, "solid"), ("BaSO4(s)", -1362.18, "solid"),
            ("Na2SO4(s)", -1270.12, "solid"), ("CuSO4(s)", -661.86, "solid"),
            ("MgSO4(s)", -1170.66, "solid"), ("ZnSO4(s)", -891.58, "solid"),
            # Nitrates
            ("NaNO3(s)", -367.04, "solid"), ("KNO3(s)", -394.93, "solid"),
            ("Ca(NO3)2(s)", -743.34, "solid"), ("NH4NO3(s)", -189.47, "solid"),
            # Sulfides
            ("FeS(s)", -100.42, "solid"), ("ZnS(s)", -201.29, "solid"),
            ("PbS(s)", -98.73, "solid"), ("CuS(s)", -53.59, "solid"),
            ("CdS(s)", -156.50, "solid"),
            # Ammonium salts
            ("NH4Cl(s)", -202.95, "solid"), ("(NH4)2SO4(s)", -902.28, "solid"),
            # Other inorganic
            ("PCl3(g)", -267.77, "gas"), ("PCl5(g)", -289.23, "gas"),
            ("SF6(g)", -1105.14, "gas"), ("BF3(g)", -1120.33, "gas"),
            ("CS2(l)", 65.27, "liquid"), ("HCN(g)", 124.60, "gas"),
            # Aqueous ions
            ("H+(aq)", 0.0, "aqueous"), ("OH-(aq)", -157.24, "aqueous"),
            ("Na+(aq)", -261.88, "aqueous"), ("K+(aq)", -283.27, "aqueous"),
            ("Ag+(aq)", 77.11, "aqueous"), ("Ca2+(aq)", -553.58, "aqueous"),
            ("Mg2+(aq)", -456.01, "aqueous"),
            ("Fe2+(aq)", -78.90, "aqueous"), ("Fe3+(aq)", -4.65, "aqueous"),
            ("Cu2+(aq)", 65.49, "aqueous"), ("Zn2+(aq)", -147.19, "aqueous"),
            ("Al3+(aq)", -485.00, "aqueous"),
            ("Cl-(aq)", -131.22, "aqueous"), ("Br-(aq)", -103.96, "aqueous"),
            ("I-(aq)", -51.59, "aqueous"),
            ("SO4^2-(aq)", -744.53, "aqueous"), ("CO3^2-(aq)", -527.81, "aqueous"),
            ("NO3-(aq)", -111.25, "aqueous"), ("NH4+(aq)", -79.31, "aqueous"),
            # Organic
            ("C6H12O6(s, glucose)", -910.44, "solid"),
            ("C12H22O11(s, sucrose)", -1551.78, "solid"),
            ("C2H5OH(l)", -174.78, "liquid"), ("CH3OH(l)", -166.27, "liquid"),
            ("C6H6(l)", 124.45, "liquid"),
            ("HCOOH(l)", -361.40, "liquid"), ("CH3COOH(l)", -389.85, "liquid"),
            ("C6H12O6(s)", -910.44, "solid"),
        ]
        for entry in compounds:
            self._db[entry[0]] = (entry[1], entry[2])

        self._aliases = {
            "water": "H2O(l)", "water vapor": "H2O(g)", "steam": "H2O(g)",
            "ammonia": "NH3(g)", "methane": "CH4(g)", "ethane": "C2H6(g)",
            "ethene": "C2H4(g)", "ethyne": "C2H2(g)",
            "carbon monoxide": "CO(g)", "carbon dioxide": "CO2(g)",
            "sulfur dioxide": "SO2(g)", "sulfur trioxide": "SO3(g)",
            "nitric oxide": "NO(g)", "nitrogen dioxide": "NO2(g)",
            "table salt": "NaCl(s)", "salt": "NaCl(s)",
            "lime": "CaO(s)", "quicklime": "CaO(s)",
            "iron oxide": "Fe2O3(s)", "rust": "Fe2O3(s)",
            "quartz": "SiO2(s)", "sand": "SiO2(s)",
            "corundum": "Al2O3(s)", "alumina": "Al2O3(s)",
            "calcite": "CaCO3(s)", "limestone": "CaCO3(s)",
            "ethanol": "C2H5OH(l)", "methanol": "CH3OH(l)",
            "benzene": "C6H6(l)",
            "formic acid": "HCOOH(l)", "acetic acid": "CH3COOH(l)",
            "glucose": "C6H12O6(s, glucose)",
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
            f"Species '{species}' not found in the standard Gibbs free energy database."
        )

    def _run_text(self, species_str: str) -> dict:
        return self._run_base(species_str.strip())

    def _lookup(self, key: str):
        if key in self._db:
            dg, state = self._db[key]
            return {
                "species": key,
                "delta_gf_kj_per_mol": dg,
                "unit": "kJ/mol",
                "state": state,
                "temperature_k": 298.15,
                "pressure": "1 bar (standard pressure)",
                "note": f"Standard Gibbs free energy of formation of {key} ({state}) at 298.15 K.",
            }

        import re
        stripped = re.sub(r'\([a-z]+\)$', '', key)
        if stripped in self._db:
            dg, state = self._db[stripped]
            return {
                "species": key,
                "delta_gf_kj_per_mol": dg,
                "unit": "kJ/mol",
                "state": state,
                "temperature_k": 298.15,
                "pressure": "1 bar (standard pressure)",
            }
        return None
