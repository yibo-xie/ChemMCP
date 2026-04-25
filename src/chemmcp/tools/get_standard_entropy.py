import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GetStandardEntropy(BaseTool):
    """
    查询标准熵 S°（298.15 K，标准状态）。
    内置常见化合物的热力学数据库。
    """
    __version__ = "0.1.0"
    name = "GetStandardEntropy"
    func_name = "get_standard_entropy"
    description = "Query standard absolute entropy (S°) for chemical species at 298 K and 1 bar."
    implementation_description = "Uses a built-in thermodynamic database of standard entropies (S° in J/(mol·K)) for common elements and compounds at standard conditions (298.15 K, 100 kPa/1 bar)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Entropy", "Standard Data", "Thermochemistry"]
    required_envs = []

    code_input_sig = [
        ("species", "str", "N/A", "Chemical species name or formula, e.g., 'H2O', 'CO2', 'NH3', 'Fe2O3'."),
    ]

    text_input_sig = [
        ("species_str", "str", "N/A", "Chemical species name or formula string."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with species, standard_s (J/(mol·K)), unit, state, temperature."),
    ]

    examples = [
        {
            "code_input": {"species": "H2O"},
            "text_input": {"species_str": "H2O"},
            "output": {
                "result": {
                    "species": "H2O",
                    "standard_s_j_per_mol_k": 69.91,
                    "unit": "J/(mol·K)",
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
                    "standard_s_j_per_mol_k": 213.79,
                    "unit": "J/(mol·K)",
                    "state": "gas",
                    "temperature_k": 298.15,
                }
            },
        },
        {
            "code_input": {"species": "Fe2O3"},
            "text_input": {"species_str": "Fe2O3"},
            "output": {
                "result": {
                    "species": "Fe2O3",
                    "standard_s_j_per_mol_k": 87.40,
                    "unit": "J/(mol·K)",
                    "state": "solid",
                    "temperature_k": 298.15,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize standard entropy database."""
        # S° values in J/(mol·K) at 298.15 K, 1 bar
        # Source: NIST / CRC Handbook
        self._db = {
            # ── Elements ──
            "H2(g)":     (130.68,   "gas"),
            "O2(g)":     (205.15,   "gas"),
            "N2(g)":     (191.61,   "gas"),
            "C(graphite)": (5.74,   "solid"),
            "C(s, graphite)": (5.74,"solid"),
            "S(s)":      (31.80,    "solid"),
            "Cl2(g)":    (223.07,   "gas"),
            "Br2(l)":    (152.23,   "liquid"),
            "I2(s)":     (116.14,   "solid"),
            "F2(g)":     (202.79,   "gas"),
            "Na(s)":     (51.21,    "solid"),
            "Fe(s)":     (27.28,    "solid"),
            "Zn(s)":     (41.63,    "solid"),
            "Cu(s)":     (33.15,    "solid"),
            "Ag(s)":     (42.55,    "solid"),
            "Ca(s)":     (41.59,    "solid"),
            "Mg(s)":     (32.68,    "solid"),
            "Al(s)":     (28.33,    "solid"),
            "Si(s)":     (18.83,    "solid"),

            # ── Inorganic compounds ──
            # Water
            "H2O(l)":    (69.91,    "liquid"),
            "H2O(g)":    (188.83,   "gas"),
            "H2O":       (69.91,    "liquid"),

            # Hydrogen compounds
            "HCl(g)":    (186.91,   "gas"),
            "HBr(g)":    (198.70,   "gas"),
            "HI(g)":     (206.59,   "gas"),
            "H2S(g)":    (205.79,   "gas"),
            "HF(g)":     (173.78,   "gas"),
            "NH3(g)":    (192.45,   "gas"),
            "NH3":       (192.45,   "gas"),
            "CH4(g)":    (186.26,   "gas"),
            "CH4":       (186.26,   "gas"),
            "C2H6(g)":   (229.60,   "gas"),
            "C2H4(g)":   (219.32,   "gas"),
            "C2H2(g)":   (200.94,   "gas"),

            # Oxides
            "CO(g)":     (197.67,   "gas"),
            "CO":        (197.67,   "gas"),
            "CO2(g)":    (213.79,   "gas"),
            "CO2":       (213.79,   "gas"),
            "SO2(g)":    (248.21,   "gas"),
            "SO2":       (248.21,   "gas"),
            "SO3(g)":    (256.76,   "gas"),
            "NO(g)":     (210.76,   "gas"),
            "NO":        (210.76,   "gas"),
            "NO2(g)":    (240.06,   "gas"),
            "N2O(g)":    (220.00,   "gas"),
            "N2O4(g)":   (304.38,   "gas"),
            "P4O10(s)":  (228.86,   "solid"),
            "Na2O(s)":   (75.06,    "solid"),
            "MgO(s)":    (26.94,    "solid"),
            "MgO":       (26.94,    "solid"),
            "CaO(s)":    (39.75,    "solid"),
            "CaO":       (39.75,    "solid"),
            "Al2O3(s)":  (50.92,    "solid"),
            "Al2O3":     (50.92,    "solid"),
            "Fe2O3(s)":  (87.40,    "solid"),
            "Fe2O3":     (87.40,    "solid"),
            "Fe3O4(s)":  (146.44,   "solid"),
            "FeO(s)":    (60.75,    "solid"),
            "CuO(s)":    (42.63,    "solid"),
            "Cu2O(s)":   (93.14,    "solid"),
            "ZnO(s)":    (43.64,    "solid"),
            "Ag2O(s)":   (121.30,   "solid"),
            "PbO2(s)":   (67.00,    "solid"),
            "PbO(s)":    (66.50,    "solid"),
            "SiO2(s)":   (41.84,    "solid"),
            "SiO2":      (41.84,    "solid"),
            "Cr2O3(s)":  (81.20,    "solid"),
            "MnO2(s)":   (53.05,    "solid"),
            "MnO(s)":    (59.71,    "solid"),
            "TiO2(s)":   (50.62,    "solid"),
            "BaO(s)":    (70.42,    "solid"),

            # Peroxide
            "H2O2(l)":   (109.60,   "liquid"),
            "H2O2(g)":   (232.99,   "gas"),
            "Na2O2(s)":  (94.82,    "solid"),

            # Acids
            "H2SO4(l)": (156.90,   "liquid"),
            "HNO3(l)":  (155.60,   "liquid"),
            "H3PO4(s)": (110.50,   "solid"),

            # Hydroxides
            "NaOH(s)":  (64.46,    "solid"),
            "KOH(s)":   (79.32,    "solid"),
            "Ca(OH)2(s)": (83.39,  "solid"),
            "Al(OH)3(s)": (88.00,   "solid"),
            "Fe(OH)3(s)": (106.70,  "solid"),
            "Mg(OH)2(s)": (63.18,   "solid"),

            # Halides
            "NaCl(s)":  (72.13,    "solid"),
            "NaCl":     (72.13,    "solid"),
            "KCl(s)":   (82.59,    "solid"),
            "AgCl(s)":  (96.20,    "solid"),
            "AgBr(s)":  (107.11,   "solid"),
            "AgI(s)":   (115.10,   "solid"),
            "CaCl2(s)": (104.60,   "solid"),
            "MgCl2(s)": (89.62,    "solid"),
            "AlCl3(s)": (109.29,   "solid"),
            "FeCl2(s)": (118.00,   "solid"),
            "FeCl3(s)": (142.30,   "solid"),
            "CuCl2(s)": (108.00,   "solid"),
            "ZnCl2(s)": (111.46,   "solid"),

            # Carbonates
            "CaCO3(s)": (92.90,    "solid"),
            "CaCO3":    (92.90,    "solid"),
            "Na2CO3(s)": (134.98,   "solid"),
            "NaHCO3(s)": (101.70,   "solid"),
            "BaCO3(s)": (112.10,   "solid"),
            "MgCO3(s)": (65.69,    "solid"),

            # Sulfates
            "CaSO4(s)": (106.70,   "solid"),
            "BaSO4(s)": (132.16,   "solid"),
            "Na2SO4(s)": (149.58,   "solid"),
            "CuSO4(s)": (113.00,   "solid"),
            "MgSO4(s)": (126.40,   "solid"),
            "ZnSO4(s)": (119.65,   "solid"),

            # Nitrates
            "NaNO3(s)": (116.52,   "solid"),
            "KNO3(s)":  (151.04,   "solid"),
            "Ca(NO3)2(s)": (193.30,  "solid"),
            "NH4NO3(s)": (151.08,   "solid"),

            # Sulfides
            "FeS(s)":    (60.29,    "solid"),
            "ZnS(s)":    (57.70,    "solid"),
            "PbS(s)":    (91.30,    "solid"),
            "CuS(s)":    (66.53,    "solid"),
            "CdS(s)":    (72.40,    "solid"),

            # Ammonium salts
            "NH4Cl(s)": (94.56,    "solid"),
            "(NH4)2SO4(s)": (220.08,  "solid"),

            # Other inorganic
            "PCl3(g)":   (311.67,   "gas"),
            "PCl5(g)":   (352.88,   "gas"),
            "SF6(g)":   (291.54,   "gas"),
            "BF3(g)":   (254.12,   "gas"),
            "CS2(l)":   (151.34,   "liquid"),
            "HCN(g)":   (201.78,   "gas"),

            # Ions (aqueous)
            "H+(aq)":   (0.0,      "aqueous"),
            "OH-(aq)":  (-10.75,   "aqueous"),
            "Na+(aq)":  (58.41,    "aqueous"),
            "K+(aq)":   (102.50,   "aqueous"),
            "Ag+(aq)":  (72.68,    "aqueous"),
            "Ca2+(aq)": (-53.10,   "aqueous"),
            "Mg2+(aq)": (-138.10,  "aqueous"),
            "Fe2+(aq)": (-137.70,  "aqueous"),
            "Fe3+(aq)": (-315.90,  "aqueous"),
            "Cu2+(aq)": (-98.00,   "aqueous"),
            "Zn2+(aq)": (-112.10,  "aqueous"),
            "Al3+(aq)": (-321.70,  "aqueous"),
            "Cl-(aq)":  (56.48,    "aqueous"),
            "Br-(aq)":  (111.30,   "aqueous"),
            "I-(aq)":   (111.30,   "aqueous"),
            "SO4^2-(aq)": (20.10,   "aqueous"),
            "CO3^2-(aq)": (-56.90,   "aqueous"),
            "NO3-(aq)": (146.40,   "aqueous"),
            "NH4+(aq)": (112.84,   "aqueous"),

            # Organic
            "C6H12O6(s, glucose)": (212.13, "solid"),
            "C12H22O11(s, sucrose)": (360.24, "solid"),
            "C2H5OH(l)": (160.70,   "liquid"),
            "CH3OH(l)": (126.70,   "liquid"),
            "C6H6(l)":  (173.26,   "liquid"),  # benzene
            "HCOOH(l)": (129.00,   "liquid"),
            "CH3COOH(l)": (159.83,  "liquid"),
            "C6H12O6(s)": (212.13,   "solid"),
        }

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
            "benzene": "C6H6(l)", "glucose": "C6H12O6(s, glucose)",
            "formic acid": "HCOOH(l)", "acetic acid": "CH3COOH(l)",
        }

    def _run_base(self, species: str) -> dict:
        """Core logic: look up S° for given species."""
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
            f"Species '{species}' not found in the standard entropy database. "
            f"The database contains ~150 common elements and compounds."
        )

    def _run_text(self, species_str: str) -> dict:
        return self._run_base(species_str.strip())

    def _lookup(self, key: str) -> Optional[dict]:
        if key in self._db:
            s_val, state = self._db[key]
            return {
                "species": key,
                "standard_s_j_per_mol_k": s_val,
                "unit": "J/(mol·K)",
                "state": state,
                "temperature_k": 298.15,
                "pressure": "1 bar (standard pressure)",
                "note": f"Standard absolute entropy of {key} ({state}) at 298.15 K.",
            }

        import re
        stripped = re.sub(r'\([a-z]+\)$', '', key)
        if stripped in self._db:
            s_val, state = self._db[stripped]
            return {
                "species": key,
                "standard_s_j_per_mol_k": s_val,
                "unit": "J/(mol·K)",
                "state": state,
                "temperature_k": 298.15,
                "pressure": "1 bar (standard pressure)",
            }

        return None
