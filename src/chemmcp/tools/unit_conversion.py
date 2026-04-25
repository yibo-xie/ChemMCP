import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

UNIT_CATEGORIES = {
    "length": {"m":(1,"m"),"cm":(0.01,"m"),"mm":(0.001,"m"),"km":(1000,"m"),"nm":(1e-9,"m"),"pm":(1e-12,"m"),"angstrom":(1e-10,"m"),"in":(0.0254,"m"),"ft":(0.3048,"m")},
    "mass": {"kg":(1,"kg"),"g":(0.001,"kg"),"mg":(1e-6,"kg"),"ug":(1e-9,"kg"),"amu":(1.66053906660e-27,"kg"),"lb":(0.45359237,"kg")},
    "temperature": {"K":(None,"K"),"C":(None,"\u00b0C"),"F":(None,"\u00b0F")},
    "volume": {"L":(0.001,"m\u00b3"),"mL":(1e-6,"m\u00b3"),"m3":(1,"m\u00b3"),"cm3":(1e-6,"m\u00b3")},
    "pressure": {"Pa":(1,"Pa"),"kPa":(1000,"Pa"),"MPa":(1e6,"Pa"),"atm":(101325,"Pa"),"bar":(1e5,"Pa"),"torr":(133.322368421,"Pa"),"mmHg":(133.322368421,"Pa"),"psi":(6894.757293168,"Pa")},
    "energy": {"J":(1,"J"),"kJ":(1000,"J"),"cal":(4.184,"J"),"kcal":(4184,"J"),"eV":(1.602176634e-19,"J"),"MeV":(1.602176634e-13,"J"),"kWh":(3.6e6,"J")},
    "concentration": {"M":(1,"mol/L"),"mM":(0.001,"mol/L"),"uM":(1e-6,"mol/L"),"nM":(1e-9,"mol/L")},
    "time": {"s":(1,"s"),"min":(60,"s"),"h":(3600,"s"),"day":(86400,"s"),"year":(31557600,"s")},
}

def _find_unit(unit_str):
    for cat, units in UNIT_CATEGORIES.items():
        if unit_str in units:
            return cat, units
    return None, None


@ChemMCPManager.register_tool
class UnitConversion(BaseTool):
    __version__ = "0.1.0"
    name = "UnitConversion"
    func_name = "convert_unit"
    description = "Convert between common chemistry units: length, mass, temperature, volume, pressure, energy, concentration, time."
    implementation_description = "Uses built-in conversion factors for chemistry-relevant units with special handling for temperature."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Unit Conversion", "Chemistry", "SI Units"]
    required_envs = []

    code_input_sig = [
        ("value", "float", "N/A", "Numerical value to convert."),
        ("from_unit", "str", "N/A", "Source unit (e.g., 'atm', 'cal', 'eV', 'nm', 'L')."),
        ("to_unit", "str", "N/A", "Target unit (e.g., 'Pa', 'J', 'MeV', 'pm', 'mL')."),
    ]
    text_input_sig = [
        ("conversion_str", "str", "N/A", "Space-separated: value from_unit to_unit (e.g., '1 atm Pa' or '25 C K')."),
    ]
    output_sig = [
        ("input_value", "float", "Original input value."),
        ("input_unit", "str", "Original unit."),
        ("output_value", "float", "Converted result value."),
        ("output_unit", "str", "Target unit."),
        ("category", "str", "Unit category (length, mass, etc.)."),
        ("formula", "str", "Conversion formula or method used."),
    ]
    examples = [
        {
            "code_input": {"value": 1.0, "from_unit": "atm", "to_unit": "Pa"},
            "text_input": {"conversion_str": "1 atm Pa"},
            "output": {"input_value": 1.0, "input_unit": "atm", "output_value": 101325.0, "output_unit": "Pa", "category": "pressure", "formula": "multiply by 101325"},
        },
        {
            "code_input": {"value": 25.0, "from_unit": "C", "to_unit": "K"},
            "text_input": {"conversion_str": "25 C K"},
            "output": {"input_value": 25.0, "input_unit": "\u00b0C", "output_value": 298.15, "output_unit": "K", "category": "temperature", "formula": "K = \u00b0C + 273.15"},
        },
        {
            "code_input": {"value": 1.0, "from_unit": "eV", "to_unit": "kJ/mol"},
            "text_input": {"conversion_str": "1 eV kJ/mol"},
            "output": {"input_value": 1.0, "input_unit": "eV", "output_value": 96.485, "output_unit": "kJ/mol", "category": "energy", "formula": "1 eV/particle = 96.485 kJ/mol"},
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _convert_temp(self, value, from_u, to_u):
        """Special temperature conversion."""
        # Convert to Kelvin first
        if from_u == "K":
            k_val = value
        elif from_u == "C" or from_u == "\u00b0C":
            k_val = value + 273.15
        elif from_u == "F" or from_u == "\u00b0F":
            k_val = (value - 32) * 5.0 / 9.0 + 273.15
        else:
            raise ChemMCPError(f"Unknown temperature unit: {from_u}")
        # Convert from Kelvin to target
        if to_u == "K":
            return round(k_val, 6)
        elif to_u == "C" or to_u == "\u00b0C":
            return round(k_val - 273.15, 6)
        elif to_u == "F" or to_u == "\u00b0F":
            return round((k_val - 273.15) * 9.0 / 5.0 + 32, 6)
        else:
            raise ChemMCPError(f"Unknown temperature unit: {to_u}")

    def _run_base(self, value: float, from_unit: str, to_unit: str) -> dict:
        fu = from_unit.strip()
        tu = to_unit.strip()

        # Special: eV to kJ/mol (per-particle to per-mole)
        if fu.lower() == "ev" and tu.lower() == "kj/mol":
            kj_mol = value * 96.48533212
            return {"input_value": value, "input_unit": fu, "output_value": round(kj_mol, 6), "output_unit": tu, "category": "energy", "formula": "1 eV = 96.485 kJ/mol"}
        if fu.lower() == "kj/mol" and tu.lower() == "ev":
            ev_val = value / 96.48533212
            return {"input_value": value, "input_unit": fu, "output_value": round(ev_val, 10), "output_unit": tu, "category": "energy", "formula": "1 kJ/mol = 0.01036 eV"}

        cat_f, units_f = _find_unit(fu)
        cat_t, units_t = _find_unit(tu)

        if cat_f is None:
            raise ChemMCPError(f"Unknown unit '{fu}'. Categories: {list(UNIT_CATEGORIES.keys())}")
        if cat_t is None:
            raise ChemMCPError(f"Unknown unit '{tu}'. Categories: {list(UNIT_CATEGORIES.keys())}")

        # Temperature special case
        if cat_f == "temperature":
            if cat_t != "temperature":
                raise ChemMCPError(f"Cannot convert between temperature and non-temperature units.")
            result = self._convert_temp(value, fu, tu)
            return {"input_value": value, "input_unit": fu, "output_value": result, "output_unit": tu, "category": "temperature", "formula": "Temperature conversion formula"}

        # Cross-category conversion via SI
        factor_f, si_f = units_f[fu]
        factor_t, si_t = units_t[tu]
        si_value = value * factor_f
        result = si_value / factor_t

        return {
            "input_value": value,
            "input_unit": fu,
            "output_value": round(result, 10),
            "output_unit": tu,
            "category": cat_f,
            "formula": f"{value} {fu} \u00d7 {factor_f} / {factor_t} = {round(result, 10)} {tu}",
        }

    def _run_text(self, conversion_str: str) -> dict:
        parts = conversion_str.strip().split()
        if len(parts) < 3:
            raise ChemMCPError("Need at least 3 values: <value> <from_unit> <to_unit>")
        val = float(parts[0])
        fu = parts[1]
        tu = " ".join(parts[2:])  # handle multi-word like "kJ/mol"
        return self._run_base(val, fu, tu)
