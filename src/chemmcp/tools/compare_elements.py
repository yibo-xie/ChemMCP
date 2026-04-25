import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element, PERIODIC_TABLE

logger = logging.getLogger(__name__)

SUPPORTED_COMPARE_PROPERTIES = {
    "atomic_weight": {"label": "Atomic Weight (u)", "key": "atomic_weight"},
    "electronegativity": {"label": "Pauling Electronegativity", "key": "en_pauling"},
    "atomic_number": {"label": "Atomic Number", "key": "atomic_number"},
    "ionization_energy_1": {"label": "First Ionization Energy (kJ/mol)", "special": True},
    "electron_affinity": {"label": "Electron Affinity (kJ/mol)", "special": True},
    "period": {"label": "Period", "key": "period"},
    "group": {"label": "Group", "key": "group"},
}

# IE1 data for comparison
_IE1 = {"H":1312,"He":2372,"Li":520,"Be":900,"B":801,"C":1086,"N":1402,"O":1314,"F":1681,"Ne":2081,
        "Na":496,"Mg":738,"Al":577,"Si":787,"P":1012,"S":1000,"Cl":1251,"Ar":1521,"K":419,"Ca":590,
        "Sc":633,"Ti":658,"V":650,"Cr":653,"Mn":717,"Fe":762,"Co":760,"Ni":737,"Cu":745,"Zn":906,
        "Ga":579,"Ge":762,"As":947,"Se":941,"Br":1140,"Kr":1351,"Rb":403,"Sr":550,"Y":600,"Zr":640,
        "Nb":652,"Mo":685,"Tc":702,"Ru":711,"Rh":720,"Pd":804,"Ag":731,"Cd":867,"In":558,"Sn":709,
        "Sb":834,"Te":869,"I":1008,"Xe":1170,"Cs":376,"Ba":503,"La":538,"Ce":534,"Pr":527,"Nd":533,
        "Pm":538,"Sm":544,"Eu":547,"Gd":593,"Tb":565,"Dy":573,"Ho":581,"Er":589,"Tm":597,"Yb":603,
        "Lu":524,"Hf":659,"Ta":761,"W":770,"Re":760,"Os":840,"Ir":880,"Pt":870,"Au":890,"Hg":1007,
        "Tl":589,"Pb":716,"Bi":703}

_EA = {"H":72.8,"Li":59.6,"B":26.7,"C":121.9,"N":-6.8,"O":141.0,"F":328.0,"Na":52.9,"Al":42.5,
       "Si":133.6,"P":72.0,"S":200.4,"Cl":349.0,"K":48.4,"Ca":2.37,"Sc":18.1,"Ti":7.6,"V":50.6,
       "Cr":64.3,"Fe":14.5,"Co":63.9,"Ni":112.0,"Cu":118.4,"Ga":28.9,"Ge":119.0,"As":78.2,"Se":195.0,
       "Br":324.5,"Rb":46.9,"Sr":5.03,"Y":29.6,"Zr":41.1,"Nb":86.2,"Mo":71.9,"Tc":53.0,"Ru":101.3,
       "Rh":109.7,"Pd":53.7,"Ag":125.6,"In":28.9,"Sn":107.3,"Sb":101.2,"Te":190.2,"I":295.2,"Cs":45.5,
       "Ba":13.95,"La":48.5,"Au":222.8,"Pt":205.3,"Tl":19.2,"Pb":35.2,"Bi":91.2}


@ChemMCPManager.register_tool
class CompareElements(BaseTool):
    __version__ = "0.1.0"
    name = "CompareElements"
    func_name = 'compare_elements'
    description = "Compare properties of multiple elements side by side."
    implementation_description = "Accepts a list of element symbols and a property name, returns a comparison table with values sorted by rank. Supports atomic weight, electronegativity, atomic number, first ionization energy, electron affinity, period, and group."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Comparison", "Element Properties"]
    required_envs = []

    code_input_sig = [
        ('elements', 'list', 'N/A', 'List of element symbols to compare (e.g., [\"Na\", \"K\", \"Li\"])'),
        ('property', 'str', 'electronegativity', 'Property to compare: atomic_weight, electronegativity, atomic_number, ionization_energy_1, electron_affinity, period, group'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Comma-separated element symbols and property, e.g., \"Na,K,Li electronegativity\"'),
    ]
    output_sig = [
        ('comparison', 'dict', 'Ranked comparison table with element values and ranking'),
        ('property', 'str', 'Property that was compared'),
        ('trend_note', 'str', 'Brief trend analysis'),
    ]
    
    examples = [
        {'code_input': {'elements': ['Li', 'Na', 'K'], 'property': 'electronegativity'}, 'text_input': {'query': 'compare Li Na K electronegativity'}, 'output': {'comparison': [...], 'property': 'electronegativity', 'trend_note': '...'}},
    ]
    def _run_base(self, elements: list, property: str = "electronegativity") -> dict:
        if property not in SUPPORTED_COMPARE_PROPERTIES:
            avail = list(SUPPORTED_COMPARE_PROPERTIES.keys())
            raise ChemMCPInputError(f"Unsupported property: '{property}'. Available: {avail}")

        results = {}
        for el in elements:
            data = get_element(el)
            if data is None:
                results[el] = None
                continue
            sym = data["symbol"]
            prop_info = SUPPORTED_COMPARE_PROPERTIES[property]
            if prop_info.get("special"):
                if property == "ionization_energy_1":
                    val = _IE1.get(sym)
                elif property == "electron_affinity":
                    val = _EA.get(sym)
                else:
                    val = None
            else:
                val = data.get(prop_info["key"])
            results[sym] = val

        # Sort by value (None last)
        valid = {k: v for k, v in results.items() if v is not None}
        sorted_items = sorted(valid.items(), key=lambda x: x[1], reverse=True)

        label = SUPPORTED_COMPARE_PROPERTIES[property]["label"]
        return {
            "property": property,
            "property_label": label,
            "comparison": {sym: round(val, 3) if isinstance(val, float) else val for sym, val in results.items()},
            "ranking": [{"rank": i+1, "element": sym, "value": round(val, 3) if isinstance(val, float) else val} for i, (sym, val) in enumerate(sorted_items)],
            "trend_note": self._generate_trend(property, sorted_items),
        }

    def _generate_trend(self, prop, sorted_items):
        if len(sorted_items) < 2:
            return "Insufficient data for trend analysis."
        highest = sorted_items[0]
        lowest = sorted_items[-1]
        return f"Highest: {highest[0]} ({highest[1]}), Lowest: {lowest[0]} ({lowest[1]}). Difference: {abs(highest[1]-lowest[1]) if isinstance(highest[1],(int,float)) else 'N/A'}"


if __name__ == "__main__":
    run_mcp_server()
