import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_POURBAIX_DATA = {
    "H+/H2":      {"E0": 0.0,     "n": 2, "reaction": "2H+ + 2e- <=> H2(g)",       "pH_dependent": -0.0591},
    "O2/H2O":     {"E0": 1.229,   "n": 4, "reaction": "O2 + 4H+ + 4e- <=> 2H2O",   "pH_dependent": -0.0591},
    "Fe2+/Fe":    {"E0": -0.44,   "n": 2, "reaction": "Fe2+ + 2e- <=> Fe(s)",      "pH_dependent": 0.0},
    "Fe3+/Fe2+":  {"E0": 0.771,   "n": 1, "reaction": "Fe3+ + e- <=> Fe2+",         "pH_dependent": 0.0},
    "Cu2+/Cu":    {"E0": 0.337,   "n": 2, "reaction": "Cu2+ + 2e- <=> Cu(s)",       "pH_dependent": 0.0},
    "Zn2+/Zn":    {"E0": -0.763,  "n": 2, "reaction": "Zn2+ + 2e- <=> Zn(s)",       "pH_dependent": 0.0},
    "Al3+/Al":    {"E0": -1.662,  "n": 3, "reaction": "Al3+ + 3e- <=> Al(s)",       "pH_dependent": 0.0},
    "Ni2+/Ni":    {"E0": -0.25,   "n": 2, "reaction": "Ni2+ + 2e- <=> Ni(s)",       "pH_dependent": 0.0},
    "Ag+/Ag":     {"E0": 0.7996,  "n": 1, "reaction": "Ag+ + e- <=> Ag(s)",         "pH_dependent": 0.0},
    "Au3+/Au":    {"E0": 1.498,   "n": 3, "reaction": "Au3+ + 3e- <=> Au(s)",       "pH_dependent": 0.0},
    "Pb2+/Pb":    {"E0": -0.126,  "n": 2, "reaction": "Pb2+ + 2e- <=> Pb(s)",       "pH_dependent": 0.0},
    "Cr3+/Cr":    {"E0": -0.744,  "n": 3, "reaction": "Cr3+ + 3e- <=> Cr(s)",       "pH_dependent": 0.0},
    "MnO4-/Mn2+": {"E0": 1.51,    "n": 5, "reaction": "MnO4- + 8H+ + 5e- <=> Mn2+ + 4H2O", "pH_dependent": -0.0944},
    "Co2+/Co":    {"E0": -0.28,   "n": 2, "reaction": "Co2+ + 2e- <=> Co(s)",       "pH_dependent": 0.0},
    "Mg2+/Mg":    {"E0": -2.37,   "n": 2, "reaction": "Mg2+ + 2e- <=> Mg(s)",       "pH_dependent": 0.0},
}


@ChemMCPManager.register_tool
class PourbaixDiagramLookup(BaseTool):
    """Pourbaix (E-pH) diagram lookup and analysis tool."""
    __version__ = "0.1.0"
    name = "PourbaixDiagramLookup"
    func_name = "lookup_pourbaix"
    description = "Look up Pourbaix diagram data for common redox couples, calculate Nernst potentials at given pH, and determine species stability regions."
    implementation_description = "Built-in database of ~15 redox couples with standard potentials and pH dependencies. Computes Nernst equation E = E0 - (0.05916/n)*pH."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Pourbaix", "Redox", "Thermodynamics"]
    required_envs = []

    code_input_sig = [
        ("couple_name", "str", "N/A", "Redox couple name (e.g., 'Fe3+/Fe2+', 'Cu2+/Cu'). Use 'list' for all."),
        ("ph", "float", "7.0", "pH value. Default 7.0."),
        ("temperature_k", "float", "298.15", "Temperature (K). Default 298.15."),
        ("activity_oxidized", "float", "1.0", "Activity of oxidized species. Default 1.0."),
        ("activity_reduced", "float", "1.0", "Activity of reduced species. Default 1.0."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: couple_name [pH] [T_K] [a_ox] [a_red]."),
    ]

    output_sig = [
        ("couple_info", "dict", "Couple info: E0, n, reaction, type."),
        ("nernst_potential_v", "float", "Nernst potential at given conditions (V vs SHE)."),
        ("water_stability_window", "dict", "Water window: E_H2 and E_O2 (V vs SHE)."),
        ("stability_analysis", "str", "Analysis of position relative to water stability window."),
        ("available_couples", "list", "All available couples (when couple_name='list')."),
    ]

    examples = [
        {
            "code_input": {
                "couple_name": "Fe3+/Fe2+",
                "ph": 7.0,
                "temperature_k": 298.15,
                "activity_oxidized": 1.0,
                "activity_reduced": 1.0,
            },
            "text_input": {
                "input_params": "Fe3+/Fe2+ 7.0",
            },
            "output": {
                "couple_info": {"E0": 0.771, "n": 1, "reaction": "Fe3+ + e- <=> Fe2+"},
                "nernst_potential_v": 0.771,
                "water_stability_window": {"E_H2": -0.414, "E_O2": 0.815},
                "stability_analysis": "within water stability window",
                "available_couples": None,
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618
        self.F = 96485.33212

    def _get_water_window(self, ph, T=298.15):
        factor = (2.303 * self.R * T) / self.F
        return {"E_H2": round(-factor * ph, 4), "E_O2": round(1.229 - factor * ph, 4)}

    def _run_base(self, couple_name, ph=7.0, temperature_k=298.15,
                   activity_oxidized=1.0, activity_reduced=1.0) -> dict:
        if couple_name.lower() == "list":
            return {
                "couple_info": None, "nernst_potential_v": None,
                "water_stability_window": self._get_water_window(ph, temperature_k),
                "stability_analysis": f"Available: {sorted(_POURBAIX_DATA.keys())}",
                "available_couples": sorted(_POURBAIX_DATA.keys()),
            }

        data = _POURBAIX_DATA.get(couple_name)
        if data is None:
            matches = [k for k in _POURBAIX_DATA if couple_name.lower() in k.lower()]
            if matches:
                couple_name = matches[0]
                data = _POURBAIX_DATA[couple_name]
            else:
                raise ChemMCPError(f"Couple '{couple_name}' not found. Available: {sorted(_POURBAIX_DATA.keys())}")

        E0 = data["E0"]; n = data["n"]; ph_dep = data.get("pH_dependent", 0.0)
        factor = (2.303 * self.R * temperature_k) / self.F
        E_nernst = E0 - (factor / n) * ph_dep * ph
        if activity_oxidized > 0 and activity_reduced > 0:
            E_nernst += (factor / n) * math.log10(activity_oxidized / activity_reduced)

        window = self._get_water_window(ph, temperature_k)
        within = window["E_H2"] <= round(E_nernst, 3) <= window["E_O2"]
        if within:
            analysis = f"{couple_name} (E={round(E_nernst, 3)}V) within water window ({window['E_H2']} to {window['E_O2']}V) at pH {ph}."
        elif E_nernst < window["E_H2"]:
            analysis = f"{couple_name} BELOW water window (E_H2={window['E_H2']}V) at pH {ph}."
        else:
            analysis = f"{couple_name} ABOVE water window (E_O2={window['E_O2']}V) at pH {ph}."

        logger.info(f"Pourbaix: {couple_name} at pH={ph}, E={round(E_nernst, 4)}V")
        return {
            "couple_info": dict(data), "nernst_potential_v": round(E_nernst, 6),
            "water_stability_window": window, "stability_analysis": analysis,
            "available_couples": None,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            kwargs = {"couple_name": parts[0]}
            if len(parts) > 1: kwargs["ph"] = float(parts[1])
            if len(parts) > 2: kwargs["temperature_k"] = float(parts[2])
            if len(parts) > 3: kwargs["activity_oxidized"] = float(parts[3])
            if len(parts) > 4: kwargs["activity_reduced"] = float(parts[4])
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}")
