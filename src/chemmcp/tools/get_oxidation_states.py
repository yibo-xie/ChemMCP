import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)

# Common oxidation states database with stability notes
OXIDATION_STATES: dict = {
    "H":  {"states": [-1, +1],   "most_common": +1, "notes": {"+1": "common in acids, water, organic compounds", "-1": "hydrides (NaH, CaH2)"}},
    "He": {"states": [0],        "most_common": 0,  "notes": {0: "noble gas, inert"}},
    "Li": {"states": [+1],      "most_common": +1, "notes": {+1: "only oxidation state, highly electropositive"}},
    "Be": {"states": [+2],      "most_common": +2, "notes": {+2: "only common state"}},
    "B":  {"states": [+3],      "most_common": +3, "notes": {+3: "borates, boric acid"}},
    "C":  {"states": [-4, -3, -2, -1, +1, +2, +3, +4], "most_common": [+2, +4, -4],
           "notes": {-4: "methane (CH4)", -2: "CO", +2: "CO", +4: "CO2, carbonates"}},
    "N":  {"states": [-3, -2, -1, +1, +2, +3, +4, +5], "most_common": [-3, +3, +5],
           "notes": {-3: "ammonia (NH3)", +3: "nitrites", +5: "nitrates"}},
    "O":  {"states": [-2, -1, -0.5, +1, +2], "most_common": -2,
           "notes": {-2: "oxides, water", -1: "peroxides (H2O2)", "+1": "OF2", "+2": "OF2 alternative"}},
    "F":  {"states": [-1],      "most_common": -1, "notes": {-1: "most electronegative element, only -1 state"}},
    "Ne": {"states": [0],       "most_common": 0,  "notes": {0: "noble gas, inert"}},
    "Na": {"states": [+1],      "most_common": +1, "notes": {+1: "only common state, alkali metal"}},
    "Mg": {"states": [+2],      "most_common": +2, "notes": {+2: "only common state"}},
    "Al": {"states": [+3],      "most_common": +3, "notes": {+3: "aluminum compounds"}},
    "Si": {"states": [-4, -3, -2, -1, +1, +2, +3, +4], "most_common": [+4, -4],
           "notes": {-4: "silane (SiH4)", +4: "silica (SiO2), silicates"}},
    "P":  {"states": [-3, -2, +1, +3, +4, +5], "most_common": [-3, +5],
           "notes": {-3: "phosphine (PH3)", +3: "phosphites", +5: "phosphates"}},
    "S":  {"states": [-2, -1, +1, +2, +3, +4, +5, +6], "most_common": [-2, +6, +4],
           "notes": {-2: "sulfides (H2S)", +4: "SO2", +6: "sulfates"}},
    "Cl": {"states": [-1, +1, +3, +4, +5, +6, +7], "most_common": [-1, +1, +3, +5, +7],
           "notes": {-1: "chlorides (HCl)", +1: "hypochlorites", +3: "chlorites", +5: "chlorates", +7: "perchlorates"}},
    "Ar": {"states": [0],       "most_common": 0,  "notes": {0: "noble gas, inert"}},
    "K":  {"states": [+1],      "most_common": +1, "notes": {+1: "alkali metal"}},
    "Ca": {"states": [+2],      "most_common": +2, "notes": {+2: "alkaline earth metal"}},
    "Sc": {"states": [+3],      "most_common": +3, "notes": {+3: "rare earth chemistry"}},
    "Ti": {"states": [+2, +3, +4], "most_common": [+4, +3],
           "notes": {+4: "titanium dioxide (TiO2)", +3: "aqueous Ti(III) complexes"}},
    "V":  {"states": [+2, +3, +4, +5], "most_common": [+5, +4, +3],
           "notes": {+5: "vanadates", +4: "VO2", +3: "V2O3"}},
    "Cr": {"states": [+1, +2, +3, +4, +5, +6], "most_common": [+3, +6],
           "notes": {+3: "Cr(III) salts (green)", +6: "chromates/dichromates (orange)"}},
    "Mn": {"states": [+2, +3, +4, +5, +6, +7], "most_common": [+2, +4, +7],
           "notes": {+2: "Mn(II) salts (pale pink)", +4: "MnO2", +7: "permanganates (purple)"}},
    "Fe": {"states": [+2, +3, +6], "most_common": [+2, +3],
           "notes": {+2: "ferrous (FeSO4, pale green)", +3: "ferric (Fe2O3, reddish-brown)", "+6": "ferrates (rare)"}},
    "Co": {"states": [+2, +3], "most_common": [+2, +3],
           "notes": {+2: "cobalt(II) (pink)", +3: "cobalt(III) (blue-green)"}},
    "Ni": {"states": [+2, +3], "most_common": +2, "notes": {+2: "nickel(II) (green)"}},
    "Cu": {"states": [+1, +2], "most_common": [+2, +1],
           "notes": {+1: "cuprous (Cu2O, red)", +2: "cupric (CuSO4, blue)"}},
    "Zn": {"states": [+2],      "most_common": +2, "notes": {+2: "only common state"}},
    "Ga": {"states": [+1, +2, +3], "most_common": +3, "notes": {+3: "gallium(III)"}},
    "Ge": {"states": [-4, +2, +4], "most_common": [+4, +2],
           "notes": {+4: "germanium dioxide", +2: "Ge(II) compounds"}},
    "As": {"states": [-3, +3, +5], "most_common": [+3, +5],
           "notes": {-3: "arsine (AsH3)", +3: "arsenous oxide", +5: "arsenate"}},
    "Se": {"states": [-2, +2, +4, +6], "most_common": [-2, +4, +6],
           "notes": {-2: "selenides", +4: "selenites (SeO3^2-)", +6: "selenates"}},
    "Br": {"states": [-1, +1, +3, +4, +5, +7], "most_common": [-1, +5],
           "notes": {-1: "bromides (HBr)", +5: "bromates"}},
    "Kr": {"states": [0, +2],   "most_common": 0,  "notes": {0: "noble gas", +2: "KrF2 (very rare)"}},
    "Rb": {"states": [+1],      "most_common": +1, "notes": {+1: "alkali metal"}},
    "Sr": {"states": [+2],      "most_common": +2, "notes": {+2: "alkaline earth metal"}},
    "Y":  {"states": [+3],      "most_common": +3, "notes": {+3: "rare earth"}},
    "Zr": {"states": [+4],      "most_common": +4, "notes": {+4: "zirconium(IV)"}},
    "Nb": {"states": [+2, +3, +4, +5], "most_common": +5, "notes": {+5: "niobium(V)"}},
    "Mo": {"states": [+2, +3, +4, +5, +6], "most_common": [+6, +5, +4],
           "notes": {+6: "molybdates", +4: "MoO2"}},
    "Tc": {"states": [+4, +7],  "most_common": +7, "notes": {+7: "pertechnetate (TcO4^-)"}},
    "Ru": {"states": [+2, +3, +4, +6, +8], "most_common": [+3, +4],
           "notes": {+3: "Ru(III)", +4: "RuO2"}},
    "Rh": {"states": [+1, +2, +3, +4, +6], "most_common": [+3], "notes": {+3: "rhodium(III)"}},
    "Pd": {"states": [+2, +4],  "most_common": [+2, +4], "notes": {+2: "Pd(II)", +4: "Pd(IV)"}},
    "Ag": {"states": [+1, +2, +3], "most_common": +1, "notes": {+1: "silver(I) (AgNO3)"}},
    "Cd": {"states": [+2],      "most_common": +2, "notes": {+2: "cadmium(II)"}},
    "In": {"states": [+1, +2, +3], "most_common": +3, "notes": {+3: "indium(III)"}},
    "Sn": {"states": [-4, +2, +4], "most_common": [+4, +2],
           "notes": {+4: "stannic (SnO2)", +2: "stannous (SnCl2)", -4: "stannide"}},
    "Sb": {"states": [-3, +3, +5], "most_common": [+3, +5],
           "notes": {+3: "antimonous (Sb2O3)", +5: "antimonic (Sb2O5)", -3: "stibine"}},
    "Te": {"states": [-2, +2, +4, +6], "most_common": [-2, +4, +6],
           "notes": {-2: "tellurides", +4: "tellurites", +6: "tellurates"}},
    "I":  {"states": [-1, +1, +3, +5, +7], "most_common": [-1, +5],
           "notes": {-1: "iodides (HI)", +5: "iodates (IO3^-)", +7: "periodates (IO4^-)"}},
    "Xe": {"states": [+2, +4, +6, +8], "most_common": [+4, +6],
           "notes": {+4: "XeF4", +6: "XeF6/XeO3", +8: "XeO4", +2: "XeF2"}},
    "Cs": {"states": [+1],      "most_common": +1, "notes": {+1: "alkali metal"}},
    "Ba": {"states": [+2],      "most_common": +2, "notes": {+2: "alkaline earth metal"}},
    "La": {"states": [+3],      "most_common": +3, "notes": {+3: "lanthanum(III)"}},
    "Ce": {"states": [+3, +4],  "most_common": [+3, +4], "notes": {+3: "Ce(III)", +4: "Ce(IV) (CeO2)"}},
    "Pr": {"states": [+3, +4],  "most_common": +3, "notes": {+3: "Pr(III)", +4: "Pr(IV) (PrO2)"}},
    "Nd": {"states": [+3],      "most_common": +3, "notes": {+3: "Nd(III)"}},
    "Pm": {"states": [+3],      "most_common": +3, "notes": {+3: "Pm(III)"}},
    "Sm": {"states": [+2, +3],  "most_common": [+3, +2], "notes": {+3: "Sm(III)", +2: "Sm(II)"}},
    "Eu": {"states": [+2, +3],  "most_common": [+2, +3], "notes": {+2: "Eu(II)", +3: "Eu(III)"}},
    "Gd": {"states": [+3],      "most_common": +3, "notes": {+3: "Gd(III)"}},
    "Tb": {"states": [+3, +4],  "most_common": +3, "notes": {+3: "Tb(III)", +4: "Tb(IV)"}},
    "Dy": {"states": [+3],      "most_common": +3, "notes": {+3: "Dy(III)"}},
    "Ho": {"states": [+3],      "most_common": +3, "notes": {+3: "Ho(III)"}},
    "Er": {"states": [+3],      "most_common": +3, "notes": {+3: "Er(III)"}},
    "Tm": {"states": [+3],      "most_common": +3, "notes": {+3: "Tm(III)"}},
    "Yb": {"states": [+2, +3],  "most_common": [+3, +2], "notes": {+3: "Yb(III)", +2: "Yb(II)"}},
    "Lu": {"states": [+3],      "most_common": +3, "notes": {+3: "Lu(III)"}},
    "Hf": {"states": [+4],      "most_common": +4, "notes": {+4: "hafnium(IV)"}},
    "Ta": {"states": [+5],      "most_common": +5, "notes": {+5: "tantalum(V)"}},
    "W":  {"states": [+2, +3, +4, +5, +6], "most_common": [+6],
           "notes": {+6: "tungstates (WO4^2-)"}},
    "Re": {"states": [-1, +1, +4, +5, +6, +7], "most_common": [+7, +4],
           "notes": {+7: "perrhenate (ReO4^-)", +4: "ReO2"}},
    "Os": {"states": [+2, +3, +4, +6, +8], "most_common": [+4, +8],
           "notes": {+8: "osmium tetroxide (OsO4)", +4: "OsO2"}},
    "Ir": {"states": [+3, +4, +6], "most_common": [+3, +4],
           "notes": {+3: "iridium(III)", +4: "IrO2"}},
    "Pt": {"states": [+2, +4],  "most_common": [+4, +2],
           "notes": {+4: "Pt(IV)", +2: "Pt(II)"}},
    "Au": {"states": [+1, +3],  "most_common": [+3, +1],
           "notes": {+1: "aurous (AuCN)", +3: "auric (HAuCl4)"}},
    "Hg": {"states": [+1, +2],  "most_common": +2, "notes": {+2: "mercuric (HgCl2)", +1: "mercurous (Hg2Cl2)"}},
    "Tl": {"states": [+1, +3],  "most_common": [+1, +3],
           "notes": {+1: "thallous (Tl2O)", +3: "thallic (Tl2O3)"}},
    "Pb": {"states": [+2, +4],  "most_common": [+2, +4],
           "notes": {+2: "plumbous (Pb(NO3)2)", +4: "plumbic (PbO2)"}},
    "Bi": {"states": [+3, +5],  "most_common": +3, "notes": {+3: "bismuth(III)", +5: "bismuthate(V)"}},
    "Po": {"states": [+2, +4, +6], "most_common": [+4, +2], "notes": {+4: "PoO2", +2: "PoCl2"}},
    "At": {"states": [+1, +3, +5, +7], "most_common": [+1, -1], "notes": {+1: "astatine(I)", -1: "astatides"}},
    "Rn": {"states": [+2],      "most_common": +2, "notes": {+2: "RnF2 (theoretical)"}},
    "Fr": {"states": [+1],      "most_common": +1, "notes": {+1: "alkali metal"}},
    "Ra": {"states": [+2],      "most_common": +2, "notes": {+2: "alkaline earth metal"}},
    "Ac": {"states": [+3],      "most_common": +3, "notes": {+3: "actinium(III)"}},
    "Th": {"states": [+4],      "most_common": +4, "notes": {+4: "thorium(IV) as ThO2"}},
    "Pa": {"states": [+3, +4, +5], "most_common": +5, "notes": {+5: "protactinyl (PaO2+)"}},
    "U":  {"states": [+3, +4, +5, +6], "most_common": [+6, +4],
           "notes": {+6: "uranyl (UO2^2+)", +4: "UO2", +3: "U(III)"}},
    "Np": {"states": [+3, +4, +5, +6], "most_common": [+5, +4],
           "notes": {+5: "neptunyl (NpO2^+)"}},
    "Pu": {"states": [+3, +4, +5, +6], "most_common": [+4, +3],
           "notes": {+4: "PuO2", +3: "Pu(III) solutions"}},
}


@ChemMCPManager.register_tool
class GetOxidationStates(BaseTool):
    __version__ = "0.1.0"
    name = "GetOxidationStates"
    func_name = 'get_oxidation_states'
    description = "Query common oxidation states of an element with stability information."
    implementation_description = "Uses a built-in database of oxidation states for all common elements. Returns all known oxidation states, the most common ones, and descriptive notes for each state including typical compounds."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "Oxidation States", "Redox Chemistry"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., Fe, Mn, Cl)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('element', 'str', 'Element symbol'),
        ('oxidation_states', 'list', 'All known oxidation states'),
        ('most_common', 'list/int', 'Most common oxidation state(s)'),
        ('state_details', 'dict', 'Details for each oxidation state with examples'),
    ]
    
    examples = [
        {'code_input': {'element': 'Fe'}, 'text_input': {'element': 'Fe'}, 'output': {'element': 'Fe', 'oxidation_states': [+2, +3, +6], 'most_common': [+2, +3], 'state_details': {...}}},
        {'code_input': {'element': 'Mn'}, 'text_input': {'element': 'Mn'}, 'output': {'element': 'Mn', 'oxidation_states': [+2,+3,+4,+5,+6,+7], 'most_common': [+2, +4, +7], 'state_details': {...}}},
    ]
    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        sym = data["symbol"]
        if sym not in OXIDATION_STATES:
            return {
                "element": sym,
                "oxidation_states": [],
                "most_common": None,
                "state_details": {},
                "note": f"Oxidation state data not available for {sym}. Data covers elements H through Pu.",
            }
        ox = OXIDATION_STATES[sym]
        return {
            "element": sym,
            "oxidation_states": ox["states"],
            "most_common": ox["most_common"],
            "state_details": ox["notes"],
        }


if __name__ == "__main__":
    run_mcp_server()
