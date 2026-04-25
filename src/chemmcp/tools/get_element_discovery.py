import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.periodic_table import get_element

logger = logging.getLogger(__name__)

# Element discovery data
DISCOVERY_DATA: dict = {
    "H":  {"discoverer": "Henry Cavendish", "year": 1766, "place": "London, UK", "name_origin": "Greek 'hydro genes' (water-former)", "etymology": "Named because it produces water when burned."},
    "He": {"discoverer": "William Ramsay, P.T. Cleve, Nils Langlet", "year": 1895, "place": "Uppsala, Sweden / London, UK", "name_origin": "Greek 'helios' (Sun)", "etymology": "First detected in the Sun's spectrum (1868) by Jules Janssen and Joseph Lockyer before being found on Earth."},
    "Li": {"discoverer": "Johan August Arfwedson", "year": 1817, "place": "Stockholm, Sweden", "name_origin": "Greek 'lithos' (stone)", "etymology": "Discovered in the mineral petalite."},
    "Be": {"discoverer": "Louis Nicolas Vauquelin", "year": 1798, "place": "Paris, France", "name_origin": "Greek 'beryllos' (beryl)", "etymology": "Found in emeralds and beryl."},
    "B":  {"discoverer": "Joseph Louis Gay-Lussac, Louis Jacques Thénard (isolated by Humphry Davy)", "year": 1808, "place": "Paris, France / London, UK", "name_origin": "Borax ('borax')", "etymology": "Derived from the Arabic 'buraq' or Persian 'burah'."},
    "C":  {"discoverer": "Known since antiquity", "year": None, "place": "Ancient civilizations", "name_origin": "Latin 'carbo' (coal/charcoal)", "etymology": "One of the first elements known to humans; charcoal, diamond, graphite are natural forms."},
    "N":  {"discoverer": "Daniel Rutherford", "year": 1772, "place": "Edinburgh, Scotland", "name_origin": "Greek 'nitron genes' (nitre/soda-forming)", "etymology": "Rutherford called it 'noxious air' when he isolated it."},
    "O":  {"discoverer": "Carl Wilhelm Scheele (1771), Antoine Lavoisier (1774, published first)", "year": 1774, "place": "Uppsala, Sweden / Paris, France", "name_origin": "Greek 'oxy genes' (acid-former)", "etymology": "Lavoisier mistakenly believed oxygen was in all acids."},
    "F":  {"discoverer": "Henri Moissan", "year": 1886, "place": "Paris, France", "name_origin": "Latin 'fluere' (flow)", "etymology": "First isolated by electrolysis of KF in liquid HF."},
    "Ne": {"discoverer": "William Ramsay, Morris Travers", "year": 1898, "place": "London, UK", "name_origin": "Greek 'neos' (new)", "etymology": "Discovered by fractional distillation of liquid air."},
    "Na": {"discoverer": "Humphry Davy", "year": 1807, "place": "London, UK", "name_origin": "English/Latin 'sodium' / Arabic 'suda' (headache remedy)", "etymology": "Isolated by electrolysis of caustic soda (NaOH)."},
    "Mg": {"discoverer": "Joseph Black (recognized as element), Humphry Davy (isolated)", "year": 1755, "place": "Edinburgh, Scotland / London, UK", "name_origin": "Magnesia (district in Thessaly, Greece)", "etymology": "Named after the mineral magnesite."},
    "Al": {"discoverer": "Hans Christian Ørsted", "year": 1825, "place": "Copenhagen, Denmark", "name_origin": "Latin 'alumen' (alum)", "etymology": "Ørsted isolated impure aluminum; Friedrich Wöhler purified it in 1827."},
    "Si": {"discoverer": "Jöns Jacob Berzelius", "year": 1824, "place": "Stockholm, Sweden", "name_origin": "Latin 'silex' (flint)", "etymology": "Berzelius prepared amorphous silicon by heating potassium with silicon tetrafluoride."},
    "P":  {"discoverer": "Hennig Brand", "year": 1669, "place": "Hamburg, Germany", "name_origin": "Greek 'phosphoros' (light-bringer)", "etymology": "Brand discovered it while trying to create the philosopher's stone from urine."},
    "S":  {"discoverer": "Known since antiquity", "year": None, "place": "Ancient civilizations", "name_origin": "Sanskrit 'sulvere' (Latin: 'sulfur')", "etymology": "Known to ancient Chinese, Greeks, and Egyptians; referred to as 'brimstone' in the Bible."},
    "Cl": {"discoverer": "Carl Wilhelm Scheele", "year": 1774, "place": "Uppsala, Sweden", "name_origin": "Greek 'chloros' (pale green)", "etymology": "Scheele produced chlorine but thought it was a compound; Davy proved it was an element in 1810."},
    "Ar": {"discoverer": "Lord Rayleigh, William Ramsay", "year": 1894, "place": "London, UK", "name_origin": "Greek 'argos' (lazy/inert)", "etymology": "Discovered by noticing that nitrogen from air was denser than nitrogen from chemical sources."},
    "K":  {"discoverer": "Humphry Davy", "year": 1807, "place": "London, UK", "name_origin": "English/Latin 'potassium' / Arabic 'qali' (alkali)", "etymology": "Isolated by electrolysis of potash (K2CO3)."},
    "Ca": {"discoverer": "Humphry Davy", "year": 1808, "place": "London, UK", "name_origin": "Latin 'calx' (lime)", "etymology": "Davy obtained calcium by electrolysis of lime (CaO)."},
    "Fe": {"discoverer": "Known since antiquity (~3500 BCE)", "year": None, "place": "Anatolia (modern Turkey)", "name_origin": "English/Germanic 'iron'", "etymology": "One of the oldest known metals; iron age began ~1200 BCE."},
    "Cu": {"discoverer": "Known since antiquity (~9000 BCE)", "year": None, "place": "Middle East", "name_origin": "Latin 'cuprum' (Cyprus island)", "etymology": "Mined extensively in Cyprus during Roman times."},
    "Zn": {"discoverer": "Known in India/China since antiquity; recognized as element by Andreas Marggraf", "year": 1746, "place": "Germany", "name_origin": "German 'Zink' (tooth-like, pointed)", "etymology": "Marggraf is credited with recognizing zinc as a distinct metal."},
    "Ag": {"discoverer": "Known since antiquity", "year": None, "place": "Ancient Anatolia", "name_origin": "English 'silver' (Proto-Germanic 'silubr')", "etymology": "One of the first seven metals known to humanity; mentioned in Genesis."},
    "Au": {"discoverer": "Known since antiquity", "year": None, "place": "Ancient Egypt/Mesopotamia", "name_origin": "Old English 'geol' (yellow)", "etymology": "Treasured since prehistoric times; found naturally in pure form."},
    "Hg": {"discoverer": "Known since antiquity (~1500 BCE)", "year": None, "place": "Ancient Egypt/China", "name_origin": "Latin 'hydrargyrum' (liquid silver)", "etymology": "The only metal that is liquid at room temperature; known to ancient Chinese and Egyptians."},
    "Pb": {"discoverer": "Known since antiquity (~7000 BCE)", "year": None, "place": "Middle East", "name_origin": "English 'lead' (Proto-Germanic 'laudan')", "etymology": "One of the oldest metals; used by Romans for plumbing (hence the symbol Pb from Latin plumbum)."},
    "U":  {"discoverer": "Martin Heinrich Klaproth", "year": 1789, "place": "Berlin, Germany", "name_origin": "Planet Uranus", "etymology": "Named after the planet Uranus, which had been discovered just 8 years earlier (1781)."},
}


@ChemMCPManager.register_tool
class GetElementDiscovery(BaseTool):
    __version__ = "0.1.0"
    name = "GetElementDiscovery"
    func_name = 'get_element_discovery'
    description = "Get discovery history of an element including discoverer, year, place of discovery, name origin, and etymology."
    implementation_description = "Uses a built-in historical database covering elements discovered from antiquity through the modern era. Returns discoverer(s), year, location, naming origin, and interesting historical notes."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Table", "History of Chemistry", "Element Discovery"]
    required_envs = []

    code_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol (e.g., O, U, He)'),
    ]
    text_input_sig = [
        ('element', 'str', 'N/A', 'Element symbol'),
    ]
    output_sig = [
        ('element', 'str', 'Element symbol'),
        ('discoverer', 'str', 'Discoverer(s)'),
        ('year', 'int/str', 'Year of discovery'),
        ('place', 'str', 'Place of discovery'),
        ('name_origin', 'str', 'Origin of the element name'),
        ('etymology', 'str', 'Etymological details'),
    ]
    
        
    examples = [
        {'code_input': {'element': 'O'}, 'text_input': {'element': 'O'}, 'output': {'element': 'O', 'discoverer': 'Carl Wilhelm Scheele', 'year': 1774, 'place': 'Sweden', 'name_origin': 'Greek', 'etymology': '...'}},
        {'code_input': {'element': 'U'}, 'text_input': {'element': 'U'}, 'output': {'element': 'U', 'discoverer': 'Martin Heinrich Klaproth', 'year': 1789, 'place': 'Germany', 'name_origin': 'Uranus', 'etymology': '...'}},
    ]
    def _run_base(self, element: str) -> dict:
        data = get_element(element)
        if data is None:
            raise ChemMCPInputError(f"Element not found: {element}")
        sym = data["symbol"]
        if sym not in DISCOVERY_DATA:
            return {
                "element": sym,
                "note": f"Discovery history not available for {sym} in this database. Data covers ~30 historically significant elements.",
            }
        d = DISCOVERY_DATA[sym]
        return {
            "element": sym,
            "discoverer": d["discoverer"],
            "year": d["year"],
            "place": d["place"],
            "name_origin": d["name_origin"],
            "etymology": d["etymology"],
        }


if __name__ == "__main__":
    run_mcp_server()
