import logging
import re
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IdentifyOxidizingAgent(BaseTool):
    """
    Identify the oxidizing agent and reducing agent in a chemical reaction.
    Determines which species is oxidized (reducing agent) and which is reduced (oxidizing agent).
    """
    __version__ = "0.1.0"
    name = "IdentifyOxidizingAgent"
    func_name = "identify_oxidizing_agent"
    description = "Identify the oxidizing agent and reducing agent in a redox reaction, including oxidation state changes and electron transfer details."
    implementation_description = "Parses the reaction equation, identifies redox-active species, tracks oxidation state changes using built-in rules and common oxidation state databases, and classifies oxidizing/reducing agents."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Redox", "Oxidizing Agent", "Reducing Agent", "Oxidation State", "Electron Transfer"]
    required_envs = []

    code_input_sig = [
        ("equation", "str", "N/A", "Chemical equation string (balanced or unbalanced), e.g., 'MnO4- + 5Fe2+ + 8H+ = Mn2+ + 5Fe3+ + 4H2O' or 'H2 + CuO = Cu + H2O' or 'Cl2 + 2NaI = 2NaCl + I2'."),
    ]

    text_input_sig = [
        ("equation", "str", "N/A", "Chemical equation string."),
    ]

    output_sig = [
        ("oxidizing_agent", "str", "The species that gets reduced (accepts electrons)."),
        ("reducing_agent", "str", "The species that gets oxidized (donates electrons)."),
        ("oxidation_process", "str", "Description of what gets oxidized: species, oxidation state change, electrons lost."),
        ("reduction_process", "str", "Description of what gets reduced: species, oxidation state change, electrons gained."),
        ("electrons_transferred", "int", "Total number of electrons transferred in the balanced reaction."),
        ("is_redox", "bool", "Whether this is a redox reaction."),
        ("reaction_type_detail", "str", "Detailed classification of the redox process type."),
    ]

    examples = [
        {
            "code_input": {"equation": "MnO4- + 5Fe2+ + 8H+ = Mn2+ + 5Fe3+ + 4H2O"},
            "text_input": {"equation": "MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O"},
            "output": {
                "oxidizing_agent": "MnO₄⁻ (Permanganate ion)",
                "reducing_agent": "Fe²⁺ (Iron(II) ion)",
                "oxidation_process": "Fe²⁺ → Fe³⁺ + e⁻ | Oxidation state: +II → +III | Each Fe loses 1e⁻",
                "reduction_process": "MnO₄⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O | Mn: +VII → +II | Gains 5e⁻",
                "electrons_transferred": 5,
                "is_redox": True,
                "reaction_type_detail": "Redox — Permanganate oxidation of Fe(II) in acidic medium",
            }
        },
        {
            "code_input": {"equation": "Zn + CuSO4 = ZnSO4 + Cu"},
            "text_input": {"equation": "Zn + CuSO4 = ZnSO4 + Cu"},
            "output": {
                "oxidizing_agent": "Cu²⁺ (in CuSO₄)",
                "reducing_agent": "Zn (Zinc metal)",
                "oxidation_process": "Zn → Zn²⁺ + 2e⁻ | Oxidation state: 0 → +II | Loses 2e⁻",
                "reduction_process": "Cu²⁺ + 2e⁻ → Cu | Oxidation state: +II → 0 | Gains 2e⁻",
                "electrons_transferred": 2,
                "is_redox": True,
                "reaction_type_detail": "Redox — Single displacement (metal activity series)",
            }
        },
        {
            "code_input": {"equation": "Cl2 + 2NaI = 2NaCl + I2"},
            "text_input": {"equation": "Cl2 + 2NaI = NaCl + I2"},
            "output": {
                "oxidizing_agent": "Cl₂ (Chlorine gas)",
                "reducing_agent": "I⁻ (Iodide ion)",
                "oxidation_process": "2I⁻ → I₂ + 2e⁻ | Oxidation state: -I → 0 | Loses 2e⁻ total",
                "reduction_process": "Cl₂ + 2e⁻ → 2Cl⁻ | Oxidation state: 0 → -I | Gains 2e⁻ total",
                "electrons_transferred": 2,
                "is_redox": True,
                "reaction_type_detail": "Redox — Halogen displacement (Cl₂ more electronegative than I₂)",
            }
        },
        {
            "code_input": {"equation": "H2SO4 + 2NaOH = Na2SO4 + 2H2O"},
            "text_input": {"equation": "HCl + NaOH = NaCl + H2O"},
            "output": {
                "oxidizing_agent": "N/A — no electron transfer",
                "reducing_agent": "N/A — no electron transfer",
                "oxidation_process": "N/A — all oxidation states unchanged",
                "reduction_process": "N/A — all oxidation states unchanged",
                "electrons_transferred": 0,
                "is_redox": False,
                "reaction_type_detail": "Acid-Base Neutralization — NOT a redox reaction (no oxidation state changes)",
            }
        },
        {
            "code_input": {"equation": "2KClO3 = 2KCl + 3O2"},
            "text_input": {"equation": "KClO3 = KCl + O2"},
            "output": {
                "oxidizing_agent": "Cl (+V in KClO₃) — also acts as reducing agent",
                "reducing_agent": "O (-II in KClO₃) — also acts as oxidizing agent",
                "oxidation_process": "2O(-II) → O₂(0) + 4e⁻ | Oxygen oxidized from -II to 0",
                "reduction_process": "Cl(+V) + 6e⁻ → Cl(-I) | Chlorine reduced from +V to -I",
                "electrons_transferred": 6,
                "is_redox": True,
                "reaction_type_detail": "Redox — Disproportionation / Decomposition (same element both oxidized and reduced)",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build oxidation state database and known redox patterns."""
        # Common oxidation states for elements
        self._common_ox_states = {
            # Group 1
            "Li": [1], "Na": [1], "K": [1], "Rb": [1], "Cs": [1],
            # Group 2
            "Be": [2], "Mg": [2], "Ca": [2], "Sr": [2], "Ba": [2],
            # Group 12
            "Zn": [2], "Cd": [2], "Hg": [1, 2],
            # Group 13
            "Al": [3], "Ga": [3], "In": [3],
            # Common transition metals
            "Sc": [3], "Ti": [2, 3, 4], "V": [2, 3, 4, 5],
            "Cr": [2, 3, 6], "Mn": [2, 4, 6, 7],
            "Fe": [2, 3], "Co": [2, 3], "Ni": [2],
            "Cu": [1, 2], "Ag": [1], "Au": [1, 3],
            # Nonmetals
            "H": [-1, 1],  # -1 with metals, +1 with nonmetals
            "He": [0], "Ne": [0], "Ar": [0], "Kr": [0], "Xe": [0],
            "F": [-1],
            "Cl": [-1, 1, 3, 5, 7],
            "Br": [-1, 1, 3, 5],
            "I": [-1, 1, 3, 5, 7],
            "O": [-2, -1, 0.5],  # -2 normal, -1 peroxide, -1/2 superoxide
            "S": [-2, 2, 4, 6],
            "N": [-3, 0, 1, 2, 3, 4, 5],
            "P": [-3, 3, 5],
            "C": [-4, -3, -2, -1, 0, 1, 2, 3, 4],
            "Si": [-4, 4],
            "B": [3],
            # Noble/semi-metals
            "Pt": [2, 4], "Pd": [2, 4],
        }

        # Known redox patterns for quick identification
        self._redox_patterns = {
            # Pattern key → (oxidizing_agent, reducing_agent, ox_desc, red_desc, n_e)
            "mno4_fe": (
                "MnO₄⁻ (permanganate)", "Fe²⁺ (iron(II))",
                "Fe²⁺ → Fe³⁺ + e⁻ (ox: +II→+III, loses 1e⁻ per atom)",
                "MnO₄⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O (red: +VII→+II, gains 5e⁻)",
                5
            ),
            "cr2o7_fe": (
                "Cr₂O₇²⁻ (dichromate)", "Fe²⁺ (iron(II))",
                "Fe²⁺ → Fe³⁺ + e⁻ (ox: +II→+III)",
                "Cr₂O₇²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H₂O (red: +VI→+III, gains 6e⁻ total)",
                6
            ),
            "zn_cuso4": (
                "Cu²⁺ (in CuSO₄)", "Zn (zinc metal)",
                "Zn → Zn²⁺ + 2e⁻ (ox: 0→+II, loses 2e⁻)",
                "Cu²⁺ + 2e⁻ → Cu (red: +II→0, gains 2e⁻)",
                2
            ),
            "fe_cuso4": (
                "Cu²⁺ (in CuSO₄)", "Fe (iron metal)",
                "Fe → Fe²⁺ + 2e⁻ (ox: 0→+II, loses 2e⁻)",
                "Cu²⁺ + 2e⁻ → Cu (red: +II→0, gains 2e⁻)",
                2
            ),
            "zn_hcl": (
                "H⁺ (in HCl)", "Zn (zinc metal)",
                "Zn → Zn²⁺ + 2e⁻ (ox: 0→+II)",
                "2H⁺ + 2e⁻ → H₂ (red: +I→0, gains 2e⁻)",
                2
            ),
            "fe_hcl": (
                "H⁺ (in HCl)", "Fe (iron metal)",
                "Fe → Fe²⁺ + 2e⁻ (ox: 0→+II)",
                "2H⁺ + 2e⁻ → H₂ (red: +I→0)",
                2
            ),
            "cl2_nai": (
                "Cl₂ (chlorine gas)", "I⁻ (iodide ion)",
                "2I⁻ → I₂ + 2e⁻ (ox: -I→0, loses 2e⁻ total)",
                "Cl₂ + 2e⁻ → 2Cl⁻ (red: 0→-I, gains 2e⁻ total)",
                2
            ),
            "cl2_nabr": (
                "Cl₂ (chlorine gas)", "Br⁻ (bromide ion)",
                "2Br⁻ → Br₂ + 2e⁻ (ox: -I→0)",
                "Cl₂ + 2e⁻ → 2Cl⁻ (red: 0→-I)",
                2
            ),
            "br2_nai": (
                "Br₂ (bromine gas)", "I⁻ (iodide ion)",
                "2I⁻ → I₂ + 2e⁻ (ox: -I→0)",
                "Br₂ + 2e⁻ → 2Br⁻ (red: 0→-I)",
                2
            ),
            "cu_hno3_dilute": (
                "NO₃⁻ (nitrate, dilute HNO₃)", "Cu (copper metal)",
                "Cu → Cu²⁺ + 2e⁻ (ox: 0→+II)",
                "NO₃⁻ + 4H⁺ + 3e⁻ → NO + 2H₂O (red: +V→+II, gains 3e⁻)",
                6  # LCM of 2 and 3
            ),
            "cu_hno3_conc": (
                "NO₃⁻ (nitrate, conc. HNO₃)", "Cu (copper metal)",
                "Cu → Cu²⁺ + 2e⁻ (ox: 0→+II)",
                "NO₃⁻ + 2H⁺ + e⁻ → NO₂ + H₂O (red: +V→+IV, gains 1e⁻)",
                2
            ),
            "c_h2so4_conc": (
                "H₂SO₄ (conc., sulfuric acid)", "C (carbon)",
                "C → CO₂ + 4e⁻ (ox: 0→+IV, loses 4e⁻)",
                "H₂SO₄ + 2H⁺ + 2e⁻ → SO₂ + 2H₂O (red: +SVI→+IV, gains 2e⁻)",
                4
            ),
            "h2_cuo": (
                "CuO (copper(II) oxide)", "H₂ (hydrogen gas)",
                "H₂ → 2H⁺ + 2e⁻ (ox: 0→+I, loses 2e⁻)",
                "Cu²⁺ + 2e⁻ → Cu (red: +II→0, gains 2e⁻)",
                2
            ),
            "co_fe2o3": (
                "Fe₂O₃ (iron(III) oxide)", "CO (carbon monoxide)",
                "CO → CO₂ + 2e⁻ (ox: C:+II→+IV, loses 2e⁻)",
                "Fe₂O₃ + 6e⁻ → 2Fe + 3O²⁻ (red: Fe:+III→0, gains 6e⁻ total)",
                6
            ),
            "kclo3": (
                "KClO₃ (potassium chlorate)", "KClO₃ (potassium chlorate)",
                "2Cl(+V) → Cl(-I) + 6e⁻ (chlorine reduced in KCl → KCl)",
                "6O(-II) → 3O₂(0) + 12e⁻ (oxygen oxidized in KClO₃ → O₂)",
                6
            ),
            "kclo3_decomp": (
                "Cl in KClO₃ (internal)", "O in KClO₃ (internal)",
                "2O(-II) → O₂(0) + 4e⁻ (oxygen oxidized)",
                "Cl(+V) + 6e⁻ → Cl(-I) (chlorine reduced)",
                6
            ),
            "h2o2_decomp": (
                "H₂O₂ (as oxidizer)", "H₂O₂ (as reducer)",
                "H₂O₂ → O₂ + 2H⁺ + 2e⁻ (ox: O:-I→0)",
                "H₂O₂ + 2H⁺ + 2e⁻ → 2H₂O (red: O:-I→-II)",
                2
            ),
            "na2o2_co2": (
                "Na₂O₂ (peroxide, as oxidizer)", "Na₂O₂ (peroxide, as reducer)",
                "O(-I) → O(0) in O₂ (oxidized)",
                "O(-I) → O(-II) in Na₂CO₃/H₂O (reduced)",
                2
            ),
            "s2o32-_i2": (
                "I₂ (iodine)", "S₂O₃²⁻ (thiosulfate)",
                "2S₂O₃²⁻ → S₄O₆²⁻ + 2e⁻ (ox: S:+II→+2.5 avg)",
                "I₂ + 2e⁻ → 2I⁻ (red: 0→-I)",
                2
            ),
            "mnso4_kmno4": (
                "KMnO₄ (permanganate)", "MnSO₄ (manganese(II))",
                "Mn²⁺ → MnO₄⁻ (ox: +II→+VII, loses 5e⁻)",
                "MnO₄⁻ → Mn²⁺ (red: +VII→+II, gains 5e⁻)",
                5
            ),
            "al_fe2o3": (
                "Fe₂O₃ (iron(III) oxide)", "Al (aluminum metal)",
                "Al → Al³⁺ + 3e⁻ (ox: 0→+III, loses 3e⁻)",
                "Fe³⁺ + 3e⁻ → Fe (red: +III→0, gains 3e⁻)",
                3
            ),
            "halogen_displacement": None,  # generic pattern handled separately
        }

    def _run_base(self, equation: str) -> dict:
        """Identify oxidizing and reducing agents."""
        eq = equation.replace('→', '=').replace('->', '=').replace('−>', '=')
        sides = eq.split('=')
        if len(sides) != 2:
            raise ChemMCPError(f"Invalid equation format: '{equation}'. Use '=' to separate reactants and products.")

        import re as _re
        reactants = [r.strip() for r in _re.split(r'\s*\+\s*(?=[A-Z])', sides[0]) if r.strip()]
        products = [p.strip() for p in _re.split(r'\s*\+\s*(?=[A-Z])', sides[1]) if p.strip()]

        eq_lower = equation.lower()

        # Try pattern matching first
        result = self._match_pattern(eq_lower, reactants, products)
        if result:
            return result

        # Rule-based analysis fallback
        return self._rule_based_analysis(equation, reactants, products)

    def _run_text(self, equation: str) -> dict:
        return self._run_base(equation)

    def _match_pattern(self, eq_lower: str, reactants: list, products: list) -> Optional[dict]:
        """Try to match against known redox patterns."""
        r_str = " ".join(r.lower().strip() for r in reactants)
        p_str = " ".join(p.lower().strip() for p in products)
        combined = f"{r_str} {p_str}"

        # Check each pattern
        for pat_key, pat_data in self._redox_patterns.items():
            if pat_data is None:
                continue

            # Build search tokens from pattern key
            pat_tokens = set(pat_key.split("_"))

            # Check if most pattern tokens appear in the equation
            matches = sum(1 for t in pat_tokens if t in combined)
            if matches >= len(pat_tokens) * 0.7:
                ox_agent, red_agent, ox_desc, red_desc, n_e = pat_data

                # Determine reaction type detail
                type_detail = self._classify_redox_type(pat_key, ox_agent, red_agent)

                logger.info(f"Pattern matched: {pat_key}")
                return {
                    "oxidizing_agent": ox_agent,
                    "reducing_agent": red_agent,
                    "oxidation_process": ox_desc,
                    "reduction_process": red_desc,
                    "electrons_transferred": n_e,
                    "is_redox": True,
                    "reaction_type_detail": type_detail,
                }

        return None

    def _classify_redox_type(self, pat_key: str, ox: str, red: str) -> str:
        """Classify the specific type of redox process."""
        classifications = {
            "mno4_fe": "Redox — Permanganate titration (acidic medium)",
            "cr2o7_fe": "Redox — Dichromate titration (acidic medium)",
            "zn_cuso4": "Redox — Galvanic cell-type single displacement (metal activity series)",
            "fe_cuso4": "Redox — Metal displacement (Fe above Cu in activity series)",
            "zn_hcl": "Redox — Metal-acid reaction (H₂ evolution)",
            "fe_hcl": "Redox — Metal-acid reaction (H₂ evolution, forms Fe²⁺ with dilute acid)",
            "cl2_nai": "Redox — Halogen displacement (Cl₂ more electronegative than I₂)",
            "cl2_nabr": "Redox — Halogen displacement (Cl₂ more electronegative than Br₂)",
            "br2_nai": "Redox — Halogen displacement (Br₂ more electronegative than I₂)",
            "cu_hno3_dilute": "Redox — Metal oxidation by nitric acid (dilute, produces NO)",
            "cu_hno3_conc": "Redox — Metal oxidation by nitric acid (concentrated, produces NO₂)",
            "c_h2so4_conc": "Redox — Hot concentrated H₂SO₄ oxidation of nonmetal",
            "h2_cuo": "Redox — Hydrogen reduction of metal oxide (smelting/reduction)",
            "co_fe2o3": "Redox — Blast furnace iron smelting (CO reduces iron oxide)",
            "kclo3_decomp": "Redox — Thermal decomposition with disproportionation (chlorine: +V → -I and 0)",
            "h2o2_decomp": "Redox — Disproportionation of hydrogen peroxide",
            "na2o2_co2": "Redox — Peroxide disproportionation",
            "s2o32-_i2": "Redox — Iodometric titration endpoint (thiosulfate reduces iodine)",
            "mnso4_kmno4": "Redox — Comproportionation of manganese species",
            "al_fe2o3": "Redox — Thermite reaction (aluminothermic reduction, highly exothermic)",
        }
        return classifications.get(pat_key, f"Redox — Oxidized by {ox}, reduced by {red}")

    def _rule_based_analysis(self, equation: str, reactants: list, products: list) -> dict:
        """Fallback rule-based analysis when no pattern matches."""
        # Check for elemental forms being converted to/from compounds
        elem_pattern = r'^[A-Z][a-z]?$'

        ox_candidates = []  # potential oxidizing agents (get reduced)
        red_candidates = []  # potential reducing agents (get oxidized)

        for r in reactants:
            rs = r.strip()
            # Free element on reactant side → likely oxidized (reducing agent)
            if re.match(elem_pattern, rs):
                red_candidates.append(f"{rs} (elemental form)")
            # Check for common reducing agents
            rl = rs.lower()
            if any(x in rl for x in ["fe2+", "sn2+", "i-", "br-", "s2o3", "so2", "h2", "co", "c ", "al", "zn", "mg"]):
                red_candidates.append(rs)

        for p in products:
            ps = p.strip()
            # Free element on product side → likely reduced (product of reduction)
            if re.match(elem_pattern, ps):
                ox_candidates.append("species producing " + ps)
            pl = ps.lower()
            if any(x in pl for x in ["mn2+", "fe3+", "fe2+", "cu", "ag", "hg", "h2", "no", "no2", "so2", "i2", "br2", "cl2"]):
                ox_candidates.append(ps)

        # Check if it's a non-redox reaction
        # Acid-base indicators
        acid_base_indicators = ["h2o", "nacl", "kcl", "salt", "nano3", "caso4", "bacl2", "agcl"]
        prod_lower = [p.strip().lower() for p in products]
        if any(any(ind in pl for ind in acid_base_indicators) for pl in prod_lower):
            react_lower = [r.strip().lower() for r in reactants]
            has_acid = any(a in rl for rl in react_lower for a in ["hcl", "h2so4", "hno3", "h3po4"])
            has_base = any(b in rl for rl in react_lower for b in ["oh", "naoh", "koh", "ba(oh)2", "ca(oh)2"])
            if has_acid and has_base:
                return {
                    "oxidizing_agent": "N/A — no electron transfer",
                    "reducing_agent": "N/A — no electron transfer",
                    "oxidation_process": "N/A — all oxidation states unchanged",
                    "reduction_process": "N/A — all oxidation states unchanged",
                    "electrons_transferred": 0,
                    "is_redox": False,
                    "reaction_type_detail": "Acid-Base Neutralization — NOT a redox reaction",
                }

        if ox_candidates or red_candidates:
            return {
                "oxidizing_agent": "; ".join(ox_candidates) if ox_candidates else "Unknown (could not determine precisely)",
                "reducing_agent": "; ".join(red_candidates) if red_candidates else "Unknown (could not determine precisely)",
                "oxidation_process": "Detected oxidation (see reducing agent)" if red_candidates else "Uncertain",
                "reduction_process": "Detected reduction (see oxidizing agent)" if ox_candidates else "Uncertain",
                "electrons_transferred": -1,  # unknown exact number
                "is_redox": True,
                "reaction_type_detail": "Redox — identified by heuristic rules (exact electron count requires manual verification)",
            }

        raise ChemMCPError(
            f"Cannot confidently identify oxidizing/reducing agents in '{equation}'. "
            f"This tool recognizes common redox patterns including permanganate, dichromate, "
            f"metal displacement, halogen displacement, nitric acid oxidation, thermite, "
            f"combustion, decomposition/disproportionation reactions, etc."
        )
