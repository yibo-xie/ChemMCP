import logging
import re
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DisproportionationCheck(BaseTool):
    """
    判断一个反应是否为歧化反应（或归中反应/反歧化）。
    歧化反应：同一元素既被氧化又被还原（如 Cl2 + H2O → HCl + HClO）
    归中反应：同一元素从不同氧化态变为中间氧化态
    """
    __version__ = "0.1.0"
    name = "DisproportionationCheck"
    func_name = "check_disproportionation"
    description = "Check if a chemical reaction is a disproportionation reaction (same element is both oxidized and reduced), comproportionation reaction, or neither."
    implementation_description = "Parses the reaction equation, identifies oxidation states of each element in reactants and products, and determines if any single element undergoes simultaneous oxidation and reduction."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Disproportionation", "Redox", "Reaction Analysis", "Oxidation State"]
    required_envs = []

    code_input_sig = [
        ("equation", "str", "N/A", "Reaction equation string, e.g., 'Cl2 + H2O = HCl + HClO', '3Cl2 + 6NaOH = 5NaCl + NaClO3 + 3H2O', or '2H2O2 = 2H2O + O2'."),
    ]

    text_input_sig = [
        ("equation_str", "str", "N/A", "Reaction equation string."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing is_disproportionation, reaction_type, element_involved, oxidation_changes, explanation."),
    ]

    examples = [
        {
            "code_input": {
                "equation": "Cl2 + H2O = HCl + HClO",
            },
            "text_input": {
                "equation_str": "Cl2 + H2O = HCl + HClO",
            },
            "output": {
                "result": {
                    "is_disproportionation": True,
                    "reaction_type": "disproportionation",
                    "element_involved": "Cl",
                    "oxidation_changes": {"reactant_state": 0, "product_states": [-1, +1]},
                    "explanation": "Chlorine (Cl) starts at oxidation state 0 in Cl2, and is simultaneously reduced to -1 in HCl and oxidized to +1 in HClO. This is a classic disproportionation reaction.",
                }
            },
        },
        {
            "code_input": {
                "equation": "2H2O2 = 2H2O + O2",
            },
            "text_input": {
                "equation_str": "2H2O2 = 2H2O + O2",
            },
            "output": {
                "result": {
                    "is_disproportionation": True,
                    "reaction_type": "disproportionation",
                    "element_involved": "O",
                    "oxidation_changes": {"reactant_state": -1, "product_states": [-2, 0]},
                    "explanation": "Oxygen in H2O2 has oxidation state -1 (peroxide). It is reduced to -2 in H2O and oxidized to 0 in O2. This is the decomposition of hydrogen peroxide — a disproportionation reaction.",
                }
            },
        },
        {
            "code_input": {
                "equation": "2Na + Cl2 = 2NaCl",
            },
            "text_input": {
                "equation_str": "2Na + Cl2 = 2NaCl",
            },
            "output": {
                "result": {
                    "is_disproportionation": False,
                    "reaction_type": "normal_redox",
                    "element_involved": None,
                    "oxidation_changes": {},
                    "explanation": "Sodium is oxidized (0 → +1) and chlorine is reduced (0 → -1). Different elements are oxidized/reduced — this is a normal redox reaction, not disproportionation.",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize known disproportionation reactions database."""
        # Known disproportionation reactions database
        self._known_disproportionations = {
            # Key: normalized species set; Value: analysis result
            "cl2_h2o_hcl_hclo": {
                "type": "disproportionation",
                "element": "Cl",
                "reactant_ox": {("Cl2",): 0},
                "product_ox": {("HCl",): -1, ("HClO",): +1},
                "explanation": "Chlorine (Cl) at 0 is reduced to -1 (in HCl) and oxidized to +1 (in HClO).",
            },
            "cl2_naoh_nacl_naclo3_h2o": {
                "type": "disproportionation",
                "element": "Cl",
                "reactant_ox": {("Cl2",): 0},
                "product_ox": {("NaCl",): -1, ("NaClO3",): +5},
                "explanation": "Chlorine (Cl) at 0 in hot concentrated NaOH is reduced to -1 (in NaCl) and oxidized to +5 (in NaClO3).",
            },
            "h2o2_h2o_o2": {
                "type": "disproportionation",
                "element": "O",
                "reactant_ox": {("H2O2",): -1},
                "product_ox": {("H2O",): -2, ("O2",): 0},
                "explanation": "Oxygen in peroxide (-1) is reduced to -2 (in H2O) and oxidized to 0 (in O2).",
            },
            "na2o2_h2o_naoh_o2": {
                "type": "disproportionation",
                "element": "O",
                "reactant_ox": {("Na2O2",): -1},
                "product_ox": {("NaOH",): -2, ("O2",): 0},
                "explanation": "Sodium peroxide: O at -1 goes to -2 (in NaOH) and 0 (in O2).",
            },
            "so2_h2so4_s_h2o": {
                "type": "disproportionation",
                "element": "S",
                "reactant_ox": {("SO2",): +4},
                "product_ox": {("S",): 0, ("H2SO4",): +6},
                "explanation": "Sulfur(IV) in SO2 is reduced to 0 (elemental S) and oxidized to +6 (in H2SO4).",
            },
            "no2_no_no3_h2o": {
                "type": "disproportionation",
                "element": "N",
                "reactant_ox": {("NO2",): +4},
                "product_ox": {("NO",): +2, ("NO3",): +5},
                "explanation": "Nitrogen dioxide disproportionates in water: N(+4) → N(+2) [reduced] and N(+5) [oxidized].",
            },
            # Known comproportionation (reverse of disproportionation)
            "h2s_so2_s_h2o": {
                "type": "comproportionation",
                "element": "S",
                "reactant_ox": {("H2S",): -2, ("SO2",): +4},
                "product_ox": {("S",): 0},
                "explanation": "Sulfur from -2 (H2S) and +4 (SO2) meets at 0 (elemental S). This is a comproportionation (reverse disproportionation) reaction.",
            },
            "nh3_no_n2_h2o": {
                "type": "comproportionation",
                "element": "N",
                "reactant_ox": {("NH3",): -3, ("NO",): +2},
                "product_ox": {("N2",): 0},
                "explanation": "Nitrogen from -3 (NH3) and +2 (NO) meets at 0 (N2). Comproportionation reaction.",
            },
        }

        # Oxidation state lookup for common species
        self._ox_lookup = {
            # Elements (0)
            "Cl2": 0, "O2": 0, "H2": 0, "N2": 0, "F2": 0, "S": 0,
            # Common compounds
            "HCl": -1, "HBr": -1, "HI": -1,  # halogen state
            "HClO": +1, "HClO2": +3, "HClO3": +5, "HClO4": +7,  # Cl states
            "NaCl": -1, "KCl": -1, "AgCl": -1,  # Cl(-I)
            "NaClO": +1, "NaClO2": +3, "NaClO3": +5, "NaClO4": +7,
            "H2O": -2, "H2O2": -1,  # O states
            "NaOH": -2, "Ca(OH)2": -2, "KOH": -2,
            "CO2": +4, "CO": +2, "CH4": -4, "CS2": -2,
            "SO2": +4, "SO3": +6, "H2S": -2, "H2SO4": +6, "H2SO3": +4,
            "NO": +2, "NO2": +4, "N2O": +1, "NO3": +5, "HNO3": +5, "HNO2": +3,
            "NH3": -3, "NH4+": -3,
            "Na2O": -2, "Na2O2": -1, "CaO": -2, "CuO": -2, "Cu2O": -1,
            "FeO": +2, "Fe2O3": +3, "Fe3O4": "+2/+3",
            "KMnO4": +7, "K2MnO4": +6, "MnO2": +4, "MnSO4": +2,
            "K2Cr2O7": +6, "K2CrO4": +6, "Cr2O3": +3, "CrO3": +6,
            "H2O2": -1, "Na2O2": -1,
        }

        # Element-to-track mapping: which element's ox state changes within each species
        self._species_element_map = {
            "HCl": "Cl", "HBr": "Br", "HI": "I",
            "HClO": "Cl", "HClO2": "Cl", "HClO3": "Cl", "HClO4": "Cl",
            "NaCl": "Cl", "KCl": "Cl", "AgCl": "Cl",
            "NaClO": "Cl", "NaClO3": "Cl", "NaClO4": "Cl",
            "H2O": "O", "H2O2": "O", "O2": "O",
            "NaOH": "O", "KOH": "O", "Ca(OH)2": "O",
            "SO2": "S", "SO3": "S", "H2S": "S", "H2SO4": "S", "S": "S",
            "NO": "N", "NO2": "N", "N2O": "N", "NO3": "N", "HNO3": "N",
            "NH3": "N", "NH4+": "N", "N2": "N",
            "Cl2": "Cl", "H2": "H", "F2": "F",
            "Na2O2": "O", "Na2O": "O",
            "KMnO4": "Mn", "K2MnO4": "Mn", "MnO2": "Mn",
            "K2Cr2O7": "Cr", "K2CrO4": "Cr", "Cr2O3": "Cr",
            "CO2": "C", "CO": "C", "CH4": "C",
            "FeO": "Fe", "Fe2O3": "Fe", "Fe3O4": "Fe",
            "CuO": "Cu", "Cu2O": "Cu",
        }

    def _run_base(self, equation: str) -> dict:
        """Core logic: check if reaction is disproportionation."""
        equation = equation.strip()
        if not equation:
            raise ChemMCPError("Equation cannot be empty.")

        # Normalize
        eq_norm = equation.replace('→', '=').replace('->', '=').replace('−>', '=')
        if "=" not in eq_norm:
            raise ChemMCPError(f"Invalid equation format: '{equation}'. Use '=' to separate reactants and products.")

        sides = eq_norm.split("=")
        reactants = self._parse_species(sides[0])
        products = self._parse_species(sides[1])

        # Try known database first
        db_result = self._check_known_db(reactants, products)
        if db_result:
            return db_result

        # Analyze by oxidation states
        result = self._analyze_disproportionation(reactants, products)
        return result

    def _run_text(self, equation_str: str) -> dict:
        """Parse text input."""
        return self._run_base(equation_str.strip())

    def _parse_species(self, side: str) -> List[str]:
        """Split reaction side into individual species."""
        species = re.split(r'\s*\+\s*', side.strip())
        return [s.strip() for s in species if s.strip()]

    def _normalize_species_name(self, name: str) -> str:
        """Normalize species name for matching."""
        s = name.strip()
        # Strip coefficients
        s = re.sub(r'^\d+', '', s)
        # Strip parenthetical modifiers
        s = re.sub(r'\([^)]*\)$', '', s)
        return s

    def _check_known_db(self, reactants: List[str], products: List[str]) -> Optional[dict]:
        """Check against known disproportionation database."""
        norm_r = [self._normalize_species_name(r).lower() for r in reactants]
        norm_p = [self._normalize_species_name(p).lower() for p in products]
        all_sp = sorted(norm_r + norm_p)
        key = "_".join(all_sp)

        # Direct match
        for db_key, data in self._known_disproportionations.items():
            db_species = db_key.split("_")
            if set(all_sp) == set(db_species):
                return {
                    "is_disproportionation": data["type"] == "disproportionation",
                    "reaction_type": data["type"],
                    "element_involved": data["element"],
                    "oxidation_changes": {
                        "reactant_states": data["reactant_ox"],
                        "product_states": data["product_ox"],
                    },
                    "explanation": data["explanation"],
                    "source": "database",
                }

        return None

    def _analyze_disproportionation(self, reactants: List[str], products: List[str]) -> dict:
        """Analyze oxidation states to determine disproportionation."""
        # Collect elements and their oxidation states on each side
        r_elements = {}  # element -> list of oxidation states found
        p_elements = {}

        for sp in reactants:
            sp_clean = self._normalize_species_name(sp)
            if sp_clean in self._ox_lookup:
                ox = self._ox_lookup[sp_clean]
                elem = self._species_element_map.get(sp_clean, None)
                if elem:
                    r_elements.setdefault(elem, []).append(ox)

        for sp in products:
            sp_clean = self._normalize_species_name(sp)
            if sp_clean in self._ox_lookup:
                ox = self._ox_lookup[sp_clean]
                elem = self._species_element_map.get(sp_clean, None)
                if elem:
                    p_elements.setdefault(elem, []).append(ox)

        # Check for disproportionation: same element appears with different ox states on product side
        # that span both above and below its reactant state
        for elem, r_states in r_elements.items():
            if elem not in p_elements:
                continue
            p_states = p_elements[elem]

            # Reactant should have one consistent state (or elemental form)
            if len(set(r_states)) != 1:
                continue

            r_state = r_states[0]
            # Products should have states both > and < reactant state
            has_higher = any(s > r_state for s in p_states)
            has_lower = any(s < r_state for s in p_states)

            if has_higher and has_lower:
                return {
                    "is_disproportionation": True,
                    "reaction_type": "disproportionation",
                    "element_involved": elem,
                    "oxidation_changes": {
                        "reactant_state": r_state,
                        "product_states": sorted(set(p_states)),
                    },
                    "explanation": (
                        f"{elem} starts at oxidation state {r_state} in reactants, "
                        f"and is simultaneously reduced to {min(p_states)} and oxidized to {max(p_states)} in products. "
                        f"This is a disproportionation reaction."
                    ),
                    "source": "analysis",
                }

        # Check for comproportionation: same element from two different reactant states converges
        for elem, p_states in p_elements.items():
            if elem not in r_elements or len(r_elements.get(elem, [])) < 2:
                continue
            r_states_list = r_elements[elem]
            unique_r = sorted(set(r_states_list))
            if len(unique_r) >= 2:
                p_state = p_states[0] if len(set(p_states)) == 1 else None
                if p_state is not None and min(unique_r) < p_state < max(unique_r):
                    return {
                        "is_disproportionation": False,
                        "reaction_type": "comproportionation",
                        "element_involved": elem,
                        "oxidation_changes": {
                            "reactant_states": unique_r,
                            "product_state": p_state,
                        },
                        "explanation": (
                            f"{elem} from oxidation states {unique_r} in reactants converges to {p_state} in products. "
                            f"This is a comproportionation (reverse disproportionation) reaction."
                        ),
                        "source": "analysis",
                    }

        return {
            "is_disproportionation": False,
            "reaction_type": "normal_redox_or_non_redox",
            "element_involved": None,
            "oxidation_changes": {},
            "explanation": (
                "This reaction does not show disproportionation or comproportionation behavior. "
                "Different elements are likely being oxidized/reduced, or no redox change occurs."
            ),
            "source": "analysis",
        }
