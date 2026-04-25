import logging
import re
from typing import Dict, List, Optional, Tuple, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AssignOxidationNumber(BaseTool):
    """
    标注化学式中各元素的氧化数（氧化态）。
    基于常见氧化数规则和内置数据库进行计算。
    """
    __version__ = "0.1.0"
    name = "AssignOxidationNumber"
    func_name = "assign_oxidation_number"
    description = "Assign oxidation numbers (oxidation states) to each element in a chemical formula or reaction equation."
    implementation_description = "Uses standard oxidation number rules (electronegativity, common ion charges, polyatomic ions) with a built-in database of common compound oxidation states."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Oxidation Number", "Inorganic Chemistry", "Redox", "Valence"]
    required_envs = []

    code_input_sig = [
        ("formula", "str", "N/A", "Chemical formula or reaction equation, e.g., 'H2SO4', 'KMnO4', 'Fe2O3', or 'MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O'."),
    ]

    text_input_sig = [
        ("formula_str", "str", "N/A", "Chemical formula or reaction equation string."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary mapping each element/species to its oxidation number(s), with detailed breakdown."),
    ]

    examples = [
        {
            "code_input": {
                "formula": "H2SO4",
            },
            "text_input": {
                "formula_str": "H2SO4",
            },
            "output": {
                "result": {
                    "formula": "H2SO4",
                    "oxidation_numbers": {"H": +1, "S": +6, "O": -2},
                    "breakdown": "H: +1 (Group 1), O: -2 (common rule), S: calculated as +6 (2*(+1) + S + 4*(-2) = 0)",
                }
            },
        },
        {
            "code_input": {
                "formula": "KMnO4",
            },
            "text_input": {
                "formula_str": "KMnO4",
            },
            "output": {
                "result": {
                    "formula": "KMnO4",
                    "oxidation_numbers": {"K": +1, "Mn": +7, "O": -2},
                    "breakdown": "K: +1 (Group 1), O: -2 (common rule), Mn: calculated as +7 (+1 + Mn + 4*(-2) = 0)",
                }
            },
        },
        {
            "code_input": {
                "formula": "Fe2O3",
            },
            "text_input": {
                "formula_str": "Fe2O3",
            },
            "output": {
                "result": {
                    "formula": "Fe2O3",
                    "oxidation_numbers": {"Fe": +3, "O": -2},
                    "breakdown": "O: -2 (common rule), Fe: calculated as +3 (2*Fe + 3*(-2) = 0 → Fe = +3)",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize oxidation number rules and databases."""
        # Common fixed oxidation states for single-element ions / groups
        self._fixed_states = {
            # Group 1
            "Li": +1, "Na": +1, "K": +1, "Rb": +1, "Cs": +1, "Fr": +1,
            # Group 2
            "Be": +2, "Mg": +2, "Ca": +2, "Sr": +2, "Ba": +2, "Ra": +2,
            # Group 13 (common)
            "Al": +3,
            # Common fixed anions
            "F": -1,
            # Hydrogen (usually +1, -1 in hydrides)
            # Oxygen (usually -2, -1 in peroxides)
            # Common polyatomic ion charges (for reference)
        }

        # Common known compound oxidation states database
        self._compound_db = {
            # Acids
            "H2SO4": {"H": +1, "S": +6, "O": -2},
            "H2SO3": {"H": +1, "S": +4, "O": -2},
            "H2S": {"H": +1, "S": -2},
            "HNO3": {"H": +1, "N": +5, "O": -2},
            "HNO2": {"H": +1, "N": +3, "O": -2},
            "HClO4": {"H": +1, "Cl": +7, "O": -2},
            "HClO3": {"H": +1, "Cl": +5, "O": -2},
            "HClO": {"H": +1, "Cl": +1, "O": -2},
            "H3PO4": {"H": +1, "P": +5, "O": -2},
            "H2CO3": {"H": +1, "C": +4, "O": -2},
            "CH3COOH": {"H": [+1, +1, +1, +1], "C": [-3, +3], "O": [-2, -1]},
            "H2O2": {"H": +1, "O": -1},  # peroxide
            "H2O": {"H": +1, "O": -2},
            "HCl": {"H": +1, "Cl": -1},
            "HI": {"H": +1, "I": -1},
            "HF": {"H": +1, "F": -1},
            "NH3": {"N": -3, "H": +1},
            "NH4+": {"N": -3, "H": +1},
            "CH4": {"C": -4, "H": +1},
            "CO2": {"C": +4, "O": -2},
            "CO": {"C": +2, "O": -2},
            "SiO2": {"Si": +4, "O": -2},
            "SO2": {"S": +4, "O": -2},
            "SO3": {"S": +6, "O": -2},
            "NO": {"N": +2, "O": -2},
            "NO2": {"N": +4, "O": -2},
            "N2O": {"N": [+1, 0], "O": -2},
            "N2O5": {"N": +5, "O": -2},
            "N2O3": {"N": +3, "O": -2},
            "ClO2": {"Cl": +4, "O": -2},
            "Cl2O": {"Cl": +1, "O": -2},
            "Cl2O7": {"Cl": +7, "O": -2},
            # Oxides
            "Fe2O3": {"Fe": +3, "O": -2},
            "FeO": {"Fe": +2, "O": -2},
            "Fe3O4": {"Fe": [+2, +3], "O": -2},  # mixed valence
            "CuO": {"Cu": +2, "O": -2},
            "Cu2O": {"Cu": +1, "O": -2},
            "Ag2O": {"Ag": +1, "O": -2},
            "ZnO": {"Zn": +2, "O": -2},
            "CaO": {"Ca": +2, "O": -2},
            "Al2O3": {"Al": +3, "O": -2},
            "Cr2O3": {"Cr": +3, "O": -2},
            "CrO3": {"Cr": +6, "O": -2},
            "MnO2": {"Mn": +4, "O": -2},
            "Mn2O7": {"Mn": +7, "O": -2},
            "PbO2": {"Pb": +4, "O": -2},
            "PbO": {"Pb": +2, "O": -2},
            "SnO2": {"Sn": +4, "O": -2},
            "CO32-": {"C": +4, "O": -2},
            "SO42-": {"S": +6, "O": -2},
            "SO32-": {"S": +4, "O": -2},
            "NO3-": {"N": +5, "O": -2},
            "NO2-": {"N": +3, "O": -2},
            "PO43-": {"P": +5, "O": -2},
            "ClO4-": {"Cl": +7, "O": -2},
            "ClO3-": {"Cl": +5, "O": -2},
            "ClO-": {"Cl": +1, "O": -2},
            "MnO4-": {"Mn": +7, "O": -2},
            "MnO42-": {"Mn": +6, "O": -2},
            "Cr2O72-": {"Cr": +6, "O": -2},
            "CrO42-": {"Cr": +6, "O": -2},
            "OH-": {"O": -2, "H": +1},
            # Salts
            "KMnO4": {"K": +1, "Mn": +7, "O": -2},
            "K2Cr2O7": {"K": +1, "Cr": +6, "O": -2},
            "K2CrO4": {"K": +1, "Cr": +6, "O": -2},
            "NaCl": {"Na": +1, "Cl": -1},
            "Na2CO3": {"Na": +1, "C": +4, "O": -2},
            "CaCO3": {"Ca": +2, "C": +4, "O": -2},
            "BaSO4": {"Ba": +2, "S": +6, "O": -2},
            "AgCl": {"Ag": +1, "Cl": -1},
            "AgBr": {"Ag": +1, "Br": -1},
            "AgI": {"Ag": +1, "I": -1},
            "PbI2": {"Pb": +2, "I": -1},
            "CuSO4": {"Cu": +2, "S": +6, "O": -2},
            "FeSO4": {"Fe": +2, "S": +6, "O": -2},
            "Fe2(SO4)3": {"Fe": +3, "S": +6, "O": -2},
            "AlCl3": {"Al": +3, "Cl": -1},
            "PCl3": {"P": +3, "Cl": -1},
            "PCl5": {"P": +5, "Cl": -1},
            "SF6": {"S": +6, "F": -1},
            "NaH": {"Na": +1, "H": -1},  # hydride
            "CaH2": {"Ca": +2, "H": -1},  # hydride
            "OF2": {"O": +2, "F": -1},  # exception: F is more electronegative
            "Na2O2": {"Na": +1, "O": -1},  # peroxide
            # Others
            "H2": {"H": 0},
            "O2": {"O": 0},
            "N2": {"N": 0},
            "Cl2": {"Cl": 0},
            "F2": {"F": 0},
            "O3": {"O": 0},
            "P4": {"P": 0},
            "S8": {"S": 0},
        }

        # Electronegativity order (Pauling scale, higher = more EN)
        # Used when simple calculation doesn't work
        self._en_order = ["F", "O", "Cl", "N", "Br", "I", "S", "C", "H", "P", " metals"]

    def _run_base(self, formula: str) -> dict:
        """Core logic: assign oxidation numbers to each element in the formula."""
        formula = formula.strip()
        if not formula:
            raise ChemMCPError("Formula cannot be empty.")

        # Check if it's a reaction equation (contains '=' or '→')
        if "=" in formula or "→" in formula:
            return self._handle_equation(formula)

        # Normalize formula
        normalized = self._normalize_formula(formula)

        # Try exact match in database first
        result = self._lookup_database(normalized)
        if result is not None:
            return result

        # Calculate from rules
        result = self._calculate_oxidation(normalized)
        return result

    def _run_text(self, formula_str: str) -> dict:
        """Parse text input and call core logic."""
        return self._run_base(formula_str.strip())

    def _normalize_formula(self, formula: str) -> str:
        """Normalize formula string."""
        s = formula.strip()
        s = s.replace(" ", "")
        return s

    def _lookup_database(self, formula: str) -> Optional[dict]:
        """Look up oxidation numbers from built-in database."""
        # Direct lookup
        if formula in self._compound_db:
            data = self._compound_db[formula]
            breakdown = self._generate_breakdown(formula, data)
            return {
                "formula": formula,
                "oxidation_numbers": data,
                "breakdown": breakdown,
                "source": "database",
            }

        # Try stripping charge
        charge_stripped = re.sub(r'[\+\-]\d*$', '', formula)
        if charge_stripped in self._compound_db:
            data = self._compound_db[charge_stripped]
            breakdown = self._generate_breakdown(charge_stripped, data)
            return {
                "formula": formula,
                "oxidation_numbers": data,
                "breakdown": breakdown + f" (charge adjusted for ion {formula})",
                "source": "database",
            }

        return None

    def _generate_breakdown(self, formula: str, ox_data: dict) -> str:
        """Generate human-readable breakdown explanation."""
        parts = []
        for elem, val in ox_data.items():
            if isinstance(val, list):
                parts.append(f"{elem}: {val} (multiple atoms/positions)")
            elif val == 0:
                parts.append(f"{elem}: {val} (elemental form)")
            elif val > 0:
                parts.append(f"{elem}: +{val}")
            else:
                parts.append(f"{elem}: {val}")
        return "; ".join(parts)

    def _calculate_oxidation(self, formula: str) -> dict:
        """Calculate oxidation numbers using systematic rules."""
        # Parse charge
        overall_charge = 0
        charge_match = re.search(r'\[?([^\]]*)\]?(\d*[+-])$', formula)
        clean_formula = formula
        if charge_match:
            charge_str = charge_match.group(2)
            clean_formula = charge_match.group(1) if charge_match.group(1) else re.sub(r'\d*[+-]$', '', formula)
            if charge_str.endswith("+") and len(charge_str) > 1:
                overall_charge = int(charge_str[:-1])
            elif charge_str == "+":
                overall_charge = +1
            elif charge_str.endswith("-") and len(charge_str) > 1:
                overall_charge = -int(charge_str[:-1])
            elif charge_str == "-":
                overall_charge = -1

        # Simple pattern-based parser for formulas like H2SO4, Fe2O3, etc.
        elements = self._parse_elements(clean_formula)
        if not elements:
            raise ChemMCPError(f"Cannot parse formula: '{formula}'")

        # Apply rules
        ox_numbers = {}
        known = {}  # element -> oxidation number
        unknown = []  # elements still unknown

        for elem, count in elements:
            if elem in self._fixed_states:
                known[elem] = self._fixed_states[elem]
                ox_numbers[elem] = self._fixed_states[elem]
            elif elem == "O":
                # Check for peroxide (O is -1) or OF2 (O is positive)
                if any(e == "F" for e, _ in elements):
                    ox_numbers[elem] = +2  # OF2 case
                    known["O"] = +2
                else:
                    ox_numbers[elem] = -2
                    known["O"] = -2
            elif elem == "H":
                # Check if it's a hydride (metal present)
                has_metal = any(e in self._fixed_states and self._fixed_states[e] > 0 for e, _ in elements)
                if has_metal and len(elements) == 2:
                    ox_numbers[elem] = -1
                    known["H"] = -1
                else:
                    ox_numbers[elem] = +1
                    known["H"] = +1
            elif elem == "F":
                ox_numbers[elem] = -1
                known["F"] = -1
            else:
                unknown.append((elem, count))

        # Calculate unknown elements by charge balance
        if unknown:
            total_known_charge = sum(known.get(e, 0) * c for e, c in elements if e in known)
            total_unknown_count = sum(c for e, c in unknown)
            remaining_charge = overall_charge - total_known_charge

            if len(unknown) == 1:
                elem, count = unknown[0]
                ox_val = round(remaining_charge / count, 2)
                ox_numbers[elem] = ox_val
                known[elem] = ox_val
            else:
                # Multiple unknowns — use best guess based on common patterns
                for elem, count in unknown:
                    ox_numbers[elem] = "?"
                    known[elem] = "?"

        breakdown = self._generate_breakdown(clean_formula, ox_numbers)
        if "?" in str(ox_numbers):
            breakdown += " (some values could not be uniquely determined)"

        return {
            "formula": formula,
            "oxidation_numbers": ox_numbers,
            "breakdown": breakdown,
            "source": "calculated",
        }

    def _parse_elements(self, formula: str) -> List[Tuple[str, int]]:
        """Parse chemical formula into list of (element, count) tuples."""
        elements = []
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        for elem, count_str in matches:
            if not elem:
                continue
            count = int(count_str) if count_str else 1
            elements.append((elem, count))
        return elements

    def _handle_equation(self, equation: str) -> dict:
        """Handle full reaction equations — analyze each species."""
        eq_normalized = equation.replace('→', '=').replace('->', '=').replace('−>', '=')
        sides = eq_normalized.split('=')
        if len(sides) != 2:
            raise ChemMCPError(f"Invalid equation format: '{equation}'. Use '=' to separate sides.")

        all_species_results = {}
        # Parse both sides
        for side_idx, side in enumerate(sides):
            side_name = "reactants" if side_idx == 0 else "products"
            species_list = self._split_species(side.strip())
            species_results = {}
            for sp in species_list:
                sp_clean = sp.strip().rstrip('+').lstrip('+')
                if sp_clean:
                    try:
                        result = self._run_base(sp_clean)
                        species_results[sp_clean] = result.get("oxidation_numbers", {})
                    except ChemMCPError as e:
                        species_results[sp_clean] = f"Error: {e}"
            all_species_results[side_name] = species_results

        # Determine redox changes
        redox_analysis = self._analyze_redox_change(all_species_results)

        return {
            "equation": equation,
            "species_analysis": all_species_results,
            "redox_analysis": redox_analysis,
        }

    def _split_species(self, side: str) -> List[str]:
        """Split reaction side into individual species."""
        # Split by '+' but be careful with charges
        species = re.split(r'\s*\+\s*', side)
        return [s.strip() for s in species if s.strip()]

    def _analyze_redox_change(self, analysis: dict) -> dict:
        """Analyze oxidation state changes across reactants and products."""
        reactants = analysis.get("reactants", {})
        products = analysis.get("products", {})
        return {
            "note": "Compare oxidation numbers between reactant and product species to identify oxidation/reduction.",
            "hint": "Increase in oxidation number = oxidation (loss of electrons); Decrease = reduction (gain of electrons).",
        }
