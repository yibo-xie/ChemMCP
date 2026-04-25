import logging
import re
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IdentifyReactionType(BaseTool):
    """
    Identify the type of chemical reaction (combination, decomposition, single replacement,
    double displacement, redox, combustion, acid-base, precipitation, etc.).
    """
    __version__ = "0.1.0"
    name = "IdentifyReactionType"
    func_name = "identify_reaction_type"
    description = "Identify the type(s) of a chemical reaction from its equation string."
    implementation_description = "Uses pattern matching against known reaction type characteristics: reactant/product count patterns, oxidation state changes, acid-base pairs, precipitate formation rules, combustion patterns, etc."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Reaction Classification", "Reaction Types", "Chemical Analysis", "Pattern Matching"]
    required_envs = []

    code_input_sig = [
        ("equation", "str", "N/A", "Chemical equation string (balanced or unbalanced), e.g., 'MnO4- + 5Fe2+ + 8H+ = Mn2+ + 5Fe3+ + 4H2O' or 'H2+O2=H2O'."),
    ]

    text_input_sig = [
        ("equation", "str", "N/A", "Chemical equation string."),
    ]

    output_sig = [
        ("reaction_type", "str", "Primary reaction type (e.g., 'Redox', 'Combination', 'Acid-Base')."),
        ("sub_type", "str", "Sub-type if applicable (e.g., 'Single Displacement', 'Oxidation-Reduction (Disproportionation)')."),
        ("description", "str", "Brief explanation of why this reaction is classified as such."),
        ("oxidation_info", "str", "If redox: describes what is oxidized and reduced (N/A for non-redox)."),
    ]

    examples = [
        {
            "code_input": {"equation": "MnO4- + 5Fe2+ + 8H+ = Mn2+ + 5Fe3+ + 4H2O"},
            "text_input": {"equation": "MnO4- + 5Fe2+ + 8H+ = Mn2+ + 5Fe3+ + 4H2O"},
            "output": {
                "reaction_type": "Redox",
                "sub_type": "Oxidation-Reduction (Electron Transfer)",
                "description": "This is a redox reaction where Fe²⁺ is oxidized to Fe³⁺ (loses 1e⁻ per atom) and MnO₄⁻ is reduced to Mn²⁺ (gains 5e⁻). Permanganate acts as oxidizing agent in acidic medium.",
                "oxidation_info": "Oxidation: Fe²⁺ → Fe³⁺ + e⁻ | Reduction: MnO₄⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O",
            }
        },
        {
            "code_input": {"equation": "H2SO4 + 2NaOH = Na2SO4 + 2H2O"},
            "text_input": {"equation": "H2SO4 + NaOH = Na2SO4 + H2O"},
            "output": {
                "reaction_type": "Acid-Base",
                "sub_type": "Neutralization (Double Displacement)",
                "description": "A classic acid-base neutralization reaction where H₂SO₄ (acid) reacts with NaOH (base) to form salt (Na₂SO₄) and water.",
                "oxidation_info": "N/A — no electron transfer occurs; all oxidation states remain unchanged.",
            }
        },
        {
            "code_input": {"equation": "CH4 + 2O2 = CO2 + 2H2O"},
            "text_input": {"equation": "CH4 + 2O2 = CO2 + 2H2O"},
            "output": {
                "reaction_type": "Combustion / Redox",
                "sub_type": "Hydrocarbon Combustion (Complete)",
                "description": "Complete combustion of methane: carbon is oxidized from -4 to +4, oxygen is reduced from 0 to -2. Highly exothermic.",
                "oxidation_info": "Oxidation: C(-IV) → C(IV) + 8e⁻ | Reduction: O₂(0) → O²(-II) + 2e⁻ per O atom",
            }
        },
        {
            "code_input": {"equation": "Zn + CuSO4 = ZnSO4 + Cu"},
            "text_input": {"equation": "Zn + CuSO4 = ZnSO4 + Cu"},
            "output": {
                "reaction_type": "Redox / Single Replacement",
                "sub_type": "Single Displacement (Metal Displacement)",
                "description": "Zinc displaces copper from copper sulfate because Zn is more reactive than Cu (higher activity series position). Zn is oxidized, Cu²⁺ is reduced.",
                "oxidation_info": "Oxidation: Zn → Zn²⁺ + 2e⁻ | Reduction: Cu²⁺ + 2e⁻ → Cu",
            }
        },
        {
            "code_input": {"equation": "AgNO3 + NaCl = AgCl + NaNO3"},
            "text_input": {"equation": "AgNO3 + NaCl = AgCl + NaNO3"},
            "output": {
                "reaction_type": "Double Displacement / Precipitation",
                "sub_type": "Precipitation Reaction",
                "description": "Ag⁺ ions combine with Cl⁻ ions to form insoluble silver chloride (AgCl) precipitate. A classic metathesis/precipitation reaction.",
                "oxidation_info": "N/A — ion exchange without change in oxidation states.",
            }
        },
        {
            "code_input": {"equation": "2KClO3 = 2KCl + 3O2"},
            "text_input": {"equation": "2KClO3 = 2KCl + O2"},
            "output": {
                "reaction_type": "Decomposition / Redox",
                "sub_type": "Thermal Decomposition (Disproportionation-like / Redox)",
                "description": "Potassium chlorate decomposes upon heating to produce potassium chloride and oxygen gas. Chlorine undergoes disproportionation (+V → -I and 0).",
                "oxidation_info": "Oxidation: 2Cl(+V) → Cl₂(0) + 10e⁻ (in O2) | Reduction: Cl(+V) + 6e⁻ → Cl(-I) (in KCl)",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        # Common oxidation state reference (for redox detection)
        self._ox_states = {
            # Common fixed states
            "H": [1, -1],  # +1 with nonmetals, -1 with metals
            "O": [-2, -1, 0],  # -2 usually, -1 in peroxides, 0 in O2
            "F": [-1],
            "Cl": [-1, 1, 3, 5, 7],  # varies widely
            "Br": [-1],
            "I": [-1],
            "Na": [1], "K": [1], "Li": [1], "Rb": [1], "Cs": [1],
            "Ca": [2], "Mg": [2], "Ba": [2], "Sr": [2], "Zn": [2],
            "Al": [3],
            "Ag": [1],
        }

        # Solubility rules (simplified)
        self._insoluble = frozenset([
            "agcl", "agbr", "agi", "ag2co3", "ag2so4",  # Ag salts
            "baso4", "caso4(partial)", "pbso4",  # sulfates
            "caco3", "baco3", "srco3",  # carbonates
            "cu(oh)2", "fe(oh)2", "fe(oh)3", "mg(oh)2", "al(oh)3",  # hydroxides
        ])

    def _run_base(self, equation: str) -> dict:
        """Identify reaction type from equation string."""
        eq = equation.replace('→', '=').replace('->', '=')
        sides = eq.split('=')
        if len(sides) != 2:
            raise ChemMCPError(f"Invalid equation format: '{equation}'")

        # Smart split that handles charged species (Fe2+, H+, etc.)
        import re as _re
        reactants = _re.split(r'\s*\+\s*(?=[A-Z])', sides[0])
        reactants = [r.strip() for r in reactants if r.strip()]
        products = _re.split(r'\s*\+\s*(?=[A-Z])', sides[1])
        products = [p.strip() for p in products if p.strip()]

        r_count = len(reactants)
        p_count = len(products)
        all_species = reactants + products

        # Build analysis
        results = []
        sub_types = []
        ox_info = "N/A"

        # ── Pattern-based classification ──

        eq_lower = equation.lower()

        # 1. Combustion check: hydrocarbon/Organic + O2 → CO2 + H2O
        has_o2_reactant = any("o2" in r.lower().strip() and len(r.strip()) <= 3 for r in reactants)
        has_co2_prod = any("co2" in p.lower().strip() for p in products)
        has_h2o_prod = any(re.match(r'^\d*H2O$', p.strip(), re.IGNORECASE) for p in products)
        organic_reactant = any(
            re.search(r'[Cc][Hh]', r) or
            any(c.isalpha() and c not in 'HCOoNnSsPpFfClBriKkNaamMggaAlLizZnnFFeeCCuu'
                for c in r if c.isalpha())
            for r in reactants
        )

        if has_o2_reactant and (has_co2_prod or has_h2o_prod):
            results.append("Combustion")
            sub_types.append("Hydrocarbon Combustion (Complete)" if has_co2_prod else "Combustion")

        # 2. Acid-Base check
        acids = []
        bases = []
        for r in reactants:
            rl = r.lower().strip()
            if any(a in rl for a in ["h2so4", "hcl", "hno3", "h3po4", "h2co3", "hclo4", "hf", "hbr", "hi", "h2s", "ch3coOH".lower()]):
                acids.append(r)
            if any(b in rl for b in ["naoh", "koh", "ba(oh)2", "ca(oh)2", "nh3", "nh4oh"]):
                bases.append(r)
        has_water_prod = any(re.match(r'^\d*H2O$', p.strip(), re.IGNORECASE) for p in products)
        has_salt_prod = any(
            any(sym in p.lower() for sym in ["no3", "so4", "cl", "co3", "po4"])
            and not any(e in p.lower() for e in ["h+", "oh-", "h2"])
            for p in products
        )
        if acids and bases:
            results.append("Acid-Base")
            sub_types.append("Neutralization")
            if r_count == 2 and p_count == 2:
                sub_types.append("Double Displacement")

        # 3. Precipitation check
        if r_count >= 2 and p_count >= 2 and not acids:
            prod_lower = [p.lower().strip() for p in products]
            for insol in self._insoluble:
                if any(insol in pl for pl in prod_lower):
                    if "Precipitation" not in results:
                        results.insert(0, "Precipitation")
                        sub_types.append("Precipitation Reaction")
                    break

        # 4. Combination (synthesis): A + B → AB
        if r_count >= 2 and p_count == 1 and not results:
            results.append("Combination")
            sub_types.append("Synthesis / Combination")

        # 5. Decomposition: AB → A + B
        if r_count == 1 and p_count >= 2 and not results:
            results.append("Decomposition")
            sub_types.append("Decomposition")

        # 6. Single replacement: A + BC → AC + B (element + compound → new compound + element)
        if r_count == 2 and p_count == 2:
            is_single_repl = False
            # Check if one reactant is a pure element and one product is a pure element
            pure_elem_react = None
            pure_elem_prod = None
            for r in reactants:
                if re.match(r'^[A-Z][a-z]?$', r.strip()):
                    pure_elem_react = r.strip()
            for p in products:
                if re.match(r'^[A-Z][a-z]?$', p.strip()):
                    pure_elem_prod = p.strip()
            # Also check: element from pure reactant appears in a product compound
            # and element from pure product appeared in a reactant compound
            if pure_elem_react and pure_elem_prod and pure_elem_react != pure_elem_prod:
                # Verify: pure reactant element is in a product compound
                in_prod_compound = any(pure_elem_react in p for p in products if not re.match(r'^[A-Z][a-z]?$', p.strip()))
                # And pure product element was in a reactant compound
                in_react_compound = any(pure_elem_prod in r for r in reactants if not re.match(r'^[A-Z][a-z]?$', r.strip()))
                if in_prod_compound and in_react_compound:
                    is_single_repl = True

            # Fallback: original element uniqueness check
            if not is_single_repl:
                react_elems = set()
                for r in reactants:
                    elems = re.findall(r'[A-Z][a-z]?', r)
                    react_elems.update(elems)
                prod_elems = set()
                for p in products:
                    elems = re.findall(r'[A-Z][a-z]?', p)
                    prod_elems.update(elems)
                unique_r = react_elems - prod_elems
                unique_p = prod_elems - react_elems
                if len(unique_r) == 1 and len(unique_p) == 1:
                    is_single_repl = True

            if is_single_repl:
                if "Single Replacement" not in sub_types:
                    results.append("Single Replacement")
                    sub_types.append("Single Displacement (Metal Displacement)")

        # 7. Double displacement: AB + CD → AD + CB
        if r_count == 2 and p_count == 2 and "Single Replacement" not in str(sub_types) and "Acid-Base" not in results:
            if not results:
                results.append("Double Displacement")
                sub_types.append("Metathesis / Double Displacement")

        # 8. Redox detection via oxidation state changes (check after specific types)
        redox_detected, ox_detail = self._detect_redox(reactants, products)
        if redox_detected:
            if not any("Redox" in r for r in results):
                results.append("Redox / " + results[0] if results else "Redox")
            ox_info = ox_detail

        # Determine primary type
        primary = results[0] if results else "Unknown"
        subtype_str = "; ".join(dict.fromkeys(sub_types)) if sub_types else "General"

        # Generate description
        desc = self._generate_description(primary, subtype_str, reactants, products, ox_info)

        return {
            "reaction_type": primary,
            "sub_type": subtype_str,
            "description": desc,
            "oxidation_info": ox_info,
        }

    def _run_text(self, equation: str) -> dict:
        return self._run_base(equation)

    def _detect_redox(self, reactants: list, products: list) -> tuple:
        """Detect if oxidation states change between reactants and products."""
        # Simple heuristic: look for elemental forms (oxidation number 0) becoming compounds
        # or vice versa, plus known redox species
        redox_indicators = [
            ("mno4-", "mn2+"), ("cr2o72-", "cr3+"), ("fe2+", "fe3+"), ("fe3+", "fe2+"),
            ("sn2+", "sn4+"), ("cu", "cu2+"), ("cu2+", "cu"), ("zn", "zn2+"),
            ("zn2+", "zn"), ("al", "al3+"), ("na", "na+"), ("mg", "mg2+"),
            ("o2", None), ("cl2", None), ("h2", None), ("f2", None),
            ("no3-", "no"), ("no3-", "no2"), ("so4^2-", "so2"), ("so2", "so4^2-"),
            ("clo-", "cl-"), ("cl2", "cl-"), ("i2", "i-"),
            ("c2o42-", "co2"), ("h2o2", "h2o"), ("h2o2", "o2"),
            ("kclo3", "kcl"), ("kmno4", "mn2+"), ("k2cr2o7", "cr3+"),
        ]

        r_str = " ".join(r.lower() for r in reactants)
        p_str = " ".join(p.lower() for p in products)

        ox_changes = []
        for ox_sp, red_sp in redox_indicators:
            if ox_sp in r_str:
                if red_sp and red_sp in p_str:
                    ox_changes.append(f"{ox_sp} → {red_sp}")
                elif red_sp is None:
                    # Element being consumed (reduced or oxidized)
                    ox_changes.append(f"{ox_sp} (elemental) → compound")

        if ox_changes:
            detail = f"Oxidation/reduction detected: {' | '.join(ox_changes)}"
            return True, detail

        # Check for free elements
        elem_pattern = r'^[A-Z][a-z]?$'
        for r in reactants:
            if re.match(elem_pattern, r.strip()):
                for p in products:
                    if not re.match(elem_pattern, p.strip()) and len(p.strip()) > 2:
                        return True, f"Element {r} is oxidized/reduced to form compound"

        return False, "N/A"

    def _generate_description(self, primary: str, sub_type: str,
                               reactants: list, products: list, ox_info: str) -> str:
        """Generate human-readable description of the classification."""
        r_str = " + ".join(reactants)
        p_str = " + ".join(products)

        descs = {
            "Redox": f"This is a redox (oxidation-reduction) reaction involving electron transfer between species.",
            "Combustion": f"This is a combustion reaction where a fuel reacts rapidly with oxygen, releasing energy as heat and light.",
            "Acid-Base": f"This is an acid-base neutralization reaction producing salt and water.",
            "Precipitation": f"This is a precipitation reaction where an insoluble solid (precipitate) forms from aqueous solutions.",
            "Combination": f"This is a combination (synthesis) reaction where two or more substances combine to form a single product.",
            "Decomposition": f"This is a decomposition reaction where a single compound breaks down into two or simpler substances.",
            "Single Replacement": f"This is a single replacement (displacement) reaction where a more reactive element replaces a less reactive one in a compound.",
            "Double Displacement": f"This is a double displacement (metathesis) reaction where cations and anions exchange partners.",
        }

        base_desc = descs.get(primary, f"This reaction is classified as '{primary}'.")
        detail = f" Reaction: {r_str} → {p_str}"

        return base_desc + detail
