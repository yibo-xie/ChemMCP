"""
从 IUPAC 名称解析配合物结构工具
Parse an IUPAC coordination compound name back to structural information.
"""
import logging
import re
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ComplexFromName(BaseTool):
    """
    从配合物的IUPAC名称反向解析出结构信息：
    金属离子、配体列表、氧化态、电荷、结构式。
    """
    __version__ = "0.1.0"
    name = "ComplexFromName"
    func_name = "complex_from_name"
    description = "Parse an IUPAC coordination compound name to extract metal ion, ligands, oxidation state, charge, and structural formula."
    implementation_description = "Reverse of IUPAC naming: parses ligand prefixes (di-, tri-, tetra-, bis-, tris-...), ligand names, metal name (including -ate → element mapping), and Roman numeral oxidation state."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Nomenclature", "IUPAC", "Parsing"]
    required_envs = []

    code_input_sig = [
        ("iupac_name", "str", "N/A", "The IUPAC name of the coordination compound, e.g., 'pentaamminechlorocobalt(III)' or 'hexacyanoferrate(II)'."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "The IUPAC name string (same as code_input)."),
    ]

    output_sig = [
        ("metal_ion", "str", "Metal element symbol and oxidation state."),
        ("metal_element", "str", "Metal element symbol only."),
        ("oxidation_state", "int", "Oxidation state of the metal."),
        ("ligands", "list", "List of ligand dictionaries with name, count, formula."),
        ("coordination_number", "int", "Total coordination number (sum of denticity × count)."),
        ("complex_charge", "int", "Estimated overall charge of the complex."),
        ("structural_formula", "str", "Structural formula string like [Co(NH3)5Cl]??."),
        ("is_anionic", "bool", "Whether the complex is anionic."),
        ("parse_confidence", "str", "Confidence level of the parse (high/medium/low)."),
    ]

    examples = [
        {
            "code_input": {
                "iupac_name": "pentaamminechlorocobalt(III)",
            },
            "text_input": {
                "query": "pentaamminechlorocobalt(III)"
            },
            "output": {
                "metal_ion": "Co3+",
                "metal_element": "Co",
                "oxidation_state": 3,
                "ligands": [
                    {"name": "ammine", "count": 5, "formula": "NH3"},
                    {"name": "chloro", "count": 1, "formula": "Cl"},
                ],
                "coordination_number": 6,
                "complex_charge": 2,
                "structural_formula": "[Co(NH3)5Cl]2+",
                "is_anionic": False,
                "parse_confidence": "high",
            }
        },
        {
            "code_input": {
                "iupac_name": "hexacyanoferrate(II)",
            },
            "text_input": {
                "query": "hexacyanoferrate(II)"
            },
            "output": {
                "metal_ion": "Fe2+",
                "metal_element": "Fe",
                "oxidation_state": 2,
                "ligands": [
                    {"name": "cyano", "count": 6, "formula": "CN"},
                ],
                "coordination_number": 6,
                "complex_charge": -4,
                "structural_formula": "[Fe(CN)6]4-",
                "is_anionic": True,
                "parse_confidence": "high",
            }
        },
        {
            "code_input": {
                "iupac_name": "tris(ethylenediamine)cobalt(III)",
            },
            "text_input": {
                "query": "tris(ethylenediamine)cobalt(III)"
            },
            "output": {
                "metal_ion": "Co3+",
                "metal_element": "Co",
                "oxidation_state": 3,
                "ligands": [
                    {"name": "ethylenediamine", "count": 3, "formula": "en"},
                ],
                "coordination_number": 6,
                "complex_charge": 3,
                "structural_formula": "[Co(en)3]3+",
                "is_anionic": False,
                "parse_confidence": "high",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize reverse lookup databases."""
        # Ligand pattern → (canonical_name, formula, denticity)
        self._ligand_patterns = [
            (r"aqua|aqu[aá]", ("aqua", "H2O", 1)),
            (r"ammin[ea]", ("ammine", "NH3", 1)),
            (r"carbonyl", ("carbonyl", "CO", 1)),
            (r"nitrosyl", ("nitrosyl", "NO", 1)),
            (r"fluor(?:o|ido)", ("fluoro", "F", 1)),
            (r"chlor(?:o|ido)", ("chloro", "Cl", 1)),
            (r"brom(?:o|ido)", ("bromo", "Br", 1)),
            (r"iod(?:o|ido)", ("iodo", "I", 1)),
            (r"hydrox(?:o|ide)", ("hydroxo", "OH", 1)),
            (r"cyan(?:o|o|ide|ido)", ("cyano", "CN", 1)),
            (r"thiocyanat(?:o)?", ("thiocyanato", "SCN", 1)),
            (r"isothiocyanat(?:o)?", ("isothiocyanato", "NCS", 1)),
            (r"\bnitro(?!\s*rito)", ("nitro", "NO2", 1)),
            (r"nitrit(?:o-[no])|(?:nitrito)", ("nitrito", "ONO", 1)),
            (r"oxalat(?:o)?", ("oxalato", "C2O4", 2)),
            (r"carbonat(?:o)?", ("carbonato", "CO3", 2)),
            (r"sulfat(?:o)?", ("sulfato", "SO4", 2)),
            (r"acetylacetonat(?:o)?|acac", ("acetylacetonato", "acac", 2)),
            (r"ethylenediamine|en\b", ("ethylenediamine", "en", 2)),
            (r"phenanthrolin(?:e)?|phen", ("phenanthroline", "phen", 2)),
            (r"bipyridin(?:e)?|bpy", ("bipyridine", "bpy", 2)),
            (r"ethylenediaminetetraacetat(?:o)?|edta", ("ethylenediaminetetraacetato", "EDTA", 6)),
        ]

        # -ate suffix → element symbol mapping (reverse of naming)
        self._ate_to_element = {
            "aluminate": "Al",
            "antimonate": "Sb",
            "arsenate": "As",
            "bismuthate": "Bi",
            "chromate": "Cr",
            "ferrate": "Fe",
            "indate": "In",
            "manganate": "Mn",
            "nickelate": "Ni",
            "plumbate": "Pb",
            "platinate": "Pt",
            "stannate": "Sn",
            "tantalate": "Ta",
            "titanate": "Ti",
            "tungstate": "W",
            "uranate": "U",
            "vanadate": "V",
            "zincate": "Zn",
            "aurate": "Au",
            "argentate": "Ag",
            "cuprate": "Cu",
            "cadmate": "Cd",
            "mercurate": "Hg",
            "gallate": "Ga",
            "germanate": "Ge",
            "stibate": "Sb",
            "plumbate": "Pb",
            "cobaltate": "Co",
            "palladate": "Pd",
            "osmate": "Os",
            "iridate": "Ir",
            "rhodate": "Rh",
            "ruthenate": "Ru",
        }

        # Direct metal name → element (for cationic complexes)
        self._metal_names = {
            "cobalt": "Co", "iron": "Fe", "copper": "Cu", "silver": "Ag",
            "gold": "Au", "zinc": "Zn", "nickel": "Ni", "chromium": "Cr",
            "manganese": "Mn", "vanadium": "V", "titanium": "Ti",
            "platinum": "Pt", "palladium": "Pd", "mercury": "Hg",
            "aluminum": "Al", "tin": "Sn", "lead": "Pb", "calcium": "Ca",
            "magnesium": "Mg", "barium": "Ba", "strontium": "Sr",
            "scandium": "Sc", "rhodium": "Rh", "iridium": "Ir",
            "osmium": "Os", "ruthenium": "Ru", "technetium": "Tc",
            "rhenium": "Re", "tungsten": "W", "uranium": "U",
            "gallium": "Ga", "germanium": "Ge", "indium": "In",
            "thallium": "Tl", "bismuth": "Bi", "antimony": "Sb",
            "arsenic": "As", "cadmium": "Cd", "boron": "B",
            "silicon": "Si", "phosphorus": "P",
        }

        # Prefix patterns
        self._prefix_pattern = r"^(?:(bis|tris|tetrakis|pentakis|hexakis|heptakis|octakis|nonakis|decakis)\(([^)]+)\)|(?:(mono)?)(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?(.+?)$"

        # Roman numerals
        self._roman_map = {
            "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5,
            "VI": 6, "VII": 7, "VIII": 8, "IX": 9, "X": 10,
            "XI": 11, "XII": 12,
        }

        # Ligand charges for charge estimation
        self._ligand_charges = {
            "aqua": 0, "ammine": 0, "carbonyl": 0, "nitrosyl": 0,
            "fluoro": -1, "chloro": -1, "bromo": -1, "iodo": -1,
            "hydroxo": -1, "cyano": -1, "thiocyanato": -1,
            "isothiocyanato": -1, "nitro": -1, "nitrito": -1,
            "oxalato": -2, "carbonato": -2, "sulfato": -2,
            "acetylacetonato": -1, "ethylenediamine": 0,
            "phenanthroline": 0, "bipyridine": 0,
            "ethylenediaminetetraacetato": -4,
        }

    def _parse_roman(self, roman_str: str) -> int:
        """Parse Roman numeral to integer."""
        r = roman_str.strip().upper()
        if r in self._roman_map:
            return self._roman_map[r]
        raise ChemMCPError(f"Cannot parse Roman numeral: '{roman_str}'")

    def _extract_ligand_block(self, name_without_metal: str) -> list:
        """Extract ligand prefix + name pairs from the name before the metal."""
        ligands = []
        remaining = name_without_metal

        while remaining:
            remaining = remaining.strip()
            if not remaining:
                break

            # Try special prefix with parentheses: tris(ethylenediamine)
            m_special = re.match(
                r"^(bis|tris|tetrakis|pentakis|hexakis|heptakis|octakis)\(([^)]+)\)(.*)",
                remaining, re.IGNORECASE
            )
            if m_special:
                prefix_word = m_special.group(1)
                lg_content = m_special.group(2)
                remaining = m_special.group(3)

                count_map = {"bis": 2, "tris": 3, "tetrakis": 4, "pentakis": 5, "hexakis": 6}
                count = count_map.get(prefix_word.lower(), 1)

                canon, formula, dent = self._match_ligand(lg_content)
                ligands.append({"name": canon, "count": count, "formula": formula, "denticity": dent})
                continue

            # Try regular Greek prefix + ligand name
            # Match greedily but stop at known ligand boundaries
            m_normal = re.match(
                r"^((?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?([a-z][a-z\-]*)(.*)",
                remaining, re.IGNORECASE
            )
            if m_normal:
                prefix_str = m_normal.group(1) or ""
                lg_name_candidate = m_normal.group(2)
                rest = m_normal.group(3) or ""

                # Check if lg_name_candidate contains multiple ligands — split at known boundaries
                best_split = self._split_ligand_name(lg_name_candidate)
                if best_split:
                    lg_name_candidate = best_split[0]
                    rest = best_split[1] + rest

                remaining = rest

                count_map_n = {
                    "di": 2, "tri": 3, "tetra": 4, "penta": 5, "hexa": 6,
                    "hepta": 7, "octa": 8, "nona": 9, "deca": 10,
                }
                count = count_map_n.get(prefix_str.lower(), 1) if prefix_str else 1

                canon, formula, dent = self._match_ligand(lg_name_candidate)
                ligands.append({"name": canon, "count": count, "formula": formula, "denticity": dent})
                continue

            # Fallback: consume one word-like token
            m_fallback = re.match(r"^([a-zA-Z]+\-?[a-zA-Z]*)(.*)", remaining)
            if m_fallback:
                token = m_fallback.group(1)
                remaining = m_fallback.group(2)
                canon, formula, dent = self._match_ligand(token)
                ligands.append({"name": canon, "count": 1, "formula": formula, "denticity": dent})
                continue

            break

        return ligands

    def _match_ligand(self, name: str) -> tuple:
        """Match a ligand name against known patterns."""
        name_lower = name.lower().strip()
        for pattern, result in self._ligand_patterns:
            if re.search(pattern, name_lower):
                return result
        # Return as-is if unknown
        return (name_lower, name_lower, 1)

    def _split_ligand_name(self, candidate: str):
        """Try to split a greedy ligand name match into (first_ligand, remainder).
        Returns None if no split is needed."""
        c = candidate.lower()
        # Sort patterns by length descending to match longest first
        sorted_patterns = sorted(self._ligand_patterns, key=lambda p: len(p[0]), reverse=True)
        for pattern, (canon, formula, dent) in sorted_patterns:
            m = re.match(pattern, c)
            if m and m.end() < len(c):  # partial match — something follows
                return (canon, c[m.end():])
        return None

    def _run_base(self, iupac_name: str) -> dict:
        """Parse IUPAC name into structural components."""
        name = iupac_name.strip()

        # Step 1: Extract oxidation state from Roman numerals in parentheses at end
        ox_match = re.search(r"\(([IVXLCDM]+)\)\s*$", name)
        if not ox_match:
            raise ChemMCPError(
                f"Cannot find oxidation state (Roman numerals in parentheses) at end of name: '{iupac_name}'"
            )
        oxidation_state = self._parse_roman(ox_match.group(1))
        name_before_ox = name[:ox_match.start()].strip()

        # Step 2: Detect if anionic (-ate suffix on metal)
        is_anionic = bool(re.search(r"ate$", name_before_ox, re.IGNORECASE))

        # Step 3: Find where metal name starts (scan backwards for -ate or known metal name)
        metal_part, ligand_part = self._split_metal_ligands(name_before_ox, is_anionic)

        # Step 4: Identify metal element
        metal_lower = metal_part.lower().strip()
        if is_anionic:
            metal_element = self._ate_to_element.get(metal_lower)
            if not metal_element:
                # Try partial match
                for ate_name, elem in self._ate_to_element.items():
                    if ate_name in metal_lower or metal_lower in ate_name:
                        metal_element = elem
                        break
            if not metal_element:
                metal_element = metal_part.capitalize()
                logger.warning(f"Unknown -ate metal: {metal_part}, using as-is")
        else:
            metal_element = self._metal_names.get(metal_lower, metal_part.capitalize())

        # Step 5: Parse ligands
        ligands = self._extract_ligand_block(ligand_part)

        # Step 6: Calculate coordination number
        coord_num = sum(lg["count"] * lg["denticity"] for lg in ligands)

        # Step 7: Estimate complex charge
        total_ligand_charge = sum(
            self._ligand_charges.get(lg["name"], 0) * lg["count"] for lg in ligands
        )
        complex_charge = oxidation_state + total_ligand_charge

        # Step 8: Build structural formula
        formula_parts = []
        for lg in ligands:
            if lg["count"] > 1:
                formula_parts.append(f"({lg['formula']}){lg['count']}")
            else:
                formula_parts.append(lg["formula"])

        charge_str = ""
        if complex_charge > 0:
            charge_str = "+" + str(complex_charge) if complex_charge > 1 else "+"
        elif complex_charge < 0:
            charge_str = str(complex_charge)

        structural_formula = f"[{metal_element}({''.join(formula_parts)})]{charge_str}"

        logger.info(f"Parsed '{iupac_name}' → {structural_formula}")

        return {
            "metal_ion": f"{metal_element}{charge_str.replace('+', '+').replace('-', '-')}" if charge_str else f"{metal_element}{oxidation_state:+d}",
            "metal_element": metal_element,
            "oxidation_state": oxidation_state,
            "ligands": [{"name": lg["name"], "count": lg["count"], "formula": lg["formula"]} for lg in ligands],
            "coordination_number": coord_num,
            "complex_charge": complex_charge,
            "structural_formula": structural_formula,
            "is_anionic": is_anionic,
            "parse_confidence": "high" if metal_element and ligands else "low",
        }

    def _split_metal_ligands(self, name_before_ox: str, is_anionic: bool) -> tuple:
        """Split the name into metal part and ligand part."""
        if is_anionic:
            # Metal ends with -ate; find it
            ate_match = re.search(r"([a-zA-Z]+)ate$", name_before_ox, re.IGNORECASE)
            if ate_match:
                return ate_match.group(0), name_before_ox[:ate_match.start()]
            # Fallback: last word
            parts = name_before_ox.rsplit(None, 1)
            return parts[-1], parts[0] if len(parts) > 1 else ""
        else:
            # Cationic: metal name is the last word before (Roman)
            # Known metal names
            for metal_name in sorted(self._metal_names.keys(), key=len, reverse=True):
                pattern = re.compile(rf"({re.escape(metal_name)})$", re.IGNORECASE)
                m = pattern.search(name_before_ox)
                if m:
                    return m.group(1), name_before_ox[:m.start()]
            # Fallback
            parts = name_before_ox.rsplit(None, 1)
            return parts[-1], parts[0] if len(parts) > 1 else ""

    def _run_text(self, query: str) -> dict:
        """Text interface delegates to base."""
        return self._run_base(query)
