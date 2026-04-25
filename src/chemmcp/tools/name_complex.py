"""
配合物 IUPAC 命名工具
IUPAC naming for coordination compounds.
"""
import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NameComplex(BaseTool):
    """
    按照IUPAC规则对配合物进行命名。
    输入金属离子、配体列表、配位数、氧化态，输出标准IUPAC名称。
    """
    __version__ = "0.1.0"
    name = "NameComplex"
    func_name = "name_complex"
    description = "Generate IUPAC name for a coordination compound from its structural components."
    implementation_description = "Applies IUPAC Red Book naming rules: ligands in alphabetical order → metal name (with -ate suffix if anionic complex) → oxidation state in Roman numerals. Handles common ligand names, prefixes, and special cases."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Nomenclature", "IUPAC", "Complexes"]
    required_envs = []

    code_input_sig = [
        ("metal_ion", "str", "N/A", "Metal element symbol, e.g., 'Co', 'Fe', 'Pt', 'Cu'."),
        ("ligands", "list", "N/A", "List of ligand names (strings), e.g., ['ammine', 'chloro', 'aqua']."),
        ("coordination_number", "int", "6", "Total coordination number (total number of ligand donor atoms)."),
        ("oxidation_state", "int", "N/A", "Oxidation state of the metal center (integer)."),
        ("is_cationic_complex", "bool", "True", "Whether the coordination entity is cationic/neutral (True) or anionic (False). If anionic, metal gets '-ate' suffix."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'metal_ion oxidation_state ligand1:count ligand2:count ... [anionic]', e.g., 'Co 3 ammine:5 chloro:1' or 'Fe 3 cyano:6 anionic'."),
    ]

    output_sig = [
        ("iupac_name", "str", "The IUPAC name of the coordination compound."),
        ("formula_suggestion", "str", "Suggested structural formula string."),
        ("naming_breakdown", "dict", "Step-by-step breakdown of how the name was constructed."),
    ]

    examples = [
        {
            "code_input": {
                "metal_ion": "Co",
                "ligands": ["ammine", "chloro"],
                "coordination_number": 6,
                "oxidation_state": 3,
                "is_cationic_complex": True,
            },
            "text_input": {
                "query": "Co 3 ammine:5 chloro:1"
            },
            "output": {
                "iupac_name": "pentaamminechlorocobalt(III)",
                "formula_suggestion": "[Co(NH3)5Cl]?",
                "naming_breakdown": {
                    "ligands_ordered": "ammine (×5), chloro (×1)",
                    "prefixes": "penta-, (mono omitted)",
                    "metal_part": "cobalt",
                    "oxidation": "(III)",
                    "charge_note": "Cationic complex",
                }
            }
        },
        {
            "code_input": {
                "metal_ion": "Fe",
                "ligands": ["cyano"],
                "coordination_number": 6,
                "oxidation_state": 2,
                "is_cationic_complex": False,
            },
            "text_input": {
                "query": "Fe 2 cyano:6 anionic"
            },
            "output": {
                "iupac_name": "hexacyanoferrate(II)",
                "formula_suggestion": "[Fe(CN)6]4-",
                "naming_breakdown": {
                    "ligands_ordered": "cyano (×6)",
                    "prefixes": "hexa-",
                    "metal_part": "ferrate (iron → ferrate for anionic)",
                    "oxidation": "(II)",
                    "charge_note": "Anionic complex",
                }
            }
        },
        {
            "code_input": {
                "metal_ion": "Pt",
                "ligands": ["ammine", "chloro"],
                "coordination_number": 4,
                "oxidation_state": 2,
                "is_cationic_complex": False,
            },
            "text_input": {
                "query": "Pt 2 ammine:2 chloro:2 anionic"
            },
            "output": {
                "iupac_name": "diamminedichloroplatinate(II)",
                "formula_suggestion": "[Pt(NH3)2Cl2]2-",
                "naming_breakdown": {
                    "ligands_ordered": "ammine (×2), chloro (×2)",
                    "prefixes": "di-, di-",
                    "metal_part": "platinate (platinum → platinate for anionic)",
                    "oxidation": "(II)",
                    "charge_note": "Anionic complex",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize ligand naming database and rules."""
        # Ligand name map: canonical name → (alphabetical key, formula abbreviation, denticity)
        self._ligand_db = {
            # Neutral ligands
            "aqua": ("aqua", "H2O", 1),
            "water": ("aqua", "H2O", 1),
            "ammine": ("ammine", "NH3", 1),
            "ammonia": ("ammine", "NH3", 1),
            "carbonyl": ("carbonyl", "CO", 1),
            "nitrosyl": ("nitrosyl", "NO", 1),
            # Anionic ligands (end with -o)
            "fluoro": ("fluoro", "F", 1),
            "fluorido": ("fluoro", "F", 1),
            "chloro": ("chloro", "Cl", 1),
            "chlorido": ("chloro", "Cl", 1),
            "bromo": ("bromo", "Br", 1),
            "bromido": ("bromo", "Br", 1),
            "iodo": ("iodo", "I", 1),
            "iodido": ("iodo", "I", 1),
            "hydroxo": ("hydroxo", "OH", 1),
            "hydroxide": ("hydroxo", "OH", 1),
            "cyano": ("cyano", "CN", 1),
            "cyanide": ("cyano", "CN", 1),
            "thiocyanato": ("thiocyanato", "SCN", 1),
            "thiocyanide": ("thiocyanato", "SCN", 1),
            "isothiocyanato": ("isothiocyanato", "NCS", 1),
            "nitro": ("nitro", "NO2", 1),
            "nitrito-n": ("nitrito-N", "ONO", 1),
            "nitrito-o": ("nitrito-O", "ONO", 1),
            "oxalato": ("oxalato", "C2O4", 2),
            "ox": ("oxalato", "C2O4", 2),
            "carbonato": ("carbonato", "CO3", 2),
            "sulfato": ("sulfato", "SO4", 2),
            "acetylacetonato": ("acetylacetonato", "acac", 2),
            "acac": ("acetylacetonato", "acac", 2),
            # Chelating diamines
            "ethylenediamine": ("ethylenediamine", "en", 2),
            "en": ("ethylenediamine", "en", 2),
            "1,10-phenanthroline": ("phenanthroline", "phen", 2),
            "phen": ("phenanthroline", "phen", 2),
            "2,2'-bipyridine": ("bipyridine", "bpy", 2),
            "bpy": ("bipyridine", "bpy", 2),
            "ethylenediaminetetraacetato": ("ethylenediaminetetraacetato", "EDTA", 6),
            "edta": ("ethylenediaminetetraacetato", "EDTA", 6),
        }

        # Greek prefixes for multiplicity
        self._prefixes = {
            1: "",       # mono- usually omitted
            2: "di-",
            3: "tri-",
            4: "tetra-",
            5: "penta-",
            6: "hexa-",
            7: "hepta-",
            8: "octa-",
            9: "nona-",
            10: "deca-",
        }

        # Special prefixes for ligand names already containing Greek numerals
        self._special_prefixes = {
            1: "",
            2: "bis-",
            3: "tris-",
            4: "tetrakis-",
            5: "pentakis-",
            6: "hexakis-",
        }

        # Metals that change name in anionic complexes
        self._metal_ate_names = {
            "al": "aluminate",
            "aluminum": "aluminate",
            "sb": "antimonate",
            "antimony": "antimonate",
            "as": "arsenate",
            "arsenic": "arsenate",
            "bi": "bismuthate",
            "bismuth": "bismuthate",
            "cr": "chromate",
            "chromium": "chromate",
            "fe": "ferrate",
            "iron": "ferrate",
            "in": "indate",
            "indium": "indate",
            "mn": "manganate",
            "manganese": "manganate",
            "nb": "niobate",
            "niobium": "niobate",
            "ni": "nickelate",
            "nickel": "nickelate",
            "pb": "plumbate",
            "lead": "plumbate",
            "pt": "platinate",
            "platinum": "platinate",
            "sn": "stannate",
            "tin": "stannate",
            "ta": "tantalate",
            "tantalum": "tantalate",
            "ti": "titanate",
            "titanium": "titanate",
            "w": "tungstate",
            "tungsten": "tungstate",
            "u": "uranate",
            "uranium": "uranate",
            "v": "vanadate",
            "vanadium": "vanadate",
            "zn": "zincate",
            "zinc": "zincate",
        }

        # Roman numeral conversion
        self._roman_numerals = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
            (1, "I")
        ]

    def _to_roman(self, num: int) -> str:
        """Convert integer to Roman numerals."""
        result = ""
        for value, symbol in self._roman_numerals:
            while num >= value:
                result += symbol
                num -= value
        return result

    def _needs_special_prefix(self, ligand_name: str) -> bool:
        """Check if ligand name requires bis-/tris-/etc. prefix."""
        special_keywords = [
            "ethyl", "di", "tri", "tetra", "phen", "bipy", "en ",
            "acetyl", "phenanthro", "bipyridine", "ethylenediamine",
        ]
        ln_lower = ligand_name.lower()
        return any(kw in ln_lower for kw in special_keywords)

    def _run_base(self, metal_ion: str, ligands: List[str],
                  coordination_number: int = 6, oxidation_state: int = 0,
                  is_cationic_complex: bool = True) -> dict:
        """Generate IUPAC name for a coordination compound."""
        if not ligands:
            raise ChemMCPError("At least one ligand must be provided.")

        if oxidation_state == 0:
            raise ChemMCPError("Oxidation state is required.")

        # Step 1: Resolve each ligand to canonical name and count occurrences
        resolved_ligands = []  # list of (canonical_name, formula, denticity)
        for lg in ligands:
            lg_lower = lg.lower().strip()
            if lg_lower not in self._ligand_db:
                raise ChemMCPError(
                    f"Unknown ligand: '{lg}'. Known ligands include: "
                    f"aqua, ammine, carbonyl, fluoro, chloro, bromo, iodo, "
                    f"hydroxo, cyano, thiocyanato, nitro, oxalato, "
                    f"ethylenediamine(en), phenanthroline(phen), bipyridine(bpy), EDTA, acac."
                )
            entry = self._ligand_db[lg_lower]
            resolved_ligands.append(entry)

        # Step 2: Count ligand occurrences (group by canonical alphabetical name)
        ligand_counts = {}
        for canon_name, formula, dent in resolved_ligands:
            key = canon_name  # alphabetically sorted key
            if key not in ligand_counts:
                ligand_counts[key] = {"count": 0, "formula": formula, "denticity": dent}
            ligand_counts[key]["count"] += 1

        # Step 3: Sort ligands alphabetically by canonical name (IUPAC rule)
        sorted_ligands = sorted(ligand_counts.items(), key=lambda x: x[0])

        # Step 4: Build name parts
        name_parts = []
        formula_parts = []
        ligand_detail = []

        for canon_name, info in sorted_ligands:
            count = info["count"]
            formula = info["formula"]

            # Choose prefix type
            if self._needs_special_prefix(canon_name):
                prefix = self._special_prefixes.get(count, f"{count}-")
            else:
                prefix = self._prefixes.get(count, f"{count}-")

            # Omit "mono-" prefix
            if count == 1:
                prefix = ""

            name_parts.append(f"{prefix}{canon_name}")
            if count > 1:
                formula_parts.append(f"({formula}){count}")
            else:
                formula_parts.append(formula)
            ligand_detail.append(f"{canon_name} (×{count})")

        # Step 5: Metal name (with -ate suffix if anionic complex)
        metal_lower = metal_ion.lower().strip()
        if not is_cationic_complex:
            metal_name = self._metal_ate_names.get(metal_lower, f"{metal_ion}ate")
        else:
            metal_name = metal_ion.capitalize()

        # Step 6: Oxidation state in Roman numerals
        ox_roman = self._to_roman(abs(oxidation_state))

        # Assemble full name
        full_name = f"{''.join(name_parts)}{metal_name}({ox_roman})"

        # Suggest formula
        formula_str = f"[{metal_ion}({''.join(formula_parts)})]"
        if not is_cationic_complex:
            formula_str += "?"  # charge to be determined

        logger.info(f"IUPAC named: {full_name}")

        return {
            "iupac_name": full_name,
            "formula_suggestion": formula_str,
            "naming_breakdown": {
                "ligands_ordered": ", ".join(ligand_detail),
                "prefixes_used": [p for p in name_parts],
                "metal_part": metal_name + (" (→ -ate suffix)" if not is_cationic_complex else ""),
                "oxidation_state": f"{oxidation_state} → ({ox_roman})",
                "complex_type": "cationic/neutral" if is_cationic_complex else "anionic",
            }
        }

    def _run_text(self, query: str) -> dict:
        """Parse text query: 'metal oxidation_state ligand1[:count] ligand2[:count] ... [anionic]'"""
        parts = query.strip().split()

        if len(parts) < 3:
            raise ChemMCPError(
                "Format: 'metal oxidation_state ligand1[:count] ligand2[:count] ... [anionic]'\n"
                "Example: 'Co 3 ammine:5 chloro:1'\n"
                "Example: 'Fe 2 cyano:6 anionic'"
            )

        metal_ion = parts[0]
        try:
            oxidation_state = int(parts[1])
        except ValueError:
            raise ChemMCPError(f"Oxidation state must be an integer, got: '{parts[1]}'")

        # Parse remaining tokens as ligands
        ligands = []
        is_anionic = False
        for token in parts[2:]:
            tlower = token.lower()
            if tlower in ("anionic", "anion", "negative"):
                is_anionic = True
                continue
            if ":" in token:
                lg_name, count_str = token.rsplit(":", 1)
                try:
                    count = int(count_str)
                except ValueError:
                    count = 1
                for _ in range(count):
                    ligands.append(lg_name)
            else:
                ligands.append(token)

        return self._run_base(metal_ion, ligands, oxidation_state=oxidation_state,
                              is_cationic_complex=not is_anionic)
