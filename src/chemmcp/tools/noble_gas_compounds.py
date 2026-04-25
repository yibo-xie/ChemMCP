import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NobleGasCompounds(BaseTool):
    """
    稀有气体化合物查询工具。
    覆盖 He, Ne, Ar, Kr, Xe, Rn 的化合物（尤其是Xe的各类氟化物、氧化物、含氧酸等），
    以及稀有气体化合物的合成条件、结构和性质。
    """
    __version__ = "0.1.0"
    name = "NobleGasCompounds"
    func_name = "get_noble_gas_compounds"
    description = "Query noble gas (Group 18) compounds, focusing on xenon fluorides/oxides/oxyfluorides, krypton difluoride, radon compounds, and the discovery history of noble gas chemistry."
    implementation_description = "Built-in database of known noble gas compounds including synthesis conditions, molecular structures, oxidation states, chemical properties, and historical context. Covers Neil Bartlett's 1962 breakthrough (XePtF6) that shattered the 'inert gas' dogma."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Noble Gases", "Group 18", "Xenon Compounds", "Krypton", "Inert Gas Chemistry"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'Xe', 'xenon', 'all' for all)."),
        ("property_type", "str", "compounds", "'compounds', 'synthesis', 'structure', 'history', 'trends', or 'all'."),
        ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element [property_type]'. Example: 'Xe compounds' or 'all history'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested data."),
    ]

    examples = [
        {
            "code_input": {"element": "Xe", "property_type": "compounds"},
            "text_input": {"input_params": "Xe compounds"},
            "output": {"result": {"element": "Xenon", "known_compounds": [...], "oxidation_states": [...]}}
        },
    ]

    DATABASE = {
        "He": {
            "name": "Helium",
            "known_compounds": [
                {"formula": "HeH⁺", "name": "Helium hydride ion", "notes": "Strongest known acid; detected in interstellar space; not a neutral compound"},
                {"formula": "Na₂He", "name": "Sodium helium compound", "notes": "Stable only at >113 GPa; electride-like structure; predicted then synthesized (2017)"},
            ],
            "chemistry_possible": "Extremely limited — highest IE of any element (2372 kJ/mol); only under extreme pressure or as ions",
            "clathrates": "Forms clathrate compounds (He trapped in ice/crystal lattices) but no true chemical bonds",
        },
        "Ne": {
            "name": "Neon",
            "known_compounds": [
                {"formula": "Ne(AuF)(SbF₆)₃", "name": "Neon-gold complex", "notes": "Theoretical/predicted; extremely high-pressure species"},
                {"formula": "[Ne⊂C₆₀]", "name": "Ne@C60 endohedral fullerene", "notes": "Ne atom trapped inside C60 cage; physical encapsulation not chemical bonding"},
            ],
            "chemistry_possible": "Virtually none under normal conditions — second highest IE (2081 kJ/mol)",
        },
        "Ar": {
            "name": "Argon",
            "known_compounds": [
                {"formula": "HArF", "name": "Argon fluorohydride", "notes": "First true neutral argon compound (2000); stable only below 27 K in matrix isolation"},
                {"formula": "ArF₂", "name": "Argon difluoride", "notes": "Metastable; observed only in low-T matrices"},
                {"formula": "[AuF][SbF₆]·Ar", "name": "Argon coordination compound", "notes": "Ar coordinated to cationic metal center (2019)"},
                {"formula": "Ar@C60/C70", "name": "Endohedral fullerenes", "notes": "Encapsulated Ar atoms"},
            ],
            "chemistry_limited": "Very limited — requires extreme conditions or matrix isolation",
        },
        "Kr": {
            "name": "Krypton",
            "known_compounds": [
                {"formula": "KrF₂", "name": "Krypton difluoride", "notes": "MOST IMPORTANT Kr compound! Colorless crystalline solid; powerful fluorinating agent; linear molecule (D∞h); decomposes >130°C; photolytic or electric discharge synthesis from Kr + F₂"},
                {"formula": "Kr(OTeF₅)₂", "name": "Bis(pentafluoroorthotellurate)krypton", "notes": "Analogous to Xe compounds with OTeF5 ligand"},
                {"formula": "[KrF]⁺[MF₆]⁻", "name": "KrF+ salts", "notes": "KrF⁺[AsF₆]⁻, KrF⁺[SbF₆]⁻; strong oxidizers"},
                {"formula": "Kr@C60", "name": "Endohedral fullerene", "notes": "Encapsulated Kr"},
            ],
            "oxidation_states_known": "+2 (only)",
            "key_fact": "KrF₂ is one of the strongest known fluorinating agents (stronger than F₂ itself thermodynamically for some reactions)",
        },
        "Xe": {
            "name": "Xenon",
            "known_compounds": [
                # Fluorides
                {"formula": "XeF₂", "name": "Xenon difluoride", "ox_state": +2, "geometry": "Linear (D∞h)", "color": "Colorless crystals",
                 "synthesis": "Xe + F₂ → XeF₂ (excess Xe, 400°C, 1 atm; or UV light/silent discharge)", "properties": "Stable solid; mild fluorinating agent; hydrolyzes slowly; used as etchant in microelectronics"},
                {"formula": "XeF₄", "name": "Xenon tetrafluoride", "ox_state": +4, "geometry": "Square planar (D4h)", "color": "Colorless crystals",
                 "synthesis": "Xe + 2F₂ → XeF₄ (400°C, 6 atm, 1:5 Xe:F₂ ratio)", "properties": "Reactive fluorinating agent; reacts violently with water; strong oxidizer"},
                {"formula": "XeF₆", "name": "Xenon hexafluoride", "ox_state": +6, "geometry": "Distorted octahedral (fluxional C3v)", "color": "Colorless crystals",
                 "synthesis": "Xe + 3F₂ → XeF₆ (250°C, 50 atm, excess F₂)", "properties": "Most reactive Xe fluoride; reacts violently with water; strong fluorinating/oxidizing agent"},
                # Oxides
                {"formula": "XeO₃", "name": "Xenon trioxide", "ox_state": +6, "geometry": "Trigonal pyramidal", "color": "Colorless explosive crystals",
                 "synthesis": "XeF₆ + 3H₂O → XeO₃ + 6HF (careful hydrolysis of XeF₆)", "properties": "HIGHLY EXPLOSIVE when dry; powerful oxidizer; dissolves in water to give xenic acid (H₂XeO₄)"},
                {"formula": "XeO₄", "name": "Xenon tetroxide", "ox_state": +8, "geometry": "Tetrahedral", "color": "Yellow explosive solid",
                 "synthesis": "From perxenates (Na₄XeO₆) + concentrated H₂SO₄", "properties": "EXTREMELY EXPLOSIVE; unstable above -35°C"},
                # Oxyfluorides
                {"formula": "XeOF₂", "name": "Xenon oxydifluoride", "ox_state": +4, "geometry": "See-saw (T-shaped derivative)", "color": "Colorless",
                 "synthesis": "Partial hydrolysis of XeF₄", "notes": "Less common"},
                {"formula": "XeOF₄", "name": "Xenon oxytetrafluoride", "ox_state": +6, "geometry": "Square pyramidal (C4v)", "color": "Colorless liquid",
                 "synthesis": "XeF₆ + H₂O (controlled) → XeOF₄ + 2HF", "properties": "Reacts with water to give XeO₃"},
                {"formula": "XeO₂F₂", "name": "Xenon dioxydifluoride", "ox_state": +6, "geometry": "See-saw-like", "color": "Colorless",
                 "notes": "Unstable intermediate"},
                {"formula": "XeO₃F₂", "name": "Xenon trioxydifluoride", "ox_state": +8, "geometry": "Trigonal bipyramidal", "color": "Colorless",
                 "notes": "Very unstable"},
                {"formula": "XeO₄", "name": "Xenon tetroxide", "ox_state": +8, "geometry": "Tetrahedral", "color": "Yellow",
                 "notes": "Explosive above -35°C"},
                # Other
                {"formula": "Xe[MF₆] salts", "name": "Xenon hexafluorometallate(V) salts", "ox_state": +1, "examples": "Xe⁺[PtF₆]⁻ (THE FIRST noble gas compound!), Xe⁺[RuF₆]⁻, Xe⁺[PuF₆]⁻",
                 "historical_significance": "Neil Bartlett's 1962 discovery of XePtF₆ proved noble gases CAN form compounds — ended 'inert gas' era"},
                # Perxenates
                {"formula": "Na₄XeO₆", "name": "Sodium perxenate", "ox_state": +8, "notes": "Powerful oxidizing agent; source of Xe(VIII)"},
                {"formula": "H₄XeO₆ / H₂XeO₆", "name": "Perxenic acid", "ox_state": +8, "notes": "Very weak acid; powerful oxidizer"},
            ],
            "oxidation_states": [+2, +4, +6, +8],
            "most_stable_ox_state": "+4 (XeF₄) and +6 (XeF₆) are most stable; +8 is strongly oxidizing",
            "VSEPR_notes": "XeF₂: linear (3 lone pairs), XeF₄: square planar (2 LP), XeF₆: distorted octahedron (1 LP)",
        },
        "Rn": {
            "name": "Radon",
            "known_compounds": [
                {"formula": "RnF₂", "name": "Radon difluoride", "notes": "Predicted/claimed; ionic character expected (Rn²⁺ more stable than Xe²⁺ due to lower IE); difficult to study due to radioactivity"},
                {"formula": "RnF₄", "name": "Radon tetrafluoride", "notes": "Reported but not fully characterized; should be more stable than XeF₄ based on trends"},
                {"formula": "RnO₃/RnO₄?", "name": "Radon oxides", "notes": "Predicted but unconfirmed; radioactivity makes characterization nearly impossible"},
            ],
            "study_difficulty": "All Rn isotopes radioactive (longest: ²²²Rn, t½ = 3.8 days); highly toxic (alpha emitter, lung cancer risk); compounds studied only in trace amounts",
            "predicted_chemistry": "Should be MORE reactive than Xe (lower ionization energy); RnF₂ likely ionic (Rn²⁺ + 2F⁻)",
        },
    }

    HISTORY = """
    ## The Death of "Inert Gas" Dogma

    ### Before 1962
    - Noble gases considered completely "inert" — full valence shell (ns²np⁶)
    - No known compounds; Group 18 called "zero group"
    - Linus Pauling predicted XeF₆ could exist (1933) but was ignored

    ### 1962: The Breakthrough
    **Neil Bartlett** (University of British Columbia) noticed PtF₆ was an incredibly strong oxidizing agent:
    - PtF₆ could oxidize O₂ to O₂⁺[PtF₆]⁻
    - IE(O₂ → O₂⁺) = 1175 kJ/mol
    - IE(Xe → Xe⁺) = 1170 kJ/mol — SLIGHTLY LOWER!
    - Conclusion: PtF₆ should be able to oxidize Xe!

    **Experiment:** Mixed deep-red PtF₆ vapor with Xe gas at room temperature
    **Result:** Immediate formation of yellow-orange solid — **Xe⁺[PtF₆]⁻**
    This was the FIRST EVER noble gas compound.

    ### After 1962 (Gold Rush Period)
    Within months, multiple groups reported:
    - XeF₂, XeF₄, XeF₆ (direct reaction of Xe with F₂)
    - XeO₃, oxoacids, perxenates
    - KrF₂ discovered (1963)

    ### Modern Understanding
    Noble gases are NOT inert — they are simply **less reactive**. Their chemistry is dominated by:
    1. Highly electronegative partners (F, O in high oxidation states)
    2. Strong oxidizing conditions
    3. Large polarizable atoms (Xe, Kr) where relativistic effects help
    """

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "compounds") -> dict:
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()

        if prop_type == "history":
            return {"result": {"history_text": self.HISTORY.strip()}}

        if element == "All":
            result = {}
            for sym in self.DATABASE:
                d = self.DATABASE[sym]
                if prop_type == "all":
                    result[sym] = d
                elif prop_type == "compounds":
                    result[sym] = {"name": d.get("name"), "known_compounds": d.get("known_compounds", []),
                                     "oxidation_states": d.get("oxidation_states")}
                else:
                    result[sym] = {k: v for k, v in d.items() if k != "name"}
            return {"result": result}

        if element not in self.DATABASE:
            raise ChemMCPError(f"Element '{element}' not found. Options: {list(self.DATABASE.keys()) + ['All']}")

        data = self.DATABASE[element]
        return {"result": {**{"element": element}, **data}}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            elem = parts[0] if parts else "Xe"
            prop = parts[1] if len(parts) > 1 else "compounds"
            return self._run_base(elem, prop)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'element [property_type]'")
