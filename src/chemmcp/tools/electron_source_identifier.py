"""
Electron Source Identifier (Tool #152)
识别反应中的电子给体（电子源）：分析反应中氧化数降低的物种，
确定哪个物种提供电子，分析电子流向和还原过程。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# 常见还原剂及其失电子特征
_REDUCING_AGENTS = {
    # 分子/离子 -> (典型氧化态变化, 描述)
    'Na': ('0 → +1', 'Sodium: strong reducing agent, Na → Na⁺ + e⁻'),
    'K': ('0 → +1', 'Potassium: strong reducing agent, K → K⁺ + e⁻'),
    'Li': ('0 → +1', 'Lithium: strong reducing agent, Li → Li⁺ + e⁻'),
    'Mg': ('0 → +2', 'Magnesium: Mg → Mg²⁺ + 2e⁻'),
    'Zn': ('0 → +2', 'Zinc: Zn → Zn²⁺ + 2e⁻, mild reducing agent'),
    'Fe': ('0 → +2/+3', 'Iron: Fe → Fe²⁺/Fe³⁺'),
    'Al': ('0 → +3', 'Aluminum: Al → Al³⁺ + 3e⁻, thermite reaction'),
    'Sn': ('0 → +2/+4', 'Tin: Sn → Sn²⁺/Sn⁴⁺'),
    'H2': ('0 → +1', 'Hydrogen gas: H₂ → 2H⁺ + 2e⁻'),
    'C': ('0 → +II/+IV', 'Carbon: C → CO/CO₂, reducing in metallurgy'),
    'CO': ('+II → +IV', 'Carbon monoxide: CO → CO₂ (blast furnace)'),
    'SO2': ('+IV → +VI', 'Sulfur dioxide: SO₂ → SO₃/SO₄²⁻'),
    'H2S': ('-II → 0/+IV', 'Hydrogen sulfide: S(-II) → S(0)/SO₂'),
    'Na2S2O3': ('+II → +VI/-II→0', 'Thiosulfate: complex redox behavior'),
    'H2C2O4': ('+III → +IV', 'Oxalic acid: C₂O₄²⁻ → 2CO₂ + 2e⁻'),
    'NaBH4': ('H(-I) → +I', 'Sodium borohydride: hydride donor, H⁻ → H⁺'),
    'LiAlH4': ('H(-I) → +I', 'Lithium aluminum hydride: stronger hydride donor'),
    'N2H4': ('-II → 0', 'Hydrazine: N₂H₄ → N₂ + 4e⁻ + 4H⁺'),
    'NH3': ('-III → 0/-II', 'Ammonia: can be oxidized to N₂ or NO'),
    'HSO3-': ('+IV → +VI', 'Bisulfite: SO₃²⁻ → SO₄²⁻ + 2e⁻'),
    'SO3(2-)': ('+IV → +VI', 'Sulfite: SO₃²⁻ → SO₄²⁻ + 2e⁻'),
    'As2O3': ('+III → +V', 'Arsenic trioxide: As(III) → As(V)'),
    'H2O2(reducing)': ('-I → 0', 'H₂O₂ as reductant: H₂O₂ → O₂ + 2H⁺ + 2e⁻'),
    'SnCl2': ('+II → +IV', 'Stannous chloride: Sn(II) → Sn(IV), reduces Fe³⁺→Fe²⁺, Hg²⁺→Hg'),
    'FeSO4': ('+II → +III', 'Ferrous sulfate: Fe²⁺ → Fe³⁺ + e⁻'),
    'TiCl3': ('+III → +IV', 'Titanium(III): Ti³⁺ → Ti⁴⁺ + e⁻'),
    'VCl2': ('+II → +III', 'Vanadium(II): V²⁺ → V³⁺ + e⁻'),
    'Cr(II)': ('+II → +III', 'Chromium(II): Cr²⁺ → Cr³⁺ + e⁻, strong reductant'),
    'ascorbic_acid': ('organic', 'Ascorbic acid (vitamin C): oxidized to dehydroascorbic acid'),
    'oxalic_acid': ('+III → +IV', 'Oxalic acid: (COOH)₂ → 2CO₂ + 2e⁻'),
    'formate': ('+II → +IV', 'Formate ion: HCOO⁻ → CO₂ + 2e⁻ + H⁺'),
    'glucose': ('organic', 'Glucose: aldehyde group oxidized to carboxyl (reducing sugar test)'),
    'I-': ('-I → 0', 'Iodide: 2I⁻ → I₂ + 2e⁻'),
    'Br-': ('-I → 0', 'Bromide: 2Br⁻ → Br₂ + 2e⁻ (strong oxidant needed)'),
    'metal_carbonyl': ('low → higher', 'Metal carbonyls: electron-rich metal centers'),
}

# 官能团还原特征
_FUNCTIONAL_GROUP_REDUCTION = {
    'aldehyde': {'reduction_state_change': 'C(+I) → -I', 'product': 'primary alcohol', 'source': 'reducing_agent'},
    'ketone': {'reduction_state_change': 'C(0) → -I', 'product': 'secondary alcohol', 'source': 'reducing_agent'},
    'carboxylic_acid': {'reduction_state_change': 'C(+III) → -I', 'product': 'primary alcohol (strong reduction)', 'source': 'LiAlH4'},
    'ester': {'reduction_state_change': 'C(+III)/0 → -I/-I', 'product': 'two alcohols', 'source': 'LiAlH4/NaBH4(limited)'},
    'amide': {'reduction_state_change': 'C(+III) → -I', 'product': 'amine', 'source': 'LiAlH4'},
    'epoxide': {'reduction_state_change': 'ring opens', 'product': 'alcohol', 'source': 'LiAlH4/Red-Al'},
    'alkene': {'reduction_state_change': 'C(-II) → -III', 'product': 'alkane', 'source': 'H2/Pt or dissolving metal'},
    'alkyne': {'reduction_state_change': 'C(-I) → -II/-III', 'product': 'cis/trans-alkane or alkene', 'source': 'H2/Lindlar or Na/NH3'},
    'nitro': {'reduction_state_change': 'N(+III) → -III', 'product': 'amine', 'source': 'H2/Pd, Zn/HCl, Sn/HCl, Fe/HCl'},
    'nitrile': {'reduction_state_change': 'C(+II)N(-III) → -I(N)(-III)', 'product': 'primary amine', 'source': 'LiAlH4 or H2/Raney Ni'},
    'oxime': {'reduction_state_change': 'C(+I)N(-I) → -I(N)(-III)', 'product': 'amine', 'source': 'LiAlH4 or catalytic hydrogenation'},
    'azide': {'reduction_state_change': 'N(+I/3) → -III', 'product': 'primary amine', 'source': 'H2/Pd or PPh3/H2O'},
    'imine': {'reduction_state_change': 'C(+I)N(-II) → -I(N)(-III)', 'product': 'amine', 'source': 'NaBH4 or NaBH3CN or H2/Pd'},
    'disulfide': {'reduction_state_change': 'S(-I) → -II', 'product': 'thiol (2 RSH)', 'source': 'Zn/HCl or NaBH4'},
    'sulfoxide': {'reduction_state_change': 'S(0) → -II', 'product': 'sulfide', 'source': 'PPh3/I2 or MoOPH'},
    'peroxide': {'reduction_state_change': 'O(-I) → -II', 'product': 'alcohol/water', 'source': 'PPh3 or sulfide'},
    'quinone': {'reduction_state_change': 'quinone → hydroquinone', 'product': 'hydroquinone', 'source': 'Na2S2O4 or metal hydride'},
    'αβ_unsat_carbonyl': {'reduction_state_change': 'C=C reduced', 'product': 'saturated carbonyl', 'source': 'H2/Pd or dissolving metal (1,2 vs 1,4)'},
}


@ChemMCPManager.register_tool
class ElectronSourceIdentifier(BaseTool):
    __version__ = "0.1.0"
    name = "ElectronSourceIdentifier"
    func_name = 'identify_electron_source'
    description = "Identify the electron source (electron donor/reducing agent) in a chemical reaction. Analyze reaction SMILES, chemical equations, or redox processes to determine which species donates electrons and gets oxidized."
    implementation_description = "Uses pattern matching against a comprehensive database of known reducing agents, functional group reduction signatures, and oxidation state rules. Parses reactions to identify species that increase in oxidation state (lose electrons). Supports both inorganic and organic redox reactions."
    categories = ["Reaction"]
    tags = ["Electron Transfer", "Redox", "Reduction", "Reaction Analysis", "Reducing Agent"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("reaction_input", "str", "N/A", "Reaction SMILES string, chemical equation, or description of the redox process."),
        ("analysis_mode", "str", "detailed", "Analysis mode: 'brief' for summary, 'detailed' for full analysis with electron counting."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: reaction_input [analysis_mode]. E.g., 'CH3CHO+NaBH4→CH3CH2OH detailed'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing source_species, source_reasoning, electron_flow, oxidation_state_changes, confidence, and analysis_details."),
    ]

    examples = [
        {
            "code_input": {
                "reaction_input": "CH3CHO + NaBH4 → CH3CH2OH",
                "analysis_mode": "detailed",
            },
            "text_input": {"query": "CH3CHO+NaBH4→CH3CH2OH detailed"},
            "output": {
                "result": {
                    "primary_electron_source": "NaBH4 (sodium borohydride)",
                    "source_species": ["NaBH4", "hydride (H⁻)"],
                    "source_reasoning": "NaBH4 is a hydride-donating reducing agent. The hydride ion (H⁻) donates 2 electrons (as H:⁻) to the carbonyl carbon of acetaldehyde, reducing it to ethanol. Borohydride itself is oxidized as H(-I) → H(+I).",
                    "electrons_transferred": 2,
                    "half_reaction_oxidation": "H⁻ (from BH4⁻) → H⁺ + 2e⁻  (or BH4⁻ → B(OH)3 + ...)",
                    "electron_flow_description": "Hydride from NaBH4 attacks electrophilic carbonyl carbon. Carbon oxidation state changes from +I to -I (gains 2 e⁻ equivalent).",
                    "oxidation_state_changes": [
                        {"species": "C (carbonyl)", "from": "+I", "to": "-I", "change": "-2", "electrons_gained": 2},
                        {"species": "H (from BH4⁻)", "from": "-I", "to": "+I", "change": "+2", "electrons_lost": 2},
                    ],
                    "confidence": "high",
                    "reaction_type": "reduction (hydride transfer)",
                }
            },
        },
        {
            "code_input": {
                "reaction_input": "Zn + CuSO4 → ZnSO4 + Cu",
                "analysis_mode": "brief",
            },
            "text_input": {"query": "Zn+CuSO4→ZnSO4+Cu brief"},
            "output": {
                "result": {
                    "primary_electron_source": "Zn (zinc metal)",
                    "source_species": ["Zn(s)"],
                    "source_reasoning": "Zinc metal is more easily oxidized than copper. Zn → Zn²⁺ + 2e⁻, donating electrons to reduce Cu²⁺ → Cu(s).",
                    "electrons_transferred": 2,
                    "oxidation_state_changes": [
                        {"species": "Zn", "from": "0", "to": "+2", "change": "+2"},
                        {"species": "Cu", "from": "+2", "to": "0", "change": "-2"},
                    ],
                    "confidence": "high",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize the reducing agent database."""
        self.reducing_agents = dict(_REDUCING_AGENTS)
        self.fg_reduction = dict(_FUNCTIONAL_GROUP_REDUCTION)

    def _run_base(self, reaction_input: str, analysis_mode: str = "detailed") -> dict:
        """
        Core logic: identify electron sources in a reaction.
        """
        if not reaction_input or not reaction_input.strip():
            raise ChemMCPInputError("Reaction input cannot be empty.")

        rxn = reaction_input.strip()
        mode = analysis_mode.lower() if analysis_mode else "detailed"

        # Step 1: Try to match known reducing agents
        matched_sources = self._match_reducing_agents(rxn)

        # Step 2: Detect organic reduction patterns
        organic_patterns = self._detect_organic_reduction(rxn)

        # Step 3: Parse oxidation state changes
        ox_changes = self._parse_oxidation_changes(rxn)

        # Step 4: Determine primary electron source
        primary_source = self._determine_primary_source(matched_sources, organic_patterns, ox_changes)

        # Step 5: Count electrons transferred
        e_count = self._count_electrons(ox_changes, matched_sources)

        result = {
            "result": {
                "reaction_input": rxn,
                "primary_electron_source": primary_source.get("name", "unknown"),
                "source_species": primary_source.get("species", []),
                "source_reasoning": primary_source.get("reasoning", "Could not determine electron source definitively."),
                "matched_reducing_agents": matched_sources,
                "organic_reduction_patterns": organic_patterns,
                "oxidation_state_changes": ox_changes,
                "electrons_transferred": e_count.get("total", None),
                "electrons_transferred_detail": e_count,
                "electron_flow_description": self._describe_electron_flow(primary_source, ox_changes),
                "confidence": primary_source.get("confidence", "low"),
                "analysis_mode": mode,
                "additional_notes": self._generate_notes(primary_source, matched_sources, organic_patterns),
            }
        }

        if mode == "detailed":
            result["result"]["half_reactions"] = self._propose_half_reactions(primary_source, ox_changes)

        logger.info(f"ElectronSource: identified {primary_source.get('name','?')} as electron source in: {rxn[:50]}")
        return result

    def _match_reducing_agents(self, rxn: str) -> list:
        """Match known reducing agents in the reaction string."""
        matches = []
        rxn_upper = rxn.upper()

        for agent, info in self.reducing_agents.items():
            clean_agent = agent.replace('(', '').replace(')', '').replace('+', '')
            if clean_agent.upper() in rxn_upper or agent.replace('_', '') in rxn_upper:
                matches.append({
                    "agent": agent,
                    "oxidation_change": info[0],
                    "description": info[1],
                })

        # Pattern-based matches
        patterns = [
            (r'\[H\]|\[H-\]', '[H] / [H-]', 'Hydride notation', 'Hydride transfer reagent (NaBH4, LiAlH4, etc.)'),
            (r'NaBH4|borohydride', 'NaBH4', 'Sodium borohydride', 'Mild reducing agent: aldehydes/ketones → alcohols'),
            (r'LiAlH4|LAH|lithium.aluminum', 'LiAlH4', 'Lithium aluminum hydride', 'Strong reducing: esters/amides/acids → alcohols/amines'),
            (r'Red.Al|red.aluminum|Na\[(?:CH2CH2OCH3)3\]', 'Red-Al', 'Sodium bis(methoxyethoxy)aluminum hydride', 'Alternative to LAH, milder'),
            (r'DIBAL|DIBAL.H|i-Bu2AlH', 'DIBAL-H', 'Diisobutylaluminum hydride', 'Ester→aldehyde at low T; nitrile→aldehyde'),
            (r'H2[/\s]?Pt|H2[/\s]?Pd|H2[/\s]?Ni|H2[/\s]?Raney|catalytic.*hydrogen', 'H2/metal catalyst', 'Catalytic hydrogenation', 'H₂ + Pd/Pt/Ni: reduces C=C, C≡C, NO2, CN, C=O'),
            (r'Lindlar|Lindlar.*catalyst|CaCO3.Pd', 'Lindlar catalyst', 'Pd/CaCO4 poisoned', 'Alkyne → cis-alkene only'),
            (r'Na.*NH3|Na.*liq.*ammonia|dissolving.metal', 'Na/NH3(l)', 'Dissolving metal reduction', 'Alkyne → trans-alkene; Birch reduction arenes'),
            (r'Birch.*reduc|Na.*NH3.*arene', 'Birch conditions', 'Birch reduction', 'Arene → 1,4-cyclohexadiene'),
            (r'Sn.*HCl|tin.*HCl|Fe.*HCl|iron.*HCl', 'Sn/HCl or Fe/HCl', 'Metal/acid reduction', 'Nitro → amine (classic method)'),
            (r'Zn.*HCl|zinc.*HCl', 'Zn/HCl', 'Zinc/acid reduction', 'Nitro → amine; Clemmensen reduction (C=O→CH2)'),
            (r'Clemmensen|Zn(Hg).*HCl', 'Clemmensen', 'Clemmensen reduction', 'Carbonyl → methylene under acidic conditions'),
            (r'Wolff.Kishner|NH2NH2.*KOH.*heat', 'Wolff-Kishner', 'Wolff-Kishner reduction', 'Carbonyl → methylene under basic conditions'),
            (r'Mozingo|HS(CH2)2SH.*BF3|Mozingo.*reduc', 'Mozingo', 'Mozingo modification', 'Modified Wolff-Kishner, milder'),
            (r'SnCl2|stannous|Sn\(II\)', 'SnCl2', 'Stannous chloride', 'Reduces NO2→NHOH, Fe3+→Fe2+, Hg2+→Hg'),
            (r'TiCl3|titanium.*III', 'TiCl3', 'Titanium(III) chloride', 'Reduces nitro→hydroxylamine; McFadyen-Stevens'),
            (r'CrSO4|chromous.*sulfate|Cr\(II\)', 'Cr(II)', 'Chromous sulfate', 'Very strong reductant: alkyl halides → radicals'),
            (r'NADH|NADPH|reduced.*nicotinamide', 'NADH/NADPH', 'Biological reductant', 'Biological hydride donor: reduces C=O, C=C'),
            (r'FADH2', 'FADH2', 'Reduced flavin', 'Biological reductant: electron carrier'),
            (r'ferredoxin|cytochrome.*reduced', 'Ferredoxin/Cyt red', 'Biological electron carriers', 'Chain of biological redox carriers'),
            (r'samarium.*iodide|SmI2', 'SmI2', 'Samarium(II) iodide', 'Single-electron reductant: pinacol coupling, ketone→alcohol'),
            (r'Zn.*Cu.*couple|Zn(Cu)', 'Zn-Cu couple', 'Clemmensen-type', 'Reduces α-halo ketones; Simmons-Smith precursor'),
            (r'PPh3.*reduc|triphenylphosphine.*reduc', 'PPh3', 'Triphenylphosphine', 'Reduces peroxides, disulfides, azides (Staudinger)'),
            (r'Na2S2O4|dithionite|hydrosulfite', 'Na2S2O4', 'Sodium dithionite', 'Reduces quinones→hydroquinones, azo compounds'),
            (r'thiourea.dioxide|H2NC(SH)NH2', 'Thiourea dioxide', 'Thiourea dioxide', 'Alternative dithionite-like reductant'),
            (r'ascorbic.acid|vitamin.C', 'Ascorbic acid', 'Ascorbic acid (vitamin C)', 'Mild biological reductant: antioxidant'),
            (r'glucose|reducing.sugar', 'Glucose', 'Reducing sugars', 'Aldehyde end can reduce Tollens/Fehling/Benedict reagents'),
            (r'CO.*reduc|carbon.monoxide.*reduc', 'CO', 'Carbon monoxide', 'Reducing agent in metallurgy (blast furnace)'),
            (r'C.*reduc|coke.*reduc|charcoal', 'C (coke/charcoal)', 'Carbon', 'Reducing in metallurgy: extracts metals from ores'),
            (r'H2S.*reduc|hydrogen.sulfide.*reduc', 'H2S', 'Hydrogen sulfide', 'Can reduce halogens, some metal ions'),
            (r'SO2.*reduc|sulfur.dioxide.*reduc', 'SO2', 'Sulfur dioxide', 'Reduces dichromate, permanganate, Fe(III)'),
            (r'I-|iodide', 'I⁻', 'Iodide ion', 'Reduces: peroxides, Fe3+, Cu2+, H2O2, etc.'),
            (r'organometallic|R.Mg|R.Li|R2Cu', 'Organometallic', 'Grignard/organolithium', 'Strong nucleophiles/electron donors (R:⁻ character)'),
            (r'photochem.*reduc|photoreduc', 'Photochemical', 'Photoreduction', 'Light-induced electron transfer (e.g., Ru(bpy)3²+)'),
        ]

        for pattern, name, desc, detail in patterns:
            if re.search(pattern, rxn, re.IGNORECASE):
                if not any(m["agent"] == name for m in matches):
                    matches.append({
                        "agent": name,
                        "oxidation_change": "varies",
                        "description": f"{desc}: {detail}",
                    })

        return matches

    def _detect_organic_reduction(self, rxn: str) -> list:
        """Detect organic functional group reduction patterns."""
        patterns = []
        rxn_lower = rxn.lower()

        # Aldehyde/ketone → alcohol
        if re.search(r'aldehyde|ketone|C=O|methanal|ethanal|acetone|benzaldehyde', rxn_lower):
            if re.search(r'alcohol|ol|ethanol|propanol|butanol|OH', rxn_lower):
                patterns.append({"pattern": "carbonyl (aldehyde/ketone) → alcohol", "electrons_gained": 2})

        # Carboxylic acid/ester → alcohol
        if re.search(r'carboxylic|acid|ester|COOR|COOH', rxn_lower):
            if re.search(r'alcohol|ol', rxn_lower):
                patterns.append({"pattern": "carboxylic derivative → alcohol", "electrons_gained": 4})

        # Alkene/alkyne → alkane
        if re.search(r'alkene|alkyne|=C|#C|ethylene|acetylene', rxn_lower):
            if re.search(r'alkane|ethane|propane|butane|saturated', rxn_lower):
                patterns.append({"pattern": "unsaturated → saturated (hydrogenation)", "electrons_gained": 2 or 4})

        # Nitro → amine
        if re.search(r'nitro|NO2|nitrobenzene', rxn_lower):
            if re.search(r'amine|aniline|NH2|amino', rxn_lower):
                patterns.append({"pattern": "nitro → amine", "electrons_gained": 6})

        # Nitrile → amine/aldehyde
        if re.search(r'nitrile|CN|cyanide', rxn_lower):
            if re.search(r'amine|NH2', rxn_lower):
                patterns.append({"pattern": "nitrile → amine", "electrons_gained": 4})
            elif re.search(r'aldehyde|formyl', rxn_lower):
                patterns.append({"pattern": "nitrile → aldehyde (partial reduction)", "electrons_gained": 2})

        # Epoxide → alcohol
        if re.search(r'epoxide|oxirane|ethoxyirane', rxn_lower):
            if re.search(r'alcohol|diol|OH', rxn_lower):
                patterns.append({"pattern": "epoxide → alcohol (ring opening)", "electrons_gained": 2})

        # Halide → alkane (dehalogenation)
        if re.search(r'halide|bromide|chloride|iodide|X', rxn_lower):
            if re.search(r'alkane|H|dehalogen', rxn_lower):
                patterns.append({"pattern": "halide → alkane (dehalogenation)", "electrons_gained": 2})

        # Disulfide → thiol
        if re.search(r'disulfide|S-S|RSSR', rxn_lower):
            if re.search(r'thiol|SH|mercaptan', rxn_lower):
                patterns.append({"pattern": "disulfide → thiol", "electrons_gained": 2})

        # Quinone → hydroquinone
        if re.search(r'quinone|cyclohexadienedione', rxn_lower):
            if re.search(r'hydroquinone|dihydroxybenzene', rxn_lower):
                patterns.append({"pattern": "quinone → hydroquinone", "electrons_gained": 2})

        return patterns

    def _parse_oxidation_changes(self, rxn: str) -> list:
        """Parse oxidation state changes — focus on species being OXIDIZED (losing e⁻)."""
        changes = []

        couples = [
            (r'Na\b|Na(?!\w)|sodium', r'Na[+]|NaOH|NaCl|Na2SO4', 'Na', '0', '+1', +1),
            (r'K\b|potassium', r'K[+]|KOH|KCl', 'K', '0', '+1', +1),
            (r'Zn\b|zinc', r'Zn[+2]|ZnCl2|ZnSO4', 'Zn', '0', '+2', +2),
            (r'Fe\b.*metal|Fe\(s\)|iron.*metal', r'Fe[+2]|FeSO4|FeCl2', 'Fe', '0', '+2', +2),
            (r'Mg\b|magnesium', r'Mg[+2]|MgCl2|MgSO4', 'Mg', '0', '+2', +2),
            (r'Al\b|aluminum', r'Al[+3]|AlCl3|Al2(SO4)3', 'Al', '0', '+3', +3),
            (r'SnCl2|Sn\(II\)|stannous', r'Sn[+4]|SnCl4', 'Sn', '+2', '+4', +2),
            (r'FeSO4|Fe[+2]|ferrous', r'Fe[+3]|FeCl3|Fe2(SO4)3', 'Fe', '+2', '+3', +1),
            (r'H2\b', r'H2O|H+', 'H', '0', '+1', +1),
            (r'C\b.*reduc|CO\b', r'CO2|CO3', 'C', 'varies', '+IV', 'varies'),
            (r'SO2|sulfur.dioxide', r'SO3|SO4|sulfate', 'S', '+IV', '+VI', +2),
            (r'H2S|hydrogen.sulfide', r'S|SO2', 'S', '-II', '0 or +IV', 'varies'),
            (r'I[-−]|iodide', r'I2', 'I', '-I', '0', +1),
            (r'Br[-−]|bromide', r'Br2', 'Br', '-I', '0', +1),
            (r'organometallic|R-Mg|R-Li', r'R-H|R-OH', 'C (in R)', '~ -III to ~varies', 'higher', 'varies'),
            (r'N2H4|hydrazine', r'N2|NH3', 'N', '-II', '0 or -III', 'varies'),
            (r'NaBH4|borohydride', r'B(OH)3|borate', 'H (in BH4)', '-I', '+I', +2),
            (r'LiAlH4', r'Al(OH)3|Al3+', 'H (in AlH4)', '-I', '+I', +2),
            (r'C2O4(2-)|oxalate', r'CO2', 'C', '+III', '+IV', +1),
            (r'glucose|aldose', r'carboxylic|acid|lactone', 'C (anomeric)', '+I', '+III or +IV', '+2 to +6'),
        ]

        for rp, pp, elem, f_os, t_os, e_ch in couples:
            if re.search(rp, rxn, re.IGNORECASE) and re.search(pp, rxn, re.IGNORECASE):
                changes.append({
                    "element": elem,
                    "from_oxidation_state": f_os,
                    "to_oxidation_state": t_os,
                    "change": e_ch,
                    "is_oxidation": isinstance(e_ch, (int, float)) and e_ch > 0,
                    "is_reduction": isinstance(e_ch, (int, float)) and e_ch < 0,
                })

        return changes

    def _determine_primary_source(self, matched_sources, organic_patterns, ox_changes):
        """Determine the primary electron source."""
        strong_reductants = ['Na', 'K', 'Li', 'LiAlH4', 'Cr(II)', 'CrSO4',
                             'SmI2', 'metal/acid (Zn/HCl, Fe/HCl, Sn/HCl)',
                             'dissolving metal (Na/NH3)']

        for src_info in matched_sources:
            agent = src_info.get("agent", "")
            for strong in strong_reductants:
                if strong.lower() in agent.lower():
                    return {
                        "name": agent,
                        "species": [agent],
                        "reasoning": src_info.get("description", f"{strong} is a strong reducing agent that donates electrons."),
                        "confidence": "high",
                        "type": "strong_reductant",
                    }

        if matched_sources:
            best = matched_sources[0]
            return {
                "name": best["agent"],
                "species": [best["agent"]],
                "reasoning": best.get("description", "Identified as a reducing agent."),
                "confidence": "medium-high",
                "type": "known_reductant",
            }

        oxidations = [c for c in ox_changes if c.get("is_oxidation")]
        if oxidations:
            elem = oxidations[0]["element"]
            return {
                "name": f"{elem}-containing species (oxidized)",
                "species": [f"{elem} species"],
                "reasoning": f"Oxidation state of {elem} increases from {oxidations[0]['from_oxidation_state']} to {oxidations[0]['to_oxidation_state']}, indicating it loses electrons.",
                "confidence": "medium",
                "type": "oxidation_detected",
            }

        if organic_patterns:
            return {
                "name": "reducing agent (implied)",
                "species": ["[H] / reductant"],
                "reasoning": f"Organic reduction pattern(s) detected: {[p['pattern'] for p in organic_patterns]}. A reducing agent must be present to donate these electrons.",
                "confidence": "medium-low",
                "type": "inferred_from_organic_reduction",
            }

        return {
            "name": "undetermined",
            "species": [],
            "reasoning": "Could not identify a clear electron source. Provide more specific reaction details.",
            "confidence": "low",
            "type": "unknown",
        }

    def _count_electrons(self, ox_changes, matched_sources):
        total = None
        detail = {}
        oxidations = [c for c in ox_changes if c.get("is_oxidation")]
        reductions = [c for c in ox_changes if c.get("is_reduction")]

        if oxidations and reductions:
            ox_e = sum(c["change"] for c in oxidations if isinstance(c["change"], (int, float)))
            red_e = sum(abs(c["change"]) for c in reductions if isinstance(c["change"], (int, float)))
            total = max(ox_e, red_e)
            detail = {"electrons_lost_by_source": ox_e, "electrons_gained_by_sink": red_e, "total": total}
        elif oxidations:
            total = sum(c["change"] for c in oxidations if isinstance(c["change"], (int, float)))
            detail = {"electrons_lost_by_source": total, "total": total}
        else:
            detail = {"note": "Insufficient data.", "total": None}
        return detail

    def _describe_electron_flow(self, primary_source, ox_changes):
        name = primary_source.get("name", "?")
        reasoning = primary_source.get("reasoning", "")
        parts = [f"Primary electron source: **{name}**.", reasoning]
        oxidations = [c for c in ox_changes if c.get("is_oxidation")]
        reductions = [c for c in ox_changes if c.get("is_reduction")]
        if oxidations:
            sources = [f"{c['element']} ({c['from_oxidation_state']} → {c['to_oxidation_state']})" for c in oxidations]
            parts.append(f"Species oxidized (lose e⁻): {', '.join(sources)}.")
        if reductions:
            sinks = [f"{c['element']}" for c in reductions]
            parts.append(f"Electron sink(s): {', '.join(sinks)} — these gain electrons.")
        return " ".join(parts)

    def _propose_half_reactions(self, primary_source, ox_changes):
        half_rxns = []
        name = primary_source.get("name", "")

        # Oxidation half-reaction (source)
        agent_map = {
            'NaBH4': ("BH4⁻ + 8OH⁻ → BO2⁻ + 6H2O + 8e⁻", 8, "or simplified: H⁻ → ½H2 + e⁻"),
            'LiAlH4': ("AlH4⁻ + 4OH⁻ → Al(OH)4⁻ + 4H2↑ + 4e⁻", 4, "hydride donation"),
            'H2/metal catalyst': ("H2 → 2H⁺ + 2e⁻", 2, "heterogeneous catalysis"),
            'Na/NH3(l)': ("Na → Na⁺ + e⁻", 1, "dissolving metal: dissolved e⁻(solv)"),
            'Zn/HCl': ("Zn → Zn²⁺ + 2e⁻", 2, "acidic conditions"),
            'Sn/HCl': ("Sn → Sn²⁺ + 2e⁻", 2, "or Sn²⁺ → Sn⁴⁺ + 2e⁻"),
            'Fe/HCl': ("Fe → Fe²⁺ + 2e⁻", 2, "iron in acid"),
            'Zn': ("Zn → Zn²⁺ + 2e⁻", 2, "direct oxidation"),
            'Na': ("Na → Na⁺ + e⁻", 1, "alkali metal"),
            'K': ("K → K⁺ + e⁻", 1, "alkali metal"),
            'I⁻': ("2I⁻ → I₂ + 2e⁻", 2, "halide oxidation"),
            'SO2': ("SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻", 2, "sulfite oxidation"),
            'H2C2O4': ("C2O4²⁻ → 2CO2 + 2e⁻", 2, "oxalate oxidation"),
            'H2S': ("H2S + 4H2O → SO4²⁻ + 10H⁺ + 8e⁻", 8, "full oxidation"),
            'N2H4': ("N2H4 → N2 + 4H⁺ + 4e⁻", 4, "hydrazine oxidation"),
            'glucose': ("C6H12O6 + 3H2O → 2CH3COOH + 4H⁺ + 8e⁻", 8, "simplified glucose oxidation"),
        }

        for key, val in agent_map.items():
            if key.lower() in name.lower():
                half_rxns.append({
                    "type": "oxidation (electron source)",
                    "equation": val[0],
                    "electrons": val[1],
                    "note": val[2],
                })
                break
        else:
            half_rxns.append({
                "type": "oxidation (electron source)",
                "equation": f"(Half-reaction for {name} — provide specifics)",
                "electrons": "unknown",
            })

        if ox_changes:
            for c in ox_changes:
                if c.get("is_reduction"):
                    half_rxns.append({
                        "type": "reduction (electron sink)",
                        "equation": f"{c['element']} ({c['from_oxidation_state']}) → {c['to_oxidation_state']} + {abs(c['change'])}e⁻",
                        "electrons": abs(c["change"]) if isinstance(c["change"], (int, float)) else "unknown",
                    })

        return half_rxns

    def _generate_notes(self, primary_source, matched_sources, organic_patterns):
        notes = []
        conf = primary_source.get("confidence", "")
        if conf == "low":
            notes.append("⚠️ Low confidence: Provide a more detailed reaction equation including all reactants and products.")
        if len(matched_sources) > 1:
            notes.append(f"ℹ️ Multiple potential reductants detected ({len(matched_sources)}). Strongest identified as primary source.")
        if organic_patterns and not matched_sources:
            notes.append("📌 Organic reduction detected but no specific reductant named.")
            notes.append("💡 For carbonyl reductions: NaBH4 (mild, aldehydes/ketones), LiAlH4 (strong, esters/acids too).")
            notes.append("💡 For nitro reductions: H2/Pd (clean), Fe/HCl (classic industrial), SnCl2 (stops at hydroxylamine).")
            notes.append("💡 For alkene/alkyne reductions: H2/Pd (→ alkane), Lindlar (→ cis-alkene), Na/NH3 (→ trans-alkene).")
        t = primary_source.get("type", "")
        if t == "inferred_from_organic_reduction":
            notes.append("🔬 Consider the substrate sensitivity: chemoselectivity matters when multiple reducible groups present.")
        return notes

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        rxn_input = parts[0] if parts else ""
        mode = parts[1] if len(parts) > 1 else "detailed"
        return self._run_base(rxn_input, mode)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
