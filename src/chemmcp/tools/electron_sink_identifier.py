"""
Electron Sink Identifier (Tool #151)
识别反应中的电子受体（电子阱）：分析反应中氧化数升高的物种，
确定哪个物种接受电子，分析电子流向和氧化还原过程。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# 常见氧化剂及其得电子特征
_OXIDIZING_AGENTS = {
    # 分子/离子 -> (典型氧化态变化, 描述)
    'O2': ('0 → -2', 'Oxygen: reduced to oxide/water (O₂ → O²⁻ or H₂O)'),
    'H2O2': ('-1 → -2', 'Hydrogen peroxide: each O goes from -1 to -2'),
    'KMnO4': ('Mn(+7) → +2/+4/+6', 'Permanganate: Mn(VII) reduced, color purple→colorless/brown'),
    'K2Cr2O7': ('Cr(+6) → +3', 'Dichromate: Cr(VI) → Cr(III), orange→green'),
    'MnO2': ('Mn(+4) → +2', 'Manganese dioxide: Mn(IV) → Mn(II)'),
    'CrO3': ('Cr(+6) → +3', 'Chromium trioxide: strong oxidant'),
    'HNO3': ('N(+5) → +2/+4', 'Nitric acid: NO₃⁻ → NO/NO₂ depending on concentration'),
    'H2SO4(conc)': ('S(+6) → +4', 'Concentrated sulfuric acid: SO₄²⁻ → SO₂'),
    'FeCl3': ('Fe(+3) → +2', 'Iron(III): Fe³⁺ → Fe²⁺'),
    'Ce(SO4)2': ('Ce(+4) → +3', 'Ceric ion: Ce(IV) → Ce(III)'),
    'NaOCl': ('Cl(+1) → -1', 'Bleach: ClO⁻ → Cl⁻'),
    'Cl2': ('0 → -1', 'Chlorine: Cl₂ → 2Cl⁻'),
    'Br2': ('0 → -1', 'Bromine: Br₂ → 2Br⁻'),
    'I2': ('0 → -1', 'Iodine: I₂ → 2I⁻'),
    'F2': ('0 → -1', 'Fluorine: strongest elemental oxidant'),
    'PbO2': ('Pb(+4) → +2', 'Lead dioxide: Pb(IV) → Pb(II)'),
    'NaBO3': ('O(-1/2) → -2', 'Sodium perborate: peroxide oxygen reduced'),
    'OXONE': ('O(-1/2) → -1', 'Potassium peroxymonosulfate'),
    'NBS': ('Br(+1) → -1', 'N-Bromosuccinimide: brominating/oxidizing agent'),
    'NCS': ('Cl(+1) → -1', 'N-Chlorosuccinimide'),
    'DDQ': ('quinone → hydroquinone', 'Dichlorodicyanobenzoquinone: accepts 2e⁻/2H⁺'),
    'CAN': ('Ce(+4) → +3', 'Ceric ammonium nitrate'),
    'IBX': ('I(+3) → +1/-1', '2-Iodoxybenzoic acid: hypervalent iodine oxidant'),
    'Dess-Martin': ('I(+3) → +1', 'Dess-Martin periodinane: alcohol oxidation'),
    'Swern_oxidant': ('DMSO → DMS', 'Swern oxidation: DMSO reduced by electrons'),
    'TPAP': ('Ru(+7) → +5', 'Tetra-n-propylammonium perruthenate'),
    'PCC': ('Cr(+6) → +3', 'Pyridinium chlorochromate'),
    'PDC': ('Cr(+6) → +3', 'Pyridinium dichromate'),
    'Ag2O': ('Ag(+1) → 0', 'Silver(I) oxide: Ag⁺ → Ag(s)'),
    'CuCl2': ('Cu(+2) → +1/0', 'Copper(II): Cu²⁺ → Cu⁺/Cu'),
    'Fehling': ('Cu(+2) → +1', "Fehling's solution: Cu²⁺ → Cu₂O (red ppt)"),
    'Benedict': ('Cu(+2) → +1', "Benedict's solution: similar to Fehling"),
    'Tollens': ('Ag(+1) → 0', "Tollens' reagent: silver mirror test"),
    'SeO2': ('Se(+4) → 0', 'Selenium dioxide: Se(IV) → Se(0), allylic oxidant'),
    'O3': ('0 → -2', 'Ozone: ozonolysis oxidant'),
    'mCPBA': ('O(-1) → -2', 'meta-Chloroperoxybenzoic acid: epoxidizing agent'),
    'peracid': ('O(-1) → -2', 'Peracids (RCO₃H): epoxidation/Baeyer-Villiger'),
    'NMO': ('N(-1) → -3', 'N-Methylmorpholine N-oxide: co-oxidant'),
    'quat_amine_oxide': ('N(+1) → -3', 'Amine N-oxide: [O] transfer reagent'),
    'PhI(OAc)2': ('I(+3) → +1', '(Diacetoxyiodo)benzene: hypervalent iodine'),
    'IBX': ('I(+3) → -1', '2-Iodoxybenzoic acid'),
    'HIO4': ('I(+7) → +5/+3', 'Periodic acid: glycol cleavage'),
}

# 官能团氧化特征（被氧化的基团 -> 氧化产物）
_FUNCTIONAL_GROUP_OXIDATION = {
    'alcohol_primary': {'oxidation_state_change': 'C(-I) → 0 → +II', 'product': 'aldehyde → carboxylic acid', 'sink': 'oxidizing_agent'},
    'alcohol_secondary': {'oxidation_state_change': 'C(-I) → 0', 'product': 'ketone', 'sink': 'oxidizing_agent'},
    'alcohol_tertiary': {'note': 'No α-H, resists oxidation (except C-C bond cleavage)'},
    'aldehyde': {'oxidation_state_change': 'C(+I) → +III', 'product': 'carboxylic acid', 'sink': 'oxidizing_agent'},
    'alkene': {'oxidation_state_change': 'C(-II) → -I/0/+III', 'product': 'diol/epoxide/cleaved carbonyl', 'sink': 'oxidizing_agent'},
    'alkyne': {'oxidation_state_change': 'C(-I) → +III', 'product': 'dicarboxylic acid/α-diketone', 'sink': 'oxidizing_agent'},
    'arene_side_chain': {'oxidation_state_change': 'benzylic C → C=O', 'product': 'benzoic acid (regardless of chain length)', 'sink': 'KMnO4/Na2Cr2O7'},
    'phenol': {'oxidation_state_change': 'aromatic ring', 'product': 'quinone/polymer', 'sink': 'oxidizing_agent'},
    'thiol': {'oxidation_state_change': 'S(-II) → -I/0', 'product': 'disulfide/sulfonic acid', 'sink': 'mild/strong oxidant'},
    'thioether': {'oxidation_state_change': 'S(-II) → 0/+IV/+VI', 'product': 'sulfoxide/sulfone', 'sink': 'H2O2/mCPBA/peracid'},
    'amine': {'oxidation_state_change': 'N(-III) → -II/-I/0/+I', 'product': 'hydroxylamine/nitroso/nitro', 'sink': 'oxidizing_agent'},
    'phosphine': {'oxidation_state_change': 'P(-III) → +V', 'product': 'phosphine oxide', 'sink': 'O2/H2O2'},
}


@ChemMCPManager.register_tool
class ElectronSinkIdentifier(BaseTool):
    __version__ = "0.1.0"
    name = "ElectronSinkIdentifier"
    func_name = 'identify_electron_sink'
    description = "Identify the electron sink (electron acceptor/oxidizing agent) in a chemical reaction. Analyze reaction SMILES, chemical equations, or redox processes to determine which species accepts electrons and gets reduced."
    implementation_description = "Uses pattern matching against a comprehensive database of known oxidizing agents, functional group oxidation signatures, and oxidation state rules. Parses reaction SMILES/equations to identify species that decrease in oxidation state (gain electrons). Supports both inorganic and organic redox reactions."
    categories = ["Reaction"]
    tags = ["Electron Transfer", "Redox", "Oxidation", "Reaction Analysis", "Oxidizing Agent"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("reaction_input", "str", "N/A", "Reaction SMILES string, chemical equation, or description of the redox process."),
        ("analysis_mode", "str", "detailed", "Analysis mode: 'brief' for summary, 'detailed' for full analysis with electron counting."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: reaction_input [analysis_mode]. E.g., 'CC(=O)O>[O]>CC(=O)O detailed'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing sink_species, sink_reasoning, electron_flow, oxidation_state_changes, confidence, and analysis_details."),
    ]

    examples = [
        {
            "code_input": {
                "reaction_input": "2KMnO4 + 16HCl → 2KCl + 2MnCl2 + 5Cl2 + 8H2O",
                "analysis_mode": "detailed",
            },
            "text_input": {"query": "2KMnO4+16HCl→2KCl+2MnCl2+5Cl2+8H2O detailed"},
            "output": {
                "result": {
                    "primary_electron_sink": "KMnO4 (permanganate)",
                    "sink_species": ["MnO4-"],
                    "sink_reasoning": "Mn oxidation state decreases from +7 to +2 (gains 5 e⁻ per Mn atom). In acidic medium, permanganate is a powerful oxidizing agent.",
                    "electrons_transferred_per_reaction": 10,
                    "half_reaction_reduction": "MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O",
                    "electron_flow_description": "Chloride ions (Cl⁻) donate electrons to permanganate (MnO4⁻). Each Mn(VII) accepts 5 electrons to become Mn(II).",
                    "oxidation_state_changes": [
                        {"species": "Mn", "from": "+7", "to": "+2", "change": "-5", "electrons_gained": 5},
                        {"species": "Cl", "from": "-1", "to": "0", "change": "+1", "electrons_lost": 1},
                    ],
                    "confidence": "high",
                    "reaction_type": "redox (disproportionation of HCl / reduction of Mn)",
                }
            },
        },
        {
            "code_input": {
                "reaction_input": "CH3CH2OH + [O] → CH3CHO + H2O",
                "analysis_mode": "brief",
            },
            "text_input": {"query": "CH3CH2OH+[O]→CH3CHO+H2O brief"},
            "output": {
                "result": {
                    "primary_electron_sink": "[O] (oxidizing agent)",
                    "sink_species": ["oxidizing agent"],
                    "sink_reasoning": "Primary alcohol oxidized to aldehyde; the oxidizing agent (e.g., PCC, Swern, Dess-Martin) accepts 2 electrons from the alcohol carbon.",
                    "electrons_transferred": 2,
                    "oxidation_state_changes": [
                        {"species": "C (alcohol)", "from": "-I", "to": "+I", "change": "+2"},
                    ],
                    "confidence": "high",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize the oxidizing agent database."""
        self.oxidizing_agents = dict(_OXIDIZING_AGENTS)
        self.fg_oxidation = dict(_FUNCTIONAL_GROUP_OXIDATION)

    def _run_base(self, reaction_input: str, analysis_mode: str = "detailed") -> dict:
        """
        Core logic: identify electron sinks in a reaction.
        """
        if not reaction_input or not reaction_input.strip():
            raise ChemMCPInputError("Reaction input cannot be empty.")

        rxn = reaction_input.strip()
        mode = analysis_mode.lower() if analysis_mode else "detailed"

        # Step 1: Try to match known oxidizing agents
        matched_sinks = self._match_oxidizing_agents(rxn)

        # Step 2: Detect organic oxidation patterns
        organic_patterns = self._detect_organic_oxidation(rxn)

        # Step 3: Parse oxidation state changes from equation
        ox_changes = self._parse_oxidation_changes(rxn)

        # Step 4: Determine primary electron sink
        primary_sink = self._determine_primary_sink(matched_sinks, organic_patterns, ox_changes)

        # Step 5: Count electrons transferred
        e_count = self._count_electrons(ox_changes, matched_sinks)

        # Build result
        result = {
            "result": {
                "reaction_input": rxn,
                "primary_electron_sink": primary_sink.get("name", "unknown"),
                "sink_species": primary_sink.get("species", []),
                "sink_reasoning": primary_sink.get("reasoning", "Could not determine electron sink definitively."),
                "matched_oxidizing_agents": matched_sinks,
                "organic_oxidation_patterns": organic_patterns,
                "oxidation_state_changes": ox_changes,
                "electrons_transferred": e_count.get("total", None),
                "electrons_transferred_detail": e_count,
                "electron_flow_description": self._describe_electron_flow(primary_sink, ox_changes),
                "confidence": primary_sink.get("confidence", "low"),
                "analysis_mode": mode,
                "additional_notes": self._generate_notes(primary_sink, matched_sinks, organic_patterns),
            }
        }

        if mode == "detailed":
            result["result"]["half_reactions"] = self._propose_half_reactions(primary_sink, ox_changes)

        logger.info(f"ElectronSink: identified {primary_sink.get('name','?')} as electron sink in: {rxn[:50]}")
        return result

    def _match_oxidizing_agents(self, rxn: str) -> list:
        """Match known oxidizing agents in the reaction string."""
        matches = []
        rxn_upper = rxn.upper()

        # Direct formula matches
        for agent, info in self.oxidizing_agents.items():
            if agent.upper() in rxn_upper or agent.replace('(', '').replace(')', '') in rxn_upper:
                matches.append({
                    "agent": agent,
                    "oxidation_change": info[0],
                    "description": info[1],
                })

        # Pattern-based matches
        patterns = [
            (r'\[O\]', '[O]', 'Generic oxidant notation ([O])', 'Accepts 2H + 2e⁻ typically'),
            (r'oxid|OXID', 'oxidation', 'Oxidation indicated', 'General oxidation process'),
            (r'MnO4[-−]?', 'MnO4⁻', 'Permanganate', 'Strong oxidant, Mn(VII)→lower'),
            (r'Cr2O7[^)]*', 'Cr2O7²⁻', 'Dichromate', 'Strong oxidant, Cr(VI)→Cr(III)'),
            (r'HNO3', 'HNO3', 'Nitric acid', 'Oxidizing acid'),
            (r'H2SO4.*conc|conc.*H2SO4', 'conc. H2SO4', 'Conc. sulfuric acid', 'Hot conc. H2SO4 is oxidizing'),
            (r'Ce\(?IV\)?|Ce\(?4\+\)?|Ce\(SO', 'Ce(IV)', 'Ceric ion', 'Ce(IV)→Ce(III)'),
            (r'NAD[\+]?|NADP[\+]?', 'NAD(P)+', 'Nicotinamide cofactor', 'Biological oxidant: NAD+→NADH'),
            (r'FAD|FMN', 'FAD/FMN', 'Flavin cofactor', 'Biological oxidant: FAD→FADH2'),
            (r'O2', 'O2', 'Molecular oxygen', 'Terminal biological oxidant'),
            (r'DMP|Dess.Martin', 'Dess-Martin periodinane', 'Alcohol oxidation', 'Hypervalent I(III)→I(I)'),
            (r'Swern|DMSO.*oxalyl|(COCl)2.*DMSO', 'Swern conditions', 'Swern oxidation', 'DMSO→DMS (reduced)'),
            (r'PCC', 'PCC', 'PCC', 'Cr(VI) oxidant for alcohols'),
            (r'TPCA|TPAP', 'TPAP', 'Ruthenium oxidant', 'Ru(VII)→Ru(V)'),
            (r'mCPBA|MMPP|peracid', 'peracid/mCPBA', 'Peroxide oxidant', 'Epoxidation/Baeyer-Villiger'),
            (r'SeO2', 'SeO2', 'Selenium dioxide', 'Allylic/benzylic oxidation'),
            (r'O3', 'O3', 'Ozone', 'Ozonolysis oxidant'),
            (r'IBX', 'IBX', 'IBX', 'Hypervalent iodine oxidant'),
            (r'Ag2O|Tollens', 'Ag(I)/Tollens', 'Silver-based oxidant', 'Aldehyde→carboxylic acid, Ag mirror'),
            (r'Cu\(II\)|CuCl2|CuSO4.*oxid', 'Cu(II)', 'Copper(II) oxidant', 'Cu(II)→Cu(I) or Cu(0)'),
            (r'CAN|ceric.*nitrate', 'CAN', 'Ceric ammonium nitrate', 'Ce(IV) oxidant'),
            (r'DDQ', 'DDQ', 'DDQ (quinone)', 'Dehydrogenating agent, accepts 2H'),
            (r'activated.*MnO2|MnO2', 'MnO2', 'Activated MnO2', 'Mild oxidant for allylic/benzylic alcohols'),
            (r'Baeyer|cold.*KMnO4|dilute.*KMnO4', 'Baeyer test', 'Cold dilute KMnO4', 'Alkene→diol (purple→colorless)'),
            (r'I2.*NaOH|Iodoform', 'Iodoform test', 'I2/NaOH (iodoform)', 'Methyl ketone→CHI3 (yellow ppt)'),
        ]

        for pattern, name, desc, detail in patterns:
            if re.search(pattern, rxn, re.IGNORECASE):
                # Avoid duplicates
                if not any(m["agent"] == name for m in matches):
                    matches.append({
                        "agent": name,
                        "oxidation_change": "varies",
                        "description": f"{desc}: {detail}",
                    })

        return matches

    def _detect_organic_oxidation(self, rxn: str) -> list:
        """Detect organic functional group oxidation patterns."""
        patterns = []
        rxn_lower = rxn.lower()

        # Alcohol → aldehyde/carboxylic acid
        if re.search(r'OH|alcohol|ethanol|propanol|butanol', rxn_lower):
            if re.search(r'CHO|aldehyde|methanal|ethanal', rxn_lower):
                patterns.append({"pattern": "alcohol → aldehyde", "electrons_lost": 2})
            elif re.search(r'COOH|carboxylic|acid|oic', rxn_lower):
                patterns.append({"pattern": "alcohol → carboxylic acid", "electrons_lost": 4})

        # Aldehyde → carboxylic acid
        if re.search(r'CHO|aldehyde', rxn_lower) and re.search(r'COOH|carboxylic|acid', rxn_lower):
            patterns.append({"pattern": "aldehyde → carboxylic acid", "electrons_lost": 2})

        # Alkene → diol/epoxide/cleavage
        if re.search(r'C=C|=C|alkene|ethylene|propene', rxn_lower):
            if re.search(r'OH.*OH|diol|glycol', rxn_lower):
                patterns.append({"pattern": "alkene → vicinal diol", "electrons_lost": 2})
            elif re.search(r'epox|oxirane|oxir', rxn_lower):
                patterns.append({"pattern": "alkene → epoxide", "electrons_lost": 2})
            elif re.search(r'C=O|ketone|aldehyde|cleav', rxn_lower):
                patterns.append({"pattern": "alkene → oxidative cleavage", "electrons_lost": varies})

        # Thiol → disulfide
        if re.search(r'SH|thiol', rxn_lower) and re.search(r'S-S|disulfide', rxn_lower):
            patterns.append({"pattern": "thiol → disulfide", "electrons_lost": 2})

        # Thioether → sulfoxide/sulfone
        if re.search(r'S-|sulfide|thioether', rxn_lower):
            if re.search(r'S=O|sulfoxide', rxn_lower):
                patterns.append({"pattern": "thioether → sulfoxide", "electrons_lost": 2})

        # Amine → higher oxidation state
        if re.search(r'NH2|amine', rxn_lower):
            if re.search(r'NO|nitroso', rxn_lower):
                patterns.append({"pattern": "amine → nitroso", "electrons_lost": 2})
            elif re.search(r'NO2|nitro', rxn_lower):
                patterns.append({"pattern": "amine → nitro", "electrons_lost": 6})

        # Arene side-chain oxidation
        if re.search(r'benzyl|phenyl.*CH|toluene|ethylbenz', rxn_lower) and re.search(r'COOH|benzoic', rxn_lower):
            patterns.append({"pattern": "alkylbenzene → benzoic acid", "electrons_lost": "many"})

        # SMILES-based detection
        if '=' in rxn and ('O' in rxn or '[' in rxn):
            # Check for common SMILES oxidation patterns
            if re.search(r'CCO', rxn) and re.search(r'CC=O', rxn):
                patterns.append({"pattern": "SMILES: alcohol → aldehyde detected", "electrons_lost": 2})

        return patterns

    def _parse_oxidation_changes(self, rxn: str) -> list:
        """Attempt to parse oxidation state changes from the reaction."""
        changes = []

        # Common redox couples
        redox_couples = [
            # (reactant_pattern, product_pattern, element, from_os, to_os, e_change)
            (r'MnO4[-−]?|Mn\(?VII\)?', r'Mn[+2]|MnCl2|MnSO4', 'Mn', '+7', '+2', -5),
            (r'Cr2O7[^)]*|Cr\(?VI\)?|CrO3', r'Cr[+3]|CrCl3', 'Cr', '+6', '+3', -3),
            (r'Fe[+3]|FeCl3|Fe2(SO4)3', r'Fe[+2]|FeCl2|FeSO4', 'Fe', '+3', '+2', -1),
            (r'Cu[+2]|CuCl2|CuSO4', r'Cu[+]|Cu2O|Cu', 'Cu', '+2', '+1', -1),
            (r'Ag[+]|AgNO3|Ag2O', r'Ag\(s\)|Ag', 'Ag', '+1', '0', -1),
            (r'Ce[+4]|Ce\(?IV\)?', r'Ce[+3]|Ce\(?III\)?', 'Ce', '+4', '+3', -1),
            (r'Cl[-−]?|Cl2', r'Cl2', 'Cl', '-1', '0', +1),
            (r'Br[-−]?|Br2', r'Br2', 'Br', '-1', '0', +1),
            (r'I[-−]?|I2', r'I2', 'I', '-1', '0', +1),
            (r'S[−-][−-]?', r'SO4|SO3|SO2', 'S', '-2', '+4/+6', '+6 to +8'),
            (r'N[−-]?H3?|NH3', r'N2|NO|NO2|N2O', 'N', '-3', 'varies', 'varies'),
            (r'H2', r'H2O|H+', 'H', '0', '+1', +1),
            (r'C\b.*metal|organometallic|R-M', r'C-H|C-OH|C=O', 'C', 'varies', 'varies', 'varies'),
        ]

        for rp, pp, elem, f_os, t_os, e_ch in redox_couples:
            if re.search(rp, rxn, re.IGNORECASE) and re.search(pp, rxn, re.IGNORECASE):
                changes.append({
                    "element": elem,
                    "from_oxidation_state": f_os,
                    "to_oxidation_state": t_os,
                    "change": e_ch,
                    "is_reduction": isinstance(e_ch, (int, float)) and e_ch < 0,
                    "is_oxidation": isinstance(e_ch, (int, float)) and e_ch > 0,
                })

        return changes

    def _determine_primary_sink(self, matched_agents, organic_patterns, ox_changes):
        """Determine the primary electron sink from all evidence."""
        # Priority 1: Known strong oxidizing agent explicitly present
        strong_oxidants = ['KMnO4', 'K2Cr2O7', 'CrO3', 'NaBO3', 'OXONE', 'F2', 'O3',
                          'CAN', 'mCPBA', 'Mg(MnO4)2', 'Ce(SO4)2']

        for agent_info in matched_agents:
            agent = agent_info.get("agent", "")
            for strong in strong_oxidants:
                if strong.upper() in agent.upper():
                    return {
                        "name": agent,
                        "species": [agent],
                        "reasoning": agent_info.get("description", f"{strong} is a strong oxidizing agent that accepts electrons."),
                        "confidence": "high",
                        "type": "strong_inorganic_oxidant",
                    }

        # Priority 2: Any matched oxidizing agent
        if matched_agents:
            best = matched_agents[0]
            return {
                "name": best["agent"],
                "species": [best["agent"]],
                "reasoning": best.get("description", "Identified as an oxidizing agent in the reaction."),
                "confidence": "medium-high",
                "type": "known_oxidant",
            }

        # Priority 3: Oxidation state changes showing reduction
        reductions = [c for c in ox_changes if c.get("is_reduction")]
        if reductions:
            elem = reductions[0]["element"]
            return {
                "name": f"{elem}-containing species (reduced)",
                "species": [f"{elem} species"],
                "reasoning": f"Oxidation state of {elem} decreases from {reductions[0]['from_oxidation_state']} to {reductions[0]['to_oxidation_state']}, indicating it gains electrons.",
                "confidence": "medium",
                "type": "reduction_detected",
            }

        # Priority 4: Organic oxidation implies presence of oxidant
        if organic_patterns:
            return {
                "name": "oxidizing agent (implied)",
                "species": ["[O] / oxidant"],
                "reasoning": f"Organic oxidation pattern(s) detected: {[p['pattern'] for p in organic_patterns]}. An oxidizing agent must be present to accept these electrons.",
                "confidence": "medium-low",
                "type": "inferred_from_organic_oxidation",
            }

        # Default
        return {
            "name": "undetermined",
            "species": [],
            "reasoning": "Could not identify a clear electron sink from the input. Provide more specific reaction details.",
            "confidence": "low",
            "type": "unknown",
        }

    def _count_electrons(self, ox_changes, matched_agents):
        """Estimate total electrons transferred."""
        total = None
        detail = {}

        reductions = [c for c in ox_changes if c.get("is_reduction")]
        oxidations = [c for c in ox_changes if c.get("is_oxidation")]

        if reductions and oxidations:
            # Sum up electrons
            red_e = sum(abs(c["change"]) for c in reductions if isinstance(c["change"], (int, float)))
            ox_e = sum(c["change"] for c in oxidations if isinstance(c["change"], (int, float)))
            total = max(red_e, ox_e)
            detail = {"electrons_gained_by_sink": red_e, "electrons_lost_by_source": ox_e, "total": total}
        elif reductions:
            total = sum(abs(c["change"]) for c in reductions if isinstance(c["change"], (int, float)))
            detail = {"electrons_gained_by_sink": total, "total": total}
        elif matched_agents:
            detail = {"note": "Electron count based on matched oxidant stoichiometry needed.", "total": None}
        else:
            detail = {"note": "Insufficient data to count electrons.", "total": None}

        return detail

    def _describe_electron_flow(self, primary_sink, ox_changes):
        """Generate human-readable electron flow description."""
        name = primary_sink.get("name", "?")
        reasoning = primary_sink.get("reasoning", "")

        reductions = [c for c in ox_changes if c.get("is_reduction")]
        oxidations = [c for c in ox_changes if c.get("is_oxidation")]

        parts = [f"Primary electron sink: **{name}**."]
        parts.append(reasoning)

        if oxidations:
            sources = [f"{c['element']}" for c in oxidations]
            parts.append(f"Electron source(s): {', '.join(sources)} — these species are oxidized (lose electrons).")

        if reductions:
            sinks = [f"{c['element']} ({c['from_oxidation_state']} → {c['to_oxidation_state']})" for c in reductions]
            parts.append(f"Electron sink(s): {', '.join(sinks)} — these species are reduced (gain electrons).")

        return " ".join(parts)

    def _propose_half_reactions(self, primary_sink, ox_changes):
        """Propose balanced half-reactions."""
        half_rxns = []

        # Reduction half-reaction (sink)
        name = primary_sink.get("name", "")
        type_ = primary_sink.get("type", "")

        if "permanganate" in name.lower() or "MnO4" in name:
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O  (acidic)",
                "alternative": "MnO4⁻ + 2H2O + 3e⁻ → MnO2 + 4OH⁻  (basic)",
                "electrons": 5,
            })
        elif "dichromate" in name.lower() or "Cr2O7" in name or "Cr(VI)" in name:
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "Cr2O7²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H₂O",
                "electrons": 6,
            })
        elif "Fe(III)" in name or "Fe+3" in name or "iron(III)" in name.lower():
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "Fe³⁺ + e⁻ → Fe²⁺",
                "electrons": 1,
            })
        elif "Cu(II)" in name or "copper" in name.lower():
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "Cu²⁺ + 2e⁻ → Cu(s)  or  Cu²⁺ + e⁻ → Cu⁺",
                "electrons": "1 or 2",
            })
        elif "Ag" in name or "silver" in name.lower():
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "Ag⁺ + e⁻ → Ag(s)",
                "electrons": 1,
            })
        elif "Ce(IV)" in name or "ceric" in name.lower():
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "Ce⁴⁺ + e⁻ → Ce³⁺",
                "electrons": 1,
            })
        elif "oxygen" in name.lower() or "O2" == name:
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "O2 + 4H⁺ + 4e⁻ → 2H₂O  (acidic)  or  O2 + 2H2O + 4e⁻ → 4OH⁻  (basic)",
                "electrons": 4,
            })
        elif "[O]" in name or "oxidant" in name.lower():
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": "[O] + 2H⁺ + 2e⁻ → H₂O  (generic oxidant)",
                "electrons": 2,
            })
        else:
            half_rxns.append({
                "type": "reduction (electron sink)",
                "equation": f"(Half-reaction for {name} — provide specific conditions for balancing)",
                "electrons": "unknown",
            })

        # Oxidation half-reaction (source)
        if ox_changes:
            for c in ox_changes:
                if c.get("is_oxidation"):
                    half_rxns.append({
                        "type": "oxidation (electron source)",
                        "equation": f"{c['element']} (oxidation state {c['from_oxidation_state']}) → {c['element']} (oxidation state {c['to_oxidation_state']}) + {abs(c['change'])}e⁻",
                        "electrons": abs(c["change"]) if isinstance(c["change"], (int, float)) else "unknown",
                    })

        return half_rxns

    def _generate_notes(self, primary_sink, matched_agents, organic_patterns):
        """Generate additional analytical notes."""
        notes = []
        conf = primary_sink.get("confidence", "")

        if conf == "low":
            notes.append("⚠️ Low confidence: Consider providing the reaction in SMILES format or a more detailed chemical equation.")
            notes.append("Tip: Include specific oxidant names (e.g., PCC, Swern, KMnO4) for better identification.")

        if len(matched_agents) > 1:
            notes.append(f"ℹ️ Multiple potential oxidants detected ({len(matched_agents)}). The strongest/most specific one is identified as the primary sink.")

        if organic_patterns and not matched_agents:
            notes.append("📌 Organic oxidation detected but no specific oxidant named. Common choices depend on substrate sensitivity.")

        if primary_sink.get("type") == "inferred_from_organic_oxidation":
            notes.append("💡 For alcohol oxidations: PCC (mild, stops at aldehyde), Jones (goes to acid), Dess-Martin (mild, selective).")
            notes.append("💡 For alkene oxidations: mCPBA (epoxidation), OsO4 (syn diol), KMnO4 (cold: diol; hot: cleavage), O3 (ozonolysis).")

        return notes

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        rxn_input = parts[0] if parts else ""
        mode = parts[1] if len(parts) > 1 else "detailed"
        return self._run_base(rxn_input, mode)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
