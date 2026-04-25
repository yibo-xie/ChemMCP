"""
Jones Oxidation (Tool #165)
Jones 氧化反应：CrO3/H2SO4/丙酮体系将醇氧化为羰基化合物（伯醇→酸，仲醇→酮）。
涵盖：机理（铬酸酯形成→消除）、Jones试剂制备、
伯醇过氧化为羧酸、适用范围、铬(VI)毒性/环境问题、限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_JONES_DATA = {
    'mechanism': [
        ('1', 'Chromate ester formation', 'Alcohol attacks Cr(VI) in chromic acid (H2CrO4) → chromate ester (R-O-CrO2H). This is a reversible, acid-catalyzed step.'),
        ('2', 'Elimination (rate-determining)', 'Base (water or acetone) removes α-H from chromate ester → carbonyl compound + Cr(IV). E2-like concerted elimination.'),
        ('3', 'Chromium fate', 'Cr(IV) disproportionates: 2 Cr(IV) → Cr(III) + Cr(VI), or further oxidizes substrate. Final product: green Cr(III) solution.'),
    ],
    'jones_reagent': {
        'preparation': 'Dissolve CrO3 (26.72 g, 0.267 mol) in water (23 mL), then carefully add concentrated H2SO4 (23 mL) dropwise with cooling. Dilute to 100 mL with water.',
        'concentration': '~2.67 M in CrO3 (~2.5 M H2SO4)',
        'appearance': 'Clear orange-red solution',
        'stability': 'Stable for months at RT; can be titrated for exact concentration',
        'storage': 'Store in glass bottle (NOT plastic — H2SO4 degrades many plastics)',
        'safety': 'EXTREMELY CORROSIVE and TOXIC; Cr(VI) is carcinogenic; use full PPE in fume hood',
    },
    'reaction_outcomes': {
        'primary_alcohol': {
            'product': 'CARBOXYLIC ACID (via aldehyde hydrate → further oxidation)',
            'pathway': 'RCH2OH → RCHO → RCH(OH)2 (aldehyde hydrate, favored in aqueous acid) → RCOOH',
            'note': 'Aldehyde intermediate is NOT isolated under standard Jones conditions — it gets over-oxidized to the acid',
            'to_stop_at_aldehyde': 'Use PCC, Swern, or Dess-Martin instead',
        },
        'secondary_alcohol': {
            'product': 'KETONE (stops here — no further oxidation possible)',
            'note': 'Clean conversion; one of the best methods for alcohol → ketone',
        },
        'tertiary_alcohol': {'product': 'NO REACTION (no α-H to eliminate)', 'note': 'Tertiary alcohols are not oxidized by Jones reagent'},
    },
    'scope': [
        'Secondary aliphatic alcohols → ketones: EXCELLENT — clean, high-yielding',
        'Primary aliphatic alcohols → carboxylic acids: EXCELLENT — reliable over-oxidation',
        'Benzylic alcohols: GOOD → benzoic acids (primary) or ketones (secondary)',
        'Allylic alcohols: MAY work but double bond can be affected by acidic conditions',
        'Cyclic alcohols: EXCELLENT → cyclic ketones (cyclohexanone, etc.)',
        'Steroid alcohols: CLASSIC application → steroid ketones',
        'Diols: Selective oxidation of one OH possible if sterically differentiated',
    ],
    'limitations': [
        ('Acid-sensitive groups destroyed', 'Acetals, THP ethers, Boc groups, t-Boc esters, epoxides hydrolyze in strong acid', 'Solution: Use neutral oxidation (Swern, Dess-Martin, TPAP) for acid-labile substrates'),
        ('Over-oxidation of primary alcohols', 'Primary alcohols always give carboxylic acids, NOT aldehydes', 'Solution: Need aldehyde? Use Swern, PCC, Dess-Martin, or TEMPO/NaOCl instead'),
        ('Cr(VI) toxicity and environmental concern', 'Cr(VI) compounds are carcinogenic, mutagenic, highly toxic to aquatic life', 'Solution: Minimize scale; proper waste disposal as hazardous waste; consider greener alternatives'),
        ('Oxidation of other functional groups', 'Thiols → disulfides; sulfides → sulfoxides/sulfones; some alkenes may react', 'Solution: Check compatibility; protect sensitive groups'),
        ('Not compatible with base-labile groups either', 'Strongly ACIDIC medium (H2SO4)', 'Solution: Protect base-sensitive groups before oxidation; use alternative method'),
        ('Stereochemistry at α-carbon', 'If chiral center at α-carbon, racemization can occur via enol/enolate formation', 'Solution: Jones is generally stereospecific (no enolization) but very acidic conditions may cause epimerization'),
        ('Solubility issues', 'Non-polar substrates may have limited solubility in aqueous acetone', 'Solution: Add co-solvent (dioxane, THF); use Collins oxidation (CrO3·pyridine in CH2Cl2) for non-polar substrates'),
        ('Exotherm on addition', 'Adding Jones reagent to alcohol can be exothermic', 'Solution: Add slowly with cooling; control addition rate'),
    ],
    'alternatives_to_jones': {
        'PCC/PDC (pyridinium chlorochromate)': 'Neutral conditions; primary → aldehyde possible; still uses Cr(VI)',
        'Collins oxidation (CrO3·pyridine)': 'Anhydrous CH2Cl2; milder than Jones; still Cr-based',
        'Dess-Martin periodinane': 'Neutral, mild, high-yielding; expensive',
        'Swern oxidation': 'Low-T, no metal; primary → aldehyde; foul-smelling byproduct',
        "TEMPO/NaOCl (Anelli oxidation)": "Aqueous, mild, environmentally friendlier",
        'IBX/DMP': 'High-yielding; selective',
        'TPAP/NMO': 'Catalytic Ru; mild; excellent for sensitive molecules',
        'Oppenauer oxidation': 'For secondary alcohols only; Al(OiPr)3 as catalyst',
    },
    'typical_yields': {
        'secondary_alcohol_to_ketone': '80-98%',
        'primary_alcohol_to_acid': '75-95%',
        'benzylic_alcohol_to_acid': '85-97%',
        'cyclic_alcohol_to_ketone': '85-96%',
        'steroid_oxidation': '80-95%',
    },
    'workup_procedure': (
        'After TLC shows complete consumption of starting material:\n'
        '(1) Carefully pour reaction mixture into cold saturated NaHCO3 solution (CAUTION: CO2 evolution!)\n'
        '(2) Extract with Et2O or EtOAc (×3)\n'
        '(3) Wash organic layer with water, then brine\n'
        '(4) Dry (Na2SO4 or MgSO4), filter, concentrate\n'
        '(5) Purify by column chromatography or recrystallization\n'
        'NOTE: Aqueous layer turns GREEN (Cr(III)) — this indicates successful oxidation'
    ),
}


@ChemMCPManager.register_tool
class JonesOxidation(BaseTool):
    __version__ = "0.1.0"
    name = "JonesOxidation"
    func_name = 'analyze_jones_oxidation'
    description = "Jones oxidation analysis: CrO3/H2SO4/acetone-mediated oxidation of alcohols. Covers mechanism (chromate ester formation → elimination, Cr(VI)→Cr(III)), Jones reagent preparation (composition, safety), reaction outcomes (primary→carboxylic acid via over-oxidation, secondary→ketone, tertiary→no reaction), scope (7 categories), limitations (8 items including Cr(VI) toxicity, acid sensitivity, over-oxidation), comparison with 8 alternative oxidations, typical yields, and workup procedure."
    implementation_description = "Comprehensive knowledge base covering: 3-step mechanism with chromate ester intermediate, detailed Jones reagent recipe (with safety warnings), 3 reaction outcome categories, 7 scope categories, 8 limitations with practical solutions, 8 alternative oxidation methods compared, yield benchmarks, and standardized workup procedure."
    categories = ["Reaction"]
    tags = ["Jones", "Oxidation", "Chromium", "Alcohol", "Ketone", "Carboxylic Acid", "Cr(VI)", "Organic Synthesis"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("alcohol_smiles", "str", "N/A", "SMILES or name of the alcohol substrate (primary, secondary, or tertiary)."),
        ("solvent", "str", "acetone", "Solvent (standard: acetone)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: alcohol [solvent]. E.g., 'cyclohexanol acetone' or '1-hexanol acetone'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing alcohol_analysis, mechanism_steps, jones_reagent_details, product_prediction, outcome_type, scope, limitations, alternatives_comparison, typical_yields, safety_notes, workup, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"alcohol_smiles": "cyclohexanol", "solvent": "acetone"},
            "text_input": {"query": "cyclohexanol acetone"},
            "output": {"result": {
                "reaction": "Jones oxidation: cyclohexanol → cyclohexanone",
                "alcohol_analysis": {"type": "secondary cyclic alcohol", "outcome": "ketone"},
                "product": "cyclohexanone",
                "reagents": "Jones reagent (CrO3/H2SO4) in acetone",
                "conditions": {"T": "0°C → RT", "time": "30 min - 4 h"},
                "yield": "90-97%",
                "monitoring": "Orange color disappears; green color appears (Cr(III) formed)",
            }},
        },
        {
            "code_input": {"alcohol_smiles": "1-hexanol", "solvent": "acetone"},
            "text_input": {"query": "1-hexanol acetone"},
            "output": {"result": {
                "reaction": "Jones oxidation: 1-hexanol → hexanoic acid",
                "alcohol_analysis": {"type": "primary aliphatic alcohol", "outcome": "carboxylic acid (OVER-OXIDIZED)"},
                "product": "hexanoic acid (NOT hexanal!)",
                "pathway": "1-hexanol → hexanal (not isolated) → hexanal hydrate → hexanoic acid",
                "yield": "82-95%",
                "warning": "Primary alcohols ALWAYS give carboxylic acids with Jones — use Swern/PCC/Dess-Martin for aldehydes",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_JONES_DATA)

    def _run_base(self, alcohol_smiles: str, solvent: str = "acetone") -> dict:
        if not alcohol_smiles:
            raise ChemMCPInputError("Alcohol substrate is required.")

        alc = self._analyze_alcohol(alcohol_smiles)
        prod = self._predict_product(alc)
        cond = self._optimize(alc, solvent)
        limits = self._relevant_limitations(alc)

        result = {
            "result": {
                "reaction": f"Jones oxidation: {alcohol_smiles} → {prod.get('name','?')}",
                "alcohol_analysis": alc,
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['mechanism']],
                "jones_reagent": self.data['jones_reagent'],
                "product_prediction": prod,
                "outcome_type": prod.get('outcome_type', ''),
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "alternatives": self.data['alternatives_to_jones'],
                "typical_yields": self._estimate_yield(alc),
                "safety_notes": "CR(VI) IS CARCINOGENIC — wear gloves, goggles, lab coat; work in fume hood; dispose as hazardous waste",
                "workup_procedure": self.data['workup_procedure'],
                "summary": f"Jones oxidation: {alcohol_smiles} → {prod.get('name','?')} ({prod.get('outcome_type','')}). Yield: {self._estimate_yield(alc)}. WARNING: Cr(VI) toxic!",
            }
        }
        logger.info(f"Jones oxidation: {alcohol_smiles}")
        return result

    def _analyze_alcohol(self, smi):
        s = (smi or "").strip().lower()
        alc_types = [
            ('primary aliphatic alcohol', ['1-hexanol', '1-butanol', 'ethanol', r'primary.*alcohol'], '→ carboxylic acid (over-oxidized)'),
            ('primary benzylic alcohol', ['benzyl alcohol', r'phch2oh'], '→ benzoic acid'),
            ('secondary aliphatic alcohol', ['2-butanol', 'isopropanol', r'secondary.*alcohol'], '→ ketone'),
            ('secondary cyclic alcohol', ['cyclohexanol', r'cyclic.*alcohol', r'cyclopentanol'], '→ cyclic ketone'),
            ('secondary allylic alcohol', ['cholesterol', r'crotyl.*alcohol', r'3.*ol.*sterol'], '→ enone (may have acid-catalyzed side reactions)'),
            ('tertiary alcohol', ['t-butanol', r'tertiary.*alcohol', r'tert-butanol'], '→ NO REACTION'),
        ]
        for atype, pats, outcome in alc_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "outcome": outcome}
        return {"type": "unknown_alcohol", "input": smi, "outcome": "unknown"}

    def _predict_product(self, alc):
        at = alc.get('type', '')
        outcome = alc.get('outcome', '')
        if 'primary' in at:
            return {"name": f"{alc.get('input','?')} → corresponding carboxylic acid", "class": "carboxylic acid",
                    "outcome_type": "over-oxidation (primary alcohol → carboxylic acid)", "yield": "80-95%"}
        elif 'secondary' in at:
            return {"name": f"{alc.get('input','?')} → corresponding ketone", "class": "ketone",
                    "outcome_type": "clean oxidation (secondary alcohol → ketone)", "yield": "85-98%"}
        elif 'tertiary' in at:
            return {"name": "NO REACTION — tertiary alcohols lack α-hydrogen", "class": "no reaction",
                    "outcome_type": "no oxidation", "yield": "N/A"}
        return {"name": f"oxidized product from {alc.get('input','?')}", "class": "unknown", "outcome_type": "unknown", "yield": "70-90%"}

    def _optimize(self, alc, solvent):
        at = alc.get('type', '')
        return {
            "reagent": "Jones reagent (CrO3/H2SO4 in H2O)",
            "solvent": solvent or "acetone (or dioxane/acetone mixtures)",
            "temperature": "0°C (ice bath) during addition, then RT",
            "procedure": (
                "1. Dissolve alcohol (1 eq) in acetone (0.1-0.5 M)\n"
                "2. Cool to 0°C (ice bath)\n"
                "3. Add Jones reagent dropwise until orange color persists\n"
                "4. Stir at 0°C → RT until TLC shows completion (30 min - 4 h)\n"
                "5. Quench by pouring into cold sat. NaHCO3 / add isopropanol to destroy excess oxidant"
            ),
            "stoichiometry": "~1.5 eq CrO3 per OH group (typically slight excess)",
            "monitoring": "Color change: orange (Cr(VI)) → green (Cr(III)); TLC for starting material consumption",
            "time": "30 min - 4 h depending on substrate",
        }

    def _relevant_limitations(self, alc):
        relevant = []
        at = alc.get('type', '')
        relevant.append({"issue": "Cr(VI) toxicity", "problem": "Carcinogenic, mutagenic, environmentally hazardous",
                        "solution": "Minimize scale; proper PPE and disposal; consider greener alternatives"})
        if 'primary' in at:
            relevant.append({"issue": "Over-oxidation", "problem": "Primary alcohols give carboxylic acids, not aldehydes",
                           "Solution": "Need aldehyde? Use Swern, Dess-Martin, PCC, or TEMPO/NaOCl"})
        relevant.append({"issue": "Acidic conditions", "problem": "Strong acid (H2SO4) destroys acid-labile protecting groups",
                        "Solution": "Protect acid-sensitive groups first; use neutral oxidation alternatives"})
        if 'allylic' in at.lower():
            relevant.append({"issue": "Double bond effects", "problem": "Acidic conditions may cause alkene isomerization/cyclization",
                           "Solution": "Consider milder oxidation (Dess-Martin, Swern)"})
        return relevant

    def _estimate_yield(self, alc):
        at = alc.get('type', '')
        if 'secondary' in at and 'cyclic' in at: return "85-97%"
        if 'secondary' in at: return "85-98%"
        if 'primary' in at and 'benzyl' in at: return "85-97%"
        if 'primary' in at: return "80-95%"
        if 'tertiary' in at: return "N/A (no reaction)"
        return "75-92%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        alc = parts[0] if parts else ""
        solv = parts[1] if len(parts) > 1 else "acetone"
        return self._run_base(alc, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
