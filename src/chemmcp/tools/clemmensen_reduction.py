"""
Clemmensen Reduction (Tool #170)
Clemmensen 还原反应：锌汞齐/浓盐酸将羰基还原为亚甲基。
涵盖：机理（锌表面介导的碳正离子或碳负离子途径）、
锌汞齐制备（Zn + HgCl2）、浓盐酸水溶液条件、
底物范围（特别适合芳基酮）、与Wolff-Kishner的对比、限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_CLEMMENSEN_DATA = {
    'mechanism': {
        'overview': 'The mechanism remains debated; two main pathways are proposed — both occur on the surface of zinc metal:',
        'carbene_pathway': (
            '1. Carbonyl oxygen is protonated (conc. HCl) → oxocarbenium ion\n'
            '2. Zn transfers 2e− to carbonyl carbon → radical anion on Zn surface\n'
            '3. Further protonation/electron transfer → zinc carbenoid species (:CR2 on Zn surface)\n'
            '4. Carbenoid protonated → alkane product (CH2)\n'
            'This pathway explains why acid-sensitive groups are destroyed (strongly acidic medium).'
        ),
        'carbanion_pathway': (
            '1. Carbonyl forms zinc-chelated complex on Zn surface\n'
            '2. Electron transfer from Zn → carbanion-like intermediate\n'
            '3. Protonation by HCl → alcohol intermediate (on surface)\n'
            '4. Further reduction of alcohol → alkane\n'
            'This pathway involves a zinc-bound alcohol/carbocation intermediate.'
        ),
        'consensus': 'Both pathways may operate depending on substrate and exact conditions. The key point: reduction occurs on the Zn surface in strongly acidic media.',
    },
    'zinc_amalgam_preparation': {
        'method': 'Treat granular Zn with aqueous HgCl2 (mercuric chloride) solution',
        'procedure': (
            "1. Add granular zinc (30 mesh, 10 g) to a flask\n"
            "2. Prepare solution of HgCl2 (0.5 g) in water (5 mL) + concentrated HCl (1 mL)\n"
            "3. Add HgCl2 solution to Zn, swirl for 5-10 min\n"
            "4. Decant liquid, wash Zn amalgam with water, then acetone, then ether\n"
            "5. Use immediately or store under inert atmosphere (amalgam degrades over time)"
        ),
        'appearance': 'Zn becomes slightly grayish; surface Hg layer visible',
        'role_of_Hg': 'Hg increases the overpotential of H2 evolution on Zn surface, allowing carbonyl reduction to compete with H2 formation. Also improves Zn surface reactivity.',
        'alternatives': ('Activated Zn (without Hg): possible but less effective; '
                         'Zn dust activated by HCl washing: sometimes works; '
                         'Ultra-dilute conditions can use plain Zn dust'),
    },
    'reaction_conditions': {
        'acid': 'Concentrated hydrochloric acid (conc. HCl, ~12 M)',
        'amount': 'Large excess (often as solvent or co-solvent)',
        'solvent': ['Water', 'Ethanol', 'Toluene/Acetic acid mixtures', 'Dioxane/water'],
        'temperature': 'Reflux (for aqueous/EtOH: ~78-110°C)',
        'time': 'Several hours to overnight (sometimes days for stubborn substrates)',
        'zinc_amalgam': 'Freshly prepared Zn(Hg); large excess (10-100 eq relative to substrate)',
        'setup': 'Reflux condenser (to prevent HCl evaporation); addition funnel for adding more Zn/HCl if needed',
        'modifications': {
            'original_Clemmensen': 'Zn(Hg) / conc. HCl(aq), reflux — classic conditions',
            'Yamamoto_Kishi modification': 'Zn(Hg) / TMSCl / HCl / THF — milder, anhydrous; works for acid-sensitive substrates!',
            'ethanolic_HCl': 'HCl gas dissolved in ethanol — milder than aq. conc. HCl',
            'dilute_method': 'Very dilute conditions with long reaction time — gentler',
        },
    },
    'scope': [
        'Aryl alkyl ketones: EXCELLENT — Ar-CO-R → Ar-CH2-R (the CLASSIC Clemmensen application)',
        'Diaryl ketones: GOOD — Ar-CO-Ar → Ar-CH2-Ar',
        'Aliphatic ketones: MODERATE — works but W-K may be better for base-stable aliphatic substrates',
        'Cyclic ketones: GOOD — cyclohexanone → cyclohexane; cyclopentanone → cyclopentane',
        'Steroidal ketones: EXCELLENT — traditional steroid C=O → CH2 method (before W-K became popular)',
        'α,β-Unsaturated ketones: CAN work — may saturate C=C as well (full reduction to saturated alkane)',
        'Aldehydes: WORKABLE but may be tricky — aldehydes can polymerize under acidic conditions',
        'Keto-acids: CAN reduce C=O while preserving COOH (sometimes)',
    ],
    'limitations': [
        ('Strongly ACIDIC conditions', 'Conc. HCl destroys ACID-LABILE groups: acetals, THP ethers, Boc, t-Boc esters, epoxides, glycosidic bonds', 'Solution: Use Wolff-Kishner (basic) for acid-sensitive substrates; or Yamamoto-Kishi modified Clemmensen'),
        ('Also not compatible with many base-sensitive groups', 'Although acidic, some groups that are only stable in neutral pH range also suffer', 'Solution: Evaluate each functional group carefully'),
        ('Does NOT reduce carboxylic acids or esters', 'Only aldehydes and ketones (and some α,β-unsaturated systems)', 'Solution: Need multi-step sequence for full acid → methylene conversion'),
        ('Zinc amalgam preparation required', 'Extra step before reaction; HgCl2 is TOXIC', 'Solution: Prepare fresh each time; use proper Hg waste disposal'),
        ('Mercury environmental concerns', 'Hg is a heavy metal pollutant; must be properly disposed', 'Solution: Collect all Hg-containing waste as hazardous waste; consider alternatives'),
        ('Long reaction time', 'Some substrates require 12-48 hours of reflux', 'Solution: Add fresh portions of Zn/HCl; ensure efficient reflux'),
        ('Functional group interference', 'Some groups get reduced (C=C in enones), others destroyed (OH may dehydrate)', 'Solution: Check compatibility; protect sensitive groups'),
        ('Acid-catalyzed side reactions', 'Dehydration of alcohols, rearrangement of carbocations, pinacol rearrangement possible', 'Solution: Consider substrate-specific side reactions; milder alternatives'),
    ],
    'vs_wolff_kishner': {
        'clemmensen_better_for': [
            'Acid-stable, base-sensitive substrates (opposite of W-K)',
            'Aryl ketones (traditional strength)',
            'When hydrazine availability is limited',
            'Substrates that decompose at high temperature (W-K needs ~200°C)',
        ],
        'wolff_kishner_better_for': [
            'Base-stable, acid-sensitive substrates (opposite of Clemmensen)',
            'When mercury/zinc waste is undesirable',
            'Highly acid-labile molecules',
        ],
        'neither_works_for': [
            'Substances sensitive to BOTH strong acid AND strong base',
            '→ Use thioacetal/desulfurization (Raney Ni), or other alternatives',
        ],
    },
    'typical_yields': {
        'aryl_alkyl_ketone': '75-95%',
        'diaryl_ketone': '70-90%',
        'aliphatic_ketone': '60-85%',
        'cyclic_ketone': '80-93%',
        'steroid_ketone': '75-92%',
        'αβ_unsaturated_ketone': '50-80% (may be messy)',
        'aldehyde': '50-75% (can be low due to polymerization)',
    },
    'workup_procedure': (
        "After TLC/GC shows completion:\n"
        "(1. Cool reaction mixture to RT\n"
        "(2. Filter through Celite® to remove Zn residue (rinse well with solvent!)\n"
        "(3. Separate organic layer (if biphasic); extract aqueous layer with Et2O/EtOAc (×3)\n"
        "(4. Wash combined organics with water, sat. NaHCO3 (CAUTION: CO2!), then brine\n"
        "(5. Dry (Na2SO4 or MgSO4), filter, concentrate\n"
        "(6. Purify by column chromatography or distillation\n"
        "NOTE: Zn/Hg waste must be collected as HEAVY METAL HAZARDOUS WASTE!"
    ),
}


@ChemMCPManager.register_tool
class ClemmensenReduction(BaseTool):
    __version__ = "0.1.0"
    name = "ClemmensenReduction"
    func_name = 'analyze_clemmensen_reduction'
    description = "Clemmensen reduction analysis: zinc amalgam/concentrated HCl-mediated reduction of carbonyl groups (aldehydes/ketones) to methylene. Covers mechanism (two pathways: carbene-like on Zn surface vs carbanion; both emphasize Zn-surface chemistry), zinc amalgam preparation (Zn + HgCl2 procedure with Hg's role), detailed condition parameters (conc. HCl, reflux, modifications including Yamamoto-Kishi), scope (8 categories, excellent for aryl ketones), limitations (8 items including acid sensitivity, Hg concerns), comparison with Wolff-Kishner (selection guide), typical yields, and workup procedure."
    implementation_description = "Comprehensive knowledge base covering: dual mechanistic pathways (carbene vs carbanion on Zn surface), complete Zn(Hg) amalgam preparation protocol, 4 condition variants (original, Yamamoto-Kishi, ethanolic, dilute), 8 scope categories, 8 limitations with solutions, W-K vs Clemmensen decision tree, yield benchmarks, and standardized workup with hazardous waste note."
    categories = ["Reaction"]
    tags = ["Clemmensen", "Reduction", "Carbonyl", "Methylene", "Zinc", "Mercury", "Acidic Conditions", "Zinc Amalgam"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("carbonyl_smiles", "str", "N/A", "SMILES or name of the carbonyl compound (aldehyde or ketone)."),
        ("modification", "str", "classic", "Condition modification: 'classic' (conc. HCl aq.), 'Yamamoto-Kishi', or 'ethanolic'."),
        ("solvent", "str", "aqueous/HCl", "Solvent system."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: carbonyl [modification] [solvent]. E.g., 'acetophenone classic aqueous/HCl'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing carbonyl_analysis, predicted_product, mechanism_details, zinc_amalgam_protocol, optimal_conditions, scope, limitations, vs_wolff_kishner_comparison, typical_yields, safety_notes, workup, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"carbonyl_smiles": "acetophenone", "modification": "classic", "solvent": "aqueous/HCl"},
            "text_input": {"query": "acetophenone classic aqueous/HCl"},
            "output": {"result": {
                "reaction": "Clemmensen: acetophenone (PhCOMe) → ethylbenzene (PhEt)",
                "carbonyl_analysis": {"type": "aromatic ketone (aryl alkyl)", "product_type": "alkylarene"},
                "predicted_product": "ethylbenzene (PhCH2CH3)",
                "transformation": "C=O → CH2 (carbonyl reduced to methylene)",
                "conditions": {"reagent": "Zn(Hg) + conc. HCl (aq)", "T": "reflux", "time": "6-24 h"},
                "yield": "82-94%",
                "key_advantage": "Acidic conditions — base-sensitive groups survive; excellent for aryl ketones",
            }},
        },
        {
            "code_input": {"carbonyl_smiles": "fluorenone", "modification": "classic", "solvent": "aqueous/HCl"},
            "text_input": {"query": "fluorenone classic aqueous/HCl"},
            "output": {"result": {
                "reaction": "Clemmensen: fluorenone → fluorene",
                "carbonyl_analysis": {"type": "diaryl cyclic ketone", "product_type": "polycyclic arene"},
                "predicted_product": "fluorene",
                "yield": "85-95%",
                "note": "Classic textbook example of Clemmensen reduction — fluorenone to fluorene is nearly quantitative",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_CLEMMENSEN_DATA)

    def _run_base(self, carbonyl_smiles: str, modification: str = "classic", solvent: str = "aqueous/HCl") -> dict:
        if not carbonyl_smiles:
            raise ChemMCPInputError("Carbonyl compound is required.")

        carb = self._analyze_carbonyl(carbonyl_smiles)
        prod = self._predict_product(carb)
        cond = self._optimize(carb, modification, solvent)
        limits = self._relevant_limitations(carb)

        result = {
            "result": {
                "reaction": f"Clemmensen: {carbonyl_smiles} → {prod.get('name','?')}",
                "carbonyl_analysis": carb,
                "predicted_product": prod,
                "mechanism_details": self.data['mechanism'],
                "zinc_amalgam_protocol": self.data['zinc_amalgam_preparation'],
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "vs_wolff_kishner": self.data['vs_wolff_kishner'],
                "typical_yields": self._estimate_yield(carb),
                "safety_notes": "⚠️ CONC. HCl = SEVERE CORROSIVE! HgCl2 = HIGHLY TOXIC! Heavy metal waste! Full PPE + fume hood!",
                "workup_procedure": self.data['workup_procedure'],
                "summary": f"Clemmensen: {prod.get('name','?')}. Yield: {self._estimate_yield(carb)}. ⚠️ STRONG ACID + Hg TOXICITY!",
            }
        }
        logger.info(f"Clemmensen: {carbonyl_smiles}")
        return result

    def _analyze_carbonyl(self, smi):
        s = (smi or "").strip().lower()
        carb_types = [
            ('aromatic ketone (aryl alkyl)', ['acetophenone', r'phenyl.*ketone', r'phcor', r'propiophenone'], 'Ar-CO-R → Ar-CH2-R'),
            ('diaryl ketone', ['benzophenone', r'diphenyl_ketone', r'phcop h', r'fluorenone'], 'Ar-CO-Ar → Ar-CH2-Ar'),
            ('aliphatic ketone', ['2-butanone', r'aliphatic.*ketone', r'cyclohexanone'], 'R-CO-R\' → R-CH2-R\''),
            ('aldehyde', ['benzaldehyde', r'butyraldehyde', r'hexanal', r'aldehyde'], 'R-CHO → R-CH3'),
            ('cyclic ketone', ['cyclohexanone', r'cyclopentanone', r'cyclic.*ketone', r'fluorenone'], 'Ring C=O → CH2'),
            ('α,β-unsaturated ketone', ['benzylideneacetone', r'enone', r'αβ.*unsat'], 'May fully reduce C=C and C=O'),
            ('steroidal ketone', ['cholestenone', r'steroid.*ketone'], 'Steroid C=O → CH2'),
        ]
        for ctype, pats, note in carb_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": ctype, "input": smi, "note": note}
        return {"type": "unknown_carbonyl", "input": smi}

    def _predict_product(self, carb):
        ct = carb.get('type', '')
        inp = carb.get('input', '?')
        if 'aromatic ketone' in ct:
            return {"name": f"{inp} → corresponding alkylarene (Ar-CH2-R)", "transformation": "C=O → CH2",
                    "yield": "75-95%"}
        elif 'diaryl' in ct:
            return {"name": f"{inp} → diarylmethane (Ar-CH2-Ar\')", "transformation": "C=O → CH2",
                    "yield": "70-90%"}
        elif 'aliphatic ketone' in ct:
            return {"name": f"{inp} → corresponding alkane (R-CH2-R\')", "transformation": "C=O → CH2",
                    "yield": "60-85%"}
        elif 'aldehyde' in ct:
            return {"name": f"{inp} → corresponding alkane (R-CH3)", "transformation": "C=O → CH2",
                    "yield": "50-75%"}
        elif 'cyclic' in ct.lower():
            return {"name": f"{inp} → corresponding cycloalkane/arene", "transformation": "ring C=O → CH2",
                    "yield": "80-95%" if 'fluorenone' in inp.lower() else "80-93%"}
        elif 'unsaturated' in ct.lower():
            return {"name": f"{inp} → fully/partially reduced product", "transformation": "C=O → CH2; C=C may also be reduced",
                    "yield": "50-80%"}
        return {"name": f"reduced product from {inp}", "transformation": "C=O → CH2", "yield": "65-88%"}

    def _optimize(self, carb, modification, solvent):
        mod = (modification or "classic").strip().lower()
        cond = {
            "reagent": "Zinc amalgam Zn(Hg) (freshly prepared, large excess) + conc. HCl (aq)",
            "solvent": solvent or "water / ethanol / toluene-HOAc mixtures",
            "temperature": "Reflux (78-110°C depending on solvent)",
            "modification": mod,
            "zinc_amount": "10-100 eq (large excess; reaction consumes Zn)",
            "acid_amount": "Conc. HCl as co-solvent (large excess, typically 20-50 mL per mmol)",
            "time": "6-24 hours (some substrates need longer)",
            "setup": "Reflux condenser; add Zn/HCl in portions if needed",
            "procedure_classic": (
                "1. Prepare Zn(Hg) amalgam from granular Zn + HgCl2\n"
                "2. Dissolve/suspend carbonyl compound (1 eq) in solvent\n"
                "3. Add conc. HCl (20-50 mL per mmol)\n"
                "4. Add Zn(Hg) (large excess)\n"
                "5. Reflux with vigorous stirring until TLC shows completion\n"
                "6. Filter hot mixture through Celite®, extract product"
            ) if mod == 'classic' else "Follow specific modification protocol",
        }
        if 'yamamoto' in mod or 'kishi' in mod:
            cond["note"] = "Yamamoto-Kishi: Zn(Hg)/TMSCl/HCl/THF — milder, anhydrous; better for acid-sensitive substrates"
        return cond

    def _relevant_limitations(self, carb):
        relevant = []
        relevant.append({"issue": "Strongly acidic conditions", "problem": "Conc. HCl destroys acid-labile protecting groups (acetal, THP, Boc, etc.)",
                        "Solution": "Use Wolff-Kishner for acid-sensitive substrates; try Yamamoto-Kishi modification"})
        relevant.append({"issue": "Mercury toxicity/environment", "problem": "HgCl2 used in amalgam prep; Hg waste generated",
                        "Solution": "Proper hazardous waste disposal; consider W-K or other alternatives"})
        ct = carb.get('type', '')
        if 'aldehyde' in ct.lower():
            relevant.append({"issue": "Aldehyde complications", "problem": "Aldehydes can polymerize under strongly acidic conditions",
                           "Solution": "Lower T; shorter time; accept lower yield or use alternative"})
        if 'unsaturated' in ct.lower():
            relevant.append({"issue": "C=C reduction", "problem": "α,β-Unsaturated ketones may have C=C reduced too",
                           "Solution": "May be desirable (full saturation) or problematic; test small scale"})
        return relevant

    def _estimate_yield(self, carb):
        ct = carb.get('type', '')
        if 'fluorenone' in (carb.get('input') or '').lower(): return "85-95%"
        if 'aromatic ketone' in ct: return "75-95%"
        if 'cyclic' in ct.lower(): return "80-93%"
        if 'diaryl' in ct: return "70-90%"
        if 'aliphatic' in ct: return "60-85%"
        if 'aldehyde' in ct: return "50-75%"
        if 'unsaturated' in ct: return "50-80%"
        return "65-88%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        carb = parts[0] if parts else ""
        mod = parts[1] if len(parts) > 1 else "classic"
        solv = parts[2] if len(parts) > 2 else "aqueous/HCl"
        return self._run_base(carb, mod, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
