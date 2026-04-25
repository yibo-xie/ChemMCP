"""
Wolff-Kishner Reduction (Tool #169)
Wolff-Kishner 还原反应：通过腙中间体在强碱性高温条件下将羰基还原为亚甲基。
涵盖：机理（腙形成→碱催化去质子化→放N2→碳负离子→质子化）、
Huang-Minlon改良法（高沸点溶剂/回流）、条件（强碱、~200°C）、
底物范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_WOLFF_KISHNER_DATA = {
    'mechanism': [
        ('1', 'Hydrazone formation', 'Carbonyl compound (aldehyde or ketone) + hydrazine (NH2NH2) → hydrazone + H2O. Acid-catalyzed; often requires removal of water to drive equilibrium.'),
        ('2', 'Deprotonation', 'Strong base (KOH or NaOH) deprotonates the hydrazone → resonance-stabilized anion. This is usually RATE-LIMITING.'),
        ('3', 'Loss of N2', 'Anion collapses with expulsion of N2 gas → carbanion intermediate. This step is IRREVERSIBLE (driven by N2 evolution).'),
        ('4', 'Protonation', 'Carbanion is protonated by solvent → alkane product (methylene = CH2). The carbonyl carbon has been reduced from C=O to CH2.'),
    ],
    'huang_minlon_modification': {
        'description': 'The practical variant used in most labs: high-boiling solvent (HO(CH2CH2O)2H = diethylene glycol, bp 245°C) + KOH + heat (~200°C)',
        'advantages': ['One-pot procedure (hydrazone formation and decomposition in same pot)', 'High temperature achievable without autoclave', 'Water removed azeotropically during heating'],
        'solvent': 'Diethylene glycol (DEG) or triethylene glycol (bp 285°C) — high boiling point essential',
        'base': 'KOH pellets or flakes (large excess, 5-10 eq)',
        'temperature': '180-210°C',
        'hydrazine': 'Hydrazine hydrate (NH2NH2·H2O, 2-5 eq) or anhydrous NH2NH2',
    },
    'conditions_comparison': {
        'classical_Wolf_Kishner': {'T': '~200°C in sealed tube', 'solvent': 'No solvent (neat) or high-boiling ether', 'base': 'KOH or NaONHNH2', 'note': 'Original method; requires sealed tube'},
        'Huang_Minlon': {'T': '180-210°C (reflux)', 'solvent': '(HOCH2CH2)2O (diethylene glycol)', 'base': 'KOH (5-10 eq)', 'note': 'Most common modern variant'},
        'Cram_modification': {'T': '20-25°C (RT!)', 'reagent': 'NaCN + NH2NH2 in DMSO', 'note': 'Room temperature W-K! For TMS-hydrazones'},
        'Bamford_Stevens': {'variant': 'Related: tosylhydrazone + base → alkene (not alkane)', 'note': 'Gives alkenes via carbene pathway'},
    },
    'scope': [
        'Aliphatic aldehydes → alkanes: GOOD — RCHO → RCH3',
        'Aliphatic ketones → alkanes: EXCELLENT — RCOR\' → RCH2R\' (the flagship transformation)',
        'Aromatic ketones → alkylarenes: EXCELLENGTH — ArCOR → ArCH2R (very useful for arene side-chain shortening)',
        'Diketones: CAN be reduced — one or both C=O groups depending on conditions',
        'α,β-Unsaturated ketones: VARIABLE — may reduce C=O only, or also reduce C=C (conjugate reduction possible)',
        'Steroidal ketones: CLASSIC application — steroid C=O → CH2 without affecting other functionality',
        'Cyclic ketones: EXCELLENT — cyclohexanone → cyclohexane; cyclopentanone → cyclopentane',
        'Acid-sensitive substrates: EXCELLENT — BASIC conditions (opposite of Clemmensen which is acidic)',
    ],
    'limitations': [
        ('Harsh basic conditions', 'Strong base (KOH/NaOH) at 180-210°C — very harsh!', 'Solution: Substrate must be base-stable; use Clemmensen for base-sensitive compounds'),
        ('High temperature required', '180-210°C needed for reasonable reaction rates', 'Solution: Use high-boiling DEG solvent; sealed tube for lower-boiling alternatives'),
        ('Hydrazine toxicity', 'Hydrazine is HIGHLY TOXIC (carcinogenic), potentially explosive', 'Solution: Handle in fume hood with PPE; use hydrazine hydrate (safer than anhydrous); minimize exposure'),
        ('Not compatible with base-sensitive groups', 'Esters, amides, β-lactams may hydrolyze under these conditions', 'Solution: Protect or remove base-labile groups before W-K'),
        ('α,β-Unsaturated systems', 'May give mixture of C=O reduced and/or conjugate-reduced products', 'Solution: Optimize conditions; consider alternative (LiAlH4 then dehydration/hydrogenation)'),
        ('Long reaction time', 'Typically 2-8 hours at reflux temperature', 'Solution: Plan accordingly; microwave-assisted variants can shorten time'),
        ('Water removal important', 'Hydrazone formation is equilibrium-driven; water must be removed', 'Solution: Huang-Minlon: water distills off azeotropically with some DEG during heating'),
        ('Does NOT reduce carboxylic acids / esters', 'Only aldehydes and ketones are reduced', 'Solution: Convert acid to ester → reduce to alcohol → tosylate → LiAlH4 etc. for full reduction chain'),
    ],
    'vs_clemmensen': {
        'W_K_advantages': ['Basic conditions (good for acid-sensitive substrates)', 'No zinc/mercury waste', 'Generally cleaner for aromatic ketones', 'Better functional group tolerance for base-stable groups'],
        'Clemmensen_advantages': ['Faster (sometimes)', 'Works on some acid-stable substrates that are base-sensitive', 'No toxic hydrazine needed'],
        'selection_rule': 'Base-stable substrate → Wolff-Kishner; Acid-stable but base-sensitive → Clemmensen; Sensitive to BOTH → find alternative (thioacetal/desulfurization, etc.)',
    },
    'typical_yields': {
        'aliphatic_ketone_to_alkane': '70-92%',
        'aromatic_ketone_to_alkylarene': '75-95%',
        'aldehyde_to_alkane': '60-85%',
        'cyclic_ketone': '80-95%',
        'steroid_ketone': '70-90%',
        'αβ_unsaturated_ketone': '50-80% (may be messy)',
    },
    'workup_procedure': (
        "After reflux period (TLC shows completion):\n"
        "(1. Cool reaction mixture to RT\n"
        "(2. Dilute carefully with water (CAUTION: hot basic solution!)\n"
        "(3. Extract with Et2O or EtOAc (×3)\n"
        "(4. Wash organic layer with water, then brine\n"
        "(5. Dry (Na2SO4), filter, concentrate\n"
        "(6. Purify by column chromatography or distillation\n"
        "NOTE: N2 gas evolution indicates active reduction occurring!"
    ),
}


@ChemMCPManager.register_tool
class WolffKishnerReduction(BaseTool):
    __version__ = "0.1.0"
    name = "WolffKishnerReduction"
    func_name = 'analyze_wolff_kishner_reduction'
    description = "Wolff-Kishner reduction analysis: hydrazone-mediated reduction of carbonyl groups (aldehydes/ketones) to methylene under strongly basic conditions at high temperature. Covers mechanism (4 steps: hydrazone formation, base-catalyzed deprotonation, N2 loss, protonation), Huang-Minlon modification (DEG/KOH/reflux protocol), condition comparison table (classical vs Huang-Minlon vs Cram vs Bamford-Stevens), scope (8 categories including steroids and cyclic ketones), limitations (8 items including hydrazine toxicity), comparison with Clemmensen reduction, typical yields, and workup procedure."
    implementation_description = "Comprehensive knowledge base covering: 4-step mechanism emphasizing N2-driven irreversibility, detailed Huang-Minlon practical protocol, 4 condition variants compared, 8 scope categories, 8 limitations with solutions, W-K vs Clemmensen selection guide, yield benchmarks, and standardized workup."
    categories = ["Reaction"]
    tags = ["Wolff-Kishner", "Reduction", "Carbonyl", "Methylene", "Hydrazone", "Basic Conditions", "Huang-Minlon"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("carbonyl_smiles", "str", "N/A", "SMILES or name of the carbonyl compound (aldehyde or ketone)."),
        ("variant", "str", "Huang-Minlon", "Variant: 'Huang-Minlon' (standard), 'Classical', or 'Cram'."),
        ("base", "str", "KOH", "Base: 'KOH' or 'NaOH'."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: carbonyl [variant] [base]. E.g., 'acetophenone Huang-Minlon KOH'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing carbonyl_analysis, predicted_product, mechanism_steps, huang_minlon_details, optimal_conditions, scope, limitations, vs_clemmensen_comparison, typical_yields, safety_notes, workup, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"carbonyl_smiles": "acetophenone", "variant": "Huang-Minlon", "base": "KOH"},
            "text_input": {"query": "acetophenone Huang-Minlon KOH"},
            "output": {"result": {
                "reaction": "Wolff-Kishner: acetophenone (PhCOMe) → ethylbenzene (PhEt)",
                "carbonyl_analysis": {"type": "aromatic ketone (aryl alkyl)", "product_type": "alkylarene"},
                "predicted_product": "ethylbenzene (PhCH2CH3)",
                "transformation": "C=O → CH2 (carbonyl reduced to methylene)",
                "conditions": {"variant": "Huang-Minlon", "solvent": "diethylene glycol", "base": "KOH (5-10 eq)", "T": "190-210°C"},
                "yield": "80-93%",
                "key_advantage": "Basic conditions — acid-sensitive groups elsewhere in molecule survive",
            }},
        },
        {
            "code_input": {"carbonyl_smiles": "cyclohexanone", "variant": "Huang-Minlon", "base": "KOH"},
            "text_input": {"query": "cyclohexanone Huang-Minlon KOH"},
            "output": {"result": {
                "reaction": "Wolff-Kishner: cyclohexanone → cyclohexane",
                "carbonyl_analysis": {"type": "cyclic aliphatic ketone", "product_type": "cycloalkane"},
                "predicted_product": "cyclohexane",
                "yield": "88-96%",
                "note": "Cyclic ketones are excellent W-K substrates — clean conversion to cycloalkanes",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_WOLFF_KISHNER_DATA)

    def _run_base(self, carbonyl_smiles: str, variant: str = "Huang-Minlon", base: str = "KOH") -> dict:
        if not carbonyl_smiles:
            raise ChemMCPInputError("Carbonyl compound is required.")

        carb = self._analyze_carbonyl(carbonyl_smiles)
        prod = self._predict_product(carb)
        cond = self._optimize(carb, variant, base)
        limits = self._relevant_limitations(carb)

        result = {
            "result": {
                "reaction": f"Wolff-Kishner: {carbonyl_smiles} → {prod.get('name','?')}",
                "carbonyl_analysis": carb,
                "predicted_product": prod,
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['mechanism']],
                "huang_minlon_details": self.data['huang_minlon_modification'],
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "vs_clemmensen": self.data['vs_clemmensen'],
                "typical_yields": self._estimate_yield(carb),
                "safety_notes": "⚠️ HYDRAZINE IS HIGHLY TOXIC AND CARCINOGENIC! Hot concentrated KOH causes SEVERE BURNS! Full PPE required!",
                "workup_procedure": self.data['workup_procedure'],
                "summary": f"W-K reduction: {prod.get('name','?')}. Yield: {self._estimate_yield(carb)}. ⚠️ HIGH T + STRONG BASE + TOXIC N2H4!",
            }
        }
        logger.info(f"Wolff-Kishner: {carbonyl_smiles}")
        return result

    def _analyze_carbonyl(self, smi):
        s = (smi or "").strip().lower()
        carb_types = [
            ('aromatic ketone (aryl alkyl)', ['acetophenone', r'benzophenone', r'phcor', r'phenyl.*ketone', r'propiophenone'], 'Ar-CO-R → Ar-CH2-R'),
            ('aromatic ketone (diaryl)', ['benzophenone', r'diphenyl_ketone', r'phcop h'], 'Ar-CO-Ar → Ar-CH2-Ar'),
            ('aliphatic ketone', ['acetone', r'2-butanone', r'cyclohexanone', r'cyclopentanone', r'aliphatic.*ketone'], 'R-CO-R\' → R-CH2-R\''),
            ('aliphatic aldehyde', ['butyraldehyde', r'hexanal', r'valeraldehyde', r'aldehyde', r'butanal'], 'R-CHO → R-CH3'),
            ('aromatic aldehyde', ['benzaldehyde', r'phcho', r'benzaldehyde'], 'Ar-CHO → Ar-CH3'),
            ('α,β-unsaturated ketone', ['benzylideneacetone', r'chalc one', r'enone', r'αβ.*unsaturated'], 'May give mixed products'),
            ('cyclic ketone', ['cyclohexanone', r'cyclopentanone', r'cyclobutanone', r'cyclic.*ketone'], 'Cyclic → ring with CH2 instead of C=O'),
            ('steroidal ketone', ['cholestenone', r'steroid.*ketone', r'androst'], 'Steroid C=O → CH2'),
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
        elif 'aliphatic ketone' in ct:
            return {"name": f"{inp} → corresponding alkane (R-CH2-R\')", "transformation": "C=O → CH2",
                    "yield": "70-92%"}
        elif 'aldehyde' in ct:
            return {"name": f"{inp} → corresponding alkane (one fewer C atom? No: RCHO → RCH3)", "transformation": "C=O → CH2",
                    "yield": "60-85%"}
        elif 'cyclic' in ct.lower():
            return {"name": f"{inp} → corresponding cycloalkane", "transformation": "ring C=O → CH2",
                    "yield": "80-96%"}
        elif 'unsaturated' in ct.lower():
            return {"name": f"{inp} → partially reduced product (may be mixture)", "transformation": "C=O → CH2; C=C may or may not be affected",
                    "yield": "50-80%"}
        return {"name": f"reduced product from {inp}", "transformation": "C=O → CH2", "yield": "65-88%"}

    def _optimize(self, carb, variant, base):
        var = (variant or "Huang-Minlon").strip()
        cond = {
            "variant": var,
            "hydrazine_source": "Hydrazine hydrate (NH2NH2·H2O, 2-5 eq) or anhydrous NH2NH2",
            "base": f"{base} (5-10 eq; pellets or flakes)",
            "solvent": "Diethylene glycol ((HOCH2CH2)2O, bp 245°C)" if 'minlon' in var.lower() else "High-boiling solvent or neat",
            "temperature": "190-210°C (reflux)" if 'minlon' in var.lower() else "~200°C (sealed tube)",
            "procedure_huang_minlon": (
                "1. Dissolve carbonyl compound (1 eq) in diethylene glycol (10-20 mL/g)\n"
                "2. Add hydrazine hydrate (2-5 eq)\n"
                "3. Add KOH pellets (5-10 eq)\n"
                "4. Fit condenser; heat to 190-210°C (reflux)\n"
                "5. Reflux 2-8 hours (water + some DEG distills off early)\n"
                "6. Cool, dilute with water, extract with Et2O/EtOAc\n"
                "7. Purify product"
            ) if 'minlon' in var.lower() else "Follow classical sealed-tube protocol",
            "time": "2-8 hours (monitor by TLC of cooled aliquot)",
            "indicator": "N2 gas bubbling = active reduction occurring",
        }
        return cond

    def _relevant_limitations(self, carb):
        relevant = []
        relevant.append({"issue": "Harsh basic conditions", "problem": "KOH/NaOH at 190-210°C destroys many functional groups",
                        "Solution": "Ensure substrate is base-stable; protect sensitive groups"})
        relevant.append({"issue": "Hydrazine toxicity", "problem": "N2H4 is carcinogenic and toxic",
                        "Solution": "Use in fume hood with PPE; hydrazine hydrate safer than anhydrous"})
        ct = carb.get('type', '')
        if 'unsaturated' in ct.lower():
            relevant.append({"issue": "Conjugated system complications", "problem": "α,β-unsaturated ketones may give mixtures",
                           "Solution": "Consider alternative reduction strategies"})
        if 'ester' in (carb.get('input') or '').lower() or 'amide' in (carb.get('input') or '').lower():
            relevant.append({"issue": "Other carbonyls may be affected", "problem": "Esters/amides can hydrolyze under harsh basic conditions",
                           "Solution": "Protect or avoid having these groups present"})
        return relevant

    def _estimate_yield(self, carb):
        ct = carb.get('type', '')
        if 'cyclic' in ct.lower(): return "80-96%"
        if 'aromatic ketone' in ct: return "75-95%"
        if 'aliphatic ketone' in ct: return "70-92%"
        if 'aldehyde' in ct: return "60-85%"
        if 'unsaturated' in ct: return "50-80%"
        return "65-88%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        carb = parts[0] if parts else ""
        var = parts[1] if len(parts) > 1 else "Huang-Minlon"
        b = parts[2] if len(parts) > 2 else "KOH"
        return self._run_base(carb, var, b)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
