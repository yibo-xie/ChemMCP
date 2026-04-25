"""
Swern Oxidation (Tool #164)
Swern 氧化反应：DMSO/(COCl)2/Et3N 低温下将伯醇氧化为醛、仲醇氧化为酮。
涵盖：机理（DMSO活化→烷氧基硫鎓盐→叶立德消除）、低温要求、
后处理（无过氧化！）、与其他氧化方法的对比、底物范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_SWERN_DATA = {
    'mechanism': [
        ('1', 'DMSO activation (−78°C → −60°C)', 'Oxalyl chloride ((COCl)2) reacts with DMSO → chlorodimethylsulfonium ion + CO2 + CO. EXOTHERMIC — must add slowly at low T!'),
        ('2', 'Alcohol addition (−78°C)', 'Alcohol attacks activated sulfur → alkoxysulfonium salt intermediate. Color changes (orange/red) indicate this step.'),
        ('3', 'Base addition (−78°C → RT)', 'Et3N deprotonates the α-carbon of alkoxysulfonium salt → sulfur ylide intermediate'),
        ('4', 'Elimination', 'Ylide collapses via intramolecular elimination → carbonyl product (aldehyde/ketone) + dimethyl sulfide (Me2S, foul smell!)'),
    ],
    'temperature_profile': {
        'step_1_temp': '−78°C (dry ice/acetone) to −60°C',
        'step_2_temp': '−78°C (must keep cold during alcohol addition)',
        'step_3_temp': 'Start at −78°C, then allow to warm to RT over 30-60 min',
        'critical_note': 'LOW TEMPERATURE IS ESSENTIAL — higher temperatures lead to side reactions (Pummerer rearrangement, elimination)',
    },
    'reagents': {
        'DMSO': {'role': 'Oxidant (oxygen source)', 'equivalents': '1.5-3 eq (usually 2 eq)', 'note': 'Anhydrous! Water quenches the activation'},
        'oxalyl_chloride': {'role': 'Activates DMSO', 'equivalents': '1.0-1.5 eq (usually 1.2 eq)', 'note': 'HIGHLY toxic; lachrymator; handle in fume hood with PPE; forms CO/CO2 gas'},
        'triethylamine': {'role': 'Base for deprotonation', 'equivalents': '3-5 eq (usually 3-5 eq)', 'note': 'Must be anhydrous; scavenges HCl produced'},
        'solvent': {'choice': 'CH2Cl2 (dry, distilled from CaH2)', 'concentration': '0.1-0.5 M relative to alcohol'},
        'alternatives_to_oxalyl_chloride': {
            '(COBr)2 (oxalyl bromide)': 'Similar reactivity; more expensive',
            'TFAA (trifluoroacetic anhydride)': 'Swern variant — milder activation',
            'SO3·Py (sulfur trioxide-pyridine)': 'Parikh-Doering oxidation — related method',
        },
    },
    'advantages_vs_other_oxidations': {
        'vs_PCC_PDC': 'No chromium waste; cheaper; no heavy metal contamination',
        'vs_Dess_Martin': 'Much cheaper reagents; DMSO is inexpensive',
        'vs_Jones': 'No acidic conditions — acid-sensitive groups survive; NO OVER-OXIDATION of primary alcohols to acids!',
        'vs_TPAP_NMO': 'Cheaper than TPAP (ruthenium catalyst)',
        'vs_IBX_DMP': 'Gentler conditions; works on sensitive substrates',
        'key_advantage': 'Primary alcohols → ALDEHYDES (not carboxylic acids) — unique advantage over many oxidants',
    },
    'scope': [
        'Primary alcohols → aldehydes: EXCELLENT — the flagship transformation (no over-oxidation!)',
        'Secondary alcohols → ketones: EXCELLENT — clean conversion',
        'Allylic alcohols: GOOD — α,β-unsaturated carbonyls formed; double bond may isomerize slightly',
        'Benzylic alcohols: EXCELLENT — aromatic aldehydes cleanly formed',
        'Sugar/polyol substrates: GOOD — selective oxidation possible with protecting groups',
        'Acid-sensitive substrates: EXCELLENT — neutral to slightly basic workup',
        'Base-sensitive substrates: CAUTION — Et3N is used as base',
        'Thiol-sensitive substrates: CAUTION — Me2S byproduct can react with thiols/disulfides',
    ],
    'limitations': [
        ('Low temperature requirement', 'Reaction must be run at −78°C (dry ice/acetone); not all labs equipped', 'Solution: Use dry ice/acetone bath or cryocooler; consider Dess-Martin if cryogenic unavailable'),
        ('Foul-smelling byproduct (Me2S)', 'Dimethyl sulfide has extremely unpleasant odor (rotten cabbage/threshold ~0.02 ppm)', 'Solution: Work in efficient fume hood; trap Me2S with activated charcoal or CuSO4 solution during workup'),
        ('Moisture sensitivity', 'Water destroys oxalyl chloride and activated DMSO species', 'Solution: Rigorously anhydrous conditions; dried glassware; inert atmosphere (N2/Ar)'),
        ('Oxalyl chloride hazards', '(COCl)2 is highly toxic, corrosive lachrymator; releases CO/CO2 gas', 'Solution: Handle in fume hood with full PPE; use syringe addition behind shield'),
        ('Scale-up difficulty', 'Cryogenic conditions make large-scale reactions challenging', 'Solution: Use alternative oxidations (Dess-Martin, TEMPO/NaOCl) for >10 mmol scale'),
        ('Pummerer rearrangement side reaction', 'At elevated temperatures, alkoxysulfonium salt can undergo Pummerer rearrangement → α-oxygenated sulfide', 'Solution: Strict temperature control (≤−60°C during key steps)'),
        ('Elimination side reaction', 'Some substrates can undergo E2 elimination instead of oxidation', 'Solution: Optimize temperature; ensure proper stoichiometry'),
        ('Not suitable for phenols', 'Phenols are not oxidized under Swern conditions', 'Solution: Use other methods (Fremy\'s salt, IBX) for phenol oxidation'),
    ],
    'typical_yields': {
        'primary_alcohol_to_aldehyde': '75-95%',
        'secondary_alcohol_to_ketone': '80-98%',
        'allylic_alcohol': '70-90%',
        'benzylic_alcohol': '85-96%',
        'sugar_substrate': '65-88%',
    },
    'workup_procedure': (
        'After reaction reaches RT and stirs 30-60 min: '
        '(1) Dilute with water or sat. aq. NH4Cl, '
        '(2) Extract with CH2Cl2 or Et2O (×3), '
        '(3) Wash organic layer with water, then brine, '
        '(4) Dry (Na2SO4 or MgSO4), filter, concentrate, '
        '(5) Purify by column chromatography or distillation. '
        'NOTE: Product is NOT further oxidized — primary alcohols give aldehydes!'
    ),
}


@ChemMCPManager.register_tool
class SwernOxidation(BaseTool):
    __version__ = "0.1.0"
    name = "SwernOxidation"
    func_name = 'analyze_swern_oxidation'
    description = "Swern oxidation analysis: DMSO/(COCl)2/Et3N-mediated oxidation of primary alcohols to aldehydes and secondary alcohols to ketones at low temperature. Covers mechanism (4 steps: DMSO activation, alcohol addition, base-promoted elimination via ylide), critical temperature profile (−78°C), reagent details (DMSO, oxalyl chloride, Et3N equivalents), advantages vs other oxidations (no Cr, no over-oxidation, acid-compatible), scope (6 categories), limitations (8 items including odor, scale-up, hazards), typical yields, and workup procedure."
    implementation_description = "Comprehensive knowledge base covering: 4-step mechanism with temperature-critical notes, detailed reagent table (with alternatives to oxalyl chloride), comparison with 6 other oxidation methods (PCC, Dess-Martin, Jones, TPAP, IBX), 8 scope categories, 8 limitations with practical solutions, yield benchmarks, and standardized workup procedure."
    categories = ["Reaction"]
    tags = ["Swern", "Oxidation", "Alcohol", "Aldehyde", "Ketone", "DMSO", "Organic Synthesis"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("alcohol_smiles", "str", "N/A", "SMILES or name of the alcohol substrate (primary or secondary)."),
        ("temperature", "str", "-78°C", "Reaction temperature (standard: -78°C)."),
        ("base", "str", "Et3N", "Base for deprotonation."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: alcohol [temperature] [base]. E.g., '1-hexanol -78°C Et3N'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing alcohol_analysis, mechanism_steps, reagent_details, temperature_profile, product_prediction, scope, limitations, advantages_comparison, typical_yields, workup, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"alcohol_smiles": "1-hexanol", "temperature": "-78°C", "base": "Et3N"},
            "text_input": {"query": "1-hexanol -78°C Et3N"},
            "output": {"result": {
                "reaction": "Swern oxidation: 1-hexanol → hexanal",
                "alcohol_analysis": {"type": "primary aliphatic alcohol", "product_type": "aldehyde"},
                "product": "hexanal (hexanaldehyde)",
                "reagents": {"DMSO": "2 eq", "(COCl)2": "1.2 eq", "Et3N": "5 eq", "solvent": "CH2Cl2 (0.2 M)"},
                "conditions": {"T": "−78°C → RT", "time": "2-4 h total", "atmosphere": "N2"},
                "yield": "80-92%",
                "key_advantage": "NO over-oxidation to hexanoic acid — Swern stops at aldehyde!",
            }},
        },
        {
            "code_input": {"alcohol_smiles": "cholesterol", "temperature": "-78°C", "base": "Et3N"},
            "text_input": {"query": "cholesterol -78°C Et3N"},
            "output": {"result": {
                "reaction": "Swern oxidation: cholesterol (3β-ol) → cholest-4-en-3-one",
                "alcohol_analysis": {"type": "secondary allylic alcohol (steroid)", "product_type": "α,β-unsaturated ketone"},
                "product": "cholestenone (Δ⁴-3-ketosteroid)",
                "yield": "78-90%",
                "note": "Classic steroid oxidation — Swern preserves acid-sensitive steroid skeleton",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_SWERN_DATA)

    def _run_base(self, alcohol_smiles: str, temperature: str = "-78°C", base: str = "Et3N") -> dict:
        if not alcohol_smiles:
            raise ChemMCPInputError("Alcohol substrate is required.")

        alc = self._analyze_alcohol(alcohol_smiles)
        prod = self._predict_product(alc)
        cond = self._optimize(alc, temperature, base)
        limits = self._relevant_limitations(alc)

        result = {
            "result": {
                "reaction": f"Swern oxidation: {alcohol_smiles} → {prod.get('name','?')}",
                "alcohol_analysis": alc,
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['mechanism']],
                "reagent_details": self.data['reagents'],
                "temperature_profile": self.data['temperature_profile'],
                "product_prediction": prod,
                "advantages_vs_other_oxidations": self.data['advantages_vs_other_oxidations'],
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "typical_yields": self._estimate_yield(alc),
                "workup_procedure": self.data['workup_procedure'],
                "summary": f"Swern oxidation: {alcohol_smiles} → {prod.get('name','?')}. Yield: {self._estimate_yield(alc)}. Key: NO over-oxidation of primary alcohols!",
            }
        }
        logger.info(f"Swern oxidation: {alcohol_smiles}")
        return result

    def _analyze_alcohol(self, smi):
        s = (smi or "").strip().lower()
        alc_types = [
            ('primary aliphatic alcohol', ['1-hexanol', '1-butanol', 'ethanol', r'primary.*alcohol', r'ch2oh'], '→ aldehyde'),
            ('primary benzylic alcohol', ['benzyl alcohol', r'phch2oh', r'benzyl.*alcohol'], '→ aromatic aldehyde (excellent yield)'),
            ('primary allylic alcohol', ['allyl alcohol', r'ch2=chch2oh', r'cinnamyl.*alcohol'], '→ α,β-unsaturated aldehyde'),
            ('secondary aliphatic alcohol', ['2-butanol', 'isopropanol', r'chohol', r'cyclohexanol', r'secondary.*alcohol'], '→ ketone'),
            ('secondary allylic alcohol', ['cholesterol', r'3.*ol.*steroid', r'crotyl.*alcohol', r'secondary.*allylic'], '→ α,β-unsaturated ketone'),
            ('sugar/polyol', ['glucose', r'sugar', r'polyol', r'ribos'], 'Selective oxidation possible'),
        ]
        for atype, pats, prod_hint in alc_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "product_hint": prod_hint}
        return {"type": "unknown_alcohol", "input": smi, "product_hint": "unknown"}

    def _predict_product(self, alc):
        at = alc.get('type', '')
        hint = alc.get('product_hint', '')
        name_map = {
            '→ aldehyde': f"{alc.get('input', '?')} → corresponding aldehyde",
            '→ aromatic aldehyde (excellent yield)': f"{alc.get('input', '?')} → aryl aldehyde (high yield expected)",
            '→ α,β-unsaturated aldehyde': f"{alc.get('input', '?')} → enal (may have E/Z mixture)",
            '→ ketone': f"{alc.get('input', '?')} → corresponding ketone",
            '→ α,β-unsaturated ketone': f"{alc.get('input', '?')} → enone",
        }
        return {
            "name": name_map.get(hint, f"oxidized product from {alc.get('input', '?')}"),
            "class": "aldehyde" if 'aldehyde' in hint else "ketone" if 'ketone' in hint else "carbonyl compound",
            "yield": "80-95%" if 'benzyl' in at else "75-92%" if 'primary' in at else "80-98%" if 'secondary' in at else "70-90%",
        }

    def _optimize(self, alc, temperature, base):
        at = alc.get('type', '')
        return {
            "reagents": {"DMSO": "2.0 eq (anhydrous)", "(COCl)2": "1.2 eq", "base": f"{base} (5 eq)", "solvent": "CH2Cl2 (dry, 0.1-0.5 M)"},
            "temperature": temperature or "−78°C",
            "procedure": (
                "1. Add (COCl)2 (1.2 eq) to DMSO (2 eq) in CH2Cl2 at −78°C dropwise over 5 min\n"
                "2. Stir 5-10 min at −78°C (solution turns orange/yellow)\n"
                "3. Add alcohol (1 eq) in CH2Cl2 dropwise at −78°C\n"
                "4. Stir 15-30 min at −78°C\n"
                "5. Add Et3N (5 eq) dropwise at −78°C\n"
                "6. Allow to warm to RT over 30-60 min with stirring\n"
                "7. Quench with water, extract with CH2Cl2"
            ),
            "atmosphere": "N2 or Ar (anhydrous conditions essential)",
            "time": "2-4 hours total",
            "monitoring": "TLC (disappearance of starting alcohol)",
        }

    def _relevant_limitations(self, alc):
        relevant = []
        relevant.append({"issue": "Low temperature (−78°C)", "problem": "Requires cryogenic cooling; scale-up difficult",
                        "solution": "Dry ice/acetone bath; for large scale consider Dess-Martin oxidation"})
        relevant.append({"issue": "Dimethyl sulfide odor", "problem": "Me2S byproduct has extremely foul smell",
                        "solution": "Efficient fume hood; trap with activated charcoal or aq. CuSO4"})
        at = alc.get('type', '')
        if 'sugar' in at.lower() or 'polyol' in at.lower():
            relevant.append({"issue": "Selectivity", "problem": "Multiple OH groups present — which one gets oxidized?",
                           "Solution": "Protect other OH groups first; use selective conditions"})
        return relevant

    def _estimate_yield(self, alc):
        at = alc.get('type', '')
        if 'benzyl' in at: return "85-96%"
        if 'secondary' in at: return "80-98%"
        if 'primary' in at: return "75-95%"
        if 'allylic' in at: return "70-90%"
        return "75-92%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        alc = parts[0] if parts else ""
        temp = parts[1] if len(parts) > 1 else "-78°C"
        b = parts[2] if len(parts) > 2 else "Et3N"
        return self._run_base(alc, temp, b)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
