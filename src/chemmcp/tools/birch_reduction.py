"""
Birch Reduction (Tool #168)
Birch 还原反应：碱金属（Na/Li）/液氨/醇体系将芳香环还原为1,4-环己二烯。
涵盖：机理（溶剂化电子转移→自由基阴离子→质子化→再电子转移→烯醇负离子）、
反芳香性回避（1,4-还原模式）、区域选择性（EWG在α位，EDG在β位）、
条件（液氨−33°C）、范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_BIRCH_DATA = {
    'mechanism': [
        ('1', 'Solvated electron transfer', 'Alkali metal (Na or Li) dissolves in liquid NH3 → solvated electron (e−(am)). Electron adds to aromatic π-system → radical anion. This is the key reducing species!'),
        ('2', 'Protonation', 'Proton source (t-BuOH, EtOH) protonates the radical anion → cyclohexadienyl radical. Position of protonation is REGIOSELECTIVE.'),
        ('3', 'Second electron transfer', 'Another solvated electron adds to the radical → cyclohexadienyl anion'),
        ('4', 'Second protonation', 'Proton source protonates the anion → 1,4-cyclohexadiene product (non-conjugated diene). The 1,4-pattern avoids antiaromaticity in the anion intermediate.'),
    ],
    'regioselectivity_rules': {
        'summary': 'The position of substituents in the 1,4-cyclohexadiene product follows predictable patterns:',
        'electron_withdrawing_group_EWG': (
            'EWG (COOR, COR, CN, CONR2, etc.) ends up at the α-position '
            '(the sp³-hybridized carbon adjacent to one remaining double bond). '
            'Reason: EWG stabilizes the radical anion intermediate when it is at the α-position.'
        ),
        'electron_donating_group_EDG': (
            'EDG (alkyl, OR, NR2, etc.) ends up at the β-position '
            '(on one of the double bond carbons). '
            'Reason: EDG destabilizes the radical anion at the α-position; the system avoids putting negative charge near EDG.'
        ),
        'examples': {
            'Benzoic acid / methyl benzoate': '→ 2,5-cyclohexadiene carboxylic acid (COOH at α-position)',
            'Anisole (methoxybenzene)': '→ 1-methoxy-1,4-cyclohexadiene (OMe at β-position, on double bond)',
            'Toluene (methylbenzene)': '→ 1-methyl-1,4-cyclohexadiene (Me at β-position)',
            'Benzonitrile': '→ 2,5-cyclohexadienecarbonitrile (CN at α-position)',
            'Acetophenone': '→ 1-phenylethanol derivative after workup (C=O reduced to enol ether then hydrolyzed)',
        },
    },
    'reaction_conditions': {
        'metal': {'options': ['Na (sodium)', 'Li (lithium)'], 'note': 'Na is most common; Li gives slightly different selectivity; K is too reactive'},
        'solvent': {'primary': 'Liquid ammonia (NH3(l), bp −33°C)', 'alternative': 'Amines (ethylamine, methylamine) for higher boiling systems',
                    'co-solvent': 'THF or Et2O (10-20%) to improve substrate solubility'},
        'proton_source': {'options': ['t-BuOH (tert-butanol)', 'EtOH (ethanol)', 'i-PrOH', 'NH4Cl (workup)'],
                         'amount': '1-10 eq depending on substrate; t-BuOH is common (weaker acid = milder)'},
        'temperature': '−33°C (boiling point of liquid NH3) or lower (dry ice/acetone bath at −78°C possible)',
        'atmosphere': 'Inert (N2/Ar) not strictly required but recommended',
        'apparatus': '3-neck flask with dry ice condenser; liquid NH3 condensed in situ or added from cylinder',
        'color': 'DEEP BLUE color indicates solvated electrons present — this is a VISUAL indicator of active reducing conditions!',
    },
    'scope': [
        'Benzene derivatives: EXCELLENT — the classic Birch substrates',
        'Benzene with EWG (COOR, CN, COR): EXCELLENT — EWG at α-position after reduction',
        'Benzene with EDG (alkyl, OR, NR2): GOOD — EDG at β-position after reduction',
        'Fused aromatics (naphthalene, anthracene): EXCELLENT — often easier than benzene itself',
        'Heteroaromatics: SELECTIVE — pyridine, quinoline can be reduced; furan/thiophene may decompose',
        'Phenols: SPECIAL — phenolates are NOT reduced (negative charge repels electrons); free phenols ARE reduced',
        'Anilines: SPECIAL — similar to phenols; free aniline reduces, anilinate does not',
        'Conjugated enones/enals: CAN be reduced via conjugate (1,4-) reduction pathway',
        'Alkynes: Can be reduced to trans-alkenes under modified Birch conditions',
    ],
    'limitations': [
        ('Liquid ammonia handling', 'NH3(l) boils at −33°C; requires cryogenic equipment; pressure build-up; toxic fumes', 'Solution: Use proper cryogenic setup (dry ice condenser); work in well-ventilated fume hood; consider amine alternatives'),
        ('Phenols and anilines (as phenolates/anilinates)', 'Deprotonated forms carry negative charge — REPELS solvated electrons → NO REDUCTION', 'Solution: Use free (neutral) phenol/aniline; or protect as ether/amide'),
        ('Benzoic acids (as carboxylates)', 'Carboxylate anions are not reduced efficiently', 'Solution: Esterify first; reduce ester; hydrolyze back'),
        ('Functional group incompatibility', 'C-X bonds (halides), azides, some nitro groups, epoxides may be affected by solvated electrons', 'Solution: Check each group\'s compatibility; protect sensitive groups'),
        ('Over-reduction', 'Prolonged reaction time or excess metal can give fully hydrogenated (cyclohexane) products', 'Solution: Monitor carefully; quench with NH4Cl when conversion complete'),
        ('Regiochemistry for polysubstituted benzenes', 'Multiple substituents compete for α/β positions; prediction becomes complex', 'Solution: Use empirical rules; test small scale; literature precedent helps'),
        ('Substrate solubility', 'Many organic compounds have limited solubility in liquid NH3', 'Solution: Add THF/Et2O co-solvent (10-20%); use ultrasonication'),
        ('Safety concerns', 'Na metal + NH3 = explosive hazard if water/O2 present; Na is pyrophoric', 'Solution: Rigorously exclude water/air; add Na slowly; use proper safety equipment'),
    ],
    'typical_yields': {
        'benzene_ewg_substituted': '70-92%',
        'benzene_edg_substituted': '65-88%',
        'naphthalene': '80-95%',
        'anthracene': '85-97%',
        'benzoic_ester': '72-90%',
        'anisole': '60-85%',
        'heteroaromatic': '50-80%',
    },
    'workup_procedure': (
        "After TLC shows completion:\n"
        "(1) Carefully add solid NH4Cl (or saturated aq. NH4Cl) to quench excess metal — CAUTION: H2 evolution!\n"
        "(2) Allow NH3 to evaporate (use cold trap or let warm slowly)\n"
        "(3) Dilute with water, extract with Et2O or EtOAc (×3)\n"
        "(4) Wash with water, brine; dry (Na2SO4); concentrate\n"
        "(5) Purify by column chromatography or distillation\n"
        "NOTE: Blue color should disappear upon quenching (electrons consumed)"
    ),
}


@ChemMCPManager.register_tool
class BirchReduction(BaseTool):
    __version__ = "0.1.0"
    name = "BirchReduction"
    func_name = 'analyze_birch_reduction'
    description = "Birch reduction analysis: alkali metal (Na/Li)/liquid ammonia/proton source-mediated 1,4-reduction of aromatic rings to non-conjugated cyclohexadienes. Covers mechanism (4 steps: solvated e− transfer → radical anion → protonation → second e− → anion → protonation), regioselectivity rules (EWG at α-position, EDG at β-position, with examples for 6 common substituents), detailed conditions (metal choice, liquid NH3, proton sources, deep blue visual indicator), scope (9 categories including fused aromatics and heterocycles), limitations (8 items including liquid NH3 handling, safety), typical yields, and workup procedure."
    implementation_description = "Comprehensive knowledge base covering: 4-step mechanism emphasizing solvated electron chemistry, complete regioselectivity prediction table (EWG vs EDG rules with 6 examples), detailed condition parameters (with safety notes), 9 scope categories, 8 limitations with practical solutions, yield benchmarks, and standardized quench/workup protocol."
    categories = ["Reaction"]
    tags = ["Birch", "Reduction", "Aromatic", "Dissolving Metal", "Liquid Ammonia", "Solvated Electron", "1,4-Cyclohexadiene"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("aromatic_smiles", "str", "N/A", "SMILES or name of the aromatic compound to be reduced."),
        ("metal", "str", "Na", "Reducing metal: 'Na' (sodium) or 'Li' (lithium)."),
        ("proton_source", "str", "t-BuOH", "Proton source: 't-BuOH', 'EtOH', or 'i-PrOH'."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: aromatic [metal] [proton_source]. E.g., 'methyl_benzoate Na t-BuOH'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing aromatic_analysis, predicted_product_regiochemistry, mechanism_steps, conditions_details, regioselectivity_rules, scope, limitations, typical_yields, safety_notes, workup, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"aromatic_smiles": "methyl benzoate", "metal": "Na", "proton_source": "t-BuOH"},
            "text_input": {"query": "methyl_benzoate Na t-BuOH"},
            "output": {"result": {
                "reaction": "Birch reduction: methyl benzoate → methyl 2,5-cyclohexadiene-1-carboxylate",
                "aromatic_analysis": {"type": "benzene with EWG (ester)", "substituent_class": "EWG"},
                "predicted_product": "methyl 2,5-cyclohexadiene-1-carboxylate (COOMe at α-position)",
                "regioselectivity": "EWG (COOMe) ends up at α-position (sp³ carbon next to double bond)",
                "conditions": {"metal": "Na", "solvent": "liq. NH3 / THF", "proton_source": "t-BuOH", "T": "−33°C to −78°C"},
                "visual_indicator": "Deep blue color = solvated electrons present (active!)",
                "yield": "75-90%",
                "product_type": "1,4-cyclohexadiene (non-conjugated diene)",
            }},
        },
        {
            "code_input": {"aromatic_smiles": "anisole", "metal": "Na", "proton_source": "EtOH"},
            "text_input": {"query": "anisole Na EtOH"},
            "output": {"result": {
                "reaction": "Birch reduction: anisole (methoxybenzene) → 1-methoxy-1,4-cyclohexadiene",
                "aromatic_analysis": {"type": "benzene with EDG (methoxy)", "substituent_class": "EDG"},
                "predicted_product": "1-methoxy-1,4-cyclohexadiene (OMe at β-position, on double bond carbon)",
                "regioselectivity": "EDG (OMe) ends up at β-position (on double bond carbon)",
                "yield": "65-85%",
                "note": "EDG-substituted benzenes generally give slightly lower yields than EWG-substituted ones",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_BIRCH_DATA)

    def _run_base(self, aromatic_smiles: str, metal: str = "Na", proton_source: str = "t-BuOH") -> dict:
        if not aromatic_smiles:
            raise ChemMCPInputError("Aromatic substrate is required.")

        arom = self._analyze_aromatic(aromatic_smiles)
        prod = self._predict_product(arom)
        cond = self._optimize(arom, metal, proton_source)
        limits = self._relevant_limitations(arom)

        result = {
            "result": {
                "reaction": f"Birch reduction: {aromatic_smiles} → {prod.get('name','?')}",
                "aromatic_analysis": arom,
                "predicted_product": prod,
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['mechanism']],
                "regioselectivity_rules": self.data['regioselectivity_rules'],
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "typical_yields": self._estimate_yield(arom),
                "safety_notes": "⚠️ Na + liq. NH3 = EXPLOSIVE if H2O present! Pyrophoric Na! Work in fume hood with face shield!",
                "workup_procedure": self.data['workup_procedure'],
                "summary": f"Birch reduction: {prod.get('name','?')} ({prod.get('pattern','1,4-diene')}). Yield: {self._estimate_yield(arom)}. ⚠️ CRYOGENIC + PYROPHORIC!",
            }
        }
        logger.info(f"Birch reduction: {aromatic_smiles}")
        return result

    def _analyze_aromatic(self, smi):
        s = (smi or "").strip().lower()
        arom_types = [
            ('benzene with EWG (ester)', ['methyl benzoate', r'benzoate', r'cooch3.*phenyl', r'c6h5coo'], 'EWG: COOR → α-position'),
            ('benzene with EWG (ketone)', ['acetophenone', r'benzophenone', r'coc6h5', r'comethyl'], 'EWG: COR → α-position'),
            ('benzene with EWG (cyano)', ['benzonitrile', r'benzyl cyanide', r'phcn', r'cyanobenzene'], 'EWG: CN → α-position'),
            ('benzene with EDG (alkyl)', ['toluene', r'methylbenzene', r'xylene', r'c6h5ch3'], 'EDG: alkyl → β-position'),
            ('benzene with EDG (ether)', ['anisole', r'methoxybenzene', r'phoch3', r'phenetole'], 'EDG: OR → β-position'),
            ('fused aromatic (naphthalene)', ['naphthalene', r'naphth'], 'Fused arene — reduces more easily than benzene'),
            ('fused aromatic (anthracene)', ['anthracene'], 'Very easily reduced — often quantitative'),
            ('heteroaromatic (pyridine)', ['pyridine', r'pyridin'], 'N-heteroaromatic — selective reduction possible'),
            ('phenol', ['phenol', r'hydroxybenzene'], 'Free phenol reduces; phenolate does NOT'),
            ('benzoic acid', ['benzoic acid', r'carboxybenzene'], 'Free acid reduces poorly; carboxylate NOT at all'),
        ]
        for atype, pats, note in arom_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "note": note}
        return {"type": "unknown_aromatic", "input": smi}

    def _predict_product(self, arom):
        at = arom.get('type', '')
        inp = arom.get('input', '?')
        if 'EWG' in at:
            return {"name": f"{inp} → 1,4-cyclohexadiene with substituent at α-position", "pattern": "1,4-cyclohexadiene (non-conjugated)",
                    "regiochemistry": "EWG at α-position (sp³ carbon adjacent to C=C)"}
        elif 'EDG' in at:
            return {"name": f"{inp} → 1,4-cyclohexadiene with substituent at β-position", "pattern": "1,4-cyclohexadiene (non-conjugated)",
                    "regiochemistry": "EDG at β-position (on C=C carbon)"}
        elif 'fused' in at.lower():
            return {"name": f"{inp} → partially reduced fused ring (1,4-dihydro derivative)", "pattern": "partial 1,4-reduction",
                    "regiochemistry": "Follows same EWG/EDG rules per ring"}
        elif 'phenol' in at.lower():
            return {"name": f"{inp} → 1,4-cyclohexadienol (free phenol form)", "pattern": "1,4-cyclohexadiene",
                    "regiochemistry": "OH at β-position (EDG behavior)"}
        return {"name": f"reduced product from {inp}", "pattern": "1,4-cyclohexadiene (non-conjugated diene)", "regiochemistry": "depends on substituents"}

    def _optimize(self, arom, metal, proton_source):
        return {
            "metal": metal or "Na" + " (usually 2-4 eq",
            "solvent": "Liquid NH3 (anhydrous) + THF (10-20% co-solvent for solubility)",
            "proton_source": f"{proton_source} (2-10 eq; t-BuOH = mild, EtOH = stronger)",
            "temperature": "−33°C (bp of liq. NH3) to −78°C (dry ice/acetone)",
            "apparatus": "3-neck flask with dry ice condenser; add NH3 gas or condense from cylinder",
            "procedure": (
                "1. Condense NH3 (~50-100 mL per mmol substrate) into flask at −78°C\n"
                "2. Add crystal of Fe(NO3)3 (catalytic — promotes electron transfer)\n"
                "3. Add small piece of Na → DEEP BLUE color appears\n"
                "4. Add substrate in THF solution\n"
                "5. Add Na in small pieces until blue color PERSISTS (excess reducing equivalents)\n"
                "6. Add proton_source dropwise\n"
                "7. Stir until blue color fades (reaction complete) or TLC confirms\n"
                "8. Quench carefully with solid NH4Cl"
            ),
            "time": "30 min - 4 h",
            "monitoring": "Blue color intensity; TLC (disappearance of starting material)",
        }

    def _relevant_limitations(self, arom):
        relevant = []
        relevant.append({"issue": "Liquid ammonia handling", "problem": "Cryogenic (bp −33°C); pressure hazard; toxic fumes",
                        "Solution": "Proper cryogenic setup; fume hood; cold trap"})
        at = arom.get('type', '')
        if 'acid' in at.lower() or 'carboxyl' in at.lower():
            relevant.append({"issue": "Carboxylate not reduced", "problem": "Deprotonated carboxylic acid doesn't reduce",
                           "Solution": "Esterify first, reduce, then hydrolyze"})
        if 'phenol' in at.lower():
            relevant.append({"issue": "Phenol vs phenolate", "problem": "Only free phenol reduces; basic conditions give unreactive phenolate",
                           "Solution": "Ensure neutral/n acidic conditions; don't add excess base"})
        relevant.append({"issue": "Metal + NH3 safety", "problem": "Na is pyrophoric; H2 evolution on quenching",
                        "Solution": "Exclude water rigorously; quench slowly with NH4Cl; wear face shield"})
        return relevant

    def _estimate_yield(self, arom):
        at = arom.get('type', '')
        if 'anthracene' in at: return "85-97%"
        if 'naphthalene' in at: return "80-95%"
        if 'EWG' in at: return "70-92%"
        if 'EDG' in at: return "65-88%"
        if 'hetero' in at: return "50-80%"
        return "70-88%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        arom = parts[0] if parts else ""
        met = parts[1] if len(parts) > 1 else "Na"
        ps = parts[2] if len(parts) > 2 else "t-BuOH"
        return self._run_base(arom, met, ps)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
