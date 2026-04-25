"""
Sharpless Asymmetric Dihydroxylation (Tool #167)
Sharpless 不对称双羟化反应：OsO4/手性配体（DHQD/PHAL或DHQ/PHAL）对烯烃的不对称顺式双羟化。
涵盖：对映选择性（配体选择决定进攻面）、[OsO4(L)]活性物种、
机理（[3+2]环加成→锇酸酯→水解）、AD-mix-α/β商品化试剂、
范围和限制（OsO4毒性）。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_SHARPLESS_DHD_DATA = {
    'mechanism': [
        ('1', 'Catalyst-ligand complex formation', 'OsO4 + chiral ligand (DHQD- or DHQ- PHAL, ALO, etc.) → [OsO4(L)] complex. This chiral environment directs face selectivity.'),
        ('2', '[3+2] Cycloaddition (rate-determining)', '[OsO4(L)] adds across C=C bond in a concerted [3+2] manner → cyclic osmate(VI) ester. The CHIRAL LIGAND determines which alkene face is attacked (Re vs Si).'),
        ('3', 'Hydrolysis', 'Oxidative agent (K3Fe(CN)6 or NMO) hydrolyzes osmate ester → vicinal diol (cis!) + Os(VI) species. STEREOCHEMISTRY: syn addition (both OH on same face).'),
        ('4', 'Reoxidation', 'Os(VI) is reoxidized to OsO4 by co-oxidant (K3Fe(CN)6 in AD-mix; NMO in Upjohn variant). Catalytic in osmium!'),
    ],
    'ligand_systems': {
        'DHQD-PHAL (AD-mix-β)': {
            'full_name': 'Hydroquinidine 1,4-phthalazinediyl diether',
            'face_selectivity': 'Delivers oxygen from the RE face of the alkene (for most terminal and 1,2-disubstituted alkenes)',
            'commercial_form': 'AD-mix-β (contains K2OsO2(OH)4, DHQD-PHAL, K3Fe(CN)6, K2CO3)',
            'typical_ee': '88-99% for good substrates',
        },
        'DHQ-PHAL (AD-mix-α)': {
            'full_name': 'Hydroquinine 1,4-phthalazinediyl diether',
            'face_selectivity': 'Delivers oxygen from the SI face of the alkene (opposite of AD-mix-β)',
            'commercial_form': 'AD-mix-α (contains K2OsO2(OH)4, DHQ-PHAL, K3Fe(CN)6, K2CO3)',
            'typical_ee': '88-99% for good substrates',
        },
        'other_ligands': {
            'DHQD-AL / DHQ-AL': 'Acyloin-type ligands — for certain substrate classes',
            'DHQD-IND / DHQ-IND': 'Indoline-based — for specific steric environments',
            'DHQD-PYR / DHQ-PYR': 'Pyrimidine-based — alternative electronic properties',
            '(DHQD)2PHAL / (DHQ)2PHAL': 'Dimeric ligands — for enhanced selectivity',
        },
        'selection_guide': (
            "For terminal alkenes: AD-mix-α gives one enantiomer, AD-mix-β gives the other.\n"
            "For aryl-substituted alkenes: test both if configuration unknown.\n"
            "For Z-alkenes: generally higher ee than E-alkenes.\n"
            "Rule of thumb: if unsure, run small-scale tests with BOTH AD-mix-α and β."
        ),
    },
    'ad_mix_commercial': {
        'AD-mix-α': {'contents': 'K2OsO2(OH)4 (0.6 mol%), DHQ-PHAL (1.5 mol%), K3Fe(CN)6 (3 eq), K2CO3 (3 eq)',
                     'solvent': 't-BuOH/H2O (1:1)', 'T': '0°C'},
        'AD-mix-β': {'contents': 'K2OsO2(OH)4 (0.6 mol%), DHQD-PHAL (1.5 mol%), K3Fe(CN)6 (3 eq), K2CO3 (3 eq)',
                     'solvent': 't-BuOH/H2O (1:1)', 'T': '0°C'},
        'usage': 'Simply add AD-mix (1.4 g per mmol alkene) to alkene in t-BuOH/H2O at 0°C. Stir until complete.',
        'advantage': 'Pre-measured, convenient, reliable — no need to weigh toxic Os compounds!',
    },
    'scope': [
        'Terminal alkenes: EXCELLENT — 88-99% ee typical; gives primary-secondary 1,2-diols',
        'Aryl alkenes (styrenes): EXCELLENT — high ee; benzylic position gets OH',
        'cis-Alkenes (Z): VERY GOOD — often higher ee than trans counterparts',
        'trans-Alkenes (E): GOOD — moderate to high ee depending on substituents',
        '1,2-Disubstituted alkenes: GOOD — the classic AD substrate class',
        'Trisubstituted alkenes: MODERATE — works but slower; ee may vary',
        'Tetrasubstituted alkenes: POOR/SLOW — very hindered; not recommended',
        'Heteroaryl alkenes: GOOD — furan, thiophene, pyrole derivatives work',
        'Conjugated dienes: SELECTIVE — usually mono-dihydroxylates at less hindered double bond',
    ],
    'limitations': [
        ('OsO4 extreme toxicity', 'OsO4 is VOLATILE (bp 130°C), HIGHLY TOXIC (can penetrate skin, damage eyes, is carcinogenic), and EXPENSIVE', 'Solution: Use AD-mix (catalytic Os); work in fume hood with PPE; never handle pure OsO4 crystals; use NaHSO3 quench'),
        ('Over-oxidation', 'Vicinal diols can be further oxidized to α-hydroxy carbonyls under forcing conditions', 'Solution: Monitor reaction carefully; stop at diol stage; avoid excess oxidant'),
        ('Slow for tetrasubstituted alkenes', 'Very sterically hindered alkenes react very slowly', 'Solution: Use higher T (RT), longer time, more catalyst; consider alternative methods'),
        ('Cost', 'Osmium is a precious metal (~$500/g for OsO4); ligands are also expensive', 'Solution: Use catalytic conditions (0.2-2 mol% Os); AD-mix is cost-effective per reaction'),
        ('Aqueous conditions', 'Standard AD uses t-BuOH/H2O mixture — water-sensitive substrates may be problematic', 'Solution: Use anhydrous conditions with NMO as co-oxidant (Upjohn variant in acetone/H2O)'),
        ('Stereochemistry always cis', 'Only gives syn (cis) diols — cannot make anti diols directly', 'Solution: For anti diols, use epoxide opening or other methods'),
        ('E vs Z difference', 'E-alkenes typically give lower ee than Z-alkenes', 'Solution: Optimize ligand choice; accept lower ee or use different method'),
        ('Regioselectivity for unsymmetrical dienes', 'Which double bond gets dihydroxylated can be unpredictable', 'Solution: Test empirically; use directing groups'),
    ],
    'typical_outcomes': {
        'terminal_alkene': {'ee': '90-99%', 'yield': '75-95%', 'product': 'primary-secondary vicinal diol'},
        'styrene_derivative': {'ee': '88-98%', 'yield': '80-95%', 'product': 'benzylic 1,2-diol'},
        'cis_alkene_Z': {'ee': '92-99%', 'yield': '78-94%', 'product': 'cis-1,2-diol'},
        'trans_alkene_E': {'ee': '70-95%', 'yield': '72-92%', 'product': 'three/erythro diol mixture depends on ligand'},
        'trisubstituted': {'ee': '60-90%', 'yield': '60-85%', 'product': 'trisubstituted 1,2-diol'},
    },
    'workup_quench': (
        "CRITICAL: Quench excess OsO4 with SATURATED AQUEOUS NaHSO3 (sodium metabisulfite)!\n"
        "OsO4 is reduced to non-volatile Os(VI)/Os(VII) species safe for disposal.\n"
        "Then extract with EtOAc, wash with water, dry (Na2SO4), concentrate, purify.\n"
        "WARNING: Do NOT discard reaction waste without proper Os quenching!"
    ),
}


@ChemMCPManager.register_tool
class SharplessDihydroxylation(BaseTool):
    __version__ = "0.1.0"
    name = "SharplessDihydroxylation"
    func_name = 'analyze_sharpless_dihydroxylation'
    description = "Sharpless asymmetric dihydroxylation analysis: OsO4/chiral ligand-catalyzed enantioselective syn-dihydroxylation of alkenes. Covers mechanism (4 steps: catalyst-ligand complex, [3+2] cycloaddition, hydrolysis, reoxidation), ligand systems (DHQD-PHAL/AD-mix-β vs DHQ-PHAL/AD-mix-α with face-selectivity rules), commercial AD-mix details, scope (9 categories from terminal to conjugated dienes), limitations (8 items including OsO4 toxicity/handling), typical outcomes (ee/yield table by alkene type), and safety-critical quench procedure."
    implementation_description = "Comprehensive knowledge base covering: 4-step mechanism emphasizing syn addition stereochemistry, detailed ligand comparison (DHQD vs DHQ series with selection guide), AD-mix-α/β commercial reagent composition and usage, 9 scope categories, 8 limitations with solutions, outcome benchmarks table (5 alkene types), and mandatory Os quenching protocol."
    categories = ["Reaction"]
    tags = ["Sharpless", "Asymmetric Dihydroxylation", "Osmium", "Diol", "Enantioselectivity", "Alkene", "Catalysis", "AD-mix"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("alkene_smiles", "str", "N/A", "SMILES or name of the alkene substrate."),
        ("ad_mix_type", "str", "AD-mix-β", "AD-mix type: 'AD-mix-α' or 'AD-mix-β'."),
        ("solvent", "str", "t-BuOH/H2O", "Solvent system."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: alkene [ad_mix_type] [solvent]. E.g., 'styrene AD-mix-β t-BuOH/H2O'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing alkene_analysis, predicted_product, mechanism_steps, ligand_details, ad_mix_info, optimal_conditions, scope, limitations, typical_outcomes, safety_notes, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"alkene_smiles": "styrene", "ad_mix_type": "AD-mix-β", "solvent": "t-BuOH/H2O"},
            "text_input": {"query": "styrene AD-mix-β t-BuOH/H2O"},
            "output": {"result": {
                "reaction": "Sharpless AD: styrene → (R)-1-phenyl-1,2-ethanediol (with AD-mix-β)",
                "alkene_analysis": {"type": "terminal aryl alkene (styrene)", "substitution": "monosubstituted"},
                "predicted_product": "(R)-1-phenyl-1,2-ethanediol (or opposite enantiomer with AD-mix-α)",
                "stereochemistry": "syn (cis) addition — both OH on same face",
                "reagents": "AD-mix-β (K2OsO2(OH)4/DHQD-PHAL/K3Fe(CN)6/K2CO3)",
                "conditions": {"T": "0°C", "solvent": "t-BuOH/H2O (1:1)", "time": "6-48 h"},
                "yield": "82-95%",
                "ee": "90-98%",
            }},
        },
        {
            "code_input": {"alkene_smiles": "trans-β-methylstyrene", "ad_mix_type": "AD-mix-α", "solvent": "t-BuOH/H2O"},
            "text_input": {"query": "trans-beta-methylstyrene AD-mix-alpha t-BuOH/H2O"},
            "output": {"result": {
                "reaction": "Sharpless AD: (E)-β-methylstyrene → (1S,2R)-1-phenylpropane-1,2-diol (with AD-mix-α)",
                "alkene_analysis": {"type": "trans (E) 1,2-disubstituted aryl alkene", "substitution": "disubstituted E"},
                "stereochemistry": "syn (cis) diol",
                "note": "E-alkenes give lower ee than Z-alkenes typically; may need optimization",
                "yield": "75-89%",
                "ee": "72-92%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_SHARPLESS_DHD_DATA)

    def _run_base(self, alkene_smiles: str, ad_mix_type: str = "AD-mix-β", solvent: str = "t-BuOH/H2O") -> dict:
        if not alkene_smiles:
            raise ChemMCPInputError("Alkene substrate is required.")

        alk = self._analyze_alkene(alkene_smiles)
        prod = self._predict_product(alk, ad_mix_type)
        cond = self._optimize(alk, ad_mix_type, solvent)
        limits = self._relevant_limitations(alk)

        result = {
            "result": {
                "reaction": f"Sharpless AD: {alkene_smiles} → {prod.get('name','?')} ({ad_mix_type})",
                "alkene_analysis": alk,
                "predicted_product": prod,
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['mechanism']],
                "ligand_details": self.data['ligand_systems'],
                "ad_mix_info": self.data['ad_mix_commercial'],
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "typical_outcomes": self._data_lookup_outcome(alk),
                "safety_critical": "QUENCH WITH SAT. aq. NaHSO3 BEFORE WORKUP! OsO4 is extremely toxic.",
                "summary": f"Sharpless AD: {prod.get('name','?')}. Stereochemistry: syn (cis). Yield: {prod.get('yield','?')}. EE: {prod.get('ee','?')}. ⚠️ OsO4 HAZARD!",
            }
        }
        logger.info(f"Sharpless AD: {alkene_smiles}")
        return result

    def _analyze_alkene(self, smi):
        s = (smi or "").strip().lower()
        alk_types = [
            ('terminal alkene (aliphatic)', ['1-hexene', r'propene', r'butene', r'terminal.*alkene', r'=ch2'], 'Monosubstituted; excellent AD substrate'),
            ('terminal aryl alkene (styrene type)', ['styrene', r'vinylbenzene', r'phch=ch2'], 'Aryl-terminal; excellent ee expected'),
            ('cis (Z) 1,2-disubstituted', ['cis-2-butene', r'(z)-', r'cis.*alkene', r'z-alkene'], 'Z geometry; higher ee than E'),
            ('trans (E) 1,2-disubstituted', ['trans-2-butene', r'(e)-', r'trans.*alkene', r'e-alkene'], 'E geometry; moderate ee'),
            ('trisubstituted alkene', ['trimethyl ethylene', r'trisubstituted', r'(ch3)2c=chch3'], 'More hindered; slower'),
            ('tetrasubstituted alkene', ['tetramethyl ethylene', r'tetrasubstituted'], 'Very slow; not recommended'),
            ('conjugated diene', ['1,3-butadiene', r'diene', r'conjugated.*diene'], 'May give mono-dihydroxylation'),
            ('aryl internal alkene', ['β-methylstyrene', r' stilbene', r'phch=chph'], 'Internal aryl alkene'),
        ]
        for atype, pats, note in alk_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "note": note}
        return {"type": "unknown_alkene", "input": smi}

    def _predict_product(self, alk, ad_mix):
        at = alk.get('type', '')
        am = (ad_mix or "AD-mix-β").strip()
        name_map = {
            'terminal': f"{alk.get('input','?')} → chiral 1,2-diol (vicinal diol, syn addition)",
            'styrene': f"{alk.get('input','?')} → chiral 1-phenyl-1,2-ethanediol",
            'cis': f"{alk.get('input','?')} → chiral cis-1,2-diol (high ee expected)",
            'trans': f"{alk.get('input','?')} → chiral three/erythro 1,2-diol (moderate ee)",
            'trisubstituted': f"{alk.get('input','?')} → chiral trisubstituted 1,2-diol",
            'diene': f"{alk.get('input','?')} → mono-dihydroxylated product (usually at less hindered C=C)",
        }
        key = 'terminal' if 'terminal' in at else 'styrene' if 'styrene' in at else 'cis' if 'z' in at or 'cis' in at else 'trans' if 'e' in at or 'trans' in at else 'trisubstituted' if 'tri' in at else 'diene' if 'diene' in at else 'general'
        outcomes = self.data['typical_outcomes']
        outcome_key = 'terminal_alkene' if 'terminal' in at else 'styrene_derivative' if 'styrene' in at else 'cis_alkene_Z' if 'cis' in at else 'trans_alkene_E' if 'trans' in at else 'trisubstituted' if 'tri' in at else None
        oc = outcomes.get(outcome_key, {}) if outcome_key else {}
        return {
            "name": name_map.get(key, f"chiral vicinal diol from {alk.get('input','?')}"),
            "stereochemistry": "syn (cis) addition — both OH groups on same face of original double bond",
            "enantiomer_depends_on": f"{'Si-face attack with AD-mix-β' if 'β' in am else 'Re-face attack with AD-mix-α'} (verify experimentally for new substrates)",
            "yield": oc.get('yield', '70-92%'),
            "ee": oc.get('ee', '80-95%'),
        }

    def _optimize(self, alk, ad_mix, solvent):
        at = alk.get('type', '')
        return {
            "reagent": f"{ad_mix} (pre-packaged: K2OsO4(OH)2, chiral ligand, K3Fe(CN)6, K2CO3)",
            "amount": "~1.4 g AD-mix per mmol of alkene (or individual components)",
            "solvent": solvent or "t-BuOH/H2O (1:1 v/v)",
            "temperature": "0°C (ice bath) for best ee; RT possible but lower ee",
            "atmosphere": "Air OK (co-oxidant is K3Fe(CN)6 — air-stable system)",
            "concentration": "0.05-0.2 M in organic phase",
            "time": "6-48 hours (monitor by TLC; some substrates are fast, others slow)",
            "monitoring": "TLC (disappearance of starting alkene); color change (yellow → colorless/brown)",
            "quench": "SATURATED aqueous NaHSO3 (sodium metabisulfite) — stir 30 min to reduce all Os species!",
            "workup": ("After NaHSO3 quench: extract with EtOAc (×3), wash with water, dry (Na2SO4), "
                       "concentrate, purify by column chromatography"),
        }

    def _relevant_limitations(self, alk):
        relevant = []
        relevant.append({"issue": "OsO4 toxicity", "problem": "Extremely toxic, volatile, penetrates skin, expensive",
                        "Solution": "Use AD-mix (only catalytic Os); full PPE; fume hood; quench with NaHSO3"})
        at = alk.get('type', '')
        if 'tetra' in at.lower():
            relevant.append({"issue": "Very slow reaction", "problem": "Tetrasubstituted alkenes are highly hindered",
                           "Solution": "Use higher T (RT), longer time, more catalyst; consider alternative"})
        if 'e' in at.lower() or 'trans' in at.lower():
            relevant.append({"issue": "Lower ee for E-alkenes", "problem": "Trans alkenes typically give lower enantioselectivity",
                           "Solution": "Optimize ligand choice; accept lower ee or use alternative method"})
        return relevant

    def _data_lookup_outcome(self, alk):
        at = alk.get('type', '')
        outcomes = self.data['typical_outcomes']
        if 'terminal' in at: return outcomes.get('terminal_alkene', {})
        if 'styrene' in at: return outcomes.get('styrene_derivative', {})
        if 'cis' in at: return outcomes.get('cis_alkene_Z', {})
        if 'trans' in at or 'e' in at: return outcomes.get('trans_alkene_E', {})
        if 'tri' in at: return outcomes.get('trisubstituted', {})
        return {}

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        alk = parts[0] if parts else ""
        adm = parts[1] if len(parts) > 1 else "AD-mix-β"
        solv = parts[2] if len(parts) > 2 else "t-BuOH/H2O"
        return self._run_base(alk, adm, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
