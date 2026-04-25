"""
Claisen Condensation (Tool #155)
Claisen 缩合反应详解：酯烯醇化、酰基取代机理、Dieckmann 缩合、交叉 Claisen。
Detailed analysis of the Claisen condensation and its variants.
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_CLAISEN_DATA = {
    'ester_types': {
        'ethyl_acetate': {'alpha_h_pKa': 25, 'enolization': 'requires strong base', 'product': 'acetoacetate (β-ketoester)', 'example': 'CH3COOEt self-Claisen → CH3COCH2COOEt'},
        'ethyl_propionate': {'alpha_h_pKa': 25, 'product': 'ethyl 2-methyl-3-oxopentanoate'},
        'ethyl_butyrate': {'alpha_h_pKa': 25, 'product': 'ethyl 2-ethyl-3-oxohexanoate'},
        'ethyl_benzoate': {'alpha_h_pCa': None, 'no_alpha_H': True, 'note': 'No α-H — cannot form enolate; can only act as electrophile in crossed Claisen'},
        'ethyl_formate': {'alpha_h_pKa': None, 'no_alpha_H': True, 'note': 'No α-H; excellent electrophile → gives α-formyl-β-ketoester (or α-ketoaldehyde after hydrolysis/decarbonylation)'},
        'methyl_esters': {'note': 'Similar to ethyl esters; MeOH formed instead of EtOH (lower bp, easier to remove)'},
        'tert-butyl_esters': {'note': 't-Bu ester cannot enolize (no α-H on the alkoxy group matters for transesterification); useful as electrophile only'},
        'diethyl_oxalate': {'alpha_h_pKa': 13, 'very_acidic': True, 'product': 'oxaloacetic ester derivatives'},
        'ethyl_cyanoacetate': {'alpha_h_pKa': 9, 'very_acidic': True, 'product': 'cyanoacetate condensation products'},
    },

    'base_systems': {
        'NaOEt/EtOH': {'type': 'classic', 'strength': 'pKa EtOH ≈ 16', 'conditions': 'reflux in absolute EtOH', 'notes': 'Classic Claisen conditions; base must match ester alkoxy group to prevent mixed esters via transesterification', 'catalytic_vs_stoich': '1 equiv minimum (consumed to form β-ketoester salt)'},
        'NaH/THF': {'type': 'strong, non-nucleophilic', 'strength': 'pKa H2 ≈ 35', 'conditions': '0°C→RT, THF', 'notes': 'Clean enolate formation; no transesterification risk; H2 gas evolution', 'catalytic_vs_stoich': '1 equiv'},
        'NaOMe/MeOH': {'type': 'for methyl esters', 'notes': 'Use with methyl esters (matching alkoxy)'},
        'LDA/THF': {'type': 'kinetic control', 'conditions': '-78°C, THF', 'notes': 'For regioselective enolization of unsymmetrical ketones/esters', 'catalytic_vs_stoich': '1 equiv'},
        'KHMDS/THF': {'type': 'sterically hindered strong base', 'notes': 'Alternative to LDA; K+ counterion gives more reactive enolate'},
    },

    'mechanism_steps': [
        ('1', 'Enolate formation', 'Base (e.g., OEt⁻) abstracts α-proton from ester → resonance-stabilized enolate'),
        ('2', 'Nucleophilic acyl substitution', 'Enolate attacks carbonyl carbon of another ester molecule → tetrahedral intermediate'),
        ('3', 'Elimination of alkoxide', 'Tetrahedral intermediate collapses, expelling -OR (e.g., -OEt) → β-ketoester'),
        ('4', 'Deprotonation of product', 'The β-ketoester is more acidic (pKa ~11) than starting ester → deprotonated by remaining base → stable enolate salt'),
        ('5', 'Acidic workup', 'Add dilute acid → protonate enolate → neutral β-ketoester product'),
    ],

    'variants': {
        'Crossed_Claisen': {'description': 'Two different esters; one has α-H (enolate), one does not (electrophile)', 'key_requirement': 'One ester MUST lack α-H (e.g., ethyl benzoate, ethyl formate, aromatic ester)', 'example': 'EtOAc + EtOBz → benzoylacetate', 'risk': 'If both have α-H → mixture of 4 products'},
        'Dieckmann_condensation': {'description': 'Intramolecular Claisen of diesters → cyclic β-ketoester', 'ring_sizes': '5- and 6-membered rings favored (entropic advantage); 7-membered possible but slower', 'example': 'Diethyl heptanedioate (diethyl pimelate) → ethyl 2-oxocyclopentanecarboxylate (5-membered ring)', 'ring_preference': '5 > 6 >> 7 membered rings'},
        'Claisen-Schmidt': {'note': 'This is actually aldol (ketone/aldehyde + aldehyde), not a true Claisen. Named similarly but different mechanism.'},
        'Acetoacetic_ester_synthesis': {'description': 'Classic route to substituted acetones via alkylation of acetoacetate enolate', 'sequence': '(1) Claisen self-condensation of EtOAc → acetoacetic ester (2) Alkylation (3) Hydrolysis + decarboxylation → methyl ketone'},
        'Malonic_ester_synthesis': {'description': 'Analogous sequence with malonate esters → substituted acetic acids', 'relation': 'Parallel to acetoacetic ester synthesis'},
        'Thorpie-Ziegler_modification': {'description': 'Claisen of dinitriles or ketoesters', 'notes': 'Broader scope beyond simple esters'},
    },

    'scope_limitations': {
        'scope': [
            'Aliphatic esters with α-H work well',
            'Ethyl acetate is the classic substrate → ethyl acetoacetate',
            'Cyclic diesters → Dieckmann cyclization (excellent for 5/6 rings)',
            'Crossed Claisen works cleanly when one ester lacks α-H',
            'β-Ketoesters are versatile synthetic intermediates (alkylation, Knoevenagel, etc.)',
        ],
        'limitations': [
            'Esters without α-H cannot undergo self-Claisen (only as electrophiles)',
            'Both esters having α-H in crossed Claisen → messy mixture',
            'Sterically hindered α-positions (e.g., t-butyl ester α-position) react slowly',
            'Transesterification can compete if base/ester alkoxy groups do not match',
            'Requires anhydrous conditions (base consumed by any water present)',
            'Product must be acidic enough to be fully deprotonated (drives equilibrium)',
        ],
    },
}


@ChemMCPManager.register_tool
class ClaisenCondensation(BaseTool):
    __version__ = "0.1.0"
    name = "ClaisenCondensation"
    func_name = 'analyze_claisen_condensation'
    description = "Detailed analysis of Claisen condensation: ester enolate formation, nucleophilic acyl substitution mechanism, Dieckmann cyclization variant, crossed Claisen considerations, optimal conditions, scope, limitations, and product utility."
    implementation_description = "Comprehensive knowledge base covering: ester classification (ethyl acetate, benzoate, formate types), base systems (NaOEt matching rule, NaH, LDA), stepwise mechanism (5 steps from enolate to β-ketoester), variants (crossed Claisen, Dieckmann, acetoacetic ester synthesis), and practical laboratory guidance."
    categories = ["Reaction"]
    tags = ["Claisen", "Condensation", "Ester", "Enolate", "Acyl Substitution", "C-C Bond Formation", "Beta-Ketoester", "Dieckmann"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("ester1_smiles", "str", "N/A", "SMILES or name of ester 1 (must have α-H for enolate)."),
        ("ester2_smiles", "str", "", "SMILES or name of ester 2 (electrophile). Leave empty for self-Claisen."),
        ("base", "str", "NaOEt", "Base used: NaOEt, NaH, LDA, etc."),
        ("solvent", "str", "EtOH", "Solvent (should match ester alkoxy group for classic Claisen)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: ester1 [ester2] [base] [solvent]. E.g., 'CC(=O)OCC NaOEt EtOH'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing reaction_type, ester_analysis, mechanism_summary, product_structure, conditions, scope, limitations, and variants."),
    ]

    examples = [
        {
            "code_input": {"ester1_smiles": "CC(=O)OCC", "ester2_smiles": "", "base": "NaOEt", "solvent": "EtOH"},
            "text_input": {"query": "CC(=O)OCC NaOEt EtOH"},
            "output": {"result": {
                "reaction_type": "Self-Claisen condensation of ethyl acetate",
                "ester1": {"name": "ethyl acetate", "has_alpha_H": True, "alpha_h_pKa": 25, "role": "both enolate source AND electrophile"},
                "product": "ethyl acetoacetate (ethyl 3-oxobutanoate)",
                "product_smiles": "CC(=O)CC(=O)OCC",
                "product_class": "β-ketoester (1,3-dicarbonyl compound)",
                "mechanism": ["Enolate formation (NaOEt abstracts α-H)", "Nucleophilic attack on another EtOAc carbonyl", "Elimination of EtO⁻", "Deprotonation of β-ketoester (pKa~11)", "Acid workup → neutral product"],
                "optimal_conditions": {"base": "1 eq NaOEt", "solvent": "absolute EtOH", "temperature": "reflux (~78°C)", "time": "2-12 h", "workup": "dilute HCl"},
                "yield": "40-65% (classical procedure)",
                "key_requirement": "Absolute ethanol (anhydrous) — water consumes base",
                "driving_force": "Formation of stabilized β-ketoester enolate (pKa ~11 vs starting ester pKa ~25)",
                "downstream_chemistry": "Alkylation → substituted acetones; hydrolysis/decarboxylation → ketones; Knoevenagel condensations",
            }},
        },
        {
            "code_input": {"ester1_smiles": "CC(=O)OCC", "ester2_smiles": "c1ccccc1C(=O)OCC", "base": "NaOEt", "solvent": "EtOH"},
            "text_input": {"query": "CC(=O)OCC c1ccccc1C(=O)OCC NaOEt EtOH"},
            "output": {"result": {
                "reaction_type": "Crossed Claisen: ethyl acetate (enolate) + ethyl benzoate (electrophile)",
                "clean_crossed": True,
                "reason": "Ethyl benzoate has no α-H → acts exclusively as electrophile",
                "product": "ethyl benzoylacetate (ethyl 3-oxo-3-phenylpropanoate)",
                "yield": "60-80%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_CLAISEN_DATA)

    def _run_base(self, ester1_smiles: str, ester2_smiles: str = "", base: str = "NaOEt", solvent: str = "EtOH") -> dict:
        if not ester1_smiles:
            raise ChemMCPInputError("Ester 1 is required.")

        e1 = self._classify_ester(ester1_smiles)
        e2 = self._classify_ester(ester2_smiles) if ester2_smiles else None

        is_crossed = e2 is not None
        rxn_type = self._rxn_type(e1, e2, is_crossed)
        mech = self.data['mechanism_steps']
        cat_info = self._analyze_base(base)
        product = self._predict_product(e1, e2, is_crossed)
        conditions = self._optimize(e1, e2, base, solvent, is_crossed)

        result = {
            "result": {
                "reaction_type": rxn_type,
                "ester1_analysis": e1,
                "ester2_analysis": e2,
                "base_analysis": cat_info,
                "mechanism_summary": [{"step": s[0], "name": s[1], "description": s[2]} for s in mech],
                "product_prediction": product,
                "optimal_conditions": conditions,
                "scope": self.data['scope_limitations']['scope'],
                "limitations": self.data['scope_limitations']['limitations'],
                "variants": {k: {jk: jv for jk, jv in v.items() if jk != 'description'} for k, v in self.data['variants'].items()},
                "yield_expectation": self._estimate_yield(e1, e2, base, is_crossed),
                "summary": f"{rxn_type}. Product: {product.get('name','?')}. Driving force: β-ketoester stabilization.",
            }
        }
        logger.info(f"ClaisenCondensation: {rxn_type}")
        return result

    def _classify_ester(self, smi):
        s = (smi or "").strip().lower()
        patterns = {
            'ethyl_acetate': [r'cc(=o)occ', r'ethyl.acetate', r'et oac'],
            'ethyl_benzoate': [r'c1ccc.*c(=o)occ', r'ethyl.benzoate', r'et obz'],
            'ethyl_formate': [r'hc(=o)occ', r'ethyl.formate', r'et ohcoo'],
            'ethyl_propionate': [r'ccc(=o)occ', r'ethyl.propionate'],
            'generic_aliphatic_ester': [r'\(=o\)o', r'ester'],
        }
        for etype, pats in patterns.items():
            for pat in pats:
                if re.search(pat, s, re.IGNORECASE):
                    info = dict(self.data['ester_types'].get(etype, {}))
                    info['classified_as'] = etype
                    info['input'] = smi
                    return info
        return {'classified_as': 'unknown_ester', 'input': smi, 'has_alpha_H': True, 'note': f"Unknown ester '{smi}'"}

    def _rxn_type(self, e1, e2, is_crossed):
        if not is_crossed:
            return f"Self-Claisen condensation of {e1.get('input','?')}"
        n1, n2 = e1.get('input','?'), e2.get('input','?')
        if not e2.get('has_alpha_H', True) or not e1.get('has_alpha_H', True):
            return f"Crossed Claisen (clean): {n1} + {n2}"
        return f"Crossed Claisen (potentially messy): {n1} + {n2}"

    def _analyze_base(self, base):
        b = (base or "NaOEt").strip()
        for key, info in self.data['base_systems'].items():
            if key.lower() in b.lower() or b.lower() in key.lower():
                return {"selected": key, **info}
        return {"selected": b, "notes": f"Base '{b}' — ensure it matches ester alkoxy group to avoid transesterification."}

    def _predict_product(self, e1, e2, is_crossed):
        if not is_crossed:
            return {"name": f"β-ketoester from {e1.get('input','?')} self-condensation", "class": "β-ketoester (1,3-dicarbonyl)"}
        n1, n2 = e1.get('input','?'), e2.get('input','?')
        clean = (not e1.get('has_alpha_H', True)) or (e2 and not e2.get('has_alpha_H', True))
        return {"name": f"Crossed β-ketoester from {n1} + {n2}", "clean": clean}

    def _optimize(self, e1, e2, base, solv, is_crossed):
        cond = {
            "base_amount": "1 equivalent (minimum; consumed in reaction)",
            "solvent": f"{solv} (absolute/anhydrous)",
            "temperature": "reflux temperature of solvent",
            "time": "2-12 hours",
            "atmosphere": "N2/Ar (anhydrous)",
            "workup": "Cool, acidify with dilute HCl (pH ~4), extract with organic solvent",
            "critical_note": "MUST be anhydrous — water destroys base and stops reaction",
        }
        if is_crossed and e2 and not e2.get('has_alpha_H', True):
            cond["stoichiometry"] = f"{e1.get('input','?')} (1-2 eq, enolate) + {e2.get('input','?')} (1 eq, electrophile)"
        return cond

    def _estimate_yield(self, e1, e2, base, is_crossed):
        score = 50
        if is_crossed:
            if e2 and not e2.get('has_alpha_H', True): score = 70
            else: score = 30  # messy
        else:
            if e1.get('classified_as') == 'ethyl_acetate': score = 50  # classical yield
        if 'nah' in base.lower(): score += 10
        return f"{score-8}-{score+8}%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        e1 = parts[0] if parts else ""
        e2 = parts[1] if len(parts) > 1 else ""
        b = parts[2] if len(parts) > 2 else "NaOEt"
        s = parts[3] if len(parts) > 3 else "EtOH"
        return self._run_base(e1, e2, b, s)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
