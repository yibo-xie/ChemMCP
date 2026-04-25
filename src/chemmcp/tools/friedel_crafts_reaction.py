"""
Friedel-Crafts Reaction (Tool #159)
Friedel-Crafts 烷基化/酰基化反应：亲电芳香取代机理、
催化剂效应（AlCl3, FeCl3）、底物限制、取向规则、多烷基化问题。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_FC_DATA = {
    'alkylation': {
        'mechanism': [
            ('1', 'Lewis acid activation', 'R-X + AlCl3 → R⁺ (or tight ion pair) + AlCl4⁻ (carbocation formation)'),
            ('2', 'Electrophilic attack', 'Arenium ion (σ complex) formation — electrophile attacks aromatic ring'),
            ('3', 'Deprotonation', 'Base removes H⁺ from arenium ion → restores aromaticity → alkylated arene'),
            ('4', 'Catalyst regeneration', 'H-AlCl4 → HCl + AlCl3 (catalyst regenerated in principle)'),
        ],
        'electrophiles': [
            ('R-Cl + AlCl3 (1 eq)', 'Primary/secondary carbocation; may rearrange', 'Classic FC alkylation'),
            ('R-Br, R-I', 'More reactive than R-Cl but more expensive', 'Similar mechanism'),
            ('R-OH + acid (H2SO4, HF, H3PO4)', 'Protonation → R⁺(H2O) → R⁺', 'FC-type using alcohols (needs strong acid)'),
            ('Alkene + acid (HF, H2SO4, H3PO4)', 'Protonation of alkene → carbocation', 'Industrial method (e.g., cumene process)'),
            ('Epoxide + Lewis acid', 'Ring opens → β-cation electrophile', 'Gives 2-arylalcohols after workup'),
        ],
        'limitations': [
            ('Carbocation rearrangement', 'Primary R groups rearrange via hydride/alkyl shifts to more stable cations', 'Use acylation + reduction to avoid rearrangement'),
            ('Polyalkylation', 'Product is more electron-rich than starting material → further alkylation', 'Use excess arene (>10 eq); cannot fully prevent for very activated rings'),
            ('Deactivation', 'Strongly deactivated rings (nitrobenzene, pyridine) do NOT react', 'Nitro group strongly deactivates and coordinates to Lewis acid'),
            ('Meta-directors fail', 'EWG substituents direct meta but ring is too deactivated', 'FC does not work on m-directing substrates generally'),
            ('Ortho steric hindrance', 'Bulky electrophiles give mostly para product', 'Steric effects favor less hindered positions'),
            ('No reaction with -NH2', 'Amino group complexes with AlC13 → deactivation', 'Protect as -NHAc (acetamide) before FC'),
            ('Halobenzenes sluggish', 'Halogen is weakly deactivating but ortho/para directing', 'Very slow; needs more forcing conditions'),
        ],
        'orientation': {
            'o/p_directors_activated': ['—OH, —OR, —NHR, —NR2, —alkyl, —Ph', 'Strongly activated → fast reaction, polyalkylation likely'],
            'o/p_directors_mild': ['—F, —Cl, —Br, —I, —CH=CH2', 'Weakly deactivating/halogen: slow reaction'],
            'meta_directors_deactivated': ['—NO2, —CN, —SO3H, —CHO, —COR, —COOR, —COOH', 'Ring too deactivated → no FC reaction (usually)'],
            'steric_control': 'Ortho:para ratio depends on size of electrophile and existing substituent',
        },
    },

    'acylation': {
        'mechanism': [
            ('1', 'Acyl cation formation', 'RCOCl + AlCl3 (≥1 eq, often >1) → R-C≡O⁺ (acylium ion, resonance stabilized)'),
            ('2', 'Electrophilic attack', 'Arenium ion (σ complex) formation'),
            ('3', 'Deprotonation', 'Restore aromaticity → aryl ketone'),
            ('4', 'Note', 'AlCl3 COMPLEXES with ketone product → needs ≥1 eq AlCl3 (often 1.1+ eq); hydrolyzed during workup'),
        ],
        'advantages_over_alkylation': [
            'No carbocation rearrangement (acylium ion is resonance-stabilized)',
            'Monoacylation only (product is DEACTIVATED toward further acylation)',
            'Ketone product can be reduced to alkyl (Clemmensen or Wolff-Kishner) = "indirect alkylation"',
            'Better regiocontrol (no rearrangement artifacts)',
        ],
        'electrophiles': [
            ('Acid chloride (RCOCl) + AlCl3 (1.1 eq)', 'Most common; forms acylium ion cleanly', 'Standard FC acylation'),
            ('Acid anhydride ((RCO)2O) + AlCl3 (2+ eq)', 'Cheaper than some acid chlorides; second acyl group also reacts', 'Gives same ketone + carboxylic acid byproduct'),
            ('Carboxylic acid + very strong acid', 'Limited scope; possible with superacids', 'Not common'),
            ('Mixed anhydrides', 'Can be used', 'Specialty applications'),
        ],
        'limitations': [
            ('Deactivated rings unreactive', 'Same as alkylation: nitrobenzene etc. no go', 'Reduce or remove deactivator first'),
            ('Amino group interference', '—NH2 complexes AlCl3', 'Protect as amide'),
            ('Steric hindrance', 'Very bulky acyl groups react slowly', 'Use excess reagent or higher T'),
            ('Over-reduction risk', 'If reducing ketone to alkyl afterward', 'Choose appropriate reduction method'),
        ],
    },

    'catalysts': {
        'AlCl3': {'strength': 'very strong', 'stoich': 'catalytic (alkylation) or ≥1 eq (acylation)', 'scope': 'universal for FC', 'handling': 'moisture sensitive, exothermic on addition', 'notes': 'Gold standard FC catalyst'},
        'FeCl3': {'strength': 'strong', 'stoich': 'catalytic', 'scope': 'good for alkylations, milder than AlCl3', 'notes': 'Easier handling than AlCl3'},
        'BF3·OEt2': {'strength': 'moderate-strong', 'stoich': 'catalytic', 'scope': 'milder reactions, sensitive substrates', 'notes': 'Liquid, easier to handle'},
        'ZnCl2': {'strength': 'moderate', 'scope': 'milder FC, specific applications', 'notes': 'Used in Gattermann-Koch formylation with CO/HCl'},
        'H2SO4 (conc.)': {'strength': 'strong Brønsted', 'scope': 'Alcohols/alkenes as electrophiles', 'notes': 'Non-Lewis-acid FC variant'},
        'HF (anhydrous)': {'strength': 'very strong', 'scope': 'Excellent solvent + catalyst', 'safety': 'VERY hazardous (causes severe burns)', 'notes': 'Used industrially'},
        'ion_exchange_resins': {'strength': 'mild', 'scope': 'Green chemistry approach', 'notes': 'Heterogeneous, recyclable catalyst'},
        'Zeolites': {'scope': 'Shape-selective FC catalysis', 'notes': 'Industrial application (ethylation of benzene → ethylbenzene)'},
    },

    'comparison': {
        'FC_alkylation_vs_acylation': {
            'rearrangement': 'Alkylation: YES (carbocation); Acylation: NO (acylium stable)',
            'poly_substitution': 'Alkylation: YES (product activated); Acylation: NO (product deactivated)',
            'catalyst_amount': 'Alkylation: catalytic; Acylation: stoichiometric (complexes product)',
            'product_type': 'Alkyl arene vs aryl ketone',
            'indirect_route': 'Acylation → Clemmensen/WK = rearrangement-free alkylation',
        }
    },
}


@ChemMCPManager.register_tool
class FriedelCraftsReaction(BaseTool):
    __version__ = "0.1.0"
    name = "FriedelCraftsReaction"
    func_name = 'analyze_friedel_crafts_reaction'
    description = "Friedel-Crafts alkylation and acylation analysis: electrophilic aromatic substitution mechanism (σ complex), catalyst comparison (AlCl3, FeCl3, BF3, zeolites), substrate scope and limitations (deactivating groups, -NH2 issue, polyalkylation), orientation rules (ortho/para vs meta directing), and indirect alkylation strategy."
    implementation_description = "Comprehensive knowledge base covering both FC alkylation (5 electrophile types, 7 limitations, orientation rules) and FC acylation (mechanism with stoichiometry note, 4 advantages over alkylation, 4 electrophile types, 4 limitations), 8 catalyst options with comparison, and practical guidance including the acylation→reduction indirect alkylation strategy."
    categories = ["Reaction"]
    tags = ["Friedel-Crafts", "Electrophilic Aromatic Substitution", "Alkylation", "Acylation", "Aromatic", "Lewis Acid"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("arene_smiles", "str", "N/A", "SMILES or name of the arene (aromatic compound)."),
        ("electrophile_type", "str", "alkyl", "Electrophile type: 'alkyl' or 'acyl'."),
        ("electrophile_spec", "str", "CH3Cl", "Specific electrophile (e.g., CH3Cl, CH3COCl, C2H5OH)."),
        ("catalyst", "str", "AlCl3", "Lewis acid catalyst."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: arene [alkyl|acyl] [electrophile] [catalyst]. E.g., 'benzene alkyl CH3Cl AlCl3' or 'toluene acyl CH3COCl AlCl3'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing arene_analysis, electrophile_analysis, mechanism, orientation, product_prediction, limitations, fc_alkylation_vs_acylation, conditions, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"arene_smiles": "benzene", "electrophile_type": "alkyl", "electrophile_spec": "CH3CH2Cl", "catalyst": "AlCl3"},
            "text_input": {"query": "benzene alkyl CH3CH2Cl AlCl3"},
            "output": {"result": {
                "reaction": "FC Alkylation: benzene + ethyl chloride → ethylbenzene",
                "arene_analysis": {"name": "benzene", "activation": "neutral (no substituents)", "reactivity": "moderate"},
                "electrophile": {"type": "alkyl", "spec": "ethyl chloride", "carbocation": "ethyl cation (may rearrange? No, 2° is stable)"},
                "product": "ethylbenzene (C6H5CH2CH3)",
                "orientation": "No directing effect (unsubstituted) → single product",
                "polyalkylation_risk": "MODERATE — ethylbenzene is slightly activated → some diethylation possible",
                "rearrangement_risk": "LOW — ethyl cation (2°) is reasonably stable",
                "yield": "60-80%",
                "conditions": {"catalyst": "AlCl3 (catalytic)", "T": "40-80°C (gentle heating)", "solvent": "CS2 (traditional) or neat benzene", "time": "1-4 h"},
            }},
        },
        {
            "code_input": {"arene_smiles": "toluene", "electrophile_type": "acyl", "electrophile_spec": "CH3COCl", "catalyst": "AlCl3"},
            "text_input": {"query": "toluene acyl CH3COCl AlCl3"},
            "output": {"result": {
                "reaction": "FC Acylation: toluene + acetyl chloride → methyl phenyl ketones (ortho + para)",
                "arene_analysis": {"name": "toluene", "substituent": "-Me (o/p director, activating)", "reactivity": "faster than benzene"},
                "product": "ortho- and para-methylacetophenone mixture (para major due to sterics)",
                "orientation": "o/p directing (-Me); para favored for steric reasons (~60% para)",
                "polyacylation_risk": "NONE — ketone product is meta-directing/deactivating",
                "rearrangement_risk": "NONE — acylium ion is resonance-stabilized",
                "yield": "70-85%",
                "catalyst_note": "AlCl3 complexes with ketone product → use 1.1+ eq AlCl3; hydrolyze during workup",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_FC_DATA)

    def _run_base(self, arene_smiles: str, electrophile_type: str = "alkyl", electrophile_spec: str = "CH3Cl", catalyst: str = "AlCl3") -> dict:
        if not arene_smiles:
            raise ChemMCPInputError("Arene SMILES/name is required.")

        etype = electrophile_type.strip().lower()
        if etype not in ('alkyl', 'acyl'):
            etype = 'alkyl'

        arene = self._analyze_arene(arene_smiles)
        elec = self._analyze_electrophile(etype, electrophile_spec)
        cat = self._analyze_catalyst(catalyst)
        orient = self._predict_orientation(arene, etype)
        product = self._predict_product(arene, elec, etype)
        limits = self._check_limitations(arene, etype)
        cond = self._optimize(arene, elec, etype, catalyst)

        result = {
            "result": {
                "reaction": f"FC {etype.capitalize()}ation: {arene_smiles} + {electrophile_spec} → {product.get('name','?')}",
                "arene_analysis": arene,
                "electrophile_analysis": elec,
                "catalyst_analysis": cat,
                "mechanism": self.data[f'{etype}lation']['mechanism'] if f'{etype}lation' in self.data else [],
                "orientation": orient,
                "product_prediction": product,
                "applicable_limitations": limits,
                "fc_alkylation_vs_acylation": self.data['comparison']['FC_alkylation_vs_acylation'],
                "optimal_conditions": cond,
                "summary": f"FC {etype}lation: {product.get('name','?')}. Yield: {product.get('yield','50-80%')}. Key: {limits[0]['issue'] if limits else 'standard'}",
            }
        }
        logger.info(f"FriedelCrafts: {etype}lation of {arene_smiles}")
        return result

    def _analyze_arene(self, smi):
        s = (smi or "").strip().lower()
        arenes = {
            'benzene': {'patterns': ['benzene', 'c1ccccc1'], 'substituent': None, 'activation': 'neutral', 'directing': 'none', 'reactivity': 'baseline'},
            'toluene': {'patterns': ['toluene', 'c1ccccc1(C)', 'methylbenzene'], 'substituent': '-Me', 'activation': 'activating (o/p)', 'directing': 'o/p', 'reactivity': '>benzene'},
            'anisole': {'patterns': ['anisole', 'methoxybenzene', 'c1ccccc1(OC)'], 'substituent': '-OMe', 'activation': 'strongly activating (o/p)', 'directing': 'o/p', 'reactivity': '>>benzene'},
            'phenol': {'patterns': ['phenol', 'c1ccccc1(O)'], 'substituent': '-OH', 'activation': 'strongly activating (o/p)', 'directing': 'o/p', 'reactivity': '>>benzene', 'note': '-OH may complex with AlCl3; protect if needed'},
            'chlorobenzene': {'patterns': ['chlorobenzene', 'c1ccccc1(Cl)'], 'substituent': '-Cl', 'activation': 'weakly deactivating (o/p)', 'directing': 'o/p', 'reactivity': '<benzene'},
            'nitrobenzene': {'patterns': ['nitrobenzene', 'c1ccccc1([N+](=O)[O-])'], 'substituent': '-NO2', 'activation': 'strongly deactivating (m)', 'directing': 'm', 'reactivity': 'NO REACTION (too deactivated)', 'note': 'FC does NOT work on nitrobenzene'},
            'xylene': {'patterns': ['xylene', 'dimethylbenzene'], 'substituent': '2×-Me', 'activation': 'strongly activating', 'reactivity': '>>>benzene', 'note': 'Polyalkylation very likely'},
            'naphthalene': {'patterns': ['naphthalene'], 'substituent': 'fused ring', 'activation': 'more reactive than benzene', 'directing': 'α-position favored', 'reactivity': '>benzene'},
        }
        for aname, info in arenes.items():
            for pat in info.get('patterns', []):
                if pat in s or re.search(pat, s):
                    return dict(info)
        return {'name': smi, 'activation': 'unknown', 'directing': 'unknown'}

    def _analyze_electrophile(self, etype, spec):
        s = (spec or "").strip()
        return {
            "type": etype,
            "spec": spec,
            "active_species": "carbocation (R⁺)" if etype == 'alkyl' else "acylium ion (R-C≡O⁺)",
            "rearrangement_possible": etype == 'alkyl',
            "stability": "resonance-stabilized" if etype == 'acyl' else 'depends on structure (1° rearranges, 2°/3° stable)',
        }

    def _analyze_catalyst(self, cat):
        c = (cat or "AlCl3").strip()
        for key, info in self.data['catalysts'].items():
            if key.lower() in c.lower() or c.lower() in key.lower():
                return {"selected": key, **info}
        return {"selected": c, "notes": f"Catalyst '{c}' — verify suitability for this transformation"}

    def _predict_orientation(self, arene, etype):
        sub = arene.get('substituent')
        if not sub:
            return {"pattern": "monosubstitution (all positions equivalent for benzene)", "isomer_distribution": "single product"}
        d = arene.get('directing', '?')
        act = arene.get('activation', '')
        o_p_ratio = "ortho:para ≈ 30:60 (some meta)" if d == 'o/p' else "meta only"
        return {
            "directing_effect": d,
            "activation_level": act,
            "expected_isomers": f"{d}-substituted products" if d != 'none' else "single substitution product",
            "ratio_note": o_p_ratio if d == 'o/p' else "",
        }

    def _predict_product(self, arene, elec, etype):
        aname = arene.get('name', '?')
        espec = elec.get('spec', '?')
        r_group = re.sub(r'(Cl|Br|I|OH)$', '', espec).strip() or 'R'
        if etype == 'alkyl':
            name = f"{r_group}-substituted {aname}" if aname != '?' else f"alkylated arene"
        else:
            name = f"{r_group}-carbonyl {aname} (aryl ketone)"
        yield_est = "60-85%" if etype == 'acyl' else "50-80%"
        return {"name": name, "class": "alkyl arene" if etype == 'alkyl' else "aryl ketone", "yield": yield_est}

    def _check_limitations(self, arene, etype):
        relevant = []
        act = arene.get('activation', '')
        sub = arene.get('substituent')

        if 'deactivating' in act or 'NO' in (sub or ''):
            relevant.append({"issue": "Deactivated ring", "problem": f"Ring too deactivated ({sub}) for FC reaction", "solution": "Remove or reduce deactivating group first"})
        if etype == 'alkyl':
            relevant.append({"issue": "Polyalkylation", "problem": "Product more activated than starting material", "solution": "Use large excess of arene (10+ eq)"})
            relevant.append({"issue": "Rearrangement", "problem": "Carbocation may rearrange (especially 1° R groups)", "solution": "Use acylation + reduction (Clemmensen/WK) instead"})
        if sub and 'NH2' in sub:
            relevant.append({"issue": "Amino group interference", "problem": "-NH2 complexes with Lewis acid catalyst", "solution": "Protect as acetamide (-NHAc) before FC"})
        if not relevant:
            relevant.append({"issue": "Standard precautions apply", "problem": "Ensure anhydrous conditions; control exotherm", "solution": "Add catalyst slowly under cooling"})
        return relevant

    def _optimize(self, arene, elec, etype, catalyst):
        cond = {
            "catalyst_loading": "catalytic (0.1-0.2 eq)" if etype == 'alkyl' else "≥1.1 eq (complexes with ketone product)",
            "solvent": "CS2, CH2Cl2, CH3NO2, or neat arene (excess as solvent)",
            "temperature": "RT to gentle reflux (0-80°C typical)",
            "atmosphere": "anhydrous (exclude moisture — deactivates Lewis acid)",
            "addition_order": "Add catalyst to arene solution, then add electrophile slowly with cooling",
            "workup": "Pour onto ice/water carefully (EXOTHERMIC!); extract with organic solvent; wash, dry, purify",
            "safety": "HIGHLY EXOTHERMIC reaction initiation — add slowly with good cooling and stirring!",
        }
        return cond

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        arene = parts[0] if parts else ""
        etype = parts[1] if len(parts) > 1 else "alkyl"
        espec = parts[2] if len(parts) > 2 else "CH3Cl"
        cat = parts[3] if len(parts) > 3 else "AlCl3"
        return self._run_base(arene, etype, espec, cat)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
