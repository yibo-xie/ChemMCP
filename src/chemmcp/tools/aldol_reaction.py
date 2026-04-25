"""
Aldol Reaction (Tool #154)
羟醛缩合反应：底物范围、反应条件、产物预测、限制因素、优化策略。
Comprehensive analysis of aldol addition and condensation reactions.
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# 羟醛反应综合知识库
_ALDOL_DATA = {
    # 底物类型 → α-H酸性、烯醇化难易、亲电活性
    'substrate_types': {
        'aliphatic_aldehyde': {
            'alpha_h_pKa': 17,
            'enolization': 'easy',
            'electrophile_reactivity': 'high',
            'self_aldol_tendency': 'high',
            'examples': ['acetaldehyde (CH3CHO)', 'propionaldehyde (CH3CH2CHO)', 'butyraldehyde (CH3(CH2)2CHO)'],
        },
        'aromatic_aldehyde': {
            'alpha_h_pKa': None,  # 通常无α-H
            'enolization': 'none (no α-H)',
            'electrophile_reactivity': 'very high',
            'self_aldol_tendency': 'none',
            'note': 'Benzaldehyde has no α-H — excellent electrophile only for crossedaldol',
            'examples': ['benzaldehyde (C6H5CHO)', 'p-tolualdehyde', 'p-anisaldehyde'],
        },
        'aliphatic_ketone': {
            'alpha_h_pKa': 20,
            'enolization': 'moderate (needs stronger base)',
            'electrophile_reactivity': 'moderate (less than aldehyde)',
            'self_aldol_tendency': 'moderate',
            'examples': ['acetone ((CH3)2CO)', 'cyclohexanone', '2-butanone (MEK)'],
        },
        'cyclic_ketone': {
            'alpha_h_pKa': 19-20,
            'enolization': 'easy (ring strain helps)',
            'electrophile_reactivity': 'moderate',
            'self_aldol_tendency': 'good (forms fused rings)',
            'examples': ['cyclopentanone', 'cyclohexanone', 'cycloheptanone'],
        },
        'fluorinated_ketone': {
            'alpha_h_pKa': 9-12,
            'enolization': 'very easy (strong -I effect of F)',
            'electrophile_reactivity': 'reduced',
            'examples': ['1,1,1-trifluoroacetone'],
        },
        'β-ketoester': {
            'alpha_h_pKa': 11,
            'enolization': 'very easy (stabilized enolate)',
            'electrophile_reactivity': 'N/A (usually nucleophile)',
            'examples': ['ethyl acetoacetate', 'acetylacetone'],
        },
        'β-diketone': {
            'alpha_h_pKa': 9,
            'enolization': 'very easy',
            'examples': ['acetylacetone (2,4-pentanedione)'],
        },
    },

    # 催化剂体系
    'catalyst_systems': {
        'NaOH/KOH': {'type': 'strong base', 'strength': 'pKa of H2O = 15.7', 'conditions': 'aqueous or alcoholic, RT to reflux', 'selectivity': 'low (thermodynamic enolate)', 'side_reactions': 'multiple condensations, Cannizzaro (no α-H aldehydes)', 'notes': 'Simple but unselective'},
        'EtONa/NaOEt': {'type': 'alkoxide base', 'strength': 'pKa of EtOH ≈ 16', 'conditions': 'absolute EtOH, RT to reflux', 'selectivity': 'moderate', 'side_reactions': 'transesterification if esters present', 'notes': 'Classic Claisen-Schmidt conditions for crossed aldol'},
        'LDA': {'type': 'sterically hindered strong base', 'strength': 'pKa of conjugate acid ~35', 'conditions': 'THF, -78°C, anhydrous', 'selectivity': 'HIGH (kinetic enolate)', 'notes': 'Gold standard for regioselective enolate formation; kinetic vs thermodynamic control'},
        'LiHMDS': {'type': 'sterically hindered base', 'strength': 'pKa ~26', 'conditions': 'THF, -78°C', 'selectivity': 'high (kinetic)', 'notes': 'Alternative to LDA, less pyrophoric'},
        'NaH': {'type': 'strong base', 'strength': 'pKa of H2 ~35', 'conditions': 'THF/DMF, 0°C to RT', 'selectivity': 'thermodynamic enolate', 'notes': 'Forms enolate quantitatively; watch for H2 evolution'},
        'KHMDS': {'type': 'sterically hindered base', 'strength': 'pKa ~26', 'conditions': 'THF, -78°C to RT', 'selectivity': 'kinetic or thermodynamic depending on T', 'notes': 'Potassium counterion gives more reactive "naked" enolate'},
        'pyrrolidine/AcOH': {'type': 'secondary amine catalysis', 'strength': 'mild (organocatalytic)', 'conditions': 'RT, organocatalytic', 'selectivity': 'can be HIGH with chiral catalysts', 'notes': 'Proline-catalyzed asymmetric aldol is landmark reaction (List-Houk-List, 2000)'},
        'proline': {'type': 'chiral organocatalyst', 'strength': 'mild', 'conditions': 'DMSO, RT or 4°C', 'selectivity': 'excellent enantioselectivity possible', 'notes': 'Pioneering organocatalytic aldol (List, 2000); up to >99% ee reported'},
        'acid (H2SO4, TsOH)': {'type': 'Brønsted acid', 'strength': 'varies', 'conditions': 'heat, azeotropic water removal', 'mechanism': 'enol (not enolate) pathway', 'notes': 'Acid-catalyzed aldol: useful when base-sensitive groups present'},
    },

    # 溶剂效应
    'solvent_effects': {
        'protic (H2O, EtOH)': {'enolate_stability': 'solvated (less reactive)', 'rate': 'slower', 'selectivity': 'lower', 'notes': 'Hydrogen bonding stabilizes enolate, reduces reactivity and selectivity'},
        'aprotic_polar (DMF, DMSO)': {'enolate_stability': '"naked" (more reactive)', 'rate': 'faster', 'selectivity': 'higher', 'notes': 'No H-bonding to enolate; enhanced reactivity and selectivity'},
        'aprotic_nonpolar (THF, Et2O)': {'enolate_stability': 'moderate', 'rate': 'moderate', 'selectivity': 'good', 'notes': 'Standard choice for LDA-mediated aldols'},
        'CH2Cl2': {'useful_for': 'acid-catalyzed or Lewis acid promoted', 'notes': 'Inert solvent for many conditions'},
    },

    # 常见副反应
    'side_reactions': [
        ('Cannizzaro reaction', 'Strong base + aldehydes without α-H (e.g., formaldehyde, benzaldehyde) → alcohol + carboxylate', 'Avoid by using mild base or ensuring all aldehydes have α-H'),
        ('Multiple condensation', 'Product still has α-H → can undergo further aldol reactions', 'Control stoichiometry, use kinetic enolate, low T'),
        ('Dehydration to conjugated system', 'β-hydroxy carbonyl eliminates water under basic/thermal conditions', 'Desired in aldol condensation; control with T and time'),
        ('E/Z isomer mixture', 'Both E and Z enolates form → syn and anti diastereomers', 'Use chelating metal ions or bulky bases for stereocontrol'),
        ('Racemization', 'If chiral center adjacent to carbonyl, enolization causes racemization', 'Use mild conditions or alternative disconnection'),
        ('Transesterification', 'Ester-containing substrates with alkoxide bases', 'Use non-nucleophilic base like LDA'),
        ('Over-addition', 'Enolate adds to product carbonyl', 'Control stoichiometry and addition rate'),
    ],

    # 脱水条件
    'dehydration_conditions': {
        'thermal': {'condition': 'Heat (50-100°C) during or after aldol', 'mechanism': 'E1cb or E1', 'driving_force': 'conjugation stabilization (~15-20 kJ/mol)'},
        'acid_catalyzed': {'condition': 'Dilute H2SO4 or p-TsOH, heat', 'mechanism': 'E1 (via protonated OH)', 'efficiency': 'high'},
        'base_promoted': {'condition': 'Excess base, heat', 'mechanism': 'E1cb', 'efficiency': 'good for some systems'},
    },
}


@ChemMCPManager.register_tool
class AldolReaction(BaseTool):
    __version__ = "0.1.0"
    name = "AldolReaction"
    func_name = 'analyze_aldol_reaction'
    description = "Comprehensive aldol condensation reaction analysis: substrate scope (aldehydes, ketones), optimal conditions (catalyst, solvent, temperature), product prediction (β-hydroxy carbonyl / α,β-unsaturated carbonyl), limitations, side reactions, yield expectations, and optimization strategies."
    implementation_description = "Uses extensive knowledge base covering substrate types (aliphatic/aromatic aldehydes, ketones, cyclic ketones, β-dicarbonyls), catalyst systems (hydroxide, alkoxide, LDA, organocatalysts), solvent effects, dehydration pathways, side reactions (Cannizzaro, multiple condensation), and stereochemical models (Zimmerman-Traxler, Evans). Provides practical recommendations for laboratory execution."
    categories = ["Reaction"]
    tags = ["Aldol", "Condensation", "Carbonyl", "Enolate", "C-C Bond Formation", "Stereochemistry"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("substrate1_smiles", "str", "N/A", "SMILES or name of substrate 1 (enolizable carbonyl compound)."),
        ("substrate2_smiles", "str", "", "SMILES or name of substrate 2 (electrophilic carbonyl). Leave empty for self-aldol."),
        ("catalyst_type", "str", "base", "Catalyst type: 'base' (e.g., NaOH, LDA), 'acid' (e.g., H2SO4), 'organocatalytic' (e.g., proline), or specific name."),
        ("solvent", "str", "EtOH", "Solvent for the reaction."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: substrate1 [substrate2] [catalyst] [solvent]. E.g., 'CC=O Cc1ccccc1=O base EtOH'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing reaction_type, substrate_analysis, optimal_conditions, product_prediction, scope_limitations, yield_expectation, side_reactions, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"substrate1_smiles": "CC=O", "substrate2_smiles": "", "catalyst_type": "base", "solvent": "EtOH"},
            "text_input": {"query": "CC=O base EtOH"},
            "output": {"result": {
                "reaction_type": "Self-aldol of acetaldehyde",
                "substrate1": {"name": "acetaldehyde", "type": "aliphatic_aldehyde", "has_alpha_H": True, "alpha_h_pKa": 17, "role": "both enolate source AND electrophile"},
                "product_aldol": "3-hydroxybutanal (aldol adduct)",
                "product_condensation": "crotonaldehyde (but-2-enal, after dehydration)",
                "optimal_conditions": {"catalyst": "10% NaOH aq. or NaOEt/EtOH", "temperature": "0°C → RT then heat for dehydration", "time": "1-4 h (addition), 1-2 h (dehydration)"},
                "yield_expectation": "60-80% (condensation product)",
                "major_limitation": "Multiple condensations possible (product still has α-H)",
                "recommendations": ["Control stoichiometry", "Consider crossed aldol with aromatic aldehyde for cleaner product"],
            }},
        },
        {
            "code_input": {"substrate1_smiles": "CC(=O)C", "substrate2_smiles": "Cc1ccccc1=O", "catalyst_type": "base", "solvent": "EtOH"},
            "text_input": {"query": "CC(=O)C Cc1ccccc1=O base EtOH"},
            "output": {"result": {
                "reaction_type": "Crossed aldol (Claisen-Schmidt): acetone + 4-methylbenzaldehyde",
                "substrate1": {"name": "acetone", "type": "aliphatic_ketone", "role": "enolate source (nucleophile)"},
                "substrate2": {"name": "4-methylbenzaldehyde", "type": "aromatic_aldehyde", "role": "electrophile (no α-H, cannot enolize)"},
                "product_condensation": "4-methylbenzylideneacetone (mesityl oxide analog)",
                "yield_expectation": "70-85%",
                "clean_crossed": True,  # Only one component can enolize
                "optimal_conditions": {"catalyst": "10% NaOH aq./EtOH", "T": "0-5°C → RT", "stoichiometry": "acetone (2 eq) : aldehyde (1 eq)"},
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_ALDOL_DATA)

    def _run_base(self, substrate1_smiles: str, substrate2_smiles: str = "", catalyst_type: str = "base", solvent: str = "EtOH") -> dict:
        if not substrate1_smiles:
            raise ChemMCPInputError("Substrate 1 is required.")

        # Classify substrates
        sub1 = self._classify_substrate(substrate1_smiles)
        sub2 = self._classify_substrate(substrate2_smiles) if substrate2_smiles else None

        # Determine reaction type
        is_crossed = sub2 is not None
        rxn_type = self._determine_reaction_type(sub1, sub2, is_crossed)

        # Analyze catalyst
        cat_info = self._analyze_catalyst(catalyst_type)

        # Predict products
        products = self._predict_products(sub1, sub2, is_crossed)

        # Optimal conditions
        conditions = self._optimize_conditions(sub1, sub2, catalyst_type, solvent, is_crossed)

        # Scope & limitations
        scope_lim = self._analyze_scope_limitations(sub1, sub2, is_crossed)

        # Yield expectation
        yield_exp = self._estimate_yield(sub1, sub2, catalyst_type, is_crossed)

        # Side reactions
        side_rxns = self._predict_side_reactions(sub1, sub2, catalyst_type, is_crossed)

        result = {
            "result": {
                "reaction_type": rxn_type,
                "substrate1_analysis": sub1,
                "substrate2_analysis": sub2,
                "catalyst_analysis": cat_info,
                "product_prediction": products,
                "optimal_conditions": conditions,
                "scope_and_limitations": scope_lim,
                "yield_expectation": yield_exp,
                "predicted_side_reactions": side_rxns,
                "recommendations": self._generate_recommendations(sub1, sub2, catalyst_type, is_crossed),
                "stereochemical_notes": self._stereochemistry_notes(sub1, sub2, catalyst_type),
                "summary": self._build_summary(rxn_type, products, yield_exp, is_crossed),
            }
        }
        logger.info(f"AldolReaction: {rxn_type} analyzed")
        return result

    def _classify_substrate(self, smi_or_name):
        s = (smi_or_name or "").strip().lower()
        data = self.data['substrate_types']

        # Pattern matching
        patterns = {
            'aliphatic_aldehyde': [r'^c(c=o)$', r'acetaldehyde', r'propion', r'butyr', r'ch3cho', r'valeraldehyde'],
            'aromatic_aldehyde': [r'c1.*ccc.*c1.*c=O', r'benzaldehyde', r'tolualdehyde', r'anisaldehyde', r'phcho'],
            'aliphatic_ketone': [r'cc(=O)c', r'acetone', r'mek', r'butanone', r'(ch3)2co'],
            'cyclic_ketone': [r'C1CCCCC1=O', r'cyclohexanone', r'cyclopentanone', r'cycloheptanone'],
            'beta_ketoester': [r'acetoacetate', r'ethyl.acetoacetate', r'beta.ketoester'],
            'beta_diketone': [r'acetylacetone', r'pentanedione', r'beta.diketone'],
        }

        for stype, pats in patterns.items():
            for pat in pats:
                if re.search(pat, s, re.IGNORECASE):
                    info = dict(data.get(stype, {}))
                    info['classified_as'] = stype
                    info['input'] = smi_or_name
                    return info

        return {
            'classified_as': 'unknown',
            'input': smi_or_name,
            'has_alpha_H': True,  # assume yes unless known otherwise
            'note': f"Could not classify '{smi_or_name}' precisely. Treat as generic carbonyl.",
        }

    def _determine_reaction_type(self, s1, s2, is_crossed):
        if not is_crossed:
            name = s1.get('classified_as', '?').replace('_', ' ')
            return f"Self-aldol of {s1.get('input', name)}"
        n1 = s1.get('input', 'substrate 1')
        n2 = s2.get('input', 'substrate 2')
        # Check if it's Claisen-Schmidt type (aromatic aldehyde + enolizable partner)
        if s2 and s2.get('classified_as') == 'aromatic_aldehyde':
            return f"Crossed aldol (Claisen-Schmidt): {n1} + {n2}"
        if s1.get('classified_as') == 'aromatic_aldehyde':
            return f"Crossed aldol (Claisen-Schmidt): {n2} + {n1}"
        return f"Crossed aldol: {n1} + {n2}"

    def _analyze_catalyst(self, cat):
        cat_lower = (cat or "base").lower().strip()
        systems = self.data['catalyst_systems']

        for key, info in systems.items():
            if key.lower() in cat_lower or cat_lower in key.lower():
                return {"selected": key, **info}

        if cat_lower.startswith('bas'):
            return {"selected": "generic base", "type": "base", "notes": f"Base-type catalyst specified: {cat}. See NaOH/KOH or NaOEt for typical conditions."}
        elif cat_lower.startswith('aci'):
            return {"selected": "generic acid", "type": "Brønsted acid", "mechanism": "enol pathway", "notes": f"Acid-catalyzed aldol via enol intermediate."}
        return {"selected": cat, "type": "unspecified", "notes": f"Catalyst '{cat}' not in standard database. Provide specifics for detailed analysis."}

    def _predict_products(self, s1, s2, is_crossed):
        if not is_crossed:
            inp = s1.get('input', '?')
            return {
                "aldol_adduct": f"β-hydroxy carbonyl dimer from {inp}",
                "condensation_product": f"α,β-unsaturated carbonyl from {inp}",
                "functional_groups_adduct": ["carbonyl", "secondary alcohol"],
                "functional_groups_condensed": ["carbonyl", "conjugated alkene"],
            }

        n1 = s1.get('input', '?')
        n2 = s2.get('input', '?')
        return {
            "aldol_adduct": f"β-hydroxy carbonyl from {n1} (enolate) + {n2} (electrophile)",
            "condensation_product": f"α,β-unsaturated carbonyl (cross-condensation of {n1} + {n2})",
            "crossed_cleanliness": "clean" if (s1.get('has_alpha_H', True) != s2.get('has_alpha_H', True)) else "potentially messy (both can enolize)",
        }

    def _optimize_conditions(self, s1, s2, cat, solv, is_crossed):
        cat_lower = (cat or "").lower()
        base = "10% NaOH(aq)/EtOH" if 'naoh' in cat_lower or cat_lower == 'base' else cat

        cond = {
            "catalyst_loading": "10-20 mol%" if 'organocatalytic' in cat_lower else "1-2 equiv (stoichiometric)" if 'lda' in cat_lower or 'nah' in cat_lower else "10-20 mol%",
            "solvent": solv or "EtOH/H2O or THF (anhydrous for LDA)",
            "temperature": "-78°C → RT (for LDA)" if 'lda' in cat_lower else "0°C → RT (addition), then reflux (dehydration)",
            "atmosphere": "N2 or Ar (especially for LDA, Grignard-type bases)",
            "workup": "Quench with sat. NH4Cl or dilute acid; extract with EtOAc",
            "purification": "Column chromatography or distillation (if volatile)",
        }

        if is_crossed and s2 and s2.get('classified_as') == 'aromatic_aldehyde':
            cond.update({
                "stoichiometry": f"{s1.get('input','?')} (1.5-2 eq) as enolate source, {s2.get('input','?')} (1 eq) as electrophile",
                "addition_order": f"Add {s2.get('input','?')} slowly to pre-formed enolate of {s1.get('input','?')}",
                "temperature": "0-5°C (minimize side reactions)",
            })

        return cond

    def _analyze_scope_limitations(self, s1, s2, is_crossed):
        lim = []
        scope = []

        if is_crossed:
            if s1.get('has_alpha_H') and s2 and s2.get('has_alpha_H'):
                lim.append("Both substrates have α-H → 4 possible aldol products (self + crossed)")
                lim.append("Solution: Use LDA for kinetic enolate of one component only")
                scope.append("Works well with one enolizable + one non-enolizable partner (Claisen-Schmidt)")
            elif s2 and not s2.get('has_alpha_H'):
                scope.append("Clean crossed aldol: only one enolate direction possible")
                scope.append("Aromatic aldehydes are excellent electrophiles (no α-H, highly electrophilic)")

        s1_type = s1.get('classified_as', '')
        if s1_type == 'aliphatic_aldehyde':
            lim.append("Aliphatic aldehydes tend toward multiple condensations (product retains α-H)")
            scope.append("Works well with aromatic aldehydes → cinnamaldehyde derivatives")

        if s1_type == 'aliphatic_ketone':
            scope.append("Ketones are less reactive electrophiles → better as enolate sources")
            lim.append("Ketone self-aldol is slower than aldehyde (less reactive carbonyl)")

        return {"scope": scope or ["Standard aldol chemistry applies"], "limitations": lim or ["Standard precautions apply"]}

    def _estimate_yield(self, s1, s2, cat, is_crossed):
        score = 50  # baseline

        if is_crossed:
            if s2 and not s2.get('has_alpha_H'): score += 30  # clean crossed
            elif s1.get('has_alpha_H') and s2 and s2.get('has_alpha_H'): score -= 15  # messy
        else:
            if s1.get('classified_as') == 'aliphatic_aldehyde': score -= 10  # multiple condensation
            if s1.get('classified_as') == 'cyclic_ketone': score += 10  # ring formation favorable

        if 'lda' in cat.lower() or 'proline' in cat.lower(): score += 15  # selective
        if 'organocatalytic' in cat.lower(): score += 10
        if 'naoh' in cat.lower() or cat == 'base': score -= 5  # unselective

        score = max(10, min(95, score))
        return f"{score-5}-{score+5}% (estimated)"

    def _predict_side_reactions(self, s1, s2, cat, is_crossed):
        relevant = []
        for name, desc, prevention in self.data['side_reactions']:
            include = False
            if 'Cannizzaro' in name:
                if (s1.get('classified_as') in ('aromatic_aldehyde',) or
                    (s2 and s2.get('classified_as') in ('aromatic_aldehyde',))):
                    include = True
            if 'Multiple' in name:
                if not is_crossed or (s1.get('has_alpha_H') and s2 and s2.get('has_alpha_H')):
                    include = True
            if 'Dehydration' in name:
                include = True  # almost always relevant
            if 'Racemization' in name:
                include = True  # general risk
            if include:
                relevant.append({"reaction": name, "description": desc, "prevention": prevention})
        return relevant or [{"reaction": "None major expected under optimized conditions"}]

    def _stereochemistry_notes(self, s1, s2, cat):
        notes = [
            "New stereocenters: up to 2 (α and β positions of new bond)",
            "Zimmerman-Traxler model: 6-membered cyclic TS determines syn/anti diastereoselectivity",
            "E-enolate → anti aldol; Z-enolate → syn aldol (Zimmerman-Traxler)",
        ]
        if 'proline' in cat.lower() or 'organocatalytic' in cat.lower():
            notes.append("Organocatalytic (proline): Enamine mechanism gives high enantioselectivity (up to >99% ee)")
            notes.append("Anti-selectivity typical for proline-catalyzed aldols")
        if 'lda' in cat.lower():
            notes.append("LDA at -78°C: Kinetic (E)-enolate favored → anti product (via chair TS)")
        return notes

    def _generate_recommendations(self, s1, s2, cat, is_crossed):
        recs = []
        if is_crossed:
            if s2 and not s2.get('has_alpha_H'):
                recs.append("✅ Good candidate for clean crossed aldol — use enolizable component as nucleophile")
                recs.append("Add electrophilic aldehyde slowly to pre-formed enolate")
            elif s1.get('has_alpha_H') and s2 and s2.get('has_alpha_H'):
                recs.append("⚠️ Both components enolizable — consider LDA for kinetic enolate control")
                recs.append("Or use a pre-formed en equivalent (silyl enol ether, enamine)")
        else:
            recs.append("Self-aldol: control stoichiometry to minimize multiple additions")
            if s1.get('classified_as') == 'aliphatic_aldehyde':
                recs.append("Consider using aromatic aldehyde for cleaner crossed variant")

        if cat == 'base' or 'naoh' in cat.lower():
            recs.append("Upgrade to LDA/NaOEt for improved selectivity if needed")

        return recs

    def _build_summary(self, rxn_type, products, yield_exp, is_crossed):
        prod = products.get('condensation_product', products.get('aldol_adduct', '?'))
        return f"{rxn_type}. Product: {prod}. Expected yield: {yield_exp}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        s1 = parts[0] if parts else ""
        s2 = parts[1] if len(parts) > 1 else ""
        cat = parts[2] if len(parts) > 2 else "base"
        solv = parts[3] if len(parts) > 3 else "EtOH"
        return self._run_base(s1, s2, cat, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
