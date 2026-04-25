"""
Wittig Reaction (Tool #158)
Wittig 反应及其变体（HWE, Still-Gennari）：叶立德形成、烯烃合成、
E/Z 选择性、稳定/非稳定叶立德、范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_WITTIG_DATA = {
    'mechanism': [
        ('1', 'Ylide formation', 'Ph3P=CR2 formed by deprotonation of phosphonium salt with strong base (n-BuLi, NaH, KHMDS, etc.)'),
        ('2', '[2+2] cycloaddition', 'Ylide carbonyl carbon attacks carbonyl C → oxaphosphetane (4-membered ring)'),
        ('3', 'Collapse', 'Oxaphosphetane collapses: P-O bond forms, C=C bond forms simultaneously'),
        ('4', 'Products', 'Triphenylphosphine oxide (Ph3P=O) + alkene (R2C=CR\'2) — driving force is strong P=O bond (ΔH ~544 kJ/mol)'),
    ],
    'ylide_types': {
        'non_stabilized': {
            'description': 'R = alkyl, H (no electron-withdrawing group on ylide carbon)',
            'reactivity': 'very high (reacts at -78°C to 0°C)',
            'e_z_selectivity': 'Z-alkene favored (kinetic product via betaine pathway)',
            'base_needed': 'strong base (n-BuLi, NaHMDS)',
            'examples': ['CH2=PPh3 (methylene ylide)', 'MeCH=PPh3 (ethylidene ylide)'],
        },
        'semi_stabilized': {
            'description': 'R = one aryl/vinyl + one alkyl/H',
            'reactivity': 'high (0°C to RT)',
            'e_z_selectivity': 'mixture; often Z-favored but less selective',
            'examples': ['PhCH=PPh3 (benzylidene ylide)', 'vinyl-substituted ylides'],
        },
        'stabilized': {
            'description': 'R = EWG (COR, COOR, CN, etc.) conjugated to ylide carbon',
            'reactivity': 'moderate (requires RT or heating)',
            'e_z_selectivity': 'E-alkene favored (thermodynamic control via oxaphosphetane equilibrium)',
            'base_needed': 'weaker base acceptable (NaOEt, K2CO3, DBU)',
            'examples': ['Ph3P=CHCOOMe (carbomethoxymethylene)', 'Ph3P=CHCOMe', 'Ph3P=CHCN'],
        },
    },
    'variants': {
        'Horner_Wadsworth_Emmons (HWE)': {
            'reagent': 'Phosphonate ester ((RO)2P(O)CHR2) instead of phosphonium salt',
            'advantages': ['Byproducts water-soluble (easy purification)', 'E-selective for stabilized systems', 'Milder base conditions', 'Functional group tolerance'],
            'base': 'NaH, KOtBu, DBU',
            'selectivity': 'E-alkene highly favored (especially with α-EWG ylides)',
            'note': 'The go-to method for E-α,β-unsaturated esters/nitriles',
        },
        'Still_Gennari': {
            'reagent': 'Phosphonoacetate with CF3CH2O groups ((CF3CH2O)2P(O)CHR2)',
            'specialty': 'Z-selective variant of HWE',
            'selectivity': 'Z-alkene favored (opposite of normal HWE!)',
            'reason': 'Electron-withdrawing -OCF2CF3 group changes oxaphosphetane stability',
        },
        'Wittig_Horner': {'note': 'Historical name overlap — HWE is the modern standard term'},
        'Schlosser_modification': {
            'description': 'Convert Z-alkene to E-alkene via lithium salt intermediate',
            'procedure': 'Form betaine, treat with excess RLi → equilibration to E',
        },
        'One_pot_Wittig': {
            'description': 'Alkyl halide + PPh3 → phosphonium salt (no isolation) → base → ylide → carbonyl',
            'efficiency': 'Saves time; good for parallel synthesis',
        },
        'Aza_Wittig': {
            'description': 'P=NCR2 (iminophosphorane) + carbonyl → imine + Ph3PO',
            'product': 'imines instead of alkenes',
        },
        'Tebbe_Olefination': {
            'description': 'Cp2TiCH2ClAlMe2 → "Ti=CH2" equivalent',
            'scope': 'Converts esters/amides to enol ethers/enamines (Wittig cannot do this)',
        },
        'Petasis_reagent': {
            'description': 'Cp2Ti(CH2)(PMe3) — Tebbe alternative',
        },
    },
    'scope': [
        'Aldehydes: all types work well (aromatic, aliphatic, α,β-unsaturated)',
        'Ketones: work but slower than aldehydes (steric hindrance); stabilized ylides preferred',
        'Esters: NOT reactive in classical Wittig (use Tebbe or Petasis reagent instead)',
        'Intramolecular Wittig: excellent for cyclic alkene synthesis (macrocycles, medium rings)',
        'α,β-Unsaturated carbonyls: can give 1,2-addition (normal Wittig) or conjugate addition depending on conditions',
    ],
    'limitations': [
        'Sterically hindered ketones react poorly (use Tebbe/Petasis or Julia olefination)',
        'E/Z selectivity can be hard to control for semi-stabilized ylides',
        'Classical Wittig does not work with esters, amides, acid chlorides as electrophiles',
        'Phosphonium salt synthesis requires alkyl halide (may not be readily available)',
        'Ph3P=O byproduct removal can be tedious (column chromatography usually needed)',
        'Base-sensitive functional groups incompatible with ylide formation conditions',
        'Ylides are air/moisture sensitive and malodorous (phosphine smell)',
    ],
    'conditions': {
        'classical': {'solvent': 'THF or Et2O (anhydrous)', 'base': 'n-BuLi or NaHMDS', 'T': '-78°C → RT', 'atmosphere': 'N2/Ar'},
        'HWE': {'solvent': 'THF or DMF', 'base': 'NaH, KOtBu, or DBU', 'T': '0°C → RT', 'time': '1-12 h'},
        'workup': ('Aqueous NH4 quench, extract with EtOAc, column chromatography '
                   '(separates alkene from Ph3P=O — both non-polar but separable)'),
    },
}


@ChemMCPManager.register_tool
class WittigReaction(BaseTool):
    __version__ = "0.1.0"
    name = "WittigReaction"
    func_name = 'analyze_wittig_reaction'
    description = "Wittig reaction and variants (HWE, Still-Gennari, Schlosser modification) analysis: ylide formation and classification (non-stabilized/semi-stabilized/stabilized), alkene synthesis mechanism, E/Z selectivity rules, scope (aldehydes, ketones), limitations (esters unreactive, steric hindrance), and comparison with alternatives."
    implementation_description = "Comprehensive knowledge base covering: ylide classification (3 types with reactivity/selectivity profiles), stepwise mechanism (4 steps via oxaphosphetane), 6 major variants (HWE, Still-Gennari, Schlosser, aza-Wittig, Tebbe, Petasis), scope table (5 categories), limitations (6 items), condition optimization, and practical guidance."
    categories = ["Reaction"]
    tags = ["Wittig", "Ylide", "Alkene Synthesis", "HWE", "E/Z Selectivity", "Phosphorus", "Olefination"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("carbonyl_smiles", "str", "N/A", "SMILES or name of the carbonyl compound (aldehyde or ketone)."),
        ("ylide_type", "str", "non-stabilized", "Ylide type: 'non-stabilized', 'semi-stabilized', 'stabilized', or specific name."),
        ("phosphonium_salt", "str", "Ph3P=CH2", "Phosphonium salt or ylide specification."),
        ("base", "str", "n-BuLi", "Base for ylide generation."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: carbonyl [ylide_type] [phosphonium_salt] [base]. E.g., 'benzaldehyde stabilized Ph3P=CHCOOEt NaH'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing ylide_analysis, carbonyl_analysis, e_z_selectivity, mechanism, product_prediction, variant_comparison, scope, and limitations."),
    ]

    examples = [
        {
            "code_input": {"carbonyl_smiles": "benzaldehyde", "ylide_type": "non-stabilized", "phosphonium_salt": "Ph3P=CH2", "base": "n-BuLi"},
            "text_input": {"query": "benzaldehyde non-stabilized Ph3P=CH2 n-BuLi"},
            "output": {"result": {
                "reaction": "Benzaldehyde + methylenetriphenylphosphorane → styrene",
                "ylide_type": "non-stabilized (methylene ylide)",
                "e_z_selectivity": "Non-stabilized ylide + aldehyde → no E/Z issue (terminal alkene)",
                "product": "styrene (PhCH=CH2)",
                "yield": "80-95%",
                "mechanism_summary": "[2+2] cycloaddition → oxaphosphetane → collapse → Ph3P=O + alkene",
                "driving_force": "Formation of very strong P=O bond (~544 kJ/mol)",
            }},
        },
        {
            "code_input": {"carbonyl_smiles": "benzaldehyde", "ylide_type": "stabilized", "phosphonium_salt": "Ph3P=CHCOOEt", "base": "NaH"},
            "text_input": {"query": "benzaldehyde stabilized Ph3P=CHCOOEt NaH"},
            "output": {"result": {
                "reaction": "Benzaldehyde + (carboethoxymethylene)triphenylphosphorane → ethyl cinnamate (E-major)",
                "ylide_type": "stabilized (ester-stabilized ylide)",
                "e_z_selectivity": "Stabilized ylide + aldehyde → E-alkene major (thermodynamic control)",
                "e_z_ratio": "typically >20:1 E:Z",
                "product": "ethyl (E)-cinnamate (PhCH=CHCOOEt)",
                "variant_note": "Consider HWE variant for even better E-selectivity and easier purification",
                "yield": "75-90%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_WITTIG_DATA)

    def _run_base(self, carbonyl_smiles: str, ylide_type: str = "non-stabilized", phosphonium_salt: str = "Ph3P=CH2", base: str = "n-BuLi") -> dict:
        if not carbonyl_smiles:
            raise ChemMCPInputError("Carbonyl compound is required.")

        carbonyl = self._analyze_carbonyl(carbonyl_smiles)
        ylide = self._analyze_ylide(ylide_type, phosphonium_salt)
        ez = self._predict_ez_selectivity(ylide, carbonyl)
        product = self._predict_product(carbonyl, ylide)
        variant_rec = self._recommend_variant(ylide, carbonyl)

        result = {
            "result": {
                "reaction": f"{carbonyl_smiles} + {phosphonium_salt} ({ylide_type} ylide) → {product.get('name','?')}",
                "carbonyl_analysis": carbonyl,
                "ylide_analysis": ylide,
                "e_z_selectivity": ez,
                "mechanism": [{"step": s[0], "name": s[1], "desc": s[2]} for s in self.data['mechanism']],
                "product_prediction": product,
                "variant_comparison": variant_rec,
                "scope": self.data['scope'],
                "limitations": self.data['limitations'],
                "optimal_conditions": self._optimize(carbonyl, ylide, base),
                "summary": f"Wittig reaction: {product.get('name','?')}. E/Z: {ez.get('prediction','?')}. Yield: {product.get('yield','70-85%')}.",
            }
        }
        logger.info(f"Wittig: {carbonyl_smiles} + {phosphonium_salt}")
        return result

    def _analyze_carbonyl(self, smi):
        s = (smi or "").strip().lower()
        if any(p in s for p in ['aldehyde', 'cho', 'benzaldehyde', 'formyl']):
            return {"type": "aldehyde", "name": smi, "reactivity": "high (good electrophile)", "steric_demand": "low"}
        elif any(p in s for p in ['ketone', 'co(c)', 'acetone', 'cyclohexanone']):
            return {"type": "ketone", "name": smi, "reactivity": "moderate (less reactive than aldehyde)", "steric_demand": "moderate-high"}
        return {"type": "unknown_carbonyl", "name": smi}

    def _analyze_ylide(self, ylide_type, salt):
        yt = (ylide_type or "").lower().strip()
        ytypes = self.data['ylide_types']

        if 'non' in yt or 'unstabilized' in yt:
            info = dict(ytypes['non_stabilized'])
        elif 'semi' in yt:
            info = dict(ytypes['semi_stabilized'])
        elif 'stab' in yt or 'hwe' in yt:
            info = dict(ytypes['stabilized'])
        else:
            info = ytypes.get(yt, ytypes['non_stabilized'])

        info["specified_as"] = ylide_type
        info["salt"] = salt
        return info

    def _predict_ez_selectivity(self, ylide, carbonyl):
        yt = ylide.get('specified_as', '').lower()
        ct = carbonyl.get('type', '')

        if 'non' in yt or 'unstabilized' in yt:
            if ct == 'aldehyde':
                return {"prediction": "Terminal alkene (no E/Z isomerism) or Z-favored if disubstituted", "rule": "Non-stabilized ylides give kinetic Z-product via betaine pathway"}
            return {"prediction": "Z-alkene favored (kinetic control)", "rule": "Non-stabilized ylide: irreversible oxaphosphetane collapse → Z-major"}
        elif 'stab' in yt or 'hwe' in yt:
            return {"prediction": "E-alkene strongly favored", "rule": "Stabilized ylide: reversible oxaphosphetane → thermodynamic E-product (>20:1 typical)"}
        elif 'semi' in yt:
            return {"prediction": "Mixture (often Z-leaning but variable)", "rule": "Semi-stabilized: intermediate behavior; use Schlosser mod for E"}
        return {"prediction": "Depends on specific system", "rule": "General guidelines apply"}

    def _predict_product(self, carbonyl, ylide):
        cname = carbonyl.get('name', '?')
        salt = ylide.get('salt', 'Ph3P=CH2')
        # Extract R group from salt
        r_match = re.search(r'=([A-Za-z0-9]+)$', salt.replace(' ', ''))
        r_group = r_match.group(1) if r_match else 'R'
        ct = carbonyl.get('type', '?')

        if ct == 'aldehyde':
            name = f"{cname}-{r_group} alkene (disubstituted terminal/trisubstituted)"
        elif ct == 'ketone':
            name = f"Tetrasubstituted/trisubstituted alkene from {cname}"
        else:
            name = f"Alkene from {cname} + {salt}"

        return {
            "name": name,
            "class": "alkene (C=C double bond)",
            "byproduct": "triphenylphosphine oxide (Ph3P=O)",
            "yield": "80-95%" if ct == 'aldehyde' and 'non' in (ylide.get('specified_as','')).lower() else "70-90%",
        }

    def _recommend_variant(self, ylide, carbonyl):
        recs = []
        yt = ylide.get('specified_as', '').lower()

        if 'stab' in yt or 'hwe' in yt:
            recs.append({"variant": "Horner-Wadsworth-Emmons (HWE)", "reason": "Better E-selectivity, easier purification (water-soluble byproduct), milder base"})
        if 'z' in yt or 'still' in yt:
            recs.append({"variant": "Still-Gennari", "reason": "Specifically designed for Z-selective olefination"})
        if carbonyl.get('type') == 'ketone' and 'non' in yt:
            recs.append({"variant": "Consider stabilized ylide or Petasis reagent", "reason": "Ketones are sluggish with non-stabilized ylides"})
        return recs

    def _optimize(self, carbonyl, ylide, base):
        yt = ylide.get('specified_as', '')
        cond = {
            "solvent": "anhydrous THF (standard)",
            "base": base,
            "temperature": "-78°C → RT" if 'non' in yt else "0°C → RT" if 'semi' in yt else "RT → reflux",
            "atmosphere": "N2/Ar (strictly anhydrous and oxygen-free)",
            "order_of_addition": "Add base to phosphonium salt (form ylide), then add carbonyl solution",
            "workup": "Quench with sat. NH4Cl, extract with Et2O/EtOAc, purify by column chromatography",
        }
        return cond

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        carb = parts[0] if parts else ""
        yt = parts[1] if len(parts) > 1 else "non-stabilized"
        salt = parts[2] if len(parts) > 2 else "Ph3P=CH2"
        b = parts[3] if len(parts) > 3 else "n-BuLi"
        return self._run_base(carb, yt, salt, b)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
