"""
Buchwald-Hartwig Amination (Tool #163)
Buchwald-Hartwig 胺化反应：钯催化芳基卤化物与胺的C-N交叉偶联反应。
涵盖：催化循环（氧化加成、胺配位/去质子化、还原消除）、
配体重要性（BINAP、XPhos等富电子大位阻膦配体）、碱效应、
底物范围（伯/仲胺、芳胺、酰胺、N-杂环）和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_BUCHWALD_HARTWIG_DATA = {
    'catalytic_cycle': [
        ('1', 'Oxidative addition', 'Pd(0) inserts into Ar-X bond → Pd(II)(Ar)(X). Ligand controls rate and selectivity. Ar-I > Ar-Br > Ar-Cl >> Ar-OTf.'),
        ('2', 'Amine coordination & deprotonation', 'Amine coordinates to Pd(II); base deprotonates → Pd(II)(Ar)(NR²R³) complex. This step is often RATE-DETERMINING.'),
        ('3', 'Reductive elimination', 'Pd(II)(Ar)(NR²R³) eliminates Ar-NR²R³ (C-N bond formed!) + regenerates Pd(0). Sterics/electronics of ligand control this step.'),
    ],
    'ligand_importance': {
        'why_ligands_matter': (
            'Ligands are CRITICAL in B-H amination — they control: '
            '(1) oxidative addition rate (electron-rich ligands accelerate), '
            '(2) amine coordination/deprotonation (bulky ligands create open site), '
            '(3) reductive elimination (bulky ligands promote C-N bond formation), '
            '(4) prevent Pd(0) aggregation into inactive Pd black'
        ),
        'ligand_categories': {
            'bidentate_phosphines': {
                'examples': ['BINAP', 'DavePhos', 'Xantphos', 'DPPF'],
                'characteristics': 'Chelating effect stabilizes Pd; good for standard substrates',
                'limitations': 'Less effective for challenging aryl chlorides',
            },
            'bulky_biaryl_phosphines_Buchwald_type': {
                'examples': ['XPhos', 'BrettPhos', 'RuPhos', 'DavePhos', 'SPhos', 'JohnPhos'],
                'characteristics': 'STATE OF THE ART: electron-rich + very bulky; enable aryl chloride couplings; high activity at low loading',
                'note': 'Revolutionized B-H amination since ~2000',
            },
            'N_heterocyclic_carbenes_NHC': {
                'examples': ['IPr', 'IMes', 'PEPPSI series', 'SIPr'],
                'characteristics': 'Very strong σ-donors; highly active; air-stable complexes',
            },
            'dialkylbiaryl_phosphines': {
                'examples': ['CyJohnPhos', 'tBuBrettPhos', 'tBuXPhos'],
                'characteristics': 'Even bulkier than Buchwald ligands; for extremely hindered substrates',
            },
        },
    },
    'catalyst_systems': {
        'Pd2(dba)3/BINAP': {'type': 'Classic system', 'loading': 'Pd 1-3 mol%, BINAP 2-4 mol%', 'scope': 'Historically important; good for aryl bromides/iodides'},
        'Pd(OAc)2/XPhos': {'type': 'Modern Buchwald ligand', 'loading': 'Pd 0.5-1 mol%, XPhos 1-2 mol%', 'scope': 'Excellent for aryl chlorides, bromides; broad scope'},
        'Pd(OAc)2/BrettPhos': {'type': 'Bulky Buchwald ligand', 'loading': 'Pd 0.5-1 mol%, BrettPhos 1-2 mol%', 'scope': 'Best for hindered substrates; aryl chlorides'},
        'Pd(OAc)2/RuPhos': {'type': 'Buchwald ligand', 'loading': 'Pd 0.5-1 mol%, RuPhos 1-2 mol%', 'scope': 'Excellent for primary aliphatic amines'},
        'Pd2(dba)3/DavePhos': {'type': 'Buchwald ligand', 'loading': 'Pd 1-2 mol%, DavePhos 2-4 mol%', 'scope': 'Versatile; good for cyclic amines'},
        'Pd-NHC (PEPPSI-IPr)': {'type': 'NHC Pd catalyst', 'loading': 'Pd 0.5-1 mol%', 'scope': 'Air-stable; no external ligand needed'},
        'Pd(OAc)2/P(t-Bu)3': {'type': 'Bulky trialkylphosphine', 'loading': 'Pd 1-2 mol%, P(t-Bu)3 2-4 mol%', 'scope': 'Very active but air-sensitive ligand'},
        'BrettPhos precatalyst (G3)': {'type': 'Pre-formed Pd-ligand complex', 'loading': '0.5-1 mol% total', 'scope': 'Most convenient; off-the-shelf ready to use'},
    },
    'base_options': [
        ('NaOtBu (sodium tert-butoxide)', 'STRONG base; most common choice for B-H amination', 'Essential for deprotonation of coordinated amine'),
        ('K3PO4 (potassium phosphate)', 'Strong, non-nucleophilic', 'Good for base-sensitive substrates'),
        ('Cs2CO3', 'Strong, soluble in organic solvents', 'For challenging couplings'),
        ('NaOH (aqueous)', 'Strong, cheap', 'Aqueous conditions possible'),
        ('K2CO3', 'Moderate base', 'For milder conditions or sensitive substrates'),
        ('DBU', 'Strong, non-ionic organic base', 'Special applications'),
    ],
    'amine_scope': [
        ('Primary aliphatic amines (RNH2)', 'WORK WELL — give secondary arylamines', 'Use RuPhos or BrettPhos for best results'),
        ('Secondary aliphatic amines (R2NH)', 'WORK WELL — give tertiary arylamines', 'Standard conditions usually suffice'),
        ('Anilines (ArNH2)', 'WORK — give diarylamines', 'May need milder conditions; electron-poor anilines easier'),
        ('N-H heterocycles (pyrrole, indole, carbazole, etc.)', 'WORK — N-arylation of heterocycles', 'Important for pharmaceutical synthesis'),
        ('Amides (NH2COR)', 'POSSIBLE — give N-aryl amides', 'Needs optimization; weaker nucleophile'),
        ('Sulfonamides (NHSO2R)', 'WORK — give N-aryl sulfonamides', 'Acidic NH facilitates deprotonation'),
        ('Hydrazines (NHNH2)', 'WORK — give aryl hydrazines', 'Careful with over-arylation'),
        ('Ammonia (NH3)', 'CHALLENGING — gives primary arylamines', 'Requires special conditions; can over-arylate'),
        ('Amino acid esters', 'WORK — chiral amine coupling', 'No racemization under typical conditions'),
    ],
    'aryl_halide_scope': [
        ('Aryl iodides', 'Excellent reactivity', 'Standard substrate'),
        ('Aryl bromides', 'Good — most commonly used', 'Work well with modern ligands'),
        ('Aryl chlorides', 'Possible with modern ligands', 'Cheap; needs XPhos/BrettPhos-type ligands'),
        ('Aryl triflates', 'Very reactive', 'From phenols; excellent alternative'),
        ('Aryl tosylates/mesylates', 'Possible with active catalysts', 'Cheaper phenol derivatives'),
        ('Heteroaryl halides', 'Good — pyridine, pyrimidine, thiophene, etc.', 'May need optimized conditions'),
        ('Stereo-defined vinyl halides', 'Possible — gives enamines', 'Stereochemistry may be retained'),
    ],
    'solvent_systems': [
        ('Toluene or dioxane', 'Most common non-polar solvents', 'Standard choice; high boiling point'),
        ('tert-Butanol (t-BuOH)', 'Polar protic solvent', 'Excellent for many B-H reactions; dissolves inorganic bases well'),
        ('THF', 'Ethereal solvent', 'Good alternative'),
        ('DMF/toluene mixtures', 'Mixed polarity', 'For polar substrates'),
        ('Water (micellar)', 'Green chemistry approach', 'With surfactants like TPGS-750-M'),
    ],
    'limitations': [
        ('β-Hydride elimination (alkyl amines)', 'Primary alkyl amines with β-hydrogens can undergo elimination → enamine/imine byproducts', 'Solution: Use bulky ligands (BrettPhos, RuPhos); lower T; use NaOtBu as base'),
        ('Over-arylation of ammonia/primary amines', 'NH3 → PhNH2 → Ph2NH → Ph3N (multiple arylation)', 'Solution: Use large excess of NH3; protect as carbamate then deprotect'),
        ('Steric hindrance', 'Ortho-substituted aryl halides or bulky amines react slowly', 'Solution: Use BrettPhos/tBuBrettPhos; higher T (100-110°C); longer time'),
        ('Electron-rich aryl halides', 'Electron-donating groups slow oxidative addition', 'Solution: More active catalyst (Pd2(dba)3/XPhos); higher T; longer time'),
        ('Competing hydrodehalogenation', 'Ar-X + H-source → Ar-H (reduction instead of amination)', 'Solution: Optimize base/catalyst ratio; ensure proper degassing'),
        ('N-Arylation vs O-arylation (ambident nucleophiles)', 'Amino alcohols, aminophenols can give mixture of N- and O-arylated products', 'Solution: Protect OH group; choose appropriate ligand/base/solvent'),
        ('Catalyst cost', 'Buchwald ligands are expensive ($50-200/g for some)', 'Solution: Low loading (0.1-0.5 mol%) with G3/G4 precatalysts; ligand recycling'),
        ('Pyridine-type nitrogen interference', 'Substrate-bound pyridine N can coordinate Pd and poison catalyst', 'Solution: Protect as N-oxide; add extra ligand; use more catalyst'),
    ],
    'typical_yields': {
        'aryl_bromide_primary_amine': '75-95%',
        'aryl_bromide_secondary_amine': '80-97%',
        'aryl_chloride_modern_ligand': '65-92%',
        'aniline_arylation': '70-90%',
        'heterocycle_N_arylation': '60-88%',
        'amide_N_arylation': '55-80%',
        'sterically_hindered': '50-82%',
    },
}


@ChemMCPManager.register_tool
class BuchwaldHartwigAmination(BaseTool):
    __version__ = "0.1.0"
    name = "BuchwaldHartwigAmination"
    func_name = 'analyze_buchwald_hartwig_amination'
    description = "Buchwald-Hartwig amination analysis: Pd-catalyzed C-N cross-coupling of aryl halides with amines. Covers full catalytic cycle (oxidative addition, amine coordination/deprotonation, reductive elimination), critical role of ligands (BINAP, XPhos, BrettPhos, RuPhos, DavePhos, NHC categories), 8 catalyst systems, 6 base options, 9 amine types (primary/secondary aliphatic, anilines, N-heterocycles, amides, sulfonamides, hydrazines, ammonia), 7 aryl halide types, 5 solvent systems, scope, limitations (8 items with solutions), and typical yields."
    implementation_description = "Comprehensive knowledge base covering: 3-step catalytic cycle with detailed mechanistic notes, ligand importance theory (4 control mechanisms), 4 ligand categories with examples, 8 catalyst systems (classic BINAP to modern G3 precatalyst), 6 base options, 9 amine substrate categories, 7 aryl halide partner types, 5 solvent systems, 8 limitations with practical solutions, and yield benchmarks."
    categories = ["Reaction"]
    tags = ["Buchwald-Hartwig", "C-N Coupling", "Palladium", "Amination", "Cross-Coupling", "Catalysis", "Amine"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("aryl_halide_smiles", "str", "N/A", "SMILES or name of the aryl halide (iodide, bromide, chloride, triflate)."),
        ("amine_smiles", "str", "N/A", "SMILES or name of the amine coupling partner."),
        ("ligand", "str", "XPhos", "Ligand for Pd catalyst (critical parameter!)."),
        ("base", "str", "NaOtBu", "Base for amine deprotonation."),
        ("solvent", "str", "toluene", "Solvent system."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: aryl_halide amine [ligand] [base] [solvent]. E.g., '4-bromoanisole morpholine XPhos NaOtBu toluene'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing aryl_halide_analysis, amine_analysis, catalytic_cycle, mechanism_steps, ligand_recommendation, catalyst_recommendation, optimal_conditions, scope, limitations, typical_yields, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"aryl_halide_smiles": "4-bromoanisole", "amine_smiles": "morpholine", "ligand": "XPhos", "base": "NaOtBu", "solvent": "toluene"},
            "text_input": {"query": "4-bromoanisole morpholine XPhos NaOtBu toluene"},
            "output": {"result": {
                "reaction": "B-H amination: 4-BrC6H4OMe + morpholine → 4-morpholinoanisole",
                "aryl_halide_analysis": {"type": "aryl bromide", "reactivity": "good", "note": "EDG (OMe) slightly slows oxidative addition"},
                "amine_analysis": {"type": "secondary cyclic aliphatic amine (morpholine)", "reactivity": "good nucleophile"},
                "product": "4-(4-methoxyphenyl)morpholine",
                "catalyst_recommendation": "Pd(OAc)2 (1 mol%) + XPhos (2 mol%) or BrettPhos-G3 precatalyst (0.5 mol%)",
                "conditions": {"T": "100°C", "time": "4-16 h", "atmosphere": "N2"},
                "yield": "85-96%",
                "key_note": "Morpholine is an excellent B-H partner — no β-H elimination issue",
            }},
        },
        {
            "code_input": {"aryl_halide_smiles": "chlorobenzene", "amine_smiles": "piperidine", "ligand": "BrettPhos", "base": "NaOtBu", "solvent": "t-BuOH"},
            "text_input": {"query": "chlorobenzene piperidine BrettPhos NaOtBu t-BuOH"},
            "output": {"result": {
                "reaction": "B-H amination: Ph-Cl + piperidine → N-phenylpiperidine",
                "aryl_halide_analysis": {"type": "aryl chloride", "reactivity": "challenging — needs modern ligand", "note": "Cheapest aryl halide source"},
                "amine_analysis": {"type": "secondary cyclic aliphatic amine (piperidine)", "reactivity": "good"},
                "catalyst_recommendation": "Pd(OAc)2 (1 mol%) + BrettPhos (2 mol%) or BrettPhos-G3 (0.5 mol%)",
                "conditions": {"T": "100-110°C", "time": "8-24 h", "atmosphere": "N2"},
                "yield": "70-90%",
                "note": "Aryl chloride requires Buchwald-type bulky biaryl phosphine ligand",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_BUCHWALD_HARTWIG_DATA)

    def _run_base(self, aryl_halide_smiles: str, amine_smiles: str, ligand: str = "XPhos", base: str = "NaOtBu", solvent: str = "toluene") -> dict:
        if not aryl_halide_smiles:
            raise ChemMCPInputError("Aryl halide is required.")
        if not amine_smiles:
            raise ChemMCPInputError("Amine is required.")

        hal = self._analyze_aryl_halide(aryl_halide_smiles)
        amine = self._analyze_amine(amine_smiles)
        cat_rec = self._recommend_catalyst(hal, amine, ligand)
        cond = self._optimize(hal, amine, ligand, base, solvent)
        limits = self._relevant_limitations(hal, amine)

        result = {
            "result": {
                "reaction": f"B-H amination: {aryl_halide_smiles} + {amine_smiles} → arylated amine product",
                "aryl_halide_analysis": hal,
                "amine_analysis": amine,
                "catalytic_cycle": self.data['catalytic_cycle'],
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['catalytic_cycle']],
                "ligand_importance": self.data['ligand_importance'],
                "catalyst_recommendation": cat_rec,
                "optimal_conditions": cond,
                "scope": self.data['scope'] if hasattr(self.data, 'scope') else list(self.data.get('aryl_halide_scope', [])) + list(self.data.get('amine_scope', [])),
                "applicable_limitations": limits,
                "typical_yields": self._estimate_yield(hal, amine),
                "summary": f"B-H amination predicted: {self._estimate_yield(hal, amine)}. Key concern: {limits[0]['issue'] if limits else 'standard'}.",
            }
        }
        logger.info(f"B-H amination: {aryl_halide_smiles} + {amine_smiles}")
        return result

    def _analyze_aryl_halide(self, smi):
        s = (smi or "").strip().lower()
        hal_types = [
            ('aryl iodide', ['iodo', r'-i\b', r'iodo'], 'Most reactive; expensive'),
            ('aryl bromide', ['bromo', r'-br\b', r'bromo'], 'Standard substrate; good balance'),
            ('aryl chloride', ['chloro', r'-cl\b', r'chloro'], 'Least reactive; needs modern ligands; cheapest'),
            ('aryl triflate', ['triflate', r'-otf', r'otf'], 'Very reactive; from phenols'),
            ('heteroaryl halide', ['bromopyridine', 'bromothiophene', 'chloropyridine', r'pyrid.*br', r'thioph.*br'], 'Heterocyclic; may need optimization'),
        ]
        for htype, pats, note in hal_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": htype, "input": smi, "note": note}
        return {"type": "unknown_aryl_halide", "input": smi}

    def _analyze_amine(self, smi):
        s = (smi or "").strip().lower()
        amine_types = [
            ('primary aliphatic amine', ['butylamine', 'benzylamine', 'hexylamine', r'n-butylamine', r'primary.*amine', r'rnh2'], 'Gives secondary arylamine; watch for β-H elimination'),
            ('secondary aliphatic amine', ['morpholine', 'piperidine', 'pyrrolidine', r'secondary.*amine', r'r2nh'], 'Gives tertiary arylamine; generally clean'),
            ('aniline (aromatic amine)', ['aniline', r'phenylamine', r'phnh2', r'aminobenzene'], 'Gives diarylamine; electron-poor anilines are more reactive'),
            ('N-heterocycle', ['indole', 'carbazole', 'pyrrole', 'imidazole', r'nh-heterocycle'], 'N-arylation of heterocycles; important for pharma'),
            ('amide', ['acetamide', 'benzamide', r'nh2cor', r'conh2'], 'Weaker nucleophile; needs optimization'),
            ('sulfonamide', ['tosylamide', 'sulfonamide', r'nhsO2r'], 'Acidic NH facilitates deprotonation'),
            ('ammonia', ['ammonia', 'nh3'], 'Challenging — can over-arylate to di-/tri-arylamine'),
            ('amino acid ester', ['glycinate', 'alaninate', r'amino.acid'], 'No racemization expected'),
        ]
        for atype, pats, note in amine_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "note": note}
        return {"type": "unknown_amine", "input": smi}

    def _recommend_catalyst(self, hal, amine, ligand):
        ht = hal.get('type', '')
        at = amine.get('type', '')
        lg = (ligand or "XPhos").strip()

        # Aryl chloride needs modern ligand
        if 'chloride' in ht:
            if any(x in lg.lower() for x in ['brettphos', 'xphos', 'ruphos', 'tbuf']):
                return {"primary": f"Pd(OAc)2 (1 mol%) + {lg} (2 mol%) or {lg}-G3 precatalyst (0.5 mol%)",
                        "reason": f"{lg} enables aryl chloride activation"}
            return {"primary": "Pd(OAc)2 (1 mol%) + BrettPhos (2 mol%) or BrettPhos-G3 (0.5 mol%)",
                    "reason": "Aryl chlorides require Buchwald-type bulky biaryl phosphine"}

        # Primary aliphatic amine — recommend RuPhos
        if 'primary' in at and 'aliphatic' in at:
            return {"primary": "Pd2(dba)3 (1 mol%) + RuPhos (2 mol%) or Pd(OAc)2/RuPhos",
                    "reason": "RuPhos excels with primary aliphatic amines (minimizes β-H elimination)"}

        # Anilines — milder conditions
        if 'aniline' in at:
            return {"primary": "Pd2(dba)3 (1 mol%) + XPhos (2 mol%) or DavePhos",
                    "reason": "Anilines work well with standard Buchwald ligands"}

        # N-heterocycles
        if 'heterocycle' in at:
            return {"primary": "Pd(OAc)2 (1 mol%) + XPhos (2 mol%) or BrettPhos (1 mol%)",
                    "reason": "N-Heterocycle arylation works well with modern Buchwald ligands"}

        # Default for standard cases
        return {"primary": f"Pd(OAc)2 (1 mol%) + {lg} (2 mol%) or Pd2(dba)3/{lg}",
                "reason": "Standard B-H conditions for aryl bromides/iodides"}

    def _optimize(self, hal, amine, ligand, base, solvent):
        ht = hal.get('type', '')
        cond = {
            "catalyst_loading": "Pd 0.5-1 mol%" if 'chloride' in ht else "Pd 1-2 mol%",
            "ligand_loading": "2-4 mol% (relative to Pd)",
            "base": f"{base} (1.5-2.5 eq; NaOtBu is standard)",
            "solvent": solvent or "toluene or t-BuOH (both excellent)",
            "concentration": "0.1-0.5 M",
            "temperature": "100-110°C" if 'chloride' in ht else "80-100°C" if 'bromide' in ht else "70-90°C",
            "time": "4-24 hours",
            "atmosphere": "N2 or Ar (strictly oxygen-free)",
            "monitoring": "TLC or GC/MS for consumption of limiting reagent",
            "workup": ("Quench with water, extract with EtOAc (×3), wash with brine, dry (Na2SO4), "
                       "concentrate, purify by column chromatography"),
        }
        return cond

    def _relevant_limitations(self, hal, amine):
        relevant = []
        at = amine.get('type', '')

        if 'primary' in at and 'aliphatic' in at:
            relevant.append({"issue": "β-Hydride elimination", "problem": "Primary alkyl amines with β-H can form enamine/imine byproducts",
                           "Solution": "Use RuPhos or BrettPhos; lower temperature; ensure proper base"})
        if 'ammonia' in at:
            relevant.append({"issue": "Over-arylation", "problem": "NH3 → PhNH2 → Ph2NH → Ph3N (multiple substitutions)",
                           "Solution": "Use large excess of NH3; use protecting group strategy"})
        if 'heterocycle' in at:
            relevant.append({"issue": "Regioselectivity (ambident)", "problem": "Some heterocycles have multiple N sites",
                           "Solution": "Choose appropriate ligand/base; control stoichiometry"})
        ht = hal.get('type', '')
        if 'ortho' in (hal.get('input') or '').lower():
            relevant.append({"issue": "Steric hindrance", "problem": "Ortho-substituted aryl halides are slow",
                           "Solution": "Use tBuBrettPhos; higher T (110°C)"})
        relevant.append({"issue": "Ligand cost", "problem": "Buchwald ligands are expensive specialty chemicals",
                        "Solution": "Use low loading (0.1-0.5 mol%) with G3/G4 precatalysts"})
        return relevant

    def _estimate_yield(self, hal, amine):
        ht = hal.get('type', '')
        at = amine.get('type', '')

        if 'bromide' in ht and 'secondary' in at: return "80-97%"
        if 'bromide' in ht and 'primary' in at: return "75-95%"
        if 'chloride' in ht: return "65-92%"
        if 'aniline' in at: return "70-90%"
        if 'heterocycle' in at: return "60-88%"
        if 'amide' in at: return "55-80%"
        return "70-90%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        hal = parts[0] if parts else ""
        amine = parts[1] if len(parts) > 1 else ""
        lg = parts[2] if len(parts) > 2 else "XPhos"
        base = parts[3] if len(parts) > 3 else "NaOtBu"
        solv = parts[4] if len(parts) > 4 else "toluene"
        return self._run_base(hal, amine, lg, base, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
