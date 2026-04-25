"""
Heck Reaction (Tool #161)
Heck 反应（Mizoroki-Heck 反应）：钯催化烯烃与芳基/乙烯基卤化物的偶联反应。
涵盖：催化循环（氧化加成、烯烃配位/插入、β-氢消除、碱再生）、区域选择性、
立体化学、催化剂体系、碱/溶剂效应、底物范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_HECK_DATA = {
    'catalytic_cycle': [
        ('1', 'Oxidative addition', 'Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X) complex. Rate: RI > RBr >> RCl. Vinyl/Ar halides preferred.'),
        ('2', 'Alkene coordination & migratory insertion', 'Alkene coordinates to Pd(II) center, then R¹ migrates onto the LESS substituted alkene carbon (regioselectivity rule). Forms σ-alkyl-Pd complex.'),
        ('3', 'β-Hydride elimination', 'Syn-periplanar β-H eliminates → forms new C=C double bond (alkene product) + Pd-H. Stereochemistry: usually trans (E) or mixture depending on substrate.'),
        ('4', 'Base-mediated regeneration', 'Base (Et3N, NaOAc, K2CO3, etc.) removes H from Pd-H → HX salt + regenerates Pd(0). Completes catalytic cycle.'),
    ],
    'regioselectivity_rules': {
        'terminal_alkene': 'Aryl/vinyl group adds to the LESS substituted (terminal) carbon → branched product favored',
        'internal_alkene': 'More subtle; electronic and steric factors compete. Electron-withdrawing groups on alkene can reverse selectivity.',
        '1,1_disubstituted_alkene': 'Single regioisomer possible — aryl adds to the CH2 group',
        '1,2_disubstituted_alkene': 'Mixture of regioisomers common; steric bulk influences outcome',
        'electronic_control': 'EWG on alkene directs aryl to the β-position (farther from EWG); EDG can reverse this.',
    },
    'stereochemistry': {
        'general': 'β-Hydride elimination gives predominantly TRANS (E) alkene due to steric factors in the transition state',
        'exceptions': 'Cyclic alkenes give single stereoisomer; certain chiral ligands can induce asymmetry (asymmetric Heck reaction)',
        'asymmetric_heck': 'Possible with chiral phosphine ligands; used in natural product synthesis for quaternary stereocenters',
    },
    'catalyst_systems': {
        'Pd(OAc)2': {'type': 'Pd(II) pre-catalyst', 'loading': '1-5 mol%', 'scope': 'Most widely used; reduced in situ to Pd(0)', 'notes': 'Cheap, versatile'},
        'Pd(PPh3)4': {'type': 'Pd(0) pre-catalyst', 'loading': '2-5 mol%', 'scope': 'Classic catalyst; good for standard substrates', 'notes': 'Air-sensitive'},
        'Pd(dba)2 + P(o-tol)3': {'type': 'Pd(0) + monodentate phosphine', 'loading': '1-3 mol%', 'scope': 'Very active; good for hindered substrates'},
        'Pd(dppf)Cl2': {'type': 'Bidentate phosphine Pd(II)', 'loading': '2-5 mol%', 'scope': 'Good for vinyl halides; stabilizes Pd intermediates'},
        "Herrmann's catalyst": {'type': 'Pd(II) cyclopalladated ferrocenyl ligand', 'loading': '1-2 mol%', 'scope': 'Thermally stable; works at high T (>130°C); air-stable', 'notes': 'Excellent for industrial applications'},
        'Pd-NHC (PEPPSI)': {'type': 'N-heterocyclic carbene Pd', 'loading': '0.5-2 mol%', 'scope': 'Highly active, stable to air/moisture'},
        'Hermann-Beller palladacycle': {'type': 'Thermally robust Pd complex', 'loading': '0.1-1 mol%', 'scope': 'Low loading possible; high temperature tolerant'},
    },
    'base_options': [
        ('Et3N (triethylamine)', 'Most common; acts as both base AND solvent', 'Standard choice for many Heck reactions'),
        ('NaOAc (sodium acetate)', 'Mild base; good in polar aprotic solvents (DMF, NMP)', 'Common alternative'),
        ('K2CO3', 'Moderate base; inexpensive', 'Used with phase-transfer conditions or polar solvents'),
        ('Cs2CO3', 'Stronger, more soluble', 'For challenging substrates'),
        ('iPr2NEt (Hünig\'s base)', 'Sterically hindered; non-nucleophilic', 'Prevents competitive nucleophilic attack'),
        ('NaHCO3', 'Very mild', 'Aqueous conditions; green chemistry approach'),
        ('Inorganic bases (K3PO4, NaOAc)', 'Solid bases', 'Can be used under solvent-free or neat conditions'),
    ],
    'solvent_systems': [
        ('DMF (N,N-dimethylformamide)', 'Most common polar aprotic solvent', 'High boiling point (153°C); dissolves most reagents'),
        ('NMP (N-methyl-2-pyrrolidone)', 'High-boiling polar aprotic', 'Similar to DMF; slightly higher bp (202°C)'),
        ('Acetonitrile (MeCN)', 'Polar aprotic', 'Lower boiling (82°C); good for lower-T reactions'),
        ('DMAc (N,N-dimethylacetamide)', 'High-boiling polar aprotic', 'Alternative to DMF/NMP'),
        ('Water', 'Green solvent', 'Possible with water-soluble ligands/catalysts; surfactant-assisted'),
        ('Ionic liquids', 'Tunable solvents', 'Specialty applications; catalyst recycling'),
        ('Neat (solvent-free)', 'Industrial approach', 'Alkene as solvent/reagent; large-scale processes'),
        ('Et3N (as solvent/base)', 'Dual function', 'Classic Mizoroki-Heck conditions'),
    ],
    'scope': [
        'Aryl iodides: EXCELLENT — most reactive, work at 80-100°C',
        'Aryl bromides: GOOD — standard substrates; 100-140°C typically',
        'Aryl chlorides: POSSIBLE — need modern bulky ligands or high T (>140°C)',
        'Aryl triflates: VERY GOOD — from phenols; reactive like iodides',
        'Vinyl halides/triflates: EXCELLENT — stereochemistry retained (important!)',
        'Terminal alkenes: EXCELLENT — clean regioselectivity (aryl to terminal C)',
        'Acrylates/acrylonitrile: EXCELLENT — electron-poor alkenes are highly reactive',
        'Cyclic alkenes (norbornene, etc.): EXCELLENT — endo/exo selectivity interesting',
        'Heteroaryl halides: GOOD — pyridine, thiophene, furan derivatives work',
    ],
    'limitations': [
        ('Over-alkylation', 'Product alkene can undergo SECOND Heck reaction → diarylated byproduct', 'Solution: Use excess alkene (2-5 eq); monitor conversion carefully'),
        ('Regioselectivity issues', 'Internal/disubstituted alkenes give mixtures of regioisomers', 'Solution: Use directing groups on alkene; choose appropriate ligand'),
        ('Isomerization', 'Double bond can isomerize under reaction conditions (Pd-H addition/elimination)', 'Solution: Lower temperature; shorter reaction time; use appropriate ligand'),
        ('β-Hydride from alkyl halides', 'Alkyl halides can undergo competing β-hydride elimination', 'Solution: Use aryl/vinyl halides primarily; special conditions for alkyl'),
        ('Steric hindrance', 'Ortho-substituted aryl halides or tetrasubstituted alkenes are slow', 'Solution: Bulky ligands (PtBu3, XPhos); higher T; longer time'),
        ('Heteroatom poisoning', 'Amines, thiols, phosphines can coordinate Pd and deactivate catalyst', 'Solution: Protect heteroatoms; add extra ligand'),
        ('Decomposition of Pd(0)', 'Pd can form Pd black (inactive precipitate) especially at high T', 'Solution: Use stabilizing ligands; proper degassing; avoid impurities'),
        ('Cost', 'Pd catalysts expensive; ligands add cost', 'Solution: Low loading (0.1-1 mol%) with active catalysts; Pd recycling'),
    ],
    'typical_yields': {
        'aryl_iodide_terminal_alkene': '75-95%',
        'aryl_bromide_terminal_alkene': '70-92%',
        'aryl_chloride_modern_ligand': '60-85%',
        'vinyl_halide': '70-93%',
        'acrylate_ester': '80-98%',
        'heteroaryl': '55-85%',
        'ortho_substituted': '50-80%',
        'cyclic_alkene': '65-90%',
    },
}


@ChemMCPManager.register_tool
class HeckReaction(BaseTool):
    __version__ = "0.1.0"
    name = "HeckReaction"
    func_name = 'analyze_heck_reaction'
    description = "Heck (Mizoroki-Heck) reaction analysis: Pd-catalyzed coupling of alkenes with aryl/vinyl halides. Covers full catalytic cycle (oxidative addition, alkene coordination/migratory insertion, β-hydride elimination, base regeneration), regioselectivity rules (aryl to less-substituted terminus), stereochemistry (trans-favored), 7 catalyst systems (from classic Pd(OAc)2 to PEPPSI), 7 base options, 8 solvent systems, scope (9 categories), limitations (8 items with solutions), and typical yields."
    implementation_description = "Comprehensive knowledge base covering: 4-step catalytic cycle with detailed mechanistic notes, regioselectivity rules (5 categories), stereochemistry guidelines including asymmetric Heck, 7 catalyst systems, 7 base options with roles, 8 solvent systems, 9 scope categories, 8 limitations with practical solutions, and yield benchmarks."
    categories = ["Reaction"]
    tags = ["Heck", "Cross-Coupling", "Palladium", "Alkene", "C-C Bond Formation", "Catalysis", "Mizoroki-Heck"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("halide_smiles", "str", "N/A", "SMILES or name of the organic halide (aryl/vinyl iodide, bromide, chloride, or triflate)."),
        ("alkene_smiles", "str", "N/A", "SMILES or name of the alkene coupling partner."),
        ("ligand", "str", "PPh3", "Ligand for Pd catalyst (optional, leave blank for ligandless)."),
        ("base", "str", "Et3N", "Base for Pd-H elimination / HX scavenging."),
        ("solvent", "str", "DMF", "Solvent system."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: halide alkene [ligand] [base] [solvent]. E.g., 'iodobenzene styrene PPh3 Et3N DMF'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing halide_analysis, alkene_analysis, catalytic_cycle, mechanism_steps, regioselectivity_prediction, stereochemistry, catalyst_recommendation, optimal_conditions, scope, limitations, typical_yields, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"halide_smiles": "iodobenzene", "alkene_smiles": "styrene", "ligand": "PPh3", "base": "Et3N", "solvent": "DMF"},
            "text_input": {"query": "iodobenzene styrene PPh3 Et3N DMF"},
            "output": {"result": {
                "reaction": "Heck reaction: Ph-I + PhCH=CH2 → stilbene (E-major)",
                "halide_analysis": {"type": "aryl iodide", "reactivity": "excellent (most reactive)", "note": "Standard Heck substrate"},
                "alkene_analysis": {"type": "terminal alkene (styrene)", "regioselectivity": "aryl adds to β-carbon (less substituted)"},
                "product": "(E)-stilbene (major) + minor Z-isomer",
                "catalyst_recommendation": "Pd(OAc)2 (2 mol%) + PPh3 (4 mol%) or Pd(PPh3)4 (3 mol%)",
                "conditions": {"T": "100°C", "time": "4-16 h", "atmosphere": "N2"},
                "yield": "80-95%",
                "stereochemistry": "E-alkene major (trans β-hydride elimination)",
            }},
        },
        {
            "code_input": {"halide_smiles": "4-bromoacetophenone", "alkene_smiles": "n-butyl acrylate", "ligand": "P(o-tol)3", "base": "NaOAc", "solvent": "DMF"},
            "text_input": {"query": "4-bromoacetophenone butyl_acrylate P(o-tol)3 NaOAc DMF"},
            "output": {"result": {
                "reaction": "Heck reaction: 4-BrC6H4COMe + CH2=CHCOOBu → cinnamic acid derivative",
                "product": "butyl (E)-3-(4-acetylphenyl)acrylate",
                "regioselectivity": "Clean — acrylate has electron-withdrawing group, aryl adds to β-position (terminal C)",
                "yield": "82-96%",
                "note": "Electron-poor alkenes like acrylates are excellent Heck partners",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_HECK_DATA)

    def _run_base(self, halide_smiles: str, alkene_smiles: str, ligand: str = "PPh3", base: str = "Et3N", solvent: str = "DMF") -> dict:
        if not halide_smiles:
            raise ChemMCPInputError("Organic halide is required.")
        if not alkene_smiles:
            raise ChemMCPInputError("Alkene coupling partner is required.")

        hal = self._analyze_halide(halide_smiles)
        alk = self._analyze_alkene(alkene_smiles)
        regio = self._predict_regioselectivity(alk)
        stereo = self._predict_stereochemistry(alk)
        cat_rec = self._recommend_catalyst(hal, ligand)
        cond = self._optimize(hal, alk, ligand, base, solvent)
        limits = self._relevant_limitations(hal, alk)

        result = {
            "result": {
                "reaction": f"Heck reaction: {halide_smiles} + {alkene_smiles} → coupled alkene product",
                "halide_analysis": hal,
                "alkene_analysis": alk,
                "catalytic_cycle": self.data['catalytic_cycle'],
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['catalytic_cycle']],
                "regioselectivity_prediction": regio,
                "stereochemistry": stereo,
                "catalyst_recommendation": cat_rec,
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "typical_yields": self._estimate_yield(hal, alk),
                "summary": f"Heck reaction predicted: {self._estimate_yield(hal, alk)}. Regio: {regio.get('prediction','?')}. Stereo: {stereo.get('prediction','?')}. Key concern: {limits[0]['issue'] if limits else 'standard'}.",
            }
        }
        logger.info(f"Heck: {halide_smiles} + {alkene_smiles}")
        return result

    def _analyze_halide(self, smi):
        s = (smi or "").strip().lower()
        hal_types = [
            ('aryl iodide', ['iodo', r'-i\b', r'iodo'], 'Most reactive; expensive; standard Heck substrate'),
            ('aryl bromide', ['bromo', r'-br\b', r'bromo'], 'Standard substrate; good balance of cost/reactivity'),
            ('aryl chloride', ['chloro', r'-cl\b', r'chloro'], 'Least reactive; needs modern ligands/high T; cheapest'),
            ('aryl triflate', ['triflate', r'-otf', r'otf'], 'Very reactive; from phenols via Tf2O'),
            ('vinyl halide', ['vinyl.*br', r'ch=ch.*hal', r'vinyliod', r'vinylbro'], 'Stereochemistry RETAINED (cis→cis, trans→trans)'),
            ('heteroaryl halide', ['bromopyridine', 'bromothiophene', 'bromofuran', r'pyrid.*br', r'thioph.*br'], 'Heterocyclic; may need optimization'),
        ]
        for htype, pats, note in hal_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": htype, "input": smi, "note": note}
        return {"type": "unknown_halide", "input": smi}

    def _analyze_alkene(self, smi):
        s = (smi or "").strip().lower()
        alk_types = [
            ('terminal alkene', ['ethylene', 'propene', r'^\w*acrylat', r'acrylonitrile', r' styrene', r'terminal', r'ch2=ch'], 'Excellent regioselectivity; aryl adds to terminal carbon'),
            ('acrylate ester', ['acrylate', 'methyl_acrylate', 'butyl_acrylate', 'ethyl_acrylate', r'ch2=chcoo'], 'Electron-poor; very reactive; clean regioselectivity'),
            ('styrenic alkene', ['styrene', r'phch=ch2', r'phenylethene'], 'Terminal; gives stilbene derivatives'),
            ('electron-rich alkene', ['vinyl ether', 'enamine', r'chor=ch'], 'Reactive but may polymerize'),
            ('cyclic alkene', ['norbornene', 'cyclohexene', 'cyclopentene', r'cyclic'], 'Endo selectivity often observed'),
            ('internal alkene', ['2-butene', 'internal', r'ch3ch=chch3'], 'Regioselectivity may be problematic'),
            ('1,1-disubstituted', ['isobutylene', r'ch2=c(ch3)2'], 'Single regioisomer expected'),
        ]
        for atype, pats, note in alk_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "note": note}
        return {"type": "unknown_alkene", "input": smi}

    def _predict_regioselectivity(self, alk):
        at = alk.get('type', '')
        rules = self.data['regioselectivity_rules']

        if 'terminal' in at or 'acrylate' in at or 'styrenic' in at:
            return {"prediction": "Clean — aryl/vinyl group adds to LESS substituted (terminal) carbon",
                    "rule": rules.get('terminal_alkene', '')}
        elif 'cyclic' in at:
            return {"prediction": "Depends on ring size and substitution pattern; often single isomer",
                    "rule": "Cyclic alkenes have constrained geometry"}
        elif 'internal' in at or '1,1' in at:
            return {"prediction": "May give mixture of regioisomers",
                    "rule": rules.get('internal_alkene', '')}
        return {"prediction": "Follows standard Heck regioselectivity (aryl to less substituted end)",
                "rule": "General rule applies"}

    def _predict_stereochemistry(self, alk):
        at = alk.get('type', '')
        stereo = self.data['stereochemistry']

        if 'acrylate' in at:
            return {"prediction": "E (trans) alkene strongly favored", "rule": "β-Hydride elimination gives trans product; acrylates give E-selectivity >20:1"}
        elif 'terminal' in at:
            return {"prediction": "E (trans) favored for disubstituted products; no E/Z for terminal", "rule": stereo['general']}
        elif 'vinyl' in at.lower():
            return {"prediction": "Vinyl configuration largely RETAINED", "rule": "Vinyl halides transfer configuration through oxidative addition"}
        return {"prediction": "E (trans) alkene generally favored", "rule": stereo['general']}

    def _recommend_catalyst(self, hal, ligand):
        ht = hal.get('type', '')

        if 'chloride' in ht:
            return {"primary": "Pd(OAc)2 (2 mol%) + SPhos/XPhos/PtBu3 (4 mol%) or Herrmann's catalyst (1-2 mol%)",
                    "reason": "Aryl chlorides require modern bulky ligands or thermally robust palladacycles"}
        if 'triflate' in ht:
            return {"primary": "Pd(OAc)2 (1-2 mol%) + P(o-tol)3 (2-4 mol%) or Pd(PPh3)4 (2 mol%)",
                    "reason": "Triflates are very reactive; standard catalysts suffice"}
        if 'vinyl' in ht:
            return {"primary": "Pd(dppf)Cl2 (2-3 mol%) or Pd(OAc)2/PPh3",
                    "reason": "Bidentate ligands help retain vinyl stereochemistry"}

        lg = (ligand or "PPh3").strip().lower()
        if any(x in lg for x in ['sphos', 'xphos', 'ruphos']):
            return {"primary": f"Pd(OAc)2 (1 mol%) + {ligand} (2 mol%)",
                    "reason": "Modern Buchwald-type ligand — excellent activity"}
        else:
            return {"primary": f"Pd(OAc)2 (2 mol%) + {ligand} (4 mol%) or Pd(PPh3)4 (2-3 mol%)",
                    "reason": "Standard catalyst system for aryl bromides/iodides"}

    def _optimize(self, hal, alk, ligand, base, solvent):
        ht = hal.get('type', '')
        at = alk.get('type', '')
        cond = {
            "catalyst_loading": "0.5-2 mol% Pd" if 'chloride' in ht else "1-3 mol% Pd (standard)",
            "ligand_loading": "2-4 mol% (relative to Pd); ligandless also possible for simple cases",
            "base": f"{base} (2-3 eq; serves as HX scavenger)",
            "solvent": solvent or "DMF or NMP (high-boiling polar aprotic)",
            "concentration": "0.1-0.5 M",
            "temperature": "130-140°C" if 'chloride' in ht else "100-120°C" if 'bromide' in ht else "80-110°C",
            "time": "4-24 hours",
            "atmosphere": "N2 or Ar (degassed solvents recommended)",
            "alkene_equivalents": "1.5-3 eq (excess prevents over-arylation of product)",
            "monitoring": "TLC or GC/MS for consumption of limiting reagent",
            "workup": ("Dilute with water, extract with EtOAc (×3), wash with brine, dry (Na2SO4), "
                       "concentrate, purify by column chromatography"),
        }
        return cond

    def _relevant_limitations(self, hal, alk):
        relevant = []
        ht = hal.get('type', '')
        at = alk.get('type', '')

        relevant.append({"issue": "Over-alkylation", "problem": "Product alkene can undergo a second Heck reaction → diarylated byproduct",
                        "solution": "Use excess alkene (1.5-3 eq); monitor carefully; isolate product promptly"})
        if 'internal' in at or 'disubstituted' in at:
            relevant.append({"issue": "Regioselectivity", "problem": "Non-terminal alkenes can give mixtures of regioisomers",
                           "solution": "Use directing groups; optimize ligand/solvent; accept mixture"})
        if 'vinyl' not in ht and 'alkyl' not in ht:
            pass  # standard case
        if 'ortho' in (hal.get('input') or '').lower():
            relevant.append({"issue": "Steric hindrance", "problem": "Ortho-substituted aryl halides react slowly",
                           "solution": "Use bulky ligands (XPhos, PtBu3); higher T; longer time"})
        relevant.append({"issue": "Pd black formation", "problem": "Pd(0) can precipitate as inactive Pd black at high T",
                        "solution": "Use stabilizing ligands; ensure proper degassing; add fresh Pd if needed"})
        return relevant

    def _estimate_yield(self, hal, alk):
        ht = hal.get('type', '')
        at = alk.get('type', '')

        if 'iodide' in ht and ('terminal' in at or 'acrylate' in at): return "85-98%"
        if 'bromide' in ht and ('terminal' in at or 'acrylate' in at): return "80-95%"
        if 'chloride' in ht: return "60-85%"
        if 'triflate' in ht: return "78-94%"
        if 'vinyl' in ht: return "70-93%"
        if 'acrylate' in at: return "80-98%"
        if 'hetero' in ht: return "55-85%"
        return "70-90%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        hal = parts[0] if parts else ""
        alk = parts[1] if len(parts) > 1 else ""
        lg = parts[2] if len(parts) > 2 else "PPh3"
        base = parts[3] if len(parts) > 3 else "Et3N"
        solv = parts[4] if len(parts) > 4 else "DMF"
        return self._run_base(hal, alk, lg, base, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
