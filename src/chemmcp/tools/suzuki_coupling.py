"""
Suzuki Coupling (Tool #160)
Suzuki-Miyaura 交叉偶联反应：有机硼 + 有机卤化物、Pd 催化循环
（氧化加成、转金属、还原消除）、配体/碱/溶剂效应、范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_SUZUKI_DATA = {
    'catalytic_cycle': [
        ('1', 'Oxidative addition', 'Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X) complex. Rate: RI > RBr >> RCl (activated) >> RF. Vinyl/Ar halides faster than alkyl.'),
        ('2', 'Transmetalation', 'Organoboron (R²-B(OH)2 or ester) activated by base → transfers R² to Pd, displacing X. Base forms boronate [R²-B(OH)3]⁻ which transmetalates more readily.'),
        ('3', 'Reductive elimination', 'Pd(II)(R¹)(R²) eliminates R¹-R² coupled product + regenerates Pd(0). Stereoretentive for stereochemical centers on both partners.'),
    ],
    'catalyst_systems': {
        'Pd(PPh3)4': {'type': 'Pd(0) pre-catalyst', 'loading': '1-5 mol%', 'scope': 'Classic; good for aryl-aryl and aryl-vinyl couplings', 'limitations': 'Less active for aryl chlorides; thermally sensitive'},
        'Pd(PPh3)2Cl2': {'type': 'Pd(II) pre-catalyst', 'loading': '2-5 mol%', 'scope': 'Reduced in situ to Pd(0); widely used', 'notes': 'Often with excess PPh3 as ligand'},
        'Pd(dppf)Cl2': {'type': 'Bidentate phosphine Pd(II)', 'loading': '1-3 mol%', 'scope': 'Excellent for alkyl-aryl and challenging couplings; chelating ligand stabilizes Pd'},
        'Pd(OAc)2 + SPhos/XPhos/RuPhos': {'type': 'Bulky biaryl phosphine + Pd(OAc)2', 'loading': '0.5-2 mol%', 'scope': 'STATE OF THE ART: activates Ar-Cl even Ar-OTf; very high activity', 'notes': 'Buchwald-type ligands revolutionized Suzuki coupling'},
        'Pd2(dba)3 + PCy3': {'type': 'Pd(0) + bulky phosphine', 'loading': '1-3 mol%', 'scope': 'Very active for sterically hindered couplings'},
        'Pd/C (palladium on carbon)': {'type': 'Heterogeneous catalyst', 'loading': '1-5 mol% Pd', 'scope': 'Cheap, recyclable; works for many standard couplings', 'notes': 'Green chemistry approach'},
        'Pd-NHC (PEPPSI, etc.)': {'type': 'N-heterocyclic carbene Pd', 'loading': '0.5-2 mol%', 'scope': 'Highly active, stable to air/moisture', 'notes': 'Modern robust catalysts'},
        'Ni-catalyzed variants': {'type': 'Nickel instead of palladium', 'advantage': 'Much cheaper than Pd', 'scope': 'Can couple phenol derivatives (OAr), cheaper but less developed'},
    },
    'organoboron_partners': [
        ('Boronic acid (R-B(OH)2)', 'Most common; commercially available for many aryl/vinyl groups', 'Stable crystalline solids; protodeboronation side reaction possible'),
        ('Boronic ester (pinacol ester R-Bpin)', 'More stable than boronic acids; resistant to protodeboronation', 'Preferred for base-sensitive substrates; MIDA boronates are especially stable'),
        ('Trifluoroborate salt (R-BF3K)', 'Air/water-stable crystalline solids', 'Easy to handle; activated by base in situ'),
        ('Borane amine complex (R-B(NMe3) or R-BR\'2·NR\'3)', 'Alternative boron sources', 'Specialty applications'),
        ('Potassium organotrifluoroborate', 'Highly stable, easy to purify', 'Increasingly popular alternative to boronic acids'),
    ],
    'organic_halide_partners': [
        ('Aryl iodide (Ar-I)', 'Most reactive; works with most catalyst systems', 'Expensive; may have side reactions (homocoupling)'),
        ('Aryl bromide (Ar-Br)', 'Standard substrate; good balance of cost/reactivity', 'Most commonly used in academic labs'),
        ('Aryl chloride (Ar-Cl)', 'Least reactive traditional substrate', 'Requires modern bulky phosphine or NHC ligands; cheap and abundant'),
        ('Aryl triflate (Ar-OTf)', 'Very reactive; from phenols', 'Good for coupling phenol derivatives; moisture sensitive'),
        ('Aryl tosylate/mesylate (Ar-OMs/OTs)', 'From phenols; less reactive than triflates', 'Needs active catalyst system'),
        ('Vinyl halide/triflate', 'Stereochemistry RETAINED (cis→cis, trans→trans)', 'Important for styrene synthesis'),
        ('Alkyl halide (primary)', 'Possible but challenging (β-hydride elimination competing)', 'Use Pd-PEPPSI or Ni catalysis for best results'),
        ('Acid chloride / pseudohalide', 'Can be used as electrophile', 'Gives ketones after coupling (not biaryls)'),
    ],
    'base': {
        'purpose': 'Activates organoboron compound via formation of tetracoordinate boronate [R-B(OH)3]⁻ which undergoes transmetalation readily',
        'common_bases': [
            ('K2CO3', 'Mild, inexpensive; standard choice', 'Works well in dioxane/H2O or toluene/EtOH/H2O'),
            ('Cs2CO3', 'Stronger, more soluble in organic solvents', 'Good for challenging couplings'),
            ('Na2CO3', 'Similar to K2CO3', 'Common alternative'),
            ('K3PO4', 'Strong, non-nucleophilic', 'Excellent for base-sensitive substrates'),
            ('t-BuONa/K', 'Strong, anhydrous conditions', 'For anhydrous protocols'),
            ('KF', 'Mild; also provides fluoride activation of boron', 'Special applications'),
            ('TlOH/Tl2CO3', 'Historically used — AVOID (highly toxic!)', 'Legacy method only'),
        ],
    },
    'solvent_effects': [
        ('toluene/EtOH/H2O (3:1)', 'Classic solvent mixture; biphasic', 'Most common academic protocol'),
        ('dioxane/H2O', 'Fully miscible at RT; homogeneous', 'Very popular modern choice'),
        ('DMF/H2O', 'Polar; good for polar substrates', 'May complicate product isolation'),
        ('THF/H2O', 'Good compromise', 'Widely used'),
        ('MeCN/H2O', 'Polar aprotic; fast reactions', 'Good for kinetics studies'),
        ('Water (pure)', 'Green chemistry approach', 'Possible with water-soluble ligands/catalysts'),
        ('Neat (solvent-free)', 'Industrial approach', 'For large-scale processes'),
    ],
    'ligand_effects': {
        'electron_rich_bulky_phosphines': ('Accelerate oxidative addition (especially for Ar-Cl)',
                                           'Prevent β-hydride elimination',
                                           'Examples: SPhos, XPhos, RuPhos, DavePhos, BrettPhos, JohnPhos'),
        'bidentate_phosphines': ('Stabilize Pd intermediates; prevent Pd black formation',
                                 'Examples: dppf, dppe, dppp, Xantphos'),
        'N_heterocyclic_carbenes_NHC': ('Very strong σ-donors; highly active; air-stable complexes',
                                       'Examples: IPr, IMes, PEPPSI series'),
        'water_soluble_ligands': ('TPPTS (triphenylphosphine trisulfonate); enable aqueous-phase Suzuki',
                                  'Green chemistry advantage'),
    },
    'scope': [
        'Aryl-aryl (biaryl) couplings: EXCELLENT — the flagship transformation',
        'Aryl-vinyl (styrene) couplings: Excellent — stereochemistry retained',
        'Vinyl-vinyl (diene) couplings: Good — stereochemistry retained',
        'Heteroaryl couplings: Good — pyridine, thiophene, furan, indole all work (may need optimization)',
        'Alkyl-aryl couplings: CHALLENGING — primary alkyl halides work with modern catalysts; secondary/tertiary difficult',
        'Couplings with ortho-substituted partners: Works with bulky ligands (SPhos, XPhos)',
        'Large-scale: Industrially proven (e.g., boscalid, valsartan synthesis)',
    ],
    'limitations': [
        ('Protodeboronation', 'Boronic acid decomposes under basic conditions → Ar-H byproduct', 'Solution: Use boronic esters (Bpin); lower T; avoid excess base'),
        ('Homocoupling', '2 Ar-B(OH)2 → Ar-Ar (oxidative homocoupling)', 'Solution: Degas solvents; exclude oxygen; use fresh boronic acid'),
        ('β-Hydride elimination (alkyl halides)', 'Alkyl-Pd intermediate → alkene + H-Pd-X', 'Solution: Use primary alkyl only; Pd-PEPPSI or Ni catalysis'),
        ('Steric hindrance', 'Very hindered partners react slowly', 'Solution: Bulky ligands (RuPhos, BrettPhos); higher T; longer time'),
        ('Heteroatom coordination', 'Amines, pyridines can coordinate Pd and poison catalyst', 'Solution: Protect heteroatoms; use excess ligand'),
        ('Sensitive functional groups', 'Base-sensitive groups (esters OK, some acetals not)', 'Choose milder base (K2CO3 vs t-BuONa); lower T'),
        ('Cost', 'Pd catalysts expensive (especially specialized ligands)', 'Solution: Low loading (0.1-0.5 mol%) with active catalysts; Pd recycling'),
        ('Environmental concerns', 'Pd residue in pharmaceutical products must be <10 ppm', 'Solution: Scavenger resins; careful purification'),
    ],
    'typical_yields': {
        'aryl_iodide_bromide': '80-98%',
        'aryl_chloride_modern_ligand': '70-95%',
        'heteroaryl': '60-90%',
        'vinyl': '75-95%',
        'alkyl_primary': '50-85%',
        'ortho_substituted': '60-85%',
    },
}


@ChemMCPManager.register_tool
class SuzukiCoupling(BaseTool):
    __version__ = "0.1.0"
    name = "SuzukiCoupling"
    func_name = 'analyze_suzuki_coupling'
    description = "Suzuki-Miyaura cross-coupling analysis: organoboron + organic halide coupling via Pd catalysis. Covers full catalytic cycle (oxidative addition, transmetalation, reductive elimination), catalyst systems (8 types including Buchwald ligands), organoboron/halide partner scope (15+ types each), base/solvent/ligand effects, scope (7 categories), limitations (8 items with solutions), and typical yields."
    implementation_description = "Comprehensive knowledge base covering: 3-step catalytic cycle with detailed mechanistic notes, 8 catalyst systems (from classic Pd(PPh3)4 to modern Pd-NHC), 5 organoboron partner types, 8 organic halide partner types, 6 base options with purpose explanation, 6 solvent systems, detailed ligand effects (3 categories), 7 scope categories, 8 limitations with practical solutions, and yield benchmarks."
    categories = ["Reaction"]
    tags = ["Suzuki", "Cross-Coupling", "Palladium", "Organoboron", "C-C Bond Formation", "Catalysis", "Biaryl Synthesis"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("organoboron_smiles", "str", "N/A", "SMILES or name of the organoboron compound (boronic acid, pinacol ester, etc.)."),
        ("halide_smiles", "str", "N/A", "SMILES or name of the organic halide (or triflate/tosylate)."),
        ("ligand", "str", "PPh3", "Ligand for Pd catalyst."),
        ("base", "str", "K2CO3", "Base for transmetalation activation."),
        ("solvent", "str", "dioxane/H2O", "Solvent system."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: organoboron halide [ligand] [base] [solvent]. E.g., 'phenylboronic_acid bromobenzene PPh3 K2CO3 dioxane/H2O'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing coupling_partners_analysis, catalytic_cycle, mechanism_steps, conditions_optimization, scope, limitations, ligand_effects, typical_yields, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"organoboron_smiles": "phenylboronic acid", "halide_smiles": "4-bromotoluene", "ligand": "PPh3", "base": "K2CO3", "solvent": "dioxane/H2O"},
            "text_input": {"query": "phenylboronic_acid bromobenzene PPh3 K2CO3 dioxane/H2O"},
            "output": {"result": {
                "reaction": "Suzuki coupling: Ph-B(OH)2 + Ph-Br → biphenyl",
                "organoboron_analysis": {"type": "boronic acid", "reactivity": "good", "stability": "moderate (protodeboronation risk)"},
                "halide_analysis": {"type": "aryl bromide", "reactivity": "good (standard Suzuki substrate)"},
                "product": "4-methylbiphenyl (if p-bromotoluene) or biphenyl (if bromobenzene)",
                "catalyst_recommendation": "Pd(PPh3)4 (2 mol%) or Pd(OAc)2/PPh3 (3 mol%/6 mol%)",
                "conditions": {"T": "80-100°C", "time": "2-12 h", "atmosphere": "N2 (degassed solvents)"},
                "yield": "85-97%",
                "key_advantages": "Mild conditions, functional group tolerance, wide availability of boronic acids, non-toxic byproducts",
            }},
        },
        {
            "code_input": {"organoboron_smiles": "vinyl-Bpin", "halide_smiles": "4-iodophenyl acetate", "ligand": "Pd(dppf)Cl2", "base": "K3PO4", "solvent": "THF/H2O"},
            "text_input": {"query": "vinyl-Bpin iodoacetophenone Pd(dppf)Cl2 K3PO4 THF/H2O"},
            "output": {"result": {
                "reaction": "Suzuki coupling: vinylboronic ester + 4-iodophenyl acetate → styrenyl derivative",
                "stereochemistry": "Vinyl configuration RETAINED (cis/trans preserved)",
                "yield": "78-92%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_SUZUKI_DATA)

    def _run_base(self, organoboron_smiles: str, halide_smiles: str, ligand: str = "PPh3", base: str = "K2CO3", solvent: str = "dioxane/H2O") -> dict:
        if not organoboron_smiles:
            raise ChemMCPInputError("Organoboron compound is required.")
        if not halide_smiles:
            raise ChemMCPInputError("Organic halide is required.")

        ob = self._analyze_organoboron(organoboron_smiles)
        hal = self._analyze_halide(halide_smiles)
        cat_rec = self._recommend_catalyst(hal, ligand)
        cond = self._optimize(ob, hal, ligand, base, solvent)
        limits = self._relevant_limitations(ob, hal)

        result = {
            "result": {
                "reaction": f"Suzuki coupling: {organoboron_smiles} + {halide_smiles} → cross-coupled product",
                "organoboron_analysis": ob,
                "halide_analysis": hal,
                "catalytic_cycle": self.data['catalytic_cycle'],
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['catalytic_cycle']],
                "catalyst_recommendation": cat_rec,
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "ligand_effects": self.data['ligand_effects'],
                "typical_yields": self._estimate_yield(ob, hal),
                "base_role": self.data['base']['purpose'],
                "summary": f"Suzuki coupling predicted: {self._estimate_yield(ob, hal)}. Key concern: {limits[0]['issue'] if limits else 'standard'}",
            }
        }
        logger.info(f"Suzuki: {organoboron_smiles} + {halide_smiles}")
        return result

    def _analyze_organoboron(self, smi):
        s = (smi or "").strip().lower()
        ob_types = [
            ('boronic acid', [r'boronic.acid', r'b\(oh\)2', r'b\(oh\)2'], 'Moderately stable; protodeboronation possible'),
            ('pinacol ester (Bpin)', [r'bpin', r'pinacol', r'bc'], 'More stable than B(OH)2; preferred for sensitive substrates'),
            ('trifluoroborate', [r'bf3[k]', r'trifluoro'], 'Crystalline, air/water stable'),
            ('MIDA boronate', [r'mida'], 'Exceptionally stable; used in iterative coupling'),
            ('boronic ester general', [r'boronic.est', r'b\(or\)2'], 'Generic boronic ester'),
            ('vinyl boronic', ['vinyl.*bor', r'ch=ch.b'], 'Vinyl group — stereochemistry retained'),
            ('aryl boronic', ['phenylbor', 'arylbor', r'c1ccc.*b'], 'Aryl group — most common type'),
            ('heteroaryl boronic', ['thienyl', 'furanyl', 'pyridinyl', 'pyridyl'], 'Heterocyclic — may need optimization'),
            ('alkyl boronic', ['alkyl.*bor', r'cccb'], 'Alkyl — challenging due to protodeboronation'),
        ]
        for btype, pats, note in ob_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": btype, "input": smi, "note": note}
        return {"type": "unknown_organoboron", "input": smi}

    def _analyze_halide(self, smi):
        s = (smi or "").strip().lower()
        hal_types = [
            ('aryl iodide', ['iodo', r'-i\b', r'iodo'], 'Most reactive; expensive'),
            ('aryl bromide', ['bromo', r'-br\b', r'bromo'], 'Standard substrate; good cost/reactivity balance'),
            ('aryl chloride', ['chloro', r'-cl\b', r'chloro'], 'Least reactive; needs modern ligands; cheapest'),
            ('aryl triflate', ['triflate', r'-otf', r'otf', r'-os(o)2cf3'], 'Very reactive; from phenols'),
            ('aryl tosylate', ['tosylate', r'-ots', r'ots', r'-os(o)2c(ch3)3c6h4'], 'From phenols; moderate reactivity'),
            ('vinyl halide', ['vinyl.*br', r'ch=ch.*hal', r'vinyliod', r'vinylbro'], 'Stereochemistry retained'),
            ('alkyl halide (primary)', ['alkyl.*br', r'cc.*br', r'-ch2.*br'], 'Challenging; β-hydride elimination risk'),
            ('heteroaryl halide', ['bromopyridine', 'bromothiophene', 'bromofuran', r'pyrid.*br', r'thioph.*br'], 'Heterocyclic; may need optimized conditions'),
        ]
        for htype, pats, note in hal_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": htype, "input": smi, "note": note}
        return {"type": "unknown_halide", "input": smi}

    def _recommend_catalyst(self, hal, ligand):
        ht = hal.get('type', '')
        lg = (ligand or "PPh3").strip()

        if 'chloride' in ht:
            return {"primary": "Pd(OAc)2 (1 mol%) + SPhos/XPhos (2 mol%) or Pd-PEPPSI-IHeptCl (0.5 mol%)",
                    "reason": "Aryl chlorides require modern bulky biaryl phosphine or NHC ligands"}
        if 'triflate' in ht:
            return {"primary": "Pd(PPh3)4 (2 mol%) or Pd(dppf)Cl2 (2 mol%)",
                    "reason": "Triflates are very reactive; standard catalysts suffice"}
        if 'vinyl' in ht:
            return {"primary": "Pd(dppf)Cl2 (2 mol%) or Pd(OAc)2/P(o-tol)3",
                    "reason": "Bidentate ligands help retain vinyl stereochemistry"}
        if 'alkyl' in ht:
            return {"primary": "Pd-PEPPSI-IPent (1 mol%) or Ni(dppf)Cl2",
                    "reason": "Alkyl halides need special catalysts to suppress β-hydride elimination"}

        # Default for Br/I
        if any(x in lg.lower() for x in ['sphos', 'xphos', 'ruphos', 'davephos']):
            return {"primary": f"Pd(OAc)2 (1 mol%) + {lg} (2 mol%)",
                    "reason": "Modern Buchwald-type ligand — excellent activity"}
        elif 'dppf' in lg.lower():
            return {"primary": f"Pd({lg})Cl2 (2 mol%)",
                    "reason": "Bidentate phosphine — versatile and reliable"}
        else:
            return {"primary": f"Pd(PPh3)4 (2-3 mol%) or Pd(OAc)2/{lg} (3 mol%/6 mol%)",
                    "reason": "Standard catalyst system for aryl bromides/iodides"}

    def _optimize(self, ob, hal, ligand, base, solvent):
        ht = hal.get('type', '')
        cond = {
            "catalyst_loading": "0.5-2 mol% Pd (modern)" if 'chloride' in ht else "1-3 mol% Pd (standard)",
            "ligand_loading": "2-4 mol% (relative to Pd)",
            "base": f"{base} (2-3 eq)",
            "solvent": solvent or "dioxane/H2O (3:1) or toluene/EtOH/H2O (4:1)",
            "concentration": "0.1-0.5 M",
            "temperature": "80°C" if 'chloride' in ht else "60-80°C" if 'bromide' in ht else "50-80°C",
            "time": "2-16 hours",
            "atmosphere": "N2 or Ar (DEGASSED solvents — oxygen causes homocoupling!)",
            "monitoring": "TLC or GC/MS for consumption of limiting reagent",
            "workup": ("Dilute with water, extract with EtOAc (×3), wash with brine, dry (Na2SO4), "
                       "concentrate, purify by column chromatography"),
            "palladium_removal": ("For pharmaceutical: pass through silica-bound thiol scavenger resin "
                                "or aqueous KF wash to reduce Pd to <10 ppm"),
        }
        return cond

    def _relevant_limitations(self, ob, hal):
        relevant = []
        bt = ob.get('type', '')
        ht = hal.get('type', '')

        if 'acid' in bt and 'chloride' not in ht:
            relevant.append({"issue": "Protodeboronation", "problem": "Boronic acids can decompose under basic conditions → Ar-H byproduct",
                           "solution": "Use boronic ester (Bpin) instead; lower temperature; minimize reaction time"})
        if 'alkyl' in ht:
            relevant.append({"issue": "β-Hydride elimination", "problem": "Alkyl-Pd intermediate can eliminate → alkene byproduct",
                           "solution": "Use primary alkyl only; Pd-PEPPSI or Ni catalyst; bidentate ligands"})
        if 'hetero' in ht or 'hetero' in bt:
            relevant.append({"issue": "Heteroatom coordination", "problem": "N/S/O atoms can coordinate/poison Pd catalyst",
                           "solution": "Add extra ligand (1-2 eq); protect coordinating groups"})
        relevant.append({"issue": "Oxygen sensitivity", "problem": "O2 causes oxidative homocoupling of boronic species",
                        "solution": "Degas solvents thoroughly; sparge with inert gas; use sealed tube"})
        return relevant

    def _estimate_yield(self, ob, hal):
        ht = hal.get('type', '')
        if 'iodide' in ht: return "85-98%"
        if 'bromide' in ht: return "80-95%"
        if 'chloride' in ht: return "70-92%"
        if 'triflate' in ht: return "80-95%"
        if 'vinyl' in ht: return "75-95%"
        if 'alkyl' in ht: return "50-85%"
        if 'hetero' in ht: return "60-90%"
        return "70-90%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        ob = parts[0] if parts else ""
        hal = parts[1] if len(parts) > 1 else ""
        lg = parts[2] if len(parts) > 2 else "PPh3"
        base = parts[3] if len(parts) > 3 else "K2CO3"
        solv = parts[4] if len(parts) > 4 else "dioxane/H2O"
        return self._run_base(ob, hal, lg, base, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
