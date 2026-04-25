"""
Sonogashira Coupling (Tool #162)
Sonogashira 偶联反应：钯/铜共催化末端炔烃与芳基/乙烯基卤化物的偶联反应。
涵盖：催化循环（Pd氧化加成、铜乙炔化合物形成/转金属、还原消除）、
铜助催化剂作用、胺碱、催化剂体系、炔烃保护/脱保护、范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_SONOGASHIRA_DATA = {
    'catalytic_cycle': [
        ('1', 'Oxidative addition (Pd)', 'Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X). Rate: RI > RBr >> RCl.'),
        ('2', 'Copper acetylide formation (Cu)', 'Terminal alkyne deprotonated by amine base → copper(I) acetylide forms via transmetallation with CuI. Cu acts as a TRANSMETALLATION SHUTTLE.'),
        ('3', 'Transmetalation (Cu→Pd)', 'Copper acetylide transfers alkynyl group to Pd(II) center, displacing X. This is the rate-accelerating role of Cu.'),
        ('4', 'Reductive elimination', 'Pd(II)(R¹)(C≡CR²) eliminates R¹-C≡CR² coupled product + regenerates Pd(0).'),
    ],
    'copper_role': {
        'purpose': 'Accelerates transmetalation step dramatically; without Cu the reaction is much slower (Castro-Stephens conditions)',
        'mechanism': 'Forms copper(I) acetylide intermediate which transfers to Pd more readily than direct alkyne deprotonation',
        'copper_free_variants': ('Possible with special ligands (Pd-PEPPSI, Buchwald ligands) but slower; '
                                 'Cu is still recommended for most applications'),
    },
    'catalyst_systems': {
        'Pd(PPh3)2Cl2/CuI': {'type': 'Classic Sonogashira system', 'loading': 'Pd 2-5 mol%, Cu 2-10 mol%', 'scope': 'Most widely used; good for standard substrates'},
        'Pd(PPh3)4/CuI': {'type': 'Pd(0) pre-catalyst + CuI', 'loading': 'Pd 1-3 mol%, Cu 2-5 mol%', 'scope': 'Very common alternative'},
        'PdCl2(PPh3)2/CuI': {'type': 'Pd(II) pre-catalyst + CuI', 'loading': 'Pd 2-5 mol%, Cu 5 mol%', 'scope': 'Reduced in situ to Pd(0); widely used'},
        'Pd(dppf)Cl2/CuI': {'type': 'Bidentate phosphine + CuI', 'loading': 'Pd 1-3 mol%, Cu 2-5 mol%', 'scope': 'Good for challenging couplings; chelating effect'},
        'Pd(OAc)2/XPhos/CuI': {'type': 'Modern Buchwald ligand + CuI', 'loading': 'Pd 0.5-2 mol%, Cu 1-3 mol%', 'scope': 'STATE OF THE ART: activates Ar-Cl even at RT'},
        'Pd-NHC/CuI (PEPPSI)': {'type': 'NHC Pd catalyst + CuI', 'loading': 'Pd 0.5-1 mol%, Cu 1-2 mol%', 'scope': 'Highly active; air-stable complexes'},
        'Copper-only (Cadiot-Chodkiewicz)': {'type': 'No Pd — Cu only', 'scope': 'For coupling of two alkynes (not aryl halides)'},
    },
    'base_and_solvent': {
        'amines_as_base_solvent': [
            ('Et3N / THF or DMF', 'Most common; Et3N as base AND co-solvent', 'Standard academic protocol'),
            ('iPr2NEt (Hünig\'s base)', 'Stronger, non-nucleophilic base', 'Prevents competitive nucleophilic substitution on aryl halide'),
            ('Piperidine / DMF', 'Stronger base than Et3N', 'For less reactive substrates'),
            ('Morpholine', 'Moderate base', 'Alternative amine solvent'),
            ('Et3N alone (neat)', 'Solvent-free conditions', 'Industrial approach; large scale'),
        ],
        'other_bases': [
            ('K2CO3 / DMF-H2O', 'Inorganic base with water co-solvent', 'Green chemistry approach'),
            ('Cs2CO3', 'Strong, soluble in organic solvents', 'For challenging substrates'),
            ('K3PO4', 'Strong, non-nucleophilic', 'Good for base-sensitive substrates'),
        ],
    },
    'alkyne_types': [
        ('Terminal alkyne (RC≡CH)', 'Standard Sonogashira partner', 'Must have terminal H for deprotonation'),
        ('TMS-protected alkyne (RC≡CTMS)', 'Can be deprotected IN SITU with TBAF or K2CO3/MeOH', 'One-pot protection/coupling/deprotection'),
        ('TES/TIPS-protected alkyne', 'Bulkier silyl groups', 'More stable; may need stronger deprotection'),
        ('Trimethylsilylacetylene (TMSA)', 'Common building block', 'Gives terminal alkyne after deprotection'),
        ('Diaryl/ dialkyl internal alkyne', 'NOT suitable for standard Sonogashira', 'No acidic proton — cannot form copper acetylide'),
        ('Terminal conjugated alkyne', 'Phenylacetylene, propiolic acid derivatives', 'Highly reactive; electron-poor especially so'),
    ],
    'halide_partners': [
        ('Aryl iodide (Ar-I)', 'Most reactive; works under mild conditions', 'Expensive but reliable'),
        ('Aryl bromide (Ar-Br)', 'Standard substrate; good cost/reactivity balance', 'Most commonly used'),
        ('Aryl chloride (Ar-Cl)', 'Least reactive; needs modern ligands', 'Cheap and abundant; requires optimization'),
        ('Aryl triflate (Ar-OTf)', 'Very reactive; from phenols', 'Good alternative to iodides'),
        ('Vinyl halide/triflate', 'Stereochemistry RETAINED', 'Important for enyne synthesis'),
        ('Acid chloride / pseudohalide', 'Can be used as electrophile', 'Gives ynones after coupling'),
        ('Heteroaryl halide', 'Pyridine, thiophene, furan halides', 'May need optimized conditions'),
    ],
    'scope': [
        'Aryl-alkyne couplings: EXCELLENT — the flagship transformation',
        'Vinyl-alkyne (enyne) couplings: Excellent — stereochemistry retained',
        'Heteroaryl couplings: Good — pyridine, thiophene, furan all work',
        'Electron-poor alkynes (propiolic esters): Very fast, high yielding',
        'Electron-rich alkynes (dialkylaminoalkynes): Work well',
        'TMS-alkynes with in-situ deprotection: Convenient one-pot procedure',
        'Aqueous Sonogashira: Possible with surfactants or water-soluble ligands',
        'Solid-supported Sonogashira: For combinatorial chemistry',
    ],
    'limitations': [
        ('Glaser homocoupling', '2 RC≡CH → RC≡C-C≡CR (oxidative dimerization of alkyne)', 'Solution: Degas thoroughly; use CuI (not Cu(II)); control O2; add reducing agent'),
        ('Alkyne polymerization', 'Alkyne can oligomerize/polymerize under basic conditions', 'Solution: Use moderate temperature; proper stoichiometry; avoid excess base'),
        ('Over-reaction of product', 'Product alkyne can undergo SECOND coupling (diarylation)', 'Solution: Use excess alkyne (1.5-3 eq); monitor conversion'),
        ('Dehalogenation side reaction', 'Ar-X + Base → Ar-H (reductive dehalogenation)', 'Solution: Optimize catalyst loading; use appropriate base'),
        ('Steric hindrance', 'Ortho-substituted aryl halides or bulky alkynes react slowly', 'Solution: Bulky ligands (XPhos, RuPhos); higher T'),
        ('Copper removal', 'Copper residues difficult to remove from product', 'Solution: Pass through silica; aqueous NH4OH wash; metal scavenger resins'),
        ('Terminal alkyne handling', 'Some terminal alkynes are volatile/toxic/malodorous', 'Solution: Use TMS-protected version; work in fume hood'),
        ('Alkyne isomerization', 'Terminal alkyne can isomerize to internal alkyne under basic conditions', 'Solution: Lower temperature; milder base; shorter time'),
    ],
    'typical_yields': {
        'aryl_iodide_terminal_alkyne': '80-98%',
        'aryl_bromide_terminal_alkyne': '75-95%',
        'aryl_chloride_modern_ligand': '65-90%',
        'heteroaryl': '60-88%',
        'vinyl_halide': '72-93%',
        'tms_alkyne_one_pot': '70-90%',
    },
}


@ChemMCPManager.register_tool
class SonogashiraCoupling(BaseTool):
    __version__ = "0.1.0"
    name = "SonogashiraCoupling"
    func_name = 'analyze_sonogashira_coupling'
    description = "Sonogashira coupling analysis: Pd/Cu-co-catalyzed cross-coupling of terminal alkynes with aryl/vinyl halides. Covers full catalytic cycle (Pd oxidative addition, copper acetylide formation/transmetalation, reductive elimination), Cu co-catalyst role (transmetalation shuttle), 7 catalyst systems (classic Pd(PPh3)2Cl2/CuI to modern NHC), amine base/solvent effects (5 options), 6 alkyne types including TMS deprotection, 7 halide partner types, scope (8 categories), limitations (8 items with solutions), and typical yields."
    implementation_description = "Comprehensive knowledge base covering: 4-step catalytic cycle with detailed Cu role explanation, 7 catalyst systems, 5 amine base/solvent options plus 3 inorganic bases, 6 alkyne types (including TMS-protection strategies), 7 halide partner types, 8 scope categories, 8 limitations with practical solutions, and yield benchmarks."
    categories = ["Reaction"]
    tags = ["Sonogashira", "Cross-Coupling", "Palladium", "Copper", "Alkyne", "C-C Bond Formation", "Catalysis"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("halide_smiles", "str", "N/A", "SMILES or name of the organic halide (aryl/vinyl iodide, bromide, chloride, or triflate)."),
        ("alkyne_smiles", "str", "N/A", "SMILES or name of the terminal alkyne (or TMS-protected alkyne)."),
        ("ligand", "str", "PPh3", "Ligand for Pd catalyst."),
        ("base", "str", "Et3N", "Base (usually amine, also serves as solvent/co-solvent)."),
        ("solvent", "str", "THF/Et3N", "Solvent system (often contains amine base)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: halide alkyne [ligand] [base] [solvent]. E.g., 'iodobenzene phenylacetylene PPh3 Et3N THF'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing halide_analysis, alkyne_analysis, catalytic_cycle, mechanism_steps, catalyst_recommendation, optimal_conditions, scope, limitations, typical_yields, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"halide_smiles": "iodobenzene", "alkyne_smiles": "phenylacetylene", "ligand": "PPh3", "base": "Et3N", "solvent": "THF"},
            "text_input": {"query": "iodobenzene phenylacetylene PPh3 Et3N THF"},
            "output": {"result": {
                "reaction": "Sonogashira coupling: Ph-I + PhC≡CH → diphenylacetylene (tolane)",
                "halide_analysis": {"type": "aryl iodide", "reactivity": "excellent", "note": "Most reactive Sonogashira substrate"},
                "alkyne_analysis": {"type": "terminal alkyne (phenylacetylene)", "reactivity": "good"},
                "product": "diphenylacetylene (tolane)",
                "catalyst_recommendation": "Pd(PPh3)2Cl2 (2 mol%) + CuI (5 mol%)",
                "conditions": {"T": "RT-80°C", "time": "2-12 h", "atmosphere": "N2"},
                "yield": "88-98%",
                "key_advantages": "Mild conditions, high functional group tolerance, atom-economical",
            }},
        },
        {
            "code_input": {"halide_smiles": "4-bromobenzaldehyde", "alkyne_smiles": "TMS-acetylene", "ligand": "PPh3", "base": "iPr2NEt", "solvent": "DMF"},
            "text_input": {"query": "4-bromobenzaldehyde TMS-acetylene PPh3 iPr2NEt DMF"},
            "output": {"result": {
                "reaction": "Sonogashira coupling: 4-BrC6H4CHO + TMS-C≡CH → 4-ethynylbenzaldehyde (after TMS deprotection)",
                "alkyne_note": "TMS-protected alkyne — can be deprotected in situ with K2CO3/MeOH or after workup",
                "product": "4-ethynylbenzaldehyde (useful building block)",
                "yield": "75-92%",
                "deprotection": "Add K2CO3 in MeOH after coupling completes, or TBAF in THF",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_SONOGASHIRA_DATA)

    def _run_base(self, halide_smiles: str, alkyne_smiles: str, ligand: str = "PPh3", base: str = "Et3N", solvent: str = "THF/Et3N") -> dict:
        if not halide_smiles:
            raise ChemMCPInputError("Organic halide is required.")
        if not alkyne_smiles:
            raise ChemMCPInputError("Terminal alkyne is required.")

        hal = self._analyze_halide(halide_smiles)
        alk = self._analyze_alkyne(alkyne_smiles)
        cat_rec = self._recommend_catalyst(hal, ligand)
        cond = self._optimize(hal, alk, ligand, base, solvent)
        limits = self._relevant_limitations(hal, alk)

        result = {
            "result": {
                "reaction": f"Sonogashira coupling: {halide_smiles} + {alkyne_smiles} → coupled alkyne product",
                "halide_analysis": hal,
                "alkyne_analysis": alk,
                "catalytic_cycle": self.data['catalytic_cycle'],
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['catalytic_cycle']],
                "copper_role": self.data['copper_role'],
                "catalyst_recommendation": cat_rec,
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "typical_yields": self._estimate_yield(hal, alk),
                "summary": f"Sonogashira coupling predicted: {self._estimate_yield(hal, alk)}. Key concern: {limits[0]['issue'] if limits else 'standard'}.",
            }
        }
        logger.info(f"Sonogashira: {halide_smiles} + {alkyne_smiles}")
        return result

    def _analyze_halide(self, smi):
        s = (smi or "").strip().lower()
        hal_types = [
            ('aryl iodide', ['iodo', r'-i\b', r'iodo'], 'Most reactive; expensive'),
            ('aryl bromide', ['bromo', r'-br\b', r'bromo'], 'Standard substrate; good cost/reactivity'),
            ('aryl chloride', ['chloro', r'-cl\b', r'chloro'], 'Least reactive; needs modern ligands; cheapest'),
            ('aryl triflate', ['triflate', r'-otf', r'otf'], 'Very reactive; from phenols'),
            ('vinyl halide', ['vinyl.*br', r'ch=ch.*hal'], 'Stereochemistry retained'),
            ('heteroaryl halide', ['bromopyridine', 'bromothiophene', r'pyrid.*br', r'thioph.*br'], 'Heterocyclic; may need optimization'),
        ]
        for htype, pats, note in hal_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": htype, "input": smi, "note": note}
        return {"type": "unknown_halide", "input": smi}

    def _analyze_alkyne(self, smi):
        s = (smi or "").strip().lower()
        alk_types = [
            ('terminal alkyne (aromatic)', ['phenylacetylene', 'tolanyl', r'phc≡ch', r'phenylethynyl'], 'Aromatic terminal alkyne; very common Sonogashira partner'),
            ('terminal alkyne (aliphatic)', ['1-hexyne', '1-octyne', 'propargyl', r'hc≡cc'], 'Aliphatic terminal alkyne; works well'),
            ('TMS-protected alkyne', ['tms', 'trimethylsilyl', r'tms.c≡ch', r'-c≡ctms'], 'Silyl-protected; can be deprotected in situ'),
            ('conjugated alkyne (EWG)', ['propiolate', 'propionic', r'c≡cco', r'c≡ccoom'], 'Electron-poor; highly reactive'),
            ('terminal alkyne general', ['terminal.*alkyne', r'c≡ch', r'ethynyl'], 'Generic terminal alkyne'),
            ('terminal alkyne (heteroatom)', ['propargyl alcohol', 'propargylamine', r'hoch2c≡ch'], 'Contains heteroatom; may need protection'),
        ]
        for atype, pats, note in alk_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": atype, "input": smi, "note": note}
        return {"type": "unknown_alkyne", "input": smi}

    def _recommend_catalyst(self, hal, ligand):
        ht = hal.get('type', '')

        if 'chloride' in ht:
            return {"primary": "Pd(OAc)2 (1 mol%) + XPhos/SPhos (2 mol%) + CuI (2 mol%) or Pd-PEPPSI (0.5 mol%) + CuI (1 mol%)",
                    "reason": "Aryl chlorides require modern bulky biaryl phosphine or NHC ligands"}
        if 'triflate' in ht:
            return {"primary": "Pd(PPh3)2Cl2 (2 mol%) + CuI (5 mol%)",
                    "reason": "Triflates are very reactive; classic system suffices"}
        if 'vinyl' in ht:
            return {"primary": "Pd(dppf)Cl2 (2 mol%) + CuI (3 mol%)",
                    "reason": "Bidentate ligand helps retain vinyl stereochemistry"}

        lg = (ligand or "PPh3").strip().lower()
        if any(x in lg for x in ['xphos', 'sphos', 'ruphos']):
            return {"primary": f"Pd(OAc)2 (1 mol%) + {ligand} (2 mol%) + CuI (2 mol%)"}
        else:
            return {"primary": f"Pd(PPh3)2Cl2 (2-3 mol%) + CuI (5 mol%) or Pd(PPh3)4 (2 mol%) + CuI (3 mol%)",
                    "reason": "Classic Sonogashira system — robust and reliable"}

    def _optimize(self, hal, alk, ligand, base, solvent):
        ht = hal.get('type', '')
        has_tms = 'tms' in (alk.get('type') or '').lower()

        cond = {
            "catalyst_loading": "Pd 0.5-2 mol%, Cu 1-5 mol%" if 'chloride' in ht else "Pd 1-3 mol%, Cu 2-5 mol%",
            "ligand_loading": "2-4 mol% (relative to Pd)",
            "base": f"{base} (2-3 eq; often also serves as co-solvent)",
            "solvent": solvent or "THF/Et3N (2:1) or DMF/Et3N",
            "concentration": "0.05-0.2 M",
            "temperature": "RT-50°C" if 'iodide' in ht else "50-80°C" if 'bromide' in ht else "80-110°C",
            "time": "2-16 hours",
            "atmosphere": "N2 or Ar (DEGASSED — oxygen causes Glaser homocoupling!)",
            "copper_iodide": "CuI (2-5 mol%) — essential co-catalyst",
            "monitoring": "TLC or GC/MS for consumption of limiting reagent",
            "workup": ("Dilute with water, extract with EtOAc (×3), wash with brine, wash with aq. NH4OH (remove Cu), "
                       "dry (Na2SO4), concentrate, purify by column chromatography"),
        }
        if has_tms:
            cond["deprotection"] = ("After coupling: add K2CO3 (2 eq) in MeOH, stir 30 min at RT, "
                                    "or treat with TBAF (1.1 eq) in THF")
        return cond

    def _relevant_limitations(self, hal, alk):
        relevant = []
        relevant.append({"issue": "Glaser homocoupling", "problem": "Oxidative dimerization of terminal alkyne → diyne byproduct",
                        "solution": "Degas solvents thoroughly; exclude oxygen; use CuI not Cu(II); add reducing agent"})
        relevant.append({"issue": "Over-coupling", "problem": "Product alkyne can undergo second Sonogashira coupling",
                        "solution": "Use excess alkyne (1.5-3 eq); monitor carefully"})
        at = alk.get('type', '')
        if 'tms' in at.lower():
            relevant.append({"issue": "TMS deprotection", "problem": "TMS group must be removed to get free terminal alkyne",
                           "solution": "In situ deprotection with K2CO3/MeOH or post-coupling TBAF treatment"})
        ht = hal.get('type', '')
        if 'ortho' in (hal.get('input') or '').lower():
            relevant.append({"issue": "Steric hindrance", "problem": "Ortho-substituted substrates are slow",
                           "solution": "Bulky ligands; higher temperature"})
        return relevant

    def _estimate_yield(self, hal, alk):
        ht = hal.get('type', '')
        if 'iodide' in ht: return "85-98%"
        if 'bromide' in ht: return "75-95%"
        if 'chloride' in ht: return "65-90%"
        if 'triflate' in ht: return "78-94%"
        if 'vinyl' in ht: return "72-93%"
        if 'hetero' in ht: return "60-88%"
        return "70-92%"

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        hal = parts[0] if parts else ""
        alk = parts[1] if len(parts) > 1 else ""
        lg = parts[2] if len(parts) > 2 else "PPh3"
        base = parts[3] if len(parts) > 3 else "Et3N"
        solv = parts[4] if len(parts) > 4 else "THF/Et3N"
        return self._run_base(hal, alk, lg, base, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
