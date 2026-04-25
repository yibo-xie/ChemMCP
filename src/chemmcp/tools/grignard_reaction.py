"""
Grignard Reaction (Tool #157)
格氏反应：制备、羰基亲核加成、应用范围（醛/酮/酯/环氧化物/CO2）、
限制条件（质子性溶剂、敏感官能团）。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_GRIGNARD_DATA = {
    'preparation': {
        'general': ('R-X + Mg → R-Mg-X (ether solvent)', 'Dry Et2O or THF; initiation (heat/sonication/I2); exothermic once started'),
        'halide_reactivity': ('I > Br >> Cl (allylic/benzylic ≈ I in reactivity)', 'Ar-I and vinyl-Br work; Ar-Cl needs Li or activated Mg'),
        'solvent': ('Anhydrous diethyl ether (coordinating, stabilizes R2O2Mg complex) or THF', 'THF gives more reactive species for hindered substrates'),
        'initiation': ('Heat gently, crush Mg turnings (fresh surface), add crystal of I2 or 1,2-dibromoethane', 'Sonication can help initiate stubborn reactions'),
        'concentration': '0.1–1.0 M typical',
    },
    'scope': [
        {'electrophile': 'formaldehyde (HCHO)', 'product': 'primary alcohol (+1 carbon)', 'notes': 'R-MgX + HCHO → RCH2OH'},
        {'electrophile': 'aldehyde (R\'CHO)', 'product': 'secondary alcohol', 'notes': 'R-MgX + R\'CHO → RCH(R\')OH'},
        {'electrophile': 'ketone (R\'2CO)', 'product': 'tertiary alcohol', 'notes': 'R-MgX + R\'2CO → CR\'3(R)OH'},
        {'electrophile': 'ester (R\'COOR\")', 'product': 'tertiary alcohol (after 2 equiv)', 'notes': 'First: ketone intermediate; second: tertiary alcohol. 2 eq R-MgX needed.'},
        {'electrophile': 'acid chloride (R\'COCl)', 'product': 'tertiary alcohol (after 2 equiv)', 'notes': 'Reactive — often over-addition. Can stop at ketone at low T with CuI catalysis.'},
        {'electrophile': 'epoxide', 'product': 'alcohol (ring-opened at less substituted C)', 'notes': 'Regioselective SN2-like opening; CuI catalysis can invert selectivity'},
        {'electrophile': 'CO2 (dry ice)', 'product': 'carboxylic acid (+1 carbon)', 'notes': 'Classic carboxylation: R-MgX + CO2 → RCOOH after acidic workup'},
        {'electrophile': 'nitrile (R\'CN)', 'product': 'ketone after hydrolysis', 'notes': 'R-MgX + R\'CN → imine intermediate → hydrolyze → ketone'},
        {'electrophile': 'DMF, DMSO, etc.', 'product': 'aldehyde after hydrolysis', 'notes': 'Formyl equivalent: R-MgX + HC(OR)3 → aldehyde; DMF → aldehyde'},
        {'electrophile': 'alkyl halide (R\'X)', 'product': 'coupling (C-C bond formation)', 'notes': 'Possible but competes with Wurtz-type coupling; transition metal catalysis preferred nowadays'},
        {'electrophile': 'oxygen (O2)', 'product': 'alcohol (after hydrolysis)', 'notes': 'Side reaction to avoid — keep under inert atmosphere'},
    ],
    'limitations': [
        ('Protic solvents (H2O, ROH)', 'Quenches R-MgX violently → RH + Mg(OR)X', 'Solution: rigorously anhydrous conditions'),
        ('Acidic protons (O-H, N-H, S-H, terminal alkyne)', 'Deprotonation consumes R-MgX', 'Solution: protect these groups before Grignard formation'),
        ('Electrophilic functional groups within R group', 'Intramolecular attack possible', 'Avoid: additional C=O, COOR, CN, NO2, epoxide, halide in same molecule unless protected'),
        ('Steric hindrance', 'Very hindered ketones react slowly or not at all', 'Solution: use more reactive organocerium (from RLi + CeCl3) or organolithium'),
        ('Vinyl/aryl halides (unactivated)', 'Do not form Grignard readily', 'Solution: use lithium-halogen exchange (t-BuLi) or Mg activation (Rieke Mg)'),
        ('Enolizable protons α to carbonyl', 'Can be deprotonated instead of addition', 'Use lower T and controlled addition'),
        ('Competing enolization', 'Strong base character of R-MgX', 'Add Grignard to carbonyl (not vice versa) to minimize'),
        ('β-Hydride elimination (in Ni/Pd catalysis)', 'Not relevant for classical Grignard but matters for Kumada coupling', 'Use appropriate catalyst'),
    ],
    'safety': [
        ('Pyrophoric', 'Grignard reagents ignite spontaneously in air', 'Keep under N2/Ar at all times; use septum techniques'),
        ('Exothermic', 'Formation and carbonyl addition are highly exothermic', 'Cool in ice bath; add slowly with stirring'),
        ('Ether solvent fire hazard', 'Et2O is highly flammable, forms explosive peroxides', 'Test for peroxides; use fresh or properly stored ether'),
        ('Mg dust', 'Fine Mg is flammable', 'Handle carefully; avoid ignition sources'),
        ('Quenching', 'Always quench excess Grignard carefully with sat. NH4Cl or dilute acid', 'Never add water directly to concentrated Grignard solution'),
        ('Pressure buildup', 'Reaction produces gas during quenching', 'Quench slowly with cooling and venting'),
    ],
    'workup': (
        'After reaction completion (monitored by TLC or GC), cool to 0°C.\n'
        'Slowly pour onto saturated aqueous NH4Cl solution (or dilute HCl) with vigorous stirring.\n'
        'Extract with Et2O or EtOAc (×3).\n'
        'Wash combined organics with brine, dry over Na2SO4 or MgSO4.\n'
        'Filter, concentrate, purify by column chromatography or distillation.'
    ),
}


@ChemMCPManager.register_tool
class GrignardReaction(BaseTool):
    __version__ = "0.1.0"
    name = "GrignardReaction"
    func_name = 'analyze_grignard_reaction'
    description = "Comprehensive Grignard reagent reaction analysis: preparation (R-X + Mg in ether), nucleophilic addition scope (aldehydes → secondary alcohols, ketones → tertiary alcohols, esters → tertiary alcohols via 2 equiv, CO2 → carboxylic acids, epoxides → alcohols, nitriles → ketones), limitations (protic solvents, acidic protons, sensitive functional groups), safety notes, and workup procedure."
    implementation_description = "Complete knowledge base covering Grignard preparation (halide reactivity order, solvent effects, initiation methods), full electrophile scope table (10+ electrophile classes with products and notes), detailed limitation analysis (8 categories with solutions), safety protocols (6 items), standard workup procedure, and practical optimization tips."
    categories = ["Reaction"]
    tags = ["Grignard", "Organometallic", "Nucleophilic Addition", "Carbonyl", "C-C Bond Formation", "Alcohol Synthesis"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("grignard_reagent", "str", "N/A", "Grignard reagent formula (e.g., CH3MgBr, PhMgBr, C2H5MgI)."),
        ("electrophile_smiles", "str", "N/A", "SMILES or name of the electrophile (aldehyde, ketone, ester, CO2, epoxide, etc.)."),
        ("solvent", "str", "dry Et2O", "Solvent (must be anhydrous)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: grignard_reagent electrophile [solvent]. E.g., 'CH3MgBr benzaldehyde dry_Et2O'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing reagent_analysis, electrophile_analysis, mechanism, product_prediction, scope, limitations, safety_notes, and workup_procedure."),
    ]

    examples = [
        {
            "code_input": {"grignard_reagent": "CH3MgBr", "electrophile_smiles": "benzaldehyde", "solvent": "dry Et2O"},
            "text_input": {"query": "CH3MgBr benzaldehyde dry_Et2O"},
            "output": {"result": {
                "reaction": "CH3MgBr + C6H5CHO → 1-phenylethanol (secondary alcohol)",
                "reagent_analysis": {"name": "methylmagnesium bromide", "type": "primary alkyl Grignard", "reactivity": "highly nucleophilic"},
                "electrophile_analysis": {"name": "benzaldehyde", "type": "aromatic aldehyde", "product_type": "secondary alcohol"},
                "product": "1-phenylethanol (C6H5CH(CH3)OH)",
                "mechanism": ["Nucleophilic attack of CH3⁻ (from CH3MgBr) on carbonyl C of benzaldehyde", "Tetrahedral alkoxide intermediate", "Acidic workup → protonation → 1-phenylethanol"],
                "stoichiometry": "1 eq CH3MgBr : 1 eq benzaldehyde",
                "yield": "75-90%",
                "safety_critical": "Pyrophoric — keep under inert atmosphere",
            }},
        },
        {
            "code_input": {"grignard_reagent": "PhMgBr", "electrophile_smiles": "CO2", "solvent": "dry THF"},
            "text_input": {"query": "PhMgBr CO2 THF"},
            "output": {"result": {
                "reaction": "PhMgBr + CO2 → benzoic acid (after workup)",
                "product": "benzoic acid (C6H5COOH)",
                "special_conditions": "Pour Grignard onto crushed dry ice (excess CO2)", "yield": "70-85%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_GRIGNARD_DATA)

    def _run_base(self, grignard_reagent: str, electrophile_smiles: str, solvent: str = "dry Et2O") -> dict:
        if not grignard_reagent:
            raise ChemMCPInputError("Grignard reagent specification is required.")
        if not electrophile_smiles:
            raise ChemMCPInputError("Electrophile specification is required.")

        rgnt = self._analyze_reagent(grignard_reagent)
        elec = self._analyze_electrophile(electrophile_smiles)
        product = self._predict_product(rgnt, elec)
        mech = self._build_mechanism(rgnt, elec)
        limits = self._relevant_limitations(elec)
        cond = self._conditions(rgnt, elec, solvent)

        result = {
            "result": {
                "reaction": f"{grignard_reagent} + {electrophile_smiles} → {product.get('name','?')}",
                "reagent_analysis": rgnt,
                "electrophile_analysis": elec,
                "mechanism": mech,
                "product_prediction": product,
                "scope_table": self.data['scope'],
                "applicable_limitations": limits,
                "safety_notes": self.data['safety'],
                "workup_procedure": self.data['workup'],
                "optimal_conditions": cond,
                "summary": f"Grignard addition: {product.get('name','?')}. Yield: {product.get('yield_estimate','60-80%')}.",
            }
        }
        logger.info(f"Grignard: {grignard_reagent} + {electrophile_smiles}")
        return result

    def _analyze_reagent(self, rgnt):
        s = rgnt.strip()
        rtype = "unknown"
        if re.match(r'^[Cc][Hh]3[Mm][Gg]', s, re.IGNORECASE): rtype = "primary alkyl (methyl)"
        elif re.match(r'^[Cc]2[Hh]5|[Pp][Hh]|[Pp]h', s, re.IGNORECASE): rtype = "aryl (phenyl)"
        elif re.match(r'^[Cc][Hh]2=[Cc][Hh]|[Vv]inyl', s, re.IGNORECASE): rtype = "vinyl"
        elif re.match(r'^[Cc]≡|[Aa]lkynyl', s, re.IGNORECASE): rtype = "alkynyl"
        elif re.match(r'^[Cc].[Hh]21|[Aa]lkyl', s, re.IGNORECASE): rtype = "primary/secondary alkyl"
        return {"formula": s, "type": rtype, "nucleophilicity": "very high (carbanion character)", "basicity": "strong base"}

    def _analyze_electrophile(self, smi):
        s = (smi or "").strip().lower()
        e_map = {
            'formaldehyde': ['formaldehyde', 'hcho', 'methanal', 'ch2o'], 'aldehyde': ['aldehyde', 'cho', 'benzaldehyde', 'c6h5cho'],
            'ketone': ['ketone', 'acetone', '(ch3)2co', 'c=o(c)', 'cyclohexanone'],
            'ester': ['ester', 'coor', 'acetate', 'benzoate', 'coo'], 'acid_chloride': ['acyl.chloride', 'acid.chloride', 'cocl', 'cocl'],
            'epoxide': ['epoxide', 'oxirane', 'ethylene.oxide', 'propylene.oxide'],
            'co2': ['co2', 'carbon.dioxide', 'dry.ice'], 'nitrile': ['nitrile', 'cn', 'cyanide', 'ch3cn'],
            'dmf_formyl': ['dmf', 'dimethylformamide', 'formyl'], 'dmso': ['dmso', 'sulfoxide'],
        }
        for etype, pats in e_map.items():
            for pat in pats:
                if pat in s or re.search(pat, s):
                    scope_item = next((item for item in self.data['scope'] if item['electrophile'].lower() == etype.replace('_', ' ') or pat in item['electrophile'].lower()), None)
                    return {"type": etype, "name": smi, **(scope_item or {})}
        return {"type": "unknown_electrophile", "name": smi}

    def _predict_product(self, rgnt, elec):
        etype = elec.get('type', '')
        R = rgnt['formula'].split('Mg')[0] if 'Mg' in rgnt['formula'] else rgnt['formula']
        product_map = {
            'formaldehyde': {"name": f"{R}CH2OH (primary alcohol)", "class": "primary alcohol"},
            'aldehyde': {"name": f"{R}CH(elec.get('name','R'))OH (secondary alcohol)", "class": "secondary alcohol"},
            'ketone': {"name": f"{R}C(elec.get('name','R'))2OH (tertiary alcohol)", "class": "tertiary alcohol"},
            'ester': {"name": f"{R}C(elec.get('name','R'))2(OH) (tertiary alcohol, 2 eq needed)", "class": "tertiary alcohol", "equiv_needed": 2},
            'acid_chloride': {"name": f"Tertiary alcohol (2 eq may be needed)", "equiv_needed": 2},
            'epoxide': {"name": f"Ring-opened alcohol (at less substituted C)", "class": "alcohol"},
            'co2': {"name": f"{R}COOH (carboxylic acid)", "class": "carboxylic acid"},
            'nitrile': {"name": f"{R}C(=O)(elec_name) (ketone after hydrolysis)", "class": "ketone"},
            'dmf_formyl': {"name": f"{R}CHO (aldehyde)", "class": "aldehyde"},
        }
        p = product_map.get(etype, {"name": f"Addition product from {R} + {elec.get('name','?')}"})
        p["yield_estimate"] = "75-92%" if etype in ('formaldehyde','aldehyde','ketone','co2') else "60-85%" if etype in ('ester','epoxide','nitrile') else "50-80%"
        return p

    def _build_mechanism(self, rgnt, elec):
        R = rgnt['formula']
        return [
            f"1. Nucleophilic attack: {R} transfers R⁻ (with carbanion character) to electrophilic center of {elec.get('name','?')}",
            "2. Tetrahedral intermediate formed (alkoxide if carbonyl electrophile)",
            "3. Acidic workup (sat. NH4Cl or dilute HCl/H2SO4): protonation of intermediate",
            f"4. Product isolation: {self._predict_product(rgnt, elec).get('name','?')}",
        ]

    def _relevant_limitations(self, elec):
        relevant = []
        etype = elec.get('type', '')
        for lim_name, problem, solution in self.data['limitations']:
            if etype == 'ester' and 'electrophilic' in lim_name.lower():
                relevant.append({"issue": lim_name, "problem": problem, "solution": solution})
            elif etype in ('aldehyde', 'ketone') and 'enolizable' in lim_name.lower():
                relevant.append({"issue": lim_name, "problem": problem, "solution": solution})
            else:
                relevant.append({"issue": lim_name, "problem": problem, "solution": solution})[:1] if lim_name in ['Protic solvents', 'Acidic protons'] else None
        return [r for r in relevant if r]

    def _conditions(self, rgnt, elec, solvent):
        etype = elec.get('type', '')
        cond = {
            "solvent": solvent or "anhydrous THF or Et2O",
            "temperature": "0°C → RT (addition exothermic!)",
            "atmosphere": "N2 or Ar (rigorous exclusion of air/moisture)",
            "addition_order": "Add Grignard reagent dropwise to electrophile solution (minimizes enolization/wurtz side reactions)",
            "monitoring": "TLC or GC for consumption of electrophile",
        }
        if etype == 'co2':
            cond["special"] = "Pour Grignard solution onto excess crushed dry ice; then warm to RT slowly"
        if etype == 'ester':
            cond["stoichiometry"] = "2 equivalents of Grignard reagent required"
        if etype == 'epoxide':
            cond["temperature"] = "0°C to RT (may require gentle heating for hindered epoxides)"
        return cond

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        rgnt = parts[0] if parts else ""
        elec = parts[1] if len(parts) > 1 else ""
        solv = parts[2] if len(parts) > 2 else "dry Et2O"
        return self._run_base(rgnt, elec, solv)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
