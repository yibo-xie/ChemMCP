"""
Diels-Alder Reaction (Tool #156)
Diels-Aler [4+2] 环加成反应：二烯/亲二烯体要求、立体化学(endo/exo)、
区域选择性、取代基效应、逆 Diels-Alder。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_DA_DATA = {
    'diene_requirements': {
        'cisoid_conformation': {'description': 'Diene must be in s-cis (or able to adopt it) for orbital overlap', 'critical': True},
        'planarity': {'description': 'All 4 π-atoms of diene must be coplanar (or nearly so)', 'critical': True},
        'electron_rich': {'description': 'Normal electron demand: diene is electron-rich (HOMO-controlled)', 'note': 'Inverse demand: electron-poor diene + electron-rich dienophile'},
        'substituent_effect': {'description': 'EDG on diene raises HOMO → faster reaction; EWG lowers HOMO → slower'},
    },
    'dienophile_requirements': {
        'π_system': {'description': 'Must have π bond (alkene, alkyne, C=O, C≡N, etc.)', 'critical': True},
        'normal_demand': {'description': 'Electron-deficient (EWG-substituted) for normal D-A', 'examples': ['maleic anhydride', 'maleimide', 'p-benzoquinone', 'acrylate', 'nitroethylene', 'acrolein']},
        'inverse_demand': {'description': 'Electron-rich for inverse-demand D-A', 'examples': ['vinyl ether', 'enamine', 'electron-rich alkene']},
        'steric': {'description': 'Less substituted alkenes generally more reactive (less steric hindrance at π bond)'},
    },
    'stereochemistry': {
        'endo_exo': {
            'description': 'Endo vs exo approach of dienophile relative to diene',
            'endo_rule': ('Alder endo rule: Endo product is kinetically favored (secondary orbital interactions '
                         'between substituents of dienophile and diene π system stabilize the endo TS)'),
            'exo_thermodynamic': 'Exo product is often thermodynamically more stable (less steric clash)',
            'typical_ratio': 'Endo:exo ranges from >20:1 to 1:1 depending on substituents and conditions',
            'cavity_template': 'Catalytic antibodies / chiral Lewis acids can bias endo:exo ratio',
        },
        'suprafacial': {
            'description': '[4s+2s] cycloaddition — all bonding occurs on same face of each π system',
            'stereospecificity': 'STEREOSPECIFIC: cis-dienophile → cis-substituted cyclohexene; trans-dienophile → trans-substituted',
            'retention': 'Configuration of both diene and dienophile is retained in product',
        },
        'asymmetric_D-A': {
            'chiral_auxiliaries': ['Oppolzer sultam', 'acyl oxazolidinone (Evans)', 'bornyl ester'],
            'chiral_catalysts': ['Jacobsen/Cr(III) complexes', 'MacMillan imidazolidinones', 'chiral Al/B complexes'],
            'ee_range': 'Up to >99% ee reported with optimized systems',
        },
    },
    'regioselectivity': {
        'ortho_meta_para': {
            'description': 'Unsymmetrical components → regioisomeric products possible',
            'rule': ('EWG on dienophile ends up closer to more substituted end of diene (ortho-like position). '
                    'Alders rule: "maximal accumulation of double bonds"'),
            'para_rule': '1-substituted diene + 1-substituted dienophile → ortho major (not para)',
        },
        'ipso_positioning': 'Substituted patterns follow frontier molecular orbital coefficient matching',
    },
    'substituent_effects_table': [
        ('diene EDG (+OMe, +Me, +NR2)', 'Rate ↑↑ (raises HOMO)', 'Normal demand faster'),
        ('diene EWG (+COMe, +CN, +NO2)', 'Rate ↓↓ (lowers HOMO)', 'May switch to inverse demand'),
        ('dienophile EWG (+COMe, +CN, +CHO, +NO2)', 'Rate ↑↑ (lowers LUMO)', 'Normal demand faster'),
        ('dienophile EDG (+OMe, +NR2)', 'Rate ↓ (raises LUMO)', 'Normal demand slower; good for inverse demand'),
        ('dienophile = alkyne', 'Rate slower than alkene', 'Product: cyclohexadiene (aromatizable)'),
        ('dienophile = C=O', 'Possible but slow', 'Product: dihydropyran (hetero-D-A)'),
        ('dienophile = C≡N', 'Moderate reactivity', 'Product: dihydropyridine precursor'),
    ],
    'retro_diels_alder': {
        'description': 'Reverse reaction: cyclohexene → diene + dienophile',
        'driving_force': ('Aromatization (loss of benzene from bicyclic adduct), '
                          'release of stable gas (N2, CO2, CO, SO2), '
                          'formation of very stable molecule (e.g., C=O, aromatic)'),
        'applications': ['Protection strategy (D-A adduct → retro-D-A later)', 'Synthesis of reactive intermediates', 'Polymer curing', 'Generation of o-xylylenes'],
        'temperature': 'Typically requires heat (>100°C) or flash vacuum pyrolysis',
    },
    'conditions': {
        'thermal': {'T': '25-150°C (depends on reactivity)', 'solvent': 'any inert solvent (toluene, PhMe, CH2Cl2, or neat)', 'time': 'minutes to days', 'catalyst': 'none required'},
        'high_pressure': {'P': '10-15 kbar', 'effect': 'Accelerates dramatically (negative ΔV‡); allows unreactive pairs'},
        'lewis_acid': {'catalyst': 'AlCl3, Et2AlCl, BF3·OEt2, SnCl4, TiCl4, Yb(OTf)3, Sc(OTf)3', 'effect': 'Lowers LUMO of dienophile → rate acceleration; can improve endo selectivity and ee (with chiral LA)'},
        'microwave': {'effect': 'Rapid heating; can reduce reaction time from days to minutes'},
        'water': {'note': 'Hydrophobic effect in water can accelerate D-A (hydrophobic packing)'},

    },
}


@ChemMCPManager.register_tool
class DielsAlderReaction(BaseTool):
    __version__ = "0.1.0"
    name = "DielsAlderReaction"
    func_name = 'analyze_diels_alder_reaction'
    description = "Diels-Alder [4+2] cycloaddition analysis: diene and dienophile requirements (s-cis conformation, orbital symmetry), stereochemistry (endo/exo, suprafacial stereospecificity), regioselectivity (Alder rule), substituent effects, retro-Diels-Alder, optimal conditions, and scope."
    implementation_description = "Comprehensive knowledge base covering: diene conformational requirements (s-cis, planarity), dienophile electronic demands (normal/inverse), endo rule (kinetic control) vs exo (thermodynamic), suprafacial stereospecificity, regioselectivity rules (ortho/meta/para), Lewis acid catalysis, high-pressure promotion, and retro-D-A applications."
    categories = ["Reaction"]
    tags = ["Diels-Alder", "Pericyclic", "Cycloaddition", "Stereochemistry", "Regioselectivity", "Endo-Exo"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("diene_smiles", "str", "N/A", "SMILES or name of the diene component."),
        ("dienophile_smiles", "str", "", "SMILES or name of the dienophile component."),
        ("temperature_c", "float", "25", "Reaction temperature in °C."),
        ("solvent", "str", "", "Solvent (leave empty for thermal/neat conditions)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: diene_smiles [dienophile_smiles] [temp_C] [solvent]. E.g., 'C=CC=C maleic_anhydride 80 toluene'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing diene_analysis, dienophile_analysis, stereochemistry, regiochemistry, product_prediction, conditions, feasibility, and endo_exo_ratio."),
    ]

    examples = [
        {
            "code_input": {"diene_smiles": "C=CC=C", "dienophile_smiles": "maleic anhydride", "temperature_c": 80, "solvent": ""},
            "text_input": {"query": "C=CC=C maleic_anhydride 80"},
            "output": {"result": {
                "reaction": "Butadiene + maleic anhydride → cis-5-norbornene-endo-2,3-dicarboxylic anhydride",
                "diene_analysis": {"name": "1,3-butadiene", "s_cis_available": True, "planar": True, "electron_character": "moderately electron-rich"},
                "dienophile_analysis": {"name": "maleic anhydride", "type": "strong electron-deficient alkene", "excellent_dienophile": True},
                "stereochemistry": {"endo_product": "major (kinetic)", "endo_ratio": ">20:1 (endo:exo)", "approach": "[4s+2s] suprafacial — stereospecific"},
                "product": {"name": "bicyclo[2.2.1]hept-5-ene-2,3-dicarboxylic anhydride", "bicyclic": True, "endobicyclic": True},
                "feasibility": "excellent — classic D-A pair, highly favorable",
                "optimal_conditions": {"T": "80-110°C or RT (very reactive pair)", "solvent": "toluene or neat", "time": "1-24 h"},
                "yield": "85-98%",
            }},
        },
        {
            "code_input": {"diene_smiles": "c1ccc(C=C)c1", "dienophile_smiles": "acrylate", "temperature_c": 25, "solvent": ""},
            "text_input": {"query": "c1ccc(C=C)c1 acrylate 25"},
            "output": {"result": {
                "reaction": "ortho-Quinodimethane + acrylate → dihydronaphthalene derivative",
                "diene_analysis": {"name": "o-xylylene (generated in situ)", "reactive": "very reactive (strained)"},
                "feasibility": "good — o-xylylenes are highly reactive dienes",
                "yield": "70-90%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_DA_DATA)

    def _run_base(self, diene_smiles: str, dienophile_smiles: str = "", temperature_c: float = 25, solvent: str = "") -> dict:
        if not diene_smiles:
            raise ChemMCPInputError("Diene SMILES/name is required.")

        diene = self._analyze_diene(diene_smiles)
        dienophile = self._analyze_dienophile(dienophile_smiles) if dienophile_smiles else None

        T = temperature_c
        solv = solvent or "neat/toluene"

        # Stereochemical analysis
        stereo = self._analyze_stereochemistry(diene, dienophile)

        # Regiochemistry
        regio = self._analyze_regio(diene, dienophile)

        # Product prediction
        product = self._predict_product(diene, dienophile)

        # Feasibility assessment
        feas = self._assess_feasibility(diene, dienophile, T)

        # Conditions optimization
        cond = self._optimize_conditions(diene, dienophile, T, solv)

        result = {
            "result": {
                "reaction": f"{diene.get('name', diene_smiles)} + {dienophile.get('name', dienophile_smiles) if dienophile else '?'} → D-A adduct",
                "diene_analysis": diene,
                "dienophile_analysis": dienophile,
                "stereochemistry": stereo,
                "regiochemistry": regio,
                "product_prediction": product,
                "feasibility_assessment": feas,
                "optimal_conditions": cond,
                "substituent_effects": self.data['substituent_effects_table'],
                "retro_da_notes": self.data['retro_diels_alder'] if product.get('bicyclic') else None,
                "summary": f"D-A reaction: {feas['rating']}. {cond.get('notes','')}",
            }
        }
        logger.info(f"DielsAlder: {feas['rating']}")
        return result

    def _analyze_diene(self, smi):
        s = (smi or "").strip().lower()
        dienes = {
            'butadiene': {'patterns': [r'^c=cc=c$', r'butadiene', r'1,3-butadiene'], 's_cis': True, 'planar': True, 'reactivity': 'moderate', 'name': '1,3-butadiene'},
            'isoprene': {'patterns': [r'CC=CC=C', r'isoprene'], 's_cis': True, 'reactivity': 'good (EDG +Me)', 'name': 'isoprene (2-methyl-1,3-butadiene)'},
            'cyclopentadiene': {'patterns': [r'C1=CC=CC1', r'cyclopentadiene'], 's_cis': True, 'planar': True, 'reactivity': 'very high (strained, locked s-cis)', 'name': 'cyclopentadiene'},
            'ortho_quinodimethane': {'patterns': [r'c1ccc.*c=c', r'o.xylylene', r'quinodimethane'], 's_cis': True, 'reactivity': 'extremely high (transient, generated in situ)', 'name': 'o-xylylene (o-quinodimethane)'},
            'danishefsky': {'patterns': [r'danishefsky', r'silyloxy.butadiene'], 's_cis': True, 'reactivity': 'good (electron-rich)', 'hetero_D_A': True, 'name': "Danishefsky's diene (1-methoxy-3-TMS-oxybutadiene)"},
            'furan': {'patterns': [r'c1ccoc1', r'furan'], 's_cis': True, 'reactivity': 'moderate (aromaticity loss penalized)', 'reversible': True, 'name': 'furan (heterodiene)'},
            'cyclohexadiene': {'patterns': [r'C1=CC=CCC1', r'1,3-cyclohexadiene'], 's_cis': True, 'locked_s_cis': True, 'name': '1,3-cyclohexadiene'},
            'alpha_pyrone': {'patterns': [r'pyrone', r'pyran.one'], 's_cis': True, 'acts_as_both': True, 'lactone': True, 'name': 'α-pyrone (can act as diene + dienophile)'},
            'general_acyclic': {'patterns': [r'C=C.*C=C', r'diene'], 's_cis': 'must check', 'reactivity': 'depends on substitution', 'name': smi},
        }
        for dtype, info in dienes.items():
            for pat in info.get('patterns', []):
                if re.search(pat, s, re.IGNORECASE):
                    return dict(info)
        return {'name': smi, 's_cis': 'unknown — verify conformation', 'reactivity': 'unknown', 'note': f"Unrecognized diene '{smi}'"}

    def _analyze_dienophile(self, smi):
        s = (smi or "").strip().lower()
        dienophiles = {
            'maleic_anhydride': {'patterns': [r'maleic.anhydride', r'O=C1OC(=O)C=C1'], 'ewg': 'two carbonyls + anhydride', 'reactivity': 'excellent', 'name': 'maleic anhydride'},
            'maleimide': {'patterns': [r'maleimide', r'n[so].*c=c'], 'ewg': 'carbonyl + imide', 'reactivity': 'excellent', 'name': 'N-substituted maleimide'},
            'p_benzoquinone': {'patterns': [r'benzoquinone', r'quinone'], 'ewg': 'two conjugated C=O', 'reactivity': 'very good', 'name': 'p-benzoquinone'},
            'acrylate': {'patterns': [r'acrylate', r'C=CC(=O)O'], 'ewg': 'ester carbonyl', 'reactivity': 'good', 'name': 'alkyl acrylate'},
            'acrolein': {'patterns': [r'acrolein', r'C=CC=O'], 'ewg': 'aldehyde', 'reactivity': 'good', 'name': 'acrolein (propenal)'},
            'nitroethylene': {'patterns': [r'nitroethylene', r'C=C[N+](=O)[O-]'], 'ewg': 'nitro (very strong)', 'reactivity': 'excellent', 'name': 'nitroethylene'},
            'acetylene': {'patterns': [r'C#C', r'acetylene'], 'type': 'alkyne', 'reactivity': 'slower than alkene', 'name': 'acetylene'},
            'ethylene': {'patterns': [r'C=C', r'ethylene'], 'ewg': 'none', 'reactivity': 'poor (requires high T/P)', 'name': 'ethylene'},
            'singlet_oxygen': {'patterns': [r'singlet.oxygen', r'o2.*singlet'], 'type': 'diatomic', 'reactivity': 'good with dienes', 'name': 'singlet oxygen (^1O2)'},
        }
        for dtype, info in dienophiles.items():
            for pat in info.get('patterns', []):
                if re.search(pat, s, re.IGNORECASE):
                    return dict(info)
        return {'name': smi, 'reactivity': 'unknown', 'note': f"Unrecognized dienophile '{smi}'"}

    def _analyze_stereochemistry(self, diene, dieno):
        return {
            "cycloaddition_mode": "[4s+2s] suprafacial (concerted)",
            "stereospecificity": "Yes — configuration of both components is preserved in product",
            "endo_product": "Kinetically favored (Alder endo rule; secondary orbital interactions)",
            "exo_product": "Thermodynamically favored (less steric strain)",
            "predicted_endo_exo_ratio": "Highly dependent on substituents; typically endo-major for cyclic dienophiles (anhydrides, quinones)",
            "chiral_induction_possible": "Yes — via chiral auxiliaries, chiral Lewis acids, or organocatalysis",
            "relative_config": "cis-dienophile → cis-substituted in cyclohexene ring; trans → trans",
        }

    def _analyze_regio(self, diene, dieno):
        unsym_diene = any(p in (diene.get('name','')).lower() for p in ['isoprene', 'substituted', 'unsym'])
        unsym_dieno = dieno and dieno.get('ewg') and 'substituted' in (dieno.get('name','')).lower()

        if unsym_diene or unsym_dieno:
            return {
                "regioselectivity": "expected (unsymmetrical components)",
                "rule": "Alder ortho rule: maximal accumulation of double bonds in TS (EWG of dienophile ends up near more substituted diene terminus)",
                "major_regioisomer": "ortho-like (1,2-disubstituted) rather than meta (1,3)",
                "possible_minor": "meta regioisomer may form as minor product",
            }
        return {"regioselectivity": "N/A (symmetrical components → single product)"}

    def _predict_product(self, diene, dieno):
        d_name = diene.get('name','?')
        dieno_name = dieno.get('name','?') if dieno else 'unknown'
        bicyclic = 'cyclopentadiene' in d_name.lower() or 'furan' in d_name.lower()

        return {
            "name": f"Cycloadduct from {d_name} + {dieno_name}",
            "ring_system": "bicyclic (norbornene-type)" if bicyclic else "cyclohexene derivative",
            "bicyclic": bicyclic,
            "endobicyclic": bicyclic,
            "new_bonds_formed": "2 σ bonds (C-C) + 1 new π bond (if dienophile was alkene) or 2 π bonds (if alkyne)",
            "degree_of_unsaturation_change": "-1 (two π bonds consumed, one new π formed for alkene dienophile)",
        }

    def _assess_feasibility(self, diene, dieno, T):
        score = 50
        factors = []

        d_react = diene.get('reactivity', '')
        if 'very high' in d_react or 'extremely' in d_react:
            score += 30; factors.append("Very reactive diene")
        elif 'good' in d_react:
            score += 15; factors.append("Reactive diene")
        elif 'moderate' in d_react:
            score += 5

        if dieno:
            dw_react = dieno.get('reactivity', '')
            if 'excellent' in dw_react:
                score += 30; factors.append("Excellent dienophile")
            elif 'very good' in dw_react:
                score += 20
            elif 'good' in dw_react:
                score += 10
            elif 'poor' in dw_react:
                score -= 20; factors.append("Poor dienophile — needs high T or Lewis acid")

        if T > 100:
            score += 10; factors.append("Elevated T helps overcome activation barrier")
        elif T < 0:
            score -= 10

        rating = "excellent" if score >= 80 else "very good" if score >= 65 else "good" if score >= 45 else "moderate" if score >= 25 else "poor — consider Lewis acid catalyst or higher T"
        return {"rating": rating, "score": score, "factors": factors}

    def _optimize_conditions(self, diene, dieno, T, solv):
        cond = {
            "temperature": f"{T}°C (adjust based on reactivity: -78°C to 150°C typical range)",
            "solvent": solv,
            "time": "minutes (reactive pairs) to days (unreactive pairs)",
            "atmosphere": "N2/Ar (air-sensitive dienes/dienophiles should be protected)",
            "workup": "Cool, concentrate, purify by column chromatography or recrystallization",
        }
        d_react = diene.get('reactivity', '')
        if 'poor' in (dieno or {}).get('reactivity', '') or 'poor' in d_react:
            cond["catalyst_option"] = "Add Lewis acid (Et2AlCl 10 mol%, Yb(OTf)3, Sc(OTf)3) — can accelerate 10^6-10^8 fold"
            cond["pressure_option"] = "High pressure (10 kbar) as alternative"
        if 'cyclopentadiene' in (diene.get('name','')).lower():
            cond["note"] = "Cyclopentadiene dimerizes reversibly; crack dimer before use (heat to 170°C then distill)"
        return cond

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        d = parts[0] if parts else ""
        dino = parts[1] if len(parts) > 1 else ""
        t = float(parts[2]) if len(parts) > 2 else 25
        s = parts[3] if len(parts) > 3 else ""
        return self._run_base(d, dino, t, s)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
