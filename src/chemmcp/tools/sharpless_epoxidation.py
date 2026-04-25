"""
Sharpless Asymmetric Epoxidation (Tool #166)
Sharpless 不对称环氧化反应：Ti(OiPr)4/酒石酸酯/t-BuOOH 对烯丙醇的不对称环氧化。
涵盖：对映选择性规则（L/D-酒石酸酯决定绝对构型）、Ti-酒石酸络合物催化机理、
底物要求（游离烯丙醇OH必需）、预测模型、范围和限制。
"""
import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

_SHARPLESS_EPOX_DATA = {
    'enantioselectivity_rules': {
        'summary': 'The ABSOLUTE CONFIGURATION of the epoxide product is predicted by the DOUBLE BOND GEOMETRY of the allylic alcohol and the TARTRATE ESTER used:',
        'E_allylic_alcohol': {
            '(+)-DET (L-(+)-diethyl tartrate)': 'L-Epoxide (also called (2S,3S)-epoxy alcohol for standard E-allylic alcohols)',
            '(−)-DET (D-(−)-diethyl tartrate)': 'D-Epoxide ((2R,3R)-epoxy alcohol)',
            'memory_aid': 'E-allylic alcohol with (+)-DET gives L-product; think of it as E-plus → L',
        },
        'Z_allylic_alcohol': {
            '(+)-DET (L-(+)-diethyl tartrate)': 'D-Epoxide ((2R,3S)-epoxy alcohol for Z-substrates)',
            '(−)-DET (D-(−)-diethyl tartrate)': 'L-Epoxide ((2S,3R)-epoxy alcohol)',
            'memory_aid': 'Z-allylic alcohol REVERSES the rule compared to E',
        },
        'mnemonic': (
            "When the OH group is drawn to the RIGHT in the standard Fischer-like projection of the allylic alcohol, "
            "(+)-DET delivers oxygen from the BOTTOM face → predictable absolute configuration. "
            "This is one of the most reliable asymmetric reactions in organic chemistry — errors are exceedingly rare."
        ),
    },
    'mechanism': [
        ('1', 'Ti-tartrate complex formation', 'Ti(OiPr)4 + tartrate ester (DET or DIT) → chiral Ti-dikis(tartrate) dimer. This is the ACTIVE CATALYST. Stoichiometry: 1 Ti : 2 tartrate typically.'),
        ('2', 'Allylic alcohol coordination', 'The free OH of the allylic alcohol coordinates to Ti center — this is why FREE ALLYLIC OH IS ESSENTIAL! No coordination = no enantioselectivity.'),
        ('3', 'TBHP coordination & activation', 't-BuOOH coordinates to Ti and is activated as oxidant. The Ti-peroxo species transfers oxygen to alkene.'),
        ('4', 'Oxygen transfer (enantioselective)', '[3+2]-like oxygen delivery from Ti-peroxo to the Si or Re face of the coordinated alkene. Tartrate chirality determines which face is attacked.'),
        ('5', 'Product release', 'Epoxide product dissociates; catalyst regenerates with fresh TBHP.'),
    ],
    'catalyst_system': {
        'Ti(OiPr)4': {'role': 'Lewis acid / metal center', 'loading': '5-10 mol%', 'note': 'Source of Ti(IV); moisture sensitive'},
        'tartrate_ester': {'options': ['(+)-DET (diethyl L-tartrate)', '(−)-DET (diethyl D-tartrate)', '(+)-DIT (diisopropyl L-tartrate)', '(−)-DIT (diisopropyl D-tartrate)'],
                         'loading': '6-12 mol% (typically 2 eq per Ti)', 'role': 'Chiral ligand — determines face selectivity'},
        't-BuOOH (TBHP)': {'form': 'Usually as solution in water or decane', 'equivalents': '1.2-2.0 eq', 'role': 'Terminal oxidant (oxygen atom donor)'},
        'molecular_sieves': {'type': '3Å or 4Å molecular sieves', 'purpose': 'Remove water (reaction is water-sensitive); improve ee'},
    },
    'substrate_requirements': [
        ('Free allylic alcohol OH', 'ESSENTIAL — must have -CH=CH-CH2OH or similar motif with free OH', 'Protected (ether, ester) allylic alcohols DO NOT work'),
        ('Double bond position', 'Must be allylic (C=C adjacent to CH-OH) or homoallylic (one carbon further)', 'Homoallylic alcohols work but with lower ee'),
        ('Double bond geometry', 'E or Z geometry determines which enantiomer is formed', 'Critical for prediction!'),
        ('Substitution pattern', 'Mono-, di-, tri-, tetrasubstituted double bonds all work', 'More substituted = slower but often higher ee'),
    ],
    'scope': [
        'E-allylic primary alcohols: EXCELLENT — >90% ee typical, high yields',
        'Z-allylic primary alcohols: EXCELLENT — >90% ee typical (opposite configuration to E)',
        'E-allylic secondary alcohols: GOOD — works but may need optimization',
        'Trisubstituted allylic alcohols: GOOD — high stereoselectivity',
        'Tetrasubstituted allylic alchols: POSSIBLE but slow',
        'Homoallylic alcohols: WORKS — lower ee than true allylic (70-85%)',
        '2,3-Epoxy alcohols are VERSATILE intermediates: can be opened regioselectively → 1,2- or 1,3-diols, amino alcohols, etc.',
    ],
    'limitations': [
        ('Requires free allylic OH', 'Ether-protected or ester-protected allylic alcohols do NOT coordinate to Ti → no reaction or racemic', 'Solution: Must have free OH; protect other OH groups differently if needed'),
        ('Limited to allylic alcohols', 'Non-allylic alkenes cannot be epoxidized enantioselectively by this method', 'Solution: Use Jacobsen/Katsuki (Mn-salen) or Shi (fructose-derived ketone) for unfunctionalized alkenes'),
        ('Tartrate esters expensive', 'DET/DIT are costly (~$50-200/g depending on source/purity)', 'Solution: Can be recovered/recycled; use catalytic versions with molecular sieves'),
        ('TBHP hazards', 'tert-Butyl hydroperoxide is a strong oxidizer; shock-sensitive when concentrated; handle with care', 'Solution: Use commercial solutions (5-6 M in decane or 70% aq.); avoid concentrating'),
        ('Moisture sensitivity', 'Water degrades Ti-tartrate complex → lower activity and ee', 'Solution: Use activated 3Å/4Å MS; anhydrous solvents; inert atmosphere'),
        ('Low temperature required', 'Typically −20°C to −40°C for best ee', 'Solution: Dry ice/acetonitrile slurry or cryocooler'),
        ('Non-allylic alkenes in molecule', 'Other C=C bonds NOT adjacent to OH are not epoxidized selectively', 'Solution: They may react slowly — monitor selectivity'),
        ('Scale considerations', 'Large-scale use of TBHP requires safety review', 'Solution: Consider flow chemistry or alternative methods for scale-up'),
    ],
    'typical_outcomes': {
        'ee_values': '90-99% ee for standard E/Z-allylic primary alcohols',
        'yields': '70-95% (depending on substrate)',
        'cis_epoxide_selectivity': 'Generally forms cis-epoxides from acyclic substrates',
    },
}


@ChemMCPManager.register_tool
class SharplessEpoxidation(BaseTool):
    __version__ = "0.1.0"
    name = "SharplessEpoxidation"
    func_name = 'analyze_sharpless_epoxidation'
    description = "Sharpless asymmetric epoxidation analysis: Ti(OiPr)4/tartrate ester/TBHP-catalyzed enantioselective epoxidation of allylic alcohols. Covers enantioselectivity rules (E/Z geometry + DET stereochemistry → absolute config, with mnemonic), mechanism (5 steps: Ti-tartrate formation, allylic alcohol coordination, TBHP activation, enantioface-selective O-transfer), substrate requirements (free allylic OH essential), catalyst system details (Ti(OiPr)4, DET/DIT, TBHP, molecular sieves), scope (7 categories), limitations (8 items), and typical outcomes (90-99% ee)."
    implementation_description = "Comprehensive knowledge base covering: complete enantioselectivity prediction rules (E vs Z × (+)/(−)-DET matrix with mnemonic), 5-step mechanistic cycle emphasizing free OH requirement, detailed catalyst/reagent table, 4 substrate requirements, 7 scope categories, 8 limitations with solutions, and outcome benchmarks."
    categories = ["Reaction"]
    tags = ["Sharpless", "Asymmetric Epoxidation", "Titanium", "Tartrate", "Enantioselectivity", "Allylic Alcohol", "Catalysis"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("allylic_alcohol_smiles", "str", "N/A", "SMILES or name of the allylic alcohol substrate (must have free OH)."),
        ("double_bond_geometry", "str", "E", "Double bond geometry: 'E' or 'Z'."),
        ("tartrate_ester", "str", "(+)-DET", "Tartrate ester: '(+)-DET', '(-)-DET', '(+)-DIT', or '(-)-DIT'."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: allylic_alcohol [geometry] [tartrate]. E.g., 'geraniol E (+)-DET'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing substrate_analysis, predicted_configuration, mechanism_steps, catalyst_details, enantioselectivity_rules, optimal_conditions, scope, limitations, typical_outcomes, and recommendations."),
    ]

    examples = [
        {
            "code_input": {"allylic_alcohol_smiles": "geraniol", "double_bond_geometry": "E", "tartrate_ester": "(+)-DET"},
            "text_input": {"query": "geraniol E (+)-DET"},
            "output": {"result": {
                "reaction": "Sharpless epoxidation: geraniol (E-allylic alcohol) + (+)-DET → (2S,3S)-epoxygeraniol",
                "substrate_analysis": {"type": "E-allylic primary alcohol (primary OH)", "geometry": "E"},
                "predicted_product": "(2S,3S)-epoxygeraniol (L-epoxide)",
                "predicted_ee": "92-98%",
                "catalyst": "Ti(OiPr)4 (5-10 mol%) + (+)-DET (6-12 mol%) + TBHP (1.5 eq)",
                "conditions": {"T": "−20°C", "solvent": "CH2Cl2 (dry)", "additives": "3Å/4Å MS", "time": "8-24 h"},
                "yield": "80-92%",
                "key_rule": "E-allylic alcohol + (+)-DET → L-epoxide (2S,3S)",
            }},
        },
        {
            "code_input": {"allylic_alcohol_smiles": "cis-2-hexen-1-ol", "double_bond_geometry": "Z", "tartrate_ester": "(−)-DET"},
            "text_input": {"query": "cis-2-hexen-1-ol Z (-)-DET"},
            "output": {"result": {
                "reaction": "Sharpless epoxidation: (Z)-hex-2-en-1-ol + (−)-DET → (2S,3R)-epoxyhexanol (L-epoxide for Z)",
                "substrate_analysis": {"type": "Z-allylic primary alcohol", "geometry": "Z"},
                "predicted_product": "(2S,3R)-epoxyhexanol",
                "predicted_ee": "90-97%",
                "key_rule": "Z-allylic alcohol + (−)-DET → L-epoxide (rule reverses vs E)",
                "yield": "78-93%",
            }},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.data = dict(_SHARPLESS_EPOX_DATA)

    def _run_base(self, allylic_alcohol_smiles: str, double_bond_geometry: str = "E", tartrate_ester: str = "(+)-DET") -> dict:
        if not allylic_alcohol_smiles:
            raise ChemMCPInputError("Allylic alcohol substrate is required.")

        sub = self._analyze_substrate(allylic_alcohol_smiles, double_bond_geometry)
        config = self._predict_config(sub, tartrate_ester)
        cond = self._optimize(sub, tartrate_ester)
        limits = self._relevant_limitations(sub)

        result = {
            "result": {
                "reaction": f"Sharpless epoxidation: {allylic_alcohol_smiles} ({double_bond_geometry}-allylic alcohol) + {tartrate_ester} → {config.get('product','?')}",
                "substrate_analysis": sub,
                "predicted_configuration": config,
                "mechanism_steps": [{"step": s[0], "name": s[1], "detail": s[2]} for s in self.data['mechanism']],
                "catalyst_details": self.data['catalyst_system'],
                "enantioselectivity_rules": self.data['enantioselectivity_rules'],
                "optimal_conditions": cond,
                "scope": self.data['scope'],
                "applicable_limitations": limits,
                "typical_outcomes": self.data['typical_outcomes'],
                "summary": f"Sharpless epoxidation: {config.get('product','?')}. Predicted ee: {self.data['typical_outcomes']['ee_values']}. Key: FREE ALLYLIC OH REQUIRED!",
            }
        }
        logger.info(f"Sharpless epoxidation: {allylic_alcohol_smiles}")
        return result

    def _analyze_substrate(self, smi, geometry):
        s = (smi or "").strip().lower()
        geo = (geometry or "E").strip().upper()
        sub_types = [
            ('E-allylic primary alcohol', ['geraniol', r'e.*allylic', r'trans-2-alken-1-ol', r'(e)-'], 'Primary OH on allylic position; E double bond'),
            ('Z-allylic primary alcohol', ['cis-2-hexen-1-ol', r'z.*allylic', r'cis-2-alken-1-ol', r'(z)-'], 'Primary OH on allylic position; Z double bond'),
            ('homoallylic alcohol', ['homoallylic', r'pent-4-en-1-ol', r'but-3-en-1-ol'], 'OH one carbon removed from allylic position (lower ee)'),
            ('secondary allylic alcohol', ['1-propen-2-ol', r'secondary.*allylic'], 'Secondary OH at allylic position'),
            ('trisubstituted allylic alcohol', ['α-methyl geraniol', r'trisubstituted.*allylic'], 'More substituted double bond'),
        ]
        for stype, pats, note in sub_types:
            for pat in pats:
                if pat in s or re.search(pat, s):
                    return {"type": stype, "input": smi, "geometry": geo, "note": note}
        return {"type": "unknown_allylic_alcohol", "input": smi, "geometry": geo}

    def _predict_config(self, sub, tartrate):
        geo = sub.get('geometry', 'E')
        te = (tartrate or "(+)-DET").strip().lower()
        stype = sub.get('type', '')

        is_e = 'e' in geo.lower() or 'e-allylic' in stype.lower() or 'trans' in stype.lower()
        is_pos = '+det' in te or 'positive' in te or '(+' in (tartrate or '')

        if is_e and is_pos:
            return {"product": "L-Epoxide (2S,3S-configuration for standard substrates)", "rule": "E-allylic + (+)-DET → L-epoxide"}
        elif is_e and not is_pos:
            return {"product": "D-Epoxide (2R,3R-configuration for standard substrates)", "rule": "E-allylic + (−)-DET → D-epoxide"}
        elif not is_e and is_pos:
            return {"product": "D-Epoxide (2R,3S-configuration for Z-substrates)", "rule": "Z-allylic + (+)-DET → D-epoxide (REVERSED vs E)"}
        elif not is_e and not is_pos:
            return {"product": "L-Epoxide (2S,3R-configuration for Z-substrates)", "rule": "Z-allylic + (−)-DET → L-epoxide (REVERSED vs E)"}
        return {"product": "Configuration depends on geometry and tartrate choice", "rule": "Apply Sharpless rules"}

    def _optimize(self, sub, tartrate):
        return {
            "catalyst": f"Ti(OiPr)4 (5-10 mol%) + {tartrate or '(+)-DET'} (6-12 mol%)",
            "oxidant": "TBHP (1.2-2.0 eq; as 5-6 M solution in decane or 70% aq.)",
            "solvent": "CH2Cl2 (dry, distilled from CaH2) or toluene",
            "additives": "Activated 3Å or 4Å molecular sieves (~50 g/mmol) — ESSENTIAL for high ee",
            "temperature": "−20°C to −40°C (dry ice/MeCN slurry works well)",
            "atmosphere": "N2 or Ar (anhydrous, oxygen-free)",
            "concentration": "0.05-0.2 M",
            "time": "8-24 hours (monitor by TLC)",
            "workup": ("Quench with sat. aq. Na2SO3/NaHCO3 (cautiously!), filter off MS, extract with EtOAc, "
                       "wash with brine, dry (Na2SO4), concentrate, purify by column chromatography"),
        }

    def _relevant_limitations(self, sub):
        relevant = []
        stype = sub.get('type', '')
        relevant.append({"issue": "Free allylic OH required", "problem": "Protected allylic alcohols (ethers, esters) don't coordinate to Ti → no enantioselectivity",
                        "Solution": "Ensure OH is free; protect other functional groups differently"})
        if 'homoallylic' in stype.lower():
            relevant.append({"issue": "Lower ee expected", "problem": "Homoallylic alcohols give reduced enantioselectivity (70-85% ee)",
                           "Solution": "Accept lower ee or consider alternative method"})
        relevant.append({"issue": "TBHP safety", "problem": "tert-Butyl hydroperoxide is a strong oxidizer; hazardous at scale",
                        "Solution": "Use dilute solutions; proper PPE; avoid concentration/evaporation"})
        return relevant

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if parts else ""
        geo = parts[1] if len(parts) > 1 else "E"
        te = parts[2] if len(parts) > 2 else "(+)-DET"
        return self._run_base(sub, geo, te)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
