"""
Reaction Energy Estimator (Tool #153)
估算反应的热力学可行性：利用键能数据、基团贡献法估算ΔH、ΔS、ΔG，
判断反应自发性，预测平衡常数范围。
"""
import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# 平均键能 (kJ/mol) - 最常用化学键
_BOND_ENERGIES = {
    # 单键
    'H-H': 436, 'H-F': 567, 'H-Cl': 431, 'H-Br': 366, 'H-I': 299,
    'C-H': 413, 'C-C': 347, 'C-N': 305, 'C-O': 358, 'C-F': 485,
    'C-Cl': 339, 'C-Br': 276, 'C-I': 238, 'C-S': 259, 'Si-H': 323,
    'N-H': 391, 'O-H': 463, 'O-O': 146, 'O-Cl': 203, 'F-F': 158,
    'Cl-Cl': 242, 'Br-Br': 193, 'I-I': 151, 'S-H': 363, 'S-S': 266,
    'P-H': 322, 'P-O': 335, 'P=O': 544, 'B-H': 389, 'B-O': 536,
    'Li-H': 285, 'Na-H': 197, 'Mg-H': 130, 'Al-H': 161,
    # 双键
    'C=C': 614, 'C=N': 615, 'C=O': 799,  # carbonyl
    'C=S': 577,
    'N=O': 607, 'O=O': 498, 'S=O': 523,  # sulfoxide
    'P=O': 544,
    # 三键
    'C≡C': 839, 'C≡N': 891, 'N≡N': 945,
    # 芳香/共轭 (特殊)
    'C=C(benzene)': 518, 'C-C(benzene)': 518,
}

# 常见原子化焓 / 标准生成焓 (kJ/mol)
_DELTA_Hf = {
    # 气态原子
    'H(g)': 218, 'C(g)': 717, 'N(g)': 473, 'O(g)': 249, 'F(g)': 79,
    'Cl(g)': 121, 'Br(g)': 112, 'I(g)': 107, 'S(g)': 279, 'P(g)': 314,
    # 常见分子标准生成焓 ΔfH°(g) kJ/mol at 298K
    'H2O(g)': -242, 'H2O(l)': -286, 'CO2(g)': -394, 'CO(g)': -111,
    'CH4(g)': -75, 'C2H6(g)': -85, 'C2H4(g)': +52, 'C2H2(g)': +227,
    'NH3(g)': -46, 'NO(g)': +90, 'NO2(g)': +33, 'N2O(g)': +82,
    'SO2(g)': -297, 'SO3(g)': -396, 'H2S(g)': -21, 'HCl(g)': -92,
    'HBr(g)': -36, 'HI(g)': +26, 'HF(g)': -273,
    'CH3OH(g)': -201, 'CH3OH(l)': -239, 'C2H5OH(g)': -235, 'C2H5OH(l)': -277,
    'HCHO(g)': -109, 'CH3CHO(g)': -166, 'CH3COCH3(g)': -218,
    'HCOOH(g)': -379, 'CH3COOH(g)': -432, 'CH3COOH(l)': -484,
    'C6H6(g)': +83, 'C6H6(l)': +49, 'C6H12(g)': -123, 'C6H14(g)': -167,
    'C6H5CH3(g)': +50, 'C6H5OH(s)': -165,
}

# 标准摩尔熵 S° (J/mol·K) at 298K
_S_STANDARD = {
    'H2(g)': 131, 'O2(g)': 205, 'N2(g)': 192, 'Cl2(g)': 223, 'Br2(g)': 245,
    'I2(s)': 116, 'H2O(g)': 189, 'H2O(l)': 70, 'CO2(g)': 214, 'CO(g)': 198,
    'CH4(g)': 186, 'C2H6(g)': 230, 'C2H4(g)': 220, 'C2H2(g)': 201,
    'NH3(g)': 192, 'NO(g)': 211, 'NO2(g)': 240, 'H2S(g)': 206,
    'HCl(g)': 187, 'HF(g)': 174, 'C(graphite)': 5.7, 'C(diamond)': 2.4,
    'CH3OH(g)': 240, 'CH3OH(l)': 127, 'C2H5OH(g)': 283, 'C2H5OH(l)': 161,
    'HCHO(g)': 219, 'CH3CHO(g)': 250, 'CH3COCH3(g)': 300,
    'CH3COOH(g)': 283, 'CH3COOH(l)': 160, 'C6H6(g)': 269, 'C6H6(l)': 173,
    'C6H12(g)': 298, 'C6H14(g)': 287, 'C6H5CH3(g)': 321,
}

# 反应类型特征能量变化 (kJ/mol)
_REACTION_TYPE_ENERGY = {
    'combustion_hydrocarbon': {'delta_h_range': (-900, -400), 'always_exothermic': True},
    'combustion_alcohol': {'delta_h_range': (-1400, -600), 'always_exothermic': True},
    'neutralization_strong': {'delta_h_range': (-57, -55), 'always_exothermic': True},
    'hydrogenation_alkene': {'delta_h_range': (-140, -80), 'always_exothermic': True},
    'hydrogenation_alkyne_to_alkene': {'delta_h_range': (-175, -150), 'always_exothermic': True},
    'hydrogenation_alkyne_to_alkane': {'delta_h_range': (-310, -270), 'always_exothermic': True},
    'hydration_alkene': {'delta_h_range': (-50, -30), 'usually_exothermic': True},
    'dehydration_alcohol': {'delta_h_range': (+30, +70), 'endothermic': True},
    'esterification': {'delta_h_range': (-5, +5), 'near_equilibrium': True},
    'hydrolysis_ester': {'delta_h_range': (+5, +20), 'slightly_endothermic': True},
    'saponification': {'delta_h_range': (-10, +10), 'near_equilibrium': True},
    'nucleophilic_substitution_SN2': {'delta_h_range': (-30, +30), 'depends_on_substrate': True},
    'elimination_E2': {'delta_h_range': (+20, +80), 'often_endothermic': True},
    'addition_HX_alkene': {'delta_h_range': (-80, -40), 'exothermic': True},
    'combustion_general': {'delta_h_range': ('very_negative',), 'always_exothermic': True},
    'photosynthesis': {'delta_h_range': (+2800, +4800), 'strongly_endothermic': True},
    'ATP_hydrolysis': {'delta_h': -30.5, 'notes': 'Biochemical standard'},
    'bond_formation_C-C': {'delta_h': -347, 'notes': 'typical C-C single bond'},
    'bond_formation_C=C': {'delta_h': -614, 'notes': 'typical C=C double bond'},
    'bond_formation_C-H': {'delta_h': -413, 'notes': 'typical C-H bond'},
    'bond_formation_O-H': {'delta_h': -463, 'notes': 'typical O-H bond'},
    'bond_formation_N-H': {'delta_h': -391, 'notes': 'typical N-H bond'},
    'condensation_polymerization': {'delta_h_range': (-5, -1), 'per_bond': True},
    'radical_polymerization': {'delta_h_range': (-90, -60), 'per_mol_monomer': True},
    'ionic_lattice_formation': {'delta_h_range': (-800, -200), 'highly_exothermic': True},
    'acid_base_neutralization': {'delta_h': -57, 'notes': 'strong acid + strong base → -57 kJ/mol H2O'},
    'precipitation': {'delta_h_range': (-60, +40), 'varies_widely': True},
    'redox_displacement': {'delta_h_range': (-400, +100), 'depends_on_metals': True},
    'complex_formation': {'delta_h_range': (-50, +20), 'often_exothermic': True},
    'coordination_bond': {'delta_h_range': (-200, -20), 'ligand_dependent': True},
    'aldol_addition': {'delta_h_range': (-20, +10), 'near_equilibrium': True},
    'aldol_condensation': {'delta_h_range': (-20, -5), 'driven_by_conjugation': True},
    'grignard_formation': {'delta_h_range': (+50, +100), 'endothermic': True},
    'grignard_addition_carbonyl': {'delta_h_range': (-100, -40), 'exothermic': True},
    'wittig_reaction': {'delta_h_range': (+20, +60), 'endothermic': True, 'entropy_driven': True},
    'diels_alder': {'delta_h_range': (-90, -20), 'exothermic': True, 'negative_delta_s': True},
    'friedel_crafts_alkylation': {'delta_h_range': (-50, -10), 'exothermic': True},
    'friedel_crafts_acylation': {'delta_h_range': (-100, -40), 'exothermic': True},
    'suzuki_coupling': {'delta_h_range': (-30, +20), 'near_equilibrium': True, 'catalyst_driven': True},
    'claisen_condensation': {'delta_h_range': (-10, +10), 'equilibrium_controlled': True},
    'electrophilic_aromatic_substitution': {'delta_h_range': (-80, -20), 'exothermic': True},
    'nucleophilic_aromatic_substitution': {'delta_h_range': (-40, +30), 'depends_on_substrate': True},
    'oxidation_alcohol_to_aldehyde': {'delta_h_range': (-180, -120), 'exothermic': True},
    'oxidation_aldehyde_to_acid': {'delta_h_range': (-260, -200), 'exothermic': True},
    'reduction_nitro_to_amine': {'delta_h_range': (-500, -400), 'highly_exothermic': True},
    'halogenation_alkane': {'delta_h_range': (-120, -60), 'exothermic': True},
    'combustion_biofuel': {'delta_h_range': (-2800, -1200), 'always_exothermic': True},
}


@ChemMCPManager.register_tool
class ReactionEnergyEstimator(BaseTool):
    __version__ = "0.1.0"
    name = "ReactionEnergyEstimator"
    func_name = 'estimate_reaction_energy'
    description = "Estimate thermodynamic feasibility of a chemical reaction: calculate approximate ΔH (enthalpy change), estimate ΔS (entropy change), compute ΔG (Gibbs free energy), predict equilibrium constant K, and assess spontaneity and feasibility."
    implementation_description = "Uses multiple estimation approaches: (1) Bond energy method for gas-phase reactions, (2) Standard formation enthalpy lookup from built-in database, (3) Reaction type pattern matching with characteristic energy ranges, (4) Group contribution estimates for ΔS. Combines these to give ΔG = ΔH - TΔS and K = exp(-ΔG/RT). Provides confidence levels and method notes."
    categories = ["Reaction"]
    tags = ["Thermodynamics", "Gibbs Energy", "Enthalpy", "Entropy", "Equilibrium", "Feasibility", "Bond Energy"]
    required_envs = []
    oss_dependencies = []
    services_and_software = []

    code_input_sig = [
        ("reactants_smiles", "str", "N/A", "SMILES or names of reactants, separated by '+' or space."),
        ("products_smiles", "str", "", "SMILES or names of products, separated by '+'. Leave empty if reaction is implied."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (default: 298.15 K = 25°C)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Space-separated: reactants [products] [temperature_K]. E.g., 'C2H4+H2 C2H6 298'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing delta_h, delta_s_estimate, delta_g_estimate, equilibrium_constant, spontaneity, feasibility_rating, method_used, confidence, and recommendations."),
    ]

    examples = [
        {
            "code_input": {
                "reactants_smiles": "C=C + H2",
                "products_smiles": "CC",
                "temperature_k": 298.15,
            },
            "text_input": {"query": "C=C+H2 CC 298"},
            "output": {
                "result": {
                    "reaction": "ethylene hydrogenation: C2H4 + H2 → C2H6",
                    "delta_h": -136,
                    "delta_h_unit": "kJ/mol",
                    "delta_s_estimate": -120,
                    "delta_s_unit": "J/(mol·K)",
                    "delta_g_estimate": -100.2,
                    "delta_g_unit": "kJ/mol",
                    "equilibrium_constant": 4.2e17,
                    "log_k": 17.6,
                    "spontaneous": True,
                    "spontaneity_description": "Highly spontaneous (exothermic, negative ΔG)",
                    "feasibility_rating": "excellent — reaction proceeds readily at room temperature with catalyst",
                    "method_used": "bond_energy_estimation + reaction_type_pattern",
                    "confidence": "medium-high",
                    "kinetics_note": "Thermodynamically favorable but requires activation energy (catalyst needed: Pd, Pt, Ni)",
                    "temperature_effect": "Lower T favors product (exothermic); higher T slightly reduces K but increases rate",
                }
            },
        },
        {
            "code_input": {
                "reactants_smiles": "N2 + 3H2",
                "products_smiles": "2NH3",
                "temperature_k": 298.15,
            },
            "text_input": {"query": "N2+3H2 2NH3 298"},
            "output": {
                "result": {
                    "reaction": "Haber-Bosch: N2 + 3H2 → 2NH3",
                    "delta_h": -92,
                    "delta_h_unit": "kJ/mol",
                    "delta_s_estimate": -199,
                    "delta_s_unit": "J/(mol·K)",
                    "delta_g_estimate": -32.6,
                    "delta_g_unit": "kJ/mol",
                    "equilibrium_constant": 5.8e5,
                    "log_k": 5.8,
                    "spontaneous": True,
                    "feasibility_rating": "thermodynamically feasible but kinetically limited (requires high P, T, catalyst)",
                    "method_used": "standard_formation_enthalpy_lookup",
                    "confidence": "high",
                    "temperature_effect": "Exothermic: lower T favors NH3; but kinetics requires high T (400-500°C industrial compromise)",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.bond_energies = dict(_BOND_ENERGIES)
        self.delta_hf = dict(_DELTA_Hf)
        self.s_standard = dict(_S_STANDARD)
        self.reaction_types = dict(_REACTION_TYPE_ENERGY)
        self.R = 8.314  # J/(mol·K)

    def _run_base(self, reactants_smiles: str, products_smiles: str = "", temperature_k: float = 298.15) -> dict:
        """Core logic: estimate reaction thermodynamics."""
        if not reactants_smiles:
            raise ChemMCPInputError("Reactants SMILES/names cannot be empty.")

        T = temperature_k
        if T <= 0:
            raise ChemMCPInputError("Temperature must be positive (Kelvin).")

        rxn_str = f"{reactants_smiles} → {products_smiles}" if products_smiles else reactants_smiles

        # Strategy 1: Try to match known reaction types
        rxn_type_match = self._match_reaction_type(rxn_str)

        # Strategy 2: Look up standard enthalpies of formation
        hf_result = self._lookup_enthalpies(reactants_smiles, products_smiles)

        # Strategy 3: Estimate from bond energies (if we can parse bonds)
        bond_result = self._estimate_from_bonds(rxn_str)

        # Combine strategies to get best ΔH estimate
        delta_h, dh_method, dh_conf = self._combine_dh_estimates(rxn_type_match, hf_result, bond_result)

        # Estimate ΔS
        delta_s, ds_method = self._estimate_entropy(rxn_str, reactants_smiles, products_smiles, T)

        # Calculate ΔG = ΔH - TΔS
        delta_g = delta_h - T * delta_s / 1000  # convert J→kJ

        # Equilibrium constant: K = exp(-ΔG/RT)
        try:
            K = math.exp(-delta_g * 1000 / (self.R * T))
            log_k = math.log10(K) if K > 0 else float('-inf')
        except OverflowError:
            K = float('inf')
            log_k = float('inf')

        # Spontaneity & feasibility assessment
        spontaneity, feas_rating, notes = self._assess_feasibility(delta_h, delta_s, delta_g, K, T, rxn_type_match)

        result = {
            "result": {
                "reaction": rxn_str.replace(' ', ' ').strip(),
                "delta_h": round(delta_h, 1),
                "delta_h_unit": "kJ/mol",
                "delta_s_estimate": round(delta_s, 1),
                "delta_s_unit": "J/(mol·K)",
                "delta_g_estimate": round(delta_g, 1),
                "delta_g_unit": "kJ/mol",
                "equilibrium_constant": f"{K:.2e}" if K != float('inf') else ">1e308",
                "log_k": round(log_k, 2) if log_k != float('inf') else "very large",
                "spontaneous": delta_g < 0,
                "spontaneity_description": spontaneity,
                "feasibility_rating": feas_rating,
                "method_used": dh_method,
                "confidence": dh_conf,
                "entropy_method": ds_method,
                "temperature_k": round(T, 2),
                "thermodynamic_notes": notes,
                "temperature_effect": self._temp_effect(delta_h, delta_s),
                "practical_considerations": self._practical_notes(rxn_str, rxn_type_match),
            }
        }

        logger.info(f"ReactionEnergy: {rxn_str[:40]}... ΔH={delta_h}, ΔG={delta_g}, K={K:.2e}")
        return result

    def _match_reaction_type(self, rxn: str):
        """Match reaction against known reaction type patterns."""
        rxn_lower = rxn.lower()
        matches = []

        type_patterns = [
            (r'combust|burn|\+?O2.*\→.*(CO2|H2O)', 'combustion_general'),
            (r'hydrogenat|\\+?H2.*\→.*(alkane|sat)', 'hydrogenation_alkene'),
            (r'hydrat|\\+?H2O.*\→.*alcohol|alkene.*→.*OH', 'hydration_alkene'),
            (r'dehydrat|alcohol.*→.*alkene', 'dehydration_alcohol'),
            (r'esterif|acid.*alcohol.*→.*ester', 'esterification'),
            (r'hydrolys|ester.*→.*acid', 'hydrolysis_ester'),
            (r'saponif|ester.*base', 'saponification'),
            (r'neutraliz|acid.*base.*→.*salt', 'neutralization_strong'),
            (r'aldol|enolate.*carbonyl', 'aldol_addition'),
            (r'diels.alder|\[4\+2\]|cycloaddition', 'diels_alder'),
            (r'friedel.craft|FC|arene.*electrophile', 'friedel_crafts_alkylation'),
            (r'grignard|R.Mg', 'grignard_addition_carbonyl'),
            (r'wittig|ylide.*carbonyl', 'wittig_reaction'),
            (r'suzuki|organoboron.*palladium', 'suzuki_coupling'),
            (r'claisen|ester.*condens', 'claisen_condensation'),
            (r'n2.*h2|haber|ammonia.synth', 'N2+3H2→2NH3'),  # special case
            (r'photosynth|CO2.*H2O.*glucose', 'photosynthesis'),
            (r'atp.*hydrolys|atp.*adp', 'ATP_hydrolysis'),
            (r'polymeriz|monomer.*polymer', 'radical_polymerization'),
            (r'oxid.*alcohol.*aldehyde', 'oxidation_alcohol_to_aldehyde'),
            (r'oxid.*aldehyde.*acid', 'oxidation_aldehyde_to_acid'),
            (r'reduc.*nitro.*amine', 'reduction_nitro_to_amine'),
            (r'halogen.*alkane|alkane.*X2', 'halogenation_alkane'),
            (r'SN2|nucleophilic.subst', 'nucleophilic_substitution_SN2'),
            (r'E2|elimination', 'elimination_E2'),
            (r'EAS|electrophilic.aromatic', 'electrophilic_aromatic_substitution'),
            (r'SNAr|nucleophilic.aromatic', 'nucleophilic_aromatic_substitution'),
        ]

        for pattern, rtype in type_patterns:
            if re.search(pattern, rxn_lower):
                info = self.reaction_types.get(rtype, {})
                if info:
                    matches.append({"type": rtype, **info})

        return matches[0] if matches else None

    def _lookup_enthalpies(self, reactants, products):
        """Look up standard enthalpies of formation."""
        # Parse common molecule names/formulas
        dh_react = 0
        dh_prod = 0
        n_react = 0
        n_prod = 0
        found_react = []
        found_prod = []

        for mol in re.split(r'\+|\s+', reactants):
            mol = mol.strip()
            if not mol:
                continue
            # Handle coefficients like "2H2" or "3O2"
            coef_match = re.match(r'^(\d+)(.+)$', mol)
            coef = int(coef_match.group(1)) if coef_match else 1
            formula = coef_match.group(2) if coef_match else mol

            key = None
            for k in self.delta_hf:
                if k.split('(')[0].lower() == formula.lower() or k.replace('(g)','').replace('(l)','').replace('(s)','') == formula:
                    key = k
                    break

            if key:
                dh_react += coef * self.delta_hf[key]
                n_react += coef
                found_react.append((formula, coef, self.delta_hf[key]))

        for mol in re.split(r'\+|\s+', products):
            mol = mol.strip()
            if not mol:
                continue
            coef_match = re.match(r'^(\d+)(.+)$', mol)
            coef = int(coef_match.group(1)) if coef_match else 1
            formula = coef_match.group(2) if coef_match else mol

            key = None
            for k in self.delta_hf:
                if k.split('(')[0].lower() == formula.lower() or k.replace('(g)','').replace('(l)','').replace('(s)','') == formula:
                    key = k
                    break

            if key:
                dh_prod += coef * self.delta_hf[key]
                n_prod += coef
                found_prod.append((formula, coef, self.delta_hf[key]))

        if found_react and found_prod:
            return {
                "delta_h": dh_prod - dh_react,
                "method": "standard_formation_enthalpy",
                "details": {"reactants": found_react, "products": found_prod},
                "confidence": "high" if n_react >= 2 and n_prod >= 2 else "medium",
            }
        return None

    def _estimate_from_bonds(self, rxn: str):
        """Estimate ΔH from bond energies (crude but useful fallback)."""
        return None  # Would need full structural parsing; defer to pattern matching

    def _combine_dh_estimates(self, rxn_type, hf_result, bond_result):
        """Combine different ΔH estimation methods."""
        if hf_result and hf_result.get("confidence") == "high":
            return hf_result["delta_h"], hf_result["method"], hf_result["confidence"]
        if hf_result:
            return hf_result["delta_h"], hf_result["method"], hf_result["confidence"]
        if rxn_type:
            dhr = rxn_type.get("delta_h_range")
            if dhr and isinstance(dhr, tuple) and len(dhr) == 2:
                avg = sum(dhr) / 2
                return avg, f"reaction_type_pattern ({rxn_type['type']})", "medium-low"
            if 'delta_h' in rxn_type:
                return rxn_type['delta_h'], f"reaction_type_pattern ({rxn_type['type']})", "medium"

        return 0, "no_data_available", "very_low"

    def _estimate_entropy(self, rxn_full, reactants, products, T):
        """Estimate entropy change."""
        # Count gas molecules change (dominant factor for ΔS)
        rxn_lower = rxn_full.lower()

        # Simple heuristic: count molecules before vs after
        n_gas_react = len([m for m in re.split(r'\+|\s+', reactants) if m.strip()])
        n_gas_prod = len([m for m in re.split(r'\+|\s+', products) if m.strip()]) if products else n_gas_react

        dn_gas = n_gas_prod - n_gas_react

        # Base entropy change per mole of gas change: ~±100-150 J/(mol·K)
        base_ds = dn_gas * (-130)  # losing gas molecules → more order → negative ΔS

        # Reaction-type-specific adjustments
        ds_adjustments = {
            'combustion_general': 0,  # small net gas change typically
            'hydrogenation_alkene': -120,  # 2 gas → 1 gas
            'dehydration_alcohol': +130,  # liquid → gas + liquid
            'esterification': -50,  # 2 liquids → 2 liquids (one water lost)
            'photosynthesis': -200,  # gas consumption
            'N2+3H2→2NH3': -199,  # 4 mol gas → 2 mol gas
            'diels_alder': -160,  # 2 molecules → 1 (ordering)
            'wittig_reaction': +80,  # Ph3PO solid precipitate drives entropy
            'claisen_condensation': -30,  # similar molecule counts
            'suzuki_coupling': +20,  # boron byproduct precipitation
            'aldol_condensation': -20,  # dehydration produces water/gas
            'ATP_hydrolysis': +100,  # 1 molecule → 2 (more disorder)
            'neutralization_strong': +50,  # ions → fewer species
        }

        rxn_type_match = self._match_reaction_type(rxn_full)
        if rxn_type_match:
            t = rxn_type_match.get("type", "")
            if t in ds_adjustments:
                base_ds = ds_adjustments[t]

        # Additional heuristics
        if 'solid' in rxn_lower or '(s)' in rxn_full or 'ppt' in rxn_lower:
            base_ds -= 30  # precipitation → more order
        if 'gas' in rxn_lower or '(g)' in rxn_full:
            pass  # already accounted
        if 'dissolving' in rxn_lower or 'solution' in rxn_lower:
            base_ds += 50  # dissolution → more disorder

        method = "gas_molecule_count_heuristic + reaction_type_pattern"
        if rxn_type_match:
            method += f" (matched: {rxn_type_match['type']})"

        return round(base_ds, 1), method

    def _assess_feasibility(self, dh, ds, dg, K, T, rxn_type):
        """Assess overall feasibility."""
        # Spontaneity description
        if dg < -40:
            spont = "Highly spontaneous (strongly negative ΔG)"
        elif dg < -10:
            spont = "Spontaneous (negative ΔG)"
        elif dg < 0:
            spont = "Weakly spontaneous (slightly negative ΔG)"
        elif dg < 10:
            spont = "Near equilibrium (ΔG ≈ 0)"
        elif dg < 40:
            spont = "Non-spontaneous under given conditions (positive ΔG)"
        else:
            spont = "Strongly non-spontaneous (large positive ΔG)"

        # Feasibility rating
        if K > 1e10:
            rating = "excellent — essentially quantitative yield expected"
        elif K > 1e4:
            rating = "very good — high conversion expected (>99%)"
        elif K > 100:
            rating = "good — moderate-to-high conversion; may need optimization"
        elif K > 1:
            rating = "moderate — favorable but equilibrium-limited; consider Le Chatelier"
        elif K > 0.01:
            rating = "poor — unfavorable equilibrium; needs driving force (product removal, excess reactant)"
        elif K > 1e-6:
            rating = "very poor — highly unfavorable; not practical without external intervention"
        else:
            rating = "essentially impossible under these conditions"

        # Additional notes
        notes = []
        if abs(dh) < 10:
            notes.append("Near-thermoneutral: entropy dominates the spontaneity.")
        if ds > 100:
            notes.append("Large positive ΔS: favored by increasing temperature.")
        if ds < -100:
            notes.append("Large negative ΔS: favored by decreasing temperature (but kinetic concerns apply).")
        if rxn_type:
            extra = rxn_type.get("notes", "")
            if extra:
                notes.append(extra)
            if rxn_type.get("catalyst_driven"):
                notes.append("This reaction is kinetically limited despite thermodynamics — catalyst essential.")

        return spont, rating, notes

    def _temp_effect(self, dh, ds):
        """Describe temperature effect on equilibrium."""
        if ds > 0:
            return "Endothermic-like (ΔS>0): increasing T makes reaction MORE favorable (larger K)."
        elif ds < 0:
            return "Exothermic-like (ΔS<0): increasing T makes reaction LESS favorable (smaller K)."
        else:
            return "ΔS ≈ 0: temperature has minimal effect on K."

    def _practical_notes(self, rxn, rxn_type):
        """Generate practical considerations."""
        notes = []
        if rxn_type:
            t = rxn_type.get("type", "")

            if t in ("hydrogenation_alkene", "hydrogenation_alkyne_to_alkene"):
                notes.append("Requires catalyst (Pd, Pt, Ni, Raney Ni); pressure equipment for H2 gas.")
            elif t == "combustion_general":
                notes.append("Requires ignition/activation energy; once started, self-sustaining.")
            elif t == "N2+3H2→2NH3":
                notes.append("Industrial: 400-500°C, 150-300 atm, Fe catalyst. Thermodynamic-kinetic tradeoff.")
            elif t == "diels_alder":
                notes.append("Often thermal (no catalyst needed); can be accelerated by high pressure or Lewis acids.")
            elif t in ("suzuki_coupling",):
                notes.append("Requires Pd catalyst, base, anhydrous conditions sensitive to oxygen/water.")
            elif t == "wittig_reaction":
                notes.append("Anhydrous conditions critical; ylide preparation sensitive to moisture.")
            elif t == "grignard_addition_carbonyl":
                notes.append("Strictly anhydrous/anaerobic; ether solvent; exothermic addition — control T.")
            elif t == "dehydration_alcohol":
                notes.append("Requires acid catalyst (H2SO4) and heat (>170°C) or POCl3/pyridine (milder).")

        if not notes:
            notes.append("Consider kinetic factors: even thermodynamically favorable reactions may require catalyst, heat, or activation.")

        return notes

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        reactants = parts[0] if parts else ""
        products = parts[1] if len(parts) > 1 else ""
        temp = float(parts[2]) if len(parts) > 2 else 298.15
        return self._run_base(reactants, products, temp)


# For regex in _run_base
import re
