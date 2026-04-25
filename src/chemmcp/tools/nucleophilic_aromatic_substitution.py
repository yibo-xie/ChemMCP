"""
Nucleophilic Aromatic Substitution (Tool #123)
芳香亲核取代机理（SNAr 加成-消除、苯炔中间体）。
Provides SNAr mechanism analysis via Meisenheimer complex (addition-elimination)
or benzyne intermediate pathway.
"""
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


# EWG strength for SNAr activation (ortho/para to leaving group)
EWG_STRENGTH = {
    'NO2': ('very strong', 'nitro group — powerful -M and -I'),
    'CN':  ('strong', 'cyano group — -M and -I'),
    'COR': ('strong', 'acyl group — -M and -I'),
    'CHO': ('strong', 'formyl group — -M and -I'),
    'COOR':('moderate', 'ester/carboxyl — -M and -I'),
    'CF3': ('moderate', 'trifluoromethyl — strong -I'),
    'SO3R':('strong', 'sulfonyl — -M and -I'),
    'F':   ('weak activation but best LG', 'fluorine — strong inductive withdrawal + excellent LG'),
    'Cl':  ('weak LG', 'chlorine — moderate LG'),
    'Br':  ('good LG', 'bromine — good LG'),
    'I':   ('excellent LG', 'iodine — excellent LG but poor EWG'),
}

LEAVING_GROUP_QUALITY_SNAr = {'F': 5, 'Cl': 4, 'Br': 3, 'I': 2, 'OTs': 6, 'OSO2R': 6}


@ChemMCPManager.register_tool
class NucleophilicAromaticSubstitution(BaseTool):
    __version__ = "0.1.0"
    name = "NucleophilicAromaticSubstitution"
    func_name = 'explain_snar_mechanism'
    description = "Explain nucleophilic aromatic substitution (SNAr) mechanism: addition-elimination via Meisenheimer complex for activated aryl halides, or benzyne intermediate for unactivated substrates under extreme conditions. Covers ortho/para nitro activation, leaving group effects, and the amination of chlorobenzene."
    implementation_description = "Analyzes aromatic substrate for electron-withdrawing groups (especially nitro groups ortho/para to leaving group), determines whether SNAr (addition-elimination) or benzyne pathway operates, provides stepwise mechanism with Meisenheimer complex resonance stabilization, predicts regiochemistry, and evaluates reaction feasibility."
    categories = ["Reaction"]
    tags = ["SNAr", "Nucleophilic Aromatic Substitution", "Meisenheimer Complex", "Benzyne", "Nitro Activation"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('substrate_smiles', 'str', 'N/A', 'SMILES string of the aromatic substrate with leaving group.'),
        ('nucleophile', 'str', 'OH-', 'Nucleophile. Options: OH-, NH3, CN-, CH3O-, N3-, etc.'),
        ('mechanism_type', 'str', 'auto', 'Mechanism: auto (auto-detect), SNAr, or benzyne.'),
        ('temperature_c', 'float', '25.0', 'Temperature in °C (benzyne needs >300°C or very strong base).'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: substrate_smiles nucleophile [mechanism_type] [temperature]. E.g., "O=c1c(cc(cc1)[N+](=O)[O-])Cl OH- SNAr 50".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing mechanism_pathway, activation_analysis, mechanism_steps, meisenheimer_complex_analysis, product_prediction, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'substrate_smiles': 'O=c1c(cc(cc1)[N+](=O)[O-])Cl',
                'nucleophile': 'OH-',
                'mechanism_type': 'SNAr',
                'temperature_c': 50.0,
            },
            'text_input': {'query': 'O=c1c(cc(cc1)[N+](=O)[O-])Cl OH- SNAr'},
            'output': {
                'result': {
                    'substrate': 'm-chloronitrobenzene (activated)',
                    'mechanism_pathway': 'SNAr (addition-elimination) via Meisenheimer complex',
                    'activation': 'nitro group at meta position provides moderate activation',
                    'steps': [
                        {'step': 1, 'description': 'Nucleophilic addition: OH⁻ attacks C bearing Cl → Meisenheimer complex'},
                        {'step': 2, 'description': 'Elimination: Cl⁻ departs → restored aromaticity → m-nitrophenol'},
                    ],
                    'product': 'm-nitrophenol',
                    'favorability': 'good — activated by NO₂',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'c1ccccc1Cl',
                'nucleophile': 'NH2-',
                'mechanism_type': 'benzyne',
                'temperature_c': '350',
            },
            'text_input': {'query': 'c1ccccc1Cl NH2- benzyne 350'},
            'output': {
                'result': {
                    'substrate': 'chlorobenzene (unactivated)',
                    'mechanism_pathway': 'Benzyne intermediate (extreme conditions)',
                    'steps': [
                        {'step': 1, 'description': 'Strong base (NH₂⁻) abstracts ortho proton → benzyne + HCl'},
                        {'step': 2, 'description': 'Nucleophile adds to benzyne → anionic intermediate'},
                        {'step': 3, 'description': 'Protonation → substituted arene (mixture of ortho products)'},
                    ],
                    'regiochemistry': 'mixture of ortho-substituted products (no directing effect)',
                    'favorability': 'requires extreme conditions (350°C or NaNH₂/NH₃)',
                }
            },
        },
        {
            'code_input': {
                'substrate_smiles': 'O=c1c([N+](=O)[O-])c(cc(c1)[N+](=O)[O-])Cl',
                'nucleophile': 'NH3',
                'mechanism_type': 'auto',
                'temperature_c': 25.0,
            },
            'text_input': {'query': 'O=c1c([N+](=O)[O-])c(cc(c1)[N+](=O)[O-])Cl NH3'},
            'output': {
                'result': {
                    'substrate': '2,4-dinitrochlorobenzene (highly activated)',
                    'mechanism_pathway': 'SNAr — extremely facile',
                    'activation': 'two nitro groups (ortho + para) provide maximum activation',
                    'product': '2,4-dinitroaniline',
                    'favorability': 'excellent — reacts even with weak nucleophiles at RT',
                }
            },
        },
    ]

    def _run_base(self, substrate_smiles: str, nucleophile: str = 'OH-', mechanism_type: str = 'auto', temperature_c: float = 25.0) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(substrate_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Analyze substrate
        substrate_analysis = self._analyze_substrate(mol)

        # 2. Determine mechanism pathway
        pathway = self._determine_pathway(mechanism_type, substrate_analysis, temperature_c)

        # 3. Build mechanism steps
        steps = self._build_mechanism_steps(substrate_analysis, nucleophile, pathway)

        # 4. Meisenheimer/benzyne analysis
        intermediate_analysis = self._analyze_intermediate(substrate_analysis, pathway)

        # 5. Product prediction
        product = self._predict_product(substrate_analysis, nucleophile, pathway)

        # 6. Favorability
        favorability = self._evaluate_favorability(substrate_analysis, pathway, temperature_c)

        result = {
            'result': {
                'substrate_smiles': substrate_smiles,
                'substrate_analysis': substrate_analysis,
                'nucleophile': nucleophile,
                'mechanism_pathway': pathway['name'],
                'pathway_determination': pathway.get('reason', ''),
                'mechanism_steps': steps,
                'intermediate_analysis': intermediate_analysis,
                'product_prediction': product,
                'temperature': f'{temperature_c}°C',
                'favorability': favorability,
                'summary': self._build_summary(substrate_analysis, pathway, product, favorability),
            }
        }

        logger.info(f"SNAr: {substrate_smiles} + {nucleophile} → {pathway['name']} → {favorability}")
        return result

    def _analyze_substrate(self, mol):
        """Analyze aromatic substrate for SNAr suitability."""
        # Find aromatic ring(s)
        aromatic_atoms = [a for a in mol.GetAtoms() if a.GetIsAromatic()]
        if not aromatic_atoms:
            return {'has_aromatic': False}

        # Find leaving group attached to ring
        lg_info = None
        lg_carbon_idx = None

        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in ('F', 'Cl', 'Br', 'I'):
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIsAromatic():
                        lg_info = {
                            'atom_sym': sym,
                            'atom_idx': atom.GetIdx(),
                            'quality': LEAVING_GROUP_QUALITY_SNAr.get(sym, 1),
                            'name': f'{sym}ide',
                        }
                        lg_carbon_idx = neighbor.GetIdx()
                        break

        # Find EWGs on the ring (especially ortho/para to LG)
        ewgs = []
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and lg_carbon_idx:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIsAromatic() or neighbor.GetIdx() == lg_carbon_idx:
                        continue
                    nsym = neighbor.GetSymbol()

                    # Check for nitro groups
                    if nsym == 'N':
                        for nn in neighbor.GetNeighbors():
                            if nn.GetAtomicNum() == 8:
                                bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                                if bond and bond.GetBondTypeAsDouble() >= 2.0:
                                    pos = self._get_position_relative(atom.GetIdx(), lg_carbon_idx, len(aromatic_atoms))
                                    ewgs.append({'type': 'NO2', 'position': pos, 'strength': 'very strong'})

                    # Check for carbonyl-containing groups
                    elif nsym == 'C':
                        for nn in neighbor.GetNeighbors():
                            if nn.GetAtomicNum() == 8:
                                b = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                                if b and b.GetBondTypeAsDouble() >= 2.0:
                                    pos = self._get_position_relative(atom.GetIdx(), lg_carbon_idx, len(aromatic_atoms))
                                    ewgs.append({'type': 'COR', 'position': pos, 'strength': 'strong'})

                    # CF3
                    elif nsym == 'F' and neighbor.GetIdx() != (lg_info.get('atom_idx') if lg_info else -1):
                        n_f = sum(1 for nn in neighbor.GetNeighbors() if nn.GetSymbol() == 'F')
                        if n_f >= 2:
                            pos = self._get_position_relative(atom.GetIdx(), lg_carbon_idx, len(aromatic_atoms))
                            ewgs.append({'type': 'CF3', 'position': pos, 'strength': 'moderate'})

        # Count activating EWGs at ortho/para positions
        n_ortho_para_ewg = sum(1 for e in ewgs if e['position'] in ('ortho', 'para'))
        has_ortho_nitro = any(e['type'] == 'NO2' and e['position'] == 'ortho' for e in ewgs)
        has_para_nitro = any(e['type'] == 'NO2' and e['position'] == 'para' for e in ewgs)

        return {
            'has_aromatic': True,
            'leaving_group': lg_info,
            'lg_carbon_index': lg_carbon_idx,
            'ewgs': ewgs,
            'n_ewgs_total': len(ewgs),
            'n_ortho_para_ewg': n_ortho_para_ewg,
            'has_ortho_nitro': has_ortho_nitro,
            'has_para_nitro': has_para_nitro,
            'activation_level': self._classify_activation(n_ortho_para_ewg, has_ortho_nitro),
        }

    def _get_position_relative(self, atom_idx, lg_idx, n_aromatic_atoms):
        """Determine if substituent is ortho/meta/para to LG."""
        # Simplified: use distance in ring atoms
        diff = abs(atom_idx - lg_idx)
        if n_aromatic_atoms > 0:
            min_diff = min(diff, n_aromatic_atoms - diff)
            if min_diff <= 1: return 'ortho'
            elif min_diff <= 2: return 'meta'
            else: return 'para'
        return 'unknown'

    def _classify_activation(self, n_op_ewg, has_ortho_no2):
        """Classify overall SNAr activation level."""
        if has_ortho_no2 and n_op_ewg >= 2: return 'extremely high (reacts at RT)'
        if has_ortho_no2: return 'very high (reacts under mild heating)'
        if n_op_ewg >= 2: return 'high (needs moderate heat)'
        if n_op_ewg == 1: return 'moderate (needs forcing conditions)'
        return 'low/unactivated (benzyne pathway needed)'

    def _determine_pathway(self, requested, substrate_analysis, temp_c):
        """Determine which mechanism operates."""
        if requested != 'auto':
            name = 'SNAr (Addition-Elimination)' if requested == 'SNAr' else 'Benzyne Intermediate'
            return {'name': name, 'reason': f'User-specified: {requested}'}

        activation = substrate_analysis.get('activation_level', '')
        lg = substrate_analysis.get('leaving_group', {})

        if 'extremely' in activation or 'very high' in activation:
            return {'name': 'SNAr (Addition-Elimination)', 'reason': f'Strong activation ({activation})'}
        if 'high' in activation or 'moderate' in activation:
            return {'name': 'SNAr (Addition-Elimination)', 'reason': f'Moderate-strong activation ({activation})'}

        # Unactivated — check conditions for benzyne
        if temp_c >= 200 or 'NaNH' in str(temp_c) or temp_c >= 100:
            return {'name': 'Benzyne Intermediate', 'reason': 'Unactivated substrate under extreme conditions → benzyne pathway'}

        return {
            'name': 'SNAr unlikely / Benzyne required',
            'reason': f'Low activation ({activation}). Need stronger EWGs or extreme conditions (NaNH₂/NH₃, >200°C).'
        }

    def _build_mechanism_steps(self, substrate, nu, pathway):
        """Build stepwise mechanism."""
        path_name = pathway.get('name', '')
        nu_display = nu.replace('-', '⁻') if nu.endswith('-') else nu

        if 'SNAr' in path_name:
            lg = substrate.get('leaving_group', {})
            lg_name = lg.get('name', 'X')

            steps = [
                {
                    'step': 1,
                    'name': 'Nucleophilic Addition (Rate-Determining)',
                    'equation': f'Ar-X + Nu⁻ → Meisenheimer complex (σ-complex)',
                    'details': (
                        f"The nucleophile ({nu_display}) attacks the carbon bearing the leaving group "
                        f"({lg_name}). The aromatic π system is disrupted as a new C-Nu bond forms, "
                        f"producing a cyclohexadienyl anion intermediate (Meisenheimer complex). "
                        f"The negative charge is delocalized onto ortho/para positions, especially "
                        f"stabilized by electron-withdrawing groups (nitro groups). This step is slow "
                        f"and rate-determining."
                    ),
                    'rate_determining': True,
                },
                {
                    'step': 2,
                    'name': 'Elimination of Leaving Group (Fast)',
                    'equation': f'Meisenheimer complex → Ar-Nu + X⁻',
                    'details': (
                        f"The leaving group ({lg_name}) departs with its bonding electrons, "
                        f"the π system reforms, and aromaticity is restored. This fast step releases "
                        f"the substitution product."
                    ),
                    'rate_determining': False,
                },
            ]

        else:  # Benzyne
            steps = [
                {
                    'step': 1,
                    'name': 'Deprotonation (Extreme Conditions)',
                    'equation': 'Ar-X + Strong Base (e.g., NH₂⁻) → Benzyne + HX + X⁻',
                    'details': (
                        "An extremely strong base (usually NaNH₂ in liquid NH₃ at -33°C or at 350°C) "
                        "abstracts a proton ortho to the leaving group. The resulting carbanion expels X⁻ "
                        "to form a highly reactive benzyne intermediate (a benzene ring with a triple bond)."
                    ),
                    'rate_determining': True,
                },
                {
                    'step': 2,
                    'name': 'Nucleophilic Addition to Benzyne',
                    'equation': 'Benzyne + Nu⁻ → Anionic intermediate',
                    'details': (
                        f"The nucleophile adds to one end of the benzyne triple bond. The unsymmetrical "
                        f"distribution of electron density means both ends can be attacked, giving a mixture "
                        f"of regioisomers (unless the benzyne is asymmetrically substituted)."
                    ),
                    'rate_determining': False,
                },
                {
                    'step': 3,
                    'name': 'Protonation',
                    'equation': 'Anionic intermediate + NH₃ → Product',
                    'details': (
                        "The anionic intermediate is protonated by solvent (NH₃) to give the final "
                        "substituted arene product."
                    ),
                    'rate_determining': False,
                },
            ]

        return steps

    def _analyze_intermediate(self, substrate, pathway):
        """Analyze key intermediate."""
        if 'SNAr' in pathway.get('name', ''):

            n_ewg = substrate.get('n_ortho_para_ewg', 0)
            ewgs = [e['type'] for e in substrate.get('ewgs', [])]

            return {
                'intermediate_type': 'Meisenheimer complex (σ-complex)',
                'structure': 'cyclohexadienyl anion — sp³ carbon at site of attack, 5 remaining sp² carbons forming conjugated diene',
                'charge_delocalization': (
                    f"Negative charge delocalized over ortho and para positions relative to attack site. "
                    f"{'Nitro groups provide exceptional resonance stabilization.' if 'NO2' in ewgs else ''}"
                ),
                'stabilizing_factors': [
                    f'{n_ewg} electron-withdrawing group(s) at ortho/para position(s)' if n_ewg > 0 else 'No EWG stabilization',
                    'Resonance delocalization of negative charge',
                    'Inductive withdrawal by electronegative substituents',
                ],
                'key_evidence': 'SNAr reactions show second-order kinetics: rate = k[ArX][Nu⁻]',
            }

        else:
            return {
                'intermediate_type': 'Benzyne',
                'structure': 'Benzene ring with a formal triple bond between two adjacent carbons (distorted geometry)',
                'bond_lengths': 'C≡C ~120 pm (between normal single 154 pm and double 134 pm)',
                'reactivity': 'Extremely reactive — undergoes rapid [2+2] and [4+2] cycloadditions',
                'spectroscopic_signature': 'IR stretch at ~2100 cm⁻¹ (C≡C stretch)',
                'regiochemistry': 'Nu adds to either end of triple bond → mixture unless unsymmetrically substituted',
            }

    def _predict_product(self, substrate, nu, pathway):
        """Predict product."""
        lg = substrate.get('leaving_group', {})
        nu_clean = nu.upper().replace('-', '')

        if 'SNAr' in pathway.get('name', ''):
            prod_map = {
                'OH': 'phenol derivative', 'NH3': 'aniline', 'NH2': 'aniline',
                'CN': 'aryl nitrile', 'CH3O': 'aryl methyl ether', 'N3': 'aryl azide',
                'SH': 'thiophenol', 'SR': 'thioether',
            }
            prod_type = prod_map.get(nu_clean, 'substituted arene')
            return {
                'product_type': prod_type,
                'description': f'LG ({lg.get("name", "X")}) replaced by {nu}. Aromaticity preserved.',
                'regiochemistry': 'Specific — substitution occurs at the carbon bearing the leaving group.',
            }

        else:
            prod_map = {
                'OH': 'phenol', 'NH3': 'aniline', 'NH2': 'aniline',
                'CN': 'benzonitrile', 'CH3O': 'anisole',
            }
            return {
                'product_type': prod_map.get(nu_clean, 'substituted arene'),
                'description': (
                    'Substituted arene via benzyne intermediate. '
                    'Usually a mixture of regioisomers (original site + adjacent site). '
                    'For monosubstituted benzenes: ortho-substituted product predominates.'
                ),
                'regiochemistry': 'Non-specific — mixture of possible products',
            }

    def _evaluate_favorability(self, substrate, pathway, temp_c):
        """Score favorability."""
        score = 0
        activation = substrate.get('activation_level', '')
        lg = substrate.get('leaving_group', {})

        if 'extremely' in activation: score += 5
        elif 'very high' in activation: score += 4
        elif 'high' in activation: score += 3
        elif 'moderate' in activation: score += 2
        elif 'low' in activation: score -= 2

        # Leaving group quality
        q = lg.get('quality', 0)
        if q >= 5: score += 2
        elif q >= 3: score += 1

        # Temperature
        if 'benzyne' in pathway.get('name', ''):
            if temp_c >= 300: score += 2
            elif temp_c >= 100: score += 1
            else: score -= 3

        if score >= 7: return 'excellent — proceeds readily'
        elif score >= 5: return 'good — mild conditions sufficient'
        elif score >= 3: return 'moderate — needs heating'
        elif score >= 1: return 'possible but slow'
        return 'unfavorable — requires activated substrate or extreme conditions'

    def _build_summary(self, substrate, pathway, product, favorability):
        path = pathway.get('name', '?')
        prod = product.get('product_type', '?')
        return f"Nucleophilic aromatic substitution via {path}. Product: {prod}. Favorability: {favorability}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        sub = parts[0] if len(parts) > 0 else ''
        nu = parts[1] if len(parts) > 1 else 'OH-'
        mech = parts[2] if len(parts) > 2 else 'auto'
        temp = float(parts[3]) if len(parts) > 3 else 25.0
        return self._run_base(sub, nu, mech, temp)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
