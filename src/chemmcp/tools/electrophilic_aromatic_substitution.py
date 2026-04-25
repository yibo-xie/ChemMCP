"""
Electrophilic Aromatic Substitution (Tool #122)
芳香亲电取代通用机理（σ络合物）：硝化、磺化、卤代、傅-克烷基化/酰基化。
Provides EAS mechanism analysis with σ-complex (arenium ion) intermediate,
orientation effects, and activating/deactivating group analysis.
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


# Electrophile database for EAS
ELECTROPHILES = {
    'NO2+': {
        'name': 'Nitration (NO₂⁺)', 'reagent': 'HNO₃/H₂SO₄',
        'product': 'nitroarene', 'conditions': 'mixed acid, 50-60°C',
        'electrophile_source': 'HNO₃ + 2H₂SO₄ → NO₂⁺ + H₃O⁺ + 2HSO₄⁻',
    },
    'SO3': {
        'name': 'Sulfonation (SO₃)', 'reagent': f'H₂SO₄ (conc.) / SO₃ / fuming H₂SO₄',
        'product': 'sulfonic acid', 'conditions': 'conc. H₂SO₄, heat',
        'electrophile_source': '2H₂SO₄ ⇌ SO₃ + H₃O⁺ + HSO₄⁻',
    },
    'Br2': {
        'name': 'Bromination (Br₂)', 'reagent': 'Br₂/FeBr₃ or Br₂/Fe',
        'product': 'bromoarene', 'conditions': 'FeBr₃ or AlBr₃ catalyst, RT',
        'electrophile_source': 'Br₂ + FeBr₃ → Br⁺FeBr₄⁻ (polarized complex)',
    },
    'Cl2': {
        'name': 'Chlorination (Cl₂)', 'reagent': 'Cl₂/FeCl₃ or Cl₂/Fe',
        'product': 'chloroarene', 'conditions': 'FeCl₃ or AlCl₃ catalyst, RT',
        'electrophile_source': 'Cl₂ + FeCl₃ → Cl⁺FeCl₄⁻',
    },
    'R+': {
        'name': 'Friedel-Crafts Alkylation (R⁺)', 'reagent': 'R-X/AlCl₃',
        'product': 'alkylbenzene', 'conditions': 'AlCl₃ catalyst, < RT to prevent rearrangement',
        'electrophile_source': 'R-X + AlCl₃ → R⁺(AlCl₃X)⁻ (carbocation)',
    },
    'RCO+': {
        'name': 'Friedel-Crafts Acylation (RCO⁺)', 'reagent': 'RCOCl/AlCl₃',
        'product': 'aryl ketone', 'conditions': 'AlCl₃ catalyst, RT-reflux',
        'electrophile_source': 'RCOCl + AlCl₃ → R-C⁺=O (AlCl₃Cl)⁻ (acylium ion)',
    },
}

# Substituent effects on EAS
SUBSTITUENT_EFFECTS = {
    # (effect, orientation, strength) — effect: activating/deactivating
    '-NH2': ('strongly activating', 'ortho/para', '+M dominates'),
    '-NHR': ('strongly activating', 'ortho/para', '+M dominates'),
    '-NR2': ('strongly activating', 'ortho/para', '+M dominates'),
    '-OH': ('strongly activating', 'ortho/para', '+M > -I'),
    '-OR': ('strongly activating', 'ortho/para', '+M > -I'),
    '-NHCOR': ('moderately activating', 'ortho/para', '+M (resonance with N)'),
    '-OCOR': ('weakly activating', 'ortho/para', '+M weaker'),
    '-R':  ('weakly activating', 'ortho/para', '+I (hyperconjugation)'),
    '-Ph': ('weakly activating', 'ortho/para', '+I (hyperconjugation)'),
    '-X(F,Cl,Br,I)': ('weakly deactivating', 'ortho/para', '+M < -I (halogen exception)'),
    '-CH=CHR': ('weakly activating', 'ortho/para', '+M (conjugation)'),
    '-CHO': ('moderately deactivating', 'meta', '-M and -I (electron withdrawal)'),
    '-COR': ('moderately deactivating', 'meta', '-M and -I'),
    '-COOH': ('moderately deactivating', 'meta', '-M and -I'),
    '-COOR': ('moderately deactivating', 'meta', '-M and -I'),
    '-CONH2': ('moderately deactivating', 'meta', '-M'),
    '-CF3': ('strongly deactivating', 'meta', '-I only'),
    '-CN': ('strongly deactivating', 'meta', '-M and -I'),
    '-NO2': ('strongly deactivating', 'meta', '-M and -I'),
    '-SO3H': ('strongly deactivating', 'meta', '-M and -I'),
    '-NR3+': ('strongly deactivating', 'meta', '+I (positive charge)'),
    '-CCl3': ('strongly deactivating', 'meta', '-I only'),
}


@ChemMCPManager.register_tool
class ElectrophilicAromaticSubstitution(BaseTool):
    __version__ = "0.1.0"
    name = "ElectrophilicAromaticSubstitution"
    func_name = 'explain_eas_mechanism'
    description = "Explain electrophilic aromatic substitution (EAS) mechanism: σ-complex (arenium ion) formation, resonance stabilization, orientation (ortho/meta/para), and rate effects of substituents. Covers nitration, sulfonation, halogenation, Friedel-Crafts alkylation/acylation."
    implementation_description = "Analyzes the aromatic substrate for existing substituents and their directing effects, matches with specified electrophile type, provides complete three-step EAS mechanism (electrophilic attack → arenium ion → deprotonation), predicts major/minor products based on orientation rules, and evaluates reaction feasibility."
    categories = ["Reaction"]
    tags = ["EAS", "Electrophilic Aromatic Substitution", "Arenium Ion", "Sigma Complex", "Orientation", "Friedel-Crafts"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('aromatic_smiles', 'str', 'N/A', 'SMILES string of the aromatic substrate.'),
        ('electrophile_type', 'str', 'NO2+', 'Type of EAS electrophile. Options: NO2+, SO3, Br2, Cl2, R+, RCO+.'),
        ('temperature_c', 'float', '25.0', 'Reaction temperature in °C.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: aromatic_smiles electrophile_type [temperature_c]. E.g., "c1ccccc1 NO2+ 55".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing substituent_analysis, orientation_prediction, mechanism_steps, arenium_ion_resonance, product_prediction, rate_effect, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'aromatic_smiles': 'c1ccccc1',
                'electrophile_type': 'NO2+',
                'temperature_c': 55.0,
            },
            'text_input': {'query': 'c1ccccc1 NO2+ 55'},
            'output': {
                'result': {
                    'substrate': 'benzene',
                    'reaction': 'Nitration',
                    'mechanism_type': 'EAS via σ-complex (arenium ion)',
                    'steps': [
                        {'step': 1, 'description': 'Electrophilic attack: π electrons attack NO₂⁺ → arenium ion (σ-complex)'},
                        {'step': 2, 'description': 'Deprotonation: base removes H⁺ → restores aromaticity → nitrobenzene'},
                    ],
                    'orientation': 'all positions equivalent (monosubstituted benzene)',
                    'product': 'nitrobenzene (single product)',
                    'favorability': 'excellent — benzene is highly reactive toward EAS',
                }
            },
        },
        {
            'code_input': {
                'aromatic_smiles': 'Cc1ccccc1',
                'electrophile_type': 'NO2+',
                'temperature_c': 55.0,
            },
            'text_input': {'query': 'Cc1ccccc1 NO2+'},
            'output': {
                'result': {
                    'substrate': 'toluene (methylbenzene)',
                    'existing_substituent': '-CH₃ (alkyl)',
                    'substituent_effect': 'weakly activating, ortho/para-directing',
                    'major_products': ['o-nitrotoluene (~58%)', 'p-nitrotoluene (~38%)', 'm-nitrotoluene (~4%)'],
                    'orientation_ratio': 'ortho:para ≈ 3:2 (with some steric hindrance at ortho)',
                    'rate_vs_benzene': '~25× faster than benzene',
                    'favorability': 'excellent — activated ring',
                }
            },
        },
        {
            'code_input': {
                'aromatic_smiles': 'O=c1ccccc1',
                'electrophile_type': 'NO2+',
                'temperature_c': 55.0,
            },
            'text_input': {'query': 'O=c1ccccc1 NO2+'},
            'output': {
                'result': {
                    'substrate': 'benzaldehyde',
                    'existing_substituent': '-CHO (formyl)',
                    'substituent_effect': 'moderately deactivating, meta-directing',
                    'major_product': 'm-nitrobenzaldehyde (>80%)',
                    'rate_vs_benzene': '~100× slower than benzene',
                    'favorability': 'possible — deactivated ring requires harsher conditions',
                }
            },
        },
    ]

    def _run_base(self, aromatic_smiles: str, electrophile_type: str = 'NO2+', temperature_c: float = 25.0) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(aromatic_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol = Chem.MolFromSmiles(aromatic_smiles)
        if mol is None:
            raise ChemMCPInputError("Cannot parse SMILES string.")

        # 1. Analyze aromatic system
        aromatic_analysis = self._analyze_aromatic(mol)

        # 2. Get electrophile info
        e_info = self._get_electrophile(electrophile_type)

        # 3. Identify substituents and their effects
        sub_effects = self._identify_substituents(mol)

        # 4. Predict orientation
        orientation = self._predict_orientation(sub_effects)

        # 5. Build mechanism steps
        steps = self._build_mechanism_steps(aromatic_analysis, e_info, sub_effects)

        # 6. Arenium ion resonance
        resonance = self._describe_arenium_resonance(sub_effects, e_info)

        # 7. Product prediction
        products = self._predict_products(aromatic_analysis, sub_effects, e_info, orientation)

        # 8. Rate effect
        rate = self._evaluate_rate(sub_effects)

        # 9. Favorability
        favorability = self._evaluate_favorability(aromatic_analysis, sub_effects, e_info, temperature_c)

        result = {
            'result': {
                'substrate_smiles': aromatic_smiles,
                'aromatic_analysis': aromatic_analysis,
                'electrophile_info': e_info,
                'substituent_effects': sub_effects,
                'orientation_prediction': orientation,
                'mechanism_steps': steps,
                'arenium_ion_resonance': resonance,
                'product_prediction': products,
                'rate_effect': rate,
                'temperature': f'{temperature_c}°C',
                'favorability': favorability,
                'summary': self._build_summary(aromatic_analysis, e_info, products, favorability),
            }
        }

        logger.info(f"EAS: {aromatic_smiles} + {electrophile_type} → {favorability}")
        return result

    def _analyze_aromatic(self, mol):
        """Analyze aromatic system."""
        has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
        n_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        n_rings = mol.GetRingInfo().NumRings()

        return {
            'has_aromatic_ring': has_aromatic,
            'n_aromatic_atoms': n_aromatic_atoms,
            'n_rings': n_rings,
            'is_benzene_core': n_aromatic_atoms == 6,
        }

    def _get_electrophile(self, e_type):
        """Get electrophile information."""
        e_type_clean = e_type.upper().replace(' ', '')

        # Direct match
        if e_type in ELECTROPHILES:
            return {**ELECTROPHILES[e_type], 'key': e_type}

        # Fuzzy match
        fuzzy = {
            'NITRATION': 'NO2+', 'NO2': 'NO2+', 'NITRO': 'NO2+',
            'SULFONATION': 'SO3', 'SO3H': 'SO3', 'SULFONIC': 'SO3',
            'BROMINATION': 'BR2', 'BR': 'BR2',
            'CHLORINATION': 'CL2', 'CL': 'CL2',
            'ALKYLATION': 'R+', 'FRIEDEL': 'R+', 'FCALKYL': 'R+',
            'ACYLATION': 'RCO+', 'FCACYL': 'RCO+', 'ACYL': 'RCO+',
        }
        for key, real in fuzzy.items():
            if key in e_type_clean:
                return {**ELECTROPHILES[real], 'key': real, 'matched_from': e_type}

        available = ', '.join(ELECTROPHILES.keys())
        return {'name': e_type, 'error': f'Unknown electrophile. Available: {available}'}

    def _identify_substituents(self, mol):
        """Identify substituents on aromatic ring(s)."""
        substituents = []

        for atom in mol.GetAtoms():
            if not atom.GetIsAromatic():
                continue

            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    continue

                sym = neighbor.GetSymbol()
                # Build a simple representation of the substituent
                sub_desc = self._describe_substituent(neighbor, mol)
                effect_key = self._match_substituent_to_effect(sub_desc.get('label', ''))

                substituents.append({
                    'atom_idx': neighbor.GetIdx(),
                    'symbol': sym,
                    'attachment_position': atom.GetIdx(),
                    'description': sub_desc,
                    'effect': SUBSTITUENT_EFFECTS.get(effect_key, ('unknown', 'unknown', 'unknown')),
                    'effect_key': effect_key,
                })

        return substituents

    def _describe_substituent(self, atom, mol):
        """Describe a non-aromatic substituent attached to ring."""
        sym = atom.GetSymbol()

        if sym == 'C':
            # Count neighbors to determine if it's CH3, CHO, COOH, etc.
            neighbors = atom.GetNeighbors()
            n_neighbors = len(neighbors)

            has_double_bond_o = any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()) and
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
                for n in neighbors if n.GetAtomicNum() == 8
            )

            if has_double_bond_o:
                return {'type': 'carbonyl-containing', 'label': '-COR/-CHO'}
            if n_neighbors <= 1:
                return {'type': 'alkyl (likely methyl)', 'label': '-R'}
            return {'type': 'carbon-based', 'label': '-R'}

        elif sym == 'O':
            neighbors = atom.GetNeighbors()
            if any(n.GetAtomicNum() == 1 for n in neighbors):  # O-H
                return {'type': 'hydroxyl', 'label': '-OH'}
            return {'type': 'ether/ester oxygen', 'label': '-OR'}

        elif sym == 'N':
            return {'type': 'nitrogen-containing', 'label': '-NR2'}

        elif sym in ('F', 'Cl', 'Br', 'I'):
            return {'type': 'halogen', 'label': f'-{sym}'}

        elif sym == 'S':
            return {'type': 'sulfur-containing', 'label': '-SR'}

        return {'type': 'other', 'label': f'-{sym}X'}

    def _match_substituent_to_effect(self, label):
        """Match substituent label to known effects."""
        label_upper = label.upper().replace(' ', '')
        for key in SUBSTITUENT_EFFECTS:
            key_clean = key.upper().replace('-', '').replace('(', '').replace(')', '')
            label_clean = label_upper.replace('-', '').replace('(', '').replace(')', '')
            if key_clean in label_clean or label_clean in key_clean:
                return key
        return label

    def _predict_orientation(self, sub_effects):
        """Predict orientation based on substituents."""
        if not sub_effects:
            return {
                'direction': 'all positions equivalent',
                'reason': 'No substituents (or benzene) — all 6 positions equivalent.',
                'ortho_ratio': None,
                'meta_ratio': None,
                'para_ratio': None,
            }

        # Find strongest activator/deactivator
        strength_order = [
            'strongly activating', 'moderately activating', 'weakly activating',
            'weakly deactivating', 'moderately deactivating', 'strongly deactivating'
        ]

        best = min(sub_effects, key=lambda s: strength_order.index(s['effect'][0]) if s['effect'][0] in strength_order else 99)
        effect = best['effect']

        direction = effect[1]  # ortho/para or meta
        strength = effect[0]

        # Estimate ratios
        if 'ortho/para' in direction:
            if 'strongly' in strength:
                return {'direction': 'ortho/para', 'strength': strength, 'ortho_ratio': '~50%', 'meta_ratio': '<5%', 'para_ratio': '~45%'}
            elif 'weakly' in strength:
                if 'deactivating' in strength:
                    return {'direction': 'ortho/para (halogen exception)', 'strength': strength, 'ortho_ratio': '~12%', 'meta_ratio': '<1%', 'para_ratio': '~87%'}
                return {'direction': 'ortho/para', 'strength': strength, 'ortho_ratio': '~45%', 'meta_ratio': '~5%', 'para_ratio': '~50%'}

        elif 'meta' in direction:
            if 'strongly' in strength:
                return {'direction': 'meta', 'strength': strength, 'ortho_ratio': '~trace', 'meta_ratio': '>95%', 'para_ratio': '~trace'}
            return {'direction': 'meta', 'strength': strength, 'ortho_ratio': '~7%', 'meta_ratio': '~85%', 'para_ratio': '~8%'}

        return {'direction': 'unknown', 'reason': 'Cannot determine orientation.'}

    def _build_mechanism_steps(self, aromatic, e_info, sub_effects):
        """Build EAS mechanism steps."""
        e_name = e_info.get('name', 'E⁺')
        steps = [
            {
                'step': 1,
                'name': 'Electrophilic Attack (Rate-Determining Step)',
                'equation': f'Ar-H + {e_name.split("(")[0].strip()} → Arenium ion (σ-complex)',
                'details': (
                    f"The electron-rich aromatic π-system attacks the electrophile ({e_name}). "
                    f"A C-E bond forms as one π bond breaks, producing a delocalized "
                    f"arenium ion (σ-complex) with a positive charge distributed over the ortho and para positions "
                    f"(3 resonance structures). The aromaticity is temporarily lost — this is the high-energy, "
                    f"rate-determining step."
                ),
                'rate_determining': True,
            },
            {
                'step': 2,
                'name': 'Deprotonation (Fast)',
                'equation': 'Arenium ion + :B → Substituted arene + BH⁺',
                'details': (
                    "A weak base (e.g., HSO₄⁻, FeBr₄⁻, or H₂O) removes the proton from the sp³ carbon. "
                    "The electron pair from the C-H bond reforms the π bond, restoring aromatic stability. "
                    "This step is fast and exergonic because it regains aromatic stabilization energy (~152 kJ/mol)."
                ),
                'rate_determining': False,
            },
        ]

        # Add special notes for specific reactions
        key = e_info.get('key', '')
        if key == 'R+':
            steps.append({
                'step': 3,
                'note': 'Friedel-Crafts Alkylation Caveat',
                'details': (
                    "⚠ Carbocation intermediates may undergo rearrangement (hydride/alkyl shift) "
                    "to form more stable carbocations, leading to unexpected products. "
                    "Polyalkylation is common due to alkyl group activation. "
                    "Consider FC acylation + reduction instead for unambiguous products."
                ),
            })
        elif key == 'SO3':
            steps.append({
                'step': 3,
                'note': 'Reversibility of Sulfonation',
                'details': (
                    "Sulfonation is reversible at elevated temperature. "
                    "This can be used strategically to block a position (sulfonation → substitution → desulfonation)."
                ),
            })

        return steps

    def _describe_arenium_resonance(self, sub_effects, e_info):
        """Describe arenium ion resonance structures."""
        if not sub_effects:
            return {
                'n_resonance_structures': 3,
                'charge_distribution': 'positive charge on ortho (2 positions) and para (1 position)',
                'resonance_description': (
                    "For monosubstituted benzene: the positive charge of the arenium ion "
                    "is delocalized over the two ortho positions and the para position "
                    "(3 resonance contributors). The meta position never bears positive charge."
                ),
            }

        best = sub_effects[0]
        direction = best['effect'][1]

        if 'ortho/para' in direction:
            return {
                'n_resonance_structures': 3,
                'charge_distribution': 'positive charge on ortho and para positions',
                'stabilization': (
                    f"The {best.get('description', {}).get('label', '?')} substituent can donate electron density "
                    f"by resonance into the ring, stabilizing the arenium ion when the electrophile attacks "
                    f"ortho or para (where the positive charge lands on the carbon bearing the substituent in one "
                    f"resonance structure). This lowers the activation energy for ortho/para attack."
                ),
            }
        else:
            return {
                'n_resonance_structures': 3,
                'charge_distribution': 'positive charge on ortho and para positions',
                'destabilization': (
                    f"The {best.get('description', {}).get('label', '?')} substituent withdraws electron density. "
                    f"Meta attack avoids placing positive charge directly on the substituted carbon in any "
                    f"resonance structure, making meta substitution energetically favored despite overall deactivation."
                ),
            }

    def _predict_products(self, aromatic, sub_effects, e_info, orientation):
        """Predict products."""
        product_name = e_info.get('product', 'substituted arene')
        direction = orientation.get('direction', 'unknown')

        if not sub_effects:
            return {
                'major_product': f'{product_name} (single product, all positions equivalent)',
                'isomer_distribution': 'N/A (symmetric)',
            }

        o_r = orientation.get('ortho_ratio', '?')
        m_r = orientation.get('meta_ratio', '?')
        p_r = orientation.get('para_ratio', '?')

        return {
            'major_product_direction': direction,
            'isomer_distribution': f'ortho: {o_r}, meta: {m_r}, para: {p_r}',
            'product_type': product_name,
        }

    def _evaluate_rate(self, sub_effects):
        """Evaluate relative rate vs benzene."""
        if not sub_effects:
            return {'relative_rate': '1× (reference = benzene)', 'classification': 'unsubstituted'}

        strength = sub_effects[0]['effect'][0]
        rate_map = {
            'strongly activating': '1000–10000× faster than benzene',
            'moderately activating': '10–100× faster than benzene',
            'weakly activating': '2–5× faster than benzene',
            'weakly deactivating': '0.03–0.1× (slower than benzene; halogen exception)',
            'moderately deactivating': '10⁻³–10⁻²× (much slower)',
            'strongly deactivating': '<10⁻⁶× (very slow; may require extreme conditions)',
        }
        return {
            'relative_rate': rate_map.get(strength, 'unknown'),
            'classification': strength,
        }

    def _evaluate_favorability(self, aromatic, sub_effects, e_info, temp_c):
        """Score favorability."""
        score = 5  # baseline: EAS works on most aromatics

        if not sub_effects:
            score += 3  # benzene itself is fine
        else:
            strength = sub_effects[0]['effect'][0]
            if 'activating' in strength: score += 2
            elif 'deactivating' in strength:
                score -= 2
                if 'strongly' in strength: score -= 2

        # Temperature helps for deactivated rings
        if temp_c > 50: score += 1
        if temp_c > 100: score += 1

        # Special case: strongly deactivated rings need more forcing conditions
        if sub_effects and 'strongly deactivating' in sub_effects[0]['effect'][0]:
            if temp_c < 50: score -= 2

        # FC alkylation limitations
        if e_info.get('key') == 'R+':
            score -= 1  # potential rearrangement issues

        if score >= 8: return 'excellent'
        elif score >= 6: return 'good'
        elif score >= 4: return 'moderate'
        elif score >= 2: return 'possible but may need forcing conditions'
        return 'difficult — consider alternative routes'

    def _build_summary(self, aromatic, e_info, products, favorability):
        rxn_name = e_info.get('name', 'EAS')
        prod = products.get('major_product', '?')
        return f"EAS ({rxn_name}). Product: {prod}. Favorability: {favorability}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        arom = parts[0] if len(parts) > 0 else ''
        e_type = parts[1] if len(parts) > 1 else 'NO2+'
        temp = float(parts[2]) if len(parts) > 2 else 25.0
        return self._run_base(arom, e_type, temp)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
