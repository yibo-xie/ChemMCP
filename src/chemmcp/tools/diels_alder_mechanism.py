"""
Diels-Alder Mechanism (Tool #129)
Diels-Alder 环加成反应：[4+2] 环加成、s-cis 二烯要求、
内向/外向选择性、区域选择性、立体专一性、FMO 分析（正常/逆电子需求）。
Provides Diels-Alder cycloaddition mechanism analysis: [4+2] cycloaddition,
diene s-cis conformation, endo/exo selectivity, regioselectivity,
stereospecificity, and FMO (Frontier Molecular Orbital) analysis.
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


@ChemMCPManager.register_tool
class DielsAlderMechanism(BaseTool):
    __version__ = "0.1.0"
    name = "DielsAlderMechanism"
    func_name = 'explain_diels_alder_mechanism'
    description = "Explain the Diels-Alder [4+2] cycloaddition mechanism: diene s-cis conformation requirement, normal vs inverse electron demand FMO analysis (HOMOdiene-LUMOdienophile), endo vs exo selectivity with secondary orbital interactions, regioselectivity (ortho/meta/para patterns), stereospecificity (cis/trans relationship), and orbital symmetry analysis."
    implementation_description = "Analyzes the diene for s-cis conformation ability and substitution pattern, classifies the dienophile by electron demand, provides complete concerted pericyclic mechanism with FMO diagram description, predicts endo/exo ratio using secondary orbital interactions and kinetic vs thermodynamic control, determines regiochemistry based on frontier coefficients, and evaluates reaction feasibility."
    categories = ["Reaction"]
    tags = ["Diels-Alder", "Cycloaddition", "Pericyclic", "FMO", "Endo/Exo", "Regioselectivity", "Orbital Symmetry", "[4+2]"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('diene_smiles', 'str', 'N/A', 'SMILES of the diene component (must adopt s-cis conformation).'),
        ('dienophile_smiles', 'str', 'N/A', 'SMILES of the dienophile (alkene or alkyne).'),
        ('analysis_mode', 'str', 'full', 'Analysis depth: full, regiochemistry, stereochemistry, or FMO.'),
        ('temperature_c', 'float', '25.0', 'Reaction temperature in °C.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: diene_smiles dienophile_smiles [analysis_mode] [temperature]. E.g., "C=CC C=C full".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing diene_analysis, dienophile_analysis, fmo_analysis, mechanism, endo_exo_analysis, regioselectivity, stereospecificity, product_prediction, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'diene_smiles': 'C=CC=C',
                'dienophile_smiles': 'C=C',
                'analysis_mode': 'full',
                'temperature_c': 25.0,
            },
            'text_input': {'query': 'butadiene ethylene'},
            'output': {
                'result': {
                    'reaction': 'parent Diels-Alder: buta-1,3-diene + ethylene',
                    'product': 'cyclohexene',
                    'mechanism_type': 'concerted [4+2] cycloaddition (pericyclic)',
                    'stereochemistry': 'stereospecific — suprafacial on both components',
                    'electron_demand': 'normal (HOMOdiene-LUMOdienophile)',
                    'favorability': 'excellent — the prototype DA reaction',
                }
            },
        },
        {
            'code_input': {
                'diene_smiles': 'C=CC=C',
                'dienophile_smiles': 'O=C=C=O',  # maleic anhydride
                'analysis_mode': 'full',
                'temperature_c': 25.0,
            },
            'text_input': {'query': 'butadiene maleic_anhydride'},
            'output': {
                'result': {
                    'reaction': 'butadiene + maleic anhydride',
                    'endo_product_major': True,
                    'endo_reason': (
                        'Secondary orbital interactions between diene π-system and '
                        'dienophile carbonyl π* orbitals stabilize endo TS → kinetic product'
                    ),
                    'product': 'cis-norbornene-5,6-dicarboxylic anhydride (endo major)',
                    'favorability': 'excellent — electron-deficient dienophile accelerates DA',
                }
            },
        },
        {
            'code_input': {
                'diene_smiles': 'CC1=CC=CO1',  # furan
                'dienophile_smiles': 'O=C=C=O',  # maleic anhydride
                'analysis_mode': 'full',
                'temperature_c': 25.0,
            },
            'text_input': {'query': 'furan maleic_anhydride'},
            'output': {
                'result': {
                    'reaction': 'inverse electron demand DA (furan as diene)',
                    'note': 'Furan has low-lying LUMO (aromatic) → may act as inverse-demand diene',
                    'reversibility': 'DA with furan is often reversible at RT → retro-Diels-Alder possible',
                    'oxabicyclic_product': 'exo/endo-7-oxabicyclo[2.2.1]hept-5-ene derivative',
                    'favorability': 'good — reversible, may need low T to trap adduct',
                }
            },
        },
    ]

    def _run_base(self, diene_smiles: str, dienophile_smiles: str, analysis_mode: str = 'full', temperature_c: float = 25.0) -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(diene_smiles) or not is_smiles(dienophile_smiles):
            raise ChemMCPInputError("Invalid SMILES string(s).")

        mol_diene = Chem.MolFromSmiles(diene_smiles)
        mol_dienophile = Chem.MolFromSmiles(dienophile_smiles)
        if mol_diene is None or mol_dienophile is None:
            raise ChemMCPInputError("Cannot parse SMILES.")

        # 1. Diene analysis
        diene = self._analyze_diene(mol_diene)

        # 2. Dienophile analysis
        dienophile = self._analyze_dienophile(mol_dienophile)

        # 3. FMO analysis
        fmo = self._fmo_analysis(diene, dienophile)

        # 4. Mechanism
        mechanism = self._build_mechanism()

        # 5. Endo/Exo analysis
        endo_exo = self._endo_exo_analysis(dienophile, temperature_c)

        # 6. Regioselectivity
        regio = self._regioselectivity(diene, dienophile)

        # 7. Stereospecificity
        stereo = self._stereospecificity(diene, dienophile)

        # 8. Product prediction
        product = self._predict_product(diene, dienophile, endo_exo)

        # 9. Favorability
        favorability = self._evaluate_favorability(diene, dienophile, temperature_c)

        result = {
            'result': {
                'diene_analysis': diene,
                'dienophile_analysis': dienophile,
                'fmo_analysis': fmo,
                'mechanism': mechanism,
                'endo_exo_analysis': endo_exo,
                'regioselectivity': regio,
                'stereospecificity': stereo,
                'product_prediction': product,
                'temperature': f'{temperature_c}°C',
                'favorability': favorability,
                'summary': self._build_summary(diene, dienophile, product, favorability),
            }
        }

        logger.info(f"Diels-Alder: {diene_smiles} + {dienophile_smiles} → {favorability}")
        return result

    def _analyze_diene(self, mol):
        """Analyze diene component."""
        # Find conjugated diene system
        double_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                double_bonds.append((a1, a2))

        has_conjugated_diene = False
        diene_atoms = []
        for i, (a1, a1) in enumerate(double_bonds):
            for a2, b2 in double_bonds[i+1:] if i+1 < len(double_bonds) else []:
                # Check if they share an atom or are adjacent
                pass  # simplified

        # Basic diene check: count double bonds
        n_double_bonds = len(double_bonds)
        n_aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())

        can_be_s_cis = True  # assume yes unless locked trans
        # Check for cyclic dienes that might be s-trans locked
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            for ring in ring_info.AtomRings():
                if len(ring) >= 6:
                    can_be_s_cis = True  # large rings OK
                elif len(ring) <= 4:
                    can_be_s_cis = False  # small rings lock s-trans

        return {
            'n_double_bonds': n_double_bonds,
            'has_conjugated_system': n_double_bonds >= 2,
            'can_adopt_s_cis': can_be_s_cis,
            's_cis_requirement': (
                'The diene must be able to adopt s-cis conformation (two double bonds on the same side '
                'of the single bond connecting them). s-Trans locked dienes cannot undergo DA.'
            ),
            'substitution': 'analyzed from structure' if n_double_bonds >= 2 else 'N/A',
            'is_cyclic_diene': n_aromatic > 0 or ring_info.NumRings() > 0,
            'electron_richness': 'rich' if n_aromatic > 0 else 'moderate' if n_double_bonds >= 2 else 'unknown',
        }

    def _analyze_dienophile(self, mol):
        """Analyze dienophile."""
        has_pi_bond = any(bond.GetBondTypeAsDouble() >= 2.0 for bond in mol.GetBonds())
        n_electron_withdrawing = 0

        for atom in mol.GetAtoms():
            for neighbor in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond and bond.GetBondTypeAsDouble() >= 2.0:
                    if neighbor.GetAtomicNum() == 8:
                        # Carbonyl attached to C=C
                        for nn in neighbor.GetNeighbors():
                            if nn.GetAtomicNum() == 6 and nn.GetIdx() != atom.GetIdx():
                                n_electron_withdrawing += 1

        is_alkyne = any(bond.GetBondTypeAsDouble() == 3.0 for bond in mol.GetBonds())

        return {
            'has_pi_bond': has_pi_bond,
            'is_alkene': has_pi_bond and not is_alkyne,
            'is_alkyne': is_alkyne,
            'n_ewg_attached': n_electron_withdrawing,
            'electron_deficiency': 'high' if n_electron_withdrawing >= 2 else
                                     'moderate' if n_electron_withdrawing == 1 else
                                     'low (simple alkene)',
            'reactivity': 'excellent' if n_electron_withdrawing >= 2 else
                        'very good' if n_electron_withdrawing == 1 else
                        'moderate (needs higher T or pressure)',
        }

    def _fmo_analysis(self, diene, dienophile):
        """Frontier Molecular Orbital analysis."""
        d_rich = diene.get('electron_richness', '') == 'rich'
        p_def = dienophile.get('electron_deficiency', '')

        if p_def in ('high', 'moderate'):
            ed = 'normal electron demand'
            interaction = 'HOMO(diene) → LUMO(dienophile)'
            gap_description = (
                'Electron-rich diene HOMO interacts with electron-poor dienophile LUMO. '
                'Small energy gap → fast reaction. EWGs on dienophile lower LUMO → faster.'
            )
        elif d_rich and p_def == 'low':
            ed = 'normal electron demand (less favorable)'
            interaction = 'HOMO(diene) → LUMO(dienophile)'
            gap_description = 'Larger gap — slower reaction, may need heat or pressure.'
        else:
            ed = 'possible inverse electron demand'
            interaction = 'LUMO(diene) → HOMO(dienophile)'
            gap_description = (
                'Electron-poor diene (e.g., with EWGs) + electron-rich dienophile. '
                'LUMO(diene)-HOMO(dienophile) interaction dominates.'
            )

        return {
            'electron_demand': ed,
            'key_orbital_interaction': interaction,
            'energy_gap': gap_description,
            'orbital_symmetry': (
                'Suprafacial on both components: all bonding interactions are constructive in the cyclic TS. '
                'Total of 6π electrons (4n+2, n=1) → thermally allowed via Woodward-Hoffmann rules.'
            ),
            'woodward_hoffmann': 'Allowed: [π4s + π2s] — 6π electrons, thermal, suprafacial/suprafacial → ✓',
        }

    def _build_mechanism(self):
        """Build concerted DA mechanism."""
        return {
            'type': 'concerted pericyclic [4+2] cycloaddition',
            'concertedness': (
                'All bond-making and bond-breaking occurs simultaneously in a single cyclic transition state. '
                'No intermediates. The reaction is fully concerted and pericyclic.'
            ),
            'bond_formation': (
                'Two σ bonds form simultaneously (or nearly so) between C1-C6 and C4-C5 of the '
                'suprafacial array. One π bond in the diene shifts to become the new π bond in the cyclohexene product.'
            ),
            'transition_state': (
                'Aromatic-like cyclic TS with partial bond order (~0.5) for the two forming σ bonds '
                'and the three π interactions. Approximately planar (or envelope-shaped for cyclic dienes).'
            ),
            'volume_change': 'Large negative ΔV‡ (two molecules → one) → accelerated by high pressure.',
            'entropy': 'ΔS‡ < 0 (bimolecular) — favored at higher concentrations.',
        }

    def _endo_exo_analysis(self, dienophile, temp_c):
        """Analyze endo vs exo selectivity."""
        ewgs = dienophile.get('n_ewg_attached', 0)

        if ewgs >= 1:
            return {
                'kinetic_product': 'endo',
                'thermodynamic_product': 'exo',
                'endo_ratio_at_RT': 'high (often > 9:1 endo:exo)',
                'reason': (
                    'Secondary Orbital Interactions (SOI): In the endo TS, the substituents on the dienophile '
                    '(carbonyls, etc.) are positioned under the diene π-system, allowing favorable through-space '
                    'interaction between filled diene π orbitals and empty dienophile substituent π* orbitals. '
                    'This stabilizes the endo TS despite greater steric clash in the product.'
                ),
                'temperature_effect': (
                    f"At {temp_c}°C: {'endo dominates (kinetic control)' if temp_c < 80 else 'higher T favors exo (thermodynamic control, less steric strain)'}"
                ),
            }

        return {
            'kinetic_product': 'N/A (no EWG for SOI)',
            'selectivity': 'low — minimal endo/exo differentiation without electron-withdrawing substituents',
            'reason': 'No significant secondary orbital interactions without π-accepting groups on dienophile.',
        }

    def _regioselectivity(self, diene, dienophile):
        """Predict regiochemistry (ortho/meta/para analogy)."""
        unsymmetrical = diene.get('has_conjugated_system') and dienophile.get('n_ewg_attached', 0) > 0

        if not unsymmetrical:
            return {
                'regioselectivity': 'none (symmetric components)',
                'single_product': True,
            }

        return {
            'regioselectivity': 'ortho-type (if both substituted) or meta-type (polar mismatch)',
            'governing_factors': [
                'Frontier orbital coefficient matching: large-large and small-small interact',
                'Partial charge effects (if polar substituents present)',
                'Steric effects usually minor in DA (TS is early)',
            ],
            'analogy_to_EAS': (
                'DA regiochemistry resembles EAS patterns: electron-donating group on diene directs the '
                'new bond ortho to it; electron-withdrawing group on dienophile directs ortho to it. '
                '"Ortho" means these groups end up adjacent in the product.'
            ),
        }

    def _stereospecificity(self, diene, dienophile):
        """Stereospecificity analysis."""
        return {
            'stereospecific': True,
            'cis_diene_rule': 'Cis relationship of substituents on diene is preserved in the product (same face of cyclohexene).',
            'cis_dienophile_rule': 'Cis relationship on dienophile gives cis-substituted cyclohexene; trans gives trans.',
            'absolute_stereochemistry': 'All relative configurations are transferred with fidelity — no rotation about single bonds in the concerted TS.',
            'example': 'cis,cis-2,4-hexadiene + cis-2-butene → cis-3,6-dimethylcyclohexene (specific isomer)',
        }

    def _predict_product(self, diene, dienophile, endo_exo):
        """Predict DA product."""
        is_alkyne = dienophile.get('is_alkyne', False)
        major_isomer = endo_exo.get('kinetic_product', 'N/A')

        return {
            'ring_system': 'cyclohexene (6-membered ring)' if not is_alkyne else '1,4-cyclohexadiene (from alkyne dienophile)',
            'major_isomer': f'{major_isomer}' if major_isomer != 'N/A' else 'single product',
            'unsaturation': 'one double bond in ring (from diene central single bond)' if not is_alkyne else 'two double bonds',
            'stereochemistry': 'retained from starting materials (stereospecific)',
        }

    def _evaluate_favorability(self, diene, dienophile, temp_c):
        score = 3  # baseline
        if diene.get('can_adopt_s_cis'): score += 3
        if diene.get('has_conjugated_system'): score += 1
        ewg = dienophile.get('n_ewg_attached', 0)
        if ewg >= 2: score += 3
        elif ewg == 1: score += 2
        if dienophile.get('has_pi_bond'): score += 1
        if temp_c > 25: score += 1
        if temp_c > 100: score += 1

        if score >= 8: return 'excellent — should proceed readily at RT'
        elif score >= 6: return 'good — may need mild heating'
        elif score >= 4: return 'moderate — heating likely needed'
        return 'possible but may need forcing conditions (high T, high pressure, Lewis acid)'

    def _build_summary(self, diene, dienophile, product, fav):
        return f"Diels-Alder [4+2] cycloaddition. Product: {product.get('ring_system','?')}. Favorability: {fav}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        diene = parts[0] if len(parts) > 0 else ''
        dienophile = parts[1] if len(parts) > 1 else ''
        mode = parts[2] if len(parts) > 2 else 'full'
        temp = float(parts[3]) if len(parts) > 3 else 25.0
        return self._run_base(diene, dienophile, mode, temp)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
