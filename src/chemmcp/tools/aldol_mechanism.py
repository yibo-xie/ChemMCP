"""
Aldol Mechanism (Tool #126)
羟醛缩合反应机理：碱催化/酸催化烯醇(盐)形成、亲核加成、
α,β-不饱和羰基化合物脱水。
Provides aldol condensation mechanism analysis: enolate formation (base-catalyzed)
or enol formation (acid-catalyzed), nucleophilic addition to carbonyl, and dehydration.
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
class AldolMechanism(BaseTool):
    __version__ = "0.1.0"
    name = "AldolMechanism"
    func_name = 'explain_aldol_mechanism'
    description = "Explain the aldol addition and aldol condensation mechanisms: enolate/enol formation (acid or base catalysis), nucleophilic attack on a carbonyl, β-hydroxy carbonyl product (aldol), and dehydration to α,β-unsaturated carbonyl compound. Covers crossed aldol selectivity and stereochemistry."
    implementation_description = "Analyzes both carbonyl components for α-hydrogen availability, determines enolizable vs non-enolizable partners, provides complete stepwise mechanism for base-catalyzed (enolate) or acid-catalyzed (enol) pathways, predicts crossed vs self-aldol outcomes, analyzes dehydration thermodynamics, and evaluates syn/anti diastereoselectivity."
    categories = ["Reaction"]
    tags = ["Aldol", "Condensation", "Enolate", "Carbonyl", "Dehydration", "Stereochemistry"]
    required_envs = []
    oss_dependencies = [
        ("RDKit", "https://github.com/rdkit/rdkit", "BSD"),
    ]
    services_and_software = []

    code_input_sig = [
        ('component1_smiles', 'str', 'N/A', 'SMILES of component 1 (enolizable aldehyde/ketone).'),
        ('component2_smiles', 'str', '', 'SMILES of component 2 (electrophilic aldehyde/ketone). Leave empty for self-aldol.'),
        ('catalyst', 'str', 'base', 'Catalyst type: base (e.g., OH-, LDA, NaOEt) or acid (e.g., H3O+, H+).'),
        ('conditions', 'str', 'RT then heat', 'Reaction conditions: RT, heat, etc.'),
    ]
    text_input_sig = [
        ('query', 'str', 'N/A', 'Space-separated: component1_smiles [component2_smiles] [catalyst] [conditions]. E.g., "CC=O CC=O base".'),
    ]
    output_sig = [
        ('result', 'dict', 'Dictionary containing enolate_analysis, mechanism_steps, aldol_product, dehydration_analysis, stereochemistry, crossed_aldol_analysis, and favorability.'),
    ]
    examples = [
        {
            'code_input': {
                'component1_smiles': 'CC=O',
                'component2_smiles': '',
                'catalyst': 'base',
                'conditions': 'RT then heat',
            },
            'text_input': {'query': 'CC=O'},
            'output': {
                'result': {
                    'reaction_type': 'self-aldol of acetaldehyde',
                    'catalysis': 'base-catalyzed (enolate pathway)',
                    'enolate_source': 'acetaldehyde → enolate (nucleophile)',
                    'electrophile': 'acetaldehyde (another molecule)',
                    'aldol_product': '3-hydroxybutanal (aldol)',
                    'dehydration_product': 'crotonaldehyde (but-2-enal)',
                    'mechanism_steps': [
                        {'step': 1, 'description': 'Base abstracts α-H → enolate'},
                        {'step': 2, 'description': 'Enolate attacks carbonyl C of another acetaldehyde'},
                        {'step': 3, 'description': 'Protonation → β-hydroxy aldehyde (aldol)'},
                        {'step': 4, 'description': 'Heat-induced dehydration → α,β-unsaturated aldehyde'},
                    ],
                    'favorability': 'excellent — classic aldol reaction',
                }
            },
        },
        {
            'code_input': {
                'component1_smiles': 'CC(=O)C',
                'component2_smiles': 'CC=O',
                'catalyst': 'base',
                'conditions': 'OH-/heat',
            },
            'text_input': {'query': 'CC(=O)C CC=O base'},
            'output': {
                'result': {
                    'reaction_type': 'crossed aldol (acetone + acetaldehyde)',
                    'enolate_source': 'acetone (more acidic α-H, pKa ~20) → enolate',
                    'electrophile': 'acetaldehyde (more reactive electrophile than ketone)',
                    'product': '4-hydroxy-4-methylpentan-2-one (aldol)',
                    'dehydration_product': 'mesityl oxide derivative',
                    'selectivity': 'Acetone enolate prefers acetaldehyde (ketone is less reactive electrophile)',
                    'favorability': 'good — controlled crossed aldol',
                }
            },
        },
    ]

    def _run_base(self, component1_smiles: str, component2_smiles: str = '', catalyst: str = 'base', conditions: str = 'RT then heat') -> dict:
        if not _RDKIT_AVAILABLE:
            raise ChemMCPInputError("RDKit is not available.")
        if not is_smiles(component1_smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")

        mol1 = Chem.MolFromSmiles(component1_smiles)
        if mol1 is None:
            raise ChemMCPInputError("Cannot parse SMILES string for component 1.")

        mol2 = None
        if component2_smiles:
            if not is_smiles(component2_smiles):
                raise ChemMCPInputError("Component 2 SMILES is invalid.")
            mol2 = Chem.MolFromSmiles(component2_smiles)

        # 1. Analyze both components
        comp1_analysis = self._analyze_carbonyl_component(mol1, component1_smiles)
        comp2_analysis = self._analyze_carbonyl_component(mol2, component2_smiles) if mol2 else None

        # 2. Determine reaction type (self vs crossed)
        is_crossed = mol2 is not None

        # 3. Enolate analysis
        enolate_info = self._analyze_enolate(comp1_analysis, catalyst)

        # 4. Build mechanism steps
        steps = self._build_mechanism(comp1_analysis, comp2_analysis, catalyst, is_crossed)

        # 5. Product prediction
        products = self._predict_products(comp1_analysis, comp2_analysis, is_crossed)

        # 6. Dehydration analysis
        dehydration = self._analyze_dehydration(products, conditions)

        # 7. Stereochemistry
        stereo = self._analyze_stereochemistry(comp1_analysis, comp2_analysis, catalyst)

        # 8. Crossed aldol considerations
        crossed_analysis = self._analyze_crossed(comp1_analysis, comp2_analysis, is_crossed)

        # 9. Favorability
        favorability = self._evaluate_favorability(comp1_analysis, comp2_analysis, catalyst, is_crossed)

        result = {
            'result': {
                'component1': comp1_analysis,
                'component2': comp2_analysis,
                'reaction_type': f'crossed aldol' if is_crossed else 'self-aldol',
                'catalyst': catalyst,
                'enolate_analysis': enolate_info,
                'mechanism_steps': steps,
                'product_prediction': products,
                'dehydration_analysis': dehydration,
                'stereochemistry': stereo,
                'crossed_aldol_analysis': crossed_analysis,
                'favorability': favorability,
                'summary': self._build_summary(comp1_analysis, comp2_analysis, products, favorability, is_crossed),
            }
        }

        logger.info(f"Aldol: {component1_smiles} + {component2_smiles or 'self'} ({catalyst}) → {favorability}")
        return result

    def _analyze_carbonyl_component(self, mol, smiles_str):
        """Analyze a carbonyl component."""
        if mol is None:
            return None

        has_carbonyl = False
        has_alpha_h = False
        alpha_carbons = []
        c_type = 'unknown'

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() == 2.0 and neighbor.GetAtomicNum() == 8:
                        has_carbonyl = True
                        n_non_h = sum(1 for n in atom.GetNeighbors()
                                      if n.GetAtomicNum() not in (1,) and n.GetIdx() != neighbor.GetIdx())
                        c_type = 'aldehyde' if n_non_h <= 1 else 'ketone'

                        # Check alpha positions for H
                        for n in atom.GetNeighbors():
                            if n.GetAtomicNum() == 6 and n.GetTotalNumHs() > 0:
                                has_alpha_h = True
                                n_alpha_h = n.GetTotalNumHs()
                                alpha_carbons.append({
                                    'idx': n.GetIdx(),
                                    'n_alpha_h': n_alpha_h,
                                    'is_primary': sum(1 for nn in n.GetNeighbors() if nn.GetAtomicNum() == 6) <= 1,
                                })

        return {
            'smiles': smiles_str,
            'has_carbonyl': has_carbonyl,
            'carbonyl_type': c_type,
            'has_alpha_hydrogen': has_alpha_h,
            'alpha_carbons': alpha_carbons,
            'can_form_enolate': has_alpha_h,
            'reactivity_as_electrophile': 'high' if c_type == 'aldehyde' else 'moderate' if c_type == 'ketone' else 'low',
            'alpha_h_acidity': 'pKa ~17 (aldehyde)' if c_type == 'aldehyde' else 'pKa ~20 (ketone)' if c_type == 'ketone' else 'N/A',
        }

    def _analyze_enolate(self, comp_analysis, catalyst):
        """Analyze enolate formation."""
        if not comp_analysis or not comp_analysis.get('can_form_enolate'):
            return {'can_form_enolate': False, 'note': 'No α-H available — cannot form enolate.'}

        cat_lower = catalyst.lower()
        if cat_lower.startswith('bas'):
            return {
                'pathway': 'enolate (base-catalyzed)',
                'mechanism': (
                    'Strong base (OH⁻, RO⁻, LDA, etc.) deprotonates α-carbon → '
                    'resonance-stabilized enolate anion (nucleophilic at α-carbon AND oxygen)'
                ),
                'kinetic_vs_thermodynamic': (
                    'Low T, sterically hindered base (LDA) → kinetic enolate (less substituted, O-bound). '
                    'Higher T, smaller base → thermodynamic enolate (more substituted).'
                ) if 'LDA' in catalyst or 'lda' in catalyst.lower() else 'Standard base-catalyzed enolate.',
            }
        else:
            return {
                'pathway': 'enol (acid-catalyzed)',
                'mechanism': (
                    'Acid protonates carbonyl oxygen → tautomerization via enol '
                    '(neutral nucleophile; less reactive than enolate but avoids polycondensation)'
                ),
            }

    def _build_mechanism(self, comp1, comp2, catalyst, is_crossed):
        """Build stepwise mechanism."""
        cat_lower = catalyst.lower()
        is_base = cat_lower.startswith('bas')
        donor = comp1
        acceptor = comp2 if is_crossed else comp1

        if is_base:
            steps = [
                {
                    'step': 1, 'name': 'Enolate Formation',
                    'equation': f'{donor.get("carbonyl_type","?")} + B⁻ → Enolate + BH',
                    'details': f"Base abstracts α-proton from {donor.get('carbonyl_type','?')} to form resonance-stabilized enolate.",
                },
                {
                    'step': 2, 'name': 'Nucleophilic Addition (Aldol Addition)',
                    'equation': 'Enolate + Carbonyl → Alkoxide intermediate',
                    'details': (
                        f"Enolate α-carbon attacks the carbonyl carbon of "
                        f"{acceptor.get('carbonyl_type','electrophile') if acceptor else 'another molecule of donor'} "
                        f"→ tetrahedral alkoxide intermediate."
                    ),
                },
                {
                    'step': 3, 'name': 'Protonation',
                    'equation': 'Alkoxide + H⁺ → β-hydroxy carbonyl (aldol product)',
                    'details': 'The alkoxide is protonated to give the neutral β-hydroxy aldehyde/ketone (aldol).',
                },
            ]
        else:
            steps = [
                {
                    'step': 1, 'name': 'Carbonyl Protonation',
                    'equation': 'C=O + H⁺ → C⁺-OH (activated carbonyl)',
                    'details': 'Acid activates carbonyl toward nucleophilic attack.',
                },
                {
                    'step': 2, 'name': 'Enol Formation (Tautomerization)',
                    'equation': f'{donor.get("carbonyl_type","?")} ⇌ Enol',
                    'details': 'Acid-catalyzed keto-enol tautomerization generates enol nucleophile.',
                },
                {
                    'step': 3, 'name': 'Nucleophilic Attack of Enol',
                    'equation': 'Enol + activated C=O → Oxocarbenium intermediate',
                    'details': 'Enol attacks protonated carbonyl.',
                },
                {
                    'step': 4, 'name': 'Deprotonation',
                    'equation': 'Intermediate → β-hydroxy carbonyl',
                    'details': 'Deprotonation yields the aldol product.',
                },
            ]

        # Dehydration step (separate, usually requires heat)
        steps.append({
            'step': len(steps) + 1, 'name': 'Dehydration (Aldol Condensation)',
            'equation': 'β-hydroxy carbonyl → α,β-unsaturated carbonyl + H₂O',
            'details': (
                'Under heating or basic conditions: elimination of water via E1cb mechanism '
                '(conjugate base) gives the conjugated α,β-unsaturated carbonyl compound. '
                'Driven by conjugation stabilization (~15-20 kJ/mol gain).'
            ),
            'requires_heat': True,
        })

        return steps

    def _predict_products(self, comp1, comp2, is_crossed):
        """Predict aldol and condensation products."""
        c1_type = comp1.get('carbonyl_type', '?') if comp1 else '?'
        c2_type = comp2.get('carbonyl_type', '?') if comp2 else c1_type

        if is_crossed:
            aldol_desc = f"β-hydroxy carbonyl from {c1_type} (enolate) + {c2_type} (electrophile)"
            cond_desc = f"α,β-unsaturated {c1_type}/{c2_type} cross-condensation product"
        else:
            aldol_desc = f"β-hydroxy-{c1_type} dimer (self-aldol)"
            cond_desc = f"α,β-unsaturated {c1_type} (self-condensation)"

        return {
            'aldol_product': {
                'type': 'β-hydroxy aldehyde/ketone',
                'description': aldol_desc,
                'functional_groups': ['carbonyl', 'secondary alcohol'],
            },
            'condensation_product': {
                'type': 'α,β-unsaturated carbonyl',
                'description': cond_desc,
                'functional_groups': ['carbonyl', 'alkene (conjugated)'],
                'driving_force': 'Conjugation (π-π) stabilization + entropy gain (loss of H₂O)',
            },
        }

    def _analyze_dehydration(self, products, conditions):
        """Analyze dehydration feasibility."""
        has_heat = 'heat' in conditions.lower() or 'Δ' in conditions or 'reflux' in conditions.lower()

        return {
            'spontaneous_at_RT': False,
            'favored_by_heat': True,
            'mechanism': 'E1cb (conjugate base elimination) under basic conditions; E1 under acid',
            'thermodynamic_drive': 'Conjugation stabilizes α,β-unsaturated product by ~15-20 kJ/mol',
            'will_dehydrate': has_heat or 'condensation' in conditions.lower(),
            'conditions_needed': 'Heat (50-100°C) or prolonged reaction time',
        }

    def _analyze_stereochemistry(self, comp1, comp2, catalyst):
        """Analyze stereoselectivity."""
        return {
            'new_stereocenters': 'Up to 2 new chiral centers (α and β carbons)',
            'diastereomers': 'syn (Z-enolate) and anti (E-enolate) possible',
            'zimmerman_traxler_model': (
                'A six-membered cyclic transition state determines stereochemistry: '
                'the enolate and carbonyl approach in a chair-like TS, with substituents '
                'preferring equatorial positions → predictable syn/anti ratio'
            ),
            'evans_model': 'Chelated (metal) control gives syn selectivity; non-chelated gives anti',
            'practical_note': (
                'Without chiral auxiliary or catalyst: racemic mixture of each diastereomer. '
                'With chiral catalyst (proline, etc.): enantioselective aldol possible.'
            ),
        }

    def _analyze_crossed(self, comp1, comp2, is_crossed):
        """Analyze crossed aldol issues."""
        if not is_crossed or not comp2:
            return {'not_applicable': True}

        c1_enol = comp1.get('can_form_enolate', False)
        c2_enol = comp2.get('can_form_enolate', False)
        c1_react = comp1.get('reactivity_as_electrophile', '')
        c2_react = comp2.get('reactivity_as_electrophile', '')

        issues = []
        if c1_enol and c2_enol:
            issues.append("Both components can enolize → mixture of 4 products (self + crossed)")

        strategy = []
        if c1_enol and not c2_enol:
            strategy.append(f"Use {comp1.get('carbonyl_type','?')} as enolate source (has α-H)")
            strategy.append(f"Use {comp2.get('carbonyl_type','?')} as exclusive electrophile (no α-H)")
            strategy.append("Clean crossed aldol — only one enolate direction possible")
        elif c2_enol and not c1_enol:
            strategy.append(f"Use {comp2.get('carbonyl_type','?')} as enolate source")
            strategy.append(f"Use {comp1.get('carbonyl_type','?')} as exclusive electrophile")
        else:
            strategy.append("Consider: use LDA for kinetic enolate of one component only")
            strategy.append("Or use pre-formed enolate (stoichiometric base at low T)")

        return {
            'potential_issue': issues[0] if issues else None,
            'recommended_strategy': strategy,
            'clean_crossed_possible': (c1_enol != c2_enol),
        }

    def _evaluate_favorability(self, comp1, comp2, catalyst, is_crossed):
        score = 3
        if comp1 and comp1.get('can_form_enolate'): score += 3
        if comp1 and comp1.get('carbonyl_type') == 'aldehyde': score += 1
        if comp2:
            if comp2.get('carbonyl_type') == 'aldehyde': score += 1
            if is_crossed and not comp2.get('can_form_enolate'): score += 2  # clean crossed
            if is_crossed and comp2.get('can_form_enolate'): score -= 1  # messy
        if 'LDA' in catalyst or 'lda' in catalyst.lower(): score += 1

        if score >= 7: return 'excellent'
        elif score >= 5: return 'good'
        elif score >= 3: return 'moderate'
        return 'possible but may need optimization'

    def _build_summary(self, comp1, comp2, products, fav, is_crossed):
        rxn = "crossed aldol" if is_crossed else "self-aldol"
        aldol = products.get('aldol_product', {}).get('description', '?')
        cond = products.get('condensation_product', {}).get('description', '?')
        return f"{rxn}: {aldol}. Condensation: {cond}. Favorability: {fav}."

    def _run_text(self, query: str) -> dict:
        parts = query.strip().split()
        c1 = parts[0] if len(parts) > 0 else ''
        c2 = parts[1] if len(parts) > 1 else ''
        cat = parts[2] if len(parts) > 2 else 'base'
        cond = parts[3] if len(parts) > 3 else 'RT then heat'
        return self._run_base(c1, c2, cat, cond)


if __name__ == "__main__":
    from ..utils.mcp_app import run_mcp_server
    run_mcp_server()
