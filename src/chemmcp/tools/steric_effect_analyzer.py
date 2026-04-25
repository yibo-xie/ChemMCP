import math
import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class StericEffectAnalyzer(BaseTool):
    """
    位阻效应分析工具 - 分析分子中的位阻效应对反应活性、选择性和构象稳定性的影响。
    包含A值（环己烷）、锥角（配体）、Tolman参数等定量数据。
    """
    __version__ = "0.1.0"
    name             = "StericEffectAnalyzer"
    func_name        = "analyze_steric_effect"
    description      = "Analyze steric (steric hindrance) effects on reaction rates, selectivity, and conformational stability. Includes A-values for cyclohexane, cone angles for ligands, and Tolman steric parameters."
    implementation_description = "Knowledge-based system using experimental steric parameters: A-values (cyclohexane conformational preferences), cone angles (ligand bulk), Tolman parameters, Taft Es values, and buried volume (%Vbur) to predict steric effects on reactions."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Steric Effects", "Conformational Analysis", "Ligand Design", "Cyclohexane", "A-Values"]
    required_envs    = []

    code_input_sig   = [
        ("molecule", "str", "N/A", "Molecule or substituent name to analyze: e.g., 'tert-butylcyclohexane', 'PPh3 ligand', 'isopropyl group', 'mesityl', or general terms like 'cyclohexane A-values'."),
        ("reaction_context", "str", "", "Optional reaction context: 'SN2', 'E2', 'catalysis', 'conformation', 'ligand binding'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Query string. Example: 'tert-butyl SN2' or 'PPh3 cone angle' or 'cyclohexane equatorial preference'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed steric effect analysis with quantitative parameters, predicted impact on reaction, and comparison data."),
    ]

    examples         = [
        {
            "code_input": {"molecule": "tert-butyl cyclohexane", "reaction_context": "conformation"},
            "text_input": {"input_params": "tert-butyl cyclohexane"},
            "output": {"result": "A-value analysis... strong preference for equatorial..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_database()

    def _build_database(self):
        """Build comprehensive steric parameter database."""
        
        # === A-Values for Cyclohexane (kcal/mol) ===
        # ΔGdeg for equilibrium: axial ↔ equatorial (positive = prefers equatorial)
        self.a_values = {
            # Very large groups (always equatorial)
            "tert-butyl": {"a_value": 4.9, "category": "very_large", "always_equatorial": True,
                           "description": "Extremely bulky; essentially locks the ring in one conformation with t-Bu equatorial. Used as a 'conformational lock'."},
            "triphenylmethyl (trityl)": {"a_value": "~5.5", "category": "very_large", "always_equatorial": True,
                                         "description": "Even larger than t-Bu; three phenyl rings create massive steric clash in axial position."},
            "trimethylsilyl (TMS)": {"a_value": 2.5, "category": "large", "always_equatorial": False,
                                      "description": "Large but not as extreme as t-Bu. Si-C bond longer than C-C reduces some 1,3-diaxial interactions."},
            
            # Large groups (strongly prefer equatorial)
            "isopropyl": {"a_value": 2.21, "category": "large", "always_equatorial": False,
                          "description": "Two methyl groups create significant 1,3-diaxial repulsion when axial."},
            "phenyl": {"a_value": 3.00, "category": "large", "always_equatorial": False,
                        "description": "Phenyl is larger than isopropyl due to flat shape and π-system extending into axial space."},
            "cyclohexyl": {"a_value": "~3.0", "category": "large", "always_equatorial": False,
                             "description": "Second ring creates substantial steric demand."},
            "carboxylic acid / ester": {"a_value": 0.70, "category": "medium", "always_equatorial": False,
                                          "description": "COOH/COOR has moderate size; can be axial in some cases."},
            "CN (cyano)": {"a_value": 0.23, "category": "small", "always_equatorial": False,
                            "description": "Linear group — small cross-section despite moderate length."},
            
            # Medium groups
            "methyl": {"a_value": 1.80, "category": "medium", "always_equatorial": False,
                       "description": "The reference medium-sized group. ~95% equatorial at room temperature."},
            "ethyl": {"a_value": 1.79, "category": "medium", "always_equatorial": False,
                      "description": "Nearly identical to methyl (gauche interactions similar)."},
            "F": {"a_value": 0.25, "category": "small", "always_equatorial": False,
                   "description": "Very small — F is small despite high electronegativity. Bond length matters more than atomic size here."},
            "Cl": {"a_value": 0.55, "category": "small-medium", "always_equatorial": False,
                    "description": "Larger than F due to larger atomic radius and longer bond."},
            "Br": {"a_value": 0.48, "category": "small-medium", "always_equatorial": False,
                    "description": "Similar to Cl; longer bond reduces 1,3-diaxial somewhat."},
            "I": {"a_value": 0.43, "category": "small-medium", "always_equatorial": False,
                   "description": "Longest bond among halogens → smallest effective size in 1,3-diaxial context."},
            "OH": {"a_value": 0.87, "category": "medium-small", "always_equatorial": False,
                    "description": "Moderate — H-bonding can override steric preference in protic solvents."},
            "OCH3": {"a_value": 0.60, "category": "small-medium", "always_equatorial": False,
                      "description": "Smaller than expected because C-O bond rotation allows CH3 to avoid some clashes."},
            "NH2": {"a_value": 1.2, "category": "medium", "always_equatorial": False,
                     "description": "N is smaller than O (less 1,3-diaxial repulsion), but NH2 pyramidal shape adds complexity."},
            "NO2": {"a_value": 1.15, "category": "medium", "always_equatorial": False,
                     "description": "Planar group — moderate steric demand from two oxygen atoms."},
            "CH=CH2 (vinyl)": {"a_value": 1.35, "category": "medium", "always_equatorial": False,
                                 "description": "sp² carbon creates a planar group of moderate bulk."},
            "C≡CH (ethynyl)": {"a_value": 0.50, "category": "small", "always_equatorial": False,
                                  "description": "Linear — minimal steric profile despite length."},
            "COCH3 (acetyl)": {"a_value": 0.70, "category": "medium-small", "always_equatorial": False,
                                 "description": "Similar to COOH/ester."},
        }

        # === Cone Angles (degrees) for Ligands ===
        # Tolman cone angle: measure of ligand bulk in metal complexes
        self.cone_angles = {
            # Very small ligands
            "H": {"cone_angle": 0, "category": "tiny", "notes": "Reference point"},
            "CO (carbonyl)": {"cone_angle": 95, "category": "small", "notes": "Small linear ligand; common in organometallics"},
            "PMe3 (trimethylphosphine)": {"cone_angle": 118, "category": "small-medium",
                                           "notes": "Small phosphine; very basic, good σ-donor"},
            "P(OMe)3 (trimethyl phosphite)": {"cone_angle": 107, "category": "small",
                                                "notes": "Small π-acceptor ligand"},
            
            # Medium ligands
            "PPh3 (triphenylphosphine)": {"cone_angle": 145, "category": "large",
                                           "notes": "Most common phosphine ligand; moderately bulky"},
            "AsPh3 (triphenylarsine)": {"cone_angle": 145, "category": "large",
                                         "notes": "Similar to PPh3 but softer donor"},
            "PEt3 (triethylphosphine)": {"cone_angle": 132, "category": "medium",
                                          "notes": "More flexible than PPh3; alkyl chains can rotate away"},
            "PCy3 (tricyclohexylphosphine)": {"cone_angle": 170, "category": "large",
                                                "notes": "Bulky, strong electron donor; widely used in catalysis"},
            
            # Large ligands
            "P(o-tol)3 (tri-o-tolylphosphine)": {"cone_angle": 194, "category": "very_large",
                                                   "notes": "Ortho-methyl groups add significant bulk near metal center"},
            "P(mes)3 (trimesitylphosphine)": {"cone_angle": 212, "category": "very_large",
                                               "notes": "Extremely bulky; rarely binds due to steric congestion"},
            "P(t-Bu)3 (tri-tert-butylphosphine)": {"cone_angle": 182, "category": "very_large",
                                                      "notes": "Very bulky electron-rich phosphine; used in Buchwald-Hartwig, Pd-catalyzed couplings"},
            "P(i-Pr)3 (triisopropylphosphine)": {"cone_angle": 160, "category": "large",
                                                    "notes": "Popular bulky electron-donor for cross-coupling catalysts"},
            "NHC (N-heterocyclic carbene, IMes)": {"cone_angle": 185, "category": "very_large",
                                                     "notes": "Bulky NHC ligand; very strong σ-donor"},
            "NHC (SIPr)": {"cone_angle": 192, "category": "very_large",
                             "notes": "Even bulkier NHC with isopropyl substituents"},
            "Cp (cyclopentadienyl)": {"cone_angle": "~120 (variable)", "category": "medium",
                                        "notes": "η⁵-coordination; effective size depends on substitution"},
        }

        # === Taft Steric Parameters (Es) ===
        # More negative = more sterically hindered
        self.taft_es = {
            "H": 0.00, "Me": 0.00, "Et": -0.07, "i-Pr": -0.47, "t-Bu": -1.54,
            "CH2Ph": -0.38, "CH2SiMe3": -0.41, "CF3": -0.90,
            "F": 0.27, "Cl": 0.24, "Br": 0.27, "I": 0.37,
            "Ome": 0.99, "NEt2": -0.07,
        }

        # === Reaction-specific steric effects ===
        self.reaction_effects = {
            "SN2": {
                "description": "SN2 reactions are extremely sensitive to steric hindrance at the electrophilic carbon.",
                "rate_order": "CH3-X > primary > secondary >>> tertiary (essentially no SN3)",
                "relative_rates": {"methyl": 1.0, "primary": 1.0, "secondary": 0.01-0.001, "tertiary": "~10^-6"},
                "explanation": "Backside attack requires clear access to the σ* orbital. Bulky groups block this trajectory.",
                "neighboring_group_effect": "β-branching also slows SN2 (Thorpe-Ingold effect: angle compression favors cyclization over intermolecular attack)",
            },
            "E2": {
                "description": "E2 elimination is less sensitive than SN2 but still affected by substrate sterics.",
                "regioselectivity": "Bulky bases favor Hofmann (less substituted alkene) over Zaitsev product due to difficulty accessing more hindered β-hydrogen.",
                "base_sterics": "t-BuOK (bulky) → Hofmann; EtONa (smaller) → Zaitsev",
                "explanation": "Base must abstract a β-hydrogen; bulky base cannot easily reach the more hindered (but thermodynamically favored) β-position.",
            },
            "catalysis": {
                "description": "In transition metal catalysis, ligand sterics control selectivity and activity.",
                "effects": [
                    "Bulky ligands → favor monomolecular oxidative addition (dissociative pathway)",
                    "Bulky ligands → increase selectivity (block certain approach trajectories)",
                    "Too bulky → may inhibit reactivity entirely (cannot form active complex)",
                    "Optimal steric bulk → balance between activity and selectivity",
                ],
                "examples": "Buchwald ligands (SPhos, XPhos): tuned cone angles for optimal coupling efficiency",
            },
            "conformation": {
                "description": "Sterics determine preferred conformations of acyclic and cyclic molecules.",
                "principles": [
                    "Cyclohexane: equatorial preference governed by A-values",
                    "Acyclic: staggered > eclipsed; anti > gauche (for most cases)",
                    "Butane: anti favored by ~0.9 kcal/mol over gauche (steric repulsion between methyls)",
                    "Geminal dimethyl effect: stabilizes gauche conformation (Thorpe-Ingold)",
                ],
            },
        }

    def _run_base(self, molecule: str, reaction_context: str = "") -> str:
        """Analyze steric effects."""
        mol_lower = molecule.lower().strip()
        ctx_lower = reaction_context.lower() if reaction_context else ""
        
        parts = [f"## Steric Effect Analysis: `{molecule}`\n"]
        if reaction_context:
            parts.append(f"**Reaction Context:** {reaction_context}\n")

        # Search A-values
        a_match = None
        for key, data in self.a_values.items():
            if key in mol_lower or mol_lower in key or any(w in mol_lower for w in key.replace('(', '').replace(')', '').split()):
                a_match = (key, data)
                break
        
        # Search cone angles
        ca_match = None
        for key, data in self.cone_angles.items():
            if key in mol_lower or mol_lower in key:
                ca_match = (key, data)
                break

        # Search Taft Es
        es_match = None
        for key, val in self.taft_es.items():
            if key.lower() in mol_lower or mol_lower in key.lower():
                es_match = (key, val)
                break

        found_any = a_match or ca_match or es_match

        if not found_any:
            # General overview mode
            parts.append("### ⚠️ No Direct Match — General Reference\n")
            parts += self._generate_reference_overview(ctx_lower)
            return "\n".join(parts)

        # Display matches
        if a_match:
            key, data = a_match
            parts.append(f"### 📐 Cyclohexane A-Value Analysis: **{key}**\n")
            parts.append(f"| Parameter | Value |")
            parts.append(f"|---|---|")
            parts.append(f"| **A-Value (ΔGdeg)** | **{data['a_value']} kcal/mol** |")
            parts.append(f"| **Category** | {data['category'].replace('_', ' ').title()} |")
            parts.append(f"| **Always Equatorial?** | {'Yes ✅ (conformational lock)' if data['always_equatorial'] else 'No (equilibrium mixture)'} |")
            parts.append("")
            parts.append(f"> {data['description']}\n")
            
            # Calculate equilibrium composition
            av = float(str(data['a_value']).replace('~', ''))
            keq = math.exp(av / (0.001987 * 298))  # Keq = exp(-ΔG/RT), but ΔG = -RT ln K for eq→ax
            # Actually: ΔGdeg = -RT ln K where K = [eq]/[ax]
            pct_eq = 100 * av / (av + math.exp(-av / (0.001987 * 298)) * 100) if av > 0 else 50
            # Simplified: use Boltzmann distribution
            try:
                pct_eq = 100 * math.exp(av / (0.001987 * 298)) / (1 + math.exp(av / (0.001987 * 298)))
            except:
                pct_eq = 95 if av > 2 else 80 if av > 1 else 60
            
            parts.append(f"#### 📊 Equilibrium Composition at 25degC\n")
            parts.append(f"- **Equatorial:** ~{pct_eq:.1f}%")
            parts.append(f"- **Axial:** ~{100-pct_eq:.1f}%\n")

        if ca_match:
            key, data = ca_match
            parts.append(f"### 🔺 Cone Angle Analysis: **{key}**\n")
            parts.append(f"| Parameter | Value |")
            parts.append(f"|---|---|")
            parts.append(f"| **Cone Angle** | **{data['cone_angle']}deg** |")
            parts.append(f"| **Category** | {data['category'].replace('_', ' ').title()} |")
            parts.append(f"| **Notes** | {data['notes']} |\n")
            
            # Compare with other ligands
            parts.append("#### Comparison with Other Ligands\n")
            parts.append("| Ligand | Cone Angle | Category |")
            parts.append("|---|---|---|")
            def _cone_sort_key(item):
                try:
                    val = item[1]['cone_angle']
                    if isinstance(val, (int, float)):
                        return float(val)
                    s = str(val).split('(')[0].strip().replace('~','')
                    return float(s) if s.replace('.','').replace('-','').isdigit() else 999
                except:
                    return 999
            sorted_ca = sorted(self.cone_angles.items(), key=_cone_sort_key)
            for k, d in sorted_ca:
                marker = " ← current" if k == key else ""
                _ca = d['cone_angle']
                _ca_num = float(str(_ca).split('(')[0].strip().replace('~','')) if not isinstance(_ca, (int, float)) else _ca
                cat_emoji = "🔴" if _ca_num > 180 else "🟠" if _ca_num > 140 else "🟢" if _ca_num < 120 else "🟡"
                parts.append(f"| {k} | {d['cone_angle']}deg | {cat_emoji} {d['category']} |{marker}")

        if es_match:
            key, val = es_match
            parts.append(f"\n### 📏 Taft Steric Parameter (Es): **{key}**\n")
            parts.append(f"- **Es = {val:+.2f}")
            interpretation = "more sterically hindered" if val < -0.3 else "similar to methyl" if abs(val) < 0.2 else "less hindered than methyl"
            parts.append(f"- **Interpretation:** {interpretation}")

        # Reaction context analysis
        if ctx_lower:
            parts.append("\n### ⚗️ Steric Effects in Context: `{reaction_context}`\n")
            if ctx_lower in self.reaction_effects:
                rex = self.reaction_effects[ctx_lower]
                parts.append(rex["description"] + "\n")
                if "explanation" in rex:
                    parts.append(f"**Mechanism:** {rex['explanation']}\n")
                if "rate_order" in rex:
                    parts.append(f"**Rate Order:** {rex['rate_order']}\n")
                if "relative_rates" in rex:
                    parts.append("**Relative Rates:**\n")
                    for k, v in rex["relative_rates"].items():
                        parts.append(f"  - {k}: {v}\n")
                if "regioselectivity" in rex:
                    parts.append(f"**Regioselectivity:** {rex['regioselectivity']}\n")
                if "effects" in rex:
                    for e in rex["effects"]:
                        parts.append(f"- {e}\n")
            elif "sn2" in ctx_lower:
                parts.append("SN2 reactions require backside attack → severe steric penalty for hindered substrates.\n")
                parts.append("- Methyl and primary substrates react fastest\n")
                parts.append("- Secondary: 10-100× slower\n")
                parts.append("- Tertiary: essentially no SN2 (E2 dominates instead)\n")
            elif "elimination" in ctx_lower or "e2" in ctx_lower:
                parts.append("E2 elimination: bulky bases favor less substituted (Hofmann) alkene products.\n")
                parts.append("- Small base (EtO-, HO-) → Zaitsev product (more stable alkene)\n")
                parts.append("- Bulky base (t-BuO-, LDA) → Hofmann product (less hindered β-H abstraction)\n")

        return "\n".join(parts)

    def _generate_reference_overview(self, ctx):
        """Generate reference tables."""
        parts = ["### 📚 Cyclohexane A-Value Reference Table\n"]
        parts.append("| Substituent | A-Value (kcal/mol) | % Equatorial at 25degC | Category |")
        parts.append("|---|---|---|---|")
        sorted_av = sorted(self.a_values.items(), key=lambda x: float(str(x[1]['a_value']).replace('~','')), reverse=True)
        for k, d in sorted_av:
            av = float(str(d['a_value']).replace('~',''))
            try:
                pct = 100 * math.exp(av / (0.592)) / (1 + math.exp(av / (0.592)))
            except:
                pct = 99 if av > 3 else 95 if av > 2 else 75 if av > 1 else 60
            lock = " 🔒" if d.get("always_equatorial") else ""
            parts.append(f"| {k} | {d['a_value']} | ~{pct:.0f}% | {d['category']}{lock} |")

        parts.append("\n### 🔺 Ligand Cone Angle Reference\n")
        parts.append("| Ligand | Cone Angle (deg) | Use Case |")
        parts.append("|---|---|---|")
        def _cone_sort_key(item):
            try:
                val = item[1]['cone_angle']
                if isinstance(val, (int, float)):
                    return float(val)
                s = str(val).split('(')[0].strip().replace('~','')
                return float(s) if s.replace('.','').replace('-','').isdigit() else 999
            except:
                return 999
        sorted_ca = sorted(self.cone_angles.items(), key=_cone_sort_key)
        for k, d in sorted_ca:
            parts.append(f"| {k} | {d['cone_angle']} | {d['notes']} |")

        return parts

    def _run_text(self, input_params: str) -> str:
        input_params = input_params.strip()
        if not input_params:
            raise ChemMCPError("Please provide a molecule or query. Example: 'tert-butyl cyclohexane', 'PPh3 ligand'")
        parts = input_params.split()
        molecule = parts[0]
        context = " ".join(parts[1:]) if len(parts) > 1 else ""
        return self._run_base(molecule, context)
