import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HyperconjugationExplainer(BaseTool):
    """
    超共轭效应解释工具 - 解释超共轭效应对碳正离子、自由基、烯烃等体系稳定性的影响。
    分析α-C-H键数目、σ-π/σ-p轨道相互作用、稳定性排序。
    """
    __version__ = "0.1.0"
    name             = "HyperconjugationExplainer"
    func_name        = "explain_hyperconjugation"
    description      = "Explain hyperconjugation effects on molecular stability: carbocations, free radicals, alkenes, and other systems with α-C-H bonds interacting with adjacent empty/partially-filled orbitals."
    implementation_description = "Knowledge-based analysis of hyperconjugation (σC-H → π* or σC-H → p orbital donation) covering carbocation stability series, radical stability, alkene stability (Zaitsev vs Hofmann), NMR coupling constants, and conformational preferences."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Hyperconjugation", "Carbocation Stability", "Radical Stability", "Molecular Orbital Theory", "Conformational Analysis"]
    required_envs    = []

    code_input_sig   = [
        ("molecule", "str", "N/A", "Molecule or system to analyze: e.g., 'tert-butyl cation', 'propene', 'isobutylene', 'ethyl radical', 'toluene', or general terms like 'carbocation stability'."),
        ("question", "str", "", "Specific question (optional): e.g., 'why is tertiary more stable?', 'number of alpha hydrogens', 'NMR effect')."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Molecule query. Example: 'tert-butyl cation' or 'isobutylene stability'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed hyperconjugation analysis with orbital interaction diagrams, stability ranking, α-H count, and quantitative estimates."),
    ]

    examples         = [
        {
            "code_input": {"molecule": "tert-butyl cation", "question": ""},
            "text_input": {"input_params": "tert-butyl cation"},
            "output": {"result": "9 α-C-H bonds hyperconjugate..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_database()

    def _build_database(self):
        """Build hyperconjugation knowledge database."""
        
        # Carbocation stability data with hyperconjugation analysis
        self.carbocations = {
            "methyl cation": {
                "formula": "CH3+",
                "alpha_h_count": 0,
                "alpha_c_count": 0,
                "hyperconjugative_structures": 0,
                "stability_rank": 7,
                "relative_energy_kcal": "~0 (reference)",
                "description": "No α-C-H bonds available for hyperconjugation. The positive charge is fully localized on carbon. Extremely unstable — only observed in gas phase or superacid media.",
                "orbital_picture": "Empty sp² hybrid orbital on C. No filled σ(C-H) orbitals of appropriate symmetry to overlap.",
                "geometry": "Trigonal planar (sp²), ∠HCH ≈ 120deg",
            },
            "ethyl cation": {
                "formula": "CH3-CH2+",
                "alpha_h_count": 3,
                "alpha_c_count": 1,
                "hyperconjugative_structures": 3,
                "stability_rank": 6,
                "relative_energy_kcal": "+1 to +5 above t-Bu+ (very unstable)",
                "description": "Primary carbocation with 3 α-C-H bonds from the methyl group. Each C-H σ bond can donate electron density into the empty p orbital. Still very unstable — rarely observed without stabilization.",
                "orbital_picture": "Three σ(C-H) bonds of CH3 can overlap with empty p orbital on CH2+. No bond rotation needed (free rotation).",
                "geometry": "The classical structure has the empty p orbital perpendicular to the C-C bond plane; non-classical bridged structures also proposed.",
            },
            "isopropyl cation": {
                "formula": "(CH3)2CH+",
                "alpha_h_count": 6,
                "alpha_c_count": 2,
                "hyperconjugative_structures": 6,
                "stability_rank": 4,
                "relative_energy_kcal": "~5-12 above t-Bu+",
                "description": "Secondary carbocation with 6 α-C-H bonds (two methyl groups). Significantly more stable than primary due to doubled hyperconjugation. Can be observed under some conditions.",
                "orbital_picture": "Six σ(C-H) bonds can interact with the empty p orbital. Maximum overlap when C-H bonds are parallel to the p orbital (eclipsed conformations favored).",
                "note": "Preferred conformation: one C-H bond from each CH3 eclipsed with the empty p orbital for maximum overlap",
            },
            "tert-butyl cation": {
                "formula": "(CH3)3C+",
                "alpha_h_count": 9,
                "alpha_c_count": 3,
                "hyperconjugative_structures": 9,
                "stability_rank": 1,
                "relative_energy_kcal": "Most stable simple alkyl carbocation",
                "description": "Tertiary carbocation with **9 α-C-H bonds** — maximum hyperconjugation among simple alkyl carbocations. This explains why tertiary > secondary > primary > methyl in stability. Observed in stable ion pair chemistry.",
                "orbital_picture": "Nine σ(C-H) bonds from three methyl groups can donate into the empty p orbital. In the idealized geometry, one C-H bond from each CH3 is aligned parallel to the p orbital axis.",
                "geometry": "Approximately planar at cationic center; methyl groups rotate freely but prefer eclipsed orientation relative to empty p orbital",
                "key_evidence": "NMR shows equivalent methyl groups; IR/X-ray in superacids confirms structure",
            },
            "allyl cation": {
                "formula": "CH2=CH-CH2+",
                "alpha_h_count": 4,  # 2 on C3 + 2 on C1 (vinylic)
                "stability_rank": 2,
                "relative_energy_kcal": "More stable than t-Bu+ (resonance > hyperconjugation)",
                "description": "**Resonance-stabilized** — this is primarily a resonance effect (π delocalization), not just hyperconjugation. But the terminal CH2+ group still benefits from 2 α-C-H hyperconjugation. Overall more stable than t-Bu+ because π-conjugation >> σ-hyperconjugation.",
                "special_note": "Resonance dominates here; hyperconjugation is a secondary stabilizing factor",
            },
            "benzyl cation": {
                "formula": "Ph-CH2+",
                "alpha_h_count": 2,
                "stability_rank": 3,
                "relative_energy_kcal": "Comparable to allyl cation, slightly less stable",
                "description": "Resonance with benzene ring provides massive stabilization (charge delocalized over ortho/para positions + benzylic carbon). Plus 2 α-C-H hyperconjugation. One of the most stable carbocations.",
                "special_note": "Ring resonance is primary effect; hyperconjugation from benzylic CH2 is secondary",
            },
        }

        # Radical stability data
        self.radicals = {
            "methyl radical": {"alpha_h": 0, "rank": 7, "stability": "Least stable"},
            "primary radical (•CH2R)": {"alpha_h": 3, "rank": 5, "stability": "Unstable"},
            "secondary radical (•CHR2)": {"alpha_h": 6, "rank": 3, "stability": "Moderately stable"},
            "tertiary radical (•CR3)": {"alpha_h": 9, "rank": 2, "stability": "Stable"},
            "allyl radical (•CH2-CH=CH2)": {"alpha_h": 2, "rank": 1, "stability": "Very stable (resonance)"},
            "benzyl radical (Ph-CH2•)": {"alpha_h": 2, "rank": 1, "stability": "Very stable (resonance)"},
        }

        # Alkene stability (hyperconjugation with π bond)
        self.alkenes = {
            "ethene": {"alpha_h_on_sp2": 0, "substitution": "unsubstituted", "stability": "Reference (0)", "heat_of_hydrogenation": "-32.8 kcal/mol (most negative = least stable)"},
            "propene (mono-substituted)": {"alpha_h_on_sp2": 3, "substitution": "monosubstituted", "stability": "+0.0 (reference substituted)", "heat_of_hydrogenation": "-30.1 kcal/mol"},
            "1-butene (terminal)": {"alpha_h_on_sp2": 3, "substitution": "monosubstituted", "stability": "Similar to propene", "heat_of_hydrogenation": "-30.3 kcal/mol"},
            "cis-2-butene (disubstituted, cis)": {"alpha_h_on_sp2": 6, "substitution": "disubstituted (cis)", "stability": "+2.7 kcal/mol more stable than mono", "heat_of_hydrogenation": "-28.6 kcal/mol"},
            "trans-2-butene (disubstituted, trans)": {"alpha_h_on_sp2": 6, "substitution": "disubstituted (trans)", "stability": "+4.0 kcal/mol more stable than mono", "heat_of_hydrogenation": "-27.6 kcal/mol"},
            "tetramethylethylene (tetra-substituted)": {"alpha_h_on_sp2": 12, "substitution": "tetrasubstituted", "stability": "Most stable alkene", "heat_of_hydrogenation": "-26.9 kcal/mol"},
            "isobutylene (gem-disubstituted)": {"alpha_h_on_sp2": 6, "substitution": "geminal disubstituted", "stability": "Between mono and trans-di", "heat_of_hydrogenation": "-28.4 kcal/mol"},
        }

        # Special topics
        self.special_topics = {
            "anomeric_effect": {
                "title": "Anomeric Effect (Special Case of Hyperconjugation)",
                "description": "In saturated heterocycles (pyranoses, acetals), an electronegative substituent at the anomeric center prefers the axial position over equatorial — contrary to steric expectations. Caused by hyperconjugation: σ*C-O donates into σ*C-OR antibonding orbital (nO → σ* interaction also contributes).",
                "example": "α-D-glucopyranose: OMe group is axial despite steric preference for equatorial",
                "energy": "~3-5 kcal/mol stabilization for axial orientation",
            },
            "gauche_effect": {
                "title": "Gauche Effect",
                "description": "In 1,2-difluoroethane and similar molecules, the gauche conformation is more stable than anti — opposite to typical steric preferences. Hyperconjugation: σC-H → σ*C-F donation stabilizes gauche arrangement where more C-H bonds are antiperiplanar to C-F bonds.",
                "example": "F-CH2-CH2-F: gauche preferred by ~0.5-1 kcal/mol over anti",
            },
            "nmr_coupling": {
                "title": "NMR Vicinal Coupling Constants (³JHH)",
                "description": "Hyperconjugation affects J-coupling magnitudes. More hyperconjugative electron donation into a C-H σ* orbital lengthens and weakens that C-H bond, reducing coupling constant. Karplus relationship modified by substituent effects.",
                "typical_values": "sp³-sp³: 6-8 Hz (typical); electron-withdrawing groups increase J slightly; electron-donating groups can decrease it",
            },
        }

    def _find_system(self, query: str):
        """Find system in database."""
        q = query.lower().strip()
        
        # Check carbocations
        for key, data in self.carbocations.items():
            if q in key or key in q or any(word in q for word in key.split()):
                return ("carbocation", key, data)
        
        # Check radicals
        for key, data in self.radicals.items():
            if q in key or key in q or "radical" in q:
                return ("radical", key, data)
        
        # Check alkenes
        for key, data in self.alkenes.items():
            if q in key or key in q or "alkene" in q or "olefin" in q:
                return ("alkene", key, data)
        
        # Check special topics
        for key, data in self.special_topics.items():
            if q in key or key.replace("_", " ") in q:
                return ("special", key, data)

        return None

    def _run_base(self, molecule: str, question: str = "") -> str:
        """Explain hyperconjugation effects."""
        result = self._find_system(molecule)
        
        parts = [f"## Hyperconjugation Analysis: `{molecule}`\n"]
        if question:
            parts.append(f"**Question:** {question}\n")

        if result is None:
            # General overview / no specific match
            parts += self._generate_general_overview()
            return "\n".join(parts)

        sys_type, sys_name, data = result
        parts.append(f"**System Type:** {sys_type.replace('_', ' ').title()}\n")

        if sys_type == "carbocation":
            parts += self._format_carbocation(sys_name, data)
        elif sys_type == "radical":
            parts += self._format_radical(sys_name, data)
        elif sys_type == "alkene":
            parts += self._format_alkene(sys_name, data)
        elif sys_type == "special":
            parts += self._format_special(sys_name, data)

        # Always add general principles summary
        parts += self._principles_summary()

        return "\n".join(parts)

    def _format_carbocation(self, name, data):
        """Format carbocation analysis."""
        parts = []
        parts.append(f"### ⚡ {name.title()}: {data['formula']}\n")
        parts.append(f"| Property | Value |")
        parts.append(f"|---|---|")
        parts.append(f"| **α-C-H Count** | **{data['alpha_h_count']}** |")
        parts.append(f"| **Hyperconjugative Structures** | **{data['hyperconjugative_structures']}** |")
        parts.append(f"| **Stability Rank** | #{data['stability_rank']} (among common carbocations) |")
        parts.append(f"| **Relative Energy** | {data['relative_energy_kcal']} |")
        parts.append("")
        parts.append(f"#### 📖 Description\n{data['description']}\n")
        parts.append(f"#### 🔬 Orbital Picture\n{data.get('orbital_picture', 'N/A')}\n")
        if data.get('geometry'):
            parts.append(f"**Geometry:** {data['geometry']}\n")
        if data.get('key_evidence'):
            parts.append(f"**Experimental Evidence:** {data['key_evidence']}\n")
        if data.get('special_note'):
            parts.append(f"> 💡 **Note:** {data['special_note']}\n")

        # Show all carbocations ranking table
        parts.append("\n### 📊 Complete Carbocation Stability Ranking\n")
        parts.append("| Rank | Carbocation | α-C-H | Hyp. Struct. | Relative Stability |")
        parts.append("|---|---|---|---|---|")
        ranked = sorted(self.carbocations.items(), key=lambda x: x[1]["stability_rank"])
        for r, (k, d) in enumerate(ranked, 1):
            marker = " ← current" if k == name else ""
            parts.append(f"| {r} | {k} ({d['formula']}) | {d['alpha_h_count']} | {d.get('hyperconjugative_structures', 'N/A')} | {'★★★ Most Stable' if d['stability_rank']==1 else '★★' if d['stability_rank']<=3 else '★'} |{marker}")
        
        return parts

    def _format_radical(self, name, data):
        """Format radical analysis."""
        parts = []
        parts.append(f"### ⚛️ {name.title()}\n")
        parts.append(f"- **α-C-H Bonds Available:** {data['alpha_h']}")
        parts.append(f"- **Stability Rank:** #{data['rank']}")
        parts.append(f"- **Stability:** {data['stability']}")
        parts.append("")
        parts.append("Radicals follow the same trend as carbocations: **tertiary > secondary > primary > methyl**, driven by hyperconjugation of α-C-H σ bonds into the half-filled p orbital.\n")
        
        # Full ranking
        parts.append("| Radical | α-C-H | Rank | Stability |")
        parts.append("|---|---|---|---|")
        ranked = sorted(self.radicals.items(), key=lambda x: x[1]["rank"])
        for k, d in ranked:
            marker = " ← current" if k.split()[0] in name else ""
            parts.append(f"| {k} | {d['alpha_h']} | #{d['rank']} | {d['stability']} |{marker}")
        return parts

    def _format_alkene(self, name, data):
        """Format alkene stability analysis."""
        parts = []
        parts.append(f"### 🔄 {name.title()}\n")
        parts.append(f"- **Substitution:** {data['substitution']}")
        parts.append(f"- **α-C-H on sp² Carbons:** {data['alpha_h_on_sp2']}")
        parts.append(f"- **Heat of Hydrogenation:** {data['heat_of_hydrogenation']}")
        parts.append(f"- **Stability:** {data['stability']}")
        parts.append("")
        parts.append("**Why more substituted alkenes are more stable:**\n")
        parts.append("- Each alkyl group attached to sp² carbon can hyperconjugate: σ(C-H) → π*(C=C)\n")
        parts.append("- More alkyl groups = more hyperconjugation = greater stabilization\n")
        parts.append("- Also: alkyl groups are electron-donating (+I), which enriches the π bond electron density\n")
        
        # Full comparison table
        parts.append("\n### 📊 Alkene Stability Comparison (Heat of Hydrogenation)\n")
        parts.append("| Alkene | Substitution | ΔHhydrog (kcal/mol) | Relative Stability |")
        parts.append("|---|---|---|---|")
        def _alkene_sort_key(item):
            try:
                hoh = item[1]['heat_of_hydrogenation']
                val = hoh.split('(')[1].split()[0].replace(')','')
                return float(val)
            except (ValueError, IndexError):
                return 0
        sorted_alkenes = sorted(self.alkenes.items(), key=_alkene_sort_key)
        for k, d in sorted_alkenes:
            marker = " ← current" if k == name else ""
            parts.append(f"| {k} | {d['substitution']} | {d['heat_of_hydrogenation']} | {d['stability']} |{marker}")
        parts.append("\n> *More negative ΔHhydrog = less stable alkene (more heat released upon hydrogenation)*\n")
        return parts

    def _format_special(self, name, data):
        """Format special topic."""
        parts = []
        parts.append(f"### 🌟 {data['title']}\n")
        parts.append(f"{data['description']}\n")
        if data.get('example'):
            parts.append(f"**Example:** {data['example']}\n")
        if data.get('energy'):
            parts.append(f"**Magnitude:** {data['energy']}\n")
        return parts

    def _generate_general_overview(self):
        """Generate general hyperconjugation overview."""
        return [
            "### 📚 Hyperconjugation: General Overview\n",
            "**Definition:** Hyperconjugation is the stabilizing interaction between a filled σ-bonding orbital (usually C-H) and an adjacent empty or partially-filled orbital (p orbital, π*, or σ*).\n",
            "#### Core Concept\n",
            "```\n     H             H\n     |             |\n  H-C-C:  ↔  H-C=C-H  (no-bond resonance form)\n     |    ↑        ||\n     R    p         R\n",
            "  (σC-H donates electron density into empty p orbital)\n",
            "```\n",
            "\n#### Key Rules\n",
            "| Rule | Explanation |",
            "|---|---|",
            "| **More α-C-H = more stable** | Each C-H bond can contribute ~5-10 kcal/mol of stabilization |",
            "| **Order of stability (cations)** | (CH3)3C+ > (CH3)2CH+ > CH3CH2+ > CH3+ |",
            "| **Order of stability (radicals)** | (CH3)3C• > (CH3)2CH• > CH3CH2• > CH3• |",
            "| **Alkenes:** tetra > tri > di > mono > unsubstituted | Same principle: σ(C-H) → π*(C=C) |",
            "| **Bond lengthening** | C-H bonds involved in hyperconjugation are slightly elongated |",
            "| **Conformational preference** | Eclipsed/syn-periplanar arrangements maximize overlap |",
            "",
            "#### Quantitative Impact (Approximate)\n",
            "- Each α-C-H hyperconjugative interaction: **~5-12 kcal/mol** stabilization for carbocations\n",
            "- For radicals: **~4-8 kcal/mol** per α-C-H\n",
            "- For alkenes: **~2-3 kcal/mol** per additional alkyl substituent\n",
            "",
            "### Available Systems to Analyze\n",
            "**Carbocations:** " + ", ".join(self.carbocations.keys()) + "\n",
            "**Alkenes:** " + ", ".join(list(self.alkenes.keys())[:5]) + ", ...\n",
            "**Special Topics:** " + ", ".join(self.special_topics.keys()) + "\n",
        ]

    def _principles_summary(self):
        """Return principles summary section."""
        return [
            "\n---\n",
            "### 🧠 Summary: How Hyperconjugation Works\n",
            "1. **Orbital Requirement:** A filled σ orbital (donor) adjacent to an empty/partially-filled orbital (acceptor)\n",
            "2. **Geometry Matters:** Maximum overlap when σ bond is **parallel** (or anti-periplanar) to the acceptor orbital\n",
            "3. **Not Resonance:** Unlike resonance, hyperconjugation involves σ bonds (not π bonds) and is generally weaker\n",
            "4. **Evidence:** \n",
            "   - Bond lengthening of participating C-H bonds (electron diffraction)\n",
            "   - NMR coupling constant changes\n",
            "   - Conformational preferences (anomeric effect, gauche effect)\n",
            "   - Isotope effects (C-D bonds are shorter/stronger → less hyperconjugation)\n",
            "5. **Practical Consequences:**\n",
            "   - Zaitsev's rule (more substituted alkene favored in elimination)\n",
            "   - Markovnikov's rule (more stable carbocation intermediate)\n",
            "   - Carbocation rearrangements (hydride/alkyl shifts toward more hyperconjugation)\n",
            "   - Toluene acidity (benzylic C-H weakened by hyperconjugation with ring π system)\n",
        ]

    def _run_text(self, input_params: str) -> str:
        input_params = input_params.strip()
        if not input_params:
            raise ChemMCPError("Please provide a molecule or system name. Example: 'tert-butyl cation', 'propene'")
        return self._run_base(input_params)
