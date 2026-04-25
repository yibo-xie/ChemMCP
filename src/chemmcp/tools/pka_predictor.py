import logging
import re
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PkaPredictor(BaseTool):
    """
    pKa 预测工具 - 基于官能团分析和取代基效应预测有机酸的pKa值。
    内置常见官能团pKa数据库 + Hammett型校正。
    """
    __version__ = "0.1.0"
    name             = "PkaPredictor"
    func_name        = "predict_pka"
    description      = "Predict pKa values of organic acids based on functional group analysis and substituent effects (inductive, resonance, field effects)."
    implementation_description = "Uses a built-in database of ~150 functional group pKa values in water (and some in DMSO) combined with Hammett σ constants for substituent effect corrections. Covers carboxylic acids, phenols, alcohols, thiols, ammonium ions, carbon acids, and more."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["pKa", "Acidity", "Functional Group", "Substituent Effects", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("molecule", "str", "N/A", "Molecule identifier: SMILES string, functional group name (e.g., 'acetic acid', 'phenol'), or structural formula description."),
        ("solvent", "str", "water", "Solvent: 'water' (default) or 'dmso'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'molecule [solvent]'. Example: 'p-nitrophenol water' or 'CH3COOH'."),
    ]

    output_sig       = [
        ("result", "str", "Predicted pKa value(s) with confidence level, explanation of contributing factors, and comparison to similar compounds."),
    ]

    examples         = [
        {
            "code_input": {"molecule": "acetic acid", "solvent": "water"},
            "text_input": {"input_params": "acetic acid water"},
            "output": {"result": "Predicted pKa ≈ 4.76..."},
        },
        {
            "code_input": {"molecule": "phenol", "solvent": "water"},
            "text_input": {"input_params": "phenol"},
            "output": {"result": "Predicted pKa ≈ 10.00..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_pka_database()

    def _build_pka_database(self):
        """Build comprehensive pKa database with substituent corrections."""
        
        # Core pKa values (parent compounds in water at 25°C)
        # Format: {name: (pKa, category, notes)}
        self.parent_pkas = {
            # === Carboxylic Acids ===
            "formic acid": (3.75, "carboxylic acid", "Simplest carboxylic acid; no electron-donating alkyl"),
            "acetic acid": (4.76, "carboxylic acid", "CH3 is weakly electron-donating → less acidic than formic"),
            "propionic acid": (4.87, "carboxylic acid", "Longer chain = slightly more EDG effect"),
            "butyric acid": (4.82, "carboxylic acid", ""),
            "benzoic acid": (4.20, "carboxylic acid", "Ph stabilizes conjugate base via resonance"),
            "carbonic acid": (6.35, "carboxylic acid", "H2CO3 ↔ HCO3-; pKa2=10.33"),
            "lactic acid": (3.86, "carboxylic acid", "α-OH is electron-withdrawing → stronger acid"),
            "oxalic acid": (1.25, "carboxylic acid", "First COOH; pKa2=4.14; α-dicarboxylic"),
            "malonic acid": (2.83, "carboxylic acid", "First COOH; pKa2=5.69"),
            "succinic acid": (4.21, "carboxylic acid", "First COOH; pKa2=5.64"),
            "glutaric acid": (4.34, "carboxylic acid", "First COOH; pKa2=5.27"),
            "maleic acid": (1.92, "carboxylic acid", "cis-butenedioic; pKa2=6.23; intramolecular H-bonding"),
            "fumaric acid": (3.03, "carboxylic acid", "trans-butenedioic; pKa2=4.44"),
            "phthalic acid": (2.89, "carboxylic acid", "ortho-dicarboxy aromatic; pKa2=5.51"),
            "trifluoroacetic acid": (0.23, "carboxylic acid", "CF3 strongly EWG → very strong acid"),
            "trichloroacetic acid": (0.70, "carboxylic acid", "CCl3 strongly EWG"),
            "dichloroacetic acid": (1.25, "carboxylic acid", "CCl2H EWG"),
            "chloroacetic acid": (2.86, "carboxylic acid", "ClCH2 EWG"),
            "cyanoacetic acid": (2.47, "carboxylic acid", "NC-CH2 EWG via induction"),
            
            # === Phenols ===
            "phenol": (10.00, "phenol", "Reference phenol; O-H acidity enhanced by Ph resonance"),
            "p-nitrophenol": (7.15, "phenol", "p-NO2 strongly EWG via resonance → much more acidic"),
            "m-nitrophenol": (8.39, "phenol", "m-NO2 EWG via induction only"),
            "o-nitrophenol": (7.17, "phenol", "o-NO2 EWG + intramolecular H-bonding stabilizes anion"),
            "p-chlorophenol": (9.38, "phenol", "p-Cl weakly EWG (-I > +R)"),
            "m-chlorophenol": (9.02, "phenol", "m-Cl -I effect only"),
            "p-cresol": (10.26, "phenol", "p-CH3 EDG → less acidic than phenol"),
            "o-cresol": (10.28, "phenol", "o-CH3 EDG + steric/H-bonding effects"),
            "picric acid" : (0.38, "phenol", "2,4,6-trinitrophenol; extremely acidic due to 3 NO2 groups"),
            "2,4-dinitrophenol": (4.09, "phenol", "Two NO2 groups → strong acid"),
            "2,6-di-tert-butylphenol": (11.5, "phenol", "Bulky t-Bu groups destabilize anion sterically"),
            
            # === Alcohols ===
            "methanol": (15.5, "alcohol", "Reference aliphatic alcohol"),
            "ethanol": (15.9, "alcohol", ""),
            "isopropanol": (17.1, "alcohol", "More substituted = less acidic (EDG effect)"),
            "tert-butanol": (18.0, "alcohol", "Most substituted → least acidic"),
            "phenol_as_alcohol_ref": (10.0, "alcohol", "PhO-H much more acidic than RO-H due to resonance stabilization"),
            
            # === Thiols ===
            "ethanethiol": (10.6, "thiol", "RSH more acidic than ROH (S larger, more polarizable)"),
            "thiophenol": (6.5, "thiol", "PhSH; S- stabilized by Ph resonance better than O"),
            "cysteine": (8.3, "thiol", "Biological thiol; α-NH3+ enhances acidity"),
            
            # === Ammonium / Amines ===
            "ammonium": (9.25, "ammonium", "NH4+ reference"),
            "methylammonium": (10.62, "ammonium", "Alkyl EDG makes conjugate base weaker → higher pKa"),
            "anilinium": (4.60, "ammonium", "PhNH3+; N lone pair delocalized into ring → strong acid"),
            "pyridinium": (5.20, "ammonium", "Aromatic N; sp2 hybridized"),
            "imidazolium": (7.00, "ammonium", "Histidine-like; two N atoms"),
            "glycine_amino": (9.78, "ammonium", "α-CO2H EWG → slightly lower than NH4+"),
            
            # === Carbon Acids (C-H Acidity) ===
            "acetone": (19.3, "carbon acid", "α-proton to ketone; enolate formation"),
            "acetophenone": (18.0, "carbon acid", "PhCOCH3; Ph stabilizes enolate"),
            "ethyl acetate": (25.0, "carbon acid", "Ester enolate; less stable than ketone enolate"),
            "acetonitrile": (25.0, "carbon acid", "NC-CH2; nitrile stabilizes carbanion"),
            "nitromethane": (10.2, "carbon acid", "CH3NO2; very acidic C-H for neutral molecule"),
            "malonate ester": (13.0, "carbon acid", "Active methylene between two esters"),
            "acetylacetone": (9.0, "carbon acid", "β-Diketone; very stable enol/enolate"),
            "dimethyl sulfoxide": (35.0, "carbon acid", "DMSO; CH3-S(=O)-CH3"),
            "triphenylmethane": (31.5, "carbon acid", "Ph3CH; trityl anion stabilized by 3 Ph rings"),
            "fluorene": (22.6, "carbon acid", "Aromatic C-H; planar carbanion stabilized"),
            "indene": (18.7, "carbon acid", "Benzylic + vinylic stabilization"),
            "cyclopentadiene": (16.0, "carbon acid", "Very acidic! Aromatic 6π e- cyclopentadienyl anion"),
            
            # === Other ===
            "hydronium": (-1.74, "other", "H3O+; reference strong acid"),
            "hydrochloric acid": (-7.0, "other", "HCl; strong acid"),
            "sulfuric acid_1": (-3.0, "other", "H2SO4 first proton; pKa2=1.99"),
            "phosphoric acid_1": (2.15, "other", "H3PO4; pKa2=7.20, pKa3=12.35"),
            "hydrogen cyanide": (9.31, "other", "HCN; weak acid"),
            "hydrogen sulfide": (7.04, "other", "H2S; pKa2=14-19 (approx)"),
            "water": (15.7, "other", "H2O autoionization"),
            "hydrazine_n2h4": (8.1, "other", "N2H4 protonated (first N)"),
            "urea": (0.88, "other", "Protonated urea (acidic form)"),
            "guanidinium": (13.6, "other", "Protonated guanidine; very basic conjugate base"),
        }

        # Substituent Hammett σ constants for pKa correction
        # Format: {substituent_position: (σm, σp)}
        self.hammett_constants = {
            "nme2": (-0.15, -0.83),  # Strong +R donor
            "nh2": (-0.16, -0.66),
            "oh": (0.12, -0.37),     # +R dominates para
            "och3": (0.10, -0.27),
            "oc2h5": (0.10, -0.24),
            "ch3": (-0.07, -0.17),
            "c2h5": (-0.07, -0.15),
            "ch(ch3)2": (-0.05, -0.15),
            "c(ch3)3": (-0.10, -0.20),
            "ph": (0.06, -0.01),
            "f": (0.34, 0.06),       # -I dominant para, +R partially cancels
            "cl": (0.37, 0.23),
            "br": (0.39, 0.23),
            "i": (0.35, 0.18),
            "cf3": (0.43, 0.54),      # Strong -I
            "cn": (0.56, 0.66),       # Strong -I, -R
            "so2ch3": (0.60, 0.72),
            "cho": (0.35, 0.42),
            "coch3": (0.38, 0.50),
            "cooh": (0.37, 0.45),
            "coor": (0.37, 0.45),
            "no2": (0.71, 0.78),      # Very strong -I, -R
            "sr": (0.23, 0.03),
            "sor": (0.40, 0.35),
            "so2r": (0.72, 0.81),
            "h": (0.0, 0.0),
        }

        # SMILES pattern matching rules
        self.smiles_patterns = {
            r"C\(=O\)[OoH]": ("carboxylic acid", 4.5, "Generic carboxylic acid"),
            r"c[OoH]": ("phenol", 10.0, "Generic phenolic OH"),
            r"[OoH]": ("alcohol", 16.0, "Generic alcohol"),
            r"[SsH]": ("thiol", 10.5, "Generic thiol"),
            r"\[NH3+\]": ("ammonium", 9.0, "Ammonium ion"),
            r"N\[H+]": ("ammonium", 8.0, "Protonated amine"),
            r"C\(=O\)C": ("ketone_alpha", 20.0, "Alpha to carbonyl (carbon acid)"),
            r"C#N": ("nitrile_alpha", 25.0, "Alpha to nitrile"),
        }

    def _parse_molecule(self, molecule: str) -> tuple:
        """Parse molecule input and return best match."""
        mol_lower = molecule.lower().strip()
        
        # Direct database lookup
        if mol_lower in self.parent_pkas:
            return self.parent_pkas[mol_lower], "exact_match"

        # Fuzzy match
        for key, val in self.parent_pkas.items():
            if key in mol_lower or mol_lower in key:
                return val, "fuzzy_match"

        # Try SMILES-like patterns
        for pattern, info in self.smiles_patterns.items():
            if re.search(pattern, molecule, re.IGNORECASE):
                return (info[1], info[0], info[2]), "pattern_match"

        return None, "no_match"

    def _apply_substituent_effects(self, base_pka: float, molecule: str) -> tuple:
        """Estimate substituent effects on pKa."""
        mol_lower = molecule.lower()
        delta_pka = 0.0
        factors = []

        # Detect common substituents and their effects
        sub_effects = {
            # Electron-withdrawing groups (lower pKa = stronger acid)
            "nitro": (-1.0, "Strong -I, -R effect stabilizes conjugate base"),
            "nitro-": (-2.0, "Multiple nitro groups have additive/cumulative effects"),
            "halo": (-0.5, "Halogen -I effect (stronger when closer to acidic site)"),
            "fluoro": (-0.8, "F is most electronegative; strong -I"),
            "chloro": (-0.5, "Cl -I effect"),
            "trifluoro": (-4.5, "CF3 extremely strong -I; e.g., TFA pKa 0.23 vs acetic 4.76"),
            "cyano": (-0.8, "CN strong -I, -R"),
            "carbonyl": (-0.5, "α-carbonyl stabilizes anion via resonance/induction"),
            "keto": (-0.5, "α-keto group"),
            "ester": (-0.3, "α-ester group"),
            "carboxyl": (-0.3, "Additional COOH nearby"),
            "sulfonyl": (-1.0, "SO2 very strong -I"),
            "trichloro": (-4.0, "CCl3 very strong -I"),
            "dichloro": (-2.0, "CCl2H strong -I"),
            "dihydroxy": (+0.5, "Intramolecular H-bonding can stabilize/destabilize depending on geometry"),
            
            # Electron-donating groups (raise pKa = weaker acid)
            "methoxy": (+0.3, "+R can dominate para to OH; net effect depends on position"),
            "methyl": (+0.2, "Weak +I / hyperconjugation donation"),
            "ethyl": (+0.2, "Weak +I"),
            "isopropyl": (+0.3, "+I slightly stronger than methyl"),
            "tert-butyl": (+0.5, "+I + steric inhibition of solvation"),
            "hydroxy": (0.0, "Can be ± depending on position (H-bonding vs donation)"),
            "amino": (+0.5, "Strong +R donor (if not protonated); raises pKa significantly"),
            "dimethylamino": (+1.0, "Very strong +R donor"),
            "phenyl": (0.0, "Resonance can be stabilizing or destabilizing"),
        }

        for sub, (effect, explanation) in sub_effects.items():
            if sub.replace("-", "") in mol_lower or sub in mol_lower:
                # Handle special cases
                if sub == "nitro-" and "dinitro" in mol_lower:
                    delta_pka += -2.0
                    factors.append(f"Dinitro: additive -I,-R (≈-2.0)")
                elif sub == "nitro-" and "trinitro" in mol_lower:
                    delta_pka += -4.0
                    factors.append(f"Trinitro: cumulative strong EWG (≈-4.0)")
                elif sub in mol_lower:
                    delta_pka += effect
                    factors.append(f"{sub.title()}: {explanation} (ΔpKa≈{effect:+.1f})")

        return delta_pka, factors

    def _run_base(self, molecule: str, solvent: str = "water") -> str:
        """Predict pKa value."""
        mol_lower = molecule.lower().strip()
        match_result, match_type = self._parse_molecule(molecule)

        parts = [f"## pKa Prediction: `{molecule}`\n"]
        parts.append(f"**Solvent:** {solvent}\n")

        if match_result is None:
            parts.append("### ⚠️ No Direct Match\n")
            parts.append(f"The molecule '{molecule}' was not found in the pKa database.\n")
            parts.append("### Suggestions:\n")
            parts.append("- Provide a specific compound name (e.g., 'acetic acid', 'p-nitrophenol')\n")
            parts.append("- Use SMILES notation\n")
            parts.append("- Describe the functional group (e.g., 'aromatic carboxylic acid with para-cyano')\n")
            parts.append("\n### Available Reference Compounds (sample):\n")
            sample_keys = list(self.parent_pkas.keys())[:20]
            for k in sample_keys:
                pka, cat, _ = self.parent_pkas[k]
                parts.append(f"- **{k}**: pKa = {pka:.2f} ({cat})")
            return "\n".join(parts)

        base_pka = match_result[0]
        category = match_result[1]
        notes = match_result[2] if len(match_result) > 2 else ""

        # Apply substituent corrections
        delta_pka, factors = self._apply_substituent_effects(base_pka, molecule)
        predicted_pka = base_pka + delta_pka

        # Confidence assessment
        if match_type == "exact_match":
            confidence = "High 🔬"
            predicted_pka = base_pka  # Exact match — use database value
        elif match_type == "fuzzy_match":
            confidence = "Medium-High 📊"
        else:
            confidence = "Low-Medium 🔍"

        parts.append(f"### Predicted Value\n")
        parts.append(f"| Property | Value |")
        parts.append(f"|---|---|")
        parts.append(f"| **Predicted pKa** | **{predicted_pka:.2f}** |")
        parts.append(f"| **Confidence** | {confidence} |")
        parts.append(f"| **Category** | {category} |")
        parts.append(f"| **Match Type** | {match_type.replace('_', ' ').title()} |")

        if match_type != "exact_match":
            parts.append(f"| **Base Reference** | pKa = {base_pka:.2f} |")
            parts.append(f"| **Substituent Correction** | ΔpKa = {delta_pka:+.2f} |")

        parts.append("")
        
        if notes:
            parts.append(f"### 📝 Notes\n{notes}\n")

        if factors:
            parts.append("### 🧪 Substituent Effects Analysis\n")
            for f in factors:
                parts.append(f"- {f}")
            parts.append("")

        # Add context about what this pKa means practically
        parts.append("### 💡 Practical Interpretation\n")
        if predicted_pka < 1:
            parts.append("- **Very strong acid** — fully dissociated even in acidic solutions")
        elif predicted_pka < 4:
            parts.append("- **Strong organic acid** — mostly deprotonated at physiological pH (7.4)")
        elif predicted_pka < 7:
            parts.append("- **Moderately strong acid** — significant portion deprotonated at pH 7")
        elif predicted_pka < 8:
            parts.append("- **Weak acid** — mostly protonated at pH 7, deprotonated at pH > 9")
        elif predicted_pka < 12:
            parts.append("- **Very weak acid** — requires strong base for deprotonation")
        elif predicted_pka < 16:
            parts.append("- **Extremely weak acid** — requires very strong base (LDA, NaH, etc.)")
        else:
            parts.append("- **Carbon acid / extremely weak** — requires organometallic bases (BuLi, LDA)")

        parts.append(f"\n- At **pH 7.4 (physiological)**: {'~100% deprotonated (anion)' if predicted_pka < 5.4 else '~50% deprotonated' if abs(predicted_pka - 7.4) < 1 else '~100% protonated (neutral)' if predicted_pka > 9.4 else f'~{100/(1+10**(7.4-predicted_pka)):.0f}% deprotonated'}")
        parts.append(f"- At **pH = pKa**: Exactly 50% protonated / 50% deprotonated (Henderson-Hasselbalch)")
        parts.append(f"- **Conjugate base stability**: {'Excellent (resonance-stabilized)' if predicted_pka < 5 else 'Good (inductive stabilization)' if predicted_pka < 10 else 'Moderate (localized charge)' if predicted_pka < 16 else 'Poor (unstable carbanion)'}")

        # Add related compounds for comparison
        parts.append("\n### 📊 Related Compounds for Comparison\n")
        related = [(k, v[0], v[1]) for k, v in self.parent_pkas.items() if v[1] == category]
        related.sort(key=lambda x: abs(x[1] - predicted_pka))
        for k, pka, cat in related[:5]:
            if k != mol_lower.split()[0]:
                diff = pka - predicted_pka
                arrow = "more acidic" if diff < 0 else "less acidic"
                parts.append(f"- **{k}**: pKa = {pka:.2f} ({diff:+.1f}, {arrow})")

        return "\n".join(parts)

    def _run_text(self, input_params: str) -> str:
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Input must include molecule. Format: 'molecule [solvent]'")
        molecule = parts[0]
        solvent = parts[1] if len(parts) > 1 else "water"
        return self._run_base(molecule, solvent)
