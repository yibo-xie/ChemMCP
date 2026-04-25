import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Substrate classification rules for competition analysis
# Based on organic chemistry principles (March's Advanced Organic Chemistry, Clayden)

# Steric classification of alkyl halides / substrates
_SUBSTRATE_CLASSIFICATION = {
    # Primary substrates
    "methyl": {
        "class": "primary", "steric": "very low",
        "carbocation_stability": "very poor (CH3+ is extremely unstable)",
        "sn2_feasibility": "★ Excellent — minimal steric hindrance",
        "sn1_feasibility": "✗ Essentially impossible — CH3+ too unstable",
        "e2_feasibility": "✗ No β-hydrogens for methyl",
        "e1_feasibility": "✗ Impossible",
    },
    "primary_unchanged": {
        "class": "primary", "steric": "low",
        "carbocation_stability": "poor (1° carbocation)",
        "sn2_feasibility": "★ Excellent — good SN2 substrate",
        "sn1_feasibility": "✗ Very unlikely — 1° carbocation too unstable",
        "e2_feasibility": "Possible with strong base, but SN2 usually dominates",
        "e1_feasibility": "✗ Very unlikely",
    },
    "primary_benzylic": {
        "class": "primary (benzylic)", "steric": "low",
        "carbocation_stability": "good (resonance-stabilized benzylic cation)",
        "sn2_feasibility": "● Good — but SN1 competes in protic solvents",
        "sn1_feasibility": "● Possible — especially in polar protic solvents",
        "e2_feasibility": "Possible with strong base",
        "e1_feasibility": "● Possible — benzylic cation is stable enough",
    },
    "primary_allylic": {
        "class": "primary (allylic)", "steric": "low",
        "carbocation_stability": "good (resonance-stabilized allylic cation)",
        "sn2_feasibility": "● Good — SN2 and SN1 both possible",
        "sn1_feasibility": "● Possible — allylic cation stabilized by resonance",
        "e2_feasibility": "Possible",
        "e1_feasibility": "● Possible",
    },
    "neopentyl": {
        "class": "primary (neopentyl)", "steric": "high (β-branching)",
        "carbocation_stability": "poor (would rearrange to tert-butyl)",
        "sn2_feasibility": "◐ Very slow — severe steric hindrance (β-carbons block backside attack)",
        "sn1_feasibility": "◐ May proceed via rearrangement to tert-cation then capture",
        "e2_feasibility": "◐ Possible if base can access β-H",
        "e1_feasibility": "Via rearrangement",
    },
    "secondary": {
        "class": "secondary", "steric": "moderate",
        "carbocation_stability": "moderate (2° carbocation)",
        "sn2_feasibility": "● Good — moderate rate",
        "sn1_feasibility": "● Competes with SN2 in protic solvents",
        "e2_feasibility": "● Significant pathway with strong bases",
        "e1_feasibility": "● Possible at elevated temperature in protic solvents",
    },
    "secondary_benzylic": {
        "class": "secondary (benzylic)", "steric": "moderate",
        "carbocation_stability": "very good (resonance + 2°)",
        "sn2_feasibility": "● Works but SN1 often dominates",
        "sn1_feasibility": "★ Favored in protic solvents",
        "e2_feasibility": "● Possible",
        "e1_feasibility": "★ Likely in protic solvent + heat",
    },
    "tertiary": {
        "class": "tertiary", "steric": "high",
        "carbocation_stability": "good (3° carbocation)",
        "sn2_feasibility": "✗ Essentially impossible — backside blocked",
        "sn1_feasibility": "★ Dominant pathway in protic solvents",
        "e2_feasibility": "★ Dominant with strong base",
        "e1_feasibility": "★ Competes with E2 in protic solvents",
    },
    "tertiary_benzylic": {
        "class": "tertiary (benzylic/phenylethyl)", "steric": "high",
        "carbocation_stability": "excellent (3° + resonance)",
        "sn2_feasibility": "✗ Impossible",
        "sn1_feasibility": "★★ Extremely fast — excellent substrate for SN1/E1",
        "e2_feasibility": "● Occurs but SN1/E1 dominate",
        "e1_feasibility": "★★ Very fast",
    },
}

# Reagent strength classification
_REAGENT_CLASSIFICATION = {
    # Strong nucleophiles / weak bases (favor SN2)
    "I-", "Br-", "RS-", "CN-", "N3-", "CH3COO-",
    # Strong nucleophiles / strong bases (SN2 or E2)
    "HO-", "RO-", "CH3O-",
    # Strong bases / poor nucleophiles (favor E2)
    "t-BuO-", "LDA", "NH2-",
    # Weak bases / poor nucleophiles (favor SN1/E1)
    "H2O", "ROH", "CH3COOH", "carboxylic acids",
}

_BASE_STRENGTH = {
    "very_strong_base": ["t-BuO-", "LDA", "NH2-", "H-", "Ph3C-K"],
    "strong_base": ["HO-", "CH3O-", "EtO-", "NaH", "KH"],
    "moderate_base": ["CN-", "acetylide", "RMgX (acts as base)"],
    "weak_base": ["I-", "Br-", "Cl-", "CH3COO-", "N3-", "RS-", "H2O", "ROH"],
}

_NUCLEOPHILE_STRENGTH = {
    "excellent_nu": ["I-", "Br-", "RS-", "SePh-", "CN-", "N3-"],
    "good_nu": ["HO-", "RO-", "CH3O-", "CH3S-", "PhS-", "NH3", "RNH2"],
    "poor_nu": ["t-BuO-", "H2O", "ROH", "carboxylates"],
}

_SOLVENT_EFFECTS_COMPETITION = {
    "protic_polar": {
        "effect_on_sn2": "Slows down nucleophile (H-bonding solvation)",
        "effect_on_sn1": "Promotes (stabilizes ions/carbocations)",
        "effect_on_e2": "May slow strong base slightly",
        "effect_on_e1": "Promotes (stabilizes carbocation intermediate)",
        "overall": "Favors unimolecular pathways (SN1/E1) over bimolecular (SN2/E2)",
    },
    "aprotic_polar": {
        "effect_on_sn2": "★★ Dramatically accelerates ('naked' anions)",
        "effect_on_sn1": "Less effect (carbocation already stabilized)",
        "effect_on_e2": "Accelerates (base more reactive)",
        "effect_on_e1": "Minor effect",
        "overall": "Favors bimolecular pathways (SN2/E2), especially SN2",
    },
    "nonpolar": {
        "effect_on_sn2": "Poor — doesn't dissolve ionic species",
        "effect_on_sn1": "Doesn't promote ionization well",
        "effect_on_e2": "Limited dissolution of ionic reagents",
        "effect_on_e1": "Limited",
        "overall": "Generally poor for all ionic mechanisms. Radical or pericyclic pathways may dominate.",
    },
}


@ChemMCPManager.register_tool
class CompetitionAnalyzer(BaseTool):
    """
    分析竞争反应（SN1 vs SN2 vs E1 vs E2）的工具。
    基于底物结构、试剂强度、溶剂和温度，综合分析哪种反应机理占主导地位。
    """
    __version__      = "0.1.0"
    name             = "CompetitionAnalyzer"
    func_name        = "analyze_competition"
    description      = "Analyze competing reaction pathways (SN1 vs SN2 vs E1 vs E2) based on substrate structure, reagent, solvent, and temperature."
    implementation_description = "Uses rule-based expert system incorporating substrate classification (9 substrate types), reagent strength/nucleophilicity categories, solvent effects, and temperature dependence to predict dominant mechanism(s)."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Competition Analysis", "Reaction Mechanisms", "Substitution", "Elimination"]
    required_envs    = []

    code_input_sig   = [
        ("substrate_type", "str", "N/A", "Substrate type: 'methyl', 'primary', 'primary_benzylic', 'primary_allylic', 'neopentyl', 'secondary', 'secondary_benzylic', 'tertiary', 'tertiary_benzylic'."),
        ("reagent", "str", "NaOH", "Reagent/base/nucleophile (e.g., 'NaOH', 'KOtBu', 'KI', 'H2O', 'CH3COOH')."),
        ("solvent", "str", "ethanol", "Solvent: 'water', 'methanol', 'ethanol', 'acetone', 'DMSO', 'DMF', 'THF'."),
        ("temperature_c", "float", "25.0", "Temperature in °C (default: 25°C)."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'substrate_type [reagent] [solvent] [temp_C]'. Example: 'secondary NaOH ethanol 60'"),
    ]

    output_sig       = [
        ("result", "str", "Comprehensive analysis predicting the dominant mechanism with reasoning and product distribution."),
    ]

    examples         = [
        {
            "code_input": {
                "substrate_type": "secondary",
                "reagent": "NaOH",
                "solvent": "ethanol",
                "temperature_c": 60.0,
            },
            "text_input": {"input_text": "secondary NaOH ethanol 60"},
            "output": {
                "result": """## Competition Analysis: Secondary Alkyl Halide + NaOH

### Substrate: **Secondary Alkyl Halide** (R₂CH-X)

| Property | Assessment |
|----------|-----------|
| Steric Profile | Moderate — backside partially accessible |
| Carbocation Stability | Moderate (2° cation) |
| SN2 Feasibility | ● Good |
| SN1 Feasibility | ● Competes in protic solvent |
| E2 Feasibility | ● Significant with strong base |
| E1 Feasibility | ● Possible at elevated T |

### Reagent: **NaOH** → provides HO⁻ (strong base + good nucleophile)

### Solvent: **Ethanol** (protic, polar)

### Temperature: **60°C**

---

## 🏆 Predicted Outcome

### Primary Pathway: **Mixture of SN2 + E2**

| Mechanism | Likelihood | Reasoning |
|-----------|-----------|-----------|
| **SN2** | ████████░░ ~40% | HO⁻ is a good nucleophile; secondary substrate accessible; but protic solvent slows it down |
| **E2** | ████████░░ ~40% | HO⁻ is a strong base; 60°C favors elimination; secondary substrates give significant E2 |
| **SN1** | ██░░░░░░░░ ~15% | Protic solvent promotes some ionization; 60°C helps; but HO⁻ is not a great leaving group for SN1 |
| **E1** | █░░░░░░░░░ ~5% | Minor pathway; requires full ionization first |

### Key Factors:
1. **Solvent effect**: Ethanol (protic) slows SN2 somewhat but promotes ionization pathways
2. **Base strength**: HO⁻ is both strong base AND good nucleophile → competition between SN2 and E2
3. **Temperature**: 60°C shifts equilibrium toward elimination (more entropy gain in E2 than SN2)
4. **Substrate**: Secondary is the crossover point — no single mechanism dominates

### Expected Products:
- **SN2 product**: R₂CH-OH (substitution, possibly with inversion)
- **E2 products**: Mixture of alkenes (Zaitsev + Hofmann depending on base sterics)
  - With HO⁻ (non-bulky): Zaitsev (more substituted alkene) favored
  - Ratio depends on specific substrate structure

### How to Favor Each Pathway:
- **→ Pure SN2**: Use NaI/acetone (aprotic polar), RT, good nucleophile (I⁻, CN⁻)
- **→ Pure E2**: Use KOtBu/t-BuOH (bulky strong base), heat (80°C+)
- **→ SN1**: Use H₂O or ROH as solvent/nucleophile, heat, no strong base"""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, substrate_type: str, reagent: str = "NaOH", solvent: str = "ethanol", temperature_c: float = 25.0) -> str:
        """Core logic: analyze competition."""
        sub = substrate_type.strip().lower()
        reag = reagent.strip()
        solv = solvent.strip().lower()

        lines = []
        display_sub = substrate_type.replace("_", " ").title()
        lines.append(f"## Competition Analysis: {display_sub} + {reag}\n")

        # Get substrate data
        sub_data = self._match_substrate(sub)
        if sub_data is None:
            lines.append(f"⚠️ Substrate type '{substrate_type}' not recognized.")
            lines.append("**Available:** " + ", ".join(sorted(_SUBSTRATE_CLASSIFICATION.keys())))
            return "\n".join(lines)

        # Substrate assessment table
        lines.append(f"### Substrate: **{display_sub}**\n")
        lines.append("| Property | Assessment |")
        lines.append("|----------|-----------|")
        for key in ["class", "steric", "carbocation_stability", "sn2_feasibility", "sn1_feasibility", "e2_feasibility", "e1_feasibility"]:
            val = sub_data.get(key, "N/A")
            label = key.replace("_", " ").title()
            lines.append(f"| {label} | {val} |")
        lines.append("")

        # Reagent analysis
        lines.append(f"### Reagent: **{reag}**\n")
        reagent_analysis = self._analyze_reagent(reag)
        lines.append(reagent_analysis)
        lines.append("")

        # Solvent analysis
        lines.append(f"### Solvent: **{solvent.title()}**\n")
        solvent_analysis = self._classify_solvent(solv)
        lines.append(solvent_analysis)
        lines.append("")

        # Temperature factor
        lines.append(f"### Temperature: **{temperature_c}°C**\n")
        temp_factor = self._temp_factor(temperature_c)
        lines.append(temp_factor)
        lines.append("")

        # Predict outcome
        lines.append("---\n")
        lines.append("## 🏆 Predicted Outcome\n")
        predictions = self._predict(sub_data, reag, solv, temperature_c)
        lines.append(predictions)
        lines.append("")

        # How to favor each pathway
        lines.append("### How to Favor Each Pathway:\n")
        favor_tips = self._favor_tips(sub, reag, solv)
        lines.append(favor_tips)

        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        sub = parts[0] if parts else "secondary"
        reag = parts[1] if len(parts) > 1 else "NaOH"
        solv = parts[2] if len(parts) > 2 else "ethanol"
        temp = float(parts[3]) if len(parts) > 3 else 25.0
        return self._run_base(sub, reag, solv, temp)

    def _match_substrate(self, sub):
        exact = _SUBSTRATE_CLASSIFICATION.get(sub)
        if exact:
            return exact
        for k, v in _SUBSTRATE_CLASSIFICATION.items():
            if sub in k or k.startswith(sub):
                return v
        aliases = {
            "1°": "primary_unchanged", "pri": "primary_unchanged", "1": "primary_unchanged",
            "2°": "secondary", "sec": "secondary", "2": "secondary",
            "3°": "tertiary", "tert": "tertiary", "3": "tertiary",
            "benzyl": "primary_benzylic", "allyl": "primary_allylic",
        }
        return _SUBSTRATE_CLASSIFICATION.get(aliases.get(sub))

    def _analyze_reagent(self, reag):
        reag_lower = reag.lower().replace(" ", "")
        # Determine what species the reagent provides
        species_map = {
            "naoh": "HO-", "koh": "ho-", "naoch3": "CH3O-", "ko t-bu": "t-BuO-", "kotbu": "t-BuO-",
            "t-buok": "t-BuO-", "lda": "LDA", "nah": "H-", "ki": "I-", "nai": "I-",
            "kbr": "br-", "nabr": "br-", "kcn": "cn-", "nan3": "N3-",
            "h2o": "H2O", "ch3oh": "CH3OH", "etoh": "C2H5OH", "ch3cooh": "CH3COOH",
            "ch3coona": "CH3COO-", "nach3coo": "CH3COO-",
        }
        active = species_map.get(reag_lower, reag)

        base_cat = "Unknown"
        nu_cat = "Unknown"
        for cat, members in _BASE_STRENGTH.items():
            if active in members or any(active.lower() in m.lower() for m in members):
                base_cat = cat
                break
        for cat, members in _NUCLEOPHILE_STRENGTH.items():
            if active in members or any(active.lower() in m.lower() for m in members):
                nu_cat = cat
                break

        lines = []
        lines.append(f"| Active Species | **{active}** |")
        lines.append(f"| Base Strength | {base_cat.replace('_', ' ').title()} |")
        lines.append(f"| Nucleophilicity | {nu_cat.replace('_', ' ').title()} |")
        return "\n".join(lines)

    def _classify_solvent(self, solv):
        mapping = {
            "water": "protic_polar", "methanol": "protic_polar", "ethanol": "protic_polar",
            "isopropanol": "protic_polar", "t-butanol": "protic_polar",
            "acetic_acid": "protic_polar", "formic_acid": "protic_polar",
            "acetone": "aprotic_polar", "dmso": "aprotic_polar", "dmf": "aprotic_polar",
            "acetonitrile": "aprotic_polar", "thf": "aprotic_polar", "ether": "aprotic_polar",
            "dcm": "aprotic_nonpolar", "chloroform": "aprotic_nonpolar",
            "hexane": "nonpolar", "toluene": "nonpolar", "benzene": "nonpolar",
        }
        key = mapping.get(solv, "unknown")
        data = _SOLVENT_EFFECTS_COMPETITION.get(key, {})
        if not data:
            return f"⚠️ Solvent '{solv}' not classified.\n"

        lines = []
        lines.append("| Effect on Mechanism | Impact |")
        lines.append("|---------------------|--------|")
        for mech, effect in data.items():
            if mech != "overall":
                mech_label = mech.replace("effect_on_", "").upper().replace("_", " ")
                lines.append(f"| {mech_label} | {effect} |")
        lines.append(f"\n**Overall:** {data.get('overall', 'N/A')}")
        return "\n".join(lines)

    def _temp_factor(self, temp_c):
        if temp_c < 0:
            return "❄️ **Low temperature:** Favors kinetic control (SN2 over E2, substitution over elimination). Reduces all rates but SN2 least affected.\n"
        elif temp_c <= 40:
            return "🌡️ **Room temperature to moderate:** Standard conditions. All mechanisms possible depending on other factors.\n"
        elif temp_c <= 80:
            return "🔥 **Elevated temperature:** Favors elimination (E1/E2) over substitution due to positive ΔS‡ for elimination. Also promotes SN1/E1 (ionization needs thermal energy).\n"
        else:
            return "🔥🔥 **High temperature:** Strongly favors elimination (E2 > E1 > SN1 > SN2). At very high T, E2 dominates for most substrates.\n"

    def _predict(self, sub_data, reag, solv, temp_c):
        scores = {"SN2": 0, "SN1": 0, "E2": 0, "E1": 0}
        reasons = {k: [] for k in scores}

        # Substrate factors
        sn2_score = sub_data["sn2_feasibility"].count("★") * 20 + sub_data["sn2_feasibility"].count("●") * 10
        sn1_score = sub_data["sn1_feasibility"].count("★") * 20 + sub_data["sn1_feasibility"].count("●") * 10
        e2_score = sub_data["e2_feasibility"].count("★") * 20 + sub_data["e2_feasibility"].count("●") * 10
        e1_score = sub_data["e1_feasibility"].count("★") * 20 + sub_data["e1_feasibility"].count("●") * 10

        scores["SN2"] += sn2_score
        scores["SN1"] += sn1_score
        scores["E2"] += e2_score
        scores["E1"] += e1_score

        # Reagent effects
        reag_lower = reag.lower().replace(" ", "")
        if any(x in reag_lower for x in ["t-buo", "t-buok", "kotbu", "lda"]):
            scores["E2"] += 30
            scores["SN2"] -= 20
            reasons["E2"].append("Bulky strong base → E2 favored, SN2 suppressed")
        elif any(x in reag_lower for x in ["naoh", "koh", "naoch3", "ko ch3"]):
            scores["SN2"] += 15
            scores["E2"] += 15
            reasons["SN2"].append("Strong base/good nucleophile → SN2 competitive")
            reasons["E2"].append("Strong base → E2 competitive")
        elif any(x in reag_lower for x in ["i-", "ki", "nai", "br-", "kbr", "nabr", "cn-", "kcn"]):
            scores["SN2"] += 25
            scores["E2"] -= 10
            reasons["SN2"].append("Good nucleophile/weak base → SN2 strongly favored")
        elif any(x in reag_lower for x in ["h2o", "roh", "ch3oh", "etoh", "cooh"]):
            scores["SN1"] += 20
            scores["E1"] += 15
            scores["SN2"] -= 15
            reasons["SN1"].append("Weak nucleophile → unimolecular pathways favored")
            reasons["E1"].append("Weak base/poor nucleophile → E1 possible")

        # Solvent effects
        if solv in ("water", "methanol", "ethanol", "isopropanol", "t-butanol"):
            scores["SN1"] += 15
            scores["E1"] += 10
            scores["SN2"] -= 10
            reasons["SN1"].append("Protic solvent stabilizes carbocation → SN1 promoted")
        elif solv in ("dmso", "dmf", "acetone", "acetonitrile"):
            scores["SN2"] += 20
            scores["E2"] += 10
            reasons["SN2"].append("Polar aprotic → naked anions → SN2 accelerated")
        elif solv in ("dcm", "chloroform", "hexane", "toluene"):
            scores["SN1"] -= 10
            scores["SN2"] -= 10
            reasons["SN1"].append("Nonpolar/poorly solvating → ionic mechanisms slowed")

        # Temperature
        if temp_c > 50:
            scores["E2"] += 15
            scores["E1"] += 10
            scores["SN2"] -= 5
            reasons["E2"].append("Higher T favors elimination (ΔS‡ > 0)")
        elif temp_c < 10:
            scores["SN2"] += 10
            scores["E2"] -= 10
            reasons["SN2"].append("Low T favors substitution (lower Ea for SN2)")

        # Normalize to percentages
        total = sum(max(0, v) for v in scores.values())
        if total == 0:
            total = 1

        lines = []
        lines.append("| Mechanism | Likelihood | Key Reasons |")
        lines.append("|-----------|-----------|-------------|")

        sorted_mechs = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        for mech, score in sorted_mechs:
            pct = max(0, score) / total * 100
            bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
            reason_txt = "; ".join(reasons[mech][:2]) if reasons[mech] else "See substrate analysis"
            lines.append(f"| **{mech}** | {bar} ~{pct:.0f}% | {reason_txt} |")

        lines.append("")
        lines.append("### Expected Products:")
        lines.append(self._product_prediction(sorted_mechs, sub_data))
        return "\n".join(lines)

    def _product_prediction(self, sorted_mechs, sub_data):
        top = sorted_mechs[0][0]
        cls = sub_data.get("class", "")

        products = []
        if top == "SN2":
            products.append("- **SN2 product**: Substituted product with **inverted stereochemistry** (Walden inversion)")
        elif top == "SN1":
            products.append("- **SN1 product**: Substituted product with **racemization** (+ partial inversion via ion pair)")
            products.append("  - If chiral center: racemic mixture ± slight inversion")
        elif top == "E2":
            zaitsev = "- **E2 product**: **Zaitsev alkene** (more substituted, thermodynamically favored)" if "bulky" not in str(sorted_mechs) else "- **E2 product**: **Hofmann alkene** (less substituted, from bulky base)"
            products.append(zaitsev)
        elif top == "E1":
            products.append("- **E1 product**: Mixture of alkenes (**Zaitsev major**) via carbocation intermediate")
            products.append("  - Possible rearrangement products if carbocation can rearrange")

        # Check for mixture
        if len([m for m, s in sorted_mechs if s > 10]) >= 2:
            products.append("\n⚠️ **Multiple mechanisms compete** — expect product mixture. Use the tips below to direct selectivity.")

        return "\n".join(products)

    def _favor_tips(self, sub, reag, solv):
        tips = []
        sub_data = self._match_substance_for_tips(sub)
        if sub_data:
            if "primary" in sub_data.get("class", "") and "benzylic" not in sub and "allylic" not in sub:
                tips.append("- **→ Pure SN2**: Use NaI/acetone or KCN/DMSO, RT, polar aprotic solvent")
            elif "tertiary" in sub_data.get("class", ""):
                tips.append("- **→ Pure E2**: Use KOtBu/t-BuOH, 80°C+ (bulky strong base forces elimination)")
                tips.append("- **→ SN1/E1**: Use H2O or EtOH as solvent/nucleophile, heat, weak base")
            elif "secondary" in sub_data.get("class", ""):
                tips.append("- **→ SN2**: Use NaN3/DMSO or NaI/acetone, RT (polar aprotic, good Nu)")
                tips.append("- **→ E2**: Use KOtBu/t-BuOH, heat (bulky base)")
                tips.append("- **→ SN1**: Use AgNO3/EtOH or H2O/EtOH, heat (promote ionization)")
        else:
            tips = [
                "- **→ Favor SN2**: Polar aprotic solvent (DMSO, acetone), good nucleophile (I-, CN-), RT",
                "- **→ Favor SN2**: Concentrated conditions, less steric hindrance",
                "- **→ Favor E2**: Bulky strong base (KOtBu, LDA), high temperature, aprotic solvent",
                "- **→ Favor E1/SN1**: Protic solvent (EtOH, H2O), heat, weak base/nucleophile, dilute",
            ]
        return "\n".join(tips)

    def _match_substance_for_tips(self, sub):
        return _SUBSTRATE_CLASSIFICATION.get(sub)
