import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Database of reactions where kinetic vs thermodynamic control is relevant
_KINETIC_THERMO_DB = {
    "conjugated_addition_1_2_vs_1_4": {
        "name": "1,2- vs 1,4-Conjugate Addition (Michael Addition)",
        "description": "Addition of nucleophiles to α,β-unsaturated carbonyl compounds.",
        "kinetic_product": {
            "name": "1,2-Addition (direct addition)",
            "description": "Nucleophile adds directly to the carbonyl carbon.",
            "conditions": "Low temperature (-78°C to 0°C), hard nucleophiles (organolithiums, Grignards), kinetic control.",
            "stereochemistry": "Depends on nucleophile and substrate chirality.",
            "example": "CH3Li + CH2=CH-COCH3 (at -78°C) → CH3C(OLi)(CH3)CH=CH2 after workup: 1,2-adduct",
            "reversibility": "Often irreversible under kinetic conditions.",
        },
        "thermodynamic_product": {
            "name": "1,4-Addition (conjugate / Michael addition)",
            "description": "Nucleophile adds to the β-carbon via conjugate addition mechanism.",
            "conditions": "Higher temperature (RT to reflux), soft nucleophiles (enolates, thiolates, Cu-catalyzed), thermodynamic control or equilibration possible.",
            "stereochemistry": "Trans relationship between new substituents typical for cyclic systems.",
            "example": "CuI-catalyzed Me2CuLi + CH2=CH-COCH3 (RT) → CH3CH2COCH3 (via 1,4-addition then tautomerization)",
            "reversibility": "Can be reversible under basic conditions; conjugate product more stable (conjugated enolate).",
        },
        "how_to_favor_kinetic": [
            "Use LOW temperature (-78°C to 0°C)",
            "Use HARD nucleophile (RLi > RMgX > organozinc/enolate)",
            "Use non-polar solvent (THF, ether) at low T",
            "Avoid Lewis acids or catalysts that promote conjugation",
            "Quench rapidly before equilibration",
        ],
        "how_to_favor_thermodynamic": [
            "Use HIGHER temperature (RT to reflux)",
            "Use SOFT nucleophile (enolates, Gilman reagents R2CuLi, thiolates)",
            "Add Cu(I) salt as catalyst (promotes 1,4-addition)",
            "Use polar coordinating solvents",
            "Allow time for equilibration (or add base to catalyze reversibility)",
        ],
        "key_factor": "Hard nucleophiles prefer the harder carbonyl carbon (1,2); soft nucleophiles prefer softer β-carbon (1,4). Temperature controls reversibility and equilibration rate.",
    },
    "enolate_formation_kinetic_vs_thermodynamic": {
        "name": "Kinetic vs Thermodynamic Enolate Formation",
        "description": "Deprotonation of unsymmetrical ketones can give two different regioisomeric enolates.",
        "kinetic_product": {
            "name": "Kinetic enolate (less substituted)",
            "description": "Deprotonation at the less hindered α-carbon (faster deprotonation, lower Ea).",
            "conditions": "Strong sterically hindered base (LDA), LOW temperature (-78°C), aprotic solvent (THF), short reaction time.",
            "stereochemistry": "(E)-enolate typically favored kinetically for acyclic systems.",
            "example": "2-methylcyclohexanone + LDA/THF/-78°C → kinetic enolate (less substituted C=C after alkylation)",
            "reversibility": "Irreversible under these conditions (no equilibrium).",
        },
        "thermodynamic_product": {
            "name": "Thermodynamic enolate (more substituted)",
            "description": "Deprotonation at the more substituted α-carbon (more stable enolate, lower G).",
            "conditions": "Weaker base, HIGHER temperature (0°C to RT), protic solvent or longer reaction time allowing equilibration.",
            "stereochemistry": "(Z)-enolate often favored at thermodynamic equilibrium for acyclic systems.",
            "example": "2-methylcyclohexanone + NaH/THF/0°C→RT (or NaOEt/EtOH/reflux) → thermodynamic enolate (more substituted)",
            "reversibility": "Reversible — enolate can reprotonate and redeprotonate at either position.",
        },
        "how_to_favor_kinetic": [
            "LDA (or KHMDS, LiTMP) in THF at **-78°C** (dry ice/acetone)",
            "Short reaction time (< 30 min before adding electrophile)",
            "Excess LDA (1.1-2.0 eq relative to ketone)",
            "Aprotic anhydrous conditions (strict!)",
            "Add electrophile immediately after enolate formation",
        ],
        "how_to_favor_thermodynamic": [
            "Use **weaker bulkier base** (t-BuOK in t-BuOH) or use **NaH/ROH** system",
            "Higher temperature (**0°C to RT or reflux**)",
            "Protic solvent (t-BuOH, EtOH-HOP) allows proton exchange/equilibration",
            "Longer reaction time (hours) to reach equilibrium",
            "Use 18-crown-6 with K+ to enhance equilibration rate",
        ],
        "key_factor": "Kinetic enolate forms faster (less hindered H). Thermodynamic enolate is more stable (more substituted C=C, lower ΔG). The barrier to interconversion determines which one you get.",
    },
    "alkylation_of_indole_C2_vs_C3": {
        "name": "Indole Alkylation: C2 vs C3",
        "description": "Electrophilic aromatic substitution / alkylation of indole can occur at C2 or C3 position.",
        "kinetic_product": {
            "name": "C3-substituted indole",
            "description": "Attack at C3 (position 3 of indole ring) is kinetically favored due to higher electron density and less steric hindrance.",
            "conditions": "Mild electrophiles, lower temperature, standard EAS conditions.",
            "example": "Indole + alkyl halide/Lewis acid (RT) → 3-alkylindole (major kinetic product)",
        },
        "thermodynamic_product": {
            "name": "C2-substituted indole",
            "description": "Under equilibrating conditions, C2-substituted product can be favored due to greater stability (aromaticity preservation in benzenoid form).",
            "conditions": "Reversible conditions, higher T, acid catalysis allowing rearrangement.",
            "example": "Indole + certain electrophiles under reversible conditions → 2-alkylindole (after rearrangement)",
        },
        "how_to_favor_kinetic": ["Milder conditions", "Non-reversible electrophiles", "Lower temperature"],
        "how_to_favor_thermodynamic": ["Acid catalysis with reversible intermediate", "Higher temperature", "Longer reaction times"],
        "key_factor": "C3 has higher electron density (kinetically favored), but C2 substitution may be thermodynamically more stable in some cases.",
    },
    "diels_allder_endo_vs_exo": {
        "name": "Diels-Alder: endo vs exo Selectivity",
        "description": "Diels-Alder cycloaddition can give endo or exo stereoisomers.",
        "kinetic_product": {
            "name": "endo product",
            "description": "Favored by secondary orbital interactions in the transition state (lower activation energy despite being sterically more crowded).",
            "conditions": "Normal electron-demand DA, lower temperature, kinetic control.",
            "example": "Cyclopentadiene + maleic anhydride (RT, short time) → endo adduct (kinetic)",
        },
        "thermodynamic_product": {
            "name": "exo product",
            "description": "Less sterically crowded, more thermodynamically stable.",
            "conditions": "Higher temperature, longer reaction time allowing retro-Diels-Alder and re-addition (equilibration).",
            "example": "Cyclopentadiene + maleic anhydride (reflux, long time) → exo adduct (thermodynamic, after equilibration)",
        },
        "how_to_favor_kinetic": ["Lower temperature (0°C to RT)", "Shorter reaction time", "Normal electron demand"],
        "how_to_favor_thermodynamic": ["Higher temperature (reflux in xylene/toluene)", "Longer reaction time (allows retro-DA equilibrium)", "High pressure (if volume change favors exo)"],
        "key_factor": "endo rule: kinetic preference due to secondary orbital interactions. At high T, retro-DA allows equilibration to more stable exo.",
    },
    "sulfonation_of_naphthalene_alpha_vs_beta": {
        "name": "Naphthalene Sulfonation: α vs β",
        "description": "Sulfonation of naphthalene gives α- or β-naphthalenesulfonic acid depending on temperature.",
        "kinetic_product": {
            "name": "α-naphthalenesulfonic acid",
            "description": "Formed faster at low temperature due to more favorable intermediate (arenium ion stabilized by more resonance forms for attack at C1).",
            "conditions": "Low temperature (< 80°C).",
        },
        "thermodynamic_product": {
            "name": "β-naphthalenesulfonic acid",
            "description": "More stable product (less steric hindrance between SO3H group and peri-H at C8). Favored at high T where reaction is reversible.",
            "conditions": "High temperature (> 160°C). Reversible sulfonation allows equilibration.",
        },
        "how_to_favor_kinetic": ["Low temperature sulfonation (< 80°C)"],
        "how_to_favor_thermodynamic": ["High temperature sulfonation (> 160°C)", "Prolonged heating"],
        "key_factor": "Classic textbook example of kinetic vs thermodynamic control. α is kinetic (faster formation), β is thermodynamic (more stable, reversible at high T).",
    },
    "bromination_of_phenol_ortho_vs_para_bromination": {
        "name": "Phenol Bromination (ortho vs para ratio)",
        "description": "Electrophilic bromination of phenol gives ortho and para products whose ratio depends on conditions (solvent, polarity, T).",
        "kinetic_product": "Ortho-bromophenol (often kinetically favored due to intramolecular H-bonding stabilization of ortho-TS in nonpolar solvent).",
        "thermodynamic_product": "Para-bromophenol (generally more stable, less steric crowding).",
        "how_to_favor_kinetic": ["Nonpolar solvent (CS2, CCl4)", "Low temperature", "No added base"],
        "how_to_favor_thermodynamic": ["Polar solvent (water)", "Higher temperature", "Presence of base (deprotonates phenol → phenoxide which strongly favors para)"],
        "key_factor": "Solvent polarity dramatically changes the ortho/para ratio through hydrogen bonding effects and ionization.",
    },
    "aldol_addition_vs_condensation": {
        "name": "Aldol: Addition vs Condensation (Dehydration)",
        "description": "Aldol reaction can give β-hydroxy carbonyl (addition) or α,β-unsaturated carbonyl (condensation/dehydration).",
        "kinetic_product": {
            "name": "β-hydroxy carbonyl (aldol adduct)",
            "description": "Direct addition product before dehydration.",
            "conditions": "Low temperature, mild conditions, irreversible aldol (e.g., Li-enolate additions).",
        },
        "thermodynamic_product": {
            "name": "α,β-unsaturated carbonyl (condensation product)",
            "description": "Dehydrated conjugated enone — more stable due to extended conjugation.",
            "conditions": "Higher temperature, acidic or basic workup, reversible aldol conditions.",
        },
        "how_to_favor_kinetic": ["Low temperature (-78°C to 0°C)", "Irreversible enolate (Li, Zr, B enolates)", "Short reaction time", "Mild aqueous workup (avoid acid/base that promotes dehydration)"],
        "how_to_favor_thermodynamic": ["Heat (reflux)", "Acid or base catalysis during workup", "Longer reaction time", "Removal of water (drives dehydration forward)"],
        "key_factor": "The adduct is kinetic; dehydration gives the more stable conjugated system (thermodynamic). Control depends on reversibility and elimination conditions.",
    },
}


@ChemMCPManager.register_tool
class KineticVsThermodynamic(BaseTool):
    """
    判断动力学控制与热力学控制产物的工具。
    内置 7 种经典反应类型的动力学/热力学控制数据，包含产物结构、条件、选择性和控制方法。
    """
    __version__      = "0.1.0"
    name             = "KineticVsThermodynamic"
    func_name        = "analyze_kinetic_vs_thermodynamic"
    description      = "Determine whether kinetic or thermodynamic control dominates for a given reaction, predict the major product, and provide guidance on how to favor each outcome."
    implementation_description = "Uses embedded database of 7 classic kinetic vs thermodynamic control scenarios (conjugate addition, enolate formation, indole alkylation, Diels-Alder stereochemistry, naphthalene sulfonation, phenol bromination, aldol) with full condition analysis."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Kinetic Control", "Thermodynamic Control", "Reaction Mechanisms", "Selectivity"]
    required_envs    = []

    code_input_sig   = [
        ("reaction_type", "str", "N/A", "Type of reaction: 'conjugate_addition', 'enolate', 'indole', 'diels_alder', 'naphthalene_sulfonation', 'phenol_bromination', 'aldol'."),
        ("temperature_c", "float", "25.0", "Reaction temperature in °C."),
        ("conditions", "str", "standard", "Special conditions: 'low_T', 'high_T', 'hard_nucleophile', 'soft_nucleophile', 'reversible', 'irreversible', 'standard'."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'reaction_type [temperature] [conditions]'. Example: 'enolate -78 irreversible'"),
    ]

    output_sig       = [
        ("result", "str", "Analysis of kinetic vs thermodynamic control with predicted major product and how to favor each pathway."),
    ]

    examples         = [
        {
            "code_input": {"reaction_type": "enolate", "temperature_c": -78.0, "conditions": "irreversible"},
            "text_input": {"input_text": "enolate -78 irreversible"},
            "output": {
                "result": """## Kinetic vs Thermodynamic Analysis: Enolate Formation

### Reaction: Kinetic vs Thermodynamic Enolate Formation

---

### 🎯 Prediction: **KINETIC CONTROL** ⚡

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| Temperature | -78°C | ★ Very low → kinetic regime |
| Conditions | Irreversible | No equilibration possible |
| Regime | **KINETIC CONTROL** | Product determined by relative rates |

---

### Expected Major Product: **Kinetic Enolate** (Less Substituted)

**Structure:** Deprotonation occurs at the **less hindered α-carbon**

| Property | Detail |
|----------|--------|
| **Product** | Less substituted enolate → less substituted alkylated product |
| **Why?** | Lower activation energy for removing the less sterically hindered proton |
| **Stereochemistry** | Typically (E)-enolate for acyclic systems |
| **Reversibility** | ❌ Irreversible at -78°C — no equilibration |

---

### How It Works:
```
     O                            O(-)
    //  + LDA ( -78°C )   →       ||
    R-CH2-C-CH3         (kinetic)  R-CH=C-CH3   ← LESS substituted (kinetic)
          |                          |
         H  (removed first)        R-group

vs.

     O                            O(-)
    //  + Base (RT/reflux)  →      ||
    R-CH2-C-CH3        (thermo.)   R-CH=C-CH3   ← MORE substituted (thermo)
          |                       |
         H  (removed later)      R-group (more stable C=C)
```

---

### How to Switch to Thermodynamic Product:

| Change | Effect |
|--------|--------|
| Raise T to **0°C → RT** | Allows enolate equilibration |
| Switch base to **t-BuOK/t-BuOH** | Weaker, bulkier base + protic solvent |
| Use **NaH/EtOH then reflux** | Protic solvent enables proton exchange |
| Add **18-crown-6/K+** | Enhances equilibration rate |
| Use **longer reaction time** | Hours instead of minutes |

⚠️ **Key Insight:** The kinetic product forms faster but is NOT the most stable enolate. To get the thermodynamic product, you MUST allow equilibration (higher T, reversible conditions, weaker base).

---

*Reference: House, Modern Synthetic Methods; Smith, March's Advanced Organic Chemistry*"""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, temperature_c: float = 25.0, conditions: str = "standard") -> str:
        """Core logic: analyze kinetic vs thermodynamic control."""
        rt = reaction_type.strip().lower()
        cond = conditions.strip().lower() if conditions else "standard"

        lines = []
        display_name = rt.replace("_", " ").title()
        lines.append(f"## Kinetic vs Thermodynamic Analysis: {display_name}\n")

        # Match reaction type
        data = self._match_reaction(rt)
        if data is None:
            lines.append(f"⚠️ Reaction type '{reaction_type}' not found.")
            lines.append("\n**Available reaction types:**\n")
            for k in sorted(_KINETIC_THERMO_DB.keys()):
                rd = _KINETIC_THERMO_DB[k]
                lines.append(f"- `{k}` — {rd['name']}")
            return "\n".join(lines)

        # Determine control regime
        is_kinetic = self._determine_regime(temperature_c, cond)

        lines.append(f"### Reaction: {data['name']}\n")
        lines.append("---\n")

        if is_kinetic:
            lines.append("### 🎯 Prediction: **KINETIC CONTROL** ⚡\n")
        else:
            lines.append("### 🎯 Prediction: **THERMODYNAMIC CONTROL** 🌡️\n")

        # Parameters table
        lines.append("| Parameter | Value | Interpretation |")
        lines.append("|-----------|-------|----------------|")
        regime_text = "**KINETIC CONTROL**" if is_kinetic else "**THERMODYNAMIC CONTROL**"
        interp = "Product determined by relative rates (faster pathway wins)" if is_kinetic else "Product determined by relative stability (most stable product wins)"
        lines.append(f"| Temperature | {temperature_c}°C | {'★ Low → kinetic' if temperature_c < 0 else ('● Moderate' if temperature_c < 50 else '🔥 High → thermodynamic')} |")
        lines.append(f"| Conditions | {cond.title()} | {interp} |")
        lines.append(f"| Regime | {regime_text} | {'Rate-controlled' if is_kinetic else 'Equilibrium-controlled'} |")
        lines.append("")

        # Expected product
        prod = data["kinetic_product"] if is_kinetic else data["thermodynamic_product"]
        lines.append("### Expected Major Product: **{}**\n".format(prod["name"]))
        lines.append(prod["description"] + "\n")
        if prod.get("example"):
            lines.append(f"**Example:** {prod['example']}\n")
        if prod.get("stereochemistry"):
            lines.append(f"**Stereochemistry:** {prod['stereochemistry']}\n")
        if prod.get("reversibility"):
            lines.append(f"**Reversibility:** {prod['reversibility']}\n")

        # Key factor
        lines.append("---\n### 💡 Key Factor\n")
        lines.append(data["key_factor"] + "\n")

        # How to favor each
        lines.append("### 🔄 How to Switch Control:\n")
        if is_kinetic:
            lines.append("**To get Thermodynamic Product instead:**\n")
            tips = data["how_to_favor_thermodynamic"]
        else:
            lines.append("**To get Kinetic Product instead:**\n")
            tips = data["how_to_favor_kinetic"]

        for tip in tips:
            lines.append(f"- {tip}")
        lines.append("")
        lines.append("\n---\n*Analysis based on standard organic chemistry principles.*")
        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        rt = parts[0] if parts else "enolate"
        temp = float(parts[1]) if len(parts) > 1 else 25.0
        cond = parts[2] if len(parts) > 2 else "standard"
        return self._run_base(rt, temp, cond)

    def _match_reaction(self, rt):
        exact = _KINETIC_THERMO_DB.get(rt)
        if exact:
            return exact
        keywords = {
            "conjugate": "conjugated_addition_1_2_vs_1_4", "michael": "conjugated_addition_1_2_vs_1_4",
            "1,2-1,4": "conjugated_addition_1_2_vs_1_4", "14": "conjugated_addition_1_2_vs_1_4",
            "enolate": "enolate_formation_kinetic_vs_thermodynamic",
            "indole": "alkylation_of_indole_C2_vs_C3",
            "da": "diels_allder_endo_vs_exo", "diels-alder": "diels_allder_endo_vs_exo", "cycloaddition": "diels_allder_endo_vs_exo",
            "naphthalene": "sulfonation_of_naphthalene_alpha_vs_beta", "sulfonation": "sulfonation_of_naphthalene_alpha_vs_beta",
            "phenol": "bromination_of_phenol_ortho_vs_para_bromination", "bromination": "bromination_of_phenol_ortho_vs_para_bromination",
            "aldol": "aldol_addition_vs_condensation",
        }
        for kw, val in keywords.items():
            if kw in rt:
                return _KINETIC_THERMO_DB.get(val)
        return None

    def _determine_regime(self, temp_c, cond):
        # Determine based on temperature and explicit conditions
        if cond in ("irreversible", "low_t", "kinetic", "hard_nucleophile"):
            return True
        elif cond in ("reversible", "high_t", "thermodynamic", "soft_nucleophile"):
            return False
        elif cond == "standard":
            # Default: low T → kinetic, high T → thermodynamic
            if temp_c < 0:
                return True
            elif temp_c > 60:
                return False
            else:
                # Intermediate — default to kinetic for most cases
                return True
        return temp_c < 40
