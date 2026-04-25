import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Temperature database for common reaction types
# Format: (temp_min_C, temp_max_C, optimal_C, reasoning, safety_notes)
_TEMP_DB = {
    "SN2": {
        "range": (0, 100), "optimal": (25, 50),
        "reason": "SN2 is fast at room temperature. Elevated temperatures increase rate but may promote elimination side reactions.",
        "safety": "Generally safe at RT. Use reflux condenser if >50°C.",
        "notes": "For highly hindered substrates, lower temps favor SN2 over E2.",
    },
    "SN1": {
        "range": (25, 80), "optimal": (40, 60),
        "reason": "SN1 requires heat to promote ionization and carbocation formation. Too hot → elimination dominates.",
        "safety": "Use reflux apparatus. Protic solvent often used — check BP of solvent.",
        "notes": "Reflux in ethanol/water/acetic acid is typical.",
    },
    "E2": {
        "range": (0, 150), "optimal": (50, 80),
        "reason": "E2 benefits from heat (higher E2/SN2 ratio). Bulky bases like t-BuOK are often used at elevated temperature.",
        "safety": "t-BuOK/t-BuOH systems can reach 80-100°C. Use oil bath with good ventilation.",
        "notes": "Higher temperature favors more substituted (Zaitsev) alkene; lower T favors Hofmann product with bulky bases.",
    },
    "E1": {
        "range": (40, 100), "optimal": (60, 80),
        "reason": "E1 needs sufficient heat for ionization step but excessive heat causes decomposition.",
        "safety": "Reflux conditions typical. Acidic media may be present.",
        "notes": "Often concurrent with SN1 — same conditions apply.",
    },
    "Grignard_formation": {
        "range": (25, 70), "optimal": (35, 50),
        "reason": "Initiation may require gentle warming (35-40°C). Once started, exothermic — control with cooling.",
        "safety": "⚠️ HIGHLY EXOTHERMIC once initiated! Must have ice bath ready. Under N2/Ar atmosphere.",
        "notes": "Start at RT, warm gently to initiate, then control the exotherm. Ether/THF reflux (~35-70°C) is typical.",
    },
    "Grignard_reaction": {
        "range": (-20, 67), "optimal": (0, 25),
        "reason": "Grignard addition to carbonyls is usually done at 0°C to RT to control selectivity and avoid side reactions.",
        "safety": "Anhydrous conditions mandatory. Exothermic upon addition — add Grignard slowly.",
        "notes": "Sensitive substrates: -78°C (dry ice/acetone). Routine additions: 0°C → RT.",
    },
    "organolithium": {
        "range": (-78, 25), "optimal": (-78, 0),
        "reason": "Organolithium reagents are more reactive and less selective than Grignards. Low temperatures improve selectivity.",
        "safety": "⚠️ Pyrophoric! Strictly inert atmosphere. Cryogenic handling required.",
        "notes": "n-BuLi typically added at -78°C for deprotonations. Additions to carbonyls: -78°C → RT slowly.",
    },
    "reduction_NaBH4": {
        "range": (0, 65), "optimal": (0, 25),
        "reason": "NaBH4 reduces aldehydes/ketones rapidly at 0°C to RT. Esters require reflux in MeOH.",
        "safety": "H2 gas evolution! Work in well-ventilated area. No flames.",
        "notes": "MeOH or EtOH as solvent (0°C→RT for C=O; reflux for esters).",
    },
    "reduction_LiAlH4": {
        "range": (0, 67), "optimal": (0, 25),
        "reason": "LiAlH4 reductions are typically done at 0°C then warmed to RT. Reflux in ether/THF for stubborn substrates.",
        "safety": "☠️ EXTREMELY DANGEROUS with water! Quench carefully with ethyl acetate then aqueous Rochelle's salt. H2 evolution.",
        "notes": "Standard: 0°C addition, stir at RT overnight. Reflux only when necessary.",
    },
    "oxidation_PCC": {
        "range": (0, 50), "optimal": (20, 25),
        "reason": "PCC oxidations proceed well at room temperature. Heating can cause over-oxidation or chromium byproducts.",
        "safety": "Cr(VI) compound — toxic, carcinogen. Use gloves, fume hood. Proper disposal required.",
        "notes": "DCM as solvent, RT, 1-12h. Monitor by TLC.",
    },
    "oxidation_Swern": {
        "range": (-78, 25), "optimal": (-78, -60),
        "reason": "Swern oxidation requires low temperature (-78°C) to form the chlorosulfonium intermediate, then warm to RT.",
        "safety": "CO and CO2 gas evolution! Must be in fume hood. Low-temp handling required.",
        "notes": "-78°C (dry ice/acetone): add oxalyl chloride to DMSO, then alcohol, then Et3N. Warm to RT slowly.",
    },
    "oxidation_Jones": {
        "range": (0, 25), "optimal": (0, 20),
        "reason": "Jones oxidation is done at 0-20°C to prevent over-oxidation of intermediates.",
        "safety": "Strongly acidic Cr(VI) solution — corrosive, toxic, oxidizing. Full PPE required.",
        "notes": "Add Jones reagent dropwise to acetone solution at 0°C. Quench with iPrOH.",
    },
    "oxidation_Dess-Martin": {
        "range": (0, 50), "optimal": (20, 25),
        "reason": "Dess-Martin periodinane works at RT. Mild, selective for alcohols → aldehydes without over-oxidation.",
        "safety": "Explosive when dry! Keep wet. Periodinane compound — moderate toxicity.",
        "notes": "DCM, RT, 30min-2h. Very clean reaction — workup is simple Na2S2O3/NaHCO3 wash.",
    },
    "Wittig": {
        "range": (0, 67), "optimal": (0, 25),
        "reason": "Wittig reactions typically run from 0°C to RT. Non-stabilized ylides need lower temps; stabilized ylides tolerate heat.",
        "safety": "Generate ylide under N2. Ph3P=O byproduct can be tedious to remove.",
        "notes": "Non-stabilized: 0°C→RT in THF. Stabilized: reflux in benzene/toluene possible.",
    },
    "Diels-Alder": {
        "range": (25, 200), "optimal": (80, 150),
        "reason": "Diels-Alder is often reversible. Higher temperatures drive endothermic reactions forward but may cause retro-Diels-Alder.",
        "safety": "High-temperature reflux common. Sealed tube for very high T. Pressure buildup risk.",
        "notes": "Cyclopentadiene + maleic anhydride: RT (very fast). Less reactive pairs: reflux in toluene/xylene (110-140°C).",
    },
    "Friedel-Crafts": {
        "range": (25, 80), "optimal": (40, 60),
        "reason": "FC alkylation/acylation needs mild heating to activate Lewis acid catalyst. Too hot → polyalkylation/decomposition.",
        "safety": "Lewis acids (AlCl3, FeCl3) are moisture-sensitive and corrosive. HCl gas evolution possible.",
        "notes": "Anhydrous DCM or nitrobenzene as solvent. AlCl3 (1.1-2 eq). 0°C addition, then RT-reflux.",
    },
    "nitration": {
        "range": (0, 100), "optimal": (0, 50),
        "reason": "Nitrations are highly exothermic. Control temperature to prevent multiple nitrations or decomposition.",
        "safety": "☠️ Mixed acid (HNO3/H2SO4) is HIGHLY corrosive and oxidizing. Can cause severe burns. Run in fume hood behind blast shield.",
        "notes": "Add substrate to cold (0-10°C) mixed acid slowly. Aromatic nitration: 0-50°C. Aliphatic: stricter T control needed.",
    },
    "suzuki_coupling": {
        "range": (50, 120), "optimal": (80, 100),
        "reason": "Pd-catalyzed cross-coupling requires elevated temperature for transmetalation and reductive elimination steps.",
        "safety": "Pd compounds are toxic. Boronic acids generally safe. Use oil bath with condenser.",
        "notes": "Pd(PPh3)4 or Pd(dppf)Cl2. Base (K2CO3, Cs2CO3). Toluene/EtOH/H2O or dioxane/H2O. 80-100°C standard.",
    },
    "Heck_reaction": {
        "range": (80, 140), "optimal": (100, 120),
        "reason": "Heck coupling requires higher temperatures than Suzuki due to the need for migratory insertion into alkene.",
        "safety": "Pd catalyst, high temperature. Use sealed tube or oil bath with condenser.",
        "notes": "Pd(OAc)2 with phosphine ligand. DMF or acetonitrile. 100-120°C, 12-24h.",
    },
    "amide_coupling": {
        "range": (0, 50), "optimal": (20, 25),
        "reason": "Most amide couplings (EDC, HATU, DCC) proceed at RT. Some sluggish couplings may need gentle warming.",
        "safety": "Coupling reagents can be sensitizers. DCC causes allergic reactions in some people.",
        "notes": "DMF or DCM, RT, 2-24h. HATU/DIPEA: faster, more expensive. EDC/HOBt: economical.",
    },
    "peptide_synthesis": {
        "range": (0, 50), "optimal": (20, 25),
        "reason": "SPPS (solid-phase peptide synthesis) uses RT couplings. Deprotection (piperidine) also at RT.",
        "safety": "DMF exposure concerns. Piperidine is corrosive and malodorous. TFA for cleavage is highly corrosive.",
        "notes": "RT coupling (30min-2h per residue). Fmoc deprotection: 20% piperidine/DMF, RT, 5-20min.",
    },
    "hydrogenation": {
        "range": (25, 80), "optimal": (25, 50),
        "reason": "Room temperature hydrogenation is standard. Slightly elevated T may be needed for stubborn substrates.",
        "safety": "⚠️ H2 gas — FLAMMABLE! Use Parr shaker or balloon setup away from ignition sources. Catalyst is pyrophoric when dry.",
        "notes": "Pd/C (5-10%), H2 balloon, RT, 1-24h. EtOAc, EtOH, or MeOH as solvent.",
    },
    "hydroboration": {
        "range": (0, 80), "optimal": (0, 25),
        "reason": "BH3 additions are typically done at 0°C then allowed to warm to RT. Oxidation step with H2O2/NaOH at 0°C.",
        "safety": "BH3·THF is FLAMMABLE and reacts violently with water. H2O2 is strong oxidizer. Handle with care.",
        "notes": "0°C addition of BH3·THF, then RT 1-4h. Then cool to 0°C, add NaOH then H2O2 slowly.",
    },
    "enolate_formation": {
        "range": (-78, 25), "optimal": (-78, 0),
        "reason": "LDA/KHMDS deprotonations require cryogenic temperatures for kinetic control and regioselectivity.",
        "safety": "Pyrophoric bases (LDA, KHMDS). Cryogenic handling. Strictly anhydrous.",
        "notes": "LDA/THF at -78°C (dry ice/acetone) for kinetic enolates. Warmer temps give thermodynamic enolates.",
    },
    "alkylation_enolate": {
        "range": (-78, 67), "optimal": (-78, 0),
        "reason": "Enolate alkylation at low T for kinetic control. Alkyl halide added after enolate formation.",
        "safety": "Depends on base and electrophile. Low-temperature handling.",
        "notes": "Form enolate at -78°C, add alkyl halide, warm slowly to RT.",
    },
    "halogenation_alkene": {
        "range": (-20, 25), "optimal": (0, 20),
        "reason": "Halogen addition to alkenes is rapid at or below RT. Higher temperatures can cause allylic rearrangement or substitution.",
        "safety": "Br2 is highly toxic and corrosive. Cl2 is a toxic gas. Use in fume hood with proper PPE.",
        "notes": "Br2 in DCM at 0°C → RT. Anti-addition stereochemistry. Watch color discharge.",
    },
    "bromination_aromatic": {
        "range": (25, 80), "optimal": (40, 60),
        "reason": "Electrophilic aromatic bromination often requires heat unless activated by Lewis acid (FeBr3).",
        "safety": "Br2 is hazardous. FeBr3 is moisture-sensitive and corrosive.",
        "notes": "FeBr3 catalyst, Br2 in DCM or acetic acid. RT for activated rings, reflux for deactivated.",
    },
    "esterification": {
        "range": (25, 120), "optimal": (60, 80),
        "reason": "Fischer esterification is equilibrium-limited. Heat drives equilibrium toward ester (Le Chatelier). Remove water if possible.",
        "safety": "Conc. H2SO4 catalyst is corrosive. Reflux apparatus needed.",
        "notes": "Acid-catalyzed (H2SO4 or p-TsOH). Reflux with Dean-Stark trap to remove water.",
    },
    "hydrolysis_ester": {
        "range": (25, 100), "optimal": (60, 80),
        "reason": "Ester hydrolysis (basic or acidic) benefits from heat. Saponification: MeOH/H2O reflux.",
        "safety": "NaOH/MeOH reflux is caustic. Acidic hydrolysis uses conc. HCl or H2SO4.",
        "notes": "Basic: NaOH/MeOH-H2O, reflux 1-6h. Acidic: 6M HCl, reflux 2-12h.",
    },
    "protection_alcohol_TBS": {
        "range": (0, 50), "optimal": (0, 25),
        "reason": "TBS protection (TBSCl, imidazole, DMF) proceeds readily at RT.",
        "safety": "TBSCl releases HCl — use base (imidazole, Et3N). DMF penetrates skin.",
        "notes": "DMF or DCM, RT, 2-12h. For sensitive substrates: 0°C.",
    },
    "deprotection_TBS": {
        "range": (0, 50), "optimal": (20, 25),
        "reason": "TBS deprotection with TBAF or acid (HCl/AcOH) works at RT.",
        "safety": "TBAF is basic and corrosive on contact with moisture. HF source concern.",
        "notes": "TBAF/THF, RT, 30min. Or AcOH/THF/H2O (3:1:1), RT, 2-12h.",
    },
}


@ChemMCPManager.register_tool
class TemperatureAdvisor(BaseTool):
    """
    建议反应温度范围的工具。
    内置 35+ 种常见反应类型的温度数据，包含最佳温度范围、安全注意事项和特殊条件说明。
    """
    __version__      = "0.1.0"
    name             = "TemperatureAdvisor"
    func_name        = "advise_temperature"
    description      = "Recommend appropriate temperature ranges for chemical reactions, including optimal range, safety warnings, and special condition notes."
    implementation_description = "Uses embedded database of 35+ reaction types with temperature ranges (min/max/optimal), detailed reasoning, safety notes, and special handling instructions."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Temperature", "Reaction Conditions", "Lab Safety", "Organic Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("reaction_type", "str", "N/A", "Type of reaction (e.g., 'SN2', 'Grignard', 'oxidation', 'suzuki_coupling', 'Diels-Alder')."),
        ("special_conditions", "str", "None", "Special conditions: 'cryogenic', 'large_scale', 'sensitive_substrate', 'pressure'. Use 'None' for default."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'reaction_type [special_condition]'. Example: 'Grignard_formation large_scale'"),
    ]

    output_sig       = [
        ("result", "str", "Detailed temperature recommendation with range, reasoning, safety, and notes."),
    ]

    examples         = [
        {
            "code_input": {"reaction_type": "Grignard_formation", "special_conditions": "None"},
            "text_input": {"input_text": "Grignard_formation"},
            "output": {
                "result": """## Temperature Recommendation: Grignard Formation

### 🌡️ Recommended Range
**0 – 70 °C** | **Optimal: 35 – 50 °C**

### Reasoning
Initiation may require gentle warming (35-40°C). Once started, the reaction is **highly exothermic** — must be controlled with cooling.

### ⚠️ Safety Notes
- ⚠️ **HIGHLY EXOTHERMIC once initiated!** Have ice bath ready before starting.
- Operate under **N2/Ar inert atmosphere** at all times.
- Use **reflux condenser** if temperature exceeds solvent BP.

### Special Instructions
1. Start at room temperature
2. If no initiation, warm gently to 35-40°C (oil bath)
3. Once initiation observed (cloudiness, gentle reflux), **remove heat immediately**
4. Apply cooling bath if needed to maintain 35-50°C
5. Add remaining alkyl halide slowly to control exotherm

### Solvent Boiling Points Reference
- Diethyl ether: 34.6°C (reflux = ~35°C)
- THF: 66°C (reflux = ~66°C)

### Scale Considerations
On large scale (>50 mmol): The exotherm is MORE difficult to control. Consider:
- Lower initial concentration
- Slower addition rate
- More aggressive cooling capacity"""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, reaction_type: str, special_conditions: str = "None") -> str:
        """Core logic: advise temperature for given reaction type."""
        rt = reaction_type.strip()
        sc = special_conditions.strip() if special_conditions and special_conditions.upper() != "NONE" else None

        lines = []
        lines.append(f"## Temperature Recommendation: {rt}\n")

        # Match reaction type
        match_key = self._match_reaction(rt)
        data = _TEMP_DB.get(match_key)

        if data is None:
            lines.append(f"⚠️ Reaction type '{rt}' not found in database.")
            lines.append(f"\n**Available reaction types:**\n")
            for k in sorted(_TEMP_DB.keys()):
                lines.append(f"- `{k}`")
            return "\n".join(lines)

        tmin, tmax = data["range"]
        opt_lo, opt_hi = data["optimal"]

        lines.append("### 🌡️ Recommended Range\n")
        lines.append(f"**{tmin} – {tmax} °C** | **Optimal: {opt_lo} – {opt_hi} °C**\n")

        lines.append("### Reasoning\n")
        lines.append(data["reason"] + "\n")

        lines.append("### ⚠️ Safety Notes\n")
        lines.append(data["safety"] + "\n")

        lines.append("### 📝 Special Instructions\n")
        lines.append(data["notes"] + "\n")

        # Special conditions adjustments
        if sc:
            lines.append("---\n### 🔧 Special Condition Adjustment: `{sc}`\n")
            adjustment = self._get_adjustment(data, sc)
            lines.append(adjustment)

        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        rt = parts[0] if parts else "SN2"
        sc = parts[1] if len(parts) > 1 else "None"
        return self._run_base(rt, sc)

    def _match_reaction(self, rt):
        rt_lower = rt.lower().replace(" ", "_").replace("-", "_")
        exact = _TEMP_DB.get(rt_lower)
        if exact:
            return rt_lower
        for key in _TEMP_DB:
            if rt_lower in key or key in rt_lower:
                return key
        keywords = {
            "sn2": "SN2", "sn1": "SN1", "e2": "E2", "e1": "E1",
            "grignard": "Grignard_formation", "organolithium": "organolithium",
            "reduction": "reduction_NaBH4", "oxidation": "oxidation_PCC",
            "wittig": "Wittig", "diels": "Diels-Alder",
            "friedel": "Friedel-Crafts", "nitration": "nitration",
            "suzuki": "suzuki_coupling", "heck": "Heck_reaction",
            "hydrogenation": "hydrogenation", "hydroboration": "hydroboration",
            "enolate": "enolate_formation", "peptide": "peptide_synthesis",
            "amide": "amide_coupling", "halogenation": "halogenation_alkene",
            "esterification": "esterification", "hydrolysis": "hydrolysis_ester",
            "protection": "protection_alcohol_TBS", "deprotection": "deprotection_TBS",
        }
        for kw, val in keywords.items():
            if kw in rt_lower:
                return val
        return rt_lower

    def _get_adjustment(self, data, condition):
        cl = condition.lower()
        tmin, tmax = data["range"]
        opt_lo, opt_hi = data["optimal"]

        if cl == "cryogenic":
            return """**Cryogenic Protocol:**
- Use dry ice/acetone bath (**-78°C**) for initiation/sensitive steps
- Allow slow warm-up to recommended range
- Ensure proper cryogenic PPE (insulated gloves, face shield)
- Check for CO₂ asphyxiation risk in confined spaces"""
        elif cl == "large_scale":
            return f"""**Large-Scale Adjustments (>50 mmol):**
- Expected ΔT during exotherm is **significantly higher** than small scale
- Reduce initial concentration by **2-3×**
- Addition rate: **10× slower** than lab scale
- Enhanced cooling: **ice-salt bath** (-10 to -5°C) or recirculating chiller
- Temperature monitoring: **Multiple thermocouples** (top, bottom, center)
- Emergency quench protocol: **pre-chilled quenching solution ready**"""
        elif cl == "sensitive_substrate":
            return f"""**Sensitive Substrate Protocol:**
- Start **{max(tmin, -20)}°C** (lower end of range or below)
- Monitor by **TLC every 15-30 minutes**
- If conversion stalls, increase by **5-10°C increments**
- Consider **slow addition** of reagent to minimize local overheating
- Protect from **light/moisture/air** as appropriate"""
        elif cl == "pressure":
            return """**Pressure Reaction Adjustments:**
- Use **sealed tube** or **autoclave** for temperatures above solvent BP
- Calculate expected pressure: ~1 atm per 40°C above ambient for many solvents
- **Never exceed rated pressure** of glassware/equipment
- Use **safety shield** behind blast shield
- Consider **microwave reactor** for precise T/P control"""
        else:
            return f"⚠️ Unknown special condition: '{condition}'. Using standard recommendations."
