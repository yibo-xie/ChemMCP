import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive deprotection conditions database
_DEPROTECTION_DB = {
    # Silyl ethers
    "TBS": {
        "full_name": "tert-Butyldimethylsilyl (TBS/TBDMS)",
        "protected_functional_group": "Alcohol (1°, 2°, 3°)",
        "methods": [
            {
                "name": "TBAF (Tetra-n-butylammonium fluoride)",
                "reagents": "TBAF (1.0M in THF), THF as solvent",
                "conditions": "RT, 30 min – 2 h. Use 1.5–2.0 eq TBAF per silyl group.",
                "mechanism": "Fluoride-induced desilylation (Si-F bond formation is thermodynamically very favorable, K ~10¹⁸). Pentacoordinate silicate intermediate.",
                "pros": ["Fast and reliable", "Works for all silyl groups", "Mild conditions (RT)", "Chemoselective for silyl groups"],
                "cons": ["TBAF is basic (can epimerize stereocenters)", "Moisture-sensitive (contains water)", "Expensive", "Fluoride waste disposal concerns"],
                "compatibility": "❌ Base-labile groups (esters may cleave). ✅ Acetals, benzyl ethers, carbamates OK.",
                "safety": "TBAF is corrosive. HF source concern — handle in fume hood with gloves.",
            },
            {
                "name": "Acetic Acid / Water / THF",
                "reagents": "AcOH/THF/H₂O (3:1:1 v/v/v)",
                "conditions": "RT, 2–12 h. May require gentle heating (40°C) for hindered substrates.",
                "mechanism": "Acid-catalyzed hydrolysis of Si-O bond. Protonation of oxygen → nucleophilic attack by water.",
                "pros": ["Very cheap", "No special reagents needed", "Mild acid conditions", "Compatible with base-sensitive groups"],
                "cons": ["Slower than TBAF", "May not work for very hindered TBS groups", "Not compatible with acid-labile groups (acetals, Boc, THP)", "Long reaction times possible"],
                "compatibility": "✅ Base-labile groups (esters). ❌ Acetals, Boc, THP, other acid-labile PGs.",
                "safety": "Acetic acid is corrosive but low hazard. Standard PPE sufficient.",
            },
            {
                "name": "HF-Pyridine",
                "reagents": "HF·pyridine (70% HF in pyridine complex)",
                "conditions": "0°C → RT, 30 min – 3 h in pyridine or THF.",
                "mechanism": "Same fluoride mechanism as TBAF but with 'naked' HF (more reactive).",
                "pros": ["Very powerful — works when TBAF fails", "Fast even for hindered substrates"],
                "cons": "☠️ **EXTREMELY HAZARDOUS** — HF causes severe burns, systemic toxicity, can be fatal. Requires specialized plastic equipment and calcium gluconate gel on-hand.",
                "compatibility": "Similar to TBAF but harsher. Avoid acid-labile groups.",
                "safety": "🚨 **HF IS EXTREMELY DANGEROUS**. Only use if absolutely necessary. Must have Ca-gluconate gel, specialized training, plastic ware.",
            },
            {
                "name": "HCl in MeOH or dioxane",
                "reagents": "HCl (g) in MeOH, or 4M HCl in dioxane",
                "conditions": "RT, 1–4 h. Concentration dependent.",
                "mechanism": "Acid-catalyzed hydrolysis (stronger than AcOH).",
                "pros": ["Fast", "Cheap", "Volatile (easy to remove)"],
                "cons": ["Strong acid — many PGs not compatible", "May cause elimination side reactions", "HCl gas handling required"],
                "compatibility": "Only use when no other acid-labile groups present.",
                "safety": "Corrosive. Use in fume hood.",
            },
        ],
    },
    "TES": {
        "full_name": "Triethylsilyl (TES)",
        "protected_functional_group": "Alcohol",
        "methods": [
            {
                "name": "TBAF/THF",
                "reagents": "TBAF (1.0M in THF)",
                "conditions": "RT, 5–30 min (faster than TBS due to less steric hindrance).",
                "mechanism": "Fluoride-induced desilylation.",
                "pros": ["Faster than TBS removal", "Same protocol as TBS"],
                "cons": ["Less selective between TES/TBS than desired sometimes", "Basic conditions"],
                "compatibility": "❌ Base-labile groups. ✅ Most others OK.",
                "safety": "Standard TBAF precautions.",
            },
            {
                "name": "AcOH/H₂O/THF",
                "reagents": "AcOH/THF/H₂O (3:1:1)",
                "conditions": "RT, 1–6 h (faster than TBS under same conditions).",
                "mechanism": "Acid hydrolysis.",
                "pros": ["Easier to remove than TBS with acid", "Cheap"],
                "cons": ["Still not compatible with acid-labile groups"],
                "compatibility": "❌ Acetals, Boc, THP.",
                "safety": "Low hazard.",
            },
        ],
    },
    "TMS": {
        "full_name": "Trimethylsilyl (TMS)",
        "protected_functional_group": "Alcohol (often temporary protection)",
        "methods": [
            {
                "name": "K₂CO₃ in MeOH",
                "reagents": "K₂CO₃ (catalytic to 1 eq), MeOH",
                "conditions": "RT, 5–30 min. Often complete during aqueous workup.",
                "mechanism": "Base-catalyzed methanolysis / hydrolysis.",
                "pros": ["Extremely mild", "Buffered (K2CO3 is weak base)", "Often happens spontaneously during workup", "Cheap"],
                "cons": ["Too labile for most synthetic uses (only for temporary protection)"],
                "compatibility": "✅ Very chemoselective — only removes TMS among common silyl groups.",
                "safety": "Very safe.",
            },
            {
                "name": "Acetic Acid / Water",
                "reagents": "80% AcOH in water, or AcOH/THF/H₂O",
                "conditions": "RT, 5–15 min.",
                "mechanism": "Acid hydrolysis.",
                "pros": ["Fast", "Cheap"],
                "cons": ["Even dilute acid removes it — limits utility"],
                "safety": "Safe.",
            },
        ],
    },
    # Amine protecting groups
    "Boc": {
        "full_name": "tert-Butoxycarbonyl (Boc)",
        "protected_functional_group": "Amine (primary, secondary)",
        "methods": [
            {
                "name": "TFA in DCM",
                "reagents": "TFA/DCM (1:1 to 4:1 v/v), optionally with scavenger (water, TIS, anisole)",
                "conditions": "RT, 30 min – 2 h. Monitor by TLC (disappearance of starting material).",
                "mechanism": "Acid-catalyzed elimination: protonation of carbonyl O → tert-butyl cation elimination → CO₂ loss → free amine (as TFA salt).",
                "pros": ["★ Standard method — extremely well-precedented", "Fast (minutes to hours)", "Volatile (TFA and byproducts evaporate)", "Scavengers prevent t-butyl cation alkylation"],
                "cons": ["TFA is highly corrosive and has terrible odor", "Acid-labile groups (acetals, Bn? no Bn is fine, THP, TBS) will also cleave", "TFA salts of amines need basification for free amine", "Isobutylene gas evolution (toxic, flammable)"],
                "compatibility": "❌ TBS, THP, acetals, BOM. ✅ Cbz, Fmoc, Bn, methyl esters, benzyl esters.",
                "safety": "TFA causes severe burns, respiratory damage. Strong fume hood mandatory. Isobutylene is flammable gas.",
            },
            {
                "name": "HCl in dioxane (4M)",
                "reagents": "4 M HCl in 1,4-dioxane",
                "conditions": "RT, 1–3 h. Or reflux for stubborn cases.",
                "mechanism": "Same acidolysis as TFA. Amine obtained as HCl salt.",
                "pros": ["Amine precipitates as HCl salt (easy isolation)", "No evaporation needed (filter solid)", "Clean product", "Dioxane easy to remove"],
                "cons": ["HCl/dioxane solution is carcinogenic (!)", "Dioxane forms peroxides", "Not suitable for acid-labile substrates", "HCl salt needs neutralization for free amine"],
                "compatibility": "Same as TFA — acid-labile groups incompatible.",
                "safety": "⚠️ Dioxane is a suspected carcinogen and forms explosive peroxides. HCl is corrosive.",
            },
            {
                "name": "Me₃SiI (TMS iodide)",
                "reagents": "TMSI in DCM or chloroform",
                "conditions": "0°C → RT, 30 min – 2 h.",
                "mechanism": "Silyl iodide-mediated cleavage: silylation of carbonyl oxygen → iodide attack → elimination.",
                "pros": ["Neutral conditions (no strong acid!)", "Fast", "Can be chemoselective in presence of some acid-labile groups"],
                "cons": ["TMSI is EXPENSIVE and moisture-sensitive", "Highly electrophilic (reacts with many functional groups)", "Iodide byproducts can be hard to remove", "Not commonly used (specialized method)"],
                "compatibility": "More compatible than TFA with some acid-sensitive groups. Still reactive toward nucleophiles.",
                "safety": "TMSI releases HI on contact with moisture. Corrosive. Handle under N2.",
            },
        ],
    },
    "Cbz": {
        "full_name": "Benzyloxycarbonyl (Cbz/Z)",
        "protected_functional_group": "Amine",
        "methods": [
            {
                "name": "Hydrogenolysis (H₂, Pd/C)",
                "reagents": "10% Pd/C (5-20 wt%), H₂ balloon (1 atm) or Parr shaker (30-50 psi), EtOH or EtOAc or MeOH as solvent",
                "conditions": "RT, 1–24 h. Monitor by TLC (UV-active Cbz group disappears).",
                "mechanism": "Pd-catalyzed hydrogenolytic cleavage of benzylic O-C bond: adsorption on Pd surface → H₂ addition → toluene + CO₂ + free amine.",
                "pros": ["Clean byproducts (toluene + CO₂ + amine)", "Catalyst can be filtered off", "Well-established (peptide chemistry standard)", "Mild conditions (RT, atmospheric pressure)"],
                "cons": ["NOT compatible with: C=C (reduction), C≡C (reduction), NO₂ (→ NH₂), aryl halides (dehalogenation), Bn ether (cleaved), benzyl ester (cleaved)", "Pd/C is PYROPHORIC WHEN DRY", "Requires H₂ gas (flammable)", "Slow for sterically hindered substrates"],
                "compatibility": "❌ Alkenes, alkynes, nitro groups, benzyl ethers/esters, other hydrogenolysis-sensitive groups. ✅ TBS, acetals, Boc, Fmoc, t-butyl esters.",
                "safety": "🔥 Keep Pd/C WET! H₂ is FLAMMABLE. Work away from ignition sources. Filter through Celite.",
            },
            {
                "name": "HBr in Acetic Acid (33% or 48%)",
                "reagents": "33% or 48% HBr in acetic acid",
                "conditions": "RT, 30 min – 2 h. Or 0°C for sensitive substrates.",
                "mechanism": "Strong acidolysis: protonation → benzylic cation → cleavage → toluene + CO₂ + amine·HBr salt.",
                "pros": ["Fast", "No H₂/Pd needed (compatible with reducible groups!)","Complete conversion"],
                "cons": ["☠️ VERY HARSH — strongly acidic, oxidizing", "HBr/AcOH is HIGHLY CORROSIVE", "Not compatible with acid-labile groups at all", "May cause racemization of chiral centers", "HBr fumes are toxic and corrosive"],
                "compatibility": "❌ Almost all acid-labile groups. Only use when no other option.",
                "safety": "🚨 HBr/AcOH causes SEVERE burns. Use full face shield, heavy gloves, fume hood behind blast shield.",
            },
            {
                "name": "BCl₃ in DCM",
                "reagents": "BCl₃ (1.0 M in DCM or neat), DCM as solvent",
                "conditions": "-78°C → -30°C or -78°C → RT, 30 min – 2 h. Quench with MeOH at -78°C.",
                "mechanism": "Lewis acid-mediated cleavage: BCl₃ coordinates to carbonyl O → activates toward nucleophilic attack by chloride.",
                "pros": ["Anhydrous conditions (compatible with some groups that can't survive aqueous acid)", "Can be performed at low temperature", "Alternative when H₂/HBr not suitable"],
                "cons": ["BCl₃ reacts VIOLENTLY with water (HCl gas!)", "Expensive", "Low-temperature cryogenic handling needed", "Over-reaction possible"],
                "compatibility": "Better than protic acid for some sensitive substrates. Still a Lewis acid — will affect other Lewis-basic sites.",
                "safety": "BCl₃ releases HCl on quench. Cryogenic handling. Strict anhydrous technique required.",
            },
        ],
    },
    "Fmoc": {
        "full_name": "9-Fluorenylmethoxycarbonyl (Fmoc)",
        "protected_functional_group": "Amine (especially in SPPS)",
        "methods": [
            {
                "name": "Piperidine (20%) in DMF",
                "reagents": "20% piperidine in DMF (v/v)",
                "conditions": "RT, 5–20 min (standard SPPS protocol). For solution-phase: 20–60 min.",
                "mechanism": "β-Elimination: piperidine deprotonates C9 of fluorenyl → dibenzofulvene + CO₂ + free amine (as piperidine salt, which is freely soluble).",
                "pros": ["★ THE standard SPPS deprotection method", "Very fast (5-20 min)", "Mild (no acid, no H₂)", "Dibenzofulvene byproduct is easily removed (reacts with piperidine)", "UV monitoring possible (Fmoc is UV-active)"],
                "cons": ["DMF is a reproductive hazard (use with care)", "Piperidine is corrosive and foul-smelling", "Base-labile groups (esters) may be affected over long exposure", "Dibenzofulvene can alkylate nucleophiles (add scavenger like piperazine)"],
                "compatibility": "❌ Base-labile groups (esters via elimination, Fmoc itself). ✅ Acid-labile groups (Boc, acetals, TBS), Cbz, Bn.",
                "safety": "Piperidine is corrosive. DMF requires gloves (penetrates skin). Good ventilation needed.",
            },
            {
                "name": "DBU (1,8-Diazabicyclo[5.4.0]undec-7-ene) in DMF",
                "reagents": "2% DBU in DMF",
                "conditions": "RT, 1–5 min (very fast!).",
                "mechanism": "Same β-elimination as piperidine but DBU is a stronger, non-nucleophilic base.",
                "pros": ["Extremely fast (minutes)", "Lower loading needed (2% vs 20%)", "Non-nucleophilic (less side reactions than piperidine)"],
                "cons": ["DBU is MORE BASIC → more side reactions with base-sensitive groups", "Aspartimide formation risk (in peptides)", "More expensive than piperidine", "Fast = harder to control for partial deprotection"],
                "compatibility": "More aggressive than piperidine. Avoid with base-sensitive substrates.",
                "safety": "DBU is corrosive. Standard base precautions.",
            },
        ],
    },
    # Ester protecting groups
    "methyl_ester": {
        "full_name": "Methyl Ester (-COOMe)",
        "protected_functional_group": "Carboxylic Acid",
        "methods": [
            {
                "name": "LiOH Saponification",
                "reagents": "LiOH·H₂O (1–5 eq), THF/H₂O (3:1) or MeOH/H₂O",
                "conditions": "RT, 1–12 h. For hindered esters: 40°C or LiOH/THF/H₂O/H₂O₂.",
                "mechanism": "Nucleophilic acyl substitution: OH⁻ attacks carbonyl → tetrahedral intermediate → MeO⁻ leaving → carboxylate.",
                "pros": ["★ Standard method — works for almost all methyl esters", "Mild (RT)", "Li⁺ is beneficial for solubility of organic intermediates", "Clean (LiOH byproducts are water-soluble)"],
                "cons": ["Base-labile groups NOT compatible (enolates, β-ketoesters may decarboxylate)", "Epimerization possible at α-stereocenters", "Over-hydrolysis of other esters if present"],
                "compatibility": "❌ Other esters (non-selective), base-labile groups. ✅ Ethers, acetals, silyl ethers (briefly), carbamates.",
                "safety": "LiOH is caustic. Standard base precautions.",
            },
            {
                "name": "NaOH/MeOH (Zemplén conditions)",
                "reagents": "NaOH (1–2 eq), MeOH or MeOH/H₂O",
                "conditions": "Reflux, 1–6 h. Or RT for activated esters.",
                "mechanism": "Same saponification as LiOH. Methanol as solvent/nucleophile.",
                "pros": ["Very cheap", "Simple setup (reflux in round-bottom flask)", "Well-established"],
                "cons": ["Harsher than LiOH (reflux vs RT)", "NaOH is less soluble in organic/aqueous mixtures than LiOH", "May cause more epimerization"],
                "compatibility": "Same as LiOH — base-labile groups problematic.",
                "safety": "NaOH is caustic. MeOH is flammable. Reflux precautions.",
            },
            {
                "name": "Trimethyltin hydroxide (for base-sensitive substrates)",
                "reagents": "Me₃SnOH, wet DCM or benzene",
                "conditions": "Reflux, 2–12 h.",
                "mechanism": "Organometallic hydrolysis: Sn coordinates to carbonyl, facilitates hydrolysis under near-neutral pH.",
                "pros": ["Near-neutral pH (compatible with base-sensitive groups!)","Chemoselective for esters in presence of amides"],
                "cons": "☠️ **ORGANOTIN COMPOUND — HIGHLY TOXIC**. Tin residues hard to remove. Environmental hazard.",
                "compatibility": "✅ More compatible with base-sensitive groups than LiOH/NaOH.",
                "safety": "☠️ Organotin compounds are NEUROTOXINS. Avoid unless absolutely necessary. Special disposal required.",
            },
        ],
    },
    # Benzyl (ether)
    "Bn": {
        "full_name": "Benzyl Ether (-CH₂Ph)",
        "protected_functional_group": "Alcohol",
        "methods": [
            {
                "name": "Hydrogenolysis (H₂, Pd/C)",
                "reagents": "10% Pd/C (10-20 wt%), H₂ (1 atm or 30-50 psi), EtOH/EtOAc/MeOH",
                "conditions": "RT, 2–24 h. Pearlman's catalyst (Pd(OH)₂/C) is also excellent.",
                "mechanism": "Pd-catalyzed hydrogenolysis of benzylic C-O bond → toluene + alcohol.",
                "pros": ["Clean (toluene is volatile)", "Catalyst filtered off", "Very reliable", "Mild (RT)"],
                "cons": ["Same limitations as Cbz-H₂: no alkenes, nitro, etc.", "Pd/C pyrophoric when dry", "Benzylidene acetals also cleave (both Ph-CH bonds break)"],
                "compatibility": "❌ C=C, C≡C, NO₂, benzyl esters, Cbz, N-Cbz. ✅ TBS, Fmoc, Boc, acetals, PMB (slower), esters.",
                "safety": "Keep Pd/C WET! H₂ flammable.",
            },
            {
                "name": "Lewis Acid (BCl₃ or TMSI)",
                "reagents": "BCl₃ (1M in DCM) or TMSI in DCM/chloroform",
                "conditions": "-78°C → RT, 1–3 h. Quench carefully with MeOH.",
                "mechanism": "Lewis acid coordinates to benzylic oxygen → weakens C-O bond → nucleophilic attack by halide.",
                "pros": ["No H₂ needed (compatible with reducible groups!)", "Can be chemoselective (different rates for different benzyl-type groups)"],
                "cons": ["Harsh Lewis acids — many functional groups affected", "Cryogenic handling", "Expensive reagents", "Byproduct removal tedious"],
                "compatibility": "Check each substrate individually. Generally affects Lewis-basic sites.",
                "safety": "BCl₃ + H₂O → HCl. TMSI + H₂O → HI. Both dangerous.",
            },
        ],
    },
    # Acetal
    "acetal": {
        "full_name": "Acetal / Ketal (e.g., 1,3-dioxolane from ethylene glycol)",
        "protected_functional_group": "Aldehyde / Ketone",
        "methods": [
            {
                "name": "Aqueous Acid (p-TsOH or PPTS)",
                "reagents": "p-TsOH (cat.) or PPTS (cat.) in wet acetone or THF/H₂O",
                "conditions": "RT, 1–6 h. Can accelerate with gentle heating (40°C).",
                "mechanism": "Acid-catalyzed hydrolysis: protonation → hemiacetal → carbonyl + glycol.",
                "pros": ["Very simple setup", "Mild acid (PPTS is very mild)", "Cheap", "Works for most acetals"],
                "cons": ["Not compatible with other acid-labile groups", "Water-sensitive substrates may have issues", "Acetals from ketones (ketals) are slower to hydrolyze"],
                "compatibility": "❌ TBS (may cleave), THP, Boc, BOM. ✅ Bn, Cbz, Fmoc, esters (briefly), silyl ethers (short exposure).",
                "safety": "p-TsOH is corrosive but low hazard. PPTS is milder.",
            },
            {
                "name": "Acetic Acid / Water",
                "reagents": "80% AcOH in water, or AcOH/THF/H₂O (4:1:1)",
                "conditions": "RT, 2–12 h. Or 40°C for ketals.",
                "mechanism": "Acid-catalyzed hydrolysis (weaker acid than p-TsOH → slower but milder).",
                "pros": ["Milder than p-TsOH", "Selective deprotection possible (rate differences between different acetals)", "Cheap"],
                "cons": ["Slow for ketals", "Still acid-labile group incompatibility", "Long reaction times"],
                "compatibility": "Slightly more tolerant than p-TsOH but still incompatible with very acid-labile groups.",
                "safety": "Acetic acid — low hazard.",
            },
        ],
    },
}


@ChemMCPManager.register_tool
class DeprotectionConditions(BaseTool):
    """
    给出脱保护条件的工具。
    内置全面的脱保护条件数据库，涵盖 TBS、TES、TMS、Boc、Cbz、Fmoc、甲基酯、苄基、缩醛等常见保护基团。
    每种方法包含试剂、条件、机理、优缺点、兼容性和安全信息。
    """
    __version__      = "0.1.0"
    name             = "DeprotectionConditions"
    func_name        = "get_deprotection_conditions"
    description      = "Provide detailed deprotection conditions for common protecting groups, including multiple methods with reagents, mechanisms, compatibility, and safety data."
    implementation_description = "Uses embedded database of 9+ protecting groups with 20+ deprotection methods, each containing reagents, conditions, mechanism, pros/cons, compatibility matrix, and safety warnings."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Deprotection", "Protecting Groups", "Organic Synthesis", "Functional Groups"]
    required_envs    = []

    code_input_sig   = [
        ("protecting_group", "str", "N/A", "Name of the protecting group (e.g., 'TBS', 'Boc', 'Cbz', 'Fmoc', 'methyl_ester', 'Bn', 'acetal')."),
        ("method_preference", "str", "all", "Preferred method type: 'mild', 'fast', 'cheap', 'anhydrous', 'neutral', 'all' for all methods."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'protecting_group [preference]'. Example: 'Boc mild' or 'TBS all'"),
    ]

    output_sig       = [
        ("result", "str", "Detailed deprotection conditions with reagents, procedure, mechanism, compatibility, and safety."),
    ]

    examples         = [
        {
            "code_input": {"protecting_group": "Boc", "method_preference": "standard"},
            "text_input": {"input_text": "Boc"},
            "output": {
                "result": """## Deprotection Conditions: Boc (tert-Butoxycarbonyl)

**Protected Functional Group:** Amine (primary, secondary)

---

### Method #1 (Standard): TFA in DCM ⭐

| Item | Detail |
|------|--------|
| **Reagents** | TFA/DCM (1:1 to 4:1 v/v); optional scavenger: H₂O, TIS, anisole |
| **Conditions** | RT, 30 min – 2 h |
| **Workup** | Evaporate TFA (co-evaporate with toluene 2-3×), dissolve residue, basify (NaHCO₃), extract |

**Mechanism:** Acid-catalyzed elimination:
1. Protonation of carbonyl oxygen
2. Elimination of tert-butyl cation (+ isobutylene gas)
3. Decarboxylation (CO₂ release)
4. Free amine as TFA salt

**Compatibility:**
- ✅ Compatible with: Cbz, Fmoc, Bn, methyl/benzyl esters, TMS
- ❌ Incompatible with: TBS, THP, acetals, BOM, other acid-labile PGs

**Safety:** ⚠️ TFA causes severe burns. Use fume hood. Isobutylene gas is flammable.

---

### Method #2 (Alternative): HCl/Dioxane (4M)

| Item | Detail |
|------|--------|
| **Reagents** | 4 M HCl in dioxane |
| **Conditions** | RT, 1–3 h |
| **Workup** | Filter precipitated amine·HCl salt, wash with dioxane |

**Note:** Amine isolated directly as crystalline HCl salt — convenient!

**Safety:** ⚠️ Dioxane is suspected carcinogen. Forms peroxides."""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, protecting_group: str, method_preference: str = "all") -> str:
        """Core logic: provide deprotection conditions."""
        pg = protecting_group.strip()
        pref = method_preference.strip().lower() if method_preference else "all"

        lines = []
        lines.append(f"## Deprotection Conditions: {pg}\n")

        # Match PG name
        pg_key = self._match_pg(pg)
        data = _DEPROTECTION_DB.get(pg_key)

        if data is None:
            lines.append(f"⚠️ Protecting group '{pg}' not found in database.")
            lines.append(f"\n**Available protecting groups:**\n")
            for k in sorted(_DEPROTECTION_DB.keys()):
                d = _DEPROTECTION_DB[k]
                lines.append(f"- `{k}` ({d['full_name']}) — {d['protected_functional_group']}")
            return "\n".join(lines)

        lines.append(f"**Protected Functional Group:** {data['protected_functional_group']}\n")
        lines.append("---\n")

        methods = data.get("methods", [])
        filtered_methods = self._filter_methods(methods, pref)

        for i, method in enumerate(filtered_methods):
            is_std = (i == 0 and pref == "all") or (pref == "standard")
            title = f"### Method #{i+1}{' (⭐ Standard)' if is_std else ''}: {method['name']}\n"
            lines.append(title)

            lines.append("| Item | Detail |")
            lines.append("|------|--------|")
            lines.append(f"| **Reagents** | {method['reagents']} |")
            lines.append(f"| **Conditions** | {method['conditions']} |")
            lines.append("")

            if method.get("workup"):
                lines.append(f"| **Workup** | {method['workup']} |")
                lines.append("")

            lines.append("**Mechanism:** " + method.get("mechanism", "N/A").split(".")[0] + ".\n")

            lines.append("**Compatibility:**")
            lines.append(method.get("compatibility", "N/A"))
            lines.append("")

            lines.append("✅ **Pros:**")
            for p in method["pros"]:
                lines.append(f"- {p}")
            lines.append("\n❌ **Cons:**")
            for c in method["cons"]:
                lines.append(f"- {c}")
            lines.append(f"\n⚠️ **Safety:** {method['safety']}\n")
            lines.append("---\n")

        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        pg = parts[0] if parts else "TBS"
        pref = parts[1] if len(parts) > 1 else "all"
        return self._run_base(pg, pref)

    def _match_pg(self, pg):
        pg_lower = pg.lower()
        mapping = {
            "tbs": "TBS", "tbdms": "TBS", "tbdms": "TBS", "tert-butyldimethylsilyl": "TBS",
            "tes": "TES", "triethylsilyl": "TES",
            "tms": "TMS", "trimethylsilyl": "TMS",
            "boc": "Boc", "t-boc": "Boc", "tert-butoxycarbonyl": "Boc", "t-butoxycarbonyl": "Boc",
            "cbz": "Cbz", "z": "Cbz", "benzyloxycarbonyl": "Cbz",
            "fmoc": "Fmoc", "9-fluorenylmethoxycarbonyl": "Fmoc",
            "methyl_ester": "methyl_ester", "methyl ester": "methyl_ester", "ome": "methyl_ester", "coome": "methyl_ester",
            "bn": "Bn", "benzyl": "Bn", "benzyl ether": "Bn",
            "acetal": "acetal", "ketal": "acetal", "1,3-dioxolane": "acetal",
        }
        exact = mapping.get(pg_lower)
        if exact:
            return exact
        for key in _DEPROTECTION_DB:
            if pg_lower in key.lower() or key.lower() in pg_lower:
                return key
        return pg

    def _filter_methods(self, methods, pref):
        if pref == "all" or pref == "standard":
            return methods
        scored = []
        for m in methods:
            score = 50
            name_lower = m["name"].lower()
            if pref == "mild":
                if any(x in name_lower for x in ["k2co3", "piperidine", "acoh", "acetic"]):
                    score += 30
                elif any(x in name_lower for x in ["hbr", "bcl3", "hf", "tmsi"]):
                    score -= 20
            elif pref == "fast":
                if any(x in m["conditions"].lower() for x in ["min", "5 min", "30 min", "rt"]):
                    score += 20
                if any(x in m["conditions"].lower() for x in ["12 h", "24 h"]):
                    score -= 10
            elif pref == "cheap":
                if any(x in name_lower for x in ["k2co3", "naoh", "tfa", "acetic", "piperidine", "hcl"]):
                    score += 25
                elif any(x in name_lower for x in ["tbaf", "tmsi", "bcl3", "tin"]):
                    score -= 20
            elif pref == "anhydrous":
                if "aqueous" not in name_lower and "water" not in name_lower.lower():
                    score += 20
            elif pref == "neutral":
                if any(x in name_lower for x in ["tmsi", "tin"]):
                    score += 30
                elif any(x in name_lower for x in ["tfa", "hcl", "naoh", "lioh", "ptsa"]):
                    score -= 10
            scored.append((score, m))
        scored.sort(key=lambda x: x[0], reverse=True)
        return [m for s, m in scored]
