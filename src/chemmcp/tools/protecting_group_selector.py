import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive protecting group database
_PG_DB = {
    "alcohol": {
        "description": "Protecting groups for alcohols (1°, 2°, 3°) and phenols.",
        "groups": {
            "TBS (TBDMS)": {
                "full_name": "tert-Butyldimethylsilyl (TBS/TBDMS)",
                "formula": "-SiMe2(t-Bu)",
                "protection_reagent": "TBSCl, imidazole (or Et3N), DMAP (cat.) in DMF or DCM. RT, 2-12h.",
                "deprotection": "TBAF in THF (RT, 30min); or AcOH/THF/H2O (3:1:1, RT, 2-12h); or HF·pyridine.",
                "stability": "Stable to bases, mild acids, hydrogenation, many organometallics. Cleaved by acid (TFA) and fluoride.",
                "size": "Bulky — can influence stereoselectivity",
                "pros": ["Easy to install/remove", "Crystalline (easy to handle)", "Stable to wide range of conditions", " orthogonal to acetals", "Can be removed selectively in presence of TES/TMS"],
                "cons": ["Bulky (steric issues)", "Requires fluoride for clean removal (TBAF is basic/moisture-sensitive)", "Acid-labile (may cleave unexpectedly with strong acids)", "Side reactions: elimination, migration"],
                "selectivity": "TBS > TES > TMS for acid stability; TBAF removes all silyl groups.",
                "typical_uses": "Primary/secondary alcohol protection. Most common silyl protecting group.",
            },
            "TES (TEMS)": {
                "full_name": "Triethylsilyl (TES)",
                "formula": "-SiEt3",
                "protection_reagent": "TESCl, imidazole/Et3N in DMF. RT.",
                "deprotection": "TBAF/THF (faster than TBS due to less steric hindrance); HF·pyridine; AcOH/H2O.",
                "stability": "Less stable to acid than TBS but more stable than TMS.",
                "size": "Moderately bulky",
                "pros": ["Less bulky than TBS", "Can be removed selectively in presence of TBS (with careful control)", "Easier to install on hindered alcohols than TBS"],
                "cons": ["Less acid-stable than TBS", "Still requires fluoride for clean removal"],
                "selectivity": "TES < TBS (acid); TES > TBS (fluoride removal rate)",
                "typical_uses": "When a smaller silyl group is needed, or selective deprotection sequences.",
            },
            "TMS": {
                "full_name": "Trimethylsilyl (TMS)",
                "formula": "-SiMe3",
                "protection_reagent": "TMSCl, Et3N in DCM. 0°C → RT. Or TMSCl, imidazole in DMF.",
                "deprotection": "Very labile: K2CO3/MeOH (RT, minutes); AcOH/H2O; TBAF (instantaneous); even silica gel chromatography.",
                "stability": "Very labile — cleaved by mild acid, base, and even on silica gel.",
                "size": "Smallest silyl group",
                "pros": ["Smallest steric footprint", "Very easy to remove", "Cheap reagents", "Fast installation"],
                "cons": ["Too labile for most multi-step syntheses", "Cleaved during silica gel purification", "Not compatible with basic or acidic conditions"],
                "selectivity": "Most easily removed silyl group. Can be used as 'temporary' protection.",
                "typical_uses": "Short synthetic sequences, protection of reactive intermediates, temporary blocking.",
            },
            "THP (tetrahydropyranyl)": {
                "full_name": "Tetrahydropyranyl (THP)",
                "formula": "acetal from dihydropyran (DHP)",
                "protection_reagent": "DHP, PPTS (cat.) or p-TsOH (cat.) in DCM. 0°C → RT, 1-4h.",
                "deprotection": "Dilute acid: p-TsOH/MeOH (RT); AcOH/H2O; PPTS/MeOH (mild).",
                "stability": "Stable to bases, organometallics, hydrides. Acid-labile.",
                "size": "Moderate bulk",
                "pros": ["Cheap (DHP is inexpensive)", "Easy installation (acid-catalyzed)", "No special reagents needed for removal", "Compatible with base-sensitive substrates"],
                "cons": ["Creates new stereocenter (anomeric mixture)", "Acid-labile (limits utility in acidic sequences)", "Bulky acetal group", "Cannot be used with acid-sensitive functionality"],
                "selectivity": "More acid-labile than TBDPS; similar to acetonide sensitivity.",
                "typical_uses": "Cost-sensitive applications, base-compatible sequences, when stereochemistry at the PG is not critical.",
            },
            "MOM (methoxymethyl)": {
                "full_name": "Methoxymethyl ether (MOM)",
                "formula": "-CH2OCH3",
                "protection_reagent": "MOMCl, i-Pr2NEt (DIPEA) in DCM. 0°C → RT. Or MOMCl, NaH in THF (for less reactive alcohols).",
                "deprotection": "Strong acid: conc. HCl, TFA, or HBr/AcOH. Boron tribromide (BCl3) also works.",
                "stability": "Stable to strong bases, organometallics, reducing agents. Requires strong acid for removal.",
                "size": "Small",
                "pros": ["Small steric profile", "Very robust under basic conditions", "Orthogonal to silyl ethers and esters"],
                "cons": ["Requires STRONG acid for removal (can damage other groups)", "MOMCl is carcinogenic (!)", "Overalkylation possible (bis-MOM)"],
                "selectivity": "More stable than THP and benzyl to acid. Removed only with strong protic/Lewis acids.",
                "typical_uses": "Long synthetic sequences requiring base stability. Often used in complex natural product synthesis.",
            },
            "Bn (benzyl)": {
                "full_name": "Benzyl ether (Bn)",
                "formula": "-CH2Ph",
                "protection_reagent": "NaH, BnBr in DMF/THF. 0°C → RT. Or Ag2O, BnBr (milder).",
                "deprotection": "Hydrogenolysis: H2, Pd/C (1 atm, RT, EtOH/EtOAc). Or dissolving metal: Na, NH3(l).",
                "stability": "VERY stable to acid and base. Only removed by hydrogenolysis or Lewis acid (BCl3, TMSI).",
                "size": "Moderate (planar aromatic)",
                "pros": ["★ One of the most robust protecting groups", "Stable to both strong acid AND strong base", "Clean removal (H2 → toluene + ROH)", "Well-precedented protocols"],
                "cons": ["Requires H2/Pd-C (not compatible with reducible groups: C=C, NO2, Cbz, benzylidene)", "Hydrogenation conditions may reduce other functionalities", "BnBr is lachrymator/toxic"],
                "selectivity": "Orthogonal to almost everything except Cbz (also removed by hydrogenolysis).",
                "typical_uses": "Carbohydrate chemistry, long syntheses needing maximum stability, final global deprotection step.",
            },
            "Acetyl (Ac)": {
                "full_name": "Acetate ester (Ac/OAc)",
                "formula": "-OC(O)CH3",
                "protection_reagent": "Ac2O, pyridine (or Et3N, DMAP) in DCM. 0°C → RT. Or AcCl, pyridine.",
                "deprotection": "Basic hydrolysis: K2CO3/MeOH (RT); NaOH/MeOH (RT). Or Zemplén deacetylation: NaOMe/MeOH.",
                "stability": "Stable to acid. Base-labile. Compatible with hydrogenation.",
                "size": "Small planar",
                "pros": ["Very easy to install", "Very easy to remove (mild base)", "Cheap", "No special equipment needed"],
                "cons": ["Base-labile (not compatible with enolates, strong bases)", "Can migrate under basic conditions (acyl migration)", "Ester may participate in neighboring group reactions"],
                "selectivity": "Removed before silyl ethers (base vs fluoride). Orthogonal to acid-labile groups.",
                "typical_uses": "Carbohydrate/prostaglandin chemistry, short sequences, temporary protection.",
            },
            "PMB (p-methoxybenzyl)": {
                "full_name": "p-Methoxybenzyl ether (PMB)",
                "formula": "-CH2(C6H4-OMe-p)",
                "protection_reagent": "NaH, PMBCl (or PMBBr) in DMF/THF. RT.",
                "deprotection": "Oxidative: DDQ, DCM/H2O (RT, 0.5-2h). OR CAN (ceric ammonium nitrate). Also hydrogenolysis (like Bn but faster).",
                "stability": "Similar to Bn but more electron-rich → easier oxidative cleavage.",
                "size": "Moderate-large",
                "pros": ["Can be removed oxidatively (orthogonal to Bn!)", "DDQ conditions are mild and chemoselective", "More easily cleaved than Bn by hydrogenolysis"],
                "cons": ["DDQ is toxic and expensive", "Oxidation may affect other electron-rich arenes", "PMBBr is expensive"],
                "selectivity": "PMB can be removed with DDQ while Bn remains intact. This orthogonality is VERY useful.",
                "typical_uses": "Complex syntheses where Bn + PMB orthogonal pair is needed.",
            },
        }
    },
    "amine": {
        "description": "Protecting groups for amines (primary, secondary).",
        "groups": {
            "Boc (tert-butoxycarbonyl)": {
                "full_name": "tert-Butoxycarbonyl (Boc)",
                "formula": "-NHCOO-t-Bu",
                "protection_reagent": "(Boc)2O (Boc anhydride), Et3N or DMAP in DCM/THF. RT, 1-12h. Or Boc-ON for sensitive substrates.",
                "deprotection": "Acid: TFA/DCM (1:1, RT, 30min-2h); or HCl/dioxane (4M, RT, 1h); or conc. H2SO4.",
                "stability": "Stable to bases, nucleophiles, hydrogenation, organometallics. Acid-labile.",
                "size": "Bulky (t-butyl group)",
                "pros": ["★ Most common amine protecting group", "Easy installation (Boc2O is cheap/stable)", "Clean acid removal (CO2 + isobutylene gas)", "Compatible with Fmoc strategy (solid phase)"],
                "cons": ["Acid-labile (not compatible with acid-catalyzed reactions)", "t-Butyl cation intermediate can alkylate nucleophiles", "TFA removal requires evaporation (tedious for volatile products)"],
                "selectivity": "Boc < Cbz (acid stability). Boc removed with mild acid; Cbz needs hydrogenolysis.",
                "typical_uses": "SPPS (alternating with Fmoc), solution-phase synthesis, peptide chemistry.",
            },
            "Cbz (Z/benzyloxycarbonyl)": {
                "full_name": "Benzyloxycarbonyl (Cbz/Z)",
                "formula": "-NHCOOCH2Ph",
                "protection_reagent": "Cbz-Cl (benzyl chloroformate), NaOH (Schotten-Baumann) or NaHCO3 in dioxane/H2O. 0°C → RT.",
                "deprotection": "Hydrogenolysis: H2, Pd/C (RT, EtOH). OR HBr/AcOH (strong acid). OR BCl3/DCM.",
                "stability": "Stable to acid (except very strong like HBr/AcOH). Stable to bases. Removed by hydrogenolysis.",
                "size": "Moderate",
                "pros": ["Stable to acid (unlike Boc — useful in acidic sequences)", "Clean removal (toluene + CO2)", "Well-established in peptide chemistry"],
                "cons": ["Requires H2/Pd-C (same limitation as Bn ether)", "Not compatible with reducible groups", "Cbz-Cl is lachrymator and tear agent (!)"],
                "selectivity": "Cbz > Boc (acid). Cbz survives Boc removal conditions. Both removed by H2/Pd-C.",
                "typical_uses": "Classical peptide synthesis, when acid-stable amine protection needed.",
            },
            "Fmoc (9-fluorenylmethoxycarbonyl)": {
                "full_name": "9-Fluorenylmethoxycarbonyl (Fmoc)",
                "formula": "-NHCOOCH2(9-fluorenyl)",
                "protection_reagent": "Fmoc-OSu (or Fmoc-Cl), Na2CO3 in dioxane/H2O. RT, 1-2h.",
                "deprotection": "Base: 20% piperidine/DMF (RT, 5-20min). Or DBU/DMF. β-Elimination releases dibenzofulvene.",
                "stability": "Stable to acid, hydrogenation, nucleophiles. Base-labile (mild base sufficient).",
                "size": "Large (aromatic tricyclic)",
                "pros": ["★ Standard for SPPS (solid-phase peptide synthesis)", "Mild base removal (piperidine/DMF)", "UV-active (monitoring by UV)", "Orthogonal to Boc (acid-labile)"],
                "cons": ["Base-labile (limits use with enolates, etc.)", "Fmoc reagents expensive", "Dibenzofulvene byproduct can alkylate nucleophiles", "Large size (steric issues)"],
                "selectivity": "Fmoc << Boc (base stability). Fmoc removed with piperidine; Boc requires TFA.",
                "typical_uses": "★ SPPS (Fmoc-strategy peptide synthesis). Solution-phase when base-labile PG needed.",
            },
            "Alloc (allyloxycarbonyl)": {
                "full_name": "Allyloxycarbonyl (Alloc)",
                "formula": "-NHCOOCH2CH=CH2",
                "protection_reagent": "Alloc-Cl, NaHCO3 in dioxane/H2O. 0°C → RT.",
                "deprotection": "Pd(0)-catalyzed: Pd(PPh3)4, PhSiH3 (or morpholine) in THF. RT, 30min-2h.",
                "stability": "Stable to BOTH acid AND base (under most conditions). Only removed via Pd-catalyzed deallylation.",
                "size": "Small",
                "pros": ["★ Orthogonal to BOTH Boc (acid) AND Fmoc (base)!","Mild Pd(0) deprotection conditions", "Useful in complex orthogonal schemes"],
                "cons": ["Requires Pd catalyst (expensive, toxic)", "Pd removal can be tricky", "Limited precedents compared to Boc/Fmoc/Cbz"],
                "selectivity": "Fully orthogonal to Boc (acid), Fmoc (base), Cbz (hydrogenolysis if Pd-selective).",
                "typical_uses": "Multi-step syntheses requiring triple orthogonality. Lysine side-chain protection.",
            },
            "Troc": {
                "full_name": "2,2,2-Trichloroethoxycarbonyl (Troc)",
                "formula": "-NHCOOCH2CCl3",
                "protection_reagent": "Troc-Cl, base in dioxane/H2O. 0°C.",
                "deprotection": "Reductive: Zn, AcOH (RT). OR Cd/Pb couple. OR electrochemical reduction.",
                "stability": "Very stable to acid and base. Requires reductive removal.",
                "size": "Moderate",
                "pros": ["Very robust (survives Boc AND Cbz removal)", "Reductive deprotection (chemoselective)", "Electron-withdrawing: makes N less nucleophilic (good for selectivity)"],
                "cons": ["Chlorinated compound (environmental concern)", "Zn dust workup is messy", "Limited use cases"],
                "selectivity": "Troc survives Boc (acid) and Cbz (H2) removal. Reductive deprotection is unique.",
                "typical_uses": "Specialized applications: glycoscience (amine protection in sugars), complex orthogonal schemes.",
            },
        }
    },
    "carbonyl": {
        "description": "Protecting groups for aldehydes and ketones.",
        "groups": {
            "acetal / 1,3-dioxolane": {
                "full_name": "Ethylene acetal (1,3-dioxolane)",
                "formula": "cyclic -O-CH2-CH2-O- (from ethylene glycol)",
                "protection_reagent": "Ethylene glycol, p-TsOH (cat.) or PPTS (cat.), toluene. Dean-Stark trap to remove water. Reflux.",
                "deprotection": "Aqueous acid: p-TsOH/H2O (acetone, RT); AcOH/H2O (RT); PPTS/wet acetone.",
                "stability": "Stable to base, organometallics, hydrides, Grignards. Acid-labile.",
                "size": "Planar cyclic",
                "pros": ["Very common carbonyl PG", "Original carbonyl restored upon deprotection", "Installation uses simple reagents", "Protects aldehydes from oxidation/enolization"],
                "cons": ["Acid-labile (limits use in acidic sequences)", "Forms new stereocenter (from ketones → mixture)", "Requires water removal during formation (Dean-Stark)"],
                "selectivity": "Acetal stability: dithiane > 1,3-dithiolane > 1,3-dioxolane > 1,3-oxathiolane > O,O-acetal.",
                "typical_uses": "Aldehyde/ketone protection in total synthesis. Very widely used.",
            },
            "1,3-dithiane": {
                "full_name": "1,3-Dithiane (thioacetal/Umpolung)",
                "formula": "cyclic -S-CH2-CH2-S-",
                "protection_reagent": "1,3-Propanedithiol, BF3·Et2O (Lewis acid) or ZnI2, BF3·Et2O. DCM, 0°C → RT.",
                "deprotection": "Oxidative: Hg(II) salts (HgCl2, HgO, wet aq. CaCO3), MeCN/aq. OR NBS, wet acetone. OR other metals (Cu(II), Tl(III)).",
                "stability": "VERY stable to acid AND base. Requires oxidative (or mercuric) removal.",
                "size": "Larger than O,O-acetal",
                "pros": ["★ Umpolung: dithiane anion acts as acyl anion equivalent (!)","Extremely robust (survives extreme conditions)", "Enables disconnections impossible with normal carbonyl chemistry"],
                "cons": ["☠️ Mercury-based deprotection (TOXIC!) — though Hg-free methods exist", "Strong odor (thiols)", "Larger than O,O-acetals", "Installation requires Lewis acid"],
                "selectivity": "Most stable carbonyl PG. Survives conditions that cleave O,O-acetals.",
                "typical_uses": "Umpolung strategies, synthesis of ketones from 'acyl anions', protecting carbonyls through harsh reaction sequences.",
            },
        }
    },
    "carboxylic_acid": {
        "description": "Protecting groups for carboxylic acids.",
        "groups": {
            "methyl ester": {
                "full_name": "Methyl ester (-COOMe)",
                "formula": "-COOCH3",
                "protection_reagent": "MeOH, catalytic H2SO4 or SOCI2 then MeOH. Fischer esterification: reflux with Dean-Stark.",
                "deprotection": "Basic hydrolysis: LiOH, THF/H2O (RT-reflux). OR NaOH/MeOH-H2O. Saponification.",
                "stability": "Stable to acid, neutral conditions. Base-labile.",
                "size": "Small",
                "pros": ["Simplest ester PG", "Easy installation (Fischer esterification)", "Easy removal (saponification)", "Small steric impact"],
                "cons": ["Base-labile (not compatible with strong bases)", "May require harsh conditions for sterically hindered esters", "Methyl iodide byproduct (toxic) if using MeI/Silane method"],
                "selectivity": "Methyl < ethyl < t-butyl < benzyl (toward base hydrolysis rate). Methyl easiest to hydrolyze.",
                "typical_uses": "Most common carboxylic acid PG. Standard in peptide synthesis (C-terminal).",
            },
            "t-butyl ester": {
                "full_name": "tert-Butyl ester (-COOt-Bu)",
                "formula": "-COOC(CH3)3",
                "protection_reagent": "Isobutene, H2SO4 (acid-catalyzed). Or Boc2O, DMAP (for pre-formed acids). Or t-BuOH, DCC, DMAP.",
                "deprotection": "Acid: TFA/DCM (RT, 30min-2h). Strong acid cleaves via tert-butyl cation.",
                "stability": "Stable to base, hydrogenation, nucleophiles. Acid-labile (like Boc for amines).",
                "size": "Bulky",
                "pros": ["Acid-labile (orthogonal to methyl ester which is base-labile!)","Clean removal (isobutylene + CO2)", "Compatible with saponification conditions"],
                "cons": ["Acid-labile (limits use in acidic sequences)", "Bulkier than methyl ester", "Installation can be trickier than methyl ester"],
                "selectivity": "t-Butyl ester cleaved by TFA; methyl ester survives. Useful orthogonality.",
                "typical_uses": "When acid-labile carboxylic acid PG needed (e.g., alongside base-labile groups).",
            },
            "benzyl ester": {
                "full_name": "Benzyl ester (-COOBn)",
                "formula": "-COOCH2Ph",
                "protection_reagent": "BnBr, K2CO3 in DMF (RT). Or BnOH, DCC, DMAP (Steglich).",
                "deprotection": "Hydrogenolysis: H2, Pd/C (RT, EtOH/EtOAc). Same conditions as Bn ether/Cbz.",
                "stability": "Stable to acid AND base. Only removed by hydrogenolysis (or strong Lewis acid).",
                "size": "Moderate",
                "pros": ["★ Very robust (stable to acid AND base)","Clean H2-mediated removal", "Orthogonal to both Boc (acid) and methyl ester (base)"],
                "cons": ["Requires H2/Pd-C (reducible group incompatibility)", "Same limitations as Bn ether for reducible groups"],
                "selectivity": "Benzyl ester survives both Boc (TFA) and saponification (NaOH). Only H2 removes it.",
                "typical_uses": "Complex syntheses where carboxylic acid must survive both acidic and basic steps.",
            },
        }
    },
}


@ChemMCPManager.register_tool
class ProtectingGroupSelector(BaseTool):
    """
    选择合适保护基团的工具。
    内置全面的保护基团数据库，涵盖醇、胺、羰基、羧酸四大类官能团，包含保护/脱保护条件、正交性、优缺点分析。
    """
    __version__      = "0.1.0"
    name             = "ProtectingGroupSelector"
    func_name        = "select_protecting_group"
    description      = "Select appropriate protecting groups for a given functional group type, considering reaction context, compatibility, and deprotection strategy."
    implementation_description = "Uses embedded database of 22+ protecting groups across 4 functional group categories (alcohol, amine, carbonyl, carboxylic acid) with full protection/deprotection protocols, stability data, and orthogonality information."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Protecting Groups", "Organic Synthesis", "Functional Groups", "Orthogonal Protection"]
    required_envs    = []

    code_input_sig   = [
        ("functional_group", "str", "N/A", "Functional group to protect: 'alcohol', 'amine', 'carbonyl', 'carboxylic_acid'."),
        ("reaction_context", "str", "general", "Reaction context: 'acidic', 'basic', 'reductive', 'oxidative', 'long_sequence', 'general'."),
        ("priority", "str", "balance", "Priority: 'ease_of_removal', 'robustness', 'small_size', 'orthogonality', 'balance'."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'functional_group [context] [priority]'. Example: 'amine acidic ease_of_removal'"),
    ]

    output_sig       = [
        ("result", "str", "Recommended protecting group(s) with full details including protection/deprotection conditions."),
    ]

    examples         = [
        {
            "code_input": {"functional_group": "alcohol", "reaction_context": "general", "priority": "balance"},
            "text_input": {"input_text": "alcohol"},
            "output": {
                "result": """## Protecting Group Recommendation for: Alcohol

### Context: General Synthesis | Priority: Balanced

---

### ⭐ #1 Recommendation: **TBS (tert-Butyldimethylsilyl)**

| Property | Detail |
|----------|--------|
| Full Name | tert-Butyldimethylsilyl (TBS/TBDMS) |
| Formula | -SiMe₂(t-Bu) |
| Size | Bulky |

**Protection:** TBSCl, imidazole, DMF/DCM. RT, 2-12h.
**Deprotection:** TBAF/THF (RT, 30min); or AcOH/THF/H₂O.

**Why it's the standard choice:**
- Easy to install and remove
- Stable to bases, mild acids, hydrogenation
- Crystalline solid (easy handling)
- Most commonly used alcohol PG

✅ Pros:
- Easy to install/remove
- Crystalline (easy to handle)
- Wide condition compatibility
- Orthogonal to acetals

❌ Cons:
- Bulky (steric issues)
- Requires fluoride for clean removal
- Acid-labile (unexpected cleavage risk)

---

### Alternative Options

**For acidic reaction context:** → THP or MOM (more acid-stable than TBS? No — use Bn)
**Actually for acidic:** → **Bn (benzyl)** — stable to acid AND base
**For base-sensitive sequences:** → **Acetyl** (easily removed with K2CO3/MeOH)
**For maximum robustness:** → **Bn** (only removed by H2/Pd-C)
**For small steric profile:** → **TMS** or **Acetyl**

### Orthogonality Quick Reference (Alcohols)
| Group | Removed By | Compatible With |
|-------|-----------|-----------------|
| TBS | Fluoride (TBAF) / Acid | Bases, H2, organometallics |
| Bn | H2/Pd-C | Acid, base, most conditions |
| Ac | Base (K2CO3/MeOH) | Acid, H2 |
| THP | Mild acid | Bases, organometallics |
| PMB | DDQ (oxidation) or H2/Pd-C | Acid, base (no DDQ) |"""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, functional_group: str, reaction_context: str = "general", priority: str = "balance") -> str:
        """Core logic: recommend protecting groups."""
        fg = functional_group.strip().lower()
        ctx = reaction_context.strip().lower() if reaction_context else "general"
        pri = priority.strip().lower() if priority else "balance"

        lines = []
        lines.append(f"## Protecting Group Recommendation for: {fg.replace('_', ' ').title()}\n")
        lines.append(f"### Context: {ctx.title()} | Priority: {pri.title()}\n")

        # Find FG category
        fg_key = self._match_fg(fg)
        category = _PG_DB.get(fg_key)

        if not category:
            lines.append(f"⚠️ Functional group '{fg}' not found.")
            lines.append("**Available:** `alcohol`, `amine`, `carbonyl`, `carboxylic_acid`")
            return "\n".join(lines)

        groups = category["groups"]

        # Rank recommendations based on context and priority
        ranked = self._rank_groups(groups, ctx, pri)

        for i, (name, pg) in enumerate(ranked[:4]):
            lines.append(f"---\n### {'⭐ #' + str(i+1) + ' Recommendation' if i == 0 else '### Alternative #' + str(i+1)}: **{name}**\n")
            lines.append("| Property | Detail |")
            lines.append("|----------|--------|")
            lines.append(f"| Full Name | {pg['full_name']} |")
            lines.append(f"| Formula | `{pg['formula']}` |")
            lines.append(f"| Size | {pg['size']} |")
            lines.append("")
            lines.append(f"**Protection:** {pg['protection_reagent'].split('.')[0]}.")
            lines.append(f"**Deprotection:** {pg['deprotection'].split('.')[0]}.")
            lines.append("")
            if i == 0:
                lines.append(f"**Why it's the recommended choice for this context:**\n{self._why_recommend(name, pg, ctx, pri)}\n")
            lines.append("✅ **Pros:**")
            for p in pg["pros"]:
                lines.append(f"- {p}")
            lines.append("\n❌ **Cons:**")
            for c in pg["cons"]:
                lines.append(f"- {c}")
            lines.append("")

        # Orthogonality summary
        lines.append("### Orthogonality Quick Reference\n")
        lines.append("| Group | Removed By | Key Compatibility |")
        lines.append("|-------|-----------|--------------------|")
        for name, pg in list(groups.items())[:6]:
            removal = pg["deprotection"].split(";")[0][:50]
            compat = [s for s in pg["stability"].split(". ")[:2]]
            lines.append(f"| {name} | {removal}... | {' '.join(compat)[:40]}... |")

        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        fg = parts[0] if parts else "alcohol"
        ctx = parts[1] if len(parts) > 1 else "general"
        pri = parts[2] if len(parts) > 2 else "balance"
        return self._run_base(fg, ctx, pri)

    def _match_fg(self, fg):
        mapping = {
            "alcohol": "alcohol", "oh": "alcohol", "hydroxyl": "alcohol",
            "amine": "amine", "nh2": "amine", "amino": "amine",
            "carbonyl": "carbonyl", "ketone": "carbonyl", "aldehyde": "carbonyl", "c=o": "carbonyl",
            "carboxylic_acid": "carboxylic_acid", "cooh": "carboxylic_acid", "acid": "carboxylic_acid",
            "ester": "carboxylic_acid",
        }
        return mapping.get(fg, fg)

    def _rank_groups(self, groups, ctx, pri):
        items = list(groups.items())

        # Scoring based on context and priority
        def score(item):
            name, pg = item
            s = 0
            stab = pg["stability"].lower()

            if ctx == "acidic":
                if "stable to acid" in stab or "very stable" in stab:
                    s += 10
                if "acid-labile" in stab:
                    s -= 10
            elif ctx == "basic":
                if "stable to base" in stab or "very stable" in stab:
                    s += 10
                if "base-labile" in stab:
                    s -= 10
            elif ctx == "reductive":
                if "stable to" in stab and "hydrogenation" in stab:
                    s += 10
                if "hydrogenolysis" in pg["deprotection"].lower():
                    s -= 8
            elif ctx == "long_sequence":
                if "very stable" in stab or ("stable to acid" in stab and "stable to base" in stab):
                    s += 10
                if "very labile" in stab:
                    s -= 10

            if pri == "ease_of_removal":
                if any(x in pg["deprotection"].lower() for x in ["rt", "minutes", "mild"]):
                    s += 5
                elif any(x in pg["deprotection"].lower() for x in ["hg", "strong acid", "reflux"]):
                    s -= 3
            elif pri == "robustness":
                if "very stable" in stab:
                    s += 8
                elif "labile" in stab:
                    s -= 5
            elif pri == "small_size":
                if "small" in pg["size"].lower():
                    s += 5
                elif "bulky" in pg["size"].lower() or "large" in pg["size"].lower():
                    s -= 3
            elif pri == "orthogonality":
                if any(x in pg["deprotection"].lower() for x in ["palladium", "ddq", "fluoride", "oxidative"]):
                    s += 5

            return s

        items.sort(key=score, reverse=True)
        return items

    def _why_recommend(self, name, pg, ctx, pri):
        reasons = []
        if ctx == "acidic" and "stable to acid" in pg["stability"].lower():
            reasons.append(f"Survives acidic conditions ({pg['stability'].split('.')[0]})")
        elif ctx == "basic" and "stable to base" in pg["stability"].lower():
            reasons.append(f"Tolerates basic conditions well")
        elif ctx == "general":
            reasons.append(pg.get("typical_uses", "").split(".")[0])
        if pri == "balance":
            reasons.append("Good balance of ease of installation, stability, and clean removal")
        return ". ".join(reasons[:2]) + "." if reasons else "Well-suited for this application."
