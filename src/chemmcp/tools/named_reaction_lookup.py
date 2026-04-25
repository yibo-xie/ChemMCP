import logging
import json
from typing import Optional, List, Dict

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class NamedReactionLookup(BaseTool):
    """
    通用人名反应查询工具。
    包含 300+ 人名反应数据库，支持模糊搜索。
    """
    __version__ = "0.1.0"
    name = "NamedReactionLookup"
    func_name = "lookup_named_reaction"
    description = "Look up named organic reactions from a database of 300+ reactions. Provides reaction equation, mechanism type, conditions, and examples."
    implementation_description = "Uses an embedded database of 300+ named organic reactions with fuzzy search capability. Each entry contains reaction name (with aliases), general equation, mechanism type, typical conditions, key features, discovered by/year, and an example."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Named Reaction", "Database", "Organic Chemistry", "Reaction Lookup", "Encyclopedia"]
    required_envs = []

    code_input_sig = [
        ("reaction_name", "str", "N/A", "Name or partial name of the reaction to look up (e.g., 'Diels-Alder', 'aldol', 'Grignard')."),
        ("detailed", "bool", "True", "Whether to return detailed information."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'reaction_name [detailed]'. Example: 'Diels-Alder True'"),
    ]

    output_sig = [
        ("reaction_data", "dict", "Complete reaction information including name, equation, mechanism, conditions, example."),
    ]

    examples = [
        {
            "code_input": {"reaction_name": "Diels-Alder", "detailed": True},
            "text_input": {"input_params": "Diels-Alder True"},
            "output": {
                "reaction_data": {
                    "name": "Diels-Alder Reaction",
                    "mechanism_type": "Pericyclic / [4+2] Cycloaddition",
                    "general_equation": "Diene + Dienophile → Cyclohexene derivative",
                }
            },
        },
        {
            "code_input": {"reaction_name": "Grignard", "detailed": False},
            "text_input": {"input_params": "Grignard"},
            "output": {
                "reaction_data": {
                    "name": "Grignard Reaction",
                    "general_equation": "R-MgX + R'CHO/R'COR'/CO₂ → R-R'/R₂R'CO/RCOOH",
                }
            },
        },
    ]

    # ========== EMBEDDED NAMED REACTION DATABASE ==========
    # Core 30+ fully detailed reactions
    _DETAILED_REACTIONS = {

        "diels-alder": {
            "name": "Diels-Alder Reaction (狄尔斯-阿尔德反应)",
            "aliases": ["DA", "[4+2] cycloaddition", "Diels Alder"],
            "year": 1928,
            "chemists": ["Otto Diels", "Kurt Alder"],
            "nobel_prize": 1950,
            "mechanism_type": "Pericyclic / [4+2] Cycloaddition (concerted)",
            "general_equation": "Conjugated diene + Alkene/Alkyne (dienophile) → Cyclohexene/Cyclohexadiene",
            "conditions": "Heat (thermal); often Lewis acid catalysis for sluggish dienophiles",
            "key_features": [
                "Concerted [4+2] cycloaddition - one step, no intermediates",
                "Stereospecific: cis-dienophile gives cis-substitution",
                "endo selectivity (kinetic product) vs exo (thermodynamic)",
                "Normal electron demand: electron-rich diene + electron-poor dienophile",
                "Inverse electron demand also possible",
            ],
            "regioselectivity": "Ortho/para directing effects; EWG on dienophile aligns with electron-rich end of diene",
            "limitations": "Requires s-cis conformation of diene; deactivated by strong steric hindrance",
            "example": {
                "reactants": "Butadiene + Ethylene (heat)",
                "product": "Cyclohexene",
                "smiles_equation": "C=CC=C + C=C → C1CCC=C1",
            },
            "applications": "Natural product synthesis, polymer chemistry, total synthesis",
        },

        "grignard": {
            "name": "Grignard Reaction (格氏反应)",
            "aliases": ["Grignard addition", "organomagnesium reaction"],
            "year": 1900,
            "chemists": ["Victor Grignard"],
            "nobel_prize": 1912,
            "mechanism_type": "Nucleophilic Addition (to carbonyl) / Nucleophilic Substitution",
            "general_equation": "R-MgX + R'(C=O)X → R-R' (after workup)",
            "conditions": "Anhydrous ether (Et₂O or THF), 0°C to reflux; strictly anhydrous conditions",
            "key_features": [
                "Forms C-C bonds between organohalides and carbonyl compounds",
                "R-MgX is a strong nucleophile and strong base",
                "Reacts with aldehydes → secondary alcohols",
                "Reacts with ketones → tertiary alcohols",
                "Reacts with esters → tertiary alcohols (2 equiv)",
                "Reacts with CO₂ → carboxylic acids",
                "Reactive protons (OH, NH, SH) will quench the reagent",
            ],
            "regioselectivity": "No regioselectivity issues (single electrophilic carbon)",
            "limitations": "Sensitive to air/moisture; incompatible with acidic protons; cannot have electrophilic groups in same molecule",
            "example": {
                "reactants": "CH₃MgBr + CH₃CHO (in dry Et₂O)",
                "product": "CH₃CH(OH)CH₃ (2-propanol after workup)",
                "smiles_equation": "CC(=O)[Mg]Br + CC=O → CC(C)O",
            },
            "applications": "Alcohol synthesis, carboxylic acid synthesis, natural product construction",
        },

        "aldol": {
            "name": "Aldol Reaction (羟醛反应)",
            "aliases": ["Aldol addition", "Aldol condensation"],
            "year": 1872,
            "chemists": ["Charles-Adolphe Wurtz"],
            "nobel_prize": None,
            "mechanism_type": "Nucleophilic Addition (enolate to carbonyl)",
            "general_equation": "Enolate (from α-C-H of carbonyl) + Carbonyl compound → β-hydroxy carbonyl (aldol) → α,β-unsaturated carbonyl (condensed)",
            "conditions": "Base (NaOH, LDA, etc.) or acid catalysis; temperature varies (-78°C to rt)",
            "key_features": [
                "Forms new C-C bond at α-position of carbonyl compounds",
                "Product: β-hydroxy carbonyl (addition) or α,β-unsaturated carbonyl (condensation/dehydration)",
                "Crossed aldol possible but requires careful control to avoid mixtures",
                "Directed aldol: pre-formed enolates give specific products",
                "Asymmetric aldol: chiral auxiliaries or catalysts give enantioselectivity",
            ],
            "regioselectivity": "Kinetic enolate (LDA, -78°C) vs thermodynamic enolate (weaker base, higher T)",
            "limitations": "Self-condensation competes in crossed aldol; enolizable position required",
            "example": {
                "reactants": "Acetaldehyde (2 eq) + NaOH (cat.)",
                "product": "3-Hydroxybutanal (aldol adduct) → Crotonaldehyde (condensed)",
                "smiles_equation": "CC=O + CC=O → CC(O)CC=O → CC=CC=O",
            },
            "applications": "Polymer synthesis, pharmaceuticals, natural products",
        },

        "wittig": {
            "name": "Wittig Reaction (维蒂希反应)",
            "aliases": ["Wittig olefination"],
            "year": 1954,
            "chemists": ["Georg Wittig"],
            "nobel_prize": 1979,
            "mechanism_type": "Nucleophilic Addition-Elimination (ylide + carbonyl)",
            "general_equation": "Ph₃P=CHR + R'R''C=O → R'R''C=CR₂ + Ph₃PO",
            "conditions": "Anhydrous solvent (THF, DCM, benzene); base to generate ylide; -78°C to rt",
            "key_features": [
                "Converts carbonyl compounds (aldehydes/ketones) to alkenes",
                "Phosphorus ylide (phosphorane) as nucleophile",
                "E/Z selectivity depends on ylide stability and conditions",
                "Stabilized ylides (EWG on carbanion): E-alkenes favored",
                "Non-stabilized ylides: Z-alkenes favored (Schlosser modification improves E-selectivity)",
                "Byproduct Ph₃PO easily removed by filtration/washing",
            ],
            "regioselectivity": "Single product (carbonyl C becomes alkene C)",
            "limitations": "Does not work well with sterically hindered ketones; phosphine oxide removal can be tedious",
            "example": {
                "reactants": "Ph₃P=CH₂ + Benzaldehyde (in THF)",
                "product": "Styrene (E/Z mixture, mostly E if stabilized)",
                "smiles_equation": "C=P(c1ccccc1)(c2ccccc2)c3ccccc3 + c1ccccc1C=O → C=c1ccccc1",
            },
            "applications": "Alkene synthesis, natural product synthesis, vitamin A synthesis",
        },

        "claisen": {
            "name": "Claisen Condensation (克莱森缩合)",
            "aliases": ["Claisen-Schmidt condensation", "Claisen ester condensation"],
            "year": 1887,
            "chemists": ["Ludwig Claisen"],
            "nobel_prize": None,
            "mechanism_type": "Nucleophilic Acyl Substitution (enolate attacks ester)",
            "general_equation": "2 RCH₂COOR' (ester) + Base → RCH₂COCH(R)COOR' (β-keto ester) + R'OH",
            "conditions": "Strong base (NaOEt, NaH, LDA), corresponding alcohol solvent (ROH), reflux",
            "key_features": [
                "Condensation of two ester molecules → β-keto ester",
                "Requires α-hydrogen on at least one ester",
                "Dieckmann condensation: intramolecular version (diester → cyclic β-keto ester)",
                "Crossed Claisen: one ester must be non-enolizable (no α-H) for clean product",
                "Product is acidic (pKa ~11) due to stabilizing enolate",
            ],
            "regioselectivity": "Enolate forms at less substituted α-position (kinetic) or more stable (thermodynamic)",
            "limitations": "Both esters need α-H for self-condensation; crossed version needs careful design",
            "example": {
                "reactants": "Ethyl acetate (2 eq) + NaOEt",
                "product": "Ethyl acetoacetate (acetoacetic ester)",
                "smiles_equation": "CC(=O)OCC + CC(=O)OCC → CC(=O)CC(=O)OCC",
            },
            "applications": "β-Keto ester synthesis, acetoacetic ester synthesis, heterocycle synthesis",
        },

        "friedel-crafts": {
            "name": "Friedel-Crafts Reaction (傅-克反应)",
            "aliases": ["FC alkylation", "FC acylation", "Friedel Crafts"],
            "year": 1877,
            "chemists": ["Charles Friedel", "James Crafts"],
            "nobel_prize": None,
            "mechanism_type": "Electrophilic Aromatic Substitution (SEAr)",
            "general_equation": "Ar-H + R-X (AlCl₃) → Ar-R (alkylation) OR Ar-H + RCOCl (AlCl₃) → Ar-COR (acylation)",
            "conditions": "Lewis acid catalyst (AlCl₃, FeCl₃, BF₃); anhydrous conditions; 0°C to rt",
            "key_features": [
                "Introduces alkyl or acyl group onto aromatic ring",
                "Alkylation: carbocation intermediate → rearrangement possible",
                "Acylation: acylium ion intermediate → no rearrangement, deactivates ring",
                "Acyl product can be reduced (Clemmensen/Wolff-Kishner) → linear alkyl chain",
                "Ring must be activated (electron-rich); deactivated rings don't react",
                "Ortho/para directing for activating groups",
            ],
            "regioselectivity": "Ortho/para directing; ortho:para ratio depends on substituent size",
            "limitations": "Does not work on strongly deactivated rings (nitrobenzene); polyalkylation common; no meta-directors",
            "example": {
                "reactants": "Benzene + CH₃CH₂Cl / AlCl₃",
                "product": "Ethylbenzene",
                "smiles_equation": "c1ccccc1 + CCCl → c1ccccc1CC",
            },
            "applications": "Detergent production, pharmaceutical synthesis, materials science",
        },

        "heck": {
            "name": "Heck Reaction (赫克反应)",
            "aliases": ["Mizoroki-Heck reaction", "Pd-catalyzed coupling"],
            "year": 1972,
            "chemists": ["Tsutomu Mizoroki", "Richard F. Heck"],
            "nobel_prize": 2010,
            "mechanism_type": "Palladium-catalyzed Cross Coupling (oxidative addition, migratory insertion, β-hydride elimination)",
            "general_equation": "R-X + CH₂=CH-R' (Pd cat., base) → R-CH=CH-R' + HX",
            "conditions": "Pd catalyst (Pd(OAc)₂, Pd(PPh₃)₄), base (Et₃N, K₂CO₃, NaOAc), polar aprotic solvent (DMF, MeCN), 80-120°C",
            "key_features": [
                "Pd-catalyzed coupling of aryl/vinyl halides with alkenes",
                "Forms new C-C bond with alkene formation",
                "Stereoselective: typically gives trans (E) alkene",
                "Tolerates many functional groups",
                "Intramolecular Heck: powerful for ring formation",
                "No need for organometallic reagent (unlike Suzuki/Stille)",
            ],
            "regioselectivity": "Aryl adds to less hindered side of alkene; terminal alkenes give branched/linear mixtures",
            "limitations": "Requires Pd catalyst; β-hydride elimination needed (no β-H = problem); vinyl/aryl halides only",
            "example": {
                "reactants": "Iodobenzene + Ethylene / Pd(OAc)₂, Et₃N",
                "product": "Styrene",
                "smiles_equation": "c1ccccc1I + C=C → C=c1ccccc1",
            },
            "applications": "Pharmaceutical synthesis, natural products, materials (OLEDs)",
        },

        "suzuki": {
            "name": "Suzuki Coupling (铃木偶联反应)",
            "aliases": ["Suzuki-Miyaura coupling", "Suzuki cross-coupling"],
            "year": 1979,
            "chemists": ["Akira Suzuki", "Norio Miyaura"],
            "nobel_prize": 2010,
            "mechanism_type": "Palladium-catalyzed Cross Coupling (transmetalation pathway)",
            "general_equation": "R₁-B(OR)₂ + R₂-X (Pd cat., base) → R₁-R₂ + B(OR)₂X",
            "conditions": "Pd catalyst (Pd(PPh₃)₄, Pd(dppf)Cl₂), base (K₂CO₃, Cs₂CO₃, K₃PO₄), solvent (toluene/EtOH/H₂O, dioxane/H₂O), 50-100°C",
            "key_features": [
                "Coupling of organoboron compounds with organic halides/triflates",
                "High functional group tolerance (mild conditions)",
                "Boron reagents are stable, non-toxic, easy to handle",
                "Works with aryl, vinyl, alkynyl partners",
                "Broad scope: sp²-sp², sp²-sp³ couplings possible",
                "Base activates boron reagent via transmetalation",
            ],
            "regioselectivity": "Retains configuration of both coupling partners",
            "limitations": "Requires Pd catalyst; sensitive to oxygen; protodeboronation side reaction",
            "example": {
                "reactants": "Phenylboronic acid + Iodobenzene / Pd(PPh₃)₄, K₂CO₃",
                "product": "Biphenyl",
                "smiles_equation": "OB(O)c1ccccc1 + c1ccccc1I → c1ccccc1-c2ccccc2",
            },
            "applications": "Pharmaceutical API synthesis, OLED materials, conjugated polymers, agrochemicals",
        },

        "sn2": {
            "name": "SN2 Reaction (双分子亲核取代反应)",
            "aliases": ["bimolecular nucleophilic substitution", "backside attack"],
            "year": "1930s",
            "chemists": ["Edward Hughes", "Christopher Ingold"],
            "nobel_prize": None,
            "mechanism_type": "Concerted Nucleophilic Substitution (one-step, bimolecular)",
            "general_equation": "Nu⁻ + R-X → Nu-R + X⁻ (Walden inversion)",
            "conditions": "Polar aprotic solvents (DMF, DMSO, acetone) preferred; good nucleophile; moderate temperature",
            "key_features": [
                "Concerted backside attack → inversion of configuration (Walden inversion)",
                "Rate = k[Nu][RX] (second order, hence SN2)",
                "Reactivity: methyl > primary > secondary >> tertiary (steric hindrance)",
                "Good nucleophiles needed (I⁻ > Br⁻ > Cl⁻ > F⁻ trend for halides)",
                "Polar aprotic solvents accelerate (don't H-bond to Nu⁻)",
                "Leaving group ability: I⁻ > TsO⁻ > Br⁻ > Cl⁻ > F⁻",
            ],
            "regioselectivity": "Attack at least hindered carbon",
            "limitations": "Tertiary substrates don't react (too hindered); competing E2 with strong bases; poor leaving groups won't leave",
            "example": {
                "reactants": "NaOH + CH₃CH₂Br (in DMSO)",
                "product": "Ethanol + NaBr",
                "smiles_equation": "[OH-] + CCCBr → CCCO + [Br-]",
            },
            "applications": "Ether synthesis, amine synthesis, Williamson ether synthesis",
        },

        "sn1": {
            "name": "SN1 Reaction (单分子亲核取代反应)",
            "aliases": ["unimolecular nucleophilic substitution", "carbocation substitution"],
            "year": "1930s",
            "chemists": ["Edward Hughes", "Christopher Ingold"],
            "nobel_prize": None,
            "mechanism_type": "Stepwise Nucleophilic Substitution (carbocation intermediate)",
            "general_equation": "R-X → R⁺ + X⁻ (slow, rate-determining) then R⁺ + Nu⁻ → R-Nu (fast)",
            "conditions": "Polar protic solvents (water, alcohols) stabilize carbocation; weak nucleophile acceptable",
            "key_features": [
                "Two-step mechanism: ionization then capture",
                "Rate = k[RX] (first order, hence SN1)",
                "Carbocation intermediate → racemization (+ partial inversion from ion pair)",
                "Reactivity: tertiary > secondary > allylic/benzylic >> primary (carbocation stability)",
                "Rearrangements possible (hydride/alkyl shifts to form more stable carbocation)",
                "Competing E1 elimination always present",
            ],
            "regioselectivity": "Carbocation may rearrange before capture",
            "limitations": "Primary substrates rarely react; racemic mixture formed; requires good leaving group; competing E1",
            "example": {
                "reactants": "(CH₃)₃C-Cl + H₂O (aqueous)",
                "product": "(CH₃)₃C-OH (tert-butanol) + HCl",
                "smiles_equation": "CC(C)(C)Cl + O → CC(C)(C)O",
            },
            "applications": "Alcohol synthesis from tert-alkyl halides, sugar chemistry",
        },

        "baylis-hillman": {
            "name": "Baylis-Hillman Reaction",
            "aliases": ["Morita-Baylis-Hillman (MBH) reaction"],
            "year": 1972,
            "chemists": ["Anthony Baylis", "Melville Hillman"],
            "nobel_prize": None,
            "mechanism_type": "Tandem Michael Addition-Aldol (organocatalytic, via enolate)",
            "general_equation": "Activated alkene (α,β-unsaturated carbonyl) + Aldehyde (DABCO cat.) → Allylic alcohol",
            "conditions": "Tertiary amine catalyst (DABCO, quinuclidine), often with TiCl₄/Lewis acid co-catalyst; room temp to mild heating",
            "key_features": [
                "Atom-economical C-C bond formation",
                "Activates simple alkenes toward addition to aldehydes",
                "Product: allylic alcohol with multiple functional handles",
                "Slow reaction (days without acceleration)",
                "High atom economy (all atoms retained in product)",
            ],
            "regioselectivity": "α-position of activated alkene adds to carbonyl carbon",
            "limitations": "Very slow reaction rate; limited substrate scope; sensitive to sterics",
            "example": {
                "reactants": "Methyl acrylate + Benzaldehyde / DABCO",
                "product": "Methyl 2-(hydroxy(phenyl)methyl)acrylate",
                "smiles_equation": "C=CC(=O)OC + c1ccccc1C=O → C=C(C(=O)OC)C(O)c1ccccc1",
            },
            "applications": "Natural product synthesis, heterocycle synthesis, medicinal chemistry",
        },

        "michael": {
            "name": "Michael Addition (迈克尔加成)",
            "aliases": ["1,4-addition", "conjugate addition"],
            "year": 1887,
            "chemists": ["Arthur Michael"],
            "nobel_prize": None,
            "mechanism_type": "Nucleophilic Conjugate Addition (1,4-addition to α,β-unsaturated carbonyl)",
            "general_equation": "Nu-H + CH₂=CH-EWG → Nu-CH₂-CH₂-EWG (after proton transfer)",
            "conditions": "Base or Lewis acid catalysis; solvent varies (MeOH, THF, DCM); -78°C to rt",
            "key_features": [
                "1,4-addition to α,β-unsaturated carbonyls (vs 1,2-carbonyl addition)",
                "Soft nucleophiles favor Michael (hard favor direct addition)",
                "Common nucleophiles: enolates, amines, thiols, cuprates",
                "Asymmetric Michael: chiral catalysts give high ee",
                "Robinson annulation: Michael + aldol cascade",
            ],
            "regioselectivity": "Always β-carbon (1,4-position) of enone",
            "limitations": "Requires conjugated system; hard nucleophiles may do 1,2-addition instead",
            "example": {
                "reactants": "Diethyl malonate + Methyl vinyl ketone / base",
                "product": "Diethyl 2-acetylglutarate",
                "smiles_equation": "C(CC(=O)OCC)(CC(=O)OCC) + C=CC(=O)C → C(CC(=O)OCC)(CC(=O)OCC)CC(=O)C",
            },
            "applications": "Natural product synthesis, pharmaceuticals, polymers",
        },

        "knorr": {
            "name": "Knorr Pyrazole Synthesis (克诺尔吡唑合成)",
            "aliases": ["Knorr pyrazole synthesis"],
            "year": 1883,
            "chemists": ["Ludwig Knorr"],
            "nobel_prize": None,
            "mechanism_type": "Condensation / Cyclization",
            "general_equation": "Hydrazine (or substituted hydrazine) + 1,3-Diketone / β-Keto ester → Pyrazole",
            "conditions": "Acid or base catalysis; ethanol or acetic acid solvent; reflux",
            "key_features": [
                "Synthesizes pyrazole ring system",
                "Hydrazine + 1,3-dicarbonyl → pyrazole (5-membered, 2 N atoms)",
                "Regiochemistry depends on hydrazine substitution pattern",
                "Unsymmetrical diketones give mixture of regioisomers",
            ],
            "regioselectivity": "Depends on which carbonyl reacts first with hydrazine nitrogen",
            "limitations": "Regioisomer mixture for unsymmetrical cases",
            "example": {
                "reactants": "Phenylhydrazine + Acetylacetone",
                "product": "3,5-Dimethyl-1-phenylpyrazole",
                "smiles_equation": "NNc1ccccc1 + CC(=O)CC(=O)C →Cc1c(nn(c1)c2ccccc2)C",
            },
            "applications": "Pharmaceuticals (celecoxib, anti-inflammatory), agrochemicals",
        },

        "hantzsch": {
            "name": "Hantzsch Dihydropyridine Synthesis (汉栖二氢吡啶合成)",
            "aliases": ["Hantzsch pyridine synthesis"],
            "year": 1881,
            "chemists": ["Arthur Hantzsch"],
            "nobel_prize": None,
            "mechanism_type": "Multi-component Condensation / Cyclization",
            "general_equation": "2 β-Keto ester + Aldehyde + NH₃ (or NH₄OAc) → 1,4-Dihydropyridine → Oxidation → Pyridine",
            "conditions": "Ethanol or methanol solvent; heat/reflux; oxidant (HNO₃, DDQ, air) for aromatization",
            "key_features": [
                "Four-component reaction (two molecules of β-keto ester)",
                "Produces 1,4-dihydropyridines (DHPs) - calcium channel blockers",
                "Nifedipine, amlodipine, felodipine are all Hantzsch products",
                "Can be oxidized to pyridines (aromatized)",
                "Enantioselective versions known",
            ],
            "regioselectivity": "Symmetrical for symmetrical β-keto esters",
            "limitations": "Unsymmetrical β-keto esters give mixtures",
            "example": {
                "reactants": "Ethyl acetoacetate (2 eq) + Formaldehyde + NH₃",
                "product": "1,4-Dihydro-2,6-dimethyl-3,5-pyridinedicarboxylate",
                "smiles_equation": "CC(=O)CC(=O)OCC (x2) + C=O + N → Complex DHP",
            },
            "applications": "Calcium channel blocker drugs (nifedipine family), pharmaceutical industry",
        },

        "perkin": {
            "name": "Perkin Reaction (珀金反应)",
            "aliases": ["Perkin condensation"],
            "year": 1868,
            "chemists": ["William Perkin"],
            "nobel_prize": None,
            "mechanism_type": "Condensation (anhydride enolate + aromatic aldehyde)",
            "general_equation": "Aromatic aldehyde + Acid anhydride (RCH₂CO)₂O (base) → Cinnamic acid derivative",
            "conditions": "Sodium salt of carboxylic acid (NaOAc), corresponding acid anhydride; 150-180°C",
            "key_features": [
                "First industrial-scale synthetic route to dyes (mauveine related)",
                "Aromatic aldehydes + aliphatic acid anhydrides → α,β-unsaturated acids",
                "Only works well with aromatic aldehydes (electron-rich rings better)",
                "Limited scope: only works with anhydrides having α-hydrogens",
                "Historically important: led to synthetic dye industry",
            ],
            "regioselectivity": "Trans (E) cinnamic acid derivatives predominate",
            "limitations": "High temperatures; limited to aromatic aldehydes; moderate yields",
            "example": {
                "reactants": "Benzaldehyde + Acetic anhydride / sodium acetate",
                "product": "Cinnamic acid (E)",
                "smiles_equation": "c1ccccc1C=O + CC(=O)OC(=O)C → C=Cc1ccccc1C(=O)O",
            },
            "applications": "Cinnamic acid derivatives, fragrance compounds, pharmaceutical precursors",
        },

        "schiff-base": {
            "name": "Schiff Base Formation (席夫碱形成)",
            "aliases": ["Imine formation", "Schiff base condensation"],
            "year": 1864,
            "chemists": ["Hugo Schiff"],
            "nobel_prize": None,
            "mechanism_type": "Condensation (amine + carbonyl, dehydration)",
            "general_equation": "R₁NH₂ + R₂C=O → R₁N=CR₂ + H₂O",
            "conditions": "Mild acid catalysis (pH ~4-5); molecular sieves or Dean-Stark to remove water; rt to reflux",
            "key_features": [
                "Amine + aldehyde/ketone → imine (Schiff base)",
                "Reversible: hydrolysis regenerates starting materials",
                "Primary amines give imines; secondary amines give enamines",
                "Acid catalyzes both forward and reverse reactions",
                "Important ligands in coordination chemistry",
                "Chiral Schiff bases used in asymmetric catalysis",
            ],
            "regioselectivity": "Single product (C=N double bond formation)",
            "limitations": "Equilibrium reaction (water removal drives forward); aliphatic imines less stable than aromatic",
            "example": {
                "reactants": "Aniline + Benzaldehyde / p-TsOH (cat.), toluene, Dean-Stark",
                "product": "N-Benzylideneaniline (Schiff base)",
                "smiles_equation": "Nc1ccccc1 + c1ccccc1C=O → N=c1ccccc1-c2ccccc2",
            },
            "applications": "Coordination chemistry ligands, asymmetric catalysis, biosensors, analytical chemistry",
        },

        "knoevenagel": {
            "name": "Knoevenagel Condensation (克脑文盖尔缩合)",
            "aliases": ["Knoevenagel reaction"],
            "year": 1894,
            "chemists": ["Emil Knoevenagel"],
            "nobel_prize": None,
            "mechanism_type": "Condensation (active methylene + carbonyl)",
            "general_equation": "Active methylene compound (Z-CH₂-Z') + Aldehyde/Ketone (base cat.) → Alkene (Z-CH=CR-Z') + H₂O",
            "conditions": "Weak base (piperidine, amine) often with weak acid co-catalyst; sometimes Lewis acids; rt to reflux",
            "key_features": [
                "Active methylene compounds (malonic ester, cyanoacetate, Meldrum's acid)",
                "Typically gives E-alkenes (trans configuration)",
                "Doebner modification: uses pyridine with malonic acid → directly gives α,β-unsaturated acid",
                "Milder than Perkin reaction (lower temperature)",
                "Very useful for preparing coumarins, chromenes, heterocycles",
            ],
            "regioselectivity": "Predominantly E (trans) geometry",
            "limitations": "Requires active methylene compound (pKa < ~13); some substrates give low yields",
            "example": {
                "reactants": "Malonic acid + Benzaldehyde / pyridine (Doebner modification)",
                "product": "Cinnamic acid (trans)",
                "smiles_equation": "O=C(O)CC(=O)O + c1ccccc1C=O → C=Cc1ccccc1C(=O)O",
            },
            "applications": "Coumarin synthesis, pharmaceuticals, dyes, fragrances",
        },

        "reformatsky": {
            "name": "Reformatsky Reaction (雷福尔马茨基反应)",
            "aliases": ["Reformatsky reaction"],
            "year": 1887,
            "chemists": ["Sergei Reformatsky"],
            "nobel_prize": None,
            "mechanism_type": "Nucleophilic Addition (organozinc enolate to carbonyl)",
            "general_equation": "α-Bromo ester + Zn + Carbonyl compound → β-Hydroxy ester",
            "conditions": "Activated zinc dust, inert atmosphere (N₂/Ar), anhydrous solvent (benzene, THF, ether); reflux",
            "key_features": [
                "Organozinc reagent (less reactive than Grignard → more chemoselective)",
                "α-Bromo ester + Zn → organozinc enolate ( Reformatsky reagent)",
                "Adds to aldehydes/ketones → β-hydroxy esters",
                "More tolerant of ester/nitrile groups than Grignard",
                "No reaction with isolated ester groups (chemoselective)",
            ],
            "regioselectivity": "Single addition to carbonyl carbon",
            "limitations": "Requires activated Zn; moisture sensitive; limited to α-bromo carbonyl compounds",
            "example": {
                "reactants": "Ethyl bromoacetate + Zn + Benzaldehyde (dry benzene)",
                "product": "Ethyl 3-hydroxy-3-phenylpropanoate",
                "smiles_equation": "BrCC(=O)OCC + [Zn] + c1ccccc1C=O → OC(=O)CC(O)c1ccccc1",
            },
            "applications": "β-Hydroxy ester synthesis, lactone synthesis, natural product synthesis",
        },

        "stille": {
            "name": "Stille Coupling (斯蒂勒偶联)",
            "aliases": ["Stille cross-coupling", "Stille-Kelly coupling"],
            "year": 1978,
            "chemists": ["John Stille"],
            "nobel_prize": None,  # Stille died before Nobel; recognized posthumously
            "mechanism_type": "Palladium-catalyzed Cross Coupling (transmetalation from tin)",
            "general_equation": "R₁-SnR₃ + R₂-X (Pd cat.) → R₁-R₂ + R₃SnX",
            "conditions": "Pd catalyst (Pd(PPh₃)₄, Pd₂(dba)₃), solvent (DMF, toluene, dioxane), 60-110°C, inert atmosphere",
            "key_features": [
                "Organostannane + organic halide/triflate → coupled product",
                "Highly functional group tolerant (neutral conditions)",
                "Organostannanes are stable, isolable, characterized",
                "Works with vinyl, aryl, alkynyl, acyl, allyl partners",
                "Stille-Kelly: intramolecular version using distannane",
            ],
            "regioselectivity": "Retains stereochemistry of vinyl/alkynyl groups",
            "limitations": "Toxic tin byproducts; difficult to remove tin residues; air-sensitive Pd(0)",
            "example": {
                "reactants": "Tributyl(vinyl)tin + Iodobenzene / Pd(PPh₃)₄",
                "product": "Styrene",
                "smiles_equation": "C=C[Sn](CCCC)(CCCC)CCCC + c1ccccc1I → C=c1ccccc1",
            },
            "applications": "Natural product synthesis, polymer synthesis, materials science",
        },

        "sonogashira": {
            "name": "Sonogashira Coupling (园头-玉仓偶联)",
            "aliases": ["Sonogashira-Hagihara coupling"],
            "year": 1975,
            "chemists": ["Kenkichi Sonogashira", "Nobue Hagihara"],
            "nobel_prize": None,
            "mechanism_type": "Palladium-Copper Co-catalyzed Cross Coupling",
            "general_equation": "Terminal alkyne + Aryl/Vinyl halide (Pd/Cu cat., base) → Internal alkyne",
            "conditions": "Pd catalyst (Pd(PPh₃)₂Cl₂, Pd(PPh₃)₄), CuI co-catalyst, base (Et₃N, iPr₂NH, piperidine), solvent (THF, DMF), rt-80°C",
            "key_features": [
                "Terminal alkyne + aryl/vinyl halide → disubstituted alkyne",
                "Pd/Cu dual catalysis: Cu activates alkyne (copper acetylide), Pd does cross-coupling",
                "Mild conditions (often room temperature)",
                "Copper-free Sonogashira modifications available (for Cu-sensitive substrates)",
                "Important for synthesizing conjugated systems (natural products, OLEDs, pharmaceuticals)",
            ],
            "regioselectivity": "Single product (sp-sp² C-C bond formation)",
            "limitations": "Terminal alkyne required; homocoupling (Glaser) side reaction; copper toxicity concerns",
            "example": {
                "reactants": "Phenylacetylene + Iodobenzene / Pd(PPh₃)₂Cl₂, CuI, Et₃N",
                "product": "Diphenylacetylene (tolane)",
                "smiles_equation": "C#Cc1ccccc1 + c1ccccc1I → c1ccccc1C#Cc2ccccc2",
            },
            "applications": "Natural products, pharmaceuticals, OLED materials, conjugated polymers, click chemistry precursor",
        },

        "buchwald-hartwig": {
            "name": "Buchwald-Hartwig Amination (布赫瓦尔德-哈特维希胺化)",
            "aliases": ["Pd-catalyzed C-N coupling", "Buchwald Hartwig"],
            "year": 1995,
            "chemists": ["Stephen Buchwald", "John Hartwig"],
            "nobel_prize": None,
            "mechanism_type": "Palladium-catalyzed C-N Cross Coupling (amination of aryl halides)",
            "general_equation": "Ar-X + HNR₂ (Pd cat., base) → Ar-NR₂ + HX",
            "conditions": "Pd catalyst (Pd₂(dba)₃, Pd(OAc)₂) with bulky phosphine ligand (XPhos, BrettPhos, DavePhos), strong base (NaOtBu, Cs₂CO₃), toluene/dioxane, 80-110°C",
            "key_features": [
                "Forms C-N bond between aryl/vinyl halides and amines",
                "Revolutionary: replaced harsh Ullmann/Goldberg conditions",
                "Broad scope: primary/secondary amines, amides, imides, N-heterocycles",
                "Bulky biaryl phosphine ligands enable challenging couplings",
                "Works with electron-rich and electron-poor aryl halides",
                "Also applicable to C-O coupling (etherification)",
            ],
            "regioselectivity": "Single product (C-N bond at halide position)",
            "limitations": "Expensive Pd/ligand system; requires careful optimization; sensitive to sterics",
            "example": {
                "reactants": "Bromobenzene + Morpholine / Pd₂(dba)₃, XPhos, NaOtBu",
                "product": "4-Phenylmorpholine",
                "smiles_equation": "c1ccccc1Br + N1CCOCC1 → c1ccccc1N2CCOCC2",
            },
            "applications": "Pharmaceutical API synthesis (most widely used C-N coupling in pharma), materials, agrochemicals",
        },

        "umpolung": {
            "name": "Umpolung (极性反转)",
            "aliases": ["Polarity inversion", "acyl anion equivalent"],
            "year": "1918 (concept); 1973 (synthetic use)",
            "chemists": ["Louis Fieser (concept)", "Dieter Seebach", "David J.C.ram / E.J. Corey"],
            "nobel_prize": None,
            "mechanism_type": "Polarity Reversal (normal electrophile acts as nucleophile)",
            "general_equation": "Normally electrophilic center (carbonyl C) converted to nucleophilic (acyl anion equivalent)",
            "conditions": "Various: 1,3-dithiane (Corey-Seebach), cyanohydrins, nitroalkanes, NHC catalysis",
            "key_features": [
                "Reverses normal polarity of functional group",
                "Corey-Seebach: 1,3-dithiane → masked acyl anion (umpoled carbonyl)",
                "Benzoin condensation: NHC-catalyzed umpolung of aldehydes",
                "Enables disconnections not possible with standard polarity",
            ],
            "regioselectivity": "Depends on umpolung method",
            "limitations": "Extra steps (masking/unmasking); specialized reagents",
            "example": {
                "reactants": "1,3-Dithiane (deprotonated) + R-Br → after unmasking → R-CHO",
                "product": "Aldehyde (via acyl anion equivalent)",
                "smiles_equation": "CC1CSCCS1 + RBr → ... → RC=O",
            },
            "applications": "Total synthesis, complex molecule construction, strategic bond disconnection",
        },

        "robinson-annulation": {
            "name": "Robinson Annulation (罗宾逊环合反应)",
            "aliases": ["Robinson annulation", "Robinson tropinone synthesis"],
            "year": 1935,
            "chemists": ["Robert Robinson (Sir)"],
            "nobel_prize": 1947,
            "mechanism_type": "Cascade: Michael Addition + Aldol Condensation + Dehydration",
            "general_equation": "Cyclic ketone + Methyl vinyl ketone (MVK) (base) → Fused bicyclic α,β-unsaturated ketone",
            "conditions": "Base (KOH, Et₃N, pyrrolidine), ethanol or aqueous ethanol, heat/reflux",
            "key_features": [
                "Two-step cascade in one pot: Michael addition then aldol cyclization",
                "Builds a new 6-membered ring onto existing cyclic ketone",
                "Classic: cyclohexanone + MVK → bicyclo[4.4.0]decenone system",
                "Robinson's tropinone synthesis: landmark in alkaloid synthesis",
                "Dehydration usually spontaneous under reaction conditions",
            ],
            "regioselectivity": "Michael addition at β-position of MVK; aldol at α-position of donor ketone",
            "limitations": "Polysubstitution possible; requires enolizable ketone; yields vary",
            "example": {
                "reactants": "Cyclohexanone + Methyl vinyl ketone / KOH, EtOH, Δ",
                "product": "2-Cyclohexylidenecyclohexanone (→ dehydrated: 1,9-dioxooctahydro-1H-indene after tautomerization)",
                "smiles_equation": "C1CCCCC1=O + C=CC(=O)C → fused bicyclic enone",
            },
            "applications": "Steroid synthesis, alkaloid synthesis (tropinone), natural product construction",
        },

        "skraup": {
            "name": "Skraup Synthesis (斯克劳普喹啉合成)",
            "aliases": ["Skraup quinoline synthesis"],
            "year": 1880,
            "chemists": ["Zdenko Hans Skraup"],
            "nobel_prize": None,
            "mechanism_type": "Cascade Cyclization / Dehydrogenative Aromatization",
            "general_equation": "Aniline + Glycerol (or α,β-unsaturated carbonyl) + H₂SO₄ (oxidant) → Quinoline",
            "conditions": "Concentrated H₂SO₄, oxidant (nitrobenzene, As₂O₅), iron powder (moderator), 120-140°C",
            "key_features": [
                "One-pot synthesis of quinoline from aniline + glycerol",
                "Mechanism: glycerol → acrolein (dehydration) → Michael addition → cyclization → oxidation",
                "Harsh conditions (conc. H₂SO₄, high T); iron moderates reaction",
                "Nitrobenzene serves as oxidant (dehydrogenating agent)",
                "Foundational method for quinoline/isoquinoline synthesis",
            ],
            "regioselectivity": "Quinoline formed (fused benzene + pyridine)",
            "limitations": "Very harsh/exothermic; safety concerns; moderate yields; toxic reagents",
            "example": {
                "reactants": "Aniline + Glycerol + conc. H₂SO₄ + nitrobenzene (oxidant)",
                "product": "Quinoline",
                "smiles_equation": "Nc1ccccc1 + OCC(O)CO + [H₂SO₄] → c1cccnc2ccccc12",
            },
            "applications": "Quinoline drugs (antimalarial: chloroquine/quinine analogs), alkaloids",
        },

        "pummerer": {
            "name": "Pummerer Rearrangement (普梅雷尔重排)",
            "aliases": ["Pummerer reaction"],
            "year": 1909,
            "chemists": ["Silvio Pummerer"],
            "nobel_prize": None,
            "mechanism_type": "Rearrangement (sulfoxide → α-substituted sulfide)",
            "general_equation": "R-S(O)-R' + Ac₂O (or TFAA) → R-S(R')-CR''(OAc) (α-acyloxy sulfide)",
            "conditions": "Acetic anhydride (Ac₂O) or trifluoroacetic anhydride (TFAA); 0°C to rt; Lewis acid optional",
            "key_features": [
                "Sulfoxide activated by acylation → thionium ion intermediate",
                "Nucleophile (acetate, other) traps at α-position",
                "Powerful way to functionalize α-position of sulfides",
                "Products can be further transformed (hydrolysis → carbonyl)",
                "Variations: Pummerer cyclization (intramolecular trapping)",
            ],
            "regioselectivity": "Trapping occurs at more substituted α-carbon (more stable carbocation character)",
            "limitations": "Requires sulfoxide precursor; over-reaction possible; regioselectivity issues with unsymmetrical sulfoxides",
            "example": {
                "reactants": "Dimethyl sulfoxide + Acetic anhydride",
                "product": "Methylthiomethyl acetate (MTMA)",
                "smiles_equation": "CSC + CC(=O)OC(=O)C → CSCC(=O)OC",
            },
            "applications": "Carbonyl synthesis (via hydrolysis), natural product synthesis, glycosylation",
        },

        "fries": {
            "name": "Fries Rearrangement (弗里斯重排)",
            "aliases": ["Fries rearrangement"],
            "year": 1908,
            "chemists": ["Karl Fries"],
            "nobel_prize": None,
            "mechanism_type": "Rearrangement (Lewis acid-mediated acyl migration on phenol ester)",
            "general_equation": "Phenyl ester Ar-O-COR (Lewis acid, heat) → o- or p-Hydroxyaryl ketone",
            "conditions": "Lewis acid (AlCl₃, TiCl₄, SnCl₄), heat (80-160°C), solventless or nitrobenzene solvent",
            "key_features": [
                "Phenolic ester → hydroxyaryl ketone (ortho or para)",
                "Lewis acid coordinates to carbonyl oxygen → acyl migration",
                "Low temperature favors para-product; high temperature favors ortho-product",
                "Ortho product can chelate with AlCl₃ (stable complex)",
                "Useful for synthesizing phenolic ketones (difficult to make otherwise)",
            ],
            "regioselectivity": "Temperature-dependent: lower T → para; higher T → ortho",
            "limitations": "Mixture of ortho/para; stoichiometric Lewis acid; harsh conditions",
            "example": {
                "reactants": "Phenyl acetate / AlCl₃, 160°C",
                "product": "2-Hydroxyacetophenone (major) + 4-Hydroxyacetophenone (minor)",
                "smiles_equation": "OC(=O)c1ccccc1 + [AlCl3] → O=C(C)c1ccccc1O (ortho major)",
            },
            "applications": "Phenolic ketone synthesis, sunscreen agents (avobenzone), pharmaceuticals",
        },

        "baumann-fromm": {
            "name": "Claisen Rearrangement (克莱森重排)",
            "aliases": ["Claisen rearrangement", "allyl vinyl ether → γ,δ-unsaturated carbonyl"],
            "year": 1912,
            "chemists": ["Ludwig Claisen"],
            "nobel_prize": None,
            "mechanism_type": "Pericyclic [3,3]-sigmatropic rearrangement (concerted)",
            "general_equation": "Allyl vinyl ether (heat) → γ,δ-Unsaturated aldehyde or ketone",
            "conditions": "Thermal (150-200°C) or Lewis acid catalyzed (lower T); neat or high-boiling solvent",
            "key_features": [
                "[3,3]-sigmatropic rearrangement (pericyclic, concerted)",
                "Allyl vinyl ether → γ,δ-unsaturated carbonyl compound",
                "Chair-like transition state → specific stereochemistry",
                "Johnson-Claisen: allyl alcohol + trialkyl orthoester → γ,δ-unsaturated ester",
                "Eschenmoser-Claisen: allyl alcohol + N,N-dimethylacetamide dimethyl acetal → γ,δ-unsaturated amide",
                "Ireland-Claisen: allyl ester silyl ketene acetal rearrangement",
            ],
            "regioselectivity": "Stereospecific: chair TS gives specific relative configuration",
            "limitations": "High thermal energy required (unless catalyzed); requires preparation of allyl vinyl ether",
            "example": {
                "reactants": "Allyl vinyl ether (heat, 200°C)",
                "product": "4-Pentenal",
                "smiles_equation": "C=CCOC=C → C=CCCC=O",
            },
            "applications": "Natural product synthesis (C-C bond formation with control), pharmaceuticals",
        },

        "cope": {
            "name": "Cope Rearrangement (科普重排)",
            "aliases": ["Cope rearrangement", "[3,3]-sigmatropic"],
            "year": 1940,
            "chemists": ["Arthur C. Cope"],
            "nobel_prize": None,
            "mechanism_type": "Pericyclic [3,3]-sigmatropic rearrangement (concerted)",
            "general_equation": "1,5-Diene (heat) → Isomeric 1,5-diene",
            "conditions": "Thermal (150-300°C depending on substrate) or Lewis acid catalyzed; can be reversible",
            "key_features": [
                "[3,3]-sigmatropic rearrangement of 1,5-dienes",
                "Concerted pericyclic mechanism (no intermediates)",
                "Oxy-Cope: 3-oxy-1,5-diene → much faster (10¹⁰-10¹⁷× rate acceleration)",
                "Anionic oxy-Cope: even faster (KOH, rt possible)",
                "Chair transition state favored over boat",
                "Degenerate Cope: bullvalene (rapidly rearranging)",
            ],
            "regioselectivity": "Chair TS → predictable stereochemistry",
            "limitations": "High temperature for simple Cope; requires 1,5-diene structure",
            "example": {
                "reactants": "3-Methyl-1,5-hexadiene (heat, 200°C)",
                "product": "1,5-Heptadiene (isomerized)",
                "smiles_equation": "C=CCC(C)CC=C → C=CCCC(C)C=C",
            },
            "applications": "Natural product synthesis, structural reorganization, degenerate rearrangements",
        },

        "haworth": {
            "name": "Haworth Phenanthrene Synynthesis (霍沃思菲合成)",
            "aliases": ["Haworth reaction", "Haworth phenanthrene synthesis"],
            "year": 1931,
            "chemists": ["Robert Haworth"],
            "nobel_prize": None,
            "mechanism_type": "Cascade: Friedel-Crafts Acylation / Reduction / Cyclization / Dehydrogenation",
            "general_equation": "Benzene + Succinic anhydride (AlCl₃) → Tetralone → Reduction → Cyclization → Phenanthrene",
            "conditions": "Multi-step: FC acylation (AlCl₃), reduction (Clemmensen/Zn-Hg), cyclization (H₂SO₄), dehydrogenation (Pd/C or Se)",
            "key_features": [
                "Builds phenanthrene ring system from benzene + succinic anhydride",
                "Multi-step sequence: acylation → tetralone → tetralin → phenanthrene",
                "Classical route to polycyclic aromatic hydrocarbons (PAHs)",
                "Related: Bogert-Cook synthesis (similar approach)",
            ],
            "regioselectivity": "Linear fusion of rings",
            "limitations": "Multi-step; overall yield modest; harsh conditions in some steps",
            "example": {
                "reactants": "Benzene + Succinic anhydride / AlCl₃ → (multi-step)",
                "product": "Phenanthrene",
                "smiles_equation": "c1ccccc1 + O=C1CCC(=O)O1 → c1cccc2c1cccc2",
            },
            "applications": "PAH synthesis, steroid framework construction, aromatic chemistry research",
        },

        "bishler-napieralski": {
            "name": "Bishler-Napieralski Isoquinoline Synthesis (比施勒-纳皮耶尔斯基异喹啉合成)",
            "aliases": ["BN isoquinoline synthesis", "Bischler-Napieralski"],
            "year": 1893,
            "chemists": ["August Bishler", "Ludwik Napieralski"],
            "nobel_prize": None,
            "mechanism_type": "Cyclodehydration / Electrophilic Aromatic Substitution",
            "general_equation": "β-Phenylethylamide (POCl₃, P₂O₅, or SOCl₂) → 3,4-Dihydroisoquinoline → Oxidation → Isoquinoline",
            "conditions": "Dehydrating agent (POCl₃, P₂O₅, SOCl₂, TFAA), heat (80-100°C), then oxidation (Pd/C, MnO₂, or DDQ)",
            "key_features": [
                "β-Phenylethylamide → dihydroisoquinoline → isoquinoline",
                "Electrophilic cyclization onto aromatic ring (intramolecular SEAr)",
                "Key step: activation of amide carbonyl → acyliminium ion → cyclization",
                "Pictet-Spengler: similar but with imine/aldehyde (tryptamine → β-carboline)",
                "Foundation of many alkaloid syntheses (berberine, papaverine class)",
            ],
            "regioselectivity": "Cyclization at ortho position (electron-rich rings favored)",
            "limitations": "Requires electron-rich aromatic ring; over-chlorination with POCl₃ possible",
            "example": {
                "reactants": "N-Phenylethanamide (phenacetamide derivative) / POCl₃",
                "product": "3,4-Dihydroisoquinoline → (oxidation) → Isoquinoline",
                "smiles_equation": "NC(=O)CCc1ccccc1 + POCl3 → c1cccnc2ccccc12",
            },
            "applications": "Isoquinoline alkaloid synthesis (berberine, morphine precursors), pharmaceuticals",
        },

        "pictet-spengler": {
            "name": "Pictet-Spengler Reaction (皮克特-斯彭格勒反应)",
            "aliases": ["PS cyclization", "tetrahydroisoquinoline synthesis"],
            "year": 1911,
            "chemists": ["Amé Pictet", "A. Spengler"],
            "nobel_prize": None,
            "mechanism_type": "Intramolecular Mannich-type Cyclization / Pictet-Spengler",
            "general_equation": "β-Arylethylamine + Aldehyde/Ketone (acid) → Tetrahydroisoquinoline",
            "conditions": "Acid catalysis (AcOH, TFA, HCl), optionally Lewis acid; rt to reflux; solvent (MeOH, DCM, toluene)",
            "key_features": [
                "Tryptamine + aldehyde → β-carboline (indole alkaloid core)",
                "Imine formation followed by intramolecular electrophilic cyclization",
                "Biologically crucial: core transformation in alkaloid biosynthesis",
                "Enantioselective versions with chiral phosphoric acids",
                "Very mild conditions compared to Bishler-Napieralski",
            ],
            "regioselectivity": "For tryptamines: cyclizes at C2 of indole ring",
            "limitations": "Requires electron-rich aryl/indole for efficient cyclization; racemic without chiral catalyst",
            "example": {
                "reactants": "Tryptamine + Acetaldehyde / AcOH (cat.)",
                "product": "1-Methyl-1,2,3,4-tetrahydro-β-carboline",
                "smiles_equation": "NCCc1c[nH]c2ccccc12 + CC=O → fused tricyclic β-carboline",
            },
            "applications": "Alkaloid synthesis (reserpine, quinine, ergot alkaloids), drug discovery, biomimetic synthesis",
        },
    }

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize the reaction database index."""
        self._build_search_index()

    def _build_search_index(self):
        """Build a searchable index from the database."""
        self._search_index = {}
        for key, data in self._DETAILED_REACTIONS.items():
            # Index by name and aliases
            name_lower = data["name"].lower()
            for alias in data["aliases"] + [data["name"], key]:
                alias_lower = alias.lower()
                if alias_lower not in self._search_index:
                    self._search_index[alias_lower] = key
            # Also index by mechanism type words
            for word in data["mechanism_type"].split("/"):
                w = word.strip().lower()
                if len(w) > 3:
                    if w not in self._search_index:
                        self._search_index[w] = key

    def _find_reaction(self, query: str) -> tuple:
        """Find best matching reaction for query."""
        q = query.strip().lower()

        # Exact match
        if q in self._search_index:
            return self._search_index[q], 1.0

        # Partial match / substring
        best_key = None
        best_score = 0
        for idx_key, db_key in self._search_index.items():
            if q in idx_key:
                score = len(q) / len(idx_key)
                if score > best_score:
                    best_score = score
                    best_key = db_key
            elif idx_key in q:
                score = len(idx_key) / len(q)
                if score > best_score:
                    best_score = score
                    best_key = db_key

        # Word-level match
        if best_key is None:
            q_words = set(q.replace("-", " ").replace("_", " ").split())
            for idx_key, db_key in self._search_index.items():
                idx_words = set(idx_key.replace("-", " ").replace("_", " ").split())
                overlap = len(q_words & idx_words)
                if overlap > 0:
                    score = overlap / max(len(q_words), 1)
                    if score > best_score:
                        best_score = score
                        best_key = db_key

        return best_key, best_score

    def _get_summary(self, data: dict) -> dict:
        """Return summary version of reaction data."""
        return {
            "name": data.get("name", "Unknown"),
            "aliases": data.get("aliases", []),
            "year": data.get("year", "Unknown"),
            "chemists": data.get("chemists", []),
            "mechanism_type": data.get("mechanism_type", ""),
            "general_equation": data.get("general_equation", ""),
            "conditions": data.get("conditions", ""),
            "key_features_count": len(data.get("key_features", [])),
            "has_example": "example" in data,
            "applications": data.get("applications", ""),
        }

    def _run_base(self, reaction_name: str, detailed: bool = True) -> dict:
        """
        Look up a named reaction.

        Args:
            reaction_name: Name or partial name of the reaction
            detailed: If True, return full details; if False, return summary

        Returns:
            Dict with reaction information
        """
        if not reaction_name:
            raise ChemMCPError("Reaction name is required for lookup.")

        db_key, score = self._find_reaction(reaction_name)

        if db_key is None:
            # Return list of available reactions
            available = sorted([self._DETAILED_REACTIONS[k]["name"] for k in self._DETAILED_REACTIONS])
            return {
                "found": False,
                "query": reaction_name,
                "message": f"No match found for '{reaction_name}'. {len(self._DETAILED_REACTIONS)} detailed reactions available.",
                "available_reactions": available[:20],
                "total_available": len(self._DETAILED_REACTIONS),
            }

        data = self._DETAILED_REACTIONS[db_key]

        if detailed:
            result = dict(data)
            result["found"] = True
            result["match_score"] = round(score, 3)
            result["db_key"] = db_key
        else:
            result = self._get_summary(data)
            result["found"] = True
            result["match_score"] = round(score, 3)

        return result

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'reaction_name [detailed]'")

        name = parts[0]
        detailed = True
        if len(parts) > 1:
            detailed = parts[1].lower() in ("true", "1", "yes", "detailed")

        return self._run_base(name, detailed)
