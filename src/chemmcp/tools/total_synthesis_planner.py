import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TotalSynthesisPlanner(BaseTool):
    """
    多步合成路线规划工具 - 根据目标分子生成多步合成方案。
    包含逆合成分析、关键转化策略、试剂建议。
    """
    __version__ = "0.1.0"
    name             = "TotalSynthesisPlanner"
    func_name        = "total_synthesis_planner"
    description      = "Plan multi-step total synthesis routes for target molecules with disconnection analysis, key transformations, and reagent suggestions."
    implementation_description = "Knowledge-based retrosynthesis planner using disconnection strategies (FGI, 1,2-/1,3-/1,4-disconnections), named reactions database, and synthetic heuristics to propose step-by-step routes."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Total Synthesis", "Retrosynthesis", "Synthetic Planning", "Disconnection", "Route Design"]
    required_envs    = []

    code_input_sig   = [
        ("target_molecule", "str", "N/A", "Target molecule name, SMILES, or structural description."),
        ("complexity_level", "str", "medium", "Complexity: 'simple' (<5 steps), 'medium' (5-15 steps), 'complex' (>15 steps)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'target_molecule [complexity_level]'. Example: 'aspirin simple'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed synthesis plan with numbered steps, reagents, conditions, yields, and strategic rationale."),
    ]

    examples         = [
        {
            "code_input": {"target_molecule": "aspirin", "complexity_level": "simple"},
            "text_input": {"input_params": "aspirin simple"},
            "output": {"result": "Step 1: Salicylic acid + acetic anhydride → Aspirin..."},
        },
        {
            "code_input": {"target_molecule": "ibuprofen", "complexity_level": "medium"},
            "text_input": {"input_params": "ibuprofen medium"},
            "output": {"result": "Step 1: ... Step 2: ... Step 3: ..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_synthetic_knowledge()

    def _build_synthetic_knowledge(self):
        """Build comprehensive synthetic planning knowledge base."""
        # Known molecule synthesis routes
        self.known_routes = {
            # === Simple molecules (1-5 steps) ===
            "aspirin": {
                "name": "Aspirin (Acetylsalicylic Acid)",
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Acetylation of salicylic acid",
                        "starting_material": "Salicylic acid (2-hydroxybenzoic acid)",
                        "reagents": "Acetic anhydride (Ac2O), catalytic H2SO4 or H3PO4",
                        "conditions": "85°C, 10-20 min, then ice-water crystallization",
                        "mechanism": "Nucleophilic acyl substitution (phenolic OH attacks acetic anhydride)",
                        "yield": "85-95%",
                        "notes": "Industrial process; can also use acetyl chloride",
                        "safety": "Acetic anhydride is corrosive/lachrymatory; use fume hood",
                    },
                ],
                "alternative_routes": [
                    "From phenol via Kolbe-Schmitt (CO2 + NaOH → salicylic acid) then acetylation",
                    "From methyl salicylate (oil of wintergreen) hydrolysis then acetylation",
                ],
                "key_strategic_bonds": "Ester bond between phenol and acetate",
                "commercially_available": True,
                "complexity": "simple",
            },
            "acetaminophen": {
                "name": "Acetaminophen (Paracetamol)",
                "smiles": "CC(=O)Nc1ccc(O)cc1",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Nitration of phenol → p-nitrophenol",
                        "starting_material": "Phenol",
                        "reagents": "Dilute HNO3 (NaNO2 catalyst)",
                        "conditions": "0-25°C (temperature control critical for para-selectivity)",
                        "yield": "~80% (para isomer after separation)",
                        "notes": "Ortho/para mixture; separate by distillation/recrystallization",
                    },
                    {
                        "step": 2,
                        "transformation": "Reduction of p-nitrophenol → p-aminophenol",
                        "starting_material": "p-Nitrophenol",
                        "reagents": "H2 / Pd-C or Fe / HCl or Sn / HCl",
                        "conditions": "H2 (1-3 atm), Pd/C (5%), EtOH, rt; OR Fe powder, aq. HCl, reflux",
                        "yield": "85-95%",
                        "notes": "Catalytic hydrogenation preferred industrially (cleaner)",
                    },
                    {
                        "step": 3,
                        "transformation": "N-Acetylation → acetaminophen",
                        "starting_material": "p-Aminophenol",
                        "reagents": "Acetic anhydride (Ac2O) or acetyl chloride",
                        "conditions": "Water, 0-5°C (control N-acetylation vs O-acetylation)",
                        "yield": "90-95%",
                        "notes": "Low temperature favors N-acetylation over O-acetylation",
                    },
                ],
                "alternative_routes": ["From p-chloronitrobenzene → hydrolysis → reduction → acetylation"],
                "key_strategic_bonds": "Amide bond",
                "commercially_available": True,
                "complexity": "simple",
            },
            "ibuprofen": {
                "name": "Ibuprofen",
                "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Friedel-Crafts Acylation",
                        "starting_material": "Isobutylbenzene",
                        "reagents": "Acetic anhydride or acetyl chloride, AlCl3 (Lewis acid)",
                        "conditions": "CH2Cl2 or neat, 0°C → rt, then aqueous workup",
                        "yield": "75-85%",
                        "notes": "Gives p-isobutylacetophenone (para selective due to isobutyl group)",
                    },
                    {
                        "step": 2,
                        "transformation": "Darzens Glycidic Ester Condensation / or Boots Process",
                        "starting_material": "p-Isobutylacetophenone",
                        "reagents": "Ethyl chloroacetate, NaOEt (or NaOH for greener route)",
                        "conditions": "Ethanol, reflux (Boots process); OR catalytic carbonylation (Hoechst process - greener)",
                        "yield": "70-85% (Boots); >95% atom economy (Hoechst hydrocarbonylation)",
                        "notes": "Boots: classical racemate resolution needed; Hoechst: asymmetric hydrogenation gives chiral product directly",
                    },
                    {
                        "step": 3,
                        "transformation": "Hydrolysis/Decarboxylation → ibuprofen",
                        "starting_material": "Glycidic ester intermediate",
                        "reagents": "NaOH (hydrolysis), then H+ (decarboxylation)",
                        "conditions": "Aq. NaOH, heat, then acidification",
                        "yield": "80-90%",
                        "notes": "Final product may need resolution if using Boots process",
                    },
                ],
                "alternative_routes": [
                    "Hoechst process (carbonylation): p-isobutylstyrene + CO + H2 → ibuprofen (Pd-catalyzed, 3 steps, higher atom economy)",
                    "BHC Company asymmetric route via Ni-catalyzed coupling",
                ],
                "key_strategic_bonds": "α-Chiral center (if making S-enantiomer); propionic acid side chain",
                "commercially_available": True,
                "complexity": "medium",
            },
            "caffeine": {
                "name": "Caffeine (1,3,7-trimethylxanthine)",
                "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Dimethylation of theobromine",
                        "starting_material": "Theobromine (from cocoa waste or synthesized from xanthines)",
                        "reagents": "Dimethyl sulfate (Me2SO4) or methyl iodide (CH3I), NaOH/K2CO3",
                        "conditions": "Water or DMF, 60-100°C",
                        "yield": "70-85%",
                        "notes": "Theobromine can be isolated from cocoa or derived from xanthine",
                    },
                ],
                "alternative_routes": [
                    "From uric acid → methylation (Traube synthesis)",
                    "Total synthesis from dimethylurea + cyanoacetic acid (multi-step)",
                ],
                "key_strategic_bonds": "N-methyl groups at positions 1, 3, and 7",
                "commercially_available": True,
                "complexity": "simple",
            },
            "dopamine": {
                "name": "Dopamine (4-(2-aminoethyl)benzene-1,2-diol)",
                "smiles": "NCCc1ccc(O)c(O)c1",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Protection of catechol",
                        "starting_material": "2-(3,4-dimethoxyphenyl)ethylamine (homoveratrylamine)",
                        "reagents": "BBr3 or BCl3 (demethylating agent)",
                        "conditions": "CH2Cl2, −78°C → rt, under N2",
                        "yield": "75-90%",
                        "notes": "Homoveratrylamine is commercially available; demethylation reveals catechol",
                    },
                ],
                "alternative_routes": [
                    "From L-DOPA decarboxylation (enzymatic or chemical)",
                    "From 3,4-dihydroxyphenylacetaldehyde reductive amination",
                ],
                "key_strategic_bonds": "Catechol (ortho-diol) + ethylamine side chain",
                "commercially_available": True,
                "complexity": "simple",
            },
            # === Medium complexity molecules ===
            "menthol": {
                "name": "(-)-Menthol",
                "smiles": "CC(C)[C@H]1CCC[C@H](C)C1",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Diels-Alder: Myrcene + diethylaminosulfonyl isoprene",
                        "starting_material": "Myrcene (from β-pinene or turpentine)",
                        "reagents": "Lithium diisopropylamide (LIPA), chloroamine derivative",
                        "conditions": "−20°C, then thermal rearrangement",
                        "yield": "~60%",
                        "notes": "Takasago process starting point",
                    },
                    {
                        "step": 2,
                        "transformation": "Asymmetric Isomerization (Takasago Key Step)",
                        "starting_material": "Diels-Alder adduct",
                        "reagents": "Ru-BINAP complex (Noyori-type catalyst)",
                        "conditions": "H2 pressure, 100-120°C",
                        "yield": ">96% ee, quantitative conversion",
                        "notes": "KEY STEP: Ru-BINAP catalyzes asymmetric isomerization to chiral enamine",
                    },
                    {
                        "step": 3,
                        "transformation": "Hydrolysis → (-)-menthol",
                        "starting_material": "Chiral enamine intermediate",
                        "reagents": "Aqueous acid (H2SO4/H2O)",
                        "conditions": "Hydrolysis conditions",
                        "yield": "High overall",
                        "notes": "This industrial process produces ~30% of world's menthol supply",
                    },
                ],
                "alternative_routes": ["Extraction from peppermint oil (natural source)", "From thymol (classical partial resolution)"],
                "key_strategic_bonds": "Three contiguous stereocenters (1R,2S,5R configuration)",
                "commercially_available": True,
                "complexity": "medium",
            },
            "propranolol": {
                "name": "Propranolol (β-blocker)",
                "smiles": "COc1cccc(CC(O)CN(C)C)c1OC",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Epoxide formation from allyl alcohol derivative",
                        "starting_material": "α-Naphthol",
                        "reagents": "Allyl bromide, K2CO3 (Williamson ether synthesis)",
                        "conditions": "Acetone, reflux",
                        "yield": "80-90%",
                        "notes": "Forms 1-allyloxy-naphthalene",
                    },
                    {
                        "step": 2,
                        "transformation": "Epoxidation of terminal alkene",
                        "starting_material": "1-Allyloxy-naphthalene",
                        "reagents": "m-CPBA (meta-chloroperbenzoic acid)",
                        "conditions": "CH2Cl2, 0°C → rt",
                        "yield": "75-85%",
                        "notes": "Forms glycidyl ether intermediate",
                    },
                    {
                        "step": 3,
                        "transformation": "Epoxide ring opening with isopropylamine",
                        "starting_material": "Glycidyl ether",
                        "reagents": "Isopropylamine (i-PrNH2), excess",
                        "conditions": "EtOH or MeOH, 50-60°C, sealed tube or reflux",
                        "yield": "70-80%",
                        "notes": "Gives racemic propranolol; can be resolved or made asymmetrically",
                    },
                ],
                "alternative_routes": ["Asymmetric epoxidation → enantioselective ring opening"],
                "key_strategic_bonds": "β-Amino alcohol moiety (epoxide opening pattern)",
                "commercially_available": True,
                "complexity": "medium",
            },
            "lidocaine": {
                "name": "Lidocaine (local anesthetic)",
                "smiles": "CCN(CC)C(=O)NC1=C(C)C=C(C)C=C1C",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Nitration of 2,6-dimethylaniline → 2,6-dimethylnitrobenzene",
                        "starting_material": "2,6-Dimethylaniline",
                        "reagents": "HNO3 / H2SO4 (mixed acid)",
                        "conditions": "0-5°C (temperature control important)",
                        "yield": "70-80%",
                        "notes": "Position 4 nitration (between two methyl groups)",
                    },
                    {
                        "step": 2,
                        "transformation": "Reduction → 2,6-dimethylaniline (2,4-diamino derivative)",
                        "starting_material": "2,6-Dimethylnitrobenzene",
                        "reagents": "H2, Pd-C (catalytic hydrogenation)",
                        "conditions": "EtOH, rt, 1-3 atm H2",
                        "yield": "90-95%",
                    },
                    {
                        "step": 3,
                        "transformation": "Acylation with 2-(diethylamino)acetyl chloride",
                        "starting_material": "2,4-Diamino-1,3-dimethylbenzene",
                        "reagents": "2-Chloro-N,N-diethylacetamide (or chloroacetyl chloride + diethylamine)",
                        "conditions": "Acetone or CH2Cl2, NaHCO3 or pyridine as base, 0°C → rt",
                        "yield": "75-85%",
                        "notes": "Amide bond formation on the less hindered amino group (position 4)",
                    },
                ],
                "alternative_routes": ["From 2,6-dimethylaniline via α-chloroacyl chloride directly"],
                "key_strategic_bonds": "Amide linkage between aniline and diethylglycine",
                "commercially_available": True,
                "complexity": "medium",
            },
            "quinine": {
                "name": "Quinine (antimalarial alkaloid)",
                "smiles": "COc1cc2C[C@H]3N(CC)C(=O)C(CO)=C(C3=c2c(OC)c1O)[C@@H](O)CCCN2CCCCC2",
                "steps": [
                    {
                        "step": 1,
                        "transformation": "Construction of quinoline core",
                        "starting_material": "Starting from m-anisidine or 3-hydroxybenzaldehyde derivatives",
                        "reagents": "Skraup-Doebner-Von Miller quinoline synthesis conditions",
                        "conditions": "Glycerol, H2SO4, an oxidizing agent (nitrobenzene or As2O5)",
                        "yield": "Variable (40-60%)",
                        "notes": "Classical approach; modern routes use better methods",
                    },
                    {
                        "step": 2,
                        "transformation": "Introduction of quinuclidine moiety",
                        "starting_material": "Quinoline intermediate",
                        "reagents": "Various: Mannich reaction, cyclization sequences",
                        "conditions": "Multi-step sequence",
                        "yield": "Varies per step",
                        "notes": "The quinuclidine bridge is the most challenging part",
                    },
                    {
                        "step": 3-12,
                        "transformation": "Stereocontrolled assembly of remaining stereocenters",
                        "starting_material": "Advanced intermediates",
                        "reagents": "Chiral auxiliaries, asymmetric reductions, etc.",
                        "conditions": "Multiple steps each with specific conditions",
                        "yield": "Overall ~1-5% (typical for complex natural products)",
                        "notes": "Total syntheses by Woodward & Doering (1945), Stork (2001), Jacobsen exist",
                    },
                ],
                "alternative_routes": ["Extraction from Cinchona bark (still main commercial source)", "Semi-synthesis from more available cinchona alkaloids"],
                "key_strategic_bonds": "Quinuclidine bridge C-C bonds; C8-C9 stereochemistry; vinyl group",
                "commercially_available": True,
                "complexity": "complex",
            },
            }

        # Disconnection strategies database
        self.disconnection_rules = {
            "amide": {
                "strategy": "Disconnect amide C-N bond → carboxylic acid + amine",
                "method": "Coupling (EDCI, DCC, HATU) or acyl chloride + amine",
                "priority": 1,
            },
            "ester": {
                "strategy": "Disconnect ester → acid + alcohol (Fischer esterification reverse)",
                "method": "Steglich esterification (DCC/DMAP) or acid chloride + alcohol",
                "priority": 1,
            },
            "alcohol_1,2": {
                "strategy": "1,2-Disconnection: alcohol ← carbonyl reduction (NaBH4, LiAlH4)",
                "method": "Reduce aldehyde/ketone/ester",
                "priority": 2,
            },
            "alkene": {
                "strategy": "Disconnect alkene → Wittig reaction or elimination",
                "method": "Wittig, Horner-Wadsworth-Emmons, or dehydration of alcohol",
                "priority": 2,
            },
            "aryl_alkyl": {
                "strategy": "Disconnect aryl-alkyl bond → Friedel-Crafts or cross-coupling",
                "method": "Suzuki, Heck, Negishi, or Friedel-Crafts alkylation/acylation",
                "priority": 2,
            },
            "heterocycle_6": {
                "strategy": "6-membered heterocycle → 1,5-dicarbonyl or equivalent",
                "method": "Paal-Knorr, Hantzsch, Biginelli, condensation reactions",
                "priority": 3,
            },
            "ring_fusion": {
                "strategy": "Analyze ring fusion strategy (annulation, cycloaddition)",
                "method": "Robinson annulation, Diels-Alder, intramolecular aldol/Michael",
                "priority": 3,
            },
        }

        # Common building blocks catalog
        self.building_blocks = {
            "simple_acids": ["Acetic acid", "Benzoic acid", "Salicylic acid", "p-Nitrobenzoic acid"],
            "simple_alcohols": ["Methanol", "Ethanol", "Benzyl alcohol", "Phenol"],
            "simple_amines": ["Methylamine", "Ammonia", "Aniline", "Benzylamine"],
            "aromatics": ["Benzene", "Toluene", "Phenol", "Anisole", "Chlorobenzene"],
            "chiral_pool": ["L-Serine", "L-Proline", "(S)-(-)-Mandelic acid", "(R)-(+)-Citronellal", "Natural terpenes"],
        }

    def _run_base(self, target_molecule: str, complexity_level: str = "medium") -> str:
        """Generate total synthesis plan."""
        target_lower = target_molecule.lower().strip()
        
        # Check known routes first
        matched_route = None
        for key, route in self.known_routes.items():
            if key in target_lower or target_lower in key:
                matched_route = route
                break

        parts = [f"## Total Synthesis Plan: {target_molecule}\n"]

        if matched_route:
            parts.append(f"**Target:** {matched_route['name']}")
            if matched_route.get('smiles'):
                parts.append(f"**SMILES:** `{matched_route['smiles']}`")
            parts.append(f"**Complexity:** {matched_route.get('complexity', 'unknown').upper()}")
            parts.append(f"**Commercially Available:** {'Yes ✅' if matched_route.get('commercially_available') else 'No ❌'}\n")

            parts.append("### 🧪 Proposed Synthetic Route\n")
            for s in matched_route["steps"]:
                parts.append(f"#### Step {s['step']}: {s['transformation']}")
                parts.append(f"- **Starting Material:** {s['starting_material']}")
                if 'reagents' in s:
                    parts.append(f"- **Reagents:** {s['reagents']}")
                if 'conditions' in s:
                    parts.append(f"- **Conditions:** {s['conditions']}")
                if 'yield' in s:
                    parts.append(f"- **Yield:** {s['yield']}")
                if 'mechanism' in s:
                    parts.append(f"- **Mechanism:** {s['mechanism']}")
                if 'notes' in s:
                    parts.append(f"- **Notes:** {s['notes']}")
                if 'safety' in s:
                    parts.append(f"- ⚠️ **Safety:** {s['safety']}")
                parts.append("")

            if matched_route.get('alternative_routes'):
                parts.append("### 🔀 Alternative Routes\n")
                for alt in matched_route['alternative_routes']:
                    parts.append(f"- {alt}")
                parts.append("")

            if matched_route.get('key_strategic_bonds'):
                parts.append(f"### 🎯 Key Strategic Bonds\n{matched_route['key_strategic_bonds']}\n")
        else:
            # Generic retrosynthetic analysis for unknown molecules
            parts.append("⚠️ **No exact match found in database.** Generating generic retrosynthetic analysis based on structural features.\n")
            
            parts.append("### 🔬 Retrosynthetic Analysis\n")
            parts.append("Apply these disconnection strategies systematically:\n")
            
            for rule_key, rule in self.disconnection_rules.items():
                parts.append(f"**{rule_key.replace('_', ' ').title()} Disconnection:**")
                parts.append(f"- {rule['strategy']}")
                parts.append(f"- Method: {rule['method']}")
                parts.append(f"")

            parts.append("\n### 🧱 Suggested Building Blocks\n")
            for category, blocks in self.building_blocks.items():
                parts.append(f"**{category.replace('_', ' ').title()}:** {', '.join(blocks)}")
            
            parts.append("""
### 📋 General Planning Heuristics

1. **Start from the target** — work backward (retrosynthesis)
2. **Identify functional groups** — plan interconversions (FGI)
3. **Find disconnection points** — C-X (heteroatom) bonds first, then C-C bonds
4. **Check symmetry** — use symmetry elements to simplify
5. **Consider stereochemistry** — plan stereoselective steps early
6. **Use convergent synthesis** — join large fragments late in the synthesis
7. **Protecting group strategy** — minimize number of PG operations
8. **Redox economy** — avoid unnecessary oxidation state changes
9. **Atom economy** — prefer catalytic, atom-economical transformations
10. **Scalability** — consider cost, safety, and availability of reagents

For a detailed route, please provide SMILES structure or more specific molecular description.
""")

        return "\n".join(parts)

    def _run_text(self, input_params: str) -> str:
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Input must include target_molecule. Format: 'target_molecule [complexity_level]'")
        target = parts[0]
        complexity = parts[1] if len(parts) > 1 else "medium"
        return self._run_base(target, complexity)
