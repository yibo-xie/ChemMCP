import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ResonanceStructureGenerator(BaseTool):
    """
    共振结构式生成工具 - 为常见官能团和分子生成主要共振结构式。
    包含电子推动箭头描述、相对稳定性排序、贡献度分析。
    """
    __version__ = "0.1.0"
    name             = "ResonanceStructureGenerator"
    func_name        = "generate_resonance_structures"
    description      = "Generate all major resonance structures for a molecule with electron-pushing arrow descriptions, stability ranking, and contribution analysis."
    implementation_description = "Pattern-based resonance structure generator covering carboxylate, nitro, enolate, benzene derivatives, amide, allylic/benzylic systems, carbocations, radicals, and more. Uses electron-pushing formalism with curved arrow notation."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Resonance", "Electron Delocalization", "Curved Arrows", "Structural Chemistry", "Molecular Orbital"]
    required_envs    = []

    code_input_sig   = [
        ("molecule", "str", "N/A", "Molecule: name (e.g., 'carboxylate', 'nitrobenzene', 'acetone enolate'), SMILES, or formula."),
        ("show_curved_arrows", "bool", "True", "Whether to include curved arrow (electron pushing) descriptions."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Molecule query. Example: 'nitromethane' or 'phenolate anion'."),
    ]

    output_sig       = [
        ("result", "str", "All major resonance structures with ASCII/descriptive representation, electron-pushing arrows, stability ranking, and contribution analysis."),
    ]

    examples         = [
        {
            "code_input": {"molecule": "carboxylate", "show_curved_arrows": True},
            "text_input": {"input_params": "carboxylate"},
            "output": {"result": "Resonance structures for carboxylate ion..."},
        },
        {
            "code_input": {"molecule": "nitro group", "show_curved_arrows": True},
            "text_input": {"input_params": "nitro"},
            "output": {"result": "Resonance structures for -NO2..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_database()

    def _build_database(self):
        """Build comprehensive resonance structure database."""
        # Each entry: {name: {structures: [...], total_description: str}}
        self.resonance_db = {
            "carboxylate": {
                "aliases": ["carboxylate", "RCOO-", "acetate", "formate", "benzoate", "-COO-", "carboxylic acid anion"],
                "description": "The carboxylate anion is one of the most important resonance-stabilized systems in organic chemistry. The negative charge is delocalized equally over both oxygen atoms.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "  O⁻\n   ||\n R-C\n   |\n   O",
                        "charge_distribution": "C=O double bond, O⁻ single bond (negative on single-bonded O)",
                        "contribution": "~50%",
                        "stability": "Major contributor — no charge separation beyond the anion itself",
                        "formal_charges": {"O_single": -1, "O_double": 0, "C": 0},
                    },
                    {
                        "number": 2,
                        "representation": "   O\n   ||\n R-C\n   |\n  O⁻",
                        "charge_distribution": "C=O double bond to other oxygen, O⁻ on former carbonyl O",
                        "contribution": "~50%",
                        "stability": "Equivalent to Structure 1 — degenerate pair",
                        "formal_charges": {"O_double": -1, "O_single": 0, "C": 0},
                        "curved_arrows": "Curved arrow from O⁻ lone pair → C-O π* (or equivalently: C=O π electrons move to O)",
                    },
                ],
                "key_features": [
                    "Both C-O bonds are identical in length (~1.26 Å, between single (1.43) and double (1.23) Å)",
                    "Both oxygens are chemically equivalent (NMR shows one signal)",
                    "Negative charge is delocalized over two electronegative atoms → very stable",
                    "This explains why carboxylic acids are much more acidic than alcohols (pKa ~5 vs ~16-18)",
                ],
                "hybrid": "Real structure is a hybrid with -½ charge on each oxygen and partial double bond character (bond order ~1.5) for both C-O bonds",
            },

            "nitro": {
                "aliases": ["nitro", "NO2", "nitro group", "nitromethane", "nitrobenzene", "-NO2"],
                "description": "The nitro group has extensive resonance delocalization involving N=O π bonding and formal positive charge on nitrogen.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "     O\n     ║\n R-N⁺\n     \\/\n      O⁻",
                        "charge_distribution": "N⁺ with one N=O double bond, one N-O⁻ single bond",
                        "contribution": "~50%",
                        "stability": "Major contributor — equivalent to structure 2",
                        "formal_charges": {"N": +1, "O_double": 0, "O_single": -1},
                    },
                    {
                        "number": 2,
                        "representation": "    ⁻O\n    //\n R-N⁺\n     ║\n     O",
                        "charge_distribution": "Other oxygen now bears the negative charge",
                        "contribution": "~50%",
                        "stability": "Degenerate with structure 1",
                        "formal_charges": {"N": +1, "O_single": -1, "O_double": 0},
                        "curved_arrows": "π electrons of N=O move to O; lone pair on O⁻ forms new π bond to N",
                    },
                ],
                "key_features": [
                    "Both N-O bonds are equal length (~1.22 Å, between N-O single 1.47 and N=O double 1.15)",
                    "Nitrogen bears a full +1 formal charge in all resonance forms",
                    "Strong -I and -R effects make NO2 extremely electron-withdrawing",
                    "Nitro group is planar (sp² nitrogen)",
                ],
                "hybrid": "N with partial + charge, each O with partial -½ charge, both N-O bonds have bond order ~1.5",
            },

            "enolate": {
                "aliases": ["enolate", "enolate ion", "ketone enolate", "acetylacetone enolate", "beta-dicarbonyl enolate", "1,3-dicarbonyl"],
                "description": "Enolates are among the most important nucleophiles in organic chemistry. The negative charge is delocalized over carbon and oxygen.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "  O:\n  ||\n -C-C=C⁻\n    |",
                        "charge_distribution": "Carbanion form — negative charge on α-carbon, C=O intact",
                        "contribution": "~15-30% (depends on system)",
                        "stability": "Less stable (negative charge on less electronegative C)",
                        "formal_charges": {"C_alpha": -1, "O": 0},
                    },
                    {
                        "number": 2,
                        "representation": "  O⁻\n  |\n  C=C-C\n    |",
                        "charge_distribution": "Oxyanion form — negative charge on oxygen, C=C double bond",
                        "contribution": "~70-85% (depends on system)",
                        "stability": "More stable (negative charge on more electronegative O)",
                        "formal_charges": {"O": -1, "C_alpha": 0},
                        "curved_arrows": "C=C π electrons move toward O; C=O π electrons move to O (or: Cα lone pair forms π bond, C=O π breaks putting electrons on O)",
                    },
                ],
                "key_features": [
                    "The oxyanion form dominates (more stable due to electronegativity)",
                    "But the carbanion character is crucial for reactivity at carbon (C-alkylation vs O-alkylation)",
                    "Hard electrophiles (Li⁺, hard metal cations) favor O-coordination → more C-character",
                    "Soft electrophiles (Pd, alkyl halides) react more at carbon",
                    "In β-diketones/β-ketoesters, enolate is even more stabilized (two carbonyls)",
                ],
                "hybrid": "Negative charge mostly on O but with significant delocalization to Cα; C-O bond has partial single, C=C partial double character",
            },

            "amide": {
                "aliases": ["amide", "peptide bond", "RCONH2", "RCONHR", "RCONR2", "acetamide", "DMF"],
                "description": "Amide resonance is fundamental to protein structure and reactivity. The C-N bond has partial double bond character.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "   O\n   ║\n R-C\n   |\n  N-R2",
                        "charge_distribution": "Neutral form: C=O double bond, C-N single bond, N lone pair",
                        "contribution": "~40-60% (less than you'd think!)",
                        "stability": "Minor contributor despite looking 'normal'",
                        "formal_charges": {},
                    },
                    {
                        "number": 2,
                        "representation": "  O⁻\n  |\n R-C\n  ║\n  N⁺-R2",
                        "charge_distribution": "Zwitterionic form: C-O⁻ single, C=N⁺ double, N bears + charge",
                        "contribution": "~40-60%",
                        "stability": "Surprisingly important contributor!",
                        "formal_charges": {"O": -1, "N": +1},
                        "curved_arrows": "N lone pair donates into C=O π* (forming C=N); C=O π electrons move to O",
                    },
                ],
                "key_features": [
                    "C-N bond length: ~1.33 Å (vs 1.47 for typical C-N single, 1.27 for C=N double) → partial double bond",
                    "Amide N is sp² hybridized (planar), not sp³",
                    "Restricted rotation about C-N bond (barrier ~15-20 kcal/mol) → cis/trans isomerism in peptides",
                    "N is much less basic than amines (pKa of conjugate acid ≈ 0-1 vs ~10 for amines)",
                    "C=O IR stretch: ~1680 cm⁻¹ (lower than typical 1715 for ketones, indicating reduced bond order)",
                    "This resonance explains protein secondary structure (β-sheet, α-helix H-bonding patterns)",
                ],
                "hybrid": "C-N partial double bond (bond order ~1.3-1.4), C=O slightly weakened (bond order < 2), N slightly positive, O slightly negative",
            },

            "phenolate": {
                "aliases": ["phenolate", "phenoxide", "phenol anion", "PhO-", "cresolate"],
                "description": "Phenoxide anion — the reason phenols (pKa ~10) are far more acidic than alcohols (pKa ~16-18). Negative charge delocalized into aromatic ring.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "  O⁻ attached directly to benzene ring\n  (negative charge localized on oxygen initially)",
                        "contribution": "~20%",
                        "stability": "Least contributing — charge localized on one atom",
                        "formal_charges": {"O": -1},
                    },
                    {
                        "number": 2,
                        "representation": "  O (neutral)\n  ║\n  ortho-C of benzene ring (with negative charge)",
                        "contribution": "~20%",
                        "stability": "Negative charge on carbon (less favorable than on O)",
                        "curved_arrows": "O⁻ lone pair → ortho C (forming C=O-like π bond in ring)",
                        "formal_charges": {"O": 0, "ortho_C": -1},
                    },
                    {
                        "number": 3,
                        "representation": "  O (neutral)\n  benzene ring with negative charge at para position",
                        "contribution": "~20%",
                        "stability": "Same as ortho — degenerate for unsubstituted phenol",
                        "curved_arrows": "O⁻ lone pair → para C via conjugated pathway",
                        "formal_charges": {"O": 0, "para_C": -1},
                    },
                    {
                        "number": 4,
                        "representation": "  O (neutral)\n  benzene ring with negative charge at other ortho position",
                        "contribution": "~20%",
                        "stability": "Degenerate with ortho contributors",
                        "formal_charges": {"O": 0, "other_ortho_C": -1},
                    },
                    {
                        "number": 5,
                        "representation": "  Kekule structure flipped (ring double bonds shifted)\n  + O⁻ / ortho / para variants",
                        "contribution": "~20% (total for flipped Kekule set)",
                        "stability": "Additional structures from the other Kekule form",
                    },
                ],
                "key_features": [
                    "4 major resonance structures place negative charge on ortho (2×) and para (1×) positions",
                    "Total of 4 significant resonance contributors (plus Kekule variants)",
                    "Charge delocalization over 4 atoms (O + 3 ring carbons) → great stabilization",
                    "Explains why p-nitrophenol (pKa 7.15) is much more acidic than phenol (pKa 10.00): additional resonance structures where NO2 accepts the negative charge",
                    "Ortho/para-directing nature of O⁻ in EAS follows from these resonance forms",
                ],
                "hybrid": "Negative charge spread over O and ortho/para ring positions; C-O bond has partial double bond character",
            },

            "allyl_cation": {
                "aliases": ["allyl cation", "allylic carbocation", "CH2=CH-CH2+", "allyl carbocation"],
                "description": "The simplest resonance-stabilized carbocation. Positive charge delocalized over two terminal carbons.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "  CH2=CH-CH2⁺  (+ on terminal C3)",
                        "contribution": "50%",
                        "stability": "Equal contributor",
                        "formal_charges": {"C3": +1},
                    },
                    {
                        "number": 2,
                        "representation": "  ⁺CH2-CH=CH2  (+ on terminal C1)",
                        "contribution": "50%",
                        "stability": "Degenerate with structure 1",
                        "formal_charges": {"C1": +1},
                        "curved_arrows": "π electrons of C=C shift toward C2, forming new π bond; C2-C3 σ electrons (bonding pair) move to C3 as new π bond to C1",
                    },
                ],
                "key_features": [
                    "Both terminal C-C bonds are equal (~1.38 Å, between single 1.54 and double 1.34)",
                    "Central carbon is sp hybridized (linear geometry, ∠C-C-C ≈ 120deg actually due to vacant p orbital)",
                    "Positive charge distributed equally over two terminal carbons",
                    "Nucleophilic attack can occur at either terminal carbon → mixture of products unless unsymmetrical substitution",
                    "Allylic rearrangements common (SN1' reactions)",
                ],
                "hybrid": "Symmetric cation with +½ charge on each terminal carbon, partial double bond character throughout C1-C2-C3 framework",
            },

            "benzyl_cation": {
                "aliases": ["benzyl cation", "PhCH2+", "benzylic carbocation"],
                "description": "Benzylic carbocation stabilized by resonance with the benzene ring. More stable than allyl cation.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "  Ph-CH2⁺  (positive on benzylic carbon)",
                        "contribution": "~20%",
                        "formal_charges": {"CH2": +1},
                    },
                    {
                        "number": 2-5,
                        "representation": "  Ortho/para positions of benzene bear positive charge\n  (4 structures: 2 ortho + 1 para × 2 Kekule forms)",
                        "contribution": "~80% combined (4 structures × ~20% each)",
                        "curved_arrows": "Ring π electrons donate into empty p orbital of CH2+",
                        "formal_charges": {"ring_C": +1, "CH2": 0},
                    },
                ],
                "key_features": [
                    "Positive charge delocalized over benzylic carbon + ortho/para ring positions (7 atoms total!)",
                    "Much more stable than simple primary carbocation (comparable to tertiary)",
                    "Benzyl cation stability: PhCH2+ > allyl CH2=CH-CH2+ > t-Bu+ > i-Pr+ > primary",
                    "Explains why benzylic SN1 reactions proceed readily",
                ],
            },

            "allyl_radical": {
                "aliases": ["allyl radical", "allylic radical", "CH2=CH-CH2•"],
                "description": "Resonance-stabilized radical. Unpaired electron delocalized over two terminal carbons.",
                "structures": [
                    {
                        "number": 1,
                        "representation": "  CH2=CH-CH2•  (radical on C3)",
                        "contribution": "50%",
                    },
                    {
                        "number": 2,
                        "representation": "  •CH2-CH=CH2  (radical on C1)",
                        "contribution": "50%",
                        "curved_arrows": "Single electron (fishhook arrows): one electron from π bond moves to C3, one electron from C2-C3 bond moves to C1",
                    },
                ],
                "key_features": [
                    "EPR spectroscopy shows spin density at both terminal carbons",
                    "Radical reactions occur at both termini",
                    "Allylic bromination (NBS) proceeds via this intermediate",
                ],
            },

            "benzene": {
                "aliases": ["benzene", "C6H6", "aromatic", "aromatic ring"],
                "description": "The archetype of aromatic resonance. Two equivalent Kekule structures plus Dewar/ionic forms (minor).",
                "structures": [
                    {
                        "number": 1,
                        "representation": "   Kekule A: alternating single-double bonds\n    C1=C2-C3=C4-C5=C6-C1 (with C1-C6 single)",
                        "contribution": "Major (but not exclusive)",
                    },
                    {
                        "number": 2,
                        "representation": "   Kekule B: shifted pattern\n    C1-C2=C3-C4=C5-C6=C1 (with C1-C6 double)",
                        "contribution": "Equivalent to Kekule A",
                        "curved_arrows": "All three π bonds shift by one position simultaneously",
                    },
                ],
                "key_features": [
                    "All C-C bonds equal (1.397 Å, between single 1.54 and double 1.34)",
                    "Fully delocalized 6π electron system (Hückel's rule: 4n+2, n=1)",
                    "Resonance energy: ~36 kcal/mol (152 kJ/mol)",
                    "Dewar benzene structures contribute negligibly (<5%)",
                    "Aromatic compounds resist addition reactions (would destroy aromaticity)",
                ],
                "hybrid": "Perfect hexagon with uniform bond order 1.5, circular π cloud above and below ring plane",
            },
        }

    def _find_entry(self, query: str):
        """Find resonance entry by name or alias."""
        q = query.lower().strip()
        
        for key, entry in self.resonance_db.items():
            if key == q:
                return entry, key
            for alias in entry.get("aliases", []):
                if q == alias.lower() or q in alias.lower() or alias.lower() in q:
                    return entry, key
        
        return None, None

    def _run_base(self, molecule: str, show_curved_arrows: bool = True) -> str:
        """Generate resonance structures."""
        entry, key = self._find_entry(molecule)
        
        if entry is None:
            available = ", ".join(self.resonance_db.keys())
            return f"## Resonance Structure Generator: `{molecule}`\n\n### ⚠️ Not Found\n\nNo resonance data for '{molecule}'.\n\n**Available systems:** {available}\n\nTry: carboxylate, nitro, enolate, amide, phenolate, allyl_cation, benzyl_cation, allyl_radical, benzene"

        parts = [f"## Resonance Structures: {entry.get('description', '').split('.')[0]}\n"]
        parts.append(f"**System:** `{key}`\n")
        parts.append(f"### 📖 Overview\n{entry['description']}\n")

        parts.append("### 🔬 Resonance Structures\n")
        for struct in entry["structures"]:
            parts.append(f"#### Structure {struct['number']}")
            
            # Show representation
            repr_text = struct.get("representation", "")
            if repr_text:
                parts.append(f"```\n{repr_text}\n```")
            
            parts.append(f"- **Charge Distribution:** {struct['charge_distribution']}")
            parts.append(f"- **Contribution:** {struct['contribution']}")
            parts.append(f"- **Stability:** {struct['stability']}")
            
            fc = struct.get("formal_charges", {})
            if fc:
                fc_str = ", ".join(f"{k}: {'+' if v > 0 else ''}{v}" for k, v in fc.items())
                parts.append(f"- **Formal Charges:** {fc_str}")
            
            if show_curved_arrows and "curved_arrows" in struct:
                parts.append(f"- **🏹 Curved Arrow Electron Pushing:** {struct['curved_arrows']}")
            
            parts.append("")

        if entry.get("key_features"):
            parts.append("### ✨ Key Features\n")
            for feature in entry["key_features"]:
                parts.append(f"- {feature}")
            parts.append("")

        if entry.get("hybrid"):
            parts.append(f"### 🔄 True Hybrid Structure\n> {entry['hybrid']}\n")

        # General resonance rules reminder
        parts.append("""
### 📐 Resonance Rules Applied

1. **Structures must be valid Lewis structures** (octet rule, correct electron count)
2. **Only electrons move** — nuclei never change position
3. **Preserve total number of unpaired electrons** (fishhook arrows for radicals)
4. **Major contributors have:** maximum covalent bonds, minimum formal charges, negative charge on electronegative atoms
5. **All structures contribute to the real hybrid** — the actual molecule is none of the individual forms
6. **More resonance structures = greater stabilization** (but only if they are significant/equivalent)
""")

        return "\n".join(parts)

    def _run_text(self, input_params: str) -> str:
        input_params = input_params.strip()
        if not input_params:
            raise ChemMCPError("Please provide a molecule name. Example: 'carboxylate', 'nitro', 'enolate'")
        return self._run_base(input_params)
