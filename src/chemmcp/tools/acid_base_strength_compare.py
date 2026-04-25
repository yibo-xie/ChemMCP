import logging
from typing import List, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AcidBaseStrengthCompare(BaseTool):
    """
    比较酸碱强度（基于结构/数据）。
    提供常见酸碱的pKa/pKb数据，并基于结构因素解释相对强弱。
    """
    __version__ = "0.1.0"
    name = "AcidBaseStrengthCompare"
    func_name = "compare_acid_base_strength"
    description = "Compare acid or base strength among multiple species, with pKa/pKb values and structural explanations (inductive effect, resonance, atomic size, hybridization, etc.)."
    implementation_description = "Uses built-in database of pKa/pKb values for common acids and bases. Ranks them by strength and provides structural reasoning for the observed order."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Acid-Base", "pKa", "pKb", "Strength", "Structure", "Inductive Effect"]
    required_envs = []

    code_input_sig = [
        ("species_list", "list", "N/A", "List of acid/base names or formulas to compare, e.g., ['HCl', 'CH3COOH', 'HF', 'H2CO3']."),
        ("compare_type", "str", "auto", "Type of comparison: 'acid' (compare acidity), 'base' (compare basicity), or 'auto' (auto-detect)."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Semicolon-separated list: 'HCl;CH3COOH;HF;H3PO4 [type]'. Type optional: acid/base/auto."),
    ]

    output_sig = [
        ("ranking", "list", "Ranked list of species from strongest to weakest, with pKa/pKb values."),
        ("structural_analysis", "str", "Explanation of strength differences based on molecular structure factors."),
        ("compare_type_used", "str", "The actual comparison type used (acid or base)."),
    ]

    examples = [
        {
            "code_input": {
                "species_list": ["HCl", "CH3COOH", "HF", "H2CO3"],
                "compare_type": "acid",
            },
            "text_input": {
                "input_str": "HCl;CH3COOH;HF;H2CO3 acid",
            },
            "output": {
                "ranking": [
                    {"name": "HCl", "pKa": -7.0, "strength": "strong"},
                    {"name": "HF", "pKa": 3.17, "strength": "weak"},
                    {"name": "CH3COOH", "pKa": 4.76, "strength": "weak"},
                    {"name": "H2CO3", "pKa": 6.35, "strength": "weak"},
                ],
                "structural_analysis": "HCl is a strong acid (complete dissociation). HF has higher pKa than expected due to strong H-F bond and hydrogen bonding. CH3COOH is a typical weak carboxylic acid. H2CO3 is very weak due to instability.",
                "compare_type_used": "acid",
            },
        },
        {
            "code_input": {
                "species_list": ["NaOH", "NH3", "aniline", "pyridine"],
                "compare_type": "base",
            },
            "text_input": {
                "input_str": "NaOH;NH3;aniline;pyridine base",
            },
            "output": {
                "ranking": [
                    {"name": "NaOH", "pKb": -1.0, "strength": "strong"},
                    {"name": "pyridine", "pKb": 8.75, "strength": "weak"},
                    {"name": "NH3", "pKb": 4.75, "strength": "weak"},
                    {"name": "aniline", "pKb": 9.38, "strength": "very weak"},
                ],
                "structural_analysis": "NaOH is a strong base. NH3 is a weak base (lone pair on N). Pyridine's N lone pair is in sp2 orbital (more s-character, less available). Aniline's lone pair is delocalized into the benzene ring, making it much weaker.",
                "compare_type_used": "base",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize pKa/pKb database."""
        # Common acids: {name: {"pKa": value, "type": "monoprotic|polyprotic", "notes": "..."}}
        self._acids = {
            # Strong mineral acids
            "HCl":      {"pKa": -7.0,   "type": "strong",     "notes": "Strong mineral acid, complete dissociation."},
            "HBr":      {"pKa": -9.0,   "type": "strong",     "notes": "Stronger than HCl due to weaker H-Br bond."},
            "HI":       {"pKa": -10.0,  "type": "strong",     "notes": "Strongest hydrogen halide acid."},
            "H2SO4":    {"pKa": -3.0,   "type": "strong",     "notes": "First dissociation strong; pKa2=1.99."},
            "HNO3":     {"pKa": -1.4,   "type": "strong",     "notes": "Strong oxoacid, resonance-stabilized conjugate base."},
            "HClO4":    {"pKa": -10.0,  "type": "strong",     "notes": "Strongest common oxoacid."},
            "HClO3":    {"pKa": -1.0,   "type": "strong",     "notes": "Strong oxoacid."},

            # Weak inorganic acids
            "HF":       {"pKa": 3.17,   "type": "weak",       "notes": "Weak due to strong H-F bond and H-bonding stabilization of HF."},
            "H3PO4":    {"pKa": 2.12,   "type": "polyprotic", "notes": "Triprotic; pKa2=7.21, pKa3=12.38."},
            "H2CO3":    {"pKa": 6.35,   "type": "polyprotic", "notes": "Diprotic; pKa2=10.33. Actually CO2(aq) equilibrium."},
            "H2S":      {"pKa": 7.02,   "type": "polyprotic", "notes": "Diprotic; pKa2≈13-19 (uncertain). Foul smell."},
            "HCN":      {"pKa": 9.31,   "type": "weak",       "notes": "Very weak; CN- stabilized by resonance."},
            "H2O":      {"pKa": 15.7,   "type": "weak",       "notes": "Very weakly amphoteric."},
            "H3BO3":    {"pKa": 9.24,   "type": "weak",       "notes": "Lewis acid behavior; accepts OH-."},

            # Carboxylic acids
            "HCOOH":    {"pKa": 3.75,   "type": "weak",       "notes": "Formic acid; no electron-donating alkyl group."},
            "CH3COOH":  {"pKa": 4.76,   "type": "weak",       "notes": "Acetic acid; methyl group slightly electron-donating."},
            "C6H5COOH": {"pKa": 4.20,   "type": "weak",       "notes": "Benzoic acid; phenyl ring mildly electron-withdrawing by induction."},

            # Phenols
            "phenol":   {"pKa": 10.0,   "type": "weak",       "notes": "Phenoxide stabilized by resonance into ring."},
            "p-nitrophenol": {"pKa": 7.15, "type": "weak",     "notes": "Nitro group strongly electron-withdrawing, strengthens acidity."},

            # Alcohols (extremely weak acids)
            "CH3CH2OH": {"pKa": 16.0,   "type": "very_weak",  "notes": "Ethanol; alkoxide is a strong base."},
        }

        # Common bases: {name: {"pKb": value, "conjugate_acid_pKa": ..., "notes": "..."}}
        self._bases = {
            "NaOH":     {"pKb": -1.0,   "type": "strong",     "notes": "Strong base, complete dissociation."},
            "KOH":      {"pKb": -1.0,   "type": "strong",     "notes": "Strong base."},
            "Ba(OH)2":  {"pKb": -1.0,   "type": "strong",     "notes": "Strong base."},
            "Ca(OH)2":  {"pKb": -1.0,   "type": "strong",     "notes": "Strong but limited solubility."},
            "NH3":      {"pKb": 4.75,   "type": "weak",       "notes": "Ammonia; lone pair on sp3 N."},
            "methylamine": {"pKb": 3.36, "type": "weak",      "notes": "Electron-donating methyl increases basicity vs NH3."},
            "ethylamine":  {"pKb": 3.25, "type": "weak",      "notes": "Similar to methylamine."},
            "aniline":   {"pKb": 9.38,   "type": "weak",       "notes": "Lone pair delocalized into benzene ring → much weaker."},
            "pyridine":  {"pKb": 8.75,   "type": "weak",       "notes": "sp2 N; lone pair in orbital with more s-character."},
            "trimethylamine": {"pKb": 4.20, "type": "weak",    "notes": "Three EDG groups, but steric hindrance in water."},
            "hydrazine": {"pKb": 5.77,   "type": "weak",       "notes": "N2H4; less basic than NH3 due to electron-withdrawing N."},
            "carbonate": {"pKb": 3.67,   "type": "weak",       "notes": "CO3^2- as base; Kb from Kw/Ka2 of H2CO3."},
            "bicarbonate": {"pKb": 7.65, "type": "weak",       "notes": "HCO3- as base; Kb from Kw/Ka1 of H2CO3."},
        }

        self._aliases = {
            "acetic acid": "CH3COOH", "acetate": "CH3COOH", "hac": "CH3COOH",
            "formic acid": "HCOOH", "methanoic acid": "HCOOH",
            "benzoic acid": "C6H5COOH",
            "hydrofluoric acid": "HF",
            "hydrocyanic acid": "HCN", "prussic acid": "HCN",
            "phosphoric acid": "H3PO4",
            "carbonic acid": "H2CO3",
            "hydrosulfuric acid": "H2S",
            "ammonia": "NH3",
            "sodium hydroxide": "NaOH", "lye": "NaOH",
            "potassium hydroxide": "KOH",
            "methyl amine": "methylamine", "ch3nh2": "methylamine",
            "phenol": "phenol", "c6h5oh": "phenol",
            "boric acid": "H3BO3",
            "perchloric acid": "HClO4",
            "nitric acid": "HNO3",
            "sulfuric acid": "H2SO4",
        }

    def _run_base(self, species_list: List[str], compare_type: str = "auto") -> dict:
        """Core logic: compare strengths."""
        if not species_list:
            raise ChemMCPError("Species list cannot be empty.")

        ctype = compare_type.lower() if compare_type else "auto"

        # Resolve names
        resolved = []
        for s in species_list:
            r = self._resolve(s)
            resolved.append((s, r))

        # Auto-detect type
        if ctype == "auto":
            acid_count = sum(1 for _, r in resolved if r in self._acids)
            base_count = sum(1 for _, r in resolved if r in self._bases)
            if acid_count >= base_count:
                ctype = "acid"
            else:
                ctype = "base"

        # Build ranking
        if ctype == "acid":
            ranking_data = []
            for original, r in resolved:
                if r in self._acids:
                    info = self._acids[r].copy()
                    info["name"] = original
                    ranking_data.append(info)
                elif r in self._bases:
                    # Convert base to its conjugate acid
                    binfo = self._bases[r]
                    conj_pka = 14.0 - binfo["pKb"] if binfo["pKb"] > 0 else None
                    ranking_data.append({
                        "name": original,
                        "pKa": conj_pka,
                        "type": "conjugate_acid_of_base",
                        "notes": f"Conjugate acid of base with pKb={binfo['pKb']}",
                    })
                else:
                    ranking_data.append({
                        "name": original,
                        "pKa": None,
                        "type": "unknown",
                        "notes": f"No data found for '{original}'",
                    })

            # Sort by pKa ascending (lower pKa = stronger acid); unknowns at end
            known = [x for x in ranking_data if x.get("pKa") is not None]
            unknown = [x for x in ranking_data if x.get("pKa") is None]
            known.sort(key=lambda x: x["pKa"])
            ranking = []
            for item in known:
                pka = item.get("pKa", 0)
                if pka < 0:
                    strength = "strong"
                elif pka < 4:
                    strength = "moderately weak"
                elif pka < 8:
                    strength = "weak"
                else:
                    strength = "very weak"
                ranking.append({
                    "name": item["name"],
                    "pKa": round(pka, 2),
                    "strength": strength,
                })
            for item in unknown:
                ranking.append({"name": item["name"], "pKa": None, "strength": "unknown"})

            analysis = self._analyze_acid_structure(ranking)

        else:  # base
            ranking_data = []
            for original, r in resolved:
                if r in self._bases:
                    info = self._bases[r].copy()
                    info["name"] = original
                    ranking_data.append(info)
                elif r in self._acids:
                    ainfo = self._acids[r]
                    conj_pkb = 14.0 - ainfo["pKa"] if ainfo["pKa"] > 0 else None
                    ranking_data.append({
                        "name": original,
                        "pKb": conj_pkb,
                        "type": "conjugate_base_of_acid",
                        "notes": f"Conjugate base of acid with pKa={ainfo['pKa']}",
                    })
                else:
                    ranking_data.append({
                        "name": original,
                        "pKb": None,
                        "type": "unknown",
                        "notes": f"No data found for '{original}'",
                    })

            known = [x for x in ranking_data if x.get("pKb") is not None]
            unknown = [x for x in ranking_data if x.get("pKb") is None]
            known.sort(key=lambda x: x["pKb"])
            ranking = []
            for item in known:
                pkb = item.get("pKb", 0)
                if pkb < 0:
                    strength = "strong"
                elif pkb < 4:
                    strength = "moderately weak"
                elif pkb < 8:
                    strength = "weak"
                else:
                    strength = "very weak"
                ranking.append({
                    "name": item["name"],
                    "pKb": round(pkb, 2),
                    "strength": strength,
                })
            for item in unknown:
                ranking.append({"name": item["name"], "pKb": None, "strength": "unknown"})

            analysis = self._analyze_base_structure(ranking)

        logger.info(f"AcidBaseStrengthCompare: compared {len(species)} species as {ctype}")
        return {
            "ranking": ranking,
            "structural_analysis": analysis,
            "compare_type_used": ctype,
        }

    def _analyze_acid_structure(self, ranking: list) -> str:
        """Generate structural explanation for acid strength differences."""
        parts = []
        for i, item in enumerate(ranking):
            name = item["name"]
            pka = item.get("pKa")
            if pka is None:
                continue
            reasons = []
            if name in ("HF",):
                reasons.append("strong H-F bond (high bond energy)")
                reasons.append("hydrogen bonding stabilizes undissociated HF")
            elif name in ("HCl", "HBr", "HI"):
                reasons.append(f"bond strength decreases down group (H-I < H-Br < H-Cl)")
            elif name == "CH3COOH":
                reasons.append("resonance stabilization of acetate conjugate base")
                reasons.append("electron-donating methyl slightly destabilizes conjugate base")
            elif name == "HCOOH":
                reasons.append("resonance stabilization of formate; no EDG to weaken it")
            elif name == "phenol":
                reasons.append("phenoxide ion stabilized by resonance delocalization into aromatic ring")
            elif name in ("H3PO4", "H2CO3", "H2S"):
                reasons.append("polyprotic oxoacid; stability of oxoanion matters")
            elif name == "HCN":
                reasons.append("CN⁻ is very stable (resonance), making HCN reluctant to donate H⁺")
            if reasons:
                parts.append(f"  • {name} (pKa={pka}): {'; '.join(reasons)}")

        if not parts:
            return "Ranking based on tabulated pKa values. Structural factors include bond polarity, conjugate base stability (resonance, inductive effects), solvation energy, and atomic size."
        return "Structural analysis:\n" + "\n".join(parts)

    def _analyze_base_structure(self, ranking: list) -> str:
        """Generate structural explanation for base strength differences."""
        parts = []
        for i, item in enumerate(ranking):
            name = item["name"]
            pkb = item.get("pKb")
            if pkb is None:
                continue
            reasons = []
            if name in ("NaOH", "KOH", "Ba(OH)2", "Ca(OH)2"):
                reasons.append("ionic hydroxide — complete dissociation in water")
            elif name == "NH3":
                reasons.append("lone pair on sp³ nitrogen available for protonation")
            elif name in ("methylamine", "ethylamine"):
                reasons.append("+I (inductive) effect of alkyl group increases e⁻ density on N")
            elif name == "aniline":
                reasons.append("lone pair delocalized into π-system of benzene ring (resonance)")
                reasons.append("much weaker than aliphatic amines despite having N atom")
            elif name == "pyridine":
                reasons.append("lone pair in sp² orbital (33% s-char) — held more tightly")
                reasons.append("less available for protonation vs sp³ N in NH3")
            elif name == "trimethylamine":
                reasons.append("three +I groups enhance basicity, but steric inhibition of solvation in water")
            elif name in ("carbonate", "bicarbonate"):
                reasons.append("anionic base — basicity derived from hydrolysis")

            if reasons:
                parts.append(f"  • {name} (pKb={pkb}): {'; '.join(reasons)}")

        if not parts:
            return "Ranking based on tabulated pKb values. Key factors: availability of lone pair, electron-donating/withdrawing substituents, resonance delocalization, hybridization (s-character), solvation effects."
        return "Structural analysis:\n" + "\n".join(parts)

    def _run_text(self, input_str: str) -> dict:
        """Parse semicolon-separated input."""
        try:
            parts = input_str.strip().split(";")
            species_list = [s.strip() for s in parts[0].split(";") if s.strip()] if ";" in input_str else [s.strip() for s in parts if s.strip()]
            # Re-parse properly
            tokens = input_str.strip().split()
            species_list = []
            compare_type = "auto"
            for tok in tokens:
                if tok.lower() in ("acid", "base", "auto"):
                    compare_type = tok.lower()
                else:
                    species_list.append(tok.replace(";", ""))
            if not species_list:
                raise ValueError("No species provided.")
            return self._run_base(species_list, compare_type)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'sp1;sp2;... [acid|base|auto]'")

    def _resolve(self, name: str) -> str:
        """Resolve name/alias to canonical key."""
        n = name.strip()
        if n in self._acids or n in self._bases:
            return n
        nl = n.lower()
        if nl in self._aliases:
            return self._aliases[nl]
        for k in list(self._acids.keys()) + list(self._bases.keys()):
            if k.lower() == nl:
                return k
        return n
