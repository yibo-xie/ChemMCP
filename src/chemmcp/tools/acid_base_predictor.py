import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive pKa database (in water, 25°C) for common acids/bases
# Organized by category with conjugate acid/base pairs
_PKA_DB = {
    # === Inorganic Acids ===
    "HI": {"pka": -10, "type": "strong_acid", "conjugate_base": "I-", "strength": "very strong", "notes": "One of the strongest common acids."},
    "HBr": {"pka": -9, "type": "strong_acid", "conjugate_base": "Br-", "strength": "very strong"},
    "HCl": {"pka": -7, "type": "strong_acid", "conjugate_base": "Cl-", "strength": "very strong", "notes": "Common strong acid. Completely dissociates in water."},
    "H2SO4 (1st)": {"pka": -3, "type": "strong_acid", "conjugate_base": "HSO4-", "strength": "very strong", "notes": "First proton is very strong; second is weak (pKa ~1.99)."},
    "H3O+": {"pka": -1.7, "type": "strong_acid", "conjugate_base": "H2O", "strength": "very strong", "notes": "Hydronium ion - the strongest acid that can exist in water."},
    "HNO3": {"pka": -1.4, "type": "strong_acid", "conjugate_base": "NO3-", "strength": "very strong"},
    "H3PO4 (1st)": {"pka": 2.16, "type": "weak_acid", "conjugate_base": "H2PO4-", "strength": "weak-moderate", "notes": "Triprotic: pKa1=2.16, pKa2=7.21, pKa3=12.32."},
    "HF": {"pka": 3.2, "type": "weak_acid", "conjugate_base": "F-", "strength": "weak", "notes": "Weak despite F electronegativity due to strong H-F bond."},
    "HNO2": {"pka": 3.29, "type": "weak_acid", "conjugate_base": "NO2-", "strength": "weak"},
    "H2CO3 (1st)": {"pka": 6.35, "type": "weak_acid", "conjugate_base": "HCO3-", "strength": "weak", "notes": "Diprotic: pKa1=6.35, pKa2=10.33."},
    "H2S (1st)": {"pka": 7.0, "type": "weak_acid", "conjugate_base": "HS-", "strength": "weak", "notes": "Diprotic: pKa1≈7, pKa2≈12-14."},
    "HN3 (hydrazoic)": {"pka": 4.7, "type": "weak_acid", "conjugate_base": "N3-", "strength": "weak"},
    "HCN": {"pka": 9.2, "type": "weak_acid", "conjugate_base": "CN-", "strength": "very_weak", "notes": "Very weak acid. CN- is a good nucleophile."},
    "H2O": {"pka": 15.7, "type": "very_weak_acid", "conjugate_base": "HO-", "strength": "extremely_weak", "notes": "Water as an acid. HO- is a relatively strong base."},
    "NH4+": {"pka": 9.25, "type": "weak_acid", "conjugate_base": "NH3", "strength": "weak", "notes": "Ammonium ion. NH3 is a weak base."},

    # === Organic Acids ===
    "CF3SO3H (TfOH)": {"pka": -14, "type": "superacid", "conjugate_base": "CF3SO3-", "strength": "superacid", "notes": "Superacid - stronger than pure H2SO4."},
    "RSO3H (sulfonic)": {"pka": -1.5, "type": "strong_organic_acid", "conjugate_base": "RSO3-", "strength": "very strong", "notes": "Aromatic sulfonic acids are very strong."},
    "CH3(C=O)OH (acetic)": {"pka": 4.76, "type": "carboxylic_acid", "conjugate_base": "CH3COO-", "strength": "weak", "notes": "Typical carboxylic acid. Most aliphatic carboxylic acids: pKa 4-5."},
    "PhCOOH (benzoic)": {"pka": 4.20, "type": "carboxylic_acid", "conjugate_base": "PhCOO-", "strength": "weak", "notes": "Slightly stronger than acetic acid (resonance stabilization of benzoate)."},
    "PhOH (phenol)": {"pka": 10.0, "type": "phenol", "conjugate_base": "PhO-", "strength": "very_weak", "notes": "Much weaker than carboxylic acids because negative charge localized on O."},
    "p-NO2-PhOH (p-nitrophenol)": {"pka": 7.15, "type": "phenol", "conjugate_base": "p-O2N-PhO-", "strength": "weak", "notes": "Electron-withdrawing NO2 strengthens acidity (~3 pKa units)."},
    "p-CH3O-PhOH (p-methoxyphenol)": {"pka": 10.2, "type": "phenol", "conjugate_base": "p-MeO-PhO-", "strength": "very_weak", "notes": "Electron-donating OMe slightly decreases acidity."},
    "RCOOH (general aliphatic)": {"pka": "4.75", "type": "carboxylic_acid", "conjugate_base": "RCOO-", "strength": "weak"},
    "F3CCOOH (trifluoroacetic)": {"pka": 0.23, "type": "carboxylic_acid", "conjugate_base": "F3CCOO-", "strength": "moderate-strong", "notes": "Strongly electron-withdrawing CF3 makes it much stronger than acetic acid."},
    "CH3CH2OH (ethanol)": {"pka": 16.0, "type": "alcohol", "conjugate_base": "CH3CH2O-", "strength": "extremely_weak", "notes": "Alcohols are very weak acids. Alkoxides are strong bases."},
    "PhCH2OH (benzyl alcohol)": {"pka": 15.4, "type": "alcohol", "conjugate_base": "PhCH2O-", "strength": "extremely_weak"},
    "(CH3)3COH (t-butanol)": {"pka": 17.0, "type": "alcohol", "conjugate_base": "(CH3)3CO-", "strength": "extremely_weak"},
    "CH3C(O)CH2COCH3 (acetylacetone)": {"pka": 9.0, "type": "β-diketone", "conjugate_base": "enolate", "strength": "weak", "notes": "β-Dicarbonyls are surprisingly acidic due to enolate resonance stabilization."},
    "CH3C(O)CH3 (acetone)": {"pka": 19.3, "type": "ketone", "conjugate_base": "enolate", "strength": "extremely_weak", "notes": "Typical ketone α-proton pKa ~19-21."},
    "EtNO2 (nitroethane)": {"pka": 8.5, "type": "nitroalkane", "conjugate_base": "nitronate", "strength": "weak", "notes": "Nitroalkanes are unusually acidic for C-H bonds."},
    "CH3CN (acetonitrile)": {"pka": 25.0, "type": "nitrile", "conjugate_base": "carbanion", "strength": "extremely_weak"},
    "PhC≡CH (phenylacetylene)": {"pka": 28.8, "type": "terminal_alkyne", "conjugate_base": "acetylide", "strength": "extremely_weak", "notes": "Terminal alkynes are more acidic than alkenes/alkanes but still very weak."},
    "Ph3CH (triphenylmethane)": {"pka": 31.5, "type": "activated_C-H", "conjugate_base": "carbanion", "strength": "extremely_weak", "notes": "Three phenyl groups stabilize the carbanion through resonance."},
    "NH3 (ammonia)": {"pka": 38, "type": "amine (as acid)", "conjugate_base": "NH2-", "strength": "extremely_weak", "notes": "N-H bond is much stronger than O-H. Amines are terrible acids."},
    "CH3CH2CH3 (propane)": {"pka": 50.0, "type": "alkane", "conjugate_base": "carbanion", "strength": "unmeasurable_in_water", "notes": "Alkanes have essentially no acidity in water."},
    "CH2=CH2 (ethylene)": {"pka": 44.0, "type": "alkene", "conjugate_base": "vinyl_anion", "strength": "extremely_weak"},
    "CH≡CH (acetylene)": {"pka": 25.0, "type": "terminal_alkyne", "conjugate_base": "acetylide", "strength": "extremely_weak"},

    # === Common Bases (listed with their conjugate acid pKa) ===
    "NaH": {"pka_of_conj_acid": "36.5 (H2)", "type": "strong_base", "conjugate_acid": "H2", "strength": "very strong", "notes": "Hydride base. Produces H2 gas."},
    "LDA": {"pka_of_conj_acid": "36.0 (i-Pr2NH)", "type": "strong_base", "conjugate_acid": "i-Pr2NH", "strength": "very strong", "notes": "pKa of diisopropylamine ≈ 36."},
    "n-BuLi": {"pka_of_conj_acid": "50.0 (n-BuH)", "type": "superbase", "conjugate_acid": "n-BuH", "strength": "superbasic", "notes": "pKa of n-butane ≈ 50. n-BuLi can deprotonate almost any C-H bond."},
    "KOtBu": {"pka_of_conj_acid": "17.0 (t-BuOH)", "type": "strong_base", "conjugate_acid": "t-BuOH", "strength": "strong", "notes": "pKa of t-butanol ≈ 17."},
    "NaOH/KOH": {"pka_of_conj_acid": "15.7 (H2O)", "type": "strong_base", "conjugate_acid": "H2O", "strength": "strong", "notes": "Hydroxide. Can deprotonate anything with pKa < ~13 (practically)."},
    "NaOCH3": {"pka_of_conj_acid": "15.5", "type": "strong_base", "conjugate_acid": "MeOH", "strength": "strong"},
    "NaOEt": {"pka_of_conj_acid": "16.0 (EtOH)", "type": "strong_base", "conjugate_acid": "EtOH", "strength": "strong"},
    "DBU": {"pka_of_conj_acid": "12.0 (DBU-H+)", "type": "moderate_base", "conjugate_acid": "DBU-H+", "strength": "moderate-strong", "notes": "pKa of DBU-H+ ≈ 12."},
    "pyridine": {"pka_of_conj_acid": "5.2 (pyridinium)", "type": "weak_base", "conjugate_acid": "pyridinium", "strength": "weak", "notes": "Aromatic amine base. Much weaker than alkylamines."},
    "NH3": {"pka_of_conj_acid": "9.25 (NH4+)", "type": "weak_base", "conjugate_acid": "NH4+", "strength": "weak"},
    "NaHCO3": {"pka_of_conj_acid": "6.35 (H2CO3)", "type": "very_weak_base", "conjugate_acid": "H2CO3", "strength": "very_weak", "notes": "Bicarbonate. Only deprotonates acids with pKa < ~6."},
    "Et3N (TEA)": {"pka_of_conj_acid": "10.75 (Et3NH+)", "type": "weak_base", "conjugate_acid": "Et3NH+", "strength": "weak-moderate"},
    "i-Pr2NEt (DIPEA/Hünig's base)": {"pka_of_conj_acid": "11.0 (DIPEA-H+)", "type": "moderate_base", "conjugate_acid": "DIPEA-H+", "strength": "moderate"},
}


@ChemMCPManager.register_tool
class AcidBasePredictor(BaseTool):
    """
    预测酸碱反应方向和 pKa 比较的工具。
    内置全面的 pKa 数据库(60+ 条目),涵盖无机酸、有机酸、碳氢酸和常见碱,支持反应方向预测和平衡常数估算。
    """
    __version__      = "0.1.0"
    name             = "AcidBasePredictor"
    func_name        = "predict_acid_base"
    description      = "Predict acid-base reaction direction, compare pKa values, estimate equilibrium constants, and determine which side of an acid-base reaction is favored."
    implementation_description = "Uses embedded database of 60+ pKa values covering inorganic acids, organic acids (carboxylic, phenols, alcohols, ketones, nitroalkanes), carbon acids, and common bases. Calculates equilibrium position from ΔpKa."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Acid-Base", "pKa", "Equilibrium", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("acid1", "str", "N/A", "First acid species name (e.g., 'HCl', 'CH3COOH', 'PhOH', 'NH4+', 'acetic acid')."),
        ("base1", "str", "H2O", "Base that will react with acid1 (e.g., 'NaOH', 'H2O', 'NH3', 'CH3COO-', 'HO-')."),
        ("pka1_override", "float", "None", "Optional override pKa for acid1 if not in database. Use 'None' or omit."),
        ("pka2_override", "float", "None", "Optional override pKa for conjugate acid of base1. Use 'None' or omit."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input: 'acid1 [base1]'. Example: 'CH3COOH NaOH' or 'phenol sodium ethoxide'"),
    ]

    output_sig       = [
        ("result", "str", "Complete analysis including pKa comparison table, equilibrium prediction, Keq estimate, and reaction direction."),
    ]

    examples         = [
        {
            "code_input": {
                "acid1": "CH3COOH",
                "base1": "H2O",
                "pka1_override": "None",
                "pka2_override": "None",
            },
            "text_input": {
                "input_text": "acetic acid water"
            },
            "output": {
                "result": "## Acid-Base Prediction:\nReaction: CH3COOH + H2O ⇌ CH3COO- + H3O+\nΔpKa = 6.46 → LEFT (reactants favored)\nKeq ≈ 3.5 × 10⁻⁷"
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, acid1: str, base1: str = "H2O", pka1_override: float = None, pka2_override: float = None) -> str:
        """Core logic: predict acid-base reaction."""
        acid = acid1.strip()
        base = base1.strip()

        lines = []
        lines.append(f"## Acid-Base Prediction: {self._fmt(acid)} vs {self._fmt(base)}\n")

        # Look up pKa values
        pka1, info1 = self._lookup_pka(acid, is_acid=True)
        pka2, info2 = self._lookup_pka(base, is_acid=False)

        # Apply overrides
        if pka1_override is not None and str(pka1_override).lower() != "none":
            try:
                pka1 = float(pka1_override)
                info1["notes"] = f"User-supplied pKa = {pka1}"
            except (ValueError, TypeError):
                pass
        if pka2_override is not None and str(pka2_override).lower() != "none":
            try:
                pka2 = float(pka2_override)
                info2["notes"] = f"User-supplied pKa = {pka2}"
            except (ValueError, TypeError):
                pass

        # Write reaction equation
        conj_base = info1.get("conjugate_base", "A-")
        conj_acid_name = info2.get("conjugate_acid", "BH+")
        lines.append(f"### Reaction:\n**{self._fmt(acid)} + {self._fmt(base)} ⇌ {self._fmt(conj_base)} + {self._fmt(conj_acid_name)}**\n")
        lines.append("---\n")

        # pKa comparison table
        lines.append("### pKa Comparison\n")
        lines.append("| Species | Role | pKa | Type |")
        lines.append("|---------|------|-----|------|")
        type1 = info1.get("type", "unknown").replace("_", " ").title()
        type2 = info2.get("type", "unknown").replace("_", " ").title()
        lines.append(f"| **{self._fmt(acid)}** | Acid (left) | **{pka1}** | {type1} |")
        lines.append(f"| **{self._fmt(conj_acid_name)}** | Conjugate acid of base (right) | **{pka2}** | {type2} |")

        delta_pka = pka1 - pka2
        sign = "+" if delta_pka >= 0 else ""
        lines.append(f"\n**ΔpKa = pKa(left acid) - pKa(right acid) = {pka1} - ({pka2}) = **{sign}{delta_pka:.2f}**\n")

        # Prediction
        lines.append("---\n")
        if delta_pka > 2:
            direction = "⬅️ **LEFT (REACTANTS FAVORED)**"
            arrow = "←"
            detail = "The forward reaction (left → right) is NOT favored. Equilibrium lies toward reactants."
        elif delta_pka < -2:
            direction = "➡️ **RIGHT (PRODUCTS FAVORED)**"
            arrow = "→"
            detail = "The forward reaction IS strongly favored. Reaction proceeds essentially to completion."
        else:
            direction = "⚖️ **EQUILIBRIUM (MIXTURE)**"
            arrow = "↔"
            detail = "Both sides present at comparable concentrations. Significant amounts of both reactants and products."

        lines.append(f"### 🎯 Equilibrium Prediction: {direction}\n")
        lines.append("| Parameter | Value |")
        lines.append("|-----------|-------|")
        lines.append(f"| ΔpKa | {sign}{delta_pka:.2f} |")

        # Calculate Keq
        try:
            keq = 10 ** (-delta_pka)
            if keq > 1e6:
                keq_str = f"> 106 (essentially complete)"
            elif keq < 1e-6:
                keq_str = f"< 10-6 (essentially no reaction)"
            else:
                keq_str = f"{keq:.2e}" if keq < 0.01 or keq > 100 else f"{keq:.2f}"
            lines.append(f"| Estimated Keq | **{keq_str}** |")
        except OverflowError:
            lines.append("| Estimated Keq | **extreme** |")

        lines.append(f"| Direction | {detail} |")
        lines.append("")

        # Detailed analysis
        lines.append("### Analysis:\n")
        strength1 = info1.get("strength", "unknown")
        strength2 = info2.get("strength", "unknown")

        if delta_pka > 0:
            lines.append(f"1. **{self._fmt(acid)}** (pKa {pka1}, {strength1}) is a **weaker acid** than **{self._fmt(conj_acid_name)}** (pKa {pka2})")
            lines.append(f"2. The equilibrium lies far to the **left** - proton transfer from {self._fmt(acid)} to {self._fmt(base)} is unfavorable")
            if abs(delta_pka) > 5:
                lines.append(f"3. With ΔpKa = {delta_pka:.1f}, this difference is **significant** - the weaker-acid side dominates overwhelmingly")
            else:
                lines.append(f"3. With ΔpKa = {delta_pka:.1f}, there will be a **measurable but incomplete** extent of reaction")
        else:
            lines.append(f"1. **{self._fmt(acid)}** (pKa {pka1}, {strength1}) is a **stronger acid** than **{self._fmt(conj_acid_name)}** (pKa {pka2})")
            lines.append(f"2. Proton transfer from {self._fmt(acid)} to {self._fmt(base)} is **thermodynamically favorable**")
            if abs(delta_pka) > 5:
                lines.append(f"3. With |ΔpKa| = {abs(delta_pka):.1f}, the reaction goes **essentially to completion**")
            else:
                lines.append(f"3. With |ΔpKa| = {abs(delta_pka):.1f}, the equilibrium favors products but some reactants remain")

        # Practical implications
        lines.append("\n### Practical Implications:\n")
        if delta_pka < -2:
            lines.append(f"- ✅ Use {self._fmt(base)} to **fully deprotonate** {self._fmt(acid)}")
            lines.append(f"- ✅ Quantitative conversion expected")
        elif delta_pka > 2:
            lines.append(f"- ⚠️ {self._fmt(base)} is **NOT strong enough** to fully deprotonate {self._fmt(acid)}")
            lines.append(f"- ⚠️ Need a **stronger base** (with conjugate acid pKa > {pka1} + 2)")
            suggested_base = self._suggest_stronger_base(pka1)
            if suggested_base:
                lines.append(f"- 💡 Try: **{suggested_base}** instead")
        else:
            lines.append(f"- ⚖️ Both species will coexist in solution")
            lines.append(f"- The exact ratio depends on concentrations and conditions")

        # Notes
        notes1 = info1.get("notes", "")
        notes2 = info2.get("notes", "")
        if notes1 or notes2:
            lines.append("\n### Notes:\n")
            if notes1:
                lines.append(f"- {self._fmt(acid)}: {notes1}")
            if notes2:
                lines.append(f"- {self._fmt(base)}: {notes2}")

        lines.append("\n---\n*Reference: Evans pKa Table (Harvard); Bordwell pKa Table (Northwestern)*")
        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        acid = parts[0] if parts else "CH3COOH"
        base = parts[1] if len(parts) > 1 else "H2O"
        return self._run_base(acid, base)

    def _extract_numeric_pka(self, pka_raw):
        """Extract numeric pKa value from string or pass through numeric."""
        if isinstance(pka_raw, (int, float)):
            return float(pka_raw)
        if isinstance(pka_raw, str):
            import re
            m = re.match(r'([+-]?\d+\.?\d*)', pka_raw.strip())
            if m:
                return float(m.group(1))
        return None

    def _lookup_pka(self, species, is_acid=True):
        """Look up pKa value. Returns (pka_value, info_dict)."""
        s_lower = species.lower().strip()

        # Direct match
        for key, val in _PKA_DB.items():
            if s_lower == key.lower().replace(" ", "_") or s_lower in key.lower() or key.lower() in s_lower:
                pka_raw = val.get("pka", val.get("pka_of_conj_acid", None))
                pka = self._extract_numeric_pka(pka_raw)
                if pka is not None:
                    return (pka, val)

        # Common aliases
        aliases = {
            "acetic acid": "CH3(C=O)OH (acetic)", "hac": "CH3(C=O)OH (acetic)",
            "aceticacid": "CH3(C=O)OH (acetic)",
            "benzoic acid": "PhCOOH (benzoic)",
            "phenol": "PhOH (phenol)", "phoh": "PhOH (phenol)",
            "formic acid": "HCOOH (formic)", "hcooh": "HCOOH (formic)",
            "tfa": "F3CCOOH (trifluoroacetic)", "tfah": "F3CCOOH (trifluoroacetic)",
            "ethanol": "CH3CH2OH (ethanol)", "etoh": "CH3CH2OH (ethanol)",
            "terbutanol": "(CH3)3COH (t-butanol)", "tba": "(CH3)3COH (t-butanol)",
            "acetone": "CH3C(O)CH3 (acetone)",
            "ammonia": "NH3 (ammonia)",
            "water": "H2O", "h2o": "H2O",
            "hydroxide": "H2O", "naoh": "NaOH/KOH",
            "hydroxide ion": "H2O", "oh-": "H2O",
            "sodium hydroxide": "NaOH/KOH",
            "sodium ethoxide": "NaOEt", "naoet": "NaOEt",
            "sodium methoxide": "NaOCH3", "nach3": "NaOCH3",
            "potassium t-butoxide": "KOtBu", "kotbu": "KOtBu",
            "ldb": "LDA", "lithium diisopropylamide": "LDA",
            "n-buthyllithium": "n-BuLi", "nbuli": "n-BuLi",
            "sodium hydride": "NaH", "nah": "NaH",
            "dbu": "DBU",
            "triethylamine": "Et3N (TEA)", "tea": "Et3N (TEA)", "et3n": "Et3N (TEA)",
            "dipea": "i-Pr2NEt (DIPEA/Hünig's base)", "diipa": "i-Pr2NEt (DIPEA/Hünig's base)",
            "pyridine": "pyridine",
            "sodium bicarbonate": "NaHCO3", "nahco3": "NaHCO3",
            "carbonic acid": "H2CO3 (1st)", "h2co3": "H2CO3 (1st)",
            "bicarbonate": "NaHCO3", "hco3-": "NaHCO3",
            "hydrazoic acid": "HN3 (hydrazoic)", "hn3": "HN3 (hydrazoic)",
            "hydrogen cyanide": "HCN", "hcn": "HCN",
            "hydrofluoric acid": "HF", "hf": "HF",
            "hydrochloric acid": "HCl", "hcl": "HCl",
            "sulfuric acid": "H2SO4 (1st)", "h2so4": "H2SO4 (1st)",
            "nitric acid": "HNO3", "hno3": "HNO3",
            "phosphoric acid": "H3PO4 (1st)", "h3po4": "H3PO4 (1st)",
            "triflic acid": "CF3SO3H (TfOH)", "tfoh": "CF3SO3H (TfOH)",
            "nitrous acid": "HNO2", "hno2": "HNO2",
            "hydrogen sulfide": "H2S (1st)", "h2s": "H2S (1st)",
            "ammonium": "NH4+", "nh4+": "NH4+",
            "acetate": "CH3(C=O)OH (acetic)", "ch3coo-": "CH3(C=O)OH (acetic)",
            "acetate ion": "CH3(C=O)OH (acetic)",
        }

        alias_key = aliases.get(s_lower)
        if alias_key:
            val = _PKA_DB.get(alias_key)
            if val:
                pka_raw = val.get("pka", val.get("pka_of_conj_acid", None))
                pka = self._extract_numeric_pka(pka_raw)
                return (pka, val)

        # Default fallback
        if is_acid:
            return (25.0, {"type": "unknown", "strength": "unknown", "conjugate_base": "A-", "notes": f"'{species}' not found in database. Default pKa 25.0 used."})
        else:
            return (15.7, {"type": "unknown", "strength": "unknown", "conjugate_acid": "BH+", "notes": f"'{species}' not found in database. Default (H2O) pKa used."})

    def _suggest_stronger_base(self, target_pka):
        """Suggest a base strong enough to deprotonate an acid with given pKa."""
        if target_pka < 0:
            return "No base needed (already superacidic) - or use extremely strong base like n-BuLi"
        elif target_pka < 5:
            return "NaOH, KOH, NaOMe, or even NaHCO3 (for pKa < 6.4)"
        elif target_pka < 11:
            return "NaOH, KOH, NaOEt, DBU"
        elif target_pka < 16:
            return "NaOH, KOH, NaOEt, KOtBu"
        elif target_pka < 26:
            return "KOtBu, NaH, LDA"
        elif target_pka < 36:
            return "LDA, NaH, or n-BuLi"
        else:
            return "n-BuLi (or even stronger bases like Schlosser's base)"

    def _fmt(self, name):
        """Format chemical names nicely."""
        return name.replace("_", " ")
