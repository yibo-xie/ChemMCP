import logging
import re
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PredictProducts(BaseTool):
    """
    Predict reaction products given reactants.
    Uses reaction type classification and chemical knowledge to predict likely products.
    """
    __version__ = "0.1.0"
    name = "PredictProducts"
    func_name = "predict_products"
    description = "Predict the products of a chemical reaction given the reactants and optional conditions."
    implementation_description = "Uses pattern matching against common reaction types (combustion, acid-base, precipitation, single/double displacement, redox) with built-in databases of known reactions and activity series rules."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Reaction Prediction", "Products", "Chemical Synthesis", "Pattern Matching"]
    required_envs = []

    code_input_sig = [
        ("reactants", "str", "N/A", "Reactant species separated by '+', e.g., 'NaCl + AgNO3' or 'CH4 + O2' or 'Zn + H2SO4'."),
        ("conditions", "str", "room temperature", "Optional reaction conditions: 'heat', 'catalyst', 'light', 'electrolysis', 'acidic', 'basic', 'dilute', 'concentrated', etc."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query string: 'reactants [conditions]', e.g., 'NaCl + AgNO3' or 'Cu + HNO3 concentrated'."),
    ]

    output_sig = [
        ("predicted_products", "str", "Predicted product formula(s)."),
        ("balanced_equation", "str", "The predicted balanced equation."),
        ("reaction_type", "str", "The type of reaction that produces these products."),
        ("confidence", "str", "Confidence level: 'high', 'medium', or 'low'."),
        ("notes", "str", "Additional notes about side reactions, conditions required, or caveats."),
    ]

    examples = [
        {
            "code_input": {
                "reactants": "H2 + O2",
                "conditions": "spark/ignition",
            },
            "text_input": {"query": "H2 + O2 spark"},
            "output": {
                "predicted_products": "H2O",
                "balanced_equation": "2H2 + O2 = 2H2O",
                "reaction_type": "Combination / Combustion",
                "confidence": "high",
                "notes": "Requires ignition energy (activation barrier ~400 kJ/mol). Produces water; if limited O2, may also produce H2O2.",
            }
        },
        {
            "code_input": {
                "reactants": "NaCl + AgNO3",
                "conditions": "aqueous",
            },
            "text_input": {"query": "NaCl + AgNO3 aqueous"},
            "output": {
                "predicted_products": "AgCl + NaNO3",
                "balanced_equation": "NaCl + AgNO3 = AgCl↓ + NaNO3",
                "reaction_type": "Double Displacement / Precipitation",
                "confidence": "high",
                "notes": "AgCl is a white precipitate (Ksp ≈ 1.8×10⁻¹⁰). Reaction goes to completion in aqueous solution.",
            }
        },
        {
            "code_input": {
                "reactants": "Zn + H2SO4",
                "conditions": "dilute",
            },
            "text_input": {"query": "Zn + H2SO4 dilute"},
            "output": {
                "predicted_products": "ZnSO4 + H2",
                "balanced_equation": "Zn + H2SO4(dil) = ZnSO4 + H2↑",
                "reaction_type": "Single Replacement / Redox",
                "confidence": "high",
                "notes": "Zinc is above hydrogen in the activity series. With concentrated H2SO4, SO2 is produced instead of H2.",
            }
        },
        {
            "code_input": {
                "reactants": "C6H12O6",
                "conditions": "fermentation (yeast)",
            },
            "text_input": {"query": "C6H12O6 fermentation yeast"},
            "output": {
                "predicted_products": "C2H5OH + CO2",
                "balanced_equation": "C6H12O6 = 2C2H5OH + 2CO2",
                "reaction_type": "Decomposition / Fermentation",
                "confidence": "high",
                "notes": "Alcoholic fermentation by yeast. Anaerobic conditions. Temperature ~25-35°C optimal.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build reaction prediction database."""
        # Metal activity series (most reactive → least reactive)
        self._activity_series = [
            "K", "Na", "Ca", "Mg", "Al", "Zn", "Fe", "Ni", "Sn", "Pb",
            "H", "Cu", "Hg", "Ag", "Pt", "Au",
        ]
        self._activity_set = set(self._activity_series)

        # Known reaction patterns: (reactant_pattern_key) → prediction
        self._reaction_db = self._build_reaction_db()

    def _build_reaction_db(self) -> dict:
        """Build comprehensive reaction database."""
        db = {}

        # ── Combustion ──
        db["ch4_o2"] = {
            "products": ["CO2", "H2O"],
            "equation": "CH4 + 2O2 = CO2 + 2H2O",
            "type": "Combustion / Redox",
            "confidence": "high",
            "notes": "Complete combustion. Blue flame. ΔH° = -890 kJ/mol.",
        }
        db["c2h5oh_o2"] = {
            "products": ["CO2", "H2O"],
            "equation": "C2H5OH + 3O2 = 2CO2 + 3H2O",
            "type": "Combustion / Redox",
            "confidence": "high",
            "notes": "Ethanol combustion. Complete combustion.",
        }
        db["h2_o2"] = {
            "products": ["H2O"],
            "equation": "2H2 + O2 = 2H2O",
            "type": "Combination / Combustion",
            "confidence": "high",
            "notes": "Highly exothermic. Produces water.",
        }

        # ── Acid-Base Neutralization ──
        db["hcl_naoh"] = {
            "products": ["NaCl", "H2O"],
            "equation": "HCl + NaOH = NaCl + H2O",
            "type": "Acid-Base Neutralization",
            "confidence": "high",
            "notes": "Strong acid + strong base. pH neutral product.",
        }
        db["h2so4_naoh"] = {
            "products": ["Na2SO4", "H2O"],
            "equation": "H2SO4 + 2NaOH = Na2SO4 + 2H2O",
            "type": "Acid-Base Neutralization",
            "confidence": "high",
            "notes": "Diprotic acid neutralization.",
        }
        db["hcl_caco3"] = {
            "products": ["CaCl2", "H2O", "CO2"],
            "equation": "2HCl + CaCO3 = CaCl2 + H2O + CO2↑",
            "type": "Acid-Carbonate Reaction",
            "confidence": "high",
            "notes": "Effervescence due to CO2 gas. Carbonate dissolves in acid.",
        }
        db["naoh_co2"] = {
            "products": ["Na2CO3", "H2O"],
            "equation": "2NaOH + CO2 = Na2CO3 + H2O",
            "type": "Acid-Gas Reaction",
            "confidence": "high",
            "notes": "CO2 absorption by NaOH solution.",
        }

        # ── Precipitation / Double Displacement ──
        db["nacl_agno3"] = {
            "products": ["AgCl", "NaNO3"],
            "equation": "NaCl + AgNO3 = AgCl↓ + NaNO3",
            "type": "Precipitation / Double Displacement",
            "confidence": "high",
            "notes": "White precipitate of AgCl. Ksp = 1.8×10⁻¹⁰. Light-sensitive.",
        }
        db["bacl2_na2so4"] = {
            "products": ["BaSO4", "NaCl"],
            "equation": "BaCl2 + Na2SO4 = BaSO4↓ + 2NaCl",
            "type": "Precipitation / Double Displacement",
            "confidence": "high",
            "notes": "White precipitate of BaSO4. Used in X-ray imaging. Very insoluble (Ksp = 1.1×10⁻¹⁰).",
        }
        db["agno3_nacl"] = db["nacl_agno3"]

        # ── Single Replacement (Metal Activity Series) ──
        db["zn_h2so4_dilute"] = {
            "products": ["ZnSO4", "H2"],
            "equation": "Zn + H2SO4(dil) = ZnSO4 + H2↑",
            "type": "Single Replacement / Redox",
            "confidence": "high",
            "notes": "Zn displaces H from dilute acid. Bubbles of H2 observed.",
        }
        db["fe_hcl"] = {
            "products": ["FeCl2", "H2"],
            "equation": "Fe + 2HCl = FeCl2 + H2↑",
            "type": "Single Replacement / Redox",
            "confidence": "high",
            "notes": "Iron(II) chloride forms (not FeCl3 with dilute HCl).",
        }
        db["zn_cuso4"] = {
            "products": ["ZnSO4", "Cu"],
            "equation": "Zn + CuSO4 = ZnSO4 + Cu↓",
            "type": "Single Replacement / Redox",
            "confidence": "high",
            "notes": "Zn is more reactive than Cu. Blue solution fades as Cu deposits (red-brown solid).",
        }
        db["fe_cuso4"] = {
            "products": ["FeSO4", "Cu"],
            "equation": "Fe + CuSO4 = FeSO4 + Cu↓",
            "type": "Single Replacement / Redox",
            "confidence": "high",
            "notes": "Iron displaces copper. Color changes from blue to pale green.",
        }
        db["cl2_nabr"] = {
            "products": ["Br2", "NaCl"],
            "equation": "Cl2 + 2NaBr = 2NaCl + Br2",
            "type": "Halogen Single Replacement / Redox",
            "confidence": "high",
            "notes": "Cl₂ is more reactive than Br₂. Orange color appears (Br₂ in water).",
        }
        db["f2_h2o"] = {
            "products": ["HF", "O2"],
            "equation": "2F2 + 2H2O = 4HF + O2",
            "type": "Single Replacement / Redox (Violent)",
            "confidence": "high",
            "notes": "Extremely violent! Fluorine reacts explosively with water.",
        }

        # ── Thermal Decomposition ──
        db["kclo3_heat"] = {
            "products": ["KCl", "O2"],
            "equation": "2KClO3 = 2KCl + 3O2↑",
            "type": "Thermal Decomposition / Redox",
            "confidence": "high",
            "notes": "Requires MnO2 catalyst and heating (~200°C). Oxygen production for lab use.",
        }
        db["caco3_heat"] = {
            "products": ["CaO", "CO2"],
            "equation": "CaCO3 = CaO + CO2↑",
            "type": "Thermal Decomposition",
            "confidence": "high",
            "notes": "Limestone decomposition. Requires >840°C. Quicklime production.",
        }
        db["cu(oh)2_heat"] = {
            "products": ["CuO", "H2O"],
            "equation": "Cu(OH)2 = CuO + H2O",
            "type": "Thermal Decomposition",
            "confidence": "high",
            "notes": "Blue precipitate turns black upon heating.",
        }
        db["h2o2_decomp"] = {
            "products": ["H2O", "O2"],
            "equation": "2H2O2 = 2H2O + O2↑",
            "type": "Decomposition / Disproportionation",
            "confidence": "high",
            "notes": "Catalyzed by MnO2 or catalase enzyme. Spontaneous at room temp (slow).",
        }

        # ── Specific Redox Reactions ──
        db["mno4_fe2+_h+_acidic"] = {
            "products": ["Mn2+", "Fe3+", "H2O"],
            "equation": "MnO4^- + 5Fe^{2+} + 8H^+ = Mn^{2+} + 5Fe^{3+} + 4H2O",
            "type": "Redox (Permanganate Titration)",
            "confidence": "high",
            "notes": "Classic titration reaction. Purple MnO4⁻ becomes colorless Mn²⁺ as it reacts.",
        }
        db["cu_hno3_dilute"] = {
            "products": ["Cu(NO3)2", "NO", "H2O"],
            "equation": "3Cu + 8HNO3(dil) = 3Cu(NO3)2 + 2NO↑ + 4H2O",
            "type": "Redox / Oxidation by Nitric Acid",
            "confidence": "high",
            "notes": "Colorless NO gas (turns brown on air contact → NO2). Dilute HNO3.",
        }
        db["cu_hno3_concentrated"] = {
            "products": ["Cu(NO3)2", "NO2", "H2O"],
            "equation": "Cu + 4HNO3(conc) = Cu(NO3)2 + 2NO2↑ + 2H2O",
            "type": "Redox / Oxidation by Nitric Acid",
            "confidence": "high",
            "notes": "Brown NO2 gas evolved. Concentrated HNO3 acts as oxidizing agent.",
        }
        db["c_concentrated_h2so4"] = {
            "products": ["CO2", "SO2", "H2O"],
            "equation": "C + 2H2SO4(conc) = CO2↑ + 2SO2↑ + 2H2O",
            "type": "Redox / Hot Concentrated H2SO4 Oxidation",
            "confidence": "high",
            "notes": "Hot conc. H2SO4 oxidizes carbon. Both CO2 and SO2 gases produced.",
        }

        # ── Organic Reactions ──
        db["c6h12o6_fermentation"] = {
            "products": ["C2H5OH", "CO2"],
            "equation": "C6H12O6 = 2C2H5OH + 2CO2",
            "type": "Fermentation / Enzymatic Decomposition",
            "confidence": "high",
            "notes": "Alcoholic fermentation by yeast. Anaerobic. Optimal T = 25-35°C.",
        }
        db["ch3cooh_c2h5oh"] = {
            "products": ["CH3COOC2H5", "H2O"],
            "equation": "CH3COOH + C2H5OH ⇌ CH3COOC2H5 + H2O",
            "type": "Esterification (Condensation)",
            "confidence": "high",
            "notes": "Fischer esterification. Requires acid catalyst (H2SO4) and heat. Equilibrium reaction.",
        }

        # ── Synthesis / Combination ──
        db["nh3_hcl"] = {
            "products": ["NH4Cl"],
            "equation": "NH3 + HCl = NH4Cl",
            "type": "Combination (Acid-Base)",
            "confidence": "high",
            "notes": "Forms white solid ammonium chloride smoke. Exothermic.",
        }
        db["caO_h2o"] = {
            "products": ["Ca(OH)2"],
            "equation": "CaO + H2O = Ca(OH)2",
            "type": "Combination / Hydration",
            "confidence": "high",
            "notes": "Slaking of lime. Highly exothermic. Forms slaked lime.",
        }
        db["n2_h2"] = {
            "products": ["NH3"],
            "equation": "N2 + 3H2 ⇌ 2NH3",
            "type": "Combination / Haber Process",
            "confidence": "high",
            "notes": "Haber process conditions: 400-500°C, 150-250 atm, Fe catalyst. Equilibrium reaction.",
        }
        db["n2_o2"] = {
            "products": ["NO"],
            "equation": "N2 + O2 = 2NO",
            "type": "Combination (High-Temperature)",
            "confidence": "high",
            "notes": "Lightning or electric arc (>2000°C). Endothermic. NO formation in atmosphere.",
        }

        return db

    def _normalize_output(self, entry: dict) -> dict:
        """Normalize output dict to match output_sig keys.
        Internal DB uses 'products' but output_sig promises 'predicted_products'.
        Also normalizes 'equation' → 'balanced_equation', 'type' → 'reaction_type'."""
        result = dict(entry)
        # Rename keys to match output_sig
        if "products" in result and "predicted_products" not in result:
            result["predicted_products"] = ", ".join(result["products"]) if isinstance(result["products"], list) else result["products"]
        if "type" in result and "reaction_type" not in result:
            result["reaction_type"] = result.pop("type")
        if "equation" in result and "balanced_equation" not in result:
            result["balanced_equation"] = result.pop("equation")
        return result

    def _run_base(self, reactants: str, conditions: str = "room temperature") -> dict:
        """Predict products from reactants string."""
        r_str = reactants.strip()
        cond_str = (conditions or "room temperature").strip().lower()

        # Normalize reactants key
        norm_key = self._normalize_reactants(r_str)
        cond_key = cond_str.replace(" ", "_")

        # Try direct lookup with conditions
        lookup_key = f"{norm_key}_{cond_key}"
        if lookup_key in self._reaction_db:
            entry = self._reaction_db[lookup_key]
            logger.info(f"Direct match: {lookup_key}")
            return self._normalize_output(entry)

        # Try without conditions
        if norm_key in self._reaction_db:
            entry = self._reaction_db[norm_key]
            logger.info(f"Match (no condition): {norm_key}")
            return self._normalize_output(entry)

        # Fuzzy matching
        best_match, best_score = None, 0
        for key, entry in self._reaction_db.items():
            score = self._fuzzy_match_score(norm_key, cond_key, key)
            if score > best_score:
                best_match, best_score = (key, entry), score

        if best_match and best_score >= 0.5:
            logger.info(f"Fuzzy match ({best_score:.2f}): {best_match[0]}")
            result = dict(best_match[1])
            result["confidence"] = "medium" if best_score < 0.9 else result.get("confidence", "medium")
            result["notes"] += " [Fuzzy matched — please verify]"
            return self._normalize_output(result)

        # Rule-based fallback prediction
        return self._rule_based_predict(r_str, cond_str)

    def _run_text(self, query: str) -> dict:
        """Parse text query."""
        parts = query.strip().rsplit(None, 1)
        if len(parts) >= 2 and parts[-1].lower() in (
            "heat", "dilute", "concentrated", "aqueous", "spark", "light",
            "catalyst", "electrolysis", "acidic", "basic", "fermentation",
            "yeast", "ignition", "conc", "dry"
        ):
            return self._run_base(" ".join(parts[:-1]), parts[-1])
        return self._run_base(query)

    def _normalize_reactants(self, r_str: str) -> str:
        """Normalize reactants string to a lookup key."""
        r_str = r_str.lower().replace(" ", "").replace("+", "_")
        # Remove stoichiometric coefficients
        import re
        r_str = re.sub(r'^\d+', '', r_str)
        r_str = re.sub(r'_\d+', '_', r_str)
        return r_str.strip("_")

    def _fuzzy_match_score(self, norm_r: str, cond: str, db_key: str) -> float:
        """Calculate fuzzy match score between input and database key."""
        # Split db_key into reactants part and condition part
        if "_" in db_key:
            parts = db_key.rsplit("_", 1)
            db_r, db_cond = parts[0], parts[1]
        else:
            db_r, db_cond = db_key, ""

        # Compare reactants
        r_tokens = set(norm_r.split("_"))
        db_tokens = set(db_r.split("_"))

        overlap = len(r_tokens & db_tokens)
        total = max(len(r_tokens), len(db_tokens))
        r_score = overlap / total if total > 0 else 0

        # Condition bonus
        c_score = 1.0 if not cond or cond == db_cond else (0.5 if cond in db_cond or db_cond in cond else 0)

        return r_score * 0.8 + c_score * 0.2

    def _rule_based_predict(self, r_str: str, cond: str) -> dict:
        """Fallback rule-based prediction when no direct match found."""
        reactants = [r.strip() for r in r_str.split('+')]

        # Check for metal + acid
        for r in reactants:
            if r.strip() in self._activity_set:
                idx = self._activity_series.index(r.strip())
                if idx < self._activity_series.index("H"):
                    # Metal can displace H
                    for a in reactants:
                        al = a.strip().lower()
                        if any(acid in al for acid in ["hcl", "h2so4", "hno3"]):
                            salt = f"{r.strip()}2" if any(x in al for x in ["hcl"]) else \
                                    f"{r.strip()}SO4" if "h2so4" in al else \
                                    f"{r.strip()}(NO3)2" if "hno3" in al else f"{r.strip()}_salt"
                            return {
                                "predicted_products": f"{salt} + H2",
                                "balanced_equation": f"2{r.strip()} + 2{a.strip()} = {salt} + H2↑",
                                "reaction_type": "Single Replacement / Redox",
                                "confidence": "low",
                                "notes": f"Predicted based on activity series ({r.strip()} is above H). Please verify experimentally.",
                            }

        raise ChemMCPError(
            f"Cannot confidently predict products for: '{r_str}' (conditions: '{cond}'). "
            f"This tool supports common reaction patterns including combustion, acid-base, "
            f"precipitation, single/double displacement, thermal decomposition, and specific redox reactions. "
            f"For more complex or organic reactions, consider using ForwardSynthesis or Retrosynthesis tools."
        )
