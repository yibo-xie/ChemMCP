import logging
from typing import Dict, List, Optional, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GasEvolutionPrediction(BaseTool):
    """
    预测气体逸出反应。
    基于常见气体生成模式（酸+金属→H2、酸+碳酸盐→CO2、过氧化物分解→O2等）判断反应是否会产生气体。
    """
    __version__ = "0.1.0"
    name = "GasEvolutionPrediction"
    func_name = "predict_gas_evolution"
    description = "Predict whether a chemical reaction will produce gas (evolution). Identifies the gas product and provides the balanced equation."
    implementation_description = "Uses a built-in pattern database of common gas-producing reactions: acid + metal → H2, acid + carbonate/bicarbonate/sulfite/sulfide → CO2/SO2/H2S, peroxide decomposition → O2, ammonia reactions, etc."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Gas Evolution", "Reaction Prediction", "Effervescence", "Chemical Reaction"]
    required_envs = []

    code_input_sig = [
        ("reactants", "list", "N/A", "List of reactant species, e.g., ['Zn', 'HCl'], ['Na2CO3', 'HCl'], ['H2O2', 'MnO2']."),
        ("conditions", "str", "room_temperature", "Optional conditions: 'room_temperature', 'heating', 'electrolysis', 'catalyst'. Default is room temperature."),
    ]

    text_input_sig = [
        ("input_str", "str", "N/A", "Space or comma-separated reactants, e.g., 'Zn + HCl', 'Na2CO3 HCl', 'H2O2 MnO2 catalyst'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with will_evolve_gas, gas_products, predicted_equation, gas_type, explanation."),
    ]

    examples = [
        {
            "code_input": {
                "reactants": ["Zn", "HCl"],
                "conditions": "room_temperature",
            },
            "text_input": {
                "input_str": "Zn HCl",
            },
            "output": {
                "result": {
                    "will_evolve_gas": True,
                    "gas_products": ["H2"],
                    "predicted_equation": "Zn + 2HCl = ZnCl2 + H2↑",
                    "gas_type": "hydrogen",
                    "explanation": "Active metal (Zn) reacts with non-oxidizing acid (HCl) to produce hydrogen gas via single displacement: Zn + 2H⁺ → Zn²⁺ + H₂↑.",
                }
            },
        },
        {
            "code_input": {
                "reactants": ["Na2CO3", "HCl"],
                "conditions": "room_temperature",
            },
            "text_input": {
                "input_str": "Na2CO3 HCl",
            },
            "output": {
                "result": {
                    "will_evolve_gas": True,
                    "gas_products": ["CO2"],
                    "predicted_equation": "Na2CO3 + 2HCl = 2NaCl + H2O + CO2↑",
                    "gas_type": "carbon_dioxide",
                    "explanation": "Carbonate reacts with acid to produce carbon dioxide gas: CO₃²⁻ + 2H⁺ → H₂O + CO₂↑.",
                }
            },
        },
        {
            "code_input": {
                "reactants": ["NaOH", "HCl"],
                "conditions": "room_temperature",
            },
            "text_input": {
                "input_str": "NaOH HCl",
            },
            "output": {
                "result": {
                    "will_evolve_gas": False,
                    "gas_products": [],
                    "predicted_equation": "NaOH + HCl = NaCl + H2O",
                    "gas_type": None,
                    "explanation": "Acid-base neutralization produces salt and water — no gas evolution expected under normal conditions.",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize gas evolution pattern database."""
        # Pattern database: each entry defines a class of gas-producing reactions
        self._patterns = [
            # ── Hydrogen production (acid + active metal) ──
            {
                "type": "hydrogen",
                "gas": "H2",
                "reactant_patterns": [
                    (["metal_active"], ["acid_non_oxidizing"]),
                ],
                "active_metals": {"K", "Ca", "Na", "Mg", "Al", "Zn", "Fe", "Ni", "Sn", "Pb"},
                "non_oxidizing_acids": {"HCl", "HBr", "HI", "dilute H2SO4", "CH3COOH"},
                "equation_template": "{metal} + {acid} = {salt} + H2↑",
                "explanation_template": "Active metal {metal} reacts with non-oxidizing acid {acid} to produce hydrogen gas via single displacement.",
            },
            # ── Carbon dioxide (acid + carbonate/bicarbonate) ──
            {
                "type": "carbon_dioxide",
                "gas": "CO2",
                "reactant_patterns": [
                    (["carbonate"], ["acid"]),
                    (["bicarbonate"], ["acid"]),
                ],
                "carbonates": {"Na2CO3", "K2CO3", "CaCO3", "BaCO3", "MgCO3", "ZnCO3", "(NH4)2CO3", "FeCO3", "CuCO3", "Ag2CO3", "CO3^2-"},
                "bicarbonates": {"NaHCO3", "KHCO3", "Ca(HCO3)2", "NH4HCO3", "HCO3^-"},
                "acids": {"HCl", "H2SO4", "HNO3", "HBr", "HI", "CH3COOH", "H3PO4", "H+"},
                "equation_template": "{carbonate} + {acid} = {salt} + H2O + CO2↑",
                "explanation_template": "{carbonate} reacts with acid {acid} to produce carbon dioxide gas: CO₃²⁻ + 2H⁺ → H₂O + CO₂↑.",
            },
            # ── Sulfur dioxide (acid + sulfite) ──
            {
                "type": "sulfur_dioxide",
                "gas": "SO2",
                "reactant_patterns": [
                    (["sulfite"], ["acid"]),
                ],
                "sulfites": {"Na2SO3", "K2SO3", "BaSO3", "CaSO3", "SO3^2-", "HSO3^-"},
                "acids": {"HCl", "H2SO4", "HNO3", "H+"},
                "equation_template": "{sulfite} + {acid} = {salt} + H2O + SO2↑",
                "explanation_template": "Sulfite reacts with acid to produce sulfur dioxide gas: SO₃²⁻ + 2H⁺ → H₂O + SO₂↑.",
            },
            # ── Hydrogen sulfide (acid + sulfide) ──
            {
                "type": "hydrogen_sulfide",
                "gas": "H2S",
                "reactant_patterns": [
                    (["sulfide"], ["non_oxidizing_acid"]),
                ],
                "sulfides": {"Na2S", "K2S", "FeS", "ZnS", "PbS", "CdS", "CuS", "S^2-", "HS^-"},
                "non_oxidizing_acids": {"HCl", "dilute H2SO4", "HBr"},
                "equation_template": "{sulfide} + {acid} = {salt} + H2S↑",
                "explanation_template": "Metal sulfide reacts with non-oxidizing acid to produce hydrogen sulfide gas (rotten egg smell).",
            },
            # ── Oxygen (peroxide decomposition, chlorate decomposition) ──
            {
                "type": "oxygen",
                "gas": "O2",
                "reactant_patterns": [
                    (["peroxide"], ["catalyst_or_heat"]),
                    (["chlorate"], ["heat_or_catalyst"]),
                    (["potassium_permanganate"], ["heat"]),
                ],
                "peroxides": {"H2O2", "Na2O2", "BaO2"},
                "catalysts": {"MnO2", "Fe2O3", "catalyst"},
                "chlorates": {"KClO3", "NaClO3"},
                "equation_templates": {
                    "peroxide": "2H2O2 --(MnO2)--> 2H2O + O2↑",
                    "chlorate": "2KClO3 --(heat/MnO2)--> 2KCl + 3O2↑",
                    "kmno4": "2KMnO4 --(heat)--> K2MnO4 + MnO2 + O2↑",
                },
                "explanation_template": "Decomposition of {reactant} releases oxygen gas.",
            },
            # ── Nitrogen oxides / Ammonia ──
            {
                "type": "nitrogen_oxide",
                "gas": ["NO", "NO2"],
                "reactant_patterns": [
                    (["copper"], ["dilute_HNO3"]),   # NO
                    (["copper"], ["concentrated_HNO3"]),  # NO2
                    (["nitrate"], ["heat_conc_acid"]),
                ],
                "equation_templates": {
                    "cu_dilute_hno3": "3Cu + 8HNO3(dilute) = 3Cu(NO3)2 + 2NO↑ + 4H2O",
                    "cu_conc_hno3": "Cu + 4HNO3(conc) = Cu(NO3)2 + 2NO2↑ + 2H2O",
                },
                "explanation_template": "Copper reacts with nitric acid to produce nitrogen oxide gases.",
            },
            # ── Ammonia ──
            {
                "type": "ammonia",
                "gas": "NH3",
                "reactant_patterns": [
                    (["ammonium_salt"], ["base_heated"]),
                ],
                "ammonium_salts": {"NH4Cl", "NH4NO3", "(NH4)2SO4", "NH4HCO3", "CH3COONH4"},
                "bases": {"NaOH", "KOH", "Ca(OH)2", "CaO"},
                "equation_template": "2NH4Cl + Ca(OH)2 --(heat)--> 2NH3↑ + CaCl2 + 2H2O",
                "explanation_template": "Ammonium salt reacts with base upon heating to release ammonia gas (pungent smell).",
            },
            # ── Chlorine (oxidizing acid + chloride / electrolysis) ──
            {
                "type": "chlorine",
                "gas": "Cl2",
                "reactant_patterns": [
                    (["KMnO4"], ["concentrated_HCl"]),
                    (["chloride"], ["strong_oxidizer"]),
                ],
                "equation_template": "2KMnO4 + 16HCl(conc) = 2KCl + 2MnCl2 + 5Cl2↑ + 8H2O",
                "explanation_template": "Strong oxidizer reacts with concentrated hydrochloric acid to produce chlorine gas.",
            },
        ]

        # Species classification lookup
        self._classify = {}
        for p in self._patterns:
            for k, v in p.items():
                if k.endswith("s") and isinstance(v, set) and not k.startswith(("equation", "explanation")):
                    for item in v:
                        key = k.rstrip("s") if not k.endswith("acids") else "acid"
                        self._classify[item.lower()] = key

    def _run_base(self, reactants: List[str], conditions: str = "room_temperature") -> dict:
        """Core logic: predict gas evolution."""
        if not reactants:
            raise ChemMCPError("Reactants list cannot be empty.")

        conditions_lower = conditions.lower().strip() if conditions else "room_temperature"

        # Store original reactant names for result formatting
        self._last_reactants = [r.strip() for r in reactants]

        # Classify each reactant
        classified = [self._classify_species(r.strip()) for r in reactants]
        logger.info(f"Classified reactants: {classified}")

        # Check against patterns
        for pattern in self._patterns:
            match = self._check_pattern(classified, pattern, conditions_lower)
            if match:
                return match

        return {
            "will_evolve_gas": False,
            "gas_products": [],
            "predicted_equation": None,
            "gas_type": None,
            "explanation": (
                f"No known gas-producing pattern matches the combination of {', '.join(reactants)} "
                f"under '{conditions}' conditions. This does not guarantee no gas is produced — "
                f"only that it doesn't match common textbook patterns."
            ),
        }

    def _run_text(self, input_str: str) -> dict:
        """Parse text input."""
        s = input_str.strip()
        # Split by + or space/comma
        parts = re.split(r'[\s+,]+', s)
        reactants = [p.strip() for p in parts if p.strip()]

        # Check for condition keywords
        conditions = "room_temperature"
        cond_keywords = {"heating", "heat", "hot", "catalyst", "electrolysis", "electrolyze"}
        filtered_reactants = []
        for r in reactants:
            if r.lower() in cond_keywords:
                conditions = r.lower()
            else:
                filtered_reactants.append(r)

        return self._run_base(filtered_reactants, conditions)

    def _classify_species(self, species: str) -> str:
        """Classify a chemical species into a category."""
        s = species.lower().strip()

        # Direct lookup
        if s in self._classify:
            return self._classify[s]

        # Pattern matching
        import re

        # Active metals
        if s in {"k", "ca", "na", "mg", "al", "zn", "fe", "ni", "sn", "pb"}:
            return "metal_active"

        # Copper
        if s == "cu":
            return "copper"

        # Acids
        if s in {"hcl", "hbr", "hi"}:
            return "acid_non_oxidizing"
        if s in {"h2so4", "hno3", "ch3cooh", "h3po4"}:
            return "acid"
        if "dilute" in s and "hno3" in s:
            return "dilute_HNO3"
        if "conc" in s and "hno3" in s:
            return "concentrated_HNO3"
        if "conc" in s and "hcl" in s:
            return "concentrated_HCL"
        if "dilute" in s and "h2so4" in s:
            return "acid_non_oxidizing"

        # Carbonates
        if "co3" in s or "carbonate" in s:
            if "hco3" in s or "bicarbon" in s:
                return "bicarbonate"
            return "carbonate"

        # Sulfites
        if "so3" in s and "sulfite" in s.lower():
            return "sulfite"

        # Sulfides
        if s.endswith("s") and len(s) <= 4 and s[0].isupper():
            return "sulfide"

        # Peroxides
        if "o2" in s or "peroxide" in s.lower():
            return "peroxide"

        # Chlorates
        if "clo3" in s:
            return "chlorate"

        # Permanganate
        if "mno4" in s:
            return "potassium_permanganate"

        # Bases
        if s in {"naoh", "koh", "ca(oh)2", "cao", "ba(oh)2"}:
            return "base_heated" if "heat" in "" else "base"

        # Ammonium salts
        if "nh4" in s:
            return "ammonium_salt"

        # Catalysts
        if s in {"mno2", "fe2o3", "catalyst"}:
            return "catalyst_or_heat"

        # Strong oxidizers
        if "kmno4" in s:
            return "strong_oxidizer"

        # Nitrates
        if "no3" in s:
            return "nitrate"

        return "unknown"

    def _check_pattern(self, classified: List[str], pattern: dict, conditions: str) -> Optional[dict]:
        """Check if classified reactants match a gas evolution pattern."""
        required_categories = set()

        # Build required category set from pattern
        has_carbonate = any(c in pattern.get("carbonates", set()) for c in [""])  # placeholder
        # Simpler approach: check if our classified list contains needed types

        # Special handling by gas type
        ptype = pattern.get("type", "")

        # Hydrogen: need active metal + non-oxidizing acid
        if ptype == "hydrogen":
            has_metal = any("metal" in c for c in classified)
            has_acid = any("acid" in c for c in classified)
            has_oxidizing = any(x in classified for x in ["concentrated_HNO3", "concentrated_HCL"])
            if has_metal and has_acid and not has_oxidizing:
                metal = next((r for r, c in zip(self._last_reactants, classified) if "metal" in c), "M")
                acid = next((r for r, c in zip(self._last_reactants, classified) if "acid" in c), "HA")
                return self._build_result(pattern, metal=metal, acid=acid)

        # CO2: need carbonate/bicarbonate + acid
        if ptype == "carbon_dioxide":
            has_c = any(c in ["carbonate", "bicarbonate"] for c in classified)
            has_a = any(c.startswith("acid") or c == "acid_non_oxidizing" for c in classified)
            if has_c and has_a:
                carb = next((r for r, c in zip(self._last_reactants, classified) if c in ["carbonate", "bicarbonate"]), "M2CO3")
                acid = next((r for r, c in zip(self._last_reactants, classified) if "acid" in c), "HA")
                return self._build_result(pattern, carbonate=carb, acid=acid)

        # SO2: sulfite + acid
        if ptype == "sulfur_dioxide":
            has_s = any(c == "sulfite" for c in classified)
            has_a = any(c.startswith("acid") for c in classified)
            if has_s and has_a:
                sulf = next((r for r, c in zip(self._last_reactants, classified) if c == "sulfite"), "SO3^2-")
                acid = next((r for r, c in zip(self._last_reactants, classified) if "acid" in c), "HA")
                return self._build_result(pattern, sulfite=sulf, acid=acid)

        # H2S: sulfide + non-oxidizing acid
        if ptype == "hydrogen_sulfide":
            has_s = any(c == "sulfide" for c in classified)
            has_a = any(c in ["acid_non_oxidizing"] for c in classified)
            if has_s and has_a:
                sulf = next((r for r, c in zip(self._last_reactants, classified) if c == "sulfide"), "MS")
                acid = next((r for r, c in zip(self._last_reactants, classified) if "acid" in c), "HA")
                return self._build_result(pattern, sulfide=sulf, acid=acid)

        # O2: peroxide/chlorate/permanganate + catalyst/heat
        if ptype == "oxygen":
            has_peroxide = any(c == "peroxide" for c in classified)
            has_chlorate = any(c == "chlorate" for c in classified)
            has_kmno4 = any(c == "potassium_permanganate" for c in classified)
            has_trigger = any(c in ["catalyst_or_heat", "catalyst"] for c in classified) or "heat" in conditions
            if (has_peroxide or has_chlorate or has_kmno4) and has_trigger:
                react = next((r for r, c in zip(self._last_reactants, classified) if c in ["peroxide", "chlorate", "potassium_permanganate"]), "X")
                return self._build_result(pattern, reactant=react)

        # NO/NO2: copper + nitric acid
        if ptype == "nitrogen_oxide":
            has_cu = any(c == "copper" for c in classified)
            has_dil_hno3 = any(c == "dilute_HNO3" for c in classified)
            has_conc_hno3 = any(c == "concentrated_HNO3" for c in classified)
            if has_cu and has_dil_hno3:
                return {
                    "will_evolve_gas": True,
                    "gas_products": ["NO"],
                    "predicted_equation": "3Cu + 8HNO3(dilute) = 3Cu(NO3)2 + 2NO↑ + 4H2O",
                    "gas_type": "nitric_oxide_NO",
                    "explanation": "Copper reacts with dilute nitric acid to produce colorless NO gas (turns brown in air as NO → NO2).",
                }
            if has_cu and has_conc_hno3:
                return {
                    "will_evolve_gas": True,
                    "gas_products": ["NO2"],
                    "predicted_equation": "Cu + 4HNO3(conc) = Cu(NO3)2 + 2NO2↑ + 2H2O",
                    "gas_type": "nitrogen_dioxide_NO2",
                    "explanation": "Copper reacts with concentrated nitric acid to produce brown NO2 gas.",
                }

        # NH3: ammonium salt + base + heat
        if ptype == "ammonia":
            has_nh4 = any(c == "ammonium_salt" for c in classified)
            has_base = any("base" in c for c in classified)
            is_heated = "heat" in conditions or "heating" in conditions
            if has_nh4 and has_base and is_heated:
                nh4 = next((r for r, c in zip(self._last_reactants, classified) if c == "ammonium_salt"), "NH4X")
                return self._build_result(pattern)

        # Cl2: strong oxidizer + concentrated HCl
        if ptype == "chlorine":
            has_ox = any(c in ["strong_oxidizer", "potassium_permanganate"] for c in classified)
            has_conc_hcl = any(c in ["concentrated_HCL"] for c in classified)
            if has_ox and has_conc_hcl:
                return self._build_result(pattern)

        return None

    def _build_result(self, pattern: dict, **kwargs) -> dict:
        """Build result dict from matched pattern."""
        gas = pattern.get("gas", "")
        ptype = pattern.get("type", "")

        # Get equation
        eq_templates = pattern.get("equation_templates", {})
        if ptype == "oxygen":
            reactant_val = kwargs.get("reactant", "X")
            if "h2o2" in str(reactant_val).lower() or "na2o2" in str(reactant_val).lower():
                eq = eq_templates.get("peroxide", "{react} decomposes to H2O + O2")
            elif "clo3" in str(reactant_val).lower():
                eq = eq_templates.get("chlorate", "{react} decomposes to KCl + O2")
            elif "mno4" in str(reactant_val).lower():
                eq = eq_templates.get("kmno4", "{react} decomposes")
            else:
                eq = eq_templates.get("peroxide", "{react} decomposes")
        elif ptype == "nitrogen_oxide":
            eq = "Nitrogen oxide gas produced from nitric acid reaction"
        elif ptype == "hydrogen":
            metal = kwargs.get("metal", "M")
            acid = kwargs.get("acid", "HA")
            salt = f"{metal}({acid.lstrip('H') if acid.startswith('H') and acid != 'HF' and acid != 'HI' else 'Cl'}" if acid.startswith('H') else f"{metal}{acid}"
            eq = f"{metal} + {acid} = {salt} + H2↑"
        elif ptype == "carbon_dioxide":
            carb = kwargs.get("carbonate", "M2CO3")
            acid = kwargs.get("acid", "HA")
            eq = f"{carb} + {acid} = salt + H2O + CO2↑"
        elif ptype == "ammonia":
            eq = pattern.get("equation_template", "Ammonia gas released upon heating")
        elif ptype == "chlorine":
            eq = pattern.get("equation_template", "Cl2 gas produced")
        else:
            eq = pattern.get("equation_template", "Gas-producing reaction detected.")

        # Format equation safely
        try:
            eq = eq.format(**kwargs)
        except (KeyError, IndexError):
            pass  # keep as-is if format fails

        expl = pattern.get("explanation_template", "Gas-producing reaction detected.")
        try:
            expl = expl.format(**kwargs)
        except (KeyError, IndexError):
            pass

        return {
            "will_evolve_gas": True,
            "gas_products": [gas] if isinstance(gas, str) else gas,
            "predicted_equation": eq,
            "gas_type": ptype,
            "explanation": expl,
        }
