import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Digestion method database
DIGESTION_METHODS = {
    "wet_acid": {
        "name": "Wet Acid Digestion (湿法消解)",
        "description": "Uses strong acids (HNO3, HCl, H2SO4, HClO4) at elevated temperature to decompose organic matrix.",
        "suitable_for": ["biological", "food", "soil", "sludge", "plant_tissue", "water_sediment"],
        "unsuitable_for": ["volatile_elements", "high_organic_load"],
        "typical_temperature": "120-250°C",
        "time_range": "1-8 hours",
        "acids": ["HNO3", "HCl", "H2SO4", "HClO4"],
        "pros": ["Low equipment cost", "Suitable for most sample types", "Good for multi-element analysis"],
        "cons": ["Time-consuming", "Uses large amounts of acid", "Risk of contamination", "Safety hazards with HClO4"],
        "elements_at_risk": ["As", "Sb", "Se", "Hg", "Pb (volatile)", "Sn"],
        "safety_notes": ["Use fume hood at all times", "Never add HClO4 to organic material directly", "Start with HNO3 pre-digestion"],
    },
    "microwave": {
        "name": "Microwave-Assisted Digestion (微波消解)",
        "description": "Uses closed-vessel microwave heating with high pressure and temperature for rapid digestion.",
        "suitable_for": ["biological", "food", "soil", "geological", "pharmaceutical", "petroleum", "environmental"],
        "unsuitable_for": ["peroxide-sensitive_samples", "very_large_sample_amounts"],
        "typical_temperature": "180-260°C",
        "time_range": "15-60 minutes",
        "acids": ["HNO3", "HF", "H2O2", "HCl"],
        "pros": ["Very fast", "Low contamination risk", "Closed system prevents volatile loss", "Reproducible", "Small acid volume"],
        "cons": ["High equipment cost", "Limited sample size per vessel", "Requires specialized vessels", "Pressure safety concerns"],
        "elements_at_risk": ["B, Si (need HF for silicates)"],
        "safety_notes": ["Do not exceed vessel pressure rating", "Allow cooling before opening", "Check vessel integrity regularly"],
    },
    "dry_ashing": {
        "name": "Dry Ashing (干灰化法)",
        "description": "Organic matter is destroyed by heating in a muffle furnace at high temperature.",
        "suitable_for": ["food", "plant_material", "biological_tissue", "animal_feed"],
        "unsuitable_for": ["volatile_elements", "halide-containing_samples"],
        "typical_temperature": "450-550°C",
        "time_range": "4-24 hours",
        "acids": [],
        "pros": ["No acid reagents needed", "Large sample amounts possible", "Simple equipment", "Low blank values"],
        "cons": ["Volatile element loss (Hg, As, Se, Pb, Cd)", "Long time required", "Risk of refractory compound formation", "Not suitable for all elements"],
        "elements_at_risk": ["Hg", "As", "Se", "Pb", "Cd", "Cr(VI) reduction to Cr(III)"],
        "safety_notes": ["Control temperature carefully (<550°C to avoid volatilization losses)", "Use ashing aids (Mg(NO3)2) for retention of volatile elements"],
    },
}


@ChemMCPManager.register_tool
class DigestionProtocolSelector(BaseTool):
    """
    根据样品基质推荐消解方法（湿法、微波、干灰化）。
    基于样品类型、目标元素和实验条件进行智能推荐。
    """
    __version__ = "0.1.0"
    name = "DigestionProtocolSelector"
    func_name = "select_digestion_protocol"
    description = "Recommend the optimal digestion protocol (wet acid, microwave, or dry ashing) based on sample matrix and target analytes."
    implementation_description = "Uses a rule-based decision tree matching sample type, target elements, and constraints against a database of digestion methods."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Sample Preparation", "Digestion", "Analytical Chemistry", "Laboratory"]
    required_envs = []

    code_input_sig = [
        ("sample_type", "str", "N/A", "Type of sample matrix (e.g., 'biological', 'food', 'soil', 'water', 'geological', 'plant')."),
        ("target_elements", "str", "N/A", "Target elements to analyze, comma-separated (e.g., 'Fe,Zn,Cu,Pb,Cd,Hg')."),
        ("sample_mass_g", "float", "0.5", "Approximate sample mass in grams."),
        ("priority", "str", "balanced", "Priority: 'speed', 'completeness', 'safety', 'cost', or 'balanced'."),
        ("equipment_available", "str", "all", "Equipment available: 'all', 'microwave_only', 'basic' (hotplate only), or 'furnace_only'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'sample_type target_elements [sample_mass_g] [priority] [equipment_available]'. Elements are comma-separated without spaces."),
    ]

    output_sig = [
        ("recommended_method", "str", "Recommended digestion method name."),
        ("method_details", "dict", "Full details of the recommended method."),
        ("rationale", "str", "Explanation of why this method was recommended."),
        ("protocol_steps", "list", "Step-by-step procedure."),
        ("warnings", "list", "Specific warnings for this combination."),
        ("alternatives", "list", "Alternative methods if primary recommendation not feasible."),
    ]

    examples = [
        {
            "code_input": {
                "sample_type": "food",
                "target_elements": "Pb,Cd,As,Hg",
            },
            "text_input": {
                "input_params": "food Pb,Cd,As,Hg",
            },
            "output": {
                "recommended_method": "Microwave-Assisted Digestion",
                "rationale": "Food samples with volatile target elements (Hg, As) require closed-vessel microwave digestion to prevent loss.",
                "warnings": ["Hg and As are volatile - closed vessel essential.", "Consider adding gold stabilizer for Hg analysis."],
            },
        },
        {
            "code_input": {
                "sample_type": "plant",
                "target_elements": "K,Ca,Mg,Fe,Mn,Zn",
                "sample_mass_g": 1.0,
                "equipment_available": "basic",
            },
            "text_input": {
                "input_params": "plant K,Ca,Mg,Fe,Mn,Zn 1.0 balanced basic",
            },
            "output": {
                "recommended_method": "Wet Acid Digestion",
                "rationale": "Non-volatile elements with basic equipment available; wet acid digestion is cost-effective and suitable.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # Volatile elements that need closed-vessel methods
    VOLATILE_ELEMENTS = {"HG", "AS", "SB", "SE", "PB", "CD", "SN", "BI"}

    def _run_base(
        self,
        sample_type: str,
        target_elements: str,
        sample_mass_g: float = 0.5,
        priority: str = "balanced",
        equipment_available: str = "all",
    ) -> dict:
        """Core logic: select best digestion protocol."""
        if not sample_type or not target_elements:
            raise ChemMCPError("Both sample_type and target_elements are required.")

        # Parse target elements
        elements = set(e.strip().upper() for e in target_elements.split(","))
        has_volatile = bool(elements & self.VOLATILE_ELEMENTS)
        has_si = "SI" in elements  # Silicates need HF

        sample_key = sample_type.lower().strip()

        # Decision logic
        method_key = None
        rationale_parts = []
        warnings = []

        # Check for volatile elements → prefer microwave (closed vessel)
        if has_volatile:
            warnings.append(f"Volatile elements detected ({', '.join(elements & self.VOLATILE_ELEMENTS)}). Closed-vessel method strongly recommended.")
            if equipment_available in ("all", "microwave_only"):
                method_key = "microwave"
                rationale_parts.append("Target contains volatile elements requiring closed-vessel conditions.")
            elif equipment_available == "basic":
                method_key = "wet_acid"
                rationale_parts.append("Volatile elements present but only hotplate available — use wet digestion with reflux condenser.")
                warnings.append("WARNING: Volatile element loss possible with open-vessel wet digestion!")

        # Check for silicates → need HF (microwave)
        elif has_si:
            method_key = "microwave"
            rationale_parts.append("Silicate-containing samples require HF, best handled in microwave closed vessels.")
            if equipment_available != "microwave" and equipment_available != "all":
                warnings.append("HF requires specialized equipment (microwave/PTEE vessels). Silicate dissolution may be incomplete.")

        # Equipment constraint
        if method_key is None:
            if equipment_available == "microwave_only":
                method_key = "microwave"
                rationale_parts.append("Only microwave equipment available.")
            elif equipment_available == "furnace_only":
                method_key = "dry_ashing"
                rationale_parts.append("Only muffle furnace available.")
                if has_volatile:
                    warnings.append("CRITICAL: Dry ashing will cause loss of volatile elements!")
            elif equipment_available == "basic":
                method_key = "wet_acid"
                rationale_parts.append("Basic equipment (hotplate) available; wet acid digestion selected.")
            else:
                # All equipment available — smart selection based on priority
                if priority == "speed":
                    method_key = "microwave"
                    rationale_parts.append("Speed priority: microwave is fastest.")
                elif priority == "cost":
                    method_key = "dry_ashing" if sample_key in ("food", "plant_material", "animal_feed") else "wet_acid"
                    rationale_parts.append("Cost priority: lowest-cost method selected.")
                elif priority == "safety":
                    method_key = "microwave"
                    rationale_parts.append("Safety priority: closed-vessel microwave minimizes exposure.")
                else:
                    # Balanced default
                    method_key = "microwave"
                    rationale_parts.append("Balanced selection: microwave offers best overall performance.")

        method = DIGESTION_METHODS[method_key]
        protocol = self._generate_protocol(method_key, sample_type, elements)

        # Alternatives
        alternatives = []
        for k, v in DIGESTION_METHODS.items():
            if k != method_key:
                alt_score = self._score_alternative(k, sample_key, elements, equipment_available)
                if alt_score > 0.3:
                    alternatives.append({"method": v["name"], "key": k, "suitability_score": round(alt_score, 2)})

        logger.info(f"Digestion selected: {method_key} for {sample_type}, targets={target_elements}")
        return {
            "recommended_method": method["name"],
            "method_details": {k: v for k, v in method.items() if k != "name"},
            "rationale": " ".join(rationale_parts),
            "protocol_steps": protocol,
            "warnings": warnings,
            "alternatives": sorted(alternatives, key=lambda x: x["suitability_score"], reverse=True),
        }

    def _generate_protocol(self, method_key, sample_type, elements) -> list:
        steps = []
        if method_key == "microwave":
            steps = [
                f"Weigh approximately {0.25}-{0.5} g of {sample_type} sample into PTFE digestion vessel.",
                "Add 6 mL concentrated HNO3 and allow pre-reaction for 15 min (if organic).",
                "Add 2 mL H2O2 (30%) as auxiliary oxidant.",
                "Seal vessels and place in rotor according to manufacturer instructions.",
                "Run microwave program: ramp to 180°C over 10 min, hold for 20 min.",
                "Cool to <50°C before opening vessels.",
                "Transfer digestate to volumetric flask, dilute to mark with deionized water.",
                "If solution is cloudy, filter through 0.45 μm membrane.",
            ]
        elif method_key == "wet_acid":
            steps = [
                f"Weigh {0.5}-{1.0} g of {sample_type} sample into Erlenmeyer flask or beaker.",
                "Add 10 mL concentrated HNO3, cover with watch glass.",
                "Heat on hotplate at low temperature (~80°C) until initial reaction subsides.",
                "Increase temperature to 150°C and continue heating until brown fumes cease.",
                "If necessary, add 2-3 mL HClO4 cautiously (in fume hood!) to complete digestion.",
                "Continue heating until solution becomes clear/colorless or white residue forms.",
                "Cool, add 5 mL deionized water, warm to dissolve salts.",
                "Transfer quantitatively to volumetric flask and dilute to mark.",
            ]
        elif method_key == "dry_ashing":
            steps = [
                f"Weigh {2.0}-{5.0} g of {sample_type} sample into porcelain crucible.",
                "Dry in oven at 105°C overnight to remove moisture.",
                "Place crucible in cold muffle furnace, then program: ramp to 300°C over 1 h, hold 30 min.",
                "Ramp to 450-500°C over 1 h, hold 4-8 h until ash is white/gray.",
                "Cool in desiccator.",
                "Dissolve ash in 5-10 mL 1:1 HNO3 with gentle warming.",
                "Transfer to volumetric flask, dilute to mark with deionized water.",
            ]
        return steps

    def _score_alternative(self, key, sample_type, elements, equipment) -> float:
        score = 1.0
        method = DIGESTION_METHODS[key]
        if sample_type not in method["suitable_for"]:
            score -= 0.4
        if key == "dry_ashing" and (elements & self.VOLATILE_ELEMENTS):
            score -= 0.5
        if key == "microwave" and equipment == "furnace_only":
            score -= 0.8
        if key == "wet_acid" and equipment == "furnace_only":
            score -= 0.9
        return max(0.0, score)

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            sample_type = parts[0]
            elements = parts[1]
            mass = float(parts[2]) if len(parts) > 2 else 0.5
            priority = parts[3] if len(parts) > 3 else "balanced"
            equip = parts[4] if len(parts) > 4 else "all"
            return self._run_base(sample_type, elements, mass, priority, equip)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'sample_type elements [mass] [priority] [equipment]'")
