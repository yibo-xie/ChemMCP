import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Homogenization method database
HOMOGENIZATION_METHODS = {
    "blender": {
        "name": "High-Speed Blender / Homogenizer",
        "suitable_for": ["fresh_fruits", "vegetables", "meat", "soft_tissue", "plant_material", "food_products"],
        "particle_size_um": "100-500",
        "sample_amount_g": "10-500 g (wet weight)",
        "time_range": "30 sec - 5 min",
        "temperature_control": "Ice bath or short bursts to avoid heating; cryogenic for heat-sensitive analytes",
        "pros": ["Fast", "Large capacity", "Inexpensive", "Easy to use", "Good for high-moisture samples"],
        "cons": ["Heat generation", "Air incorporation", "Cross-contamination risk between samples", "Not suitable for dry/hard materials"],
        "equipment_needed": "Laboratory blender (e.g., Waring, IKA), stainless steel or disposable containers",
        "protocol": [
            "1. Cut sample into ~2 cm pieces (if solid).",
            "2. Weigh representative portion into blender container.",
            "3. For heat-sensitive analytes: pre-chill container, work quickly, use pulse mode.",
            "4. Blend at high speed in 15-30 second pulses with 30 sec cooling intervals.",
            "5. Check homogeneity — no visible large particles should remain.",
            "6. Transfer homogenate to labeled container; mix thoroughly before subsampling.",
            "7. Clean equipment thoroughly between samples to prevent carryover.",
        ],
        "notes": "For trace analysis: dedicate a blender to low-blank work or use disposable liners.",
    },
    "ball_mill": {
        "name": "Ball Mill / Planetary Mill",
        "suitable_for": ["dry_plant_material", "seeds", "grains", "soil", "sediment", "dried_food", "minerals", "hard_samples"],
        "particle_size_um": "< 100 (fine grinding possible to <20 μm)",
        "sample_amount_g": "1-50 g (depending on jar size)",
        "time_range": "2-30 min (depends on hardness and target fineness)",
        "temperature_control": "Cooling pauses; cryogenic milling with liquid N2 for thermolabile/volatile analytes",
        "pros": ["Very fine and uniform particle size", "Reproducible", "Enclosed system (minimal contamination/loss)", "Suitable for hard/dry materials"],
        "cons": ["Small batch size", "Longer time per sample", "Equipment cost", "Risk of cross-contamination if not cleaned properly"],
        "equipment_needed": "Planetary ball mill (e.g., Retsch MM series), grinding jars (agate, zirconia, steel), balls",
        "protocol": [
            "1. Pre-dry sample if moisture >10% (air-dry or oven-dry at ≤40°C for sensitive analytes).",
            "2. Break up large pieces to <1 cm.",
            "3. Weigh sample into grinding jar with appropriate grinding balls (typically 10:1 ball-to-sample mass ratio).",
            "4. For cryogenic milling: add liquid N2, allow to equilibrate, then mill.",
            "5. Set milling program: typically 200-400 rpm, 5-30 min with reverse rotation intervals.",
            "6. Allow jar to cool before opening (especially after cryogenic milling).",
            "7. Collect powder with clean spatula; store in airtight container.",
            "8. Clean jar and balls thoroughly between samples.",
        ],
        "notes": "Choose jar material based on analysis: agate (trace metals), zirconia (general hard), PTFE (avoid contamination).",
    },
    "mortar_pestle_cryogenic": {
        "name": "Mortar & Pestle (Cryogenic Grinding)",
        "suitable_for": ["animal_tissue", "frozen_tissue", "plant_tissue", "small_amounts", "heat_sensitive_analytes"],
        "particle_size_um": "50-200 (variable)",
        "sample_amount_g": "0.5-20 g",
        "time_range": "5-20 min (including liquid N2 addition)",
        "temperature_control": "Continuous liquid N2 addition maintains sample at <-150°C",
        "pros": ["Excellent for heat-sensitive/volatile analytes", "Low cost", "No electrical equipment needed", "Small amounts feasible", "Minimal thermal degradation"],
        "cons": ["Labor intensive", "Inconsistent particle size", "Liquid N2 handling safety concerns", "Not scalable", "Sample loss possible during transfer"],
        "equipment_needed": "Cryogenic mortar and pestle (pre-chilled, usually stainless steel or agate), liquid N2 dewar, PPE (cryo-gloves, face shield)",
        "protocol": [
            "1. Pre-chill mortar and pestle by adding small amount of liquid N2 and letting it evaporate.",
            "2. Cut frozen tissue/sample into small pieces (~0.5 cm) on dry ice.",
            "3. Transfer pieces to chilled mortar, add liquid N2 to cover sample.",
            "4. Wait until N2 stops boiling vigorously (sample is fully frozen).",
            "5. Grind firmly in circular motion, applying downward pressure.",
            "6. Add more liquid N2 if sample thaws during grinding (it becomes sticky/pasty).",
            "7. Continue grinding until uniform fine powder is achieved.",
            "8. Transfer powder to pre-chilled, labeled container using cold spatula.",
            "9. Store immediately at -80°C or proceed to extraction.",
        ],
        "notes": "Essential for enzyme activity assays, labile metabolites, and volatile compound analysis. Always wear cryogenic PPE.",
    },
    "rotor_stator_homogenizer": {
        "name": "Rotor-Stator Homogenizer (e.g., Polytron, Ultra-Turrax)",
        "suitable_for": ["biological_tissue", "cell_disruption", "microbial_cultures", "emulsions", "suspensions", "soft_solid_mixtures"],
        "particle_size_um": "1-50 μm (cell disruption level)",
        "sample_amount_g": "0.5-100 g (or mL for liquids)",
        "time_range": "30 sec - 3 min",
        "temperature_control": "Ice bath essential; intermittent operation to prevent heating",
        "pros": ["Very fine homogenization", "Reproducible speed settings", "Good for cell lysis", "Fast", "Various probe sizes available"],
        "cons": ["Heat generation (requires cooling)", "Probe cleaning between samples", "Limited sample volume per probe size", "Aerosol generation (biohazard risk)"],
        "equipment_needed": "Rotor-stator homogenizer with interchangeable probes, ice bath, tubes/beakers",
        "protocol": [
            "1. Select appropriate probe size for sample volume (probe should be immersed 2-3x its diameter).",
            "2. Place sample in tube/beaker on ice.",
            "3. Immerse probe, start at low speed, then increase gradually.",
            "4. Homogenize in 10-30 second bursts with 30 second cooling on ice.",
            "5. Move probe slowly to ensure all sample is processed.",
            "6. Total processing time: typically 1-3 min depending on tissue type.",
            "7. Rinse probe thoroughly with solvent (water, buffer, or extraction solvent) between samples.",
            "8. Keep homogenate on ice until subsampling or extraction.",
        ],
        "notes": "For RNA work: use RNase-free conditions, process quickly, keep everything cold. Aerosol containment recommended for biohazardous samples.",
    },
    "bead_beater": {
        "name": "Bead Beating / Bead Mill",
        "suitable_for": ["microbial_cells", "tough_plant_material", "spores", "soil_slurry", "biofilm", "lysis_for_omics"],
        "particle_size_um": "Cell disruption level (complete lysis achievable)",
        "sample_amount_g": "0.05-5 g (or mL)",
        "time_range": "30 sec - 10 min",
        "temperature_control": "Cooling blocks; short bursts with cooling intervals; cryogenic beads available",
        "pros": ["Most effective cell disruption method", "Fast", "High-throughput compatible", "Reproducible", "Works on very tough samples"],
        "cons": ["Heat generation", "Potential sample heating", "Bead carryover contamination risk", "Shearing of DNA/RNA (for genomics applications)"],
        "equipment_needed": "Bead beater (e.g., MP Biomedicals FastPrep, Qiagen TissueLyser), lysing matrix tubes with beads",
        "protocol": [
            "1. Select appropriate bead type and size: glass (0.5 mm for bacteria), zirconia/silica (for tough tissues), steel (very tough).",
            "2. Weigh sample into bead beating tube containing beads.",
            "3. Add extraction buffer or appropriate medium (usually 1-2 mL per 100 mg sample).",
            "4. Cap tightly; place tube in adapter.",
            "5. Run beating program: typically 4.0-6.5 m/s for 30-60 sec (adjust based on sample toughness).",
            "6. Cool on ice for 2-3 min; repeat if needed (most samples need 1-3 cycles).",
            "7. Centrifuge to pellet debris; collect supernatant for analysis.",
            "8. For downstream applications requiring intact DNA: reduce speed/time to minimize shearing.",
        ],
        "notes": "The gold standard for microbial lysis and environmental sample homogenization. Optimize bead type, size, and speed for each sample type.",
    },
}


@ChemMCPManager.register_tool
class HomogenizationProtocol(BaseTool):
    """
    固体样品均质化方案设计。
    根据样品类型、目标分析物和样品量推荐合适的均质化方法和操作方案。
    """
    __version__ = "0.1.0"
    name = "HomogenizationProtocol"
    func_name = "design_homogenization_protocol"
    description = "Design an optimal solid sample homogenization protocol based on sample type, target analytes, and sample amount."
    implementation_description = "Uses a database of homogenization methods with detailed protocols, matching sample properties to the most appropriate technique."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Homogenization", "Sample Preparation", "Solid Sample", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("sample_type", "str", "N/A", "Type of solid sample (e.g., 'plant_tissue', 'animal_tissue', 'soil', 'food', 'seeds', 'biological_cells')."),
        ("sample_state", "str", "fresh", "Physical state of sample: 'fresh', 'frozen', 'dried', 'lyophilized'."),
        ("target_analytes", "str", "general", "Target analyte class: 'general', 'volatile', 'thermolabile', 'metals', 'DNA_RNA', 'proteins', 'pesticides'."),
        ("sample_amount_g", "float", "5.0", "Approximate sample amount in grams."),
        ("available_equipment", "str", "standard", "Equipment available: 'basic' (blender only), 'standard', 'advanced' (has ball mill, bead beater, etc.)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'sample_type [sample_state] [target_analytes] [sample_amount_g] [available_equipment]'"),
    ]

    output_sig = [
        ("recommended_method", "dict", "Recommended homogenization method with full details."),
        ("protocol_steps", "list", "Step-by-step protocol."),
        ("rationale", "str", "Why this method was selected."),
        ("alternatives", "list", "Alternative methods if primary not available."),
        ("critical_points", "list", "Critical points to watch for success."),
    ]

    examples = [
        {
            "code_input": {
                "sample_type": "animal_tissue",
                "sample_state": "frozen",
                "target_analytes": "thermolabile",
                "sample_amount_g": 2.0,
            },
            "text_input": {
                "input_params": "animal_tissue frozen thermolabile 2.0",
            },
            "output": {
                "recommended_method": {"name": "Mortar & Pestle (Cryogenic Grinding)"},
                "rationale": "Cryogenic grinding prevents thermal degradation of labile compounds in frozen tissue.",
            },
        },
        {
            "code_input": {
                "sample_type": "soil",
                "sample_state": "dried",
                "target_analytes": "metals",
                "sample_amount_g": 10.0,
            },
            "text_input": {
                "input_params": "soil dried metals 10.0",
            },
            "output": {
                "recommended_method": {"name": "Ball Mill / Planetary Mill"},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        sample_type: str,
        sample_state: str = "fresh",
        target_analytes: str = "general",
        sample_amount_g: float = 5.0,
        available_equipment: str = "standard",
    ) -> dict:
        """Core logic: design homogenization protocol."""
        stype = sample_type.lower().strip().replace(" ", "_")
        state = sample_state.lower().strip()
        analytes = target_analytes.lower().strip()
        equip = available_equipment.lower().strip()

        # Method selection logic
        method_key = self._select_method(stype, state, analytes, sample_amount_g, equip)
        method = HOMOGENIZATION_METHODS[method_key]

        # Build rationale
        rationale = self._build_rationale(method_key, stype, state, analytes)

        # Alternatives
        alternatives = []
        for key, m in HOMOGENIZATION_METHODS.items():
            if key != method_key:
                alt_score = self._score_alternative(key, stype, state, analytes, equip)
                if alt_score >= 0.25:
                    alternatives.append({
                        "method": m["name"],
                        "key": key,
                        "score": round(alt_score, 2),
                        "best_for": ", ".join(m["suitable_for"][:3]),
                    })
        alternatives.sort(key=lambda x: x["score"], reverse=True)

        # Critical points
        critical = self._critical_points(method_key, state, analytes)

        logger.info(f"Homogenization protocol designed: {method['name']} for {sample_type} ({state})")
        return {
            "recommended_method": {k: v for k, v in method.items() if k != "protocol"},
            "protocol_steps": method["protocol"],
            "rationale": rationale,
            "alternatives": alternatives[:3],
            "critical_points": critical,
        }

    def _select_method(self, stype, state, analytes, amount, equip):
        """Select best homogenization method."""

        # Cryogenic priority for volatile/thermolabile/frozen tissue
        if analytes in ("volatile", "thermolabile") and stype in ("animal_tissue", "plant_tissue"):
            return "mortar_pestle_cryogenic"

        # Cell disruption / omics → bead beater
        if analytes in ("dna_rna", "proteins", "microbial_cells") or stype == "microbial_cells":
            return "bead_beater"

        # Biological soft tissue → rotor-stator
        if stype in ("animal_tissue", "biological_tissue") and state == "fresh":
            if equip in ("standard", "advanced"):
                return "rotor_stator_homogenizer"
            else:
                return "blender"

        # Dry/hard materials → ball mill
        if state in ("dried", "lyophilized") or stype in ("soil", "sediment", "seeds", "grains", "minerals"):
            if equip == "advanced":
                return "ball_mill"
            elif equip == "basic" and stype in ("soil", "sediment"):
                return "mortar_pestle_cryogenic"  # Can grind dry too
            else:
                return "ball_mill"

        # Food / plant fresh → blender
        if stype in ("fresh_fruits", "vegetables", "food_products", "plant_material") and state == "fresh":
            return "blender"

        # Default selection based on equipment
        if equip == "basic":
            return "blender"
        elif equip == "advanced":
            return "ball_mill" if state in ("dried", "lyophilized") else "rotor_stator_homogenizer"
        else:
            return "blender"

    def _build_rationale(self, key, stype, state, analytes):
        m = HOMOGENIZATION_METHODS[key]
        parts = [f"{m['name']} selected because:"]
        if any(s.replace("_", "") in stype.replace("_", "") for s in m["suitable_for"]):
            parts.append(f"It is well-suited for {stype.replace('_', ' ')}.")
        if state == "frozen" and key == "mortar_pestle_cryogenic":
            parts.append("Cryogenic grinding maintains sample integrity for frozen samples.")
        if state == "dried" and key == "ball_mill":
            parts.append("Ball mill achieves fine, uniform particle size for dry materials.")
        if analytes == "thermolabile" and "cryogenic" in key:
            parts.append("Cryogenic conditions prevent thermal degradation of labile analytes.")
        if analytes in ("dna_rna", "proteins") and key == "bead_beater":
            parts.append("Bead beating provides complete cell disruption for biomolecule extraction.")
        return " ".join(parts)

    def _score_alternative(self, key, stype, state, analytes, equip):
        score = 0.0
        m = HOMOGENIZATION_METHODS[key]
        if any(s.replace("_", "") in stype.replace("_", "") or stype.replace("_", "") in s.replace("_", "") for s in m["suitable_for"]):
            score += 2.0
        if state == "frozen" and "cryogenic" in key:
            score += 1.5
        elif state == "dried" and key == "ball_mill":
            score += 1.5
        if analytes in ("dna_rna", "proteins") and key == "bead_beater":
            score += 2.0
        if equip == "basic" and key == "blender":
            score += 1.5
        if equip == "advanced" and key in ("ball_mill", "bead_beater"):
            score += 1.0
        return max(0.0, min(score / 4.0, 1.0))

    def _critical_points(self, key, state, analytes):
        points = []
        common = [
            "Clean all equipment thoroughly between samples to prevent cross-contamination.",
            "Document exact homogenization parameters (time, speed, temperature) for reproducibility.",
        ]
        specific = {
            "blender": [
                "Use pulse mode to minimize heat buildup — overheating degrades many analytes.",
                "Do not overfill (>50% capacity) — this reduces homogenization efficiency.",
            ],
            "ball_mill": [
                "Choose grinding jar material compatible with your analysis (agate for metals, zirconia for general).",
                "Do not open jar immediately after milling — pressure buildup and/or extreme cold hazard.",
            ],
            "mortar_pestle_cryogenic": [
                "ALWAYS wear cryogenic PPE (face shield, insulated gloves, lab coat).",
                "Add liquid N2 frequently — sample must remain brittle/crumbly, not thaw and become gummy.",
                "Work in a well-ventilated area — N2 displacement can cause O2 deficiency.",
            ],
            "rotor_stator_homogenizer": [
                "Keep sample on ice AT ALL TIMES — rotor-stator generates significant heat.",
                "Ensure probe is fully immersed to avoid air incorporation and aerosol generation.",
                "Start at low speed to prevent sample splashing.",
            ],
            "bead_beater": [
                "Optimize bead type and beating speed — too aggressive shears DNA/RNA.",
                "Allow sufficient cooling between cycles — heat denatures proteins/DNA.",
                "Check for bead carryover in supernatant (centrifuge adequately).",
            ],
        }
        points.extend(common)
        points.extend(specific.get(key, []))
        return points

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            stype = parts[0]
            state = parts[1] if len(parts) > 1 else "fresh"
            analytes = parts[2] if len(parts) > 2 else "general"
            amount = float(parts[3]) if len(parts) > 3 else 5.0
            equip = parts[4] if len(parts) > 4 else "standard"
            return self._run_base(stype, state, analytes, amount, equip)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'sample_type [state] [analytes] [amount_g] [equipment]'")
