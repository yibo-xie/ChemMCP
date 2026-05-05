import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Matrix effect database and mitigation strategies
MATRIX_EFFECT_DATA = {
    "blood_plasma": {
        "name": "Blood / Plasma",
        "major_interferences": ["proteins", "phospholipids", "salts", "lipids", "endogenous_compounds"],
        "matrix_effect_severity": "High",
        "recommended_strategies": [
            {"strategy": "Protein precipitation (PPT)", "details": "Add 3 vol ACN or MeOH, vortex, centrifuge at 10000g for 10 min. Simple but limited cleanup.", "effectiveness": "Moderate (50-70% ME reduction)"},
            {"strategy": "SPE (HLB or mixed-mode)", "details": "Use HLB cartridge for broad-spectrum cleanup. Best for LC-MS analysis.", "effectiveness": "High (80-95% ME reduction)"},
            {"strategy": "Dilute-and-shoot", "details": "Dilute sample 5-20x with mobile phase. Fast but may compromise sensitivity.", "effectiveness": "Low-Moderate (30-50% ME reduction)"},
        ],
        "internal_standard": "Stable isotope-labeled internal standard (SIL-IS) strongly recommended.",
        "calibration": "Matrix-matched calibration or standard addition method.",
    },
    "urine": {
        "name": "Urine",
        "major_interferences": ["salts", "urea", "creatinine", "metabolites", "pH_variability"],
        "matrix_effect_severity": "Moderate-High",
        "recommended_strategies": [
            {"strategy": "Dilution + pH adjustment", "details": "Dilute 2-10x, adjust pH to optimize extraction recovery.", "effectiveness": "Moderate"},
            {"strategy": "LLE or SPE", "details": "LLE for nonpolar analytes; SPE (C18/HLB) for broader range.", "effectiveness": "High"},
            {"strategy": "Enzymatic hydrolysis", "details": "For conjugated metabolites: β-glucuronidase/sulfatase treatment before extraction.", "effectiveness": "N/A (needed for total analyte measurement)"},
        ],
        "internal_standard": "SIL-IS recommended; structural analog acceptable if SIL unavailable.",
        "calibration": "Matrix-matched calibration in pooled blank urine preferred.",
    },
    "food": {
        "name": "Food (solid/semi-solid)",
        "major_interferences": ["pigments", "fats/oils", "proteins", "carbohydrates", "fiber", "water"],
        "matrix_effect_severity": "Very High",
        "recommended_strategies": [
            {"strategy": "QuEChERS extraction", "details": "ACN extraction + MgSO4 salting-out + PSA/C18/GCB cleanup. Standard for pesticide multi-residue.", "effectiveness": "High (80-90%)"},
            {"strategy": "SPE after extraction", "details": "GCB for pigment removal, Florisil for lipid cleanup post-extraction.", "effectiveness": "Very High (90-98%)"},
            {"strategy": "Freezing lipid removal", "details": "Freeze extract at -20°C, centrifuge to remove precipitated lipids.", "effectiveness": "Moderate for lipids only"},
        ],
        "internal_standard": "SIL-IS essential for accurate quantification.",
        "calibration": "Matrix-matched calibration mandatory. Use same matrix type if possible.",
    },
    "soil_sediment": {
        "name": "Soil / Sediment",
        "major_interferences": ["humic_substances", "clay minerals", "organic_matter", "heavy_metals", "particulates"],
        "matrix_effect_severity": "Very High",
        "recommended_strategies": [
            {"strategy": "ASE/PLE + cleanup", "details": "Accelerated solvent extraction with acetone/hexane or DCM/acetone, followed by Florisil/SPE cleanup.", "effectiveness": "High"},
            {"strategy": "Soxhlet extraction + GPC", "details": "Soxhlet extraction with appropriate solvent, gel permeation chromatography for bulk matrix removal.", "effectiveness": "Very High"},
            {"strategy": "Ultrasonic assisted extraction", "details": "Ultrasonication with solvent mixture, centrifugation/filtration, SPE cleanup.", "effectiveness": "Moderate-High"},
        ],
        "internal_standard": "Surrogate standards added before extraction to monitor recovery.",
        "calibration": "Matrix-matched calibration or standard addition recommended.",
    },
    "water_environmental": {
        "name": "Environmental Water (surface/ground/wastewater)",
        "major_interferences": ["dissolved_organic_matter", "suspended_solids", "humic_acids", "variable_ionic_strength", "microorganisms"],
        "matrix_effect_severity": "Low-Moderate",
        "recommended_strategies": [
            {"strategy": "Direct SPE (large volume)", "details": "Pass 0.5-2 L water through SPE cartridge (HLB/C18). Concentrates analytes and removes most matrix.", "effectiveness": "High"},
            {"strategy": "LLE (liquid-liquid)", "details": "Traditional LLE with DCM or ethyl acetate. Good for non-polar compounds.", "effectiveness": "Moderate-High"},
            {"strategy": "Filtration (0.45μm) + SPE", "details": "Pre-filter to remove particulates, then SPE. Essential for turbid samples.", "effectiveness": "High"},
        ],
        "internal_standard": "Added before extraction to correct for recovery losses.",
        "calibration": "Solvent calibration often acceptable for clean waters; matrix-matched for wastewater.",
    },
    "tissue_biological": {
        "name": "Biological Tissue (liver, brain, fat, etc.)",
        "major_interferences": ["high_lipid_content", "proteins", "cellular_debris", "endogenous_compounds"],
        "matrix_effect_severity": "Very High",
        "recommended_strategies": [
            {"strategy": "Homogenization + LLE", "details": "Homogenize tissue, extract with ethyl acetate/hexane, freeze to remove lipids.", "effectiveness": "Moderate-High"},
            {"strategy": "SPE (hybrid/mixed-mode)", "details": "Post-extraction cleanup with mixed-mode SPE for comprehensive lipid/protein removal.", "effectiveness": "High-Very High"},
            {"strategy": "Matrix solid-phase dispersion (MSPD)", "details": "Blend tissue with C18/bonded silica directly, elute analytes. Integrated extraction+cleanup.", "effectiveness": "High"},
        ],
        "internal_standard": "SIL-IS mandatory. Add before homogenization.",
        "calibration": "Matrix-matched calibration using blank tissue of same species/type.",
    },
}


@ChemMCPManager.register_tool
class MatrixMatchingAdvisor(BaseTool):
    """
    基质匹配建议：减少基质效应干扰。
    根据样品类型、检测方法和目标分析物提供基质效应评估和解决方案。
    """
    __version__ = "0.1.0"
    name = "MatrixMatchingAdvisor"
    func_name = "advise_matrix_matching"
    description = "Provide matrix effect assessment and matching strategies to reduce matrix interference in quantitative analysis."
    implementation_description = "Uses a rule-based database of common matrices with interference profiles and validated mitigation strategies."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Matrix Effect", "Sample Preparation", "LC-MS", "Quantitative Analysis", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("sample_matrix", "str", "N/A", "Type of sample matrix (e.g., 'blood_plasma', 'urine', 'food', 'soil_sediment', 'water_environmental', 'tissue_biological')."),
        ("detection_method", "str", "LC-MS", "Detection method: 'LC-MS', 'LC-MS/MS', 'GC-MS', 'HPLC-UV', 'general'."),
        ("analyte_polarity", "str", "moderate", "Analyte polarity: 'nonpolar', 'moderate', 'polar', 'ionic'."),
        ("expected_concentration_range", "str", "trace", "Expected concentration level: 'trace_ppb', 'low_ppm', 'moderate', 'high'."),
        ("available_equipment", "str", "standard", "Equipment available: 'basic', 'standard', 'advanced' (has UPLC-MS/MS, QuEChERS kits, etc.)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'sample_matrix [detection_method] [analyte_polarity] [conc_level] [equipment]'"),
    ]

    output_sig = [
        ("matrix_assessment", "dict", "Assessment of matrix effect severity and major interferences."),
        ("recommended_strategy", "dict", "Primary recommended strategy with details."),
        ("calibration_approach", "str", "Recommended calibration approach."),
        ("internal_standard_advice", "str", "Internal standard recommendations."),
        ("alternative_strategies", "list", "Backup strategies if primary not feasible."),
    ]

    examples = [
        {
            "code_input": {
                "sample_matrix": "food",
                "detection_method": "LC-MS/MS",
                "analyte_polarity": "moderate",
                "expected_concentration_range": "trace_ppb",
            },
            "text_input": {
                "input_params": "food LC-MS/SS moderate trace_ppb",
            },
            "output": {
                "matrix_assessment": {"severity": "Very High"},
                "recommended_strategy": {"strategy": "QuEChERS extraction"},
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        sample_matrix: str,
        detection_method: str = "LC-MS",
        analyte_polarity: str = "moderate",
        expected_concentration_range: str = "trace_ppb",
        available_equipment: str = "standard",
    ) -> dict:
        """Core logic: advise on matrix matching."""
        matrix_key = sample_matrix.lower().strip().replace(" ", "_")
        det = detection_method.upper().strip()
        polarity = analyte_polarity.lower().strip()
        conc = expected_concentration_range.lower().strip()

        # Look up matrix data
        if matrix_key not in MATRIX_EFFECT_DATA:
            # Fuzzy match
            matched = None
            for k in MATRIX_EFFECT_DATA:
                if matrix_key in k or k in matrix_key:
                    matched = k
                    break
            if matched is None:
                raise ChemMCPError(
                    f"Unknown matrix type: '{sample_matrix}'. "
                    f"Available: {', '.join(MATRIX_EFFECT_DATA.keys())}"
                )
            matrix_key = matched

        mdata = MATRIX_EFFECT_DATA[matrix_key]

        # Select best strategy based on equipment and constraints
        strategies = mdata["recommended_strategies"]
        primary_idx = self._select_best_strategy(strategies, det, polarity, conc, available_equipment)
        primary = strategies[primary_idx]
        alternatives = [s for i, s in enumerate(strategies) if i != primary_idx]

        # Calibration advice
        cal_advice = self._calibration_advice(mdata["calibration"], det, conc)

        logger.info(f"Matrix matching advised for {sample_matrix}: {primary['strategy']}")
        return {
            "matrix_assessment": {
                "matrix_type": mdata["name"],
                "severity": mdata["matrix_effect_severity"],
                "major_interferences": mdata["major_interferences"],
                "risk_level": self._assess_risk(det, mdata["matrix_effect_severity"]),
            },
            "recommended_strategy": primary,
            "calibration_approach": cal_advice,
            "internal_standard_advice": mdata["internal_standard"],
            "alternative_strategies": alternatives,
            "general_workflow": self._build_workflow(matrix_key, primary, det),
        }

    def _select_best_strategy(self, strategies, det, polarity, conc, equipment):
        """Select the best strategy index based on constraints."""
        scores = []
        for s in strategies:
            score = 0.0
            name = s["strategy"].lower()
            eff = s["effectiveness"]

            if "very high" in eff.lower():
                score += 3.0
            elif "high" in eff.lower():
                score += 2.0
            elif "moderate" in eff.lower():
                score += 1.0

            if det in ("LC-MS", "LC-MS/MS") and "spe" in name:
                score += 1.5
            if det == "GC-MS" and ("lle" in name or "quchers" in name):
                score += 1.5
            if conc == "trace_ppb" and "high" in eff.lower():
                score += 1.0
            if equipment == "basic" and ("ppt" in name or "dilut" in name):
                score += 1.5
            if equipment == "advanced" and ("spe" in name or "quchers" in name):
                score += 1.0

            scores.append(score)

        max_idx = 0
        max_score = scores[0] if scores else 0
        for i, sc in enumerate(scores):
            if sc > max_score:
                max_score = sc
                max_idx = i
        return max_idx

    def _assess_risk(self, det, severity):
        """Assess overall risk level based on detection method and matrix severity."""
        risk_map = {
            ("LC-MS", "Very High"): "CRITICAL — Severe ion suppression/enhancement likely",
            ("LC-MS", "High"): "HIGH — Significant matrix effects expected",
            ("LC-MS", "Moderate-High"): "MODERATE-HIGH — May need correction",
            ("LC-MS", "Moderate"): "LOW-MODERATE — Usually manageable",
            ("LC-MS", "Low-Moderate"): "LOW — Minimal concern",
            ("LC-MS/MS", "Very High"): "HIGH — MRM reduces but doesn't eliminate matrix effects",
            ("LC-MS/MS", "High"): "MODERATE-HIGH — MRM helps but matrix-matching still needed for quant",
            ("GC-MS", "Very High"): "MODERATE — Less affected than LC-MS but co-elution possible",
            ("HPLC-UV", "Very High"): "LOW-MODERATE — UV less prone to matrix effects than MS",
        }
        key = (det, severity)
        return risk_map.get(key, f"UNKNOWN — Assess empirically for {det} with {severity} matrix")

    def _calibration_advice(self, base_cal, det, conc):
        if det in ("LC-MS", "LC-MS/MS"):
            if conc == "trace_ppb":
                return f"{base_cal} For trace analysis, consider also: post-column infusion test to map matrix effect regions, and use matrix factor evaluation per FDA/EMA guidelines."
            return base_cal
        elif det == "GC-MS":
            if "solvent" in base_cal.lower():
                return base_cal
            return f"{base_cal} GC-MS generally has lower matrix sensitivity; solvent calibration may be adequate for clean extracts."
        return base_cal

    def _build_workflow(self, matrix_key, strategy, det):
        steps = [
            f"1. Sample collection and storage: Follow SOP for {matrix_key.replace('_', ' ')} samples.",
            f"2. Sample preparation: Apply selected strategy ({strategy['strategy']}).",
            f"3. Internal standard addition: Add IS before extraction if possible.",
            f"4. Extraction/cleanup: Follow protocol details for {strategy['strategy']}.",
            f"5. Reconstitution/dilution: Prepare in injection-compatible solvent.",
            f"6. Analysis by {det}: Include QC samples (blank, matrix blank, spiked matrix).",
            f"7. Data review: Check IS response consistency, peak shape, retention time stability.",
        ]
        return steps

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            matrix = parts[0]
            det = parts[1] if len(parts) > 1 else "LC-MS"
            pol = parts[2] if len(parts) > 2 else "moderate"
            conc = parts[3] if len(parts) > 3 else "trace_ppb"
            equip = parts[4] if len(parts) > 4 else "standard"
            return self._run_base(matrix, det, pol, conc, equip)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'matrix [det] [polarity] [conc] [equipment]'")
