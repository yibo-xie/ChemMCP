"""
Ion Suppression Checker - Evaluates potential ion suppression/enhancement effects
in LC-MS analysis based on analyte properties, matrix type, and chromatography.
"""

import logging
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Matrix type risk data: (matrix_name, base_risk_score, typical_suppressants, notes)
_MATRIX_RISK_DATA: Dict[str, Dict] = {
    "plasma": {
        "base_risk": 0.75,
        "suppressants": ["phospholipids", "salts", "proteins", "triglycerides", "lysophospholipids"],
        "typical_suppression_pct": "30-70%",
        "notes": "Most challenging biological matrix — high phospholipid content",
    },
    "serum": {
        "base_risk": 0.72,
        "suppressants": ["phospholipids", "salts", "albumin", "fatty acids"],
        "typical_suppression_pct": "25-65%",
        "notes": "Similar to plasma but slightly lower protein content",
    },
    "urine": {
        "base_risk": 0.45,
        "suppressants": ["salts", "urea", "creatinine", "metabolites"],
        "typical_suppression_pct": "10-40%",
        "notes": "Variable composition; salt concentration varies with hydration",
    },
    "tissue_homogenate": {
        "base_risk": 0.80,
        "suppressants": ["lipids", "phospholipids", "proteins", "cell debris", "pigments"],
        "typical_suppression_pct": "35-80%",
        "notes": "Very complex matrix; requires extensive sample prep",
    },
    "feces": {
        "base_risk": 0.85,
        "suppressants": ["complex organics", "bile acids", "undigested matter", "microflora products"],
        "typical_suppression_pct": "40-90%",
        "notes": "Extremely challenging; highly variable composition",
    },
    "soil_extract": {
        "base_risk": 0.70,
        "suppressants": ["humic substances", "organic matter", "metal ions", "particulates"],
        "typical_suppression_pct": "20-60%",
        "notes": "Humic acids cause strong suppression in ESI",
    },
    "water_environmental": {
        "base_risk": 0.25,
        "suppressants": ["dissolved organic matter", "particulates"],
        "typical_suppression_pct": "5-25%",
        "notes": "Relatively clean; depends on source (wastewater > drinking water)",
    },
    "food": {
        "base_risk": 0.78,
        "suppressants": ["sugars", "fats", "proteins", "pigments", "additives"],
        "typical_suppression_pct": "30-75%",
        "notes": "Highly variable by food type; fatty foods worst for ESI",
    },
    "cell_culture": {
        "base_risk": 0.55,
        "suppressants": ["media components", "salts", "amino acids", "vitamins"],
        "typical_suppression_pct": "15-50%",
        "notes": "Depends on media composition (e.g., FBS content)",
    },
    "csf": {
        "base_risk": 0.35,
        "suppressants": ["proteins", "salts (low)"],
        "typical_suppression_pct": "5-30%",
        "notes": "Relatively clean compared to plasma; low protein",
    },
    "saliva": {
        "base_risk": 0.40,
        "suppressants": ["proteins", "mucins", "bacteria metabolites"],
        "typical_suppression_pct": "10-35%",
        "notes": "Moderate complexity; mucins can cause issues",
    },
    "milk": {
        "base_risk": 0.73,
        "suppressants": ["fats (high)", "proteins (casein), lactose, calcium"],
        "typical_suppression_pct": "25-65%",
        "notes": "High fat and protein content; defatting critical",
    },
}

# Chromatography method impact on ion suppression
_CHROMATOGRAPHY_IMPACT = {
    "reverse_phase":       {"suppression_reduction": 0.4,  "notes": "RP separates many matrix components; good retention of non-polar analytes"},
    "hilic":               {"suppression_reduction": 0.3,  "notes": "HILIC elutes salts early — can co-elute with polar analytes"},
    "ion_exchange":        {"suppression_reduction": 0.5,  "notes": "Good separation from neutral suppressants"},
    "mixed_mode":          {"suppression_reduction": 0.55, "notes": "Best for complex matrices; dual mechanism separation"},
    "normal_phase":        {"suppression_reduction": 0.35, "notes": "Less common for bioanalysis"},
    "no_chromatography":   {"suppression_reduction": 0.0,  "notes": "Flow injection — maximum suppression risk!"},
}

# Ionization mode susceptibility
_IONIZATION_SUSCEPTIBILITY = {
    "ESI+":  {"base_susceptibility": 0.85, "worst_coeluters": "phospholipids, TFA, amines", "notes": "ESI+ is most prone to competition for droplet charge"},
    "ESI-":  {"base_susceptibility": 0.65, "worst_coeluters": "fatty acids, phospholipids, phenols", "notes": "ESI- less susceptible than ESI+ but still affected"},
    "APCI+": {"base_susceptibility": 0.35, "worst_coeluters": "less sensitive", "notes": "APCI less prone to suppression (gas-phase ionization)"},
    "APCI-": {"base_susceptibility": 0.30, "worst_coeluters": "less sensitive", "notes": "APCI- most robust against matrix effects"},
}


@ChemMCPManager.register_tool
class IonSuppressionChecker(BaseTool):
    """
    离子抑制效应评估器 — 评估 LC-MS 分析中潜在的离子抑制/增强效应。
    
    基于分析物性质、基质类型、电离模式和色谱方法，评估离子抑制风险，
    并提供缓解策略建议。
    """
    __version__      = "0.1.0"
    name             = "IonSuppressionChecker"
    func_name        = "check_ion_suppression"
    description      = "Evaluate potential ion suppression/enhancement effects in LC-MS analysis based on analyte properties, matrix type, ionization mode, and chromatographic method."
    implementation_description = "Uses a multi-factor risk scoring model combining matrix complexity, ionization mode susceptibility, chromatographic separation effectiveness, and analyte-specific properties. Provides quantitative risk score (0-1) and actionable mitigation strategies."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Mass Spectrometry", "Ion Suppression", "Matrix Effects", "LC-MS", "Bioanalysis"]
    required_envs    = []

    code_input_sig   = [
        ("analyte_properties", "str", "N/A", "Analyte info: compound name or dict-like string with 'logP', 'pKa', 'polarity' if available."),
        ("matrix_type", "str", "plasma", "Sample matrix: 'plasma', 'serum', 'urine', 'tissue_homogenate', 'feces', 'soil_extract', 'water_environmental', 'food', 'cell_culture', 'csf', 'saliva', 'milk'."),
        ("ionization_mode", "str", "ESI+", "Ionization: 'ESI+', 'ESI-', 'APCI+', 'APCI-'."),
        ("chromatography_method", "str", "reverse_phase", "LC method: 'reverse_phase', 'hilic', 'ion_exchange', 'mixed_mode', 'normal_phase', 'no_chromatography'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'analyte [matrix_type] [ionization_mode] [chromatography_method]'. Example: 'Ciprofloxacin plasma ESI+ reverse_phase'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with risk_score (0-1), risk_level, contributing_factors (with individual scores), mitigation_strategies (prioritized), matrix_effects_summary, and experimental validation recommendations."),
    ]

    examples         = [
        {
            "code_input": {
                "analyte_properties": "Ciprofloxacin (logP=0.28, pKa=6.09/8.74)",
                "matrix_type": "plasma",
                "ionization_mode": "ESI+",
                "chromatography_method": "reverse_phase",
            },
            "text_input": {
                "input_params": "Ciprofloxacin plasma ESI+ reverse_phase"
            },
            "output": {
                "result": {
                    "risk_score": 0.68,
                    "risk_level": "moderate-high",
                    "estimated_suppression_pct": "30-50%",
                    "contributing_factors": [
                        {"factor": "Matrix complexity (plasma)", "score": 0.75, "weight": "high"},
                        {"factor": "ESI+ ionization", "score": 0.85, "weight": "high"},
                        {"factor": "RP chromatography", "score": 0.40, "weight": "medium"},
                        {"factor": "Polar analyte (low logP)", "score": 0.55, "weight": "medium"},
                    ],
                    "mitigation_strategies": [
                        {"strategy": "Phospholipid removal SPE", "priority": 1, "expected_improvement": "40-60% reduction"},
                        {"strategy": "Chromatographic separation optimization", "priority": 2, "expected_improvement": "20-30% reduction"},
                        {"strategy": "Isotope-labeled internal standard", "priority": 3, "expected_improvement": "Compensates residual effect"},
                    ],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, analyte_properties: str, matrix_type: str = "plasma",
                  ionization_mode: str = "ESI+", chromatography_method: str = "reverse_phase") -> dict:
        """Core logic."""
        # Validate inputs
        if matrix_type not in _MATRIX_RISK_DATA:
            raise ChemMCPError(f"Unknown matrix_type '{matrix_type}'. Options: {list(_MATRIX_RISK_DATA.keys())}")
        if ionization_mode not in _IONIZATION_SUSCEPTIBILITY:
            raise ChemMCPError(f"Unknown ionization_mode '{ionization_mode}'. Options: {list(_IONIZATION_SUSCEPTIBILITY.keys())}")
        if chromatography_method not in _CHROMATOGRAPHY_IMPACT:
            raise ChemMCPError(f"Unknown chromatography_method '{chromatography_method}'. Options: {list(_CHROMATOGRAPHY_IMPACT.keys())}")

        # Get base risk factors
        matrix_data = _MATRIX_RISK_DATA[matrix_type]
        ion_data = _IONIZATION_SUSCEPTIBILITY[ionization_mode]
        chrom_data = _CHROMATOGRAPHY_IMPACT[chromatography_method]

        # Parse analyte properties
        log_p = self._extract_logp(analyte_properties)
        pka_info = self._extract_pka(analyte_properties)
        is_polar = log_p is not None and log_p < 1.0
        is_basic = pka_info is not None and any(5 <= pk <= 11 for pk in pka_info)
        is_acidic = pka_info is not None and any(2 <= pk <= 6 for pk in (pka_info or []))

        # Calculate contributing factor scores
        factors = []

        # Factor 1: Matrix complexity
        matrix_score = matrix_data["base_risk"]
        factors.append({
            "factor": f"Matrix complexity ({matrix_type})",
            "raw_score": round(matrix_score, 3),
            "weight": "high",
            "details": matrix_data["notes"],
            "known_suppressants": matrix_data["suppressants"],
        })

        # Factor 2: Ionization mode susceptibility
        ion_score = ion_data["base_susceptibility"]
        factors.append({
            "factor": f"Ionization mode ({ionization_mode})",
            "raw_score": round(ion_score, 3),
            "weight": "high",
            "details": ion_data["notes"],
            "worst_coeluters": ion_data["worst_coeluters"],
        })

        # Factor 3: Chromatography mitigation
        chrom_risk = 1.0 - chrom_data["suppression_reduction"]
        factors.append({
            "factor": f"Chromatography ({chromatography_method})",
            "raw_score": round(chrom_risk, 3),
            "weight": "medium",
            "details": chrom_data["notes"],
        })

        # Factor 4: Analyte polarity (polar analytes more prone to early elution → co-elution with matrix)
        polarity_score = 0.5
        if is_polar:
            polarity_score = 0.7  # polar compounds often elute near void volume
        elif log_p is not None and log_p > 3:
            polarity_score = 0.3  # well-retained, less suppression
        factors.append({
            "factor": f"Analyte polarity (logP≈{log_p if log_p else 'unknown'})",
            "raw_score": round(polarity_score, 3),
            "weight": "medium",
            "details": "Polar analytes may co-elute with matrix components at the solvent front",
        })

        # Factor 5: Acid/base character
        ab_score = 0.4
        if is_basic and ionization_mode == "ESI+":
            ab_score = 0.65  # basic compounds compete with matrix amines
        elif is_acidic and ionization_mode == "ESI-":
            ab_score = 0.60  # acidic compounds compete with matrix acids
        factors.append({
            "factor": "Acid-base character / charge competition",
            "raw_score": round(ab_score, 3),
            "weight": "medium-low",
            "details": "Ionizable analytes can experience charge competition in ESI",
        })

        # Calculate weighted overall risk score
        weights = {"high": 0.3, "medium": 0.18, "medium-low": 0.12}
        total_risk = sum(f["raw_score"] * weights.get(f["weight"], 0.1) for f in factors)
        total_risk = min(1.0, max(0.0, total_risk))

        # Determine risk level
        if total_risk >= 0.75:
            risk_level = "critical"
        elif total_risk >= 0.55:
            risk_level = "moderate-high"
        elif total_risk >= 0.35:
            risk_level = "moderate"
        elif total_risk >= 0.18:
            risk_level = "low"
        else:
            risk_level = "minimal"

        # Estimate suppression percentage
        est_suppression = self._estimate_suppression(total_risk)

        # Generate mitigation strategies
        mitigations = self._generate_mitigations(total_risk, matrix_type, ionization_mode, chromatography_method, is_polar)

        return {
            "result": {
                "analyte": analyte_properties,
                "matrix_type": matrix_type,
                "ionization_mode": ionization_mode,
                "chromatography_method": chromatography_method,
                "risk_score": round(total_risk, 3),
                "risk_level": risk_level,
                "estimated_suppression_range": est_suppression,
                "contributing_factors": factors,
                "mitigation_strategies": mitigations,
                "matrix_effects_summary": (
                    f"{matrix_data['notes']}. "
                    f"Typical suppression in this matrix: {matrix_data['typical_suppression_pct']}. "
                    f"Primary suppressants: {', '.join(matrix_data['suppressants'][:5])}."
                ),
                "validation_recommendation": self._validation_rec(risk_level, matrix_type),
                "notes": (
                    "Ion Suppression Notes:\n"
                    "• Ion suppression is the #1 source of inaccuracy in LC-MS bioanalysis\n"
                    "• Always use stable isotope-labeled internal standards (SIL-IS) when available\n"
                    "• Post-column infusion experiment is the gold standard for mapping suppression zones\n"
                    "• Matrix-matched calibration is essential when SIL-IS is unavailable\n"
                    "• APCI can be used as alternative when ESI suppression is severe\n"
                    "• Sample clean-up (SPE, PPT, LLE) is the most effective prevention strategy"
                ),
            }
        }

    def _extract_logp(self, props: str) -> Optional[float]:
        """Try to extract logP from analyte properties string."""
        import re
        match = re.search(r'log[Pp]\s*[:=]\s*([+-]?\d+\.?\d*)', props)
        if match:
            return float(match.group(1))
        return None

    def _extract_pka(self, props: str) -> Optional[List[float]]:
        """Try to extract pKa values from properties string."""
        import re
        matches = re.findall(r'p[Ka]\s*[:=]\s*([\d.]+)', props)
        return [float(m) for m in matches] if matches else None

    def _estimate_suppression(self, risk: float) -> str:
        """Convert risk score to estimated suppression range."""
        if risk >= 0.75:
            return "50-90%"
        elif risk >= 0.55:
            return "30-60%"
        elif risk >= 0.35:
            return "15-40%"
        elif risk >= 0.18:
            return "5-20%"
        else:
            return "<10%"

    def _generate_mitigations(self, risk: float, matrix: str, mode: str, chrom: str, is_polar: bool) -> List[dict]:
        """Generate prioritized mitigation strategies."""
        strategies = []

        # Strategy selection based on risk profile
        strategies.append({"strategy": "Stable isotope-labeled internal standard (SIL-IS)", "priority": 1,
                          "effectiveness": 95, "effort": "high", "expected_improvement": "Corrects for residual suppression"})
        
        if matrix in ("plasma", "serum", "tissue_homogenate", "milk"):
            strategies.append({"strategy": "Phospholipid removal (HybridSPE/PPT+)", "priority": 2,
                              "effectiveness": 85, "effort": "medium", "expected_improvement": "40-70% suppression reduction"})
            strategies.append({"strategy": "Protein precipitation + dilution", "priority": 3,
                              "effectiveness": 60, "effort": "low", "expected_improvement": "20-40% reduction"})

        if mode.startswith("ESI") and risk > 0.5:
            alt_mode = "APCI+" if mode == "ESI+" else "APCI-"
            strategies.append({"strategy": f"Switch to {alt_mode} ionization", "priority": 4,
                              "effectiveness": 75, "effort": "low", "expected_improvement": "Significantly less suppression-prone"})

        if chrom == "no_chromatography":
            strategies.append({"strategy": "Add LC separation (even short gradient)", "priority": 2,
                              "effectiveness": 90, "effort": "medium", "expected_improvement": "Dramatic improvement over flow injection"})
        elif is_polar and chrom == "reverse_phase":
            strategies.append({"strategy": "Use HILIC or polar-embedded RP column", "priority": 3,
                              "effectiveness": 70, "effort": "medium", "expected_improvement": "Better retention of polar analytes"})

        strategies.append({"strategy": "Optimize chromatographic separation (shift RT away from matrix region)", "priority": 5,
                          "effectiveness": 70, "effort": "medium", "expected_improvement": "Reduce co-elution"})
        strategies.append({"strategy": "Increase sample dilution factor", "priority": 6,
                          "effectiveness": 50, "effort": "very_low", "expected_improvement": "Simple but reduces sensitivity"})

        # Sort by priority
        strategies.sort(key=lambda x: x["priority"])
        return strategies

    def _validation_rec(self, risk_level: str, matrix: str) -> str:
        """Generate validation recommendation."""
        if risk_level in ("critical", "moderate-high"):
            return (
                "STRONGLY RECOMMENDED: Perform post-column infusion experiment to map suppression zones. "
                "Compare calibration slopes in neat solvent vs matrix extract. "
                "Use Matuszewski validation: matrix factor = peak area(matrix/neat) should be 0.85-1.15."
            )
        elif risk_level == "moderate":
            return (
                "RECOMMENDED: Check matrix factor at LLOQ and ULOQ levels. "
                "Consider post-column infusion if developing regulated bioanalysis method."
            )
        else:
            return (
                "Optional: Standard validation with matrix-matched calibrators should be sufficient. "
                "Monitor for batch-to-batch variability."
            )

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            analyte = parts[0]
            matrix = parts[1] if len(parts) > 1 else "plasma"
            mode = parts[2] if len(parts) > 2 else "ESI+"
            chrom = parts[3] if len(parts) > 3 else "reverse_phase"
            return self._run_base(analyte, matrix, mode, chrom)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'analyte [matrix] [mode] [chromatography]'")
