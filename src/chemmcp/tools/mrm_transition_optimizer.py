"""
MRM Transition Optimizer - Optimizes selection of precursor → product ion transitions
for Multiple Reaction Monitoring (MRM) LC-MS/MS quantitative analysis.
"""

import logging
import math
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Common fragment ions for small molecule MRM optimization
# (fragment_formula, typical_mz_range, specificity_notes)
_COMMON_PRODUCT_IONS = [
    ("[M+H-H₂O]⁺", "precursor-18", "Dehydration product — common but not specific"),
    ("[M+H-NH₃]⁺", "precursor-17", "Ammonia loss — amines/amides"),
    ("[M+H-CO]⁺",  "precursor-28", "Carbonyl loss — moderate specificity"),
    ("[M+H-CO₂]⁺", "precursor-44", "Decarboxylation — carboxylic acids"),
    ("[M+H-HCOOH]⁺", "precursor-46", "Formic acid loss — ESI+ common"),
    ("[M+H-CH₃COOH]⁺", "precursor-60", "Acetic acid loss — esters/acids"),
    ("[M+H-C₂H₅OH]⁺", "precursor-46", "Ethanol loss"),
    ("[M+H-C₃H₆O]⁺", "precursor-58", "Acetone/propionaldehyde loss"),
    ("[M+H-C₄H₈O₂]⁺", "precursor-88", "Butyric acid / ester loss"),
    ("Quan ion", "highest_intensity", "Quantifier transition — most intense product"),
    ("Qual ion", "second_highest", "Qualifier transition — confirmatory"),
]

# CE optimization reference data: compound_type → (CE_slope, CE_intercept)
# Based on empirical formula: CE = slope * m/z + intercept
_CE_PARAMS = {
    "small_molecule":      {"slope": 0.038, "intercept": -2.0,  "range": (5, 60)},
    "peptide":             {"slope": 0.044, "intercept": -4.5,  "range": (15, 70)},
    "lipid":               {"slope": 0.035, "intercept": -1.0,  "range": (10, 55)},
    "nucleotide":          {"slope": 0.040, "intercept": -3.0,  "range": (10, 65)},
    "carbohydrate":        {"slope": 0.032, "intercept": 0.0,   "range": (5, 45)},
    "halogenated_compound":{"slope": 0.036, "intercept": -1.5,  "range": (8, 50)},
}

# Instrument-specific adjustments
_INSTRUMENT_ADJUST = {
    "triple_quad":   {"ce_factor": 1.0,  "dwell_min_ms": 5,   "notes": "Standard triple quadrupole"},
    "qtrap":         {"ce_factor": 1.05, "dwell_min_ms": 5,   "notes": "QTrap with enhanced sensitivity"},
    "qtof":          {"ce_factor": 0.85, "dwell_min_ms": 10,  "notes": "Q-TOF (less optimal for MRM)"},
    "orbitrap":      {"ce_factor": 0.7,  "dwell_min_ms": 20,  "notes": "Orbitrap (PRM mode recommended)"},
}


@ChemMCPManager.register_tool
class MRmTransitionOptimizer(BaseTool):
    """
    MRM 离子对优化选择器 — 为 LC-MS/MS 定量分析优化 MRM 跃迁参数。
    
    根据前体离子、候选产物离子和化合物类型，推荐最优的 MRM 跃迁对，
    包括碰撞能量（CE）、去簇电压等参数建议。
    """
    __version__      = "0.1.0"
    name             = "MRmTransitionOptimizer"
    func_name        = "optimize_mrm_transitions"
    description      = "Optimize MRM (Multiple Reaction Monitoring) transition selection for LC-MS/MS quantitative method development."
    implementation_description = "Uses empirical collision energy formulas (CE = slope × m/z + intercept) calibrated by compound class and instrument type. Ranks transitions by predicted specificity and sensitivity. Provides dwell time scheduling recommendations."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Mass Spectrometry", "MRM", "LC-MS/MS", "Quantitation", "Method Development"]
    required_envs    = []

    code_input_sig   = [
        ("compound_name_or_smiles", "str", "N/A", "Compound identifier: name or SMILES string."),
        ("precursor_mz", "float", "N/A", "Precursor ion m/z value."),
        ("product_ions", "list", "None", "List of candidate product ion m/z values (optional; will suggest if not provided)."),
        ("collision_energies_to_test", "list", "None", "List of CE values to test in eV (optional; auto-generates if not provided)."),
        ("compound_type", "str", "small_molecule", "Compound category: 'small_molecule', 'peptide', 'lipid', 'nucleotide', 'carbohydrate', 'halogenated_compound'."),
        ("instrument_type", "str", "triple_quad", "Instrument type: 'triple_quad', 'qtrap', 'qtof', 'orbitrap'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'compound precursor_mz [product_mz1,mz2,...] [compound_type] [instrument_type]'. Example: 'Cocaine 304.154 182.117,198.091 small_molecule triple_quad'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with optimized_transitions list (ranked), recommended_CE_curve, quantifier_transition, qualifier_transitions, dwell_time_schedule, interference_risk_assessment, and method_recommendations."),
    ]

    examples         = [
        {
            "code_input": {
                "compound_name_or_smiles": "Cocaine",
                "precursor_mz": 304.154,
                "product_ions": [182.117, 198.091, 150.068, 82.049],
                "compound_type": "small_molecule",
                "instrument_type": "triple_quad",
            },
            "text_input": {
                "input_params": "Cocaine 304.154 182.117,198.091 small_molecule triple_quad"
            },
            "output": {
                "result": {
                    "compound": "Cocaine",
                    "precursor_mz": 304.154,
                    "quantifier_transition": {"precursor": 304.154, "product": 182.117, "recommended_CE": 25, "specificity_score": 92},
                    "qualifier_transitions": [
                        {"precursor": 304.154, "product": 198.091, "recommended_CE": 22, "specificity_score": 78},
                        {"precursor": 304.154, "product": 150.068, "recommended_CE": 30, "specificity_score": 65},
                    ],
                    "ce_optimization_curve": [{"CE": 10, "expected_response": 30}, {"CE": 20, "expected_response": 78}, {"CE": 28, "expected_response": 95}, {"CE": 40, "expected_response": 72}, {"CE": 50, "expected_response": 45}],
                    "method_recommendations": ["Use 304→182 as quantifier with CE ~25 eV", "Add 304→198 as primary qualifier", "Dwell time ≥ 10 ms per transition"],
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, compound_name_or_smiles: str, precursor_mz: float,
                  product_ions: Optional[List[float]] = None,
                  collision_energies_to_test: Optional[List[float]] = None,
                  compound_type: str = "small_molecule",
                  instrument_type: str = "triple_quad") -> dict:
        """Core logic."""
        if precursor_mz <= 0:
            raise ChemMCPError("Precursor m/z must be positive.")
        if compound_type not in _CE_PARAMS:
            raise ChemMCPError(f"Unknown compound_type '{compound_type}'. Options: {list(_CE_PARAMS.keys())}")
        if instrument_type not in _INSTRUMENT_ADJUST:
            raise ChemMCPError(f"Unknown instrument_type '{instrument_type}'. Options: {list(_INSTRUMENT_ADJUST.keys())}")

        # Get CE parameters
        ce_p = _CE_PARAMS[compound_type]
        inst_adj = _INSTRUMENT_ADJUST[instrument_type]

        # Generate suggested product ions if none provided
        if not product_ions:
            product_ions = self._suggest_product_ions(precursor_mz, compound_type)

        # Generate CE test range if not provided
        if not collision_energies_to_test:
            lo, hi = ce_p["range"]
            step = max(2, (hi - lo) // 12)
            collision_energies_to_test = list(range(lo, hi + 1, step))

        # Calculate optimal CE for each transition
        transitions = []
        for prod_mz in product_ions:
            if prod_mz >= precursor_mz:
                logger.warning(f"Product ion m/z {prod_mz} >= precursor {precursor_mz}, skipping")
                continue

            opt_ce = self._calc_optimal_ce(precursor_mz, prod_mz, ce_p, inst_adj)
            
            # Specificity scoring
            spec_score = self._score_specificity(precursor_mz, prod_mz, compound_type)

            # Generate response curve
            response_curve = self._generate_ce_curve(opt_ce, collision_energies_to_test)

            transitions.append({
                "precursor_mz": round(precursor_mz, 4),
                "product_mz": round(prod_mz, 4),
                "mass_transition": round(precursor_mz - prod_mz, 4),
                "recommended_CE_eV": round(opt_ce, 1),
                "CE_range_eV": (round(max(5, opt_ce - 10), 1), round(opt_ce + 8, 1)),
                "specificity_score": spec_score,
                "sensitivity_rank": 0,  # fill after sorting
                "response_curve": response_curve,
            })

        # Sort by specificity score descending
        transitions.sort(key=lambda x: x["specificity_score"], reverse=True)
        for i, t in enumerate(transitions):
            t["sensitivity_rank"] = i + 1

        # Select quantifier (best overall) and qualifiers
        quantifier = transitions[0] if transitions else None
        qualifiers = transitions[1:4] if len(transitions) > 1 else []

        # Dwell time recommendation
        n_trans = len(transitions)
        cycle_time_suggested = max(0.5, n_trans * inst_adj["dwell_min_ms"] / 1000)

        # Interference assessment
        interference = self._assess_interference(transitions, compound_type)

        return {
            "result": {
                "compound": compound_name_or_smiles,
                "precursor_mz": round(precursor_mz, 4),
                "compound_type": compound_type,
                "instrument_type": instrument_type,
                "optimized_transitions": transitions[:10],
                "quantifier_transition": quantifier,
                "qualifier_transitions": qualifiers,
                "total_transitions_evaluated": len(transitions),
                "ce_optimization_curve_for_quantifier": quantifier["response_curve"] if quantifier else [],
                "dwell_time_schedule": {
                    "min_dwell_ms": inst_adj["dwell_min_ms"],
                    "recommended_dwell_ms": min(100, max(inst_adj["dwell_min_ms"], int(500 / n_trans))),
                    "estimated_cycle_time_s": round(cycle_time_suggested, 3),
                },
                "interference_risk_assessment": interference,
                "method_recommendations": self._generate_method_recommendations(quantifier, qualifiers, inst_adj, compound_type),
                "notes": (
                    "MRM Optimization Notes:\n"
                    "• Always verify transitions with authentic standards\n"
                    "• Use at least 2 transitions per analyte (1 quan + 1+ qual)\n"
                    f"• Instrument: {inst_adj['notes']}\n"
                    "• Optimize DP (declustering potential) and CXP (cell exit potential) separately\n"
                    "• Consider matrix-matched calibration for complex samples"
                ),
            }
        }

    def _calc_optimal_ce(self, prec_mz: float, prod_mz: float, ce_params: dict, inst_adj: dict) -> float:
        """Calculate optimal CE using linear formula adjusted for instrument."""
        base_ce = ce_params["slope"] * prod_mz + ce_params["intercept"]
        adj_ce = base_ce * inst_adj["ce_factor"]
        lo, hi = ce_params["range"]
        return max(lo, min(hi, adj_ce))

    def _suggest_product_ions(self, prec_mz: float, comp_type: str) -> List[float]:
        """Suggest likely product ions based on neutral losses."""
        suggestions = []
        # Common neutral losses to try
        losses = [18.0106, 17.0027, 27.9949, 43.9898, 46.0055, 58.0419, 60.0211,
                  80.0, 98.0, 120.0,  # generic
                  prec_mz * 0.4, prec_mz * 0.5, prec_mz * 0.6, prec_mz * 0.7]  # fraction-based
        for loss in losses:
            prod = prec_mz - loss
            if prod > 50 and prod < prec_mz - 10:
                suggestions.append(round(prod, 1))
        return sorted(set(suggestions))[:8]

    def _score_specificity(self, prec_mz: float, prod_mz: float, comp_type: str) -> int:
        """Score transition specificity (0-100)."""
        score = 50
        
        mass_diff = prec_mz - prod_mz
        
        # Larger mass difference generally more specific
        if mass_diff > 200:
            score += 25
        elif mass_diff > 100:
            score += 18
        elif mass_diff > 50:
            score += 10
        else:
            score -= 10

        # Product ion in low m/z region is less specific (more interference)
        if prod_mz < 100:
            score -= 15
        elif prod_mz > 300:
            score += 10

        # Very characteristic fragments get bonus
        if abs(mass_diff - 43.9898) < 1:  # CO2 loss
            score += 5
        if abs(mass_diff - 18.0106) < 1:  # H2O loss
            score -= 5  # very common, less specific

        return max(0, min(100, score))

    def _generate_ce_curve(self, opt_ce: float, ce_list: List[float]) -> List[dict]:
        """Generate theoretical response vs CE curve (Gaussian-like around optimum)."""
        curve = []
        for ce in ce_list:
            # Gaussian response centered on optimum
            sigma = max(8.0, opt_ce * 0.3)
            response = 100 * math.exp(-((ce - opt_ce) ** 2) / (2 * sigma ** 2))
            # Add some noise floor
            response = max(5, response)
            curve.append({"CE_eV": ce, "expected_relative_response_pct": round(response, 1)})
        return curve

    def _assess_interference(self, transitions: list, comp_type: str) -> dict:
        """Assess risk of isobaric interference."""
        risks = []
        for t in transitions:
            if t["product_mz"] < 100:
                risks.append(f"m/z {t['product_mz']:.1f}: HIGH interference risk (low m/z region)")
            elif t["specificity_score"] < 40:
                risks.append(f"Transition {t['precursor_mz']:.1f}→{t['product_mz']:.1f}: MODERATE specificity concern")

        risk_level = "low"
        if any("HIGH" in r for r in risks):
            risk_level = "high"
        elif any("MODERATE" in r for r in risks):
            risk_level = "medium"

        return {"risk_level": risk_level, "concerns": risks, "mitigation": "Use high-resolution MS/MS confirmation or alternative transitions"}

    def _generate_method_recommendations(self, quantifier: dict, qualifiers: list, inst_adj: dict, comp_type: str) -> List[str]:
        recs = []
        if quantifier:
            recs.append(
                f"Primary quantifier: {quantifier['precursor_mz']}→{quantifier['product_mz']} "
                f"at CE={quantifier['recommended_CE_eV']} eV"
            )
        if qualifiers:
            q_strs = [f"{q['precursor_mz']}→{q['product_mz']}" for q in qualifiers]
            recs.append(f"Qualifiers: {', '.join(q_strs)}")
        recs.append(f"Dwell time ≥ {inst_adj['dwell_min_ms']} ms per transition")
        recs.append(f"Schedule CE optimization experiments ±10 eV around recommended values")
        if comp_type == "small_molecule":
            recs.append("Consider using both positive and negative mode screening")
        return recs

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            compound = parts[0]
            prec_mz = float(parts[1])
            prods = None
            comp_type = "small_molecule"
            inst = "triple_quad"

            if len(parts) > 2 and parts[2].lower() != "none":
                prods = [float(x) for x in parts[2].split(",")]
            if len(parts) > 3:
                comp_type = parts[3]
            if len(parts) > 4:
                inst = parts[4]

            return self._run_base(compound, prec_mz, prods, None, comp_type, inst)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'compound precursor_mz [prod1,prod2,...] [type] [instrument]'")
