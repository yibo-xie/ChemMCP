"""
Collision Energy Optimizer - Suggests optimal collision energy (CE)
for MS/MS experiments based on precursor m/z, compound type, and instrument.
"""

import logging
import math
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# CE optimization parameters: (compound_type, slope, intercept, min_CE, max_CE, optimal_range_center)
# Formula: CE = slope * precursor_mz + intercept
# Based on published empirical data from various instrument platforms
_CE_DATABASE: Dict[str, Dict[str, float]] = {
    "small_molecule":      {"slope": 0.038, "intercept": -2.0,  "min_ce": 5,   "max_ce": 60,  "optimal_ratio": 0.42},
    "peptide":             {"slope": 0.044, "intercept": -4.5,  "min_ce": 15,  "max_ce": 70,  "optimal_ratio": 0.45},
    "lipid":               {"slope": 0.035, "intercept": -1.0,  "min_ce": 10,  "max_ce": 55,  "optimal_ratio": 0.40},
    "nucleotide":          {"slope": 0.040, "intercept": -3.0,  "min_ce": 10,  "max_ce": 65,  "optimal_ratio": 0.43},
    "carbohydrate":        {"slope": 0.032, "intercept": 0.0,   "min_ce": 5,   "max_ce": 45,  "optimal_ratio": 0.38},
    "halogenated_compound":{"slope": 0.036, "intercept": -1.5,  "min_ce": 8,   "max_ce": 50,  "optimal_ratio": 0.41},
    "steroid":             {"slope": 0.039, "intercept": -3.0,  "min_ce": 10,  "max_ce": 55,  "optimal_ratio": 0.43},
    "drug_like":           {"slope": 0.037, "intercept": -2.5,  "min_ce": 8,   "max_ce": 58,  "optimal_ratio": 0.42},
    "metabolite":          {"slope": 0.036, "intercept": -2.0,  "min_ce": 8,   "max_ce": 55,  "optimal_ratio": 0.41},
    "polymer":             {"slope": 0.030, "intercept": 0.0,   "min_ce": 5,   "max_ce": 40,  "optimal_ratio": 0.36},
}

# Instrument-specific CE adjustment factors
_INSTRUMENT_CE_ADJUST = {
    "triple_quad":   {"factor": 1.00, "ce_step": 2,  "notes": "AB Sciex/Agilent/Waters triple quadrupole"},
    "qtrap":         {"factor": 1.05, "ce_step": 2,  "notes": "AB Sciex QTrap — slightly higher CE"},
    "tsq":           {"factor": 0.95, "ce_step": 2,  "notes": "Thermo TSQ series"},
    "xevo_tqd":      {"factor": 0.98, "ce_step": 2,  "notes": "Waters Xevo TQD"},
    "qtof":          {"factor": 0.85, "ce_step": 5,  "notes": "Q-TOF (lower CE for in-source-like fragmentation)"},
    "orbitrap":      {"factor": 0.70, "ce_step": 10, "notes": "Orbitrap HCD (normalized CE scale)"},
    "fticr":         {"factor": 0.65, "ce_step": 10, "notes": "FT-ICR CID/IRMPD"},
}


@ChemMCPManager.register_tool
class CollisionEnergyOptimizer(BaseTool):
    """
    碰撞能量优化建议器 — 为 MS/MS 实验推荐最优碰撞能量（CE）。
    
    基于前体离子 m/z、化合物类型和仪器平台，使用经验公式计算推荐碰撞能量，
    并提供优化曲线和实验建议。
    """
    __version__      = "0.1.0"
    name             = "CollisionEnergyOptimizer"
    func_name        = "optimize_collision_energy"
    description      = "Suggest optimal collision energy (CE) for MS/MS experiments based on precursor m/z, compound type, and instrument platform."
    implementation_description = "Uses the linear CE formula (CE = slope × m/z + intercept) calibrated for different compound classes and instrument types. Generates theoretical breakdown curves to help plan CE ramping experiments."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Mass Spectrometry", "Collision Energy", "MS/MS", "Method Development", "CID"]
    required_envs    = []

    code_input_sig   = [
        ("precursor_mz", "float", "N/A", "Precursor ion m/z value."),
        ("ionization_mode", "str", "ESI+", "Ionization mode: 'ESI+', 'ESI-', 'APCI+', 'APCI-', 'EI'."),
        ("compound_type", "str", "small_molecule", "Compound category: 'small_molecule', 'peptide', 'lipid', 'nucleotide', 'carbohydrate', 'halogenated_compound', 'steroid', 'drug_like', 'metabolite', 'polymer'."),
        ("instrument_type", "str", "triple_quad", "Instrument: 'triple_quad', 'qtrap', 'tsq', 'xevo_tqd', 'qtof', 'orbitrap', 'fticr'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'precursor_mz [ionization_mode] [compound_type] [instrument_type]'. Example: '286.14 ESI+ small_molecule triple_quad'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with recommended_CE, CE_range, optimization_curve_points (CE vs relative intensity), rationale, and experimental protocol suggestions."),
    ]

    examples         = [
        {
            "code_input": {
                "precursor_mz": 286.1438,
                "ionization_mode": "ESI+",
                "compound_type": "small_molecule",
                "instrument_type": "triple_quad",
            },
            "text_input": {
                "input_params": "286.14 ESI+ small_molecule triple_quad"
            },
            "output": {
                "result": {
                    "precursor_mz": 286.1438,
                    "recommended_CE_eV": 25,
                    "CE_range_eV": (15, 33),
                    "formula_used": "CE = 0.038 × m/z + (-2.0)",
                    "optimization_curve": [
                        {"CE": 10, "relative_intensity_pct": 20},
                        {"CE": 17, "relative_intensity_pct": 65},
                        {"CE": 25, "relative_intensity_pct": 100},
                        {"CE": 33, "relative_intensity_pct": 75},
                        {"CE": 45, "relative_intensity_pct": 30},
                    ],
                    "rationale": "Small molecule on triple quad; CE scales linearly with product ion m/z",
                }
            },
        },
        {
            "code_input": {
                "precursor_mz": 532.3050,
                "ionization_mode": "ESI+",
                "compound_type": "peptide",
                "instrument_type": "orbitrap",
            },
            "text_input": {
                "input_params": "532.31 ESI+ peptide orbitrap"
            },
            "output": {
                "result": {
                    "recommended_CE_eV": 24,
                    "instrument_note": "Orbitrap uses normalized CE (NCE); value shown as equivalent eV",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, precursor_mz: float, ionization_mode: str = "ESI+",
                  compound_type: str = "small_molecule",
                  instrument_type: str = "triple_quad") -> dict:
        """Core logic."""
        if precursor_mz <= 0:
            raise ChemMCPError("Precursor m/z must be positive.")
        if compound_type not in _CE_DATABASE:
            raise ChemMCPError(f"Unknown compound_type '{compound_type}'. Options: {list(_CE_DATABASE.keys())}")
        if instrument_type not in _INSTRUMENT_CE_ADJUST:
            raise ChemMCPError(f"Unknown instrument_type '{instrument_type}'. Options: {list(_INSTRUMENT_CE_ADJUST.keys())}")

        # Get parameters
        ce_params = _CE_DATABASE[compound_type]
        inst_adj = _INSTRUMENT_CE_ADJUST[instrument_type]

        # Calculate base CE
        base_ce = ce_params["slope"] * precursor_mz + ce_params["intercept"]
        adj_ce = base_ce * inst_adj["factor"]

        # Clamp to valid range
        opt_ce = max(ce_params["min_ce"], min(ce_params["max_ce"], round(adj_ce, 1)))

        # Define optimization range (±40% of optimum or ±10 eV, whichever is larger)
        margin = max(opt_ce * 0.4, 10)
        ce_lo = max(ce_params["min_ce"], round(opt_ce - margin))
        ce_hi = min(ce_params["max_ce"], round(opt_ce + margin))

        # Generate optimization curve points
        curve = self._generate_breakdown_curve(opt_ce, ce_lo, ce_hi, inst_adj["ce_step"])

        # Rationale
        rationale = self._build_rationale(precursor_mz, compound_type, instrument_type, opt_ce, ce_params)

        # Experimental protocol suggestions
        protocol = self._suggest_protocol(opt_ce, ce_lo, ce_hi, instrument_type, compound_type)

        return {
            "result": {
                "precursor_mz": round(precursor_mz, 4),
                "ionization_mode": ionization_mode.upper(),
                "compound_type": compound_type,
                "instrument_type": instrument_type,
                "recommended_CE_eV": opt_ce,
                "CE_range_eV": (ce_lo, ce_hi),
                "base_formula": f"CE = {ce_params['slope']} × m/z + ({ce_params['intercept']})",
                "instrument_adjustment_factor": inst_adj["factor"],
                "optimization_curve": curve,
                "rationale": rationale,
                "experimental_protocol": protocol,
                "notes": (
                    "Collision Energy Optimization Notes:\n"
                    "• CE formula is empirical — always verify with standards\n"
                    "• For MRM: optimize each transition individually\n"
                    "• For PRM/DIA: use CE ramping across the m/z range\n"
                    f"• {inst_adj['notes']}\n"
                    "• Higher CE → more low-m/z fragments; lower CE → more precursor survival\n"
                    "• Consider scheduled MRM with narrow CE windows for best sensitivity"
                ),
            }
        }

    def _generate_breakdown_curve(self, opt_ce: float, lo: int, hi: int, step: int) -> List[dict]:
        """Generate theoretical breakdown curve (Gaussian response around optimum)."""
        curve = []
        sigma = max(8.0, opt_ce * 0.3)
        ce = lo
        while ce <= hi:
            resp = 100 * math.exp(-((ce - opt_ce) ** 2) / (2 * sigma ** 2))
            resp = max(3, resp)  # noise floor
            curve.append({"CE_eV": ce, "relative_intensity_pct": round(resp, 1)})
            ce += step
        return curve

    def _build_rationale(self, prec_mz: float, comp_type: str, inst: str, opt_ce: float, params: dict) -> str:
        """Build explanation for recommended CE."""
        parts = [
            f"Base CE = {params['slope']} × {prec_mz:.1f} + ({params['intercept']}) = {params['slope'] * prec_mz + params['intercept']:.1f} eV",
        ]
        
        if comp_type == "peptide":
            parts.append("Peptides require higher CE due to multiple backbone cleavage sites")
        elif comp_type == "lipid":
            parts.append("Lipids fragment at lower CE due to labile ester/ether linkages")
        elif comp_type == "carbohydrate":
            parts.append("Carbohydrates need lower CE to avoid complete fragmentation to non-informative ions")

        if inst in ("orbitrap", "fticr"):
            parts.append(f"{inst} uses normalized collision energy; convert using instrument-specific calibration")

        parts.append(f"Final recommendation: {opt_ce} eV (adjusted for {inst})")
        return ". ".join(parts)

    def _suggest_protocol(self, opt_ce: int, lo: int, hi: int, inst: str, comp_type: str) -> dict:
        """Suggest experimental protocol."""
        n_points = min(11, max(5, (hi - lo) // max(2, (hi - lo) // 11)))
        step = max(2, (hi - lo) // (n_points - 1))

        protocol = {
            "step_1_screening": f"Coarse screen: test CE values from {lo} to {hi} eV in steps of {step} eV",
            "step_2_fine_tune": f"Fine optimization: ±5 eV around best response (~{opt_ce} eV) in 2 eV steps",
            "step_3_verify": "Verify with authentic standard at 3 concentration levels",
            "estimated_time_min": n_points * 2,  # ~2 min per point including equilibration
        }

        if inst == "triple_quad":
            protocol["step_4_mrm"] = "For MRM: repeat optimization for each product ion transition separately"

        return protocol

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            prec_mz = float(parts[0])
            mode = parts[1] if len(parts) > 1 else "ESI+"
            comp_type = parts[2] if len(parts) > 2 else "small_molecule"
            inst = parts[3] if len(parts) > 3 else "triple_quad"
            return self._run_base(prec_mz, mode, comp_type, inst)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'precursor_mz [mode] [compound_type] [instrument]'")
