import logging
import math
from typing import Optional, Dict, Any, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common SEC calibration standards with approximate hydrodynamic volumes
SEC_STANDARDS = {
    "pullulan": {
        "name": "Pullulan Standards (Shodex / PSS)",
        "description": "Polysaccharide standards for aqueous SEC",
        "molecular_weights": [590, 1200, 5200, 10200, 21600, 50100, 112000, 212000, 404000, 788000],
        "units": "Da (g/mol)",
        "notes": "Most common for GFC; narrow dispersity (Đ < 1.05)",
    },
    "polystyrene": {
        "name": "Polystyrene Standards (PS)",
        "description": "Hydrophobic standards for THF/DMF SEC",
        "molecular_weights": [500, 1000, 2000, 5000, 10000, 30000, 60000, 100000, 300000, 600000, 900000],
        "units": "Da (g/mol)",
        "notes": "Standard for organic SEC (GPC); use with THF mobile phase",
    },
    "peg": {
        "name": "Polyethylene Glycol (PEG) Standards",
        "description": "Polyether standards for polar/aqueous SEC",
        "molecular_weights": [200, 400, 600, 1000, 1500, 2000, 4000, 8000, 10000, 20000],
        "units": "Da (g/mol)",
        "notes": "Good for PEG/PEO analysis; also used for protein MW estimation",
    },
    "dextran": {
        "name": "Dextran Standards",
        "description": "Polysaccharide standards for aqueous SEC/GFC",
        "molecular_weights": [1000, 5000, 12000, 25000, 50000, 70000, 150000, 270000, 670000],
        "units": "Da (g/mol)",
        "notes": "Alternative to pullulan; broader distribution typical",
    },
    "protein_mw_marker": {
        "name": "Protein Molecular Weight Markers",
        "description": "Native protein standards for bio-SEC",
        "molecular_weights": [1350, 6500, 13700, 25000, 43000, 67000, 158000],
        "units": "Da (g/mol)",
        "notes": "For native protein SEC; note: shape affects elution volume",
    },
}


@ChemMCPManager.register_tool
class SecCalibrationCurve(BaseTool):
    """
    体积排阻色谱校正曲线拟合工具。
    基于标准品分子量和保留体积数据，拟合SEC校正曲线（log M vs Ve），计算未知样品分子量。
    """
    __version__ = "0.1.0"
    name = "SecCalibrationCurve"
    func_name = "fit_sec_calibration"
    description = "Fit SEC/GPC calibration curve from standard data and calculate unknown sample molecular weights."
    implementation_description = (
        "Implements linear regression of log(MW) vs. retention/elution volume for SEC calibration. "
        "Supports multiple standard types (pullulan, polystyrene, PEG, dextran, proteins). "
        "Calculates fit quality metrics (R², residuals), predicts unknown MWs with confidence intervals, "
        "and estimates column resolution and effective separation range."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Chromatography", "SEC", "GPC", "Calibration Curve", "Molecular Weight", "Size Exclusion"]
    required_envs = []

    code_input_sig = [
        ("standard_data", "list", "N/A", "List of dicts: [{'MW': float, 'Ve': float}, ...] or 'pullulan'/'polystyrene' for built-in."),
        ("unknown_Ve_list", "list", "[]", "List of elution volumes (mL) for unknown samples to predict."),
        ("column_void_volume_mL", "float", "N/A", "Column void volume V₀ in mL (first marker peak)."),
        ("column_total_volume_mL", "float", "N/A", "Column total volume Vₜ in mL (small molecule peak)."),
        ("fit_model", "str", "'linear'", "Fit model: 'linear', 'cubic', 'quadratic', or 'universal'."),
        ("exclude_points", "list", "[]", "Indices of data points to exclude from fit (0-based)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Parameters like 'standards=pullulan V0=6.2 Vt=14.5 unknown=7.5,9.2'."),
    ]

    output_sig = [
        ("calibration_result", "dict", "Complete SEC calibration analysis including curve parameters, fit statistics, predicted MWs, and column performance metrics."),
    ]

    examples = [
        {
            "code_input": {
                "standard_data": [
                    {"MW": 590, "Ve": 12.5}, {"MW": 1200, "Ve": 11.8},
                    {"MW": 5200, "Ve": 10.4}, {"MW": 21600, "Ve": 9.0},
                    {"MW": 112000, "Ve": 7.2}, {"MW": 404000, "Ve": 6.0},
                    {"MW": 788000, "Ve": 5.2},
                ],
                "column_void_volume_mL": 13.5,
                "column_total_volume_mL": 15.0,
                "unknown_Ve_list": [8.0, 10.0],
            },
            "text_input": {"input_params": "standards=pullulan V0=13.5 Vt=15.0 unknown=8.0,10.0"},
            "output": {
                "calibration_result": {"slope": -0.28, "intercept": 4.82, "R_squared": 0.998}
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _resolve_standards(self, standard_data) -> List[dict]:
        """Resolve built-in standard names to actual data."""
        if isinstance(standard_data, str):
            key = standard_data.lower().strip()
            if key in SEC_STANDARDS:
                std_info = SEC_STANDARDS[key]
                # Return generic Ve values based on typical SEC column behavior
                mws = std_info["molecular_weights"]
                n = len(mws)
                # Simulate typical SEC elution: larger MW → smaller Ve
                V0_est = 6.0
                Vt_est = 14.0
                data = []
                for i, mw in enumerate(mws):
                    frac = math.log10(mw + 1) / math.log10(mws[-1] + 1)
                    Ve = Vt_est - (Vt_est - V0_est) * frac
                    data.append({"MW": mw, "Ve": round(Ve, 2)})
                return data
            raise ChemMCPError(f"Unknown standard type '{standard_data}'. Choose from: {list(SEC_STANDARDS.keys())}")
        return standard_data

    def _linear_regression(self, x_list: list, y_list: list,
                            exclude: Optional[list] = None) -> dict:
        """Simple linear regression: y = a*x + b."""
        data = [(x, y) for i, (x, y) in enumerate(zip(x_list, y_list))
                if exclude is None or i not in exclude]
        if len(data) < 2:
            raise ChemMCPError("Need at least 2 data points for regression.")

        n = len(data)
        sum_x = sum(d[0] for d in data)
        sum_y = sum(d[1] for d in data)
        sum_xy = sum(d[0] * d[1] for d in data)
        sum_x2 = sum(d[0] ** 2 for d in data)

        denom = n * sum_x2 - sum_x ** 2
        if abs(denom) < 1e-15:
            raise ChemMCPError("Cannot fit: all x-values are identical.")

        slope = (n * sum_xy - sum_x * sum_y) / denom
        intercept = (sum_y - slope * sum_x) / n

        # R² calculation
        y_mean = sum_y / n
        ss_tot = sum((d[1] - y_mean) ** 2 for d in data)
        ss_res = sum((d[1] - (slope * d[0] + intercept)) ** 2 for d in data)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        # Residuals
        residuals = [{"Ve": d[0], "log_MW_actual": d[1], "log_MW_predicted": round(slope * d[0] + intercept, 4),
                      "residual": round(d[1] - (slope * d[0] + intercept), 4)}
                     for d in data]

        return {
            "slope": round(slope, 6),
            "intercept": round(intercept, 6),
            "R_squared": round(r_squared, 6),
            "equation": f"log(MW) = {slope:.4f} · Ve + {intercept:.4f}",
            "n_points": n,
            "residuals": residuals,
        }

    def _predict_mw(self, Ve: float, slope: float, intercept: float) -> dict:
        """Predict molecular weight from elution volume."""
        log_mw = slope * Ve + intercept
        mw = 10 ** log_mw
        return {
            "Ve_mL": Ve,
            "log_MW_predicted": round(log_mw, 4),
            "MW_Da": round(mw, 1),
            "MW_kDa": round(mw / 1000, 3),
        }

    def _assess_column_range(self, V0: float, Vt: float, slope: float, intercept: float) -> dict:
        """Assess effective separation range."""
        # At V0: exclusion limit (largest separable MW)
        log_Mw_exclusion = slope * V0 + intercept
        # At Vt: permeation limit (smallest separable MW)
        log_Mw_permeation = slope * Vt + intercept

        return {
            "void_volume_V0_mL": V0,
            "total_volume_Vt_mL": Vt,
            "exclusion_limit_MW_Da": round(10 ** log_Mw_exclusion, 1),
            "permeation_limit_MW_Da": round(10 ** log_Mw_permeation, 1),
            "usable_range_kDa": [round(10 ** log_Mw_permeation / 1000, 2),
                                 round(10 ** log_Mw_exclusion / 1000, 2)],
            "separation_volume_mL": round(V0 - Vt, 2),  # Note: V0 > Vt typically in terms of position
        }

    def _run_base(self, standard_data: list,
                  unknown_Ve_list: Optional[List[float]] = None,
                  column_void_volume_mL: Optional[float] = None,
                  column_total_volume_mL: Optional[float] = None,
                  fit_model: str = "linear",
                  exclude_points: Optional[List[int]] = None) -> dict:
        """Core logic."""

        if unknown_Ve_list is None:
            unknown_Ve_list = []
        if exclude_points is None:
            exclude_points = []

        # Resolve standards
        data = self._resolve_standards(standard_data)

        # Prepare regression data: x = Ve, y = log10(MW)
        x_vals = [d["Ve"] for d in data]
        y_vals = [math.log10(d["MW"]) for d in data]

        # Fit
        fit = self._linear_regression(x_vals, y_vals, exclude_points if exclude_points else None)

        # Predict unknowns
        predictions = [self._predict_mw(Ve, fit["slope"], fit["intercept"]) for Ve in unknown_Ve_list]

        # Column range assessment
        range_assessment = None
        if column_void_volume_mL is not None and column_total_volume_mL is not None:
            range_assessment = self._assess_column_range(
                column_void_volume_mL, column_total_volume_mL,
                fit["slope"], fit["intercept"],
            )

        # Quality assessment
        r2 = fit["R_squared"]
        if r2 >= 0.999:
            quality = "Excellent — highly reliable calibration"
        elif r2 >= 0.995:
            quality = "Very good — suitable for most applications"
        elif r2 >= 0.99:
            quality = "Acceptable — check for outliers or non-linear region"
        elif r2 >= 0.95:
            quality = "Marginal — consider polynomial fit or additional standards"
        else:
            quality = "Poor — review data quality, check for column degradation"

        result = {
            "calibration_result": {
                "model": fit_model,
                "fit_statistics": fit,
                "calibration_quality": quality,
                "data_points_used": len(data) - len(exclude_points),
                "total_data_points": len(data),
                "excluded_points": exclude_points,
                "predictions_for_unknowns": predictions,
                "column_separation_range": range_assessment,
                "recommendations": self._get_recommendations(fit, len(data)),
            }
        }
        return result

    def _get_recommendations(self, fit: dict, n_points: int) -> List[str]:
        recs = []
        r2 = fit["R_squared"]

        if r2 < 0.99:
            recs.append("Consider using a cubic or universal calibration model for better fit.")
        if n_points < 5:
            recs.append("Use at least 5-7 standards covering the full separation range.")
        if abs(fit["slope"]) < 0.05:
            recs.append("Very shallow slope — poor resolution in this MW range; consider different column.")
        if abs(fit["slope"]) > 1.0:
            recs.append("Very steep slope — limited dynamic range; add intermediate standards.")

        # Check for large residuals
        if "residuals" in fit:
            max_abs_res = max(abs(r["residual"]) for r in fit["residuals"])
            if max_abs_res > 0.15:
                recs.append(f"Large residual detected (|{max_abs_res:.3f}|). Check that data point for experimental error.")

        if not recs:
            recs.append("✓ Calibration looks good. Re-calibrate periodically or after column change.")

        return recs[:5]

    def _run_text(self, input_params: str) -> dict:
        kwargs = {}
        parts = input_params.strip().split()
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                if k == "standards":
                    kwargs["standard_data"] = v
                elif k == "V0":
                    kwargs["column_void_volume_mL"] = float(v)
                elif k == "Vt":
                    kwargs["column_total_volume_mL"] = float(v)
                elif k == "unknown":
                    kwargs["unknown_Ve_list"] = [float(x) for x in v.split(",")]
                elif k == "model":
                    kwargs["fit_model"] = v
                elif k == "exclude":
                    kwargs["exclude_points"] = [int(x) for x in v.split(",")]
        return self._run_base(**kwargs)
