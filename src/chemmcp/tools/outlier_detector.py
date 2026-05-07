import json
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ---- Dixon Q critical values table (for common α and n=3..30) ----
# Format: DIXON_Q_CRITICAL[alpha][n] = critical value
_DIXON_Q_CRITICAL: Dict[float, Dict[int, float]] = {
    0.10: {
        3: 0.886, 4: 0.679, 5: 0.557, 6: 0.482, 7: 0.434,
        8: 0.479, 9: 0.441, 10: 0.409, 11: 0.517, 12: 0.490,
        13: 0.467, 14: 0.449, 15: 0.435, 16: 0.423, 17: 0.413,
        18: 0.404, 19: 0.396, 20: 0.389, 21: 0.383, 22: 0.378,
        23: 0.373, 24: 0.370, 25: 0.366, 26: 0.363, 27: 0.360,
        28: 0.358, 29: 0.355, 30: 0.353,
    },
    0.05: {
        3: 0.941, 4: 0.765, 5: 0.642, 6: 0.560, 7: 0.507,
        8: 0.554, 9: 0.512, 10: 0.477, 11: 0.576, 12: 0.546,
        13: 0.521, 14: 0.501, 15: 0.485, 16: 0.472, 17: 0.460,
        18: 0.450, 19: 0.441, 20: 0.434, 21: 0.427, 22: 0.421,
        23: 0.416, 24: 0.411, 25: 0.407, 26: 0.403, 27: 0.399,
        28: 0.396, 29: 0.393, 30: 0.390,
    },
    0.01: {
        3: 0.988, 4: 0.889, 5: 0.780, 6: 0.698, 7: 0.637,
        8: 0.683, 9: 0.635, 10: 0.597, 11: 0.679, 12: 0.642,
        13: 0.615, 14: 0.592, 15: 0.572, 16: 0.555, 17: 0.540,
        18: 0.527, 19: 0.516, 20: 0.505, 21: 0.496, 22: 0.488,
        23: 0.480, 24: 0.473, 25: 0.467, 26: 0.461, 27: 0.456,
        28: 0.451, 29: 0.447, 30: 0.443,
    },
}


def _grubbs_critical(n: int, alpha: float) -> float:
    """
    Approximate Grubbs critical value using the formula:
    G_crit ≈ t_(α/(2n), n-2) × √((n-1)/(n-2))
    where t is the Student's t-distribution critical value (approximated).
    Uses the Stephens approximation formula.
    """
    from math import sqrt as _sqrt
    from math import sqrt
    if n < 3:
        return float('inf')
    # Approximate using formula from ISO 5725 / standard tables
    # G = (n-1)/sqrt(n) * sqrt(t²/(n-2+t²))
    # We approximate t-critical using inverse of CDF approximation
    p = 1 - alpha / (2 * n)  # two-sided equivalent for max outlier
    # Use approximate t-value via normal approximation for large n
    # For small n, use lookup-like approximations
    if alpha == 0.05:
        _t_approx = {3:1.153,4:1.463,5:1.672,6:1.822,7:1.938,8:2.032,
                     9:2.110,10:2.176,11:2.232,12:2.280,13:2.322,
                     14:2.360,15:2.393,16:2.423,17:2.451,18:2.476,
                     19:2.499,20:2.520,22:2.558,25:2.605,30:2.663,
                     35:2.711,40:2.750,50:2.812}
    elif alpha == 0.01:
        _t_approx = {3:1.155,4:1.496,5:1.749,6:1.944,7:2.097,8:2.220,
                     9:2.323,10:2.410,11:2.485,12:2.550,13:2.607,
                     14:2.658,15:2.703,16:2.743,17:2.779,18:2.812,
                     19:2.842,20:2.869,22:2.916,25:2.974,30:3.044}
    elif alpha == 0.10:
        _t_approx = {3:1.148,4:1.425,5:1.602,6:1.729,7:1.827,8:1.906,
                     9:1.972,10:2.027,11:2.075,12:2.116,13:2.152,
                     14:2.184,15:2.212,16:2.237,17:2.260,18:2.280,
                     19:2.299,20:2.316,22:2.346,25:2.385,30:2.433}
    else:
        # Approximate via normal distribution
        import math as _m
        z = _sqrt(2) * _erfinv(1 - alpha / (2 * n)) if n > 2 else 3.0
        return _sqrt((n - 1) ** 2 / n) * z / _sqrt(n - 2 + z**2)

    if n in _t_approx:
        t_val = _t_approx[n]
    else:
        # Linear interpolation for intermediate n
        keys = sorted(k for k in _t_approx if k < n)
        if keys:
            k1 = keys[-1]
            k2 = min((k for k in _t_approx if k > k1), default=k1)
            if k2 > k1:
                frac = (n - k1) / (k2 - k1)
                t_val = _t_approx[k1] + frac * (_t_approx[k2] - _t_approx[k1])
            else:
                t_val = _t_approx[k1]
        else:
            t_val = list(_t_approx.values())[-1]

    g = ((n - 1) / _sqrt(n)) * _sqrt(t_val ** 2 / (n - 2 + t_val ** 2))
    return g


def _erfinv(x: float) -> float:
    """Approximate inverse error function."""
    import math as _m
    a = [0.000000000000000015249067145857658684, 0.00000000000000495614179399820814667,
         0.00000000000167389250722874903875, 0.00000000013994121059807494,
         0.00000000538078556214343031, 0.00000012486056894513651,
         0.00000173411061668058963, 0.00001552334358640781,
         0.00008844190098916769, 0.00034008264799896471,
         0.00085047020215134475, 0.0013337298837011033,
         0.0012505592878743723, 0.00042849203027709052]
    w = -_m.log((1 - x) * (1 + x), 2)
    if w < 1:
        t = w * (a[0] + w * (a[1] + w * a[2]))
    elif w < 4:
        t = w * (a[3] + w * (a[4] + w * (a[5] + w * a[6])))
    else:
        t = w * (a[7] + w * (a[8] + w * (a[9] + w * (a[10] + w * (a[11] + w * (a[12] + w * a[13]))))))
    p = _m.sqrt(t)
    result = p - (a[0] + a[1]*t + a[2]*t*t) / (1 + 2*a[0]*p + (2*a[1]+a[0]*a[0])*t)
    return result if x >= 0 else -result


@ChemMCPManager.register_tool
class OutlierDetector(BaseTool):
    """
    异常值检验工具。
    支持 Grubbs 检验和 Dixon Q 检验，支持顺序剔除多个异常值。
    """
    __version__ = "0.1.0"
    name = "OutlierDetector"
    func_name = "detect_outliers"
    description = "Detect outliers using Grubbs' test or Dixon's Q test with support for sequential removal."
    implementation_description = "Implements Grubbs' test (G-statistic based on mean and SD) and Dixon's Q test (gap/range ratio). Supports sequential removal of multiple outliers; Dixon only allows α∈{0.1,0.05,0.01} and sample size 3–30."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Outlier", "Grubbs", "Dixon Q", "Statistics", "QA/QC"]
    required_envs = []

    code_input_sig = [
        ("data", "list", "N/A", "List of numerical values."),
        ("method", "str", "grubbs", "Test method: 'grubbs' or 'dixon'."),
        ("alpha", "float", "0.05", "Significance level (Dixon: 0.01, 0.05, or 0.10)."),
        ("tail", "str", "two-sided", "Tail type: 'two-sided', 'high', or 'low'."),
        ("max_outliers", "int", "1", "Maximum number of outliers to detect sequentially."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters."),
    ]

    output_sig = [
        ("method", "str", "Test method used."),
        ("alpha", "float", "Significance level used."),
        ("tail", "str", "Tail type used."),
        ("sample_size", "int", "Sample size."),
        ("detected", "list", "List of detected outliers with details (value, index, direction, statistic, threshold)."),
        ("outliers", "list", "List of outlier values."),
        ("clean_data", "list", "Data with outliers removed."),
        ("interpretation", "str", "Natural language interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "data": [12.5, 12.7, 12.8, 12.6, 12.7, 12.5, 12.4, 18.3, 12.6, 12.7],
                "method": "grubbs",
                "alpha": 0.05,
                "max_outliers": 2,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"outliers": [18.3], "clean_data": [12.5, 12.7, 12.8, 12.6, 12.7, 12.5, 12.4, 12.6, 12.7]},
        },
        {
            "code_input": {
                "data": [0.225, 0.227, 0.230, 0.228, 0.226, 0.295],
                "method": "dixon",
                "alpha": 0.05,
            },
            "text_input": {"params_str": "see code input"},
            "output": {"outliers": [0.295]},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    @staticmethod
    def _mean(data: List[float]) -> float:
        return sum(data) / len(data)

    @staticmethod
    def _std(data: List[float], ddof: int = 1) -> float:
        n = len(data)
        m = sum(data) / n
        return math.sqrt(sum((x - m) ** 2 for x in data) / (n - ddof))

    def _grubbs_test(self, data: List[float], alpha: float, tail: str) -> dict:
        """Perform single Grubbs test on data."""
        n = len(data)
        if n < 3:
            return {"is_outlier": False, "statistic": 0, "critical_value": float('inf'), "outlier_value": None, "index": None}

        mean = self._mean(data)
        std = self._std(data)

        if std < 1e-15:
            return {"is_outlier": False, "statistic": 0, "critical_value": float('inf'), "outlier_value": None, "index": None}

        if tail == "high":
            idx = data.index(max(data))
            G = (max(data) - mean) / std
            val = max(data)
        elif tail == "low":
            idx = data.index(min(data))
            G = (mean - min(data)) / std
            val = min(data)
        else:  # two-sided
            if abs(max(data) - mean) >= abs(min(data) - mean):
                idx = data.index(max(data))
                G = abs(max(data) - mean) / std
                val = max(data)
            else:
                idx = data.index(min(data))
                G = abs(mean - min(data)) / std
                val = min(data)

        crit = _grubbs_critical(n, alpha)
        is_outl = G > crit

        return {
            "is_outlier": is_outl,
            "statistic": round(G, 6),
            "critical_value": round(crit, 6),
            "outlier_value": val,
            "index": idx,
            "direction": "high" if val == max(data) else "low" if val == min(data) else "unknown",
        }

    def _dixon_test(self, data: List[float], alpha: float, tail: str) -> dict:
        """Perform single Dixon Q test on data."""
        n = len(data)
        if n < 3 or n > 30:
            raise ChemMCPError(f"Dixon Q test requires sample size between 3 and 30. Got n={n}.")
        if alpha not in (0.01, 0.05, 0.10):
            raise ChemMCPError(f"Dixon Q test only supports α ∈ {{0.01, 0.05, 0.10}}. Got α={alpha}.")

        sorted_data = sorted(data)
        range_val = sorted_data[-1] - sorted_data[0]

        if range_val < 1e-15:
            return {"is_outlier": False, "statistic": 0, "critical_value": _DIXON_Q_CRITICAL[alpha][n],
                    "outlier_value": None, "index": None}

        # Choose appropriate Q formula based on n
        if tail == "high":
            # Test maximum value
            if n <= 7:
                Q = (sorted_data[-1] - sorted_data[-2]) / range_val
            elif n <= 10:
                Q = (sorted_data[-1] - sorted_data[-2]) / (sorted_data[-1] - sorted_data[1])
            elif n <= 13:
                Q = (sorted_data[-1] - sorted_data[-2]) / (sorted_data[-1] - sorted_data[2])
            else:
                Q = (sorted_data[-1] - sorted_data[-3]) / (sorted_data[-1] - sorted_data[2])
            val = sorted_data[-1]
            idx = data.index(val)
        elif tail == "low":
            if n <= 7:
                Q = (sorted_data[1] - sorted_data[0]) / range_val
            elif n <= 10:
                Q = (sorted_data[1] - sorted_data[0]) / (sorted_data[-2] - sorted_data[0])
            elif n <= 13:
                Q = (sorted_data[1] - sorted_data[0]) / (sorted_data[-3] - sorted_data[0])
            else:
                Q = (sorted_data[2] - sorted_data[0]) / (sorted_data[-3] - sorted_data[0])
            val = sorted_data[0]
            idx = data.index(val)
        else:  # two-sided: test both ends, use larger Q
            # High end
            if n <= 7:
                Q_high = (sorted_data[-1] - sorted_data[-2]) / range_val
                Q_low = (sorted_data[1] - sorted_data[0]) / range_val
            elif n <= 10:
                Q_high = (sorted_data[-1] - sorted_data[-2]) / (sorted_data[-1] - sorted_data[1])
                Q_low = (sorted_data[1] - sorted_data[0]) / (sorted_data[-2] - sorted_data[0])
            elif n <= 13:
                Q_high = (sorted_data[-1] - sorted_data[-2]) / (sorted_data[-1] - sorted_data[2])
                Q_low = (sorted_data[1] - sorted_data[0]) / (sorted_data[-3] - sorted_data[0])
            else:
                Q_high = (sorted_data[-1] - sorted_data[-3]) / (sorted_data[-1] - sorted_data[2])
                Q_low = (sorted_data[2] - sorted_data[0]) / (sorted_data[-3] - sorted_data[0])

            if Q_high >= Q_low:
                Q = Q_high
                val = sorted_data[-1]
                idx = data.index(val)
            else:
                Q = Q_low
                val = sorted_data[0]
                idx = data.index(val)

        crit = _DIXON_Q_CRITICAL[alpha].get(n, 0.5)
        is_outl = Q > crit

        return {
            "is_outlier": is_outl,
            "statistic": round(Q, 6),
            "critical_value": round(crit, 6),
            "outlier_value": val,
            "index": idx,
            "direction": "high" if val == max(data) else "low" if val == min(data) else "unknown",
        }

    def _run_base(
        self,
        data: List[float],
        method: str = "grubbs",
        alpha: float = 0.05,
        tail: str = "two-sided",
        max_outliers: int = 1,
    ) -> dict:
        """Core logic: detect outliers with sequential removal."""
        if not data:
            raise ChemMCPError("Data list is empty.")
        if len(data) < 3:
            raise ChemMCPError(f"Need at least 3 data points. Got {len(data)}.")

        meth = method.lower().strip()
        tail_type = tail.lower().strip()
        max_out = max(1, min(max_outliers, len(data) - 2))

        working_data = list(data)
        original_indices = list(range(len(data)))
        all_detected: List[Dict[str, Any]] = []
        step = 0

        while step < max_out and len(working_data) >= 3:
            step += 1
            if meth == "grubbs":
                result = self._grubbs_test(working_data, alpha, tail_type)
            elif meth == "dixon":
                if len(working_data) > 30:
                    break  # Dixon doesn't support n > 30
                result = self._dixon_test(working_data, alpha, tail_type)
            else:
                raise ChemMCPError(f"Unknown method '{method}'. Use 'grubbs' or 'dixon'.")

            if not result["is_outlier"]:
                break

            outl_val = result["outlier_value"]
            orig_idx = original_indices[working_data.index(outl_val)] if outl_val in working_data else -1

            detection = {
                "step": step,
                "value": outl_val,
                "original_index": orig_idx,
                "direction": result.get("direction", "unknown"),
                "statistic": result["statistic"],
                "threshold": result["critical_value"],
            }
            all_detected.append(detection)

            # Remove outlier for next iteration
            if outl_val in working_data:
                oi = working_data.index(outl_val)
                working_data.pop(oi)
                original_indices.pop(oi)

        outliers_list = [d["value"] for d in all_detected]

        # Build interpretation
        if all_detected:
            details = "; ".join(
                f"x[{d['original_index']}]={d['value']} ({d['direction']}, G={d['statistic']:.4f}>{d['threshold']:.4f})"
                for d in all_detected
            )
            interpretation = (
                f"{meth.capitalize()} test (α={alpha}, {tail_type}): "
                f"detected {len(all_detected)} outlier(s): {details}. "
                f"These values are statistically inconsistent with the rest at the {alpha} significance level."
            )
        else:
            interpretation = (
                f"{meth.capitalize()} test (α={alpha}, {tail_type}): "
                f"no outliers detected. All values are statistically consistent."
            )

        logger.info(f"Outlier detection ({meth}, α={alpha}): found {len(all_detected)} outlier(s)")
        return {
            "method": meth,
            "alpha": alpha,
            "tail": tail_type,
            "sample_size": len(data),
            "detected": all_detected,
            "outliers": outliers_list,
            "clean_data": working_data,
            "interpretation": interpretation,
        }

    def _run_text(self, params_str: str) -> dict:
        try:
            kwargs = json.loads(params_str.strip())
        except json.JSONDecodeError:
            raise ChemMCPError("Invalid JSON input.")
        return self._run_base(**kwargs)
