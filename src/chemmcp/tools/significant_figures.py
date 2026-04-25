import logging
import math
from decimal import Decimal, ROUND_HALF_UP

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

SUPERSCRIPT_MAP = {
    "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3",
    "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077",
    "8": "\u2078", "9": "\u2079", "-": "\u207b",
}


@ChemMCPManager.register_tool
class SignificantFigures(BaseTool):
    __version__ = "0.1.0"
    name = "SignificantFigures"
    func_name = "handle_significant_figures"
    description = "Handle significant figures: count sig figs, round to N sig figs, convert to scientific notation."
    implementation_description = "Implements IUPAC/chemistry significant figure conventions for counting, rounding, and arithmetic."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Significant Figures", "Precision", "Data Analysis"]
    required_envs = []

    code_input_sig = [
        ("operation", "str", "N/A", "Operation: 'count', 'round', 'scientific', 'add', 'subtract', 'multiply', 'divide'."),
        ("value_a", "str", "N/A", "First value as string."),
        ("value_b", "str", "", "Second value (for arithmetic)."),
        ("n_sig_figs", "int", "4", "Number of significant figures for rounding."),
    ]
    text_input_sig = [
        ("params_str", "str", "N/A", "Space-separated: operation value_a [value_b] [n_sig_figs]."),
    ]
    output_sig = [
        ("operation", "str", "Operation performed."),
        ("result", "str", "Result value as string."),
        ("sig_figs", "int", "Number of significant figures in result."),
        ("explanation", "str", "Step-by-step explanation."),
    ]
    examples = [
        {
            "code_input": {"operation": "count", "value_a": "0.004500", "value_b": "", "n_sig_figs": 4},
            "text_input": {"params_str": "count 0.004500"},
            "output": {"operation": "count", "result": "0.004500", "sig_figs": 4,
                       "explanation": "Leading zeros not significant; trailing zeros after decimal are significant."},
        },
        {
            "code_input": {"operation": "round", "value_a": "3.14159265", "value_b": "", "n_sig_figs": 4},
            "text_input": {"params_str": "round 3.14159265 4"},
            "output": {"operation": "round", "result": "3.142", "sig_figs": 4,
                       "explanation": "Rounded to 4 significant figures."},
        },
        {
            "code_input": {"operation": "multiply", "value_a": "3.5", "value_b": "4.20", "n_sig_figs": 4},
            "text_input": {"params_str": "multiply 3.5 4.20"},
            "output": {"operation": "multiply", "result": "14.7", "sig_figs": 2,
                       "explanation": "3.5(2 sf) x 4.20(3 sf) -> rounded to 2 sf."},
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _count_sig_figs(self, value_str):
        s = value_str.strip()
        try:
            d = Decimal(s)
        except Exception:
            raise ChemMCPError(f"Invalid number format: '{s}'")
        if "E" in s.upper():
            mantissa = s.upper().split("E")[0]
            if "." in mantissa:
                mantissa = mantissa.rstrip("0").rstrip(".")
            else:
                mantissa = mantissa.lstrip("0")
            return len(mantissa.replace(".", ""))
        if s.startswith("-") or s.startswith("+"):
            s = s[1:]
        if Decimal(s) == 0:
            return 1
        parts = s.split(".")
        raw_int = parts[0].lstrip("0")
        int_part = raw_int if raw_int else ""  # empty if number < 1
        frac_part = parts[1].rstrip("0") if len(parts) > 1 and parts[1] else ""
        combined = int_part + frac_part
        return len(combined) if combined else 1

    def _round_to_sig_figs(self, value_str, n):
        d = Decimal(value_str.strip())
        if d == 0:
            return "0" + ("." + "0" * (n-1)) if n > 1 else "0"
        sign = -1 if d < 0 else 1
        d_abs = abs(d)
        if d_abs == 0:
            return "0"
        magnitude = int(math.floor(math.log10(float(d_abs))))
        factor = Decimal(10) ** (magnitude - n + 1)
        rounded = (d_abs / factor).quantize(Decimal(1), rounding=ROUND_HALF_UP) * factor
        result_str = f"{rounded:f}"
        if "." in result_str:
            int_r, frac_r = result_str.split(".")
            if frac_r.rstrip("0"):
                result_str = int_r + "." + frac_r.rstrip("0")
            else:
                result_str = int_r
        if sign < 0:
            result_str = "-" + result_str
        return result_str

    def _to_scientific(self, value_str):
        d = Decimal(value_str.strip())
        if d == 0:
            return "0 x 10^0"
        exp = int(math.floor(math.log10(abs(float(d)))))
        coeff = float(d) / (10 ** exp)
        sf = self._count_sig_figs(value_str)
        coeff_rounded = round(coeff, max(sf - 1, 1))
        superscript = "".join(SUPERSCRIPT_MAP.get(c, c) for c in str(abs(exp)))
        return f"{coeff_rounded} x 10{superscript}"

    def _run_base(self, operation, value_a, value_b="", n_sig_figs=4):
        op = operation.strip().lower()
        if op == "count":
            n = self._count_sig_figs(value_a)
            return {"operation": "count", "result": value_a, "sig_figs": n,
                    "explanation": f"'{value_a}' has {n} significant figure(s)."}
        elif op == "round":
            r = self._round_to_sig_figs(value_a, n_sig_figs)
            return {"operation": "round", "result": r, "sig_figs": min(n_sig_figs, self._count_sig_figs(r)),
                    "explanation": f"'{value_a}' rounded to {n_sig_figs} sf = '{r}'."}
        elif op == "scientific":
            sci = self._to_scientific(value_a)
            return {"operation": "scientific", "result": sci, "sig_figs": self._count_sig_figs(value_a),
                    "explanation": f"'{value_a}' in scientific notation: {sci}."}
        elif op in ("add", "subtract"):
            if not value_b:
                raise ChemMCPError(f"'{op}' requires two values.")
            da = Decimal(value_a)
            db = Decimal(value_b)
            raw_result = da + db if op == "add" else da - db
            dec_a = self._decimal_places(value_a)
            dec_b = self._decimal_places(value_b)
            min_dec = min(dec_a, dec_b)
            result_str = f"{float(raw_result):.{min_dec}f}" if min_dec >= 0 else str(int(round(float(raw_result))))
            return {"operation": op, "result": result_str, "sig_figs": self._count_sig_figs(result_str),
                    "explanation": f"{value_a} ({dec_a} dp) {op} {value_b} ({dec_b} dp) = {result_str}."}
        elif op in ("multiply", "divide"):
            if not value_b:
                raise ChemMCPError(f"'{op}' requires two values.")
            da = float(value_a)
            db = float(value_b)
            sf_a = self._count_sig_figs(value_a)
            sf_b = self._count_sig_figs(value_b)
            raw = da * db if op == "multiply" else (da / db if db != 0 else float("inf"))
            min_sf = min(sf_a, sf_b)
            result_str = self._round_to_sig_figs(str(raw), min_sf)
            return {"operation": op, "result": result_str, "sig_figs": min_sf,
                    "explanation": f"{value_a}({sf_a} sf) {op} {value_b}({sf_b} sf) = {result_str} ({min_sf} sf)."}
        else:
            raise ChemMCPError(f"Unknown operation '{op}'.")

    def _decimal_places(self, value_str):
        s = value_str.strip()
        if "." in s:
            return len(s.split(".")[1].rstrip("0") or "0")
        return 0

    def _run_text(self, params_str):
        parts = params_str.strip().split()
        if len(parts) < 2:
            raise ChemMCPError("Need at least 2 params.")
        kw = {"operation": parts[0], "value_a": parts[1]}
        if len(parts) > 2:
            try:
                kw["n_sig_figs"] = int(parts[2])
            except ValueError:
                kw["value_b"] = parts[2]
        if len(parts) > 3:
            try:
                kw["n_sig_figs"] = int(parts[3])
            except ValueError:
                pass
        return self._run_base(**kw)
