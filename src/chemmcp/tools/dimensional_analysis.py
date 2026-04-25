import logging
import re
from fractions import Fraction

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

BASE_DIMENSIONS = {"L": "Length", "M": "Mass", "T": "Time", "I": "Electric Current",
                    "O": "Temperature", "N": "Amount of Substance", "J": "Luminous Intensity"}

QUANTITY_DIMENSIONS = {
    "length": {"L": 1}, "area": {"L": 2}, "volume": {"L": 3},
    "velocity": {"L": 1, "T": -1}, "speed": {"L": 1, "T": -1},
    "acceleration": {"L": 1, "T": -2},
    "mass": {"M": 1}, "density": {"M": 1, "L": -3},
    "force": {"M": 1, "L": 1, "T": -2},
    "pressure": {"M": 1, "L": -1, "T": -2}, "stress": {"M": 1, "L": -1, "T": -2},
    "energy": {"M": 1, "L": 2, "T": -2}, "work": {"M": 1, "L": 2, "T": -2},
    "heat": {"M": 1, "L": 2, "T": -2}, "power": {"M": 1, "L": 2, "T": -3},
    "momentum": {"M": 1, "L": 1, "T": -1},
    "frequency": {"T": -1}, "angular_velocity": {"T": -1},
    "temperature": {"O": 1},
    "electric_charge": {"I": 1, "T": 1}, "current": {"I": 1},
    "voltage": {"M": 1, "L": 2, "T": -3, "I": -1},
    "resistance": {"M": 1, "L": 2, "T": -3, "I": -2},
    "distance": {"L": 1}, "time_q": {"T": 1},
    "amount_of_substance": {"N": 1}, "molar_mass": {"M": 1, "N": -1},
    "molarity": {"N": 1, "L": -3}, "concentration": {"N": 1, "L": -3},
    "rate_constant_1": {"T": -1}, "decay_constant": {"T": -1},
    "half_life_dim": {"T": 1}, "activity_bq": {"T": -1},
    "angle": {}, "solid_angle": {}, "dimensionless": {},
}


@ChemMCPManager.register_tool
class DimensionalAnalysis(BaseTool):
    __version__ = "0.1.0"
    name = "DimensionalAnalysis"
    func_name = "analyze_dimensions"
    description = "Perform dimensional analysis: query dimensions, check equation consistency, derive unknown dimensions."
    implementation_description = "Uses built-in database of physical quantities with MLTO system dimensions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Dimensional Analysis", "Units", "Physics", "Chemistry"]
    required_envs = []

    code_input_sig = [
        ("operation", "str", "N/A", "Operation: 'query', 'check', or 'derive'."),
        ("quantity", "str", "", "Quantity name for 'query'."),
        ("left_side", "str", "", "Space-separated left-side quantities for 'check'."),
        ("right_side", "str", "", "Space-separated right-side quantities for 'check'."),
        ("expression", "str", "", "Expression string for 'derive', e.g., 'force * velocity / area'."),
    ]
    text_input_sig = [
        ("query_str", "str", "N/A", "Query string, e.g., 'query pressure' or 'check energy force distance'."),
    ]
    output_sig = [
        ("operation", "str", "Operation performed."),
        ("result", "dict", "Detailed analysis result."),
    ]
    examples = [
        {
            "code_input": {"operation": "query", "quantity": "pressure", "left_side": "", "right_side": "", "expression": ""},
            "text_input": {"query_str": "query pressure"},
            "output": {
                "operation": "query",
                "result": {
                    "quantity": "pressure",
                    "dimensional_formula": "M.L^-1.T^-2",
                    "is_dimensionless": False,
                },
            }
        },
        {
            "code_input": {"operation": "check", "quantity": "", "left_side": "energy", "right_side": "force distance", "expression": ""},
            "text_input": {"query_str": "check energy force distance"},
            "output": {
                "operation": "check",
                "result": {
                    "consistent": True,
                    "message": "Equation is dimensionally consistent.",
                },
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _format_dims(self, dims):
        parts = []
        superscript_map = {"0":"0","1":"","2":"2","3":"3","4":"4","5":"5","6":"6","7":"7","8":"8","9":"9","-":"-"}
        for d in ["M", "L", "T", "I", "O", "N", "J"]:
            if d in dims and dims[d] != 0:
                exp = dims[d]
                if exp == 1:
                    parts.append(d)
                else:
                    es = "".join(superscript_map.get(c, str(c)) for c in str(exp))
                    parts.append(f"{d}^{es}")
        return ".".join(parts) if parts else "dimensionless"

    def _multiply_dims(self, dim_list):
        result = {}
        for dims in dim_list:
            for d, exp in dims.items():
                result[d] = result.get(d, 0) + exp
        return {d: e for d, e in result.items() if e != 0}

    def _find_quantity(self, name):
        key = name.strip().lower()
        if key in QUANTITY_DIMENSIONS:
            return QUANTITY_DIMENSIONS[key]
        matches = [k for k in QUANTITY_DIMENSIONS if key in k or k in key]
        if len(matches) == 1:
            return QUANTITY_DIMENSIONS[matches[0]]
        elif len(matches) > 1:
            raise ChemMCPError(f"Ambiguous quantity '{name}'. Matches: {matches}")
        raise ChemMCPError(f"Unknown quantity '{name}'.")

    def _run_base(self, operation, quantity="", left_side="", right_side="", expression=""):
        op = operation.strip().lower()
        if op == "query":
            if not quantity:
                raise ChemMCPError("'query' requires a 'quantity' parameter.")
            dims = self._find_quantity(quantity)
            fmt = self._format_dims(dims)
            return {"operation": "query", "result": {
                "quantity": quantity,
                "dimensional_formula": fmt,
                "base_dims": dims,
                "is_dimensionless": len(dims) == 0,
            }}
        elif op == "check":
            if not left_side or not right_side:
                raise ChemMCPError("'check' requires both 'left_side' and 'right_side'.")
            left_dims = self._multiply_dims([self._find_quantity(q) for q in left_side.split()])
            right_dims = self._multiply_dims([self._find_quantity(q) for q in right_side.split()])
            consistent = left_dims == right_dims
            return {"operation": "check", "result": {
                "consistent": consistent,
                "message": "Consistent" if consistent else f"Not consistent: LHS={self._format_dims(left_dims)} != RHS={self._format_dims(right_dims)}",
                "left_dims": self._format_dims(left_dims),
                "right_dims": self._format_dims(right_dims),
            }}
        elif op == "derive":
            if not expression:
                raise ChemMCPError("'derive' requires an 'expression' parameter.")
            tokens = expression.replace("*", " * ").replace("/", " / ").split()
            current_dims = {}
            op_type = "*"
            for token in tokens:
                if token == "*":
                    op_type = "*"
                elif token == "/":
                    op_type = "/"
                else:
                    try:
                        qdims = self._find_quantity(token)
                        if op_type == "*":
                            current_dims = self._multiply_dims([current_dims, qdims]) if current_dims else dict(qdims)
                        else:
                            result = dict(current_dims)
                            for d, e in qdims.items():
                                result[d] = result.get(d, 0) - e
                            current_dims = {d: e for d, e in result.items() if e != 0}
                    except ChemMCPError:
                        pass
            return {"operation": "derive", "result": {
                "expression": expression,
                "dimensional_formula": self._format_dims(current_dims),
                "base_dims": current_dims,
                "is_dimensionless": len(current_dims) == 0,
            }}
        else:
            raise ChemMCPError(f"Unknown operation '{op}'.")

    def _run_text(self, query_str):
        parts = query_str.strip().split()
        if not parts:
            raise ChemMCPError("Empty query.")
        op = parts[0].lower()
        if op == "query" and len(parts) >= 2:
            return self._run_base("query", quantity=" ".join(parts[1:]))
        elif op == "check" and len(parts) >= 4:
            rest = " ".join(parts[1:])
            if "=" in rest:
                lhs, rhs = rest.split("=", 1)
                return self._run_base("check", left_side=lhs.strip(), right_side=rhs.strip())
            else:
                all_q = rest.split()
                mid = len(all_q) // 2
                return self._run_base("check", left_side=" ".join(all_q[:mid]), right_side=" ".join(all_q[mid:]))
        elif op == "derive" and len(parts) >= 2:
            return self._run_base("derive", expression=" ".join(parts[1:]))
        else:
            return self._run_base("query", quantity=query_str.strip())
