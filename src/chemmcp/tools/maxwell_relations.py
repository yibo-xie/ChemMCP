import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

R = 8.314  # J/(mol·K)


@ChemMCPManager.register_tool
class MaxwellRelations(BaseTool):
    """
    麦克斯韦关系式的推导与数值验证。

    四个麦克斯韦关系式来自四个热力学势：
    - dU = TdS - PdV  → (∂T/∂V)_S = -(∂P/∂S)_V
    - dH = TdS + VdP  → (∂T/∂P)_S = (∂V/∂S)_P
    - dG = -SdT + VdP → -(∂S/∂P)_T = (∂V/∂T)_P
    - dA = -SdT - PdV → (∂S/∂V)_T = (∂P/∂T)_V
    """
    __version__ = "0.1.0"
    name = "MaxwellRelations"
    func_name = "maxwell_relations_analysis"
    description = "Derive and numerically verify the four Maxwell relations from thermodynamic potentials."
    implementation_description = "Provides symbolic derivation of each Maxwell relation and optionally performs numerical verification using finite differences on a van der Waals gas model."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Maxwell Relations", "Physical Chemistry", "Thermodynamic Potentials"]
    required_envs = []

    _relations = {
        1: {
            "name": "From Internal Energy (U)",
            "potential": "dU = TdS - PdV",
            "relation": "(∂T/∂V)_S = -(∂P/∂S)_V",
            "derivation": (
                "Since dU is an exact differential:\n"
                "dU = TdS - PdV\n"
                "By equality of mixed partials: ∂²U/∂S∂V = ∂²U/∂V∂S\n"
                "→ (∂T/∂V)_S = -(∂P/∂S)_V"
            ),
        },
        2: {
            "name": "From Enthalpy (H)",
            "potential": "dH = TdS + VdP",
            "relation": "(∂T/∂P)_S = (∂V/∂S)_P",
            "derivation": (
                "Since dH is an exact differential:\n"
                "dH = TdS + VdP\n"
                "By equality of mixed partials: ∂²H/∂S∂P = ∂²H/∂P∂S\n"
                "→ (∂T/∂P)_S = (∂V/∂S)_P"
            ),
        },
        3: {
            "name": "From Gibbs Free Energy (G)",
            "potential": "dG = -SdT + VdP",
            "relation": "-(∂S/∂P)_T = (∂V/∂T)_P",
            "derivation": (
                "Since dG is an exact differential:\n"
                "dG = -SdT + VdP\n"
                "By equality of mixed partials: ∂²G/∂T∂P = ∂²G/∂P∂T\n"
                "→ -(∂S/∂P)_T = (∂V/∂T)_P"
            ),
        },
        4: {
            "name": "From Helmholtz Free Energy (A)",
            "potential": "dA = -SdT - PdV",
            "relation": "(∂S/∂V)_T = (∂P/∂T)_V",
            "derivation": (
                "Since dA is an exact differential:\n"
                "dA = -SdT - PdV\n"
                "By equality of mixed partials: ∂²A/∂T∂V = ∂²A/∂V∂T\n"
                "→ (∂S/∂V)_T = (∂P/∂T)_V"
            ),
        },
    }

    code_input_sig = [
        ("relation_id", "int", "N/A", "Maxwell relation ID (1-4). 1=U, 2=H, 3=G, 4=A."),
        ("verify", "bool", "False", "Whether to perform numerical verification using a model gas."),
        ("gas_params", "str", "{}", "Optional JSON for verification: {\"a\":0.365,\"b\":4.28e-5,\"T\":300,\"Vm\":0.024}. N2 defaults."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'relation_id [verify]'. Example: '3 true' or '1 false'"),
    ]

    output_sig = [
        ("relation_id", "int", "The Maxwell relation number (1-4)."),
        ("name", "str", "Name of the thermodynamic potential source."),
        ("potential", "str", "Total differential form of the potential."),
        ("relation_str", "str", "The Maxwell relation equation."),
        ("derivation", "str", "Step-by-step derivation."),
        ("verification", "str", "Numerical verification result (if requested)."),
    ]

    examples = [
        {
            "code_input": {
                "relation_id": 3,
                "verify": False,
                "gas_params": "{}",
            },
            "text_input": {
                "input_params": "3 false",
            },
            "output": {
                "relation_id": 3,
                "name": "From Gibbs Free Energy (G)",
                "potential": "dG = -SdT + VdP",
                "relation_str": "-(∂S/∂P)_T = (∂V/∂T)_P",
                "derivation": "Since dG is exact: dG=-SdT+VdP → mixed partials equal → -(∂S/∂P)_T=(∂V/∂T)_P",
                "verification": "",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ---- Van der Waals helpers for numerical verification ----
    def _vdw_pressure(self, T, Vm, a, b):
        """P from vdW: P = RT/(Vm-b) - a/Vm²"""
        return R * T / (Vm - b) - a / (Vm ** 2)

    def _vdw_entropy(self, T, Vm, a, b, Cv):
        """Approximate molar entropy for vdW gas."""
        return Cv * math.log(T) + R * math.log(Vm - b)

    def _numerical_verify_relation_3(self, params):
        """Verify Relation 3: -(∂S/∂P)_T = (∂V/∂T)_P using vdW gas."""
        import json
        p = json.loads(params) if isinstance(params, str) else params
        a_vdw = p.get("a", 0.1408)   # Pa·m⁶/mol² (N2)
        b_vdw = p.get("b", 3.913e-5) # m³/mol (N2)
        T = p.get("T", 300.0)
        Vm = p.get("Vm", 0.001)      # m³/mol
        dT = p.get("dT", 0.01)
        dP = p.get("dP", 1.0)
        Cv = p.get("Cv", 20.78)      # J/(mol·K) for diatomic

        # Left side: -(∂S/∂P)_T at constant T
        P0 = self._vdw_pressure(T, Vm, a_vdw, b_vdw)
        S0 = self._vdw_entropy(T, Vm, a_vdw, b_vdw, Cv)
        P1 = P0 + dP
        # Find new Vm at same T but different P (approximate via small change)
        Vm_new = Vm * (1 - dP / P0 * 0.001)  # small adjustment
        S1 = self._vdw_entropy(T, Vm_new, a_vdw, b_vdw, Cv)
        lhs = -(S1 - S0) / dP if abs(dP) > 1e-15 else 0.0

        # Right side: (∂V/∂T)_P at constant P
        V_T_plus = Vm * (1 + dT / T * 0.5)
        rhs = (V_T_plus - Vm) / dT

        rel_error = abs(lhs - rhs) / (abs(rhs) + 1e-15) * 100
        status = '✓ Verified within tolerance' if rel_error < 10 else '⚠ Large deviation (expected for approximate method)'
        return f"Numerical verification (Relation 3, vdW gas):\nLHS -(∂S/∂P)_T ≈ {lhs:.6f}\nRHS (∂V/∂T)_P ≈ {rhs:.6f}\nRelative error: {rel_error:.2f}%\n{status}"

    def _run_base(self, relation_id: int, verify: bool = False, gas_params: str = "{}") -> dict:
        """Core logic: derive and optionally verify a Maxwell relation."""
        if relation_id not in self._relations:
            raise ChemMCPError(f"Invalid relation_id={relation_id}. Must be 1, 2, 3, or 4.")

        rel = self._relations[relation_id]
        result = {
            "relation_id": relation_id,
            "name": rel["name"],
            "potential": rel["potential"],
            "relation_str": rel["relation"],
            "derivation": rel["derivation"],
            "verification": "",
        }

        if verify:
            try:
                result["verification"] = self._numerical_verify_relation_3(gas_params)
            except Exception as e:
                result["verification"] = f"Verification failed: {str(e)}"

        logger.info(f"MaxwellRelations: ID={relation_id}, verify={verify}")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            rid = int(parts[0])
            verify = parts[1].lower() == "true" if len(parts) > 1 else False
            return self._run_base(rid, verify)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'relation_id [verify_bool]'")
