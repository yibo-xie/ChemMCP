import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PhaseRuleAnalyzer(BaseTool):
    """
    吉布斯相律分析工具。
    计算系统的自由度 F = C - P + 2（或 C - P + 1 for condensed systems, or C - P + N with extra constraints）。
    """
    __version__ = "0.1.0"
    name = "PhaseRuleAnalyzer"
    func_name = "analyze_phase_rule"
    description = "Gibbs phase rule analysis: determine degrees of freedom (variance) for a given system."
    implementation_description = "Applies Gibbs phase rule: F = C - P + N (where N=2 for T&P variables, N=1 for condensed). Handles special constraints (azeotrope, reaction equilibrium, etc.)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Phase Rule", "Physical Chemistry", "Equilibrium"]
    required_envs    = []

    code_input_sig = [
        ("n_components", "int", "N/A", "Number of components C (chemically independent constituents)."),
        ("n_phases", "int", "N/A", "Number of phases P present in equilibrium."),
        ("n_variables", "int", "2", "Number of intensive variables N (default 2 for T and P; use 1 for condensed system at constant pressure)."),
        ("n_constraints", "int", "0", "Number of additional independent constraints (e.g., azeotropic condition, reaction equilibrium, stoichiometric relation)."),
        ("system_name", "str", "", "Optional name/description of the system."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: n_components n_phases [n_variables n_constraints system_name]."),
    ]

    output_sig = [
        ("degrees_of_freedom", "int", "Degrees of freedom F (variance) of the system."),
        ("phase_rule_equation", "str", "The phase rule equation used."),
        ("interpretation", "str", "Physical interpretation of what F means for this system."),
        ("system_type", "str", "Classification: invariant, univariant, bivariant, etc."),
        ("examples_scenarios", "list", "Example scenarios consistent with this F value."),
    ]

    examples = [
        {
            "code_input": {
                "n_components": 1,
                "n_phases": 3,
                "n_variables": 2,
                "n_constraints": 0,
                "system_name": "Water triple point",
            },
            "text_input": {
                "input_params": "1 3 2 0 Water_triple_point"
            },
            "output": {
                "degrees_of_freedom": 0,
                "phase_rule_equation": "F = C - P + 2 = 1 - 3 + 2 = 0",
                "interpretation": "Invariant system: both temperature and pressure are fixed at the triple point.",
                "system_type": "Invariant (F=0)",
                "examples_scenarios": ["Triple point of pure substance"],
            }
        },
        {
            "code_input": {
                "n_components": 2,
                "n_phases": 2,
                "n_variables": 2,
                "n_constraints": 1,
                "system_name": "Azeotropic binary mixture",
            },
            "text_input": {
                "input_params": "2 2 2 1 Azeotropic_binary_mixture"
            },
            "output": {
                "degrees_of_freedom": 1,
                "phase_rule_equation": "F = C - P + 2 - R' = 2 - 2 + 2 - 1 = 1",
                "interpretation": "Univariant: fixing one variable (T or composition) determines the rest.",
                "system_type": "Univariant (F=1)",
                "examples_scenarios": ["Boiling point curve of an azeotrope"],
            }
        },
        {
            "code_input": {
                "n_components": 3,
                "n_phases": 2,
                "n_variables": 1,
                "n_constraints": 1,
                "system_name": "Ternary liquid-liquid extraction",
            },
            "text_input": {
                "input_params": "3 2 1 1 Ternary_LLE"
            },
            "output": {
                "degrees_of_freedom": 1,
                "phase_rule_equation": "F = C - P + 1 - R' = 3 - 2 + 1 - 1 = 1",
                "interpretation": "Condensed ternary LLE system with one constraint.",
                "system_type": "Univariant (F=1)",
                "examples_scenarios": ["Binodal curve on ternary diagram"],
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        n_components: int,
        n_phases: int,
        n_variables: int = 2,
        n_constraints: int = 0,
        system_name: str = "",
    ) -> dict:
        if n_components < 1:
            raise ChemMCPError("Number of components must be >= 1.")
        if n_phases < 1:
            raise ChemMCPError("Number of phases must be >= 1.")
        if n_variables < 1:
            raise ChemMCPError("Number of variables must be >= 1.")
        if n_constraints < 0:
            raise ChemMCPError("Constraints cannot be negative.")

        # Gibbs phase rule: F = C - P + N - R'
        F = n_components - n_phases + n_variables - n_constraints

        if F < 0:
            raise ChemMCPError(
                f"Negative degrees of freedom (F={F}): system is over-constrained. "
                f"Check component count, phase count, and constraints."
            )

        # Classification
        if F == 0:
            stype = "Invariant (F=0)"
            interp = (
                f"The system is invariant: all intensive properties are fixed. "
                f"No variables can be changed independently without changing the number of phases."
            )
            scenarios = [
                f"{'(' + system_name + ') ' if system_name else ''}All T, P, and compositions are fixed",
            ]
        elif F == 1:
            stype = "Univariant (F=1)"
            interp = (
                f"The system is univariant: exactly one intensive variable can be changed independently. "
                f"All other variables adjust accordingly along an equilibrium line/curve."
            )
            scenarios = [
                f"Coexistence curve (e.g., melting, vaporization)",
                f"Eutectic or peritectic line",
                f"Azeotropic boiling point vs pressure",
            ]
        elif F == 2:
            stype = "Bivariant (F=2)"
            interp = (
                f"The system is bivariant: two independent variables can be changed freely. "
                f"The system exists as a region (area) on a phase diagram."
            )
            scenarios = [
                f"Single-phase region on binary T-x diagram",
                f"Two-phase area on ternary diagram",
            ]
        elif F == 3:
            stype = "Trivariant (F=3)"
            interp = "Three independent variables; system has large flexibility."
            scenarios = [f"Single-phase ternary mixture at constant P"]
        else:
            stype = f"Multivariant (F={F})"
            interp = f"The system has {F} independent degrees of freedom."
            scenarios = [f"Complex multi-component system"]

        eq_str = f"F = C - P + {n_variables}"
        if n_constraints > 0:
            eq_str += f" - {n_constraints}"
        eq_str += f" = {n_components} - {n_phases} + {n_variables}"
        if n_constraints > 0:
            eq_str += f" - {n_constraints}"
        eq_str += f" = {F}"

        return {
            "degrees_of_freedom": F,
            "phase_rule_equation": eq_str,
            "interpretation": interp,
            "system_type": stype,
            "examples_scenarios": scenarios,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            c = int(parts[0])
            p = int(parts[1])
            n_var = int(parts[2]) if len(parts) > 2 else 2
            n_con = int(parts[3]) if len(parts) > 3 else 0
            name = parts[4] if len(parts) > 4 else ""
            return self._run_base(c, p, n_var, n_con, name)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
