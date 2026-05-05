import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IonTransportNumber(BaseTool):
    """Ion transport number (transference number) calculation tool."""
    __version__ = "0.1.0"
    name = "IonTransportNumber"
    func_name = "calculate_transport_number"
    description = "Calculate ion transport numbers from mobility, concentration, conductivity, or Hittorf method data."
    implementation_description = "t_i = c_i*u_i / sum(c_j*u_j); Kohlrausch law; Hittorf method analysis."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Ion Transport", "Conductivity"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "'mobility', 'conductivity', or 'hittorf'."),
        ("cation_mobility", "float", "None", "Cation mobility u+ (m2/(V*s)). For 'mobility'."),
        ("anion_mobility", "float", "None", "Anion mobility u-. For 'mobility'."),
        ("cation_concentration", "float", "None", "Cation concentration (mol/m3). For 'mobility'."),
        ("anion_concentration", "float", "None", "Anion concentration (mol/m3). For 'mobility'."),
        ("cation_molar_conductivity", "float", "None", "Cation molar conductivity lambda+ (S*m2/mol). For 'conductivity'."),
        ("anion_molar_conductivity", "float", "None", "Anion molar conductivity lambda-. For 'conductivity'."),
        ("initial_cation_conc", "float", "None", "Initial cation conc near cathode (mol/L). For 'hittorf'."),
        ("final_cation_conc", "float", "None", "Final cation conc near cathode (mol/L). For 'hittorf'."),
        ("charge_passed_coulomb", "float", "None", "Total charge passed Q (C). For 'hittorf'."),
        ("volume_liter", "float", "1.0", "Electrolyte volume (L). For 'hittorf'. Default 1.0."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Mode-specific params space-separated."),
    ]

    output_sig = [
        ("transport_number_cation", "float", "Cation transport number t+."),
        ("transport_number_anion", "float", "Anion transport number t-."),
        ("method_used", "str", "Description of calculation method."),
        ("details", "dict", "Additional details."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "mobility",
                "cation_mobility": 7.62e-8,
                "anion_mobility": 7.91e-8,
                "cation_concentration": 100.0,
                "anion_concentration": 100.0,
                "cation_molar_conductivity": None,
                "anion_molar_conductivity": None,
                "initial_cation_conc": None,
                "final_cation_conc": None,
                "charge_passed_coulomb": None,
                "volume_liter": 1.0,
            },
            "text_input": {
                "input_params": "mobility 7.62e-8 7.91e-8 100 100",
            },
            "output": {
                "transport_number_cation": 0.4906,
                "transport_number_anion": 0.5094,
                "method_used": "Mobility-based",
                "details": {},
            }
        },
        {
            "code_input": {
                "mode": "conductivity",
                "cation_mobility": None,
                "anion_mobility": None,
                "cation_concentration": None,
                "anion_concentration": None,
                "cation_molar_conductivity": 0.00735,
                "anion_molar_conductivity": 0.00763,
                "initial_cation_conc": None,
                "final_cation_conc": None,
                "charge_passed_coulomb": None,
                "volume_liter": 1.0,
            },
            "text_input": {
                "input_params": "conductivity 0.00735 0.00763",
            },
            "output": {
                "transport_number_cation": 0.4906,
                "transport_number_anion": 0.5094,
                "method_used": "Molar conductivity-based (Kohlrausch)",
                "details": {},
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.F = 96485.33212

    def _run_base(self, mode, **kwargs) -> dict:
        mode = mode.lower().strip()

        if mode == "mobility":
            u_plus = kwargs.get("cation_mobility"); u_minus = kwargs.get("anion_mobility")
            c_plus = kwargs.get("cation_concentration"); c_minus = kwargs.get("anion_concentration")
            if any(v is None for v in [u_plus, u_minus, c_plus, c_minus]):
                raise ChemMCPError("'mobility' needs cation/anion mobility and concentration.")
            cond_plus = c_plus * u_plus; cond_minus = c_minus * u_minus; total = cond_plus + cond_minus
            if total <= 0: raise ChemMCPError("Total conductivity must be positive.")
            t_plus = cond_plus / total; t_minus = cond_minus / total
            return {"transport_number_cation": round(t_plus, 6), "transport_number_anion": round(t_minus, 6),
                    "method_used": "Mobility-based: t+ = c+*u+/(c+*u+ + c-*u-)", "details": {}}

        elif mode == "conductivity":
            lp = kwargs.get("cation_molar_conductivity"); lm = kwargs.get("anion_molar_conductivity")
            if any(v is None for v in [lp, lm]):
                raise ChemMCPError("'conductivity' needs both molar conductivities.")
            total = lp + lm
            if total <= 0: raise ChemMCPError("Total must be positive.")
            return {"transport_number_cation": round(lp/total, 6), "transport_number_anion": round(lm/total, 6),
                    "method_used": "Kohlrausch: t+ = lambda+/(lambda+ + lambda-)", "details": {}}

        elif mode == "hittorf":
            ci = kwargs.get("initial_cation_conc"); cf = kwargs.get("final_cation_conc")
            Q = kwargs.get("charge_passed_coulomb"); V = kwargs.get("volume_liter", 1.0)
            if any(v is None for v in [ci, cf, Q]): raise ChemMCPError("'hittorf' needs initial/final conc and Q.")
            delta_n = abs(ci - cf) * V; n_tot = Q / self.F
            if n_tot <= 0: raise ChemMCPError("Q must yield positive mole equivalent.")
            t_plus = delta_n / n_tot
            return {"transport_number_cation": round(t_plus, 6), "transport_number_anion": round(1-t_plus, 6),
                    "method_used": f"Hittorf: t+ = |dc|*V*F/Q", "details": {"delta_moles": delta_n, "total_moles": n_tot}}
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use: mobility, conductivity, hittorf.")

    def _run_text(self, input_params: str) -> dict:
        try:
            p = input_params.split(); mode = p[0]
            if mode == "mobility":
                return self._run_base(mode, cation_mobility=float(p[1]), anion_mobility=float(p[2]),
                    cation_concentration=float(p[3]), anion_concentration=float(p[4]))
            elif mode == "conductivity":
                return self._run_base(mode, cation_molar_conductivity=float(p[1]), anion_molar_conductivity=float(p[2]))
            elif mode == "hittorf":
                return self._run_base(mode, initial_cation_conc=float(p[1]), final_cation_conc=float(p[2]),
                    charge_passed_coulomb=float(p[3]))
            else: raise ValueError(mode)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {str(e)}")
