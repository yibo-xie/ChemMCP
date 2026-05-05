import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ButlerVolmerKinetics(BaseTool):
    """Butler-Volmer equation electrode kinetics tool."""
    __version__ = "0.1.0"
    name = "ButlerVolmerKinetics"
    func_name = "butler_volmer"
    description = "Compute electrode kinetics using the full Butler-Volmer equation with forward/inverse solving and Tafel analysis."
    implementation_description = "Full BV: j = j0*[exp(alpha_a*n*F*eta/RT) - exp(-alpha_c*n*F*eta/RT)]. Newton-Raphson inverse solver included."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Kinetics", "Butler-Volmer", "Tafel"]
    required_envs = []

    code_input_sig = [
        ("calculation_mode", "str", "N/A", "'forward': j from eta; 'inverse': eta from j; 'full': full analysis."),
        ("overpotential_v", "float", "None", "Overpotential eta (V). For 'forward'/'full' modes."),
        ("current_density_ma_cm2", "float", "None", "Current density j (mA/cm2). For 'inverse' mode."),
        ("exchange_current_density_ma_cm2", "float", "N/A", "Exchange current density j0 (mA/cm2)."),
        ("anodic_coefficient", "float", "0.5", "Anodic transfer coeff alpha_a. Default 0.5."),
        ("cathodic_coefficient", "float", "0.5", "Cathodic transfer coeff alpha_c. Default 0.5."),
        ("electrons_transferred", "int", "1", "Electrons n. Default 1."),
        ("temperature_k", "float", "298.15", "Temperature (K). Default 298.15."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: mode j0 [eta_or_j] [alpha_a] [alpha_c] [n] [T]."),
    ]

    output_sig = [
        ("current_density_ma_cm2", "float", "Current density j (mA/cm2)."),
        ("overpotential_v", "float", "Overpotential eta (V)."),
        ("anodic_current_density", "float", "Anodic component j_a (mA/cm2)."),
        ("cathodic_current_density", "float", "Cathodic component j_c (mA/cm2)."),
        ("tafel_slope_anodic_mv", "float", "Anodic Tafel slope b_a (mV/dec)."),
        ("tafel_slope_cathodic_mv", "float", "Cathodic Tafel slope b_c (mV/dec)."),
        ("regime", "str", "Kinetic regime classification."),
        ("analysis_summary", "str", "Human-readable summary."),
    ]

    examples = [
        {
            "code_input": {
                "calculation_mode": "forward",
                "overpotential_v": 0.05,
                "current_density_ma_cm2": None,
                "exchange_current_density_ma_cm2": 0.001,
                "anodic_coefficient": 0.5,
                "cathodic_coefficient": 0.5,
                "electrons_transferred": 1,
                "temperature_k": 298.15,
            },
            "text_input": {
                "input_params": "forward 0.001 0.05",
            },
            "output": {
                "current_density_ma_cm2": 1.906,
                "overpotential_v": 0.05,
                "anodic_current_density": 1.953,
                "cathodic_current_density": 0.047,
                "tafel_slope_anodic_mv": 118.32,
                "tafel_slope_cathodic_mv": 118.32,
                "regime": "tafel_anodic",
                "analysis_summary": "At eta=+50mV: j=1.91 mA/cm2, anodic branch dominates",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.F = 96485.33212; self.R = 8.314462618

    def _compute_forward(self, eta, j0, aa, ac, n, T):
        frt = self.F / (self.R * T)
        ja = j0 * math.exp(aa * n * frt * eta)
        jc = j0 * math.exp(-ac * n * frt * eta)
        return ja - jc, ja, jc

    def _compute_inverse(self, jt, j0, aa, ac, n, T, max_iter=100, tol=1e-10):
        frt = self.F / (self.R * T)
        if abs(jt) > 0:
            if jt > 0: eta = (self.R*T/(aa*n*self.F))*math.log(jt/j0+1)
            else: eta = (self.R*T/(ac*n*self.F))*math.log(1-jt/j0)
        else: eta = 0.0
        for _ in range(max_iter):
            j, ja, jc = self._compute_forward(eta, j0, aa, ac, n, T)
            res = j - jt
            djd = ja*aa*n*frt + jc*ac*n*frt
            if abs(djd) < 1e-30: break
            enew = eta - res/djd
            if abs(enew-eta) < tol: eta=enew; break
            eta = enew
        return eta

    def _determine_regime(self, eta, ja, jc):
        if abs(eta) < 0.01: return "linear (ohmic)"
        elif eta > 0.05 and ja > 10*jc: return "tafel_anodic"
        elif eta < -0.05 and jc > 10*ja: return "tafel_cathodic"
        else: return "mixed"

    def _run_base(self, calculation_mode, exchange_current_density_ma_cm2,
                   overpotential_v=None, current_density_ma_cm2=None,
                   anodic_coefficient=0.5, cathodic_coefficient=0.5,
                   electrons_transferred=1, temperature_k=298.15) -> dict:
        mode = calculation_mode.lower().strip(); j0 = exchange_current_density_ma_cm2
        n = electrons_transferred; T = temperature_k; aa = anodic_coefficient; ac = cathodic_coefficient
        if j0 <= 0: raise ChemMCPError("j0 must be positive.")
        if T <= 0: raise ChemMCPError("T must be positive.")

        ba = (math.log(10)*self.R*T)/(aa*n*self.F)*1000
        bc = (math.log(10)*self.R*T)/(ac*n*self.F)*1000

        if mode in ("forward", "full"):
            if overpotential_v is None: raise ChemMCPError("'forward' needs overpotential_v.")
            eta = overpotential_v; j, ja, jc = self._compute_forward(eta, j0, aa, ac, n, T)
        elif mode == "inverse":
            if current_density_ma_cm2 is None: raise ChemMCPError("'inverse' needs current_density_ma_cm2.")
            eta = self._compute_inverse(current_density_ma_cm2, j0, aa, ac, n, T)
            j, ja, jc = self._compute_forward(eta, j0, aa, ac, n, T)
        else: raise ChemMCPError(f"Unknown mode '{mode}'.")

        regime = self._determine_regime(eta, ja, jc)
        parts = []
        if mode == "inverse": parts.append(f"For j={current_density_ma_cm2}: eta={round(abs(eta)*1000,1)}mV ({'anodic' if eta>0 else 'cathodic'})")
        parts.append(f"j={round(j,4)} (ja={round(ja,4)}, jc={round(jc,4)})")
        parts.append(f"Tafel: ba={round(ba,1)}, bc={round(bc,1)} mV/dec")
        parts.append(f"Regime: {regime}")

        logger.info(f"Butler-Volmer [{mode}]: eta={eta}, j={j}, regime={regime}")
        return {"current_density_ma_cm2": round(j,6), "overpotential_v": round(eta,6),
                "anodic_current_density": round(ja,6), "cathodic_current_density": round(jc,6),
                "tafel_slope_anodic_mv": round(ba,2), "tafel_slope_cathodic_mv": round(bc,2),
                "regime": regime, "analysis_summary": "; ".join(parts)}

    def _run_text(self, input_params: str) -> dict:
        try:
            p = input_params.split(); kw = {"calculation_mode": p[0], "exchange_current_density_ma_cm2": float(p[1])}
            if len(p)>2:
                v=float(p[2]); kw["overpotential_v" if p[0] in ("forward","full") else "current_density_ma_cm2"]=v
            if len(p)>3: kw["anodic_coefficient"]=float(p[3])
            if len(p)>4: kw["cathodic_coefficient"]=float(p[4])
            if len(p)>5: kw["electrons_transferred"]=int(p[5])
            if len(p)>6: kw["temperature_k"]=float(p[6])
            return self._run_base(**kw)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {str(e)}")
