import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FreundlichIsotherm(BaseTool):
    """Freundlich isotherm adsorption model tool."""
    __version__ = "0.1.0"
    name = "FreundlichIsotherm"
    func_name = "freundlich_isotherm"
    description = "Freundlich isotherm: fit K_f and n from data, calculate uptake, assess adsorption intensity."
    implementation_description = "q = K_f * P^(1/n). Log-log linearization with least-squares fitting. Intensity interpretation via 1/n value."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Adsorption", "Freundlich", "Surface Chemistry", "Isotherm"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "'calculate': q from P; 'fit': fit K_f,n from data."),
        ("pressure", "float", "None", "Pressure P (for 'calculate')."),
        ("kf", "float", "None", "Freundlich constant K_f (for 'calculate')."),
        ("n_exponent", "float", "None", "Exponent n, typically >1 (for 'calculate')."),
        ("experimental_pressures", "list", "None", "Pressures for 'fit'."),
        ("experimental_uptakes", "list", "None", "Uptakes for 'fit'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Calculate: 'calculate P K_f n'. Fit: 'fit P1,P2,... q1,q2,...'."),
    ]

    output_sig = [
        ("uptake_q", "float", "Calculated q. None for fit."),
        ("fitted_Kf", "float", "Fitted K_f. Only for fit."),
        ("fitted_n", "float", "Fitted n. Only for fit."),
        ("one_over_n", "float", "Heterogeneity factor 1/n."),
        ("r_squared", "float", "R2 of log-log fit. Only for fit."),
        ("intensity_interpretation", "str", "Interpretation of 1/n."),
        ("analysis_summary", "str", "Summary."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "calculate",
                "pressure": 0.3,
                "kf": 5.0,
                "n_exponent": 2.0,
                "experimental_pressures": None,
                "experimental_uptakes": None,
            },
            "text_input": {
                "input_params": "calculate 0.3 5.0 2.0",
            },
            "output": {
                "uptake_q": 2.7385,
                "fitted_Kf": None,
                "fitted_n": None,
                "one_over_n": 0.5,
                "r_squared": None,
                "intensity_interpretation": "favorable adsorption",
                "analysis_summary": "At P=0.3, Kf=5.0, n=2.0: q=2.74",
            }
        },
        {
            "code_input": {
                "mode": "fit",
                "pressure": None,
                "kf": None,
                "n_exponent": None,
                "experimental_pressures": [0.05, 0.1, 0.2, 0.5, 1.0],
                "experimental_uptakes": [0.56, 1.0, 1.79, 3.54, 5.62],
            },
            "text_input": {
                "input_params": "fit 0.05,0.1,0.2,0.5,1.0 0.56,1.0,1.79,3.54,5.62",
            },
            "output": {
                "uptake_q": None,
                "fitted_Kf": 5.62,
                "fitted_n": 2.0,
                "one_over_n": 0.5,
                "r_squared": 0.999,
                "intensity_interpretation": "favorable adsorption over entire range",
                "analysis_summary": "Freundlich fit: Kf=5.62, n=2.00, R2=0.999",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self): pass

    def _loglog_reg(self, xs, ys):
        lx=[math.log(x) for x in xs if x>0]; ly=[math.log(y) for y in ys if y>0]
        if len(lx)<2: raise ChemMCPError("Need >=2 positive points.")
        n=len(lx); sx=sum(lx); sy=sum(ly); sxy=sum(x*y for x,y in zip(lx,ly))
        sx2=sum(x*x for x in lx); D=n*sx2-sx*sx
        if abs(D)<1e-30: raise ChemMCPError("Cannot fit.")
        slope=(n*sxy-sx*sy)/D; intercept=(sy-slope*sx)/n
        Kf=math.exp(intercept); nv=1/slope if abs(slope)>1e-30 else float("inf")
        ym=sy/n; sst=sum((y-ym)**2 for y in ly); ssr=sum((y-(intercept+slope*x))**2 for x,y in zip(lx,ly))
        return Kf,nv,slope,(1-ssr/sst if sst>0 else 1.0)

    def _interpret(self, invn):
        if invn<0.1: return f"1/n={invn:.2f}: nearly irreversible"
        elif invn<0.25: return f"1/n={invn:.2f}: strongly favorable"
        elif invn<0.5: return f"1/n={invn:.2f}: moderately favorable"
        elif invn<1.0: return f"1/n={invn:.2f}: favorable"
        elif abs(invn-1)<0.05: return "1/n~1: linear (C-type)"
        else: return f"1/n={invn:.2f}: unfavorable"

    def _run_base(self, mode, pressure=None, kf=None, n_exponent=None,
                   experimental_pressures=None, experimental_uptakes=None) -> dict:
        mode=mode.lower().strip()
        if mode=="calculate":
            if any(v is None for v in [pressure,kf,n_exponent]): raise ChemMCPError("'calculate' needs P, kf, n.")
            if pressure<=0 or n_exponent<=0: raise ChemMCPError("P>0, n>0 required.")
            q=kf*(pressure**(1.0/n_exponent))
            invn=1.0/n_exponent
            return {"uptake_q":round(q,6),"fitted_Kf":None,"fitted_n":None,
                    "one_over_n":round(invn,6),"r_squared":None,
                    "intensity_interpretation":self._interpret(invn),
                    "analysis_summary":f"At P={pressure}, Kf={kf}, n={n_exponent}: q={round(q,4)}"}
        elif mode=="fit":
            if experimental_pressures is None or experimental_uptakes is None:
                raise ChemMCPError("'fit' needs data lists.")
            if len(experimental_pressures)!=len(experimental_uptakes):
                raise ChemMCPError("Lists equal length required.")
            Kf,nv,invn,r2=self._loglog_reg(experimental_pressures,experimental_uptakes)
            return {"uptake_q":None,"fitted_Kf":round(Kf,6),"fitted_n":round(nv,6),
                    "one_over_n":round(invn,6),"r_squared":round(r2,6),
                    "intensity_interpretation":self._interpret(invn),
                    "analysis_summary":f"Freundlich fit: Kf={round(Kf,4)}, n={round(nv,3)} (1/n={round(invn,3)}), R2={round(r2,4)}"}
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use calculate/fit.")

    def _run_text(self,s:str)->dict:
        try:
            p=s.split(); m=p[0]
            if m=="calculate": return self._run_base(m,float(p[1]),float(p[2]),float(p[3]))
            elif m=="fit":
                pp=[float(x) for x in p[1].split(",")]; qu=[float(x) for x in p[2].split(",")]
                return self._run_base(m,experimental_pressures=pp,experimental_uptakes=qu)
            else: raise ValueError(m)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")
