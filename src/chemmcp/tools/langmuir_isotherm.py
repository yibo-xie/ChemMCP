import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class LangmuirIsotherm(BaseTool):
    """Langmuir isotherm adsorption model tool."""
    __version__ = "0.1.0"
    name = "LangmuirIsotherm"
    func_name = "langmuir_isotherm"
    description = "Langmuir isotherm: calculate coverage, fit K/qmax from data, Lineweaver-Burk linearization."
    implementation_description = "theta = K*P/(1+KP); q = qmax*K*P/(1+KP). Linear regression on 1/q vs 1/P for fitting."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Adsorption", "Surface Chemistry", "Langmuir", "Isotherm"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "'calculate': q from P; 'fit': fit K,qmax from data; 'linearize': LB plot data."),
        ("pressure_or_data", "float", "None", "Pressure P (for 'calculate')."),
        ("equilibrium_constant_k", "float", "None", "Langmuir constant K (for 'calculate')."),
        ("max_adsorption_qmax", "float", "None", "Max capacity qmax (for 'calculate')."),
        ("experimental_pressures", "list", "None", "Pressures [P1,P2,...] (for 'fit'/'linearize')."),
        ("experimental_uptakes", "list", "None", "Uptakes [q1,q2,...] (for 'fit'/'linearize')."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Calculate: 'calculate P K qmax'. Fit: 'fit P1,P2,... q1,q2,...'."),
    ]

    output_sig = [
        ("uptake_q", "float", "Adsorption amount q."),
        ("surface_coverage_theta", "float", "Fractional coverage theta."),
        ("fitted_K", "float", "Fitted K. Only for 'fit'."),
        ("fitted_qmax", "float", "Fitted qmax. Only for 'fit'."),
        ("r_squared", "float", "R2 of fit. Only for 'fit'/'linearize'."),
        ("linearized_data", "list", "LB coordinates [(1/P,1/q)]. Only for 'linearize'."),
        ("analysis_summary", "str", "Summary."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "calculate",
                "pressure_or_data": 0.5,
                "equilibrium_constant_k": 2.0,
                "max_adsorption_qmax": 10.0,
                "experimental_pressures": None,
                "experimental_uptakes": None,
            },
            "text_input": {
                "input_params": "calculate 0.5 2.0 10.0",
            },
            "output": {
                "uptake_q": 9.0909,
                "surface_coverage_theta": 0.9091,
                "fitted_K": None,
                "fitted_qmax": None,
                "r_squared": None,
                "linearized_data": None,
                "analysis_summary": "At P=0.5 atm, K=2.0: theta=0.909, q=9.09",
            }
        },
        {
            "code_input": {
                "mode": "fit",
                "pressure_or_data": None,
                "equilibrium_constant_k": None,
                "max_adsorption_qmax": None,
                "experimental_pressures": [0.05, 0.1, 0.2, 0.5, 1.0],
                "experimental_uptakes": [1.67, 3.0, 5.0, 7.14, 8.33],
            },
            "text_input": {
                "input_params": "fit 0.05,0.1,0.2,0.5,1.0 1.67,3.0,5.0,7.14,8.33",
            },
            "output": {
                "uptake_q": None,
                "surface_coverage_theta": None,
                "fitted_K": 2.0,
                "fitted_qmax": 10.01,
                "r_squared": 0.999,
                "linearized_data": None,
                "analysis_summary": "Langmuir fit: K=2.00, qmax=10.01, R2=0.999",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self): pass

    def _linreg(self, xs, ys):
        n=len(xs)
        sx=sum(xs); sy=sum(ys); sxy=sum(x*y for x,y in zip(xs,ys)); sx2=sum(x*x for x in xs)
        D=n*sx2-sx*sx
        if abs(D)<1e-30: raise ChemMCPError("Cannot fit.")
        b=(n*sxy-sx*sy)/D; a=(sy-b*sx)/n
        ym=sy/n
        sst=sum((y-ym)**2 for y in ys); ssr=sum((y-(a+b*x))**2 for x,y in zip(xs,ys))
        r2=1-ssr/sst if sst>0 else 1.0
        return a,b,r2

    def _run_base(self, mode, pressure_or_data=None, equilibrium_constant_k=None,
                   max_adsorption_qmax=None, experimental_pressures=None,
                   experimental_uptakes=None) -> dict:
        mode=mode.lower().strip()
        if mode=="calculate":
            if any(v is None for v in [pressure_or_data,equilibrium_constant_k,max_adsorption_qmax]):
                raise ChemMCPError("'calculate' needs P, K, qmax.")
            P=pressure_or_data; K=equilibrium_constant_k; qm=max_adsorption_qmax
            if K<0 or qm<=0: raise ChemMCPError("K>=0, qmax>0 required.")
            theta=(K*P)/(1+K*P); q=qm*theta
            return {"uptake_q":round(q,6), "surface_coverage_theta":round(theta,6),
                    "fitted_K":None,"fitted_qmax":None,"r_squared":None,"linearized_data":None,
                    "analysis_summary":f"At P={P}, K={K}: theta={round(theta,4)}, q={round(q,4)}"}
        elif mode in ("fit","linearize"):
            if experimental_pressures is None or experimental_uptakes is None:
                raise ChemMCPError(f"'{mode}' needs data lists.")
            if len(experimental_pressures)!=len(experimental_uptakes):
                raise ChemMCPError("Lists must have equal length.")
            inv_p=[1/p if p!=0 else float("inf") for p in experimental_pressures]
            inv_q=[1/q if q!=0 else float("inf") for q in experimental_uptakes]
            valid=[(ip,iq) for ip,iq in zip(inv_p,inv_q) if ip!=float("inf") and iq!=float("inf")]
            if len(valid)<2: raise ChemMCPError("Need >=2 valid points.")
            ipl,iql=zip(*valid); a,b,r2=self._linreg(list(ipl),list(iql))
            qm_fit=1/a if abs(a)>1e-30 else float("inf")
            K_fit=a/b if abs(b)>1e-30 else float("inf")
            ldata=[(ip,iq) for ip,iq in zip(inv_p,inv_q)] if mode=="linearize" else None
            return {"uptake_q":None,"surface_coverage_theta":None,
                    "fitted_K":round(K_fit,6),"fitted_qmax":round(qm_fit,6),
                    "r_squared":round(r2,6),"linearized_data":ldata,
                    "analysis_summary":f"Langmuir {mode}: K={round(K_fit,4)}, qmax={round(qm_fit,4)}, R2={round(r2,4)}"}
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use calculate/fit/linearize.")

    def _run_text(self, s:str)->dict:
        try:
            p=s.split(); m=p[0]
            if m=="calculate": return self._run_base(m,float(p[1]),float(p[2]),float(p[3]))
            elif m in ("fit","linearize"):
                pp=[float(x) for x in p[1].split(",")]; qu=[float(x) for x in p[2].split(",")]
                return self._run_base(m,experimental_pressures=pp,experimental_uptakes=qu)
            else: raise ValueError(m)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")
