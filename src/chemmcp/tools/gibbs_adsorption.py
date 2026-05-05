import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class GibbsAdsorption(BaseTool):
    """Gibbs adsorption equation - surface excess concentration tool."""
    __version__ = "0.1.0"
    name = "GibbsAdsorption"
    func_name = "gibbs_adsorption"
    description = "Calculate surface excess concentration Gamma using Gibbs adsorption equation from surface tension vs concentration data."
    implementation_description = "Gamma = -(c/RT)*(dgamma/dc). Numerical differentiation of gamma-c data. Molecular area from A=1/(Gamma*N_A)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Surface Chemistry", "Gibbs Adsorption", "Surface Excess"]
    required_envs = []

    R_GAS = 8.314462618
    AVOGADRO = 6.02214076e23

    code_input_sig = [
        ("mode", "str", "N/A", "'calculate': Gamma from dg/dc; 'from_data': from (c,gamma) data; 'molecular_area': area from Gamma."),
        ("concentration", "float", "None", "Solute concentration c (mol/L). For 'calculate'."),
        ("dgammadc", "float", "None", "Derivative dgamma/dc. For 'calculate'."),
        ("temperature_k", "float", "298.15", "Temperature (K). Default 298.15."),
        ("concentrations_data", "list", "None", "Concentrations [c1,c2,...] (mol/L). For 'from_data'."),
        ("surface_tensions_data", "list", "None", "Surface tensions [g1,g2,...] (N/m). For 'from_data'."),
        ("target_concentration", "float", "None", "Target c for Gamma in 'from_data'. If None, returns all."),
        ("surface_excess", "float", "None", "Surface excess Gamma (mol/m2). For 'molecular_area'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Calculate: 'calculate c dgammadc [T]'. Data: 'from_data c_list g_list [target_c]'. Area: 'molecular_area Gamma'."),
    ]

    output_sig = [
        ("surface_excess_mol_m2", "float", "Surface excess Gamma (mol/m2)."),
        ("surface_excess_mol_cm2", "float", "Gamma (mol/cm2)."),
        ("molecular_area_nm2", "float", "Area per molecule (nm2)."),
        ("dgammadc_value", "float", "dgamma/dc used/computed."),
        ("details", "dict", "Additional details."),
        ("analysis_summary", "str", "Summary."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "calculate",
                "concentration": 0.1,
                "dgammadc": -0.08,
                "temperature_k": 298.15,
                "concentrations_data": None,
                "surface_tensions_data": None,
                "target_concentration": None,
                "surface_excess": None,
            },
            "text_input": {
                "input_params": "calculate 0.1 -0.08",
            },
            "output": {
                "surface_excess_mol_m2": 3.23e-06,
                "surface_excess_mol_cm2": 3.23e-10,
                "molecular_area_nm2": 0.514,
                "dgammadc_value": -0.08,
                "details": {},
                "analysis_summary": "At c=0.1 mol/L, dgdc=-0.08: Gamma=3.23e-6 mol/m2, A=0.51 nm2",
            }
        },
        {
            "code_input": {
                "mode": "from_data",
                "concentration": None,
                "dgammadc": None,
                "temperature_k": 298.15,
                "concentrations_data": [0.001, 0.005, 0.01, 0.05, 0.1, 0.5],
                "surface_tensions_data": [0.0728, 0.0715, 0.0698, 0.0630, 0.0560, 0.0450],
                "target_concentration": 0.05,
                "surface_excess": None,
            },
            "text_input": {
                "input_params": "from_data 0.001,0.005,0.01,0.05,0.1,0.5 0.0728,0.0715,0.0698,0.063,0.056,0.045 0.05",
            },
            "output": {
                "surface_excess_mol_m2": 2.89e-06,
                "surface_excess_mol_cm2": 2.89e-10,
                "molecular_area_nm2": 0.575,
                "dgammadc_value": -0.142,
                "details": {},
                "analysis_summary": "At c=0.050 mol/L: Gamma=2.89e-6 mol/m2, A=0.58 nm2",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self): pass

    def _mol_area(self, G):
        if G is None or G<=0: return None
        return 1/(G*self.AVOGADRO)*1e18  # m2 -> nm2

    def _num_deriv(self, xs, ys, xt):
        n=len(xs)
        if xt is not None:
            idx=min(range(n), key=lambda i: abs(xs[i]-xt))
            if idx==0: return (ys[1]-ys[0])/(xs[1]-xs[0]) if xs[1]!=xs[0] else None
            elif idx==n-1: return (ys[-1]-ys[-2])/(xs[-1]-xs[-2]) if xs[-1]!=xs[-2] else None
            else: return (ys[idx+1]-ys[idx-1])/(xs[idx+1]-xs[idx-1]) if xs[idx+1]!=xs[idx-1] else None
        else:
            derivs=[]
            for i in range(1,n-1):
                dx=xs[i+1]-xs[i-1]
                if dx!=0: derivs.append((xs[i],(ys[i+1]-ys[i-1])/dx))
            return derivs

    def _run_base(self, mode, concentration=None, dgammadc=None, temperature_k=298.15,
                   concentrations_data=None, surface_tensions_data=None,
                   target_concentration=None, surface_excess=None) -> dict:
        mode=mode.lower().strip(); T=temperature_k; RT=self.R_GAS*T

        if mode=="calculate":
            if concentration is None or dgammadc is None:
                raise ChemMCPError("'calculate' needs c and dgammadc.")
            if concentration<0: raise ChemMCPError("c>=0 required.")
            Gamma=-(concentration/RT)*dgammadc
            area=self._mol_area(Gamma)
            return {"surface_excess_mol_m2":round(Gamma,10),"surface_excess_mol_cm2":round(Gamma*1e-4,12),
                    "molecular_area_nm2":round(area,4) if area else None,
                    "dgammadc_value":dgammadc,"details":{},
                    "analysis_summary":f"c={concentration}, dgdc={dgammadc}, T={T}K: Gamma={Gamma:.3e} mol/m2, A={area:.2f} nm2" if area else f"Gamma={Gamma:.3e}"}

        elif mode=="from_data":
            if concentrations_data is None or surface_tensions_data is None:
                raise ChemMCPError("'from_data' needs data lists.")
            if len(concentrations_data)!=len(surface_tensions_data):
                raise ChemMCPError("Lists equal length required.")
            if len(concentrations_data)<3: raise ChemMCPError("Need >=3 points.")
            cl=concentrations_data; gl=surface_tensions_data
            if target_concentration is not None:
                d=self._num_deriv(cl,gl,target_concentration)
                if d is None: raise ChemMCPError(f"Cannot compute derivative at c={target_concentration}.")
                Gamma=-(target_concentration/RT)*d; area=self._mol_area(Gamma)
                return {"surface_excess_mol_m2":round(Gamma,10),"surface_excess_mol_cm2":round(Gamma*1e-4,12),
                        "molecular_area_nm2":round(area,4) if area else None,"dgammadc_value":round(d,6),
                        "details":{},"analysis_summary":f"At c={target_concentration}: dgdc={d:.4f}, Gamma={Gamma:.3e}, A={area:.2f}" if area else f"Gamma={Gamma:.3e}"}
            else:
                results=[]
                for i in range(1,len(cl)-1):
                    dx=cl[i+1]-cl[i-1]
                    if dx==0: continue
                    d=(gl[i+1]-gl[i-1])/dx; Gi=-(cl[i]/RT)*d; ai=self._mol_area(Gi)
                    results.append({"concentration":cl[i],"surface_excess_mol_m2":round(Gi,10),
                        "molecular_area_nm2":round(ai,4) if ai else None,"dgammadc":round(d,6)})
                return {"surface_excess_mol_m2":None,"surface_excess_mol_cm2":None,
                        "molecular_area_nm2":None,"dgammadc_value":None,
                        "details":{"results":results},"analysis_summary":f"{len(results)} interior points computed"}

        elif mode=="molecular_area":
            if surface_excess is None: raise ChemMCPError("'molecular_area' needs Gamma.")
            area=self._mol_area(surface_excess)
            return {"surface_excess_mol_m2":surface_excess,"surface_excess_mol_cm2":round(surface_excess*1e-4,12),
                    "molecular_area_nm2":round(area,4) if area else None,"dgammadc_value":None,
                    "details":{},"analysis_summary":f"Gamma={surface_excess:.3e} -> A={area:.2f} nm2" if area else f"Invalid Gamma"}
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use calculate/from_data/molecular_area.")

    def _run_text(self,s:str)->dict:
        try:
            p=s.split(); m=p[0]
            if m=="calculate":
                kw={"mode":m,"concentration":float(p[1]),"dgammadc":float(p[2])}
                if len(p)>3: kw["temperature_k"]=float(p[3])
                return self._run_base(**kw)
            elif m=="from_data":
                cc=[float(x) for x in p[1].split(",")]; gg=[float(x) for x in p[2].split(",")]
                kw={"mode":m,"concentrations_data":cc,"surface_tensions_data":gg}
                if len(p)>3: kw["target_concentration"]=float(p[3])
                return self._run_base(**kw)
            elif m=="molecular_area": return self._run_base(m,surface_excess=float(p[1]))
            else: raise ValueError(m)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")
