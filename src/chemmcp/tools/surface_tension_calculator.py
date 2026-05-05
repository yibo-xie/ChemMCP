import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SurfaceTensionCalculator(BaseTool):
    """Surface tension and surface energy calculation tool."""
    __version__ = "0.1.0"
    name = "SurfaceTensionCalculator"
    func_name = "surface_tension"
    description = "Calculate surface tension, Laplace pressure, capillary rise, and temperature-dependent surface properties."
    implementation_description = "Eotvos rule, Young-Laplace equation (dP=2g/R), Jurin's law for capillary rise: h=2g*cos(q)/(rgr)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Surface Chemistry", "Surface Tension", "Capillarity"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "'temperature', 'laplace', 'capillary', or 'surface_energy'."),
        ("temperature_k", "float", "298.15", "Temperature (K). Default 298.15."),
        ("surface_tension_n_m", "float", "None", "Surface tension gamma (N/m). For laplace/capillary/SE."),
        ("critical_temperature_k", "float", "None", "Critical temp Tc (K). For 'temperature'."),
        ("eutvos_constant_k", "float", "2.1e-7", "Eotvos constant k. Default 2.1e-7."),
        ("molar_mass_g_mol", "float", "None", "Molar mass M (g/mol). For 'temperature'."),
        ("density_kg_m3", "float", "None", "Density rho (kg/m3). For 'temperature'/'capillary'."),
        ("radius_m", "float", "None", "Radius r (m). For 'laplace'/'capillary'."),
        ("contact_angle_deg", "float", "0", "Contact angle theta (deg). For 'capillary'. Default 0."),
        ("gravity", "float", "9.80665", "Gravity g (m/s2). Default 9.80665."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Mode-specific params."),
    ]

    output_sig = [
        ("surface_tension_n_m", "float", "Surface tension gamma (N/m)."),
        ("surface_energy_j_m2", "float", "Surface energy (J/m2)."),
        ("laplace_pressure_pa", "float", "Laplace pressure dP (Pa). Only laplace mode."),
        ("capillary_rise_m", "float", "Capillary rise h (m). Only capillary mode."),
        ("details", "dict", "Extra details."),
        ("analysis_summary", "str", "Summary."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "laplace",
                "temperature_k": 298.15,
                "surface_tension_n_m": 0.0728,
                "critical_temperature_k": None,
                "eutvos_constant_k": 2.1e-7,
                "molar_mass_g_mol": None,
                "density_kg_m3": None,
                "radius_m": 1e-5,
                "contact_angle_deg": 0,
                "gravity": 9.80665,
            },
            "text_input": {
                "input_params": "laplace 0.0728 1e-5",
            },
            "output": {
                "surface_tension_n_m": 0.0728,
                "surface_energy_j_m2": 0.0728,
                "laplace_pressure_pa": 14560.0,
                "capillary_rise_m": None,
                "details": {},
                "analysis_summary": "Spherical bubble (r=10um): dP=14.56 kPa",
            }
        },
        {
            "code_input": {
                "mode": "capillary",
                "temperature_k": 298.15,
                "surface_tension_n_m": 0.0728,
                "critical_temperature_k": None,
                "eutvos_constant_k": 2.1e-7,
                "molar_mass_g_mol": None,
                "density_kg_m3": 998.0,
                "radius_m": 5e-5,
                "contact_angle_deg": 0,
                "gravity": 9.80665,
            },
            "text_input": {
                "input_params": "capillary 0.0728 5e-5 998 0",
            },
            "output": {
                "surface_tension_n_m": 0.0728,
                "surface_energy_j_m2": 0.0728,
                "laplace_pressure_pa": None,
                "capillary_rise_m": 0.297,
                "details": {},
                "analysis_summary": "Water in glass capillary (r=50um): h~29.7 cm",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self): pass

    def _run_base(self, mode, **kw) -> dict:
        mode=mode.lower().strip(); T=kw.get("temperature_k",298.15); g=kw.get("gravity",9.80665)

        if mode=="temperature":
            Tc=kw.get("critical_temperature_k"); M=kw.get("molar_mass_g_mol"); rho=kw.get("density_kg_m3")
            if any(v is None for v in [Tc,M,rho]): raise ChemMCPError("'temperature' needs Tc, M, rho.")
            if T>=Tc: raise ChemMCPError(f"T={T} must be < Tc={Tc}.")
            ke=kw.get("eutvos_constant_k",2.1e-7); Vm=(M/1000)/rho
            gamma=ke*(Tc-T)/(Vm**(2/3))
            return {"surface_tension_n_m":round(gamma,6),"surface_energy_j_m2":round(gamma,6),
                    "laplace_pressure_pa":None,"capillary_rise_m":None,
                    "details":{"method":"Eotvos"},"analysis_summary":f"At {T}K: gamma={round(gamma*1000,3)} mN/m"}

        elif mode=="laplace":
            gamma=kw.get("surface_tension_n_m"); R=kw.get("radius_m")
            if gamma is None or R is None: raise ChemMCPError("'laplace' needs gamma and R.")
            dP=2*gamma/R
            return {"surface_tension_n_m":round(gamma,6),"surface_energy_j_m2":round(gamma,6),
                    "laplace_pressure_pa":round(dP,4),"capillary_rise_m":None,
                    "details":{},"analysis_summary":f"Interface (R={R*1e6:.4g}mm): dP={round(dP,4)} Pa ({round(dP/1000,4)} kPa)"}

        elif mode=="capillary":
            gamma=kw.get("surface_tension_n_m"); r=kw.get("radius_m")
            rho=kw.get("density_kg_m3"); th=kw.get("contact_angle_deg",0)
            if any(v is None for v in [gamma,r,rho]): raise ChemMCPError("'capillary' needs gamma, r, rho.")
            thr=math.radians(th); h=(2*gamma*math.cos(thr))/(rho*g*r)
            return {"surface_tension_n_m":round(gamma,6),"surface_energy_j_m2":round(gamma,6),
                    "laplace_pressure_pa":None,"capillary_rise_m":round(h,6),
                    "details":{},"analysis_summary":f"Capillary (r={r*1e6:.1f}um, th={th}deg): h={round(h*100,2)} cm"}

        elif mode=="surface_energy":
            gamma=kw.get("surface_tension_n_m")
            if gamma is None: raise ChemMCPError("'surface_energy' needs gamma.")
            return {"surface_tension_n_m":round(gamma,6),"surface_energy_j_m2":round(gamma,6),
                    "laplace_pressure_pa":None,"capillary_rise_m":None,
                    "details":{},"analysis_summary":f"gamma={gamma} N/m -> SE={gamma} J/m2"}
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use temperature/lapillary/capillary/surface_energy.")

    def _run_text(self,s:str)->dict:
        try:
            p=s.split(); m=p[0]; kw={"mode":m}
            if m=="temperature":
                kw["critical_temperature_k"]=float(p[1]); kw["molar_mass_g_mol"]=float(p[2]); kw["density_kg_m3"]=float(p[3])
                if len(p)>4: kw["temperature_k"]=float(p[4])
            elif m=="laplace":
                kw["surface_tension_n_m"]=float(p[1]); kw["radius_m"]=float(p[2])
            elif m=="capillary":
                kw["surface_tension_n_m"]=float(p[1]); kw["radius_m"]=float(p[2]); kw["density_kg_m3"]=float(p[3])
                if len(p)>4: kw["contact_angle_deg"]=float(p[4])
            elif m=="surface_energy": kw["surface_tension_n_m"]=float(p[1])
            else: raise ValueError(m)
            return self._run_base(**kw)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")
