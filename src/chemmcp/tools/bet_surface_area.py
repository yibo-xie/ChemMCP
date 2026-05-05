import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BETSurfaceArea(BaseTool):
    """BET specific surface area calculation tool."""
    __version__ = "0.1.0"
    name = "BETSurfaceArea"
    func_name = "bet_surface_area"
    description = "Calculate specific surface area using BET method from N2 adsorption isotherm data."
    implementation_description = "Multi-point BET: linearizes as P/[V(P0-P)] vs P/P0 to obtain Vm and C, then SBET=Vm*NA*sigma/Vmolar."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Adsorption", "Surface Area", "BET", "Materials Science"]
    required_envs = []

    AVOGADRO = 6.02214076e23
    MOLAR_VOLUME_STP = 22414.0

    code_input_sig = [
        ("relative_pressures", "list", "N/A", "Relative pressures P/P0 (typically 0.05-0.35)."),
        ("adsorbed_volumes", "list", "N/A", "Adsorbed volumes (cm3 STP/g) at each P/P0."),
        ("cross_sectional_area", "float", "0.162", "Cross-sectional area sigma (nm2). Default N2 at 77K."),
        ("sample_mass_g", "float", "1.0", "Sample mass (g). Default 1.0."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Format: 'pp0_list(comma) vol_list(comma) [sigma]'."),
    ]

    output_sig = [
        ("monolayer_capacity_cm3_g", "float", "Monolayer capacity Vm (cm3 STP/g)."),
        ("bet_constant_C", "float", "BET constant C."),
        ("specific_surface_area_m2_g", "float", "BET surface area SBET (m2/g)."),
        ("correlation_coefficient_r2", "float", "R2 of BET linear fit."),
        ("bet_plot_data", "list", "BET plot data [(P/P0, P/(V(P0-P))), ...]."),
        ("analysis_summary", "str", "Summary with validity check."),
    ]

    examples = [
        {
            "code_input": {
                "relative_pressures": [0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
                "adsorbed_volumes": [120.0, 130.0, 142.0, 155.0, 170.0, 188.0],
                "cross_sectional_area": 0.162,
                "sample_mass_g": 1.0,
            },
            "text_input": {
                "input_params": "0.05,0.1,0.15,0.2,0.25,0.3 120,130,142,155,170,188",
            },
            "output": {
                "monolayer_capacity_cm3_g": 108.5,
                "bet_constant_C": 95.3,
                "specific_surface_area_m2_g": 472.1,
                "correlation_coefficient_r2": 0.999,
                "bet_plot_data": [[0.05, 0.00044], [0.10, 0.00085]],
                "analysis_summary": "BET: Vm=108.5, C=95.3, S=472.1 m2/g, R2=0.999 Valid",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self): pass

    def _linreg(self, xs, ys):
        n=len(xs); sx=sum(xs); sy=sum(ys); sxy=sum(x*y for x,y in zip(xs,ys))
        sx2=sum(x*x for x in xs); D=n*sx2-sx*sx
        if abs(D)<1e-30: raise ChemMCPError("Cannot regress.")
        b=(n*sxy-sx*sy)/D; a=(sy-b*sx)/n
        ym=sy/n; sst=sum((y-ym)**2 for y in ys); ssr=sum((y-(a+b*x))**2 for x,y in zip(xs,ys))
        return a,b,(1-ssr/sst if sst>0 else 1.0)

    def _run_base(self, relative_pressures, adsorbed_volumes,
                   cross_sectional_area=0.162, sample_mass_g=1.0) -> dict:
        if len(relative_pressures)!=len(adsorbed_volumes):
            raise ChemMCPError("Lists must have equal length.")
        if len(relative_pressures)<3:
            raise ChemMCPError("Need >=3 points.")

        bx=[]; by=[]
        for pp0,v in zip(relative_pressures,adsorbed_volumes):
            if pp0<=0 or pp0>=1 or v<=0: continue
            bx.append(pp0); by.append(pp0/(v*(1-pp0)))
        if len(bx)<3: raise ChemMCPError("Need >=3 valid points.")

        a,b,r2=self._linreg(bx,by)
        if abs(a+b)<1e-30: raise ChemMCPError("Invalid fit.")
        Vm=1/(a+b); C=b/a+1 if abs(a)>1e-30 else float("inf")
        sig_m2=cross_sectional_area*1e-18
        SBET=(Vm*self.AVOGADRO*sig_m2)/self.MOLAR_VOLUME_STP/sample_mass_g

        pdata=[[round(x,6),round(y,6)] for x,y in zip(bx,by)]
        valid=Vm>0 and C>0
        st="Valid" if valid else "Invalid"
        summary=f"BET: Vm={round(Vm,2)}, C={round(C,1)}, S={round(SBET,1)} m2/g, R2={round(r2,4)} {st}"
        logger.info(summary)
        return {"monolayer_capacity_cm3_g":round(Vm,4), "bet_constant_C":round(C,4),
                "specific_surface_area_m2_g":round(SBET,4), "correlation_coefficient_r2":round(r2,6),
                "bet_plot_data":pdata, "analysis_summary":summary}

    def _run_text(self, s:str)->dict:
        try:
            p=s.split()
            pp=[float(x) for x in p[0].split(",")]; vv=[float(x) for x in p[1].split(",")]
            kw={"relative_pressures":pp,"adsorbed_volumes":vv}
            if len(p)>2: kw["cross_sectional_area"]=float(p[2])
            return self._run_base(**kw)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")
