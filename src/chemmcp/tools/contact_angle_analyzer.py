import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ContactAngleAnalyzer(BaseTool):
    """Contact angle and wettability analysis tool."""
    __version__ = "0.1.0"
    name = "ContactAngleAnalyzer"
    func_name = "analyze_contact_angle"
    description = "Analyze contact angle, wettability, solid SFE (Owens-Wendt), adhesion work, and spreading coefficient."
    implementation_description = "Young equation, Owens-Wendt two-liquid method, Dupre: W_SL=gamma_LV*(1+cos(theta)), S=gamma*(cos-1)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Surface Chemistry", "Contact Angle", "Wettability", "Surface Energy"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "'analyze', 'young', 'surface_energy', or 'adhesion'."),
        ("contact_angle_deg", "float", "None", "Contact angle theta (deg). For analyze/adhesion."),
        ("liquid_surface_tension", "float", "None", "Liquid surface tension gamma_LV (N/m). For young/adhesion."),
        ("solid_surface_tension", "float", "None", "Solid surface tension gamma_SV (N/m). For young."),
        ("solid_liquid_tension", "float", "None", "Solid-liquid tension gamma_SL (N/m). For young."),
        ("liquid_dispersion_part", "float", "None", "Dispersion part of gamma_LV. For SE."),
        ("liquid_polar_part", "float", "None", "Polar part of gamma_LV. For SE."),
        ("angle_liquid1_deg", "float", "None", "Contact angle with liquid 1 (non-polar). For SE."),
        ("angle_liquid2_deg", "float", "None", "Contact angle with liquid 2 (polar). For SE."),
        ("gamma_liquid1", "float", "0.0508", "Surface tension of liquid 1. Default CH2I2."),
        ("gamma_liquid2", "float", "0.0728", "Surface tension of liquid 2. Default water."),
        ("d1_fraction", "float", "1.0", "Dispersion fraction of liq1. Default 1.0."),
        ("d2_fraction", "float", "0.215", "Dispersion fraction of liq2. Default water (0.215)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Mode-specific params."),
    ]

    output_sig = [
        ("contact_angle_deg", "float", "Contact angle theta (deg)."),
        ("wettability_class", "str", "Wettability classification."),
        ("cos_theta", "float", "cos(theta) value."),
        ("work_adhesion_j_m2", "float", "Work of adhesion (J/m2). Only adhesion."),
        ("spreading_coefficient", "float", "Spreading coefficient (J/m2). Only adhesion."),
        ("solid_sfe_total", "float", "Total solid SFE (mN/m). Only SE."),
        ("solid_sfe_dispersion", "float", "Dispersive SFE component (mN/m). Only SE."),
        ("solid_sfe_polar", "float", "Polar SFE component (mN/m). Only SE."),
        ("analysis_summary", "str", "Summary."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "analyze",
                "contact_angle_deg": 75.0,
                "liquid_surface_tension": None,
                "solid_surface_tension": None,
                "solid_liquid_tension": None,
                "liquid_dispersion_part": None,
                "liquid_polar_part": None,
                "angle_liquid1_deg": None,
                "angle_liquid2_deg": None,
                "gamma_liquid1": 0.0508,
                "gamma_liquid2": 0.0728,
                "d1_fraction": 1.0,
                "d2_fraction": 0.215,
            },
            "text_input": {
                "input_params": "analyze 75",
            },
            "output": {
                "contact_angle_deg": 75.0,
                "wettability_class": "hydrophobic",
                "cos_theta": 0.2588,
                "work_adhesion_j_m2": None,
                "spreading_coefficient": None,
                "solid_sfe_total": None,
                "solid_sfe_dispersion": None,
                "solid_sfe_polar": None,
                "analysis_summary": "theta=75 deg: hydrophobic, cos=0.259",
            }
        },
        {
            "code_input": {
                "mode": "adhesion",
                "contact_angle_deg": 65.0,
                "liquid_surface_tension": 0.0728,
                "solid_surface_tension": None,
                "solid_liquid_tension": None,
                "liquid_dispersion_part": None,
                "liquid_polar_part": None,
                "angle_liquid1_deg": None,
                "angle_liquid2_deg": None,
                "gamma_liquid1": 0.0508,
                "gamma_liquid2": 0.0728,
                "d1_fraction": 1.0,
                "d2_fraction": 0.215,
            },
            "text_input": {
                "input_params": "adhesion 0.0728 65",
            },
            "output": {
                "contact_angle_deg": 65.0,
                "wettability_class": "hydrophobic",
                "cos_theta": 0.4226,
                "work_adhesion_j_m2": 0.1037,
                "spreading_coefficient": -0.0419,
                "solid_sfe_total": None,
                "solid_sfe_dispersion": None,
                "solid_sfe_polar": None,
                "analysis_summary": "theta=65 on water: Wad=103.7 mJ/m2, S=-41.9 mJ/m2",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self): pass

    def _classify(self, th):
        th = th % 360
        if th <= 5: return "superhydrophilic (theta<=5, complete wetting)"
        elif th <= 30: return "hydrophilic (5<theta<=30, high wettability)"
        elif th <= 90: return "hydrophobic (30<theta<=90, partial wetting)"
        elif th <= 150: return "strongly hydrophobic (90<theta<=150)"
        else: return "superhydrophobic (theta>150, lotus effect)"

    def _run_base(self, mode, **kw) -> dict:
        mode = mode.lower().strip()

        if mode == "analyze":
            th = kw.get("contact_angle_deg")
            if th is None: raise ChemMCPError("'analyze' needs contact_angle_deg.")
            ct = round(math.cos(math.radians(th)), 6); wc = self._classify(th)
            return {"contact_angle_deg": float(th), "wettability_class": wc, "cos_theta": ct,
                    "work_adhesion_j_m2": None, "spreading_coefficient": None,
                    "solid_sfe_total": None, "solid_sfe_dispersion": None, "solid_sfe_polar": None,
                    "analysis_summary": f"theta={th} deg: {wc}, cos={ct}"}

        elif mode == "young":
            glv = kw.get("liquid_surface_tension"); gsv = kw.get("solid_surface_tension")
            gsl = kw.get("solid_liquid_tension")
            if any(v is None for v in [glv, gsv, gsl]): raise ChemMCPError("'young' needs all tensions.")
            if glv <= 0: raise ChemMCPError("gamma_LV must be positive.")
            ct = (gsv - gsl) / glv; ct = max(-1, min(1, ct)); th = math.degrees(math.acos(ct))
            return {"contact_angle_deg": round(th, 4), "wettability_class": self._classify(th),
                    "cos_theta": round(ct, 6), "work_adhesion_j_m2": None, "spreading_coefficient": None,
                    "solid_sfe_total": None, "solid_sfe_dispersion": None, "solid_sfe_polar": None,
                    "analysis_summary": f"Young: cos=({gsv}-{gsl})/{glv}={round(ct,4)} -> theta={round(th,2)} deg"}

        elif mode == "adhesion":
            glv = kw.get("liquid_surface_tension"); th = kw.get("contact_angle_deg")
            if glv is None or th is None: raise ChemMCPError("'adhesion' needs gamma_LV and theta.")
            ct = math.cos(math.radians(th)); Wad = glv * (1 + ct); S = glv * (ct - 1)
            return {"contact_angle_deg": float(th), "wettability_class": self._classify(th),
                    "cos_theta": round(ct, 6), "work_adhesion_j_m2": round(Wad, 6),
                    "spreading_coefficient": round(S, 6),
                    "solid_sfe_total": None, "solid_sfe_dispersion": None, "solid_sfe_polar": None,
                    "analysis_summary": f"theta={th}, gLV={glv}: Wad={round(Wad*1000,2)} mJ/m2, S={round(S*1000,2)} mJ/m2"}

        elif mode == "surface_energy":
            th1 = kw.get("angle_liquid1_deg"); th2 = kw.get("angle_liquid2_deg")
            gl1 = kw.get("gamma_liquid1", 0.0508); gl2 = kw.get("gamma_liquid2", 0.0728)
            d1f = kw.get("d1_fraction", 1.0); d2f = kw.get("d2_fraction", 0.215)
            if any(v is None for v in [th1, th2]): raise ChemMCPError("'SE' needs both angles.")
            gl1d = gl1 * d1f; gl1p = gl1 * (1 - d1f); gl2d = gl2 * d2f; gl2p = gl2 * (1 - d2f)
            c1 = math.cos(math.radians(th1)); c2 = math.cos(math.radians(th2))
            A1 = math.sqrt(gl1d); B1 = math.sqrt(gl1p); A2 = math.sqrt(gl2d); B2 = math.sqrt(gl2p)
            Y1 = (1 + c1) * math.sqrt(gl1) / 2; Y2 = (1 + c2) * math.sqrt(gl2) / 2
            det = A1 * B2 - A2 * B1
            if abs(det) < 1e-20: raise ChemMCPError("Cannot solve: liquids too similar.")
            a = (Y1 * B2 - Y2 * B1) / det; b = (A1 * Y2 - A2 * Y1) / det
            gsd = max(a * a, 0); gsp = max(b * b, 0); gst = gsd + gsp
            return {"contact_angle_deg": None, "wettability_class": self._classify(th2) if th2 else None,
                    "cos_theta": None, "work_adhesion_j_m2": None, "spreading_coefficient": None,
                    "solid_sfe_total": round(gst, 4), "solid_sfe_dispersion": round(gsd, 4),
                    "solid_sfe_polar": round(gsp, 4),
                    "analysis_summary": f"Owens-Wendt: gS_d={round(gsd,2)}, gS_p={round(gsp,2)}, total={round(gst,2)} mN/m"}
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use analyze/young/adhesion/surface_energy.")

    def _run_text(self, s: str) -> dict:
        try:
            p = s.split(); m = p[0]; kw = {"mode": m}
            if m == "analyze": kw["contact_angle_deg"] = float(p[1])
            elif m == "young":
                kw["liquid_surface_tension"] = float(p[1]); kw["solid_surface_tension"] = float(p[2])
                kw["solid_liquid_tension"] = float(p[3])
            elif m == "adhesion":
                kw["liquid_surface_tension"] = float(p[1]); kw["contact_angle_deg"] = float(p[2])
            elif m == "surface_energy":
                kw["angle_liquid1_deg"] = float(p[1]); kw["angle_liquid2_deg"] = float(p[2])
                if len(p) > 3: kw["gamma_liquid1"] = float(p[3])
                if len(p) > 4: kw["gamma_liquid2"] = float(p[4])
            else: raise ValueError(m)
            return self._run_base(**kw)
        except Exception as e:
            raise ChemMCPError(f"Parse error: {e}")
