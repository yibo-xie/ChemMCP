import logging
import math
from typing import List, Optional, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Typical J-coupling ranges (Hz) based on number of bonds and structural context
# Format: (J_min_hz, J_max_hz, typical_hz, notes)
J_COUPLING_DATA = {
    # ¹H-¹H couplings
    "1J_HH": (250, 320, 280, "Direct bond — not observed in solution NMR (decoupled)"),
    "2J_HH_geminal": (-25, -10, -15, "Geminal protons; sign depends on H-C-H angle"),
    "3J_HH_vicinal_sp3": (5, 14, 7, "Vicinal, free rotation (Karplus average ~7 Hz)"),
    "3J_HH_vicinal_sp2_trans": (12, 19, 16, "Trans alkene (dihedral ≈ 180°)"),
    "3J_HH_vicinal_sp2_cis": (6, 13, 10, "Cis alkene (dihedral ≈ 0°)"),
    "3J_HH_vicinal_aromatic_ortho": (6, 9, 8, "Ortho aromatic"),
    "4J_HH_meta_aromatic": (1, 4, 2, "Meta aromatic (long-range)"),
    "5J_HH_para_aromatic": (0, 1, 0.5, "Para aromatic (very weak)"),
    # ¹³C-¹H couplings
    "1J_CH_sp3": (115, 135, 125, "C-H one-bond sp³"),
    "1J_CH_sp2": (150, 170, 160, "C-H one-bond sp² (alkene/aromatic)"),
    "1J_CH_sp": (245, 255, 250, "C-H one-bond sp (terminal alkyne)"),
    "2J_CH": (-5, 5, 1, "Two-bond C-H (geminal)"),
    "3J_CH_vicinal": (0, 10, 5, "Three-bond C-H vicinal"),
    # ¹³C-¹³C couplings (rarely observed)
    "1J_CC_sp3": (30, 55, 40, "C-C one-bond sp³"),
    "1J_CC_sp2": (50, 70, 60, "C-C one-bond sp²"),
    # Couplings to heteronuclei
    "1J_CF": (150, 360, 270, "C-F one-bond; large due to high s-character"),
    "1J_CP": (125, 280, 200, "C-P one-bond"),
    "2J_HF": (30, 65, 48, "H-C-F two-bond (geminal)"),
    "3J_HF_vicinal": (0, 25, 10, "H-C-C-F three-bond"),
}

# Karplus equation parameters for ³JHH coupling
# J(θ) = A·cos²(θ) + B·cos(θ) + C
KARPLUS_PARAMS = {
    "sp3_default": {"A": 7.0, "B": -1.0, "C": 0.5},   # Free-rotating average
    "sp2_trans": {"A": 9.5, "B": -0.3, "C": 0.8},     # Rigid trans alkene
    "sp2_cis": {"A": 7.2, "B": -1.8, "C": 1.2},       # Rigid cis alkene
    "aromatic_ortho": {"A": 8.0, "B": -0.5, "C": 0.6},  # Aromatic ortho
}


@ChemMCPManager.register_tool
class SpinSpinCoupling(BaseTool):
    """
    自旋-自旋耦合常数（NMR J耦合）计算工具。
    基于键数、杂化方式和二面角（Karplus关系）计算J耦合常数（Hz）。
    支持¹H-¹H、¹³C-¹H、¹³C-¹³C及异核耦合。
    """
    __version__ = "0.1.0"
    name = "SpinSpinCoupling"
    func_name = "calculate_spin_spin_coupling"
    description = "Calculate NMR spin-spin coupling constants (J, in Hz) based on number of bonds, hybridization, and dihedral angle using Karplus relation for vicinal couplings."
    implementation_description = (
        "Uses empirical J-coupling databases for different bond numbers (¹J–⁵J) and nuclear pairs "
        "(¹H-¹H, ¹³C-¹H, ¹³C-¹³C, heteronuclear). For ³JHH vicinal couplings, applies Karplus equation: "
        "J(θ) = A·cos²θ + B·cosθ + C with parameter sets for sp³/sp²/aromatic systems."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["NMR", "J-Coupling", "Spin-Spin Coupling", "Karplus Equation", "Spectroscopy"]
    required_envs = []

    code_input_sig = [
        ("coupled_nuclei", "str", "N/A", "Coupled nuclear pair: 'HH', 'CH', 'CC', 'HF', 'HP', etc."),
        ("num_bonds", "int", "N/A", "Number of bonds between coupled nuclei (1–5)."),
        ("hybridization", "str", "sp3", "Hybridization of intervening atoms: 'sp3', 'sp2', 'sp'."),
        ("geometry", "str", "default", "Geometry/context: 'default', 'trans', 'cis', 'aromatic_ortho', 'aromatic_meta', 'geminal'."),
        ("dihedral_angle_deg", "float", "None", "Dihedral angle θ in degrees for Karplus calculation (optional, only for 3J vicinal)."),
        ("detail_level", "str", "standard", "Detail level: 'basic', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: nuclei num_bonds [hybridization] [geometry] [dihedral_angle]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing J_hz, coupling_type, number_of_bonds, karplus_analysis, splitting_pattern."),
    ]

    examples = [
        {
            "code_input": {
                "coupled_nuclei": "HH",
                "num_bonds": 3,
                "hybridization": "sp3",
                "geometry": "default",
                "dihedral_angle_deg": None,
                "detail_level": "standard",
            },
            "text_input": {
                "input_params": "HH 3 sp3 default",
            },
            "output": {
                "result": {
                    "coupled_nuclei": "¹H-¹H",
                    "j_coupling_hz": 7.0,
                    "coupling_type": "³JHH vicinal (sp³, free rotation)",
                    "num_bonds": 3,
                    "splitting_pattern": "quartet/triplet (n+1 rule)",
                    "karplus_analysis": {
                        "dihedral_angle_used": "averaged over free rotation",
                        "average_J": 7.0,
                        "karplus_equation": "J(θ) = 7.0·cos²(θ) − 1.0·cos(θ) + 0.5",
                        "note": "Free rotation averages to ~7 Hz",
                    },
                    "typical_range_hz": (5, 14),
                }
            }
        },
        {
            "code_input": {
                "coupled_nuclei": "HH",
                "num_bonds": 3,
                "hybridization": "sp2",
                "geometry": "trans",
                "dihedral_angle_deg": 180.0,
                "detail_level": "detailed",
            },
            "text_input": {
                "input_params": "HH 3 sp2 trans 180 detailed",
            },
            "output": {
                "result": {
                    "coupled_nuclei": "¹H-¹H",
                    "j_coupling_hz": 10.6,
                    "coupling_type": "³JHH vicinal (sp², trans)",
                    "num_bonds": 3,
                    "splitting_pattern": "doublet of doublets (large J)",
                    "karplus_analysis": {
                        "dihedral_angle_deg": 180.0,
                        "J_calc": 10.6,
                        "karplus_equation": "J(θ) = 9.5·cos²(θ) − 0.3·cos(θ) + 0.8",
                        "max_possible_J": 10.6,
                        "note": "Maximum at θ=180° (trans configuration)",
                    },
                    "stereochemical_info": "Large trans J confirms E (trans) alkene geometry.",
                    "typical_range_hz": (12, 19),
                }
            }
        },
        {
            "code_input": {
                "coupled_nuclei": "CH",
                "num_bonds": 1,
                "hybridization": "sp2",
                "geometry": "default",
                "dihedral_angle_deg": None,
                "detail_level": "basic",
            },
            "text_input": {
                "input_params": "CH 1 sp2 basic",
            },
            "output": {
                "result": {
                    "coupled_nuclei": "¹³C-¹H",
                    "j_coupling_hz": 160,
                    "coupling_type": "¹JCH one-bond (sp²)",
                    "num_bonds": 1,
                    "note": "One-bond CH coupling observed in proton-coupled ¹³C spectra.",
                    "typical_range_hz": (150, 170),
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.j_db = dict(J_COUPLING_DATA)
        self.karplus_params = dict(KARPLUS_PARAMS)

    def _run_base(
        self,
        coupled_nuclei: str,
        num_bonds: int,
        hybridization: str = "sp3",
        geometry: str = "default",
        dihedral_angle_deg: Optional[float] = None,
        detail_level: str = "standard",
    ) -> dict:
        """Core logic: calculate J-coupling constant."""
        nuclei = coupled_nuclei.upper().strip()
        nb = num_bonds

        if nb < 1 or nb > 5:
            raise ChemMCPError("num_bonds must be between 1 and 5.")

        hyb = hybridization.lower().strip()
        geo = geometry.lower().strip()
        dl = detail_level.lower()

        # Build lookup key
        key = self._build_lookup_key(nuclei, nb, hyb, geo)

        if nb == 3 and nuclei == "HH" and dihedral_angle_deg is not None:
            # Use Karplus equation for explicit dihedral angle
            j_val = self._karplus_calc(dihedral_angle_deg, hyb, geo)
            karplus_data = self._full_karplus_analysis(dihedral_angle_deg, hyb, geo)
        else:
            # Look up from database
            db_entry = self.j_db.get(key)
            if db_entry:
                j_min, j_max, j_typical, notes = db_entry
                if geo == "trans" or geo == "cis":
                    j_val = j_typical
                elif nb == 3 and nuclei == "HH":
                    j_val = j_typical
                else:
                    j_val = j_typical
                karplus_data = None
            else:
                # Estimate from general trends
                j_val, notes = self._estimate_j(nuclei, nb, hyb, geo)
                karplus_data = None

        j_val = round(j_val, 1)

        result = {
            "coupled_nuclei": f"^{self._format_nucleus(nuclei[0])}-{self._format_nucleus(nuclei[1])}",
            "j_coupling_hz": j_val,
            "coupling_type": self._describe_coupling(nuclei, nb, hyb, geo),
            "num_bonds": nb,
            "splitting_pattern": self._predict_splitting(nb, j_val, nuclei),
            "typical_range_hz": self._get_range(key, nuclei, nb),
        }

        if karplus_data:
            result["karplus_analysis"] = karplus_data

        if dl == "detailed":
            result["structural_implications"] = self._structural_implications(nuclei, nb, j_val, hyb, geo)
            result["sign_information"] = self._sign_info(nuclei, nb)

        return {"result": result}

    def _build_lookup_key(self, nuclei: str, nb: int, hyb: str, geo: str) -> str:
        """Build database lookup key."""
        n_str = f"{nuclei[0]}{nuclei[1]}" if len(nuclei) == 2 else nuclei
        prefix = f"{nb}J_{n_str}"

        if nb == 2:
            return f"{prefix}_geminal"
        elif nb == 3:
            if geo == "trans" and hyb == "sp2":
                return f"{prefix}_vicinal_sp2_trans"
            elif geo == "cis" and hyb == "sp2":
                return f"{prefix}_vicinal_sp2_cis"
            elif geo == "aromatic_ortho":
                return f"{prefix}_vicinal_aromatic_ortho"
            elif geo == "aromatic_meta":
                return f"{prefix}_meta_aromatic"
            elif hyb == "sp3":
                return f"{prefix}_vicinal_sp3"
            return f"{prefix}_vicinal_{hyb}"
        elif nb == 4:
            if geo == "aromatic_meta":
                return f"{prefix}_meta_aromatic"
            return f"{prefix}"
        elif nb == 5:
            if geo == "aromatic" or geo == "para":
                return f"{prefix}_para_aromatic"
            return f"{prefix}"
        else:
            return f"{prefix}_{hyb}"

    def _karplus_calc(self, theta_deg: float, hyb: str, geo: str) -> float:
        """Calculate J using Karplus equation: J(θ) = A·cos²θ + B·cosθ + C."""
        theta_rad = math.radians(theta_deg)
        cos_t = math.cos(theta_rad)

        # Select parameter set
        if geo == "trans" and hyb == "sp2":
            params = self.karplus_params["sp2_trans"]
        elif geo == "cis" and hyb == "sp2":
            params = self.karplus_params["sp2_cis"]
        elif geo == "aromatic_ortho":
            params = self.karplus_params["aromatic_ortho"]
        else:
            params = self.karplus_params["sp3_default"]

        A, B, C = params["A"], params["B"], params["C"]
        return A * cos_t ** 2 + B * cos_t + C

    def _full_karplus_analysis(self, theta_deg: float, hyb: str, geo: str) -> dict:
        """Full Karplus analysis with curve characteristics."""
        params = self.karplus_params.get(f"{hyb}_{geo}", self.karplus_params.get("sp3_default", self.karplus_params["sp3_default"]))
        A, B, C = params["A"], params["B"], params["C"]

        theta_rad = math.radians(theta_deg)
        cos_t = math.cos(theta_rad)
        j_val = A * cos_t ** 2 + B * cos_t + C

        # Find max and min of Karplus curve
        # dJ/dθ = -2A·sinθ·cosθ - B·sinθ = -sinθ·(2A·cosθ + B) = 0
        # sinθ = 0 → θ = 0°, 180°
        # 2A·cosθ + B = 0 → cosθ = -B/(2A)
        j_at_0 = A * 1 + B * 1 + C
        j_at_180 = A * 1 - B * + C
        cos_extreme = -B / (2 * A) if abs(A) > 0.01 else 0
        cos_extreme = max(-1, min(1, cos_extreme))
        j_min = A * cos_extreme ** 2 + B * cos_extreme + C

        return {
            "dihedral_angle_deg": round(theta_deg, 1),
            "J_calc": round(j_val, 1),
            "karplus_equation": f"J(θ) = {A}·cos²(θ) + {B}·cos(θ) + {C}",
            "J_at_0deg": round(j_at_0, 1),
            "J_at_180deg": round(j_at_180, 1),
            "J_minimum": round(j_min, 1),
            "note": self._karplus_note(theta_deg, j_val),
        }

    @staticmethod
    def _karplus_note(theta: float, j: float) -> str:
        if theta >= 165 or theta <= 15:
            return f"Near extreme of Karplus curve (θ={theta:.0f}°). Large |J| indicates well-defined geometry."
        elif theta >= 60 and theta <= 120:
            return f"Near minimum of Karplus curve (θ≈90°). Small J typical for gauche/perpendicular arrangement."
        else:
            return f"Intermediate dihedral angle gives moderate J value."

    @staticmethod
    def _format_nucleus(n: str) -> str:
        mapping = {"H": "1H", "C": "13C", "F": "19F", "P": "31P", "N": "15N"}
        return mapping.get(n, n)

    @staticmethod
    def _describe_coupling(nuclei: str, nb: int, hyb: str, geo: str) -> str:
        bond_names = {1: "one-bond", 2: "two-bond (geminal)", 3: "three-bond (vicinal)",
                      4: "four-bond (long-range)", 5: "five-bond (long-range)"}
        bname = bond_names.get(nb, f"{nb}-bond")
        parts = [f"^{nb}J{nuclei} {bname}"]
        if hyb != "default":
            parts.append(f"({hyb})")
        if geo != "default":
            parts.append(f", {geo}")
        return " ".join(parts)

    @staticmethod
    def _predict_splitting(nb: int, j: float, nuclei: str) -> str:
        if nb == 1:
            return "Doublet (if observed)"
        elif nb == 2:
            return "Doublet of doubleets (or quartet if equivalent)"
        elif nb == 3:
            if j > 11:
                return "Doublet of doublets (large J, distinct splitting)"
            elif j > 3:
                return "Triplet/quartet pattern (moderate J)"
            else:
                return "Broadened singlet or small splitting (small J)"
        elif nb >= 4:
            return f"Long-range coupling ({nb}J), often appears as line broadening"
        return "Complex multiplet"

    def _get_range(self, key: str, nuclei: str, nb: int) -> Tuple[float, float]:
        entry = self.j_db.get(key)
        if entry:
            return (entry[0], entry[1])
        # Fallback ranges
        fallback = {
            1: (100, 300), 2: (-25, 20), 3: (0, 20), 4: (0, 5), 5: (0, 2)
        }
        return fallback.get(nb, (0, 50))

    def _estimate_j(self, nuclei: str, nb: int, hyb: str, geo: str) -> tuple:
        """Estimate J when no exact DB match."""
        base = {1: 150, 2: -12, 3: 7, 4: 1, 5: 0.3}.get(nb, 1)
        if "F" in nuclei:
            base *= 3
        note = f"Estimated value for ^{nb}J{nuclei}."
        return base, note

    def _structural_implications(self, nuclei: str, nb: int, j: float, hyb: str, geo: str) -> str:
        implications = []
        if nb == 3 and nuclei == "HH":
            if j > 14:
                implications.append("Large vicinal J suggests anti-periplanar arrangement (θ ≈ 180°).")
            elif j < 3:
                implications.append("Small vicinal J suggests gauche arrangement (θ ≈ 60°) or orthogonal orientation.")
            else:
                implications.append("Moderate vicinal J suggests intermediate dihedral angle.")
            if geo == "trans":
                implications.append("Trans alkene confirmed by large J (>12 Hz).")
            elif geo == "cis":
                implications.append("Cis alkene indicated by smaller J (6–13 Hz vs trans 12–19 Hz).")
        elif nb == 2 and nuclei == "HH":
            if j < -15:
                implications.append("Large geminal |J| indicates strained ring or electronegative substituent effects.")
        elif nb == 1 and "C" in nuclei and "H" in nuclei:
            if j > 200:
                implications.append("Very large ¹JCH suggests sp-hybridized carbon (alkyne).")
            elif j > 155:
                implications.append("Large ¹JCH consistent with sp² carbon (alkene/aromatic).")
            else:
                implications.append("Typical sp³ ¹JCH range (115–135 Hz).")
        return " ".join(implications) if implications else "Standard coupling."

    @staticmethod
    def _sign_info(nuclei: str, nb: int) -> str:
        signs = {
            (1, "CH"): "Positive (¹JCH > 0 always)",
            (2, "HH"): "Negative (²JHH < 0 typically, -10 to -20 Hz)",
            (3, "HH"): "Positive (³JHH > 0 typically)",
            (1, "CC"): "Positive (¹JCC > 0)",
        }
        return signs.get((nb, nuclei), "Sign depends on specific molecular context.")

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ChemMCPError("Need at least: coupled_nuclei num_bonds")

            kwargs = {
                "coupled_nuclei": parts[0],
                "num_bonds": int(parts[1]),
            }
            modes = {"basic", "standard", "detailed"}
            geos = {"default", "trans", "cis", "aromatic_ortho", "aromatic_meta", "geminal"}
            hybs = {"sp3", "sp2", "sp"}

            for p in parts[2:]:
                pl = p.lower()
                if pl in modes:
                    kwargs["detail_level"] = pl
                elif pl in geos:
                    kwargs["geometry"] = pl
                elif pl in hybs:
                    kwargs["hybridization"] = pl
                else:
                    try:
                        kwargs["dihedral_angle_deg"] = float(p)
                    except ValueError:
                        pass

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
