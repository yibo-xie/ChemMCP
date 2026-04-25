import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CouplingConstantAnalyzer(BaseTool):
    """
    偶合常数分析工具。
    基于 Karplus 方程分析偶合常数与二面角/构象的关系。
    """
    __version__ = "0.1.0"
    name = "CouplingConstantAnalyzer"
    func_name = "analyze_coupling_constant"
    description = "Analyze NMR coupling constants (J values) and their relationship to molecular conformation. Uses Karplus equation for vicinal coupling prediction."
    implementation_description = "Implements the Karplus equation J(θ) = A·cos²θ + B·cosθ + C to relate vicinal coupling constants (³JHH) to dihedral angles. Also provides reference ranges for geminal, vicinal, long-range, and heteronuclear couplings."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["NMR", "Coupling Constant", "Karplus Equation", "Conformation", "Spectroscopy"]
    required_envs = []

    code_input_sig = [
        ("dihedral_angle", "float", "N/A", "Dihedral angle in degrees between coupled protons (0-180°)."),
        ("coupling_type", "str", "vicinal", "Type of coupling: 'vicinal' (³J), 'geminal' (²J), 'long_range' (⁴J, ⁵J), or 'heteronuclear'."),
        ("system_type", "str", "general", "System type: 'general', 'allylic', 'aromatic', 'alicyclic', 'carbocycle', 'H-C-C-X'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'dihedral_angle [coupling_type] [system_type]'. Example: '60 vicinal general'"),
    ]

    output_sig = [
        ("analysis", "dict", "Complete analysis including predicted J value, conformation interpretation, and Karplus parameters."),
    ]

    examples = [
        {
            "code_input": {"dihedral_angle": 60.0, "coupling_type": "vicinal", "system_type": "general"},
            "text_input": {"input_params": "60"},
            "output": {
                "analysis": {
                    "predicted_J_hz": 2.8,
                    "conformation": "gauche",
                }
            },
        },
        {
            "code_input": {"dihedral_angle": 180.0, "coupling_type": "vicinal", "system_type": "general"},
            "text_input": {"input_params": "180 vicinal general"},
            "output": {
                "analysis": {
                    "predicted_J_hz": 12.5,
                    "conformation": "anti-periplanar (trans)",
                }
            },
        },
    ]

    # ========== KARPLUS PARAMETERS FOR DIFFERENT SYSTEMS ==========
    # Format: {system_type: (A, B, C) for J(θ) = A·cos²θ + B·cosθ + C}
    _KARPLUS_PARAMS = {
        # General H-C-C-H system (original Karplus, refined)
        "general": (7.0, -1.0, 0.0),
        "general_refined": (13.7, -0.73, 0.0),   # More accurate modern parametrization
        "classic": (8.5, -0.3, 0.0),              # Classic Karplus 1963

        # Substituent-dependent systems
        "H-C-C-X": (8.4, -0.6, 0.0),             # One proton replaced by substituent X
        "X-C-C-X": (9.5, -1.0, 0.0),             # Both protons replaced by substituents

        # Specific structural contexts
        "allylic": (10.2, -0.8, 0.0),             # H-C=C-C-H allylic coupling
        "aromatic_ortho": (8.0, -0.5, 0.0),       # Ortho aromatic H-C-C-H
        "alicyclic": (8.0, -1.0, 0.3),            # Cyclohexane-type chair conformations
        "carbocycle_6ring": (8.0, -1.0, 0.3),     # 6-membered carbocycle
        "carbocycle_5ring": (6.5, -0.8, 0.5),     # 5-membered ring
        "olefinic": (11.0, -1.4, 0.4),            # H-C=C-C= type systems
        "protein_backbone": (7.9, -1.05, 0.65),   # H-N-Cα-C-H protein backbone (φ angle)
        "dna_sugar": (9.5, -1.0, 0.3),            # DNA ribose/deoxyribose sugar pucker
        "ethane_gauche": (8.0, -1.0, 0.0),        # Ethane-like staggered
    }

    # Reference coupling constant ranges (Hz)
    _REFERENCE_RANGES = {
        "geminal_2J_HH_aliphatic": (-12, -15, "Aliphatic geminal ²JHH; depends on H-C-H bond angle"),
        "geminal_2J_HH_olefinic": (0, 3, "Olefinic geminal ²JHH (=CH₂)"),
        "geminal_2J_HH_aromatic": (1, 3, "Aromatic/heteroaromatic geminal"),
        "vicinal_3J_HH_sp3 gauche": (2, 4, "Gauche (θ≈60°): small coupling"),
        "vicinal_3J_HH_sp3_antiperiplanar": (8, 14, "Anti-periplanar (θ≈180°): large coupling"),
        "vicinal_3J_HH_sp3_syn": (2, 5, "Syn (θ≈0°): moderate coupling"),
        "vicinal_3J_HH_olefinic_cis": (6, 12, "Cis olefinic (Z-alkene)"),
        "vicinal_3J_HH_olefinic_trans": (12, 18, "Trans olefinic (E-alkene)"),
        "vicinal_3J_HH_aromatic_ortho": (6, 9, "Ortho aromatic coupling"),
        "vicinal_3J_HH_aromatic_meta": (1, 3, "Meta aromatic coupling (through-bond)"),
        "vicinal_3J_HH_aromatic_para": (0, 1, "Para aromatic coupling (very small)"),
        "long_range_4J_W_coupling": (0, 3, "W-path (allylic/benzilic) long-range coupling"),
        "long_range_4J_M_coupling": (0, 1, "M-path (meta aromatic) coupling"),
        "long_range_5J_allenic": (0, 7, "Allenic homoallylic coupling"),
        "heteronuclear_1J_CH": (125, 250, "One-bond C-H coupling (¹JCH); ~125 Hz sp³, ~150 Hz sp², ~250 Hz sp"),
        "heteronuclear_1J_CF": (150, 370, "One-bond C-F coupling (¹JCF); ~45 Hz aliphatic, ~285 Hz aryl"),
        "heteronuclear_2J_CHF": (0, 60, "Two-bond H-C-F coupling (²JHCF)"),
        "heteronuclear_1J_CC": (30, 80, "One-bond C-C coupling (direct detection only)"),
        "heteronuclear_1J_CN": (10, 20, "One-bond C-N coupling"),
        "heteronuclear_2J_HH_C=O": (0, 3, "Two-bond coupling across carbonyl (W-coupling)"),
    }

    # Conformation interpretation guide
    _CONFORMATION_GUIDE = {
        0: ("syn-periplanar (eclipsed)", "Eclipsed conformation; high energy; rarely observed in flexible molecules"),
        30: ("synclinal (gauche+)", "Gauche conformation; common in staggered rotamers"),
        60: ("gauche (±)", "Classic gauche; typical for freely rotating bonds"),
        90: ("orthogonal/clinal (90°)", "Perpendicular dihedral; minimal orbital overlap"),
        120: ("anticlinal (gauche-)", "Approaching anti; intermediate coupling"),
        150: ("anti-clinal", "Near anti-periplanar"),
        180: ("anti-periplanar (trans)", "Anti/trans; maximum orbital overlap for σ bond formation"),
    }

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _karplus_equation(self, theta_deg: float, params: tuple) -> float:
        """
        Calculate J using Karplus equation.

        J(θ) = A·cos²(θ) + B·cos(θ) + C

        Args:
            theta_deg: Dihedral angle in degrees
            params: (A, B, C) coefficients

        Returns:
            Coupling constant J in Hz
        """
        theta_rad = math.radians(theta_deg)
        cos_theta = math.cos(theta_rad)
        A, B, C = params
        J = A * (cos_theta ** 2) + B * cos_theta + C
        return round(J, 2)

    def _interpret_conformation(self, theta_deg: float) -> dict:
        """Interpret the dihedral angle in terms of molecular conformation."""
        # Find closest standard angle
        theta = theta_deg % 360
        if theta > 180:
            theta = 360 - theta  # symmetry around 180°

        best_match = None
        best_diff = 999
        for std_angle, (conf_name, desc) in self._CONFORMATION_GUIDE.items():
            diff = abs(theta - std_angle)
            if diff < best_diff:
                best_diff = diff
                best_match = (std_angle, conf_name, desc)

        return {
            "input_angle": round(theta_deg, 1),
            "effective_angle": round(theta, 1),
            "nearest_standard": best_match[0],
            "conformation_name": best_match[1],
            "description": best_match[2],
            "deviation_from_standard": round(best_diff, 1),
        }

    def _get_reference_range(self, coupling_type: str, j_value: float) -> list:
        """Find relevant reference ranges for a given J value."""
        relevant = []
        for key, (lo, hi, desc) in self._REFERENCE_RANGES.items():
            if coupling_type.lower() in key.lower() or coupling_type == "general":
                relevant.append({
                    "range_name": key,
                    "range_lo": lo,
                    "range_hi": hi,
                    "description": desc,
                    "j_in_range": lo <= abs(j_value) <= hi,
                })
        return sorted(relevant, key=lambda x: x["j_in_range"], reverse=True)

    def _run_base(self, dihedral_angle: float, coupling_type: str = "vicinal",
                  system_type: str = "general") -> dict:
        """
        Analyze coupling constant.

        Args:
            dihedral_angle: Dihedral angle in degrees
            coupling_type: Type of coupling
            system_type: Molecular system type for Karplus parameter selection

        Returns:
            Dict with full analysis
        """
        if not isinstance(dihedral_angle, (int, float)):
            raise ChemMCPError("dihedral_angle must be a number.")

        # Get Karplus parameters
        params_key = system_type if system_type in self._KARPLUS_PARAMS else "general_refined"
        params = self._KARPLUS_PARAMS.get(system_type, self._KARPLUS_PARAMS["general_refined"])

        # Calculate J value
        if coupling_type.lower() == "vicinal":
            j_predicted = self._karplus_equation(dihedral_angle, params)
            method = f"Karplus equation with {system_type} parameters (A={params[0]}, B={params[1]}, C={params[2]})"
        elif coupling_type.lower() == "geminal":
            # Geminal coupling is angle-independent (depends on bond angle)
            j_predicted = round(-12.0 + (abs(dihedral_angle - 109.5) / 100) * 5, 1)
            method = "Geminal coupling estimate based on H-C-H bond angle deviation from tetrahedral"
        elif coupling_type.lower() in ("long_range", "4j", "5j"):
            # Long-range coupling (small values)
            j_predicted = round(abs(math.sin(math.radians(dihedral_angle))) * 2.5, 1)
            method = "Long-range coupling estimate (through-space/through-bond)"
        elif coupling_type.lower() == "heteronuclear":
            j_predicted = round(145.0 + 50 * math.cos(math.radians(dihedral_angle)), 1)
            method = "Heteronuclear coupling estimate (¹JCH range)"
        else:
            j_predicted = self._karplus_equation(dihedral_angle, params)
            method = f"Karplus equation (default: {system_type})"

        # Conformation interpretation
        conform = self._interpret_conformation(dihedral_angle)

        # Reference ranges
        ref_ranges = self._get_reference_range(coupling_type, j_predicted)

        return {
            "analysis": {
                "input_dihedral_angle_deg": float(dihedral_angle),
                "coupling_type": coupling_type,
                "system_type": system_type,
                "predicted_J_hz": j_predicted,
                "method": method,
                "conformation": conform,
                "reference_ranges": ref_ranges[:5],
                "karplus_parameters_used": {
                    "A": params[0],
                    "B": params[1],
                    "C": params[2],
                    "equation": f"J({dihedral_angle}°) = {params[0]}·cos²({dihedral_angle}°) + {params[1]}·cos({dihedral_angle}°) + {params[2]} = {j_predicted} Hz",
                },
                "note": (
                    f"Vicinal coupling ³JHH follows Karplus relationship with dihedral angle.\n"
                    f"Typical ranges: syn (0°): 0-4 Hz | gauche (60°): 2-4 Hz | anti (180°): 8-14 Hz\n"
                    f"Actual J values depend on substituents, electronegativity, and bond length."
                ),
            }
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'angle [type] [system]'")

        try:
            angle = float(parts[0])
        except ValueError:
            raise ChemMCPError(f"Invalid angle: '{parts[0]}' must be a number")

        ctype = parts[1] if len(parts) > 1 else "vicinal"
        stype = parts[2] if len(parts) > 2 else "general"

        return self._run_base(angle, ctype, stype)
