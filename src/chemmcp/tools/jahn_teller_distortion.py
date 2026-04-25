import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class JahnTellerDistortion(BaseTool):
    """
    预测 Jahn-Teller 畸变的类型和程度。
    基于 Jahn-Teller 定理：非线性分子在对称性简并电子组态下会发生几何畸变以消除简并。
    覆盖八面体（Oh）和四面体（Td）配合物，针对 d^4-d^7（高/低自旋）和 d^8-d^9 组态。
    """
    __version__ = "0.1.0"
    name = "JahnTellerDistortion"
    func_name = "predict_jahn_teller"
    description = "Predict Jahn-Teller distortion type and extent for octahedral or tetrahedral coordination complexes based on d-electron configuration, spin state, and geometry."
    implementation_description = "Uses Jahn-Teller theorem rules: identifies orbitally degenerate ground states in Oh/Td geometries and predicts distortion direction (elongation/compression) and magnitude (strong/weak/none)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Coordination Chemistry", "Jahn-Teller", "Crystal Field", "d-electrons", "Distortion"]
    required_envs = []

    code_input_sig = [
        ("d_electron_count", "int", "N/A", "Number of d electrons (0-10)."),
        ("geometry", "str", "octahedral", "Coordination geometry: 'octahedral' or 'tetrahedral'."),
        ("spin_state", "str", "unknown", "Spin state: 'high_spin', 'low_spin', or 'unknown' (will analyze both)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'd_electron_count geometry spin_state'. Example: '9 octahedral high_spin'."),
    ]

    output_sig = [
        ("distortion_type", "str", "Type of distortion: 'strong_elongation', 'weak_elongation', 'strong_compression', 'weak_compression', or 'none'."),
        ("explanation", "str", "Detailed explanation of which orbitals are degenerate and how the distortion lifts the degeneracy."),
        ("electron_config", "str", "The d-orbital electron configuration (e.g., '(t2g)^6(eg)^3')."),
        ("examples", "str", "Well-known examples of complexes with this configuration."),
    ]

    examples = [
        {
            "code_input": {
                "d_electron_count": 9,
                "geometry": "octahedral",
                "spin_state": "high_spin",
            },
            "text_input": {
                "input_params": "9 octahedral high_spin"
            },
            "output": {
                "distortion_type": "strong_elongation",
                "explanation": "Cu(II) d^9: (t2g)^6(eg)^3. The eg orbital set is singly degenerate (one electron in dx2-y2, two in dz2). Elongation along z-axis lowers dz2 energy and raises dx2-y2, stabilizing the system.",
                "electron_config": "(t2g)^6(eg)^3",
                "examples": "[Cu(H2O)6]2+, [Cu(NH3)4(H2O)2]2+, CuCl2, CuO (all show characteristic elongated octahedra)",
            }
        },
        {
            "code_input": {
                "d_electron_count": 4,
                "geometry": "octahedral",
                "spin_state": "low_spin",
            },
            "text_input": {
                "input_params": "4 octahedral low_spin"
            },
            "output": {
                "distortion_type": "strong_elongation",
                "explanation": "Low-spin d^4 (e.g., Cr(II), Mn(III)): (t2g)^4. The t2g set is orbitally degenerate with one hole. Elongation splits t2g, stabilizing the configuration.",
                "electron_config": "(t2g)^4(eg)^0",
                "examples": "[Cr(CN)6]4-, [Mn(CN)6]3-",
            }
        },
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize Jahn-Teller lookup tables."""
        # Octahedral (Oh) Jahn-Teller active configurations
        # Key: (d_n, spin) -> (distortion, config, explanation, examples)
        self._oh_rules = {
            # High-spin configurations
            (1, "high"): ("none", "(t2g)^1(eg)^0", "Single electron in triply-degenerate t2g; first-order JT inactive (odd electron in T term can have Kramers degeneracy).", "Ti(III), V(IV) complexes"),
            (2, "high"): ("none", "(t2g)^2(eg)^0", "Two electrons in t2g following Hund's rule; ground term ^3T1g is orbitally degenerate but often dynamically JT active.", "V(III), Ti(II)"),
            (3, "high"): ("none", "(t2g)^3(eg)^0", "Half-filled t2g (spherically symmetric); no orbital degeneracy.", "V(II), Cr(III) - no JT distortion"),
            (4, "high"): ("weak_elongation", "(t2g)^3(eg)^1", "High-spin d^4: eg set has one electron (singly degenerate). Weak JT effect possible.", "Cr(II), Mn(III) - e.g., [Cr(H2O)6]2+, [MnF6]3-"),
            (5, "high"): ("none", "(t2g)^3(eg)^2", "High-spin d^5: half-filled t2g + half-filled eg; spherically symmetric (^6A1g). No JT.", "Mn(II), Fe(III) - no JT distortion"),
            (6, "high"): ("weak_elongation", "(t2g)^4(eg)^2", "High-spin d^6: t2g has one hole (orbitally degenerate). Weak JT effect.", "Fe(II), Co(III) HS - e.g., [Fe(H2O)6]2+"),
            (7, "high"): ("none", "(t2g)^5(eg)^2", "High-spin d^7: t2g has one missing electron from half-fill; one hole but often weak/no static JT.", "Co(II) HS - weak or dynamic JT only"),
            (8, "high"): ("none", "(t2g)^6(eg)^2", "d^8: filled t2g + 2 in eg; typically square planar rather than JT-distorted octahedral.", "Ni(II) - prefers square planar geometry"),
            (9, "high"): ("strong_elongation", "(t2g)^6(eg)^3", "d^9: eg set is singly occupied asymmetrically (dx2-y2:1, dz2:2). Strong JT elongation along z-axis.", "[Cu(H2O)6]2+, [Cu(NH3)6]2+ - classic strong JT"),
            (10, "high"): ("none", "(t2g)^6(eg)^4", "d^10: fully filled; spherically symmetric. No JT.", "Cu(I), Zn(II), Ga(III) - no JT"),

            # Low-spin configurations
            (1, "low"): ("none", "(t2g)^1(eg)^0", "Low-spin d^1: single electron in t2g; same as high-spin for d^1.", "Ti(III), V(IV)"),
            (2, "low"): ("none", "(t2g)^2(eg)^0", "Low-spin d^2: two electrons in t2g.", "V(III), Ti(II)"),
            (3, "low"): ("none", "(t2g)^3(eg)^0", "Low-spin d^3: half-filled t2g; spherically symmetric.", "Cr(III), V(II) - no JT"),
            (4, "low"): ("strong_elongation", "(t2g)^4(eg)^0", "Low-spin d^4: t2g has one hole → strongly JT active (E_g ground term).", "[Cr(CN)6]4-, [Mn(CN)6]3- - strong JT"),
            (5, "low"): ("none", "(t2g)^5(eg)^0", "Low-spin d^5: one hole in t2g from half-fill; ^2T2g term but Kramers doublet protects.", "[Fe(CN)6]3- - weak or no static JT"),
            (6, "low"): ("none", "(t2g)^6(eg)^0", "Low-spin d^6: t2g fully filled; ^1A1g ground state. No JT.", "[Co(NH3)6]3+, [Fe(CN)6]4- - no JT"),
            (7, "low"): ("weak_elongation", "(t2g)^6(eg)^1", "Low-spin d^7: one electron in eg set; JT active.", "[Co(NH3)6]2+, [Ir(NH3)6]2+ - JT active"),
            (8, "low"): ("none", "(t2g)^6(eg)^2", "Low-spin d^8: usually adopts square planar.", "[Pt(NH3)4]2+, Pd(II), Au(III) - square planar"),
            (9, "low"): ("strong_elongation", "(t2g)^6(eg)^3", "Low-spin d^9: same as high-spin d^9.", "[Cu(CN)6]5- - strong JT"),
            (10, "low"): ("none", "(t2g)^6(eg)^4", "d^10: fully filled.", "No JT"),
        }

        # Tetrahedral (Td) Jahn-Teller rules (weaker than Oh due to smaller Δ)
        self._td_rules = {
            (1, "high"): ("none", "(e)^1(t2)^0", "d^1 Td: e orbital non-degenerate occupancy.", "Weak or no JT in Td"),
            (2, "high"): ("none", "(e)^2(t2)^0", "d^2 Td: e set filled.", "No JT"),
            (3, "high"): ("weak_elongation", "(e)^2(t2)^1", "d^3 Td: t2 set singly occupied; JT possible but weak (Δ_Td small).", "[CoX4]2- complexes - weak JT"),
            (4, "high"): ("none", "(e)^2(t2)^2", "d^4 Td.", "Usually no significant JT"),
            (5, "high"): ("none", "(e)^2(t2)^3", "d^5 Td: half-filled t2.", "No JT (half-filled)"),
            (6, "high"): ("weak_elongation", "(e)^3(t2)^3", "d^6 Td: one hole in e set.", "Weak JT possible"),
            (7, "high"): ("none", "(e)^4(t2)^3", "d^7 Td: e filled, t2 has one hole.", "Weak or no JT"),
            (8, "high"): ("none", "(e)^4(t2)^4", "d^8 Td.", "No JT"),
            (9, "high"): ("weak_elongation", "(e)^4(t2)^5", "d^9 Td: one hole in t2.", "Weak JT possible"),
            (10, "high"): ("none", "(e)^4(t2)^5", "d^10 Td: nearly filled.", "No JT"),
        }

    def _run_base(self, d_electron_count: int, geometry: str = "octahedral", spin_state: str = "unknown") -> dict:
        """Core logic: predict Jahn-Teller distortion."""
        if not isinstance(d_electron_count, int) or d_electron_count < 0 or d_electron_count > 10:
            raise ChemMCPError("d_electron_count must be an integer between 0 and 10.")

        geometry = geometry.lower().strip()
        if geometry not in ("octahedral", "tetrahedral"):
            raise ChemMCPError("Geometry must be 'octahedral' or 'tetrahedral'.")

        # Normalize: accept "high_spin"/"low_spin" as well as "high"/"low"
        raw_spin = spin_state.lower().strip()
        spin_map = {"high_spin": "high", "low_spin": "low", "high": "high", "low": "low", "unknown": "unknown"}
        spin_state = spin_map.get(raw_spin, raw_spin)

        if geometry == "octahedral":
            if spin_state == "unknown":
                # Return both spin states
                results = {}
                for ss in ["high", "low"]:
                    key = (d_electron_count, ss)
                    if key in self._oh_rules:
                        rule = self._oh_rules[key]
                        results[f"{ss}_spin"] = {
                            "distortion_type": rule[0],
                            "explanation": rule[2],
                            "electron_config": rule[1],
                            "examples": rule[3],
                        }
                    else:
                        results[f"{ss}_spin"] = {"distortion_type": "unknown", "explanation": f"No rule for d^{d_electron_count} {ss} spin."}
                return results
            else:
                key = (d_electron_count, spin_state)
                if key not in self._oh_rules:
                    raise ChemMCPError(f"No Jahn-Teller rule for d^{d_electron_count} ({spin_state}_spin) in {geometry} geometry.")
                rule = self._oh_rules[key]
                return {
                    "distortion_type": rule[0],
                    "explanation": rule[2],
                    "electron_config": rule[1],
                    "examples": rule[3],
                }
        else:  # tetrahedral
            if spin_state == "unknown":
                spin_state = "high"  # tetrahedral almost always high-spin
            key = (d_electron_count, spin_state)
            if key not in self._td_rules:
                raise ChemMCPError(f"No Jahn-Teller rule for d^{d_electron_count} in {geometry} geometry.")
            rule = self._td_rules[key]
            return {
                "distortion_type": rule[0],
                "explanation": rule[2],
                "electron_config": rule[1],
                "examples": rule[3],
            }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input and call core logic."""
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need at least d_electron_count and geometry.")

            d_n = int(parts[0])
            geo = parts[1]
            spin = parts[2] if len(parts) > 2 else "unknown"

            return self._run_base(d_n, geo, spin)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'd_electron_count geometry [spin_state]'")
