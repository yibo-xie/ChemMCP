import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CalculatePH(BaseTool):
    """
    计算溶液 pH 值。
    支持强酸、强碱、弱酸、弱碱、缓冲溶液、盐类水解等多种场景。
    """
    __version__ = "0.1.0"
    name = "CalculatePH"
    func_name = "calculate_ph"
    description = "Calculate solution pH for various scenarios: strong acid/base, weak acid/base, buffer solutions, and salt hydrolysis."
    implementation_description = "Implements standard analytical chemistry formulas for each scenario with proper approximations and exact quadratic solutions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["pH", "Acid-Base", "Solution Chemistry", "Equilibrium"]
    required_envs = []

    code_input_sig = [
        ("scenario", "str", "N/A", "Scenario type: 'strong_acid', 'strong_base', 'weak_acid', 'weak_base', 'buffer', 'salt_hydrolysis'."),
        ("concentration", "float", "N/A", "Concentration of the solute (mol/L)."),
        ("Ka", "float", "None", "Acid dissociation constant (needed for weak_acid, buffer, salt of weak_base)."),
        ("Kb", "float", "None", "Base dissociation constant (needed for weak_base, salt of weak_acid)."),
        ("salt_type", "str", "None", "For salt_hydrolysis: 'strong_acid_weak_base', 'weak_acid_strong_base', 'weak_acid_weak_base'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated string: 'scenario concentration [Ka] [Kb] [salt_type]'. Example: 'weak_acid 0.1 1.8e-5'"),
    ]

    output_sig = [
        ("ph", "float", "Calculated pH value."),
        ("h_conc", "float", "[H+] concentration in mol/L."),
        ("oh_conc", "float", "[OH-] concentration in mol/L."),
        ("scenario", "str", "The scenario used for calculation."),
        ("explanation", "str", "Step-by-step explanation of the calculation."),
    ]

    examples = [
        {
            "code_input": {
                "scenario": "weak_acid",
                "concentration": 0.1,
                "Ka": 1.8e-5,
                "Kb": None,
                "salt_type": None,
            },
            "text_input": {
                "input_params": "weak_acid 0.1 1.8e-5"
            },
            "output": {
                "ph": 2.87,
                "h_conc": 0.00134,
                "oh_conc": 7.46e-12,
                "scenario": "weak_acid",
                "explanation": "Weak acid equilibrium.",
            }
        },
        {
            "code_input": {
                "scenario": "strong_acid",
                "concentration": 0.01,
                "Ka": None,
                "Kb": None,
                "salt_type": None,
            },
            "text_input": {
                "input_params": "strong_acid 0.01"
            },
            "output": {
                "ph": 2.0,
                "h_conc": 0.01,
                "oh_conc": 1.0e-12,
                "scenario": "strong_acid",
                "explanation": "Strong acid fully dissociates.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.Kw = 1.0e-14  # 水的离子积常数 (25°C)

    def _run_base(self, scenario: str, concentration: float, Ka: float = None, Kb: float = None, salt_type: str = None) -> dict:
        """核心逻辑：根据场景计算pH"""
        if concentration <= 0:
            raise ChemMCPError("Concentration must be positive.")

        scenario = scenario.lower().strip() if scenario else ""
        C = concentration

        if scenario == "strong_acid":
            return self._calc_strong_acid(C)
        elif scenario == "strong_base":
            return self._calc_strong_base(C)
        elif scenario == "weak_acid":
            if Ka is None:
                raise ChemMCPError("Ka is required for weak_acid scenario.")
            return self._calc_weak_acid(C, Ka)
        elif scenario == "weak_base":
            if Kb is None:
                raise ChemMCPError("Kb is required for weak_base scenario.")
            return self._calc_weak_base(C, Kb)
        elif scenario == "buffer":
            if Ka is None:
                raise ChemMCPError("Ka and [A-]/[HA] ratio needed for buffer scenario.")
            # For buffer, we need ratio; use default 1:1 if not specified
            return self._calc_buffer(C, Ka)
        elif scenario == "salt_hydrolysis":
            return self._calc_salt_hydrolysis(C, Ka, Kb, salt_type)
        else:
            raise ChemMCPError(
                f"Unknown scenario '{scenario}'. Supported: strong_acid, strong_base, "
                f"weak_acid, weak_base, buffer, salt_hydrolysis"
            )

    def _calc_strong_acid(self, C: float) -> dict:
        """强酸：完全解离，[H+] = C"""
        h = C
        oh = self.Kw / h
        ph = -math.log10(h)
        return {
            "ph": round(ph, 4),
            "h_conc": h,
            "oh_conc": round(oh, 15),
            "scenario": "strong_acid",
            "explanation": f"Strong acid fully dissociates: [H+] = {C} M → pH = {-math.log10(C):.4f}",
        }

    def _calc_strong_base(self, C: float) -> dict:
        """强碱：完全解离，[OH-] = C"""
        oh = C
        h = self.Kw / oh
        ph = -math.log10(h)
        poh = -math.log10(oh)
        return {
            "ph": round(ph, 4),
            "h_conc": round(h, 15),
            "oh_conc": oh,
            "scenario": "strong_base",
            "explanation": f"Strong base fully dissociates: [OH-] = {C} M, pOH = {poh:.4f}, pH = 14 - pOH = {ph:.4f}",
        }

    def _calc_weak_acid(self, C: float, Ka: float) -> dict:
        """弱酸：HA ⇌ H+ + A-, Ka = x²/(C-x), 精确二次求解"""
        disc = Ka ** 2 + 4 * Ka * C
        h = (-Ka + math.sqrt(disc)) / 2
        ph = -math.log10(h) if h > 0 else 7.0
        oh = self.Kw / h
        alpha = h / C

        explanation = (
            f"Weak acid HA ⇌ H+ + A-. Solving Ka = x²/({C}-x) exactly: "
            f"x = [H+] = {h:.6f} M, α = {alpha:.6f}, pH = {ph:.4f}"
        )
        return {
            "ph": round(ph, 4),
            "h_conc": round(h, 10),
            "oh_conc": round(oh, 15),
            "scenario": "weak_acid",
            "explanation": explanation,
        }

    def _calc_weak_base(self, C: float, Kb: float) -> dict:
        """弱碱：B + H2O ⇌ BH+ + OH-, Kb = x²/(C-x)"""
        disc = Kb ** 2 + 4 * Kb * C
        oh = (-Kb + math.sqrt(disc)) / 2
        poh = -math.log10(oh) if oh > 0 else 7.0
        ph = 14 - poh
        h = self.Kw / oh if oh > 0 else 1e-7
        alpha = oh / C

        explanation = (
            f"Weak base B + H2O ⇌ BH+ + OH-. Solving Kb = x²/({C}-x): "
            f"[OH-] = {oh:.6f} M, pOH = {poh:.4f}, pH = {ph:.4f}"
        )
        return {
            "ph": round(ph, 4),
            "h_conc": round(h, 15),
            "oh_conc": round(oh, 10),
            "scenario": "weak_base",
            "explanation": explanation,
        }

    def _calc_buffer(self, C: float, Ka: float) -> dict:
        """缓冲溶液（简化：等浓度混合）"""
        # Henderson-Hasselbalch: pH = pKa + log([A-]/[HA])
        # 默认 1:1 缓冲对
        pKa = -math.log10(Ka)
        ph = pKa + math.log10(1.0)  # 1:1 ratio
        h = 10 ** (-ph)
        oh = self.Kw / h

        explanation = (
            f"Buffer solution (1:1 ratio). Henderson-Hasselbalch: "
            f"pH = pKa + log([A-]/[HA]) = {pKa:.2f} + 0 = {ph:.4f}"
        )
        return {
            "ph": round(ph, 4),
            "h_conc": round(h, 15),
            "oh_conc": round(oh, 15),
            "scenario": "buffer",
            "explanation": explanation,
        }

    def _calc_salt_hydrolysis(self, C: float, Ka: float, Kb: float, salt_type: str) -> dict:
        """盐类水解"""
        st = (salt_type or "").lower()

        if st in ("strong_acid_weak_base", "sawb"):
            # 强酸弱碱盐（如 NH4Cl）：NH4+ + H2O ⇌ NH3 + H3O+
            # Kh = Kw/Kb(NH3)
            if Kb is None or Kb <= 0:
                raise ChemMCPError("Kb needed for strong_acid_weak_base salt hydrolysis.")
            Kh = self.Kw / Kb
            disc = Kh ** 2 + 4 * Kh * C
            h = (-Kh + math.sqrt(disc)) / 2
            ph = -math.log10(h) if h > 0 else 7.0
            oh = self.Kw / h
            explanation = f"Salt of strong acid + weak base: Kh=Kw/Kb={Kh:.2e}, [H+]={h:.6f} M, pH={ph:.4f}"

        elif st in ("weak_acid_strong_base", "wasb"):
            # 弱酸盐（如 CH3COONa）：A- + H2O ⇌ HA + OH-
            if Ka is None or Ka <= 0:
                raise ChemMCPError("Ka needed for weak_acid_strong_base salt hydrolysis.")
            Kh = self.Kw / Ka
            disc = Kh ** 2 + 4 * Kh * C
            oh = (-Kh + math.sqrt(disc)) / 2
            poh = -math.log10(oh) if oh > 0 else 7.0
            ph = 14 - poh
            h = self.Kw / oh
            explanation = f"Salt of weak acid + strong base: Kh=Kw/Ka={Kh:.2e}, [OH-]={oh:.6f} M, pH={ph:.4f}"

        elif st in ("weak_acid_weak_base", "wawb"):
            # 弱酸弱碱盐：(如 NH4CH3COO)
            if Ka is None or Kb is None:
                raise ChemMCPError("Both Ka and Kb needed for weak_acid_weak_base salt hydrolysis.")
            h = math.sqrt(self.Kw * Ka / Kb)
            ph = -math.log10(h) if h > 0 else 7.0
            oh = self.Kw / h
            explanation = f"Salt of weak acid + weak base: [H+]=sqrt(Kw*Ka/Kb)={h:.6f} M, pH={ph:.4f}"
        else:
            raise ChemMCPError(f"Unknown salt_type '{salt_type}'. Use: strong_acid_weak_base, weak_acid_strong_base, or weak_acid_weak_base")

        return {
            "ph": round(ph, 4),
            "h_conc": round(h, 15),
            "oh_conc": round(oh, 15),
            "scenario": "salt_hydrolysis",
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need at least scenario and concentration.")
            scenario = parts[0]
            C = float(parts[1])
            Ka = float(parts[2]) if len(parts) > 2 and parts[2].lower() != "none" else None
            Kb = float(parts[3]) if len(parts) > 3 and parts[3].lower() != "none" else None
            salt_type = parts[4] if len(parts) > 4 else None
            return self._run_base(scenario, C, Ka, Kb, salt_type)
        except (ValueError, IndexError) as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'scenario conc [Ka] [Kb] [salt_type]'")
