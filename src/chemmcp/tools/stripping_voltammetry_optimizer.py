import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class StrippingVoltammetryOptimizer(BaseTool):
    """
    溶出伏安法条件优化工具。
    根据分析物特性、浓度范围、电极类型，推荐最优的溶出伏安法实验参数。
    """
    __version__ = "0.1.0"
    name = "StrippingVoltammetryOptimizer"
    func_name = "optimize_stripping_voltammetry"
    description = "Optimize stripping voltammetry parameters: deposition potential, time, scan rate, stirring rate, etc."
    implementation_description = "Uses rule-based heuristics and electrochemical principles to recommend optimal ASV/CSV/SqW parameters based on analyte properties, concentration range, and electrode type."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Electrochemistry", "Stripping Voltammetry", "Analytical Chemistry", "Optimization"]
    required_envs = []

    code_input_sig = [
        ("analyte", "str", "N/A", "Analyte name or symbol (e.g., 'Pb', 'Cd', 'Cu', 'Zn')."),
        ("technique", "str", "ASV", "Technique: 'ASV' (anodic), 'CSV' (cathodic), 'SqW' (square wave), 'DPV' (differential pulse)."),
        ("electrode_type", "str", "HMDE", "Working electrode: 'HMDE' (hanging mercury drop), 'TFME' (thin film mercury), 'BiFE' (bismuth film), 'CPE' (carbon paste), 'GC' (glassy carbon)."),
        ("concentration_range", "str", "1e-9 to 1e-6", "Expected concentration range in mol/L (e.g., '1e-9 to 1e-6')."),
        ("matrix", "str", "water", "Sample matrix: 'water', 'urine', 'blood', 'soil_extract', 'food', 'industrial'."),
    ]

    text_input_sig = [
        ("params_str", "str", "N/A", "JSON string with all parameters, or space-separated: analyte [technique] [electrode_type] [concentration_range] [matrix]."),
    ]

    output_sig = [
        ("optimized_params", "dict", "Optimized parameter set."),
        ("reasoning", "list", "List of reasoning steps for each parameter recommendation."),
        ("warnings", "list", "Any warnings or caveats."),
        ("technique_notes", "str", "Notes on the selected technique."),
    ]

    examples = [
        {
            "code_input": {
                "analyte": "Pb",
                "technique": "SqW",
                "electrode_type": "BiFE",
                "concentration_range": "1e-9 to 1e-7",
                "matrix": "water",
            },
            "text_input": {"input_params": "Pb SqW BiFE 1e-9 to 1e-7 water"},
            "output": {
                "optimized_params": {"deposition_potential_V": -1.2, "deposition_time_s": 120,
                                   "scan_rate_V_s": 0.02, "stirring_rate_rpm": 600,
                                   "equilibration_time_s": 10, "pulse_amplitude_V": 0.025,
                                   "frequency_Hz": 25, "step_potential_V": 0.004},
            },
        },
    ]

    # --- Analyte database: standard potentials and optimal conditions ---
    _ANALYTE_DB: Dict[str, Dict] = {
        "Pb": {"E_half": -0.40, "atomic_mass": 207.2, "charge": 2, "category": "heavy_metal"},
        "Cd": {"E_half": -0.60, "atomic_mass": 112.4, "charge": 2, "category": "heavy_metal"},
        "Cu": {"E_half": +0.05, "atomic_mass": 63.55, "charge": 2, "category": "heavy_metal"},
        "Zn": {"E_half": -1.00, "atomic_mass": 65.38, "charge": 2, "category": "heavy_metal"},
        "Tl": {"E_half": -0.48, "atomic_mass": 204.38, "charge": 1, "category": "heavy_metal"},
        "In": {"E_half": -0.56, "atomic_mass": 114.82, "charge": 3, "category": "heavy_metal"},
        "Bi": {"E_half": -0.10, "atomic_mass": 208.98, "charge": 3, "category": "heavy_metal"},
        "Sb": {"E_half": -0.15, "atomic_mass": 121.76, "charge": 3, "category": "heavy_metal"},
        "Hg": {"E_half": +0.45, "atomic_mass": 200.59, "charge": 2, "category": "heavy_metal"},
        "Ag": {"E_half": +0.40, "atomic_mass": 107.87, "charge": 1, "category": "heavy_metal"},
        "Ni": {"E_half": -0.85, "atomic_mass": 58.69, "charge": 2, "category": "transition_metal"},
        "Co": {"E_half": -0.78, "atomic_mass": 58.93, "charge": 2, "category": "transition_metal"},
        "As": {"E_half": -0.50, "atomic_mass": 74.92, "charge": 3, "category": "metalloid"},
        "Se": {"E_half": -0.06, "atomic_mass": 78.97, "charge": 2, "category": "non_metal"},
        "CrVI": {"E_half": +0.35, "atomic_mass": 52.00, "charge": 3, "category": "transition_metal"},
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _parse_concentration(self, conc_range: str) -> tuple:
        """Parse concentration range string into (low, high) in mol/L."""
        try:
            parts = conc_range.lower().replace("to", "-").split("-")
            low = float(parts[0].strip())
            high = float(parts[1].strip()) if len(parts) > 1 else low * 10
            return low, high
        except Exception:
            return 1e-9, 1e-6

    def _run_base(
        self,
        analyte: str = "Pb",
        technique: str = "ASV",
        electrode_type: str = "HMDE",
        concentration_range: str = "1e-9 to 1e-6",
        matrix: str = "water",
    ) -> dict:
        """Core logic: optimize SV parameters based on inputs."""
        warnings_list: List[str] = []
        reasoning: List[str] = []

        # Normalize analyte
        an = analyte.strip().capitalize()
        if an not in self._ANALYTE_DB:
            # Try case-insensitive
            for key in self._ANALYTE_DB:
                if key.lower() == an.lower():
                    an = key
                    break
            else:
                warnings_list.append(f"Analyte '{analyte}' not in database; using generic heavy metal defaults.")
                an = "Pb"  # fallback

        info = self._ANALYTE_DB[an]
        E_half = info["E_half"]
        charge = info["charge"]

        # Parse concentration
        conc_low, conc_high = self._parse_concentration(concentration_range)
        log_conc_mid = (abs(__import__("math").log10(conc_low)) + abs(__import__("math").log10(conc_high))) / 2

        # --- Deposition Potential ---
        # More negative than E_half by ~0.3-0.5 V for complete deposition
        dep_pot = round(E_half - 0.4, 2)
        if dep_pot > 0.1:
            dep_pot = -0.1  # minimum negative for cathodic deposition
        reasoning.append(f"Deposition potential {dep_pot} V: set ~0.4 V more negative than E½={E_half} V for {an}.")

        # --- Deposition Time ---
        # Longer time for lower concentrations
        if log_conc_mid > 8:
            dep_time = 300  # very dilute: 5 min
        elif log_conc_mid > 7:
            dep_time = 180  # dilute: 3 min
        elif log_conc_mid > 6:
            dep_time = 120  # moderate: 2 min
        else:
            dep_time = 60   # concentrated: 1 min
        reasoning.append(f"Deposition time {dep_time} s: adjusted for concentration range (~10^{9-log_conc_mid:.0f} M).")

        # --- Scan Rate ---
        tech_upper = technique.upper()
        if tech_upper in ("SQW", "DPV"):
            scan_rate = 0.02  # typical for pulse techniques
        elif tech_upper == "ASV":
            scan_rate = 0.05  # linear sweep
        else:
            scan_rate = 0.01
        reasoning.append(f"Scan rate {scan_rate} V/s: optimized for {technique} technique.")

        # --- Stirring Rate ---
        elec = electrode_type.upper()
        if elec in ("HMDE", "TFME"):
            stir_rate = 600
        elif elec == "BiFE":
            stir_rate = 500
        elif elec in ("CPE", "GC"):
            stir_rate = 400
        else:
            stir_rate = 500
        reasoning.append(f"Stirring rate {stir_rate} rpm: matched to {electrode_type} electrode hydrodynamics.")

        # --- Equilibration Time ---
        eq_time = max(5, int(dep_time / 12))
        reasoning.append(f"Equilibration time {eq_time} s: allow solution to become quiescent after stirring stops.")

        # --- Pulse / Waveform Parameters (for SqW/DPV) ---
        params: Dict[str, Any] = {
            "deposition_potential_V": dep_pot,
            "deposition_time_s": dep_time,
            "scan_rate_V_s": scan_rate,
            "stirring_rate_rpm": stir_rate,
            "equilibration_time_s": eq_time,
        }

        if tech_upper == "SQW":
            params.update({
                "pulse_amplitude_V": 0.025,
                "frequency_Hz": 25,
                "step_potential_V": 0.004,
            })
            reasoning.append("Square wave: amplitude=25 mV, f=25 Hz, step=4 mV (standard compromise).")
        elif tech_upper == "DPV":
            params.update({
                "pulse_amplitude_V": 0.050,
                "pulse_width_ms": 50,
                "step_potential_V": 0.004,
            })
            reasoning.append("Differential pulse: amplitude=50 mV, width=50 ms, step=4 mV.")
        else:
            params.update({
                "pulse_amplitude_V": None,
                "frequency_Hz": None,
                "step_potential_V": None,
            })

        # Matrix-specific adjustments
        matrix_lower = matrix.lower()
        if matrix_lower in ("urine", "blood"):
            warnings_list.append("Complex biological matrix: consider standard addition method and sample digestion.")
            params["deposition_time_s"] = min(dep_time + 60, 600)
        elif matrix_lower in ("soil_extract", "food", "industrial"):
            warnings_list.append("Complex matrix: may need pH adjustment and masking agents.")
            if matrix_lower == "soil_extract":
                warnings_list.append("Soil extracts often contain surfactants that can interfere; consider UV digestion.")

        # Electrode-specific notes
        if elec == "BiFE":
            params["note"] = "Bismuth film electrode: environmentally friendly alternative to Hg; prepare fresh film before each batch."
        elif elec == "HMDE":
            warnings_list.append("Mercury electrode: follow proper disposal protocols for Hg waste.")

        # Technique notes
        tech_notes_map = {
            "ASV": "Anodic Stripping Voltammetry: oxidize deposited metal from electrode during scan.",
            "CSV": "Cathodic Stripping Voltammetry: reduce deposited film/compound during scan.",
            "SqW": "Square Wave Voltammetry: fast scan, excellent sensitivity, good for trace analysis.",
            "DPV": "Differential Pulse Voltammetry: good resolution of overlapping peaks.",
        }
        tech_note = tech_notes_map.get(tech_upper, f"{technique}: stripping voltammetry technique.")

        logger.info(f"Optimized SV params for {an}/{technique}/{electrode_type}: dep_pot={dep_pot}V, dep_time={dep_time}s")
        return {
            "optimized_params": params,
            "reasoning": reasoning,
            "warnings": warnings_list,
            "technique_notes": tech_note,
        }

    def _run_text(self, params_str: str) -> dict:
        """Parse text/JSON input."""
        import json
        params_str = params_str.strip()
        if params_str.startswith("{"):
            try:
                kwargs = json.loads(params_str)
            except json.JSONDecodeError:
                raise ChemMCPError(f"Invalid JSON input: {params_str}")
        else:
            parts = params_str.split()
            kwargs = {"analyte": parts[0] if len(parts) > 0 else "Pb"}
            if len(parts) > 1: kwargs["technique"] = parts[1]
            if len(parts) > 2: kwargs["electrode_type"] = parts[2]
            if len(parts) > 3: kwargs["concentration_range"] = parts[3]
            if len(parts) > 4: kwargs["matrix"] = parts[4]
        return self._run_base(**kwargs)
