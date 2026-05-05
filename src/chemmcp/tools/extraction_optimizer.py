import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common log P (partition coefficient) values for solvent-water systems
LOG_P_DATA = {
    # (solute, organic_solvent) -> logP(octanol/water) approximate
    "benzene": 2.13, "toluene": 2.73, "ethyl_acetate": 0.89,
    "dichloromethane": 1.25, "chloroform": 1.97, "diethyl_ether": 0.98,
    "hexane": 3.50, "heptane": 4.00, "phenol": 1.46, "aniline": 0.90,
    "benzoic_acid": 1.87, "naphthalene": 3.30, "caffeine": -0.07,
    "ibuprofen": 3.97, "paracetamol": 0.46, "aspirin": 1.19,
    "atrazine": 2.61, "DDT": 6.91, "lindane": 3.72,
}


@ChemMCPManager.register_tool
class ExtractionOptimizer(BaseTool):
    """
    液-液萃取条件优化：计算分配系数、萃取效率和最佳溶剂体积。
    支持单次和多次连续萃取优化。
    """
    __version__ = "0.1.0"
    name = "ExtractionOptimizer"
    func_name = "optimize_extraction"
    description = "Optimize liquid-liquid extraction conditions: calculate partition coefficient, extraction efficiency, and optimal solvent volume."
    implementation_description = "Uses distribution law D = C_org/C_aq and extraction efficiency formulas for single and multiple extractions."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Extraction", "Sample Preparation", "Partition Coefficient", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("solute_name", "str", "N/A", "Name of the solute to extract (e.g., 'benzene', 'benzoic_acid')."),
        ("aqueous_volume_ml", "float", "N/A", "Volume of aqueous phase (mL)."),
        ("organic_volume_ml", "float", "N/A", "Volume of organic solvent per extraction (mL)."),
        ("num_extractions", "int", "1", "Number of extraction steps."),
        ("partition_coefficient", "float", "N/A", "Partition coefficient Kd = C_org / C_aq (if known; otherwise auto-looked up)."),
        ("pKa", "float", "None", "pKa of solute (for pH-dependent extraction optimization)."),
        ("pH", "float", "7.0", "Current pH of aqueous phase."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'solute_name aqueous_volume_ml organic_volume_ml [num_extractions] [Kd] [pKa] [pH]'."),
    ]

    output_sig = [
        ("partition_coefficient_Kd", "float", "Partition coefficient used in calculation."),
        ("single_extraction_efficiency_pct", "float", "Extraction efficiency with one extraction (%)."),
        ("total_efficiency_pct", "float", "Total extraction efficiency after all extractions (%)."),
        ("remaining_fraction", "float", "Fraction remaining in aqueous phase."),
        ("optimal_strategy", "str", "Recommended extraction strategy."),
        ("details", "list", "Step-by-step extraction details."),
        ("ph_optimization", "dict", "pH optimization suggestions (if pKa provided)."),
    ]

    examples = [
        {
            "code_input": {
                "solute_name": "benzoic_acid",
                "aqueous_volume_ml": 100.0,
                "organic_volume_ml": 50.0,
                "num_extractions": 1,
            },
            "text_input": {
                "input_params": "benzoic_acid 100.0 50.0 1",
            },
            "output": {
                "partition_coefficient_Kd": 74.13,
                "single_extraction_efficiency_pct": 97.37,
                "total_efficiency_pct": 97.37,
                "remaining_fraction": 0.0263,
                "optimal_strategy": "Single extraction with V_org=50mL gives >97% recovery.",
                "details": [{"step": 1, "V_org_ml": 50.0, "efficiency_pct": 97.37, "remaining_frac": 0.0263}],
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        solute_name: str,
        aqueous_volume_ml: float,
        organic_volume_ml: float,
        num_extractions: int = 1,
        partition_coefficient: float = None,
        pKa: float = None,
        pH: float = 7.0,
    ) -> dict:
        """Core logic: optimize LLE conditions."""
        if aqueous_volume_ml <= 0 or organic_volume_ml <= 0:
            raise ChemMCPError("Volumes must be positive.")
        if num_extractions < 1:
            raise ChemMCPError("Number of extractions must be >= 1.")

        # Look up Kd if not provided
        if partition_coefficient is None:
            key = solute_name.lower().replace(" ", "_")
            if key in LOG_P_DATA:
                log_p = LOG_P_DATA[key]
                partition_coefficient = math.pow(10, log_p)
            else:
                raise ChemMCPError(
                    f"Partition coefficient for '{solute_name}' not found. "
                    f"Please provide partition_coefficient explicitly. "
                    f"Available: {', '.join(sorted(LOG_P_DATA.keys()))}"
                )

        # pH adjustment for ionizable compounds
        effective_kd = partition_coefficient
        ph_info = None
        if pKa is not None:
            # For acidic compounds: fraction neutral = 1 / (1 + 10^(pH-pKa))
            # Effective Kd = Kd * f_neutral
            if pH < 14:
                f_neutral = 1.0 / (1.0 + math.pow(10, pH - pKa))
                effective_kd = partition_coefficient * f_neutral
                ph_info = {
                    "pKa": pKa,
                    "current_pH": pH,
                    "neutral_fraction": round(f_neutral, 4),
                    "effective_Kd": round(effective_kd, 4),
                    "recommended_pH_range": f"< {pKa - 1.5:.1f} (acidic)" if pKa < 10 else f"> {pKa + 1.5:.1f} (basic)",
                    "note": "Lower pH increases neutral form for acidic compounds, improving extraction.",
                }

        Kd = effective_kd
        V_aq = aqueous_volume_ml
        V_org = organic_volume_ml

        # Calculate extraction details
        details = []
        remaining_frac = 1.0
        for i in range(num_extractions):
            # Fraction extracted in this step
            step_extracted = (Kd * V_org) / (Kd * V_org + V_aq)
            step_efficiency = step_extracted * 100.0
            remaining_frac *= (1.0 - step_extracted)
            details.append({
                "step": i + 1,
                "V_org_ml": V_org,
                "efficiency_pct": round(step_efficiency, 2),
                "remaining_frac": round(remaining_frac, 6),
            })

        total_efficiency = (1.0 - remaining_frac) * 100.0
        single_efficiency = ((Kd * V_org) / (Kd * V_org + V_aq)) * 100.0

        # Optimal strategy recommendation
        if total_efficiency >= 99.0:
            strategy = f"{num_extractions}-extraction protocol achieves excellent recovery ({total_efficiency:.1f}%)."
        elif total_efficiency >= 95.0:
            strategy = f"{num_extractions}-extraction protocol gives good recovery ({total_efficiency:.1f}%). Consider increasing to {num_extractions + 1} extractions or solvent volume for >99%."
        else:
            suggested_vol = V_aq / Kd * 50  # Target ~98% single extraction
            strategy = (
                f"Current recovery ({total_efficiency:.1f}%) may be insufficient. "
                f"Recommendations: 1) Increase organic volume to ~{suggested_vol:.0f} mL per extraction, "
                f"2) Increase to {min(num_extractions + 2, 5)} extraction steps, "
                f"3) Adjust pH if compound is ionizable."
            )

        logger.info(f"Extraction optimized for {solute_name}: Kd={Kd:.2f}, total_eff={total_efficiency:.1f}%")
        return {
            "partition_coefficient_Kd": round(Kd, 4),
            "single_extraction_efficiency_pct": round(single_efficiency, 2),
            "total_efficiency_pct": round(total_efficiency, 2),
            "remaining_fraction": round(remaining_frac, 6),
            "optimal_strategy": strategy,
            "details": details,
            "ph_optimization": ph_info,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            solute = parts[0]
            v_aq = float(parts[1])
            v_org = float(parts[2])
            n_ext = int(parts[3]) if len(parts) > 3 else 1
            kd = float(parts[4]) if len(parts) > 4 else None
            pka = float(parts[5]) if len(parts) > 5 else None
            ph = float(parts[6]) if len(parts) > 6 else 7.0
            return self._run_base(solute, v_aq, v_org, n_ext, kd, pka, ph)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'solute V_aq V_org [n_ext] [Kd] [pKa] [pH]'")
