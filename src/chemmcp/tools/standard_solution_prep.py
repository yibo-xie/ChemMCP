import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common molar masses (g/mol)
MOLAR_MASSES = {
    "NaCl": 58.44, "KCl": 74.55, "NaOH": 40.00, "KOH": 56.11,
    "Na2CO3": 105.99, "NaHCO3": 84.01, "CaCO3": 100.09,
    "HCl": 36.46, "H2SO4": 98.08, "HNO3": 63.01, "CH3COOH": 60.05,
    "C6H12O6": 180.16, "C12H22O11": 342.30, "Na2S2O3": 158.11,
    "EDTA": 292.24, "AgNO3": 169.87, "CuSO4": 159.61,
    "KMnO4": 158.03, "K2Cr2O7": 294.18, "FeSO4": 151.91,
    "Pb(NO3)2": 331.21, "ZnSO4": 161.47, "MgSO4": 120.37,
}


@ChemMCPManager.register_tool
class StandardSolutionPrep(BaseTool):
    """
    标准溶液配制指导：包括称量、定容步骤。
    支持常见化学品的摩尔质量自动查询，计算所需称量质量。
    """
    __version__ = "0.1.0"
    name = "StandardSolutionPrep"
    func_name = "prepare_standard_solution"
    description = "Provide guidance for standard solution preparation including weighing and volume fixation steps."
    implementation_description = "Calculates required mass using m = C * V * Mw formula with built-in molar mass database for common chemicals."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Sample Preparation", "Standard Solution", "Analytical Chemistry", "Laboratory"]
    required_envs = []

    code_input_sig = [
        ("solute", "str", "N/A", "Chemical formula or name of the solute (e.g., 'NaCl', 'K2Cr2O7')."),
        ("target_concentration", "float", "N/A", "Target concentration (in mol/L for molar, or g/L for mass concentration)."),
        ("target_volume_ml", "float", "N/A", "Target final volume in milliliters (mL)."),
        ("concentration_unit", "str", "mol/L", "Concentration unit: 'mol/L' or 'g/L'."),
        ("molar_mass", "float", "N/A", "Molar mass in g/mol (optional; if not provided, will look up from database)."),
        ("purity_percent", "float", "100.0", "Purity of the reagent as a percentage (default 100%)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'solute target_concentration target_volume_ml [concentration_unit] [molar_mass] [purity_percent]'."),
    ]

    output_sig = [
        ("solute", "str", "The chemical formula/name of the solute."),
        ("molar_mass_g_mol", "float", "Molar mass used (g/mol)."),
        ("mass_required_g", "float", "Mass of solute to weigh (g), corrected for purity."),
        ("target_volume_ml", "float", "Target final volume (mL)."),
        ("preparation_steps", "list", "Step-by-step preparation instructions."),
        ("notes", "list", "Important notes and precautions."),
    ]

    examples = [
        {
            "code_input": {
                "solute": "K2Cr2O7",
                "target_concentration": 0.02,
                "target_volume_ml": 250.0,
                "concentration_unit": "mol/L",
            },
            "text_input": {
                "input_params": "K2Cr2O7 0.02 250.0 mol/L",
            },
            "output": {
                "solute": "K2Cr2O7",
                "molar_mass_g_mol": 294.18,
                "mass_required_g": 1.4709,
                "target_volume_ml": 250.0,
                "preparation_steps": [
                    "Weigh 1.471 g of K2Cr2O7 (purity corrected) using an analytical balance.",
                    "Transfer quantitatively to a clean 250 mL volumetric flask.",
                    "Dissolve in approximately 150-200 mL of deionized water.",
                    "Swirl to complete dissolution.",
                    "Carefully add deionized water to the calibration mark.",
                    "Stopper and invert the flask at least 10 times to ensure homogeneity.",
                    "Label the flask with: K2Cr2O7, 0.02 mol/L, date prepared, preparer's name.",
                ],
                "notes": [
                    "K2Cr2O7 is toxic and a strong oxidizer — wear gloves and eye protection.",
                    "Can be dried at 140-150°C for 2h before weighing if primary standard grade is needed.",
                    "Store in an amber glass bottle away from light and organic materials.",
                ],
            },
        },
        {
            "code_input": {
                "solute": "NaCl",
                "target_concentration": 0.1,
                "target_volume_ml": 100.0,
                "concentration_unit": "mol/L",
                "purity_percent": 98.5,
            },
            "text_input": {
                "input_params": "NaCl 0.1 100.0 mol/L 98.5",
            },
            "output": {
                "solute": "NaCl",
                "molar_mass_g_mol": 58.44,
                "mass_required_g": 0.5933,
                "target_volume_ml": 100.0,
                "preparation_steps": [
                    "Weigh 0.593 g of NaCl (purity corrected) using an analytical balance.",
                    "Transfer quantitatively to a clean 100 mL volumetric flask.",
                    "Dissolve in approximately 60-80 mL of deionized water.",
                    "Swirl to complete dissolution.",
                    "Carefully add deionized water to the calibration mark.",
                    "Stopper and invert the flask at least 10 times to ensure homogeneity.",
                    "Label the flask with: NaCl, 0.1 mol/L, date prepared, preparer's name.",
                ],
                "notes": [
                    "Purity correction applied: raw mass adjusted from 0.584 g to 0.593 g (98.5% purity).",
                    "Store in a tightly sealed container to prevent moisture absorption.",
                ],
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        solute: str,
        target_concentration: float,
        target_volume_ml: float,
        concentration_unit: str = "mol/L",
        molar_mass: float = None,
        purity_percent: float = 100.0,
    ) -> dict:
        """Core logic: calculate weighing mass and generate preparation steps."""
        if target_concentration <= 0:
            raise ChemMCPError("Target concentration must be positive.")
        if target_volume_ml <= 0:
            raise ChemMCPError("Target volume must be positive.")
        if purity_percent <= 0 or purity_percent > 100:
            raise ChemMCPError("Purity percent must be between 0 and 100.")

        # Look up molar mass if not provided
        if molar_mass is None:
            # Try exact match first
            if solute in MOLAR_MASSES:
                molar_mass = MOLAR_MASSES[solute]
            else:
                raise ChemMCPError(
                    f"Molar mass for '{solute}' not found in database. "
                    f"Please provide molar_mass explicitly. Available: {', '.join(sorted(MOLAR_MASSES.keys()))}"
                )

        V_L = target_volume_ml / 1000.0  # Convert mL to L

        if concentration_unit == "mol/L":
            # m = C * V * Mw
            mass_theoretical = target_concentration * V_L * molar_mass
        elif concentration_unit == "g/L":
            mass_theoretical = target_concentration * V_L
        else:
            raise ChemMCPError(f"Unsupported concentration unit: {concentration_unit}. Use 'mol/L' or 'g/L'.")

        # Purity correction
        purity_factor = purity_percent / 100.0
        mass_required = mass_theoretical / purity_factor

        # Generate preparation steps
        steps = self._generate_steps(solute, mass_required, target_volume_ml, target_concentration, concentration_unit)
        notes = self._generate_notes(solute, purity_percent)

        logger.info(f"Standard solution prep: {solute}, {target_concentration} {concentration_unit}, {target_volume_ml}mL -> mass={mass_required:.4f}g")
        return {
            "solute": solute,
            "molar_mass_g_mol": round(molar_mass, 2),
            "mass_required_g": round(mass_required, 4),
            "target_volume_ml": target_volume_ml,
            "preparation_steps": steps,
            "notes": notes,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            solute = parts[0]
            conc = float(parts[1])
            vol = float(parts[2])
            unit = parts[3] if len(parts) > 3 else "mol/L"
            mw = float(parts[4]) if len(parts) > 4 else None
            purity = float(parts[5]) if len(parts) > 5 else 100.0
            return self._run_base(solute, conc, vol, unit, mw, purity)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'solute conc vol [unit] [molar_mass] [purity]'")

    def _generate_steps(self, solute, mass, vol_ml, conc, unit) -> list:
        mass_str = f"{mass:.3f}" if mass >= 0.01 else f"{mass:.4f}"
        steps = [
            f"Weigh {mass_str} g of {solute} (purity corrected) using an analytical balance.",
            f"Transfer quantitatively to a clean {int(vol_ml)} mL volumetric flask.",
            f"Dissolve in approximately {max(50, int(vol_ml * 0.5))}-{max(80, int(vol_ml * 0.7))} mL of deionized water.",
            "Swirl to complete dissolution.",
            "Carefully add deionized water to the calibration mark.",
            "Stopper and invert the flask at least 10 times to ensure homogeneity.",
            f"Label the flask with: {solute}, {conc} {unit}, date prepared, preparer's name.",
        ]
        return steps

    def _generate_notes(self, solute, purity) -> list:
        notes = []
        if purity < 100.0:
            notes.append(f"Purity correction applied: mass adjusted for {purity}% purity.")
        # Add chemical-specific warnings
        hazardous = {
            "K2Cr2O7": "Toxic strong oxidizer — wear gloves and eye protection. Store in amber bottle.",
            "KMnO4": "Strong oxidizer — avoid contact with organic matter. Store in amber bottle.",
            "HCl": "Corrosive — use in fume hood. Wear acid-resistant gloves.",
            "H2SO4": "Corrosive — always add acid to water, never reverse. Wear PPE.",
            "HNO3": "Corrosive oxidizer — use in fume hood. Wear full PPE.",
            "NaOH": "Corrosive — dissolves exothermically. Wear gloves and eye protection.",
            "AgNO3": "Light-sensitive — store in amber bottle. Will stain skin and clothing.",
            "Pb(NO3)2": "Toxic — avoid ingestion/inhalation. Use in fume hood.",
        }
        if solute in hazardous:
            notes.append(hazardous[solute])
        notes.append("Use Class A volumetric glassware for accurate results.")
        return notes
