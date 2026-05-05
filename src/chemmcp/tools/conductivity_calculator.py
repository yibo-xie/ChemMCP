import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ConductivityCalculator(BaseTool):
    """
    电导率和摩尔电导率计算工具。
    κ = Λm · c,  Kohlrausch: Λm = Λm° - K√c,  α = Λm/Λm°
    """
    __version__      = "0.1.0"
    name             = "ConductivityCalculator"
    func_name        = "conductivity_calculator"
    description      = "Calculate electrical conductivity, molar conductivity, degree of dissociation using Kohlrausch's law."
    implementation_description = "Implements Kohlrausch's law: Lambda_m = Lambda_m0 - K*sqrt(c). Computes conductivity kappa, molar conductivity Lambda_m, and dissociation alpha."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Conductivity", "Kohlrausch", "Electrolyte", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("calc_mode", "str", "N/A", "'from_lambda' (from ionic conductivities), 'from_kappa' (from measured kappa), or 'kohlrausch_fit'."),
        ("concentration_mol_m3", "float", "N/A", "Concentration in mol/m^3."),
        ("lambda_plus_sm2mol", "float", "0", "Cation molar ionic conductivity at infinite dilution (S·m²/mol)."),
        ("lambda_minus_sm2mol", "float", "0", "Anion molar ionic conductivity at infinite dilution (S·m²/mol)."),
        ("nu_plus", "int", "1", "Stoichiometric number of cation. Default: 1."),
        ("nu_minus", "int", "1", "Stoichiometric number of anion. Default: 1."),
        ("measured_kappa_sm", "float", "0", "Measured conductivity κ in S/m (for 'from_kappa' mode only)."),
        ("kohlrausch_K", "float", "0", "Kohlrausch constant K (for 'kohlrausch_fit' mode only)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated: 'mode;c;lambda_plus;lambda_minus[;nu_plus][;nu_minus][;kappa_or_K]'. Example: 'from_lambda;100;7.35e-3;7.63e-3'."),
    ]

    output_sig       = [
        ("calc_mode", "str", "Calculation mode used."),
        ("concentration_mol_m3", "float", "Concentration used (mol/m³)."),
        ("concentration_mol_L", "float", "Concentration in mol/L."),
        ("Lambda_m0_sm2mol", "float", "Limiting molar conductivity Λm° (S·m²/mol)."),
        ("Lambda_m_sm2mol", "float", "Molar conductivity Λm (S·m²/mol)."),
        ("kappa_Sm", "float", "Electrical conductivity κ (S/m)."),
        ("degree_of_dissociation", "float", "Degree of dissociation α = Λm/Λm°."),
        ("kohlrausch_K_value", "float", "Kohlrausch K constant if applicable."),
        ("summary", "str", "Human-readable summary."),
    ]

    examples         = [
        {
            "code_input": {
                "calc_mode": "from_lambda",
                "concentration_mol_m3": 100.0,
                "lambda_plus_sm2mol": 7.35e-3,
                "lambda_minus_sm2mol": 7.63e-3,
                "nu_plus": 1,
                "nu_minus": 1,
                "measured_kappa_sm": 0.0,
                "kohlrausch_K": 0.0,
            },
            "text_input": {
                "input_params": "from_lambda;100;7.35e-3;7.63e-3"
            },
            "output": {
                "calc_mode": "from_lambda",
                "concentration_mol_m3": 100.0,
                "concentration_mol_L": 0.10,
                "Lambda_m0_sm2mol": 0.01498,
                "Lambda_m_sm2mol": 0.01048,
                "kappa_Sm": 1.048,
                "degree_of_dissociation": 0.700,
                "kohlrausch_K_value": 0.00027,
                "summary": "KCl 0.10 M: Lambda_m0=0.0150, Lambda_m=0.0105, kappa=1.048 S/m, alpha=0.700.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, calc_mode: str, concentration_mol_m3: float,
                  lambda_plus_sm2mol: float = 0.0, lambda_minus_sm2mol: float = 0.0,
                  nu_plus: int = 1, nu_minus: int = 1,
                  measured_kappa_sm: float = 0.0, kohlrausch_K: float = 0.0) -> dict:
        """Calculate conductivity quantities."""
        if concentration_mol_m3 < 0:
            raise ChemMCPError("Concentration cannot be negative.")

        mode = calc_mode.lower().replace("-", "_")
        c = concentration_mol_m3
        c_mol_L = c / 1000.0

        valid_modes = {"from_lambda", "from_kappa", "kohlrausch_fit"}
        if mode not in valid_modes:
            raise ChemMCPError(f"Unknown calc_mode '{calc_mode}'. Use: {valid_modes}")

        # Limiting molar conductivity: Λm° = ν₊λ₊° + ν₋λ₋°
        Lambda_0 = nu_plus * lambda_plus_sm2mol + nu_minus * lambda_minus_sm2mol

        sqrt_c = math.sqrt(c) if c >= 0 else 0.0

        if mode == "from_lambda":
            # Kohlrausch: Λm = Λm° - K√c, with estimated K for typical electrolytes
            # For 1:1 electrolyte at 25°C, K ≈ 8.6×10⁻³ (when c in mol/L)
            # Convert to mol/m³ basis: K_m3 = K_L / sqrt(1000)
            K_est = 8.6e-3 / math.sqrt(1000) if Lambda_0 > 0 else 0
            Lambda_m = max(0, Lambda_0 - K_est * sqrt_c)
            kappa = Lambda_m * c
            kohlrausch_K_val = round(K_est, 8)

        elif mode == "from_kappa":
            kappa = measured_kappa_sm
            Lambda_m = kappa / c if c > 0 else 0
            kohlrausch_K_val = None

        elif mode == "kohlrausch_fit":
            Lambda_m = max(0, Lambda_0 - kohlrausch_K * sqrt_c)
            kappa = Lambda_m * c
            kohlrausch_K_val = kohlrausch_K
        else:
            Lambda_m = 0
            kappa = 0
            kohlrausch_K_val = None

        alpha = Lambda_m / Lambda_0 if Lambda_0 > 0 else 0

        summary = (
            f"Conductivity calculation ({mode}):\n"
            f"c = {c:.2f} mol/m³ ({c_mol_L:.4f} mol/L)\n"
            f"Λm° = {Lambda_0:.6f} S·m²/mol\n"
            f"Λm = {Lambda_m:.6f} S·m²/mol\n"
            f"κ = Λm × c = {kappa:.4f} S/m\n"
            f"α = Λm/Λm° = {alpha:.4f}"
        )

        return {
            "calc_mode": mode,
            "concentration_mol_m3": c,
            "concentration_mol_L": round(c_mol_L, 6),
            "Lambda_m0_sm2mol": round(Lambda_0, 8),
            "Lambda_m_sm2mol": round(Lambda_m, 8),
            "kappa_Sm": round(kappa, 6),
            "degree_of_dissociation": round(alpha, 4),
            "kohlrausch_K_value": kohlrausch_K_val,
            "summary": summary,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse semicolon-separated text input."""
        parts = input_params.strip().split(";")
        if len(parts) < 2:
            raise ChemMCPError(
                "Text input requires at least mode and concentration. "
                "Format: 'mode;c;lp;lm[;np][;nm][;kappa/K]'"
            )

        try:
            mode = parts[0].strip()
            c = float(parts[1].strip())
            lp = float(parts[2].strip()) if len(parts) > 2 else 0.0
            lm = float(parts[3].strip()) if len(parts) > 3 else 0.0
            np_ = int(parts[4].strip()) if len(parts) > 4 else 1
            nm_ = int(parts[5].strip()) if len(parts) > 5 else 1
            extra = float(parts[6].strip()) if len(parts) > 6 else 0.0
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse values from '{input_params}': {e}")

        kwargs = {
            "calc_mode": mode, "concentration_mol_m3": c,
            "lambda_plus_sm2mol": lp, "lambda_minus_sm2mol": lm,
            "nu_plus": np_, "nu_minus": nm_,
        }
        if mode == "from_kappa":
            kwargs["measured_kappa_sm"] = extra
        elif mode == "kohlrausch_fit":
            kwargs["kohlrausch_K"] = extra

        return self._run_base(**kwargs)
