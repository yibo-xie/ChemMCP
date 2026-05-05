import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class DebyeHuckelActivity(BaseTool):
    """
    Debye-Hückel理论计算离子活度系数工具。
    极限定律: log₁₀(γ±) = -A|z₊z₋|√I
    扩展公式: log₁₀(γ±) = -A|z₊z₋|√I / (1 + Ba√I)
    其中 A 和 B 是与温度和溶剂有关的常数。
    I = (1/2)Σ ci·zi² (离子强度)
    """
    __version__      = "0.1.0"
    name             = "DebyeHuckelActivity"
    func_name        = "debye_huckel_activity"
    description      = "Calculate ion activity coefficients using Debye-Hückel theory (limiting law and extended equation)."
    implementation_description = "Implements both the Debye-Hückel limiting law and extended Debye-Hückel equation for mean ionic activity coefficients. Computes ionic strength from ion concentrations, then applies temperature-dependent A and B constants."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Debye-Hückel", "Activity Coefficient", "Electrolyte", "Physical Chemistry", "Ionic Strength"]
    required_envs    = []

    code_input_sig   = [
        ("ions", "list", "N/A", "List of dicts with keys 'charge' (int) and 'concentration_m' (float, molality in mol/kg). Example: [{'charge': 1, 'concentration_m': 0.1}, {'charge': -1, 'concentration_m': 0.1}]."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin. Default: 298.15 K (25°C)."),
        ("equation_type", "str", "extended", "Equation type: 'limiting' or 'extended'. Default: 'extended'."),
        ("ion_size_param_a_angstrom", "float", "4.0", "Ion size parameter a in Ångströms (for extended equation). Default: 4.0."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Semicolon-separated string: 'z1,c1;z2,c2;...[;T][;equation_type][;a]'. Example: '+1,0.1;-1,0.1;298.15;extended;4'."),
    ]

    output_sig       = [
        ("temperature_K", "float", "Temperature used (K)."),
        ("ionic_strength_m", "float", "Ionic strength I in mol/kg."),
        ("ion_details", "list", "Details of each ion (charge, concentration, contribution to I)."),
        ("A_constant", "float", "Debye-Hückel constant A (mol/kg)^(-1/2)."),
        ("B_constant", "float", "Debye-Hückel constant B (Å⁻¹·(mol/kg)^(-1/2))."),
        ("equation_type", "str", "Which equation was used."),
        ("log_gamma", "float", "log₁₀(γ±)."),
        ("gamma_pm", "float", "Mean ionic activity coefficient γ±."),
        ("activity_coefficient_valid", "bool", "Whether the result is within the valid range of the model (I < ~0.6 mol/kg)."),
        ("summary", "str", "Human-readable summary of the calculation."),
    ]

    examples         = [
        {
            "code_input": {
                "ions": [{"charge": 1, "concentration_m": 0.01}, {"charge": -1, "concentration_m": 0.01}],
                "temperature_k": 298.15,
                "equation_type": "limiting",
                "ion_size_param_a_angstrom": 4.0,
            },
            "text_input": {
                "input_params": "+1,0.01;-1,0.01;298.15;limiting"
            },
            "output": {
                "temperature_K": 298.15,
                "ionic_strength_m": 0.01,
                "A_constant": 0.509,
                "B_constant": 3.287,
                "equation_type": "limiting",
                "log_gamma": -0.051,
                "gamma_pm": 0.889,
                "activity_coefficient_valid": True,
                "ion_details": [
                    {"charge": 1, "concentration_molkg": 0.01, "contribution_to_I": 0.01},
                    {"charge": -1, "concentration_molkg": 0.01, "contribution_to_I": 0.01}
                ],
                "summary": "0.01 M 1:1 electrolyte at 25°C: gamma_pm = 0.889 (limiting law).",
            }
        },
        {
            "code_input": {
                "ions": [{"charge": 2, "concentration_m": 0.005}, {"charge": -1, "concentration_m": 0.01}],
                "temperature_k": 298.15,
                "equation_type": "extended",
                "ion_size_param_a_angstrom": 4.0,
            },
            "text_input": {
                "input_params": "+2,0.005;-1,0.01;298.15;extended;4"
            },
            "output": {
                "temperature_K": 298.15,
                "ionic_strength_m": 0.015,
                "A_constant": 0.509,
                "B_constant": 3.287,
                "equation_type": "extended",
                "log_gamma": -0.153,
                "gamma_pm": 0.703,
                "activity_coefficient_valid": True,
                "ion_details": [
                    {"charge": 2, "concentration_molkg": 0.005, "contribution_to_I": 0.02},
                    {"charge": -1, "concentration_molkg": 0.01, "contribution_to_I": 0.01}
                ],
                "summary": "0.005 M CaCl2 at 25 C: gamma_pm = 0.703 (extended Debye-Huckel).",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _get_AB_constants(self, T: float) -> tuple:
        """Calculate Debye-Hückel A and B constants for water at given T.
        At 25°C (298.15 K): A ≈ 0.509 (mol/kg)^(-1/2), B ≈ 3.287 Å⁻¹(mol/kg)^(-1/2)
        Uses approximate temperature dependence.
        """
        # Reference values at 298.15 K for aqueous solution
        if abs(T - 298.15) < 0.01:
            return 0.509, 3.287

        # Temperature-dependent approximation (water)
        # A ∝ 1/(εT)^(3/2), B ∝ 1/(εT)^(1/2)
        # Simplified linear scaling from 298.15K reference
        T_ref = 298.15
        A_ref = 0.509
        B_ref = 3.287

        ratio = (T_ref / T) ** 1.5
        A = A_ref * ratio * (298.15 / T)
        B = B_ref * math.sqrt(T_ref / T)

        return round(A, 4), round(B, 4)

    def _calc_ionic_strength(self, ions: list) -> float:
        """Calculate ionic strength I = (1/2) Σ ci·zi²."""
        I = 0.0
        for ion in ions:
            z = ion.get("charge", 0)
            c = ion.get("concentration_m", 0.0)
            I += c * z * z
        return I / 2.0

    def _run_base(self, ions: list, temperature_k: float = 298.15,
                  equation_type: str = "extended",
                  ion_size_param_a_angstrom: float = 4.0) -> dict:
        """Calculate mean ionic activity coefficient using Debye-Hückel theory."""
        if not ions or len(ions) < 2:
            raise ChemMCPError("Must provide at least two ions (cation and anion).")
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")

        etype = equation_type.lower()
        if etype not in ("limiting", "extended"):
            raise ChemMCPError(f"equation_type must be 'limiting' or 'extended', got '{equation_type}'")

        T = temperature_k
        a_ang = ion_size_param_a_angstrom

        # Ionic strength
        I = self._calc_ionic_strength(ions)

        # Ion details
        ion_details = []
        for ion in ions:
            z = ion["charge"]
            c = ion["concentration_m"]
            ion_details.append({
                "charge": z,
                "concentration_molkg": c,
                "contribution_to_I": c * z * z,
            })

        # Get |z+ * z-|
        charges = [abs(ion["charge"]) for ion in ions]
        z_pos_max = max((ion["charge"] for ion in ions if ion["charge"] > 0), default=1)
        z_neg_max = max((-ion["charge"] for ion in ions if ion["charge"] < 0), default=1)
        zz_product = abs(z_pos_max * z_neg_max)

        # A and B constants
        A_val, B_val = self._get_AB_constants(T)

        sqrt_I = math.sqrt(I) if I >= 0 else 0.0

        # Apply equation
        if etype == "limiting":
            log_gamma = -A_val * zz_product * sqrt_I
        else:  # extended
            denom = 1.0 + B_val * a_ang * sqrt_I
            log_gamma = -A_val * zz_product * sqrt_I / denom

        gamma_pm = 10.0 ** log_gamma

        # Validity check (model reliable for I < 0.6 mol/kg approximately)
        valid = I < 0.6

        summary = (
            f"Debye-Hückel ({etype} equation) at {T} K:\n"
            f"Ionic strength I = {I:.4f} mol/kg\n"
            f"|z₊·z₋| = {zz_product}\n"
            f"A = {A_val:.4f}, B = {B_val:.4f}\n"
            f"log₁₀(γ±) = {log_gamma:.4f}\n"
            f"γ± = {gamma_pm:.4f}"
            + ("" if valid else " ⚠️ Model may be inaccurate at this ionic strength.")
        )

        return {
            "temperature_K": T,
            "ionic_strength_m": round(I, 6),
            "ion_details": ion_details,
            "A_constant": A_val,
            "B_constant": B_val,
            "equation_type": etype,
            "log_gamma": round(log_gamma, 4),
            "gamma_pm": round(gamma_pm, 4),
            "activity_coefficient_valid": valid,
            "summary": summary,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse semicolon-separated text input."""
        parts = input_params.strip().split(";")
        if len(parts) < 2:
            raise ChemMCPError(
                "Text input requires at least ion data. "
                "Format: 'z1,c1;z2,c2;...[;T][;eq_type][;a]'"
            )

        # Parse ions
        ions = []
        for p in parts[:len(parts)]:
            subparts = p.strip().split(",")
            if len(subparts) == 2:
                try:
                    z = int(subparts[0].strip().replace("+", ""))
                    if subparts[0].strip().startswith("-"):
                        z = -z
                    c = float(subparts[1].strip())
                    ions.append({"charge": z, "concentration_m": c})
                except ValueError:
                    continue

        if not ions:
            raise ChemMCPError(f"Could not parse ion data from '{input_params}'")

        T = float(parts[len(ions)].strip()) if len(parts) > len(ions) else 298.15
        eq_type = parts[len(ions) + 1].strip() if len(parts) > len(ions) + 1 else "extended"
        a = float(parts[len(ions) + 2].strip()) if len(parts) > len(ions) + 2 else 4.0

        return self._run_base(ions, T, eq_type, a)
