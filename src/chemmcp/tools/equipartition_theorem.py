import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class EquipartitionTheorem(BaseTool):
    """
    能量均分定理验证工具。
    经典能量均分定理：每个二次自由度贡献 (1/2)kT 的平均能量。
    对于 n 摩尔物质：U = (f/2) * nRT，Cv = (f/2) * R，Cp = Cv + R = ((f+2)/2) * R
    """
    __version__      = "0.1.0"
    name             = "EquipartitionTheorem"
    func_name        = "equipartition_theorem"
    description      = "Verify and apply the equipartition theorem for energy distribution among molecular degrees of freedom."
    implementation_description = "Uses classical equipartition theorem: each quadratic degree of freedom contributes (1/2)RT per mole to internal energy and (1/2)R to Cv. Supports monatomic, diatomic (rigid/non-rigid), polyatomic gases and solids."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Equipartition Theorem", "Statistical Mechanics", "Thermodynamics", "Heat Capacity", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("substance_type", "str", "N/A", "Type of substance: 'monatomic', 'diatomic_rigid', 'diatomic_nonrigid', 'polyatomic_linear', 'polyatomic_nonlinear', 'solid'."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin (K). Default: 298.15."),
        ("n_moles", "float", "1.0", "Amount of substance in moles. Default: 1.0."),
        ("custom_dof", "int", "-1", "Custom total degrees of freedom (overrides substance_type if >= 0). Default: -1 (use substance_type)."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated string: 'substance_type [temperature_k] [n_moles] [custom_dof]', e.g., 'diatomic_nonrigid 300' or 'monatomic 298.15 1'."),
    ]

    output_sig       = [
        ("substance_type", "str", "Type of substance used in calculation."),
        ("degrees_of_freedom", "int", "Total number of quadratic degrees of freedom."),
        ("internal_energy_J", "float", "Total internal energy U in Joules."),
        ("cv_molar", "float", "Molar heat capacity at constant volume Cv,m in J/(mol·K)."),
        ("cp_molar", "float", "Molar heat capacity at constant pressure Cp,m in J/(mol·K)."),
        ("gamma", "float", "Heat capacity ratio γ = Cp/Cv."),
        ("energy_per_dof_J_per_mol", "float", "Energy per degree of freedom per mole in J/(mol)."),
        ("explanation", "str", "Detailed breakdown of degrees of freedom contributions."),
    ]

    examples         = [
        {
            "code_input": {
                "substance_type": "monatomic",
                "temperature_k": 300.0,
                "n_moles": 1.0,
                "custom_dof": -1,
            },
            "text_input": {
                "input_params": "monatomic 300"
            },
            "output": {
                "substance_type": "monatomic",
                "degrees_of_freedom": 3,
                "internal_energy_J": 3741.455,
                "cv_molar": 12.472,
                "cp_molar": 20.785,
                "gamma": 1.667,
                "energy_per_dof_J_per_mol": 1247.152,
                "explanation": "Monatomic ideal gas: 3 translational DOF. U = (3/2)RT, Cv = (3/2)R, Cp = (5/2)R, γ = 5/3.",
            }
        },
        {
            "code_input": {
                "substance_type": "diatomic_nonrigid",
                "temperature_k": 298.15,
                "n_moles": 1.0,
                "custom_dof": -1,
            },
            "text_input": {
                "input_params": "diatomic_nonrigid 298.15"
            },
            "output": {
                "substance_type": "diatomic_nonrigid",
                "degrees_of_freedom": 7,
                "internal_energy_J": 8730.395,
                "cv_molar": 29.099,
                "cp_molar": 37.412,
                "gamma": 1.286,
                "energy_per_dof_J_per_mol": 1247.152,
                "explanation": "Non-rigid diatomic: 3 trans + 2 rot + 1 vib(2 quad) = 7 DOF. U = (7/2)RT, Cv = (7/2)R, Cp = (9/2)R, γ = 9/7.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618  # J/(mol·K), universal gas constant

    def _get_dof(self, substance_type: str) -> int:
        """Return total quadratic degrees of freedom for given substance type."""
        dof_map = {
            "monatomic": 3,              # 3 translational
            "diatomic_rigid": 5,         # 3 trans + 2 rotational
            "diatomic_nonrigid": 7,      # 3 trans + 2 rot + 1 vib (2 quadratic)
            "polyatomic_linear": 7,      # 3 trans + 2 rot + (3N-5)vib*... simplified as 7 for typical
            "polyatomic_nonlinear": 8,   # 3 trans + 3 rot + (3N-6)vib*... simplified as 8
            "solid": 6,                  # 3 kinetic + 3 potential (harmonic)
        }
        st = substance_type.lower().replace("-", "_").replace(" ", "_")
        if st not in dof_map:
            valid = ", ".join(dof_map.keys())
            raise ChemMCPError(f"Unknown substance_type '{substance_type}'. Valid types: {valid}")
        return dof_map[st]

    def _run_base(self, substance_type: str, temperature_k: float = 298.15, n_moles: float = 1.0, custom_dof: int = -1) -> dict:
        """Apply equipartition theorem to compute thermodynamic quantities."""
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive (in Kelvin).")
        if n_moles <= 0:
            raise ChemMCPError("Number of moles must be positive.")

        if custom_dof >= 0:
            f = custom_dof
            st_label = f"custom (dof={custom_dof})"
        else:
            f = self._get_dof(substance_type)
            st_label = substance_type

        R = self.R
        T = temperature_k
        n = n_moles

        # Core equipartition results
        U = (f / 2.0) * n * R * T          # Internal energy
        Cv_mol = (f / 2.0) * R              # Molar Cv
        Cp_mol = Cv_mol + R                 # Molar Cv + R
        gamma = Cp_mol / Cv_mol if Cv_mol > 0 else float('inf')
        energy_per_dof = (1.0 / 2.0) * R * T  # RT/2 per mole per DOF

        # Build explanation
        explanations = {
            "monatomic": "Monatomic ideal gas: 3 translational DOF. U = (3/2)RT, Cv = (3/2)R, Cp = (5/2)R, γ = 5/3.",
            "diatomic_rigid": "Rigid diatomic gas: 3 translational + 2 rotational DOF. U = (5/2)RT, Cv = (5/2)R, Cp = (7/2)R, γ = 7/5.",
            "diatomic_nonrigid": "Non-rigid diatomic: 3 trans + 2 rot + 1 vibrational (2 quadratic terms) = 7 DOF. U = (7/2)RT, Cv = (7/2)R, Cp = (9/2)R, γ = 9/7.",
            "polyatomic_linear": "Linear polyatomic: 3 trans + 2 rot + vibrational DOF. Simplified model: ~7 effective DOF at moderate T.",
            "polyatomic_nonlinear": "Non-linear polyatomic: 3 trans + 3 rot + vibrational DOF. Simplified model: ~8 effective DOF at moderate T.",
            "solid": "Solid (classical): 3 kinetic + 3 potential (harmonic oscillator) = 6 DOF. U = 3nRT, Cv = 3R (Dulong-Petit law).",
        }
        explanation = explanations.get(st_label, f"Custom: {f} quadratic degrees of freedom. U = ({f}/2)nRT.")

        return {
            "substance_type": st_label,
            "degrees_of_freedom": f,
            "internal_energy_J": round(U, 3),
            "cv_molar": round(Cv_mol, 3),
            "cp_molar": round(Cp_mol, 3),
            "gamma": round(gamma, 3),
            "energy_per_dof_J_per_mol": round(energy_per_dof, 3),
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse space-separated text input."""
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Text input requires at least substance_type. Format: 'substance_type [T] [n_moles] [custom_dof]'")

        substance_type = parts[0]
        T = float(parts[1]) if len(parts) > 1 else 298.15
        n = float(parts[2]) if len(parts) > 2 else 1.0
        dof = int(parts[3]) if len(parts) > 3 else -1

        return self._run_base(substance_type, T, n, dof)
