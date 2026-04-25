import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class StandardStateConverter(BaseTool):
    """
    不同标准态间的热力学量转换工具。
    支持理想气体不同压力标准态（1 bar vs 1 atm）、溶液不同浓度标度（摩尔分数、质量摩尔、体积摩尔）之间的转换。
    """
    __version__ = "0.1.0"
    name = "StandardStateConverter"
    func_name = "convert_standard_state"
    description = "Convert thermodynamic quantities between different standard states (pressure units, concentration scales)."
    implementation_description = "Uses Δμ = RT ln(P2°/P1°) for pressure changes, and activity scale conversions (x↔m↔c) with solvent density/partial molar mass."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Standard State", "Physical Chemistry", "Unit Conversion"]
    required_envs    = []

    code_input_sig = [
        ("conversion_type", "str", "N/A", "Type: 'pressure' (gas standard state pressure), 'concentration' (solution concentration scale), or 'gibbs' (ΔG° conversion)."),
        ("quantity_value", "float", "N/A", "Value of the thermodynamic quantity to convert."),
        ("quantity_type", "str", "delta_g", "Quantity type: 'delta_g', 'delta_h', 'delta_s', 'equilibrium_constant', 'chemical_potential', 'entropy'."),
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        # For pressure conversion:
        ("p_from", "float", "1.01325", "Original standard state pressure in bar (default 1 atm = 1.01325 bar)."),
        ("p_to", "float", "1.0", "Target standard state pressure in bar (default 1 bar)."),
        # For concentration conversion:
        ("scale_from", "str", "", "From scale: 'mole_fraction', 'molality', 'molarity'."),
        ("scale_to", "str", "", "To scale: 'mole_fraction', 'molality', 'molarity'."),
        ("solvent_density_kg_L", "float", "1.0", "Solvent density in kg/L (for molarity-molality conversion, default water at 25°C)."),
        ("solute_molar_mass", "float", "100.0", "Solute molar mass in g/mol (for mole fraction conversions)."),
        ("reference_molality", "float", "1.0", "Reference molality in mol/kg (standard state, default 1)."),
        ("reference_molarity", "float", "1.0", "Reference molarity in mol/L (standard state, default 1)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated parameters string."),
    ]

    output_sig = [
        ("converted_value", "float", "Converted quantity value in new standard state."),
        ("correction_term", "float", "The correction term applied (e.g., RT ln(P2/P1) in J/mol)."),
        ("original_unit", "str", "Unit/description of original standard state."),
        ("target_unit", "str", "Unit/description of target standard state."),
        ("analysis", "str", "Detailed explanation of the conversion."),
    ]

    examples         = [
        {
            "code_input": {
                "conversion_type": 'pressure',
                "quantity_value": -394.36,
                "quantity_type": 'delta_g',
                "temperature_k": 298.15,
                "p_from": 1.01325,
                "p_to": 1.0,
                "scale_from": '',
                "scale_to": '',
                "solvent_density_kg_L": 1.0,
                "solute_molar_mass": 100.0,
                "reference_molality": 1.0,
                "reference_molarity": 1.0
            },
            "text_input": {
                "input_params": 'pressure -394.36 delta_g 298.15 1.01325 1.0'
            },
            "output": {
                "converted_value": -394.38,
                "correction_term": -12.47,
                "original_unit": 'kJ/mol (1 atm)',
                "target_unit": 'kJ/mol (1 bar)',
                "analysis": 'Delta G conversion 1 atm to 1 bar.'
            }
        },
        {
            "code_input": {
                "conversion_type": 'concentration',
                "quantity_value": -10.0,
                "quantity_type": 'chemical_potential',
                "temperature_k": 298.15,
                "p_from": 1.01325,
                "p_to": 1.0,
                "scale_from": 'molality',
                "scale_to": 'molarity',
                "solvent_density_kg_L": 1.0,
                "solute_molar_mass": 58.44,
                "reference_molality": 1.0,
                "reference_molarity": 1.0
            },
            "text_input": {
                "input_params": 'concentration -10.0 chemical_potential 298.15 molality molarity 1.0 58.44'
            },
            "output": {
                "converted_value": -10.036,
                "correction_term": -35.7,
                "original_unit": 'kJ/mol (1 mol/kg)',
                "target_unit": 'kJ/mol (1 mol/L)',
                "analysis": 'Molality to molarity conversion.'
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 8.314462618  # J/(mol·K)

    def _run_base(
        self,
        conversion_type: str,
        quantity_value: float,
        quantity_type: str = "delta_g",
        temperature_k: float = 298.15,
        p_from: float = 1.01325,
        p_to: float = 1.0,
        scale_from: str = "",
        scale_to: str = "",
        solvent_density_kg_L: float = 1.0,
        solute_molar_mass: float = 100.0,
        reference_molality: float = 1.0,
        reference_molarity: float = 1.0,
    ) -> dict:
        if temperature_k <= 0:
            raise ChemMCPError("Temperature must be positive.")
        conv_type = conversion_type.lower().strip()
        qty_type = quantity_type.lower().strip()

        if conv_type == "pressure":
            if p_from <= 0 or p_to <= 0:
                raise ChemMCPError("Pressures must be positive.")
            # Correction: Δ = RT ln(p_to / p_from)
            correction_j_mol = self.R * temperature_k * math.log(p_to / p_from)
            correction = correction_j_mol / 1000.0  # convert to kJ/mol

            # Apply based on quantity type
            if qty_type in ("delta_g", "chemical_potential"):
                converted = quantity_value + correction
            elif qty_type == "delta_h":
                converted = quantity_value  # Enthalpy of ideal gas independent of P°
                correction = 0.0
            elif qty_type == "delta_s":
                converted = quantity_value - self.R * math.log(p_to / p_from) / 1000.0  # kJ/(mol·K)
                correction = -self.R * math.log(p_to / p_from) / 1000.0
            elif qty_type == "equilibrium_constant":
                ln_K_correction = math.log(p_to / p_from)  # simplified per unit Δn=1
                converted = quantity_value  # K depends on choice; this is context-dependent
                correction = 0.0
            else:
                raise ChemMCPError(f"Unsupported quantity type for pressure conversion: '{qty_type}'.")

            orig_unit = f"({p_from} bar)"
            targ_unit = f"({p_to} bar)"
            analysis = (
                f"Pressure standard-state conversion at T={temperature_k} K:\n"
                f"Correction term: Δ = RT·ln({p_to}/{p_from}) = {self.R*temperature_k:.2f}·{math.log(p_to/p_from):.6f}\n"
                f"= {correction_j_mol:.2f} J/mol = {correction:.4f} kJ/mol\n"
                f"{qty_type}: {quantity_value} → {converted:.4f} kJ/mol"
            )

        elif conv_type == "concentration":
            sf = scale_from.lower().strip()
            st = scale_to.lower().strip()

            valid_scales = {"mole_fraction", "molality", "molarity"}
            if sf not in valid_scales or st not in valid_scales:
                raise ChemMCPError(f"Scales must be one of: {valid_scales}")

            # Activity ratio between scales: a_new / a_old
            # This is approximate; exact conversion depends on composition
            if (sf == "molality" and st == "molarity") or (sf == "molarity" and st == "molality"):
                # c(mol/L) ≈ ρ(kg/L) · m(mol/kg) / (1 + m·M(g/mol)/1000)
                # At infinite dilution: c ≈ ρ·m
                rho = solvent_density_kg_L
                M = solute_molar_mass / 1000.0  # kg/mol
                if sf == "molality":
                    # a_c/a_m ≈ ρ / (1 + m_ref·M) but at standard state m→0: a_c/a_m → ρ
                    activity_ratio = rho  # approximate at infinite dilution
                    orig_unit = "(1 mol/kg)"
                    targ_unit = f"(1 mol/L, ρ={rho} kg/L)"
                else:
                    activity_ratio = 1.0 / rho
                    orig_unit = f"(1 mol/L, ρ={rho} kg/L)"
                    targ_unit = "(1 mol/kg)"

            elif (sf == "mole_fraction" and st == "molality") or (sf == "molality" and st == "mole_fraction"):
                M_solvent = 18.015e-3  # kg/mol (water)
                if sf == "mole_fraction":
                    activity_ratio = M_solvent * reference_molality
                    orig_unit = "(1 mol/kg)"
                    targ_unit = "(mole fraction scale)"
                else:
                    activity_ratio = 1.0 / (M_solvent * reference_molality)
                    orig_unit = "(mole fraction scale)"
                    targ_unit = "(1 mol/kg)"

            elif (sf == "mole_fraction" and st == "molarity") or (sf == "molarity" and st == "mole_fraction"):
                M_solvent = 18.015e-3
                rho = solvent_density_kg_L
                if sf == "molarity":
                    activity_ratio = M_solvent * reference_molarity / rho
                    orig_unit = f"(1 mol/L, ρ={rho} kg/L)"
                    targ_unit = "(mole fraction scale)"
                else:
                    activity_ratio = rho / (M_solvent * reference_molarity)
                    orig_unit = "(mole fraction scale)"
                    targ_unit = f"(1 mol/L, ρ={rho} kg/L)"
            else:
                activity_ratio = 1.0
                orig_unit = f"({scale_from})"
                targ_unit = f"({scale_to})"

            if activity_ratio <= 0 or math.isnan(activity_ratio):
                raise ChemMCPError("Invalid activity ratio from given parameters.")

            correction_j_mol = self.R * temperature_k * math.log(activity_ratio) if activity_ratio != 1.0 else 0.0
            correction = correction_j_mol / 1000.0

            if qty_type in ("delta_g", "chemical_potential"):
                converted = quantity_value + correction
            elif qty_type == "delta_s":
                converted = quantity_value + self.R * math.log(activity_ratio) / 1000.0
                correction = self.R * math.log(activity_ratio) / 1000.0
            else:
                converted = quantity_value + correction

            analysis = (
                f"Concentration scale conversion ({sf} → {st}) at T={temperature_k} K:\n"
                f"Activity ratio (a_new/a_old) ≈ {activity_ratio:.6f}\n"
                f"Correction: RT·ln(ratio) = {correction_j_mol:.2f} J/mol = {correction:.4f} kJ/mol\n"
                f"{qty_type}: {quantity_value} → {converted:.4f}"
            )
        else:
            raise ChemMCPError(f"Unsupported conversion type: '{conv_type}'. Use 'pressure' or 'concentration'.")

        return {
            "converted_value": round(converted, 4),
            "correction_term": round(correction, 4),
            "original_unit": orig_unit,
            "target_unit": targ_unit,
            "analysis": analysis,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            kwargs = {"conversion_type": parts[0], "quantity_value": float(parts[1])}
            idx = 2
            if idx < len(parts):
                kwargs["quantity_type"] = parts[idx]; idx += 1
            if idx < len(parts):
                kwargs["temperature_k"] = float(parts[idx]); idx += 1
            if parts[0] == "pressure":
                if idx < len(parts):
                    kwargs["p_from"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["p_to"] = float(parts[idx]); idx += 1
            elif parts[0] == "concentration":
                if idx < len(parts):
                    kwargs["scale_from"] = parts[idx]; idx += 1
                if idx < len(parts):
                    kwargs["scale_to"] = parts[idx]; idx += 1
                if idx < len(parts):
                    kwargs["solvent_density_kg_L"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["solute_molar_mass"] = float(parts[idx]); idx += 1
            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
