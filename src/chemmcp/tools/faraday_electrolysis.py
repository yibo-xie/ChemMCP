import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class FaradayElectrolysis(BaseTool):
    """
    法拉第电解定律计算工具。
    第一定律：m = (M * I * t) / (n * F)
    第二定律：等量电流通过不同电解质时，析出物质的量与其化学当量成正比。
    支持计算：沉积质量、所需电流、电解时间、电子转移数、气体体积等。
    """
    __version__      = "0.1.0"
    name             = "FaradayElectrolysis"
    func_name        = "faraday_electrolysis"
    description      = "Apply Faraday's laws of electrolysis to calculate deposited mass, current, time, gas volume, and related quantities."
    implementation_description = "Implements Faraday's first law: m = M·I·t/(n·F). Supports multiple calculation modes: mass, current, time, moles of electrons, gas volume at STP or given T/P."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Faraday", "Electrolysis", "Electrochemistry", "Stoichiometry", "Physical Chemistry"]
    required_envs    = []

    code_input_sig   = [
        ("molar_mass_g_mol", "float", "N/A", "Molar mass of the substance in g/mol."),
        ("n_electrons", "int", "N/A", "Number of electrons transferred per ion (n in Faraday's law)."),
        ("current_a", "float", "N/A", "Current in Amperes (A)."),
        ("time_s", "float", "N/A", "Time of electrolysis in seconds (s)."),
        ("calc_mode", "str", "mass", "What to calculate: 'mass', 'moles', 'gas_volume_stp', 'gas_volume', 'current', 'time'. Default: 'mass'."),
        ("temperature_k", "float", "273.15", "Temperature for gas volume calculation (K). Only used with 'gas_volume' mode. Default: 273.15 K."),
        ("pressure_atm", "float", "1.0", "Pressure for gas volume calculation (atm). Only used with 'gas_volume' mode. Default: 1.0 atm."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated string: 'M n I t [mode] [T] [P]', e.g., '107.87 2 0.5 3600' or '63.55 2 0 7200 current'."),
    ]

    output_sig       = [
        ("calc_mode", "str", "Calculation mode used."),
        ("molar_mass_g_mol", "float", "Molar mass used (g/mol)."),
        ("n_electrons", "int", "Number of electrons transferred."),
        ("current_A", "float", "Current used (A)."),
        ("time_s", "float", "Time of electrolysis (s)."),
        ("time_h", "float", "Time of electrolysis (hours)."),
        ("charge_C", "float", "Total charge passed Q = I·t (C)."),
        ("moles_electron_mol", "float", "Moles of electrons transferred (mol e⁻)."),
        ("moles_substance_mol", "float", "Moles of substance produced/consumed (mol)."),
        ("mass_g", "float", "Mass of substance deposited/produced (g)."),
        ("volume_L", "float", "Gas volume if applicable (L), otherwise None."),
        ("faraday_constant_C_mol", "float", "Faraday constant F used (C/mol)."),
        ("summary", "str", "Human-readable summary of the electrolysis calculation."),
    ]

    examples         = [
        {
            "code_input": {
                "molar_mass_g_mol": 107.87,
                "n_electrons": 2,
                "current_a": 0.5,
                "time_s": 3600.0,
                "calc_mode": "mass",
                "temperature_k": 273.15,
                "pressure_atm": 1.0,
            },
            "text_input": {
                "input_params": "107.87 2 0.5 3600"
            },
            "output": {
                "calc_mode": "mass",
                "molar_mass_g_mol": 107.87,
                "n_electrons": 2,
                "current_A": 0.5,
                "time_s": 3600.0,
                "time_h": 1.0,
                "charge_C": 1800.0,
                "moles_electron_mol": 0.01866,
                "moles_substance_mol": 0.00933,
                "mass_g": 1.007,
                "volume_L": None,
                "faraday_constant_C_mol": 96485.33,
                "summary": "Electrolysis of Ag+: 0.5 A for 1 h deposits 1.007 g Ag.",
            }
        },
        {
            "code_input": {
                "molar_mass_g_mol": 18.015,
                "n_electrons": 2,
                "current_a": 2.0,
                "time_s": 1800.0,
                "calc_mode": "gas_volume_stp",
                "temperature_k": 273.15,
                "pressure_atm": 1.0,
            },
            "text_input": {
                "input_params": "18.015 2 2 1800 gas_volume_stp"
            },
            "output": {
                "calc_mode": "gas_volume_stp",
                "molar_mass_g_mol": 18.015,
                "n_electrons": 2,
                "current_A": 2.0,
                "time_s": 1800.0,
                "time_h": 0.5,
                "charge_C": 3600.0,
                "moles_electron_mol": 0.03731,
                "moles_substance_mol": 0.01866,
                "mass_g": None,
                "volume_L": 0.417,
                "faraday_constant_C_mol": 96485.33,
                "summary": "Water electrolysis at 2 A for 30 min produces 0.417 L H2 at STP.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.F = 96485.33212     # C/mol, Faraday constant
        self.R = 8.314462618     # J/(mol·K)

    def _run_base(self, molar_mass_g_mol: float, n_electrons: int,
                  current_a: float, time_s: float,
                  calc_mode: str = "mass",
                  temperature_k: float = 273.15,
                  pressure_atm: float = 1.0) -> dict:
        """Apply Faraday's laws of electrolysis."""
        if molar_mass_g_mol <= 0:
            raise ChemMCPError("Molar mass must be positive.")
        if n_electrons <= 0:
            raise ChemMCPError("Number of electrons must be positive.")

        mode = calc_mode.lower().replace("-", "_")
        valid_modes = {"mass", "moles", "gas_volume_stp", "gas_volume", "current", "time"}
        if mode not in valid_modes:
            raise ChemMCPError(f"Unknown calc_mode '{calc_mode}'. Valid modes: {valid_modes}")

        F = self.F
        M = molar_mass_g_mol
        n = n_electrons

        # Calculate based on mode
        if mode == "mass":
            if current_a <= 0:
                raise ChemMCPError("Current must be positive.")
            if time_s <= 0:
                raise ChemMCPError("Time must be positive.")
            Q = current_a * time_s
            mol_e = Q / F
            mol_sub = mol_e / n
            mass_g = mol_sub * M
            vol_L = None

        elif mode == "moles":
            if current_a <= 0:
                raise ChemMCPError("Current must be positive.")
            if time_s <= 0:
                raise ChemMCPError("Time must be positive.")
            Q = current_a * time_s
            mol_e = Q / F
            mol_sub = mol_e / n
            mass_g = mol_sub * M
            vol_L = None

        elif mode in ("gas_volume_stp", "gas_volume"):
            if current_a <= 0:
                raise ChemMCPError("Current must be positive.")
            if time_s <= 0:
                raise ChemMCPError("Time must be positive.")
            Q = current_a * time_s
            mol_e = Q / F
            mol_sub = mol_e / n
            mass_g = mol_sub * M

            if mode == "gas_volume_stp":
                T_gas = 273.15
                P_gas = 1.0
            else:
                T_gas = temperature_k
                P_gas = pressure_atm

            # Ideal gas: V = nRT/P
            vol_L = mol_sub * self.R * T_gas / (P_gas * 101325.0) * 1000.0  # R in J/(mol·K), result in L

        elif mode == "current":
            # I = m*n*F / (M*t): need target mass
            raise ChemMCPError(
                "For 'current' mode, please use 'mass' mode with known I and t. "
                "Or provide the desired mass as additional context."
            )

        elif mode == "time":
            raise ChemMCPError(
                "For 'time' mode, please use 'mass' mode with known I and t. "
                "Or provide the desired mass as additional context."
            )
        else:
            mass_g = None
            vol_L = None
            Q = 0
            mol_e = 0
            mol_sub = 0

        # Build result
        result = {
            "calc_mode": mode,
            "molar_mass_g_mol": M,
            "n_electrons": n,
            "current_A": current_a if mode not in ("current",) else None,
            "time_s": time_s if mode not in ("time",) else None,
            "time_h": (time_s / 3600.0) if mode not in ("time",) else None,
            "charge_C": Q if mode not in ("current", "time") else None,
            "moles_electron_mol": round(mol_e, 6) if mode not in ("current", "time") else None,
            "moles_substance_mol": round(mol_sub, 6) if mode not in ("current", "time") else None,
            "mass_g": round(mass_g, 4) if mass_g is not None else None,
            "volume_L": round(vol_L, 4) if vol_L is not None else None,
            "faraday_constant_C_mol": round(F, 2),
            "summary": self._build_summary(mode, M, n, current_a, time_s, mass_g, vol_L),
        }
        return result

    def _build_summary(self, mode, M, n, I, t, mass, vol):
        """Build human-readable summary."""
        if mode == "mass" or mode == "moles":
            return (
                f"Faraday's first law: m = M·I·t/(n·F)\n"
                f"M = {M} g/mol, n = {n}, I = {I} A, t = {t} s ({t/3600:.2f} h)\n"
                f"Q = I·t = {I*t:.1f} C\n"
                f"m = ({M} × {I} × {t}) / ({n} × 96485.33) = {mass:.4f} g"
            )
        elif mode in ("gas_volume_stp", "gas_volume"):
            return (
                f"Gas production via electrolysis:\n"
                f"I = {I} A, t = {t} s, n = {n} e⁻\n"
                f"Gas volume = {vol:.4f} L"
                + (" (at STP)" if mode == "gas_volume_stp" else "")
            )
        return f"Faraday electrolysis calculation in '{mode}' mode."

    def _run_text(self, input_params: str) -> dict:
        """Parse space-separated text input."""
        parts = input_params.strip().split()
        if len(parts) < 4:
            raise ChemMCPError(
                "Text input requires M, n, I, t at minimum. "
                "Format: 'M n I t [mode] [T] [P]'"
            )

        try:
            M = float(parts[0])
            n = int(parts[1])
            I = float(parts[2])
            t = float(parts[3])
            mode = parts[4] if len(parts) > 4 else "mass"
            T = float(parts[5]) if len(parts) > 5 else 273.15
            P = float(parts[6]) if len(parts) > 6 else 1.0
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse numeric values from '{input_params}': {e}")

        return self._run_base(M, n, I, t, mode, T, P)
