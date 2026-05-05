import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class CompressibilityFactor(BaseTool):
    """
    压缩因子计算工具 (MCP #299)。
    计算真实气体的压缩因子 Z = PV/nRT，并分析其物理意义：
    - Z = 1: 理想气体行为
    - Z < 1: 分子间吸引力占主导（实际压力 < 理想压力）
    - Z > 1: 分子间排斥力（体积效应）占主导
    
    支持多种计算方式：
    - 从 P, V, n, T 直接计算
    - 通过维里方程近似：Z = 1 + B(T)/Vm + C(T)/Vm² + ...
    - 对比态原理：Z = f(Pr, Tr)
    - 范德华方程推导
    """
    __version__ = "0.1.0"
    name = "CompressibilityFactor"
    func_name = "calculate_compressibility_factor"
    description = "Calculate compressibility factor Z=PV/nRT for real gases using direct measurement, virial equation, van der Waals, or corresponding states principle."
    implementation_description = (
        "Multiple methods for Z calculation: (1) Direct from P,V,T data, "
        "(2) Virial expansion Z=1+B/Vm+C/Vm², (3) van der Waals derived, "
        "(4) Corresponding states approximation. Includes second virial coefficient B(T) for common gases."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Compressibility Factor", "Real Gas", "Virial Equation", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("method", "str", "direct", "Calculation method: 'direct', 'virial', 'vdw', 'corresponding_states'."),
        ("P", "float", "None", "Pressure in atm."),
        ("T", "float", "None", "Temperature in K."),
        ("V", "float", "None", "Volume in L."),
        ("n", "float", "1.0", "Amount in mol."),
        ("gas_species", "str", "None", "Gas name for virial coefficient lookup (e.g., 'N2', 'CO2', 'Ar', 'CH4')."),
        ("Vm", "float", "None", "Molar volume in L/mol (alternative to V+n)."),
        ("Pr", "float", "None", "Reduced pressure (for corresponding states)."),
        ("Tr", "float", "None", "Reduced temperature (for corresponding states)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query like: 'N2 at P=100atm T=300K', 'CO2 Z at critical point', 'CH4 virial T=273 Vm=0.5'."),
    ]

    output_sig = [
        ("Z", "float", "Compressibility factor (dimensionless)."),
        ("interpretation", "str", "Physical meaning of the Z value and dominant molecular interactions."),
        ("method_used", "str", "Which calculation method was used."),
        ("details", "dict", "Additional details: virial coefficients, reduced parameters, etc."),
        ("deviation_from_ideal", "str", "How much the gas deviates from ideal behavior."),
    ]

    examples = [{'code_input': {'method': 'direct', 'P': 100.0, 'T': 300.0, 'V': 2.0, 'n': 1.0, 'Pr': 'N/A', 'Tr': 'N/A', 'Vm': 'N/A', 'gas_species': 'N/A'}, 'text_input': {'query': 'gas at P=100atm T=300K V=2L n=1mol'}, 'output': {'Z': '<value>', 'details': 'N/A', 'deviation_from_ideal': 'N/A', 'interpretation': 'N/A', 'method_used': 'N/A'}}, {'code_input': {'method': 'virial', 'gas_species': 'N2', 'T': 273.15, 'Vm': 0.5, 'P': 'N/A', 'Pr': 'N/A', 'Tr': 'N/A', 'V': 'N/A', 'n': 'N/A'}, 'text_input': {'query': 'N2 virial T=273K Vm=0.5L/mol'}, 'output': {'Z': '<value>', 'details': 'N/A', 'deviation_from_ideal': 'N/A', 'interpretation': 'N/A', 'method_used': 'N/A'}}]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 0.08206  # L·atm/(K·mol)

        # Second virial coefficients B(T) in cm³/mol at various temperatures
        # Approximate values from experimental data
        # Format: {T(K): B(cm³/mol)}
        self._B_coefficients = {
            "He":   {100: 11.7, 200: 11.3, 300: 11.9, 400: 12.9, 500: 14.0},
            "H2":   {100: -14.8, 200: -7.1, 300: -1.6, 400: 2.8, 500: 6.6},
            "Ne":   {100: -5.9, 200: -0.5, 300: 3.6, 400: 6.6, 500: 8.9},
            "N2":   {100: -159, 200: -45, 300: -4.2, 400: 16.0, 500: 31.0},
            "O2":   {100: -198, 200: -56, 300: -9.0, 400: 13.0, 500: 28.0},
            "Ar":   {100: -187, 200: -35, 300: -21.5, 400: -12.0, 500: -4.5},
            "CO":   {100: -165, 200: -43, 300: -7.0, 400: 14.0, 500: 28.0},
            "CO2":  {250: -182, 300: -126, 350: -85, 400: -55, 450: -32, 500: -14},
            "CH4":  {200: -133, 250: -95, 300: -64, 400: -22, 500: 8},
            "NH3":  {250: -220, 300: -150, 350: -105, 400: -70, 500: -30},
            "H2O":  {373: -270, 473: -115, 573: -40, 673: 10},  # steam
            "C2H6":  {250: -310, 300: -210, 350: -145, 400: -95, 500: -35},
            "C2H4":  {250: -218, 300: -140, 350: -90, 400: -50, 500: -15},
            "air(approx)": {200: -38, 300: -7, 400: 12, 500: 26},
            "Kr":   {200: -48, 300: -28, 400: -14, 500: -2},
            "Xe":   {200: -135, 300: -88, 400: -55, 500: -30},
            "SF6":  {300: -180, 350: -120, 400: -75, 500: -25},
        }

        # Critical constants for reduced parameter calculation: (Tc(K), Pc(atm))
        self._critical_constants = {
            "He":     (5.19, 2.24),
            "H2":     (33.18, 12.97),  # quantum corrected
            "Ne":     (44.4, 26.9),
            "N2":     (126.2, 33.94),
            "O2":     (154.58, 49.77),
            "Ar":     (150.87, 48.0),
            "CO":     (132.91, 34.99),
            "CO2":    (304.13, 72.86),
            "NH3":    (405.4, 112.8),
            "H2O":    (647.1, 217.7),
            "CH4":    (190.56, 45.99),
            "C2H6":   (305.32, 48.72),
            "C2H4":   (282.34, 50.41),
            "C3H8":   (369.83, 42.01),
            "CH3OH":  (512.65, 80.84),
            "SO2":    (430.64, 77.81),
            "NO2":    (431.35, 101.33),
            "Cl2":    (417.0, 76.0),
            "Br2":    (588.0, 102.0),
            "benzene":(562.05, 48.95),
            "toluene":(591.75, 41.08),
            "CCl4":   (556.4, 44.98),
            "SF6":    (318.73, 37.59),
            "Xe":     (289.73, 58.40),
            "Kr":     (209.41, 54.96),
            "air(approx)":(132.45, 37.21),
        }

    def _run_base(self, method="direct", **kwargs) -> dict:
        method = method.lower().strip()
        if method == "direct":
            return self._calc_direct(**kwargs)
        elif method in ("virial", "virial_expansion"):
            return self._calc_virial(**kwargs)
        elif method in ("vdw", "vanderwaals"):
            return self._calc_z_vdw(**kwargs)
        elif method in ("corresponding_states", "cs"):
            return self._calc_corresponding_states(**kwargs)
        else:
            raise ChemMCPInputError(f"Unknown method '{method}'. Use: direct, virial, vdw, or corresponding_states.")

    def _calc_direct(self, P=None, T=None, V=None, n=1.0, **kwargs) -> dict:
        """Z = PV/(nRT) directly from measured P,V,T."""
        missing = [x for x in ["P", "T", "V"] if locals()[x] is None and kwargs.get(x) is None]
        P = P or kwargs.get("P")
        T = T or kwargs.get("T")
        V = V or kwargs.get("V")
        n = kwargs.get("n", n)

        if any(x is None for x in [P, T, V]):
            raise ChemMCPInputError(f"Direct method requires P, T, and V. Missing: {[x for x in ['P','T','V'] if eval(x) is None]}")

        R = self.R
        Z = P * V / (n * R * T)

        interpretation, deviation = self._interpret_z(Z)

        return {
            "Z": round(Z, 6),
            "P_atm": P, "T_K": T, "V_L": V, "n_mol": n,
            "formula": f"Z = ({P})({V})/(({n})×{R}×{T})",
            "interpretation": interpretation,
            "deviation_from_ideal": deviation,
            "method_used": "direct ideal gas definition Z=PV/nRT",
        }

    def _calc_virial(self, gas_species=None, T=None, Vm=None, n=None, V=None, **kwargs) -> dict:
        """Z = 1 + B(T)/Vm + C(T)/Vm² (using only B coefficient typically)."""
        species = gas_species or kwargs.get("gas_species")
        T = T or kwargs.get("T")
        Vm = Vm or kwargs.get("Vm")

        if not species:
            raise ChemMCPInputError("Virial method requires gas_species.")
        if T is None:
            raise ChemMCPInputError("Virial method requires temperature T.")
        if Vm is None and V is not None and n is not None:
            Vm = V / n
        if Vm is None:
            raise ChemMCPInputError("Virial method requires molar volume Vm (or V and n).")

        # Look up B coefficient
        species_key = None
        for key in self._B_coefficients:
            if key.lower() == species.lower() or key.replace(" ", "").lower() == species.replace(" ", "").lower():
                species_key = key
                break

        if species_key is None:
            raise ChemMCPError(
                f"No virial coefficient data for '{species}'. "
                f"Available: {list(self._B_coefficients.keys())}"
            )

        B_data = self._B_coefficients[species_key]
        temps = sorted(B_data.keys())

        # Interpolate B(T)
        if T in B_data:
            B_cm3mol = B_data[T]
        elif T < temps[0]:
            B_cm3mol = B_data[temps[0]]
        elif T > temps[-1]:
            B_cm3mol = B_data[temps[-1]]
        else:
            # Linear interpolation
            for i in range(len(temps) - 1):
                if temps[i] <= T <= temps[i+1]:
                    t1, t2 = temps[i], temps[i+1]
                    b1, b2 = B_data[t1], B_data[t2]
                    B_cm3mol = b1 + (b2 - b1) * (T - t1) / (t2 - t1)
                    break

        # Convert B from cm³/mol to L/mol
        B_Lmol = B_cm3mol / 1000.0

        # Z ≈ 1 + B/Vm (truncated after first term; sufficient for moderate pressures)
        Z = 1 + B_Lmol / Vm

        interpretation, deviation = self._interpret_z(Z)

        return {
            "Z": round(Z, 6),
            "gas_species": species_key,
            "temperature_K": T,
            "molar_volume_Vm_L_mol": Vm,
            "B_T_L_per_mol": round(B_Lmol, 6),
            "formula": f"Z = 1 + ({B_Lmol:.6f})/{Vm} = {Z:.6f}",
            "note": "Using truncated virial expansion (B term only). Accurate for low-moderate densities.",
            "details": {"B_T_cm3_per_mol": round(B_Lmol * 1000, 4)},
            "method_used": "virial equation (truncated B-term)",
            "deviation_from_ideal": round((Z - 1.0) * 100, 4)
            }

    def _calc_z_vdw(self, gas_species=None, T=None, Vm=None, **kwargs) -> dict:
        """Derive Z from van der Waals equation."""
        # Import a,b from vdw tool's constants (simplified inline)
        # Z = Vm/(Vm-b) - a/(RT*Vm)
        # Using approximate a,b values
        ab_data = {
            "He": (0.0346, 0.0238), "H2": (0.2476, 0.02661), "N2": (1.370, 0.0387),
            "O2": (1.382, 0.03186), "CO2": (3.640, 0.04267), "NH3": (4.225, 0.03707),
            "CH4": (2.303, 0.04310), "H2O": (5.536, 0.03049), "Ar": (1.355, 0.0320),
        }
        species = gas_species or kwargs.get("gas_species", "N2")

        found = None
        for key in ab_data:
            if key.lower() == species.lower():
                found = key
                break
        if not found:
            found = "N2"

        a, b = ab_data[found]
        R = self.R

        if T is None:
            T = 298.15
        if Vm is None:
            Vm = R * T / 1.0  # ~24 L/mol at 1 atm

        Z = Vm / (Vm - b) - a / (R * T * Vm)

        interpretation, deviation = self._interpret_z(Z)

        return {
            "Z": round(Z, 6),
            "gas_species": found,
            "a_atm_L2_mol2": a, "b_L_mol": b,
            "T_K": T, "Vm_L_mol": Vm,
            "formula": f"Z = {Vm}/({Vm}-{b}) - ({a})/({R}×{T}×{Vm})",
            "interpretation": interpretation,
            "deviation_from_ideal": deviation,
            "method_used": "van der Waals equation Z = Vm/(Vm-b) - a/(RT×Vm)",
        }

    def _calc_corresponding_states(self, Pr=None, Tr=None, gas_species=None, P=None, T=None, **kwargs) -> dict:
        """Approximate Z from generalized compressibility chart (analytical approximation)."""
        # Look up critical constants
        species = gas_species or kwargs.get("gas_species", "N2")
        tc_data = self._critical_constants.get(species) or self._critical_constants.get("N2")
        Tc, Pc = tc_data

        if Pr is None and P is not None:
            Pr = P / Pc
        if Tr is None and T is not None:
            Tr = T / Tc

        if Pr is None or Tr is None:
            raise ChemMCPInputError("Corresponding states method requires either (Pr, Tr) or (gas_species, P, T).")

        # Analytical approximation of generalized compressibility chart
        # Based on correlation: Z = 1 + (Pr/Tr)(B0 + ω*B1) where ω is acentric factor
        # Simplified: use Newton-Raphson inspired approximation
        # For many gases, a reasonable approximation at moderate conditions:

        # Simplified Hougen-Watson type approximation
        B0 = 0.083 - 0.422 / Tr**1.6
        B1 = 0.139 - 0.172 / Tr**4.2
        # Assume acentric factor ω ≈ 0 for simplicity (noble gas-like)
        omega = 0.0  # simplified

        Z_approx = 1 + (Pr / Tr) * (B0 + omega * B1)

        # Clamp to physically reasonable range
        Z_approx = max(0.1, min(Z_approx, 5.0))

        interpretation, deviation = self._interpret_z(Z_approx)

        return {
            "Z": round(Z_approx, 4),
            "gas_species": species,
            "Tc_K": Tc, "Pc_atm": Pc,
            "reduced_temperature_Tr": round(Tr, 4),
            "reduced_pressure_Pr": round(Pr, 4),
            "B0": round(B0, 4), "B1": round(B1, 4),
            "acentric_factor_omega_assumed": omega,
            "note": "Simplified correlation (ω≈0). For accurate results, use actual compressibility charts or equations of state.",
            "details": {
                "reduced_temperature_Tr": round(Tr, 4),
                "reduced_pressure_Pr": round(Pr, 4),
                "method": "generalized compressibility chart (Hougen-Watson approximation)",
            },
            "method_used": "corresponding states principle (generalized compressibility chart)"
        }

    @staticmethod
    def _interpret_z(Z: float) -> tuple:
        """Return human-readable interpretation of Z value."""
        if abs(Z - 1.0) < 0.02:
            interp = (
                f"Z = {Z:.4f} ≈ 1.00: The gas behaves nearly as an ideal gas under these conditions. "
                f"Molecular attractions and excluded volume effects approximately cancel out."
            )
            dev = "Negligible deviation (<2%) from ideal behavior."
        elif Z < 0.99:
            if Z < 0.5:
                interp = (
                    f"Z = {Z:.4f} << 1: Strongly non-ideal! Attractive intermolecular forces dominate significantly. "
                    f"The gas is much more compressible than an ideal gas. This occurs at low temperatures "
                    f"and/or high pressures near the condensation region."
                )
                dev = f"Large negative deviation ({(Z-1)*100:.1f}%). Attractive forces strongly reduce pressure."
            else:
                interp = (
                    f"Z = {Z:.4f} < 1: Intermolecular attractive forces dominate over repulsive (size) effects. "
                    f"The real gas exerts less pressure than predicted by ideal gas law because molecules "
                    f"are pulled together by intermolecular forces."
                )
                dev = f"Moderate negative deviation ({(Z-1)*100:.1f}%). Attractions reduce effective pressure."
        else:  # Z > 1.01
            if Z > 2.0:
                interp = (
                    f"Z = {Z:.4f} >> 1: Strongly non-ideal! Molecular excluded volume (repulsive) effects dominate. "
                    f"At high pressure, the finite size of molecules means less free volume than ideal. "
                    f"The gas is harder to compress than an ideal gas."
                )
                dev = f"Large positive deviation ({(Z-1)*100:.1f}%). Excluded volume dominates."
            else:
                interp = (
                    f"Z = {Z:.4f} > 1: Repulsive/excluded volume effects dominate over attractions. "
                    f"At these conditions (typically high P or very small molecules like He/H₂ at low T), "
                    f"the finite molecular size reduces available volume, increasing effective pressure."
                )
                dev = f"Moderate positive deviation ({(Z-1)*100:.1f}%). Excluded volume increases effective pressure."

        return interp, dev

    def _run_text(self, query: str) -> dict:
        q = query.strip()
        # Simple parsing
        import re

        # Detect method
        method = "direct"
        if "virial" in q.lower():
            method = "virial"
        elif "vdw" in q.lower() or "van der" in q.lower():
            method = "vdw"
        elif "corresponding" in q.lower() or "reduced" in q.lower():
            method = "corresponding_states"

        # Extract numbers
        P_match = re.search(r'P[=\s]?(\d+\.?\d*)', q)
        T_match = re.search(r'T[=\s]?(\d+\.?\d*)\s*K?', q)
        V_match = re.search(r'V[=\s]?(\d+\.?\d*)', q)
        Vm_match = re.search(r'Vm[=\s]?(\d+\.?\d*)', q)

        kwargs = {}
        if P_match: kwargs["P"] = float(P_match.group(1))
        if T_match: kwargs["T"] = float(T_match.group(1))
        if V_match: kwargs["V"] = float(V_match.group(1))
        if Vm_match: kwargs["Vm"] = float(Vm_match.group(1))

        # Try to find gas species
        for gas in list(self._B_coefficients.keys())[:20]:
            if gas.lower() in q.lower():
                kwargs["gas_species"] = gas
                break

        return self._run_base(method=method, **kwargs)

