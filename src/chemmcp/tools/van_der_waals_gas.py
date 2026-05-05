import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class VanDerWaalsGas(BaseTool):
    """
    范德华方程真实气体计算工具 (MCP #298)。
    使用范德华状态方程: (P + a(n/V)²)(V - nb) = nRT
    计算真实气体的压力/体积修正，并与理想气体对比。
    
    包含常见气体的范德华常数 a, b 数据库：
    - a: 分子间吸引力参数 (atm·L²/mol²)
    - b: 分子体积排斥参数 (L/mol)
    
    支持功能：
    - 给定 V, T, n 计算 P（或反之）
    - 压缩因子 Z = PV/nRT 计算
    - 与理想气体偏差百分比
    - 临界常数计算
    """
    __version__ = "0.1.0"
    name = "VanDerWaalsGas"
    func_name = "calculate_vdw_gas"
    description = "van der Waals real gas equation calculations: (P + an²/V²)(V-nb) = nRT. Includes a,b constants for common gases and ideal gas comparison."
    implementation_description = (
        "Implements the van der Waals equation of state with built-in a,b constants for 30+ gases. "
        "Calculates corrected pressure/volume, compressibility factor Z, deviation from ideal gas behavior, "
        "and critical constants (Tc=8a/(27Rb), Pc=a/(27b²), Vc=3nb)."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["van der Waals", "Real Gas", "Equation of State", "Physical Chemistry", "Compressibility"]
    required_envs = []

    code_input_sig = [
        ("gas_species", "str", "N/A", "Gas species name or formula for a,b constant lookup (e.g., 'H2O', 'CO2', 'He', 'N2', 'CH4')."),
        ("T", "float", "N/A", "Temperature in Kelvin."),
        ("V", "float", "None", "Volume in L (if solving for pressure)."),
        ("n", "float", "1.0", "Amount in mol (default 1 mol)."),
        ("P", "float", "None", "Pressure in atm (if solving for volume)."),
        ("mode", "str", "pressure", "Mode: 'pressure' (given V,T→P), 'volume' (given P,T→V), 'compare' (ideal vs vdw comparison)."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query like: 'CO2 at T=300K V=1L n=1mol', 'He at STP compare', 'NH3 critical constants'."),
    ]

    output_sig = [
        ("gas_species", "str", "The gas being analyzed."),
        ("vdw_pressure", "float", "Pressure from van der Waals equation (atm)."),
        ("ideal_pressure", "float", "Pressure from ideal gas law (atm)."),
        ("compressibility_factor_Z", "float", "Z = PV/(nRT). Z=1 is ideal; Z<1 attractive dominant; Z>1 repulsive dominant."),
        ("deviation_percent", "float", "Percentage deviation from ideal gas: (P_vdw - P_ideal)/P_ideal × 100%."),
        ("constants_a_b", "dict", "van der Waals a (atm·L²/mol²) and b (L/mol) constants."),
        ("critical_constants", "dict", "Critical temperature Tc(K), pressure Pc(atm), volume Vc(L/mol) if applicable."),
        ("explanation", "str", "Physical interpretation of deviations."),
    ]

    examples = [{'code_input': {'gas_species': 'CO2', 'T': 300.0, 'V': 1.0, 'n': 1.0, 'mode': 'pressure', 'P': 'N/A'}, 'text_input': {'query': 'CO2 T=300K V=1L n=1mol'}, 'output': {'gas_species': 'CO2', 'compressibility_factor_Z': 'N/A', 'constants_a_b': 'N/A', 'critical_constants': 'N/A', 'deviation_percent': 'N/A', 'explanation': 'N/A', 'ideal_pressure': 'N/A', 'vdw_pressure': 'N/A'}}, {'code_input': {'gas_species': 'He', 'T': 273.15, 'V': 22.4, 'n': 1.0, 'P': 'N/A', 'mode': 'N/A'}, 'text_input': {'query': 'He at STP compare'}, 'output': {'gas_species': 'He', 'compressibility_factor_Z': 'N/A', 'constants_a_b': 'N/A', 'critical_constants': 'N/A', 'deviation_percent': 'N/A', 'explanation': 'N/A', 'ideal_pressure': 'N/A', 'vdw_pressure': 'N/A'}}]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Build van der Waals constants database."""
        self.R = 0.08206  # L·atm/(K·mol)
        # (a in atm·L²/mol², b in L/mol)
        # Source: CRC Handbook / standard physical chemistry data
        self._ab_constants = {
            # Noble gases
            "He":      (0.0346, 0.0238),
            "Ne":      (0.211, 0.0171),
            "Ar":      (1.355, 0.0320),
            "Kr":      (2.325, 0.0396),
            "Xe":      (4.194, 0.0511),
            # Diatomic gases
            "H2":      (0.2476, 0.02661),
            "N2":      (1.370, 0.0387),
            "O2":      (1.382, 0.03186),
            "F2":      (1.171, 0.0296),
            "Cl2":     (6.579, 0.05622),
            "Br2":     (9.75, 0.0591),
            "I2":      (16.9, 0.0724),
            "CO":      (1.472, 0.0395),
            "NO":      (1.358, 0.0279),
            # Other inorganic gases
            "H2O":     (5.536, 0.03049),   # water vapor (steam)
            "H2S":     (4.490, 0.04287),
            "NH3":     (4.225, 0.03707),
            "HF":      (9.433, 0.0739),
            "HCl":      (3.716, 0.04081),
            "HBr":      (4.523, 0.0443),
            "HI":       (6.311, 0.0530),
            "CO2":     (3.640, 0.04267),
            "N2O":     (3.832, 0.04415),
            "SO2":     (6.803, 0.05636),
            "NO2":      (5.284, 0.04424),
            "O3":       (4.546, 0.0443),
            "Hg":       (8.200, 0.0170),
            # Light hydrocarbons
            "CH4":     (2.303, 0.04310),
            "C2H2":    (4.490, 0.05120),
            "C2H4":    (4.612, 0.05826),
            "C2H6":    (5.562, 0.06380),
            "C3H8":    (8.779, 0.08445),
            "C4H10(iso)":(12.87, 0.1142),
            # Other organic vapors
            "benzene": (18.78, 0.1197),
            "C6H6":    (18.78, 0.1197),
            "toluene": (20.59, 0.1385),
            "CCl4":    (19.82, 0.1268),
            "CHCl3":   (15.37, 0.1022),
            "CH3OH":    (9.706, 0.06702),
            "C2H5OH":   (12.18, 0.08407),
            "acetone":  (13.91, 0.0994),
            "(CH3)2CO":(13.91, 0.0994),
            "ether":    (17.06, 0.1224),
            "CS2":      (11.62, 0.0769),
            # Others
            "air(approx)": (1.36, 0.0365),
            "Neon":    (0.211, 0.0171),
            "krypton": (2.325, 0.0396),
            "xenon":   (4.194, 0.0511),
            "PH3":      (4.632, 0.05152),
            "AsH3":     (4.976, 0.0558),
            "SF6":      (5.520, 0.08660),
            "UF6":     (13.74, 0.1118),
            "CCl2F2":  (10.78, 0.0998),
            "CHClF2":   (6.184, 0.0818),
            "N2H4":     (5.313, 0.0461),
            "HCN":      (11.31, 0.0727),
            "acetylene":(4.490, 0.05120),
            "ethylene": (4.612, 0.05826),
            "ethane":   (5.562, 0.06380),
            "methane":  (2.303, 0.04310),
            "water vapor": (5.536, 0.03049),
            "ammonia":  (4.225, 0.03707),
            "carbon dioxide": (3.640, 0.04267),
        }

    def _run_base(self, gas_species: str, T: float, V: float = None, n: float = 1.0,
                 P: float = None, mode: str = "pressure") -> dict:
        """
        Main calculation using van der Waals equation.
        Mode 'pressure': given V, T, n → calculate P
        Mode 'volume': given P, T, n → solve for V (iterative)
        """
        R = self.R

        # Look up a, b constants
        species = gas_species.strip()
        if species not in self._ab_constants:
            # Case-insensitive lookup
            found = None
            for key in self._ab_constants:
                if key.lower() == species.lower() or key.replace(" ", "").lower() == species.replace(" ", "").lower():
                    found = key
                    break
            if not found:
                avail = sorted(self._ab_constants.keys())
                raise ChemMCPError(
                    f"Unknown gas species '{species}'. Available: {avail}\n"
                    f"You can also provide custom a,b values via the tool's advanced mode."
                )
            species = found

        a, b = self._ab_constants[species]

        # Critical constants
        Tc = 8 * a / (27 * R * b) if b > 0 else 0
        Pc = a / (27 * b**2) if b > 0 else 0
        Vc_mol = 3 * b  # molar critical volume

        mode = mode.lower().strip()

        if mode in ("pressure", "p"):
            if V is None:
                raise ChemMCPInputError("Volume V must be provided for pressure calculation mode.")
            if T <= 0:
                raise ChemMCPInputError("Temperature T must be > 0 K.")
            if V <= n * b:
                raise ChemMCPInputError(f"Volume V={V} L must be > nb={n*b:.4f} L (excluded volume).")

            # P_ideal = nRT/V
            P_ideal = n * R * T / V

            # P_vdw = nRT/(V-nb) - a(n/V)²
            P_vdw = n * R * T / (V - n * b) - a * (n / V) ** 2

            # Compressibility factor from vdw
            Z = P_vdw * V / (n * R * T)

            # Deviation
            if P_ideal != 0:
                dev_pct = (P_vdw - P_ideal) / abs(P_ideal) * 100
            else:
                dev_pct = 0

            # Physical interpretation
            if Z < 0.99:
                interpretation = (
                    f"Z = {Z:.4f} < 1: Attractive forces dominate. "
                    f"The a(n/V)² correction term reduces pressure compared to ideal gas. "
                    f"This is typical for {species} at relatively low T where intermolecular "
                    f"attractions are significant."
                )
            elif Z > 1.01:
                interpretation = (
                    f"Z = {Z:.4f} > 1: Repulsive (excluded volume) effects dominate. "
                    f"The molecular size (b term) increases effective pressure. "
                    f"This is typical for small molecules like He/H₂ at high pressure or low T."
                )
            else:
                interpretation = (
                    f"Z = {Z:.4f} ≈ 1: Gas behaves nearly ideally under these conditions. "
                    f"Both attractive and repulsive corrections are small."
                )

            return {
                "gas_species": species,
                "temperature_K": T,
                "volume_L": V,
                "amount_mol": n,
                "vdw_pressure_atm": round(P_vdw, 4),
                "ideal_pressure_atm": round(P_ideal, 4),
                "compressibility_factor_Z": round(Z, 6),
                "deviation_percent": round(dev_pct, 4),
                "constants_a_b": {
                    "a_atm_L2_mol-2": round(a, 4),
                    "b_L_mol-1": round(b, 6)
                },
                "critical_constants": {
                    "Tc_K": round(Tc, 2),
                    "Pc_atm": round(Pc, 2),
                    "Vc_L_per_mol": round(Vc_mol, 4),
                },
                "reduced_temperature_Tr": round(T / Tc, 4) if Tc > 0 else None,
                "interpretation": interpretation,
                "equation_used": f"(P + a(n/V)²)(V - nb) = nRT\nP = nRT/(V-nb) - a(n/V)²\nP = ({n})({R})({T})/({V}-{n}{b}) - ({a})({n}/{V})²",
            }

        elif mode in ("volume", "v"):
            if P is None:
                raise ChemMCPInputError("Pressure P must be provided for volume calculation mode.")
            # Solve van der Waals equation for V iteratively (Newton-Raphson)
            # Start with ideal gas estimate
            V_est = n * R * T / P
            for iteration in range(50):
                if V_est <= n * b:
                    V_est = n * b * 1.5
                f = (P + a * (n / V_est)**2) * (V_est - n * b) - n * R * T
                # df/dV
                df_dv = P + a * n**2 / V_est**2 - 2 * a * n**2 * (V_est - n * b) / V_est**3
                if abs(df_dv) < 1e-15:
                    break
                V_new = V_est - f / df_dv
                if abs(V_new - V_est) < 1e-8:
                    V_est = V_new
                    break
                V_est = V_new

            P_ideal_calc = n * R * T / V_est
            Z = P * V_est / (n * R * T)
            P_vdw_check = n * R * T / (V_est - n * b) - a * (n / V_est) ** 2

            return {
                "gas_species": species,
                "temperature_K": T,
                "pressure_atm": P,
                "volume_L": V_est,
                "vdw_volume_L": round(V_est, 6),
                "ideal_volume_L": round(P_ideal_calc, 4),
                "amount_mol": n,
                "vdw_pressure_atm": round(P_vdw_check, 4),
                "ideal_pressure_atm": round(P_ideal_calc, 4),
                "compressibility_factor_Z": round(Z, 6),
                "constants_a_b": {
                    "a_atm_L2_mol-2": round(a, 4),
                    "b_L_mol-1": round(b, 6)
                },
                "critical_constants": {
                    "Tc_K": round(Tc, 2),
                    "Pc_atm": round(Pc, 2),
                    "Vc_L_per_mol": round(Vc_mol, 4),
                },
                "iterations_used": iteration + 1,
            }

        elif mode == "compare":
            # Compare ideal vs vdw across a range of conditions
            if V is None:
                V = n * R * T / 1.0  # assume ~1 atm
            result = self._run_base(gas_species, T, V, n, None, "pressure")
            result["mode"] = "comparison"
            return result

        else:
            raise ChemMCPInputError(f"Unknown mode '{mode}'. Use: 'pressure', 'volume', or 'compare'.")

    def _run_text(self, query: str) -> dict:
        q = query.lower().strip()
        # Parse natural language queries
        import re

        # Extract gas name
        species = None
        for gas in list(self._ab_constants.keys())[:30]:
            if gas.lower() in q or gas.replace(" ", "").lower() in q.replace(" ", ""):
                species = gas
                break
        if species is None:
            raise ChemMCPInputError(
                f"Could not identify gas species in query '{q}'. "
                f"Specify gas name explicitly, e.g., 'CO2', 'He', 'NH3', 'CH4'."
            )

        # Extract T
        T_match = re.search(r'T[=\s]?(\d+\.?\d*)\s*K?', q)
        T_val = float(T_match.group(1)) if T_match else 298.15

        # Extract V
        V_match = re.search(r'V[=\s]?(\d+\.?\d*)\s*L?', q)
        V_val = float(V_match.group(1)) if V_match else None

        # Extract n
        n_match = re.search(r'n[=\s]?(\d+\.?\d*)\s*mol?', q)
        n_val = float(n_match.group(1)) if n_match else 1.0

        # Detect mode
        if "compare" in q or "comparison" in q:
            mode = "compare"
        elif "volume" in q or "solve.*v" in q:
            mode = "volume"
        else:
            mode = "pressure"

        if V_val is None and mode != "volume":
            # Default to 1 mol at ~1 atm → V ≈ RT/P
            V_val = self.R * T_val / 1.0

        return self._run_base(species, T_val, V_val, n_val, mode=mode)

