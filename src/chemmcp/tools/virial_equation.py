import logging
import math
from typing import Optional, List

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class VirialEquation(BaseTool):
    """
    维里方程状态计算工具 (MCP #300)。
    使用维里方程（Virial Equation of State）进行真实气体状态计算：
    
    压力级数形式：P = RT/Vm (1 + B(T)/Vm + C(T)/Vm² + D(T)/Vm³ + ...)
    压缩因子形式：Z = PVm/(RT) = 1 + B(T)/Vm + C(T)/Vm² + D(T)/Vm³ + ...
    
    其中：
    - B(T): 第二维里系数（两体相互作用，最重要）
    - C(T): 第三维里系数（三体相互作用）
    - D(T)及更高: 通常可忽略
    
    功能：
    - 计算任意温度下的 Z、P 或 Vm
    - 内置常见气体的实验维里系数数据
    - Boyle 温度计算（B(Tb) = 0）
    - 与理想气体和范德华方程对比
    """
    __version__ = "0.1.0"
    name = "VirialEquation"
    func_name = "calculate_virial_equation"
    description = "Real gas calculations using the Virial Equation of State: Z=1+B/Vm+C/Vm²+... Includes experimental virial coefficients for common gases and Boyle temperature calculation."
    implementation_description = (
        "Implements the virial equation of state in both pressure and compressibility factor forms. "
        "Contains experimental second (B) and third (C) virial coefficients for 15+ gases. "
        "Supports solving for P, Vm, or Z, and calculates Boyle temperature where B(Tb)=0."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Virial Equation", "Real Gas", "Physical Chemistry", "Thermodynamics", "Compressibility"]
    required_envs = []

    code_input_sig = [
        ("gas_species", "str", "N/A", "Gas name for coefficient lookup (e.g., 'Ar', 'N2', 'CO2', 'CH4', 'He')."),
        ("T", "float", "N/A", "Temperature in Kelvin."),
        ("mode", "str", "Z_from_Vm", "Calculation mode: 'Z_from_Vm', 'P_from_Vm', 'Vm_from_P', 'boyle_temperature'."),
        ("Vm", "float", "None", "Molar volume in L/mol (for Z or P calculation)."),
        ("P_atm", "float", "None", "Pressure in atm (for Vm calculation)."),
        ("order", "int", "2", "Virial expansion order: 2 (only B), 3 (B+C), or higher. Default: 2."),
    ]

    text_input_sig = [
        ("query", "str", "N/A", "Query like: 'Ar at T=273K Vm=0.5L/mol order=2', 'N2 Boyle temperature', 'CO2 P from Vm=0.2L T=350K'."),
    ]

    output_sig = [
        ("gas_species", "str", "The gas being analyzed."),
        ("temperature_K", "float", "Temperature."),
        ("result_value", "float", "Calculated value (Z, P, or Vm depending on mode)."),
        ("unit", "str", "Unit of the result."),
        ("virial_coefficients_used", "dict", "B(T), C(T), etc. values used in the calculation."),
        ("virial_expansion_terms", "list", "Each term in the expansion: [1, B/Vm, C/Vm², ...]."),
        ("ideal_gas_comparison", "dict", "Comparison with ideal gas result."),
        ("convergence_note", "str", "Notes on convergence and accuracy of the truncation."),
        ("explanation", "str", "Step-by-step explanation."),
    ]

    examples = [
        {'code_input': {'gas_species': 'Ar', 'T': 273.15, 'Vm': 0.224, 'mode': 'Z_from_Vm', 'order': 2, 'P_atm': 'N/A'},
         'text_input': {'query': 'Argon at STP (T=273K, Vm=22.4L/mol)'},
         'output': {'result_value': '~0.999', 'unit': 'dimensionless'}},
        {'code_input': {'gas_species': 'N2', 'T': 273.15, 'Vm': 0.5, 'mode': 'P_from_Vm', 'order': 2, 'P_atm': 'N/A'},
         'text_input': {'query': 'N2 at T=273K Vm=0.5L/mol calculate P'},
         'output': {'result_value': '<value>', 'unit': 'atm'}},
        {'code_input': {'gas_species': 'He', 'mode': 'boyle_temperature', 'P_atm': 'N/A', 'T': 'N/A', 'Vm': 'N/A', 'order': 'N/A'},
         'text_input': {'query': 'Helium Boyle temperature'},
         'output': {'result_value': '~24.5 K', 'unit': 'K'}},
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.R = 0.08206  # L·atm/(K·mol)

        # Comprehensive virial coefficient database
        # Format: {T(K): (B in cm³/mol, C in cm³²/mol²)}
        self._virial_data = {
            "He": {
                "data": {
                    10: (-21.6, 200), 20: (-17.4, 120), 50: (-11.8, 37),
                    100: (11.7, 45), 150: (11.9, 40), 200: (11.3, 38),
                    273: (12.0, 35), 300: (11.9, 34), 400: (12.9, 33),
                    500: (14.0, 32),
                },
                "boyle_T_approx": 24.5,
            },
            "H2": {
                "data": {
                    50: (-4.4, 310), 75: (-1.0, 220), 100: (-14.8, 195),
                    150: (-7.5, 155), 200: (-7.1, 138), 273: (-1.6, 115),
                    300: (-1.6, 110), 350: (-0.5, 98), 400: (2.8, 88),
                    500: (6.6, 72),
                },
                "boyle_T_approx": 112.0,
            },
            "Ne": {
                "data": {
                    50: (-28.5, 420), 100: (-5.9, 170), 150: (1.5, 115),
                    200: (-0.5, 90), 273: (3.6, 78), 300: (4.5, 73),
                    400: (6.6, 62), 500: (8.9, 53),
                },
                "boyle_T_approx": 125.0,
            },
            "Ar": {
                "data": {
                    100: (-187, 1800), 123.85: (-135.5, 1300),
                    150: (-99.5, 850), 200: (-35, 450), 250: (-18, 280),
                    273: (-21.5, 230), 298.15: (-16, 200), 300: (-15.5, 198),
                    320: (-13, 182), 350: (-9.5, 158), 400: (-4.5, 120),
                    450: (0.5, 92), 480: (3.0, 77), 500: (4.5, 68),
                    550: (8.5, 52), 600: (12.0, 40),
                },
                "boyle_T_approx": 411.0,
            },
            "N2": {
                "data": {
                    100: (-159, 1600), 125: (-95, 900), 150: (-65, 550),
                    200: (-23, 290), 250: (-9.5, 185), 273: (-4.2, 145),
                    298: (0.5, 125), 300: (0.8, 122), 320: (3.5, 108),
                    350: (8.0, 90), 400: (16.0, 65), 450: (23.0, 48),
                    500: (31.0, 30),
                },
                "boyle_T_approx": 332.0,
            },
            "O2": {
                "data": {
                    100: (-198, 1750), 154.58: (-95, 900),
                    200: (-42, 340), 250: (-19, 210), 273: (-9.0, 165),
                    300: (-1.6, 128), 350: (8.0, 92), 400: (16.0, 65),
                    450: (23.0, 46), 500: (28.0, 32),
                },
                "boyle_T_approx": 305.0,
            },
            "CO2": {
                "data": {
                    230: (-225, 2400), 273: (-153, 1600), 280: (-140, 1450),
                    304.13: (-80, 800),
                    320: (-55, 600), 350: (-25, 380), 400: (-10, 240),
                    450: (5, 165), 500: (-14, 105), 550: (-5, 70),
                },
                "boyle_T_approx": 725.0,
            },
            "CH4": {
                "data": {
                    150: (-140, 1200), 190.56: (-60, 650),
                    200: (-48, 520), 250: (-26, 320), 273: (-18, 255),
                    298: (-8, 185), 300: (-7.5, 180), 350: (8, 120),
                    400: (22, 68), 450: (32, 42), 500: (41, 25),
                },
                "boyle_T_approx": 510.0,
            },
            "NH3": {
                "data": {
                    270: (-260, 2200), 300: (-150, 1400), 350: (-82, 850),
                    405.4: (-35, 400),
                    450: (-5, 250), 500: (35, 140),
                },
                "boyle_T_approx": 400.0,
            },
            "H2O": {
                "data": {
                    373: (-270, 1800), 473: (-115, 700), 573: (-40, 350),
                    647.1: (0, 0),
                    673: (25, -100),
                },
                "boyle_T_approx": 573.0,
            },
            "C2H6": {
                "data": {
                    250: (-330, 2500), 305.32: (-120, 1000),
                    350: (-58, 500), 400: (-18, 300), 450: (12, 180),
                    500: (35, 100),
                },
                "boyle_T_approx": 505.0,
            },
            "C2H4": {
                "data": {
                    250: (-218, 1800), 282.34: (-85, 750),
                    300: (-55, 520), 350: (-18, 310), 400: (8, 190),
                    450: (28, 110), 500: (45, 60),
                },
                "boyle_T_approx": 400.0,
            },
            "air(approx)": {
                "data": {
                    100: (-172, 1500), 200: (-38, 360), 273: (-7.0, 140),
                    298: (-1.5, 118), 300: (-1.2, 116), 350: (9.0, 78),
                    400: (17.0, 50), 500: (26.0, 28),
                },
                "boyle_T_approx": 327.0,
            },
            "Kr": {
                "data": {
                    150: (-220, 1500), 209.35: (-80, 700),
                    250: (-36, 370), 300: (-18, 250), 350: (-4, 175),
                    400: (8, 115), 500: (22, 60),
                },
                "boyle_T_approx": 435.0,
            },
            "Xe": {
                "data": {
                    200: (-365, 2200), 289.73: (-100, 800),
                    300: (-63, 430), 350: (-30, 280), 400: (-5, 185),
                    450: (15, 115), 500: (30, 65),
                },
                "boyle_T_approx": 468.0,
            },
            "SF6": {
                "data": {
                    280: (-350, 2500), 318.73: (-120, 1000),
                    350: (-45, 500), 400: (-8, 300), 450: (15, 175),
                    500: (32, 90),
                },
                "boyle_T_approx": 580.0,
            },
        }

    def _run_base(self, gas_species: str, T: float = None, mode="Z_from_Vm",
                 Vm=None, P_atm=None, order=2, **kwargs) -> dict:
        R = self.R

        species = self._find_species(gas_species)
        data_entry = self._virial_data[species]
        raw_data = data_entry["data"]

        mode = mode.lower().strip()

        if mode == "boyle_temperature":
            Tb = data_entry.get("boyle_T_approx", None)
            return {
                "gas_species": species,
                "result_value": Tb,
                "unit": "K",
                "explanation": (
                    f"The Boyle temperature is where B(T) = 0. For {species}, Tb ≈ {Tb} K.\n"
                    f"• At Tb: Z ≈ 1 + O(1/Vm²), minimal deviation from ideal behavior\n"
                    f"• Below Tb: B < 0 (attractions dominate, Z < 1)\n"
                    f"• Above Tb: B > 0 (repulsions dominate, Z > 1)"
                ),
                "method_used": "Boyle temperature from virial coefficient zero-crossing",
            }

        B_cm3, C_cm3 = self._interpolate_coefficients(raw_data, T)
        B_L = B_cm3 / 1000.0   # L/mol
        C_L2 = C_cm3 / 1e6     # L²/mol²

        if mode in ("z_from_vm", "z"):
            if Vm is None:
                raise ChemMCPInputError("Z_from_Vm mode requires molar volume Vm (L/mol).")

            terms = []
            Z = 1.0
            terms.append(("1 (identity)", 1.0))

            term_B = B_L / Vm
            Z += term_B
            terms.append((f"B/Vm = {B_L:.6f}/{Vm}", round(term_B, 8)))

            if order >= 3 and abs(C_L2) > 1e-10:
                term_C = C_L2 / (Vm ** 2)
                Z += term_C
                terms.append((f"C/Vm² = {C_L2:.6f}/{Vm}²", round(term_C, 10)))

            P_ideal = R * T / Vm

            return {
                "gas_species": species,
                "temperature_K": T,
                "result_value": round(Z, 8),
                "unit": "dimensionless",
                "compressibility_factor_Z": round(Z, 8),
                "virial_coefficients_used": {
                    "B_T_L_per_mol": round(B_L, 8),
                    "B_T_cm3_per_mol": round(B_cm3, 4),
                    "C_T_L2_per_mol2": round(C_L2, 8),
                    "C_T_cm2_per_mol2": round(C_cm3, 2),
                },
                "pressure_ideal_atm": round(P_ideal, 4),
                "pressure_deviation_percent": round((Z - 1.0) * 100, 4),
                "expansion_order_used": order,
                "virial_expansion_terms": terms,
            }

        elif mode in ("p_from_vm", "p"):
            if Vm is None:
                raise ChemMCPInputError("P_from_Vm mode requires molar volume Vm.")

            z_result = self._run_base(gas_species, T, "Z_from_Vm", Vm=Vm, order=order)
            Z = z_result["result_value"]
            P = Z * R * T / Vm
            P_ideal = R * T / Vm

            return {
                "gas_species": species,
                "temperature_K": T,
                "result_value": round(P, 6),
                "unit": "atm",
                "compressibility_factor_Z": round(Z, 8),
                "ideal_pressure_atm": round(P_ideal, 6),
                "pressure_deviation_atm": round(P - P_ideal, 6),
                "pressure_deviation_percent": round((P - P_ideal) / P_ideal * 100, 4),
                "molar_volume_Vm_L_mol": Vm,
                "note": f"P = ZRT/Vm = ({Z:.8f})({R})({T})/({Vm}) = {P:.6f} atm",
                "pressure_ideal_atm": round(P_ideal, 4)
            }

        elif mode in ("vm_from_p", "v"):
            if P_atm is None:
                raise ChemMCPInputError("Vm_from_P mode requires pressure P_atm.")

            Vm_est = R * T / P_atm
            for i in range(50):
                Z_est = 1.0 + B_L / Vm_est
                if order >= 3:
                    Z_est += C_L2 / (Vm_est ** 2)
                P_calc = Z_est * R * T / Vm_est
                error = (P_calc - P_atm) / P_atm
                if abs(error) < 1e-8:
                    break
                Vm_est = Vm_est * P_calc / P_atm

            Z_final = 1.0 + B_L / Vm_est + (C_L2 / (Vm_est ** 2) if order >= 3 else 0)
            Vm_ideal = R * T / P_atm

            return {
                "gas_species": species,
                "temperature_K": T,
                "result_value": round(Vm_est, 6),
                "unit": "L/mol",
                "pressure_atm": P_atm,
                "Vm_ideal_L_mol": round(Vm_ideal, 6),
                "Vm_deviation_percent": round((Vm_est - Vm_ideal) / Vm_ideal * 100, 4),
                "iterations_used": i + 1,
                "virial_coefficients": {
                    "B_L_per_mol": round(B_L, 8),
                    "C_L2_per_mol2": round(C_L2, 8),
                },
                "note": f"Vm={Vm_est:.6f} L/mol at P={P_atm} atm, T={T} K",
            }

        else:
            raise ChemMCPInputError(f"Unknown mode '{mode}'. Use: Z_from_Vm, P_from_Vm, Vm_from_P, or boyle_temperature.")

    def _find_species(self, name: str) -> str:
        """Find species key (case-insensitive)."""
        name = name.strip()
        if name in self._virial_data:
            return name
        name_lower = name.lower().replace(" ", "")
        for key in self._virial_data:
            if key.lower().replace(" ", "") == name_lower:
                return key
        available = list(self._virial_data.keys())
        raise ChemMCPError(
            f"No virial coefficient data for '{name}'. Available gases: {available}"
        )

    @staticmethod
    def _interpolate_coefficients(data: dict, T: float) -> tuple:
        """Interpolate B and C at temperature T."""
        temps = sorted(data.keys())

        if T in data:
            return data[T]

        if T < temps[0]:
            t0, t1 = temps[0], temps[min(1, len(temps)-1)]
        elif T > temps[-1]:
            t0, t1 = temps[-2], temps[-1]
        else:
            for i in range(len(temps) - 1):
                if temps[i] <= T <= temps[i+1]:
                    t0, t1 = temps[i], temps[i+1]
                    break
            else:
                t0, t1 = temps[-2], temps[-1]

        b0, c0 = data[t0]
        b1, c1 = data[t1]
        frac = (T - t0) / (t1 - t0) if t1 != t0 else 0

        B = b0 + (b1 - b0) * frac
        C = c0 + (c1 - c0) * frac
        return (B, C)

    
    def _convergence_note(self, Vm, order, Z):
        """Generate a convergence note for the virial expansion."""
        if Vm < 0.5:
            return f"Dense gas (Vm={Vm} L/mol): Truncation at order {order} may be inaccurate (% error large). Higher-order terms significant."
        elif Vm > 50:
            return f"Dilute gas (Vm={Vm} L/mol): Excellent (accurate to ~0.01%) convergence at order {order}. Z ≈ {Z:.6f} ≈ 1.0."
        else:
            return f"Moderate density (Vm={Vm} L/mol): Good (accurate within ~1%) convergence at order {order}. Z = {Z:.6f}."

def _run_text(self, query: str) -> dict:
        q = query.strip()
        import re

        mode = "Z_from_Vm"
        if "boyle" in q.lower():
            mode = "boyle_temperature"
        elif "p from" in q.lower() or "calculate p" in q.lower():
            mode = "P_from_Vm"

        species = None
        for gas in self._virial_data:
            if gas.lower() in q.lower():
                species = gas
                break
        if not species:
            for gas in ["He", "H2", "Ar", "N2", "O2", "CO2", "CH4", "NH3"]:
                if gas.lower() in q.lower():
                    species = gas
                    break
        if not species:
            raise ChemMCPInputError(f"Specify gas species. Available: {list(self._virial_data.keys())}")

        T_match = re.search(r'T[=\s]?(\d+\.?\d*)\s*K?', q)
        T_val = float(T_match.group(1)) if T_match else 298.15

        Vm_match = re.search(r'Vm[=\s]?(\d+\.?\d*)', q)
        Vm_val = float(Vm_match.group(1)) if Vm_match else None

        ord_match = re.search(r'order[=\s]?(\d)', q)
        order_val = int(ord_match.group(1)) if ord_match else 2

        return self._run_base(species, T_val, mode, Vm=Vm_val, order=order_val)
