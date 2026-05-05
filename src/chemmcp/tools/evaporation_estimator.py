"""
Evaporation Estimator — 蒸发浓缩时间和条件估算
基于蒸发速率模型估算浓缩时间，推荐最佳蒸发条件
"""
import logging
import math
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见溶剂物理常数 ────────────────────────────────────────────────
SOLVENT_DATA: Dict[str, dict] = {
    "water": {
        "bp_C": 100.0, "enthalpy_vaporization_kJ_mol": 40.65,
        "density_g_mL": 1.0, "surface_tension_mN_m": 72.8,
        "vapor_pressure_kPa_25C": 3.17,
        "viscosity_cP_25C": 0.89,
        "safety": "Non-flammable; safe for rotary evaporation.",
    },
    "methanol": {
        "bp_C": 64.7, "enthalpy_vaporization_kJ_mol": 35.21,
        "density_g_mL": 0.791, "surface_tension_mN_m": 22.6,
        "vapor_pressure_kPa_25C": 16.9,
        "viscosity_cP_25C": 0.54,
        "safety": "FLAMMABLE; TOXIC; use cold trap with liquid N₂.",
    },
    "ethanol": {
        "bp_C": 78.4, "enthalpy_vaporization_kJ_mol": 38.56,
        "density_g_mL": 0.789, "surface_tension_mN_m": 22.1,
        "vapor_pressure_kPa_25C": 7.87,
        "viscosity_cP_25C": 1.07,
        "safety": "FLAMMABLE; use proper ventilation.",
    },
    "acetone": {
        "bp_C": 56.0, "enthalpy_vaporization_kJ_mol": 29.1,
        "density_g_mL": 0.784, "surface_tension_mN_m": 23.7,
        "vapor_pressure_kPa_25C": 30.6,
        "viscosity_cP_25C": 0.30,
        "safety": "HIGHLY FLAMMABLE; evaporates very quickly.",
    },
    "dichloromethane": {
        "bp_C": 39.8, "enthalpy_vaporization_kJ_mol": 28.0,
        "density_g_mL": 1.33, "surface_tension_mN_m": 28.1,
        "vapor_pressure_kPa_25C": 47.3,
        "viscosity_cP_25C": 0.41,
        "safety": "TOXIC (suspected carcinogen); use fume hood; cold trap required.",
    },
    "ethyl_acetate": {
        "bp_C": 77.1, "enthalpy_vaporization_kJ_mol": 32.0,
        "density_g_mL": 0.902, "surface_tension_mN_m": 23.9,
        "vapor_pressure_kPa_25C": 12.6,
        "viscosity_cP_25C": 0.43,
        "safety": "FLAMMABLE; irritant.",
    },
    "hexane": {
        "bp_C": 68.7, "enthalpy_vaporization_kJ_mol": 28.85,
        "density_g_mL": 0.655, "surface_tension_mN_m": 18.4,
        "vapor_pressure_kPa_25C": 17.0,
        "viscosity_cP_25C": 0.30,
        "safety": "FLAMMABLE; neurotoxic; use in fume hood.",
    },
    "toluene": {
        "bp_C": 110.6, "enthalpy_vaporization_kJ_mol": 33.2,
        "density_g_mL": 0.867, "surface_tension_mN_m": 27.9,
        "vapor_pressure_kPa_25C": 3.80,
        "viscosity_cP_25C": 0.55,
        "safety": "FLAMMABLE; TOXIC (reproductive hazard); use fume hood.",
    },
    "chloroform": {
        "bp_C": 61.2, "enthalpy_vaporization_kJ_mol": 29.24,
        "density_g_mL": 1.48, "surface_tension_mN_m": 26.5,
        "vapor_pressure_kPa_25C": 26.0,
        "viscosity_cP_25C": 0.54,
        "safety": "TOXIC (carcinogen); decomposes to phosgene; avoid light/heat.",
    },
    "thf": {
        "bp_C": 66.0, "enthalpy_vaporization_kJ_mol": 29.3,
        "density_g_mL": 0.889, "surface_tension_mN_m": 26.4,
        "vapor_pressure_kPa_25C": 19.0,
        "viscosity_cP_25C": 0.46,
        "safety": "FORMS PEROXIDES over time; FLAMMABLE; check for peroxides before evaporation.",
    },
    "dmf": {
        "bp_C": 153.0, "enthalpy_vaporization_kJ_mol": 47.0,
        "density_g_mL": 0.944, "surface_tension_mN_m": 37.1,
        "vapor_pressure_kPa_25C": 0.36,
        "viscosity_cP_25C": 0.80,
        "safety": "High boiling — needs high vacuum; reproductive toxin.",
    },
    "dmso": {
        "bp_C": 189.0, "enthalpy_vaporization_kJ_mol": 52.8,
        "density_g_mL": 1.099, "surface_tension_mN_m": 43.0,
        "vapor_pressure_kPa_25C": 0.06,
        "viscosity_cP_25C": 1.99,
        "safety": "Very high boiling; difficult to remove by rotavap; consider lyophilization or high-vacuum pump.",
    },
}


@ChemMCPManager.register_tool
class EvaporationEstimator(BaseTool):
    """
    蒸发浓缩估算器：根据溶剂、体积、温度、压力条件估算蒸发时间，
    推荐旋转蒸发仪等设备的最佳操作参数。
    """
    __version__ = "0.1.0"
    name = "EvaporationEstimator"
    func_name = "estimate_evaporation"
    description = "Estimate evaporation/concentration time and recommend optimal conditions for rotary evaporation, nitrogen blow-down, or vacuum concentration."
    implementation_description = "Uses solvent physical properties (boiling point, enthalpy of vaporization, vapor pressure) combined with operating conditions (temperature, pressure, surface area) to estimate evaporation rate and time via a semi-empirical model."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Evaporation", "Rotary Evaporator", "Concentration", "Sample Preparation", "Solvent"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "estimate", "Mode: 'estimate', 'solvent_info', 'condition_recommendation'."),
        ("initial_volume_mL", "float", "100", "Initial solution volume in mL."),
        ("target_volume_mL", "float", "1", "Target/final volume in mL (for concentration)."),
        ("solvent", "str", "water", "Solvent name (key from built-in database)."),
        ("temperature_C", "float", "40", "Bath temperature in °C (rotavap water bath)."),
        ("pressure_mbar", "float", "200", "System pressure in mbar (vacuum)."),
        ("flask_size_mL", "float", "500", "Rotary evaporator flask size in mL."),
        ("rotation_rpm", "float", "100", "Rotation speed in RPM."),
        ("coolant_temp_C", "float", "-10", "Condenser coolant temperature in °C."),
        ("method", "str", "rotavap", "Method: 'rotavap', 'n2_blow', 'vacuum_centrifuge', 'air_dry'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., 'estimate 500 50 methanol 40 200' or 'solvent_info acetone'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with estimated time, recommended conditions, safety notes, and step-by-step protocol."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "estimate",
                "initial_volume_mL": 500,
                "target_volume_mL": 50,
                "solvent": "methanol",
                "temperature_C": 40,
                "pressure_mbar": 200,
                "flask_size_mL": 1000,
                "rotation_rpm": 100,
                "coolant_temp_C": -10,
                "method": "rotavap",
            },
            "text_input": {
                "input_params": "estimate 500 50 methanol 40 200",
            },
            "output": {
                "result": {
                    "mode": "estimate",
                    "note": "Estimated time based on semi-empirical model.",
                }
            }
        },
        {
            "code_input": {
                "mode": "solvent_info",
                "initial_volume_mL": 0,
                "target_volume_mL": 0,
                "solvent": "dichloromethane",
                "temperature_C": 35,
                "pressure_mbar": 300,
                "flask_size_mL": 500,
                "rotation_rpm": 80,
                "coolant_temp_C": 0,
                "method": "rotavap",
            },
            "text_input": {
                "input_params": "solvent_info dichloromethane",
            },
            "output": {
                "result": {
                    "mode": "solvent_info",
                    "solvent": "dichloromethane",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _get_solvent(self, name: str) -> dict:
        key = name.lower().strip().replace(" ", "_")
        if key not in SOLVENT_DATA:
            available = ", ".join(sorted(SOLVENT_DATA.keys()))
            raise ChemMCPError(
                f"Unknown solvent '{name}'. Available solvents: {available}"
            )
        return SOLVENT_DATA[key]

    def _estimate_boiling_point(self, bp_normal: float, pressure_mbar: float) -> float:
        """
        Estimate boiling point at reduced pressure using Clausius-Clapeyron approximation.
        Simplified: T₂ ≈ T₁ / (1 - T₁·R·ln(P₂/P₁)/ΔHvap)
        """
        if pressure_mbar >= 1013.25:
            return bp_normal
        if pressure_mbar <= 0:
            return -100.0  # essentially no boiling

        T1 = bp_normal + 273.15  # K
        P1 = 1013.25  # mbar (1 atm)
        P2 = max(pressure_mbar, 0.1)
        dH = 45000  # J/mol typical value

        try:
            R_gas = 8.314
            ln_ratio = math.log(P2 / P1)
            denom = 1 - (T1 * R_gas * ln_ratio) / dH
            if denom <= 0:
                return bp_normal * 0.4  # rough lower bound
            T2_K = T1 / denom
            return min(T2_K - 273.15, bp_normal)  # shouldn't exceed normal BP
        except (ValueError, ZeroDivisionError):
            return bp_normal * 0.5

    def _calc_evap_rate(self, solvent: dict, temp_C: float, pressure_mbar: float,
                        flask_mL: float, rotation_rpm: float) -> float:
        """
        Semi-empirical evaporation rate estimation (mL/min).
        Base rate scaled by: ΔT (superheat), vacuum factor, surface area, rotation.
        """
        bp_est = self._estimate_boiling_point(solvent["bp_C"], pressure_mbar)
        delta_T = temp_C - bp_est

        if delta_T < -20:
            return 0.01  # minimal evaporation if too cold
        elif delta_T < 0:
            delta_T = 1  # some evaporation even below estimated BP

        # Base rate: higher vapor pressure → faster evaporation
        vp = solvent.get("vapor_pressure_kPa_25C", 10)

        # Vacuum enhancement factor
        vacuum_factor = max(0.1, 1013.25 / max(pressure_mbar, 1))

        # Surface area factor (larger flask → thinner film at same volume)
        area_factor = (flask_mL / 500) ** 0.5

        # Rotation factor (more rpm → better film renewal)
        rotation_factor = 0.5 + (rotation_rpm / 200) ** 0.5

        # Temperature driving force
        temp_factor = max(0, 1 + delta_T / 30)

        # Base evaporation rate coefficient (empirically calibrated)
        base_rate = vp * 0.02  # mL/min at reference conditions

        rate = base_rate * vacuum_factor * area_factor * rotation_factor * temp_factor

        # Cap at physically reasonable maximum
        return min(max(rate, 0.001), 500)

    def _run_base(self, mode: str = "estimate", initial_volume_mL: float = 100.0,
                  target_volume_mL: float = 1.0, solvent: str = "water",
                  temperature_C: float = 40.0, pressure_mbar: float = 200.0,
                  flask_size_mL: float = 500.0, rotation_rpm: float = 100.0,
                  coolant_temp_C: float = -10.0, method: str = "rotavap") -> dict:

        m = mode.lower().strip()

        if m == "solvent_info":
            sdata = self._get_solvent(solvent)
            return {"result": {
                "mode": "solvent_info",
                "solvent": solvent,
                **sdata,
                "estimated_bp_at_pressure": self._estimate_boiling_point(
                    sdata["bp_C"], pressure_mbar),
                "recommended_method": self._recommend_method(sdata),
            }}

        elif m == "condition_recommendation":
            sdata = self._get_solvent(solvent)
            return {"result": {
                "mode": "condition_recommendation",
                "solvent": solvent,
                **self._recommend_conditions(sdata, initial_volume_mL, target_volume_mL),
            }}

        elif m == "estimate":
            sdata = self._get_solvent(solvent)
            bp_est = self._estimate_boiling_point(sdata["bp_C"], pressure_mbar)

            volume_to_remove = max(0, initial_volume_mL - target_volume_mL)
            if volume_to_remove <= 0:
                raise ChemMCPError(
                    f"Target volume ({target_volume_mL} mL) must be less than "
                    f"initial volume ({initial_volume_mL} mL)."
                )

            # Use average rate (rate decreases as volume decreases)
            rate_initial = self._calc_evap_rate(sdata, temperature_C, pressure_mbar,
                                                 flask_size_mL, rotation_rpm)

            # Approximate: use geometric mean of initial and final rates
            vol_mid = (initial_volume_mL + target_volume_mL) / 2
            # Simulate decreasing rate with volume
            avg_rate = rate_initial * 0.6  # empirical reduction factor

            est_time_min = volume_to_remove / avg_rate if avg_rate > 0 else float('inf')

            condenser_ok = coolant_temp_C <= (bp_est - 10)

            return {"result": {
                "mode": "estimate",
                "method": method,
                "solvent": solvent,
                "initial_volume_mL": initial_volume_mL,
                "target_volume_mL": target_volume_mL,
                "volume_to_remove_mL": round(volume_to_remove, 2),
                "concentration_factor": round(initial_volume_mL / target_volume_mL, 1),
                "normal_boiling_point_C": sdata["bp_C"],
                "estimated_boiling_point_at_pressure_C": round(bp_est, 1),
                "bath_temperature_C": temperature_C,
                "pressure_mbar": pressure_mbar,
                "delta_T_superheat_C": round(temperature_C - bp_est, 1),
                "estimated_evap_rate_mL_min_initial": round(rate_initial, 3),
                "estimated_evap_rate_mL_min_avg": round(avg_rate, 3),
                "estimated_time_min": round(est_time_min, 1),
                "estimated_time_hr": round(est_time_min / 60, 2),
                "condenser_efficiency": "OK" if condenser_ok else "WARNING: Coolant too warm!",
                "safety_note": sdata.get("safety", ""),
                "recommendations": self._build_recommendations(
                    sdata, temperature_C, pressure_mbar, bp_est, est_time_min),
                "protocol_steps": self._generate_protocol(
                    solvent, sdata, initial_volume_mL, target_volume_mL,
                    temperature_C, pressure_mbar, est_time_min),
            }}
        else:
            raise ChemMCPError(
                f"Unknown mode '{mode}'. Use: 'estimate', 'solvent_info', "
                "or 'condition_recommendation'."
            )

    @staticmethod
    def _recommend_method(sdata: dict) -> str:
        bp = sdata["bp_C"]
        vp = sdata.get("vapor_pressure_kPa_25C", 10)
        if bp > 150:
            return "High-vacuum rotary evaporation or lyophilization"
        elif bp > 100:
            return "Rotary evaporation with good vacuum (<100 mbar)"
        elif vp > 25:
            return "Rotary evaporation (easy removal) or N₂ blow-down"
        else:
            return "Standard rotary evaporation"

    def _recommend_conditions(self, sdata: dict, V_i: float, V_f: float) -> dict:
        bp = sdata["bp_C"]
        rec_temp = min(bp - 10, 45) if bp < 80 else min(bp - 20, 60)
        rec_pressure = max(200, int(bp * 2)) if bp < 80 else max(100, int(bp))
        return {
            "recommended_bath_temp_C": max(rec_temp, 25),
            "recommended_pressure_mbar": rec_pressure,
            "recommended_rotation_rpm": 80,
            "recommended_coolant": "ice/water" if bp < 60 else "ice/salt (-5°C)" if bp < 90 else "dry ice/acetone",
            "estimated_concentration_ratio": round(V_i / max(V_f, 0.1), 1),
        }

    def _build_recommendations(self, sdata: dict, temp: float, pressure: float,
                                bp_est: float, est_time: float) -> List[str]:
        recs = []
        delta_T = temp - bp_est
        if delta_T > 30:
            recs.append(f"⚠ Bath temp much higher than BP at set pressure (ΔT={delta_T:.0f}°C). Risk of bumping.")
        elif delta_T < 5:
            recs.append(f"ΔT only {delta_T:.1f}°C — evaporation will be slow. Consider lowering pressure or raising bath temp slightly.")
        else:
            recs.append(f"✓ Good ΔT (~{delta_T:.0f}°C) for steady evaporation.")

        if pressure > 500:
            recs.append("Consider stronger vacuum for faster evaporation.")
        if est_time > 120:
            recs.append(f"Long run expected (~{est_time/60:.0f}h). Plan accordingly.")
        if est_time < 2:
            recs.append("Fast evaporation expected — monitor closely to avoid drying completely.")

        # Solvent-specific
        if sdata["bp_C"] < 40:
            recs.append("Low-boiling solvent: ensure excellent condenser cooling (dry ice/acetone recommended).")
        if "peroxides" in sdata.get("safety", "").lower():
            recs.append("⚠ Check this solvent for peroxide formation before concentrating!")
        return recs

    def _generate_protocol(self, solvent: str, sdata: dict, V_i: float, V_f: float,
                           temp: float, pressure: float, est_time: float) -> List[str]:
        steps = [
            f"1. Assemble rotary evaporator with appropriate size flask (≥{max(V_i*2, 250):.0f} mL).",
            f"2. Pour {V_i:.0f} mL of {solvent} solution into flask (fill ≤½ capacity).",
            f"3. Attach flask to rotavap bump trap; engage clip securely.",
            f"4. Set water bath to {temp}°C; turn on circulation.",
            f"5. Start rotation at ~80-100 RPM.",
            f"6. Gradually apply vacuum to {pressure} mbar (watch for bumping).",
            f"7. Monitor evaporation; estimated time: {est_time:.0f} min.",
            f"8. When volume reaches ~{V_f:.0f} mL, release vacuum and stop rotation.",
            f"9. If needed, flush with N₂ before removing flask.",
            f"10. Clean flask promptly (especially for {solvent}).",
        ]
        return steps

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            mode = parts[0]
            if mode == "estimate":
                Vi = float(parts[1]) if len(parts) > 1 else 100
                Vf = float(parts[2]) if len(parts) > 2 else 1
                sol = parts[3] if len(parts) > 3 else "water"
                T = float(parts[4]) if len(parts) > 4 else 40
                P = float(parts[5]) if len(parts) > 5 else 200
                return self._run_base(mode, initial_volume_mL=Vi, target_volume_mL=Vf,
                                       solvent=sol, temperature_C=T, pressure_mbar=P)
            elif mode == "solvent_info":
                sol = parts[1] if len(parts) > 1 else "water"
                return self._run_base(mode, solvent=sol)
            elif mode == "condition_recommendation":
                sol = parts[1] if len(parts) > 1 else "water"
                Vi = float(parts[2]) if len(parts) > 2 else 100
                Vf = float(parts[3]) if len(parts) > 3 else 1
                return self._run_base(mode, solvent=sol, initial_volume_mL=Vi,
                                       target_volume_mL=Vf)
            else:
                raise ValueError(f"Unknown text mode: {mode}")
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
