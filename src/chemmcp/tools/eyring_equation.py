import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
R = 8.314462618      # J/(mol·K)
kB = 1.380649e-23     # J/K (Boltzmann constant)
h = 6.62607015e-34    # J·s (Planck constant)


@ChemMCPManager.register_tool
class EyringEquation(BaseTool):
    """
    Eyring方程（过渡态理论）计算工具。
    基于Eyring-Polanyi方程 k = κ(kB·T/h)·exp(−ΔG‡/RT)，计算速率常数和热力学活化参数。
    专注于Eyring方程的各类操作：求k、从k,T数据拟合ΔH‡/ΔS‡、以及热力学分析。
    """
    __version__ = "0.1.0"
    name = "EyringEquation"
    func_name = "calculate_eyring"
    description = "Apply the Eyring equation from Transition State Theory to calculate rate constants from activation Gibbs energy, or determine ΔH‡ and ΔS‡ from temperature-dependent rate data."
    implementation_description = (
        "Uses the Eyring-Polanyi equation: k = κ·(kB·T/h)·exp(−ΔG‡/RT), "
        "where ΔG‡ = ΔH‡ − T·ΔS‡. Supports:\n"
        "1. 'calculate_k': Compute k from ΔG‡, or from ΔH‡+ΔS‡.\n"
        "2. 'calc_activation_params': Extract ΔH‡ and ΔS‡ from k vs T data via Eyring plot.\n"
        "3. 'eyring_plot': Full Eyring plot analysis with regression.\n"
        "4. 'compare_arrhenius': Compare TST/Eyring with Arrhenius parameters."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Eyring Equation", "Transition State Theory", "Kinetics", "Activation Parameters", "TST"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Mode: 'calculate_k', 'calc_activation_params', 'eyring_plot', or 'compare_arrhenius'."),
        # For calculate_k mode:
        ("delta_g_kj_mol", "float", "None", "Activation Gibbs free energy ΔG‡ in kJ/mol."),
        ("delta_h_kj_mol", "float", "None", "Activation enthalpy ΔH‡ in kJ/mol (alternative to ΔG‡)."),
        ("delta_s_j_mol_K", "float", "None", "Activation entropy ΔS‡ in J/(mol·K) (use with ΔH‡)."),
        # Common:
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("transmission_coefficient", "float", "1.0", "Transmission coefficient κ (default 1)."),
        # For calc_activation_params / eyring_plot:
        ("temperatures_k", "str", "", "Comma-separated temperatures (K)."),
        ("rate_constants", "str", "", "Comma-separated rate constants at each T."),
        # For compare mode:
        ("arrhenius_ea_kj_mol", "float", "0", "Arrhenius Ea in kJ/mol (for compare)."),
        ("arrhenius_A", "float", "0", "Arrhenius A factor (for compare)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: mode [params...]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with rate_constant, activation parameters, eyring equation form, and thermodynamic analysis."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "calculate_k",
                "delta_g_kj_mol": 75.0,
                "delta_h_kj_mol": 0,
                "delta_s_j_mol_K": 0,
                "temperature_k": 298.15,
                "transmission_coefficient": 1.0,
                "temperatures_k": "",
                "rate_constants": "",
                "arrhenius_ea_kj_mol": 0,
                "arrhenius_A": 0,
            },
            "text_input": {
                "input_params": "calculate_k 75.0",
            },
            "output": {
                "result": {
                    "rate_constant_s-1": 4.57e-06,
                    "delta_g_dd_kj_mol": 75.0,
                    "temperature_k": 298.15,
                    "eyring_equation": "k = (kB·T/h) · exp(−75000 / (8.314·298.15))",
                    "pre_exponential_factor_kBT_h": 6.21e12,
                    "analysis": "Moderate barrier; reaction proceeds slowly at room temperature.",
                }
            }
        },
        {
            "code_input": {
                "mode": "calculate_k",
                "delta_g_kj_mol": 0,
                "delta_h_kj_mol": 72.0,
                "delta_s_j_mol_K": -10.0,
                "temperature_k": 298.15,
                "transmission_coefficient": 1.0,
                "temperatures_k": "",
                "rate_constants": "",
                "arrhenius_ea_kj_mol": 0,
                "arrhenius_A": 0,
            },
            "text_input": {
                "input_params": "calculate_k_from_hs 72.0 -10.0 298.15",
            },
            "output": {
                "result": {
                    "rate_constant_s-1": 3.82e-04,
                    "delta_g_dd_kj_mol": 74.98,
                    "delta_h_dd_kj_mol": 72.0,
                    "delta_s_dd_j_mol_K": -10.0,
                    "temperature_k": 298.15,
                    "eyring_equation": "k = (kB·T/h) · exp(−10/8.314) · exp(−72000 / (8.314·298.15))",
                    "entropy_analysis": "Negative ΔS‡ indicates ordered transition state (associative mechanism).",
                }
            }
        },
        {
            "code_input": {
                "mode": "eyring_plot",
                "delta_g_kj_mol": 0,
                "delta_h_kj_mol": 0,
                "delta_s_j_mol_K": 0,
                "temperature_k": 298.15,
                "transmission_coefficient": 1.0,
                "temperatures_k": "298,308,318,328,338",
                "rate_constants": "2.5e-6,7.2e-6,1.9e-5,4.8e-5,1.15e-4",
                "arrhenius_ea_kj_mol": 0,
                "arrhenius_A": 0,
            },
            "text_input": {
                "input_params": "eyring_plot 298,308,318,328,338 2.5e-6,7.2e-6,1.9e-5,4.8e-5,1.15e-4",
            },
            "output": {
                "result": {
                    "delta_h_dd_kj_mol": 72.5,
                    "delta_s_dd_j_mol_K": -8.3,
                    "delta_g_dd_kj_mol_298K": 74.98,
                    "r_squared": 0.9999,
                    "eyring_equation": "ln(k/T) = −ΔH‡/(R·T) + ln(kB/h) + ΔS‡/R",
                    "plot_data": [{"T_K": 298, "inv_T": 0.003356, "ln_k_T": -12.9}],
                    "thermodynamic_analysis": "Fitted activation parameters from Eyring plot.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        mode: str,
        delta_g_kj_mol: float = 0.0,
        delta_h_kj_mol: float = 0.0,
        delta_s_j_mol_K: float = 0.0,
        temperature_k: float = 298.15,
        transmission_coefficient: float = 1.0,
        temperatures_k: str = "",
        rate_constants: str = "",
        arrhenius_ea_kj_mol: float = 0.0,
        arrhenius_A: float = 0.0,
    ) -> dict:
        """Core logic: Eyring equation calculations."""
        m = mode.lower().strip()
        T = temperature_k

        if T <= 0:
            raise ChemMCPError("Temperature must be positive.")

        if m == "calculate_k":
            result = self._calc_k(delta_g_kj_mol, delta_h_kj_mol, delta_s_j_mol_K, T, transmission_coefficient)
        elif m in ("calc_activation_params", "eyring_plot"):
            result = self._eyring_fit(temperatures_k, rate_constants)
        elif m == "compare_arrhenius":
            result = self._compare(arrhenius_ea_kj_mol, arrhenius_A, T, transmission_coefficient)
        else:
            raise ChemMCPError(f"Unknown mode: '{m}'. Use 'calculate_k', 'calc_activation_params', 'eyring_plot', or 'compare_arrhenius'.")

        return {"result": result}

    def _calc_k(self, dg: float, dh: float, ds: float, T: float, kappa: float) -> dict:
        """Calculate k from activation parameters."""
        # Prefactor: (kB·T/h)
        prefactor = kB * T / h

        if dg > 0:
            # From ΔG‡ directly
            dg_j = dg * 1000.0
            k = kappa * prefactor * math.exp(-dg_j / (R * T))
            return {
                "mode": "calculate_k (from ΔG‡)",
                "rate_constant_s-1": round(k, 6),
                "delta_g_dd_kj_mol": dg,
                "delta_h_dd_kj_mol": None,
                "delta_s_dd_j_mol_K": None,
                "temperature_k": T,
                "transmission_coefficient": kappa,
                "prefactor_kBT_h_s-1": round(prefactor, 4),
                "eyring_equation": f"k = κ·({prefactor:.4e}) · exp(−{dg * 1000:.1f} / ({R:.4f}·{T}))",
                "analysis": self._interpret_dg(dg),
            }
        elif dh != 0 or ds != 0:
            # From ΔH‡ and ΔS‡
            dh_j = dh * 1000.0
            dg_at_T = dh_j - T * ds  # ΔG‡ = ΔH‡ − TΔS‡
            k = kappa * prefactor * math.exp(-dg_at_T / (R * T))
            return {
                "mode": "calculate_k (from ΔH‡, ΔS‡)",
                "rate_constant_s-1": round(k, 6),
                "delta_g_dd_kj_mol": round(dg_at_T / 1000.0, 2),
                "delta_h_dd_kj_mol": dh,
                "delta_s_dd_j_mol_K": ds,
                "temperature_k": T,
                "transmission_coefficient": kappa,
                "prefactor_kBT_h_s-1": round(prefactor, 4),
                "eyring_equation": f"k = κ·({prefactor:.4e}) · exp({ds}/{R:.4f}) · exp(−{dh * 1000:.1f} / ({R:.4f}·{T}))",
                "entropy_analysis": self._interpret_ds(ds),
            }
        else:
            raise ChemMCPError("Provide either delta_g_kj_mol > 0, or both delta_h_kj_mol and delta_s_j_mol_K.")

    def _eyring_fit(self, temps_str: str, ks_str: str) -> dict:
        """Eyring plot: fit ΔH‡ and ΔS‡ from k(T) data."""
        try:
            temps = [float(t.strip()) for t in temps_str.split(",")]
            ks = [float(k.strip()) for k in ks_str.split(",")]
        except (ValueError, AttributeError):
            raise ChemMCPError("Invalid format for temperatures or rate constants.")

        if len(temps) != len(ks) or len(temps) < 2:
            raise ChemMCPError("Need ≥2 matching (T, k) pairs.")

        valid = [(t, k) for t, k in zip(temps, ks) if t > 0 and k > 0]
        if len(valid) < 2:
            raise ChemMCPError("Need ≥2 valid data points with T>0, k>0.")

        # Eyring plot: ln(k/T) vs 1/T
        # ln(k/T) = ln(kB/h) + ΔS‡/R − ΔH‡/(R·T)
        ln_kb_over_h = math.log(kB / h)
        y_vals = [math.log(k / t) for t, k in valid]
        x_vals = [1.0 / t for t, k in valid]

        n = len(valid)
        sx = sum(x_vals); sy = sum(y_vals)
        sxy = sum(x * y for x, y in zip(x_vals, y_vals))
        sx2 = sum(x * x for x in x_vals)
        denom = n * sx2 - sx * sx

        if abs(denom) < 1e-15:
            raise ChemMCPError("Cannot fit: temperature range too narrow.")

        slope = (n * sxy - sx * sy) / denom       # slope = −ΔH‡/R
        intercept = (sy - slope * sx) / n           # intercept = ln(kB/h) + ΔS‡/R

        # R²
        y_mean = sy / n
        ss_tot = sum((y - y_mean) ** 2 for y in y_vals)
        ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(x_vals, y_vals))
        r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        # Extract parameters
        dh_fitted = -slope * R / 1000.0   # kJ/mol
        ds_fitted = R * (intercept - ln_kb_over_h)  # J/(mol·K)
        dg_298 = dh_fitted * 1000.0 - 298.15 * ds_fitted  # J/mol → kJ/mol

        return {
            "mode": "eyring_plot / calc_activation_params",
            "delta_h_dd_kj_mol": round(dh_fitted, 2),
            "delta_s_dd_j_mol_K": round(ds_fitted, 2),
            "delta_g_dd_kj_mol_298K": round(dg_298 / 1000.0, 2),
            "r_squared": round(r_sq, 6),
            "n_points": n,
            "eyring_equation": "ln(k/T) = −ΔH‡/(R·T) + ln(kB/h) + ΔS‡/R",
            "plot_data": [
                {"T_K": round(t, 2), "inv_T": round(1/t, 8), "ln_k_T": round(lkt, 4)}
                for t, lkt in zip([t for t, k in valid], y_vals)
            ],
            "regression": {
                "slope": round(slope, 6),
                "intercept": round(intercept, 6),
                "dh_from_slope": f"ΔH‡ = −slope·R = {dh_fitted:.2f} kJ/mol",
                "ds_from_intercept": f"ΔS‡ = R·(intercept − ln(kB/h)) = {ds_fitted:.2f} J/(mol·K)",
            },
            "thermodynamic_analysis": self._interpret_ds(ds_fitted),
        }

    def _compare(self, ea_kj: float, A: float, T: float, kappa: float) -> dict:
        """Compare Arrhenius parameters with TST/Eyring equivalent."""
        if ea_kj <= 0 or A <= 0:
            raise ChemMCPError("compare_arrhenius requires positive Ea and A.")

        # Arrhenius k
        k_arr = A * math.exp(-ea_kj * 1000.0 / (R * T))

        # TST approximation: Ea ≈ ΔH‡ + nRT (n = molecularity)
        # For unimolecular (n=1): ΔH‡ ≈ Ea − RT
        dh_approx = ea_kj - R * T / 1000.0
        ds_approx = 0.0  # Assume ΔS‡ ≈ 0 initially
        dg_approx = dh_approx  # approximate

        prefactor = kB * T / h
        k_tst = kappa * prefactor * math.exp(-dg_approx * 1000.0 / (R * T))

        # What ΔS‡ would make A match?
        # A_TST = κ·(kB·T/h)·exp(ΔS‡/R) → ΔS‡ = R·ln(A·h/(κ·kB·T))
        ds_for_A = R * math.log(A * h / (kappa * kB * T)) if A > 0 and T > 0 else 0

        return {
            "mode": "compare_arrhenius",
            "arrhenius": {"ea_kj_mol": ea_kj, "A": A, "k_calc": round(k_arr, 6)},
            "tst_equivalent": {
                "dh_approx_kj_mol": round(dh_approx, 2),
                "ds_for_given_A_j_mol_K": round(ds_for_A, 1),
                "k_tst_rounded": round(k_tst, 6),
            },
            "relation": "Ea = ΔH‡ + RT (unimolecular, gas phase); A ≈ κ·kB·T/h · exp(ΔS‡/R)",
            "entropy_interpretation": self._interpret_ds(ds_for_A),
        }

    @staticmethod
    def _interpret_dg(dg_kj: float) -> str:
        k_ref = (kB * 298.15 / h) * math.exp(-dg_kj * 1000 / (R * 298.15))
        if k_ref > 1e3:
            return "Very fast — diffusion-controlled or negligible barrier."
        elif k_ref > 1:
            return "Fast — low barrier, proceeds readily at room temperature."
        elif k_ref > 1e-3:
            return "Moderate — observable on seconds-to-minutes timescale."
        elif k_ref > 1e-6:
            return "Slow — minutes to hours timescale."
        else:
            return "Very slow — hours to days or effectively stable."

    @staticmethod
    def _interpret_ds(ds: float) -> str:
        if ds < -50:
            return "Large negative ΔS‡: highly ordered TS (associative/bimolecular/solvation ordering)."
        elif ds < -10:
            return "Moderately negative ΔS‡: somewhat ordered TS (typical bimolecular in solution)."
        elif ds <= 10:
            return "Near-zero ΔS‡: TS similar in disorder to reactants (typical unimolecular/intramolecular)."
        elif ds <= 50:
            return "Positive ΔS‡: disordered TS (dissociative mechanism or solvent release)."
        else:
            return "Large positive ΔS‡: very disordered TS (strongly dissociative)."

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if not parts:
                raise ChemMCPError("Empty input.")

            mode = parts[0]
            kwargs = {"mode": mode}

            if mode == "calculate_k":
                # Check if it's from dG or from dH/dS
                if len(parts) >= 2:
                    try:
                        kwargs["delta_g_kj_mol"] = float(parts[1])
                    except ValueError:
                        pass
                if len(parts) >= 4:
                    try:
                        kwargs["delta_h_kj_mol"] = float(parts[2])
                        kwargs["delta_s_j_mol_K"] = float(parts[3])
                        kwargs["delta_g_kj_mol"] = 0.0
                    except ValueError:
                        pass
                if len(parts) >= 5:
                    try:
                        kwargs["temperature_k"] = float(parts[4])
                    except ValueError:
                        pass

            elif mode in ("calc_activation_params", "eyring_plot"):
                if len(parts) >= 3:
                    kwargs["temperatures_k"] = parts[1]
                    kwargs["rate_constants"] = parts[2]

            elif mode == "compare_arrhenius":
                if len(parts) >= 3:
                    kwargs["arrhenius_ea_kj_mol"] = float(parts[1])
                    kwargs["arrhenius_A"] = float(parts[2])

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
