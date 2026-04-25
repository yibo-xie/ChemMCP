import logging
import math

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TransitionStateTheory(BaseTool):
    """
    过渡态理论（TST/ACT）计算速率常数工具。
    基于Eyring方程，从活化吉布斯自由能、活化焓和活化熵计算速率常数。
    支持不同温度下的k值计算和作图数据。
    """
    __version__ = "0.1.0"
    name = "TransitionStateTheory"
    func_name = "calculate_tst_rate_constant"
    description = "Calculate rate constants using Transition State Theory (Eyring equation) from activation parameters (ΔG‡, ΔH‡, ΔS‡)."
    implementation_description = "Uses Eyring-Polanyi equation: k = (kB·T/h)·exp(−ΔG‡/RT) = κ·(kB·T/h)·exp(ΔS‡/R)·exp(−ΔH‡/RT), where κ is transmission coefficient (default 1)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Kinetics", "Transition State Theory", "Eyring Equation", "Activation Parameters", "TST"]
    required_envs    = []

    code_input_sig = [
        ("calculation_mode", "str", "N/A", "Mode: 'from_gibbs' (from ΔG‡), 'from_enthalpy_entropy' (from ΔH‡ and ΔS‡), 'eyring_plot' (fit ΔH‡ and ΔS‡ from k vs T data), or 'compare' (compare TST with Arrhenius)."),
        # For from_gibbs mode:
        ("delta_g_double_dagger", "float", "N/A", "Activation Gibbs free energy ΔG‡ in kJ/mol (for from_gibbs mode)."),
        # For from_enthalpy_entropy mode:
        ("delta_h_double_dagger", "float", "N/A", "Activation enthalpy ΔH‡ in kJ/mol (for from_enthalpy_entropy mode)."),
        ("delta_s_double_dagger", "float", "N/A", "Activation entropy ΔS‡ in J/(mol·K) (for from_enthalpy_entropy mode)."),
        # Common parameters:
        ("temperature_k", "float", "298.15", "Temperature in Kelvin."),
        ("transmission_coefficient", "float", "1.0", "Transmission coefficient κ (default 1, typically 0.5–2)."),
        # For eyring_plot mode:
        ("temperatures_k", "str", "", "Temperatures for Eyring plot, comma-separated."),
        ("rate_constants", "str", "", "Rate constants at each temperature, comma-separated."),
        # For compare mode:
        ("arrhenius_ea", "float", "", "Arrhenius activation energy Ea in kJ/mol (for compare mode)."),
        ("arrhenius_a", "float", "", "Arrhenius pre-exponential factor A (for compare mode)."),
        # Units
        ("k_unit_output", "str", "s⁻¹", "Desired output unit for k: 's⁻¹', 'min⁻¹', etc."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated parameters string."),
    ]

    output_sig = [
        ("rate_constant", "float", "Calculated rate constant k."),
        ("k_unit", "str", "Unit of rate constant."),
        ("delta_g_dd", "float", "ΔG‡ at given temperature (kJ/mol)."),
        ("delta_h_dd", "float", "ΔH‡ (kJ/mol), if applicable."),
        ("delta_s_dd", "float", "ΔS‡ (J/(mol·K)), if applicable."),
        ("eyring_equation", "str", "The Eyring equation used."),
        ("analysis", "str", "Detailed analysis including physical interpretation of activation parameters."),
        ("plot_data", "list", "Data points for Eyring plot or temperature dependence curve."),
        ("comparison", "dict", "Comparison data (for compare mode)."),
    ]

    examples         = [
        {
            "code_input": {
                "calculation_mode": 'from_gibbs',
                "delta_g_double_dagger": 75.0,
                "delta_h_double_dagger": 0.0,
                "delta_s_double_dagger": 0.0,
                "temperature_k": 298.15,
                "transmission_coefficient": 1.0,
                "temperatures_k": '',
                "rate_constants": '',
                "arrhenius_ea": 0.0,
                "arrhenius_a": 0.0,
                "k_unit_output": 's^-1'
            },
            "text_input": {
                "input_params": 'from_gibbs 75.0 298.15'
            },
            "output": {
                "rate_constant": 4.57e-06,
                "k_unit": 's^-1',
                "delta_g_dd": 75.0,
                "delta_h_dd": 0,
                "delta_s_dd": 0,
                "eyring_equation": 'k = (kB*T/h) * exp(-dG/RT)',
                "analysis": 'Slow reaction at room temp.',
                "plot_data": [],
                "comparison": 0
            }
        },
        {
            "code_input": {
                "calculation_mode": 'from_enthalpy_entropy',
                "delta_g_double_dagger": 0.0,
                "delta_h_double_dagger": 72.0,
                "delta_s_double_dagger": -10.0,
                "temperature_k": 298.15,
                "transmission_coefficient": 1.0,
                "temperatures_k": '',
                "rate_constants": '',
                "arrhenius_ea": 0.0,
                "arrhenius_a": 0.0,
                "k_unit_output": 's^-1'
            },
            "text_input": {
                "input_params": 'from_enthalpy_entropy 72.0 -10.0 298.15'
            },
            "output": {
                "rate_constant": 0.000382,
                "k_unit": 's^-1',
                "delta_g_dd": 74.98,
                "delta_h_dd": 72.0,
                "delta_s_dd": -10.0,
                "eyring_equation": 'k = kappa*(kB*T/h)*exp(dS/R)*exp(-dH/RT)',
                "analysis": 'Negative dS suggests ordered TS.',
                "plot_data": [],
                "comparison": 0
            }
        },
        {
            "code_input": {
                "calculation_mode": 'eyring_plot',
                "delta_g_double_dagger": 0.0,
                "delta_h_double_dagger": 0.0,
                "delta_s_double_dagger": 0.0,
                "temperature_k": 298.15,
                "transmission_coefficient": 1.0,
                "temperatures_k": '298,308,318,328,338',
                "rate_constants": '2.5e-6,7.2e-6,1.9e-5,4.8e-5,1.15e-4',
                "arrhenius_ea": 0.0,
                "arrhenius_a": 0.0,
                "k_unit_output": 's^-1'
            },
            "text_input": {
                "input_params": 'eyring_plot 298,308,318,328,338 2.5e-6,7.2e-6,1.9e-5,4.8e-5,1.15e-4'
            },
            "output": {
                "rate_constant": 0.0,
                "k_unit": '',
                "delta_g_dd": 74.98,
                "delta_h_dd": 72.5,
                "delta_s_dd": -8.3,
                "eyring_equation": 'ln(k/T) = -dH/(R*T) + (ln(kB/h) + dS/R)',
                "analysis": 'Eyring plot yields activation parameters.',
                "plot_data": [{'t_k': 298, 'ln_k_T': -12.9}],
                "comparison": 0
            }
        }
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.kB = 1.380649e-23   # J/K (Boltzmann constant)
        self.h = 6.62607015e-34   # J·s (Planck constant)
        self.R = 8.314462618      # J/(mol·K)
        self.NA = 6.02214076e23   # mol⁻¹ (Avogadro's number)

    def _calc_k_from_dg(self, dg_kj_mol: float, T: float, kappa: float = 1.0) -> float:
        """k = κ·(kB·T/h)·exp(−ΔG‡/(RT))"""
        dg_j_mol = dg_kj_mol * 1000.0
        exponent = -dg_j_mol / (self.R * T)
        k = kappa * (self.kB * T / self.h) * math.exp(exponent)
        return k

    def _calc_k_from_dh_ds(self, dh_kj_mol: float, ds_j_mol_K: float, T: float, kappa: float = 1.0) -> tuple:
        """Returns (k, ΔG‡_at_T)"""
        dh_j_mol = dh_kj_mol * 1000.0
        dg_at_T = dh_j_mol - T * ds_j_mol_K  # ΔG‡ = ΔH‡ − TΔS‡
        k = self._calc_k_from_dg(dg_at_T / 1000.0, T, kappa)
        return k, dg_at_T / 1000.0

    def _linear_regression(self, x_vals, y_vals):
        n = len(x_vals)
        sx = sum(x_vals); sy = sum(y_vals)
        sxy = sum(xi*yi for xi, yi in zip(x_vals, y_vals))
        sx2 = sum(xi*xi for xi in x_vals)
        denom = n * sx2 - sx*sx
        if abs(denom) < 1e-15:
            return 0.0, sy/n, 1.0
        slope = (n * sxy - sx * sy) / denom
        intercept = (sy - slope * sx) / n
        y_mean = sy / n
        ss_tot = sum((yi-y_mean)**2 for yi in y_vals)
        ss_res = sum((yi-(slope*xi+intercept))**2 for xi, yi in zip(x_vals, y_vals))
        r_sq = 1 - ss_res/ss_tot if ss_tot > 0 else 1.0
        return slope, intercept, max(0, r_sq)

    def _run_base(
        self,
        calculation_mode: str,
        delta_g_double_dagger: float = 0.0,
        delta_h_double_dagger: float = 0.0,
        delta_s_double_dagger: float = 0.0,
        temperature_k: float = 298.15,
        transmission_coefficient: float = 1.0,
        temperatures_k: str = "",
        rate_constants: str = "",
        arrhenius_ea: float = 0.0,
        arrhenius_a: float = 0.0,
        k_unit_output: str = "s⁻¹",
    ) -> dict:
        mode = calculation_mode.lower().strip()
        T = temperature_k

        if T <= 0:
            raise ChemMCPError("Temperature must be positive.")

        if mode == "from_gibbs":
            k_val = self._calc_k_from_dg(delta_g_double_dagger, T, transmission_coefficient)

            eq_str = f"k = κ·(k_B·T/h) · exp(−ΔG‡/RT)"
            interp = (
                f"ΔG‡ = {delta_g_double_dagger} kJ/mol at {T} K\n"
                f"Rate constant k = {k_val:.4e} s⁻¹\n\n"
                f"Interpretation:\n"
                + self._interpret_dg(delta_g_double_dagger, T)
            )
            plot_data = []

        elif mode == "from_enthalpy_entropy":
            k_val, dg_at_T = self._calc_k_from_dh_ds(delta_h_double_dagger, delta_s_double_dagger, T, transmission_coefficient)

            eq_str = f"k = κ·(k_B·T/h) · exp(ΔS‡/R) · exp(−ΔH‡/RT)"
            entropy_interp = self._interpret_ds(delta_s_double_dagger)
            interp = (
                f"ΔH‡ = {delta_h_double_dagger} kJ/mol\n"
                f"ΔS‡ = {delta_s_double_dagger} J/(mol·K)\n"
                f"ΔG‡({T} K) = ΔH‡ − TΔS‡ = {dg_at_T:.2f} kJ/mol\n"
                f"k = {k_val:.4e} s⁻¹\n\n"
                f"{entropy_interp}"
            )
            plot_data = []

        elif mode == "eyring_plot":
            temps = [float(t.strip()) for t in temperatures_k.split(",")]
            ks = [float(k.strip()) for k in rate_constants.split(",")]
            if len(temps) != len(ks):
                raise ChemMCPError("temperatures_k and rate_constants must have same length.")

            # Eyring plot: ln(k/T) vs 1/T
            # ln(k/T) = ln(kB/h) + ΔS‡/R − ΔH‡/(R·T)
            ln_k_over_T = [math.log(k / t) for k, t in zip(ks, temps)]
            inv_T = [1.0 / t for t in temps]
            slope, intercept, r_sq = self._linear_regression(inv_T, ln_k_over_T)

            # slope = −ΔH‡/R → ΔH‡ = −slope·R
            dh_fitted = -slope * self.R / 1000.0  # kJ/mol
            # intercept = ln(kB/h) + ΔS‡/R → ΔS‡ = R·(intercept − ln(kB/h))
            ln_kb_over_h = math.log(self.kB / self.h)
            ds_fitted = self.R * (intercept - ln_kb_over_h)  # J/(mol·K)

            eq_str = "ln(k/T) = −ΔH‡/(R·T) + (ln(kB/h) + ΔS‡/R)"
            k_val = 0.0  # Not a single-k calculation
            interp = (
                f"Eyring Plot Analysis ({len(temps)} points):\n"
                f"Slope = {slope:.4f} → ΔH‡ = {dh_fitted:.2f} kJ/mol\n"
                f"Intercept = {intercept:.4f} → ΔS‡ = {ds_fitted:.2f} J/(mol·K)\n"
                f"R² = {r_sq:.6f}\n\n"
                f"ΔG‡(298.15 K) ≈ {dh_fitted - 298.15 * ds_fitted / 1000:.2f} kJ/mol\n"
                + self._interpret_ds(ds_fitted)
            )
            plot_data = [{"t_k": t, "inv_t": round(1/t, 8), "ln_k_T": round(lkt, 6)} for t, lkt in zip(temps, ln_k_over_T)]

        elif mode == "compare":
            if arrhenius_ea <= 0 or arrhenius_a <= 0:
                raise ChemMCPError("Compare mode requires valid arrhenius_ea (>0) and arrhenius_a (>0).")

            # Arrhenius: k_arr = A·exp(−Ea/RT)
            k_arr = arrhenius_a * math.exp(-arrhenius_ea * 1000.0 / (self.R * T))

            # TST approximation: Ea ≈ ΔH‡ + RT (for gas phase, unimolecular)
            dh_approx = arrhenius_ea - self.R * T / 1000.0
            k_tst, dg_tst = self._calc_k_from_dh_ds(dh_approx, 0.0, T, transmission_coefficient)

            eq_str = "Comparing Arrhenius: k=A·exp(−Ea/RT) vs TST: k=κ(kBT/h)exp(−ΔG‡/RT)"
            k_val = k_arr
            interp = (
                f"At T = {T} K:\n"
                f"Arrhenius: k = {k_arr:.4e} (Ea={arrhenius_ea} kJ/mol, A={arrhenius_a:.4e})\n"
                f"TST approx (ΔH‡≈Ea−RT={dh_approx:.2f} kJ/mol): k = {k_tst:.4e}\n"
                f"Relation: Ea = ΔH‡ + RT (for unimolecular gas-phase reactions)\n"
                f"For solution: Ea ≈ ΔH‡ + RT also holds approximately."
            )
            plot_data = []
            comparison = {
                "arrhenius_k": k_arr,
                "tst_k": k_tst,
                "ea_kj_mol": arrhenius_ea,
                "dh_approx_kj_mol": round(dh_approx, 2),
            }

        else:
            raise ChemMCPError(
                f"Unsupported mode: '{mode}'. Use 'from_gibbs', 'from_enthalpy_entropy', 'eyring_plot', or 'compare'."
            )

        result = {
            "rate_constant": round(k_val, 6) if isinstance(k_val, (int, float)) else k_val,
            "k_unit": k_unit_output,
            "eyring_equation": eq_str,
            "analysis": interp,
            "plot_data": plot_data[:8],
        }

        if mode == "from_gibbs":
            result["delta_g_dd"] = delta_g_double_dagger
        elif mode == "from_enthalpy_entropy":
            result["delta_g_dd"] = round(dg_at_T, 2)
            result["delta_h_dd"] = delta_h_double_dagger
            result["delta_s_dd"] = delta_s_double_dagger
        elif mode == "eyring_plot":
            result["delta_h_dd"] = round(dh_fitted, 2)
            result["delta_s_dd"] = round(ds_fitted, 2)
            result["delta_g_dd"] = round(dh_fitted - 298.15 * ds_fitted / 1000, 2)
        elif mode == "compare":
            result.update(comparison)

        return result

    def _interpret_dg(self, dg_kj: float, T: float) -> str:
        """Provide physical interpretation of ΔG‡."""
        # Reference: typical ΔG‡ values at room temp
        k_ref = self._calc_k_from_dg(dg_kj, T)
        if k_ref > 1e3:
            return "Very fast: diffusion-controlled or barrierless process."
        elif k_ref > 1:
            return "Fast reaction: low barrier, proceeds readily at room temperature."
        elif k_ref > 1e-3:
            return "Moderate rate: observable on seconds-to-minutes timescale."
        elif k_ref > 1e-6:
            return "Slow reaction: minutes to hours timescale."
        elif k_ref > 1e-9:
            return "Very slow: hours to days timescale."
        else:
            return "Extremely slow: effectively stable under these conditions."

    def _interpret_ds(self, ds: float) -> str:
        """Interpret activation entropy."""
        if ds < -50:
            return (
                "Large negative ΔS‡: highly ordered transition state.\n"
                "Suggests associative mechanism (bond formation dominates),\n"
                "bimolecular step, or significant solvation ordering."
            )
        elif ds < -10:
            return (
                "Moderately negative ΔS‡: somewhat ordered transition state.\n"
                "Typical of bimolecular reactions in solution."
            )
        elif ds <= 10:
            return (
                "Near-zero ΔS‡: transition state similar in disorder to reactants.\n"
                "Typical of unimolecular reactions or intramolecular processes."
            )
        elif ds <= 50:
            return (
                "Positive ΔS‡: disordered transition state.\n"
                "Suggests dissociative mechanism (bond breaking dominates)\n"
                "or release of solvent molecules upon reaching TS."
            )
        else:
            return (
                "Large positive ΔS‡: very disordered transition state.\n"
                "Strongly suggests dissociative mechanism or\n"
                "significant increase in degrees of freedom."
            )

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            mode = parts[0]
            kwargs = {"calculation_mode": mode}
            idx = 1

            if mode == "from_gibbs":
                kwargs["delta_g_double_dagger"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["temperature_k"] = float(parts[idx]); idx += 1
            elif mode == "from_enthalpy_entropy":
                kwargs["delta_h_double_dagger"] = float(parts[idx]); idx += 1
                kwargs["delta_s_double_dagger"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["temperature_k"] = float(parts[idx]); idx += 1
            elif mode == "eyring_plot":
                kwargs["temperatures_k"] = parts[idx]; idx += 1
                kwargs["rate_constants"] = parts[idx]; idx += 1
            elif mode == "compare":
                kwargs["arrhenius_ea"] = float(parts[idx]); idx += 1
                kwargs["arrhenius_a"] = float(parts[idx]); idx += 1
                if idx < len(parts):
                    kwargs["temperature_k"] = float(parts[idx]); idx += 1

            if idx < len(parts) and "temperature_k" not in kwargs:
                try:
                    kwargs["temperature_k"] = float(parts[idx]); idx += 1
                except ValueError:
                    pass

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
