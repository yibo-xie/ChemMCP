import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Gas constant
R = 8.314  # J/(mol·K)


@ChemMCPManager.register_tool
class EntropyCalculator(BaseTool):
    """
    计算系统熵变，支持多种过程类型。
    """
    __version__ = "0.1.0"
    name = "EntropyCalculator"
    func_name = "calculate_entropy_change"
    description = "Calculate entropy change for various thermodynamic processes: reaction, phase transition, heating/cooling, and ideal gas mixing."
    implementation_description = "Supports four modes: 'reaction' (ΔS° from standard entropies), 'phase_transition' (ΔS=ΔH_trans/T), 'heating' (Cp·ln(T2/T1)), 'mixing' (ideal: -R·Σxi·ln(xi))."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Entropy", "Physical Chemistry", "Process Analysis"]
    required_envs = []

    code_input_sig = [
        ("mode", "str", "N/A", "Process mode: 'reaction', 'phase_transition', 'heating', or 'mixing'."),
        ("params_str", "str", "N/A", "JSON string of mode-specific parameters. See description for each mode."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: mode params_json. Example: 'heating \"{\\\"cp\\\":29.1,\\\"t1\\\":300,\\\"t2\\\":400}\"'"),
    ]

    output_sig = [
        ("delta_s", "float", "Entropy change ΔS in J/(mol·K)."),
        ("explanation", "str", "Step-by-step calculation explanation."),
    ]

    examples = [
        {
            "code_input": {
                "mode": "phase_transition",
                "params_str": '{"delta_h_trans": 40700.0, "temperature_k": 373.15}',
            },
            "text_input": {
                "input_params": "phase_transition {\"delta_h_trans\":40700.0,\"temperature_k\":373.15}",
            },
            "output": {
                "delta_s": 109.06,
                "explanation": "Phase transition at T=373.15 K: ΔS = ΔH_trans/T = 40700.0/373.15 = 109.06 J/(mol·K)",
            },
        },
        {
            "code_input": {
                "mode": "heating",
                "params_str": '{"cp": 29.1, "t1": 300.0, "t2": 400.0}',
            },
            "text_input": {
                "input_params": "heating {\"cp\":29.1,\"t1\":300,\"t2\":400}",
            },
            "output": {
                "delta_s": 8.4166,
                "explanation": "Heating: ΔS = Cp·ln(T2/T1) = 29.1·ln(400/300) = 8.42 J/(mol·K)",
            },
        },
    ]

    # Standard molar entropies S° (J/(mol·K)) at 298 K
    _standard_entropies = {
        "H2(g)": 130.7, "O2(g)": 205.2, "N2(g)": 191.6,
        "CO(g)": 197.7, "CO2(g)": 213.8, "H2O(g)": 188.8,
        "H2O(l)": 69.9, "NH3(g)": 192.5, "CH4(g)": 186.3,
        "C2H6(g)": 229.6, "C2H4(g)": 219.3, "C2H2(g)": 200.9,
        "HCl(g)": 186.9, "SO2(g)": 248.2, "NO(g)": 210.8,
        "Fe(s)": 27.3, "NaCl(s)": 72.1, "CaCO3(s)": 92.9,
        "C(graphite)": 5.74, "C(diamond)": 2.38,
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, mode: str, params_str: str) -> dict:
        """Core logic: calculate entropy change for different process types."""
        import json

        try:
            params = json.loads(params_str)
        except json.JSONDecodeError:
            raise ChemMCPError(f"Invalid JSON in params_str: {params_str}")

        mode = mode.lower().strip()

        if mode == "reaction":
            return self._calc_reaction_entropy(params)
        elif mode == "phase_transition":
            return self._calc_phase_transition(params)
        elif mode == "heating":
            return self._calc_heating(params)
        elif mode == "mixing":
            return self._calc_mixing(params)
        else:
            raise ChemMCPError(f"Unknown mode '{mode}'. Use 'reaction', 'phase_transition', 'heating', or 'mixing'.")

    def _calc_reaction_entropy(self, params):
        species = params.get("species", {})
        delta_s = 0.0
        details = []
        for sp, coeff in species.items():
            sp_key = sp.strip()
            if sp_key not in self._standard_entropies:
                raise ChemMCPError(f"Unknown substance '{sp_key}' in entropy database.")
            s_val = self._standard_entropies[sp_key]
            contribution = coeff * s_val
            delta_s += contribution
            details.append(f"  {coeff}×{sp_key}: {coeff} × {s_val} = {contribution:.2f}")
        explanation = f"Reaction entropy (ΣS°(products) - ΣS°(reactants)):\n" + "\n".join(details) + f"\nΔS° = {delta_s:.2f} J/(mol·K)"
        return {"delta_s": round(delta_s, 4), "explanation": explanation}

    def _calc_phase_transition(self, params):
        dh = params.get("delta_h_trans")
        t = params.get("temperature_k")
        if dh is None or t is None:
            raise ChemMCPError("Phase transition requires 'delta_h_trans' (J/mol) and 'temperature_k' (K).")
        if t <= 0:
            raise ChemMCPError("Temperature must be positive in Kelvin.")
        delta_s = dh / t
        explanation = f"Phase transition at T={t} K: ΔS = ΔH_trans/T = {dh}/{t} = {delta_s:.4f} J/(mol·K)"
        return {"delta_s": round(delta_s, 4), "explanation": explanation}

    def _calc_heating(self, params):
        cp = params.get("cp")
        t1 = params.get("t1")
        t2 = params.get("t2")
        if any(v is None for v in [cp, t1, t2]):
            raise ChemMCPError("Heating requires 'cp', 't1', and 't2'.")
        if t1 <= 0 or t2 <= 0:
            raise ChemMCPError("Temperatures must be positive in Kelvin.")
        delta_s = cp * math.log(t2 / t1)
        explanation = f"Heating: ΔS = Cp·ln(T2/T1) = {cp}·ln({t2}/{t1}) = {delta_s:.4f} J/(mol·K)"
        return {"delta_s": round(delta_s, 4), "explanation": explanation}

    def _calc_mixing(self, params):
        mole_fractions = params.get("mole_fractions")
        if not mole_fractions or not isinstance(mole_fractions, dict):
            raise ChemMCPError("Mixing requires 'mole_fractions' as a dict of component names to mole fractions.")
        total = sum(mole_fractions.values())
        if abs(total - 1.0) > 1e-6:
            raise ChemMCPError(f"Mole fractions must sum to 1.0, got {total}.")
        delta_s = 0.0
        details = []
        for comp, xi in mole_fractions.items():
            if xi <= 0:
                raise ChemMCPError(f"Mole fraction must be positive, got {xi} for '{comp}'.")
            term = -R * xi * math.log(xi)
            delta_s += term
            details.append(f"  -R·x_{comp}·ln(x_{comp}) = -8.314×{xi}×ln({xi}) = {term:.4f}")
        explanation = f"Ideal gas mixing:\n" + "\n".join(details) + f"\nΔS_mix = {delta_s:.4f} J/(mol·K)"
        return {"delta_s": round(delta_s, 4), "explanation": explanation}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split(None, 1)
            mode = parts[0]
            params_str = parts[1] if len(parts) > 1 else "{}"
            return self._run_base(mode, params_str)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'mode params_json'")
